//! Tadpole contig mode (`Tadpole.process2`): canonical k-mer counting,
//! multi-pass seeding with decreasing depth thresholds, bidirectional
//! greedy extension with ownership, and deterministic longest-first
//! output; bubbles are popped by `super::bubble`.

use super::bubble::pop_bubbles;
use super::{is_junction, read_records, AssembleOptions, AssembleStats};
use crate::libs::asm::refine::{argmax2, number_to_base, second_highest_position};
use crate::libs::asm::table::{base_code, base_defined, Kmer, KmerFnvHasher, RefineTable};
use anyhow::Result;
use pgr::libs::kmer::key;
use pgr::libs::nt::rev_comp;
use std::collections::HashSet;
use std::io::Write;

/// Result codes from `extendToRight` (ShaveObject).
pub(crate) const DEAD_END: i32 = 1;
pub(crate) const LOOP: i32 = 7;
pub(crate) const BAD_OWNER: i32 = 11;
pub(crate) const BAD_SEED: i32 = 12;
pub(crate) const F_BRANCH: i32 = 17;
pub(crate) const B_BRANCH: i32 = 18;
pub(crate) const D_BRANCH: i32 = 19;

/// One assembled contig.
#[derive(Clone)]
pub(crate) struct Contig {
    pub(crate) bases: Vec<u8>,
    pub(crate) id: usize,
    pub(crate) coverage: f32,
    pub(crate) min_cov: usize,
    pub(crate) max_cov: usize,
    pub(crate) left_code: i32,
    pub(crate) right_code: i32,
    pub(crate) left_ratio: f32,
    pub(crate) right_ratio: f32,
    pub(crate) used: bool,
    pub(crate) associate: bool,
    pub(crate) flipped: bool,
    pub(crate) left_edges: Vec<EdgeRef>,
    pub(crate) right_edges: Vec<EdgeRef>,
}

/// Directed edge between two contigs (assemble.Edge).
#[derive(Clone)]
pub(crate) struct Edge {
    pub(crate) origin: usize,
    pub(crate) destination: usize,
    pub(crate) length: usize,
    /// bit 0: source connects on its right; bit 1: destination on its right.
    pub(crate) orientation: u8,
    pub(crate) depth: u32,
    pub(crate) bases: Vec<u8>,
}

impl Edge {
    pub(crate) fn dest_right(&self) -> bool {
        self.orientation & 2 == 2
    }

    fn flip_source(&mut self) {
        self.bases = rev_comp(&self.bases).collect();
        self.orientation ^= 1;
    }

    pub(crate) fn flip_dest(&mut self) {
        self.orientation ^= 2;
    }
}

/// Assembles reads into contigs via the k-mer graph (tadpole contigMode).
///
/// Mirrors `Tadpole.process2(contigMode)`: canonical k-mer counting, then
/// multi-pass seeding with decreasing depth thresholds, bidirectional greedy
/// extension with ownership, and deterministic longest-first output.
pub fn assemble<W: Write>(
    infiles: &[String],
    out: &mut W,
    opts: &AssembleOptions,
) -> Result<AssembleStats> {
    anyhow::ensure!(
        opts.k >= 1,
        "k-mer length must be at least 1, got {}",
        opts.k
    );
    anyhow::ensure!(
        opts.k <= key::Kmer::MAX_K,
        "k-mer length must be at most {} (the k-mer key limit), got {}",
        key::Kmer::MAX_K,
        opts.k
    );
    anyhow::ensure!(
        opts.min_count_seed >= 1,
        "min-count-seed must be at least 1, got {} (0 treats every k-mer as solid and erases error filtering)",
        opts.min_count_seed
    );

    // Read + canonicalize + phred-convert in one pass (one record buffer).
    let reads = read_records(infiles)?;

    // Pass 2: count k-mers from the canonicalized (phred) qualities.
    let table = RefineTable::build_threaded(&reads, opts.k, opts.min_prob, opts.parallel);

    // Pass 3: multi-pass seeding and contig building (BuildThread.run).
    let mut claimed: HashSet<Kmer, std::hash::BuildHasherDefault<KmerFnvHasher>> =
        HashSet::default();
    let mut contigs: Vec<Contig> = Vec::new();
    let mut id_counter = 0usize;
    for i in (1..opts.contig_passes).rev() {
        let threshold = pass_threshold(opts, i);
        scan_table(
            &table,
            threshold,
            opts,
            &mut claimed,
            &mut contigs,
            &mut id_counter,
        );
    }
    scan_table(
        &table,
        opts.min_count_seed,
        opts,
        &mut claimed,
        &mut contigs,
        &mut id_counter,
    );

    // Contig graph + bubble popping (Tadpole.processContigs/popBubbles);
    // with --no-bubbles the pre-pop contigs are kept and only sorted and
    // renumbered.
    if opts.pop_bubbles {
        process_contigs(&mut contigs, &table, opts);
        pop_bubbles(&mut contigs, opts);
    } else {
        finalize_contigs(&mut contigs);
    }

    let mut stats = AssembleStats {
        reads_in: reads.len() as u64,
        ..AssembleStats::default()
    };
    let min_contig_len = opts.resolved_min_contig_len();
    for c in &contigs {
        if c.bases.len() >= min_contig_len {
            write_contig(out, c)?;
            stats.contigs_built += 1;
            stats.bases_built += c.bases.len() as u64;
            stats.longest_contig = stats.longest_contig.max(c.bases.len());
        }
    }
    Ok(stats)
}

/// Seeding threshold for pass `i` (Java `minCountSeedCurrent` formula).
pub(crate) fn pass_threshold(opts: &AssembleOptions, i: usize) -> usize {
    let t = (opts.min_count_seed as f64 * opts.contig_pass_mult.powi(i as i32) * 0.92 - 0.25)
        .floor() as i64;
    (opts.min_count_seed as i64 + i as i64)
        .max(t)
        .min(i32::MAX as i64) as usize
}

/// One seeding scan over all table k-mers (BuildThread.processNextTable).
#[allow(clippy::too_many_arguments)]
pub(crate) fn scan_table(
    table: &RefineTable,
    threshold: usize,
    opts: &AssembleOptions,
    claimed: &mut HashSet<Kmer, std::hash::BuildHasherDefault<KmerFnvHasher>>,
    contigs: &mut Vec<Contig>,
    id_counter: &mut usize,
) {
    // Deterministic scan order by canonical k-mer sequence (the BBTools
    // hash-table cell order is memory-dependent and not portable). The
    // sorted snapshot is cached in the table, so all 16 seeding passes
    // iterate it linearly instead of re-sorting the HashMap each pass.
    let entries = table.sorted_entries();
    for (kmer, count) in entries.iter() {
        if *count < threshold as u32 {
            continue;
        }
        if claimed.contains(kmer) {
            continue;
        }
        claimed.insert(*kmer);
        if let Some(c) = make_contig(kmer, table, opts, claimed) {
            let mut c = c;
            c.id = *id_counter;
            *id_counter += 1;
            contigs.push(c);
        }
    }
}

/// Builds one contig from a claimed seed (Tadpole2.makeContig).
fn make_contig(
    seed: &Kmer,
    table: &RefineTable,
    opts: &AssembleOptions,
    claimed: &mut HashSet<Kmer, std::hash::BuildHasherDefault<KmerFnvHasher>>,
) -> Option<Contig> {
    let k = opts.k;
    // `base_at(0)` is the 3' end (last base pushed); rebuild 5'->3'.
    let mut bb: Vec<u8> = (0..k)
        .map(|i| number_to_base(seed.base_at(k - 1 - i)))
        .collect();
    debug_assert_eq!(bb.len(), k);

    let (right_status, mut right_ratio) = extend_to_right(&mut bb, table, opts, claimed);
    match right_status {
        DEAD_END | LOOP => {}
        BAD_SEED => return None,
        _ => {
            if bb.len() == k {
                // A branch or ownership failure at the seed rejects the contig.
                return None;
            }
            match right_status {
                BAD_OWNER => return None,
                F_BRANCH | D_BRANCH => {
                    right_ratio = calc_ratio(&right_counts_of(bb.as_slice(), table, opts))
                }
                B_BRANCH => right_ratio = calc_ratio(&left_counts_of(bb.as_slice(), table, opts)),
                _ => return None,
            }
        }
    }

    // Extend the left end by reverse-complementing and extending right.
    let mut rc: Vec<u8> = rev_comp(&bb).collect();
    let (left_status, mut left_ratio) = extend_to_right(&mut rc, table, opts, claimed);
    match left_status {
        DEAD_END | LOOP => {}
        BAD_SEED => return None,
        _ => match left_status {
            BAD_OWNER => return None,
            F_BRANCH | D_BRANCH => {
                left_ratio = calc_ratio(&right_counts_of(rc.as_slice(), table, opts))
            }
            B_BRANCH => left_ratio = calc_ratio(&left_counts_of(rc.as_slice(), table, opts)),
            _ => return None,
        },
    }
    bb = rev_comp(&rc).collect();

    // With bubble popping enabled (the default), BBTools keeps every contig
    // of at least k+minExtension internally; the minContigLen filter applies
    // only at output time (short contigs still anchor graph edges).
    if bb.len() >= k + opts.min_extension {
        let (coverage, min_cov, max_cov) = calc_coverage(&bb, table, k);
        if coverage < opts.min_coverage {
            return None;
        }
        // Canonical orientation (Contig.canonical + rcomp).
        let (bases, left_code, right_code, left_ratio, right_ratio) = if canonical(&bb) {
            (bb, left_status, right_status, left_ratio, right_ratio)
        } else {
            (
                rev_comp(&bb).collect(),
                right_status,
                left_status,
                right_ratio,
                left_ratio,
            )
        };
        Some(Contig {
            bases,
            id: 0,
            coverage,
            min_cov,
            max_cov,
            left_code,
            right_code,
            left_ratio,
            right_ratio,
            used: false,
            associate: false,
            flipped: false,
            left_edges: Vec::new(),
            right_edges: Vec::new(),
        })
    } else {
        None
    }
}

/// Counts of the four right/left extensions of a k-mer at `bb`'s 3'/5' end.
fn right_counts_of(bb: &[u8], table: &RefineTable, opts: &AssembleOptions) -> [u32; 4] {
    let k = opts.k;
    let mut kmer = Kmer::new(k);
    for &b in &bb[bb.len() - k..] {
        kmer.push_right(base_code(b));
    }
    table.fill_right_counts(&kmer)
}

fn left_counts_of(bb: &[u8], table: &RefineTable, opts: &AssembleOptions) -> [u32; 4] {
    let k = opts.k;
    let mut kmer = Kmer::new(k);
    for &b in &bb[bb.len() - k..] {
        kmer.push_right(base_code(b));
    }
    table.fill_left_counts(&kmer)
}

/// `extendToRight` (contig mode): bidirectional-aware greedy extension.
///
/// Returns the exit status and, for branch exits, the branch ratio.
fn extend_to_right(
    bb: &mut Vec<u8>,
    table: &RefineTable,
    opts: &AssembleOptions,
    claimed: &mut HashSet<Kmer, std::hash::BuildHasherDefault<KmerFnvHasher>>,
) -> (i32, f32) {
    let k = opts.k;
    if bb.len() < k {
        return (BAD_SEED, 0.0);
    }
    // Rightmost k-mer of the current sequence.
    let mut kmer = Kmer::new(k);
    let mut len = 0usize;
    for &b in &bb[bb.len() - k..] {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
    }
    if len < k {
        return (BAD_SEED, 0.0);
    }
    if table.get_count(&kmer) < opts.min_count_seed as u32 {
        return (BAD_SEED, 0.0);
    }

    let mut left = table.fill_left_counts(&kmer);
    let mut left_max_pos = argmax2(&left, &mut 0);
    let mut left_max = left[left_max_pos];
    let left_second_pos = second_highest_position(&left);
    let left_second = left[left_second_pos];

    let mut right = table.fill_right_counts(&kmer);
    let mut right_max_pos = argmax2(&right, &mut 0);
    let mut right_max = right[right_max_pos];
    let right_second_pos = second_highest_position(&right);
    let right_second = right[right_second_pos];

    if right_max < opts.min_count_extend as u32 {
        return (DEAD_END, 0.0);
    }
    if is_junction(
        right_max,
        right_second,
        opts.branch_mult1,
        opts.branch_mult2,
        opts.branch_lower_const,
        opts.min_count_extend,
    ) {
        let d = is_junction(
            left_max,
            left_second,
            opts.branch_mult1,
            opts.branch_mult2,
            opts.branch_lower_const,
            opts.min_count_extend,
        );
        return if d {
            (D_BRANCH, calc_ratio(&right))
        } else {
            (F_BRANCH, calc_ratio(&right))
        };
    }
    if is_junction(
        left_max,
        left_second,
        opts.branch_mult1,
        opts.branch_mult2,
        opts.branch_lower_const,
        opts.min_count_extend,
    ) {
        return (B_BRANCH, calc_ratio(&left));
    }

    // The seed was claimed by the caller (single-threaded ownership).
    let max_len = 1_000_000_000usize;
    while bb.len() < max_len {
        let b = right_max_pos as u8;
        let evicted = kmer.base_at(k - 1);
        kmer.push_right(b);

        left = table.fill_left_counts(&kmer);
        left_max_pos = argmax2(&left, &mut 0);
        left_max = left[left_max_pos];
        let left_second_pos = second_highest_position(&left);
        let left_second = left[left_second_pos];

        right = table.fill_right_counts(&kmer);
        right_max_pos = argmax2(&right, &mut 0);
        right_max = right[right_max_pos];
        let right_second_pos = second_highest_position(&right);
        let right_second = right[right_second_pos];

        let fbranch = is_junction(
            right_max,
            right_second,
            opts.branch_mult1,
            opts.branch_mult2,
            opts.branch_lower_const,
            opts.min_count_extend,
        );
        let bbranch = is_junction(
            left_max,
            left_second,
            opts.branch_mult1,
            opts.branch_mult2,
            opts.branch_lower_const,
            opts.min_count_extend,
        );
        let hbranch = left_max_pos != evicted as usize && opts.branch_mult1 > 0.0;
        if bbranch || hbranch {
            let ratio = if fbranch {
                calc_ratio(&right)
            } else {
                calc_ratio(&left)
            };
            return if fbranch {
                (D_BRANCH, ratio)
            } else {
                (B_BRANCH, ratio)
            };
        }

        bb.push(number_to_base(b));

        // Loop detection / ownership claim (single-thread id=0).
        let canonical = kmer.canonical();
        if claimed.contains(&canonical) {
            return if fbranch {
                (F_BRANCH, calc_ratio(&right))
            } else {
                (LOOP, 0.0)
            };
        }
        claimed.insert(canonical);

        if fbranch {
            return (F_BRANCH, calc_ratio(&right));
        }
        if right_max < opts.min_count_extend as u32 {
            return (DEAD_END, 0.0);
        }
    }
    (BAD_OWNER, 0.0)
}

/// `KmerTableSet.calcCoverage`: mean/min/max canonical k-mer counts.
pub(crate) fn calc_coverage(bases: &[u8], table: &RefineTable, k: usize) -> (f32, usize, usize) {
    if bases.len() < k {
        return (0.0, 0, 0);
    }
    let mut kmer = Kmer::new(k);
    let mut len = 0usize;
    let mut sum = 0u64;
    let mut max = 0usize;
    let mut min = usize::MAX;
    let mut kmers = 0usize;
    for &b in bases {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
        if len >= k {
            let count = table.get_count(&kmer) as usize;
            sum += count as u64;
            max = max.max(count);
            min = min.min(count);
            kmers += 1;
        }
    }
    if sum == 0 {
        (0.0, 0, 0)
    } else {
        (sum as f32 / kmers as f32, min, max)
    }
}

/// Mean/min/max coverage from per-k-mer canonical counts already collected
/// in output sequence order (unitig walk variant of `calc_coverage`).
pub(crate) fn cov_from_counts(counts: &[u32]) -> (f32, usize, usize) {
    if counts.is_empty() {
        return (0.0, 0, 0);
    }
    let mut sum = 0u64;
    let mut min = usize::MAX;
    let mut max = 0usize;
    for &c in counts {
        sum += c as u64;
        min = min.min(c as usize);
        max = max.max(c as usize);
    }
    (sum as f32 / counts.len() as f32, min, max)
}

/// `Contig.calcScalarsFast`: gc fraction plus dimer-based hh/caga.
pub(crate) fn calc_scalars(bases: &[u8]) -> (f32, f32, f32) {
    if bases.len() < 2 {
        return (0.0, 0.0, 0.0);
    }
    let mut counts = [0u64; 16];
    let mut prev_bad = 8u8; // "N" so the first dimer is skipped
    let mut prev_val = 0u8;
    let mut at_sum = 0u64;
    let mut gc_sum = 0u64;
    for &b in bases {
        let gcbit = b >> 1;
        at_sum += (!gcbit & 1) as u64;
        gc_sum += (gcbit & !(b >> 3) & 1) as u64;
        let mut val = (b & 6) >> 1;
        val ^= (val & 2) >> 1;
        let bad = b & 8;
        if (prev_bad | bad) == 0 {
            counts[((prev_val << 2) | val) as usize] += 1;
        }
        prev_val = val;
        prev_bad = bad;
    }
    let aa = counts[0b0000];
    let tt = counts[0b1111];
    let at = counts[0b0011];
    let ta = counts[0b1100];
    let cc = counts[0b0101];
    let gg = counts[0b1010];
    let cg = counts[0b0110];
    let gc = counts[0b1001];
    let ac = counts[0b0001];
    let tg = counts[0b1110];
    let ag = counts[0b0010];
    let ct = counts[0b0111];
    let tc = counts[0b1101];
    let ga = counts[0b1000];
    let gt = counts[0b1011];
    let ca = counts[0b0100];
    let hh = (aa + cc + gg + tt) as f32 / (aa + tt + at + ta + cc + gg + cg + gc).max(1) as f32;
    let caga = 0.5
        * (1.0
            + (ca as i64 + tg as i64 - ga as i64 - tc as i64) as f32
                / (ac + ag + ca + ga + tc + tg + ct + gt).max(1) as f32);
    let gc_frac = gc_sum as f32 / (at_sum + gc_sum).max(1) as f32;
    (gc_frac, hh, caga)
}

/// A contig is canonical iff its sequence <= its reverse complement.
pub(crate) fn canonical(bases: &[u8]) -> bool {
    let n = bases.len();
    for i in 0..n {
        let a = bases[i];
        let b = complement(bases[n - 1 - i]);
        if a < b {
            return true;
        }
        if b < a {
            return false;
        }
    }
    true
}

fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b,
    }
}

/// `calcRatio`: highest / second-highest count, 99 when no second branch.
fn calc_ratio(counts: &[u32; 4]) -> f32 {
    let mut a = 0u32;
    let mut b = 0u32;
    for &x in counts {
        if x > a {
            b = a;
            a = x;
        } else if x > b {
            b = x;
        }
    }
    if b < 1 {
        99.0
    } else {
        a as f32 / b as f32
    }
}

/// `ContigLengthComparator` (descending): length, coverage, sequence, id.
pub(crate) fn contig_cmp(a: &Contig, b: &Contig) -> std::cmp::Ordering {
    match a.bases.len().cmp(&b.bases.len()).reverse() {
        std::cmp::Ordering::Equal => {}
        x => return x,
    }
    if a.coverage != b.coverage {
        return a.coverage.partial_cmp(&b.coverage).unwrap().reverse();
    }
    match a.bases.cmp(&b.bases).reverse() {
        std::cmp::Ordering::Equal => {}
        x => return x,
    }
    a.id.cmp(&b.id).reverse()
}

/// Writes one contig in FASTA (attributes as header comment, 70-column wrap).
fn write_contig<W: Write>(w: &mut W, c: &Contig) -> Result<()> {
    let (gc, hh, caga) = calc_scalars(&c.bases);
    writeln!(
        w,
        ">contig_{} len={},cov={},gc={},min={},max={},hh={},caga={}",
        c.id,
        c.bases.len(),
        fmt_fixed(c.coverage as f64, 1),
        fmt_fixed(gc as f64, 3),
        c.min_cov,
        c.max_cov,
        fmt_fixed(hh as f64, 3),
        fmt_fixed(caga as f64, 3),
    )?;
    for chunk in c.bases.chunks(70) {
        w.write_all(chunk)?;
        w.write_all(b"\n")?;
    }
    Ok(())
}

/// `ByteBuilder.append(double, decimals)`: half-up fixed-point formatting.
pub(crate) fn fmt_fixed(x: f64, decimals: usize) -> String {
    if x == x.trunc() {
        return format!("{}", x as i64);
    }
    if decimals < 1 {
        return format!("{}", (x + 0.5) as i64);
    }
    let neg = x < 0.0;
    let x = x.abs();
    let inv = 10f64.powi(-(decimals as i32));
    let x = x + 0.5 * inv;
    let upper = x as i64;
    let lower = ((x - upper as f64) * 10f64.powi(decimals as i32)) as i64;
    format!(
        "{}{}.{:0width$}",
        if neg { "-" } else { "" },
        upper,
        lower,
        width = decimals
    )
}

/*--------------------------------------------------------------------*/
/*  Contig graph and bubble popping (Tadpole.processContigs)          */
/*--------------------------------------------------------------------*/

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

pub(crate) type EdgeRef = Rc<RefCell<Edge>>;

impl Contig {
    pub(crate) fn left_kmer(&self, k: usize) -> Kmer {
        let mut kmer = Kmer::new(k);
        for &b in &self.bases[..k] {
            kmer.push_right(base_code(b));
        }
        kmer
    }

    pub(crate) fn right_kmer(&self, k: usize) -> Kmer {
        let mut kmer = Kmer::new(k);
        let n = self.bases.len();
        for &b in &self.bases[n - k..] {
            kmer.push_right(base_code(b));
        }
        kmer
    }

    pub(crate) fn left_forward_branch(&self) -> bool {
        self.left_code == F_BRANCH
    }

    pub(crate) fn right_forward_branch(&self) -> bool {
        self.right_code == F_BRANCH
    }

    fn add_left_edge(&mut self, e: EdgeRef) {
        let (dest, orient, depth, len) = {
            let eb = e.borrow();
            (eb.destination, eb.orientation, eb.depth, eb.length)
        };
        if let Some(old) = self.get_left_edge(dest, Some(orient)) {
            let mut ob = old.borrow_mut();
            if depth >= ob.depth && (ob.depth == 1 || ob.length == len) {
                ob.bases = e.borrow().bases.clone();
                ob.length = len;
                ob.depth += depth;
                return;
            }
        }
        self.left_edges.push(e);
    }

    fn add_right_edge(&mut self, e: EdgeRef) {
        let (dest, orient, depth, len) = {
            let eb = e.borrow();
            (eb.destination, eb.orientation, eb.depth, eb.length)
        };
        if let Some(old) = self.get_right_edge(dest, Some(orient)) {
            let mut ob = old.borrow_mut();
            if depth >= ob.depth && (ob.depth == 1 || ob.length == len) {
                ob.bases = e.borrow().bases.clone();
                ob.length = len;
                ob.depth += depth;
                return;
            }
        }
        self.right_edges.push(e);
    }

    fn get_left_edge(&self, dest: usize, orientation: Option<u8>) -> Option<EdgeRef> {
        self.left_edges
            .iter()
            .find(|e| {
                let e = e.borrow();
                e.destination == dest
                    && (orientation.is_none() || orientation == Some(e.orientation))
            })
            .cloned()
    }

    pub(crate) fn get_right_edge(&self, dest: usize, orientation: Option<u8>) -> Option<EdgeRef> {
        self.right_edges
            .iter()
            .find(|e| {
                let e = e.borrow();
                e.destination == dest
                    && (orientation.is_none() || orientation == Some(e.orientation))
            })
            .cloned()
    }

    fn remove_edges_to(&mut self, dest: usize) {
        self.left_edges.retain(|e| e.borrow().destination != dest);
        self.right_edges.retain(|e| e.borrow().destination != dest);
    }

    pub(crate) fn flip(&mut self, inbound: Option<&[EdgeRef]>) {
        self.flipped = !self.flipped;
        self.bases = rev_comp(&self.bases).collect();
        std::mem::swap(&mut self.left_code, &mut self.right_code);
        std::mem::swap(&mut self.left_ratio, &mut self.right_ratio);
        std::mem::swap(&mut self.left_edges, &mut self.right_edges);
        for e in &self.left_edges {
            e.borrow_mut().flip_source();
        }
        for e in &self.right_edges {
            e.borrow_mut().flip_source();
        }
        if let Some(inbound) = inbound {
            for e in inbound {
                e.borrow_mut().flip_dest();
            }
        }
    }

    pub(crate) fn renumber(&mut self, new_id: usize, inbound: Option<&[EdgeRef]>) {
        if self.id == new_id {
            return;
        }
        for e in &self.left_edges {
            e.borrow_mut().origin = new_id;
        }
        for e in &self.right_edges {
            e.borrow_mut().origin = new_id;
        }
        if let Some(inbound) = inbound {
            for e in inbound {
                e.borrow_mut().destination = new_id;
            }
        }
        self.id = new_id;
    }
}

/// Clears a contig's edges and detaches them from live sources
/// (Contig.removeAllEdges); `inbound` is the dest-map entry for `id`.
pub(crate) fn remove_all_edges(id: usize, inbound: Option<&[EdgeRef]>, contigs: &mut [Contig]) {
    contigs[id].left_edges.clear();
    contigs[id].right_edges.clear();
    if let Some(inbound) = inbound {
        for e in inbound {
            let (dest, origin) = {
                let eb = e.borrow();
                (eb.destination, eb.origin)
            };
            if dest == id && origin != id {
                let source = &mut contigs[origin];
                if !source.used && !source.associate {
                    source.remove_edges_to(id);
                }
            }
        }
    }
}

pub(crate) fn set_used(id: usize, inbound: Option<&[EdgeRef]>, contigs: &mut [Contig]) {
    contigs[id].used = true;
    remove_all_edges(id, inbound, contigs);
}

pub(crate) fn set_associate(id: usize, inbound: Option<&[EdgeRef]>, contigs: &mut [Contig]) {
    contigs[id].associate = true;
    remove_all_edges(id, inbound, contigs);
}

/// Builds the contig end-kmer ownership map and edges
/// (Tadpole.initializeContigs + ProcessContigThread).
fn process_contigs(contigs: &mut [Contig], table: &RefineTable, opts: &AssembleOptions) {
    let k = opts.k;
    let mut end_claims: HashMap<Kmer, usize> = HashMap::new();
    for (i, c) in contigs.iter().enumerate() {
        end_claims.entry(c.left_kmer(k).canonical()).or_insert(i);
        end_claims.entry(c.right_kmer(k).canonical()).or_insert(i);
    }
    for i in 0..contigs.len() {
        process_contig_left(i, contigs, table, opts, &end_claims);
        process_contig_right(i, contigs, table, opts, &end_claims);
    }
}

fn process_contig_left(
    c_id: usize,
    contigs: &mut [Contig],
    table: &RefineTable,
    opts: &AssembleOptions,
    end_claims: &HashMap<Kmer, usize>,
) {
    if contigs[c_id].left_code == DEAD_END {
        return;
    }
    let k = opts.k;
    let kmer0 = contigs[c_id].left_kmer(k);
    let left = table.fill_left_counts(&kmer0);
    let left_max_pos = argmax2(&left, &mut 0);
    let left_max = left[left_max_pos];
    let mut edges_to_add: Vec<EdgeRef> = Vec::new();
    for x in 0..4u8 {
        let count = left[x as usize];
        if count > 0
            && is_junction(
                left_max,
                count,
                opts.branch_mult1,
                opts.branch_mult2,
                opts.branch_lower_const,
                opts.min_count_extend,
            )
        {
            let mut kmer = kmer0;
            kmer.push_left(x);
            // Tadpole1 (k <= 31) walks the left edge in reverse-complement
            // space (`processContigLeft` swaps kmer/rkmer into `exploreRight`);
            // Tadpole2 (k > 31) walks it in forward space.
            if opts.k <= 31 {
                kmer = kmer.rc();
            }
            let mut bb = vec![number_to_base(x)];
            let (target, last_length, last_orientation) =
                explore_right(&kmer, table, opts, end_claims, contigs, &mut bb);
            if let Some(target) = target {
                edges_to_add.push(Rc::new(RefCell::new(Edge {
                    origin: c_id,
                    destination: target,
                    length: last_length,
                    orientation: last_orientation,
                    depth: count,
                    bases: bb,
                })));
            }
        }
    }
    for e in edges_to_add {
        contigs[c_id].add_left_edge(e);
    }
}

fn process_contig_right(
    c_id: usize,
    contigs: &mut [Contig],
    table: &RefineTable,
    opts: &AssembleOptions,
    end_claims: &HashMap<Kmer, usize>,
) {
    if contigs[c_id].right_code == DEAD_END {
        return;
    }
    let k = opts.k;
    let kmer0 = contigs[c_id].right_kmer(k);
    let right = table.fill_right_counts(&kmer0);
    let right_max_pos = argmax2(&right, &mut 0);
    let right_max = right[right_max_pos];
    let mut edges_to_add: Vec<EdgeRef> = Vec::new();
    for x in 0..4u8 {
        let count = right[x as usize];
        if count > 0
            && is_junction(
                right_max,
                count,
                opts.branch_mult1,
                opts.branch_mult2,
                opts.branch_lower_const,
                opts.min_count_extend,
            )
        {
            let mut kmer = kmer0;
            kmer.push_right(x);
            let mut bb = vec![number_to_base(x)];
            let (target, last_length, mut last_orientation) =
                explore_right(&kmer, table, opts, end_claims, contigs, &mut bb);
            if let Some(target) = target {
                last_orientation |= 1;
                edges_to_add.push(Rc::new(RefCell::new(Edge {
                    origin: c_id,
                    destination: target,
                    length: last_length,
                    orientation: last_orientation,
                    depth: count,
                    bases: bb,
                })));
            }
        }
    }
    for e in edges_to_add {
        contigs[c_id].add_right_edge(e);
    }
}

/// `ProcessContigThread.exploreRight`: walks from an end k-mer to the next
/// contig end; returns (destination contig, path length, destination-side
/// orientation bit).
fn explore_right(
    kmer0: &Kmer,
    table: &RefineTable,
    opts: &AssembleOptions,
    end_claims: &HashMap<Kmer, usize>,
    contigs: &[Contig],
    bb: &mut Vec<u8>,
) -> (Option<usize>, usize, u8) {
    let k = opts.k;
    let mut kmer = *kmer0;
    let mut length = 1usize;
    let mut owner: Option<usize> = None;
    while length < 500 {
        owner = end_claims.get(&kmer.canonical()).copied();
        if owner.is_some() {
            break;
        }
        let left = table.fill_left_counts(&kmer);
        let left_max_pos = argmax2(&left, &mut 0);
        let left_max = left[left_max_pos];
        let left_second_pos = second_highest_position(&left);
        let left_second = left[left_second_pos];
        if is_junction(
            left_max,
            left_second,
            opts.branch_mult1,
            opts.branch_mult2,
            opts.branch_lower_const,
            opts.min_count_extend,
        ) {
            return (None, length, 0);
        }
        let right = table.fill_right_counts(&kmer);
        let right_max_pos = argmax2(&right, &mut 0);
        let right_max = right[right_max_pos];
        let right_second_pos = second_highest_position(&right);
        let right_second = right[right_second_pos];
        if right_max < opts.min_count_extend as u32 {
            return (None, length, 0);
        }
        if is_junction(
            right_max,
            right_second,
            opts.branch_mult1,
            opts.branch_mult2,
            opts.branch_lower_const,
            opts.min_count_extend,
        ) {
            return (None, length, 0);
        }
        bb.push(number_to_base(right_max_pos as u8));
        kmer.push_right(right_max_pos as u8);
        length += 1;
    }
    if let Some(owner) = owner {
        // Orientation: 0 if the destination's left k-mer matches, 2 if its
        // right k-mer matches (canonical comparison, like Java Kmer.equals).
        let dest = &contigs[owner];
        let mut temp = dest.left_kmer(k);
        let orientation = if kmer_eq(&temp, &kmer) {
            0
        } else {
            temp = dest.right_kmer(k);
            if kmer_eq(&temp, &kmer) {
                2
            } else {
                debug_assert!(false, "exploreRight destination mismatch");
                return (None, length, 0);
            }
        };
        (Some(owner), length, orientation)
    } else {
        (None, length, 0)
    }
}

pub(crate) fn kmer_eq(a: &Kmer, b: &Kmer) -> bool {
    a.canonical().cmp_bases(&b.canonical()) == std::cmp::Ordering::Equal
}

/// Deterministic longest-first sort and renumbering for the no-bubbles path
/// (bubble popping performs the same step while also renumbering edges).
fn finalize_contigs(contigs: &mut [Contig]) {
    contigs.sort_by(contig_cmp);
    for (new_id, c) in contigs.iter_mut().enumerate() {
        c.renumber(new_id, None);
    }
}

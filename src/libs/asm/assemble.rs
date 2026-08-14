//! Tadpole-compatible contig assembly (contigMode).

use crate::libs::asm::dfa::VertexStates;
use crate::libs::asm::refine::{
    argmax2, base_code, base_defined, number_to_base, second_highest_position, Kmer, KmerFnvHasher,
    RefineTable,
};
use anyhow::Result;
use pgr::libs::fmt::seq::{SeqReader, SeqRecord};
use pgr::libs::fq::qual::{from_phred, to_phred};
use pgr::libs::kmer::key;
use pgr::libs::nt::rev_comp;
use std::collections::HashSet;
use std::io::Write;

/// Result codes from `extendToRight` (ShaveObject).
const DEAD_END: i32 = 1;
const LOOP: i32 = 7;
const BAD_OWNER: i32 = 11;
const BAD_SEED: i32 = 12;
const F_BRANCH: i32 = 17;
const B_BRANCH: i32 = 18;
const D_BRANCH: i32 = 19;

/// Assembly options with tadpole.sh defaults.
#[derive(Debug, Clone)]
pub struct AssembleOptions {
    /// K-mer length.
    pub k: usize,
    /// Ignore k-mers below this error-free probability.
    pub min_prob: f32,
    /// Minimum k-mer depth to seed a contig.
    pub min_count_seed: usize,
    /// Minimum k-mer depth to continue an extension.
    pub min_count_extend: usize,
    /// Minimum added bases past the seed for a contig to be kept.
    pub min_extension: usize,
    /// Minimum output length; `0` keeps everything (the `asm unitig`
    /// default mirrors bcalm's lossless vertex decomposition — no length
    /// cleanup unless requested; `asm contig` passes tadpole's
    /// `mincontiglen` default explicitly).
    pub min_contig_len: usize,
    /// Minimum k-mer coverage for a contig.
    pub min_coverage: f32,
    /// Branch ratio at high depth (branchmult1).
    pub branch_mult1: f32,
    /// Branch ratio at low depth (branchmult2).
    pub branch_mult2: f32,
    /// Second-highest depth considered "low" (branchlower).
    pub branch_lower_const: usize,
    /// Number of seeding passes (contigPasses).
    pub contig_passes: usize,
    /// Seeding pass multiplier (contigPassMult).
    pub contig_pass_mult: f64,
    /// Merge parallel paths in the contig graph (Tadpole popbubbles).
    pub pop_bubbles: bool,
    /// Append `L:` links to unitig FASTA headers (BCALM format).
    pub emit_links: bool,
    /// Emit a GFA graph instead of FASTA.
    pub emit_gfa: bool,
    /// Emit every k-mer abundance in the FASTA header (`ab:Z:`, BCALM
    /// `-all-abundance-counts`).
    pub all_abundance_counts: bool,
    /// Experimental: classify vertices once (DFA state) and walk unitigs
    /// from the state table instead of re-scanning extension buckets.
    pub use_dfa: bool,
    /// Experimental: use pgr's FastK-style super-mer two-stage counter
    /// instead of the direct emission + sort path (no quality gating).
    pub use_supermer: bool,
    /// Minimizer length for `--supermer` (None = pgr default).
    pub supermer_m: Option<usize>,
    /// Worker threads for the whole k-mer pipeline (counting + DFA
    /// classification); `0` uses the rayon global pool (all cores). The
    /// walk stays deterministic single-threaded.
    pub parallel: usize,
}

impl Default for AssembleOptions {
    fn default() -> Self {
        Self {
            k: 31,
            min_prob: 0.5,
            min_count_seed: 3,
            min_count_extend: 2,
            min_extension: 2,
            min_contig_len: 0,
            min_coverage: 1.0,
            branch_mult1: 20.0,
            branch_mult2: 3.0,
            branch_lower_const: 3,
            contig_passes: 16,
            contig_pass_mult: 1.7,
            pop_bubbles: true,
            emit_links: false,
            emit_gfa: false,
            all_abundance_counts: false,
            use_dfa: false,
            use_supermer: false,
            supermer_m: None,
            parallel: 0,
        }
    }
}

impl AssembleOptions {
    fn resolved_min_contig_len(&self) -> usize {
        self.min_contig_len
    }
}

/// `Tadpole.isJunction(max, second)`: depth-ratio branch detection.
fn is_junction(max: u32, second: u32, opts: &AssembleOptions) -> bool {
    if second < 1
        || (second as f32) * opts.branch_mult1 < max as f32
        || (second <= opts.branch_lower_const as u32
            && (max as f32)
                >= (opts.min_count_extend as f32).max(second as f32 * opts.branch_mult2))
    {
        return false;
    }
    true
}

/// One assembled contig.
#[derive(Clone)]
struct Contig {
    bases: Vec<u8>,
    id: usize,
    coverage: f32,
    min_cov: usize,
    max_cov: usize,
    left_code: i32,
    right_code: i32,
    left_ratio: f32,
    right_ratio: f32,
    used: bool,
    associate: bool,
    flipped: bool,
    left_edges: Vec<EdgeRef>,
    right_edges: Vec<EdgeRef>,
}

/// Directed edge between two contigs (assemble.Edge).
#[derive(Clone)]
struct Edge {
    origin: usize,
    destination: usize,
    length: usize,
    /// bit 0: source connects on its right; bit 1: destination on its right.
    orientation: u8,
    depth: u32,
    bases: Vec<u8>,
}

impl Edge {
    fn dest_right(&self) -> bool {
        self.orientation & 2 == 2
    }

    fn flip_source(&mut self) {
        self.bases = rev_comp(&self.bases).collect();
        self.orientation ^= 1;
    }

    fn flip_dest(&mut self) {
        self.orientation ^= 2;
    }
}

/// Assembly statistics.
#[derive(Debug, Default, Clone)]
pub struct AssembleStats {
    pub reads_in: u64,
    pub contigs_built: u64,
    pub bases_built: u64,
    pub longest_contig: usize,
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

/// Reads all records from any number of files sequentially, canonicalizing
/// qualities like BBTools (shared by the contig and unitig modes).
/// Reads and converts records to `(sequence, phred quality)` pairs in one
/// streaming pass (one `SeqRecord` buffer alive at a time). The previous
/// version collected every `SeqRecord` first, then copied seq + phred into
/// a second structure — the two full copies together were ~0.6 GB on G37
/// full. Pairing is irrelevant for assembly (BCALM semantics): every record
/// from every file contributes its k-mers in order.
pub(crate) fn read_records(infiles: &[String]) -> Result<Vec<(Vec<u8>, Vec<u8>)>> {
    anyhow::ensure!(!infiles.is_empty(), "at least one input file is required");
    let mut reads: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    let mut rec = SeqRecord::new();
    for infile in infiles {
        let mut reader = SeqReader::new(infile)?;
        while reader.read_record(&mut rec)? {
            canonicalize_quality(&mut rec);
            reads.push((
                rec.sequence().to_vec(),
                to_phred(rec.sequence(), rec.quality_scores()),
            ));
        }
    }
    Ok(reads)
}

/// One maximal unitig (non-branching path; BCALM `graph3` semantics).
#[derive(Clone)]
pub(crate) struct Unitig {
    pub(crate) bases: Vec<u8>,
    pub(crate) id: usize,
    pub(crate) coverage: f32,
    pub(crate) min_cov: usize,
    pub(crate) max_cov: usize,
    /// The k-mer path closes back on itself (a circular contig).
    pub(crate) circular: bool,
    /// Per-k-mer canonical counts (`ab:Z:` vector); empty when not requested.
    pub(crate) abundances: Vec<u32>,
}

/// Assembles reads into maximal unitigs instead of seeded contigs.
///
/// BCALM-style compaction (`ograph.cpp` `graph3`): every solid k-mer
/// (count >= `min_count_seed`) compresses into its unique non-branching
/// path. A k-mer extends only while it has exactly one solid successor
/// whose own predecessor is also unique; parallel paths stay separate
/// (no bubble popping), and the result is independent of scan order.
pub fn assemble_unitigs<W: Write>(
    infiles: &[String],
    out: &mut W,
    opts: &AssembleOptions,
) -> Result<AssembleStats> {
    let (mut unitigs, mut stats) = assemble_unitigs_core(infiles, opts)?;
    let links = if opts.emit_links || opts.emit_gfa {
        compute_links(&unitigs, opts.k)
    } else {
        vec![Vec::new(); unitigs.len()]
    };
    if opts.emit_gfa {
        writeln!(out, "H\tVN:Z:1.0\tks:i:{}", opts.k)?;
    }
    let min_len = opts.resolved_min_contig_len();
    // GFA segments below min_len are dropped, so an L edge may only reference
    // a kept segment; otherwise the graph would dangle (zero-dangling policy).
    let kept: Vec<bool> = unitigs.iter().map(|u| u.bases.len() >= min_len).collect();
    for (i, u) in unitigs.iter_mut().enumerate() {
        u.id = i;
        if u.bases.len() >= min_len {
            if opts.emit_gfa {
                write!(out, "S\t{}\t{}", u.id, String::from_utf8_lossy(&u.bases))?;
                if !u.abundances.is_empty() {
                    let sum: u32 = u.abundances.iter().sum();
                    let n = u.abundances.len();
                    let mean = if n > 0 { sum as f64 / n as f64 } else { 0.0 };
                    write!(
                        out,
                        "\tLN:i:{}\tKC:i:{}\tkm:f:{}",
                        u.bases.len(),
                        sum,
                        fmt_fixed(mean, 1)
                    )?;
                }
                writeln!(out)?;
                for l in &links[i] {
                    if !kept[l.to] {
                        continue;
                    }
                    writeln!(
                        out,
                        "L\t{}\t{}\t{}\t{}\t{}M",
                        u.id,
                        if l.from_rc { '-' } else { '+' },
                        l.to,
                        if l.to_rc { '-' } else { '+' },
                        opts.k.saturating_sub(1),
                    )?;
                }
            } else if opts.emit_links {
                // Like the GFA `L` edges, the BCALM-style `L:` header entries
                // must only reference kept unitigs (zero-dangling policy).
                let kept_links: Vec<Link> =
                    links[i].iter().filter(|l| kept[l.to]).copied().collect();
                write_unitig(out, u, Some(&kept_links))?;
            } else {
                write_unitig(out, u, None)?;
            }
            stats.contigs_built += 1;
            stats.bases_built += u.bases.len() as u64;
            stats.longest_contig = stats.longest_contig.max(u.bases.len());
        }
    }
    Ok(stats)
}

/// Assembles reads into maximal unitigs and returns (id, bases) in memory,
/// longest-first, skipping the FASTA writer (pipeline composition).
pub fn assemble_unitigs_buf(
    infiles: &[String],
    opts: &AssembleOptions,
) -> Result<Vec<(usize, Vec<u8>)>> {
    let (unitigs, _) = assemble_unitigs_core(infiles, opts)?;
    let min_len = opts.resolved_min_contig_len();
    Ok(unitigs
        .into_iter()
        .enumerate()
        .filter(|(_, u)| u.bases.len() >= min_len)
        .map(|(i, u)| (i, u.bases))
        .collect())
}

/// Shared core of the unitig assemblers: builds and sorts unitigs.
pub(crate) fn assemble_unitigs_core(
    infiles: &[String],
    opts: &AssembleOptions,
) -> Result<(Vec<Unitig>, AssembleStats)> {
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
    // `--supermer` is the default for FASTA input (no quality scores, so
    // counting without quality gating is equivalent); FASTQ falls back to
    // the direct path to keep `min_prob` semantics.
    let fasta_input = if opts.use_supermer {
        let mut reader = SeqReader::new(&infiles[0])?;
        let mut rec = SeqRecord::new();
        if reader.read_record(&mut rec)? {
            rec.quality_scores().is_empty()
        } else {
            true
        }
    } else {
        false
    };
    if opts.use_supermer && !fasta_input {
        eprintln!(
            "note: --supermer requires FASTA input (no quality scores); falling back to direct counting"
        );
    }
    let use_supermer = opts.use_supermer && fasta_input;
    let (table, reads_in) = if use_supermer {
        let t0 = std::time::Instant::now();
        let reads = read_records(infiles)?;
        if std::env::var_os("ANCHR_SM_TIMING").is_some() {
            eprintln!(
                "read_records: {:.3}s ({} reads)",
                t0.elapsed().as_secs_f64(),
                reads.len()
            );
        }
        let n = reads.len() as u64;
        (
            RefineTable::build_supermer(reads, opts.k, opts.supermer_m)?,
            n,
        )
    } else if opts.parallel > 0 {
        let t0 = std::time::Instant::now();
        let (table, n) =
            RefineTable::build_streamed(infiles, opts.k, opts.min_prob, opts.parallel)?;
        if std::env::var_os("ANCHR_SM_TIMING").is_some() {
            eprintln!(
                "streamed count: {:.3}s ({} reads)",
                t0.elapsed().as_secs_f64(),
                n
            );
        }
        (table, n)
    } else {
        // `parallel == 0` (e.g. olc defaults): keep the historical in-memory
        // path on the rayon global pool.
        let reads = read_records(infiles)?;
        let n = reads.len() as u64;
        (
            RefineTable::build_threaded(&reads, opts.k, opts.min_prob, 0),
            n,
        )
    };

    let t_walk = std::time::Instant::now();
    let mut unitigs = if opts.use_dfa {
        let states = VertexStates::classify(&table, opts.min_count_seed as u32, opts.parallel);
        build_unitigs_dfa(&table, opts, &states)
    } else {
        build_unitigs(&table, opts)
    };
    unitigs.sort_by(unitig_cmp);
    if std::env::var_os("ANCHR_SM_TIMING").is_some() {
        eprintln!(
            "walk+build: {:.3}s ({} unitigs)",
            t_walk.elapsed().as_secs_f64(),
            unitigs.len()
        );
    }
    let stats = AssembleStats {
        reads_in,
        ..AssembleStats::default()
    };
    Ok((unitigs, stats))
}

/// One directed unitig link. `from_rc` selects the BCALM FASTA prefix:
/// `false` = out-neighbor (`L:+:` from the owning unitig's right end),
/// `true` = in-neighbor (`L:-:` from its left end). `to_rc` is the target
/// strand.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Link {
    pub(crate) to: usize,
    pub(crate) from_rc: bool,
    pub(crate) to_rc: bool,
}

/// Computes links between unitigs sharing an endpoint (k-1)-mer, matching
/// bcalm's `LinkTigs.cpp` (`link_unitigs_pass`) exactly: every unitig emits
/// in-neighbors from its begin (`L:-:` entries) and out-neighbors from its
/// end (`L:+:` entries), per the four orientation cases each, plus the
/// palindrome exception when `(k-1)` is even and the extremity (k-1)-mer is
/// a palindrome (all orientations admitted).
///
/// bcalm's 8-pass disk partitioning (`is_in_pass`) is only a memory-bounding
/// device: both ends of a shared canonical (k-1)-mer hash to the same pass,
/// so a single in-memory pass has identical semantics.
pub(crate) fn compute_links(unitigs: &[Unitig], k: usize) -> Vec<Vec<Link>> {
    if k < 2 || unitigs.is_empty() {
        return vec![Vec::new(); unitigs.len()];
    }
    let km1 = k - 1;
    let palindrome_exception = km1.is_multiple_of(2);
    // Endpoint (k-1)-mers indexed by canonical form.
    // Entry: (unitig, right end?, rc flag) — the rc flag mirrors bcalm's
    // `ExtremityInfo(!sameOrientation)` (extremities stored canonicalized).
    let mut idx: HashMap<Kmer, Vec<(usize, bool, bool)>> = HashMap::new();
    for (i, u) in unitigs.iter().enumerate() {
        for right in [false, true] {
            let actual = end_kmer1(&u.bases, k, right);
            let canon = actual.canonical();
            let same = canon.cmp_bases(&actual) == std::cmp::Ordering::Equal;
            idx.entry(canon).or_default().push((i, right, !same));
        }
    }
    let mut links = vec![Vec::new(); unitigs.len()];
    for (i, u) in unitigs.iter().enumerate() {
        // in-neighbors from the begin (k-1)-mer (bcalm `L:-:`)
        let begin = end_kmer1(&u.bases, k, false);
        let begin_same = begin.canonical().cmp_bases(&begin) == std::cmp::Ordering::Equal;
        let begin_pal =
            palindrome_exception && begin.cmp_bases(&begin.rc()) == std::cmp::Ordering::Equal;
        if let Some(cands) = idx.get(&begin.canonical()) {
            for &(j, j_right, j_rc) in cands {
                let valid = (begin_same && j_right && !j_rc)
                    || (begin_same && !j_right && j_rc)
                    || (!begin_same && j_right && j_rc)
                    || (!begin_same && !j_right && !j_rc);
                if valid || begin_pal {
                    links[i].push(Link {
                        to: j,
                        from_rc: true,  // L:-:
                        to_rc: j_right, // e_in.pos == UNITIG_END
                    });
                }
            }
        }
        // out-neighbors from the end (k-1)-mer (bcalm `L:+:`)
        let end = end_kmer1(&u.bases, k, true);
        let end_same = end.canonical().cmp_bases(&end) == std::cmp::Ordering::Equal;
        let end_pal = palindrome_exception && end.cmp_bases(&end.rc()) == std::cmp::Ordering::Equal;
        if let Some(cands) = idx.get(&end.canonical()) {
            for &(j, j_right, j_rc) in cands {
                let valid = (end_same && !j_right && !j_rc)
                    || (end_same && j_right && j_rc)
                    || (!end_same && !j_right && j_rc)
                    || (!end_same && j_right && !j_rc);
                if valid || end_pal {
                    links[i].push(Link {
                        to: j,
                        from_rc: false, // L:+:
                        to_rc: j_right, // e_out.pos == UNITIG_END
                    });
                }
            }
        }
    }
    for l in &mut links {
        l.sort_unstable();
        l.dedup();
    }
    links
}

/// The (k-1)-mer at a unitig end as a `Kmer` (actual sequence).
fn end_kmer1(bases: &[u8], k: usize, right: bool) -> Kmer {
    let mut km = Kmer::new(k - 1);
    if right {
        for &b in &bases[bases.len() - (k - 1)..] {
            km.push_right(base_code(b));
        }
    } else {
        for &b in &bases[..k - 1] {
            km.push_right(base_code(b));
        }
    }
    km
}

/// Compresses every solid k-mer into its maximal unitig (order-independent).
fn build_unitigs(table: &RefineTable, opts: &AssembleOptions) -> Vec<Unitig> {
    let k = opts.k;
    let threshold = opts.min_count_seed as u32;
    let want_abundances = opts.all_abundance_counts || opts.emit_gfa;
    let mut visited: HashSet<Kmer, std::hash::BuildHasherDefault<KmerFnvHasher>> =
        HashSet::default();
    let mut unitigs = Vec::new();
    for (seed, seed_count) in table.solid_entries(threshold).iter() {
        if visited.contains(seed) {
            continue;
        }
        // `base_at(0)` is the 3' end (last base pushed); rebuild 5'->3'.
        let mut bb: Vec<u8> = (0..k)
            .map(|i| number_to_base(seed.base_at(k - 1 - i)))
            .collect();
        visited.insert(*seed);
        // Extend right while the path stays non-branching.
        let mut circular = false;
        let mut kmer = rightmost_kmer(&bb, k);
        let mut right_counts: Vec<u32> = Vec::new();
        right_counts.push(*seed_count);
        while let Some(b) = unique_solid_out(&kmer, table, threshold) {
            let mut next = kmer;
            next.push_right(b);
            let canon = next.canonical();
            if unique_solid_in(&next, table, threshold) != 1 {
                break;
            }
            if visited.contains(&canon) {
                circular = true;
                break;
            }
            bb.push(number_to_base(b));
            right_counts.push(table.get_count(&canon));
            visited.insert(canon);
            kmer = next;
        }
        // Extend left by reverse-complementing and extending right.
        let mut rc: Vec<u8> = rev_comp(&bb).collect();
        // Canonical counts of the RC walk in RC order; the final sequence
        // is revcomp(rc), so its per-k-mer counts are the reverse (or, after
        // the canonical flip below, the same) list. This accumulates
        // coverage during the walk and removes the `calc_coverage` second
        // pass (and the `ab:Z` second pass when requested).
        let mut left_counts: Vec<u32> = Vec::with_capacity(rc.len() - k + 1);
        left_counts.push(*seed_count); // first RC k-mer = rc(seed)
        let mut rkmer = rightmost_kmer(&rc, k);
        while let Some(b) = unique_solid_out(&rkmer, table, threshold) {
            let mut next = rkmer;
            next.push_right(b);
            let canon = next.canonical();
            if unique_solid_in(&next, table, threshold) != 1 {
                break;
            }
            if visited.contains(&canon) {
                circular = true;
                break;
            }
            rc.push(number_to_base(b));
            left_counts.push(table.get_count(&canon));
            visited.insert(canon);
            rkmer = next;
        }
        bb = rev_comp(&rc).collect();
        // Canonical orientation, like the contig mode.
        let keep = canonical(&bb);
        if !keep {
            bb = rev_comp(&bb).collect();
        }
        // The RC walk covers the left prefix (positions 0..l-1 of the final
        // sequence, reversed) plus a duplicate of the last right k-mer; the
        // right walk covers positions l..l+r. Reassemble the output-order
        // counts from the two lists (drop `left_counts[0]`, which is the
        // same canonical k-mer as the last right entry).
        let mut counts: Vec<u32> = left_counts[1..].iter().rev().copied().collect();
        counts.extend_from_slice(&right_counts);
        let counts = if keep {
            counts
        } else {
            counts.into_iter().rev().collect()
        };
        let (coverage, min_cov, max_cov) = cov_from_counts(&counts);
        let abundances = if want_abundances { counts } else { Vec::new() };
        unitigs.push(Unitig {
            bases: bb,
            id: 0,
            coverage,
            min_cov,
            max_cov,
            circular,
            abundances,
        });
    }
    unitigs
}

/// DFA-state variant of `build_unitigs`: the classification pass computes
/// every solid vertex's in/out degree once (parallelizable), then the walk
/// uses O(1) state lookups instead of re-scanning the four extension
/// buckets per step. The walk order/decisions are identical to
/// `build_unitigs` (same seed scan, same visited/circular handling).
fn build_unitigs_dfa(
    _table: &RefineTable,
    opts: &AssembleOptions,
    states: &VertexStates,
) -> Vec<Unitig> {
    let k = opts.k;
    let want_abundances = opts.all_abundance_counts || opts.emit_gfa;
    let entries = states.entries();
    let mut visited = vec![0u8; entries.len()];
    let mut unitigs = Vec::new();
    for (seed, seed_count) in entries.iter() {
        let si = states.idx(seed).expect("seed in classification states");
        if visited[si] == 1 {
            continue;
        }
        // `base_at(0)` is the 3' end (last base pushed); rebuild 5'->3'.
        let mut bb: Vec<u8> = (0..k)
            .map(|i| number_to_base(seed.base_at(k - 1 - i)))
            .collect();
        visited[si] = 1;
        // Extend right while the path stays non-branching.
        let mut circular = false;
        let mut kmer = rightmost_kmer(&bb, k);
        let mut right_counts: Vec<u32> = Vec::new();
        right_counts.push(*seed_count);
        while let Some(b) = states.out_base(&kmer) {
            let mut next = kmer;
            next.push_right(b);
            let canon = next.canonical();
            if states.in_count(&next) != 1 {
                break;
            }
            let ci = states.idx(&canon).expect("canon in classification states");
            if visited[ci] == 1 {
                circular = true;
                break;
            }
            bb.push(number_to_base(b));
            right_counts.push(states.count_canonical(&canon));
            visited[ci] = 1;
            kmer = next;
        }
        // Extend left by reverse-complementing and extending right.
        let mut rc: Vec<u8> = rev_comp(&bb).collect();
        let mut left_counts: Vec<u32> = Vec::with_capacity(rc.len() - k + 1);
        left_counts.push(*seed_count); // first RC k-mer = rc(seed)
        let mut rkmer = rightmost_kmer(&rc, k);
        while let Some(b) = states.out_base(&rkmer) {
            let mut next = rkmer;
            next.push_right(b);
            let canon = next.canonical();
            if states.in_count(&next) != 1 {
                break;
            }
            let ci = states.idx(&canon).expect("canon in classification states");
            if visited[ci] == 1 {
                circular = true;
                break;
            }
            rc.push(number_to_base(b));
            left_counts.push(states.count_canonical(&canon));
            visited[ci] = 1;
            rkmer = next;
        }
        bb = rev_comp(&rc).collect();
        // Canonical orientation, like the contig mode.
        let keep = canonical(&bb);
        if !keep {
            bb = rev_comp(&bb).collect();
        }
        let mut counts: Vec<u32> = left_counts[1..].iter().rev().copied().collect();
        counts.extend_from_slice(&right_counts);
        let counts = if keep {
            counts
        } else {
            counts.into_iter().rev().collect()
        };
        let (coverage, min_cov, max_cov) = cov_from_counts(&counts);
        let abundances = if want_abundances { counts } else { Vec::new() };
        unitigs.push(Unitig {
            bases: bb,
            id: 0,
            coverage,
            min_cov,
            max_cov,
            circular,
            abundances,
        });
    }
    unitigs
}

/// The rightmost k-mer of `bb` (5'->3' order, pushed right).
fn rightmost_kmer(bb: &[u8], k: usize) -> Kmer {
    let mut kmer = Kmer::new(k);
    for &b in &bb[bb.len() - k..] {
        kmer.push_right(base_code(b));
    }
    kmer
}

/// The single solid successor of `kmer`, or None at a branch or dead end.
fn unique_solid_out(kmer: &Kmer, table: &RefineTable, threshold: u32) -> Option<u8> {
    let counts = table.fill_right_counts(kmer);
    let mut out = None;
    for b in 0..4u8 {
        if counts[b as usize] >= threshold {
            if out.is_some() {
                return None;
            }
            out = Some(b);
        }
    }
    out
}

/// Number of solid predecessors of `kmer`.
fn unique_solid_in(kmer: &Kmer, table: &RefineTable, threshold: u32) -> usize {
    let counts = table.fill_left_counts(kmer);
    (0..4).filter(|&b| counts[b] >= threshold).count()
}

/// Descending length / coverage / sequence / id order (shared shape with
/// `contig_cmp`, without the branch-code fields).
fn unitig_cmp(a: &Unitig, b: &Unitig) -> std::cmp::Ordering {
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

/// Writes one unitig in FASTA with the contig-mode header fields (no
/// left/right branch codes: unitigs have none).
fn write_unitig<W: Write>(w: &mut W, u: &Unitig, links: Option<&[Link]>) -> Result<()> {
    let (gc, hh, caga) = calc_scalars(&u.bases);
    write!(
        w,
        ">unitig_{},len={},cov={},gc={},min={},max={},hh={},caga={}",
        u.id,
        u.bases.len(),
        fmt_fixed(u.coverage as f64, 1),
        fmt_fixed(gc as f64, 3),
        u.min_cov,
        u.max_cov,
        fmt_fixed(hh as f64, 3),
        fmt_fixed(caga as f64, 3),
    )?;
    if !u.abundances.is_empty() {
        write!(w, ",ab:Z:")?;
        for (i, c) in u.abundances.iter().enumerate() {
            if i > 0 {
                write!(w, " ")?;
            }
            write!(w, "{c}")?;
        }
    }
    if let Some(links) = links {
        for l in links {
            write!(
                w,
                " L:{}:{}:{}",
                if l.from_rc { '-' } else { '+' },
                l.to,
                if l.to_rc { '-' } else { '+' },
            )?;
        }
    }
    if u.circular {
        write!(w, ",circular")?;
    }
    writeln!(w)?;
    for chunk in u.bases.chunks(70) {
        w.write_all(chunk)?;
        w.write_all(b"\n")?;
    }
    Ok(())
}

/// Seeding threshold for pass `i` (Java `minCountSeedCurrent` formula).
fn pass_threshold(opts: &AssembleOptions, i: usize) -> usize {
    let t = (opts.min_count_seed as f64 * opts.contig_pass_mult.powi(i as i32) * 0.92 - 0.25)
        .floor() as i64;
    (opts.min_count_seed as i64 + i as i64)
        .max(t)
        .min(i32::MAX as i64) as usize
}

/// One seeding scan over all table k-mers (BuildThread.processNextTable).
#[allow(clippy::too_many_arguments)]
fn scan_table(
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
    if is_junction(right_max, right_second, opts) {
        let d = is_junction(left_max, left_second, opts);
        return if d {
            (D_BRANCH, calc_ratio(&right))
        } else {
            (F_BRANCH, calc_ratio(&right))
        };
    }
    if is_junction(left_max, left_second, opts) {
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

        let fbranch = is_junction(right_max, right_second, opts);
        let bbranch = is_junction(left_max, left_second, opts);
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
fn calc_coverage(bases: &[u8], table: &RefineTable, k: usize) -> (f32, usize, usize) {
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
fn cov_from_counts(counts: &[u32]) -> (f32, usize, usize) {
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
fn calc_scalars(bases: &[u8]) -> (f32, f32, f32) {
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
fn canonical(bases: &[u8]) -> bool {
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
fn contig_cmp(a: &Contig, b: &Contig) -> std::cmp::Ordering {
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

/// Writes one contig in FASTA (SHORT_NAMES header, 70-column wrap).
fn write_contig<W: Write>(w: &mut W, c: &Contig) -> Result<()> {
    let (gc, hh, caga) = calc_scalars(&c.bases);
    writeln!(
        w,
        ">contig_{},len={},cov={},gc={},min={},max={},hh={},caga={}",
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
fn fmt_fixed(x: f64, decimals: usize) -> String {
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

/// Applies the BBTools phred round-trip to a record's quality scores.
fn canonicalize_quality(rec: &mut SeqRecord) {
    if rec.quality_scores().is_empty() {
        return;
    }
    let seq = rec.sequence().to_vec();
    let raw = rec.quality_scores().to_vec();
    let phred = to_phred(&seq, &raw);
    rec.set_quality(from_phred(&phred));
}

/*--------------------------------------------------------------------*/
/*  Contig graph and bubble popping (Tadpole.processContigs)          */
/*--------------------------------------------------------------------*/

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

type EdgeRef = Rc<RefCell<Edge>>;

impl Contig {
    fn left_kmer(&self, k: usize) -> Kmer {
        let mut kmer = Kmer::new(k);
        for &b in &self.bases[..k] {
            kmer.push_right(base_code(b));
        }
        kmer
    }

    fn right_kmer(&self, k: usize) -> Kmer {
        let mut kmer = Kmer::new(k);
        let n = self.bases.len();
        for &b in &self.bases[n - k..] {
            kmer.push_right(base_code(b));
        }
        kmer
    }

    fn left_forward_branch(&self) -> bool {
        self.left_code == F_BRANCH
    }

    fn right_forward_branch(&self) -> bool {
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

    fn get_right_edge(&self, dest: usize, orientation: Option<u8>) -> Option<EdgeRef> {
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

    fn flip(&mut self, inbound: Option<&[EdgeRef]>) {
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

    fn renumber(&mut self, new_id: usize, inbound: Option<&[EdgeRef]>) {
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
fn remove_all_edges(id: usize, inbound: Option<&[EdgeRef]>, contigs: &mut [Contig]) {
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

fn set_used(id: usize, inbound: Option<&[EdgeRef]>, contigs: &mut [Contig]) {
    contigs[id].used = true;
    remove_all_edges(id, inbound, contigs);
}

fn set_associate(id: usize, inbound: Option<&[EdgeRef]>, contigs: &mut [Contig]) {
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
        if count > 0 && is_junction(left_max, count, opts) {
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
        if count > 0 && is_junction(right_max, count, opts) {
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
        if is_junction(left_max, left_second, opts) {
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
        if is_junction(right_max, right_second, opts) {
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

fn kmer_eq(a: &Kmer, b: &Kmer) -> bool {
    a.canonical().cmp_bases(&b.canonical()) == std::cmp::Ordering::Equal
}

/// BubblePopper over the contig graph (assemble.BubblePopper).
struct BubblePopper {
    contigs: Vec<Contig>,
    dest_map: HashMap<usize, Vec<EdgeRef>>,
    k: usize,
    min_len: usize,
    center: usize,
    dest: usize,
    last_mutual_dest: i64,
    last_mutual_dest_orientation: i64,
    expansions: usize,
    contigs_absorbed: usize,
}

impl BubblePopper {
    fn dest_to_edge_map(&self) -> HashMap<usize, Vec<EdgeRef>> {
        let mut map: HashMap<usize, Vec<EdgeRef>> = HashMap::new();
        for c in &self.contigs {
            if c.used || c.associate {
                continue;
            }
            for e in &c.left_edges {
                map.entry(e.borrow().destination)
                    .or_default()
                    .push(e.clone());
            }
            for e in &c.right_edges {
                map.entry(e.borrow().destination)
                    .or_default()
                    .push(e.clone());
            }
        }
        map
    }

    fn expand(&mut self, center_id: usize) -> usize {
        self.center = center_id;
        let mut count = 0;
        while self.expand_right_simple() {
            count += 1;
        }
        while self.contigs[center_id].right_forward_branch() && self.expand_right() {
            count += 1;
            while self.expand_right_simple() {
                count += 1;
            }
        }
        let left_ok = {
            let c = &self.contigs[center_id];
            (c.left_code != LOOP && c.left_code != DEAD_END && !c.left_edges.is_empty())
                || c.left_forward_branch()
        };
        if left_ok {
            let inbound = self.dest_map.get(&center_id).cloned();
            self.contigs[center_id].flip(inbound.as_deref());
            while self.expand_right_simple() {
                count += 1;
            }
            while self.contigs[center_id].right_forward_branch() && self.expand_right() {
                count += 1;
                while self.expand_right_simple() {
                    count += 1;
                }
            }
        }
        count
    }

    fn expand_right_simple(&mut self) -> bool {
        let center_id = self.center;
        let outbound = self.contigs[center_id].right_edges.clone();
        if outbound.is_empty() || self.contigs[center_id].right_code == LOOP || outbound.len() > 1 {
            return false;
        }
        let left_edge = outbound[0].clone();
        let dest_id = left_edge.borrow().destination;
        let dest_right = left_edge.borrow().dest_right();
        if self.contigs[dest_id].used || dest_id == center_id {
            return false;
        }
        let (outbound_right, right_code) = {
            let d = &self.contigs[dest_id];
            if dest_right {
                (d.right_edges.clone(), d.right_code)
            } else {
                (d.left_edges.clone(), d.left_code)
            }
        };
        if right_code == LOOP {
            return false;
        }
        if !outbound_right.is_empty() {
            if outbound_right.len() > 1 {
                return false;
            }
            if outbound_right[0].borrow().destination != center_id {
                return false;
            }
        }
        if self.count_inbound(center_id, true) > 1 {
            return false;
        }
        if self.count_inbound(dest_id, dest_right) > 1 {
            return false;
        }
        if dest_right {
            let inbound = self.dest_map.get(&dest_id).cloned();
            self.contigs[dest_id].flip(inbound.as_deref());
        }
        self.merge(center_id, dest_id, left_edge)
    }

    fn count_inbound(&self, id: usize, dest_right: bool) -> usize {
        self.dest_map
            .get(&id)
            .map(|v| {
                v.iter()
                    .filter(|e| e.borrow().dest_right() == dest_right)
                    .count()
            })
            .unwrap_or(0)
    }

    fn merge(&mut self, left_id: usize, right_id: usize, left_edge: EdgeRef) -> bool {
        let k = self.k;
        let original_left_len = self.contigs[left_id].bases.len();
        let mut bb: Vec<u8> = self.contigs[left_id].bases.clone();
        {
            let eb = left_edge.borrow();
            if eb.bases.len() > 1 {
                bb.extend_from_slice(&eb.bases[..eb.bases.len() - 1]);
            }
        }
        bb.extend_from_slice(&self.contigs[right_id].bases);
        self.contigs[left_id].bases = bb;
        self.contigs[left_id].right_edges.clear();
        let right_right = self.contigs[right_id].right_edges.clone();
        if right_right.is_empty() {
            self.contigs[left_id].right_edges = Vec::new();
        } else {
            for e in &right_right {
                e.borrow_mut().origin = left_id;
            }
            self.contigs[left_id].right_edges = right_right;
        }
        self.redirect_edges(right_id, left_id, true);
        let inbound_right = self.dest_map.get(&right_id).cloned();
        set_used(right_id, inbound_right.as_deref(), &mut self.contigs);
        let right_len = self.contigs[right_id].bases.len();
        let (right_max_cov, right_min_cov, right_code, right_ratio, right_coverage) = {
            let r = &self.contigs[right_id];
            (
                r.max_cov,
                r.min_cov,
                r.right_code,
                r.right_ratio,
                r.coverage,
            )
        };
        {
            let left = &mut self.contigs[left_id];
            left.max_cov = left.max_cov.max(right_max_cov);
            left.min_cov = left.min_cov.min(right_min_cov);
            left.right_code = right_code;
            left.right_ratio = right_ratio;
            let coverage_sum = left.coverage as f64 * (original_left_len - k + 1) as f64
                + right_coverage as f64 * (right_len - k + 1) as f64;
            left.coverage = (coverage_sum / (left.bases.len() - k + 1) as f64) as f32;
        }
        if self.is_loop(left_id) {
            self.contigs[left_id].left_code = LOOP;
            self.contigs[left_id].right_code = LOOP;
            let inbound = self.dest_map.get(&left_id).cloned();
            remove_all_edges(left_id, inbound.as_deref(), &mut self.contigs);
        }
        self.expansions += 1;
        self.contigs_absorbed += 1;
        true
    }

    fn redirect_edges(&mut self, from: usize, to: usize, dest_right: bool) {
        if from == to {
            return;
        }
        let Some(inbound_from) = self.dest_map.remove(&from) else {
            return;
        };
        let mut inbound_to = self.dest_map.get(&to).cloned().unwrap_or_default();
        for e in &inbound_from {
            if e.borrow().dest_right() == dest_right {
                e.borrow_mut().destination = to;
                inbound_to.push(e.clone());
            }
        }
        if inbound_to.is_empty() {
            self.dest_map.remove(&to);
        } else {
            self.dest_map.insert(to, inbound_to);
        }
    }

    fn is_loop(&self, id: usize) -> bool {
        let c = &self.contigs[id];
        if c.left_code == LOOP && c.right_code == LOOP {
            return true;
        }
        if c.left_edges.len() != 1 || c.right_edges.len() != 1 {
            return false;
        }
        for e in &c.left_edges {
            let e = e.borrow();
            if e.destination != id || !e.dest_right() {
                return false;
            }
        }
        for e in &c.right_edges {
            let e = e.borrow();
            if e.destination != id || e.dest_right() {
                return false;
            }
        }
        if let Some(inbound) = self.dest_map.get(&id) {
            for e in inbound {
                if e.borrow().origin != id {
                    return false;
                }
            }
        }
        true
    }

    fn expand_right(&mut self) -> bool {
        let center_id = self.center;
        self.dest = usize::MAX;
        self.last_mutual_dest = -1;
        self.last_mutual_dest_orientation = -1;
        if !self.contigs[center_id].right_forward_branch()
            || self.contigs[center_id].right_edges.is_empty()
        {
            return false;
        }
        let outbound = self.contigs[center_id].right_edges.clone();
        let Some(left_mid_edge) = self.find_representative_mid_edge(&outbound) else {
            return false;
        };
        let mid_id = left_mid_edge.borrow().destination;
        if self.contigs[mid_id].bases.len() < self.min_len {
            return false;
        }
        let mutual_dest = self.find_mutual_dest(&outbound);
        let mutual_dest_orientation = self.last_mutual_dest_orientation;
        let mutual_dest_right = (mutual_dest_orientation & 2) == 2;
        if mutual_dest < 0 || mutual_dest_orientation < 0 {
            return false;
        }
        let dest_id = mutual_dest as usize;
        if self.contigs[dest_id].used || dest_id == center_id {
            return false;
        }
        if mutual_dest_right && !self.contigs[dest_id].right_forward_branch() {
            return false;
        }
        if !mutual_dest_right && !self.contigs[dest_id].left_forward_branch() {
            return false;
        }
        let dest_outbound = {
            let d = &self.contigs[dest_id];
            if mutual_dest_right {
                d.right_edges.clone()
            } else {
                d.left_edges.clone()
            }
        };
        if dest_outbound.is_empty() {
            return false;
        }
        let mutual_dest2 = self.find_mutual_dest(&dest_outbound);
        if mutual_dest2 < 0 || mutual_dest2 as usize != center_id {
            return false;
        }
        let Some(mid_nodes) = self.fetch_mid_nodes(&outbound, true) else {
            return false;
        };
        if !self.mid_nodes_concur(&mid_nodes) {
            return false;
        }
        if mutual_dest_right {
            let inbound = self.dest_map.get(&dest_id).cloned();
            self.contigs[dest_id].flip(inbound.as_deref());
        }
        let right_mid_edge = self.contigs[mid_id].get_right_edge(dest_id, Some(1));
        let Some(right_mid_edge) = right_mid_edge else {
            return false;
        };
        self.dest = dest_id;
        self.pop(
            center_id,
            dest_id,
            mid_id,
            left_mid_edge,
            right_mid_edge,
            &mid_nodes,
        )
    }

    fn find_representative_mid_edge(&self, edges: &[EdgeRef]) -> Option<EdgeRef> {
        let mut mid_edge: Option<EdgeRef> = None;
        let mut mid_len = 0usize;
        for e in edges {
            let c = &self.contigs[e.borrow().destination];
            let clen = c.bases.len();
            match &mid_edge {
                None => {
                    mid_edge = Some(e.clone());
                    mid_len = clen;
                }
                Some(me) => {
                    let me_depth = me.borrow().depth;
                    let e_depth = e.borrow().depth;
                    if clen >= self.min_len
                        && (mid_len < self.min_len
                            || e_depth > me_depth
                            || (e_depth == me_depth && clen > mid_len))
                    {
                        mid_edge = Some(e.clone());
                        mid_len = clen;
                    }
                }
            }
        }
        mid_edge
    }

    fn find_mutual_dest(&mut self, edges: &[EdgeRef]) -> i64 {
        self.last_mutual_dest = -2;
        self.last_mutual_dest_orientation = -1;
        for e in edges {
            let mid_id = e.borrow().destination;
            if mid_id == self.center {
                return -1;
            }
            let outbound = {
                let mid = &self.contigs[mid_id];
                if e.borrow().dest_right() {
                    mid.left_edges.clone()
                } else {
                    mid.right_edges.clone()
                }
            };
            for o in &outbound {
                let ob = o.borrow();
                if self.last_mutual_dest < 0 {
                    self.last_mutual_dest = ob.destination as i64;
                    self.last_mutual_dest_orientation = (ob.orientation & 2) as i64;
                } else if self.last_mutual_dest != ob.destination as i64
                    || self.last_mutual_dest_orientation != (ob.orientation & 2) as i64
                {
                    return -1;
                }
            }
        }
        self.last_mutual_dest
    }

    fn fetch_mid_nodes(
        &mut self,
        outbound: &[EdgeRef],
        flip_as_needed: bool,
    ) -> Option<Vec<usize>> {
        let mut mid_nodes: Vec<usize> = Vec::new();
        for e in outbound {
            let mid_id = e.borrow().destination;
            if mid_nodes.contains(&mid_id) {
                return None;
            }
            if self.contigs[mid_id].used {
                return None;
            }
            mid_nodes.push(mid_id);
            if flip_as_needed && e.borrow().dest_right() {
                let inbound = self.dest_map.get(&mid_id).cloned();
                self.contigs[mid_id].flip(inbound.as_deref());
            }
        }
        Some(mid_nodes)
    }

    fn mid_nodes_concur(&self, mid_nodes: &[usize]) -> bool {
        let center_id = self.center;
        let dest_id = self.dest;
        let mut left_dest: i64 = -1;
        let mut right_dest: i64 = -1;
        for &mid_id in mid_nodes {
            let c = &self.contigs[mid_id];
            if c.left_edges.is_empty() || c.right_edges.is_empty() {
                return false;
            }
            for e in &c.left_edges {
                let eb = e.borrow();
                if left_dest < 0 {
                    left_dest = eb.destination as i64;
                } else if left_dest != eb.destination as i64 {
                    return false;
                }
                if eb.origin == eb.destination {
                    return false;
                }
            }
            for e in &c.right_edges {
                let eb = e.borrow();
                if right_dest < 0 {
                    right_dest = eb.destination as i64;
                } else if right_dest != eb.destination as i64 {
                    return false;
                }
                if eb.origin == eb.destination {
                    return false;
                }
            }
            let incoming = self.dest_map.get(&mid_id);
            let Some(incoming) = incoming else {
                return false;
            };
            for e in incoming {
                let origin = e.borrow().origin;
                if origin != center_id && origin != dest_id {
                    return false;
                }
            }
        }
        if left_dest >= 0 && left_dest as usize != center_id {
            return false;
        }
        if right_dest >= 0 && right_dest as usize != dest_id {
            return false;
        }
        left_dest >= 0 && right_dest >= 0
    }

    fn pop(
        &mut self,
        left_id: usize,
        right_id: usize,
        mid_id: usize,
        left_mid_edge: EdgeRef,
        right_mid_edge: EdgeRef,
        mid_nodes: &[usize],
    ) -> bool {
        let k = self.k;
        let original_left_len = self.contigs[left_id].bases.len();
        let mut bb: Vec<u8> = self.contigs[left_id].bases.clone();
        {
            let eb = left_mid_edge.borrow();
            if eb.bases.len() > 1 {
                bb.extend_from_slice(&eb.bases[..eb.bases.len() - 1]);
            }
        }
        {
            let mid = &self.contigs[mid_id];
            let lim = mid.bases.len() - k + 1;
            if k - 1 < lim {
                bb.extend_from_slice(&mid.bases[k - 1..lim]);
            }
        }
        {
            let eb = right_mid_edge.borrow();
            if eb.bases.len() > 1 {
                bb.extend_from_slice(&eb.bases[..eb.bases.len() - 1]);
            }
        }
        bb.extend_from_slice(&self.contigs[right_id].bases);
        self.contigs[left_id].bases = bb;
        self.contigs[left_id].right_edges.clear();
        let right_right = self.contigs[right_id].right_edges.clone();
        if right_right.is_empty() {
            self.contigs[left_id].right_edges = Vec::new();
        } else {
            for e in &right_right {
                e.borrow_mut().origin = left_id;
            }
            self.contigs[left_id].right_edges = right_right;
        }
        self.redirect_edges(right_id, left_id, true);
        let inbound_right = self.dest_map.get(&right_id).cloned();
        set_used(right_id, inbound_right.as_deref(), &mut self.contigs);
        for &c in mid_nodes {
            let inbound = self.dest_map.get(&c).cloned();
            if c == mid_id {
                set_used(c, inbound.as_deref(), &mut self.contigs);
            } else {
                set_associate(c, inbound.as_deref(), &mut self.contigs);
            }
        }
        let right_len = self.contigs[right_id].bases.len();
        let (right_max_cov, right_min_cov, right_code, right_ratio, right_coverage) = {
            let r = &self.contigs[right_id];
            (
                r.max_cov,
                r.min_cov,
                r.right_code,
                r.right_ratio,
                r.coverage,
            )
        };
        let (mid_max_cov, mid_min_cov) = {
            let m = &self.contigs[mid_id];
            (m.max_cov, m.min_cov)
        };
        {
            let left = &mut self.contigs[left_id];
            left.max_cov = left.max_cov.max(right_max_cov).max(mid_max_cov);
            left.min_cov = left.min_cov.min(right_min_cov).min(mid_min_cov);
            left.right_code = right_code;
            left.right_ratio = right_ratio;
            let coverage_sum = left.coverage as f64 * (original_left_len - k + 1) as f64
                + right_coverage as f64 * (right_len - k + 1) as f64;
            left.coverage = (coverage_sum / (left.bases.len() - k + 1) as f64) as f32;
        }
        if self.is_loop(left_id) {
            self.contigs[left_id].left_code = LOOP;
            self.contigs[left_id].right_code = LOOP;
            let inbound = self.dest_map.get(&left_id).cloned();
            remove_all_edges(left_id, inbound.as_deref(), &mut self.contigs);
        }
        self.expansions += 1;
        self.contigs_absorbed += 1 + mid_nodes.len();
        true
    }

    fn remove_dead_edges(&self, c: &mut Contig) {
        c.left_edges.retain(|e| {
            let d = e.borrow().destination;
            let dc = &self.contigs[d];
            !(dc.used || dc.associate)
        });
        c.right_edges.retain(|e| {
            let d = e.borrow().destination;
            let dc = &self.contigs[d];
            !(dc.used || dc.associate)
        });
    }
}

/// `Tadpole.popBubbles`: one bubble-popping pass, then deterministic sort and
/// renumbering.
fn pop_bubbles(contigs: &mut Vec<Contig>, opts: &AssembleOptions) {
    let dest_map = {
        let mut map: HashMap<usize, Vec<EdgeRef>> = HashMap::new();
        for c in contigs.iter() {
            if c.used || c.associate {
                continue;
            }
            for e in &c.left_edges {
                map.entry(e.borrow().destination)
                    .or_default()
                    .push(e.clone());
            }
            for e in &c.right_edges {
                map.entry(e.borrow().destination)
                    .or_default()
                    .push(e.clone());
            }
        }
        map
    };
    let mut bp = BubblePopper {
        contigs: std::mem::take(contigs),
        dest_map,
        k: opts.k,
        min_len: 2 * opts.k - 1,
        center: 0,
        dest: usize::MAX,
        last_mutual_dest: -1,
        last_mutual_dest_orientation: -1,
        expansions: 0,
        contigs_absorbed: 0,
    };
    for i in 0..bp.contigs.len() {
        let c = &bp.contigs[i];
        if !c.used && (c.left_forward_branch() || c.right_forward_branch()) {
            bp.expand(i);
        }
    }
    let dest_map2 = bp.dest_to_edge_map();
    let mut temp: Vec<Contig> = Vec::new();
    for i in 0..bp.contigs.len() {
        if bp.contigs[i].used {
            continue;
        }
        let mut c = bp.contigs[i].clone();
        bp.remove_dead_edges(&mut c);
        temp.push(c);
    }
    temp.sort_by(contig_cmp);
    for (new_id, c) in temp.iter_mut().enumerate() {
        let inbound = dest_map2.get(&c.id).cloned();
        c.renumber(new_id, inbound.as_deref());
    }
    *contigs = temp;
}

/// Deterministic longest-first sort and renumbering for the no-bubbles path
/// (bubble popping performs the same step while also renumbering edges).
fn finalize_contigs(contigs: &mut [Contig]) {
    contigs.sort_by(contig_cmp);
    for (new_id, c) in contigs.iter_mut().enumerate() {
        c.renumber(new_id, None);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unitig(bases: &[u8]) -> Unitig {
        Unitig {
            bases: bases.to_vec(),
            id: 0,
            coverage: 0.0,
            min_cov: 0,
            max_cov: 0,
            circular: false,
            abundances: Vec::new(),
        }
    }

    fn assert_link(links: &[Link], to: usize, from_rc: bool, to_rc: bool) {
        assert!(
            links
                .iter()
                .any(|l| l.to == to && l.from_rc == from_rc && l.to_rc == to_rc),
            "missing link to {to} ({from_rc},{to_rc}): {links:?}"
        );
    }

    /// BCALM LinkTigs orientation rules (verified against bcalm v2.2.3 on
    /// G37 full: 1482/1482 unitigs with per-unitig `L:` sets identical).
    #[test]
    fn links_directions_branch_and_rc() {
        // S is a 30 bp random fragment (k = 31 -> k-1 = 30).
        let s: Vec<u8> = b"GCTAAAGACAATTACATAACATACACGTCAG"[..30].to_vec();
        assert_eq!(s.len(), 30);
        let poly_a: Vec<u8> = b"A".repeat(50);
        // Random filler fragments (poly-C/poly-G would share a canonical
        // (k-1)-mer: each is the other's reverse complement).
        let x1: Vec<u8> = b"TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAA".to_vec();
        let x2: Vec<u8> = b"CCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAG".to_vec();
        let s_rc: Vec<u8> = rev_comp(&s).collect();

        // Branch: U0's right end and U1/U2's left ends all share S.
        let uts = vec![
            unitig(&[&poly_a[..], &s[..]].concat()),
            unitig(&[&s[..], &x1[..]].concat()),
            unitig(&[&s[..], &x2[..]].concat()),
        ];
        let links = compute_links(&uts, 31);
        // U0 emits out-neighbors to U1/U2 (`L:+:`); U1/U2 emit the
        // in-neighbor from U0 (`L:-:0:-`).
        assert_eq!(links[0].len(), 2);
        assert_link(&links[0], 1, false, false);
        assert_link(&links[0], 2, false, false);
        assert_eq!(
            links[1],
            vec![Link {
                to: 0,
                from_rc: true,
                to_rc: true
            }]
        );
        assert_eq!(
            links[2],
            vec![Link {
                to: 0,
                from_rc: true,
                to_rc: true
            }]
        );

        // Reverse: U0's right end is rc(S), U1's left end is S. The two
        // share a canonical (k-1)-mer, but the actual extremities differ
        // (rc(S) != S), so bcalm's orientation cases emit no link.
        let uts = vec![
            unitig(&[&poly_a[..], &s_rc[..]].concat()),
            unitig(&[&s[..], &x1[..]].concat()),
        ];
        let links = compute_links(&uts, 31);
        assert!(links[0].is_empty() && links[1].is_empty());

        // 3'-3': both right ends are S (same actual (k-1)-mer) — also no
        // link under bcalm's out-neighbor cases (end-end with rc=false).
        let uts = vec![
            unitig(&[&poly_a[..], &s[..]].concat()),
            unitig(&[&x2[..], &s[..]].concat()),
        ];
        let links = compute_links(&uts, 31);
        assert!(links[0].is_empty() && links[1].is_empty());

        // Self-link: a unitig whose begin and end share the same (k-1)-mer
        // (e.g. a poly-C run) links to itself in both directions — bcalm
        // emits `L:-:<id>:-` and `L:+:<id>:+` (verified on G37 unitig 178).
        let uts = vec![unitig(&b"C".repeat(60))];
        let links = compute_links(&uts, 31);
        assert_eq!(
            links[0],
            vec![
                Link {
                    to: 0,
                    from_rc: false,
                    to_rc: false
                },
                Link {
                    to: 0,
                    from_rc: true,
                    to_rc: true
                },
            ]
        );
    }
}

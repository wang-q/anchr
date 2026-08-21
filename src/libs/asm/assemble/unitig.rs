//! BCALM-style maximal unitig construction (graph3 semantics): the walk
//! engine (DFA states or successor index), link computation between
//! unitigs, and the unitig-mode entry points.

use super::contig::{calc_scalars, canonical, cov_from_counts, fmt_fixed};
use super::{read_records, AssembleOptions, AssembleStats};
use crate::libs::asm::dfa::VertexStates;
use crate::libs::asm::refine::number_to_base;
use crate::libs::asm::table::{base_code, Kmer, KmerFnvHasher, RefineTable};
use anyhow::Result;
use pgr::libs::fmt::seq::{SeqReader, SeqRecord};
use pgr::libs::kmer::key;
use pgr::libs::nt::rev_comp;
use std::collections::{HashMap, HashSet};
use std::io::Write;

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
    // Guard before the `--supermer` FASTA probe (`infiles[0]`) below; all
    // CLI callers check this first, but the empty case must not panic here.
    anyhow::ensure!(!infiles.is_empty(), "at least one input file is required");
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

    let unitigs = assemble_unitigs_from_table(&table, opts);
    let stats = AssembleStats {
        reads_in,
        ..AssembleStats::default()
    };
    Ok((unitigs, stats))
}

/// Classifies solid k-mers and walks maximal unitigs from a prebuilt count
/// table (the shared tail of [`assemble_unitigs_core`]). multik reuses one
/// reads-only table per k across every master's pass 0 instead of
/// recounting the reads for each master.
pub(crate) fn assemble_unitigs_from_table(
    table: &RefineTable,
    opts: &AssembleOptions,
) -> Vec<Unitig> {
    let t_walk = std::time::Instant::now();
    let mut unitigs = if opts.use_dfa {
        let states = VertexStates::classify(table, opts.min_count_seed as u32, opts.parallel);
        build_unitigs_dfa(table, opts, &states)
    } else {
        build_unitigs(table, opts)
    };
    unitigs.sort_by(unitig_cmp);
    if std::env::var_os("ANCHR_SM_TIMING").is_some() {
        eprintln!(
            "walk+build: {:.3}s ({} unitigs)",
            t_walk.elapsed().as_secs_f64(),
            unitigs.len()
        );
    }
    unitigs
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
        // sequence, reversed) plus a duplicate of the seed k-mer; the right
        // walk covers the seed onward (positions l..l+r). Reassemble the
        // output-order counts from the two lists (drop `left_counts[0]`,
        // which is the same canonical k-mer as the first right entry).
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
    for (si, (seed, seed_count)) in entries.iter().enumerate() {
        // `si` is the seed's own index in the parallel states/counts arrays
        // (entries are unique and the index maps key -> enumeration order),
        // so no HashMap lookup is needed here.
        if visited[si] == 1 {
            continue;
        }
        // `base_at(0)` is the 3' end (last base pushed); rebuild 5'->3'.
        let mut bb: Vec<u8> = (0..k)
            .map(|i| number_to_base(seed.base_at(k - 1 - i)))
            .collect();
        visited[si] = 1;
        // Extend right while the path stays non-branching. The window keeps
        // the forward k-mer and its reverse complement in lockstep (two
        // packed shifts per step); classification prelinked each vertex's
        // unique continuation, so a step is plain array reads plus a byte
        // compare for the strand.
        let mut circular = false;
        let mut fw = rightmost_kmer(&bb, k);
        let mut rcm = fw.rc();
        let mut idx = states
            .canon_idx_pair(&fw, &rcm)
            .expect("seed in classification states");
        let mut right_counts: Vec<u32> = Vec::new();
        right_counts.push(*seed_count);
        while let Some((b, ci)) =
            states.step(fw.cmp_bases(&rcm) != std::cmp::Ordering::Greater, idx)
        {
            fw.push_right(b);
            rcm.push_left(3 - b);
            if states.in_count_at(fw.cmp_bases(&rcm) != std::cmp::Ordering::Greater, ci) != 1 {
                break;
            }
            if visited[ci] == 1 {
                circular = true;
                break;
            }
            bb.push(number_to_base(b));
            right_counts.push(states.count_at(ci));
            visited[ci] = 1;
            idx = ci;
        }
        // Extend left by reverse-complementing and extending right.
        let mut rc: Vec<u8> = rev_comp(&bb).collect();
        let mut left_counts: Vec<u32> = Vec::with_capacity(rc.len() - k + 1);
        left_counts.push(*seed_count); // first RC k-mer = rc(seed)
                                       // Same rolling fw/rc pair as the right extension (the walked strand
                                       // is the reverse complement here).
        let mut fw = rightmost_kmer(&rc, k);
        let mut rcm = fw.rc();
        let mut idx = states
            .canon_idx_pair(&fw, &rcm)
            .expect("seed in classification states");
        while let Some((b, ci)) =
            states.step(fw.cmp_bases(&rcm) != std::cmp::Ordering::Greater, idx)
        {
            fw.push_right(b);
            rcm.push_left(3 - b);
            if states.in_count_at(fw.cmp_bases(&rcm) != std::cmp::Ordering::Greater, ci) != 1 {
                break;
            }
            if visited[ci] == 1 {
                circular = true;
                break;
            }
            rc.push(number_to_base(b));
            left_counts.push(states.count_at(ci));
            visited[ci] = 1;
            idx = ci;
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
        ">unitig_{} len={},cov={},gc={},min={},max={},hh={},caga={}",
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

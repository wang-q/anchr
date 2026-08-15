//! Multi-k iterative assembly: unitig graph cross-round validation.
//!
//! MetaMDBG-style iteration (see notes/design/asm-multik.md): pass 0
//! assembles maximal unitigs at the first k, then every later k validates
//! the previous unitig graph — each link's bridge k-mer (the current-k
//! window covering the shared (k_prev-1)-mer junction) must be solid in
//! reads + previous unitigs, and every internal k-mer of a long-enough
//! unitig must stay solid (chimeric-unitig cleanup). Surviving links are
//! compacted into longer unitigs, so validated connections grow
//! monotonically instead of being re-split by a larger k.

use super::assemble::{
    assemble_unitigs_core, compute_links, read_records, AssembleOptions, Link, Unitig,
};
use super::refine::{base_code, Kmer as TdKmer, RefineTable};
use anyhow::Result;
use pgr::libs::kmer::key::Kmer;
use pgr::libs::nt::rev_comp;
use std::collections::HashSet;

/// Options for [`assemble_multik`].
#[derive(Debug, Clone)]
pub struct MultikOptions {
    /// Increasing k-mer lengths; `ks[0]` is the master k (the skeleton,
    /// assembled by pass 0) and each later k is a slave k that validates the
    /// graph. Run one master per invocation and merge masters' outputs in
    /// the template (bash parallel). Empty means auto-derive from the
    /// read-length N50 (see [`auto_ks`]).
    pub ks: Vec<usize>,
    /// Solid k-mer threshold for pass 0 (same as `asm unitig`).
    pub min_count_seed: usize,
    /// Solid k-mer threshold for cross-round validation (metaMDBG `>= 2`).
    pub min_count_extend: usize,
    /// Worker threads for counting; `0` = rayon global pool.
    pub parallel: usize,
}

impl Default for MultikOptions {
    fn default() -> Self {
        Self {
            ks: Vec::new(),
            min_count_seed: 3,
            min_count_extend: 2,
            parallel: 0,
        }
    }
}

/// One multi-k output unitig (a validated chain).
#[derive(Debug, Clone)]
pub struct MultikUnitig {
    pub bases: Vec<u8>,
    pub coverage: f32,
}

/// Assembles reads into long unitigs by iterating over increasing k.
pub fn assemble_multik(infiles: &[String], opts: &MultikOptions) -> Result<Vec<MultikUnitig>> {
    let mut ks = if opts.ks.is_empty() {
        auto_ks(read_n50(&read_records(infiles)?))
    } else {
        opts.ks.clone()
    };
    anyhow::ensure!(!ks.is_empty(), "at least one k-mer length is required");
    for &k in &ks {
        anyhow::ensure!(
            (1..=Kmer::MAX_K).contains(&k),
            "k-mer length must be in 1..={}, got {k}",
            Kmer::MAX_K
        );
    }
    ks.sort_unstable();
    ks.dedup();

    let k0 = ks[0];
    // Single-master-skeleton iteration: `ks[0]` builds the graph and every
    // later k validates it. One master per invocation — the template drives
    // several masters in parallel (bash) and merges their outputs.
    let mut chains = assemble_one(infiles, k0, &ks[1..], opts)?;
    chains.sort_by_key(|u| std::cmp::Reverse(u.bases.len()));
    Ok(chains)
}

/// Runs the single-skeleton iteration: pass 0 builds maximal unitigs at
/// `k0`, every k in `later_ks` validates the graph (bridge k-mers, internal
/// solidity, reads bridges), and validated chains are compacted into long
/// unitigs.
fn assemble_one(
    infiles: &[String],
    k0: usize,
    later_ks: &[usize],
    opts: &MultikOptions,
) -> Result<Vec<MultikUnitig>> {
    // Pass 0: maximal unitigs at the first k (BCALM graph3 semantics).
    let assemble_opts = AssembleOptions {
        k: k0,
        min_count_seed: opts.min_count_seed,
        parallel: opts.parallel,
        // FASTA input: FastK-style super-mer two-stage counting (the `asm
        // unitig` default); FASTQ falls back to the direct path inside
        // `assemble_unitigs_core`.
        use_supermer: true,
        // Same fast paths as the `asm unitig` command: adaptive minimizer
        // length and the DFA-state walk engine.
        supermer_m: Some((12).min((5).max(k0 / 4))),
        use_dfa: true,
        ..AssembleOptions::default()
    };
    let (mut unitigs, _) = assemble_unitigs_core(infiles, &assemble_opts)?;
    let mut links = compute_links(&unitigs, k0);
    // Repeat fragments (a repeated element's shared k-mers fan into four or
    // more flanking unitigs at the master k) connect to four or more
    // partners. Snapshot that branch status at pass 0: a chain through such
    // a node picks one of several genomic contexts and joins distant loci,
    // so its links never participate in recompaction. The flag propagates
    // through every unitig reindexing below.
    let mut branch: Vec<bool> = links
        .iter()
        .map(|ls| ls.iter().map(|l| l.to).collect::<HashSet<_>>().len() >= 4)
        .collect();

    // Iterative rounds: validate the graph with each larger k. Low-abundance
    // branching unitigs are pruned every round but CARRIED into the next
    // round's graph (megahit bubble re-feeding + metaMDBG unitig feedback):
    // strain/polymorphic sequences keep participating in assembly instead of
    // being emitted as dropped fragments only.
    let mut carried: Vec<Unitig> = Vec::new();
    let mut carried_branch: Vec<bool> = Vec::new();
    for &k in later_ks {
        unitigs.append(&mut carried);
        branch.append(&mut carried_branch);
        let t_round = std::time::Instant::now();
        let t0 = std::time::Instant::now();
        let table = count_at(&unitigs, infiles, k, opts.parallel)?;
        let t_count = t0.elapsed().as_secs_f64();
        let threshold = opts.min_count_extend as u32;
        // 1. Cross-round link validation (solveEdges): the bridge k-mer
        // covering the junction must be solid at the current k.
        for (i, ls) in links.iter_mut().enumerate() {
            ls.retain(|l| {
                let u = &unitigs[i];
                let v = &unitigs[l.to];
                // Short unitigs (shorter than the current k-1 window) cannot
                // provide a full bridge k-mer; skip their validation and let
                // the final compaction decide via actual extremity matching
                // (the link is guaranteed to share a (k0-1)-mer).
                if u.bases.len() < k - 1 || v.bases.len() < k - 1 {
                    return true;
                }
                bridge_kmer(u, v, l, k, k0).is_some_and(|km| table.get_count(&km) >= threshold)
            });
        }
        // 2. Chimeric-unitig cleanup (removeUnsupportedUnitigs): every
        // internal current-k k-mer of a long-enough unitig must be solid.
        remove_unsupported(&mut unitigs, &mut links, &mut branch, &table, k, threshold)?;
        // 2.5 Reads-bridge validation: every surviving link must have reads
        // fully covering a probe spanning the junction. Chimeric links (two
        // distant regions joined by a shared k-mer) have no bridging reads
        // and are dropped BEFORE recompaction, so the per-round merge cannot
        // fix them into the main path (prevents relocation misassemblies).
        bridge_filter(&unitigs, &mut links, infiles, k0, 30, 2)?;
        // 3. Recompact unique chains so the main path grows between rounds
        // (metaMDBG recompacts after every abundance-removal round). No
        // abundance pruning here — that is deferred to the final filter, so
        // single-genome coverage fluctuation never drops real content.
        recompact_graph(&mut unitigs, &mut links, &mut branch, k0);
        // 4. Prune low-abundance branching/isolated unitigs and carry them
        // into the next round (the final round's carry becomes output).
        let (kept, kept_branch) = progressive_filter(&mut unitigs, &mut links, &mut branch, k0);
        carried = kept;
        carried_branch = kept_branch;
        if std::env::var_os("ANCHR_MULTIK_TIMING").is_some() {
            let n = unitigs.len();
            let bp: usize = unitigs.iter().map(|u| u.bases.len()).sum();
            let edges: usize = links.iter().map(|l| l.len()).sum();
            eprintln!(
                "round k={k}: {n} unitigs, {bp} bp, {edges} edges, count {t_count:.3}s graph {:.3}s total {:.3}s",
                t_round.elapsed().as_secs_f64() - t_count,
                t_round.elapsed().as_secs_f64()
            );
        }
    }

    // Split unitigs at internal positions that have no reads support (the
    // abundance filter's recompaction can fuse chimeric links into a single
    // unitig — the source of G37 relocations). Every 100-mer window of a
    // unitig must occur in the reads; an unsupported window is a chimeric
    // junction and the unitig is cut there.
    split_by_bridge(&mut unitigs, &mut links, &mut branch, infiles, k0, 30, 1)?;
    // Re-verify the links recomputed by the split: the new extremities may
    // join distant regions, so every surviving link needs bridging reads.
    bridge_filter(&unitigs, &mut links, infiles, k0, 30, 2)?;

    // Megahit-style cleaning on the final (compacted) unitigs: drop short
    // low-depth tips and disconnect weak links (depth-proportional, see
    // megahit tip_remover / weak_link_remover). Doing this per round was
    // too aggressive on the k0=21 fragments (G37 longest 52.8k -> 32.6k);
    // after compaction the tips are real and few.
    tip_remover(&mut unitigs, &mut links, &mut branch, k0 * 2, 20.0);
    weak_link_remover(&mut unitigs, &mut links, 0.05);

    // Final compaction: merge validated chains into long unitigs.
    let mut chains = merge_chains(&unitigs, &links, &branch, k0)?;
    chains.extend(carried.into_iter().map(|u| MultikUnitig {
        bases: u.bases,
        coverage: u.coverage,
    }));
    chains.sort_by_key(|u| std::cmp::Reverse(u.bases.len()));
    Ok(chains)
}

/// Read-length N50 of the input records (used to derive the k sequence).
fn read_n50(reads: &[(Vec<u8>, Vec<u8>)]) -> usize {
    let mut lens: Vec<usize> = reads.iter().map(|(s, _)| s.len()).collect();
    lens.sort_unstable();
    let total: usize = lens.iter().sum();
    let mut acc = 0usize;
    for &l in lens.iter().rev() {
        acc += l;
        if acc * 2 >= total {
            return l;
        }
    }
    0
}

/// Derives an increasing k sequence from the read-length N50:
/// `k_max = min(0.8 * N50, 128)`, starting at `clamp(N50/10, 21, 31)` with
/// steps of `clamp(N50/100, 20, 30)`. The starting k is ~1/3 of the read
/// length (`clamp(N50/3, 31, 51)`): pass 0 builds the graph skeleton at the
/// first k, and a too-small skeleton (21-mer on 150 bp reads) fragments the
/// assembly at low-complexity junctions that larger k's could resolve
/// (MG1655: N50 21K at k0=21 vs 59K at k0=51). Short reads (150 bp) yield
/// 50/70/90/110; long reads (>= 10 kb) cap at 51/81/111 — mirroring
/// metaMDBG's `computeLastK` (last k-min-mer spans ~2× N50).
fn auto_ks(n50: usize) -> Vec<usize> {
    if n50 == 0 {
        return Vec::new();
    }
    let k_max = (n50 * 8 / 10).clamp(31, Kmer::MAX_K);
    let k_min = (n50 / 3).clamp(31, 51);
    let step = (n50 / 100).clamp(20, 30);
    let mut ks = Vec::new();
    let mut k = k_min;
    while k <= k_max {
        ks.push(k);
        k += step;
    }
    ks
}

/// Counts current-k k-mers over reads plus the previous unitigs (no quality
/// gating: unitigs carry no phred scores; the supermer path handles both).
fn count_at(
    unitigs: &[Unitig],
    infiles: &[String],
    k: usize,
    _parallel: usize,
) -> Result<RefineTable> {
    let mut reads = read_records(infiles)?;
    for u in unitigs {
        reads.push((u.bases.clone(), Vec::new()));
    }
    let table = RefineTable::build_supermer(reads, k, None)?;
    Ok(table)
}

/// Encodes a base slice into the assembly k-mer key (canonical lookup is
/// applied by `RefineTable::get_count`).
fn kmer_from_bases(bases: &[u8], k: usize) -> Option<TdKmer> {
    if bases.len() != k {
        return None;
    }
    let mut km = TdKmer::new(k);
    for &b in bases {
        km.push_right(base_code(b));
    }
    Some(km)
}

/// The current-k k-mer covering the `u → v` junction.
///
/// Direction is resolved by matching the actual extremity (k_prev-1)-mers,
/// not by interpreting the symbolic `to_rc` flag: `from_rc=false` (u's
/// right/3' end) puts u upstream and v (or rc(v)) downstream; `from_rc=true`
/// (u's left/5' end) puts the partner upstream and u downstream. The bridge
/// window is `upstream tail (k-1) + downstream continuation base`, which
/// covers the shared (k_prev-1)-mer and is not contained in either unitig
/// alone (so unitig self-counting cannot support it).
fn bridge_kmer(u: &Unitig, v: &Unitig, link: &Link, k: usize, k_build: usize) -> Option<TdKmer> {
    // Unitigs shorter than the previous k cannot provide the shared
    // (k_prev-1)-mer extremity plus a continuation base; treat their links
    // as unsupported (they still survive as short unitigs).
    if u.bases.len() < k_build || v.bases.len() < k_build {
        return None;
    }
    let km1 = k - 1;
    if link.from_rc {
        // u's left end is the junction source: partner upstream, u downstream.
        let u_ext = u.bases[..k_build - 1].to_vec();
        let v_begin = &v.bases[..k_build - 1];
        let v_end = &v.bases[v.bases.len() - (k_build - 1)..];
        let up: Vec<u8> = if u_ext == v_end[..] {
            v.bases.clone()
        } else if u_ext == rev_comp(v_begin).collect::<Vec<u8>>().as_slice() {
            rev_comp(&v.bases).collect()
        } else {
            return None;
        };
        let tail = up.len().saturating_sub(km1);
        let mut seq: Vec<u8> = Vec::with_capacity(k);
        seq.extend_from_slice(&up[tail..]);
        seq.push(u.bases[k_build - 1]);
        kmer_from_bases(&seq, k)
    } else {
        // u's right end is the junction source: u upstream, partner downstream.
        let u_ext = u.bases[u.bases.len() - (k_build - 1)..].to_vec();
        let v_begin = &v.bases[..k_build - 1];
        let v_end = &v.bases[v.bases.len() - (k_build - 1)..];
        let cont: u8 = if u_ext == v_begin[..] {
            v.bases[k_build - 1]
        } else if u_ext == rev_comp(v_end).collect::<Vec<u8>>().as_slice() {
            rev_comp(&v.bases).collect::<Vec<u8>>()[k_build - 1]
        } else {
            return None;
        };
        let tail = u.bases.len().saturating_sub(km1);
        let mut seq: Vec<u8> = Vec::with_capacity(k);
        seq.extend_from_slice(&u.bases[tail..]);
        seq.push(cont);
        kmer_from_bases(&seq, k)
    }
}

/// The probe spanning the `u → v` junction: `probe_half` bases on each side
/// of the shared `(k_build-1)`-mer overlap (u tail + v continuation, or the
/// reverse for a left-end link). Direction is resolved by actual extremity
/// matching (same as [`bridge_kmer`]).
fn probe_kmer(
    u: &Unitig,
    v: &Unitig,
    link: &Link,
    k_build: usize,
    probe_half: usize,
) -> Option<TdKmer> {
    let probe_len = probe_half * 2;
    let km1 = k_build - 1;
    if link.from_rc {
        // u's left end is the junction source: v upstream, u downstream.
        let u_ext = &u.bases[..km1];
        let vb = &v.bases[..km1];
        let ve = &v.bases[v.bases.len() - km1..];
        let v_dir: Vec<u8> = if u_ext == ve {
            v.bases.clone()
        } else if u_ext == rev_comp(vb).collect::<Vec<u8>>().as_slice() {
            rev_comp(&v.bases).collect()
        } else {
            return None;
        };
        if v_dir.len() < probe_half || u.bases.len() < km1 + probe_half {
            return None;
        }
        let mut seq: Vec<u8> = Vec::with_capacity(probe_len);
        seq.extend_from_slice(&v_dir[v_dir.len() - probe_half..]);
        seq.extend_from_slice(&u.bases[km1..km1 + probe_half]);
        kmer_from_bases(&seq, probe_len)
    } else {
        // u's right end is the junction source: u upstream, v downstream.
        let u_ext = &u.bases[u.bases.len() - km1..];
        let vb = &v.bases[..km1];
        let ve = &v.bases[v.bases.len() - km1..];
        let v_dir: Vec<u8> = if u_ext == vb {
            v.bases.clone()
        } else if u_ext == rev_comp(ve).collect::<Vec<u8>>().as_slice() {
            rev_comp(&v.bases).collect()
        } else {
            return None;
        };
        if u.bases.len() < probe_half || v_dir.len() < km1 + probe_half {
            return None;
        }
        let mut seq: Vec<u8> = Vec::with_capacity(probe_len);
        seq.extend_from_slice(&u.bases[u.bases.len() - probe_half..]);
        seq.extend_from_slice(&v_dir[km1..km1 + probe_half]);
        kmer_from_bases(&seq, probe_len)
    }
}

/// Reads-bridge validation: drops links whose junction-spanning probe is
/// not fully covered by at least `threshold` reads (metaMDBG
/// `computeBridgingReads` — a chimeric link joining two distant regions has
/// no reads covering the junction and is pruned). Links whose probe cannot
/// be built (short unitigs) are kept conservatively.
fn bridge_filter(
    unitigs: &[Unitig],
    links: &mut [Vec<Link>],
    infiles: &[String],
    k0: usize,
    probe_half: usize,
    threshold: u32,
) -> Result<()> {
    if unitigs.is_empty() {
        return Ok(());
    }
    let probe_len = probe_half * 2;
    let reads = read_records(infiles)?;
    let table = RefineTable::build_supermer(reads, probe_len, None)?;
    for (i, ls) in links.iter_mut().enumerate() {
        ls.retain(|l| {
            probe_kmer(&unitigs[i], &unitigs[l.to], l, k0, probe_half)
                .map(|p| table.get_count(&p) >= threshold)
                .unwrap_or(false)
        });
    }
    Ok(())
}

/// Splits unitigs at internal windows that are not supported by any read:
/// every `2*probe_half`-mer window of a unitig must occur in the reads (the
/// unitig's own sequence comes from reads, so a window with count 0 is a
/// chimeric junction — the abundance recompaction fused two distant regions).
/// Splitting keeps those junctions out of the final compaction. Links are
/// recomputed from the new extremities.
fn split_by_bridge(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    branch: &mut Vec<bool>,
    infiles: &[String],
    k0: usize,
    probe_half: usize,
    threshold: u32,
) -> Result<()> {
    if unitigs.is_empty() {
        return Ok(());
    }
    let probe_len = probe_half * 2;
    let reads = read_records(infiles)?;
    let table = RefineTable::build_supermer(reads, probe_len, None)?;
    let mut out: Vec<Unitig> = Vec::new();
    let mut out_branch: Vec<bool> = Vec::new();
    for (u, &is_branch) in unitigs.iter().zip(branch.iter()) {
        let n = u.bases.len();
        if n < probe_len {
            out.push(u.clone());
            out_branch.push(is_branch);
            continue;
        }
        // Mark windows without read support.
        let mut cut: Vec<usize> = Vec::new();
        let mut prev_cut = false;
        for i in 0..=n - probe_len {
            let ok = kmer_from_bases(&u.bases[i..i + probe_len], probe_len)
                .map(|km| table.get_count(&km) >= threshold)
                .unwrap_or(false);
            // Start a cut at the beginning of an unsupported run.
            if !ok && !prev_cut {
                cut.push(i);
            }
            prev_cut = !ok;
        }
        if cut.is_empty() {
            out.push(u.clone());
            out_branch.push(is_branch);
            continue;
        }
        // Split at the cut positions.
        let mut pieces: Vec<usize> = Vec::new();
        let mut s = 0usize;
        for &c in &cut {
            if c <= s {
                continue;
            }
            pieces.push(c - s);
            s = c;
        }
        pieces.push(n - s);
        let mut pos = 0usize;
        for len in pieces {
            if len < k0 {
                // Too short to be a unitig (cannot host a (k0-1)-mer end):
                // drop the fragment.
                pos += len;
                continue;
            }
            let mut nu = u.clone();
            nu.bases = u.bases[pos..pos + len].to_vec();
            nu.id = 0;
            out.push(nu);
            out_branch.push(is_branch);
            pos += len;
        }
    }
    *unitigs = out;
    *branch = out_branch;
    *links = compute_links(unitigs, k0);
    Ok(())
}

/// Megahit-style tip removal: short unitigs (<= `max_tip_len`) that are
/// tips (one end has no connection) with depth far below their neighbour
/// (`neighbour > depth_ratio * self`) are error tips and dropped.
fn tip_remover(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    branch: &mut Vec<bool>,
    max_tip_len: usize,
    depth_ratio: f32,
) {
    let keep: Vec<bool> = unitigs
        .iter()
        .enumerate()
        .map(|(i, u)| {
            if u.bases.len() > max_tip_len {
                return true;
            }
            let out = links[i].len();
            let in_deg = links
                .iter()
                .filter(|ls| ls.iter().any(|l| l.to == i))
                .count();
            if out + in_deg == 0 {
                return true; // isolated: handled by the abundance filter
            }
            let is_tip = (out == 0 && in_deg >= 1) || (out >= 1 && in_deg == 0);
            if !is_tip {
                return true;
            }
            // Deepest neighbour depth.
            let mut max_neighbour = 0.0f32;
            for l in links[i].iter() {
                max_neighbour = max_neighbour.max(unitigs[l.to].coverage);
            }
            for (j, ls) in links.iter().enumerate() {
                if ls.iter().any(|l| l.to == i) {
                    max_neighbour = max_neighbour.max(unitigs[j].coverage);
                }
            }
            u.coverage * depth_ratio < max_neighbour
        })
        .collect();
    retain_graph(unitigs, links, branch, &keep);
}

/// Megahit-style weak-link disconnection: at a branching unitig end
/// (out-degree >= 2), a neighbour whose depth is <= `local_ratio` of the
/// total neighbour depth is disconnected (edge dropped, node kept) — the
/// neighbour is likely a different strain sharing this region, not a real
/// continuation.
fn weak_link_remover(unitigs: &mut [Unitig], links: &mut [Vec<Link>], local_ratio: f32) {
    for (i, u) in unitigs.iter().enumerate() {
        if u.bases.is_empty() {
            continue;
        }
        // For each end (from_rc false = right, true = left) with out-degree
        // >= 2, disconnect neighbours below the proportional threshold.
        for from_rc in [false, true] {
            let mut total_depth = 0.0f32;
            let mut depths: Vec<(usize, f32)> = Vec::new();
            for (e, l) in links[i].iter().enumerate() {
                if l.from_rc != from_rc {
                    continue;
                }
                let d = unitigs[l.to].coverage;
                total_depth += d;
                depths.push((e, d));
            }
            if depths.len() <= 1 {
                continue;
            }
            for &(e, d) in &depths {
                if d <= local_ratio * total_depth {
                    // Mark for removal after the loop (retain per link).
                    links[i][e].to = usize::MAX;
                }
            }
        }
    }
    for ls in links.iter_mut() {
        ls.retain(|l| l.to != usize::MAX);
    }
}

/// Drops unitigs whose internal current-k k-mer is missing from the solid
/// table (chimeric cleanup). Unitigs shorter than `k` have no internal
/// k-mer to check and survive until the final compaction (their links were
/// already validated).
fn remove_unsupported(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    branch: &mut Vec<bool>,
    table: &RefineTable,
    k: usize,
    threshold: u32,
) -> Result<()> {
    let keep: Vec<bool> = unitigs
        .iter()
        .map(|u| {
            if u.bases.len() < k {
                return true;
            }
            let n_kmers = u.bases.len() - k + 1;
            // Tolerate a small fraction of missing internal k-mers: a single
            // genome's coverage fluctuates (a few windows below the solid
            // threshold) without making the unitig chimeric. Only unitigs
            // whose internal k-mers are largely unsupported are dropped.
            let max_missing = (n_kmers / 50).max(1);
            let mut missing = 0usize;
            for j in 0..n_kmers {
                let ok = kmer_from_bases(&u.bases[j..j + k], k)
                    .is_some_and(|km| table.get_count(&km) >= threshold);
                if !ok {
                    missing += 1;
                    if missing > max_missing {
                        break;
                    }
                }
            }
            missing <= max_missing
        })
        .collect();
    retain_graph(unitigs, links, branch, &keep);
    Ok(())
}

/// Removes dropped unitigs and remaps surviving ids and link targets.
fn retain_graph(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    branch: &mut Vec<bool>,
    keep: &[bool],
) {
    let mut remap = vec![usize::MAX; keep.len()];
    let mut n = 0usize;
    for (i, &k) in keep.iter().enumerate() {
        if k {
            remap[i] = n;
            n += 1;
        }
    }
    let mut kept_unitigs: Vec<Unitig> = Vec::with_capacity(n);
    let mut kept_links: Vec<Vec<Link>> = Vec::with_capacity(n);
    let mut kept_branch: Vec<bool> = Vec::with_capacity(n);
    for (i, (u, ls)) in unitigs.iter().zip(links.iter()).enumerate() {
        if !keep[i] {
            continue;
        }
        let mut u2 = u.clone();
        u2.id = remap[i];
        let ls2: Vec<Link> = ls
            .iter()
            .filter(|l| keep[l.to])
            .map(|l| Link {
                to: remap[l.to],
                from_rc: l.from_rc,
                to_rc: l.to_rc,
            })
            .collect();
        kept_unitigs.push(u2);
        kept_links.push(ls2);
        kept_branch.push(branch[i]);
    }
    *unitigs = kept_unitigs;
    *links = kept_links;
    *branch = kept_branch;
}

/// Progressive abundance filter (metaMDBG `removeAbundanceNoQueue`):
/// repeatedly drop unitigs below a cutoff that grows ~10% per round from
/// 1.1 up to the graph's maximum abundance, then recompact the surviving
/// graph (merge chains so the main path inherits its flanks' higher
/// abundance and is not dropped by later cutoffs). Dropped unitigs are
/// returned as independent output (metaMDBG keeps them via cutoff
/// snapshots), so low-abundance species are not lost.
fn progressive_filter(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    branch: &mut Vec<bool>,
    k_build: usize,
) -> (Vec<Unitig>, Vec<bool>) {
    if unitigs.is_empty() {
        return (Vec::new(), Vec::new());
    }
    // The cutoff climbs from 1.1 but stops at 25% of the median coverage:
    // metaMDBG's "up to the graph maximum" assumes the main path is one
    // high-abundance unitig, which holds for strain divergence but not for a
    // single genome whose k=21 unitig graph is fragmented (a repeated region
    // can push the max to 600x while the genome is 40x — climbing to the max
    // would prune every real branch into the dropped list).
    let mut covs: Vec<f32> = unitigs.iter().map(|u| u.coverage).collect();
    covs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = covs[covs.len() / 2];
    let cutoff_cap = (median * 0.25).max(1.1);
    let mut dropped: Vec<Unitig> = Vec::new();
    let mut dropped_branch: Vec<bool> = Vec::new();
    let mut t = 1.1f32;
    while t < cutoff_cap && !unitigs.is_empty() {
        // Main-path protection: unitigs on a unique chain (both ends have
        // <= 1 link) are kept regardless of coverage — coverage fluctuates
        // within one genome and a raw `cov < cutoff` rule would delete the
        // main path's lower-covered segments (max_abundance is often set by
        // a repeated region). Only branching junctions and isolated nodes
        // are pruned by abundance.
        let n = unitigs.len();
        let mut out_begin = vec![0usize; n];
        let mut out_end = vec![0usize; n];
        let mut in_begin = vec![0usize; n];
        let mut in_end = vec![0usize; n];
        for (i, ls) in links.iter().enumerate() {
            for l in ls {
                if l.from_rc {
                    out_begin[i] += 1;
                    in_end[l.to] += 1;
                } else {
                    out_end[i] += 1;
                    in_begin[l.to] += 1;
                }
            }
        }
        let keep: Vec<bool> = unitigs
            .iter()
            .enumerate()
            .map(|(i, u)| {
                let branching =
                    out_begin[i] > 1 || out_end[i] > 1 || in_begin[i] > 1 || in_end[i] > 1;
                let connected = out_begin[i] + out_end[i] + in_begin[i] + in_end[i] > 0;
                if branching || !connected {
                    u.coverage >= t
                } else {
                    true // unique chain: main path, never pruned
                }
            })
            .collect();
        for (i, &k) in keep.iter().enumerate() {
            if !k {
                let mut d = unitigs[i].clone();
                d.id = 0;
                dropped.push(d);
                dropped_branch.push(branch[i]);
            }
        }
        if keep.iter().all(|&k| k) {
            t += (t * 0.1).min(10.0);
            continue; // nothing below the current cutoff, raise it
        }
        retain_graph(unitigs, links, branch, &keep);
        // Recompact after removal so the main path grows (metaMDBG
        // recompacts every round); chimeric junctions possibly fused here
        // are cut back by `split_by_bridge` right after this filter.
        recompact_graph(unitigs, links, branch, k_build);
        t += (t * 0.1).min(10.0);
    }
    (dropped, dropped_branch)
}

/// Resolves a directed link into its oriented chain segment `(left node,
/// left rev, right node, right rev)` by matching the actual extremity
/// (k_build-1)-mers (the same matching `bridge_kmer` uses, without building
/// a k-mer). Returns `None` when the extremities do not overlap or a unitig
/// is shorter than the overlap window.
fn oriented_segment(
    i: usize,
    l: &Link,
    unitigs: &[Unitig],
    k_build: usize,
) -> Option<(usize, bool, usize, bool)> {
    let u = &unitigs[i];
    let v = &unitigs[l.to];
    // A unitig shorter than two (k-1)-mer extremities (or 90 bp at small
    // master k) has overlapping begin/end k-mers: a link through it is
    // ambiguous (the same k-mer can pair either end, folding the chain back
    // on itself — chimeric junction bridges from erroneous merged reads).
    // Such fragments stay independent output instead of being compacted
    // into the main path.
    let min_chain_len = (2 * (k_build - 1)).max(90);
    if u.bases.len() < min_chain_len || v.bases.len() < min_chain_len {
        return None;
    }
    let (li, lrev, ri, rrev) = if !l.from_rc {
        let u_ext = &u.bases[u.bases.len() - (k_build - 1)..];
        let vb = &v.bases[..k_build - 1];
        let ve = &v.bases[v.bases.len() - (k_build - 1)..];
        if u_ext == vb {
            (i, false, l.to, false)
        } else if u_ext == rev_comp(ve).collect::<Vec<u8>>().as_slice() {
            (i, false, l.to, true)
        } else {
            return None;
        }
    } else {
        let u_ext = &u.bases[..k_build - 1];
        let vb = &v.bases[..k_build - 1];
        let ve = &v.bases[v.bases.len() - (k_build - 1)..];
        if u_ext == ve {
            (l.to, false, i, false)
        } else if u_ext == rev_comp(vb).collect::<Vec<u8>>().as_slice() {
            (l.to, true, i, false)
        } else {
            return None;
        }
    };
    Some((li, lrev, ri, rrev))
}

/// In-graph recompaction: merge unique chains into longer unitigs and
/// relink the chain endpoints (metaMDBG `UnitigGraph2::recompact`). The
/// merged unitig keeps the chain head's begin (from_rc=true) links and the
/// chain tail's end (from_rc=false) links; external edges pointing into
/// the chain are redirected to the merged unitig.
fn recompact_graph(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    branch: &mut Vec<bool>,
    k_build: usize,
) {
    let n = unitigs.len();
    if n <= 1 {
        return;
    }
    let node_of = |id: usize, rev: bool| 2 * id + rev as usize;
    let mut right_of: Vec<Option<(usize, bool)>> = vec![None; 2 * n];
    let mut left_of: Vec<Option<(usize, bool)>> = vec![None; 2 * n];
    // Strict end uniqueness (SKESA's "predecessor == 1" invariant): a chain
    // segment joins only when BOTH oriented ends have exactly one link.
    // First-pass degree counting replaces first-come-first-served
    // occupancy, so a convergence node (two predecessors) is never
    // swallowed into a chain — the abundance filter, not the link iteration
    // order, decides which predecessor survives.
    let mut out_deg = vec![0usize; 2 * n];
    let mut in_deg = vec![0usize; 2 * n];
    // Deduplicate segments: compute_links emits the same junction from both
    // endpoints (u's out-link and v's in-link), so per-link degree counting
    // would double-count a unique chain segment.
    let mut seen: HashSet<(usize, usize)> = HashSet::new();
    for (i, ls) in links.iter().enumerate() {
        for l in ls {
            if branch[i] || branch[l.to] {
                continue;
            }
            let Some((li, lrev, ri, rrev)) = oriented_segment(i, l, unitigs, k_build) else {
                continue;
            };
            let ln = node_of(li, lrev);
            let rn = node_of(ri, rrev);
            if seen.insert((ln, rn)) {
                out_deg[ln] += 1;
                in_deg[rn] += 1;
            }
        }
    }
    for (i, ls) in links.iter().enumerate() {
        for l in ls {
            let Some((li, lrev, ri, rrev)) = oriented_segment(i, l, unitigs, k_build) else {
                continue;
            };
            let ln = node_of(li, lrev);
            let rn = node_of(ri, rrev);
            if out_deg[ln] != 1
                || in_deg[rn] != 1
                || right_of[ln].is_some()
                || left_of[rn].is_some()
            {
                continue;
            }
            right_of[ln] = Some((ri, rrev));
            left_of[rn] = Some((li, lrev));
        }
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        unitigs[b]
            .bases
            .len()
            .cmp(&unitigs[a].bases.len())
            .then(a.cmp(&b))
    });
    let mut visited = vec![false; n];
    // remap[old_id] = new merged unitig id (None when not yet assigned).
    let mut remap: Vec<Option<usize>> = vec![None; n];
    let mut new_unitigs: Vec<Unitig> = Vec::new();
    let mut new_branch: Vec<bool> = Vec::new();
    // Head/tail of each merged chain, per new id (for link relinking).
    let mut heads: Vec<usize> = Vec::new();
    let mut tails: Vec<usize> = Vec::new();
    for &seed in &order {
        if visited[seed] {
            continue;
        }
        let mut head = seed;
        let mut head_rev = false;
        let mut seen = vec![false; n];
        while let Some((p, prev)) = left_of[node_of(head, head_rev)] {
            if visited[p] || seen[p] {
                break;
            }
            seen[p] = true;
            head = p;
            head_rev = prev;
        }
        visited[head] = true;
        let mut seq: Vec<u8> = if head_rev {
            rev_comp(&unitigs[head].bases).collect()
        } else {
            unitigs[head].bases.clone()
        };
        let mut cov_sum = unitigs[head].coverage;
        let mut cov_n = 1usize;
        let mut chain: Vec<usize> = vec![head];
        let mut cur = head;
        let mut cur_rev = head_rev;
        while let Some((j, next_rev)) = right_of[node_of(cur, cur_rev)] {
            if visited[j] {
                break;
            }
            visited[j] = true;
            let append: Vec<u8> = if next_rev {
                rev_comp(&unitigs[j].bases).collect()
            } else {
                unitigs[j].bases.clone()
            };
            seq.extend_from_slice(&append[k_build - 1..]);
            cov_sum += unitigs[j].coverage;
            cov_n += 1;
            chain.push(j);
            cur = j;
            cur_rev = next_rev;
        }
        let new_id = new_unitigs.len();
        for &old in &chain {
            remap[old] = Some(new_id);
        }
        new_branch.push(chain.iter().any(|&old| branch[old]));
        new_unitigs.push(Unitig {
            bases: seq,
            id: new_id,
            coverage: cov_sum / cov_n as f32,
            min_cov: 0,
            max_cov: 0,
            circular: false,
            abundances: Vec::new(),
        });
        heads.push(chain[0]);
        tails.push(*chain.last().unwrap());
    }

    // Recompute links from the merged extremities: chain endpoints that
    // share a (k_build-1)-mer get linked (the next round re-validates every
    // link at a larger k), and stale edges whose extremities were absorbed
    // by the merge disappear.
    *unitigs = new_unitigs;
    *links = compute_links(unitigs, k_build);
    *branch = new_branch;
}

/// Compacts validated links into maximal chains (u + partner[(k-1)..]).
///
/// A link is traversed only when it is the unique outgoing edge of the
/// current end and the unique incoming edge of the partner end, so bubbles
/// stay separate. Each link is first resolved into an oriented chain
/// segment (`left → right`, with a per-unitig strand flag) by matching the
/// actual extremity (k_build-1)-mers; the walk then follows unique
/// chain-oriented ends (2n nodes: id * 2 + rev).
fn merge_chains(
    unitigs: &[Unitig],
    links: &[Vec<Link>],
    branch: &[bool],
    k_build: usize,
) -> Result<Vec<MultikUnitig>> {
    let n = unitigs.len();
    let node_of = |id: usize, rev: bool| 2 * id + rev as usize;
    let mut right_of: Vec<Option<(usize, bool)>> = vec![None; 2 * n];
    let mut left_of: Vec<Option<(usize, bool)>> = vec![None; 2 * n];
    // Strict end uniqueness (SKESA's "predecessor == 1" invariant), same as
    // recompact_graph: only segments whose oriented ends both have exactly
    // one link join a chain, so a convergence node is never swallowed into
    // a chain by link iteration order.
    let mut out_deg = vec![0usize; 2 * n];
    let mut in_deg = vec![0usize; 2 * n];
    // Deduplicate segments: compute_links emits the same junction from both
    // endpoints (u's out-link and v's in-link), so per-link degree counting
    // would double-count a unique chain segment.
    let mut seen: HashSet<(usize, usize)> = HashSet::new();
    for (i, ls) in links.iter().enumerate() {
        for l in ls {
            if branch[i] || branch[l.to] {
                continue;
            }
            let Some((li, lrev, ri, rrev)) = oriented_segment(i, l, unitigs, k_build) else {
                continue;
            };
            let ln = node_of(li, lrev);
            let rn = node_of(ri, rrev);
            if seen.insert((ln, rn)) {
                out_deg[ln] += 1;
                in_deg[rn] += 1;
            }
        }
    }
    for (i, ls) in links.iter().enumerate() {
        for l in ls {
            let Some((li, lrev, ri, rrev)) = oriented_segment(i, l, unitigs, k_build) else {
                continue;
            };
            let ln = node_of(li, lrev);
            let rn = node_of(ri, rrev);
            if out_deg[ln] != 1
                || in_deg[rn] != 1
                || right_of[ln].is_some()
                || left_of[rn].is_some()
            {
                continue;
            }
            right_of[ln] = Some((ri, rrev));
            left_of[rn] = Some((li, lrev));
        }
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        unitigs[b]
            .bases
            .len()
            .cmp(&unitigs[a].bases.len())
            .then(a.cmp(&b))
    });
    let mut visited = vec![false; n];
    let mut result: Vec<MultikUnitig> = Vec::new();
    for &seed in &order {
        if visited[seed] {
            continue;
        }
        // Walk back to the chain head (unique left neighbor of the begin).
        let mut head = seed;
        let mut head_rev = false;
        let mut seen = vec![false; n];
        while let Some((p, prev)) = left_of[node_of(head, head_rev)] {
            if visited[p] || seen[p] {
                break;
            }
            seen[p] = true;
            head = p;
            head_rev = prev;
        }
        visited[head] = true;
        let mut seq: Vec<u8> = if head_rev {
            rev_comp(&unitigs[head].bases).collect()
        } else {
            unitigs[head].bases.clone()
        };
        let mut cov_sum = unitigs[head].coverage;
        let mut cov_n = 1usize;
        let mut cur = head;
        let mut cur_rev = head_rev;
        while let Some((j, next_rev)) = right_of[node_of(cur, cur_rev)] {
            if visited[j] {
                break;
            }
            visited[j] = true;
            let append: Vec<u8> = if next_rev {
                rev_comp(&unitigs[j].bases).collect()
            } else {
                unitigs[j].bases.clone()
            };
            seq.extend_from_slice(&append[k_build - 1..]);
            cov_sum += unitigs[j].coverage;
            cov_n += 1;
            cur = j;
            cur_rev = next_rev;
        }
        if !seq.is_empty() {
            result.push(MultikUnitig {
                bases: seq,
                coverage: cov_sum / cov_n as f32,
            });
        }
    }
    result.sort_by_key(|u| std::cmp::Reverse(u.bases.len()));
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static SEQ: AtomicUsize = AtomicUsize::new(0);

    fn temp_fasta(seqs: &[&str]) -> String {
        let mut p = std::env::temp_dir();
        let n = SEQ.fetch_add(1, Ordering::SeqCst);
        p.push(format!("anchr_multik_{}_{}.fa", std::process::id(), n));
        let mut f = std::fs::File::create(&p).unwrap();
        for (i, s) in seqs.iter().enumerate() {
            writeln!(f, ">r{i}").unwrap();
            writeln!(f, "{s}").unwrap();
        }
        p.to_str().unwrap().to_string()
    }

    /// Builds a reads set where two unitigs share a (k1-1)-mer junction and
    /// the bridge k-mer at k2 is present (or absent) in the reads.
    fn run(seqs: &[&str], ks: &[usize]) -> Vec<MultikUnitig> {
        let f = temp_fasta(seqs);
        let opts = MultikOptions {
            ks: ks.to_vec(),
            ..Default::default()
        };
        let out = assemble_multik(std::slice::from_ref(&f), &opts).unwrap();
        let _ = std::fs::remove_file(&f);
        out
    }

    #[test]
    fn linear_chain_compacts() {
        // A 21-mer genome duplicated at high coverage with a clean linear
        // junction; k=21 then k=31 should compact the two unitigs into one.
        let g = b"ACGTACGTACGTACGTACGTTTTTACGTACGTACGTACGTACGT".to_vec();
        let mut reads: Vec<String> = Vec::new();
        for _ in 0..20 {
            reads.push(String::from_utf8(g.clone()).unwrap());
        }
        let out = run(
            &reads.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
            &[21, 31],
        );
        assert!(!out.is_empty());
        assert!(out[0].bases.len() >= g.len() - 20);
    }

    #[test]
    fn unsupported_bridge_is_dropped() {
        // u1 and u2 share a (20)-mer but no read contains the k=31 bridge;
        // pass 0 still emits the link, the k=31 round drops it, so the two
        // unitigs are not merged.
        let s = b"ACGTACGTACGTACGTACGT".to_vec(); // 20-mer shared junction
        let left = format!("AAAA{}", String::from_utf8(s.clone()).unwrap());
        let right = format!("{}CCCC", String::from_utf8(s.clone()).unwrap());
        let mut reads: Vec<String> = Vec::new();
        for _ in 0..20 {
            reads.push(left.clone());
            reads.push(right.clone());
        }
        let out = run(
            &reads.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
            &[21, 31],
        );
        // The junction k=31 bridge AAAA..S..CCCC does not appear in reads,
        // so the chain must stay split (no merged unitig spanning both).
        let merged = out
            .iter()
            .any(|u| u.bases.len() > left.len().max(right.len()) + 5);
        assert!(!merged);
    }

    #[test]
    fn auto_ks_matches_read_length() {
        // Short reads (150 bp): k0 ~1/3 of the read length.
        assert_eq!(auto_ks(150), vec![50, 70, 90, 110]);
        // Short reads (108 bp): k0 = clamp(108/3, 31, 51) = 36.
        assert_eq!(auto_ks(108), vec![36, 56, 76]);
        // Long reads (>= 10 kb): capped at 51/81/111 (128 limit).
        assert_eq!(auto_ks(15000), vec![51, 81, 111]);
        // Zero/empty input yields no ks.
        assert!(auto_ks(0).is_empty());
    }

    #[test]
    fn read_n50_is_median_heavy() {
        let reads = vec![
            (b"A".repeat(100), Vec::new()),
            (b"T".repeat(300), Vec::new()),
            (b"C".repeat(200), Vec::new()),
        ];
        // Total 600 bp; N50 = length covering half = 300.
        assert_eq!(read_n50(&reads), 300);
    }

    fn mk_unitig(bases: &str, cov: f32) -> Unitig {
        Unitig {
            bases: bases.as_bytes().to_vec(),
            id: 0,
            coverage: cov,
            min_cov: 0,
            max_cov: 0,
            circular: false,
            abundances: Vec::new(),
        }
    }

    /// The bridge k-mer covering a u→v junction equals
    /// `upstream tail (k-1) + downstream continuation` (forward strand).
    #[test]
    fn bridge_kmer_forward_strand() {
        let s = "ACGTACGTACGTACGTACGT"; // shared 20-mer
        let u = mk_unitig(&format!("{}{}", "A".repeat(40), s), 30.0);
        let v = mk_unitig(&format!("{}{}", s, "C".repeat(40)), 30.0);
        let link = Link {
            to: 1,
            from_rc: false,
            to_rc: false,
        };
        let km = bridge_kmer(&u, &v, &link, 31, 21).expect("forward bridge");
        // Expected: u[30..] (A*10 + s) + v[20] (C).
        let exp = format!("{}C", String::from_utf8_lossy(&u.bases[30..]));
        let exp_km = kmer_from_bases(exp.as_bytes(), 31).unwrap();
        assert_eq!(km, exp_km);
    }

    /// Reverse-complemented target: u's tail matches rc(v's tail), so the
    /// partner joins reverse-complemented and its continuation is rc(v)[20]
    /// (rc of the C-run = G).
    #[test]
    fn bridge_kmer_reverse_target() {
        let s = "GGATCACAGTCTACACTGCT"; // shared 20-mer (non-palindromic)
        let s_rc = "AGCAGTGTAGACTGTGATCC"; // rc(s)
        let u = mk_unitig(&format!("{}{}", "A".repeat(40), s), 30.0);
        // v stored reverse-complemented: tail = rc(s).
        let v = mk_unitig(&format!("{}{}", "C".repeat(40), s_rc), 30.0);
        let link = Link {
            to: 1,
            from_rc: false,
            to_rc: true,
        };
        let km = bridge_kmer(&u, &v, &link, 31, 21).expect("rc bridge");
        // Expected: u[30..] + rc(v)[20] = (A*10 + s) + G.
        let exp = format!("{}G", String::from_utf8_lossy(&u.bases[30..]));
        let exp_km = kmer_from_bases(exp.as_bytes(), 31).unwrap();
        assert_eq!(km, exp_km);
    }

    /// from_rc=true: u's left end is the junction source, the partner sits
    /// upstream (v tail == u head), continuation is u's own base after the
    /// shared overlap.
    #[test]
    fn bridge_kmer_left_end_source() {
        let s = "ACGTACGTACGTACGTACGT"; // shared 20-mer
        let u = mk_unitig(&format!("{}{}", s, "A".repeat(40)), 30.0);
        let v = mk_unitig(&format!("{}{}", "C".repeat(40), s), 30.0);
        let link = Link {
            to: 1,
            from_rc: true,
            to_rc: false,
        };
        let km = bridge_kmer(&u, &v, &link, 31, 21).expect("left bridge");
        // Expected: v[30..] (C*10 + s) + u[20] (A).
        let exp = format!("{}A", String::from_utf8_lossy(&v.bases[30..]));
        let exp_km = kmer_from_bases(exp.as_bytes(), 31).unwrap();
        assert_eq!(km, exp_km);
    }

    /// Progressive filter keeps unique-chain unitigs regardless of their
    /// (low) coverage and only prunes isolated nodes by abundance.
    #[test]
    fn progressive_filter_keeps_main_path() {
        // A (cov 3) -> B (cov 5): both unique-chain, both must survive.
        let a = mk_unitig(&"A".repeat(60), 3.0);
        let b = mk_unitig(&"C".repeat(60), 5.0);
        let mut unitigs = vec![a, b];
        let mut links = vec![
            vec![Link {
                to: 1,
                from_rc: false,
                to_rc: false,
            }],
            Vec::new(),
        ];
        let mut branch = vec![false; unitigs.len()];
        let (dropped, _) = progressive_filter(&mut unitigs, &mut links, &mut branch, 21);
        assert_eq!(unitigs.len(), 2, "unique chain must survive the filter");
        assert!(dropped.is_empty());
    }

    /// Isolated low-abundance unitigs are pruned and collected as dropped.
    #[test]
    fn progressive_filter_prunes_isolated() {
        let a = mk_unitig(&"A".repeat(60), 30.0);
        let iso = mk_unitig(&"C".repeat(60), 2.0);
        let mut unitigs = vec![a, iso];
        let mut links = vec![Vec::new(), Vec::new()];
        let mut branch = vec![false; unitigs.len()];
        let (dropped, _) = progressive_filter(&mut unitigs, &mut links, &mut branch, 21);
        assert_eq!(
            unitigs.len(),
            1,
            "only the high-coverage isolated unitig survives"
        );
        assert_eq!(dropped.len(), 1);
        assert_eq!(dropped[0].bases.len(), 60);
    }

    /// Recompaction merges a unique chain into one longer unitig; stale
    /// links whose extremities were absorbed by the merge disappear (the
    /// recomputed links only connect real endpoint overlaps).
    #[test]
    fn recompact_merges_chain_and_drops_stale_links() {
        let s1 = "ACGTACGTACGTACGTACGT";
        let s2 = "GATCGATCGATCGATCGATC";
        let a = mk_unitig(&format!("{}{}", "A".repeat(80), s1), 30.0);
        let b = mk_unitig(&format!("{}{}{}", s1, "G".repeat(60), s2), 30.0);
        let c = mk_unitig(&format!("{}{}", s2, "C".repeat(80)), 30.0);
        let mut unitigs = vec![a, b, c];
        let mut links = vec![
            vec![Link {
                to: 1,
                from_rc: false,
                to_rc: false,
            }],
            vec![Link {
                to: 2,
                from_rc: false,
                to_rc: false,
            }],
            Vec::new(),
        ];
        let mut branch = vec![false; unitigs.len()];
        recompact_graph(&mut unitigs, &mut links, &mut branch, 21);
        // A(100) + B[20..](80) + C[20..](80) = 260.
        assert_eq!(unitigs.len(), 1);
        assert_eq!(unitigs[0].bases.len(), 260);
    }

    /// A convergence node (two predecessors) is not swallowed into either
    /// chain: strict end uniqueness, not link iteration order, decides chain
    /// membership (SKESA's "predecessor == 1" invariant).
    #[test]
    fn merge_chains_keeps_convergence_split() {
        let s1 = "ACGTACGTACGTACGTACGT";
        // u's end == v's begin == w's end: v has two predecessors.
        let u = mk_unitig(&format!("{}{}", "A".repeat(80), s1), 30.0);
        let v = mk_unitig(&format!("{}{}", s1, "G".repeat(80)), 30.0);
        let w = mk_unitig(&format!("{}{}", "T".repeat(80), s1), 30.0);
        let unitigs = vec![u, v, w];
        let links = vec![
            vec![Link {
                to: 1,
                from_rc: false,
                to_rc: false,
            }],
            Vec::new(),
            vec![Link {
                to: 1,
                from_rc: false,
                to_rc: false,
            }],
        ];
        let branch = vec![false; unitigs.len()];
        let out = merge_chains(&unitigs, &links, &branch, 21).unwrap();
        assert_eq!(out.len(), 3, "v must not merge with either predecessor");
        let total: usize = out.iter().map(|u| u.bases.len()).sum();
        assert_eq!(total, 300, "no unitig is fused at the junction");
        drop(links);
        drop(unitigs);
    }

    /// Recompaction applies the same strict end uniqueness: the convergence
    /// node survives as its own unitig instead of being merged into the
    /// first-visited chain.
    #[test]
    fn recompact_keeps_convergence_split() {
        let s1 = "ACGTACGTACGTACGTACGT";
        let u = mk_unitig(&format!("{}{}", "A".repeat(80), s1), 30.0);
        let v = mk_unitig(&format!("{}{}", s1, "G".repeat(80)), 30.0);
        let w = mk_unitig(&format!("{}{}", "T".repeat(80), s1), 30.0);
        let mut unitigs = vec![u, v, w];
        let mut links = vec![
            vec![Link {
                to: 1,
                from_rc: false,
                to_rc: false,
            }],
            Vec::new(),
            vec![Link {
                to: 1,
                from_rc: false,
                to_rc: false,
            }],
        ];
        let mut branch = vec![false; unitigs.len()];
        recompact_graph(&mut unitigs, &mut links, &mut branch, 21);
        // v keeps both incoming links (it was not absorbed), so the merged
        // graph still has three segments and v's begin retains two edges.
        assert_eq!(unitigs.len(), 3);
        let v_in: usize = links
            .iter()
            .map(|ls| ls.iter().filter(|l| l.to == 1 && !l.from_rc).count())
            .sum();
        assert_eq!(v_in, 2, "convergence node must keep both predecessors");
    }
}

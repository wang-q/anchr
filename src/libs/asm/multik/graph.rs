//! Unitig-graph operations shared by validation rounds and finalization:
//! megahit-style tip/weak-link/bubble cleanup, unsupported-window removal,
//! progressive abundance filtering, recompaction, and the final chain
//! merge.

use super::master::{CountView, RollCanon};
use super::MultikUnitig;
use crate::libs::asm::assemble::{compute_links, Link, Unitig};
use crate::libs::asm::table::{base_code, RefineTable};
use anyhow::Result;
use pgr::libs::nt::rev_comp;
use rayon::prelude::*;
use std::collections::HashSet;

/// Megahit-style tip removal: short unitigs (<= `max_tip_len`) that are
/// tips (one end has no connection) with depth far below their neighbour
/// (`neighbour > depth_ratio * self`) are error tips and dropped.
pub(crate) fn tip_remover(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    branch: &mut Vec<bool>,
    max_tip_len: usize,
    depth_ratio: f32,
) {
    // O(n + E) degree/depth precomputation (a per-unitig scan of all links
    // was O(n * E) on large final graphs). `in_src` counts distinct source
    // unitigs per target, matching the previous `any(l.to == i)` semantics.
    let n = unitigs.len();
    let mut in_src = vec![0usize; n];
    let mut in_cov = vec![0.0f32; n];
    let mut seen = vec![false; n];
    let mut touched: Vec<usize> = Vec::new();
    for (j, ls) in links.iter().enumerate() {
        touched.clear();
        for l in ls {
            if !seen[l.to] {
                seen[l.to] = true;
                touched.push(l.to);
                in_src[l.to] += 1;
                in_cov[l.to] = in_cov[l.to].max(unitigs[j].coverage);
            }
        }
        for &t in &touched {
            seen[t] = false;
        }
    }
    let keep: Vec<bool> = unitigs
        .iter()
        .enumerate()
        .map(|(i, u)| {
            if u.bases.len() > max_tip_len {
                return true;
            }
            let out = links[i].len();
            let in_deg = in_src[i];
            if out + in_deg == 0 {
                return true; // isolated: handled by the abundance filter
            }
            let is_tip = (out == 0 && in_deg >= 1) || (out >= 1 && in_deg == 0);
            if !is_tip {
                return true;
            }
            // Deepest neighbour depth.
            let mut max_neighbour = in_cov[i];
            for l in links[i].iter() {
                max_neighbour = max_neighbour.max(unitigs[l.to].coverage);
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
pub(crate) fn weak_link_remover(unitigs: &mut [Unitig], links: &mut [Vec<Link>], local_ratio: f32) {
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

/// Banded edit-distance similarity between two unitig sequences (megahit
/// `GetSimilarity`): `1 - edit_dist / max(n, m)` when the lengths are
/// within the band, else 0. The band is `max(n,m) * (1 - min_similarity)`,
/// so a returned value below `min_similarity` is impossible.
pub(crate) fn sequence_similarity(a: &[u8], b: &[u8], min_similarity: f64) -> f64 {
    let n = a.len();
    let m = b.len();
    let max_indel = (n.max(m) as f64 * (1.0 - min_similarity)) as usize;
    if n.abs_diff(m) > max_indel || max_indel < 1 {
        return 0.0;
    }
    // Row i-1 of the banded DP; dp[j] is the edit distance between a[..i]
    // and b[..j] for |i - j| <= max_indel (entries outside the band stay
    // MAX).
    let mut dp: Vec<usize> = vec![usize::MAX; m + 1];
    for (j, item) in dp.iter_mut().enumerate().take(m.min(max_indel) + 1) {
        *item = j;
    }
    for i in 1..=n {
        let mut ndp: Vec<usize> = vec![usize::MAX; m + 1];
        if i <= max_indel {
            ndp[0] = i;
        }
        let jmin = i.saturating_sub(max_indel).max(1);
        let jmax = (i + max_indel).min(m);
        for j in jmin..=jmax {
            let sub = dp[j - 1].saturating_add(usize::from(a[i - 1] != b[j - 1]));
            let del = dp[j].saturating_add(1);
            let ins = ndp[j - 1].saturating_add(1);
            ndp[j] = sub.min(del).min(ins);
        }
        dp = ndp;
    }
    if dp[m] == usize::MAX {
        0.0
    } else {
        1.0 - dp[m] as f64 / n.max(m) as f64
    }
}

/// Megahit-style bubble merge on the final (validated) unitig graph:
/// when several unitigs diverge from one oriented end and reconverge at
/// the same partner (each middle has exactly one in- and one out-link,
/// megahit `SearchAndPopBubble`), keep the highest-coverage path and drop
/// the alternatives — but only when every alternative is length-bounded
/// (`merge_len * k`) and edit-distance-similar to the main path (megahit
/// complex bubble: `>= merge_similar`). Dropped alternatives are returned
/// as independent unitigs so variant content is not lost; the surviving
/// chain is fused by the following recompaction.
///
/// Runs after all reads-based validation, so every junction already has
/// bridging reads; branch-marked (repeat-fragment) nodes never participate.
pub(crate) fn bubble_merge(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    branch: &mut Vec<bool>,
    k_build: usize,
    merge_similar: f64,
    merge_len: usize,
) -> Vec<Unitig> {
    let n = unitigs.len();
    if n <= 2 {
        return Vec::new();
    }
    let node_of = |id: usize, rev: bool| 2 * id + rev as usize;
    let max_len = ((merge_len * k_build) as f64 / merge_similar).round() as usize;

    // Directed segment graph (the same edges recompact/merge_chains use,
    // deduplicated: compute_links emits each junction from both ends).
    let mut segs: Vec<(usize, usize)> = Vec::new();
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
                segs.push((ln, rn));
            }
        }
    }
    let mut out_deg = vec![0usize; 2 * n];
    let mut in_deg = vec![0usize; 2 * n];
    let mut outs: Vec<Vec<usize>> = vec![Vec::new(); 2 * n];
    for &(ln, rn) in &segs {
        out_deg[ln] += 1;
        in_deg[rn] += 1;
        outs[ln].push(rn);
    }

    let mut deleted = vec![false; n];
    for s in 0..2 * n {
        if out_deg[s] <= 1 || deleted[s / 2] || branch[s / 2] {
            continue;
        }
        let mut middles: Vec<(usize, bool)> = Vec::new();
        let mut right: Option<usize> = None;
        let mut ok = true;
        for &m in &outs[s] {
            let mid = m / 2;
            let mrev = m % 2 == 1;
            if branch[mid] || in_deg[m] != 1 || out_deg[m] != 1 {
                ok = false;
                break;
            }
            if unitigs[mid].bases.len() > max_len {
                ok = false;
                break;
            }
            let rn = outs[m][0];
            if right.is_none() {
                right = Some(rn);
            } else if right != Some(rn) {
                ok = false;
                break;
            }
            middles.push((mid, mrev));
        }
        let Some(r) = right else { continue };
        if !ok
            || middles.len() < 2
            || in_deg[r] != middles.len()
            || deleted[r / 2]
            || branch[r / 2]
            || s / 2 == r / 2
        {
            continue;
        }
        if middles.iter().any(|&(mid, _)| deleted[mid]) {
            continue;
        }
        // Highest-coverage middle is the main path; every other middle must
        // be length- and sequence-similar to it.
        let mut order = middles;
        order.sort_by(|&(a, ar), &(b, br)| {
            unitigs[b]
                .coverage
                .partial_cmp(&unitigs[a].coverage)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then((a, ar).cmp(&(b, br)))
        });
        let (dom, dom_rev) = order[0];
        let dom_seq: Vec<u8> = if dom_rev {
            rev_comp(&unitigs[dom].bases).collect()
        } else {
            unitigs[dom].bases.clone()
        };
        let mut all_similar = true;
        for &(mid, mrev) in &order[1..] {
            let la = (unitigs[mid].bases.len() + k_build - 1) as f64;
            let lb = (dom_seq.len() + k_build - 1) as f64;
            if lb * merge_similar > la || la * merge_similar > lb {
                all_similar = false;
                break;
            }
            let seq: Vec<u8> = if mrev {
                rev_comp(&unitigs[mid].bases).collect()
            } else {
                unitigs[mid].bases.clone()
            };
            if sequence_similarity(&dom_seq, &seq, merge_similar) < merge_similar {
                all_similar = false;
                break;
            }
        }
        if !all_similar {
            continue;
        }
        for &(mid, _) in &order[1..] {
            deleted[mid] = true;
        }
    }

    let mut variants: Vec<Unitig> = Vec::new();
    for (i, &d) in deleted.iter().enumerate() {
        if d {
            let mut u = unitigs[i].clone();
            u.id = 0;
            variants.push(u);
        }
    }
    if !variants.is_empty() {
        let keep: Vec<bool> = deleted.iter().map(|&d| !d).collect();
        retain_graph(unitigs, links, branch, &keep);
        *links = compute_links(unitigs, k_build);
    }
    variants
}

/// Drops unitigs whose internal current-k k-mer is missing from the solid
/// table (chimeric cleanup). Unitigs shorter than `k` have no internal
/// k-mer to check and survive until the final compaction (their links were
/// already validated). Windows roll one base at a time (O(total bases), not
/// O(total bases × k)).
pub(crate) fn remove_unsupported(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    branch: &mut Vec<bool>,
    table: &impl CountView,
    k: usize,
    threshold: u32,
) {
    // Per-unitig work is independent (and the k0 master's 2000+ unitigs
    // otherwise serialize inside its round), so the keep scan joins the
    // ambient pool; when the pool is busy with other masters' rounds the
    // tasks simply run inline.
    let keep: Vec<bool> = unitigs
        .par_iter()
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
            // No consecutive-run cutting here: a run of >= 2 unsupported
            // windows also occurs at real low-coverage dips and pruned long
            // unitigs wholesale (MG1655 5-group chain N50 124K -> 79.6K,
            // 0 mis either way — pure over-pruning). Chimeric joins are
            // caught by the reads-bridge validation instead.
            let mut km = RollCanon::new(k, &u.bases);
            for j in 0..n_kmers {
                if j > 0 {
                    km.push_code(base_code(u.bases[j + k - 1]));
                }
                if !table.is_solid(km.canon(), threshold) {
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
pub(crate) fn progressive_filter(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    branch: &mut Vec<bool>,
    k_build: usize,
    repeat: Option<&RefineTable>,
    repeat_k: usize,
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
        recompact_graph(unitigs, links, branch, k_build, repeat, repeat_k);
        if std::env::var_os("ANCHR_MULTIK_TRACE_DIR").is_some() {
            let dir = std::env::var("ANCHR_MULTIK_TRACE_DIR").unwrap();
            use std::io::Write;
            let path = std::path::Path::new(&dir).join(format!("prog_t{:.1}.fa", t));
            if let Ok(mut f) = std::fs::File::create(&path) {
                for (i, u) in unitigs.iter().enumerate() {
                    writeln!(f, ">u{i} len={} cov={:.1}", u.bases.len(), u.coverage).unwrap();
                    let ls: Vec<String> = links[i]
                        .iter()
                        .map(|l| format!("{}{}", l.to, if l.from_rc { "-" } else { "+" }))
                        .collect();
                    writeln!(f, "#L {}", ls.join(" ")).unwrap();
                }
            }
        }
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

/// Per-unitig short-k window statistics over the reads-only repeat table:
/// `median` (the unitig's "normal" coverage) plus `left_hi`/`right_hi`
/// (the 95th-percentile count over the first/last `SPAN` bases). A
/// junction whose junction-adjacent end elevation far exceeds the unitig's
/// own median is a repeat bridge — the shared windows carry reads from both
/// genomic copies, so recompacting that chain would fuse two distant loci
/// into a relocation chimera (DH5alpha's two relocation misassemblies join
/// loci through ~1 kb repeats whose tandem copies carry SNP variations: at
/// the 130-mer probe the shared windows are broken into distinct alleles
/// and show no elevation, while a short k (31) sees the ~2x coverage). The
/// elevated windows sit a few hundred bases from the breakpoint, so the
/// single junction probe misses them too. A high quantile is used rather
/// than the max because the max over SPAN windows is noise-dominated at low
/// coverage (≈ μ+3.9√μ, crossing 1.5× the median below ~60×): a real
/// repeat bridge keeps many windows elevated, so the 95th percentile still
/// sees the ~2× level. Windows roll one base at a time; a unitig shorter
/// than one window has no stats (0: never a base).
struct ProbeStats {
    median: u32,
    left_hi: u32,
    right_hi: u32,
}

fn probe_stats(unitigs: &[Unitig], table: &RefineTable, window_len: usize) -> Vec<ProbeStats> {
    const SPAN: usize = 2000;
    const Q: f64 = 0.95;
    let hi = |v: &mut Vec<u32>| {
        v.sort_unstable();
        let i = ((v.len() - 1) as f64 * Q) as usize;
        v[i]
    };
    unitigs
        .iter()
        .map(|u| {
            let n = u.bases.len();
            if n < window_len {
                return ProbeStats {
                    median: 0,
                    left_hi: 0,
                    right_hi: 0,
                };
            }
            let mut counts: Vec<u32> = Vec::with_capacity(n - window_len + 1);
            let mut left: Vec<u32> = Vec::new();
            let mut right: Vec<u32> = Vec::new();
            let right_cut = n.saturating_sub(SPAN);
            let mut km = RollCanon::new(window_len, &u.bases);
            for p in 0..=n - window_len {
                if p > 0 {
                    km.push_code(base_code(u.bases[p + window_len - 1]));
                }
                let c = table.get_count_canonical(km.canon());
                if p < SPAN {
                    left.push(c);
                }
                if p + window_len > right_cut {
                    right.push(c);
                }
                counts.push(c);
            }
            counts.sort_unstable();
            ProbeStats {
                median: counts[counts.len() / 2],
                left_hi: hi(&mut left),
                right_hi: hi(&mut right),
            }
        })
        .collect()
}

/// Splits unitigs at internal low-coverage seams in the short-k repeat
/// table: a narrow run of windows whose count falls far below the unitig's
/// own median is a chimeric junction — the sequence on both sides is real
/// (each locus has reads support) but the crossing windows that join two
/// distant loci are absent from the reads. Real repeats show elevated (not
/// depleted) coverage and wide low regions are genuine genomic features
/// (e.g. coverage gaps), so only narrow internal runs are cut. End regions
/// are left alone (a unitig naturally ends at a coverage drop). Links are
/// recomputed from the new extremities; the fresh ends re-enter the chain
/// merge, whose [`is_repeat_bridge`] keeps repeat ends from re-fusing.
pub(crate) fn internal_repeat_bridge_split(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    branch: &mut Vec<bool>,
    table: &RefineTable,
    repeat_k: usize,
    k_build: usize,
) {
    if unitigs.is_empty() {
        return;
    }
    // A DH5alpha relocation junction is a ~30-60 bp run of near-zero 31-mer
    // windows inside the unitig (the two loci share no crossing k-mer);
    // real content has even coverage on both sides of the seam. Isolated
    // 1-2 window dips are coverage noise and wide low regions are real
    // features, so only narrow internal runs (MIN_RUN..=MAX_RUN) are cut.
    const SPAN: usize = 2000;
    const MIN_RUN: usize = 5;
    const MAX_RUN: usize = 200;
    const GAP: usize = 3;
    const LOW_RATIO: f64 = 0.3;
    const MIN_MEDIAN: u32 = 8;
    let mut out: Vec<Unitig> = Vec::new();
    let mut out_branch: Vec<bool> = Vec::new();
    for (u, &is_branch) in unitigs.iter().zip(branch.iter()) {
        let n = u.bases.len();
        // Need an internal region on both sides of any seam, plus a window.
        if n < 2 * SPAN + repeat_k {
            out.push(u.clone());
            out_branch.push(is_branch);
            continue;
        }
        let mut counts: Vec<u32> = Vec::with_capacity(n - repeat_k + 1);
        let mut km = RollCanon::new(repeat_k, &u.bases);
        for p in 0..=n - repeat_k {
            if p > 0 {
                km.push_code(base_code(u.bases[p + repeat_k - 1]));
            }
            counts.push(table.get_count_canonical(km.canon()));
        }
        let mut sorted = counts.clone();
        sorted.sort_unstable();
        let median = sorted[sorted.len() / 2];
        if median < MIN_MEDIAN {
            out.push(u.clone());
            out_branch.push(is_branch);
            continue;
        }
        // Collect runs of depleted windows across the whole unitig (end runs
        // participate in the width logic); a window at base `p` covers
        // [p, p+repeat_k).
        let low = (LOW_RATIO * median as f64) as u32;
        let mut runs: Vec<(usize, usize)> = Vec::new(); // (start, end) inclusive
        let mut i = 0usize;
        while i < counts.len() {
            if counts[i] < low {
                let start = i;
                let mut end = i;
                while i < counts.len() && counts[i] < low {
                    end = i;
                    i += 1;
                }
                // Merge with the previous run if within GAP.
                if let Some((_, prev_end)) = runs.last_mut() {
                    if start - *prev_end - 1 <= GAP {
                        *prev_end = end;
                        continue;
                    }
                }
                runs.push((start, end));
            } else {
                i += 1;
            }
        }
        // Cut only at narrow, internal runs: a real low-coverage region is
        // wide (> MAX_RUN) and an assembly end sits within SPAN of the end.
        let mut cut: Vec<usize> = Vec::new();
        for &(start, end) in &runs {
            if end - start + 1 < MIN_RUN || end - start + 1 > MAX_RUN {
                continue;
            }
            if start < SPAN || end > n - SPAN {
                continue;
            }
            cut.push((start + end) / 2);
        }
        if cut.is_empty() {
            out.push(u.clone());
            out_branch.push(is_branch);
            continue;
        }
        if std::env::var_os("ANCHR_MULTIK_DEBUG").is_some() && n > 50_000 {
            eprintln!(
                "internal_repeat_bridge_split cut len={n} median={median} cuts={} cut_pos={:?}",
                cut.len(),
                &cut[..cut.len().min(12)]
            );
        }
        // Split at the cut positions (same as split_by_bridge).
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
            if len < k_build {
                // Too short to host a (k_build-1)-mer end: drop the fragment.
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
    *links = compute_links(unitigs, k_build);
}

/// Whether the `i → l.to` junction spans a repeat bridge: the 95th-
/// percentile short-k window count over the junction-adjacent end of
/// either unitig exceeds `REPEAT_RATIO` times that unitig's own median
/// coverage. Such chains are not compacted (recompaction would fuse two
/// distant loci into a relocation chimera — the shared windows carry reads
/// from both copies).
fn is_repeat_bridge(i: usize, l: &Link, unitigs: &[Unitig], stats: &[ProbeStats]) -> bool {
    // A repeat region has reads on both copies; a low-coverage or short
    // unitig's coverage is too noisy to be a base. 1.5: a 2-copy repeat
    // raises the junction end to ~2x the unique median, but real ratios
    // land at 1.5-1.7x (coverage noise, read sharing), so 1.8 misses them.
    const REPEAT_RATIO: f32 = 1.5;
    let si = &stats[i];
    let sj = &stats[l.to];
    // Junction-adjacent ends: `i` leaves its right (left) end when
    // `from_rc` is false (true); the link enters `to`'s left (right) end
    // when `to_rc` is false (true).
    let (im, mi) = if l.from_rc {
        (si.left_hi, si.median)
    } else {
        (si.right_hi, si.median)
    };
    let (jm, mj) = if l.to_rc {
        (sj.right_hi, sj.median)
    } else {
        (sj.left_hi, sj.median)
    };
    let hi = im as f32 >= 2.0 && mi as f32 >= 2.0 && im as f32 > REPEAT_RATIO * mi as f32;
    let hj = jm as f32 >= 2.0 && mj as f32 >= 2.0 && jm as f32 > REPEAT_RATIO * mj as f32;
    let blocked = hi || hj;
    // Debug: report every junction between long-enough unitigs checked,
    // with the coverage ratio that decided it (to trace why a relocation
    // bridge survived).
    if std::env::var_os("ANCHR_MULTIK_DEBUG").is_some()
        && unitigs[i].bases.len() >= 500
        && unitigs[l.to].bases.len() >= 500
    {
        eprintln!(
            "repeat-check {}->{} len={}x{} hi={} ({}x{}) hj={} ({}x{}) blocked={blocked}",
            i,
            l.to,
            unitigs[i].bases.len(),
            unitigs[l.to].bases.len(),
            hi,
            im,
            mi,
            hj,
            jm,
            mj
        );
    }
    blocked
}

/// In-graph recompaction: merge unique chains into longer unitigs and
/// relink the chain endpoints (metaMDBG `UnitigGraph2::recompact`). The
/// merged unitig keeps the chain head's begin (from_rc=true) links and the
/// chain tail's end (from_rc=false) links; external edges pointing into
/// the chain are redirected to the merged unitig.
/// chain. `repeat` (reads-only short-k table) gates the chain ends: a
/// junction whose shared-overlap window count far exceeds both unitigs'
/// own median (a repeat bridge) is not compacted, so recompaction cannot
/// fuse two distant loci into a relocation chimera.
pub(crate) fn recompact_graph(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    branch: &mut Vec<bool>,
    k_build: usize,
    repeat: Option<&RefineTable>,
    repeat_k: usize,
) {
    let n = unitigs.len();
    if n <= 1 {
        return;
    }
    let stats = repeat.map(|t| probe_stats(unitigs, t, repeat_k));
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
            // A repeat bridge joins two distant loci; blocking the link here
            // keeps them as separate unitigs instead of fusing a chimera.
            if let (Some(_), Some(s)) = (repeat, stats.as_deref()) {
                if is_repeat_bridge(i, l, unitigs, s) {
                    if std::env::var_os("ANCHR_MULTIK_DEBUG").is_some() {
                        eprintln!(
                            "repeat-bridge blocked {}->{} len={}x{}",
                            i,
                            l.to,
                            unitigs[i].bases.len(),
                            unitigs[l.to].bases.len()
                        );
                    }
                    continue;
                }
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
/// chain-oriented ends (2n nodes: id * 2 + rev). `repeat` (reads-only
/// short-k table) gates the chain ends: a junction whose shared-overlap
/// window count far exceeds both unitigs' own median (a repeat bridge) is
/// not compacted, so a final chain cannot fuse two distant loci into a
/// relocation chimera.
pub(crate) fn merge_chains(
    unitigs: &[Unitig],
    links: &[Vec<Link>],
    branch: &[bool],
    k_build: usize,
    repeat: Option<&RefineTable>,
    repeat_k: usize,
) -> Result<Vec<MultikUnitig>> {
    let n = unitigs.len();
    let stats = repeat.map(|t| probe_stats(unitigs, t, repeat_k));
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
            // A repeat bridge joins two distant loci; blocking the link here
            // keeps them as separate unitigs instead of fusing a chimera.
            if let (Some(_), Some(s)) = (repeat, stats.as_deref()) {
                if is_repeat_bridge(i, l, unitigs, s) {
                    if std::env::var_os("ANCHR_MULTIK_DEBUG").is_some() {
                        eprintln!(
                            "repeat-bridge blocked {}->{} len={}x{}",
                            i,
                            l.to,
                            unitigs[i].bases.len(),
                            unitigs[l.to].bases.len()
                        );
                    }
                    continue;
                }
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

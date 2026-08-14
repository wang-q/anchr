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
use super::tadpole::{base_code, Kmer as TdKmer, TadpoleTable};
use anyhow::Result;
use pgr::libs::kmer::key::Kmer;
use pgr::libs::nt::rev_comp;

/// Options for [`assemble_multik`].
#[derive(Debug, Clone)]
pub struct MultikOptions {
    /// Increasing k-mer lengths; pass 0 assembles unitigs at `ks[0]`, each
    /// later k validates the previous unitig graph. Empty means auto-derive
    /// from the read-length N50 (see [`auto_ks`]).
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

    // Pass 0: maximal unitigs at the first k (BCALM graph3 semantics).
    let assemble_opts = AssembleOptions {
        k: ks[0],
        min_count_seed: opts.min_count_seed,
        parallel: opts.parallel,
        ..AssembleOptions::default()
    };
    let (mut unitigs, _) = assemble_unitigs_core(infiles, &assemble_opts)?;
    let mut links = compute_links(&unitigs, ks[0]);

    // Iterative rounds: validate the graph with each larger k.
    let mut low_abundance: Vec<MultikUnitig> = Vec::new();
    for &k in &ks[1..] {
        let table = count_at(&unitigs, infiles, k, opts.parallel)?;
        let threshold = opts.min_count_extend as u32;
        // 1. Cross-round link validation (solveEdges): the bridge k-mer
        // covering the junction must be solid at the current k.
        for (i, ls) in links.iter_mut().enumerate() {
            ls.retain(|l| {
                bridge_kmer(&unitigs[i], &unitigs[l.to], l, k, k0)
                    .is_some_and(|km| table.get_count(&km) >= threshold)
            });
        }
        // 2. Chimeric-unitig cleanup (removeUnsupportedUnitigs): every
        // internal current-k k-mer of a long-enough unitig must be solid.
        remove_unsupported(&mut unitigs, &mut links, &table, k, threshold)?;
        // 3. Progressive abundance filter (removeAbundanceNoQueue): drop
        // the lowest-abundance unitigs so the high-abundance path can
        // recompact; dropped unitigs stay independent output (mirroring
        // metaMDBG's cutoff snapshots).
        low_abundance.extend(progressive_filter(&mut unitigs, &mut links, k0));
    }

    // Final compaction: merge validated chains into long unitigs.
    let mut chains = merge_chains(&unitigs, &links, k0)?;
    chains.extend(low_abundance);
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
/// steps of `clamp(N50/100, 20, 30)`. Short reads (108 bp) yield
/// 21/41/61/81; long reads (>= 10 kb) cap at 31/61/91/121 — mirroring
/// metaMDBG's `computeLastK` (last k-min-mer spans ~2× N50).
fn auto_ks(n50: usize) -> Vec<usize> {
    if n50 == 0 {
        return Vec::new();
    }
    let k_max = (n50 * 8 / 10).clamp(31, Kmer::MAX_K);
    let k_min = (n50 / 10).clamp(21, 31);
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
) -> Result<TadpoleTable> {
    let mut reads = read_records(infiles)?;
    for u in unitigs {
        reads.push((u.bases.clone(), Vec::new()));
    }
    TadpoleTable::build_supermer(reads, k, None)
}

/// Encodes a base slice into the assembly k-mer key (canonical lookup is
/// applied by `TadpoleTable::get_count`).
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

/// Drops unitigs whose internal current-k k-mer is missing from the solid
/// table (chimeric cleanup). Unitigs shorter than `k` have no internal
/// k-mer to check and survive until the final compaction (their links were
/// already validated).
fn remove_unsupported(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    table: &TadpoleTable,
    k: usize,
    threshold: u32,
) -> Result<()> {
    let keep: Vec<bool> = unitigs
        .iter()
        .map(|u| {
            if u.bases.len() < k {
                return true;
            }
            (0..=u.bases.len() - k).all(|j| {
                kmer_from_bases(&u.bases[j..j + k], k)
                    .is_some_and(|km| table.get_count(&km) >= threshold)
            })
        })
        .collect();
    retain_graph(unitigs, links, &keep);
    Ok(())
}

/// Removes dropped unitigs and remaps surviving ids and link targets.
fn retain_graph(unitigs: &mut Vec<Unitig>, links: &mut Vec<Vec<Link>>, keep: &[bool]) {
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
    }
    *unitigs = kept_unitigs;
    *links = kept_links;
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
    k_build: usize,
) -> Vec<MultikUnitig> {
    let max_abundance: f32 = unitigs
        .iter()
        .map(|u| u.coverage)
        .fold(0.0f32, f32::max)
        .min(10000.0);
    let mut dropped: Vec<MultikUnitig> = Vec::new();
    let mut t = 1.1f32;
    while t < max_abundance && !unitigs.is_empty() {
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
                dropped.push(MultikUnitig {
                    bases: unitigs[i].bases.clone(),
                    coverage: unitigs[i].coverage,
                });
            }
        }
        if keep.iter().all(|&k| k) {
            t += (t * 0.1).min(10.0);
            continue; // nothing below the current cutoff, raise it
        }
        retain_graph(unitigs, links, &keep);
        // Recompact after removal so merged main-path unitigs inherit the
        // higher flank abundance (metaMDBG `recompact` in the same round).
        recompact_graph(unitigs, links, k_build);
        t += (t * 0.1).min(10.0);
    }
    dropped
}

/// In-graph recompaction: merge unique chains into longer unitigs and
/// relink the chain endpoints (metaMDBG `UnitigGraph2::recompact`). The
/// merged unitig keeps the chain head's begin (from_rc=true) links and the
/// chain tail's end (from_rc=false) links; external edges pointing into
/// the chain are redirected to the merged unitig.
fn recompact_graph(unitigs: &mut Vec<Unitig>, links: &mut Vec<Vec<Link>>, k_build: usize) {
    let n = unitigs.len();
    if n <= 1 {
        return;
    }
    let node_of = |id: usize, rev: bool| 2 * id + rev as usize;
    let mut right_of: Vec<Option<(usize, bool)>> = vec![None; 2 * n];
    let mut left_of: Vec<Option<(usize, bool)>> = vec![None; 2 * n];
    for (i, ls) in links.iter().enumerate() {
        for l in ls {
            let u = &unitigs[i];
            let v = &unitigs[l.to];
            if u.bases.len() < k_build - 1 || v.bases.len() < k_build - 1 {
                continue;
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
                    continue;
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
                    continue;
                }
            };
            let ln = node_of(li, lrev);
            let rn = node_of(ri, rrev);
            if right_of[ln].is_some() || left_of[rn].is_some() {
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
    // Head/tail of each merged chain, per new id (for link relinking).
    let mut heads: Vec<usize> = Vec::new();
    let mut tails: Vec<usize> = Vec::new();
    for &seed in &order {
        if visited[seed] {
            continue;
        }
        let mut head = seed;
        let mut head_rev = false;
        while let Some((p, prev)) = left_of[node_of(head, head_rev)] {
            if visited[p] {
                break;
            }
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
    k_build: usize,
) -> Result<Vec<MultikUnitig>> {
    let n = unitigs.len();
    let node_of = |id: usize, rev: bool| 2 * id + rev as usize;
    let mut right_of: Vec<Option<(usize, bool)>> = vec![None; 2 * n];
    let mut left_of: Vec<Option<(usize, bool)>> = vec![None; 2 * n];
    for (i, ls) in links.iter().enumerate() {
        for l in ls {
            let u = &unitigs[i];
            let v = &unitigs[l.to];
            // Unitigs shorter than the shared overlap cannot join a chain.
            if u.bases.len() < k_build - 1 || v.bases.len() < k_build - 1 {
                continue;
            }
            // Resolve the oriented segment by actual extremity matching.
            let (li, lrev, ri, rrev) = if !l.from_rc {
                let u_ext = &u.bases[u.bases.len() - (k_build - 1)..];
                let vb = &v.bases[..k_build - 1];
                let ve = &v.bases[v.bases.len() - (k_build - 1)..];
                if u_ext == vb {
                    (i, false, l.to, false)
                } else if u_ext == rev_comp(ve).collect::<Vec<u8>>().as_slice() {
                    (i, false, l.to, true)
                } else {
                    continue; // extremity mismatch: skip (conservative)
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
                    continue;
                }
            };
            let ln = node_of(li, lrev);
            let rn = node_of(ri, rrev);
            if right_of[ln].is_some() || left_of[rn].is_some() {
                continue; // already occupied: bubble, not a unique chain
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
        while let Some((p, prev)) = left_of[node_of(head, head_rev)] {
            if visited[p] {
                break;
            }
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
        // Short reads (108 bp): 21/41/61/81.
        assert_eq!(auto_ks(108), vec![21, 41, 61, 81]);
        // Long reads (>= 10 kb): capped at 31/61/91/121 (128 limit).
        assert_eq!(auto_ks(15000), vec![31, 61, 91, 121]);
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
        let dropped = progressive_filter(&mut unitigs, &mut links, 21);
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
        let dropped = progressive_filter(&mut unitigs, &mut links, 21);
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
        let a = mk_unitig(&format!("{}{}", "A".repeat(40), s1), 30.0);
        let b = mk_unitig(&format!("{}{}{}", s1, "G".repeat(40), s2), 30.0);
        let c = mk_unitig(&format!("{}{}", s2, "C".repeat(40)), 30.0);
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
        recompact_graph(&mut unitigs, &mut links, 21);
        // A(60) + B[20..](60) + C[20..](40) = 160.
        assert_eq!(unitigs.len(), 1);
        assert_eq!(unitigs[0].bases.len(), 160);
    }
}

//! Bridge and probe validation: the junction bridge k-mer, reads-bridge
//! probe filters, and chimeric-junction splitting.

use super::master::RollCanon;
use crate::libs::asm::assemble::{compute_links, Link, Unitig};
use crate::libs::asm::table::{base_code, Kmer as TdKmer, RefineTable};
use pgr::libs::nt::rev_comp;

/// Encodes a base slice into the assembly k-mer key (canonical lookup is
/// applied by `RefineTable::get_count`).
pub(crate) fn kmer_from_bases(bases: &[u8], k: usize) -> Option<TdKmer> {
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
///
/// Only the required extremity windows are copied (never the whole partner
/// sequence): compacted unitigs grow long, and this runs per link per round.
pub(crate) fn bridge_kmer(
    u: &Unitig,
    v: &Unitig,
    link: &Link,
    k: usize,
    k_build: usize,
) -> Option<TdKmer> {
    // Unitigs shorter than the previous k cannot provide the shared
    // (k_prev-1)-mer extremity plus a continuation base; treat their links
    // as unsupported (they still survive as short unitigs).
    if u.bases.len() < k_build || v.bases.len() < k_build {
        return None;
    }
    let km1 = k - 1;
    let kb1 = k_build - 1;
    let v_begin = &v.bases[..kb1];
    let v_end = &v.bases[v.bases.len() - kb1..];
    let mut seq: Vec<u8> = Vec::with_capacity(k);
    if link.from_rc {
        // u's left end is the junction source: partner upstream, u downstream.
        let u_ext = &u.bases[..kb1];
        if v.bases.len() < km1 {
            return None; // partner cannot fill the (k-1) tail window
        }
        if u_ext == v_end {
            seq.extend_from_slice(&v.bases[v.bases.len() - km1..]);
        } else if u_ext == rev_comp(v_begin).collect::<Vec<u8>>().as_slice() {
            // Last (k-1) bases of rc(v) = rc of v's first (k-1) bases.
            seq.extend(rev_comp(&v.bases[..km1]));
        } else {
            return None;
        }
        seq.push(u.bases[kb1]);
    } else {
        // u's right end is the junction source: u upstream, partner downstream.
        let u_ext = &u.bases[u.bases.len() - kb1..];
        if u.bases.len() < km1 {
            return None; // u cannot fill the (k-1) tail window
        }
        seq.extend_from_slice(&u.bases[u.bases.len() - km1..]);
        if u_ext == v_begin {
            seq.push(v.bases[kb1]);
        } else if u_ext == rev_comp(v_end).collect::<Vec<u8>>().as_slice() {
            // rc(v)[k_build-1] = complement of v[len-1-(k_build-1)].
            let l = v.bases.len();
            seq.push(rev_comp(&v.bases[l - 1 - kb1..l - kb1]).next()?);
        } else {
            return None;
        }
    }
    kmer_from_bases(&seq, k)
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
/// no reads covering the junction and is pruned). `table` is the reads-only
/// probe table shared across rounds. Links whose probe cannot be built
/// (short unitigs) are kept conservatively.
pub(crate) fn bridge_filter(
    unitigs: &[Unitig],
    links: &mut [Vec<Link>],
    table: &RefineTable,
    k0: usize,
    probe_half: usize,
    threshold: u32,
) {
    if unitigs.is_empty() {
        return;
    }
    for (i, ls) in links.iter_mut().enumerate() {
        ls.retain(|l| {
            probe_kmer(&unitigs[i], &unitigs[l.to], l, k0, probe_half)
                .map(|p| table.get_count(&p) >= threshold)
                .unwrap_or(false)
        });
    }
}

/// Splits unitigs at internal windows that are not supported by any read:
/// every `2*probe_half`-mer window of a unitig must occur in the reads (the
/// unitig's own sequence comes from reads, so a window with count 0 is a
/// chimeric junction — the abundance recompaction fused two distant regions).
/// Splitting keeps those junctions out of the final compaction. Links are
/// recomputed from the new extremities. `table` is the reads-only probe
/// table shared across rounds; windows roll one base at a time.
pub(crate) fn split_by_bridge(
    unitigs: &mut Vec<Unitig>,
    links: &mut Vec<Vec<Link>>,
    branch: &mut Vec<bool>,
    table: &RefineTable,
    k0: usize,
    probe_half: usize,
    threshold: u32,
) {
    if unitigs.is_empty() {
        return;
    }
    let probe_len = probe_half * 2;
    let mut out: Vec<Unitig> = Vec::new();
    let mut out_branch: Vec<bool> = Vec::new();
    for (u, &is_branch) in unitigs.iter().zip(branch.iter()) {
        let n = u.bases.len();
        if n < probe_len {
            out.push(u.clone());
            out_branch.push(is_branch);
            continue;
        }
        // Mark windows without read support (rolling window: one push_right
        // and one lookup per position instead of an O(probe_len) encode).
        let mut cut: Vec<usize> = Vec::new();
        let mut km = RollCanon::new(probe_len, &u.bases);
        let mut ok = table.get_count_canonical(km.canon()) >= threshold;
        let mut prev_cut = !ok;
        if !ok {
            cut.push(0);
        }
        for i in 1..=n - probe_len {
            km.push_code(base_code(u.bases[i + probe_len - 1]));
            ok = table.get_count_canonical(km.canon()) >= threshold;
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
}

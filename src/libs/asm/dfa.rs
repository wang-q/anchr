//! Cuttlefish-style DFA state classification for unitig building.
//!
//! Every solid canonical k-mer gets a 4-field state (in/out degree and the
//! unique base on each side, when the degree is 1). The unitig walk then
//! needs only O(1) state lookups instead of re-scanning the four left/right
//! extension buckets at every step. The classification pass is read-only on
//! the count table and parallelizes by partitioning the sorted vertex list
//! (no CAS / locks needed), unlike cuttlefish's edge-scan update path.

use super::table::{Kmer, RefineTable};
use rayon::prelude::*;

/// Per-vertex state: solidity-adjacent in/out degrees (0/1/2+ encoded as
/// `u8` count, with 2 meaning "2 or more") plus the unique extension base
/// when a side has exactly one edge.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct VertexState {
    /// Number of distinct predecessor bases (0/1/2+).
    pub in_count: u8,
    /// Number of distinct successor bases (0/1/2+).
    pub out_count: u8,
    /// Unique predecessor base (valid iff `in_count == 1`).
    pub in_base: u8,
    /// Unique successor base (valid iff `out_count == 1`).
    pub out_base: u8,
}

/// Sentinel for "no unique continuation on this side" in the successor
/// index arrays.
const NO_SUCC: u32 = u32::MAX;

/// Classified states for every vertex in `RefineTable::sorted_entries`.
/// The unique-continuation successor of each vertex is stored as an entry
/// index (captured from the classification probes' table rows), so the
/// unitig walk steps through plain array reads — no hash index is kept.
pub(crate) struct VertexStates {
    /// Sorted solid (canonical k-mer, count) entries; kept so the walk
    /// reuses the classification pass instead of rebuilding `solid_entries`.
    entries: Vec<(Kmer, u32)>,
    states: Vec<VertexState>,
    /// Canonical counts parallel to `states` (per-unitig coverage; avoids a
    /// second `get_count` canonical+bsearch per k-mer during the walk).
    counts: Vec<u32>,
    /// Entry index of the unique solid successor (`km + out_base`), when
    /// `out_count == 1` (`NO_SUCC` otherwise).
    succ_out: Vec<u32>,
    /// Entry index of the unique solid predecessor (`in_base + km`), when
    /// `in_count == 1` (`NO_SUCC` otherwise).
    succ_in: Vec<u32>,
}

impl VertexStates {
    /// Classifies all solid k-mers (count >= `threshold`). `threads > 1`
    /// parallelizes the per-vertex pass; the result is order-independent.
    pub(crate) fn classify(table: &RefineTable, threshold: u32, threads: usize) -> Self {
        let threads = if threads == 0 {
            rayon::current_num_threads()
        } else {
            threads
        };
        let t0 = std::time::Instant::now();
        let entries = table.solid_entries(threshold);
        let mut states = vec![VertexState::default(); entries.len()];
        let mut counts = vec![0u32; entries.len()];
        let mut succ_out = vec![NO_SUCC; entries.len()];
        let mut succ_in = vec![NO_SUCC; entries.len()];
        // Row -> entry rank for solid rows (u32::MAX elsewhere): lets the
        // probes translate a found table row into the parallel entry index.
        let ranks = table.solid_row_ranks(threshold);

        let fill = |states: &mut [VertexState],
                    counts: &mut [u32],
                    succ_out: &mut [u32],
                    succ_in: &mut [u32]| {
            states
                .par_iter_mut()
                .zip(counts.par_iter_mut())
                .zip(succ_out.par_iter_mut())
                .zip(succ_in.par_iter_mut())
                .enumerate()
                .for_each(|(i, (((st, cnt), so), si))| {
                    let (km, count) = entries[i];
                    if count < threshold {
                        return;
                    }
                    *cnt = count;
                    let (in_count, in_base, in_row) = count_in(table, &km, threshold);
                    let (out_count, out_base, out_row) = count_out(table, &km, threshold);
                    *st = VertexState {
                        in_count,
                        out_count,
                        in_base,
                        out_base,
                    };
                    if in_count == 1 {
                        *si = ranks[in_row.unwrap()];
                    }
                    if out_count == 1 {
                        *so = ranks[out_row.unwrap()];
                    }
                });
        };

        if threads > 1 {
            // Ambient pool: the command wraps the whole assemble call in a
            // single rayon pool of `--parallel` threads, so classification
            // must not create a second pool (thread oversubscription).
            fill(&mut states, &mut counts, &mut succ_out, &mut succ_in);
        } else {
            states
                .iter_mut()
                .zip(counts.iter_mut())
                .zip(succ_out.iter_mut())
                .zip(succ_in.iter_mut())
                .enumerate()
                .for_each(|(i, (((st, cnt), so), si))| {
                    let (km, count) = entries[i];
                    if count < threshold {
                        return;
                    }
                    *cnt = count;
                    let (in_count, in_base, in_row) = count_in(table, &km, threshold);
                    let (out_count, out_base, out_row) = count_out(table, &km, threshold);
                    *st = VertexState {
                        in_count,
                        out_count,
                        in_base,
                        out_base,
                    };
                    if in_count == 1 {
                        *si = ranks[in_row.unwrap()];
                    }
                    if out_count == 1 {
                        *so = ranks[out_row.unwrap()];
                    }
                });
        }

        if std::env::var_os("ANCHR_DFA_TIMING").is_some() {
            eprintln!(
                "dfa classify: {:.3}s (threads={})",
                t0.elapsed().as_secs_f64(),
                threads
            );
        }
        VertexStates {
            entries,
            states,
            counts,
            succ_out,
            succ_in,
        }
    }

    /// Sorted solid entries (deterministic seed scan order).
    pub(crate) fn entries(&self) -> &[(Kmer, u32)] {
        &self.entries
    }

    /// Entry index of a rolling window pair (`fw` = walked strand, `rc`
    /// its reverse complement held in lockstep): binary search in the
    /// sorted entries (a handful of calls per unitig seed).
    pub(crate) fn canon_idx_pair(&self, fw: &Kmer, rc: &Kmer) -> Option<usize> {
        let c = if fw.cmp_bases(rc) != std::cmp::Ordering::Greater {
            fw
        } else {
            rc
        };
        self.entries
            .binary_search_by(|(km, _)| km.cmp_bases(c))
            .ok()
    }

    /// Entry index of the unique solid continuation from the oriented
    /// window at entry `idx` (`fw_is_canon`: the walked strand is the
    /// canonical one), plus its extension base; `None` at a branch or
    /// dead end. Plain array reads — the walk performs no k-mer hashing.
    pub(crate) fn step(&self, fw_is_canon: bool, idx: usize) -> Option<(u8, usize)> {
        let st = &self.states[idx];
        let succ = if fw_is_canon {
            st.out_base
        } else {
            3 - st.in_base
        };
        let si = if fw_is_canon {
            self.succ_out[idx]
        } else {
            self.succ_in[idx]
        };
        (si != NO_SUCC).then_some((succ, si as usize))
    }

    /// Oriented in-count of the window at entry `idx` (`fw_is_canon`: the
    /// window's strand as walked is the canonical one).
    pub(crate) fn in_count_at(&self, fw_is_canon: bool, idx: usize) -> u8 {
        let st = &self.states[idx];
        if fw_is_canon {
            st.in_count
        } else {
            st.out_count
        }
    }

    /// Count at a known entry index.
    pub(crate) fn count_at(&self, idx: usize) -> u32 {
        self.counts[idx]
    }
}

/// Distinct solid predecessor bases (`b + kmer[..k-1]`), the unique one,
/// and its table row (set iff exactly one is solid).
fn count_in(table: &RefineTable, km: &Kmer, threshold: u32) -> (u8, u8, Option<usize>) {
    let mut n = 0u8;
    let mut base = 0u8;
    let mut row = None;
    for b in 0..4u8 {
        let mut q = *km;
        q.push_left(b);
        if let Some((r, c)) = table.find_row(&q) {
            if c >= threshold {
                n += 1;
                base = b;
                row = Some(r);
            }
        }
    }
    (n, base, row)
}

/// Distinct solid successor bases (`kmer[1..] + b`), the unique one, and
/// its table row (set iff exactly one is solid).
fn count_out(table: &RefineTable, km: &Kmer, threshold: u32) -> (u8, u8, Option<usize>) {
    let mut n = 0u8;
    let mut base = 0u8;
    let mut row = None;
    for b in 0..4u8 {
        let mut q = *km;
        q.push_right(b);
        if let Some((r, c)) = table.find_row(&q) {
            if c >= threshold {
                n += 1;
                base = b;
                row = Some(r);
            }
        }
    }
    (n, base, row)
}

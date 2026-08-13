//! Cuttlefish-style DFA state classification for unitig building.
//!
//! Every solid canonical k-mer gets a 4-field state (in/out degree and the
//! unique base on each side, when the degree is 1). The unitig walk then
//! needs only O(1) state lookups instead of re-scanning the four left/right
//! extension buckets at every step. The classification pass is read-only on
//! the count table and parallelizes by partitioning the sorted vertex list
//! (no CAS / locks needed), unlike cuttlefish's edge-scan update path.

use super::tadpole::{Kmer, KmerFnvHasher, TadpoleTable};
use rayon::prelude::*;
use std::collections::HashMap;
use std::hash::BuildHasherDefault;

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

/// Classified states for every vertex in `TadpoleTable::sorted_entries`,
/// with a canonical-key index for the walk.
pub(crate) struct VertexStates {
    /// Sorted solid (canonical k-mer, count) entries; kept so the walk
    /// reuses the classification pass instead of rebuilding `solid_entries`.
    entries: Vec<(Kmer, u32)>,
    states: Vec<VertexState>,
    /// Canonical counts parallel to `states` (per-unitig coverage; avoids a
    /// second `get_count` canonical+bsearch per k-mer during the walk).
    counts: Vec<u32>,
    index: HashMap<Kmer, usize, BuildHasherDefault<KmerFnvHasher>>,
}

impl VertexStates {
    /// Classifies all solid k-mers (count >= `threshold`). `threads > 1`
    /// parallelizes the per-vertex pass; the result is order-independent.
    pub(crate) fn classify(table: &TadpoleTable, threshold: u32, threads: usize) -> Self {
        let threads = if threads == 0 {
            rayon::current_num_threads()
        } else {
            threads
        };
        let t0 = std::time::Instant::now();
        let entries = table.solid_entries(threshold);
        let mut states = vec![VertexState::default(); entries.len()];
        let mut counts = vec![0u32; entries.len()];

        let fill = |states: &mut [VertexState], counts: &mut [u32]| {
            states
                .par_iter_mut()
                .zip(counts.par_iter_mut())
                .enumerate()
                .for_each(|(i, (st, cnt))| {
                    let (km, count) = entries[i];
                    if count < threshold {
                        return;
                    }
                    *cnt = count;
                    let (in_count, in_base) = count_in(table, &km, threshold);
                    let (out_count, out_base) = count_out(table, &km, threshold);
                    *st = VertexState {
                        in_count,
                        out_count,
                        in_base,
                        out_base,
                    };
                });
        };

        if threads > 1 {
            // Ambient pool: the command wraps the whole assemble call in a
            // single rayon pool of `--parallel` threads, so classification
            // must not create a second pool (thread oversubscription).
            fill(&mut states, &mut counts);
        } else {
            states
                .iter_mut()
                .zip(counts.iter_mut())
                .enumerate()
                .for_each(|(i, (st, cnt))| {
                    let (km, count) = entries[i];
                    if count < threshold {
                        return;
                    }
                    *cnt = count;
                    let (in_count, in_base) = count_in(table, &km, threshold);
                    let (out_count, out_base) = count_out(table, &km, threshold);
                    *st = VertexState {
                        in_count,
                        out_count,
                        in_base,
                        out_base,
                    };
                });
        }

        let mut index =
            HashMap::with_capacity_and_hasher(entries.len(), BuildHasherDefault::default());
        for (i, (km, _)) in entries.iter().enumerate() {
            index.insert(*km, i);
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
            index,
        }
    }

    /// Sorted solid entries (deterministic seed scan order).
    pub(crate) fn entries(&self) -> &[(Kmer, u32)] {
        &self.entries
    }

    /// Index of a canonical k-mer in the parallel `states`/`counts` arrays.
    pub(crate) fn idx(&self, canon: &Kmer) -> Option<usize> {
        self.index.get(canon).copied()
    }

    /// Exact canonical count of a vertex already in canonical form (0 for
    /// absent). O(1): one index lookup + one array read.
    pub(crate) fn count_canonical(&self, canon: &Kmer) -> u32 {
        self.index.get(canon).map(|&i| self.counts[i]).unwrap_or(0)
    }

    /// Unique successor base of the *oriented* `kmer` (forward or RC), or
    /// `None` when its exiting side has 0 or >= 2 edges.
    pub(crate) fn out_base(&self, kmer: &Kmer) -> Option<u8> {
        let canon = kmer.canonical();
        let &idx = self.index.get(&canon)?;
        let st = &self.states[idx];
        if kmer.cmp_bases(&canon) == std::cmp::Ordering::Equal {
            (st.out_count == 1).then_some(st.out_base)
        } else {
            (st.in_count == 1).then_some(3 - st.in_base)
        }
    }

    /// Number of distinct predecessor bases of the *oriented* `kmer`
    /// (0/1/2+), mirroring the old `unique_solid_in` count.
    pub(crate) fn in_count(&self, kmer: &Kmer) -> u8 {
        let canon = kmer.canonical();
        let Some(&idx) = self.index.get(&canon) else {
            return 0;
        };
        let st = &self.states[idx];
        if kmer.cmp_bases(&canon) == std::cmp::Ordering::Equal {
            st.in_count
        } else {
            st.out_count
        }
    }
}

/// Distinct predecessor bases (`b + kmer[..k-1]` solid) and the unique one.
fn count_in(table: &TadpoleTable, km: &Kmer, threshold: u32) -> (u8, u8) {
    let mut n = 0u8;
    let mut base = 0u8;
    for b in 0..4u8 {
        let mut q = *km;
        q.push_left(b);
        if table.get_count(&q) >= threshold {
            n += 1;
            base = b;
        }
    }
    (n, base)
}

/// Distinct successor bases (`kmer[1..] + b` solid) and the unique one.
fn count_out(table: &TadpoleTable, km: &Kmer, threshold: u32) -> (u8, u8) {
    let mut n = 0u8;
    let mut base = 0u8;
    for b in 0..4u8 {
        let mut q = *km;
        q.push_right(b);
        if table.get_count(&q) >= threshold {
            n += 1;
            base = b;
        }
    }
    (n, base)
}

//! Layout stitching into consensus contigs (OLC stage S3).
//!
//! Overlaps are exact, so consensus is an exact stitch: walk each layout in
//! order, orient every unitig by its strand, append only the bases beyond
//! the exact overlap with the previous step, and verify the overlapping
//! prefix matches the already-stitched suffix.

use super::layout::Layout;
use super::overlap::Unitig;
use anyhow::Result;
use pgr::libs::nt::rev_comp;

/// One consensus contig.
#[derive(Debug, Clone, PartialEq)]
pub struct Contig {
    /// Consensus sequence (5' -> 3' as laid out).
    pub seq: Vec<u8>,
    /// Approximate unitig depth (`sum(unitig lengths) / contig length`).
    pub coverage: f64,
}

/// `consensus` with an approximate-containment dedup: contigs whose sequence
/// is covered by `ratio` (fraction of the shorter contig, either strand) of a
/// longer kept contig are dropped. `ratio < 1.0` merges near-duplicate
/// contigs whose boundaries differ by a few bases (multi-coverage-set
/// assembly), keeping the longest representative without losing the covered
/// sequence; `ratio = 1.0` is the exact-substring behaviour.
pub fn consensus_with_ratio(
    unitigs: &[Unitig],
    layouts: &[Layout],
    min_contig_len: usize,
    ratio: f64,
) -> Result<Vec<Contig>> {
    let mut contigs = Vec::new();
    for (ci, layout) in layouts.iter().enumerate() {
        let mut seq: Vec<u8> = Vec::new();
        let mut total = 0usize;
        for (si, step) in layout.steps.iter().enumerate() {
            let mut piece: Vec<u8> = if step.strand == '+' {
                unitigs[step.unitig].seq.clone()
            } else {
                rev_comp(&unitigs[step.unitig].seq).collect()
            };
            total += piece.len();
            if si == 0 {
                seq.append(&mut piece);
                continue;
            }
            let overlap = step.overlap_len;
            anyhow::ensure!(
                overlap <= seq.len() && overlap <= piece.len(),
                "layout contig_{} step {si}: overlap {overlap} exceeds step lengths",
                ci + 1
            );
            let start = seq.len() - overlap;
            anyhow::ensure!(
                seq[start..] == piece[..overlap],
                "layout contig_{} step {si}: overlapping bases disagree \
                 (exact overlaps must match)",
                ci + 1
            );
            seq.extend_from_slice(&piece[overlap..]);
        }
        if seq.len() >= min_contig_len {
            let coverage = total as f64 / seq.len() as f64;
            contigs.push(Contig { seq, coverage });
        }
    }
    let contigs = dedup_contained_ratio(contigs, ratio);
    if ratio < 1.0 {
        // Cross-group anchor sets re-cover the same locus with different
        // boundaries; stich boundary-differing near-duplicates into one
        // contig (exact overlap detection cannot join them).
        Ok(merge_overlapping_contigs(contigs, MIN_OVERLAP))
    } else {
        Ok(contigs)
    }
}

/// Minimum approximate overlap (bases) for stiching two contigs.
const MIN_OVERLAP: usize = 5000;

/// Seed length for locating approximate overlaps between contigs.
const SEED_LEN: usize = 31;

/// Merges boundary-differing near-duplicate contigs: a shorter contig whose
/// head aligns inside a longer one and whose tail extends past its end (or
/// the symmetric head-before case, either strand) is stitched into one
/// contig. Only single dominant high-identity alignments are merged, so
/// chimeric contigs with multi-block alignments are left untouched.
fn merge_overlapping_contigs(mut contigs: Vec<Contig>, min_overlap: usize) -> Vec<Contig> {
    contigs.sort_by_key(|c| std::cmp::Reverse(c.seq.len()));
    let mut kept: Vec<Contig> = Vec::with_capacity(contigs.len());
    for c in contigs {
        let mut merged = false;
        for k in kept.iter_mut() {
            if let Some(seq) = try_merge(&k.seq, &c.seq, min_overlap) {
                k.seq = seq;
                merged = true;
                break;
            }
        }
        if !merged {
            kept.push(c);
        }
    }
    kept
}

/// Stitches `cand` into `keeper` when they are the same locus with different
/// boundaries; returns `None` when no single dominant high-identity overlap
/// exists. Either strand is handled (the result follows `keeper`'s strand).
fn try_merge(keeper: &[u8], cand: &[u8], min_overlap: usize) -> Option<Vec<u8>> {
    if keeper.len() < 2 * SEED_LEN || cand.len() < 2 * SEED_LEN {
        return None;
    }
    let mut index: std::collections::HashMap<&[u8], Vec<usize>> = std::collections::HashMap::new();
    for (i, w) in keeper.windows(SEED_LEN).enumerate() {
        index.entry(w).or_default().push(i);
    }
    let rc: Vec<u8> = rev_comp(cand).collect();
    let direct = overlap_geometry(&index, cand);
    if let Some(seq) = merge_geometry(keeper, cand, false, direct, min_overlap) {
        return Some(seq);
    }
    let rev = overlap_geometry(&index, &rc);
    merge_geometry(keeper, cand, true, rev, min_overlap)
}

/// Dominant k-mer offset of `query` inside `keeper` (position in keeper minus
/// position in query), the number of query k-mers supporting it, and the
/// number of query k-mers present anywhere in `keeper`.
fn overlap_geometry(
    index: &std::collections::HashMap<&[u8], Vec<usize>>,
    query: &[u8],
) -> (isize, usize, usize) {
    let mut hist: std::collections::HashMap<isize, usize> = std::collections::HashMap::new();
    let mut matched = 0usize;
    for (i, w) in query.windows(SEED_LEN).enumerate() {
        if let Some(ps) = index.get(w) {
            matched += 1;
            for &p in ps {
                *hist.entry(p as isize - i as isize).or_default() += 1;
            }
        }
    }
    let (offset, hits) = hist.into_iter().max_by_key(|(_, n)| *n).unwrap_or((0, 0));
    (offset, hits, matched)
}

/// Stitches when `query` (either `cand` or its reverse complement) starts
/// inside `keeper` and extends past its end, or starts before it and shares
/// a high-identity tail overlap with `keeper`'s head.
#[allow(clippy::too_many_arguments)]
fn merge_geometry(
    keeper: &[u8],
    cand: &[u8],
    query_is_rc: bool,
    geometry: (isize, usize, usize),
    min_overlap: usize,
) -> Option<Vec<u8>> {
    let (off, hits, matched) = geometry;
    if hits < 100 || hits * 10 < matched * 8 {
        // Too short, or the query aligns to several loci (chimeric).
        return None;
    }
    let query: Vec<u8> = if query_is_rc {
        rev_comp(cand).collect()
    } else {
        cand.to_vec()
    };
    if off >= 0 {
        let off = off as usize;
        if off >= keeper.len() {
            return None;
        }
        let ov = (keeper.len() - off).min(query.len());
        if ov < min_overlap || keeper.get(off..off + SEED_LEN) != Some(&query[..SEED_LEN]) {
            // The query head must anchor inside the keeper; a chimeric head
            // aligned elsewhere drops the identity below the threshold.
            return None;
        }
        if identity(&keeper[off..off + ov], &query[..ov]) < 0.99 {
            return None;
        }
        if off + query.len() <= keeper.len() {
            // Fully contained; the dedup step above should already have
            // dropped it, but keep the longer representative regardless.
            return Some(keeper.to_vec());
        }
        let mut seq = keeper.to_vec();
        seq.extend_from_slice(&query[ov..]);
        Some(seq)
    } else {
        let head = (-off) as usize;
        if head >= query.len() {
            return None;
        }
        let ov = (query.len() - head).min(keeper.len());
        if ov < min_overlap || keeper.get(..SEED_LEN) != query.get(head..head + SEED_LEN) {
            return None;
        }
        if identity(&keeper[..ov], &query[head..head + ov]) < 0.99 {
            return None;
        }
        let mut seq = query;
        seq.extend_from_slice(&keeper[ov..]);
        Some(seq)
    }
}

/// Identity of two equally long slices.
fn identity(a: &[u8], b: &[u8]) -> f64 {
    let mm = a.iter().zip(b).filter(|(x, y)| x != y).count();
    1.0 - mm as f64 / a.len() as f64
}

/// Drops contigs whose sequence is covered by >= `ratio` of a longer kept
/// contig (either strand). Coverage accumulates every near-identical block
/// found by anchoring short seeds, so a contig split into several blocks by
/// small junction differences is still detected; multi-coverage-set unitigs
/// are exact, so the near-duplicates differ only at their boundaries and the
/// longest representative is kept.
fn dedup_contained_ratio(mut contigs: Vec<Contig>, ratio: f64) -> Vec<Contig> {
    contigs.sort_by_key(|c| std::cmp::Reverse(c.seq.len()));
    let mut kept: Vec<Contig> = Vec::with_capacity(contigs.len());
    for c in contigs {
        let rc = rev_comp(&c.seq).collect::<Vec<u8>>();
        let contained = if ratio >= 1.0 {
            // Exact substring semantics (historical behaviour).
            kept.iter()
                .any(|k| contains(&k.seq, &c.seq) || contains(&k.seq, &rc))
        } else {
            // Approximate containment: boundary-differing near-duplicates.
            kept.iter()
                .any(|k| coverage(&k.seq, &c.seq) >= ratio || coverage(&k.seq, &rc) >= ratio)
        };
        if !contained {
            kept.push(c);
        }
    }
    kept
}

/// Exact substring test.
fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    !needle.is_empty() && haystack.windows(needle.len()).any(|w| w == needle)
}

/// Fraction of `needle` covered by near-identical segments inside `haystack`
/// (either orientation is the caller's job). Anchors three seeds (head,
/// middle, tail), extends each hit while the accumulated mismatch rate stays
/// within 1%, and merges the covered intervals — a contig fully present as
/// several blocks (small junction differences) is still detected (matching
/// `anchr contained --idt 0.99` semantics).
fn coverage(haystack: &[u8], needle: &[u8]) -> f64 {
    if needle.len() < 100 {
        return if contains(haystack, needle) { 1.0 } else { 0.0 };
    }
    let seed_len = 100usize;
    let seeds = [
        0usize,
        needle.len() / 2 - seed_len / 2,
        needle.len() - seed_len,
    ];
    // Covered intervals on the needle, `[start, end)`.
    let mut intervals: Vec<(usize, usize)> = Vec::new();
    for seed_start in seeds {
        let seed = &needle[seed_start..seed_start + seed_len];
        for (p, _) in haystack
            .windows(seed_len)
            .enumerate()
            .filter(|(_, w)| *w == seed)
        {
            let mut left = 0usize;
            let mut left_mm = 0usize;
            while p as i64 - left as i64 > 0 && seed_start as i64 - left as i64 > 0 {
                if haystack[p - left - 1] == needle[seed_start - left - 1] {
                    left += 1;
                } else if left_mm * 100 <= left + left_mm {
                    left += 1;
                    left_mm += 1; // tolerate ~1% mismatches
                } else {
                    break;
                }
            }
            let mut right = 0usize;
            let mut right_mm = 0usize;
            while p + seed_len + right < haystack.len()
                && seed_start + seed_len + right < needle.len()
            {
                if haystack[p + seed_len + right] == needle[seed_start + seed_len + right] {
                    right += 1;
                } else if right_mm * 100 <= right + right_mm {
                    right += 1;
                    right_mm += 1; // tolerate ~1% mismatches
                } else {
                    break;
                }
            }
            intervals.push((seed_start - left, seed_start + seed_len + right));
        }
    }
    // Merge overlapping intervals and sum the covered length.
    intervals.sort_unstable();
    let mut covered = 0usize;
    let (mut cur_start, mut cur_end) = (0usize, 0usize);
    for (s, e) in intervals {
        if s > cur_end {
            covered += cur_end - cur_start;
            cur_start = s;
            cur_end = e;
        } else {
            cur_end = cur_end.max(e);
        }
    }
    covered += cur_end - cur_start;
    covered.min(needle.len()) as f64 / needle.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::libs::olc::layout::LayoutStep;
    use pgr::libs::nt::rev_comp;

    fn unitigs(names: &[&str], seqs: &[&str]) -> Vec<Unitig> {
        names
            .iter()
            .zip(seqs)
            .map(|(n, s)| Unitig {
                name: (*n).to_string(),
                seq: s.as_bytes().to_vec(),
            })
            .collect()
    }

    fn layout(steps: Vec<LayoutStep>) -> Layout {
        Layout { steps }
    }

    /// Forward chain stitches into the full contig.
    #[test]
    fn stitches_forward_chain() {
        let us = unitigs(
            &["u0", "u1", "u2"],
            &["AAAAAAAAACGTACGT", "ACGTACGTCCCCCCCC", "CCCCCCCCGGGGGGGG"],
        );
        let layouts = vec![layout(vec![
            LayoutStep {
                unitig: 0,
                strand: '+',
                q_start: 0,
                q_end: 16,
                overlap_len: 0,
            },
            LayoutStep {
                unitig: 1,
                strand: '+',
                q_start: 8,
                q_end: 24,
                overlap_len: 8,
            },
            LayoutStep {
                unitig: 2,
                strand: '+',
                q_start: 16,
                q_end: 32,
                overlap_len: 8,
            },
        ])];
        let contigs = consensus_with_ratio(&us, &layouts, 1, 1.0).unwrap();
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].seq, b"AAAAAAAAACGTACGTCCCCCCCCGGGGGGGG");
        // 48 unitig bases over a 32 bp contig.
        assert!((contigs[0].coverage - 1.5).abs() < 1e-9);
    }

    /// Reverse-strand step stitches via its reverse complement.
    #[test]
    fn stitches_reverse_step() {
        let u0 = "TTTTACGTAC";
        let u1 = String::from_utf8(rev_comp(b"ACGTACCCCC").collect()).unwrap();
        let us = unitigs(&["u0", "u1"], &[u0, &u1]);
        let layouts = vec![layout(vec![
            LayoutStep {
                unitig: 0,
                strand: '+',
                q_start: 0,
                q_end: 10,
                overlap_len: 0,
            },
            LayoutStep {
                unitig: 1,
                strand: '-',
                q_start: 4,
                q_end: 14,
                overlap_len: 6,
            },
        ])];
        let contigs = consensus_with_ratio(&us, &layouts, 1, 1.0).unwrap();
        assert_eq!(contigs[0].seq, b"TTTTACGTACCCCC");
    }

    /// Short layouts are dropped by the minimum contig length.
    #[test]
    fn filters_short_contigs() {
        let us = unitigs(&["u0"], &["ACGTACGT"]);
        let layouts = vec![layout(vec![LayoutStep {
            unitig: 0,
            strand: '+',
            q_start: 0,
            q_end: 8,
            overlap_len: 0,
        }])];
        assert_eq!(
            consensus_with_ratio(&us, &layouts, 9, 1.0).unwrap().len(),
            0
        );
        assert_eq!(
            consensus_with_ratio(&us, &layouts, 8, 1.0).unwrap().len(),
            1
        );
    }

    /// A disagreeing overlap is a friendly error, not a panic.
    #[test]
    fn disagreeing_overlap_errors() {
        let us = unitigs(&["u0", "u1"], &["AAAACCCC", "CCCCGGGG"]);
        let layouts = vec![layout(vec![
            LayoutStep {
                unitig: 0,
                strand: '+',
                q_start: 0,
                q_end: 8,
                overlap_len: 0,
            },
            LayoutStep {
                unitig: 1,
                strand: '+',
                q_start: 2,
                q_end: 10,
                overlap_len: 6,
            },
        ])];
        // The claimed 6 bp overlap does not match (u0 suffix "AACCCC" vs
        // u1 prefix "CCCCGG"): the stitch must fail cleanly.
        let err = consensus_with_ratio(&us, &layouts, 1, 1.0).unwrap_err();
        assert!(err.to_string().contains("disagree"), "{err}");
    }

    /// Contigs fully contained in a longer contig are dropped.
    #[test]
    fn dedups_contained_contigs() {
        let us = unitigs(&["u0", "u1"], &["AAAACCCCGGGGTTTT", "CCCCGGGG"]);
        let layouts = vec![
            layout(vec![LayoutStep {
                unitig: 0,
                strand: '+',
                q_start: 0,
                q_end: 16,
                overlap_len: 0,
            }]),
            layout(vec![LayoutStep {
                unitig: 1,
                strand: '+',
                q_start: 0,
                q_end: 8,
                overlap_len: 0,
            }]),
        ];
        let contigs = consensus_with_ratio(&us, &layouts, 1, 1.0).unwrap();
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].seq, b"AAAACCCCGGGGTTTT");
    }

    /// Reverse-complement containment also drops the shorter contig.
    #[test]
    fn dedups_rc_contained_contigs() {
        let us = unitigs(
            &["u0", "u1"],
            &["AAAATACGTACGTTTT", "CGTACGTA"], // rc(CGTACGTA) = TACGTACG
        );
        let layouts = vec![
            layout(vec![LayoutStep {
                unitig: 0,
                strand: '+',
                q_start: 0,
                q_end: 16,
                overlap_len: 0,
            }]),
            layout(vec![LayoutStep {
                unitig: 1,
                strand: '+',
                q_start: 0,
                q_end: 8,
                overlap_len: 0,
            }]),
        ];
        let contigs = consensus_with_ratio(&us, &layouts, 1, 1.0).unwrap();
        assert_eq!(contigs.len(), 1);
        assert_eq!(contigs[0].seq, b"AAAATACGTACGTTTT");
    }

    /// `--dedup-ratio < 1.0` merges near-duplicate contigs whose boundaries
    /// differ by a few bases: exact containment (ratio 1.0) keeps both, the
    /// approximate rule keeps the longer one.
    #[test]
    fn dedups_approximate_contained_contigs() {
        // long = A*95 + ACGT*25 + T*10; short = same but last two bases GG.
        let long: Vec<u8> = b"A"
            .repeat(95)
            .into_iter()
            .chain(b"ACGT".repeat(25))
            .chain(b"T".repeat(10))
            .collect();
        let short: Vec<u8> = b"A"
            .repeat(95)
            .into_iter()
            .chain(b"ACGT".repeat(25))
            .chain(b"T".repeat(8))
            .chain([b'G', b'G'])
            .collect();
        let us = unitigs(
            &["u0", "u1"],
            &[
                std::str::from_utf8(&long).unwrap(),
                std::str::from_utf8(&short).unwrap(),
            ],
        );
        let layouts = vec![
            layout(vec![LayoutStep {
                unitig: 0,
                strand: '+',
                q_start: 0,
                q_end: long.len(),
                overlap_len: 0,
            }]),
            layout(vec![LayoutStep {
                unitig: 1,
                strand: '+',
                q_start: 0,
                q_end: short.len(),
                overlap_len: 0,
            }]),
        ];
        let exact = consensus_with_ratio(&us, &layouts, 1, 1.0).unwrap();
        assert_eq!(exact.len(), 2, "not an exact substring: both kept");
        let approx = consensus_with_ratio(&us, &layouts, 1, 0.95).unwrap();
        assert_eq!(approx.len(), 1, "boundary-differing duplicate merged");
        assert_eq!(approx[0].seq.len(), long.len());
    }

    /// A contig fully covered by a longer one through two separate blocks
    /// (a small junction difference between them) is deduplicated with
    /// `--dedup-ratio < 1.0`: the blocks are summed, not just the longest.
    #[test]
    fn dedups_multi_block_contained_contigs() {
        fn pseudo(len: usize, seed: u64) -> Vec<u8> {
            let mut rng = seed;
            let mut s = Vec::with_capacity(len);
            for _ in 0..len {
                rng = rng
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                s.push(b"ACGT"[(rng >> 33) as usize % 4]);
            }
            s
        }
        // short = B1 + M' + B2 vs long = B1 + M + B2: the 5-bp middle breaks
        // the single-segment extension, but the B1/B2 blocks cover ~99.8%.
        let b1 = pseudo(200, 42);
        let b2 = pseudo(200, 43);
        let mut long = b1.clone();
        long.extend(b"CCCCC");
        long.extend(&b2);
        let mut short = b1.clone();
        short.extend(b"GGGGG");
        short.extend(&b2);
        let us = unitigs(
            &["u0", "u1"],
            &[
                std::str::from_utf8(&long).unwrap(),
                std::str::from_utf8(&short).unwrap(),
            ],
        );
        let layouts = vec![
            layout(vec![LayoutStep {
                unitig: 0,
                strand: '+',
                q_start: 0,
                q_end: long.len(),
                overlap_len: 0,
            }]),
            layout(vec![LayoutStep {
                unitig: 1,
                strand: '+',
                q_start: 0,
                q_end: short.len(),
                overlap_len: 0,
            }]),
        ];
        let exact = consensus_with_ratio(&us, &layouts, 1, 1.0).unwrap();
        assert_eq!(exact.len(), 2, "not a single exact substring: both kept");
        let approx = consensus_with_ratio(&us, &layouts, 1, 0.99).unwrap();
        assert_eq!(approx.len(), 1, "multi-block duplicate merged");
        assert_eq!(approx[0].seq.len(), long.len());
    }

    /// Distinct contigs are all kept, longest first.
    #[test]
    fn keeps_distinct_contigs() {
        let us = unitigs(&["u0", "u1"], &["AAAACCCCGGGG", "TTTTCCCCAAAA"]);
        let layouts = vec![
            layout(vec![LayoutStep {
                unitig: 0,
                strand: '+',
                q_start: 0,
                q_end: 12,
                overlap_len: 0,
            }]),
            layout(vec![LayoutStep {
                unitig: 1,
                strand: '+',
                q_start: 0,
                q_end: 12,
                overlap_len: 0,
            }]),
        ];
        let contigs = consensus_with_ratio(&us, &layouts, 1, 1.0).unwrap();
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs[0].seq, b"AAAACCCCGGGG");
        assert_eq!(contigs[1].seq, b"TTTTCCCCAAAA");
    }

    /// A candidate whose head aligns inside the keeper and whose tail
    /// extends past the keeper's end is stitched into one contig.
    #[test]
    fn merges_tail_extension_overlap() {
        let mut keeper = vec![b'A'; 5000];
        keeper.extend(std::iter::repeat(b'C').take(20000));
        keeper.extend(std::iter::repeat(b'T').take(5000));
        let mut cand = vec![b'C'; 20000];
        cand.extend(std::iter::repeat(b'T').take(5000));
        cand.extend(std::iter::repeat(b'G').take(3000));
        let merged = try_merge(&keeper, &cand, 5000).unwrap();
        assert_eq!(merged.len(), 33000);
        assert_eq!(&merged[..5000], &keeper[..5000]);
        assert_eq!(&merged[5000..25000], &cand[..20000]);
        assert_eq!(&merged[25000..30000], &keeper[25000..]);
        assert_eq!(&merged[30000..], &cand[25000..]);
    }

    /// A candidate fully inside the keeper is dropped (keeper kept).
    #[test]
    fn drops_contained_candidate() {
        let mut keeper = vec![b'A'; 1000];
        keeper.extend(std::iter::repeat(b'C').take(20000));
        keeper.extend(std::iter::repeat(b'T').take(5000));
        let cand = vec![b'C'; 20000];
        let merged = try_merge(&keeper, &cand, 5000).unwrap();
        assert_eq!(merged, keeper);
    }

    /// Reverse-complement overlap merges with the keeper's strand.
    #[test]
    fn merges_reverse_overlap() {
        let mut keeper = vec![b'A'; 5000];
        keeper.extend(std::iter::repeat(b'C').take(20000));
        keeper.extend(std::iter::repeat(b'T').take(5000));
        let mut cand = vec![b'C'; 20000];
        cand.extend(std::iter::repeat(b'T').take(5000));
        cand.extend(std::iter::repeat(b'G').take(3000));
        let rc_cand: Vec<u8> = rev_comp(&cand).collect();
        let merged = try_merge(&keeper, &rc_cand, 5000).unwrap();
        assert_eq!(merged.len(), 33000);
        assert_eq!(&merged[30000..], &b"G".repeat(3000)[..]);
    }

    /// A chimeric candidate (head aligned elsewhere) is left untouched.
    #[test]
    fn rejects_chimeric_candidate() {
        let mut keeper = vec![b'A'; 5000];
        keeper.extend(std::iter::repeat(b'C').take(20000));
        keeper.extend(std::iter::repeat(b'T').take(5000));
        let mut cand = vec![b'X'; 3000];
        cand.extend(std::iter::repeat(b'C').take(20000));
        cand.extend(std::iter::repeat(b'T').take(5000));
        cand.extend(std::iter::repeat(b'G').take(3000));
        assert!(try_merge(&keeper, &cand, 5000).is_none());
    }
}

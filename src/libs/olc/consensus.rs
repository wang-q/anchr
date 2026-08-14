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
    Ok(dedup_contained_ratio(contigs, ratio))
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
}

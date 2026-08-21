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
    // `ratio` gates approximate dedup (`cns --dedup-ratio`): <= 0 would drop
    // any contig sharing a single 31-mer with a kept contig (degenerate),
    // and > 1 silently degrades to the exact-substring path. Mirror the
    // `--merge-similar` range check in the multik driver.
    anyhow::ensure!(
        ratio > 0.0 && ratio <= 1.0,
        "dedup ratio must be in (0.0, 1.0], got {ratio}"
    );
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
    // A merge is only possible when `cand`'s head (or its reverse
    // complement) is inside the keeper, or the keeper's head (or its rc) is
    // inside `cand` (the two anchor conditions of `merge_geometry`). Those
    // O(1) set checks reject the overwhelming majority of the O(n^2) pairs
    // before the expensive seed-index/identity verification in `try_merge`.
    // Keeper sequences grow on merge, so each kept contig carries its own
    // up-to-date seed set and head seeds.
    let mut kept: Vec<Contig> = Vec::with_capacity(contigs.len());
    let mut kept_seeds: Vec<std::collections::HashSet<u64>> = Vec::new();
    let mut kept_heads: Vec<[u64; 2]> = Vec::new();
    for c in contigs {
        if c.seq.len() < min_overlap {
            // Both sides need at least `min_overlap` bases for a stitch.
            kept.push(c);
            kept_seeds.push(std::collections::HashSet::new());
            kept_heads.push([0; 2]);
            continue;
        }
        let [cand_head, cand_rc_head] = boundary_seeds(&c.seq);
        let cand_seeds = seed_set(&c.seq);
        let mut merged = false;
        for pos in 0..kept.len() {
            if kept[pos].seq.len() < min_overlap {
                continue;
            }
            let k_has =
                kept_seeds[pos].contains(&cand_head) || kept_seeds[pos].contains(&cand_rc_head);
            let c_has = cand_seeds.contains(&kept_heads[pos][0])
                || cand_seeds.contains(&kept_heads[pos][1]);
            if !k_has && !c_has {
                continue;
            }
            if let Some(seq) = try_merge(&kept[pos].seq, &c.seq, min_overlap) {
                let old_len = kept[pos].seq.len();
                let appended_start = if seq.starts_with(&kept[pos].seq) {
                    Some(old_len)
                } else if seq.ends_with(&kept[pos].seq) {
                    Some(seq.len() - old_len)
                } else {
                    None
                };
                kept[pos].seq = seq;
                if let Some(start) = appended_start {
                    extend_seed_set(&mut kept_seeds[pos], &kept[pos].seq, start);
                    if start == 0 {
                        kept_heads[pos] = boundary_seeds(&kept[pos].seq);
                    }
                }
                merged = true;
                break;
            }
        }
        if !merged {
            kept.push(c);
            kept_seeds.push(cand_seeds);
            kept_heads.push([cand_head, cand_rc_head]);
        }
    }
    kept
}

/// Packs a 31-mer into a u64 (2 bits per base; exact for 31 bases).
/// Returns [`u64::MAX`] for non-ACGT bases, which no packed seed can equal.
fn seed_u64(w: &[u8]) -> u64 {
    let mut v = 0u64;
    for &b in w {
        v = (v << 2)
            | match b {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                b'T' => 3,
                _ => return u64::MAX,
            };
    }
    v
}

/// The forward and reverse-complement head seeds of a sequence.
fn boundary_seeds(seq: &[u8]) -> [u64; 2] {
    let head = seed_u64(&seq[..SEED_LEN]);
    let rc_head = seed_u64(&rev_comp(&seq[..SEED_LEN]).collect::<Vec<_>>());
    [head, rc_head]
}

/// Distinct 31-mer seeds of a sequence. Only called on ACGT-only sequences
/// (the caller checks the boundary seeds first).
fn seed_set(seq: &[u8]) -> std::collections::HashSet<u64> {
    let mut set = std::collections::HashSet::new();
    for w in seq.windows(SEED_LEN) {
        set.insert(seed_u64(w));
    }
    set
}

/// Adds the seeds of the newly appended part (plus the junction windows)
/// after a merge. `start` is the index of the first appended base.
fn extend_seed_set(set: &mut std::collections::HashSet<u64>, seq: &[u8], start: usize) {
    let from = start.saturating_sub(SEED_LEN - 1);
    for w in seq[from..].windows(SEED_LEN) {
        set.insert(seed_u64(w));
    }
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
///
/// Normal inputs run the exact histogram: every (query window, keeper seed
/// position) pair votes for its offset. Homopolymer runs make that quadratic
/// (20k x 20k = 4e8 pairs per contig pair), so once the pair count exceeds
/// [`EXACT_WORK_CAP`] a bounded path takes over: heavy seeds (many positions
/// on both sides) contribute a piecewise-linear support function whose exact
/// breakpoints come from run pairs, light seeds contribute an exact offset
/// delta histogram, and each candidate offset is verified with an exact
/// per-window scan. The verified maximum is exact, so repetitive inputs keep
/// the same dominant-offset answer without the quadratic fan-out.
fn overlap_geometry(
    index: &std::collections::HashMap<&[u8], Vec<usize>>,
    query: &[u8],
) -> (isize, usize, usize) {
    let mut matched = 0usize;
    let mut work = 0usize;
    for w in query.windows(SEED_LEN) {
        if let Some(ps) = index.get(w) {
            matched += 1;
            work += ps.len();
        }
    }
    if work > EXACT_WORK_CAP {
        return overlap_geometry_bounded(index, query, matched);
    }
    let mut hist: std::collections::HashMap<isize, usize> = std::collections::HashMap::new();
    for (i, w) in query.windows(SEED_LEN).enumerate() {
        if let Some(ps) = index.get(w) {
            for &p in ps {
                *hist.entry(p as isize - i as isize).or_default() += 1;
            }
        }
    }
    let (offset, hits) = hist.into_iter().max_by_key(|(_, n)| *n).unwrap_or((0, 0));
    (offset, hits, matched)
}

/// Exact-histogram pair budget; above this the bounded path runs instead.
const EXACT_WORK_CAP: usize = 1_000_000;

/// Pair budget per seed above which the seed counts as heavy (run-pair
/// breakpoints instead of the plain pair loop).
const HEAVY_PAIRS: usize = 4096;

/// Candidate offsets verified by the exact per-window scan.
const MAX_CANDIDATES: usize = 32;

/// Maximal runs of consecutive window positions.
fn runs(pos: &[usize]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < pos.len() {
        let start = pos[i];
        let mut j = i + 1;
        while j < pos.len() && pos[j] == pos[j - 1] + 1 {
            j += 1;
        }
        out.push((start, pos[j - 1] + 1));
        i = j;
    }
    out
}

/// Bounded fallback for repetitive inputs. `matched` counts every query
/// window present in `index` (exact). The support of a keeper-run / query-run
/// pair is a triangle in the offset, so heavy seeds (many position pairs)
/// propose its breakpoints; light seeds accumulate an exact delta histogram.
/// Each candidate's support is then verified exactly over the full query, so
/// the returned maximum matches the exact histogram whenever the true
/// dominant offset is among the candidates.
fn overlap_geometry_bounded(
    index: &std::collections::HashMap<&[u8], Vec<usize>>,
    query: &[u8],
    matched: usize,
) -> (isize, usize, usize) {
    let mut qpos: std::collections::HashMap<&[u8], Vec<usize>> = std::collections::HashMap::new();
    for (i, w) in query.windows(SEED_LEN).enumerate() {
        if index.contains_key(w) {
            qpos.entry(w).or_default().push(i);
        }
    }
    // Heavy seeds: slope events of the support function (offset -> slope
    // delta); a run pair contributes a +1/0/-1 triangle between its four
    // breakpoints. Light seeds: exact offset delta histogram.
    let mut events: std::collections::BTreeMap<isize, i64> = std::collections::BTreeMap::new();
    let mut deltas: std::collections::HashMap<isize, usize> = std::collections::HashMap::new();
    for (w, qs) in &qpos {
        let ps = &index[w];
        if ps.len() * qs.len() > HEAVY_PAIRS {
            for (ks, ke) in runs(ps) {
                for (qr_s, qr_e) in runs(qs) {
                    *events.entry(ks as isize - qr_e as isize).or_default() += 1;
                    *events.entry(ks as isize - qr_s as isize).or_default() -= 1;
                    *events.entry(ke as isize - qr_e as isize).or_default() -= 1;
                    *events.entry(ke as isize - qr_s as isize).or_default() += 1;
                }
            }
        } else {
            for &p in ps {
                for &i in qs {
                    *deltas.entry(p as isize - i as isize).or_default() += 1;
                }
            }
        }
    }
    // Sweep the slope events to get the exact heavy-seed support at every
    // breakpoint, then keep the strongest ones.
    let mut heavy: Vec<(isize, i64)> = Vec::with_capacity(events.len());
    let mut value = 0i64;
    let mut slope = 0i64;
    let mut prev = 0isize;
    for (o, d) in events {
        value += slope * (o - prev) as i64;
        heavy.push((o, value));
        slope += d;
        prev = o;
    }
    heavy.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    heavy.truncate(MAX_CANDIDATES / 2);
    let mut light: Vec<(isize, usize)> = deltas.into_iter().filter(|(_, c)| *c >= 2).collect();
    light.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    light.truncate(MAX_CANDIDATES / 2);
    let mut candidates: Vec<isize> = Vec::with_capacity(MAX_CANDIDATES);
    candidates.extend(heavy.into_iter().map(|(o, _)| o));
    candidates.extend(light.into_iter().map(|(o, _)| o));
    candidates.sort_unstable();
    candidates.dedup();
    candidates.truncate(MAX_CANDIDATES);
    let mut best: Option<(isize, usize)> = None;
    for &o in &candidates {
        let mut support = 0usize;
        for (i, w) in query.windows(SEED_LEN).enumerate() {
            if let Some(ps) = index.get(w) {
                let pos = i as isize + o;
                if pos >= 0 && ps.binary_search(&(pos as usize)).is_ok() {
                    support += 1;
                }
            }
        }
        if best.is_none_or(|(_, b)| support > b) {
            best = Some((o, support));
        }
    }
    let (offset, hits) = best.unwrap_or((0, 0));
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

/// Identity of two equally long slices with a small-band edit-distance
/// tolerance for indels. A plain zip comparison misaligns after the first
/// indel and under-reports identity for near-duplicate contigs that differ
/// by a few indels (cross-group consensus of the same locus), which would
/// otherwise skip their merge.
fn identity(a: &[u8], b: &[u8]) -> f64 {
    debug_assert_eq!(a.len(), b.len());
    let n = a.len();
    if n == 0 {
        return 1.0;
    }
    // Banded Levenshtein distance: only cells with |i - j| <= band are
    // evaluated, so a few scattered indels do not cascade the comparison.
    // 512 covers cross-group consensus indels of up to half a kilobase
    // (e.g. a 108 bp insertion seen on G37 K160 contigs) while keeping the
    // O(n * 2*band) cost bounded.
    const BAND: usize = 512;
    let band = BAND.min(n);
    let width = 2 * band + 1;
    // dp[j] = edit distance for the current row; offsets are `j - band`
    // relative to the row index (all uninitialized cells are INF).
    let inf = usize::MAX / 4;
    let mut prev = vec![inf; width];
    let mut cur = vec![inf; width];
    prev[band] = 0; // dp[0][0]
    for (i, &ai) in a.iter().enumerate() {
        cur.fill(inf);
        let j_lo = i.saturating_sub(band);
        let j_hi = (i + band).min(n - 1);
        for (jj, &bj) in b[j_lo..=j_hi].iter().enumerate() {
            let j = j_lo + jj;
            let c = j + band - i;
            if c >= width {
                continue;
            }
            let mut best = inf;
            // substitution / match
            if i > 0 && j > 0 {
                let d = prev[c] + usize::from(ai != bj);
                best = best.min(d);
            } else if i == 0 && j == 0 {
                best = usize::from(ai != bj);
            }
            // deletion (consume a[i] only)
            if i > 0 && c + 1 < width {
                best = best.min(prev[c + 1] + 1);
            }
            // insertion (consume b[j] only)
            if j > 0 && c >= 1 {
                best = best.min(cur[c - 1] + 1);
            }
            cur[c] = best;
        }
        std::mem::swap(&mut prev, &mut cur);
    }
    let ed = prev[band + n - 1 - (n - 1)];
    let ed = if ed >= inf { n } else { ed };
    1.0 - ed as f64 / n as f64
}

/// Drops contigs whose sequence is covered by >= `ratio` of a longer kept
/// contig (either strand). Coverage accumulates every near-identical block
/// found by anchoring short seeds, so a contig split into several blocks by
/// small junction differences is still detected; multi-coverage-set unitigs
/// are exact, so the near-duplicates differ only at their boundaries and the
/// longest representative is kept.
fn dedup_contained_ratio(mut contigs: Vec<Contig>, ratio: f64) -> Vec<Contig> {
    contigs.sort_by_key(|c| std::cmp::Reverse(c.seq.len()));
    // Candidate prefilter: a shorter contig can only be contained in a
    // longer one when they share at least one 31-mer (either strand) — the
    // exact 100-mer anchors of `coverage` (and the exact-substring test for
    // short needles) imply such a shared seed. One global seed index answers
    // the candidate query in O(length) per contig instead of comparing every
    // kept pair with an O(length) window scan.
    let index = build_seed_index(&contigs);
    let mut kept: Vec<Contig> = Vec::with_capacity(contigs.len());
    let mut kept_ids: Vec<usize> = Vec::with_capacity(contigs.len());
    for (id, c) in contigs.iter().enumerate() {
        let rc = rev_comp(&c.seq).collect::<Vec<u8>>();
        let contained = if c.seq.len() < SEED_LEN || c.seq.contains(&b'N') || rc.contains(&b'N') {
            // No 31-mer prefilter applies (too short or non-ACGT): fall back
            // to the exhaustive scan.
            kept.iter()
                .any(|k| contains(&k.seq, &c.seq) || contains(&k.seq, &rc))
        } else if ratio >= 1.0 {
            // Exact substring semantics (historical behaviour).
            kept.iter()
                .any(|k| contains(&k.seq, &c.seq) || contains(&k.seq, &rc))
        } else {
            // Approximate containment: boundary-differing near-duplicates.
            let mut cands: Vec<u32> = Vec::new();
            for w in c.seq.windows(SEED_LEN) {
                push_candidates(&index, seed_u64(w), &mut cands);
                let rc_w: Vec<u8> = rev_comp(w).collect();
                push_candidates(&index, seed_u64(&rc_w), &mut cands);
            }
            cands.sort_unstable();
            cands.dedup();
            cands.into_iter().any(|k| {
                kept_ids.binary_search(&(k as usize)).is_ok()
                    && (coverage(&contigs[k as usize].seq, &c.seq) >= ratio
                        || coverage(&contigs[k as usize].seq, &rc) >= ratio)
            })
        };
        if !contained {
            kept_ids.push(id);
            kept.push(Contig {
                seq: c.seq.clone(),
                coverage: c.coverage,
            });
        }
    }
    kept
}

/// One global 31-mer seed -> contig-id index over all contigs (forward
/// strand only; non-ACGT windows are skipped, the N-fallback covers them).
fn build_seed_index(contigs: &[Contig]) -> std::collections::HashMap<u64, Vec<u32>> {
    let mut map: std::collections::HashMap<u64, Vec<u32>> = std::collections::HashMap::new();
    for (id, c) in contigs.iter().enumerate() {
        for w in c.seq.windows(SEED_LEN) {
            let s = seed_u64(w);
            if s != u64::MAX {
                map.entry(s).or_default().push(id as u32);
            }
        }
    }
    map
}

/// Appends the contig ids of a seed to `out` (no-op for non-ACGT seeds).
fn push_candidates(
    index: &std::collections::HashMap<u64, Vec<u32>>,
    seed: u64,
    out: &mut Vec<u32>,
) {
    if let Some(ids) = index.get(&seed) {
        out.extend_from_slice(ids);
    }
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
pub(crate) fn coverage(haystack: &[u8], needle: &[u8]) -> f64 {
    if needle.len() < 100 {
        return if contains(haystack, needle) { 1.0 } else { 0.0 };
    }
    // Indel-tolerant path: a dominant 31-mer offset locates `needle` inside
    // `haystack` even across large internal indels (a 120 bp insertion made
    // the fixed 100-mer anchors below miss the tail); the banded identity
    // then measures how much of `needle` is covered.
    let mut index: std::collections::HashMap<&[u8], Vec<usize>> = std::collections::HashMap::new();
    for (p, w) in haystack.windows(SEED_LEN).enumerate() {
        index.entry(w).or_default().push(p);
    }
    // Reuses the exact-offset histogram (and its bounded fallback for
    // repetitive inputs) of `overlap_geometry` so homopolymer-rich contigs
    // do not fan out quadratically.
    let (off, hits, matched) = overlap_geometry(&index, needle);
    if hits >= 100 && hits as f64 / matched.max(1) as f64 >= 0.6 {
        // A large internal indel can split the histogram into two near-equal
        // peaks (e.g. 66%/34% for a 120 bp insertion); the peak-ratio bar is
        // only a candidate filter, the >=99% identity check does the judging.
        let hay_len = haystack.len() as isize;
        let nd_len = needle.len() as isize;
        let start = off.max(0) as usize;
        let end = (off + nd_len).clamp(0, hay_len) as usize;
        if end > start {
            let a = &haystack[start..end];
            let nstart = if off < 0 { off.unsigned_abs() } else { 0 };
            let nend = nstart + (end - start);
            if nend <= needle.len() {
                let idy = identity(a, &needle[nstart..nend]);
                let covered = (end - start) as f64 * idy / nd_len as f64;
                if covered >= 0.99 {
                    return covered.min(1.0);
                }
            }
        }
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
            .chain(*b"GG")
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
        keeper.extend(std::iter::repeat_n(b'C', 20000));
        keeper.extend(std::iter::repeat_n(b'T', 5000));
        let mut cand = vec![b'C'; 20000];
        cand.extend(std::iter::repeat_n(b'T', 5000));
        cand.extend(std::iter::repeat_n(b'G', 3000));
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
        keeper.extend(std::iter::repeat_n(b'C', 20000));
        keeper.extend(std::iter::repeat_n(b'T', 5000));
        let cand = vec![b'C'; 20000];
        let merged = try_merge(&keeper, &cand, 5000).unwrap();
        assert_eq!(merged, keeper);
    }

    /// Reverse-complement overlap merges with the keeper's strand.
    #[test]
    fn merges_reverse_overlap() {
        let mut keeper = vec![b'A'; 5000];
        keeper.extend(std::iter::repeat_n(b'C', 20000));
        keeper.extend(std::iter::repeat_n(b'T', 5000));
        let mut cand = vec![b'C'; 20000];
        cand.extend(std::iter::repeat_n(b'T', 5000));
        cand.extend(std::iter::repeat_n(b'G', 3000));
        let rc_cand: Vec<u8> = rev_comp(&cand).collect();
        let merged = try_merge(&keeper, &rc_cand, 5000).unwrap();
        assert_eq!(merged.len(), 33000);
        assert_eq!(&merged[30000..], &b"G".repeat(3000)[..]);
    }

    /// A chimeric candidate (head aligned elsewhere) is left untouched.
    #[test]
    fn rejects_chimeric_candidate() {
        let mut keeper = vec![b'A'; 5000];
        keeper.extend(std::iter::repeat_n(b'C', 20000));
        keeper.extend(std::iter::repeat_n(b'T', 5000));
        let mut cand = vec![b'X'; 3000];
        cand.extend(std::iter::repeat_n(b'C', 20000));
        cand.extend(std::iter::repeat_n(b'T', 5000));
        cand.extend(std::iter::repeat_n(b'G', 3000));
        assert!(try_merge(&keeper, &cand, 5000).is_none());
    }

    /// The bounded fallback (work > [`EXACT_WORK_CAP`] pairs) must return
    /// the same maximum as the plain histogram.
    #[test]
    fn bounded_geometry_matches_exact_histogram() {
        fn exact(
            index: &std::collections::HashMap<&[u8], Vec<usize>>,
            query: &[u8],
        ) -> (isize, usize, usize) {
            let mut hist: std::collections::HashMap<isize, usize> =
                std::collections::HashMap::new();
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
        fn index_of(seq: &[u8]) -> std::collections::HashMap<&[u8], Vec<usize>> {
            let mut index: std::collections::HashMap<&[u8], Vec<usize>> =
                std::collections::HashMap::new();
            for (i, w) in seq.windows(SEED_LEN).enumerate() {
                index.entry(w).or_default().push(i);
            }
            index
        }
        fn run(c: u8, n: usize) -> Vec<u8> {
            vec![c; n]
        }
        // Equal-length runs: unique maximum, the offset must match exactly.
        let keeper = concat(&run(b'A', 100), &concat(&run(b'C', 1200), &run(b'T', 100)));
        for query in [run(b'C', 1200), concat(&run(b'C', 1200), &run(b'T', 100))] {
            let index = index_of(&keeper);
            let got = overlap_geometry(&index, &query);
            let want = exact(&index, &query);
            assert_eq!(got, want, "keeper={} query={}", keeper.len(), query.len());
        }
        // Keeper run longer than query run: the support plateaus over a range
        // of offsets (same maximum hits); the fallback may pick a different
        // plateau member than the histogram, so compare the decisive values.
        let cases: Vec<(Vec<u8>, Vec<u8>)> = vec![
            (keeper.clone(), run(b'C', 1000)),
            (
                keeper.clone(),
                concat(&run(b'X', 100), &concat(&run(b'C', 1000), &run(b'T', 100))),
            ),
        ];
        for (keeper, query) in cases {
            let index = index_of(&keeper);
            let got = overlap_geometry(&index, &query);
            let want = exact(&index, &query);
            assert_eq!(
                (got.1, got.2),
                (want.1, want.2),
                "keeper={} query={}",
                keeper.len(),
                query.len()
            );
        }
    }

    /// The first base pair must count toward the edit distance (a plain
    /// base-case `best = 0` silently treated `a[0] != b[0]` as a match,
    /// over-reporting identity for boundary-differing near-duplicates).
    #[test]
    fn identity_counts_first_base_mismatch() {
        // Differing first base, matching second: edit distance 1, not 0.
        assert!((identity(b"AC", b"GC") - 0.5).abs() < 1e-9);
        // Single-base pair: a mismatch is identity 0, not 1.
        assert!((identity(b"A", b"G") - 0.0).abs() < 1e-9);
        assert!((identity(b"A", b"A") - 1.0).abs() < 1e-9);
    }

    /// Substitution-only sequences score `1 - mismatches / len`.
    #[test]
    fn identity_scores_substitutions() {
        assert!((identity(b"ACGT", b"ACGT") - 1.0).abs() < 1e-9);
        assert!((identity(b"ACGT", b"ACGA") - 0.75).abs() < 1e-9);
        assert!((identity(b"ACGT", b"AGGT") - 0.75).abs() < 1e-9);
        // A mid-sequence indel costs 2 (delete + insert), so `ATGT` vs
        // `ACGT` is distance 1 (substitute) — the banded path finds it.
        assert!((identity(b"ATGT", b"ACGT") - 0.75).abs() < 1e-9);
    }

    /// The dedup ratio must be in (0.0, 1.0]: 0/negative would drop any
    /// contig sharing a 31-mer with a kept one, and > 1 silently degrades to
    /// exact-substring semantics (friendly error, not silent misbehaviour).
    #[test]
    fn rejects_out_of_range_dedup_ratio() {
        let us = unitigs(&["u0"], &["ACGTACGT"]);
        let layouts = vec![layout(vec![LayoutStep {
            unitig: 0,
            strand: '+',
            q_start: 0,
            q_end: 8,
            overlap_len: 0,
        }])];
        for bad in [0.0, -0.1, 1.1] {
            let err = consensus_with_ratio(&us, &layouts, 1, bad).unwrap_err();
            assert!(err.to_string().contains("dedup ratio"), "{err}");
        }
        assert!(consensus_with_ratio(&us, &layouts, 1, 1.0).is_ok());
        assert!(consensus_with_ratio(&us, &layouts, 1, 0.5).is_ok());
    }

    fn concat(a: &[u8], b: &[u8]) -> Vec<u8> {
        let mut out = a.to_vec();
        out.extend_from_slice(b);
        out
    }
}

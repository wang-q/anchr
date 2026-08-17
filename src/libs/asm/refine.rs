//! Read-level refinement along the k-mer graph.
//!
//! Per-read error correction by local reassembly through the k-mer graph,
//! conservative read extension that stops at branches, and junk/low-depth
//! read discarding. Initially ported from BBTools `tadpole.sh`
//! (correct/extend/discard modes); the implementation has since diverged
//! (long-k path, packed table, streaming counting) and only the
//! command-line compatibility notes remain. The count table itself lives
//! in [`super::table`].

use super::table::{
    base_code, base_comp_code, base_defined, canonicalize_quality, prob_error, Kmer, RefineTable,
};
use anyhow::Result;
use pgr::libs::fmt::seq::{SeqReader, SeqRecord};
use pgr::libs::fq::qual::{from_phred, to_phred};
use pgr::libs::kmer::key;
use pgr::libs::nt::rev_comp;
use std::io::Write;

/// Default k-mer length (tadpole.sh `k`).
pub const DEFAULT_K: usize = 31;

/// Options mirroring the tadpole.sh defaults used by the anchr merge flow.
#[derive(Debug, Clone)]
pub struct RefineOptions {
    /// K-mer length.
    pub k: usize,
    /// Ignore k-mers whose probability of being error-free is below this.
    pub min_prob: f32,
    /// Minimum k-mer depth to seed an extension.
    pub min_count_seed: usize,
    /// Minimum k-mer depth to continue an extension.
    pub min_count_extend: usize,
    /// Branch ratio at high depth (branchmult1).
    pub branch_mult1: f32,
    /// Branch ratio at low depth (branchmult2).
    pub branch_mult2: f32,
    /// Second-highest depth considered "low" (branchlower).
    pub branch_lower_const: usize,
    /// Error ratio multiplier (errormult1).
    pub error_mult1: f32,
    /// Alternative error ratio multiplier (errormult2).
    pub error_mult2: f32,
    /// Quality factor for the error multiplier.
    pub error_mult_q_factor: f32,
    /// Max second-highest depth for the low-depth error rule (errorlowerconst).
    pub error_lower_const: usize,
    /// Minimum depth of a k-mer to be considered correct (mincountcorrect).
    pub min_count_correct: usize,
    /// Absolute path-similarity tolerance (pathsimilarityconstant).
    pub path_similarity_constant: usize,
    /// Fractional path-similarity tolerance (pathsimilarityfraction).
    pub path_similarity_fraction: f32,
    /// K-mers to verify after an error in reassembly (errorextensionreassemble).
    pub error_extension_reassemble: usize,
    /// Do not correct bases within this distance of read ends (deadzone).
    pub dead_zone: usize,
    /// Sliding-window length for reassembly quality filtering (window).
    pub window_len: usize,
    /// Max corrections in a window (windowcount).
    pub window_count: usize,
    /// Max quality sum in a window (qualsum).
    pub window_qual_sum: usize,
    /// Undo corrections that lower k-mer coverage (eccrollback).
    pub ecc_rollback: bool,
    /// Run k-mer reassembly error correction (tadpole `ecc`; off in extend
    /// and discard-only modes, matching Java's per-mode default).
    pub ecc: bool,
    /// Require both directions to agree in the read middle (requirebidirectional).
    pub ecc_require_bidirectional: bool,
    /// Extend to the right by at most this many bases.
    pub extend_right: usize,
    /// Extend to the left by at most this many bases.
    pub extend_left: usize,
    /// Trim random trailing bases of partial extensions (extendrollback).
    pub extension_rollback: usize,
    /// Discard reads that cannot be used for assembly (tossjunk).
    pub toss_junk: bool,
    /// Discard reads containing k-mers at or below this depth (tossdepth).
    pub toss_depth: i64,
    /// Discard reads with uncorrectable errors (tossuncorrectable).
    pub toss_uncorrectable: bool,
    /// Minimum fraction of low-depth k-mers to discard a read (lowdepthfraction).
    pub low_depth_discard_fraction: f32,
    /// Only discard a pair if both reads fail (requirebothbad).
    pub require_both_bad: bool,
}

impl Default for RefineOptions {
    fn default() -> Self {
        Self {
            k: DEFAULT_K,
            min_prob: 0.5,
            min_count_seed: 3,
            min_count_extend: 2,
            branch_mult1: 20.0,
            branch_mult2: 3.0,
            branch_lower_const: 3,
            error_mult1: 16.0,
            error_mult2: 2.6,
            error_mult_q_factor: 0.002,
            error_lower_const: 4,
            min_count_correct: 3,
            path_similarity_constant: 3,
            path_similarity_fraction: 0.45,
            error_extension_reassemble: 5,
            dead_zone: 0,
            window_len: 12,
            window_count: 6,
            window_qual_sum: 80,
            ecc_rollback: true,
            ecc: false,
            ecc_require_bidirectional: true,
            extend_right: 0,
            extend_left: 0,
            extension_rollback: 3,
            toss_junk: false,
            toss_depth: -1,
            toss_uncorrectable: false,
            low_depth_discard_fraction: 0.0,
            require_both_bad: false,
        }
    }
}

/// The forward k-mers of a read (position-wise, `None` for invalid windows),
/// mirroring `KmerTableSet.fillKmers`.
fn fill_kmers(bases: &[u8], k: usize) -> Vec<Option<Kmer>> {
    let mut out = Vec::with_capacity(bases.len().saturating_sub(k - 1));
    let mut kmer = Kmer::new(k);
    let mut len = 0usize;
    let min = k - 1;
    for (i, &b) in bases.iter().enumerate() {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
        if i >= min {
            if len >= k {
                out.push(Some(kmer));
            } else {
                out.push(None);
            }
        }
    }
    out
}

/// Fill the per-window canonical k-mer counts of a read.
fn fill_counts(kmers: &[Option<Kmer>], table: &RefineTable) -> Vec<i64> {
    kmers
        .iter()
        .map(|k| {
            if let Some(kmer) = k {
                raw_count(kmer, table)
            } else {
                0
            }
        })
        .collect()
}

/// Raw table count, mirroring Java `getCount`: -1 when the k-mer is absent.
fn raw_count(kmer: &Kmer, table: &RefineTable) -> i64 {
    let c = table.get_count(kmer);
    if c == 0 {
        -1
    } else {
        c as i64
    }
}

/// `KmerTableSet.regenerateCounts`: recompute window counts starting at `ca`
/// after a base change, resetting at undefined bases (count 0 for invalid
/// windows, raw -1 for absent k-mers).
fn regenerate_counts(bases: &[u8], counts: &mut [i64], table: &RefineTable, k: usize, ca: usize) {
    let b = ca + k - 1;
    let lim = bases.len().min(b + k + 1);
    let mut kmer = Kmer::new(k);
    let mut len = 0usize;
    for (j, &base) in bases[ca..lim].iter().enumerate() {
        let i = ca + j;
        if base_defined(base) {
            kmer.push_right(base_code(base));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
        if i >= b {
            let idx = i + 1 - k;
            if len >= k {
                counts[idx] = raw_count(&kmer, table);
            } else {
                counts[idx] = 0;
            }
        }
    }
}

/// Error-correction statistics (thread-local equivalents in BBTools).
#[derive(Debug, Default, Clone)]
pub struct ErrorTracker {
    pub suspected: usize,
    pub detected_reassemble: usize,
    pub corrected_reassemble_inner: usize,
    pub corrected_reassemble_outer: usize,
    pub rollback: bool,
}

impl ErrorTracker {
    fn corrected(&self) -> usize {
        self.corrected_reassemble_inner + self.corrected_reassemble_outer
    }
}

/// `isSimilar`: two k-mer depths are similar within absolute/fractional tolerances.
fn is_similar(a: i64, b: i64, opts: &RefineOptions) -> bool {
    let min = a.min(b);
    let max = a.max(b);
    let dif = max - min;
    (dif as f32) < opts.path_similarity_constant as f32
        || (dif as f32) < (max as f32) * opts.path_similarity_fraction
}

/// `isError(high, low)` (errorPath=1).
fn is_error2(high: i64, low: i64, opts: &RefineOptions) -> bool {
    let em1 = opts.error_mult1;
    (low as f32) * em1 < high as f32
        || (low <= opts.error_lower_const as i64
            && (high as f32) >= (opts.min_count_correct as f32).max(low as f32 * opts.error_mult2))
}

/// `isError(high, low, q)` (errorPath=1, quality-weighted).
fn is_error3(high: i64, low: i64, q: u8, opts: &RefineOptions) -> bool {
    let em1 = opts.error_mult1 * (1.0 + q as f32 * opts.error_mult_q_factor);
    (low as f32) * em1 < high as f32
        || (low <= opts.error_lower_const as i64
            && (high as f32) >= (opts.min_count_correct as f32).max(low as f32 * opts.error_mult2))
}

/// `isErrorBidirectional`.
fn is_error_bidirectional(a: i64, b: i64, qa: u8, qb: u8, opts: &RefineOptions) -> bool {
    if a >= b {
        is_error3(a, b, qb, opts)
    } else {
        is_error3(b, a, qa, opts)
    }
}

/// `isSubstitution`: isolated 1bp substitution candidate.
fn is_substitution(
    ca: usize,
    error_extension: usize,
    qb: u8,
    counts: &[i64],
    k: usize,
    opts: &RefineOptions,
) -> bool {
    let cb = ca + 1;
    let a_count = counts[ca];
    let b_count = counts[cb];
    if is_error3(a_count, b_count, qb, opts)
        && similar_range(
            a_count,
            ca as isize - error_extension as isize,
            ca as isize - 1,
            counts,
            opts,
        )
        && error_range(a_count, ca + 2, ca + k, counts, opts)
    {
        let cc = ca + k;
        let cd = cc + 1;
        if cd < counts.len() {
            let c_count = counts[cc];
            let d_count = counts[cd];
            is_error2(a_count, d_count, opts) || is_error3(d_count, c_count, qb, opts)
        } else {
            true
        }
    } else {
        false
    }
}

fn similar_range(a: i64, loc1: isize, loc2: isize, counts: &[i64], opts: &RefineOptions) -> bool {
    if loc2 < 0 {
        // Java clamps loc2 to -1 and the loop body never runs (empty range).
        return true;
    }
    let lo = loc1.max(0) as usize;
    let hi = (loc2 as usize).min(counts.len() - 1);
    if lo > hi {
        return true;
    }
    counts[lo..=hi].iter().all(|&c| is_similar(a, c, opts))
}

fn error_range(a: i64, loc1: usize, loc2: usize, counts: &[i64], opts: &RefineOptions) -> bool {
    let hi = loc2.min(counts.len() - 1);
    if loc1 > hi {
        return true;
    }
    counts[loc1..=hi].iter().all(|&c| is_error2(a, c, opts))
}

/// `countErrors`: count error positions, skipping `k` after each hit.
fn count_errors(counts: &[i64], quals: Option<&[u8]>, k: usize, opts: &RefineOptions) -> usize {
    let mut possible = 0usize;
    let mut i = 1usize;
    while i < counts.len() {
        let (a, b) = (counts[i - 1], counts[i]);
        let error = match quals {
            Some(q) => is_error_bidirectional(a, b, q[i - 1], q[i + k - 1], opts),
            None => is_error_bidirectional(a, b, 20, 20, opts),
        };
        if error {
            possible += 1;
            i += k;
        } else {
            i += 1;
        }
    }
    possible
}

/// `hasErrorsFast`: sampled k-mer depth screen for likely errors.
fn has_errors_fast(kmers: &[Option<Kmer>], table: &RefineTable, opts: &RefineOptions) -> bool {
    if kmers.is_empty() {
        return false;
    }
    let incr = (opts.k / 2).clamp(1, 9);
    let mcc = opts.min_count_correct as i64;
    let mut prev = -1i64;
    let mut i = 0usize;
    while i < kmers.len() {
        let count = match &kmers[i] {
            Some(kmer) => raw_count(kmer, table),
            None => return true,
        };
        let min = count.min(prev);
        let max = count.max(prev);
        if count < mcc || (i > 0 && is_error2(max + 1, min - 1, opts)) {
            return true;
        }
        prev = count;
        i += incr;
    }
    if let Some(kmer) = kmers.last() {
        let count = match kmer {
            Some(kmer) => raw_count(kmer, table),
            None => return true,
        };
        let min = count.min(prev);
        let max = count.max(prev);
        return count < mcc || is_error2(max + 1, min - 1, opts);
    }
    false
}

/// `isJunction(max, second)` with branch-resolution thresholds.
pub(crate) fn is_junction(max: u32, second: u32, opts: &RefineOptions) -> bool {
    if second < 1
        || (second as f32) * opts.branch_mult1 < max as f32
        || (second <= opts.branch_lower_const as u32
            && (max as f32)
                >= (opts.min_count_extend as f32).max(second as f32 * opts.branch_mult2))
    {
        return false;
    }
    true
}

/// Extends a sequence to the right by at most `distance` bases, mirroring
/// `Tadpole1.extendToRight2`. Returns the number of bases added.
#[allow(clippy::too_many_arguments)]
fn extend_to_right2(
    bases: &mut Vec<u8>,
    table: &RefineTable,
    opts: &RefineOptions,
    distance: usize,
    include_junction_base: bool,
    use_left: bool,
) -> usize {
    let k = opts.k;
    let initial = bases.len();
    if initial < k {
        return 0;
    }
    // Build the rightmost k-mer.
    let mut kmer = Kmer::new(k);
    let mut rkmer = Kmer::new(k);
    let mut len = 0usize;
    for &b in &bases[initial - k..initial] {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            rkmer.push_left(base_comp_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
            rkmer.reset();
        }
    }
    if len < k {
        return 0;
    }
    let count = table.get_count(&kmer);
    if count < opts.min_count_seed as u32 {
        return 0;
    }

    let mut left_max_pos = 0usize;
    let mut left_max = opts.min_count_extend as u32;
    let mut left_second = 0u32;
    if use_left {
        let lc = table.fill_left_counts(&kmer);
        left_max_pos = argmax2(&lc, &mut left_max);
        left_second = lc[second_highest_position(&lc)];
    }

    let rc = table.fill_right_counts(&kmer);
    let mut right_max = 0u32;
    let mut right_max_pos = argmax2(&rc, &mut right_max);
    let mut right_second_pos = second_highest_position(&rc);
    let mut right_second = rc[right_second_pos];

    if right_max < opts.min_count_extend as u32 {
        return 0;
    }
    if is_junction(right_max, right_second, opts)
        || (use_left && is_junction(left_max, left_second, opts))
    {
        return 0;
    }

    let max_len = initial.saturating_add(distance);
    let mut added = 0usize;
    // Tadpole1 (k<=31) appends the junction base when the forward k-mer is
    // the canonical maximum; Tadpole2 (k>31) canonicalizes to the minimum,
    // so the condition flips.
    let canonical_is_rc = k > 31;
    while bases.len() < max_len {
        let b = right_max_pos as u8;
        let x = right_max_pos as u8;
        let x2 = 3 - x;
        let evicted = kmer.base_at(k - 1);
        kmer.push_right(x);
        rkmer.push_left(x2);

        if use_left {
            let lc = table.fill_left_counts(&kmer);
            left_max_pos = argmax2(&lc, &mut left_max);
            left_second = lc[second_highest_position(&lc)];
        }
        let rc = table.fill_right_counts(&kmer);
        right_max_pos = argmax2(&rc, &mut right_max);
        right_second_pos = second_highest_position(&rc);
        right_second = rc[right_second_pos];

        let junc_r = is_junction(right_max, right_second, opts);
        let junc_l = use_left && is_junction(left_max, left_second, opts);
        // Tadpole2 (k>31) appends the junction base when the k-mer's
        // canonical orientation is the forward one (`key()==array1` in
        // BBTools; the reverse-complement key is the other branch).
        let kmer_is_rc = kmer.cmp_bases(&rkmer).is_lt();
        if junc_r || junc_l {
            if include_junction_base
                && if canonical_is_rc {
                    kmer_is_rc
                } else {
                    !kmer_is_rc
                }
            {
                bases.push(number_to_base(b));
                added += 1;
            }
            break;
        }
        if use_left && left_max_pos != evicted as usize {
            if include_junction_base
                && if canonical_is_rc {
                    kmer_is_rc
                } else {
                    !kmer_is_rc
                }
            {
                bases.push(number_to_base(b));
                added += 1;
            }
            break;
        }
        bases.push(number_to_base(b));
        added += 1;
        if right_max < opts.min_count_extend as u32 {
            break;
        }
    }
    added
}

pub(crate) fn number_to_base(n: u8) -> u8 {
    match n {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => b'N',
    }
}

pub(crate) fn argmax2(a: &[u32; 4], max: &mut u32) -> usize {
    let mut pos = 0usize;
    *max = a[0];
    for (i, &x) in a.iter().enumerate().skip(1) {
        if x > *max {
            *max = x;
            pos = i;
        }
    }
    pos
}

pub(crate) fn second_highest_position(a: &[u32; 4]) -> usize {
    let (mut p, mut p2) = if a[0] >= a[1] { (0, 1) } else { (1, 0) };
    for i in 2..a.len() {
        let x = a[i];
        if x > a[p2] {
            if x >= a[p] {
                p2 = p;
                p = i;
            } else {
                p2 = i;
            }
        }
    }
    p2
}

/// `isJunk`: read cannot be used for assembly.
pub fn is_junk(bases: &[u8], table: &RefineTable, opts: &RefineOptions, paired: bool) -> bool {
    let k = opts.k;
    let blen = bases.len();
    if blen < k {
        return true;
    }
    let mut kmer = Kmer::new(k);
    let mut len = 0usize;
    for &b in &bases[..k] {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
    }
    if len >= k {
        let lc = table.fill_left_counts(&kmer);
        let max_pos = argmax2(&lc, &mut 0);
        if lc[max_pos] > 0 {
            return false;
        }
    }
    let mut max_depth = 0u32;
    for &b in &bases[k..] {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
        if len < k {
            continue;
        }
        {
            let depth = table.get_count(&kmer);
            if depth > max_depth {
                max_depth = depth;
                if max_depth > 1 && (!paired || max_depth > 2) {
                    return false;
                }
            }
        }
    }
    if len >= k && !paired {
        let rc = table.fill_right_counts(&kmer);
        let max_pos = argmax2(&rc, &mut 0);
        if rc[max_pos] > 0 {
            return false;
        }
    }
    true
}

/// `hasKmersAtOrBelow`: does the read have enough low-depth k-mers to toss?
pub fn has_kmers_at_or_below(
    bases: &[u8],
    table: &RefineTable,
    opts: &RefineOptions,
    too_low: u32,
    fraction: f32,
) -> bool {
    let k = opts.k;
    let blen = bases.len();
    if blen < k {
        return true;
    }
    let mut kmer = Kmer::new(k);
    let limit = ((blen - k + 1) as f32 * fraction).round().max(1.0) as usize;
    let mut valid = 0usize;
    let mut invalid = 0usize;
    let mut len = 0usize;
    for &b in bases.iter() {
        if base_defined(b) {
            kmer.push_right(base_code(b));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
        if len >= k {
            let depth = table.get_count(&kmer);
            if depth > too_low {
                valid += 1;
            } else {
                invalid += 1;
                if invalid >= limit {
                    return true;
                }
            }
        }
    }
    let limit2 = ((valid + invalid) as f32 * fraction).round().max(1.0) as usize;
    valid < 1 || invalid >= limit2
}

/// `Read.expectedErrors` (phred qualities, countUndefined=true).
pub fn expected_errors(quals: &[u8]) -> f32 {
    quals.iter().map(|&q| prob_error(q)).sum()
}

/// Error-corrects one read in place (reassemble-only path), mirroring
/// `Tadpole1.errorCorrect` + `Tadpole.reassemble`. Returns corrections applied.
pub fn error_correct(
    bases: &mut Vec<u8>,
    quals: &mut Vec<u8>,
    table: &RefineTable,
    opts: &RefineOptions,
    tracker: &mut ErrorTracker,
) -> usize {
    tracker.suspected = 0;
    tracker.detected_reassemble = 0;
    tracker.corrected_reassemble_inner = 0;
    tracker.corrected_reassemble_outer = 0;
    tracker.rollback = false;

    let kmers = fill_kmers(bases, opts.k);
    let valid = kmers.len();
    if valid < 2 {
        return 0;
    }
    let has_undefined = bases.iter().any(|&b| !base_defined(b));
    if !has_undefined && !has_errors_fast(&kmers, table, opts) {
        return 0;
    }
    let mut counts = fill_counts(&kmers, table);
    // FASTA input has no qualities; BBTools substitutes a fixed quality 20.
    let qs = if quals.is_empty() {
        None
    } else {
        Some(quals.as_slice())
    };
    let possible_errors = count_errors(&counts, qs, opts.k, opts);
    tracker.suspected = possible_errors;
    let expected = expected_errors(quals);
    let counts0 = counts.clone();
    let bases0 = bases.clone();
    let quals0 = quals.clone();

    let corrected = reassemble(
        bases,
        quals,
        table,
        opts,
        &mut counts,
        tracker,
        opts.error_extension_reassemble,
    );
    debug_assert_eq!(corrected, tracker.corrected());

    if opts.ecc_rollback && (tracker.corrected() > 0 || tracker.rollback) {
        if !tracker.rollback && tracker.corrected() > 3 {
            let mult = (0.5f32 * (0.5 + 0.01 * bases.len() as f32)).max(1.0);
            let ce = count_errors(&counts, Some(quals), opts.k, opts);
            let c1 = ce > 0 && tracker.corrected() as f32 > mult + expected;
            let c2 = tracker.corrected() as f32 > 2.5 * mult + expected;
            if c1 || c2 {
                tracker.rollback = true;
            }
        }
        if !tracker.rollback {
            for i in 0..counts.len() {
                // Java clamps both sides to 0 before the rollback comparison.
                let a = counts0[i].max(0);
                let b = counts[i].max(0);
                if b < a - 1 && !is_similar(a, b, opts) {
                    tracker.rollback = true;
                }
            }
        }
        if tracker.rollback {
            *bases = bases0;
            *quals = quals0;
            tracker.corrected_reassemble_inner = 0;
            tracker.corrected_reassemble_outer = 0;
            return 0;
        }
    }
    tracker.corrected()
}

/// `reassemble`: multi-pass local reassembly error correction.
#[allow(clippy::too_many_arguments)]
fn reassemble(
    bases: &mut [u8],
    quals: &mut [u8],
    table: &RefineTable,
    opts: &RefineOptions,
    counts: &mut Vec<i64>,
    tracker: &mut ErrorTracker,
    error_extension: usize,
) -> usize {
    if bases.len() < opts.k + 1 + opts.dead_zone {
        return 0;
    }
    let mut corrected = 0usize;
    let mut corrected_incr;
    let mut detected_incr;
    let mut uncorrected;
    let detected0 = tracker.detected_reassemble;
    corrected_incr = reassemble_pass(bases, quals, table, opts, counts, tracker, error_extension);
    corrected += corrected_incr;
    detected_incr = tracker.detected_reassemble - detected0;
    uncorrected = detected_incr.saturating_sub(corrected_incr);
    let mut passes = 1usize;
    while passes < 6 && corrected_incr > 0 && uncorrected > 0 {
        tracker.detected_reassemble -= uncorrected;
        let detected0 = tracker.detected_reassemble;
        corrected_incr =
            reassemble_pass(bases, quals, table, opts, counts, tracker, error_extension);
        corrected += corrected_incr;
        detected_incr = tracker.detected_reassemble - detected0;
        uncorrected = detected_incr.saturating_sub(corrected_incr);
        passes += 1;
    }
    corrected
}

/// `reassemble_pass`: forward + reverse passes, window filtering, consensus.
#[allow(clippy::too_many_arguments)]
fn reassemble_pass(
    bases: &mut [u8],
    quals: &mut [u8],
    table: &RefineTable,
    opts: &RefineOptions,
    counts: &mut Vec<i64>,
    tracker: &mut ErrorTracker,
    error_extension: usize,
) -> usize {
    if bases.len() < opts.k + 1 + opts.dead_zone {
        return 0;
    }
    let mut from_left = bases.to_vec();
    let mut from_right = bases.to_vec();
    let mut counts2 = counts.clone();
    reassemble_inner(
        &mut from_left,
        quals,
        table,
        opts,
        &mut counts2,
        error_extension,
    );

    from_right = rev_comp(&from_right).collect();
    let qr: Vec<u8> = quals.iter().rev().copied().collect();
    counts2 = counts.clone();
    counts2.reverse();
    reassemble_inner(
        &mut from_right,
        &qr,
        table,
        opts,
        &mut counts2,
        error_extension,
    );
    from_right = rev_comp(&from_right).collect();

    let mut corrected_inner = 0usize;
    let mut corrected_outer = 0usize;
    let mut detected_inner = 0usize;
    let mut detected_outer = 0usize;
    let mut rollback = false;
    for i in 0..bases.len() {
        let a = bases[i];
        let b = from_left[i];
        let c = from_right[i];
        if a != b || a != c {
            if b == c {
                detected_inner += 1;
            } else {
                detected_outer += 1;
                if a != b && a != c {
                    rollback = true;
                }
            }
        }
        if b == a {
            from_left[i] = 0;
        }
        if c == a {
            from_right[i] = 0;
        }
    }
    let detected = detected_inner + detected_outer;
    tracker.detected_reassemble += detected;
    if rollback || detected == 0 {
        return 0;
    }

    clear_window2(&mut from_left, quals, opts);
    // Java clears fromRight while it is in reversed orientation with the
    // reversed qualities; clearing the forward-oriented copy with reversed
    // qualities would mis-weight the two read ends.
    from_right.reverse();
    clear_window2(&mut from_right, &qr, opts);
    from_right.reverse();

    for i in 0..bases.len() {
        let a = bases[i];
        let b = from_left[i];
        let c = from_right[i];
        let mut d = a;
        if b == 0 && c == 0 {
            // nothing
        } else if b == c {
            d = b;
        } else if b == 0 {
            d = c;
        } else if c == 0 {
            d = b;
        } else if b != c {
            // keep a
        }
        if opts.ecc_require_bidirectional && b != c && i >= opts.k && i < bases.len() - opts.k {
            d = a;
        }
        if d != a {
            let mut q = if quals.is_empty() { 30 } else { quals[i] };
            if b == c {
                corrected_inner += 1;
                q = q.saturating_add(8).clamp(24, 32);
            } else {
                corrected_outer += 1;
                q = q.saturating_add(4).clamp(20, 28);
            }
            if !rollback {
                bases[i] = d;
                if !quals.is_empty() {
                    quals[i] = q;
                }
            }
        }
    }
    if rollback && corrected_inner + corrected_outer > 0 {
        tracker.rollback = true;
        return 0;
    }
    tracker.corrected_reassemble_inner += corrected_inner;
    tracker.corrected_reassemble_outer += corrected_outer;
    let corrected = corrected_inner + corrected_outer;
    if corrected > 0 {
        // Regenerate counts for all windows.
        let kmers = fill_kmers(bases, opts.k);
        *counts = fill_counts(&kmers, table);
    }
    corrected
}

/// `reassemble_inner`: per-position substitution detection and correction.
#[allow(clippy::too_many_arguments)]
fn reassemble_inner(
    bases: &mut [u8],
    quals: &[u8],
    table: &RefineTable,
    opts: &RefineOptions,
    counts: &mut [i64],
    error_extension: usize,
) -> usize {
    let k = opts.k;
    let length = bases.len();
    if length < k + 1 + opts.dead_zone {
        return 0;
    }
    let mut kmer = Kmer::new(k);
    let mut len = 0usize;
    let mut corrected = 0usize;
    let lim = length - opts.dead_zone - 1;
    for a in 0..lim {
        if base_defined(bases[a]) {
            kmer.push_right(base_code(bases[a]));
            len += 1;
        } else {
            len = 0;
            kmer.reset();
        }
        if len >= k {
            let b = a + 1;
            // len>=k implies a+1>=k, so this cannot underflow.
            let ca = a + 1 - k;
            let a_count = counts[ca];
            let qb = if quals.is_empty() { 20 } else { quals[b] };
            if is_substitution(ca, error_extension, qb, counts, k, opts) {
                let rc = table.fill_right_counts(&kmer);
                let right_max_pos = argmax2(&rc, &mut 0);
                let right_max = rc[right_max_pos];
                let right_second_pos = second_highest_position(&rc);
                let right_second = rc[right_second_pos];
                let base = bases[b];
                // Java `baseToNumber` is -1 for N, so an N never matches the
                // preferred extension and always goes through the correction.
                let num = if base_defined(base) {
                    base_code(base) as i64
                } else {
                    -1
                };
                if right_max >= opts.min_count_extend as u32 {
                    // BBTools compares the base code to the *count* here
                    // (`if(num==rightMax)`), not to the position index; the
                    // base is treated as already-correct when they coincide.
                    if num == right_max as i64 {
                    } else if (is_error3(right_max as i64, right_second as i64, qb, opts)
                        || !is_junction(right_max, right_second, opts))
                        && is_similar(a_count, right_max as i64, opts)
                    {
                        bases[b] = number_to_base(right_max_pos as u8);
                        corrected += 1;
                        // Regenerate counts for windows ca+1..=ca+k (those
                        // containing the changed base at ca+k).
                        regenerate_counts(bases, counts, table, k, ca);
                    }
                }
            }
        }
    }
    corrected
}

/// `clearWindow2`: sliding-window quality filter over correction candidates.
fn clear_window2(bb: &mut [u8], quals: &[u8], opts: &RefineOptions) -> usize {
    let len = bb.len();
    let window = opts.window_len as isize;
    let mut cleared = 0usize;
    let mut count = 0usize;
    let mut qsum = 0usize;
    for (i, prev) in (0..len as isize).zip((-window)..) {
        let b = bb[i as usize];
        if b != 0 && (quals.is_empty() || quals[i as usize] > 0) {
            count += 1;
            if !quals.is_empty() {
                qsum += quals[i as usize] as usize;
            }
            if count > opts.window_count || qsum > opts.window_qual_sum {
                let start = (i - window).max(0) as usize;
                for b in &mut bb[start..] {
                    if *b != 0 {
                        *b = 0;
                        cleared += 1;
                    }
                }
                return cleared;
            }
        }
        if prev >= 0 && bb[prev as usize] > 0 && (quals.is_empty() || quals[prev as usize] > 0) {
            count -= 1;
            if !quals.is_empty() {
                qsum -= quals[prev as usize] as usize;
            }
        }
    }
    cleared
}

/// Extends one read in place (both ends), mirroring `processRead` extension.
pub fn extend_read(
    bases: &mut Vec<u8>,
    quals: &mut Vec<u8>,
    table: &RefineTable,
    opts: &RefineOptions,
    numeric_id: u64,
) -> usize {
    let mut extension_right = 0usize;
    let mut extension_left = 0usize;
    if opts.extend_right > 0 {
        extension_right = extend_read_one_side(bases, quals, table, opts, opts.extend_right);
    }
    if opts.extend_left > 0 {
        // Reverse-complement, extend, reverse back.
        let mut rc: Vec<u8> = rev_comp(bases).collect();
        let mut rq: Vec<u8> = quals.iter().rev().copied().collect();
        extension_left = extend_read_one_side(&mut rc, &mut rq, table, opts, opts.extend_left);
        *bases = rev_comp(&rc).collect();
        *quals = rq.iter().rev().copied().collect();
    }
    let mut extension = extension_right + extension_left;
    if opts.extension_rollback > 0 {
        let mut left_mod = 0usize;
        let mut right_mod = 0usize;
        // `+1` guards against `% 0`; saturating guards against the rollback
        // value `usize::MAX` overflowing on `+1`.
        let rollback = (opts.extension_rollback as u64).saturating_add(1);
        if extension_left > 0 && extension_left < opts.extend_left {
            left_mod = extension_left.min((numeric_id % rollback) as usize);
            extension_left -= left_mod;
        }
        if extension_right > 0 && extension_right < opts.extend_right {
            right_mod = extension_right.min((numeric_id % rollback) as usize);
            extension_right -= right_mod;
        }
        if left_mod > 0 || right_mod > 0 {
            // Trim left_mod bases from the 5' end and right_mod from 3'.
            let keep_from = left_mod.min(bases.len());
            let keep_to = bases.len().saturating_sub(right_mod);
            if keep_from < keep_to {
                *bases = bases[keep_from..keep_to].to_vec();
                if !quals.is_empty() {
                    *quals = quals[keep_from..keep_to].to_vec();
                }
            } else {
                bases.clear();
                quals.clear();
            }
        }
        extension = extension_left + extension_right;
    }
    extension
}

/// Extends one read's 3' end in place, mirroring `Tadpole.extendToRight2`
/// as called by BBMerge (`extendAndMerge` / `extendRead`): right junction
/// only, no left-branch check.
pub fn extend_read_right(
    bases: &mut Vec<u8>,
    quals: &mut Vec<u8>,
    table: &RefineTable,
    opts: &RefineOptions,
    distance: usize,
) -> usize {
    // BBMerge `extendRead` calls `extendToRight2(..., false)` (no junction
    // base); the standalone `fq extend` path keeps the junction base.
    let initial = bases.len();
    if initial < opts.k {
        return 0;
    }
    let added = extend_to_right2(bases, table, opts, distance, false, false);
    if added > 0 && !quals.is_empty() {
        quals.resize(bases.len(), 30);
    }
    added
}

/// Extends one read end (3' after any RC flip) by up to `distance` bases.
fn extend_read_one_side(
    bases: &mut Vec<u8>,
    quals: &mut Vec<u8>,
    table: &RefineTable,
    opts: &RefineOptions,
    distance: usize,
) -> usize {
    let initial = bases.len();
    if initial < opts.k {
        return 0;
    }
    // BBTools never initializes the left-counts buffer for read extension
    // (`ExtendThread.leftCounts` stays null), so only the right junction is
    // considered; the left-branch check is disabled.
    let added = extend_to_right2(bases, table, opts, distance, true, false);
    if added > 0 && !quals.is_empty() {
        quals.resize(bases.len(), 30);
    }
    added
}

/// Per-read processing outcome counters (subset of tadpole.sh stats).
#[derive(Debug, Default, Clone)]
pub struct RefineStats {
    pub reads_in: u64,
    pub bases_in: u64,
    pub bases_extended: u64,
    pub reads_extended: u64,
    pub reads_corrected: u64,
    pub bases_corrected: u64,
    pub reads_detected: u64,
    pub bases_detected: u64,
    pub reads_fully_corrected: u64,
    pub reads_discarded: u64,
    pub bases_discarded: u64,
    pub rollbacks: u64,
}

/// Main per-read processing pipeline, mirroring `ExtendThread.processRead`.
#[allow(clippy::too_many_arguments)]
pub fn process_read(
    bases: &mut Vec<u8>,
    quals: &mut Vec<u8>,
    table: &RefineTable,
    opts: &RefineOptions,
    stats: &mut RefineStats,
    numeric_id: u64,
    mate: Option<usize>,
    discard_mate: &mut bool,
) -> bool {
    let initial_len = bases.len();
    let mut tracker = ErrorTracker::default();
    if opts.ecc && (opts.toss_uncorrectable || opts.ecc_rollback) {
        let corrected = error_correct(bases, quals, table, opts, &mut tracker);
        if tracker.rollback {
            stats.rollbacks += 1;
        }
        let detected = tracker.detected_reassemble;
        if detected > 0 {
            stats.reads_detected += 1;
            stats.bases_detected += detected as u64;
            if corrected > 0 {
                stats.reads_corrected += 1;
                stats.bases_corrected += corrected as u64;
            }
            if corrected == detected
                || (corrected > 0 && count_errors_from(bases, quals, table, opts) == 0)
            {
                stats.reads_fully_corrected += 1;
            } else if opts.toss_uncorrectable {
                if mate.is_some() && !opts.require_both_bad {
                    *discard_mate = true;
                }
                return true; // discard this read
            }
        }
    }

    if opts.toss_junk && is_junk(bases, table, opts, mate.is_some_and(|m| m >= opts.k)) {
        return true;
    }

    if opts.toss_depth >= 0
        && has_kmers_at_or_below(
            bases,
            table,
            opts,
            opts.toss_depth as u32,
            opts.low_depth_discard_fraction,
        )
    {
        if mate.is_some() && !opts.require_both_bad {
            *discard_mate = true;
        }
        return true;
    }

    if opts.extend_right > 0 || opts.extend_left > 0 {
        let ext = extend_read(bases, quals, table, opts, numeric_id);
        if ext > 0 {
            stats.bases_extended += ext as u64;
            stats.reads_extended += 1;
        }
    }
    stats.bases_in += initial_len as u64;
    false
}

fn count_errors_from(
    bases: &[u8],
    quals: &[u8],
    table: &RefineTable,
    opts: &RefineOptions,
) -> usize {
    let kmers = fill_kmers(bases, opts.k);
    let counts = fill_counts(&kmers, table);
    let qs = if quals.is_empty() { None } else { Some(quals) };
    count_errors(&counts, qs, opts.k, opts)
}

/// Runs the tadpole correct/extend/discard pipeline over FASTQ input
/// (1 interleaved file or 2 paired files), mirroring `tadpole.sh`.
pub fn run<W: Write>(infiles: &[String], out: &mut W, opts: &RefineOptions) -> Result<RefineStats> {
    anyhow::ensure!(
        opts.k >= 1,
        "k-mer length must be at least 1, got {}",
        opts.k
    );
    anyhow::ensure!(
        opts.k <= key::Kmer::MAX_K,
        "k-mer length must be at most {}, got {}",
        key::Kmer::MAX_K,
        opts.k
    );
    // Pass 1: read all records into memory, canonicalizing qualities.
    let mut records: Vec<SeqRecord> = Vec::new();
    let mut reader1 = SeqReader::new(&infiles[0])?;
    let mut reader2 = if infiles.len() > 1 {
        Some(SeqReader::new(&infiles[1])?)
    } else {
        None
    };
    let mut rec = SeqRecord::new();
    loop {
        if !reader1.read_record(&mut rec)? {
            break;
        }
        canonicalize_quality(&mut rec);
        records.push(rec.clone());
        if let Some(r) = reader2.as_mut() {
            if !r.read_record(&mut rec)? {
                anyhow::bail!("unpaired trailing read in {}", infiles[0]);
            }
            canonicalize_quality(&mut rec);
            records.push(rec.clone());
        } else if !reader1.read_record(&mut rec)? {
            anyhow::bail!("unpaired trailing read in {}", infiles[0]);
        } else {
            canonicalize_quality(&mut rec);
            records.push(rec.clone());
        }
    }

    // Pass 2: count k-mers from the canonicalized (phred) qualities.
    let reads: Vec<(Vec<u8>, Vec<u8>)> = records
        .iter()
        .map(|r| {
            (
                r.sequence().to_vec(),
                to_phred(r.sequence(), r.quality_scores()),
            )
        })
        .collect();
    let table = RefineTable::build(&reads, opts.k, opts.min_prob);

    // Pass 3: process pairs and write surviving reads.
    let mut stats = RefineStats {
        reads_in: records.len() as u64,
        ..Default::default()
    };
    let mut i = 0usize;
    while i < records.len() {
        let r1 = records[i].clone();
        let r2 = if i + 1 < records.len() {
            Some(records[i + 1].clone())
        } else {
            None
        };
        // BBTools assigns one numeric ID per pair (both mates share it).
        let id = (i / 2) as u64;
        let mut bases1 = r1.sequence().to_vec();
        let mut quals1 = to_phred(&bases1, r1.quality_scores());
        let mut bases2 = r2
            .as_ref()
            .map(|r| r.sequence().to_vec())
            .unwrap_or_default();
        let mut quals2 = r2
            .as_ref()
            .map(|r| to_phred(r.sequence(), r.quality_scores()))
            .unwrap_or_default();
        let mut discard_mate = false;
        let d1 = process_read(
            &mut bases1,
            &mut quals1,
            &table,
            opts,
            &mut stats,
            id,
            r2.as_ref().map(|r| r.sequence().len()),
            &mut discard_mate,
        );
        let d2 = if discard_mate {
            true
        } else if r2.is_some() {
            let mate_len = bases1.len(); // r1 length after its own processing
            process_read(
                &mut bases2,
                &mut quals2,
                &table,
                opts,
                &mut stats,
                id,
                Some(mate_len),
                &mut discard_mate,
            )
        } else {
            true
        };
        // Either read's processing may discard the other as its mate
        // (tossdepth / tossuncorrectable without requireBothBad).
        let d1 = d1 || discard_mate;
        // A pair is dropped only when both reads are discarded; otherwise
        // both are written (discarded mates keep their processed state).
        if d1 && d2 {
            stats.reads_discarded += 1 + r2.is_some() as u64;
            stats.bases_discarded += bases1.len() as u64 + bases2.len() as u64;
        } else {
            write_record(out, &r1, &bases1, &from_phred(&quals1))?;
            if let Some(r) = r2.as_ref() {
                write_record(out, r, &bases2, &from_phred(&quals2))?;
            }
        }
        i += 2;
    }
    Ok(stats)
}

fn write_record<W: Write>(w: &mut W, rec: &SeqRecord, seq: &[u8], qual: &[u8]) -> Result<()> {
    let header = if rec.comment().is_empty() {
        rec.name().to_string()
    } else {
        format!("{} {}", rec.name(), rec.comment())
    };
    if qual.is_empty() {
        pgr::libs::fmt::fq::write_fa(w, &header, seq)?;
    } else {
        pgr::libs::fmt::fq::write_fq(w, &header, seq, qual)?;
    }
    Ok(())
}

#[test]
fn junk_detects_short_read() {
    let opts = RefineOptions::default();
    let table = RefineTable::build(&[], opts.k, opts.min_prob);
    assert!(is_junk(b"ACGT".as_ref(), &table, &opts, false));
}

#[test]
fn read36_left_extension_matches_golden() {
    // Reproduce the golden `fq extend` run on the committed subset and
    // check that read 36's left extension matches BBTools.
    let infile = "tests/bbtools/Lambda/golden/ecco_sub.fq.gz";
    let mut reader = SeqReader::new(infile).unwrap();
    let mut rec = SeqRecord::new();
    let mut reads: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    while reader.read_record(&mut rec).unwrap() {
        let seq = rec.sequence().to_vec();
        let quals = to_phred(&seq, rec.quality_scores());
        reads.push((seq, quals));
    }
    let k = 62usize;
    let table = RefineTable::build(&reads, k, 0.5);
    let opts = RefineOptions {
        k,
        extend_left: 20,
        extend_right: 20,
        ..RefineOptions::default()
    };
    let r36 = &reads[35].0;
    let mut bases = r36.clone();
    let mut quals = vec![30; bases.len()];
    // The `run` pipeline assigns one numeric ID per pair; read 36 is in
    // pair 17 (0-based), which drives the extension-rollback trim.
    extend_read(&mut bases, &mut quals, &table, &opts, 17);
    let seq = String::from_utf8_lossy(&bases).into_owned();
    // Golden: input + GTGGAA on the left + GAAGGCATTAACGCCTCTGC right.
    let golden = format!(
        "GTGGAA{}{}",
        String::from_utf8_lossy(r36),
        "GAAGGCATTAACGCCTCTGC"
    );
    assert_eq!(seq, golden, "read36 extension mismatch");
}

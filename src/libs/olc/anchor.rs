//! Reliable anchor selection: per-base read coverage filtering.
//!
//! Mirrors the legacy `anchr anchors` flow (bbwrap perfectmode + basecov +
//! [lower, upper] coverage interval): reads are mapped back to the unitigs,
//! every position gets its depth, positions outside
//! `lower = max(mincov, (median - mscale*MAD)/lscale)` and
//! `upper = (median + mscale*MAD)*uscale` are excluded (low coverage =
//! errors, high coverage = repeats), and the remaining contiguous stretches
//! are the anchors — the reliable fragments fed to the OLC merge.

/// Coverage-window parameters for anchor selection (legacy defaults).
#[derive(Debug, Clone, Copy)]
pub struct AnchorOptions {
    /// Absolute floor for `lower` (legacy `mincov`).
    pub mincov: u32,
    /// Median absolute deviation multiplier (legacy `mscale`).
    pub mscale: f64,
    /// Lower-window divider (legacy `lscale`).
    pub lscale: f64,
    /// Upper-window multiplier (legacy `uscale`).
    pub uscale: f64,
    /// Minimum anchor length in bases.
    pub min_len: usize,
}

impl Default for AnchorOptions {
    fn default() -> Self {
        Self {
            mincov: 5,
            mscale: 3.0,
            lscale: 3.0,
            uscale: 2.0,
            min_len: 500,
        }
    }
}

/// One perfect-match alignment: reference index, 1-based inclusive interval.
pub type Alignment = (usize, usize, usize);

/// Accumulates per-base coverage for every reference from perfect-match
/// alignments (one sweep-line diff array per reference).
pub fn coverage_from_alignments(lens: &[usize], aligns: &[Alignment]) -> Vec<Vec<u32>> {
    let mut covs = Vec::with_capacity(lens.len());
    for &len in lens {
        covs.push(vec![0i64; len + 1]); // 1-based positions, [0] unused
    }
    for &(ri, a, b) in aligns {
        covs[ri][a] += 1;
        if b + 1 < covs[ri].len() {
            covs[ri][b + 1] -= 1;
        }
    }
    covs.into_iter()
        .map(|mut diff| {
            let mut cur = 0i64;
            for v in diff.iter_mut().skip(1) {
                cur += *v;
                *v = cur.max(0);
            }
            diff[0] = 0;
            diff.into_iter().map(|d| d as u32).collect()
        })
        .collect()
}

/// Coverage statistics for the anchor window (zero-coverage positions
/// included): the legacy formula with the caller's `AnchorOptions`.
#[derive(Debug, Clone, Copy)]
pub struct AnchorStats {
    /// Per-base depth median.
    pub median: u32,
    /// Median absolute deviation of the per-base depths.
    pub mad: u32,
    /// Lower coverage bound.
    pub lower: f64,
    /// Upper coverage bound.
    pub upper: f64,
}

/// Computes the coverage window statistics from the per-position depths of
/// all references.
pub fn anchor_stats(covs: &[Vec<u32>], opts: &AnchorOptions) -> AnchorStats {
    let mut all: Vec<u32> = Vec::new();
    for cov in covs {
        all.extend_from_slice(&cov[1..]);
    }
    if all.is_empty() {
        return AnchorStats {
            median: 0,
            mad: 0,
            lower: 0.0,
            upper: 0.0,
        };
    }
    all.sort_unstable();
    let median = all[all.len() / 2];
    let mad = {
        let mut devs: Vec<u32> = all.iter().map(|&x| x.abs_diff(median)).collect();
        devs.sort_unstable();
        devs[devs.len() / 2]
    };
    let lower = ((median as f64 - opts.mscale * mad as f64) / opts.lscale).max(opts.mincov as f64);
    let upper = (median as f64 + opts.mscale * mad as f64) * opts.uscale;
    AnchorStats {
        median,
        mad,
        lower,
        upper,
    }
}

/// Contiguous positions with `lower <= depth <= upper`, returned as
/// `(ref_index, start, end)` 1-based inclusive regions of length >= `min_len`.
/// Positions within `read_len / 2` of a contig end scale the window down
/// linearly (`lower/upper * 2*dist/read_len`, floored at `mincov`): reads
/// cannot fully cover the ends, so the expected coverage ramps up from the
/// edge (legacy `--keepedge` behavior).
pub fn anchor_regions(
    covs: &[Vec<u32>],
    opts: &AnchorOptions,
    lower: f64,
    upper: f64,
    read_len: usize,
) -> Vec<(usize, usize, usize)> {
    let half = read_len / 2;
    let mut regions = Vec::new();
    for (ri, cov) in covs.iter().enumerate() {
        let len = cov.len() - 1;
        let mut in_run = false;
        let mut start = 0usize;
        for (p, &depth) in cov.iter().enumerate().skip(1) {
            let d = depth as f64;
            let ok = if half > 0 && p <= half {
                let f = (p - 1) as f64 * 2.0 / read_len as f64;
                d >= (lower * f).max(opts.mincov as f64) && d <= upper * f
            } else if half > 0 && p > len.saturating_sub(half) {
                let f = (len - p) as f64 * 2.0 / read_len as f64;
                d >= (lower * f).max(opts.mincov as f64) && d <= upper * f
            } else {
                d >= lower && d <= upper
            };
            if ok && !in_run {
                in_run = true;
                start = p;
            } else if !ok && in_run {
                if p - start >= opts.min_len {
                    regions.push((ri, start, p - 1));
                }
                in_run = false;
            }
        }
        if in_run && cov.len() - start >= opts.min_len {
            regions.push((ri, start, cov.len() - 1));
        }
    }
    regions
}

/// Extracts the anchor sequences from the reference records.
pub fn extract_anchors(seqs: &[Vec<u8>], regions: &[(usize, usize, usize)]) -> Vec<Vec<u8>> {
    regions
        .iter()
        .map(|&(ri, a, b)| seqs[ri][a - 1..b].to_vec())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coverage_and_thresholds_basic() {
        let lens = [10usize, 8];
        let aligns = [(0usize, 1usize, 5usize), (0, 1, 5), (0, 3, 8), (1, 1, 4)];
        let covs = coverage_from_alignments(&lens, &aligns);
        assert_eq!(covs[0][1], 2);
        assert_eq!(covs[0][3], 3);
        assert_eq!(covs[0][6], 1);
        assert_eq!(covs[0][9], 0);
        assert_eq!(covs[1][1], 1);
        assert_eq!(covs[1][5], 0);
    }

    #[test]
    fn anchor_regions_filter_by_window() {
        // Depth 10 everywhere on ref 0; depth 10 except position 4-5 (depth 1).
        let mut covs = vec![vec![0u32; 11]; 1];
        covs[0][1..].fill(10);
        covs[0][4] = 1;
        covs[0][5] = 1;
        let opts = AnchorOptions {
            min_len: 2,
            ..Default::default()
        };
        // read_len = 2 -> half = 1: only the very first/last base is an edge.
        let regions = anchor_regions(&covs, &opts, 5.0, 20.0, 2);
        // Runs: 2-3 and 6-9 (base 1 and 10 are zero-width edges).
        assert_eq!(regions, vec![(0, 2, 3), (0, 6, 9)]);
    }

    /// The coverage ramp at contig ends is accepted by the scaled window
    /// (legacy `--keepedge`), instead of being clipped by the global lower.
    #[test]
    fn anchor_regions_keeps_edge_ramp() {
        // len 20; depth ramps 2..10 over the first 5 bases and 10..2 over
        // the last 5, 10 in the middle (read_len 10 -> half 5).
        let mut covs = vec![vec![0u32; 21]; 1];
        let ramp = [2u32, 4, 6, 8, 10];
        for (i, v) in ramp.iter().enumerate() {
            covs[0][i + 1] = *v;
            covs[0][20 - i] = *v;
        }
        covs[0][6..16].fill(10);
        let opts = AnchorOptions {
            min_len: 2,
            ..Default::default()
        };
        // lower = 5, upper = 20: mid-edge bases (depth 6-10) pass the scaled
        // window, the first/last (depth 2) are below mincov and stay cut.
        let regions = anchor_regions(&covs, &opts, 5.0, 20.0, 10);
        assert_eq!(regions, vec![(0, 3, 18)]);
    }
}

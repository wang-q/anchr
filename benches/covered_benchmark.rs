//! Benchmark: per-sequence coverage tiers for `anchr covered`.
//!
//! Two candidate implementations of the same workload (overlap intervals
//! stacked on one sequence, output the coverage layers):
//!
//! * `bump` — the vendored intspan-style `Coverage` (`new_len` +
//!   per-interval `bump`, then `max_tier` / `uniq_tiers`), the current
//!   anchr implementation.
//! * `sweep` — `pgr::libs::runlist` (`depth_at_least` / `depth_by_level`,
//!   one sweep over sorted start/end events), the pgr `rg coverage` core.
//!
//! Workload shapes mirror `pgr/benches/rg_merge_benchmark.rs`: disjoint
//! ranges and clusters of overlapping ranges. Both paths are asserted
//! identical (IntSpan strings, with the clamp/zero-tier adaptation on the
//! sweep side) so the numbers double as a correctness check.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use pgr::libs::ds::IntSpan;
use pgr::libs::runlist::{depth_at_least, depth_by_level};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::collections::BTreeMap;

const SEQ_LEN: u32 = 50_000;
const MIN_DEPTH: i32 = 3;
const SEED: u64 = 20260813;

// --- baseline: vendored intspan-style Coverage (kept only as the benchmark
//     reference; the production code path is `pgr::libs::runlist`) ---

#[derive(Default, Clone)]
struct Coverage {
    max: i32,
    tiers: BTreeMap<i32, IntSpan>,
}

impl Coverage {
    fn new_len(max: i32, len: i32) -> Self {
        let mut tiers: BTreeMap<i32, IntSpan> = BTreeMap::new();
        tiers.insert(-1, IntSpan::from_pair(1, len));
        tiers.insert(0, IntSpan::from_pair(1, len));
        for i in 1..=max {
            tiers.insert(i, IntSpan::new());
        }
        Self { max, tiers }
    }

    fn bump(&mut self, begin: i32, end: i32) {
        let mut tup = (begin.min(end), begin.max(end));
        if tup.0 == 0 {
            tup.0 = 1;
        }
        let mut intspan = IntSpan::from_pair(tup.0, tup.1);

        if self
            .tiers
            .get(&-1)
            .unwrap()
            .equals(self.tiers.get(&self.max).unwrap())
        {
            return;
        }

        self.tiers.entry(0).and_modify(|e| e.subtract(&intspan));

        for i in 1..=self.max {
            let intersect = self.tiers.get(&i).unwrap().intersect(&intspan);
            self.tiers.entry(i).and_modify(|e| e.merge(&intspan));
            if i + 1 > self.max {
                break;
            }
            intspan = intersect.copy();
        }
    }

    fn max_tier(&self) -> IntSpan {
        self.tiers.get(&self.max).unwrap().copy()
    }

    fn uniq_tiers(&self) -> BTreeMap<i32, IntSpan> {
        let mut tiers = self.tiers.clone();
        for i in 1..self.max {
            let intspan_next = tiers[&(i + 1)].copy();
            tiers.entry(i).and_modify(|e| e.subtract(&intspan_next));
        }
        tiers
    }
}

/// Random 1-based inclusive intervals `[start, end]` on one sequence.
/// `dense` intervals overlap heavily (high coverage depth), `sparse` ones
/// barely overlap (shallow depth).
fn make_intervals(n: usize, dense: bool, rng: &mut StdRng) -> Vec<(u32, u32)> {
    let max_len = if dense { 4_000 } else { 150 };
    (0..n)
        .map(|_| {
            let len = rng.random_range(50..=max_len);
            let start = rng.random_range(1..=SEQ_LEN - len);
            (start, start + len)
        })
        .collect()
}

/// pgr `runlist` consumes half-open `[s, e)` intervals (rg `start-end` is
/// converted with `end + 1`); the vendored `Coverage::bump` takes the same
/// 1-based inclusive `[start, end]` the overlap records carry.
fn to_half_open(intervals: &[(u32, u32)]) -> Vec<(u32, u32)> {
    intervals.iter().map(|&(s, e)| (s, e + 1)).collect()
}

// --- default output path: regions covered at least `MIN_DEPTH` ---

fn run_bump_max_tier(intervals: &[(u32, u32)]) -> IntSpan {
    let mut cov = Coverage::new_len(MIN_DEPTH, SEQ_LEN as i32);
    for &(s, e) in intervals {
        cov.bump(s as i32, e as i32);
    }
    cov.max_tier()
}

fn run_sweep_at_least(intervals: &[(u32, u32)]) -> IntSpan {
    depth_at_least(&to_half_open(intervals), MIN_DEPTH as u32)
}

// --- `--base` / `--mean` path: per-depth tiers, clamped to MIN_DEPTH,
//     with the zero-coverage tier (full length minus covered positions) ---

fn run_bump_levels(intervals: &[(u32, u32)]) -> BTreeMap<i32, IntSpan> {
    let mut cov = Coverage::new_len(MIN_DEPTH, SEQ_LEN as i32);
    for &(s, e) in intervals {
        cov.bump(s as i32, e as i32);
    }
    cov.uniq_tiers()
}

fn run_sweep_levels(intervals: &[(u32, u32)]) -> BTreeMap<i32, IntSpan> {
    let by_level = depth_by_level(&to_half_open(intervals), 1);
    let mut out: BTreeMap<i32, IntSpan> = BTreeMap::new();
    let mut covered = IntSpan::new();
    for (depth, is) in &by_level {
        let depth: i32 = depth.parse().unwrap();
        out.entry(depth.min(MIN_DEPTH)).or_default().merge(is);
        covered.merge(is);
    }
    let mut zero = IntSpan::from_pair(1, SEQ_LEN as i32);
    zero.subtract(&covered);
    out.insert(-1, IntSpan::from_pair(1, SEQ_LEN as i32));
    out.insert(0, zero);
    out
}

fn levels_string(levels: &BTreeMap<i32, IntSpan>) -> String {
    let mut parts = vec![];
    for (depth, is) in levels {
        parts.push(format!("{}:{}", depth, is));
    }
    parts.join("|")
}

fn bench_max_tier(c: &mut Criterion) {
    let mut group = c.benchmark_group("covered_max_tier");
    group.sample_size(10);
    for &n in &[200usize, 1_000] {
        for &dense in &[false, true] {
            let mut rng = StdRng::seed_from_u64(SEED ^ (n as u64) << 32 ^ u64::from(dense));
            let intervals = make_intervals(n, dense, &mut rng);
            let label = format!("{}-{}", if dense { "dense" } else { "sparse" }, n);

            assert_eq!(
                run_bump_max_tier(&intervals).to_string(),
                run_sweep_at_least(&intervals).to_string(),
                "max-tier mismatch for {label}"
            );

            group.bench_with_input(BenchmarkId::new("bump", &label), &intervals, |b, ivs| {
                b.iter(|| black_box(run_bump_max_tier(ivs)))
            });
            group.bench_with_input(BenchmarkId::new("sweep", &label), &intervals, |b, ivs| {
                b.iter(|| black_box(run_sweep_at_least(ivs)))
            });
        }
    }
    group.finish();
}

fn bench_levels(c: &mut Criterion) {
    let mut group = c.benchmark_group("covered_levels");
    group.sample_size(10);
    for &n in &[200usize, 1_000] {
        for &dense in &[false, true] {
            let mut rng = StdRng::seed_from_u64(SEED ^ (n as u64) << 32 ^ u64::from(dense));
            let intervals = make_intervals(n, dense, &mut rng);
            let label = format!("{}-{}", if dense { "dense" } else { "sparse" }, n);

            assert_eq!(
                levels_string(&run_bump_levels(&intervals)),
                levels_string(&run_sweep_levels(&intervals)),
                "levels mismatch for {label}"
            );

            group.bench_with_input(BenchmarkId::new("bump", &label), &intervals, |b, ivs| {
                b.iter(|| black_box(run_bump_levels(ivs)))
            });
            group.bench_with_input(BenchmarkId::new("sweep", &label), &intervals, |b, ivs| {
                b.iter(|| black_box(run_sweep_levels(ivs)))
            });
        }
    }
    group.finish();
}

criterion_group!(benches, bench_max_tier, bench_levels);
criterion_main!(benches);

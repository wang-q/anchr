//! Conservative read-overhang contig-end extension.
//!
//! Walks each contig end base-by-base through the reads' k-mer graph
//! (`RefineTable` right/left extension counts): a base is appended only
//! when it has a strict majority of read support (`>= min_support` reads
//! and `>= 2x` the runner-up), so junctions and repetitive contexts stop
//! the walk instead of joining distant loci. This closes small coverage
//! gaps at contig ends (megahit `local_assemble` goal) without
//! reassembling the reads.

use super::table::{base_code, Kmer, RefineTable};
use anyhow::{ensure, Result};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

/// Options for [`extend_contigs`].
pub struct ExtendOptions {
    /// Seed k-mer length for the extension walk.
    pub k: usize,
    /// Contigs shorter than this are passed through unextended (short
    /// fragments in repeat contexts join copies of the same element).
    pub min_len: usize,
    /// Maximum extension in bases per contig end.
    pub max_extend: usize,
    /// Minimum read support for each appended base.
    pub min_support: u32,
    /// Total (both ends) extension below which the extension is discarded.
    pub min_extend: usize,
}

impl Default for ExtendOptions {
    fn default() -> Self {
        Self {
            k: 31,
            max_extend: 500,
            min_len: 1000,
            min_support: 2,
            min_extend: 0,
        }
    }
}

const BASES: [u8; 4] = *b"ACGT";

/// A base is appended only when its k-mer count is at least this fraction of
/// the contig's own median count. A chimeric extension walks through a seam
/// of near-zero crossing k-mers (the two loci share no read), so the step
/// counts collapse far below the contig's coverage; real gap-closing
/// extensions stay at contig coverage throughout. Same convention as the
/// multik low-coverage seam split (`LOW_RATIO` in `asm/multik/graph.rs`).
const LOW_RATIO: f64 = 0.3;

/// A walk stops after this many consecutive steps whose best k-mer count
/// falls below `LOW_RATIO` of the contig's median. A chimeric seam is a
/// sustained low-coverage run (distant loci share no crossing read), while
/// isolated 1-2 step dips are SNP/noise and must not halt a legitimate
/// walk. Same convention as the multik low-coverage seam split (`MIN_RUN`
/// in `asm/multik/graph.rs`).
const MIN_LOW_RUN: usize = 5;

/// A walk stops after this many consecutive appended bases whose window is
/// already claimed by another input contig (`cross_contig_kmers`). A
/// relocation chimera joins a distant locus whose sequence is already
/// assembled as a different contig, so the walked k-mers beyond the junction
/// sit in that other contig's territory; once the walk stays there for
/// `MIN_FOREIGN_RUN` steps, the extension is rolled back to the last
/// self-owned base. A short isolated shared k-mer (a coincidental motif or a
/// genuine adjacent-contig overlap) does not count, keeping legitimate
/// short overhang closures intact.
const MIN_FOREIGN_RUN: usize = 5;

/// Indexes of the two most-supported extensions of a k-mer.
fn top_two(counts: &[u32; 4]) -> (usize, usize) {
    let mut order = [0usize, 1, 2, 3];
    order.sort_by(|&a, &b| counts[b].cmp(&counts[a]));
    (order[0], order[1])
}

/// Builds the seed k-mer of the first `k` bases of `seq`.
fn seed_kmer(seq: &[u8], k: usize) -> Kmer {
    let mut km = Kmer::new(k);
    for &b in seq.iter().take(k) {
        km.push_right(base_code(b));
    }
    km
}

/// Median canonical k-mer count of `seq`'s own k-mers, the extension
/// walk's coverage baseline. Counted once per contig; the bulk interior
/// dominates, so a chimeric end does not skew it.
fn contig_median_count(table: &RefineTable, seq: &[u8], k: usize) -> u32 {
    let mut counts: Vec<u32> = Vec::with_capacity(seq.len() - k + 1);
    let mut km = Kmer::new(k);
    for &b in seq.iter().take(k) {
        km.push_right(base_code(b));
    }
    counts.push(table.get_count_canonical(&km.canonical()));
    for &b in seq.iter().skip(k) {
        km.push_right(base_code(b));
        counts.push(table.get_count_canonical(&km.canonical()));
    }
    counts.sort_unstable();
    counts[counts.len() / 2]
}

/// Indexes the cross-contig ownership of every canonical k-mer over all
/// input contigs: `sole` maps a k-mer contained in **exactly one** contig to
/// that contig's index, and `multi` holds k-mers that appear in at least two
/// different contigs (conserved/duplicate elements, claimed by everyone). A
/// k-mer is *foreign to contig `i`* when it is in `multi` or owned solely by
/// a contig other than `i` — walking into it means entering sequence another
/// contig already assembled, i.e. a relocation chimera seam.
fn cross_contig_kmers(
    contigs: &[(String, Vec<u8>)],
    k: usize,
) -> (HashMap<Kmer, usize>, HashSet<Kmer>) {
    let mut sole: HashMap<Kmer, usize> = HashMap::new();
    let mut multi: HashSet<Kmer> = HashSet::new();
    let mut local: HashSet<Kmer> = HashSet::new();
    for (i, (_, seq)) in contigs.iter().enumerate() {
        if seq.len() < k {
            continue;
        }
        local.clear();
        let mut km = Kmer::new(k);
        for &b in seq.iter().take(k) {
            km.push_right(base_code(b));
        }
        local.insert(km.canonical());
        for &b in seq.iter().skip(k) {
            km.push_right(base_code(b));
            local.insert(km.canonical());
        }
        for c in local.drain() {
            if multi.contains(&c) {
                continue;
            }
            match sole.get(&c) {
                Some(&o) if o != i => {
                    sole.remove(&c);
                    multi.insert(c);
                }
                _ => {
                    sole.insert(c, i);
                }
            }
        }
    }
    (sole, multi)
}

/// True when canonical `km` lies in the claimed territory of a contig other
/// than `own.contig_id` (see [`cross_contig_kmers`]).
fn is_foreign(km: &Kmer, own: &Ownership<'_>) -> bool {
    own.multi.contains(km) || own.sole.get(km).is_some_and(|&o| o != own.contig_id)
}

/// Cross-contig ownership context for one walk: the contig being extended and
/// the sole/multi ownership indexes built by [`cross_contig_kmers`].
struct Ownership<'a> {
    contig_id: usize,
    sole: &'a HashMap<Kmer, usize>,
    multi: &'a HashSet<Kmer>,
}

/// Walks right from `seed`, appending strictly-majority bases until the
/// support drops, the best count sits below `low_threshold` for
/// `MIN_LOW_RUN` consecutive steps, or the limit is reached.
fn walk_right(
    table: &RefineTable,
    seed: Kmer,
    max_extend: usize,
    min_support: u32,
    low_threshold: u32,
    own: &Ownership<'_>,
) -> Vec<u8> {
    let mut ext = Vec::new();
    let mut cur = seed;
    let mut low_run = 0usize;
    let mut foreign_run = 0usize;
    for step in 0..max_extend {
        let counts = table.fill_right_counts(&cur);
        let (best, second) = top_two(&counts);
        if std::env::var_os("ANCHR_EXTEND_TRACE").is_some() {
            eprintln!(
                "R step{step} ext={} counts={counts:?} best={best}({}) sec={second}({})",
                String::from_utf8_lossy(&ext),
                counts[best],
                counts[second]
            );
        }
        if counts[best] < min_support || counts[best] < 2 * counts[second] {
            break;
        }
        if counts[best] < low_threshold {
            low_run += 1;
            if low_run >= MIN_LOW_RUN {
                // A sustained seam: roll back the run's pushed bases so the
                // extension never includes chimeric crossing sequence.
                ext.truncate(ext.len().saturating_sub(MIN_LOW_RUN - 1));
                break;
            }
        } else {
            low_run = 0;
        }
        ext.push(BASES[best]);
        cur.push_right(best as u8);
        // The new window sits in another contig's claimed territory a
        // sustained number of steps -> chimeric relocation into that locus.
        if is_foreign(&cur.canonical(), own) {
            foreign_run += 1;
            if foreign_run >= MIN_FOREIGN_RUN {
                // Roll back the foreign-run bases so the extension never
                // carries the distant locus sequence.
                ext.truncate(ext.len().saturating_sub(foreign_run));
                break;
            }
        } else {
            foreign_run = 0;
        }
    }
    ext
}

/// Walks left from `seed`, prepending strictly-majority bases; returns the
/// extension in genome order (the bases to attach before the contig).
fn walk_left(
    table: &RefineTable,
    seed: Kmer,
    max_extend: usize,
    min_support: u32,
    low_threshold: u32,
    own: &Ownership<'_>,
) -> Vec<u8> {
    let mut ext = Vec::new();
    let mut cur = seed;
    let mut low_run = 0usize;
    let mut foreign_run = 0usize;
    for step in 0..max_extend {
        let counts = table.fill_left_counts(&cur);
        let (best, second) = top_two(&counts);
        if std::env::var_os("ANCHR_EXTEND_TRACE").is_some() {
            eprintln!(
                "L step{step} ext={} counts={counts:?} best={best}({}) sec={second}({})",
                String::from_utf8_lossy(&ext),
                counts[best],
                counts[second]
            );
        }
        if counts[best] < min_support || counts[best] < 2 * counts[second] {
            break;
        }
        if counts[best] < low_threshold {
            low_run += 1;
            if low_run >= MIN_LOW_RUN {
                ext.truncate(ext.len().saturating_sub(MIN_LOW_RUN - 1));
                break;
            }
        } else {
            low_run = 0;
        }
        ext.push(BASES[best]);
        cur.push_left(best as u8);
        // The new window sits in another contig's claimed territory a
        // sustained number of steps -> chimeric relocation into that locus.
        if is_foreign(&cur.canonical(), own) {
            foreign_run += 1;
            if foreign_run >= MIN_FOREIGN_RUN {
                ext.truncate(ext.len().saturating_sub(foreign_run));
                break;
            }
        } else {
            foreign_run = 0;
        }
    }
    ext.reverse();
    ext
}

/// Extends each contig's ends through read-supported k-mer walks and
/// returns the extended contigs in input order (unchanged when the
/// extension is shorter than `min_extend` or the ends have no support).
pub fn extend_contigs(
    contigs: &[(String, Vec<u8>)],
    reads: Vec<(Vec<u8>, Vec<u8>)>,
    opts: &ExtendOptions,
) -> Result<Vec<(String, Vec<u8>)>> {
    ensure!(
        opts.k >= 2 && opts.k <= pgr::libs::kmer::key::Kmer::MAX_K,
        "k must be in 2..={} (the k-mer key limit)",
        pgr::libs::kmer::key::Kmer::MAX_K
    );
    if contigs.is_empty() {
        return Ok(Vec::new());
    }
    let table = RefineTable::build_supermer(reads, opts.k, None)?;
    let (sole, multi) = cross_contig_kmers(contigs, opts.k);
    let extended: Vec<(String, Vec<u8>)> = contigs
        .par_iter()
        .enumerate()
        .map(|(i, (name, seq))| {
            if seq.len() < opts.min_len || seq.len() < opts.k {
                return (name.clone(), seq.clone());
            }
            let median = contig_median_count(&table, seq, opts.k);
            let low_threshold = (LOW_RATIO * median as f64) as u32;
            let own = Ownership {
                contig_id: i,
                sole: &sole,
                multi: &multi,
            };
            let right = walk_right(
                &table,
                seed_kmer(&seq[seq.len() - opts.k..], opts.k),
                opts.max_extend,
                opts.min_support,
                low_threshold,
                &own,
            );
            let left = walk_left(
                &table,
                seed_kmer(&seq[..opts.k], opts.k),
                opts.max_extend,
                opts.min_support,
                low_threshold,
                &own,
            );
            if left.len() + right.len() < opts.min_extend {
                return (name.clone(), seq.clone());
            }
            let mut out = Vec::with_capacity(left.len() + seq.len() + right.len());
            out.extend_from_slice(&left);
            out.extend_from_slice(seq);
            out.extend_from_slice(&right);
            (name.clone(), out)
        })
        .collect();
    Ok(extended)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mk_reads(seqs: &[&[u8]]) -> Vec<(Vec<u8>, Vec<u8>)> {
        seqs.iter().map(|s| (s.to_vec(), Vec::new())).collect()
    }

    /// Deterministic random ACGT sequence (no repeats, so the extension
    /// k-mers are unique in the reads).
    fn rand_seq(rng: &mut u64, n: usize) -> Vec<u8> {
        let mut s = Vec::new();
        for _ in 0..n {
            *rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            s.push(b"ACGT"[(*rng >> 33) as usize % 4]);
        }
        s
    }

    #[test]
    fn extends_through_read_overhang() {
        // Contig = unique prefix + seed; reads = the full sequence prefix +
        // seed + 30-bp extension, repeated 10x. The right end must walk the
        // 30-bp extension and stop at the read boundary.
        let mut rng = 42u64;
        let prefix = rand_seq(&mut rng, 40);
        let seed = rand_seq(&mut rng, 11);
        let ext = rand_seq(&mut rng, 30);
        let full: Vec<u8> = [&prefix[..], &seed[..], &ext[..]].concat();
        let reads = mk_reads(&[&full[..]; 10]);
        let mut contig_seq = prefix.clone();
        contig_seq.extend_from_slice(&seed);
        let contig = vec![("c1".to_string(), contig_seq)];
        let opts = ExtendOptions {
            k: 11,
            max_extend: 60,
            min_len: 0,
            min_support: 2,
            min_extend: 0,
        };
        let out = extend_contigs(&contig, reads, &opts).unwrap();
        assert_eq!(out[0].1, full);
    }

    #[test]
    fn stops_at_ambiguous_junction() {
        // The seed's right extension splits 5/5 between two bases: no
        // strict majority, so the contig stays unchanged.
        let mut rng = 7u64;
        let prefix = rand_seq(&mut rng, 40);
        let seed = rand_seq(&mut rng, 11);
        let ext_a = rand_seq(&mut rng, 30);
        let ext_b = rand_seq(&mut rng, 30);
        let mut reads = Vec::new();
        for _ in 0..5 {
            let mut a = prefix.clone();
            a.extend_from_slice(&seed);
            a.extend_from_slice(&ext_a);
            let mut b = prefix.clone();
            b.extend_from_slice(&seed);
            b.extend_from_slice(&ext_b);
            reads.push(a);
            reads.push(b);
        }
        let mut contig_seq = prefix.clone();
        contig_seq.extend_from_slice(&seed);
        let contig = vec![("c1".to_string(), contig_seq.clone())];
        let opts = ExtendOptions {
            k: 11,
            max_extend: 60,
            min_len: 0,
            min_support: 2,
            min_extend: 0,
        };
        let read_refs: Vec<&[u8]> = reads.iter().map(|v| v.as_slice()).collect();
        let out = extend_contigs(&contig, mk_reads(&read_refs), &opts).unwrap();
        assert_eq!(out[0].1, contig_seq);
    }

    #[test]
    fn stops_at_low_coverage_seam() {
        // High-coverage extension (12 reads) followed by a low-coverage
        // bridge (2 reads). The walk must extend through the real overhang
        // at contig coverage, then stop at the seam: bridge k-mers satisfy
        // min_support=2 but sit below LOW_RATIO x median, so they must not
        // be appended.
        let mut rng = 123u64;
        let prefix = rand_seq(&mut rng, 40);
        let seed = rand_seq(&mut rng, 11);
        let high_ext = rand_seq(&mut rng, 30);
        let bridge = rand_seq(&mut rng, 20);
        let mut reads = Vec::new();
        for _ in 0..10 {
            let mut r = prefix.clone();
            r.extend_from_slice(&seed);
            r.extend_from_slice(&high_ext);
            reads.push(r);
        }
        for _ in 0..2 {
            let mut r = prefix.clone();
            r.extend_from_slice(&seed);
            r.extend_from_slice(&high_ext);
            r.extend_from_slice(&bridge);
            reads.push(r);
        }
        let mut contig_seq = prefix.clone();
        contig_seq.extend_from_slice(&seed);
        let contig = vec![("c1".to_string(), contig_seq)];
        let opts = ExtendOptions {
            k: 11,
            max_extend: 80,
            min_len: 0,
            min_support: 2,
            min_extend: 0,
        };
        let read_refs: Vec<&[u8]> = reads.iter().map(|v| v.as_slice()).collect();
        let out = extend_contigs(&contig, mk_reads(&read_refs), &opts).unwrap();
        // The right end walks the real 30-bp overhang at contig coverage and
        // stops at the seam, so the extension ends exactly at `high_ext`.
        // If the low-coverage check failed, the 2-read bridge would follow.
        assert!(
            out[0].1.ends_with(&high_ext),
            "seam was crossed into bridge"
        );
        // The bridge must not be appended after the real overhang.
        let tail: Vec<u8> = high_ext.iter().chain(bridge.iter()).copied().collect();
        assert!(
            !out[0].1.ends_with(&tail),
            "low-coverage bridge was appended"
        );
    }

    #[test]
    fn requires_min_support() {
        // A single read overhang is below min_support: no extension.
        let mut rng = 99u64;
        let prefix = rand_seq(&mut rng, 40);
        let seed = rand_seq(&mut rng, 11);
        let ext = rand_seq(&mut rng, 30);
        let full: Vec<u8> = [&prefix[..], &seed[..], &ext[..]].concat();
        let mut contig_seq = prefix.clone();
        contig_seq.extend_from_slice(&seed);
        let contig = vec![("c1".to_string(), contig_seq.clone())];
        let opts = ExtendOptions {
            k: 11,
            max_extend: 60,
            min_len: 0,
            min_support: 2,
            min_extend: 0,
        };
        let out = extend_contigs(&contig, mk_reads(&[&full]), &opts).unwrap();
        assert_eq!(out[0].1, contig_seq);
    }

    #[test]
    fn stops_when_crossing_another_contig() {
        // Two contigs: A ends in unique sequence, B owns a distant region.
        // Reads bridge A's real overhang straight into B's unique sequence
        // (solely owned by B), so the walk would join the two loci. The
        // cross-contig ownership check must stop the walk before the
        // relocation completes, never fully reproducing B's region on A's end.
        let mut rng = 5u64;
        let prefix = rand_seq(&mut rng, 30);
        let seed = rand_seq(&mut rng, 11);
        let remote = rand_seq(&mut rng, 60); // B's unique region the walk enters
        let b_filler = rand_seq(&mut rng, 40);
        let mut a_seq = prefix.clone();
        a_seq.extend_from_slice(&seed);
        let mut b_seq = b_filler.clone();
        b_seq.extend_from_slice(&remote);
        let contigs = vec![
            ("A".to_string(), a_seq.clone()),
            ("B".to_string(), b_seq.clone()),
        ];
        let mut reads: Vec<Vec<u8>> = Vec::new();
        for _ in 0..10 {
            let mut r = prefix.clone();
            r.extend_from_slice(&seed);
            r.extend_from_slice(&remote);
            reads.push(r);
        }
        // Keep B's own reads present so its region is genuinely covered.
        for _ in 0..4 {
            reads.push(b_seq.clone());
        }
        let opts = ExtendOptions {
            k: 11,
            max_extend: 100,
            min_len: 0,
            min_support: 2,
            min_extend: 0,
        };
        let read_refs: Vec<&[u8]> = reads.iter().map(|v| v.as_slice()).collect();
        let out = extend_contigs(&contigs, mk_reads(&read_refs), &opts).unwrap();
        let full: Vec<u8> = [&a_seq[..], &remote[..]].concat();
        // A must not be extended through the full B-owned region.
        assert_ne!(out[0].1, full, "crossed into another contig's territory");
        // The bulk of the remote region must not be appended.
        assert!(
            !out[0].1.ends_with(&remote[20..]),
            "remote-contig sequence was appended"
        );
    }
}

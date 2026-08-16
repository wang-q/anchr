//! Conservative read-overhang contig-end extension.
//!
//! Walks each contig end base-by-base through the reads' k-mer graph
//! (`RefineTable` right/left extension counts): a base is appended only
//! when it has a strict majority of read support (`>= min_support` reads
//! and `>= 2x` the runner-up), so junctions and repetitive contexts stop
//! the walk instead of joining distant loci. This closes small coverage
//! gaps at contig ends (megahit `local_assemble` goal) without
//! reassembling the reads.

use super::refine::{base_code, Kmer, RefineTable};
use anyhow::{ensure, Result};
use rayon::prelude::*;

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

/// Walks right from `seed`, appending strictly-majority bases until the
/// support drops or the limit is reached.
fn walk_right(table: &RefineTable, seed: Kmer, max_extend: usize, min_support: u32) -> Vec<u8> {
    let mut ext = Vec::new();
    let mut cur = seed;
    for _ in 0..max_extend {
        let counts = table.fill_right_counts(&cur);
        let (best, second) = top_two(&counts);
        if counts[best] < min_support || counts[best] < 2 * counts[second] {
            break;
        }
        ext.push(BASES[best]);
        cur.push_right(best as u8);
    }
    ext
}

/// Walks left from `seed`, prepending strictly-majority bases; returns the
/// extension in genome order (the bases to attach before the contig).
fn walk_left(table: &RefineTable, seed: Kmer, max_extend: usize, min_support: u32) -> Vec<u8> {
    let mut ext = Vec::new();
    let mut cur = seed;
    for _ in 0..max_extend {
        let counts = table.fill_left_counts(&cur);
        let (best, second) = top_two(&counts);
        if counts[best] < min_support || counts[best] < 2 * counts[second] {
            break;
        }
        ext.push(BASES[best]);
        cur.push_left(best as u8);
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
    let extended: Vec<(String, Vec<u8>)> = contigs
        .par_iter()
        .map(|(name, seq)| {
            if seq.len() < opts.min_len || seq.len() < opts.k {
                return (name.clone(), seq.clone());
            }
            let right = walk_right(
                &table,
                seed_kmer(&seq[seq.len() - opts.k..], opts.k),
                opts.max_extend,
                opts.min_support,
            );
            let left = walk_left(
                &table,
                seed_kmer(&seq[..opts.k], opts.k),
                opts.max_extend,
                opts.min_support,
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
}

//! Canonical k-mer count table shared by every assembly hot path.
//!
//! `RefineTable` wraps pgr's packed `KmerTable` with prefix-bucket lookups
//! and solid-entry snapshots; `Kmer` is the rolling 2-bit window the walks
//! use. Build paths: quality-gated direct counting, FastK-style super-mer
//! counting, and a streaming variant that never holds all reads.

use pgr::libs::fmt::seq::{SeqReader, SeqRecord};
use pgr::libs::fq::qual::{from_phred, to_phred};
use pgr::libs::kmer::key;
use rayon::prelude::*;
use std::sync::OnceLock;

/// Canonical k-mer count table (FastK 2-bit packing, forward vs
/// reverse-complement canonical form).
///
/// Storage is pgr's `KmerTable`: packed canonical keys (ascending byte
/// order, which equals `Kmer::cmp_bases` order) + parallel `u32` counts,
/// built by quality-gated emission + `pgr::kmer::count::count_keys`
/// (parallel radix sort, no HashMap). The old per-chunk `HashMap` build
/// collected every chunk map before merging (5.3 GB peak on G37 full);
/// the packed sort path keeps intermediate buffers to a single global key
/// list, so both memory and counting time drop an order of magnitude.
#[derive(Debug, Clone, Default)]
pub struct RefineTable {
    table: pgr::libs::kmer::KmerTable,
    /// 2-byte-prefix bucket offsets into `table.keys` (65537 entries, or
    /// 257 for k < 8): `get_count` binary-searches within the bucket, so a
    /// long-k query costs ~4-5 compares instead of log2(n) ≈ 20.
    prefix_index: OnceLock<Vec<u32>>,
    /// (key, count) snapshot built once on first request (the assemble
    /// passes scan it multiple times; the packed table stays the source).
    sorted: OnceLock<Vec<(Kmer, u32)>>,
}

impl RefineTable {
    /// Builds the table from reads (bases, phred qualities) with `minprob`
    /// quality filtering, mirroring `KmerTableSetU.addKmersToTable`.
    pub fn build(reads: &[(Vec<u8>, Vec<u8>)], k: usize, min_prob: f32) -> Self {
        Self::build_threaded(reads, k, min_prob, 0)
    }

    /// `build` with an explicit worker count (kept for call-site symmetry).
    ///
    /// `threads` is advisory: the parallel per-chunk emission runs on the
    /// ambient rayon pool (the caller wraps the assemble call in a single
    /// `--parallel` pool), so no private pool is ever spawned here.
    pub(crate) fn build_threaded(
        reads: &[(Vec<u8>, Vec<u8>)],
        k: usize,
        min_prob: f32,
        _threads: usize,
    ) -> Self {
        let (prob_correct, prob_correct_inv) = prob_tables();
        if reads.is_empty() {
            // Keep `k` so the empty table stays queryable: a defaulted
            // table (k=0) would panic on `get_count` (0/0 in prefix_index)
            // and on the `chunks_exact(0)` of the sorted/solid snapshots.
            return Self {
                table: pgr::libs::kmer::KmerTable {
                    k,
                    keys: Vec::new(),
                    counts: Vec::new(),
                },
                prefix_index: OnceLock::new(),
                sorted: OnceLock::new(),
            };
        }
        // Parallel per-chunk emission of quality-gated packed canonical
        // keys; each chunk is deduplicated in place by the shared pgr
        // sort+group core (`count_keys`), then the sorted per-chunk tables
        // are merged pairwise (rayon tree reduce). Only O(threads) chunk
        // buffers coexist and the duplicate-containing global key list
        // (~1 GB at G37 full) never materializes; the merged output is the
        // sorted deduplicated table regardless of merge order, so the
        // result is deterministic.
        let chunk_size = 16384usize;
        let key_bytes = k.div_ceil(4);
        let build = |reads: &[(Vec<u8>, Vec<u8>)]| -> pgr::libs::kmer::KmerTable {
            reads
                .par_chunks(chunk_size)
                .map(|chunk| {
                    let cap = chunk
                        .iter()
                        .map(|(b, _)| b.len().saturating_sub(k - 1) * key_bytes)
                        .sum();
                    let mut raw = Vec::with_capacity(cap);
                    for (bases, quals) in chunk {
                        count_read_kmers_packed(
                            &mut raw,
                            bases,
                            quals,
                            k,
                            min_prob,
                            &prob_correct,
                            &prob_correct_inv,
                        );
                    }
                    pgr::libs::kmer::count::count_keys(raw, k)
                })
                .reduce(
                    || pgr::libs::kmer::KmerTable {
                        k,
                        keys: Vec::new(),
                        counts: Vec::new(),
                    },
                    merge_tables,
                )
        };
        // Ambient pool: the caller wraps the assemble call in a single rayon
        // pool of `--parallel` threads; creating another pool here would
        // oversubscribe the machine.
        let table = build(reads);
        Self {
            table,
            prefix_index: OnceLock::new(),
            sorted: OnceLock::new(),
        }
    }

    /// FastK-style super-mer two-stage counting (pgr `kmer::supermer`).
    /// Counts every N-free k-mer without quality gating (byte-identical to
    /// the direct path on FASTA / no-quality input).
    pub(crate) fn build_supermer(
        reads: Vec<(Vec<u8>, Vec<u8>)>,
        k: usize,
        m: Option<usize>,
    ) -> anyhow::Result<Self> {
        let t0 = std::time::Instant::now();
        // Borrow the sequence buffers (slices API): the super-mer path does
        // not use qualities, and no `Vec<Vec<u8>>` is materialized.
        let seqs: Vec<&[u8]> = reads.iter().map(|(s, _)| s.as_slice()).collect();
        let table = match m {
            Some(m) => pgr::libs::kmer::supermer::build_table_slices_with_m(&seqs, k, m)?,
            None => pgr::libs::kmer::supermer::build_table_slices(&seqs, k)?,
        };
        if std::env::var_os("ANCHR_SM_TIMING").is_some() {
            eprintln!(
                "supermer count: {:.3}s (move+pack+sort+expand)",
                t0.elapsed().as_secs_f64()
            );
        }
        Ok(Self {
            table,
            prefix_index: OnceLock::new(),
            sorted: OnceLock::new(),
        })
    }

    /// Borrowing variant of [`RefineTable::build_supermer`] with the
    /// adaptive minimizer length: callers that keep their sequence buffers
    /// alive (multik reuses one reads buffer across rounds) count without
    /// moving or copying anything.
    pub(crate) fn build_supermer_slices(seqs: &[&[u8]], k: usize) -> anyhow::Result<Self> {
        let table = pgr::libs::kmer::supermer::build_table_slices(seqs, k)?;
        Ok(Self {
            table,
            prefix_index: OnceLock::new(),
            sorted: OnceLock::new(),
        })
    }

    /// Direct canonical counting over borrowed slices: byte-identical to
    /// [`RefineTable::build_supermer_slices`] but skips the super-mer
    /// stage. On unique sequence (unitigs, coverage ~1) super-mers never
    /// collapse — every window becomes an 80-byte stage-1 record that
    /// stage 2 must expand again — while the direct path sorts one
    /// `ceil(k/4)`-byte key per window and groups. High-coverage inputs
    /// (reads) still prefer the super-mer path.
    pub(crate) fn build_direct_slices(seqs: &[&[u8]], k: usize) -> anyhow::Result<Self> {
        anyhow::ensure!(
            k > 0 && k <= pgr::libs::kmer::key::Kmer::MAX_K,
            "k must be in 1..={}, got {k}",
            pgr::libs::kmer::key::Kmer::MAX_K
        );
        const CHUNK: usize = 512;
        let key_bytes = k.div_ceil(4);
        let per_chunk: Vec<Vec<u8>> = seqs
            .par_chunks(CHUNK)
            .map(|chunk| {
                let cap = chunk
                    .iter()
                    .map(|s| s.len().saturating_sub(k - 1) * key_bytes)
                    .sum();
                let mut raw = Vec::with_capacity(cap);
                for seq in chunk {
                    pgr::libs::kmer::canonical_keys(seq, k, |_, km| {
                        raw.extend_from_slice(km.to_bytes());
                    });
                }
                raw
            })
            .collect();
        let n: usize = per_chunk.iter().map(Vec::len).sum();
        let mut keys: Vec<u8> = Vec::with_capacity(n);
        for mut v in per_chunk {
            keys.append(&mut v);
        }
        let table = pgr::libs::kmer::count::count_keys(keys, k);
        Ok(Self {
            table,
            prefix_index: OnceLock::new(),
            sorted: OnceLock::new(),
        })
    }

    /// Streaming direct counter: reads `infiles` record-by-record, fans out
    /// to `threads` workers that emit packed keys in bounded chunks and
    /// merge per-worker count tables. Unlike `build_threaded` this never
    /// holds all reads (the ~0.65 GB `(seq, phred)` copy on G37 full), so
    /// peak memory stays bounded by the chunk/worker buffers. Returns the
    /// table and the number of records consumed.
    pub(crate) fn build_streamed(
        infiles: &[String],
        k: usize,
        min_prob: f32,
        threads: usize,
    ) -> anyhow::Result<(Self, u64)> {
        let threads = threads.max(1);
        let (prob_correct, prob_correct_inv) = prob_tables();
        let chunk: usize = std::env::var("ANCHR_STREAM_CHUNK")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(32768);
        let cap: usize = std::env::var("ANCHR_STREAM_CAP")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(2);
        std::thread::scope(|s| -> anyhow::Result<(Self, u64)> {
            let mut senders = Vec::with_capacity(threads);
            let mut receivers = Vec::with_capacity(threads);
            for _ in 0..threads {
                let (tx, rx) = std::sync::mpsc::sync_channel::<Vec<(Vec<u8>, Vec<u8>)>>(cap);
                senders.push(tx);
                receivers.push(rx);
            }
            let mut handles = Vec::with_capacity(threads);
            for rx in receivers {
                let (pc, pci) = (prob_correct.clone(), prob_correct_inv.clone());
                handles.push(s.spawn(move || {
                    let mut raw: Vec<u8> = Vec::new();
                    let mut table = pgr::libs::kmer::KmerTable {
                        k,
                        keys: Vec::new(),
                        counts: Vec::new(),
                    };
                    while let Ok(recs) = rx.recv() {
                        for (seq, qual) in recs {
                            count_read_kmers_packed(&mut raw, &seq, &qual, k, min_prob, &pc, &pci);
                        }
                        if !raw.is_empty() {
                            let t = count_keys_seq(std::mem::take(&mut raw), k);
                            table = merge_tables(table, t);
                        }
                    }
                    table
                }));
            }
            let mut reads_in = 0u64;
            let mut idx = 0usize;
            let mut buf: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(chunk);
            for infile in infiles {
                let mut reader = SeqReader::new(infile)?;
                let mut rec = SeqRecord::new();
                while reader.read_record(&mut rec)? {
                    canonicalize_quality(&mut rec);
                    buf.push((
                        rec.sequence().to_vec(),
                        to_phred(rec.sequence(), rec.quality_scores()),
                    ));
                    reads_in += 1;
                    if buf.len() >= chunk {
                        senders[idx % threads]
                            .send(std::mem::take(&mut buf))
                            .map_err(|_| anyhow::anyhow!("streaming count channel closed"))?;
                        idx += 1;
                    }
                }
            }
            if !buf.is_empty() {
                senders[idx % threads]
                    .send(buf)
                    .map_err(|_| anyhow::anyhow!("streaming count channel closed"))?;
            }
            drop(senders);
            let mut table = pgr::libs::kmer::KmerTable {
                k,
                keys: Vec::new(),
                counts: Vec::new(),
            };
            for h in handles {
                table = merge_tables(table, h.join().expect("streaming count worker panicked"));
            }
            Ok((
                Self {
                    table,
                    prefix_index: OnceLock::new(),
                    sorted: OnceLock::new(),
                },
                reads_in,
            ))
        })
    }

    /// Sorted-start offsets per 1-, 2- or 3-byte key prefix (bucket `p`
    /// spans `[offs[p], offs[p+1])`), built lazily in one O(n) scan. Big
    /// tables (reads) use 3-byte prefixes: 16M buckets hold well under one
    /// row each on average, so a lookup is one index read plus at most a
    /// couple of key compares instead of a ~8-probe binary search over a
    /// 64K-bucket range.
    fn prefix_index(&self) -> &[u32] {
        self.prefix_index.get_or_init(|| {
            let kb = self.table.key_bytes();
            let keys = &self.table.keys;
            let n = keys.len() / kb;
            let width = if kb == 1 {
                1
            } else if kb >= 3 && n > (1 << 20) {
                3
            } else {
                2
            };
            let entries = 256usize.pow(width as u32) + 1;
            let mut offs = vec![0u32; entries];
            let mut i = 0usize;
            for (p, slot) in offs.iter_mut().take(entries - 1).enumerate() {
                let lo = i;
                while i < n {
                    let s = i * kb;
                    let pref = match width {
                        1 => keys[s] as usize,
                        2 => ((keys[s] as usize) << 8) | keys[s + 1] as usize,
                        _ => {
                            ((keys[s] as usize) << 16)
                                | ((keys[s + 1] as usize) << 8)
                                | keys[s + 2] as usize
                        }
                    };
                    if pref > p {
                        break;
                    }
                    i += 1;
                }
                *slot = lo as u32;
            }
            offs[entries - 1] = n as u32;
            offs
        })
    }

    /// Row of canonical packed key `q` (already `kb` bytes), None when
    /// absent: prefix bucket + at most a couple of in-bucket compares.
    fn locate(&self, q: &[u8]) -> Option<usize> {
        let kb = self.table.key_bytes();
        let offs = self.prefix_index();
        let p = match offs.len() {
            257 => q[0] as usize,
            65537 => ((q[0] as usize) << 8) | q[1] as usize,
            _ => ((q[0] as usize) << 16) | ((q[1] as usize) << 8) | q[2] as usize,
        };
        let keys = &self.table.keys;
        let mut lo = offs[p] as usize;
        let mut hi = offs[p + 1] as usize;
        while lo < hi {
            let mid = (lo + hi) >> 1;
            let mid_b = &keys[mid * kb..(mid + 1) * kb];
            match mid_b.cmp(q) {
                std::cmp::Ordering::Less => lo = mid + 1,
                std::cmp::Ordering::Greater => hi = mid,
                std::cmp::Ordering::Equal => return Some(mid),
            }
        }
        None
    }

    /// Count of the canonical form of `kmer` (0 when absent).
    pub(crate) fn get_count(&self, kmer: &Kmer) -> u32 {
        let kb = self.table.key_bytes();
        let mut qbuf = [0u8; key::Kmer::MAX_K / 4];
        // Whole-byte k (k % 4 == 0): packed rc via a byte table (kb lookups
        // + half-byte canonical compare) instead of `key::Kmer::canonical`,
        // which recomputes the rc base-by-base (O(k)); at k=100 the
        // traversal's get_count canonicalization was ~30% of runtime.
        if self.table.k.is_multiple_of(4) {
            let fw = kmer.0.to_bytes();
            let mut rc = [0u8; key::Kmer::MAX_K / 4];
            for i in 0..kb {
                rc[kb - 1 - i] = REVCOMP_BYTE[fw[i] as usize];
            }
            let half = kb.div_ceil(2);
            let canon: &[u8] = if fw[..half] <= rc[..half] {
                fw
            } else {
                &rc[..kb]
            };
            qbuf[..kb].copy_from_slice(canon);
        } else {
            let canon = kmer.canonical();
            qbuf[..kb].copy_from_slice(canon.0.to_bytes());
        }
        let q = &qbuf[..kb];
        // Prefix bucket + in-bucket binary search (see `prefix_index`).
        self.locate(q).map(|r| self.table.counts[r]).unwrap_or(0)
    }

    /// Count of an already-canonicalized k-mer (0 when absent): skips the
    /// forward/rc canonicalization that `get_count` performs, for callers
    /// that roll a window and track both strands in lockstep.
    pub(crate) fn get_count_canonical(&self, kmer: &Kmer) -> u32 {
        let kb = self.table.key_bytes();
        let q = &kmer.0.to_bytes()[..kb];
        self.locate(q).map(|r| self.table.counts[r]).unwrap_or(0)
    }

    /// Row and count of an oriented k-mer's canonical form (None when
    /// absent from the table): like [`Self::get_count`] but also reports
    /// the packed-table row, so classification can link neighbouring
    /// vertices without a second lookup.
    pub(crate) fn find_row(&self, kmer: &Kmer) -> Option<(usize, u32)> {
        let kb = self.table.key_bytes();
        let mut qbuf = [0u8; key::Kmer::MAX_K / 4];
        if self.table.k.is_multiple_of(4) {
            let fw = kmer.0.to_bytes();
            let mut rc = [0u8; key::Kmer::MAX_K / 4];
            for i in 0..kb {
                rc[kb - 1 - i] = REVCOMP_BYTE[fw[i] as usize];
            }
            let half = kb.div_ceil(2);
            let canon: &[u8] = if fw[..half] <= rc[..half] {
                fw
            } else {
                &rc[..kb]
            };
            qbuf[..kb].copy_from_slice(canon);
        } else {
            let canon = kmer.canonical();
            qbuf[..kb].copy_from_slice(canon.0.to_bytes());
        }
        let q = &qbuf[..kb];
        self.locate(q).map(|r| (r, self.table.counts[r]))
    }

    /// Per-row entry rank for rows with `count >= threshold` (`u32::MAX`
    /// elsewhere), in packed-table order: `solid_entries` preserves that
    /// order, so the rank maps a found row onto the parallel entries index.
    pub(crate) fn solid_row_ranks(&self, threshold: u32) -> Vec<u32> {
        let mut ranks = vec![u32::MAX; self.table.counts.len()];
        let mut rank = 0u32;
        for (row, &c) in self.table.counts.iter().enumerate() {
            if c >= threshold {
                ranks[row] = rank;
                rank += 1;
            }
        }
        ranks
    }

    /// Deterministic (canonical k-mer) sorted snapshot of (key, count),
    /// computed once and cached for the multi-pass assemble scans.
    pub(crate) fn sorted_entries(&self) -> &[(Kmer, u32)] {
        self.sorted.get_or_init(|| {
            let kb = self.table.key_bytes();
            self.table
                .keys
                .chunks_exact(kb)
                .enumerate()
                .map(|(i, b)| {
                    let km = Kmer(pgr::libs::kmer::key::Kmer::from_bytes(self.table.k, b));
                    (km, self.table.counts[i])
                })
                .collect()
        })
    }

    /// Sorted (canonical k-mer, count) pairs with `count >= threshold`,
    /// built by scanning the packed table without materializing the
    /// below-threshold keys (the walk/classification passes only visit
    /// solid vertices).
    pub(crate) fn solid_entries(&self, threshold: u32) -> Vec<(Kmer, u32)> {
        let kb = self.table.key_bytes();
        self.table
            .keys
            .chunks_exact(kb)
            .zip(self.table.counts.iter())
            .filter(|(_, &c)| c >= threshold)
            .map(|(b, &c)| {
                (
                    Kmer(pgr::libs::kmer::key::Kmer::from_bytes(self.table.k, b)),
                    c,
                )
            })
            .collect()
    }

    /// Counts of the four right-extensions of `kmer`.
    pub(crate) fn fill_right_counts(&self, kmer: &Kmer) -> [u32; 4] {
        let mut out = [0u32; 4];
        for (i, c) in out.iter_mut().enumerate() {
            let mut x = *kmer;
            x.push_right(i as u8);
            *c = self.get_count(&x);
        }
        out
    }

    /// Counts of the four left-extensions of `kmer`.
    pub(crate) fn fill_left_counts(&self, kmer: &Kmer) -> [u32; 4] {
        let mut out = [0u32; 4];
        for (i, c) in out.iter_mut().enumerate() {
            let mut x = *kmer;
            x.push_left(i as u8);
            *c = self.get_count(&x);
        }
        out
    }
}

/// Sequential `count_keys`: sorts the packed keys with the non-parallel MSD
/// radix path and groups identical keys into counts. Used by the streamed
/// worker threads so they never spawn rayon work (single-pool pipeline).
fn count_keys_seq(mut keys: Vec<u8>, k: usize) -> pgr::libs::kmer::KmerTable {
    let key_bytes = k.div_ceil(4);
    if keys.is_empty() {
        return pgr::libs::kmer::KmerTable {
            k,
            keys,
            counts: Vec::new(),
        };
    }
    let n_keys = keys.len() / key_bytes;
    pgr::libs::ds::radix_sort::radix_sort_bytes(&mut keys, key_bytes, &mut vec![(); n_keys]);
    let mut counts: Vec<u32> = Vec::with_capacity(n_keys);
    let mut i = 0usize;
    let mut w = 0usize;
    while i < n_keys {
        let mut j = i + 1;
        while j < n_keys
            && keys[j * key_bytes..(j + 1) * key_bytes] == keys[i * key_bytes..(i + 1) * key_bytes]
        {
            j += 1;
        }
        if w != i {
            keys.copy_within(i * key_bytes..(i + 1) * key_bytes, w * key_bytes);
        }
        counts.push((j - i).min(u32::MAX as usize) as u32);
        w += 1;
        i = j;
    }
    keys.truncate(w * key_bytes);
    pgr::libs::kmer::KmerTable { k, keys, counts }
}

/// Merge two sorted, deduplicated k-mer tables into one (combining equal
/// keys' counts). Associative, so rayon tree reduction order doesn't affect
/// the result.
fn merge_tables(
    a: pgr::libs::kmer::KmerTable,
    b: pgr::libs::kmer::KmerTable,
) -> pgr::libs::kmer::KmerTable {
    debug_assert_eq!(a.k, b.k);
    let kb = a.key_bytes();
    let (na, nb) = (a.keys.len() / kb, b.keys.len() / kb);
    let mut keys = Vec::with_capacity(a.keys.len() + b.keys.len());
    let mut counts = Vec::with_capacity(na + nb);
    let (mut i, mut j) = (0usize, 0usize);
    while i < na && j < nb {
        let ka = &a.keys[i * kb..(i + 1) * kb];
        let kbk = &b.keys[j * kb..(j + 1) * kb];
        match ka.cmp(kbk) {
            std::cmp::Ordering::Less => {
                keys.extend_from_slice(ka);
                counts.push(a.counts[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                keys.extend_from_slice(kbk);
                counts.push(b.counts[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                keys.extend_from_slice(ka);
                counts.push(a.counts[i] + b.counts[j]);
                i += 1;
                j += 1;
            }
        }
    }
    while i < na {
        keys.extend_from_slice(&a.keys[i * kb..(i + 1) * kb]);
        counts.push(a.counts[i]);
        i += 1;
    }
    while j < nb {
        keys.extend_from_slice(&b.keys[j * kb..(j + 1) * kb]);
        counts.push(b.counts[j]);
        j += 1;
    }
    pgr::libs::kmer::KmerTable {
        k: a.k,
        keys,
        counts,
    }
}

/// FastK byte k-mer key (2 bits/base, bytes 5'->3') used by the assembly
/// hot paths; a thin wrapper over `key::Kmer` kept for call-site brevity.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub(crate) struct Kmer(key::Kmer);

/// FNV-1a hasher for `HashSet<Kmer>` claim/visited sets (the derived
/// `Hash` uses `RandomState`/SipHash13, which was ~7% of runtime at k=100).
#[derive(Clone, Copy)]
pub(crate) struct KmerFnvHasher(u64);

impl Default for KmerFnvHasher {
    fn default() -> Self {
        Self(0xcbf29ce484222325)
    }
}

impl std::hash::Hasher for KmerFnvHasher {
    fn write(&mut self, bytes: &[u8]) {
        let mut h = self.0;
        let mut chunks = bytes.chunks_exact(8);
        for c in &mut chunks {
            h ^= u64::from_ne_bytes(c.try_into().unwrap());
            h = h.wrapping_mul(0x100000001b3);
        }
        let rem = chunks.remainder();
        if !rem.is_empty() {
            let mut buf = [0u8; 8];
            buf[..rem.len()].copy_from_slice(rem);
            h ^= u64::from_ne_bytes(buf);
            h = h.wrapping_mul(0x100000001b3);
        }
        self.0 = h;
    }

    fn finish(&self) -> u64 {
        self.0
    }
}

/// Reverse-complement of a packed FastK byte (4 × 2-bit codes high->low:
/// `(3-c3)(3-c2)(3-c1)(3-c0)`). Exact for whole bytes; when `k % 4 != 0`
/// the boundary byte's unused bits break the byte symmetry, so the table
/// fast path is restricted to `k % 4 == 0` (e.g. k=100).
static REVCOMP_BYTE: [u8; 256] = {
    let mut t = [0u8; 256];
    let mut b = 0usize;
    while b < 256 {
        let c0 = (b >> 6) & 3;
        let c1 = (b >> 4) & 3;
        let c2 = (b >> 2) & 3;
        let c3 = b & 3;
        t[b] = (((3 - c3) << 6) | ((3 - c2) << 4) | ((3 - c1) << 2) | (3 - c0)) as u8;
        b += 1;
    }
    t
};

impl Kmer {
    /// Empty window of length `k` (`k <= key::Kmer::MAX_K`, caller-validated).
    pub(crate) fn new(k: usize) -> Self {
        Self(key::Kmer::new(k).expect("assembly k in 1..=MAX_K"))
    }

    /// Reset the window to all-zero bases.
    pub(crate) fn reset(&mut self) {
        *self = Self::new(self.0.k());
    }

    /// Base `i` as a 2-bit code (0 = 3' end, matching the legacy tadpole
    /// indexing; the underlying FastK key is indexed from the 5' end).
    pub(crate) fn base_at(&self, i: usize) -> u8 {
        self.0.base_at(self.0.k() - 1 - i)
    }

    /// Advance the window by one base: drop the 5' base, append `x`.
    pub(crate) fn push_right(&mut self, x: u8) {
        self.0.push_right(x);
    }

    /// Prepend `x` at the 5' end, dropping the 3' base.
    pub(crate) fn push_left(&mut self, x: u8) {
        self.0.push_left(x);
    }

    /// Reverse complement.
    pub(crate) fn rc(&self) -> Self {
        Self(self.0.rc())
    }

    /// Canonical key (lexicographically smaller of forward / reverse-complement).
    pub(crate) fn canonical(&self) -> Self {
        Self(self.0.canonical())
    }

    /// Lexicographic comparison of the base sequences (5' to 3').
    pub(crate) fn cmp_bases(&self, other: &Self) -> std::cmp::Ordering {
        self.0.cmp(&other.0)
    }
}

/// 2-bit base code, mirroring `AminoAcid.baseToNumber`: A=0, C=1, G=2,
/// T/U=3. Everything else maps to 0, so callers must gate on
/// [`base_defined`] first (the k-mer windows reset on N/ambiguity before
/// ever calling this).
pub(crate) fn base_code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' | b'U' | b'u' => 3,
        _ => 0,
    }
}

/// Reverse-complement code, mirroring `AminoAcid.baseToComplementNumber`:
/// A=3, C=2, G=1, T/U=0. Everything else maps to 0; callers must gate on
/// [`base_defined`] first.
pub(crate) fn base_comp_code(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 3,
        b'C' | b'c' => 2,
        b'G' | b'g' => 1,
        b'T' | b't' | b'U' | b'u' => 0,
        _ => 0,
    }
}

/// `AminoAcid.baseToNumber >= 0`: A/C/G/T/U count as defined (baseToNumber
/// is filled with -1 and only ACGTU are overwritten).
pub(crate) fn base_defined(b: u8) -> bool {
    matches!(
        b,
        b'A' | b'a' | b'C' | b'c' | b'G' | b'g' | b'T' | b't' | b'U' | b'u'
    )
}

/// `QualityTools.PROB_ERROR`: phred quality to error probability.
pub(crate) fn prob_error(q: u8) -> f32 {
    match q {
        0 => 0.75,
        1 => 0.7,
        _ => (10f64.powf(-0.1 * q as f64)) as f32,
    }
}

/// `QualityTools.PROB_CORRECT` and `PROB_CORRECT_INVERSE`, precomputed like
/// the Java arrays so the sliding `minprob` product uses exactly the same
/// float operations (multiply by a precomputed inverse, never divide).
fn prob_tables() -> (Vec<f32>, Vec<f32>) {
    let mut correct = Vec::with_capacity(128);
    let mut inverse = Vec::with_capacity(128);
    for q in 0..128u16 {
        let c = 1.0 - prob_error(q as u8);
        correct.push(c);
        inverse.push(1.0 / c);
    }
    (correct, inverse)
}

/// Emits the quality-gated canonical k-mers of one read as packed FastK
/// bytes, mirroring `KmerTableSetU.addKmersToTable` (canonical keys,
/// sliding `minprob` quality gate, N resets the window). Duplicates are
/// emitted repeatedly; `pgr::kmer::count::count_keys` groups them later.
///
/// Canonicalization uses the pgr `canonical_keys` trick: the forward window
/// and its reverse complement advance incrementally together (`push_right`
/// + `push_left(3-x)`), so the rc is never recomputed per window, and the
///   FastK mirror symmetry lets the canonical decision compare only the first
///   half of the packed bytes (`ceil(key_bytes/2)`). At long k this is the
///   dominant emission cost (k=100: `canonical()` recomputes a 25-byte rc per
///   window, 21.7% of runtime before this change).
fn count_read_kmers_packed(
    out: &mut Vec<u8>,
    bases: &[u8],
    quals: &[u8],
    k: usize,
    min_prob: f32,
    prob_correct: &[f32],
    prob_correct_inv: &[f32],
) {
    if bases.len() < k {
        return;
    }
    let min_prob2 = if min_prob > 0.0 && !quals.is_empty() {
        min_prob
    } else {
        0.0
    };
    let key_bytes = k.div_ceil(4);
    let half = key_bytes.div_ceil(2);
    let mut win = Kmer::new(k);
    let mut win_rc = Kmer::new(k); // rc of the all-zero window, advanced in lockstep
    let mut len = 0usize;
    let mut prob = 1f32;
    for (i, &b) in bases.iter().enumerate() {
        if base_defined(b) {
            let x = base_code(b);
            win.push_right(x);
            // Each new forward 3' base `x` prepends `3-x` to the rc's 5' end.
            win_rc.push_left(3 - x);
            if min_prob2 > 0.0 {
                // phred can exceed 127 on malformed (non-ASCII) quality bytes;
                // clamp so `prob_correct` (128 entries) is never indexed OOB.
                let q = (quals[i] as usize).min(127);
                prob *= prob_correct[q];
                if len >= k {
                    let oldq = (quals[i - k] as usize).min(127);
                    prob *= prob_correct_inv[oldq];
                }
            }
            len += 1;
        } else {
            len = 0;
            win.reset();
            win_rc.reset();
            prob = 1.0;
        }
        if len >= k && prob >= min_prob2 {
            let (fw, rc) = (win.0.to_bytes(), win_rc.0.to_bytes());
            out.extend_from_slice(if fw[..half] <= rc[..half] { fw } else { rc });
        }
    }
}

/// Applies the BBTools phred round-trip to a record's quality scores.
pub(crate) fn canonicalize_quality(rec: &mut SeqRecord) {
    if rec.quality_scores().is_empty() {
        return;
    }
    let seq = rec.sequence().to_vec();
    let raw = rec.quality_scores().to_vec();
    let phred = to_phred(&seq, &raw);
    rec.set_quality(from_phred(&phred));
}

#[cfg(test)]
mod tests {
    use super::*;
    use pgr::libs::nt::rev_comp;

    fn kmer_to_u128(k: &Kmer, kmer_len: usize) -> u128 {
        let mut x = 0u128;
        for i in (0..kmer_len).rev() {
            x = (x << 2) | k.base_at(i) as u128;
        }
        x
    }

    #[test]
    fn empty_reads_table_stays_queryable() {
        // A build from zero reads must still carry `k`: querying the empty
        // table (get_count / fill_* / sorted / solid snapshots) has to be
        // safe and empty instead of panicking on a defaulted k=0 table.
        let table = RefineTable::build(&[], 81, 0.5);
        let mut probe = Kmer::new(81);
        for &b in b"ACGTACGTACGTACGTACGTACGTACGTACGT" {
            probe.push_right(base_code(b));
        }
        assert_eq!(table.get_count(&probe), 0);
        assert_eq!(table.find_row(&probe), None);
        assert!(table.sorted_entries().is_empty());
        assert!(table.solid_entries(1).is_empty());
        assert!(table.solid_row_ranks(1).is_empty());
        assert!(table.fill_left_counts(&probe).iter().all(|&c| c == 0));
        assert!(table.fill_right_counts(&probe).iter().all(|&c| c == 0));
    }

    #[test]
    fn canonical_rc_is_identity() {
        // ACGT (0,1,2,3) RC is ACGT itself.
        let k = Kmer(key::Kmer::from_bases(b"ACGT", 4).unwrap());
        assert_eq!(k.rc().cmp_bases(&k), std::cmp::Ordering::Equal);
        // The canonical key is the lexicographically smaller orientation.
        let r = Kmer(key::Kmer::from_bases(b"GTCA", 4).unwrap());
        // RC(GTCA) = TGAC > GTCA, so canonical(r) must equal r itself.
        assert_eq!(r.canonical().cmp_bases(&r), std::cmp::Ordering::Equal);
        assert_eq!(r.canonical().cmp_bases(&r.rc()), std::cmp::Ordering::Less);
    }

    #[test]
    fn rolling_kmers_match_set_base_layout() {
        for k in [4usize, 31, 62, 81] {
            // Build "ACGT" repeated (truncated to k) by rolling push_right
            // and by set_base in the same orientation (base i of the window
            // occupies the low 2 bits after the window is full).
            let seq: Vec<u8> = (0..k).map(|i| b"ACGT"[i % 4]).collect();
            let direct = Kmer(key::Kmer::from_bases(&seq, k).unwrap());
            let mut rolled = Kmer::new(k);
            for &b in &seq {
                rolled.push_right(base_code(b));
            }
            assert_eq!(
                rolled.cmp_bases(&direct),
                std::cmp::Ordering::Equal,
                "k={k}"
            );

            // Rolling rc (push complements forward) must equal rc() of the
            // rolled kmer: this is the invariant extend_to_right2 relies on.
            let mut rolled_rc = Kmer::new(k);
            for &b in &seq {
                rolled_rc.push_left(base_comp_code(b));
            }
            assert_eq!(
                rolled_rc.cmp_bases(&rolled.rc()),
                std::cmp::Ordering::Equal,
                "rc-invariant k={k}"
            );

            // Canonical key must equal the lexicographically smaller of the
            // forward and rc orientations (u128 reference for k <= 62).
            if k <= 62 {
                let f = kmer_to_u128(&rolled, k);
                let r = kmer_to_u128(&rolled.rc(), k);
                assert_eq!(
                    kmer_to_u128(&rolled.canonical(), k),
                    f.min(r),
                    "canonical k={k}"
                );
            }

            // One rolling push_left (window full) must drop the oldest base
            // and prepend the new one.
            let mut rl = rolled_rc;
            rl.push_left(1);
            assert_eq!(rl.base_at(k - 1), 1, "push_left top k={k}");
            for i in 0..k - 1 {
                assert_eq!(
                    rl.base_at(i),
                    rolled_rc.base_at(i + 1),
                    "push_left shift k={k} i={i}"
                );
            }

            let mut rr = rolled;
            rr.push_right(2);
            assert_eq!(rr.base_at(0), 2, "push_right bottom k={k}");
            for i in 1..k {
                assert_eq!(
                    rr.base_at(i),
                    rolled.base_at(i - 1),
                    "push_right shift k={k} i={i}"
                );
            }
        }
    }

    #[test]
    fn counting_matches_simple_expected() {
        let reads = vec![(b"ACGTACGT".to_vec(), vec![40; 8])];
        let table = RefineTable::build(&reads, 4, 0.5);
        // 5 windows -> canonical keys ACGT(x2), CGTA(x2), GTAC(x1).
        assert_eq!(table.table.counts.iter().sum::<u32>(), 5);
        assert_eq!(table.table.keys.len(), 3); // 3 distinct packed keys (k=4 -> 1 B each)
        let mut probe = Kmer::new(4);
        for &b in b"ACGT" {
            probe.push_right(base_code(b));
        }
        assert_eq!(table.get_count(&probe), 2);
    }

    #[test]
    fn k81_table_counts_match_bruteforce() {
        // The merge phase-4 extension uses k=81 (Tadpole2 long-k path);
        // verify the k-mer table counts against a brute-force scan.
        let infile = "tests/bbtools/Lambda/golden/ext_sub.fq.gz";
        let mut reader = SeqReader::new(infile).unwrap();
        let mut rec = SeqRecord::new();
        let mut reads: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
        while reader.read_record(&mut rec).unwrap() {
            let seq = rec.sequence().to_vec();
            let quals = to_phred(&seq, rec.quality_scores());
            reads.push((seq, quals));
        }
        let k = 81usize;
        let table = RefineTable::build(&reads, k, 0.5);
        // Probe: the last 81 bases of the first read.
        let seq = &reads[0].0;
        let mut probe = Kmer::new(k);
        for &b in &seq[seq.len() - k..] {
            probe.push_right(base_code(b));
        }
        let probe_canon = probe.canonical();
        let (pc, pci) = prob_tables();
        let mut expected = 0u64;
        for (r, q) in &reads {
            let mut kk = Kmer::new(k);
            let mut len = 0usize;
            let mut prob = 1f32;
            for (i, &bb) in r.iter().enumerate() {
                if base_defined(bb) {
                    kk.push_right(base_code(bb));
                    prob *= pc[q[i] as usize];
                    if len >= k {
                        prob *= pci[q[i - k] as usize];
                    }
                    len += 1;
                } else {
                    len = 0;
                    kk.reset();
                    prob = 1.0;
                }
                if len >= k
                    && prob >= 0.5
                    && kk.canonical().cmp_bases(&probe_canon) == std::cmp::Ordering::Equal
                {
                    expected += 1;
                }
            }
        }
        assert_eq!(
            table.get_count(&probe),
            expected as u32,
            "k=81 table count mismatch"
        );
    }

    #[test]
    fn seed_kmer_count_symmetric() {
        // First Lambda read from ecco_sub.fq.gz (108 bp).
        let seq = b"AGAGATTCTTGGCGGAGAAACCATAATTGCATCTACTCGTCGCGAACCGCTTTCATCCGGCACAGTATCAAGGTATTTTATGCGCGCACGAAAAGCATC".to_vec();
        let quals = vec![40; seq.len()];
        let k = 62usize;
        let table = RefineTable::build(&[(seq.clone(), quals.clone())], k, 0.5);
        let rc: Vec<u8> = rev_comp(&seq).collect();
        for (label, s) in [("forward", &seq), ("rc", &rc)] {
            let mut kmer = Kmer::new(k);
            for &b in &s[s.len() - k..] {
                kmer.push_right(base_code(b));
            }
            eprintln!(
                "{label} tail kmer count={} words={:?} canonical={:?}",
                table.get_count(&kmer),
                kmer,
                kmer.canonical()
            );
        }
        // Directly compare the two canonical forms.
        let mut f = Kmer::new(k);
        for &b in &seq[seq.len() - k..] {
            f.push_right(base_code(b));
        }
        let mut r = Kmer::new(k);
        for &b in &rc[rc.len() - k..] {
            r.push_right(base_code(b));
        }
        eprintln!(
            "f.canonical={:?} r.canonical={:?} rc_of_f={:?}",
            f.canonical(),
            r.canonical(),
            f.rc()
        );
        // Canonical must be orientation-invariant: canonical(f) == canonical(rc(f)).
        let f_rc = f.rc();
        assert_eq!(
            f.canonical().cmp_bases(&f_rc.canonical()),
            std::cmp::Ordering::Equal,
            "canonical orientation-invariance broken"
        );

        // String-level check: rc() of a kmer must equal the reverse
        // complement of the sequence it encodes.
        for (label, s) in [("forward", &seq), ("rc", &rc)] {
            let mut kmer = Kmer::new(k);
            for &b in &s[s.len() - k..] {
                kmer.push_right(base_code(b));
            }
            let kmer_seq: Vec<u8> = (0..k).map(|i| kmer.base_at(i)).collect();
            let rc_seq: Vec<u8> = (0..k).map(|i| kmer.rc().base_at(i)).collect();
            eprintln!("{label} kmer_seq={kmer_seq:?} rc_seq={rc_seq:?}");
            // The rc of the kmer's base sequence (base_at order is 5'->3' as
            // built; verify the reversal+complement relationship explicitly).
            for i in 0..k {
                assert_eq!(
                    rc_seq[i],
                    3 - kmer_seq[k - 1 - i],
                    "{label} rc mismatch at {i}"
                );
            }
        }
    }
}

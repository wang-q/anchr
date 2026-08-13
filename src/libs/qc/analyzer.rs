//! M1 QC analyzer: statistical modules + FastQC-format text report.
//!
//! Numeric behavior follows fastqc 0.12.1 (verified against golden output
//! in `tests/qc/golden/`): integer average quality per read, `QualityCount`
//! percentile/mean over ASCII quality chars, GCModel-fraction GC histogram,
//! integer-division `%GC` in Basic Statistics.

use pgr::libs::fmt::seq::SeqRecord;
use std::collections::BTreeMap;
use std::io::Write;

use super::base_groups::make_base_groups;

// Phred bins 0..=63 (ASCII quality - offset); 64 slots keep the whole
// per-position array L1/L2-friendly (128 ASCII slots were 4× the footprint).
const QUAL_BINS: usize = 64;
const GC_BINS: usize = 101;
/// 4^7 possible 7-mers; also the column stride of the dense kmer grid.
const KMER_KEYS: usize = 1 << 14;

/// Per-base action table for the hot per-position loop. Value layout:
/// `(count_slot << 3) | pc_slot`; pc_slot 0-4 (A,C,G,T,N), 7 = skip;
/// count_slot 0-3 (A,C,G,T), 7 = no aggregate count. Replaces the 5-way
/// jump-table dispatch (indirect branches mispredict on random sequence).
static BASE_TABLE: [u8; 256] = {
    let mut t = [0x3Fu8; 256]; // 7<<3 | 7: skip both
    t[b'A' as usize] = 0;
    t[b'C' as usize] = (1 << 3) | 1;
    t[b'G' as usize] = (2 << 3) | 2;
    t[b'T' as usize] = (3 << 3) | 3;
    t[b'N' as usize] = (7 << 3) | 4;
    t
};

/// Per-position ASCII quality counts (index = ASCII char value).
#[derive(Clone, Copy)]
struct QualityCount {
    counts: [u32; QUAL_BINS],
}

impl Default for QualityCount {
    fn default() -> Self {
        QualityCount {
            counts: [0; QUAL_BINS],
        }
    }
}

impl QualityCount {
    /// Add one Phred quality value (ASCII minus the encoding offset).
    fn add(&mut self, phred: u8) {
        self.counts[(phred.min(QUAL_BINS as u8 - 1)) as usize] += 1;
    }

    fn total(&self) -> u64 {
        self.counts.iter().map(|&c| c as u64).sum()
    }

    fn mean(&self) -> f64 {
        let mut total = 0u64;
        let mut n = 0u64;
        for (i, &c) in self.counts.iter().enumerate() {
            total += i as u64 * c as u64;
            n += c as u64;
        }
        if n == 0 {
            0.0
        } else {
            total as f64 / n as f64
        }
    }

    /// fastqc `QualityCount.getPercentile`: integer `n*p/100` then the first
    /// quality whose cumulative count reaches the target (Phred value).
    fn percentile(&self, p: u64) -> f64 {
        let n = self.total();
        if n == 0 {
            return -1.0;
        }
        let target = n * p / 100;
        let mut count = 0u64;
        for (i, &c) in self.counts.iter().enumerate() {
            count += c as u64;
            if count >= target {
                return i as f64;
            }
        }
        -1.0
    }

    fn merge(&mut self, other: &QualityCount) {
        for (a, &b) in self.counts.iter_mut().zip(other.counts.iter()) {
            *a += b;
        }
    }
}

/// Per-group quality summary (mean of per-position percentiles, fastqc
/// semantics; `f64::NAN` when no position in the group had data).
#[derive(Clone, Copy)]
struct GroupQualityStats {
    mean: f64,
    median: f64,
    lq: f64,
    uq: f64,
    p10: f64,
    p90: f64,
}

impl GroupQualityStats {
    fn nan() -> Self {
        GroupQualityStats {
            mean: f64::NAN,
            median: f64::NAN,
            lq: f64::NAN,
            uq: f64::NAN,
            p10: f64::NAN,
            p90: f64::NAN,
        }
    }
}

/// Per-read-length GC model (fastqc `GCModel`): gc count -> (percentage,
/// fractional increment) pairs.
struct GcModel {
    models: Vec<Vec<(usize, f64)>>,
}

impl GcModel {
    fn new(read_length: u32) -> Self {
        let mut claiming = vec![0usize; GC_BINS];
        for pos in 0..=read_length {
            let (lp, hp) = Self::percent_range(read_length, pos);
            for p in lp..=hp {
                claiming[p] += 1;
            }
        }
        let mut models = vec![Vec::new(); read_length as usize + 1];
        for pos in 0..=read_length {
            let (lp, hp) = Self::percent_range(read_length, pos);
            for p in lp..=hp {
                models[pos as usize].push((p, 1.0 / claiming[p] as f64));
            }
        }
        GcModel { models }
    }

    fn percent_range(read_length: u32, pos: u32) -> (usize, usize) {
        let len = read_length as f64;
        let low = (pos as f64 - 0.5).clamp(0.0, len);
        let high = (pos as f64 + 0.5).clamp(0.0, len);
        let lp = ((low * 100.0) / len).round() as usize;
        let hp = ((high * 100.0) / len).round() as usize;
        (lp, hp)
    }

    fn values(&self, gc_count: u32) -> &[(usize, f64)] {
        &self.models[gc_count as usize]
    }
}

/// M1 QC statistics for one input file.
#[derive(Default)]
pub struct QcStats {
    filename: String,
    n_reads: u64,
    total_bases: u64,
    // aggregate ACGT counts (A, C, G, T); array form lets the hot loop use
    // a single table lookup + indexed increment instead of a 5-way jump
    // table (branch mispredictions were a top-2 stall source)
    base_counts: [u64; 4],
    len_min: u32,
    len_max: u32,
    // flat per-length counts (index = read length); a BTreeMap was a
    // per-read tree walk in the hot path
    len_hist: Vec<u64>,
    encoding_offset: u8,
    // per-position accumulators (index = 0-based position)
    per_base_quality: Vec<QualityCount>,
    per_base_content: Vec<[u32; 5]>, // A C G T N
    // per-read accumulators
    // integer avg ASCII quality -> reads (flat array: avg is 33..=126)
    seq_quality_hist: Vec<u64>,
    seq_qual_min: u8,
    seq_qual_max: u8,
    seq_gc_hist: Vec<f64>,
    gc_models: BTreeMap<u32, GcModel>,
    // M2: truncated-sequence counts (fastqc OverRepresentedSeqs; also feeds
    // the DuplicationLevel module)
    // (truncated length, fixed 50-byte sequence) — hash map: one hash per
    // read instead of a BTreeMap's ~log(n) memcmp comparisons per entry
    seq_counts: std::collections::HashMap<(u8, [u8; 50]), u64, FnvBuild>,
    // M2: adapter position counts (fastqc AdapterContent)
    adapter_positions: Vec<[u64; 6]>,
    longest_seq: u32,
    // M2: 7-mer position counts (fastqc KmerContent; sense strand only)
    // Dense 7-mer grid, row-major: `grid[key * kmer_stride + pos]`.
    // Replaces the per-position HashMap entry (hashbrown machinery was a
    // top-2 hotspot even with FNV); rows are contiguous so increments hit
    // sequential cache lines. When a longer read arrives the grid is
    // re-laid-out once (length steps are rare; copy cost negligible).
    kmer_grid: Vec<u32>,
    kmer_stride: usize,
    total_kmer_per_pos: Vec<u64>,
    // M3: per-tile quality (fastqc PerTileQualityScores)
    tile_ignore: bool,
    tile_split_pos: Option<usize>,
    tile_total_count: u64,
    tile_quality: std::collections::BTreeMap<u32, Vec<QualityCount>>,
}

impl QcStats {
    pub fn new(filename: &str, encoding_offset: u8) -> Self {
        QcStats {
            filename: filename.to_string(),
            encoding_offset,
            seq_gc_hist: vec![0.0; GC_BINS],
            seq_quality_hist: vec![0; 128],
            seq_qual_min: u8::MAX,
            seq_qual_max: 0,
            seq_counts: std::collections::HashMap::with_capacity_and_hasher(1 << 18, FnvBuild),
            kmer_grid: Vec::new(),
            kmer_stride: 0,
            ..Default::default()
        }
    }

    /// Consume one record (streaming; position vectors grow to max length).
    /// Consume one record; `global_index` is the 0-based index of the read
    /// in the whole file (used for the 2% kmer sampling, which must be
    /// global even when processing is split across parallel chunks).
    pub fn consume(&mut self, rec: &SeqRecord, global_index: u64) {
        self.consume_parts(
            rec.sequence(),
            rec.quality_scores(),
            rec.name(),
            global_index,
        );
    }

    /// Same as [`Self::consume`] but from borrowed record parts (zero-copy
    /// FASTQ parsing path; avoids building an owned `SeqRecord` per read).
    pub fn consume_parts(&mut self, seq: &[u8], qual: &[u8], name: &[u8], global_index: u64) {
        let len = seq.len() as u32;

        self.n_reads += 1;
        self.total_bases += len as u64;
        if self.len_hist.len() <= len as usize {
            self.len_hist.resize(len as usize + 1, 0);
        }
        self.len_hist[len as usize] += 1;
        if self.n_reads == 1 {
            self.len_min = len;
            self.len_max = len;
        } else {
            self.len_min = self.len_min.min(len);
            self.len_max = self.len_max.max(len);
        }

        if self.per_base_quality.len() < seq.len() {
            self.per_base_quality
                .resize(seq.len(), QualityCount::default());
            self.per_base_content.resize(seq.len(), [0u32; 5]);
        }
        // single pass over `seq`: per-base quality/content + aggregate
        // ACGT/GC counts + per-read average quality (was 3 passes)
        let trunc_gc = if seq.len() > 1000 {
            (seq.len() / 1000) * 1000
        } else if seq.len() > 100 {
            (seq.len() / 100) * 100
        } else {
            seq.len()
        };
        let mut qual_sum: u64 = 0;
        let mut gc_count = 0u32;
        let qual_len = qual.len();
        // Hoisted (LLVM reloads it from `self` every position otherwise),
        // and the quality path is split off: reads normally have
        // qual.len() == seq.len(), so the shared loop needs no per-position
        // `i < qual_len` compare; malformed shorter qualities take the tail.
        let offset = self.encoding_offset;
        let n = qual_len.min(seq.len());
        // per-position: quality slot + content slot + aggregate count via a
        // single table lookup (no 5-way jump table; the two guard branches
        // are ~always taken on ACGT data so they predict perfectly), and
        // the loop is unrolled 2x to halve per-position loop overhead.
        // The per_base vectors are resized to >= seq.len() above, so all
        // indices here are in bounds (checked via get_unchecked).
        let pq = &mut self.per_base_quality;
        let pc_arr = &mut self.per_base_content;
        let counts = &mut self.base_counts;
        macro_rules! process_pos {
            ($b:expr, $q:expr, $i:expr, $count_gc:expr) => {{
                unsafe { pq.get_unchecked_mut($i) }.add($q.saturating_sub(offset));
                qual_sum += $q as u64;
                let v = BASE_TABLE[$b as usize] as usize;
                if v & 7 != 7 {
                    unsafe {
                        *pc_arr
                            .get_unchecked_mut($i)
                            .get_unchecked_mut(v & 7) += 1;
                    }
                    let cs = v >> 3;
                    if cs != 7 {
                        unsafe { *counts.get_unchecked_mut(cs) += 1 };
                    }
                }
                if $count_gc && ($b == b'G' || $b == b'C') {
                    gc_count += 1;
                }
            }};
        }
        // Split at trunc_gc: only the first trunc_gc positions contribute to
        // the per-read GC count, so two loops remove the per-position
        // `i < trunc_gc` compare (fastqc truncates GC counting at 100/1000).
        let mut i = 0usize;
        let seq_p = seq.as_ptr();
        let qual_p = qual.as_ptr();
        // full positions with GC counting (first trunc_gc)
        while i + 2 <= n.min(trunc_gc) {
            process_pos!(unsafe { *seq_p.add(i) }, unsafe { *qual_p.add(i) }, i, true);
            process_pos!(
                unsafe { *seq_p.add(i + 1) },
                unsafe { *qual_p.add(i + 1) },
                i + 1,
                true
            );
            i += 2;
        }
        if i < n.min(trunc_gc) {
            process_pos!(unsafe { *seq_p.add(i) }, unsafe { *qual_p.add(i) }, i, true);
        }
        // remaining positions (beyond trunc_gc): no GC counting
        while i + 2 <= n {
            process_pos!(unsafe { *seq_p.add(i) }, unsafe { *qual_p.add(i) }, i, false);
            process_pos!(
                unsafe { *seq_p.add(i + 1) },
                unsafe { *qual_p.add(i + 1) },
                i + 1,
                false
            );
            i += 2;
        }
        if i < n {
            process_pos!(unsafe { *seq_p.add(i) }, unsafe { *qual_p.add(i) }, i, false);
        }
        for (i, &b) in seq[n..].iter().enumerate() {
            let pc = &mut self.per_base_content[n + i];
            match b {
                b'A' => {
                    pc[0] += 1;
                    self.base_counts[0] += 1;
                }
                b'C' => {
                    pc[1] += 1;
                    self.base_counts[1] += 1;
                }
                b'G' => {
                    pc[2] += 1;
                    self.base_counts[2] += 1;
                }
                b'T' => {
                    pc[3] += 1;
                    self.base_counts[3] += 1;
                }
                b'N' => pc[4] += 1,
                _ => {}
            }
            if n + i < trunc_gc && (b == b'G' || b == b'C') {
                gc_count += 1;
            }
        }

        // per-read average quality: integer division over ASCII chars
        if !qual.is_empty() {
            let avg = (qual_sum / qual.len() as u64) as i32;
            self.seq_quality_hist[avg as usize] += 1;
            self.seq_qual_min = self.seq_qual_min.min(avg as u8);
            self.seq_qual_max = self.seq_qual_max.max(avg as u8);
        }

        // per-read GC: GCModel fractions from the in-pass GC count
        let model = self
            .gc_models
            .entry(trunc_gc as u32)
            .or_insert_with(|| GcModel::new(trunc_gc as u32));
        for &(p, inc) in model.values(gc_count) {
            self.seq_gc_hist[p] += inc;
        }

        // overrep/duplication: truncate to 50 bp (fastqc default); fixed
        // [u8; 50] key avoids a per-read Vec allocation (hot path)
        let trunc = if seq.len() > 50 { &seq[..50] } else { seq };
        let mut key = [0u8; 50];
        key[..trunc.len()].copy_from_slice(trunc);
        *self.seq_counts.entry((trunc.len() as u8, key)).or_insert(0) += 1;

        // adapter content: indexOf each adapter, increment to the current
        // longest-read bound (fastqc semantics)
        if seq.len() as u32 > self.longest_seq {
            self.longest_seq = seq.len() as u32;
            let cur_max = self.longest_seq.saturating_sub(12) + 1;
            self.adapter_positions.resize(cur_max as usize, [0u64; 6]);
        }
        let cur_max = self.longest_seq.saturating_sub(12) + 1;
        // Single pass over the read locating the first occurrence of each
        // of the 6 adapters (their leading 2 bytes are pairwise distinct,
        // so one match branch suffices). Early-exits once all are found;
        // replaces 6 independent `find_subseq` scans per read.
        let ads = adapters();
        let sigs = adapter_sigs();
        let table = adapter_pair_table();
        let mut found = [None; 6];
        // bitmask of adapters still missing: single register check per
        // position (the old `found.iter().any()` re-read 6 stack slots)
        let mut mask = 0b111111u8;
        let l = seq.len();
        let mut i = 0;
        // 4 positions per iteration: two unaligned u32 loads cover the
        // overlapping pairs at i..i+3; quarters the loop-control overhead
        while i + 6 <= l && mask != 0 {
            let p0 = unsafe { read_u32_un(seq, i) };
            let p1 = unsafe { read_u32_un(seq, i + 2) };
            let a = adapter_at(table, p0 as u16);
            if a != 0xF {
                verify_adapter(a as usize, i, l, seq, sigs, ads, &mut found, &mut mask);
            }
            if mask != 0 {
                let a = adapter_at(table, (p0 >> 8) as u16);
                if a != 0xF {
                    verify_adapter(a as usize, i + 1, l, seq, sigs, ads, &mut found, &mut mask);
                }
            }
            if mask != 0 {
                let a = adapter_at(table, p1 as u16);
                if a != 0xF {
                    verify_adapter(a as usize, i + 2, l, seq, sigs, ads, &mut found, &mut mask);
                }
            }
            if mask != 0 {
                let a = adapter_at(table, (p1 >> 8) as u16);
                if a != 0xF {
                    verify_adapter(a as usize, i + 3, l, seq, sigs, ads, &mut found, &mut mask);
                }
            }
            i += 4;
        }
        // 2 positions per iteration (remaining 2-3 tail positions)
        while i + 4 <= l && mask != 0 {
            let p = unsafe { read_u32_un(seq, i) };
            let a = adapter_at(table, p as u16);
            if a != 0xF {
                verify_adapter(a as usize, i, l, seq, sigs, ads, &mut found, &mut mask);
            }
            if mask != 0 {
                let a = adapter_at(table, (p >> 8) as u16);
                if a != 0xF {
                    verify_adapter(a as usize, i + 1, l, seq, sigs, ads, &mut found, &mut mask);
                }
            }
            i += 2;
        }
        // tail: one remaining pair (odd length)
        if i + 2 <= l && mask != 0 {
            let a = adapter_at(table, u16::from_le_bytes([seq[i], seq[i + 1]]));
            if a != 0xF {
                verify_adapter(a as usize, i, l, seq, sigs, ads, &mut found, &mut mask);
            }
        }
        for (a, pos) in found.iter().enumerate() {
            if let Some(idx) = pos {
                for p in *idx..cur_max as usize {
                    self.adapter_positions[p][a] += 1;
                }
            }
        }

        // 7-mer position counts: fastqc samples every 50th read (2%) and
        // skips kmers containing Ns
        if (global_index + 1) % 50 == 0 && seq.len() >= 7 {
            let seq = if seq.len() > 500 { &seq[..500] } else { seq };
            let l = seq.len();
            if self.total_kmer_per_pos.len() < l - 6 {
                self.ensure_kmer_positions(l - 6);
            }
            // Sliding 7-mer window: O(1) per position instead of re-scanning
            // 7 bytes per window (E. coli has almost no Ns, so the old
            // `any()` scan rarely exited early; ~7x fewer byte loads).
            let mut code = 0u16;
            let mut non_acgt = 0u32;
            for &b in &seq[..7] {
                code = (code << 2) | encode_base(b);
                non_acgt += u32::from(!matches!(b, b'A' | b'C' | b'G' | b'T'));
            }
            let mut i = 0;
            while i + 7 <= l {
                if non_acgt == 0 {
                    self.kmer_grid[code as usize * self.kmer_stride + i] += 1;
                    self.total_kmer_per_pos[i] += 1;
                }
                if i + 7 == l {
                    break;
                }
                let out = seq[i];
                let inc = seq[i + 7];
                code = ((code << 2) | encode_base(inc)) & 0x3FFF;
                non_acgt += u32::from(!matches!(inc, b'A' | b'C' | b'G' | b'T'));
                non_acgt -= u32::from(!matches!(out, b'A' | b'C' | b'G' | b'T'));
                i += 1;
            }
        }

        // per-tile quality: fastqc samples all of the first 10k reads, then
        // every 10th; tile comes from the Casava 1.8+ (field 4) or legacy
        // (field 2) header
        if !self.tile_ignore && !qual.is_empty() {
            self.tile_total_count += 1;
            if self.tile_total_count > 10_000 && self.tile_total_count % 10 != 0 {
                // skip sampling
            } else if let Some(tile) = self.detect_tile(name) {
                let entry = self.tile_quality.entry(tile).or_default();
                if entry.len() < qual.len() {
                    entry.resize(qual.len(), QualityCount::default());
                }
                for (i, &q) in qual.iter().enumerate() {
                    entry[i].add(q.saturating_sub(self.encoding_offset));
                }
            }
        }
    }

    /// Encoding label from the detected offset (fastqc style). The offset
    /// is provided by the command layer via `pgr::libs::fq::qual::detect_quality_base`.
    pub fn encoding_label(&self) -> &str {
        if self.encoding_offset == 64 {
            "Solexa / Illumina 1.8"
        } else {
            "Sanger / Illumina 1.9"
        }
    }

    pub fn n_reads(&self) -> u64 {
        self.n_reads
    }

    fn detect_tile(&mut self, name: &[u8]) -> Option<u32> {
        // Split on ':' without allocating a Vec of fields (hot path: every
        // sampled read calls this). Caches the field index like fastqc:
        // Casava 1.8+ (>= 7 fields) uses field 4, legacy (>= 5) field 2.
        let pos = match self.tile_split_pos {
            Some(p) => p,
            None => {
                let colons = name.iter().filter(|&&c| c == b':').count();
                if colons >= 6 {
                    self.tile_split_pos = Some(4);
                    4
                } else if colons >= 4 {
                    self.tile_split_pos = Some(2);
                    2
                } else {
                    self.tile_ignore = true;
                    return None;
                }
            }
        };
        let mut field = 0usize;
        let mut start = 0usize;
        for (i, &c) in name.iter().enumerate() {
            if c == b':' {
                if field == pos {
                    return parse_u32_bytes(&name[start..i]);
                }
                field += 1;
                start = i + 1;
            }
        }
        if field == pos {
            return parse_u32_bytes(&name[start..]);
        }
        None
    }

    /// Merge another statistics block into `self` (parallel chunk results).
    pub fn merge(&mut self, other: &QcStats) {
        let was_empty = self.n_reads == 0;
        self.n_reads += other.n_reads;
        self.total_bases += other.total_bases;
        for (a, &b) in self.base_counts.iter_mut().zip(&other.base_counts) {
            *a += b;
        }
        if other.n_reads > 0 {
            if was_empty {
                self.len_min = other.len_min;
                self.len_max = other.len_max;
            } else {
                self.len_min = self.len_min.min(other.len_min);
                self.len_max = self.len_max.max(other.len_max);
            }
        }
        if other.len_hist.len() > self.len_hist.len() {
            self.len_hist.resize(other.len_hist.len(), 0);
        }
        for (a, &b) in self.len_hist.iter_mut().zip(&other.len_hist) {
            *a += b;
        }
        if other.per_base_quality.len() > self.per_base_quality.len() {
            self.per_base_quality
                .resize(other.per_base_quality.len(), QualityCount::default());
            self.per_base_content
                .resize(other.per_base_content.len(), [0u32; 5]);
        }
        for (a, b) in self.per_base_quality.iter_mut().zip(&other.per_base_quality) {
            a.merge(b);
        }
        for (a, b) in self.per_base_content.iter_mut().zip(&other.per_base_content) {
            for k in 0..5 {
                a[k] += b[k];
            }
        }
        for (a, &b) in self.seq_quality_hist.iter_mut().zip(&other.seq_quality_hist) {
            *a += b;
        }
        if other.seq_qual_min != u8::MAX {
            self.seq_qual_min = self.seq_qual_min.min(other.seq_qual_min);
            self.seq_qual_max = self.seq_qual_max.max(other.seq_qual_max);
        }
        for (a, &b) in self.seq_gc_hist.iter_mut().zip(&other.seq_gc_hist) {
            *a += b;
        }
        for (&key, &c) in &other.seq_counts {
            *self.seq_counts.entry(key).or_insert(0) += c;
        }
        if other.adapter_positions.len() > self.adapter_positions.len() {
            self.adapter_positions
                .resize(other.adapter_positions.len(), [0u64; 6]);
            self.longest_seq = other.longest_seq;
        }
        for (a, b) in self.adapter_positions.iter_mut().zip(&other.adapter_positions) {
            for k in 0..6 {
                a[k] += b[k];
            }
        }
        if other.total_kmer_per_pos.len() > self.total_kmer_per_pos.len() {
            self.ensure_kmer_positions(other.total_kmer_per_pos.len());
        }
        for (a, &b) in self.total_kmer_per_pos.iter_mut().zip(&other.total_kmer_per_pos) {
            *a += b;
        }
        // row-major grids: add each key's overlapping positions in place
        let s = self.kmer_stride;
        let os = other.kmer_stride;
        for key in 0..KMER_KEYS {
            for i in 0..os {
                self.kmer_grid[key * s + i] += other.kmer_grid[key * os + i];
            }
        }
    }

    /// Grow the kmer grid to `positions` columns, re-laying existing rows
    /// (the old stride only changes when a longer read arrives).
    fn ensure_kmer_positions(&mut self, positions: usize) {
        self.total_kmer_per_pos.resize(positions, 0);
        let old_stride = self.kmer_stride;
        if old_stride == positions {
            return;
        }
        let mut new_grid = vec![0u32; KMER_KEYS * positions];
        if old_stride > 0 {
            for key in 0..KMER_KEYS {
                new_grid[key * positions..key * positions + old_stride].copy_from_slice(
                    &self.kmer_grid[key * old_stride..key * old_stride + old_stride],
                );
            }
        }
        self.kmer_grid = new_grid;
        self.kmer_stride = positions;
    }

    /// Per-group quality stats: each position's percentile/mean is computed
    /// independently (positions with >100 quality chars), then averaged
    /// across the group (fastqc `PerBaseQualityScores.getPercentile/getMean`).
    fn group_quality_stats(&self) -> Vec<GroupQualityStats> {
        let groups = make_base_groups(self.per_base_quality.len() as u32);
        let mut out = Vec::with_capacity(groups.len());
        for g in &groups {
            let mut count = 0usize;
            let mut mean = 0.0;
            let mut median = 0.0;
            let mut lq = 0.0;
            let mut uq = 0.0;
            let mut p10 = 0.0;
            let mut p90 = 0.0;
            for p in (g.start - 1)..g.end {
                let q = &self.per_base_quality[p as usize];
                if q.total() > 100 {
                    count += 1;
                    mean += q.mean();
                    median += q.percentile(50);
                    lq += q.percentile(25);
                    uq += q.percentile(75);
                    p10 += q.percentile(10);
                    p90 += q.percentile(90);
                }
            }
            if count > 0 {
                let n = count as f64;
                out.push(GroupQualityStats {
                    mean: mean / n,
                    median: median / n,
                    lq: lq / n,
                    uq: uq / n,
                    p10: p10 / n,
                    p90: p90 / n,
                });
            } else {
                out.push(GroupQualityStats::nan());
            }
        }
        out
    }

    /// Aggregated per-group base counts as [A, C, G, T, N].
    fn group_contents(&self) -> Vec<[u32; 5]> {
        let groups = make_base_groups(self.per_base_content.len() as u32);
        let mut out = vec![[0u32; 5]; groups.len()];
        for (i, g) in groups.iter().enumerate() {
            for p in (g.start - 1)..g.end {
                for k in 0..5 {
                    out[i][k] += self.per_base_content[p as usize][k];
                }
            }
        }
        out
    }

    fn group_labels(&self) -> Vec<String> {
        make_base_groups(self.per_base_quality.len() as u32)
            .iter()
            .map(|g| g.to_string())
            .collect()
    }

    fn total_acgtn(&self) -> u64 {
        self.base_counts.iter().sum()
    }

    /// Render `fastqc_data.txt` (fastqc 0.12.1 compatible).
    pub fn report_txt(&self, w: &mut dyn Write) -> std::io::Result<()> {
        let offset = self.encoding_offset;

        writeln!(w, "##FastQC\t0.12.1")?;
        writeln!(w, ">>Basic Statistics\tpass")?;
        writeln!(w, "#Measure\tValue")?;
        writeln!(w, "Filename\t{}", self.filename)?;
        writeln!(w, "File type\tConventional base calls")?;
        writeln!(w, "Encoding\t{}", self.encoding_label())?;
        writeln!(w, "Total Sequences\t{}", self.n_reads)?;
        writeln!(w, "Total Bases\t{}", format_length(self.total_bases))?;
        writeln!(w, "Sequences flagged as poor quality\t0")?;
        let len_str = if self.len_min == self.len_max {
            format!("{}", self.len_min)
        } else {
            format!("{}-{}", self.len_min, self.len_max)
        };
        writeln!(w, "Sequence length\t{}", len_str)?;
        let at = self.base_counts[0] + self.base_counts[3];
        let gc = self.base_counts[2] + self.base_counts[1];
        let gc_pct = if at + gc > 0 {
            (gc * 100) / (at + gc)
        } else {
            0
        };
        writeln!(w, "%GC\t{}", gc_pct)?;
        writeln!(w, ">>END_MODULE")?;

        // Per base sequence quality
        let qualities = self.group_quality_stats();
        writeln!(
            w,
            ">>Per base sequence quality\t{}",
            self.grade_per_base_quality(&qualities)
        )?;
        writeln!(
            w,
            "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile"
        )?;
        for (i, label) in self.group_labels().iter().enumerate() {
            let q = &qualities[i];
            writeln!(
                w,
                "{}\t{}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{:.1}",
                label,
                q.mean,
                q.median,
                q.lq,
                q.uq,
                q.p10,
                q.p90,
            )?;
        }
        writeln!(w, ">>END_MODULE")?;

        // Per sequence quality scores
        writeln!(
            w,
            ">>Per sequence quality scores\t{}",
            self.grade_per_sequence_quality()
        )?;
        writeln!(w, "#Quality\tCount")?;
        if self.seq_qual_min <= self.seq_qual_max {
            for q in self.seq_qual_min..=self.seq_qual_max {
                let count = self.seq_quality_hist[q as usize];
                writeln!(w, "{}\t{:.1}", q as i32 - offset as i32, count as f64)?;
            }
        }
        writeln!(w, ">>END_MODULE")?;

        // Per base sequence content
        let contents = self.group_contents();
        writeln!(
            w,
            ">>Per base sequence content\t{}",
            self.grade_per_base_content(&contents)
        )?;
        writeln!(w, "#Base\tG\tA\tT\tC")?;
        for (i, label) in self.group_labels().iter().enumerate() {
            let [a, c, g, t, _n] = contents[i];
            let total = a + c + g + t;
            let (g, a, t, c) = if total > 0 {
                (
                    100.0 * g as f64 / total as f64,
                    100.0 * a as f64 / total as f64,
                    100.0 * t as f64 / total as f64,
                    100.0 * c as f64 / total as f64,
                )
            } else {
                (0.0, 0.0, 0.0, 0.0)
            };
            writeln!(w, "{}\t{}\t{}\t{}\t{}", label, g, a, t, c)?;
        }
        writeln!(w, ">>END_MODULE")?;

        // Per sequence GC content
        writeln!(
            w,
            ">>Per sequence GC content\t{}",
            self.grade_per_sequence_gc()
        )?;
        writeln!(w, "#GC Content\tCount")?;
        for (i, &count) in self.seq_gc_hist.iter().enumerate() {
            // fastqc outputs the GCModel-weighted read count directly
            // (sum over bins == total reads)
            writeln!(w, "{}\t{}", i, count)?;
        }
        writeln!(w, ">>END_MODULE")?;

        // Per base N content
        writeln!(
            w,
            ">>Per base N content\t{}",
            self.grade_per_base_n(&contents)
        )?;
        writeln!(w, "#Base\tN-Count")?;
        for (i, label) in self.group_labels().iter().enumerate() {
            let [a, c, g, t, n] = contents[i];
            let total = a + c + g + t + n;
            let pct = if total > 0 { 100.0 * n as f64 / total as f64 } else { 0.0 };
            writeln!(w, "{}\t{}", label, pct)?;
        }
        writeln!(w, ">>END_MODULE")?;

        // Sequence Length Distribution
        writeln!(
            w,
            ">>Sequence Length Distribution\t{}",
            self.grade_length()
        )?;
        writeln!(w, "#Length\tCount")?;
        for (len, &count) in self.len_hist.iter().enumerate().filter(|(_, &c)| c != 0) {
            writeln!(w, "{}\t{:.1}", len, count as f64)?;
        }
        writeln!(w, ">>END_MODULE")?;

        // Sequence Duplication Levels (M2)
        let (dedup_pct, slots) = self.duplication_data();
        writeln!(
            w,
            ">>Sequence Duplication Levels\t{}",
            self.grade_duplication(dedup_pct)
        )?;
        writeln!(w, "#Total Deduplicated Percentage\t{}", dedup_pct)?;
        writeln!(w, "#Duplication Level\tPercentage of total")?;
        for (i, &v) in slots.iter().enumerate() {
            writeln!(w, "{}\t{}", dup_label(i), v)?;
        }
        writeln!(w, ">>END_MODULE")?;

        // Overrepresented sequences (M2)
        let rows = self.overrep_rows();
        writeln!(
            w,
            ">>Overrepresented sequences\t{}",
            self.grade_overrep(&rows)
        )?;
        writeln!(w, "#Sequence\tCount\tPercentage\tPossible Source")?;
        for (seq, count, pct, source) in &rows {
            writeln!(
                w,
                "{}\t{}\t{}\t{}",
                String::from_utf8_lossy(seq),
                count,
                pct,
                source
            )?;
        }
        writeln!(w, ">>END_MODULE")?;

        // Adapter Content (M2)
        let adapter_groups = self.adapter_enrichment();
        let adapter_grade = self.grade_adapter(&adapter_groups);
        writeln!(w, ">>Adapter Content\t{}", adapter_grade)?;
        write!(
            w,
            "#Position\t{}",
            adapters()
                .iter()
                .map(|(n, _)| *n)
                .collect::<Vec<_>>()
                .join("\t")
        )?;
        writeln!(w)?;
        let groups = make_base_groups(self.adapter_positions.len() as u32);
        for (i, g) in groups.iter().enumerate() {
            let row = adapter_groups[i];
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                g,
                row[0],
                row[1],
                row[2],
                row[3],
                row[4],
                row[5],
            )?;
        }
        writeln!(w, ">>END_MODULE")?;

        // Kmer Content (M2; fastqc emits no section when nothing is
        // significantly enriched — as on the Lambda golden)
        let kmers = self.kmer_content();
        if !kmers.is_empty() {
            writeln!(w, ">>Kmer Content\tfail")?;
            writeln!(w, "#Sequence\tCount\tObs/Exp Max\tMax Obs/Exp Position")?;
            for (key, count, obs_exp_max, pos_group) in &kmers {
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}",
                    kmer_to_string(*key),
                    count,
                    obs_exp_max,
                    pos_group
                )?;
            }
            writeln!(w, ">>END_MODULE")?;
        }

        // Per tile sequence quality (M3; absent when no Casava-style tile
        // headers — as on the Lambda golden)
        if !self.tile_ignore && !self.tile_quality.is_empty() {
            let (grade, global_mean) = self.grade_tile();
            writeln!(w, ">>Per tile sequence quality\t{}", grade)?;
            writeln!(w, "#Tile\tBase Group\tMean Quality")?;
            for (tile, quals) in &self.tile_quality {
                let groups = make_base_groups(quals.len() as u32);
                for g in &groups {
                    let mut q = QualityCount::default();
                    for p in (g.start - 1)..g.end {
                        q.merge(&quals[p as usize]);
                    }
                    writeln!(
                        w,
                        "{}\t{}\t{}",
                        tile,
                        g,
                        q.mean() - global_mean
                    )?;
                }
            }
            writeln!(w, ">>END_MODULE")?;
        }

        Ok(())
    }

    /// Render `summary.txt`: `STATUS\tModule\tFilename`.
    pub fn report_summary(&self, w: &mut dyn Write) -> std::io::Result<()> {
        let qualities = self.group_quality_stats();
        let contents = self.group_contents();
        let (dedup_pct, _slots) = self.duplication_data();
        let overrep = self.overrep_rows();
        let rows = [
            ("Basic Statistics", "pass"),
            ("Per base sequence quality", self.grade_per_base_quality(&qualities)),
            ("Per sequence quality scores", self.grade_per_sequence_quality()),
            ("Per base sequence content", self.grade_per_base_content(&contents)),
            ("Per sequence GC content", self.grade_per_sequence_gc()),
            ("Per base N content", self.grade_per_base_n(&contents)),
            ("Sequence Length Distribution", self.grade_length()),
            ("Sequence Duplication Levels", self.grade_duplication(dedup_pct)),
            ("Overrepresented sequences", self.grade_overrep(&overrep)),
            (
                "Adapter Content",
                self.grade_adapter(&self.adapter_enrichment()),
            ),
        ];
        let mut rows = rows.to_vec();
        if !self.tile_ignore && !self.tile_quality.is_empty() {
            rows.push(("Per tile sequence quality", self.grade_tile().0));
        }
        for (name, grade) in rows {
            writeln!(w, "{}\t{}\t{}", grade.to_uppercase(), name, self.filename)?;
        }
        Ok(())
    }

    /// Render `fastqc_report.html` (tera + inline SVG; structure-compatible
    /// with fastqc, not pixel-identical).
    pub fn report_html(&self) -> anyhow::Result<String> {
        use serde_json::json;
        use tera::{Context, Tera};

        let qualities = self.group_quality_stats();
        let contents = self.group_contents();
        let labels = self.group_labels();
        let (dedup_pct, dup_slots) = self.duplication_data();
        let overrep = self.overrep_rows();
        let adapter_groups = self.adapter_enrichment();
        let adapter_names: Vec<&str> = adapters().iter().map(|(n, _)| *n).collect();

        let mut ctx = Context::new();
        ctx.insert("filename", &self.filename);
        let len_str = if self.len_min == self.len_max {
            format!("{}", self.len_min)
        } else {
            format!("{}-{}", self.len_min, self.len_max)
        };
        let at = self.base_counts[0] + self.base_counts[3];
        let gc = self.base_counts[2] + self.base_counts[1];
        let gc_pct = if at + gc > 0 { (gc * 100) / (at + gc) } else { 0 };
        let basic: Vec<serde_json::Value> = [
            ("Filename", self.filename.clone()),
            ("File type", "Conventional base calls".to_string()),
            ("Encoding", self.encoding_label().to_string()),
            ("Total Sequences", self.n_reads.to_string()),
            ("Total Bases", format_length(self.total_bases)),
            ("Sequences flagged as poor quality", "0".to_string()),
            ("Sequence length", len_str),
            ("%GC", gc_pct.to_string()),
        ]
        .iter()
        .map(|(k, v)| json!([k, v]))
        .collect();
        ctx.insert("basic", &basic);

        ctx.insert("basic_grade", "pass");
        ctx.insert(
            "pbq_grade",
            self.grade_per_base_quality(&qualities),
        );
        ctx.insert(
            "psq_grade",
            self.grade_per_sequence_quality(),
        );
        ctx.insert(
            "pbc_grade",
            self.grade_per_base_content(&contents),
        );
        ctx.insert("gc_grade", self.grade_per_sequence_gc());
        ctx.insert("n_grade", self.grade_per_base_n(&contents));
        ctx.insert("len_grade", self.grade_length());
        ctx.insert("dup_grade", self.grade_duplication(dedup_pct));
        ctx.insert("overrep_grade", self.grade_overrep(&overrep));
        ctx.insert(
            "adapter_grade",
            self.grade_adapter(&adapter_groups),
        );

        let offset = self.encoding_offset;
        let n_groups = qualities.len().max(1);
        let pbq: Vec<serde_json::Value> = qualities
            .iter()
            .enumerate()
            .map(|(i, q)| {
                json!({
                    "label": labels[i],
                    "x": i as f64 * 800.0 / n_groups as f64,
                    "mean_y": 300.0 - q.mean * 7.0,
                    "median_y": 300.0 - q.median * 7.0,
                    "mean": q.mean,
                    "median": q.median,
                    "lq": q.lq,
                    "uq": q.uq,
                    "p10": q.p10,
                    "p90": q.p90,
                })
            })
            .collect();
        ctx.insert("per_base_quality", &pbq);

        let psq: Vec<serde_json::Value> = (self.seq_qual_min..=self.seq_qual_max)
            .map(|q| json!([q as i32 - offset as i32, self.seq_quality_hist[q as usize] as f64]))
            .collect();
        ctx.insert("per_seq_quality", &psq);

        let pbc: Vec<serde_json::Value> = contents
            .iter()
            .zip(&labels)
            .map(|([a, c, g, t, _], label)| {
                let total = a + c + g + t;
                let (g, a, t, c) = if total > 0 {
                    (
                        100.0 * *g as f64 / total as f64,
                        100.0 * *a as f64 / total as f64,
                        100.0 * *t as f64 / total as f64,
                        100.0 * *c as f64 / total as f64,
                    )
                } else {
                    (0.0, 0.0, 0.0, 0.0)
                };
                json!([label, g, a, t, c])
            })
            .collect();
        ctx.insert("per_base_content", &pbc);

        let gc_rows: Vec<serde_json::Value> = self
            .seq_gc_hist
            .iter()
            .enumerate()
            .map(|(i, &v)| json!([i, (v * 180.0 / self.n_reads.max(1) as f64).min(180.0)]))
            .collect();
        ctx.insert("gc", &gc_rows);

        let n_rows: Vec<serde_json::Value> = contents
            .iter()
            .zip(&labels)
            .map(|([a, c, g, t, n], label)| {
                let total = a + c + g + t + n;
                let pct = if total > 0 {
                    100.0 * *n as f64 / total as f64
                } else {
                    0.0
                };
                json!([label, pct])
            })
            .collect();
        ctx.insert("per_base_n", &n_rows);

        let len_rows: Vec<serde_json::Value> = self
            .len_hist
            .iter()
            .enumerate()
            .filter(|(_, &c)| c != 0)
            .map(|(l, &c)| json!([l, c as f64]))
            .collect();
        ctx.insert("length", &len_rows);
        ctx.insert("dup_total", &dedup_pct);
        let dup_rows: Vec<serde_json::Value> = dup_slots
            .iter()
            .enumerate()
            .map(|(i, &v)| json!([dup_label(i), v]))
            .collect();
        ctx.insert("dup", &dup_rows);

        let overrep_rows: Vec<serde_json::Value> = overrep
            .iter()
            .map(|(s, c, p, src)| {
                json!([String::from_utf8_lossy(s), *c, *p, src])
            })
            .collect();
        ctx.insert("overrep", &overrep_rows);
        ctx.insert("adapter_names", &adapter_names);
        let adapter_rows: Vec<serde_json::Value> = adapter_groups
            .iter()
            .zip(&labels)
            .map(|(vals, label)| json!([label, vals]))
            .collect();
        ctx.insert("adapter", &adapter_rows);

        let mut tera = Tera::default();
        tera.add_raw_templates(vec![(
            "report",
            include_str!("../../../templates/qc/report.html"),
        )])?;
        Ok(tera.render("report", &ctx)?)
    }

    fn grade_per_base_quality(&self, qualities: &[GroupQualityStats]) -> &'static str {
        for q in qualities {
            if q.median < 20.0 || q.lq < 5.0 {
                return "fail";
            }
        }
        for q in qualities {
            if q.median < 25.0 || q.lq < 10.0 {
                return "warn";
            }
        }
        "pass"
    }

    fn grade_per_sequence_quality(&self) -> &'static str {
        let mut total = 0u64;
        let mut n = 0u64;
        if self.seq_qual_min <= self.seq_qual_max {
            for q in self.seq_qual_min..=self.seq_qual_max {
                let c = self.seq_quality_hist[q as usize];
                total += q as u64 * c;
                n += c;
            }
        }
        if n == 0 {
            return "pass";
        }
        let mean = total as f64 / n as f64 - self.encoding_offset as f64;
        if mean < 20.0 {
            "fail"
        } else if mean < 27.0 {
            "warn"
        } else {
            "pass"
        }
    }

    fn grade_per_base_content(&self, contents: &[[u32; 5]]) -> &'static str {
        for &[a, c, g, t, _n] in contents {
            let total = a + c + g + t;
            if total == 0 {
                continue;
            }
            let gc = 100.0 * (g + c) as f64 / total as f64;
            let dev = (gc - 50.0).abs();
            if dev > 20.0 {
                return "fail";
            }
            if dev > 10.0 {
                return "warn";
            }
        }
        "pass"
    }

    fn grade_per_base_n(&self, contents: &[[u32; 5]]) -> &'static str {
        for &[a, c, g, t, n] in contents {
            let total = a + c + g + t + n;
            if total == 0 {
                continue;
            }
            let pct = 100.0 * n as f64 / total as f64;
            if pct > 20.0 {
                return "fail";
            }
            if pct > 5.0 {
                return "warn";
            }
        }
        "pass"
    }

    /// fastqc `PerSequenceGCContent` grade: theoretical Normal(mode, stdev)
    /// distribution vs observed; deviation > 30% fail, > 15% warn.
    fn grade_per_sequence_gc(&self) -> &'static str {
        let dist = &self.seq_gc_hist;
        let total: f64 = dist.iter().sum();
        if total <= 1.0 {
            return "pass";
        }
        // mode = average of adjacent bins above 90% of the modal value
        let first_mode = dist
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0);
        let threshold = dist[first_mode] - dist[first_mode] / 10.0;
        let mut mode = 0.0;
        let mut n = 0.0;
        let mut fell_off_top = true;
        for (i, &v) in dist.iter().enumerate().skip(first_mode) {
            if v > threshold {
                mode += i as f64;
                n += 1.0;
            } else {
                fell_off_top = false;
                break;
            }
        }
        let mut fell_off_bottom = true;
        for i in (0..first_mode).rev() {
            if dist[i] > threshold {
                mode += i as f64;
                n += 1.0;
            } else {
                fell_off_bottom = false;
                break;
            }
        }
        let mode = if fell_off_bottom || fell_off_top {
            first_mode as f64
        } else {
            mode / n
        };

        let mut var = 0.0;
        for (i, &v) in dist.iter().enumerate() {
            var += (i as f64 - mode).powi(2) * v;
        }
        let stdev = (var / (total - 1.0)).sqrt();
        let mut deviation = 0.0;
        for (i, &v) in dist.iter().enumerate() {
            let pdf = normal_pdf(i as f64, mode, stdev);
            deviation += (pdf * total - v).abs();
        }
        let deviation = deviation / total * 100.0;
        if deviation > 30.0 {
            "fail"
        } else if deviation > 15.0 {
            "warn"
        } else {
            "pass"
        }
    }

    fn grade_length(&self) -> &'static str {
        "pass"
    }

    /// fastqc `DuplicationLevel`: deduplicated percentage + 16 slots.
    /// No sampling correction when `countAtUniqueLimit == total` (fastqc
    /// `getCorrectedCount` early return), which holds below the 100k unique
    /// sequence cutoff.
    fn duplication_data(&self) -> (f64, Vec<f64>) {
        let mut collated: BTreeMap<u64, u64> = BTreeMap::new();
        for &c in self.seq_counts.values() {
            *collated.entry(c).or_insert(0) += 1;
        }
        let mut slots = vec![0.0f64; 16];
        let mut dedup_total = 0.0;
        let mut raw_total = 0.0;
        for (&level, &n_seq) in &collated {
            let corr = n_seq as f64;
            dedup_total += corr;
            raw_total += corr * level as f64;
            slots[dup_slot(level)] += corr * level as f64;
        }
        for v in slots.iter_mut() {
            *v = if raw_total > 0.0 { *v / raw_total * 100.0 } else { 0.0 };
        }
        let dedup_pct = if raw_total > 0.0 {
            dedup_total / raw_total * 100.0
        } else {
            100.0
        };
        (dedup_pct, slots)
    }

    fn grade_duplication(&self, dedup_pct: f64) -> &'static str {
        if dedup_pct < 50.0 {
            "fail"
        } else if dedup_pct < 70.0 {
            "warn"
        } else {
            "pass"
        }
    }

    fn overrep_rows(&self) -> Vec<(Vec<u8>, u64, f64, String)> {
        let total = self.n_reads.max(1) as f64;
        let mut rows: Vec<(Vec<u8>, u64, f64, String)> = Vec::new();
        for (&(len, seq), &count) in &self.seq_counts {
            let seq = &seq[..len as usize];
            let pct = count as f64 / total * 100.0;
            if pct > 0.1 {
                rows.push((seq.to_vec(), count, pct, match_source(seq)));
            }
        }
        rows.sort_by(|a, b| b.1.cmp(&a.1));
        rows
    }

    fn grade_overrep(&self, rows: &[(Vec<u8>, u64, f64, String)]) -> &'static str {
        if rows.iter().any(|r| r.2 > 1.0) {
            "fail"
        } else if rows.is_empty() {
            "pass"
        } else {
            "warn"
        }
    }

    /// Per-group enrichment for each adapter: mean of per-position
    /// percentages (fastqc `AdapterContent.calculateEnrichment`).
    fn adapter_enrichment(&self) -> Vec<[f64; 6]> {
        let total = self.n_reads.max(1) as f64;
        let groups = make_base_groups(self.adapter_positions.len() as u32);
        let mut out = vec![[0.0f64; 6]; groups.len()];
        for (g, grp) in groups.iter().enumerate() {
            let span = (grp.end - grp.start + 1) as f64;
            for p in (grp.start - 1)..grp.end.min(self.adapter_positions.len() as u32) {
                for a in 0..6 {
                    out[g][a] += self.adapter_positions[p as usize][a] as f64 * 100.0 / total;
                }
            }
            for a in 0..6 {
                out[g][a] /= span;
            }
        }
        out
    }

    fn grade_adapter(&self, groups: &[[f64; 6]]) -> &'static str {
        for g in groups {
            for &v in g.iter() {
                if v > 10.0 {
                    return "fail";
                }
            }
        }
        for g in groups {
            for &v in g.iter() {
                if v > 5.0 {
                    return "warn";
                }
            }
        }
        "pass"
    }

    /// fastqc `KmerContent`: 7-mer position counts, per-group obs/exp and a
    /// binomial test (Bonferroni-corrected ×4^7); returns the significant
    /// kmers sorted by p-value (top 6).
    fn kmer_content(&self) -> Vec<(u16, u64, f64, String)> {
        let total_kmer: u64 = self.total_kmer_per_pos.iter().sum();
        if total_kmer == 0 {
            return Vec::new();
        }
        let positions = self.total_kmer_per_pos.len();
        let groups = make_base_groups(positions as u32);
        let mut hits: Vec<(u16, u64, f64, String)> = Vec::new();
        let stride = self.kmer_stride;
        for key in 0..KMER_KEYS {
            let row = key * stride;
            // row sum over all sampled positions
            let mut kmer_total = 0u64;
            for p in 0..positions {
                kmer_total += self.kmer_grid[row + p] as u64;
            }
            if kmer_total == 0 {
                continue;
            }
            let expected_prop = kmer_total as f64 / total_kmer as f64;
            let mut best: Option<(f64, f64, String)> = None; // p, obs_exp, group label
            for g in &groups {
                let mut group_count = 0u64;
                let mut group_hits = 0u64;
                for p in (g.start - 1)..g.end.min(positions as u32) {
                    group_count += self.total_kmer_per_pos[p as usize];
                    group_hits += self.kmer_grid[row + p as usize] as u64;
                }
                if group_count == 0 {
                    continue;
                }
                let predicted = expected_prop * group_count as f64;
                if group_hits as f64 <= predicted {
                    continue;
                }
                let obs_exp = group_hits as f64 / predicted;
                let p = (1.0 - binomial_cdf(group_hits, group_count, expected_prop))
                    * 16384.0; // 4^7
                if p < 0.01 && obs_exp > 5.0 {
                    let better = match &best {
                        Some((bp, _bo, _bg)) => p < *bp,
                        None => true,
                    };
                    if better {
                        best = Some((p, obs_exp, g.to_string()));
                    }
                }
            }
            if let Some((p, oe, group)) = best {
                hits.push((key as u16, kmer_total, oe, group));
                let _ = p;
            }
        }
        hits.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap_or(std::cmp::Ordering::Equal));
        hits.truncate(6);
        hits
    }

    /// fastqc-style tile grade: maximum tile mean deviation from the global
    /// mean; > 10 fail, > 5 warn. Returns (grade, global mean).
    fn grade_tile(&self) -> (&'static str, f64) {
        let mut global = QualityCount::default();
        for quals in self.tile_quality.values() {
            for q in quals {
                global.merge(q);
            }
        }
        let global_mean = global.mean();
        let mut max_dev = 0.0f64;
        for quals in self.tile_quality.values() {
            let mut q = QualityCount::default();
            for p in quals {
                q.merge(p);
            }
            max_dev = max_dev.max((q.mean() - global_mean).abs());
        }
        let grade = if max_dev > 10.0 {
            "fail"
        } else if max_dev > 5.0 {
            "warn"
        } else {
            "pass"
        };
        (grade, global_mean)
    }
}

/// fastqc duplication level slot from the observed count.
fn dup_slot(level: u64) -> usize {
    let t = level as i64 - 1;
    if t > 9999 || t < 0 {
        15
    } else if t > 4999 {
        14
    } else if t > 999 {
        13
    } else if t > 499 {
        12
    } else if t > 99 {
        11
    } else if t > 49 {
        10
    } else if t > 9 {
        9
    } else {
        t as usize
    }
}

fn dup_label(slot: usize) -> &'static str {
    match slot {
        0..=8 => ["1", "2", "3", "4", "5", "6", "7", "8", "9"][slot],
        9 => ">10",
        10 => ">50",
        11 => ">100",
        12 => ">500",
        13 => ">1k",
        14 => ">5k",
        _ => ">10k+",
    }
}

/// fastqc `Configuration/adapter_list.txt` (6 entries, longest = 12 bp),
/// embedded from the shared `data/` directory. Data from FastQC (GPL v3) —
/// see `data/FastQC_DATA_LICENSE`.
static ADAPTERS: std::sync::OnceLock<Vec<(&'static str, &'static str)>> =
    std::sync::OnceLock::new();

fn adapters() -> &'static Vec<(&'static str, &'static str)> {
    ADAPTERS.get_or_init(|| {
        include_str!("../../../data/adapter_list.txt")
            .lines()
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .filter_map(|l| {
                let parts: Vec<&str> = l.split('\t').filter(|s| !s.is_empty()).collect();
                if parts.len() >= 2 {
                    Some((parts[0], parts[1]))
                } else {
                    None
                }
            })
            .collect()
    })
}

/// Unaligned fixed-width loads for the hot adapter scan. Callers guarantee
/// `i + N <= seq.len()` before invoking (bounds are checked in the loop).
#[inline(always)]
unsafe fn read_u16_un(seq: &[u8], i: usize) -> u16 {
    debug_assert!(i + 2 <= seq.len());
    unsafe { std::ptr::read_unaligned(seq.as_ptr().add(i) as *const u16) }
}

#[inline(always)]
unsafe fn read_u32_un(seq: &[u8], i: usize) -> u32 {
    debug_assert!(i + 4 <= seq.len());
    unsafe { std::ptr::read_unaligned(seq.as_ptr().add(i) as *const u32) }
}

#[inline(always)]
unsafe fn read_u64_un(seq: &[u8], i: usize) -> u64 {
    debug_assert!(i + 8 <= seq.len());
    unsafe { std::ptr::read_unaligned(seq.as_ptr().add(i) as *const u64) }
}

/// Precomputed fixed-width signatures for the 6 built-in adapters (all
/// 12 bp in the current data): `(bytes[2..4] as u16, bytes[4..12] as u64)`.
/// After the pairwise-distinct 2-byte dispatch, these verify the remaining
/// 10 bp with two loads + compares instead of a libc memcmp. Entries with
/// other lengths are `None` and fall back to the general slice comparison.
static ADAPTER_SIGS: std::sync::OnceLock<[Option<(u16, u64)>; 6]> =
    std::sync::OnceLock::new();

fn adapter_sigs() -> &'static [Option<(u16, u64)>; 6] {
    ADAPTER_SIGS.get_or_init(|| {
        let mut sigs = [None; 6];
        for (a, (_, s)) in adapters().iter().enumerate() {
            if a >= 6 {
                break;
            }
            let b = s.as_bytes();
            if b.len() == 12 {
                sigs[a] = Some((
                    u16::from_ne_bytes([b[2], b[3]]),
                    u64::from_ne_bytes(b[4..12].try_into().unwrap()),
                ));
            }
        }
        sigs
    })
}

/// 32768-entry 4-bit lookup: pair (little-endian u16 of b0,b1) -> adapter
/// index (0xF = none), packed two per byte. 32 KB fits L1 on modern x86
/// (the 64 KB u8 version sat at the L1/L2 boundary), and the sentinel
/// encoding keeps the whole lookup to one load + shift/mask.
static ADAPTER_PAIR_TABLE: std::sync::OnceLock<[u8; 1 << 15]> =
    std::sync::OnceLock::new();

fn adapter_pair_table() -> &'static [u8; 1 << 15] {
    ADAPTER_PAIR_TABLE.get_or_init(|| {
        let mut t = [0xFFu8; 1 << 15];
        for (a, (_, s)) in adapters().iter().enumerate() {
            if a >= 6 {
                break;
            }
            let b = s.as_bytes();
            if b.len() >= 2 {
                let pair = b[0] as usize | (b[1] as usize) << 8;
                let nib = a as u8 & 0xF;
                if pair & 1 == 0 {
                    t[pair >> 1] = (t[pair >> 1] & 0xF0) | nib;
                } else {
                    t[pair >> 1] = (t[pair >> 1] & 0x0F) | (nib << 4);
                }
            }
        }
        t
    })
}

/// Look up the adapter candidate for a 2-byte pair in the packed table
/// (0xF = no adapter starts with this pair).
#[inline(always)]
fn adapter_at(table: &[u8; 1 << 15], pair: u16) -> u8 {
    let idx = pair as usize;
    let v = table[idx >> 1];
    if idx & 1 == 0 {
        v & 0xF
    } else {
        v >> 4
    }
}

/// Verify that adapter `a` really occurs at `i` (fixed-width 12 bp compare
/// via the precomputed signature; generic slice fallback for other lengths).
#[inline(always)]
fn verify_adapter(
    a: usize,
    i: usize,
    l: usize,
    seq: &[u8],
    sigs: &[Option<(u16, u64)>; 6],
    ads: &[(&str, &str)],
    found: &mut [Option<usize>; 6],
    mask: &mut u8,
) {
    if *mask & (1 << a) == 0 {
        return;
    }
    // `a` comes from the pair table (0..5 or the 0xFF sentinel, already
    // filtered), so all three indexes are in bounds; unchecked indexing
    // removes two bounds compares on the ~18 false-positive pairs per read.
    debug_assert!(a < 6);
    let matched = if let Some((tail2, body8)) = *unsafe { sigs.get_unchecked(a) } {
        i + 12 <= l
            && unsafe { read_u16_un(seq, i + 2) } == tail2
            && unsafe { read_u64_un(seq, i + 4) } == body8
    } else {
        let aseq = unsafe { ads.get_unchecked(a) }.1.as_bytes();
        i + aseq.len() <= l && &seq[i..i + aseq.len()] == aseq
    };
    if matched {
        unsafe { *found.get_unchecked_mut(a) = Some(i) };
        *mask &= !(1 << a);
    }
}

/// Full fastqc contaminant library (`data/contaminant_list.txt`, 151
/// entries), embedded at compile time via `include_str!`.
static CONTAMINANTS: std::sync::OnceLock<Vec<(&'static str, &'static str)>> =
    std::sync::OnceLock::new();

fn contaminants() -> &'static Vec<(&'static str, &'static str)> {
    CONTAMINANTS.get_or_init(|| {
        include_str!("../../../data/contaminant_list.txt")
            .lines()
            .filter(|l| !l.is_empty() && !l.starts_with('#'))
            .filter_map(|l| {
                let parts: Vec<&str> = l.split('\t').filter(|s| !s.is_empty()).collect();
                if parts.len() >= 2 {
                    Some((parts[0], parts[1]))
                } else {
                    None
                }
            })
            .collect()
    })
}

/// fastqc `Contaminant.findMatch` across the full library: allow one
/// mismatch, require a run > 20 bp for queries >= 20 bp; exact substring for
/// 8-20 bp queries. Returns the best hit (name + percent-id + length).
fn match_source(seq: &[u8]) -> String {
    let mut best: Option<(usize, usize, String)> = None; // len, id, name
    for (name, cseq) in contaminants() {
        let subjects = [cseq.as_bytes().to_vec(), revcomp(cseq.as_bytes())];
        for subject in &subjects {
            if let Some((len, id)) = contaminant_match(subject, seq) {
                let better = match &best {
                    Some((bl, bi, _bn)) => len > *bl || (len == *bl && id > *bi),
                    None => true,
                };
                if better {
                    best = Some((len, id, name.to_string()));
                }
            }
        }
    }
    match best {
        Some((len, id, name)) => format!("{name} ({id}% over {len}bp)"),
        None => "No Hit".to_string(),
    }
}

/// FastQC `Contaminant.findMatch(char[], char[], offset, direction)` over all
/// offsets: the longest >20 bp run allowing a single mismatch; percent
/// identity = (len - mismatches) * 100 / len.
fn contaminant_match(subject: &[u8], query: &[u8]) -> Option<(usize, usize)> {
    if query.is_empty() {
        return None;
    }
    if query.len() < 20 && query.len() >= 8 {
        if contains(subject, query) {
            return Some((query.len(), 100));
        }
        return None;
    }
    let mut best: Option<(usize, usize)> = None;
    let lo = 0i64 - (subject.len() as i64 - 20);
    let hi = query.len() as i64 - 20;
    for offset in lo..hi {
        let mut mismatch_count = 0usize;
        let mut start = 0usize;
        let mut end = 0usize;
        for (i, &sc) in subject.iter().enumerate() {
            let qpos = i as i64 + offset;
            if qpos < 0 {
                start = i + 1;
                continue;
            }
            if qpos >= query.len() as i64 {
                break;
            }
            if sc == query[qpos as usize] {
                end = i;
            } else {
                mismatch_count += 1;
                if mismatch_count > 1 {
                    let len = 1 + end as i64 - start as i64;
                    if len > 20 {
                        let id = ((len - (mismatch_count as i64 - 1)) * 100) / len;
                        let better = match best {
                            Some((bl, bi)) => len as usize > bl || (len as usize == bl && id as usize > bi),
                            None => true,
                        };
                        if better {
                            best = Some((len as usize, id as usize));
                        }
                    }
                    start = i + 1;
                    end = i + 1;
                }
            }
        }
        let len = 1 + end as i64 - start as i64;
        if len > 20 {
            let id = ((len - mismatch_count.min(1) as i64) * 100) / len;
            let better = match best {
                Some((bl, bi)) => len as usize > bl || (len as usize == bl && id as usize > bi),
                None => true,
            };
            if better {
                best = Some((len as usize, id as usize));
            }
        }
    }
    best
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'N' => b'N',
            other => other,
        })
        .collect()
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    find_subseq(haystack, needle).is_some()
}

fn find_subseq(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.len() < 2 || needle.len() > haystack.len() {
        return None;
    }
    // Two-byte prescreen: scan for the first byte, check the second byte
    // (plain compare), then verify the rest — avoids a memcmp on every
    // first-byte hit (~1/4 of positions) and cuts memcmp calls ~16× for
    // adapter matching (6× per read).
    let first = needle[0];
    let second = needle[1];
    let rest = &needle[2..];
    let mut i = 0;
    while i + needle.len() <= haystack.len() {
        let Some(off) = haystack[i..].iter().position(|&b| b == first) else {
            return None;
        };
        i += off;
        if i + 2 + rest.len() > haystack.len() {
            return None;
        }
        if haystack[i + 1] == second && &haystack[i + 2..i + 2 + rest.len()] == rest {
            return Some(i);
        }
        i += 1;
    }
    None
}

/// 2-bit base encoding for 7-mer keys (matches fastqc's kmer key order:
/// A=0, C=1, G=2, anything else = 3; the rolling window skips non-ACGT
/// windows before use).
#[inline(always)]
fn encode_base(b: u8) -> u16 {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        _ => 3,
    }
}

/// Parse an ASCII decimal u32 without UTF-8 validation or allocation
/// (tile ids in the read header).
#[inline]
fn parse_u32_bytes(b: &[u8]) -> Option<u32> {
    if b.is_empty() {
        return None;
    }
    let mut v = 0u32;
    for &c in b {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v.checked_mul(10)? + (c - b'0') as u32;
    }
    Some(v)
}

fn kmer_to_string(key: u16) -> String {
    const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
    let mut out = String::with_capacity(7);
    for shift in (0..7).rev() {
        out.push(BASES[((key >> (shift * 2)) & 0b11) as usize]);
    }
    out
}

/// P(X <= k) for X ~ Binomial(n, p) via the regularized incomplete beta
/// function (Numerical Recipes `betai`/`betacf`).
fn binomial_cdf(k: u64, n: u64, p: f64) -> f64 {
    if k >= n {
        return 1.0;
    }
    if p <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return 1.0;
    }
    let a = (n - k) as f64;
    let b = (k + 1) as f64;
    let x = 1.0 - p;
    let mut bt = (gammaln(a + b) - gammaln(a) - gammaln(b)
        + a * x.ln()
        + b * (1.0 - x).ln())
    .exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        bt * betacf(x, a, b) / a
    } else {
        1.0 - bt * betacf(1.0 - x, b, a) / b
    }
}

fn gammaln(xx: f64) -> f64 {
    let cof = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5,
    ];
    let mut x = xx;
    let mut y = xx;
    let mut tmp = x + 5.5;
    tmp -= (x + 0.5) * tmp.ln();
    let mut ser = 1.000000000190015;
    for (i, &c) in cof.iter().enumerate() {
        y += 1.0;
        ser += c / y;
    }
    -tmp + (2.5066282746310005 * ser / x).ln()
}

fn normal_pdf(x: f64, mean: f64, stdev: f64) -> f64 {
    if stdev <= 0.0 {
        return 0.0;
    }
    let z = (x - mean) / stdev;
    (-0.5 * z * z).exp() / (stdev * (2.0 * std::f64::consts::PI).sqrt())
}

fn betacf(x: f64, a: f64, b: f64) -> f64 {
    const MAXIT: usize = 200;
    const EPS: f64 = 3.0e-7;
    const FPMIN: f64 = 1.0e-30;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < FPMIN {
        d = FPMIN;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=MAXIT {
        let m2 = 2 * m;
        let mut aa = m as f64 * (b - m as f64) * x / ((qam + m2 as f64) * (a + m2 as f64));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m as f64) * (qab + m as f64) * x
            / ((a + m2 as f64) * (qap + m2 as f64));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < EPS {
            break;
        }
    }
    h
}

/// Minimal FNV-1a hasher: much faster than `SipHash` for the fixed 51-byte
/// `seq_counts` keys (hashbrown only needs any valid `Hasher`).
#[derive(Default)]
struct FnvHasher(u64);

impl std::hash::Hasher for FnvHasher {
    fn write(&mut self, bytes: &[u8]) {
        let mut h = self.0;
        // Hash the first 8 bytes and the last 8 bytes (2 multiplies instead
        // of 7 for a 51-byte key). The full key is still stored and
        // compared on probe, so correctness is unaffected; the suffix fold
        // keeps adapter-prefixed read families from clustering in one
        // bucket (they share the prefix but not the tail).
        if bytes.len() >= 8 {
            h ^= u64::from_ne_bytes(bytes[..8].try_into().unwrap());
            h = h.wrapping_mul(0x100000001b3);
        } else {
            for &b in bytes {
                h ^= b as u64;
                h = h.wrapping_mul(0x100000001b3);
            }
        }
        // fold the tail so different 8-byte prefixes plus distinct suffixes
        // still get distinct hashes cheaply
        if bytes.len() > 8 {
            let tail = &bytes[bytes.len() - 8..];
            h ^= u64::from_ne_bytes(tail.try_into().unwrap());
            h = h.wrapping_mul(0x100000001b3);
        }
        self.0 = h;
    }

    fn finish(&self) -> u64 {
        self.0
    }
}

#[derive(Default)]
struct FnvBuild;

impl std::hash::BuildHasher for FnvBuild {
    type Hasher = FnvHasher;
    fn build_hasher(&self) -> FnvHasher {
        FnvHasher(0xcbf29ce484222325)
    }
}

/// fastqc `BasicStats.formatLength`: bp/kbp/Mbp/Gbp, trailing `.0` dropped.
fn format_length(len: u64) -> String {
    let (value, unit) = if len >= 1_000_000_000 {
        (len as f64 / 1_000_000_000.0, " Gbp")
    } else if len >= 1_000_000 {
        (len as f64 / 1_000_000.0, " Mbp")
    } else if len >= 1_000 {
        (len as f64 / 1_000.0, " kbp")
    } else {
        (len as f64, " bp")
    };
    let raw = format!("{:.1}", value);
    let raw = raw.trim_end_matches(".0");
    format!("{raw}{unit}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_length_units() {
        assert_eq!(format_length(108), "108 bp");
        assert_eq!(format_length(216_000), "216 kbp");
        assert_eq!(format_length(2_500_000), "2.5 Mbp");
        assert_eq!(format_length(3_000_000_000), "3 Gbp");
    }

    #[test]
    fn quality_count_percentile_matches_fastqc() {
        // 4 reads: Phred 0, 2, 4, 7
        let mut q = QualityCount::default();
        for c in [0u8, 2, 4, 7] {
            q.add(c);
        }
        assert_eq!(q.percentile(50), 2.0);
        assert_eq!(q.percentile(25), 0.0);
        assert_eq!(q.percentile(90), 4.0);
        assert_eq!(q.mean(), 3.25);
    }
}

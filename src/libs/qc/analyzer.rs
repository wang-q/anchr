//! M1 QC analyzer: statistical modules + FastQC-format text report.
//!
//! Numeric behavior follows fastqc 0.12.1 (verified against golden output
//! in `tests/qc/golden/`): integer average quality per read, `QualityCount`
//! percentile/mean over ASCII quality chars, GCModel-fraction GC histogram,
//! integer-division `%GC` in Basic Statistics.

use pgr::libs::fmt::seq::SeqRecord;
use std::collections::BTreeMap;
use std::io::Write;

use super::base_groups::{make_base_groups, BaseGroup};

const QUAL_BINS: usize = 128;
const GC_BINS: usize = 101;

/// Per-position ASCII quality counts (index = ASCII char value).
#[derive(Clone, Copy)]
struct QualityCount {
    counts: [u64; QUAL_BINS],
}

impl Default for QualityCount {
    fn default() -> Self {
        QualityCount {
            counts: [0; QUAL_BINS],
        }
    }
}

impl QualityCount {
    fn add(&mut self, qc: u8) {
        self.counts[qc as usize] += 1;
    }

    fn total(&self) -> u64 {
        self.counts.iter().sum()
    }

    fn mean(&self, offset: u8) -> f64 {
        let mut total = 0u64;
        let mut n = 0u64;
        for (i, &c) in self.counts.iter().enumerate().skip(offset as usize) {
            total += (i as u64 - offset as u64) * c;
            n += c;
        }
        if n == 0 {
            0.0
        } else {
            total as f64 / n as f64
        }
    }

    /// fastqc `QualityCount.getPercentile`: integer `n*p/100` then the first
    /// quality whose cumulative count reaches the target (Phred value).
    fn percentile(&self, offset: u8, p: u64) -> f64 {
        let n = self.total();
        if n == 0 {
            return -1.0;
        }
        let target = n * p / 100;
        let mut count = 0u64;
        for (i, &c) in self.counts.iter().enumerate().skip(offset as usize) {
            count += c;
            if count >= target {
                return (i - offset as usize) as f64;
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
    g_count: u64,
    c_count: u64,
    a_count: u64,
    t_count: u64,
    len_min: u32,
    len_max: u32,
    len_hist: BTreeMap<u32, u64>,
    encoding_offset: u8,
    // per-position accumulators (index = 0-based position)
    per_base_quality: Vec<QualityCount>,
    per_base_content: Vec<[u64; 5]>, // A C G T N
    // per-read accumulators
    seq_quality_hist: BTreeMap<i32, u64>, // integer avg ASCII quality -> reads
    seq_gc_hist: Vec<f64>,
    gc_models: BTreeMap<u32, GcModel>,
    // M2: truncated-sequence counts (fastqc OverRepresentedSeqs; also feeds
    // the DuplicationLevel module)
    seq_counts: BTreeMap<Vec<u8>, u64>,
    // M2: adapter position counts (fastqc AdapterContent)
    adapter_positions: Vec<[u64; 6]>,
    longest_seq: u32,
    // M2: 7-mer position counts (fastqc KmerContent; sense strand only)
    kmer_counts: std::collections::HashMap<u16, Vec<u64>>,
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
            ..Default::default()
        }
    }

    /// Consume one record (streaming; position vectors grow to max length).
    /// Consume one record; `global_index` is the 0-based index of the read
    /// in the whole file (used for the 2% kmer sampling, which must be
    /// global even when processing is split across parallel chunks).
    pub fn consume(&mut self, rec: &SeqRecord, global_index: u64) {
        let seq = rec.sequence();
        let qual = rec.quality_scores();
        let len = seq.len() as u32;

        self.n_reads += 1;
        self.total_bases += len as u64;
        *self.len_hist.entry(len).or_insert(0) += 1;
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
            self.per_base_content.resize(seq.len(), [0u64; 5]);
        }
        for (i, &b) in seq.iter().enumerate() {
            if i < qual.len() {
                self.per_base_quality[i].add(qual[i]);
            }
            let idx = match b {
                b'A' => Some(0),
                b'C' => Some(1),
                b'G' => Some(2),
                b'T' => Some(3),
                b'N' => Some(4),
                _ => None,
            };
            if let Some(k) = idx {
                self.per_base_content[i][k] += 1;
            }
            match b {
                b'G' => self.g_count += 1,
                b'C' => self.c_count += 1,
                b'A' => self.a_count += 1,
                b'T' => self.t_count += 1,
                _ => {}
            }
        }

        // per-read average quality: integer division over ASCII chars
        if !qual.is_empty() {
            let sum: u64 = qual.iter().map(|&c| c as u64).sum();
            let avg = (sum / qual.len() as u64) as i32;
            *self.seq_quality_hist.entry(avg).or_insert(0) += 1;
        }

        // per-read GC: truncate >100bp to hundreds, then GCModel fractions
        let trunc = if seq.len() > 1000 {
            (seq.len() / 1000) * 1000
        } else if seq.len() > 100 {
            (seq.len() / 100) * 100
        } else {
            seq.len()
        };
        let gc_count = seq[..trunc].iter().filter(|&&b| b == b'G' || b == b'C').count() as u32;
        let model = self
            .gc_models
            .entry(trunc as u32)
            .or_insert_with(|| GcModel::new(trunc as u32));
        for &(p, inc) in model.values(gc_count) {
            self.seq_gc_hist[p] += inc;
        }

        // overrep/duplication: truncate to 50 bp (fastqc default)
        let trunc = if seq.len() > 50 { &seq[..50] } else { seq };
        *self.seq_counts.entry(trunc.to_vec()).or_insert(0) += 1;

        // adapter content: indexOf each adapter, increment to the current
        // longest-read bound (fastqc semantics)
        if seq.len() as u32 > self.longest_seq {
            self.longest_seq = seq.len() as u32;
            let cur_max = self.longest_seq.saturating_sub(12) + 1;
            self.adapter_positions.resize(cur_max as usize, [0u64; 6]);
        }
        let cur_max = self.longest_seq.saturating_sub(12) + 1;
        for (a, (_, aseq)) in adapters().iter().enumerate() {
            if let Some(idx) = find_subseq(seq, aseq.as_bytes()) {
                for p in idx..cur_max as usize {
                    self.adapter_positions[p][a] += 1;
                }
            }
        }

        // 7-mer position counts: fastqc samples every 50th read (2%) and
        // skips kmers containing Ns
        if (global_index + 1) % 50 == 0 && seq.len() >= 7 {
            let seq = if seq.len() > 500 { &seq[..500] } else { seq };
            if self.total_kmer_per_pos.len() < seq.len() - 6 {
                self.total_kmer_per_pos.resize(seq.len() - 6, 0);
            }
            for i in 0..=(seq.len() - 7) {
                if seq[i..i + 7].iter().any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T')) {
                    continue;
                }
                let key = encode_kmer(&seq[i..i + 7]);
                let entry = self.kmer_counts.entry(key).or_default();
                if entry.len() < self.total_kmer_per_pos.len() {
                    entry.resize(self.total_kmer_per_pos.len(), 0);
                }
                entry[i] += 1;
                self.total_kmer_per_pos[i] += 1;
            }
        }

        // per-tile quality: fastqc samples all of the first 10k reads, then
        // every 10th; tile comes from the Casava 1.8+ (field 4) or legacy
        // (field 2) header
        if !self.tile_ignore && !qual.is_empty() {
            self.tile_total_count += 1;
            if self.tile_total_count > 10_000 && self.tile_total_count % 10 != 0 {
                // skip sampling
            } else if let Some(tile) = self.detect_tile(rec.name()) {
                let entry = self.tile_quality.entry(tile).or_default();
                if entry.len() < qual.len() {
                    entry.resize(qual.len(), QualityCount::default());
                }
                for (i, &q) in qual.iter().enumerate() {
                    entry[i].add(q);
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
        let s = std::str::from_utf8(name).ok()?;
        let fields: Vec<&str> = s.split(':').collect();
        if let Some(pos) = self.tile_split_pos {
            if pos < fields.len() {
                return fields[pos].parse().ok();
            }
            self.tile_ignore = true;
            return None;
        }
        if fields.len() >= 7 {
            self.tile_split_pos = Some(4);
            return fields[4].parse().ok();
        }
        if fields.len() >= 5 {
            self.tile_split_pos = Some(2);
            return fields[2].parse().ok();
        }
        self.tile_ignore = true;
        None
    }

    /// Merge another statistics block into `self` (parallel chunk results).
    pub fn merge(&mut self, other: &QcStats) {
        let was_empty = self.n_reads == 0;
        self.n_reads += other.n_reads;
        self.total_bases += other.total_bases;
        self.g_count += other.g_count;
        self.c_count += other.c_count;
        self.a_count += other.a_count;
        self.t_count += other.t_count;
        if other.n_reads > 0 {
            if was_empty {
                self.len_min = other.len_min;
                self.len_max = other.len_max;
            } else {
                self.len_min = self.len_min.min(other.len_min);
                self.len_max = self.len_max.max(other.len_max);
            }
        }
        for (&len, &count) in &other.len_hist {
            *self.len_hist.entry(len).or_insert(0) += count;
        }
        if other.per_base_quality.len() > self.per_base_quality.len() {
            self.per_base_quality
                .resize(other.per_base_quality.len(), QualityCount::default());
            self.per_base_content
                .resize(other.per_base_content.len(), [0u64; 5]);
        }
        for (a, b) in self.per_base_quality.iter_mut().zip(&other.per_base_quality) {
            a.merge(b);
        }
        for (a, b) in self.per_base_content.iter_mut().zip(&other.per_base_content) {
            for k in 0..5 {
                a[k] += b[k];
            }
        }
        for (&q, &c) in &other.seq_quality_hist {
            *self.seq_quality_hist.entry(q).or_insert(0) += c;
        }
        for (a, &b) in self.seq_gc_hist.iter_mut().zip(&other.seq_gc_hist) {
            *a += b;
        }
        for (seq, &c) in &other.seq_counts {
            *self.seq_counts.entry(seq.clone()).or_insert(0) += c;
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
            self.total_kmer_per_pos.resize(other.total_kmer_per_pos.len(), 0);
        }
        for (a, &b) in self.total_kmer_per_pos.iter_mut().zip(&other.total_kmer_per_pos) {
            *a += b;
        }
        for (&key, pos) in &other.kmer_counts {
            let entry = self.kmer_counts.entry(key).or_default();
            if entry.len() < pos.len() {
                entry.resize(pos.len(), 0);
            }
            for (a, &b) in entry.iter_mut().zip(pos) {
                *a += b;
            }
        }
    }

    /// Per-group quality stats: each position's percentile/mean is computed
    /// independently (positions with >100 quality chars), then averaged
    /// across the group (fastqc `PerBaseQualityScores.getPercentile/getMean`).
    fn group_quality_stats(&self) -> Vec<GroupQualityStats> {
        let groups = make_base_groups(self.per_base_quality.len() as u32);
        let offset = self.encoding_offset;
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
                    mean += q.mean(offset);
                    median += q.percentile(offset, 50);
                    lq += q.percentile(offset, 25);
                    uq += q.percentile(offset, 75);
                    p10 += q.percentile(offset, 10);
                    p90 += q.percentile(offset, 90);
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
    fn group_contents(&self) -> Vec<[u64; 5]> {
        let groups = make_base_groups(self.per_base_content.len() as u32);
        let mut out = vec![[0u64; 5]; groups.len()];
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
        self.a_count + self.c_count + self.g_count + self.t_count
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
        let at = self.a_count + self.t_count;
        let gc = self.g_count + self.c_count;
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
        if let (Some(&min), Some(&max)) =
            (self.seq_quality_hist.keys().min(), self.seq_quality_hist.keys().max())
        {
            for q in min..=max {
                let count = self.seq_quality_hist.get(&q).copied().unwrap_or(0);
                writeln!(w, "{}\t{:.1}", q - offset as i32, count as f64)?;
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
        for (&len, &count) in &self.len_hist {
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
                        q.mean(self.encoding_offset) - global_mean
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
        let at = self.a_count + self.t_count;
        let gc = self.g_count + self.c_count;
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

        let psq: Vec<serde_json::Value> = self
            .seq_quality_hist
            .iter()
            .map(|(&q, &c)| json!([q - offset as i32, c as f64]))
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
            .map(|(&l, &c)| json!([l, c as f64]))
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
        for (&q, &c) in &self.seq_quality_hist {
            total += q as u64 * c;
            n += c;
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

    fn grade_per_base_content(&self, contents: &[[u64; 5]]) -> &'static str {
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

    fn grade_per_base_n(&self, contents: &[[u64; 5]]) -> &'static str {
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
        for (seq, &count) in &self.seq_counts {
            let pct = count as f64 / total * 100.0;
            if pct > 0.1 {
                rows.push((seq.clone(), count, pct, match_source(seq)));
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
        for (&key, pos) in &self.kmer_counts {
            let kmer_total: u64 = pos.iter().sum();
            let expected_prop = kmer_total as f64 / total_kmer as f64;
            let mut best: Option<(f64, f64, String)> = None; // p, obs_exp, group label
            for g in &groups {
                let mut group_count = 0u64;
                let mut group_hits = 0u64;
                for p in (g.start - 1)..g.end.min(positions as u32) {
                    group_count += self.total_kmer_per_pos[p as usize];
                    group_hits += pos[p as usize];
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
                hits.push((key, kmer_total, oe, group));
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
        let global_mean = global.mean(self.encoding_offset);
        let mut max_dev = 0.0f64;
        for quals in self.tile_quality.values() {
            let mut q = QualityCount::default();
            for p in quals {
                q.merge(p);
            }
            max_dev = max_dev.max((q.mean(self.encoding_offset) - global_mean).abs());
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
    if needle.is_empty() || needle.len() > haystack.len() {
        return None;
    }
    haystack
        .windows(needle.len())
        .position(|w| w == needle)
}

fn encode_kmer(seq: &[u8]) -> u16 {
    let mut key = 0u16;
    for &b in seq {
        key <<= 2;
        key |= match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            _ => 3,
        };
    }
    key
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
        // 4 reads: Phred 0, 2, 4, 7 (ASCII 33/35/37/40)
        let mut q = QualityCount::default();
        for c in [33u8, 35, 37, 40] {
            q.add(c);
        }
        assert_eq!(q.percentile(33, 50), 2.0);
        assert_eq!(q.percentile(33, 25), 0.0);
        assert_eq!(q.percentile(33, 90), 4.0);
        assert_eq!(q.mean(33), 3.25);
    }
}

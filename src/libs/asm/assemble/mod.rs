//! Tadpole-compatible assembly: shared options, stats, and record
//! reading. The tadpole contig mode (greedy extension + BubblePopper)
//! lives in `contig`/`bubble`; the BCALM-style unitig mode lives in
//! `unitig`.

mod bubble;
mod contig;
mod unitig;

pub(crate) use super::is_junction;
pub use contig::assemble;
pub use unitig::{assemble_unitigs, assemble_unitigs_buf};
pub(crate) use unitig::{
    assemble_unitigs_core, assemble_unitigs_from_table, compute_links, Link, Unitig,
};

use crate::libs::asm::table::canonicalize_quality;
use anyhow::Result;
use pgr::libs::fmt::seq::{SeqReader, SeqRecord};
use pgr::libs::fq::qual::to_phred;

/// Assembly options with tadpole.sh defaults.
#[derive(Debug, Clone)]
pub struct AssembleOptions {
    /// K-mer length.
    pub k: usize,
    /// Ignore k-mers below this error-free probability.
    pub min_prob: f32,
    /// Minimum k-mer depth to seed a contig.
    pub min_count_seed: usize,
    /// Minimum k-mer depth to continue an extension.
    pub min_count_extend: usize,
    /// Minimum added bases past the seed for a contig to be kept.
    pub min_extension: usize,
    /// Minimum output length; `0` keeps everything (the `asm unitig`
    /// default mirrors bcalm's lossless vertex decomposition — no length
    /// cleanup unless requested; `asm contig` passes tadpole's
    /// `mincontiglen` default explicitly).
    pub min_contig_len: usize,
    /// Minimum k-mer coverage for a contig.
    pub min_coverage: f32,
    /// Branch ratio at high depth (branchmult1).
    pub branch_mult1: f32,
    /// Branch ratio at low depth (branchmult2).
    pub branch_mult2: f32,
    /// Second-highest depth considered "low" (branchlower).
    pub branch_lower_const: usize,
    /// Number of seeding passes (contigPasses).
    pub contig_passes: usize,
    /// Seeding pass multiplier (contigPassMult).
    pub contig_pass_mult: f64,
    /// Merge parallel paths in the contig graph (Tadpole popbubbles).
    pub pop_bubbles: bool,
    /// Append `L:` links to unitig FASTA headers (BCALM format).
    pub emit_links: bool,
    /// Emit a GFA graph instead of FASTA.
    pub emit_gfa: bool,
    /// Emit every k-mer abundance in the FASTA header (`ab:Z:`, BCALM
    /// `-all-abundance-counts`).
    pub all_abundance_counts: bool,
    /// Experimental: classify vertices once (DFA state) and walk unitigs
    /// from the state table instead of re-scanning extension buckets.
    pub use_dfa: bool,
    /// Experimental: use pgr's FastK-style super-mer two-stage counter
    /// instead of the direct emission + sort path (no quality gating).
    pub use_supermer: bool,
    /// Minimizer length for `--supermer` (None = pgr default).
    pub supermer_m: Option<usize>,
    /// Worker threads for the whole k-mer pipeline (counting + DFA
    /// classification); `0` uses the rayon global pool (all cores). The
    /// walk stays deterministic single-threaded.
    pub parallel: usize,
}

impl Default for AssembleOptions {
    fn default() -> Self {
        Self {
            k: 31,
            min_prob: 0.5,
            min_count_seed: 3,
            min_count_extend: 2,
            min_extension: 2,
            min_contig_len: 0,
            min_coverage: 1.0,
            branch_mult1: 20.0,
            branch_mult2: 3.0,
            branch_lower_const: 3,
            contig_passes: 16,
            contig_pass_mult: 1.7,
            pop_bubbles: true,
            emit_links: false,
            emit_gfa: false,
            all_abundance_counts: false,
            use_dfa: false,
            use_supermer: false,
            supermer_m: None,
            parallel: 0,
        }
    }
}

impl AssembleOptions {
    fn resolved_min_contig_len(&self) -> usize {
        self.min_contig_len
    }
}

/// Assembly statistics.
#[derive(Debug, Default, Clone)]
pub struct AssembleStats {
    pub reads_in: u64,
    pub contigs_built: u64,
    pub bases_built: u64,
    pub longest_contig: usize,
}

/// Reads all records from any number of files sequentially, canonicalizing
/// qualities like BBTools (shared by the contig and unitig modes).
/// Reads and converts records to `(sequence, phred quality)` pairs in one
/// streaming pass (one `SeqRecord` buffer alive at a time). The previous
/// version collected every `SeqRecord` first, then copied seq + phred into
/// a second structure — the two full copies together were ~0.6 GB on G37
/// full. Pairing is irrelevant for assembly (BCALM semantics): every record
/// from every file contributes its k-mers in order.
pub(crate) fn read_records(infiles: &[String]) -> Result<Vec<(Vec<u8>, Vec<u8>)>> {
    anyhow::ensure!(!infiles.is_empty(), "at least one input file is required");
    let mut reads: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    let mut rec = SeqRecord::new();
    for infile in infiles {
        let mut reader = SeqReader::new(infile)?;
        while reader.read_record(&mut rec)? {
            canonicalize_quality(&mut rec);
            reads.push((
                rec.sequence().to_vec(),
                to_phred(rec.sequence(), rec.quality_scores()),
            ));
        }
    }
    Ok(reads)
}

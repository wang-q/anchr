//! FASTQ read-processing algorithms.
//!
//! Business logic migrated from pgr; Phred quality handling stays in pgr
//! (`pgr::libs::fq::qual`, `pgr::libs::fq::pairs`). The general-purpose
//! FASTA/FASTQ reader also stays in pgr (`pgr::libs::fmt::seq`); `scan`
//! here is the zero-copy FASTQ scanner used by the QC hot path.

pub mod bbnet;
pub mod clump;
pub mod merge;
pub mod norm;
pub mod overlap;
pub mod sample;
pub mod scan;
pub mod split;
pub mod trim;
pub mod trim_adapter;

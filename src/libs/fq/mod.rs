//! FASTQ read-processing algorithms.
//!
//! Business logic migrated from pgr; the FASTA/FASTQ readers and Phred
//! quality handling stay in pgr (`pgr::libs::fmt`, `pgr::libs::fq::qual`,
//! `pgr::libs::fq::pairs`).

pub mod bbnet;
pub mod clump;
pub mod merge;
pub mod norm;
pub mod overlap;
pub mod sample;
pub mod split;
pub mod trim;
pub mod trim_adapter;

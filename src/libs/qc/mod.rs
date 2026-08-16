//! Reads quality-control statistics (`anchr fq qc`).
//!
//! Design: `notes/design/fq-qc.md` (M1: statistical modules). Output format
//! follows FastQC 0.12.1 `fastqc_data.txt` / `summary.txt`; the base
//! grouping replicates the installed fastqc 0.12.1 `BaseGroup` behavior
//! (verified against golden output).

pub mod analyzer;
pub mod base_groups;

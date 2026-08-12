//! Format I/O migrated from pgr.
//!
//! FASTA/FASTQ readers stay in pgr (`pgr::libs::fmt::seq/fa/fq`); anchr only
//! hosts the SAM conversion used by `anchr sam` (ihist/to-rg).

pub mod sam;

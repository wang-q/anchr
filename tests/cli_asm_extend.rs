#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::AnchrCmd;
use std::fs;

/// Extends a contig through the reads' overhang: the contig is the first
/// half of the genome and the reads cover the full genome, so the walk
/// recovers the missing second half.
#[test]
fn command_asm_extend_recovers_overhang() {
    let dir = tempfile::tempdir().unwrap();
    let reads = dir.path().join("reads.fa");
    let contigs = dir.path().join("contigs.fa");
    let out = dir.path().join("out.fa");

    // Deterministic random genome of 300 bp; reads are 10 full-length copies.
    let mut rng = 2026u64;
    let mut genome = Vec::new();
    for _ in 0..300 {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        genome.push(b"ACGT"[(rng >> 33) as usize % 4]);
    }
    let mut fa = String::new();
    for i in 0..10 {
        fa.push_str(&format!(">r{i}\n"));
        fa.push_str(&String::from_utf8(genome.clone()).unwrap());
        fa.push('\n');
    }
    fs::write(&reads, fa).unwrap();
    // Contig = first 170 bp (a 31-mer seed near the end is covered by reads
    // that extend beyond it).
    let mut cf = String::new();
    cf.push_str(">c1\n");
    cf.push_str(std::str::from_utf8(&genome[..170]).unwrap());
    cf.push('\n');
    fs::write(&contigs, cf).unwrap();

    let (_, stderr) = AnchrCmd::new()
        .args(&[
            "asm",
            "extend",
            contigs.to_str().unwrap(),
            reads.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-len",
            "0",
            "--max-extend",
            "200",
        ])
        .run();
    assert_eq!(stderr, "");
    let data = fs::read_to_string(&out).unwrap();
    let seq: String = data.lines().filter(|l| !l.starts_with('>')).collect();
    // The extension must recover most of the missing 130 bp.
    assert!(
        seq.len() >= 250,
        "extension too short: {} bp (expected >= 250)",
        seq.len()
    );
}

/// `--min-len` guards against extending short fragments (repeat contexts
/// would join copies of the same element): shorter contigs pass through
/// unchanged.
#[test]
fn command_asm_extend_min_len_skips_short_contigs() {
    let dir = tempfile::tempdir().unwrap();
    let reads = dir.path().join("reads.fa");
    let contigs = dir.path().join("contigs.fa");
    let out = dir.path().join("out.fa");

    let mut rng = 2026u64;
    let mut genome = Vec::new();
    for _ in 0..300 {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        genome.push(b"ACGT"[(rng >> 33) as usize % 4]);
    }
    let mut fa = String::new();
    for i in 0..10 {
        fa.push_str(&format!(">r{i}\n"));
        fa.push_str(&String::from_utf8(genome.clone()).unwrap());
        fa.push('\n');
    }
    fs::write(&reads, fa).unwrap();
    // A 170 bp contig below --min-len 1000 stays unchanged even though the
    // reads would support an extension.
    let mut cf = String::new();
    cf.push_str(">c1\n");
    cf.push_str(std::str::from_utf8(&genome[..170]).unwrap());
    cf.push('\n');
    fs::write(&contigs, cf).unwrap();

    let (_, stderr) = AnchrCmd::new()
        .args(&[
            "asm",
            "extend",
            contigs.to_str().unwrap(),
            reads.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert_eq!(stderr, "");
    let data = fs::read_to_string(&out).unwrap();
    let seq: String = data.lines().filter(|l| !l.starts_with('>')).collect();
    assert_eq!(seq.len(), 170);
}

/// The help text lists the subcommand.
#[test]
fn command_asm_extend_in_help() {
    let (stdout, _) = AnchrCmd::new().args(&["asm", "--help"]).run();
    assert!(stdout.contains("extend"), "asm help must list extend");
}

/// `-o` must not overwrite any input: pointing it at a read file would
/// truncate the reads on disk, so it is rejected up front (along with the
/// contigs collision) before any file is written.
#[test]
fn command_asm_extend_outfile_not_reads() {
    let dir = tempfile::tempdir().unwrap();
    let reads = dir.path().join("reads.fa");
    let contigs = dir.path().join("contigs.fa");

    let mut rng = 2026u64;
    let mut genome = Vec::new();
    for _ in 0..300 {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        genome.push(b"ACGT"[(rng >> 33) as usize % 4]);
    }
    let mut fa = String::new();
    for i in 0..10 {
        fa.push_str(&format!(">r{i}\n"));
        fa.push_str(&String::from_utf8(genome.clone()).unwrap());
        fa.push('\n');
    }
    fs::write(&reads, fa.clone()).unwrap();
    let cf = format!(">c1\n{}\n", std::str::from_utf8(&genome[..170]).unwrap());
    fs::write(&contigs, cf).unwrap();

    let (_, stderr) = AnchrCmd::new()
        .args(&[
            "asm",
            "extend",
            contigs.to_str().unwrap(),
            reads.to_str().unwrap(),
            "-o",
            reads.to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("output file"),
        "expected an output-collision error, got: {stderr}"
    );
    // The read file must be untouched (the check runs before the writer).
    assert_eq!(fs::read_to_string(&reads).unwrap(), fa);
}

/// `--min-support 0` would append bases with zero read support (the all-zero
/// extension counts pass the `>= 2x runner-up` majority check), silently
/// extending contigs through unsupported sequence, so it is rejected.
#[test]
fn command_asm_extend_rejects_zero_min_support() {
    let dir = tempfile::tempdir().unwrap();
    let reads = dir.path().join("reads.fa");
    let contigs = dir.path().join("contigs.fa");
    let out = dir.path().join("out.fa");

    fs::write(&reads, ">r1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();
    fs::write(&contigs, ">c1\nACGTACGTACGTACGTACGTACGTACGTACGT\n").unwrap();

    let (_, stderr) = AnchrCmd::new()
        .args(&[
            "asm",
            "extend",
            contigs.to_str().unwrap(),
            reads.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--min-support",
            "0",
        ])
        .run_fail();
    assert!(
        stderr.contains("min_support must be >= 1"),
        "expected a min_support error, got: {stderr}"
    );
}

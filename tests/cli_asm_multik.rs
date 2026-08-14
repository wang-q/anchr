#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::AnchrCmd;
use std::fs;

/// Parses `unitig_<id>,len=...,cov=...` headers into (id, len, cov).
fn parse_unitigs(data: &str) -> Vec<(usize, usize, f32)> {
    data.lines()
        .filter(|l| l.starts_with('>'))
        .map(|l| {
            let head = l.trim_start_matches('>');
            let mut fields = head.split(',');
            let id: usize = fields
                .next()
                .unwrap()
                .trim_start_matches("unitig_")
                .parse()
                .unwrap();
            let len: usize = fields
                .next()
                .unwrap()
                .trim_start_matches("len=")
                .parse()
                .unwrap();
            let cov: f32 = fields
                .next()
                .unwrap()
                .trim_start_matches("cov=")
                .parse()
                .unwrap();
            (id, len, cov)
        })
        .collect()
}

/// Builds a small random synthetic genome, emits 30× identical full-length
/// reads, and checks `asm multik` compacts the reads into one full-genome
/// unitig (no repeats in the genome, so k=21 unitigs span it and k=51
/// validation keeps the junctions).
#[test]
fn command_asm_multik_compacts_synthetic_genome() {
    let dir = tempfile::tempdir().unwrap();
    let reads = dir.path().join("reads.fa");
    // Fixed-seed random 500 bp genome (no long repeats).
    let mut rng = 42u64;
    let mut genome = Vec::new();
    for _ in 0..500 {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        genome.push(b"ACGT"[(rng >> 33) as usize % 4]);
    }
    let mut fa = String::new();
    for i in 0..30 {
        fa.push_str(&format!(">r{i}\n"));
        fa.push_str(&String::from_utf8(genome.clone()).unwrap());
        fa.push('\n');
    }
    fs::write(&reads, fa).unwrap();
    let out = dir.path().join("out.fa");
    let (_, stderr) = AnchrCmd::new()
        .args(&[
            "asm",
            "multik",
            reads.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "21,51",
        ])
        .run();
    assert_eq!(stderr, "");
    let data = fs::read_to_string(&out).unwrap();
    let u = parse_unitigs(&data);
    assert!(!u.is_empty(), "expected at least one unitig");
    // A clean random genome compacts to a single full-length unitig.
    let longest = u.iter().map(|&(_, l, _)| l).max().unwrap();
    assert!(
        longest == genome.len(),
        "longest unitig {longest} != genome length {}",
        genome.len()
    );
}

/// Bubble selection: a low-coverage divergent branch is pruned by the
/// cross-round validation while the high-coverage main path stays intact.
#[test]
fn command_asm_multik_prunes_low_coverage_branch() {
    let dir = tempfile::tempdir().unwrap();
    let reads = dir.path().join("reads.fa");
    // Main path: A + M + B; divergent branch: A + X + B (same flanks, so
    // the branch bubbles at the A/M and M/B junctions).
    let mut rng = 7u64;
    let mut rand = |n: usize| {
        let mut s = Vec::new();
        for _ in 0..n {
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            s.push(b"ACGT"[(rng >> 33) as usize % 4]);
        }
        s
    };
    let a = rand(120);
    let m = rand(160);
    let x = rand(160);
    let b = rand(120);
    let main: Vec<u8> = [&a[..], &m[..], &b[..]].concat();
    let branch: Vec<u8> = [&a[..], &x[..], &b[..]].concat();
    let mut fa = String::new();
    for i in 0..30 {
        fa.push_str(&format!(">m{i}\n"));
        fa.push_str(&String::from_utf8(main.clone()).unwrap());
        fa.push('\n');
    }
    for i in 0..2 {
        fa.push_str(&format!(">b{i}\n"));
        fa.push_str(&String::from_utf8(branch.clone()).unwrap());
        fa.push('\n');
    }
    fs::write(&reads, fa).unwrap();
    let out = dir.path().join("out.fa");
    let (_, stderr) = AnchrCmd::new()
        .args(&[
            "asm",
            "multik",
            reads.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "21,51",
        ])
        .run();
    assert_eq!(stderr, "");
    let data = fs::read_to_string(&out).unwrap();
    let u = parse_unitigs(&data);
    let longest = u.iter().map(|&(_, l, _)| l).max().unwrap();
    // The main path (400 bp) must survive as the longest unitig.
    assert!(
        longest as f64 >= main.len() as f64 * 0.95,
        "main path not retained: longest unitig {longest}, main len {}",
        main.len()
    );
}

/// High-coverage strain divergence (both branches solid at the larger k)
/// is resolved by the progressive abundance filter: the low-abundance
/// branch drops out and the main path compacts into one full-length unitig,
/// while the strain branch stays as an independent unitig.
#[test]
fn command_asm_multik_resolves_high_coverage_strain_bubble() {
    let dir = tempfile::tempdir().unwrap();
    let reads = dir.path().join("reads.fa");
    let mut rng = 7u64;
    let mut rand = |n: usize| {
        let mut s = Vec::new();
        for _ in 0..n {
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            s.push(b"ACGT"[(rng >> 33) as usize % 4]);
        }
        s
    };
    let a = rand(150);
    let m = rand(200);
    let x = rand(200);
    let b = rand(150);
    let main: Vec<u8> = [&a[..], &m[..], &b[..]].concat();
    let branch: Vec<u8> = [&a[..], &x[..], &b[..]].concat();
    let mut fa = String::new();
    for i in 0..30 {
        fa.push_str(&format!(">m{i}\n"));
        fa.push_str(&String::from_utf8(main.clone()).unwrap());
        fa.push('\n');
    }
    for i in 0..5 {
        fa.push_str(&format!(">b{i}\n"));
        fa.push_str(&String::from_utf8(branch.clone()).unwrap());
        fa.push('\n');
    }
    fs::write(&reads, fa).unwrap();
    let out = dir.path().join("out.fa");
    let (_, stderr) = AnchrCmd::new()
        .args(&[
            "asm",
            "multik",
            reads.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "21,51",
        ])
        .run();
    assert_eq!(stderr, "");
    let data = fs::read_to_string(&out).unwrap();
    let u = parse_unitigs(&data);
    let longest = u.iter().map(|&(_, l, _)| l).max().unwrap();
    assert_eq!(
        longest,
        main.len(),
        "main path must compact to full length {} (got {longest})",
        main.len()
    );
    // The strain branch is retained as an independent unitig.
    let strain = u.iter().any(|&(_, l, cov)| l >= 200 && cov < 10.0);
    assert!(strain, "strain branch must survive as a separate unitig");
}

/// The help text lists the subcommand.
#[test]
fn command_asm_multik_in_help() {
    let (stdout, _) = AnchrCmd::new().args(&["asm", "--help"]).run();
    assert!(stdout.contains("multik"), "asm help must list multik");
}

/// `--kmer auto` (default) derives the k sequence from the read-length N50:
/// long reads compact into a single near-full-length unitig.
#[test]
fn command_asm_multik_auto_k_long_reads() {
    let dir = tempfile::tempdir().unwrap();
    let reads = dir.path().join("reads.fa");
    let mut rng = 2026u64;
    let mut genome = Vec::new();
    for _ in 0..20000 {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        genome.push(b"ACGT"[(rng >> 33) as usize % 4]);
    }
    let glen = genome.len();
    let rl = 3000usize;
    let mut fa = String::new();
    for i in 0..200 {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let start = (rng >> 33) as usize % (glen - rl + 1);
        fa.push_str(&format!(">r{i}\n"));
        fa.push_str(std::str::from_utf8(&genome[start..start + rl]).unwrap());
        fa.push('\n');
    }
    fs::write(&reads, fa).unwrap();
    let out = dir.path().join("out.fa");
    let (_, stderr) = AnchrCmd::new()
        .args(&[
            "asm",
            "multik",
            reads.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert_eq!(stderr, "");
    let data = fs::read_to_string(&out).unwrap();
    let u = parse_unitigs(&data);
    let longest = u.iter().map(|&(_, l, _)| l).max().unwrap();
    assert!(
        longest as f64 >= glen as f64 * 0.97,
        "auto-k long reads must compact to >=97% of the genome (got {longest}/{glen})"
    );
}

/// Full circular coverage (reads wrapping the genome origin) compacts into
/// a single unitig whose length equals the genome — the decisive "N-free
/// chromosome" end-to-end check (k-mer multiset equals the genome).
#[test]
fn command_asm_multik_full_coverage_single_contig() {
    let dir = tempfile::tempdir().unwrap();
    let reads = dir.path().join("reads.fa");
    let mut rng = 42u64;
    let mut genome = Vec::new();
    for _ in 0..20000 {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        genome.push(b"ACGT"[(rng >> 33) as usize % 4]);
    }
    let glen = genome.len();
    let rl = 3000usize;
    let mut fa = String::new();
    for i in 0..180 {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let start = (rng >> 33) as usize % (glen - rl + 1);
        fa.push_str(&format!(">r{i}\n"));
        fa.push_str(std::str::from_utf8(&genome[start..start + rl]).unwrap());
        fa.push('\n');
    }
    // Reads wrapping the origin and reads starting at the origin cover the
    // genome ends, closing the circular chromosome.
    for i in 0..20 {
        rng = rng
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let max_start = glen - 800;
        let start = if max_start > glen - rl {
            glen - rl
        } else {
            max_start + (rng >> 33) as usize % 800
        };
        let mut seq = genome[start..].to_vec();
        seq.extend_from_slice(&genome[..rl - (glen - start)]);
        fa.push_str(&format!(">w{i}\n"));
        fa.push_str(std::str::from_utf8(&seq).unwrap());
        fa.push('\n');
    }
    for i in 0..10 {
        fa.push_str(&format!(">o{i}\n"));
        fa.push_str(std::str::from_utf8(&genome[..rl]).unwrap());
        fa.push('\n');
    }
    fs::write(&reads, fa).unwrap();
    let out = dir.path().join("out.fa");
    let (_, stderr) = AnchrCmd::new()
        .args(&[
            "asm",
            "multik",
            reads.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert_eq!(stderr, "");
    let data = fs::read_to_string(&out).unwrap();
    let u = parse_unitigs(&data);
    let longest = u.iter().map(|&(_, l, _)| l).max().unwrap();
    assert_eq!(
        longest, glen,
        "full circular coverage must compact into one genome-length unitig (got {longest}/{glen})"
    );
    assert_eq!(u.len(), 1, "expected a single unitig (got {})", u.len());
}

/// Input without enough k-mers yields an empty output (no panic).
#[test]
fn command_asm_multik_empty_input_is_error() {
    let dir = tempfile::tempdir().unwrap();
    let reads = dir.path().join("empty.fa");
    fs::write(&reads, ">r1\nACGT\n").unwrap();
    let out = dir.path().join("out.fa");
    let (_, stderr) = AnchrCmd::new()
        .args(&[
            "asm",
            "multik",
            reads.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
        ])
        .run();
    assert_eq!(stderr, "");
    assert!(out.exists(), "output file must be created");
    assert!(
        fs::read_to_string(&out).unwrap().is_empty(),
        "too-short reads must produce no unitigs"
    );
}

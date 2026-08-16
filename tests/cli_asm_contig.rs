#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::AnchrCmd;
use std::io::Read;

/// Parses a FASTA into (header, sequence) pairs.
fn parse_fa(data: &[u8]) -> Vec<(String, String)> {
    let mut recs = Vec::new();
    let mut cur: Option<(String, String)> = None;
    for line in std::str::from_utf8(data).unwrap().lines() {
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(c) = cur.take() {
                recs.push(c);
            }
            cur = Some((rest.to_string(), String::new()));
        } else if let Some(c) = cur.as_mut() {
            c.1.push_str(line);
        }
    }
    if let Some(c) = cur {
        recs.push(c);
    }
    recs
}

fn read_gz(path: &str) -> Vec<u8> {
    let mut out = Vec::new();
    let mut dec = flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap());
    dec.read_to_end(&mut out).unwrap();
    out
}

/// The contig set is deterministic and statistically matches the BBTools
/// `tadpole.sh threads=1` golden. The pre-pop contig set is byte-identical;
/// bubble-popping resolutions can differ because BBTools' expand order
/// depends on its memory-dependent hash-table layout (see
/// notes/design/asm-assemble.md), so the popped output is compared by
/// sequence set rather than byte-for-byte.
#[test]
fn command_asm_contig_matches_tadpole_contig_set() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("contigs.fa");
    AnchrCmd::new()
        .args(&[
            "asm",
            "contig",
            "tests/bbtools/Lambda/R1.2k.fq.gz",
            "tests/bbtools/Lambda/R2.2k.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
        ])
        .assert()
        .success();
    let pgr = parse_fa(&std::fs::read(&out).unwrap());
    let golden = parse_fa(&read_gz(
        "tests/bbtools/Lambda/golden/tadpole_contigs31.fasta.gz",
    ));

    // Same total assembled bases as the reference.
    let pgr_bases: usize = pgr.iter().map(|(_, s)| s.len()).sum();
    let golden_bases: usize = golden.iter().map(|(_, s)| s.len()).sum();
    // Bubble-resolution differences can shift a few contigs between the
    // kept and merged sets, so allow a small total-base delta.
    assert!(
        (pgr_bases as i64 - golden_bases as i64).abs() <= 100,
        "bases {pgr_bases} vs {golden_bases}"
    );

    // Contig count within 1 of the reference (pgr is deterministic; the
    // remaining bubble-resolution differences are documented).
    assert!(
        (pgr.len() as i64 - golden.len() as i64).abs() <= 1,
        "{} vs {}",
        pgr.len(),
        golden.len()
    );

    // At least 90% of the reference contigs are present verbatim.
    let golden_seqs: std::collections::HashSet<&str> =
        golden.iter().map(|(_, s)| s.as_str()).collect();
    let shared = pgr
        .iter()
        .filter(|(_, s)| golden_seqs.contains(s.as_str()))
        .count();
    assert!(
        shared * 10 >= golden.len() * 9,
        "shared {shared}/{}",
        golden.len()
    );
}

/// Repeated runs produce byte-identical output (deterministic scan order).
#[test]
fn command_asm_contig_is_deterministic() {
    let out_dir = tempfile::tempdir().unwrap();
    let out1 = out_dir.path().join("a.fa");
    let out2 = out_dir.path().join("b.fa");
    for out in [&out1, &out2] {
        AnchrCmd::new()
            .args(&[
                "asm",
                "contig",
                "tests/bbtools/Lambda/R1.2k.fq.gz",
                "tests/bbtools/Lambda/R2.2k.fq.gz",
                "-o",
                out.to_str().unwrap(),
                "--kmer",
                "31",
            ])
            .assert()
            .success();
    }
    assert_eq!(std::fs::read(&out1).unwrap(), std::fs::read(&out2).unwrap());
}

/// `--no-bubbles` keeps the pre-pop contig set (no parallel paths merged):
/// at least as many contigs as the default bubble-popped output, and the
/// no-bubbles run is byte-identical across repeated invocations.
#[test]
fn command_asm_contig_no_bubbles_keeps_parallel_contigs() {
    let out_dir = tempfile::tempdir().unwrap();
    let default_out = out_dir.path().join("default.fa");
    let nb1 = out_dir.path().join("nb1.fa");
    let nb2 = out_dir.path().join("nb2.fa");
    AnchrCmd::new()
        .args(&[
            "asm",
            "contig",
            "tests/bbtools/Lambda/R1.2k.fq.gz",
            "tests/bbtools/Lambda/R2.2k.fq.gz",
            "-o",
            default_out.to_str().unwrap(),
            "--kmer",
            "31",
        ])
        .assert()
        .success();
    for out in [&nb1, &nb2] {
        AnchrCmd::new()
            .args(&[
                "asm",
                "contig",
                "tests/bbtools/Lambda/R1.2k.fq.gz",
                "tests/bbtools/Lambda/R2.2k.fq.gz",
                "-o",
                out.to_str().unwrap(),
                "--kmer",
                "31",
                "--no-bubbles",
            ])
            .assert()
            .success();
    }
    let default_seqs = parse_fa(&std::fs::read(&default_out).unwrap());
    let nb_seqs = parse_fa(&std::fs::read(&nb1).unwrap());
    assert!(
        nb_seqs.len() >= default_seqs.len(),
        "no-bubbles contigs {} < default {}",
        nb_seqs.len(),
        default_seqs.len()
    );
    assert_eq!(std::fs::read(&nb1).unwrap(), std::fs::read(&nb2).unwrap());
}

/// A zero k-mer length must fail cleanly instead of panicking.
#[test]
fn command_asm_contig_rejects_zero_kmer() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fa");
    AnchrCmd::new()
        .args(&[
            "asm",
            "contig",
            "tests/bbtools/Lambda/R1.2k.fq.gz",
            "tests/bbtools/Lambda/R2.2k.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "0",
        ])
        .assert()
        .failure();
}

/// A k-mer above the 128-base key limit must fail cleanly instead of
/// panicking in `Kmer::new().expect()` (zero-panic policy).
#[test]
fn command_asm_contig_rejects_kmer_above_limit() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("out.fa");
    AnchrCmd::new()
        .args(&[
            "asm",
            "contig",
            "tests/bbtools/Lambda/R1.2k.fq.gz",
            "tests/bbtools/Lambda/R2.2k.fq.gz",
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "257",
        ])
        .assert()
        .failure();
}

/// `-o` must not overwrite an input file (the writer is opened before the
/// reads are consumed).
#[test]
fn command_asm_contig_outfile_not_input() {
    let infile = "tests/bbtools/Lambda/R1.2k.fq.gz";
    AnchrCmd::new()
        .args(&["asm", "contig", infile, "-o", infile])
        .assert()
        .failure()
        .stderr(predicates::str::contains("is also an input file"));
}

/// Assembles a small synthetic repeat into contigs.
#[test]
fn command_asm_contig_small_synthetic() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let out = out_dir.path().join("out.fa");
    // Two identical 60 bp reads: the k=31 graph assembles them into a
    // single 60 bp contig (below the default 124 bp output threshold, so
    // use --min-contig-len 1).
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    std::fs::write(&infile, format!(">r1\n{seq}\n>r2\n{seq}\n")).unwrap();
    AnchrCmd::new()
        .args(&[
            "asm",
            "contig",
            infile.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
        ])
        .assert()
        .success();
    let recs = parse_fa(&std::fs::read(&out).unwrap());
    assert!(!recs.is_empty());
    assert!(recs.iter().any(|(_, s)| s.contains("ACGTACGT")));
}

/// Raising the seeding threshold drops low-depth k-mers (tadpole
/// `mincountseed`).
#[test]
fn command_asm_contig_min_count_seed() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let seq = "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGA";
    // 2 identical reads: every k-mer has count 2 (not seeded at the default
    // threshold of 3, seeded at --min-count-seed 2).
    let fa = format!(">r1\n{seq}\n>r2\n{seq}\n");
    std::fs::write(&infile, fa).unwrap();
    let default_out = out_dir.path().join("default.fa");
    AnchrCmd::new()
        .args(&[
            "asm",
            "contig",
            infile.to_str().unwrap(),
            "-o",
            default_out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
        ])
        .assert()
        .success();
    assert!(parse_fa(&std::fs::read(&default_out).unwrap()).is_empty());
    let low_out = out_dir.path().join("low.fa");
    AnchrCmd::new()
        .args(&[
            "asm",
            "contig",
            infile.to_str().unwrap(),
            "-o",
            low_out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
            "--min-count-seed",
            "2",
        ])
        .assert()
        .success();
    assert!(!parse_fa(&std::fs::read(&low_out).unwrap()).is_empty());
}

/// Raising the minimum coverage drops low-depth contigs (tadpole
/// `mincoverage`).
#[test]
fn command_asm_contig_min_coverage() {
    let out_dir = tempfile::tempdir().unwrap();
    let infile = out_dir.path().join("in.fa");
    let seq = "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGA";
    // 4 identical reads: every k-mer has count 4 (mean coverage ~4).
    let fa = format!(">r1\n{seq}\n>r2\n{seq}\n>r3\n{seq}\n>r4\n{seq}\n");
    std::fs::write(&infile, fa).unwrap();
    let default_out = out_dir.path().join("default.fa");
    AnchrCmd::new()
        .args(&[
            "asm",
            "contig",
            infile.to_str().unwrap(),
            "-o",
            default_out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
        ])
        .assert()
        .success();
    assert!(!parse_fa(&std::fs::read(&default_out).unwrap()).is_empty());
    let strict_out = out_dir.path().join("strict.fa");
    AnchrCmd::new()
        .args(&[
            "asm",
            "contig",
            infile.to_str().unwrap(),
            "-o",
            strict_out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
            "--min-coverage",
            "5",
        ])
        .assert()
        .success();
    assert!(parse_fa(&std::fs::read(&strict_out).unwrap()).is_empty());
}

/// Multiple positional files and `--list-files` give identical output;
/// unpaired reads (odd record count) are accepted like bcalm.
#[test]
fn command_asm_contig_multiple_files_and_list() {
    let out_dir = tempfile::tempdir().unwrap();
    let a = out_dir.path().join("a.fa");
    let b = out_dir.path().join("b.fa");
    let list = out_dir.path().join("files.list");
    let out_direct = out_dir.path().join("direct.fa");
    let out_list = out_dir.path().join("list.fa");
    let seq1 = "AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGA";
    let seq2 = "TGCCCAAGTTAGTTGCTCGGTAGGTCGAAACTATCCCGGACCGTAACGCACCGAAACGT";
    std::fs::write(
        &a,
        format!(">a1\n{seq1}\n>a2\n{seq1}\n>a3\n{seq2}\n>a4\n{seq2}\n"),
    )
    .unwrap();
    std::fs::write(
        &b,
        format!(">b1\n{seq1}\n>b2\n{seq1}\n>b3\n{seq2}\n>b4\n{seq2}\n"),
    )
    .unwrap();
    std::fs::write(&list, format!("{}\n{}\n", a.display(), b.display())).unwrap();
    for (out, extra) in [
        (&out_direct, Vec::<&str>::new()),
        (&out_list, vec!["--list-files"]),
    ] {
        let mut args = vec![
            "asm",
            "contig",
            if extra.is_empty() {
                a.to_str().unwrap()
            } else {
                list.to_str().unwrap()
            },
        ];
        if extra.is_empty() {
            args.push(b.to_str().unwrap());
        }
        args.extend([
            "-o",
            out.to_str().unwrap(),
            "--kmer",
            "31",
            "--min-contig-len",
            "1",
        ]);
        args.extend(extra);
        AnchrCmd::new().args(&args).assert().success();
    }
    assert_eq!(
        std::fs::read(&out_direct).unwrap(),
        std::fs::read(&out_list).unwrap(),
        "direct files and --list-files outputs differ"
    );
}

#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::AnchrCmd;
use std::fs;

/// `anchr asm anchor` keeps the well-covered unitig as an anchor and filters
/// the low-coverage one out (perfect-match read mapping + coverage window).
/// The 100-bp contig equals the read length, so the edge ramp (read_len/2
/// from each end) clips the outermost bases; the middle stays an anchor.
#[test]
fn command_asm_anchor_filters_low_coverage() {
    let dir = tempfile::tempdir().unwrap();
    let ut = dir.path().join("ut.fa");
    let reads = dir.path().join("reads.fa");
    let a100 = "A".repeat(100);
    let c100 = "C".repeat(100);
    fs::write(
        &ut,
        format!(">unitig_1,len=100,cov=5\n{a100}\n>unitig_2,len=100,cov=1\n{c100}\n"),
    )
    .unwrap();
    // 5 full-length reads of unitig_1, 1 read of unitig_2.
    let mut fa = String::new();
    for _ in 0..5 {
        fa.push_str(&format!(">r\n{a100}\n"));
    }
    fa.push_str(&format!(">r\n{c100}\n"));
    fs::write(&reads, fa).unwrap();
    let out = dir.path().join("anchors.fa");
    AnchrCmd::new()
        .args(&[
            "asm",
            "anchor",
            ut.to_str().unwrap(),
            reads.to_str().unwrap(),
            "--mincov",
            "3",
            "--min-anchor-len",
            "10",
            "-o",
            out.to_str().unwrap(),
        ])
        .assert()
        .success();
    let data = fs::read_to_string(&out).unwrap();
    let seqs: Vec<String> = data
        .lines()
        .filter(|l| !l.starts_with('>') && !l.is_empty())
        .map(|l| l.to_string())
        .collect();
    let joined = seqs.join("");
    assert!(!joined.is_empty(), "well-covered unitig yields an anchor");
    assert!(
        joined.bytes().all(|b| b == b'A'),
        "only the A unitig is anchored"
    );
}

/// The command is listed under `anchr asm --help`.
#[test]
fn command_asm_anchor_in_help() {
    let (stdout, _) = AnchrCmd::new().args(&["asm", "--help"]).run();
    assert!(stdout.contains("anchor"), "asm help must list anchor");
}

/// `--stats` is a second output: pointing it at the same path as `-o` would
/// truncate the first writer's output, so it must be rejected up front.
#[test]
fn command_asm_anchor_stats_not_outfile() {
    let dir = tempfile::tempdir().unwrap();
    let ut = dir.path().join("ut.fa");
    let reads = dir.path().join("reads.fa");
    let a100 = "A".repeat(100);
    fs::write(&ut, format!(">unitig_1,len=100,cov=5\n{a100}\n")).unwrap();
    fs::write(&reads, format!(">r1\n{a100}\n>r2\n{a100}\n>r3\n{a100}\n")).unwrap();
    let same = dir.path().join("same.tsv");
    let (_, stderr) = AnchrCmd::new()
        .args(&[
            "asm",
            "anchor",
            ut.to_str().unwrap(),
            reads.to_str().unwrap(),
            "--mincov",
            "1",
            "--min-anchor-len",
            "1",
            "-o",
            same.to_str().unwrap(),
            "--stats",
            same.to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("output files must be distinct"),
        "expected a distinct-output error, got: {stderr}"
    );
}

/// `--lscale`/`--uscale` are the coverage-window denominator/multiplier: 0
/// yields NaN/inf bounds and a silently empty anchor set, so reject it.
#[test]
fn command_asm_anchor_rejects_zero_lscale() {
    let dir = tempfile::tempdir().unwrap();
    let ut = dir.path().join("ut.fa");
    let reads = dir.path().join("reads.fa");
    fs::write(&ut, ">u1\nACGTACGTACGTACGTACGT\n").unwrap();
    fs::write(&reads, ">r1\nACGTACGTACGTACGTACGT\n").unwrap();
    let out = dir.path().join("out.fa");
    let (_, stderr) = AnchrCmd::new()
        .args(&[
            "asm",
            "anchor",
            ut.to_str().unwrap(),
            reads.to_str().unwrap(),
            "--lscale",
            "0",
            "-o",
            out.to_str().unwrap(),
        ])
        .run_fail();
    assert!(
        stderr.contains("--lscale must be > 0"),
        "expected an lscale error, got: {stderr}"
    );
}

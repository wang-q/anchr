#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::AnchrCmd;
use std::fs;

/// `anchr asm anchor` keeps the well-covered unitig as an anchor and filters
/// the low-coverage one out (perfect-match read mapping + coverage window).
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
    assert_eq!(
        joined,
        "A".repeat(100),
        "only the well-covered unitig is an anchor"
    );
}

/// The command is listed under `anchr asm --help`.
#[test]
fn command_asm_anchor_in_help() {
    let (stdout, _) = AnchrCmd::new().args(&["asm", "--help"]).run();
    assert!(stdout.contains("anchor"), "asm help must list anchor");
}

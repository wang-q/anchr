#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::AnchrCmd;
use std::io::Read;

fn read_gz(path: &str) -> Vec<u8> {
    let mut out = Vec::new();
    let mut dec = flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap());
    dec.read_to_end(&mut out).unwrap();
    out
}

fn lambda(args: &[&str], out: &str, outu: Option<&str>, ihist: Option<&str>) {
    let mut full: Vec<&str> = vec![
        "fq",
        "ec-overlap",
        "tests/bbtools/Lambda/R1.2k.fq.gz",
        "tests/bbtools/Lambda/R2.2k.fq.gz",
        "-o",
        out,
    ];
    if let Some(u) = outu {
        full.push("--outu");
        full.push(u);
    }
    if let Some(h) = ihist {
        full.push("--ihist");
        full.push(h);
    }
    full.extend_from_slice(args);
    AnchrCmd::new().args(&full).assert().success();
}

#[test]
fn command_fq_ec_overlap_rejects_outu_same_as_outfile() {
    // Data safety: --outu equal to -o would open two writers to the same path
    // and corrupt the output; it must be rejected up front.
    let (_, stderr) = AnchrCmd::new()
        .args(&[
            "fq",
            "ec-overlap",
            "tests/bbtools/Lambda/R1.2k.fq.gz",
            "tests/bbtools/Lambda/R2.2k.fq.gz",
            "-o",
            "out.fq",
            "--outu",
            "out.fq",
        ])
        .run_fail();
    assert!(stderr.contains("must be distinct"), "stderr: {stderr}");
}

#[test]
fn command_fq_ec_overlap_novector_matches_bbtools_golden() {
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("ecco.fq");
    lambda(&["--vstrict"], out.to_str().unwrap(), None, None);
    assert_eq!(
        std::fs::read(&out).unwrap(),
        read_gz("tests/bbtools/Lambda/golden/merge.novector.ecco.fq.gz")
    );
}

#[test]
fn command_fq_ec_overlap_no_mix_writes_only_corrected_pairs() {
    // `bbmerge ... ecco mix=f`: only overlapping pairs are corrected; the
    // rest are dropped when no --outu is given. Classic overlap filters
    // correct 15 of 2000 pairs on the Lambda subset (120 FASTQ lines).
    let out_dir = tempfile::tempdir().unwrap();
    let out = out_dir.path().join("ecco.fq");
    lambda(
        &["--no-mix", "--vstrict"],
        out.to_str().unwrap(),
        None,
        None,
    );
    let out = std::fs::read(&out).unwrap();
    assert_eq!(std::str::from_utf8(&out).unwrap().lines().count(), 120);
}

#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::AnchrCmd;

const SAM: &str = "\
@HD\tVN:1.6\tSO:unknown
@SQ\tSN:ut\tLN:1000
r1\t0\tut\t101\t255\t50M\t*\t0\t0\tACGT\tIIII
r2\t16\tut\t201\t255\t50M\t*\t0\t0\tACGT\tIIII
r3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII
";

/// Mapped records become 1-based inclusive ranges; unmapped records and
/// header lines are skipped.
#[test]
fn command_sam_to_rg_basic() {
    let out_dir = tempfile::tempdir().unwrap();
    let sam = out_dir.path().join("in.sam");
    std::fs::write(&sam, SAM).unwrap();
    let (stdout, _) = AnchrCmd::new()
        .args(&["sam", "to-rg", sam.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "ut:101-150\nut:201-250\n");
}

/// The range spans the full reference-consuming CIGAR (M/D/N/=/X), with
/// non-consuming operations (I/S/H/P) excluded.
#[test]
fn command_sam_to_rg_cigar_span() {
    let out_dir = tempfile::tempdir().unwrap();
    let sam = out_dir.path().join("in.sam");
    std::fs::write(
        &sam,
        "r1\t0\tut\t11\t255\t10M2I8M3D15N5=2X\t*\t0\t0\tACGT\tIIII\n",
    )
    .unwrap();
    let (stdout, _) = AnchrCmd::new()
        .args(&["sam", "to-rg", sam.to_str().unwrap()])
        .run();
    // 10 + 8 + 3 + 15 + 5 + 2 = 43 reference bases.
    assert_eq!(stdout, "ut:11-53\n");
}

/// Reads SAM from stdin.
#[test]
fn command_sam_to_rg_stdin() {
    let (stdout, _) = AnchrCmd::new()
        .args(&["sam", "to-rg", "stdin"])
        .stdin("r1\t0\tut\t5\t255\t10M\t*\t0\t0\tACGT\tIIII\n")
        .run();
    assert_eq!(stdout, "ut:5-14\n");
}

/// Malformed records are skipped by default and rejected with `--strict`.
#[test]
fn command_sam_to_rg_strict() {
    let out_dir = tempfile::tempdir().unwrap();
    let sam = out_dir.path().join("bad.sam");
    std::fs::write(
        &sam,
        "r1\t0\tut\t5\t255\t10M\t*\t0\t0\tACGT\tIIII\nonly_two_fields\n",
    )
    .unwrap();
    let (stdout, _) = AnchrCmd::new()
        .args(&["sam", "to-rg", sam.to_str().unwrap()])
        .run();
    assert_eq!(stdout, "ut:5-14\n");
    let (_, stderr) = AnchrCmd::new()
        .args(&["sam", "to-rg", sam.to_str().unwrap(), "--strict"])
        .run_fail();
    assert!(stderr.contains("malformed SAM record"), "stderr: {stderr}");
}

const SEQ50: &str = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC";

fn sam_record(name: &str, flag: u16, pos: u32, mate_pos: u32, tlen: i32) -> String {
    format!(
        "{name}\t{flag}\tut\t{pos}\t255\t50M\t=\t{mate_pos}\t{tlen}\t{SEQ50}\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
    )
}

/// Insert-size stats and histogram from a paired SAM (proper FR pairs).
#[test]
fn command_sam_ihist_basic() {
    let out_dir = tempfile::tempdir().unwrap();
    let sam = out_dir.path().join("in.sam");
    std::fs::write(
        &sam,
        format!(
            "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ut\tLN:1000\n{}\n{}\n{}\n{}\n",
            sam_record("p1", 67, 11, 31, 70),
            sam_record("p1", 147, 31, 11, -70),
            sam_record("p2", 67, 101, 121, 70),
            sam_record("p2", 147, 121, 101, -70),
        ),
    )
    .unwrap();
    let (stdout, _) = AnchrCmd::new()
        .args(&["sam", "ihist", sam.to_str().unwrap()])
        .run();
    assert_eq!(
        stdout,
        "#Mean\t70.000\n#Median\t70\n#Mode\t70\n#STDev\t0.000\n#PercentOfPairs\t100.000\n#InsertSize\tCount\n70\t2\n"
    );
}

/// Pairs are grouped by name with `/1` `/2` suffixes stripped.
#[test]
fn command_sam_ihist_name_normalization() {
    let out_dir = tempfile::tempdir().unwrap();
    let sam = out_dir.path().join("in.sam");
    std::fs::write(
        &sam,
        format!(
            "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ut\tLN:1000\n{}\n{}\n",
            sam_record("r1/1", 67, 11, 31, 70),
            sam_record("r1/2", 147, 31, 11, -70),
        ),
    )
    .unwrap();
    let (stdout, _) = AnchrCmd::new()
        .args(&["sam", "ihist", sam.to_str().unwrap()])
        .run();
    assert!(
        stdout.contains("#PercentOfPairs\t100.000"),
        "stdout: {stdout}"
    );
    assert!(stdout.contains("70\t1"), "stdout: {stdout}");
}

/// Insert sizes beyond `median ± 10*MAD` are dropped from the stats and the
/// histogram (circular-origin chimeric pairs), while PercentOfPairs still
/// counts every proper pair.
#[test]
fn command_sam_ihist_filters_outliers() {
    let out_dir = tempfile::tempdir().unwrap();
    let sam = out_dir.path().join("in.sam");
    let mut lines = String::from("@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ut\tLN:10000\n");
    for i in 0..10u32 {
        let size = 90 + i;
        lines.push_str(&format!(
            "{}\n",
            sam_record(&format!("p{i}"), 67, 11, 11 + size, size as i32,)
        ));
        lines.push_str(&format!(
            "{}\n",
            sam_record(&format!("p{i}"), 147, 11 + size, 11, -(size as i32),)
        ));
    }
    // One circular-origin chimeric pair (huge insert size).
    lines.push_str(&format!("{}\n", sam_record("px", 67, 1, 9000, 9000)));
    lines.push_str(&format!("{}\n", sam_record("px", 147, 9000, 1, -9000)));
    std::fs::write(&sam, lines).unwrap();
    let (stdout, stderr) = AnchrCmd::new()
        .args(&["sam", "ihist", sam.to_str().unwrap()])
        .run();
    // insert sizes 140-149 (span + read length), median 145, MAD 3,
    // cutoff 30: the 9000 pair is filtered.
    assert!(
        !stdout.contains("9000\t"),
        "outlier must be filtered: {stdout}\nstderr: {stderr}"
    );
    assert!(stdout.contains("#Mean\t144.500"), "stdout: {stdout}");
    assert!(stdout.contains("#Median\t144"), "stdout: {stdout}");
    assert!(
        stdout.contains("#PercentOfPairs\t100.000"),
        "stdout: {stdout}"
    );
}

/// A SAM without proper pairs yields zeroed stats and an empty histogram.
#[test]
fn command_sam_ihist_empty() {
    let out_dir = tempfile::tempdir().unwrap();
    let sam = out_dir.path().join("in.sam");
    std::fs::write(
        &sam,
        "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ut\tLN:1000\nr1\t0\tut\t11\t255\t50M\t*\t0\t0\tACGT\tIIII\n",
    )
    .unwrap();
    let (stdout, _) = AnchrCmd::new()
        .args(&["sam", "ihist", sam.to_str().unwrap()])
        .run();
    assert_eq!(
        stdout,
        "#Mean\t0.000\n#Median\t0\n#Mode\t0\n#STDev\t0.000\n#PercentOfPairs\t0.000\n#InsertSize\tCount\n"
    );
}

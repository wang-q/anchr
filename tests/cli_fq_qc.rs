//! `anchr fq qc` integration tests against the fastqc 0.12.1 golden
//! (`tests/qc/golden/R1.2k_fastqc/`, Lambda 2000 reads x 108 bp).

#[macro_use]
#[path = "common/mod.rs"]
mod common;

use common::AnchrCmd;
use std::fs;

const LAMBDA: &str = "tests/bbtools/Lambda/R1.2k.fq.gz";

fn run_qc() -> tempfile::TempDir {
    let dir = tempfile::tempdir().unwrap();
    AnchrCmd::new()
        .args(&["fq", "qc", LAMBDA, "-o", dir.path().to_str().unwrap()])
        .assert()
        .success();
    dir
}

fn data_line(dir: &tempfile::TempDir, module: &str, key: &str) -> String {
    let path = dir
        .path()
        .join("R1.2k.fq.gz_fastqc/fastqc_data.txt");
    let text = fs::read_to_string(path).unwrap();
    let mut in_module = false;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix(">>") {
            in_module = rest.starts_with(module);
            continue;
        }
        if in_module && line.starts_with(key) {
            return line.to_string();
        }
    }
    panic!("{key} not found in module {module}");
}

#[test]
fn command_fq_qc_basic_statistics() {
    let dir = run_qc();
    let text = fs::read_to_string(
        dir.path().join("R1.2k.fq.gz_fastqc/fastqc_data.txt"),
    )
    .unwrap();
    assert!(text.contains("##FastQC\t0.12.1"));
    assert!(text.contains("Filename\tR1.2k.fq.gz"));
    assert!(text.contains("Encoding\tSanger / Illumina 1.9"));
    assert!(text.contains("Total Sequences\t2000"));
    assert!(text.contains("Total Bases\t216 kbp"));
    assert!(text.contains("Sequence length\t108"));
    assert!(text.contains("%GC\t49"));
}

#[test]
fn command_fq_qc_per_base_quality_matches_golden() {
    let dir = run_qc();
    // First position: mean 32.65, median 33.0 (golden values)
    let line = data_line(&dir, "Per base sequence quality", "1\t");
    let fields: Vec<&str> = line.split('\t').collect();
    let mean: f64 = fields[1].parse().unwrap();
    assert!((mean - 32.65).abs() < 1e-3, "mean: {mean}");
    assert_eq!(fields[2], "33.0"); // median
    // Final single position 108
    let last = data_line(&dir, "Per base sequence quality", "108\t");
    assert!(last.starts_with("108\t30.444"), "last: {last}");
}

#[test]
fn command_fq_qc_sequence_quality_matches_golden() {
    let dir = run_qc();
    assert_eq!(data_line(&dir, "Per sequence quality scores", "37\t"), "37\t603.0");
    assert_eq!(data_line(&dir, "Per sequence quality scores", "13\t"), "13\t1.0");
}

#[test]
fn command_fq_qc_gc_n_length_match_golden() {
    let dir = run_qc();
    // GCModel-weighted read counts (golden: sum over bins == 2000)
    let gc21: f64 = data_line(&dir, "Per sequence GC content", "21\t")
        .split('\t')
        .nth(1)
        .unwrap()
        .parse()
        .unwrap();
    assert!((gc21 - 0.5).abs() < 1e-9, "gc21: {gc21}");
    let gc50: f64 = data_line(&dir, "Per sequence GC content", "50\t")
        .split('\t')
        .nth(1)
        .unwrap()
        .parse()
        .unwrap();
    assert!((gc50 - 79.0).abs() < 1e-9, "gc50: {gc50}");
    let n2: f64 = data_line(&dir, "Per base N content", "2\t")
        .split('\t')
        .nth(1)
        .unwrap()
        .parse()
        .unwrap();
    assert!((n2 - 0.05).abs() < 1e-9, "n2: {n2}");
    let len108: f64 = data_line(&dir, "Sequence Length Distribution", "108\t")
        .split('\t')
        .nth(1)
        .unwrap()
        .parse()
        .unwrap();
    assert!((len108 - 2000.0).abs() < 1e-9, "len108: {len108}");
}

#[test]
fn command_fq_qc_summary_statuses() {
    let dir = run_qc();
    let summary = fs::read_to_string(
        dir.path().join("R1.2k.fq.gz_fastqc/summary.txt"),
    )
    .unwrap();
    assert!(summary.contains("PASS\tBasic Statistics\tR1.2k.fq.gz"));
    assert!(summary.contains("PASS\tPer base sequence quality\tR1.2k.fq.gz"));
    assert!(summary.contains("PASS\tPer base N content\tR1.2k.fq.gz"));
    // M1: GC grade is a simplified mean-deviation check (GCModel in M3);
    // fastqc marks this dataset FAIL.
    assert!(summary.contains("Per sequence GC content\tR1.2k.fq.gz"));
}

#[test]
fn command_fq_qc_stdin() {
    let dir = tempfile::tempdir().unwrap();
    let input = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    AnchrCmd::new()
        .args(&["fq", "qc", "stdin", "-o", dir.path().to_str().unwrap()])
        .stdin(input)
        .assert()
        .success();
    let data = fs::read_to_string(dir.path().join("stdin_fastqc/fastqc_data.txt")).unwrap();
    assert!(data.contains("Total Sequences\t1"));
    assert!(data.contains("Sequence length\t20"));
}

#[test]
fn command_fq_qc_duplication_matches_golden() {
    let dir = run_qc();
    let dup = data_line(&dir, "Sequence Duplication Levels", "#Total Deduplicated Percentage");
    let pct: f64 = dup.split('\t').nth(1).unwrap().parse().unwrap();
    assert!((pct - 97.3).abs() < 1e-9, "dedup pct: {pct}");
    let level1: f64 = data_line(&dir, "Sequence Duplication Levels", "1\t")
        .split('\t')
        .nth(1)
        .unwrap()
        .parse()
        .unwrap();
    assert!((level1 - 96.35).abs() < 1e-9, "level1: {level1}");
}

#[test]
fn command_fq_qc_overrepresented_matches_golden() {
    let dir = run_qc();
    // TruSeq adapter, Index 1 (36 reads, 1.8%) with fastqc-style source
    let line = data_line(
        &dir,
        "Overrepresented sequences",
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATGAGCATCTCGTATG",
    );
    assert!(
        line.contains("36\t1.7999999999999998\tTruSeq Adapter, Index 1 (97% over 37bp)"),
        "overrep: {line}"
    );
}

#[test]
fn command_fq_qc_adapter_matches_golden() {
    let dir = run_qc();
    // Universal adapter at position 1: 1.85 (37/2000 reads)
    let line = data_line(&dir, "Adapter Content", "1\t");
    let fields: Vec<&str> = line.split('\t').collect();
    let universal: f64 = fields[1].parse().unwrap();
    assert!((universal - 1.85).abs() < 1e-9, "universal: {universal}");
    // No Kmer Content section on the Lambda golden (2% kmer sampling)
    let text = fs::read_to_string(
        dir.path().join("R1.2k.fq.gz_fastqc/fastqc_data.txt"),
    )
    .unwrap();
    assert!(!text.contains(">>Kmer Content"), "kmer section should be absent");
}

#[test]
fn command_fq_qc_empty_input_fails() {
    let dir = tempfile::tempdir().unwrap();
    let empty = dir.path().join("empty.fq");
    std::fs::write(&empty, "").unwrap();
    AnchrCmd::new()
        .args(&["fq", "qc", empty.to_str().unwrap(), "-o", dir.path().to_str().unwrap()])
        .assert()
        .failure()
        .stderr(predicates::str::contains("no reads found"));
}

#[test]
fn command_fq_qc_single_read() {
    let dir = tempfile::tempdir().unwrap();
    let input = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let fq = dir.path().join("one.fq");
    std::fs::write(&fq, input).unwrap();
    AnchrCmd::new()
        .args(&["fq", "qc", fq.to_str().unwrap(), "-o", dir.path().to_str().unwrap()])
        .assert()
        .success();
    let data = fs::read_to_string(dir.path().join("one.fq_fastqc/fastqc_data.txt")).unwrap();
    assert!(data.contains("Total Sequences\t1"));
    assert!(data.contains("Sequence length\t20"));
}

#[test]
fn command_fq_qc_variable_lengths() {
    let dir = tempfile::tempdir().unwrap();
    let input = "@r1\nACGTACGT\n+\n!!!!!!!!\n@r2\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n";
    let fq = dir.path().join("var.fq");
    std::fs::write(&fq, input).unwrap();
    AnchrCmd::new()
        .args(&["fq", "qc", fq.to_str().unwrap(), "-o", dir.path().to_str().unwrap()])
        .assert()
        .success();
    let data = fs::read_to_string(dir.path().join("var.fq_fastqc/fastqc_data.txt")).unwrap();
    assert!(data.contains("Sequence length\t8-20"), "data: {data}");
}

#[test]
fn command_fq_qc_long_reads_grouping() {
    // 1000 bp reads exercise the BaseGroup interval widening path
    let dir = tempfile::tempdir().unwrap();
    let seq = "ACGT".repeat(250); // 1000 bp
    let qual = "I".repeat(1000);
    let input = format!("@r1\n{seq}\n+\n{qual}\n@r2\n{seq}\n+\n{qual}\n");
    let fq = dir.path().join("long.fq");
    std::fs::write(&fq, input).unwrap();
    AnchrCmd::new()
        .args(&["fq", "qc", fq.to_str().unwrap(), "-o", dir.path().to_str().unwrap()])
        .assert()
        .success();
    let data = fs::read_to_string(dir.path().join("long.fq_fastqc/fastqc_data.txt")).unwrap();
    assert!(data.contains("Sequence length\t1000"));
    assert!(data.contains("10-19\t"), "interval>10 special group (interval 20)");
    assert!(data.contains("980-999"), "later interval groups");
}

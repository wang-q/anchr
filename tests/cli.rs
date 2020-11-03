use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::env;
use std::process::Command;
use tempfile::TempDir; // Run programs

#[test]
fn command_invalid() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("anchr")?;
    cmd.arg("foobar");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("which wasn't expected"));

    Ok(())
}

#[test]
fn file_doesnt_be_needed() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("anchr")?;
    cmd.arg("test").arg("tests/SKCM/meth.tsv.gz");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("which wasn't expected"));

    Ok(())
}

#[test]
fn file_doesnt_provided() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("anchr")?;
    cmd.arg("dep");
    cmd.assert().failure().stderr(predicate::str::contains(
        "The following required arguments were not provided",
    ));

    Ok(())
}

#[test]
fn file_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("anchr")?;
    cmd.arg("dep").arg("tests/file/doesnt/exist");
    cmd.assert().failure();

    Ok(())
}

#[test]
fn command_dep() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd.arg("dep").arg("check").output().unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(stdout.lines().count() > 100);
    assert!(stdout.contains("fastqc "));

    Ok(())
}

#[test]
fn command_ena() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd.arg("ena").arg("info").output().unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 196);
    assert!(stdout.contains("accession"));

    Ok(())
}

#[test]
fn command_quorum() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("quorum")
        .arg("tests/Lambda/R1.fq.gz")
        .arg("tests/Lambda/R2.fq.gz")
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(stdout.lines().count() > 150);
    assert!(stdout.contains("END_TIME"));

    Ok(())
}

#[test]
fn command_trim() -> Result<(), Box<dyn std::error::Error>> {
    let curdir = env::current_dir().unwrap();

    let tempdir = TempDir::new().unwrap();
    assert!(env::set_current_dir(&tempdir).is_ok());

    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("trim")
        .arg("tests/Lambda/R1.fq.gz")
        .arg("tests/Lambda/R2.fq.gz")
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(stdout.lines().count() > 40);
    assert!(stdout.contains("Sickle"));

    assert!(&tempdir.path().join("illumina_adapters.fa").is_file());
    assert!(&tempdir.path().join("sequencing_artifacts.fa").is_file());

    assert!(env::set_current_dir(&curdir).is_ok());
    assert!(tempdir.close().is_ok());

    Ok(())
}

#[test]
fn command_merge() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("merge")
        .arg("tests/Lambda/R1.fq.gz")
        .arg("tests/Lambda/R2.fq.gz")
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(stdout.lines().count() > 50);
    assert!(stdout.contains("Read merging"));

    Ok(())
}

#[test]
fn command_unitigs() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("unitigs")
        .arg("tests/Lambda/pe.cor.fa.gz")
        .arg("tests/Lambda/env.json")
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(stdout.lines().count() > 50);
    assert!(stdout.contains("create_k_unitigs_large_k"));

    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("unitigs")
        .arg("tests/Lambda/pe.cor.fa.gz")
        .arg("tests/Lambda/env.json")
        .arg("-u")
        .arg("tadpole")
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(stdout.lines().count() > 50);
    assert!(!stdout.contains("create_k_unitigs_large_k"));
    assert!(stdout.contains("tadpole"));

    Ok(())
}

#[test]
fn command_anchors() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("anchors")
        .arg("tests/Lambda/unitigs.fasta")
        .arg("tests/Lambda/R1.fq.gz")
        .arg("-o")
        .arg("stdout")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert!(stdout.lines().count() > 50);
    assert!(stdout.contains("bbwrap.sh"));

    Ok(())
}

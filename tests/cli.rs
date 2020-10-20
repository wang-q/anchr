use assert_cmd::prelude::*; // Add methods on commands
use itertools::Itertools;
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn command_invalid() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("Anchr")?;
    cmd.arg("foobar");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("which wasn't expected"));

    Ok(())
}

#[test]
fn file_doesnt_be_needed() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("Anchr")?;
    cmd.arg("test").arg("tests/SKCM/meth.tsv.gz");
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("which wasn't expected"));

    Ok(())
}

#[test]
fn file_doesnt_provided() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("Anchr")?;
    cmd.arg("dep");
    cmd.assert().failure().stderr(predicate::str::contains(
        "The following required arguments were not provided",
    ));

    Ok(())
}

#[test]
fn file_doesnt_exist() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("Anchr")?;
    cmd.arg("dep").arg("tests/file/doesnt/exist");
    cmd.assert().failure();

    Ok(())
}

#[test]
fn command_dep() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("Anchr")?;
    let output = cmd
        .arg("dep")
        .arg("check")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 151);
    assert!(stdout.contains("fastqc "));

    Ok(())
}

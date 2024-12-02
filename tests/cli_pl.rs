use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn command_overlap() -> anyhow::Result<()> {
    let mut bin = String::new();
    for e in &["LAshow"] {
        if let Ok(pth) = which::which(e) {
            bin = pth.to_string_lossy().to_string();
            break;
        }
    }
    if bin.is_empty() {
        return Ok(());
    } else {
        eprintln!("bin = {:#?}", bin);
    }

    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("overlap")
        .arg("tests/ovlpr/1_4.pac.fasta")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 18);
    assert!(stdout.contains("pac4745_7148"), "original names");
    assert!(stdout.contains("contained"), "relations");

    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("overlap")
        .arg("--idt")
        .arg("0.8")
        .arg("--len")
        .arg("2500")
        .arg("--serial")
        .arg("tests/ovlpr/1_4.pac.fasta")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 4);
    assert!(!stdout.contains("pac4745_7148"), "serials");

    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("overlap")
        .arg("--idt")
        .arg("0.8")
        .arg("--len")
        .arg("2500")
        .arg("--all")
        .arg("tests/ovlpr/1_4.pac.fasta")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 42);

    Ok(())
}

#[test]
fn command_overlap2() -> anyhow::Result<()> {
    let mut bin = String::new();
    for e in &["LAshow"] {
        if let Ok(pth) = which::which(e) {
            bin = pth.to_string_lossy().to_string();
            break;
        }
    }
    if bin.is_empty() {
        return Ok(());
    } else {
        eprintln!("bin = {:#?}", bin);
    }

    let tempdir = tempfile::TempDir::new().unwrap();

    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("overlap2")
        .arg("tests/ovlpr/1_4.anchor.fasta")
        .arg("tests/ovlpr/1_4.pac.fasta")
        .arg("-d")
        .arg(tempdir.path().display().to_string())
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 0);

    assert!(&tempdir.path().join("anchorLong.db").is_file());
    assert!(&tempdir.path().join("anchorLong.ovlp.tsv").is_file());

    assert!(tempdir.close().is_ok());

    Ok(())
}

#[test]
fn command_orient() -> anyhow::Result<()> {
    let mut bin = String::new();
    for e in &["LAshow"] {
        if let Ok(pth) = which::which(e) {
            bin = pth.to_string_lossy().to_string();
            break;
        }
    }
    if bin.is_empty() {
        return Ok(());
    } else {
        eprintln!("bin = {:#?}", bin);
    }

    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("orient")
        .arg("tests/ovlpr/1_4.anchor.fasta")
        .arg("tests/ovlpr/1_4.pac.fasta")
        .arg("-r")
        .arg("tests/ovlpr/1_4.2.restrict.tsv")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 24);
    assert!(stdout.contains("pac4745_7148"), "original names");

    Ok(())
}

#[test]
fn command_contained() -> anyhow::Result<()> {
    let mut bin = String::new();
    for e in &["LAshow"] {
        if let Ok(pth) = which::which(e) {
            bin = pth.to_string_lossy().to_string();
            break;
        }
    }
    if bin.is_empty() {
        return Ok(());
    } else {
        eprintln!("bin = {:#?}", bin);
    }

    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("contained")
        .arg("tests/ovlpr/contained.fasta")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("infile_0/4/0_15361"), "renamed");

    Ok(())
}

#[test]
fn command_merge() -> anyhow::Result<()> {
    let mut bin = String::new();
    for e in &["LAshow"] {
        if let Ok(pth) = which::which(e) {
            bin = pth.to_string_lossy().to_string();
            break;
        }
    }
    if bin.is_empty() {
        return Ok(());
    } else {
        eprintln!("bin = {:#?}", bin);
    }

    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("merge")
        .arg("tests/ovlpr/merge.fasta")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 2);
    assert!(stdout.contains("merge_1"), "renamed");

    Ok(())
}

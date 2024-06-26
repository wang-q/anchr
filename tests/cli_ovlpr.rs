use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs

#[test]
fn command_dazzname() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("dazzname")
        .arg("--start")
        .arg("11")
        .arg("--prefix")
        .arg("seq")
        .arg("tests/ovlpr/1_4.anchor.fasta")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 8);
    assert!(!stdout.contains("anchor148_9124"), "original names");
    assert!(stdout.contains("seq/14/0_9124"), "new names");

    Ok(())
}

#[test]
fn command_show2ovlp() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("show2ovlp")
        .arg("tests/ovlpr/1_4.show.txt")
        .arg("tests/ovlpr/1_4.replace.tsv")
        .arg("--orig")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 50);
    assert!(stdout.contains("pac7556_20928"), "original names");

    Ok(())
}

#[test]
fn command_paf2ovlp() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("paf2ovlp")
        .arg("tests/ovlpr/1_4.pac.paf")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 28);
    assert!(stdout.contains("overlap"), "overlaps");

    Ok(())
}

#[test]
fn command_covered() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("covered")
        .arg("tests/ovlpr/1_4.pac.paf.ovlp.tsv")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 8);
    assert!(stdout.contains("pac4745_7148"), "original names");
    assert!(!stdout.contains("pac4745_7148:1"), "uncovered region");

    Ok(())
}

#[test]
fn command_covered_paf() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("covered")
        .arg("tests/ovlpr/11_2.long.paf")
        .arg("--paf")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 15);
    assert!(stdout.contains("long/13141/0_10011"), "original names");
    assert!(!stdout.contains("long/13141/0_10011:1"), "uncovered region");

    Ok(())
}

#[test]
fn command_covered_longest() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("covered")
        .arg("tests/ovlpr/1_4.pac.paf.ovlp.tsv")
        .arg("--longest")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 8);
    assert!(stdout.contains("pac4745_7148"), "original names");
    assert!(!stdout.contains("pac4745_7148:1"), "uncovered region");

    Ok(())
}

#[test]
fn command_covered_base() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("covered")
        .arg("tests/ovlpr/1_4.pac.paf.ovlp.tsv")
        .arg("--base")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 98105);
    assert!(stdout.contains("pac4745_7148"), "original names");

    Ok(())
}

#[test]
fn command_covered_mean() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("covered")
        .arg("tests/ovlpr/1_4.pac.paf.ovlp.tsv")
        .arg("tests/ovlpr/1_4.pac.paf.ovlp.tsv")
        .arg("tests/ovlpr/1_4.pac.paf.ovlp.tsv")
        .arg("--mean")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 8);
    assert!(stdout.contains("pac4745_7148"), "original names");
    assert!(
        stdout.contains("pac1461_9030\t9030\t2.8"),
        "avoid duplicates"
    );

    Ok(())
}

#[test]
fn command_restrict() -> anyhow::Result<()> {
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("restrict")
        .arg("tests/ovlpr/1_4.ovlp.tsv")
        .arg("tests/ovlpr/1_4.restrict.tsv")
        .output()
        .unwrap();
    let stdout = String::from_utf8(output.stdout).unwrap();

    assert_eq!(stdout.lines().count(), 36);
    assert!(
        !stdout.contains("pac1461_9030\tpac8852_20444"),
        "no long-long overlaps"
    );

    Ok(())
}

use assert_cmd::prelude::*; // Add methods on commands
use std::env;
use std::process::Command;
use tempfile::TempDir; // Run programs

#[test]
fn command_template() -> Result<(), Box<dyn std::error::Error>> {
    let curdir = env::current_dir().unwrap();

    let tempdir = TempDir::new().unwrap();
    assert!(env::set_current_dir(&tempdir).is_ok());

    // Anchr template
    let mut cmd = Command::cargo_bin("Anchr")?;
    let output = cmd.arg("template").output().unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(stderr.lines().count(), 7);
    assert!(stderr.contains("2_trim.sh"));
    assert!(&tempdir.path().join("2_trim.sh").is_file());

    // Anchr template --fastqc
    let mut cmd = Command::cargo_bin("Anchr")?;
    let output = cmd.arg("template").arg("--fastqc").output().unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(stderr.lines().count(), 8);
    assert!(stderr.contains("2_fastqc.sh"));
    assert!(&tempdir.path().join("2_fastqc.sh").is_file());

    // Anchr template --fastqc --merge
    let mut cmd = Command::cargo_bin("Anchr")?;
    let output = cmd
        .arg("template")
        .arg("--fastqc")
        .arg("--merge")
        .output()
        .unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(stderr.lines().count(), 10);
    assert!(stderr.contains("2_merge.sh"));
    assert!(&tempdir.path().join("2_merge.sh").is_file());

    // Anchr template --quorum
    let mut cmd = Command::cargo_bin("Anchr")?;
    let output = cmd.arg("template").arg("--quorum").output().unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(stderr.lines().count(), 7);
    assert!(stderr.contains("2_quorum.sh"));
    assert!(&tempdir.path().join("2_quorum.sh").is_file());

    // cleanup
    assert!(env::set_current_dir(&curdir).is_ok());
    assert!(tempdir.close().is_ok());

    Ok(())
}

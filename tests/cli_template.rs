use assert_cmd::prelude::*; // Add methods on commands
use std::env;
use std::process::Command;
use tempfile::TempDir; // Run programs

#[test]
fn command_template() -> anyhow::Result<()> {
    let curdir = env::current_dir().unwrap();

    let tempdir = TempDir::new().unwrap();
    assert!(env::set_current_dir(&tempdir).is_ok());

    // anchr template
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd.arg("template").output().unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(stderr.lines().count(), 17);
    assert!(stderr.contains("2_trim.sh"));
    assert!(&tempdir.path().join("0_script/2_trim.sh").is_file());

    // anchr template --fastqc
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd.arg("template").arg("--fastqc").output().unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(stderr.lines().count(), 18);
    assert!(stderr.contains("2_fastqc.sh"));
    assert!(&tempdir.path().join("0_script/2_fastqc.sh").is_file());

    // anchr template --fastqc --merge
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("template")
        .arg("--fastqc")
        .arg("--merge")
        .output()
        .unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert_eq!(stderr.lines().count(), 25);
    assert!(stderr.contains("2_merge.sh"));
    assert!(&tempdir.path().join("0_script/2_merge.sh").is_file());

    // anchr template always generates 2_quorum.sh (s-filter replaced the
    // external quorum; the no_quorum fallback was removed).
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd.arg("template").output().unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert!(stderr.contains("2_quorum.sh"));
    assert!(&tempdir.path().join("0_script/2_quorum.sh").is_file());

    // anchr template --unitigger bcalm (no multik script)
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("template")
        .arg("--unitigger")
        .arg("bcalm")
        .output()
        .unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert!(stderr.contains("4_unitigs_bcalm.sh"));
    assert!(!stderr.contains("4_unitigs_multik.sh"));
    assert!(&tempdir.path().join("0_script/4_unitigs_bcalm.sh").is_file());

    let bcalm_script = std::fs::read_to_string(tempdir.path().join("0_script/4_unitigs_bcalm.sh"))?;
    assert!(bcalm_script.contains("bcalm \\"));
    assert!(bcalm_script.contains("anchr asm olc --unitigs unitigs_K*.fasta"));

    // anchr template --unitigger "multik bcalm" --merge
    let mut cmd = Command::cargo_bin("anchr")?;
    let output = cmd
        .arg("template")
        .arg("--unitigger")
        .arg("multik bcalm")
        .arg("--merge")
        .output()
        .unwrap();
    let stderr = String::from_utf8(output.stderr).unwrap();

    assert!(stderr.contains("4_unitigs_multik.sh"));
    assert!(stderr.contains("4_unitigs_bcalm.sh"));
    assert!(stderr.contains("6_unitigs_bcalm.sh"));
    assert!(&tempdir.path().join("0_script/6_unitigs_bcalm.sh").is_file());

    let master = std::fs::read_to_string(tempdir.path().join("0_script/0_master.sh"))?;
    assert!(master.contains("4_unitigs_multik"));
    assert!(master.contains("4_unitigs_bcalm"));
    assert!(master.contains("statUnitigsBcalm.md"));

    // cleanup
    assert!(env::set_current_dir(&curdir).is_ok());
    assert!(tempdir.close().is_ok());

    Ok(())
}

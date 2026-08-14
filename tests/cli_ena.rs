use anyhow::Result;
#[macro_use]
#[path = "common/mod.rs"]
mod common;
use common::AnchrCmd;
use std::fs;

fn write_temp(dir: &std::path::Path, name: &str, content: &str) -> std::path::PathBuf {
    let p = dir.join(name);
    fs::write(&p, content).unwrap();
    p
}

fn sample_json() -> &'static str {
    r#"{
  "grpA": {
    "SRX1": {
      "instrument_platform": "ILLUMINA",
      "library_layout": "PAIRED",
      "nominal_length": "150",
      "srr": ["SRR1", "SRR2"],
      "srr_info": {
        "SRR1": {"read_count": "100", "base_count": "30000"},
        "SRR2": {"read_count": "200", "base_count": "60000"}
      },
      "downloads": [
        "ftp://ftp.sra.ebi.ac.uk/vol1/run/1_1.fq.gz",
        "ftp://ftp.sra.ebi.ac.uk/vol1/run/1_2.fq.gz"
      ],
      "md5s": ["abc 1_1.fq.gz", "def 1_2.fq.gz"]
    }
  },
  "grpB": {
    "SRX2": {
      "instrument_platform": "PACBIO_SMRT",
      "library_layout": "SINGLE",
      "nominal_length": "",
      "srr": ["SRR3"],
      "srr_info": {"SRR3": {"read_count": "50", "base_count": "150000"}},
      "downloads": ["ftp://ftp.sra.ebi.ac.uk/vol1/run/2.fq.gz"],
      "md5s": ["ghi 2.fq.gz"]
    }
  }
}"#
}

#[test]
fn command_ena_manifest_writes_download_files() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let json = write_temp(dir.path(), "samples.json", sample_json());
    AnchrCmd::new()
        .current_dir(dir.path())
        .args(&["ena", "manifest", json.to_str().unwrap()])
        .assert()
        .success();

    let tsv = fs::read_to_string(dir.path().join("samples.tsv"))?;
    assert!(tsv.contains("name\tsrx\tplatform\tlayout\tilength\tsrr\tspots\tbases\n"));
    assert!(tsv.contains("grpA\tSRX1\tILLUMINA\tPAIRED\t150\tSRR1\t100\t30000\n"));
    assert!(tsv.contains("grpA\tSRX1\tILLUMINA\tPAIRED\t150\tSRR2\t200\t60000\n"));

    let ftp = fs::read_to_string(dir.path().join("samples.ftp.txt"))?;
    assert!(ftp.contains("ftp://ftp.sra.ebi.ac.uk/vol1/run/1_1.fq.gz\n"));
    assert!(ftp.contains("ftp://ftp.sra.ebi.ac.uk/vol1/run/2.fq.gz\n"));

    let md5 = fs::read_to_string(dir.path().join("samples.md5.txt"))?;
    assert!(md5.contains("abc 1_1.fq.gz\n"));
    assert!(md5.contains("ghi 2.fq.gz\n"));

    assert!(!dir.path().join("samples.ascp.sh").exists());
    Ok(())
}

#[test]
fn command_ena_manifest_applies_filters() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let json = write_temp(dir.path(), "s.json", sample_json());
    AnchrCmd::new()
        .current_dir(dir.path())
        .args(&[
            "ena",
            "manifest",
            json.to_str().unwrap(),
            "--platform",
            "illumina",
            "--layout",
            "pair",
        ])
        .assert()
        .success();

    let tsv = fs::read_to_string(dir.path().join("s.tsv"))?;
    assert!(tsv.contains("grpA\tSRX1\tILLUMINA\tPAIRED"));
    assert!(!tsv.contains("grpB"));
    Ok(())
}

#[test]
fn command_ena_manifest_writes_ascp_script() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let json = write_temp(dir.path(), "s.json", sample_json());
    AnchrCmd::new()
        .current_dir(dir.path())
        .args(&["ena", "manifest", json.to_str().unwrap(), "--ascp"])
        .assert()
        .success();

    let ascp = fs::read_to_string(dir.path().join("s.ascp.sh"))?;
    assert!(ascp.contains("era-fasp@fasp.sra.ebi.ac.uk:vol1/run/1_1.fq.gz"));
    assert!(ascp.contains("openssl dgst -md5 -r 1_1.fq.gz"));
    assert!(ascp.contains("abc"));
    Ok(())
}

#[test]
fn command_ena_manifest_rejects_invalid_json() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let json = write_temp(dir.path(), "bad.json", "not json");
    AnchrCmd::new()
        .current_dir(dir.path())
        .args(&["ena", "manifest", json.to_str().unwrap()])
        .assert()
        .failure();
    Ok(())
}

#[test]
fn command_ena_meta_requires_infile() {
    AnchrCmd::new().args(&["ena", "meta"]).assert().failure();
}

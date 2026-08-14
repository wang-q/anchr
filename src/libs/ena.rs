//! ENA metadata fetching and download-manifest generation.
//!
//! Replaces the legacy `templates/ena_info.pl` / `ena_prep.pl`: `meta`
//! queries the ENA portal filereport API for run metadata and writes JSON;
//! `manifest` turns that JSON into aria2c/md5 download lists.

use anyhow::{bail, Context, Result};
use serde_json::{json, Value};

/// ENA portal filereport endpoint.
const API: &str = "http://www.ebi.ac.uk/ena/portal/api/filereport";

/// A fetcher for one ENA accession (study/experiment → read_run JSON).
type Fetcher = dyn Fn(&str, bool) -> Result<Value>;

/// Fetches the read_run metadata for one accession via the ENA portal API.
fn fetch_read_run(accession: &str, use_sra: bool) -> Result<Value> {
    let fields = format!(
        "secondary_study_accession,secondary_sample_accession,experiment_accession,\
         run_accession,scientific_name,instrument_platform,instrument_model,\
         library_name,nominal_length,library_layout,library_source,library_selection,\
         read_count,base_count,{}",
        if use_sra {
            "sra_md5,sra_ftp"
        } else {
            "fastq_md5,fastq_ftp"
        }
    );
    let url = format!("{API}?accession={accession}&result=read_run&fields={fields}&format=json");
    let body = ureq::get(&url)
        .call()
        .with_context(|| format!("failed to fetch {url}"))?
        .into_string()?;
    serde_json::from_str(&body).with_context(|| format!("invalid JSON from {url}"))
}

/// Parses the input CSV (SRA accession, group name, ...) and returns the
/// per-group accession list (only `[DES]R\w`, `SAMN`, `PRJNA` ids).
fn parse_csv(lines: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    for line in lines.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let mut cols = trimmed.split(',');
        let Some(key) = cols.next() else { continue };
        if !key.starts_with("SRR")
            && !key.starts_with("ERR")
            && !key.starts_with("DRR")
            && !key.starts_with("SAMN")
            && !key.starts_with("PRJNA")
        {
            continue;
        }
        let name = cols.next().filter(|s| !s.is_empty()).unwrap_or(key);
        out.push((key.to_string(), name.to_string()));
    }
    out
}

/// Builds one experiment's info from a filereport run row (mirrors the
/// fields of the legacy `erx_worker`).
fn run_info(row: &Value, use_sra: bool) -> Value {
    let get = |k: &str| row.get(k).and_then(Value::as_str).unwrap_or("").to_string();
    let mut downloads = Vec::new();
    let mut md5s = Vec::new();
    let ftp = get(if use_sra { "sra_ftp" } else { "fastq_ftp" });
    let md5 = get(if use_sra { "sra_md5" } else { "fastq_md5" });
    let ftp_parts: Vec<&str> = ftp.split(';').filter(|s| !s.is_empty()).collect();
    let md5_parts: Vec<&str> = md5.split(';').filter(|s| !s.is_empty()).collect();
    for (i, path) in ftp_parts.iter().enumerate() {
        downloads.push(format!("ftp://{path}"));
        let base = path.rsplit('/').next().unwrap_or(path);
        let sum = md5_parts.get(i).copied().unwrap_or("");
        md5s.push(format!("{sum} {base}").trim().to_string());
    }
    let srr = get("run_accession");
    let mut srr_info = serde_json::Map::new();
    srr_info.insert(
        srr.clone(),
        json!({
            "read_count": row.get("read_count").and_then(Value::as_str).unwrap_or(""),
            "base_count": row.get("base_count").and_then(Value::as_str).unwrap_or(""),
        }),
    );
    json!({
        "srp": get("secondary_study_accession"),
        "srs": get("secondary_sample_accession"),
        "srx": get("experiment_accession"),
        "scientific_name": get("scientific_name"),
        "instrument_platform": get("instrument_platform"),
        "instrument_model": get("instrument_model"),
        "library_name": get("library_name"),
        "nominal_length": get("nominal_length"),
        "library_layout": get("library_layout"),
        "library_source": get("library_source"),
        "library_selection": get("library_selection"),
        "srr": [srr],
        "srr_info": srr_info,
        "downloads": downloads,
        "md5s": md5s,
    })
}

/// Queries ENA for every CSV row and aggregates the per-group JSON metadata.
fn meta_from_csv_impl(lines: &str, use_sra: bool, fetch: &Fetcher) -> Result<Value> {
    let mut master: serde_json::Map<String, Value> = serde_json::Map::new();
    for (key, name) in parse_csv(lines) {
        let rows = fetch(&key, use_sra)?;
        let runs = match rows {
            Value::Array(a) => a,
            _ => bail!("expected a JSON array from ENA for {key}"),
        };
        // Group runs by experiment accession, merging multiple runs per SRX.
        let mut srx_map: serde_json::Map<String, Value> = serde_json::Map::new();
        for row in &runs {
            let srx = row
                .get("experiment_accession")
                .and_then(Value::as_str)
                .unwrap_or("")
                .to_string();
            if srx.is_empty() {
                continue;
            }
            srx_map
                .entry(srx.clone())
                .and_modify(|info| {
                    if let Some(srr) = row.get("run_accession").and_then(Value::as_str) {
                        if let Some(a) = info["srr"].as_array_mut() {
                            a.push(Value::String(srr.to_string()));
                        }
                        if let Some(map) = info["srr_info"].as_object_mut() {
                            map.insert(
                                srr.to_string(),
                                json!({
                                    "read_count": row.get("read_count").and_then(Value::as_str).unwrap_or(""),
                                    "base_count": row.get("base_count").and_then(Value::as_str).unwrap_or(""),
                                }),
                            );
                        }
                        for (i, path) in row
                            .get(if use_sra { "sra_ftp" } else { "fastq_ftp" })
                            .and_then(Value::as_str)
                            .unwrap_or("")
                            .split(';')
                            .filter(|s| !s.is_empty())
                            .enumerate()
                        {
                            if let Some(dl) = info["downloads"].as_array_mut() {
                                dl.push(Value::String(format!("ftp://{path}")));
                            }
                            let base = path.rsplit('/').next().unwrap_or(path);
                            let sum = row
                                .get(if use_sra { "sra_md5" } else { "fastq_md5" })
                                .and_then(Value::as_str)
                                .unwrap_or("")
                                .split(';')
                                .filter(|s| !s.is_empty())
                                .nth(i)
                                .unwrap_or("");
                            if let Some(m) = info["md5s"].as_array_mut() {
                                m.push(Value::String(format!("{sum} {base}").trim().to_string()));
                            }
                        }
                    }
                })
                .or_insert_with(|| run_info(row, use_sra));
        }
        master.entry(name).or_insert_with(|| json!(srx_map));
    }
    Ok(Value::Object(master))
}

/// Fetch ENA metadata for the CSV input and write the aggregated JSON.
pub fn meta_from_csv(lines: &str, use_sra: bool) -> Result<Value> {
    meta_from_csv_impl(lines, use_sra, &fetch_read_run)
}

/// Download manifest generated from metadata JSON.
#[derive(Debug, Default)]
pub struct Manifest {
    /// Tab-separated run table (`name srx platform layout ilength srr spots bases`).
    pub tsv: String,
    /// aria2c input list (`ftp://...`).
    pub ftp: String,
    /// `md5sum --check` list (`md5 basename`).
    pub md5: String,
    /// Aspera download script (empty unless `ascp` is requested).
    pub ascp: String,
}

/// Builds the download manifest from metadata JSON, applying the optional
/// platform/layout filters and the Aspera script flag.
pub fn manifest_from_json(
    json: &Value,
    platform: Option<&str>,
    layout: Option<&str>,
    ascp: bool,
) -> Result<Manifest> {
    let mut m = Manifest::default();
    let names = json
        .as_object()
        .context("metadata JSON must be an object")?;
    let mut tsv = String::from("name\tsrx\tplatform\tlayout\tilength\tsrr\tspots\tbases\n");
    let mut sorted: Vec<(&String, &Value)> = names.iter().collect();
    sorted.sort_by(|a, b| a.0.cmp(b.0));
    for (name, srx_map) in sorted {
        let srx_obj = srx_map
            .as_object()
            .context(format!("metadata for {name} must be an object"))?;
        let mut srx_keys: Vec<&String> = srx_obj.keys().collect();
        srx_keys.sort();
        for srx in srx_keys {
            let info = &srx_obj[srx];
            let get = |k: &str| {
                info.get(k)
                    .and_then(Value::as_str)
                    .unwrap_or("")
                    .to_string()
            };
            let platform_v = get("instrument_platform");
            let layout_v = get("library_layout");
            if let Some(p) = platform {
                if !platform_v.to_lowercase().contains(&p.to_lowercase()) {
                    continue;
                }
            }
            if let Some(l) = layout {
                if !layout_v.to_lowercase().contains(&l.to_lowercase()) {
                    continue;
                }
            }
            let ilength = get("nominal_length");
            let srrs = info
                .get("srr")
                .and_then(Value::as_array)
                .cloned()
                .unwrap_or_default();
            for srr in srrs {
                let srr = srr.as_str().unwrap_or("");
                let (spots, bases) = info
                    .get("srr_info")
                    .and_then(|si| si.get(srr))
                    .map(|s| {
                        (
                            s.get("read_count")
                                .and_then(Value::as_str)
                                .unwrap_or("")
                                .to_string(),
                            s.get("base_count")
                                .and_then(Value::as_str)
                                .unwrap_or("")
                                .to_string(),
                        )
                    })
                    .unwrap_or_default();
                tsv.push_str(&format!(
                    "{name}\t{srx}\t{platform_v}\t{layout_v}\t{ilength}\t{srr}\t{spots}\t{bases}\n"
                ));
            }
            for dl in info
                .get("downloads")
                .and_then(Value::as_array)
                .unwrap_or(&Vec::new())
            {
                m.ftp.push_str(&format!("{}\n", dl.as_str().unwrap_or("")));
            }
            for md in info
                .get("md5s")
                .and_then(Value::as_array)
                .unwrap_or(&Vec::new())
            {
                m.md5.push_str(&format!("{}\n", md.as_str().unwrap_or("")));
            }
            if ascp {
                let downloads = info
                    .get("downloads")
                    .and_then(Value::as_array)
                    .cloned()
                    .unwrap_or_default();
                let md5s = info
                    .get("md5s")
                    .and_then(Value::as_array)
                    .cloned()
                    .unwrap_or_default();
                for (i, dl) in downloads.iter().enumerate() {
                    let url = dl.as_str().unwrap_or("");
                    let fn_name = url.rsplit('/').next().unwrap_or(url);
                    let md = md5s
                        .get(i)
                        .and_then(Value::as_str)
                        .unwrap_or("")
                        .to_string();
                    let md5_sum = md.split(' ').next().unwrap_or("").to_string();
                    let ascp_url = url.replacen(
                        "ftp://ftp.sra.ebi.ac.uk/",
                        "era-fasp@fasp.sra.ebi.ac.uk:",
                        1,
                    );
                    m.ascp.push_str(&format!(
                        "[ ! -e {fn_name} ] && $HOME/.aspera/connect/bin/ascp \
                         -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh \
                         -TQ -k1 -v -P33001 {ascp_url} . && \
                         if [ $(openssl dgst -md5 -r {fn_name} | cut -d' ' -f 1) != {md5_sum} ]; then \
                         echo -e '{fn_name}\\tNot OK'; rm {fn_name}; \
                         else echo -e '{fn_name}\\tOK'; fi\n"
                    ));
                }
            }
        }
    }
    m.tsv = tsv;
    Ok(m)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn csv_parsing_skips_headers_and_non_sra_ids() {
        let csv = "# comment\nSRR1,groupA\nERR2,groupB\nfoo,bar\n\nDRR3\n";
        let rows = parse_csv(csv);
        assert_eq!(
            rows,
            vec![
                ("SRR1".to_string(), "groupA".to_string()),
                ("ERR2".to_string(), "groupB".to_string()),
                ("DRR3".to_string(), "DRR3".to_string()),
            ]
        );
    }

    #[test]
    fn meta_aggregates_runs_per_experiment() {
        let fetch = |acc: &str, _sra: bool| {
            Ok(json!([
                {
                    "experiment_accession": format!("SRX{acc}"),
                    "run_accession": format!("SRR{acc}1"),
                    "secondary_study_accession": "SRP1",
                    "secondary_sample_accession": "SRS1",
                    "scientific_name": "E. coli",
                    "instrument_platform": "ILLUMINA",
                    "library_layout": "PAIRED",
                    "nominal_length": "150",
                    "read_count": "100",
                    "base_count": "30000",
                    "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/run/1.fq.gz",
                    "fastq_md5": "abc",
                }
            ]))
        };
        let meta = meta_from_csv_impl("SRR1,grp\n", false, &fetch).unwrap();
        let grp = &meta["grp"]["SRXSRR1"];
        assert_eq!(grp["scientific_name"], "E. coli");
        assert_eq!(grp["srr"], json!(["SRRSRR11"]));
        assert_eq!(
            grp["downloads"][0],
            "ftp://ftp.sra.ebi.ac.uk/vol1/run/1.fq.gz"
        );
        assert_eq!(grp["md5s"][0], "abc 1.fq.gz");
    }

    #[test]
    fn manifest_writes_tsv_ftp_md5_and_ascp() {
        let meta = json!({
            "grp": {
                "SRX1": {
                    "instrument_platform": "ILLUMINA",
                    "library_layout": "PAIRED",
                    "nominal_length": "150",
                    "srr": ["SRR1"],
                    "srr_info": {"SRR1": {"read_count": "100", "base_count": "30000"}},
                    "downloads": ["ftp://ftp.sra.ebi.ac.uk/vol1/run/1.fq.gz"],
                    "md5s": ["abc 1.fq.gz"],
                }
            }
        });
        let m = manifest_from_json(&meta, None, None, true).unwrap();
        assert!(m
            .tsv
            .contains("grp\tSRX1\tILLUMINA\tPAIRED\t150\tSRR1\t100\t30000\n"));
        assert!(m.ftp.contains("ftp://ftp.sra.ebi.ac.uk/vol1/run/1.fq.gz\n"));
        assert!(m.md5.contains("abc 1.fq.gz\n"));
        assert!(m.ascp.contains("ascp"));
        assert!(m
            .ascp
            .contains("era-fasp@fasp.sra.ebi.ac.uk:vol1/run/1.fq.gz"));
    }

    #[test]
    fn manifest_applies_platform_filter() {
        let meta = json!({
            "grp": {
                "SRX1": {
                    "instrument_platform": "PACBIO_SMRT",
                    "library_layout": "SINGLE",
                    "nominal_length": "",
                    "srr": ["SRR1"],
                    "srr_info": {"SRR1": {"read_count": "1", "base_count": "1"}},
                    "downloads": ["ftp://x/1.fq.gz"],
                    "md5s": ["a 1.fq.gz"],
                }
            }
        });
        let m = manifest_from_json(&meta, Some("illumina"), None, false).unwrap();
        assert!(m
            .tsv
            .ends_with("name\tsrx\tplatform\tlayout\tilength\tsrr\tspots\tbases\n"));
    }
}

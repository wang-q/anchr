//! `anchr fq qc` — FastQC-compatible reads quality-control statistics.
//!
//! M1 (design: notes/design/qc.md): statistical modules with
//! `fastqc_data.txt` / `summary.txt` output; per-file subdirectories named
//! `{input}_fastqc/` (fastqc layout).

use crate::libs::qc::analyzer::QcStats;
use anyhow::Context;
use clap::{Arg, ArgMatches, Command};
use pgr::libs::fmt::seq::{SeqReader, SeqRecord};
use pgr::libs::fq::qual::detect_quality_base;
use rayon::prelude::*;
use std::fs;
use std::io::Write;
use std::path::Path;

/// Build the clap subcommand for qc.
pub fn make_subcommand() -> Command {
    Command::new("qc")
        .about("Computes reads quality-control statistics (FastQC-compatible)")
        .after_help(
            r###"
Reads FASTQ file(s) and writes FastQC-compatible quality-control statistics:

* `fastqc_data.txt` — per-module data (Basic Statistics, per-base/per-sequence
  quality, base content, GC content, N content, length distribution).
* `summary.txt` — per-module pass/warn/fail.

Output goes to `<outdir>/<input>_fastqc/` for each input file (fastqc layout).
The HTML report is planned for a later milestone.

Examples:
1. Basic QC of one file:
   anchr fq qc R1.fq.gz -o qc_out
2. From stdin:
   anchr fq qc stdin -o qc_out
"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Input FASTQ file(s)"),
        )
        .arg(
            Arg::new("outdir")
                .long("outdir")
                .short('o')
                .num_args(1)
                .default_value(".")
                .help("Output directory (created if missing)"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("auto")
                .help("Worker threads (default: logical CPU count)"),
        )
        .arg(
            Arg::new("format")
                .long("format")
                .short('f')
                .num_args(1)
                .default_value("txt")
                .value_parser(["txt", "html", "both"])
                .help("Output format: txt (data + summary), html (report), both"),
        )
}

/// Execute the qc command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let outdir = args.get_one::<String>("outdir").unwrap();
    fs::create_dir_all(outdir).with_context(|| format!("failed to create {outdir}"))?;
    let parallel = crate::cmd::args::parse_parallel_auto(
        args.get_one::<String>("parallel").unwrap(),
    )?;
    let format = args.get_one::<String>("format").unwrap();

    for infile in args.get_many::<String>("infiles").unwrap() {
        let name = if infile == "stdin" {
            "stdin".to_string()
        } else {
            Path::new(infile)
                .file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| infile.clone())
        };
        run_one(infile, &name, outdir, parallel, format)?;
    }
    Ok(())
}

fn run_one(
    infile: &str,
    name: &str,
    outdir: &str,
    parallel: usize,
    format: &str,
) -> anyhow::Result<()> {
    let mut reader =
        SeqReader::new(infile).with_context(|| format!("failed to open {infile}"))?;
    let mut rec = SeqRecord::new();

    // Sample up to 200 reads for Phred encoding detection (pgr).
    let mut sample: Vec<SeqRecord> = Vec::new();
    while sample.len() < 200 && reader.read_record(&mut rec)? {
        sample.push(rec.clone());
    }
    let offset = detect_quality_base(&sample);

    let mut reads = sample;
    while reader.read_record(&mut rec)? {
        reads.push(rec.clone());
    }

    let stats = if parallel > 1 && reads.len() > parallel {
        let chunk = reads.len().div_ceil(parallel);
        reads
            .par_chunks(chunk)
            .enumerate()
            .map(|(ci, block)| {
                let mut s = QcStats::new(name, offset);
                for (li, r) in block.iter().enumerate() {
                    s.consume(r, (ci * chunk + li) as u64);
                }
                s
            })
            .reduce(
                || QcStats::new(name, offset),
                |mut a, b| {
                    a.merge(&b);
                    a
                },
            )
    } else {
        let mut s = QcStats::new(name, offset);
        for (i, r) in reads.iter().enumerate() {
            s.consume(r, i as u64);
        }
        s
    };

    if stats.n_reads() == 0 {
        anyhow::bail!("no reads found in {infile}");
    }

    let dir = Path::new(outdir).join(format!("{name}_fastqc"));
    fs::create_dir_all(&dir).with_context(|| format!("failed to create {}", dir.display()))?;

    if format == "html" || format == "both" {
        let html = stats.report_html()?;
        let html_path = dir.join("fastqc_report.html");
        let mut wh = pgr::libs::io::writer(html_path.to_str().unwrap())?;
        wh.write_all(html.as_bytes())?;
        wh.flush()?;
    }
    if format == "txt" || format == "both" {
        let data_path = dir.join("fastqc_data.txt");
        let mut w = pgr::libs::io::writer(data_path.to_str().unwrap())?;
        stats.report_txt(&mut w)?;
        w.flush()?;

        let sum_path = dir.join("summary.txt");
        let mut ws = pgr::libs::io::writer(sum_path.to_str().unwrap())?;
        stats.report_summary(&mut ws)?;
        ws.flush()?;
    }

    Ok(())
}

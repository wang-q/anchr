use anyhow::Context;
use clap::{Arg, ArgAction, ArgMatches, Command};
use std::io::Write;

use crate::libs::ena;

/// Builds the clap subcommand for ena.
pub fn make_subcommand() -> Command {
    Command::new("ena")
        .about("Queries ENA metadata and prepares download lists")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(meta::make_subcommand())
        .subcommand(manifest::make_subcommand())
}

/// Execute the ena command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    match args.subcommand() {
        Some(("meta", sub_matches)) => meta::execute(sub_matches),
        Some(("manifest", sub_matches)) => manifest::execute(sub_matches),
        _ => unreachable!(),
    }
}

mod meta {
    use super::*;

    /// Builds the clap subcommand for ena meta.
    pub fn make_subcommand() -> Command {
        Command::new("meta")
            .about("Fetches ENA run metadata into JSON")
            .after_help(
                r###"
Queries the ENA portal filereport API for every accession in the input CSV
and writes per-group run metadata as JSON. The first column is an ENA object
id (`SRR|ERR|DRR|SAMN|PRJNA`...), the second is the group name (defaults to
the id); remaining columns are ignored.

Notes:
* The JSON is the input for `anchr ena manifest`
* Each group maps to per-experiment (`SRX`) entries with run lists, download
  ftp URLs, and md5 checksums

Examples:
1. Fetch and inspect:
   anchr ena meta samples.csv -o samples.json

2. Prefer SRA archives over FASTQ files:
   anchr ena meta samples.csv --sra -o samples.json
"###,
            )
            .arg(crate::cmd::args::infile_arg_required_with_help(
                "Input CSV file (stdin for standard input)",
            ))
            .arg(crate::cmd::args::outfile_arg())
            .arg(
                Arg::new("sra")
                    .long("sra")
                    .action(ArgAction::SetTrue)
                    .help("Download SRA archives instead of FASTQ files"),
            )
    }

    /// Execute the ena meta command.
    pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
        let infile = args.get_one::<String>("infile").unwrap();
        let outfile = crate::cmd::args::get_outfile(args);
        let use_sra = args.get_flag("sra");
        let csv = if infile == "stdin" {
            std::io::read_to_string(std::io::stdin()).context("failed to read stdin")?
        } else {
            std::fs::read_to_string(infile).with_context(|| format!("failed to read {infile}"))?
        };
        let meta = ena::meta_from_csv(&csv, use_sra)?;
        let mut w =
            pgr::writer(outfile).with_context(|| format!("failed to open output {outfile}"))?;
        writeln!(w, "{}", serde_json::to_string_pretty(&meta)?)?;
        w.flush()?;
        Ok(())
    }
}

mod manifest {
    use super::*;

    /// Builds the clap subcommand for ena manifest.
    pub fn make_subcommand() -> Command {
        Command::new("manifest")
            .about("Builds download lists from metadata JSON")
            .after_help(
                r###"
Reads the JSON written by `anchr ena meta` and writes, next to the input
(basename of the JSON), the download artifacts:
* `<base>.tsv`    - tab-separated run table (name, srx, platform, layout, ...)
* `<base>.ftp.txt`- aria2c input list
* `<base>.md5.txt`- `md5sum --check` list
* `<base>.ascp.sh`- Aspera download script (only with --ascp)

Examples:
1. Generate ftp/md5 lists:
   anchr ena manifest samples.json

2. Restrict to one platform/layout:
   anchr ena manifest samples.json --platform illumina --layout pair

3. Also write an Aspera script:
   anchr ena manifest samples.json --ascp
"###,
            )
            .arg(crate::cmd::args::infile_arg_required_with_help(
                "Input JSON file written by ena meta",
            ))
            .arg(
                Arg::new("platform")
                    .long("platform")
                    .short('p')
                    .num_args(1)
                    .help("Keep runs whose platform matches (case-insensitive)"),
            )
            .arg(
                Arg::new("layout")
                    .long("layout")
                    .short('l')
                    .num_args(1)
                    .help("Keep runs whose layout matches (case-insensitive)"),
            )
            .arg(
                Arg::new("ascp")
                    .long("ascp")
                    .action(ArgAction::SetTrue)
                    .help("Also write the Aspera download script"),
            )
    }

    /// Execute the ena manifest command.
    pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
        let infile = args.get_one::<String>("infile").unwrap();
        let platform = args.get_one::<String>("platform").map(String::as_str);
        let layout = args.get_one::<String>("layout").map(String::as_str);
        let ascp = args.get_flag("ascp");
        let text =
            std::fs::read_to_string(infile).with_context(|| format!("failed to read {infile}"))?;
        let json: serde_json::Value =
            serde_json::from_str(&text).with_context(|| format!("invalid JSON in {infile}"))?;
        let m = ena::manifest_from_json(&json, platform, layout, ascp)?;
        let base = std::path::Path::new(infile)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("ena");
        write_file(&format!("{base}.tsv"), &m.tsv)?;
        write_file(&format!("{base}.ftp.txt"), &m.ftp)?;
        write_file(&format!("{base}.md5.txt"), &m.md5)?;
        if ascp {
            write_file(&format!("{base}.ascp.sh"), &m.ascp)?;
        }
        Ok(())
    }

    fn write_file(path: &str, content: &str) -> anyhow::Result<()> {
        std::fs::write(path, content).with_context(|| format!("failed to write {path}"))?;
        Ok(())
    }
}

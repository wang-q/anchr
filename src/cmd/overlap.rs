use clap::*;
use cmd_lib::*;
use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufRead, Write};
use std::path::Path;
use std::{env, fs};
use tempfile::{Builder, TempDir};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("overlap")
        .about("Detect overlaps by daligner")
        .after_help(
            r###"
This command is for small files.
All operations are running in a tempdir and no intermediate files are kept.

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}


// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let outfile = args.get_one::<String>("outfile").unwrap();

    let curdir = env::current_dir()?;
    let anchr = env::current_exe().unwrap().display().to_string();
    let tempdir = Builder::new().prefix("ovlp_").tempdir()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    run_cmd!(info "==> Paths")?;
    run_cmd!(info "    \"anchr\"   = ${anchr}")?;
    run_cmd!(info "    \"curdir\"  = ${curdir}")?;
    run_cmd!(info "    \"tempdir\" = ${tempdir_str}")?;

    //----------------------------
    // Operating
    //----------------------------
    run_cmd!(info "==> Absolute paths")?;
    // basename => abs_path
    let mut abs_infiles = vec![];
    for infile in args.get_many::<String>("infiles").unwrap() {
        if infile == "stdin" {
            abs_infiles.push("stdin".to_string());
        } else {
            let absolute = intspan::absolute_path(infile)
                .unwrap()
                .display()
                .to_string();

            abs_infiles.push(absolute.to_string());
        }
    }

    run_cmd!(info "==> Switch to tempdir")?;
    env::set_current_dir(tempdir_str)?;

    Ok(())
}

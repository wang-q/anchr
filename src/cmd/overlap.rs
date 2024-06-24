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
* This command is for small files.
* All operations are running in a tempdir and no intermediate files are kept.

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
            Arg::new("len")
                .long("len")
                .short('l')
                .num_args(1)
                .default_value("500")
                .value_parser(value_parser!(i32))
                .help("minimal length of overlaps"),
        )
        .arg(
            Arg::new("idt")
                .long("idt")
                .short('i')
                .num_args(1)
                .default_value("0.7")
                .value_parser(value_parser!(f32))
                .help("minimal identities of overlaps"),
        )
        .arg(
            Arg::new("serial")
                .long("serial")
                .action(ArgAction::SetTrue)
                .help("Serials instead of original names in outputs"),
        )
        .arg(
            Arg::new("all")
                .long("all")
                .action(ArgAction::SetTrue)
                .help("All overlaps instead of proper ones"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("8")
                .value_parser(value_parser!(i32))
                .help("Number of threads"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
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
    let min_len = *args.get_one::<i32>("len").unwrap();
    let min_idt = *args.get_one::<f32>("idt").unwrap();

    let is_serial = args.get_flag("serial");
    let is_all = args.get_flag("all");

    let parallel = *args.get_one::<i32>("parallel").unwrap();

    let outfile = args.get_one::<String>("outfile").unwrap();

    let curdir = env::current_dir()?;
    let anchr = env::current_exe().unwrap().display().to_string();
    let tempdir = Builder::new().prefix("anchr_ovlp_").tempdir()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    run_cmd!(info "==> Paths")?;
    run_cmd!(info "    \"anchr\"   = ${anchr}")?;
    run_cmd!(info "    \"curdir\"  = ${curdir}")?;
    run_cmd!(info "    \"tempdir\" = ${tempdir_str}")?;

    //----------------------------
    // Paths
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
    let abs_outfile = if outfile == "stdout" {
        outfile.to_string()
    } else {
        intspan::absolute_path(outfile)
            .unwrap()
            .display()
            .to_string()
    };

    //----------------------------
    // Operating
    //----------------------------
    let basename = "anchr_ovlp";

    run_cmd!(info "==> Switch to tempdir")?;
    env::set_current_dir(tempdir_str)?;

    run_cmd!(info "==> Preprocess reads to format them for dazzler")?;
    run_cmd!(
        ${anchr} dazzname $[abs_infiles] -o renamed.fasta
    )?;

    run_cmd!(info "==> Make the dazzler DB, each block is of size 50 MB")?;
    run_cmd!(
        fasta2DB ${basename} renamed.fasta;
        DBdust ${basename};
        DBsplit -s50 ${basename};
    )?;

    run_cmd!(info "==> Run daligner")?;
    let block_number = run_fun!(
        cat ${basename}.db |
            perl -n -e r#"/^blocks\s+=\s+(\d+)\s*$/ and print $1"#
    )
    .unwrap()
    .parse::<i32>()
    .unwrap();
    // eprintln!("block_number = {:#?}", block_number);

    run_cmd!(
        HPC.daligner ${basename} -M16 -T${parallel} -e${min_idt} -l${min_len} -s${min_len} -mdust |
            sed "s/ -vS / -S /" |
            bash
    )?;
    if block_number > 1 {
        run_cmd!(
            LAcat ${basename}.@.las > ${basename}.las
        )?;
    }

    run_cmd!(info "==> Outputs")?;
    if is_all {
        run_cmd!(
            LAshow ${basename}.db ${basename}.las > show.txt
        )?;
    } else {
        run_cmd!(
            LAshow -o ${basename}.db ${basename}.las > show.txt
        )?;
    }

    if is_serial {
        run_cmd!(
            ${anchr} show2ovlp show.txt renamed.fasta.replace.tsv -o ${abs_outfile}
        )?;
    } else {
        run_cmd!(
            ${anchr} show2ovlp show.txt renamed.fasta.replace.tsv --orig -o ${abs_outfile}
        )?;
    }

    //----------------------------
    // Done
    //----------------------------
    env::set_current_dir(&curdir)?;

    Ok(())
}

// use std::io::Read;
// fn pause() {
//     let mut stdin = std::io::stdin();
//     let mut stdout = std::io::stdout();
//
//     // We want the cursor to stay at the end of the line, so we print without a newline and flush manually.
//     write!(stdout, "Press any key to continue...").unwrap();
//     stdout.flush().unwrap();
//
//     // Read a single byte and discard
//     let _ = stdin.read(&mut [0u8]).unwrap();
// }

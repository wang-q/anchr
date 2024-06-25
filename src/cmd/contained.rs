use clap::*;
use cmd_lib::*;
use std::collections::BTreeSet;
use std::env;
use tempfile::Builder;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("contained")
        .about("Discard contained unitigs")
        .after_help(
            r###"
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
                .help("Minimal length of overlaps"),
        )
        .arg(
            Arg::new("idt")
                .long("idt")
                .short('i')
                .num_args(1)
                .default_value("0.98")
                .value_parser(value_parser!(f32))
                .help("Minimal identities of overlaps"),
        )
        .arg(
            Arg::new("ratio")
                .long("ratio")
                .short('r')
                .num_args(1)
                .default_value("0.98")
                .value_parser(value_parser!(f32))
                .help("Ratio of being contained"),
        )
        .arg(
            Arg::new("prefix")
                .long("prefix")
                .num_args(1)
                .default_value("infile")
                .help("Prefix of record names"),
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
    let min_ratio = *args.get_one::<f32>("ratio").unwrap();

    let prefix = args.get_one::<String>("prefix").unwrap();

    let parallel = *args.get_one::<i32>("parallel").unwrap();

    let outfile = args.get_one::<String>("outfile").unwrap();

    let curdir = env::current_dir()?;
    let anchr = env::current_exe().unwrap().display().to_string();
    let tempdir = Builder::new().prefix("anchr_contained_").tempdir()?;
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
    run_cmd!(info "==> Switch to tempdir")?;
    env::set_current_dir(tempdir_str)?;

    run_cmd!(info "==> Filter short unitigs then rename reads as there are duplicated names")?;
    let mut infiles = vec![];
    for (i, infile) in abs_infiles.iter().enumerate() {
        run_cmd!(
            faops filter -a ${min_len} -l 0 ${infile} stdout |
                faops dazz -p ${prefix}_${i} stdin infile.${i}.fasta
        )?;
        infiles.push(format!("infile.{}.fasta", i));
    }

    run_cmd!(info "==> `anchr overlap`")?;
    let files = infiles.clone();
    run_cmd!(
        ${anchr} overlap $[files] --len ${min_len} --idt ${min_idt} --parallel ${parallel} -o contained.ovlp.tsv
    )?;

    run_cmd!(info "==> Discard contained unitigs")?;
    let mut discards = vec![];
    let mut seen = BTreeSet::new();
    for line in &intspan::read_lines("contained.ovlp.tsv") {
        let ovlp = anchr::Overlap::new(line);
        if ovlp.is_empty() {
            continue;
        }

        // ignore self overlapping
        if ovlp.f_id == ovlp.g_id {
            continue;
        }

        // ignore poor overlaps
        if ovlp.len < min_len {
            continue;
        }
        if ovlp.idt < min_idt {
            continue;
        }

        // skip duplicated overlaps, i.e., f -> g and g -> f
        let tup = (
            ovlp.f_id.to_string().min(ovlp.g_id.to_string()),
            ovlp.f_id.to_string().max(ovlp.g_id.to_string()),
        );
        // If the set did not have this value present, true is returned.
        let not_seen = seen.insert(tup);
        if !not_seen {
            continue;
        }

        // discard contained unitigs
        if ovlp.contained == "contains" {
            discards.push(ovlp.g_id.to_string());
            continue;
        }
        if ovlp.contained == "contained" {
            discards.push(ovlp.f_id.to_string());
            continue;
        }

        // discard nearly contained unitigs
        let f_r = ovlp.len as f32 / ovlp.f_len as f32;
        let g_r = ovlp.len as f32 / ovlp.g_len as f32;

        if f_r <= min_ratio && g_r > min_ratio {
            discards.push(ovlp.g_id.to_string());
            continue;
        }
        if g_r <= min_ratio && f_r > min_ratio {
            discards.push(ovlp.f_id.to_string());
            continue;
        }
        if f_r > min_ratio && g_r > min_ratio {
            if ovlp.f_len >= ovlp.g_len {
                discards.push(ovlp.g_id.to_string());
            } else {
                discards.push(ovlp.f_id.to_string());
            }
            continue;
        }
    }
    discards.sort_unstable();
    discards.dedup();
    intspan::write_lines(
        "discard.list",
        &discards.iter().map(AsRef::as_ref).collect(),
    )?;
    // pause();

    run_cmd!(info "==> Outputs")?;
    let files = infiles.clone();
    run_cmd!(
        cat $[files] |
            faops some -i -l 0 stdin discard.list ${abs_outfile}
    )?;

    //----------------------------
    // Done
    //----------------------------
    env::set_current_dir(&curdir)?;

    Ok(())
}

// use std::io::{Read, Write};
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

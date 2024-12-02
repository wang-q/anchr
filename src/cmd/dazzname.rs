use bio::io::fasta;
use clap::*;
use std::io::{BufRead, Write};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("dazzname")
        .about("Rename FASTA records for dazz_db")
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("prefix")
                .long("prefix")
                .num_args(1)
                .default_value("read")
                .help("Prefix of record names"),
        )
        .arg(
            Arg::new("start")
                .long("start")
                .value_parser(value_parser!(usize))
                .num_args(1)
                .default_value("1")
                .help("Starting index"),
        )
        .arg(
            Arg::new("no-replace")
                .long("no-replace")
                .action(ArgAction::SetTrue)
                .help("Do not write a .replace.tsv"),
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
    let is_no_replace = args.get_flag("no-replace");

    let outfile = args.get_one::<String>("outfile").unwrap();
    let mut fa_out = fasta::Writer::new(intspan::writer(outfile));
    let mut writer_rplc = intspan::writer(&format!("{}.replace.tsv", outfile));

    let opt_prefix = args.get_one::<String>("prefix").unwrap();
    let mut opt_start = *args.get_one::<usize>("start").unwrap();

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = intspan::reader(infile);

        let fa_in = fasta::Reader::new(reader);
        for result in fa_in.records() {
            // obtain record or fail with error
            let record = result.unwrap();

            if record.is_empty() {
                continue;
            }

            let name = record.id().to_string();
            let length = record.seq().len();
            let serial = opt_start;

            let name_new = format!("{}/{}/0_{}", opt_prefix, serial, length);
            let record_new = fasta::Record::with_attrs(&name_new, None, record.seq());

            fa_out
                .write_record(&record_new)
                .expect("Write fasta file failed");

            if !is_no_replace {
                writer_rplc.write_fmt(format_args!("{}\t{}\n", name_new, name))?;
            }

            opt_start += 1;
        }
    }

    Ok(())
}

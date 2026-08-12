use clap::*;

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
    let mut fa_out = pgr::libs::fmt::fa::writer(outfile)?;

    let opt_prefix = args.get_one::<String>("prefix").unwrap();
    let mut opt_start = *args.get_one::<usize>("start").unwrap();

    let mut rplc_lines = vec![];

    //----------------------------
    // Ops
    //----------------------------
    for infile in args.get_many::<String>("infiles").unwrap() {
        let mut reader = pgr::libs::fmt::fa::reader(infile)?;
        let mut rec = pgr::libs::fmt::seq::SeqRecord::new();
        while reader.read_record(&mut rec)? {
            if rec.sequence().is_empty() {
                continue;
            }

            let name = rec.name().to_string();
            let length = rec.sequence().len();
            let serial = opt_start;

            let name_new = format!("{}/{}/0_{}", opt_prefix, serial, length);
            let record_new = pgr::libs::fmt::fa::FastaRecord::new(&name_new, rec.sequence());

            fa_out
                .write_record(&record_new)
                .expect("Write fasta file failed");

            if !is_no_replace {
                rplc_lines.push(format!("{}\t{}", name_new, name));
            }

            opt_start += 1;
        }
    }

    if !is_no_replace {
        anchr::utils::write_lines(
            &format!("{}.replace.tsv", outfile),
            &rplc_lines.iter().map(AsRef::as_ref).collect::<Vec<&str>>(),
        )?;
    }

    Ok(())
}

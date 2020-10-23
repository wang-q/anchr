use clap::*;
use itertools::Itertools;
use std::collections::HashMap;
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name("merge")
        .about("Merge Illumina PE reads with bbtools")
        .after_help(
            r#"
<R1> [R2] [Rs]

Fastq files can be gzipped
R1 and R2 are paired; or R1 is interleaved
Rs is single
"#,
        )
        .arg(
            Arg::with_name("infiles")
                .help("Sets the input file to use")
                .required(true)
                .min_values(1)
                .index(1),
        )
        .arg(
            Arg::with_name("len")
                .long("len")
                .short("l")
                .help("Filter reads less or equal to this length")
                .takes_value(true)
                .default_value("60")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("qual")
                .long("qual")
                .short("q")
                .help("Quality score for 3' end")
                .takes_value(true)
                .default_value("15")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("prefilter")
                .long("prefilter")
                .help("Prefilter=N (1 or 2) for tadpole and bbmerge")
                .takes_value(true)
                .empty_values(false),
        )
        .arg(
            Arg::with_name("ecphase")
                .long("ecphase")
                .help("Error-correct phases. Phase 2 can be skipped")
                .takes_value(true)
                .default_value("1 2 3")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("prefixm")
                .long("prefixm")
                .help("Prefix of merged reads")
                .takes_value(true)
                .default_value("M")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("prefixu")
                .long("prefixu")
                .help("Prefix of unmerged reads")
                .takes_value(true)
                .default_value("U")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("xmx")
                .long("xmx")
                .help("Set Java memory usage")
                .takes_value(true)
                .empty_values(false),
        )
        .arg(
            Arg::with_name("parallel")
                .long("parallel")
                .short("p")
                .help("Number of threads")
                .takes_value(true)
                .default_value("8")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("outfile")
                .long("outfile")
                .short("o")
                .help("Output filename. [stdout] for screen")
                .takes_value(true)
                .default_value("merge.sh")
                .empty_values(false),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> std::result::Result<(), std::io::Error> {
    let mut writer = intspan::writer(args.value_of("outfile").unwrap());

    // context from args
    let mut opt = HashMap::new();
    opt.insert("len", args.value_of("len").unwrap());
    opt.insert("qual", args.value_of("qual").unwrap());
    opt.insert(
        "prefilter",
        if args.is_present("prefilter") {
            args.value_of("prefilter").unwrap()
        } else {
            "0"
        },
    );
    opt.insert("ecphase", args.value_of("ecphase").unwrap());
    opt.insert("prefixm", args.value_of("prefixm").unwrap());
    opt.insert("prefixu", args.value_of("prefixu").unwrap());
    opt.insert(
        "xmx",
        if args.is_present("xmx") {
            args.value_of("xmx").unwrap()
        } else {
            "0"
        },
    );
    opt.insert("parallel", args.value_of("parallel").unwrap());

    let infiles = args.values_of("infiles").unwrap().collect_vec();

    let mut context = Context::new();
    context.insert("opt", &opt);
    context.insert("args", &infiles);

    // eprintln!("{:#?}", context);

    // many templates
    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("merge", include_str!("../../templates/merge.tera.sh")),
    ])
    .unwrap();

    // eprintln!("{:#?}", tera);

    let rendered = tera.render("merge", &context).unwrap();

    writer.write_all(rendered.as_ref())?;

    Ok(())
}

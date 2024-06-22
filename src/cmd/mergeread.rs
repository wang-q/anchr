use clap::*;
use itertools::Itertools;
use std::collections::HashMap;
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("mergeread")
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
            Arg::new("infiles")
                .help("Sets the input file to use")
                .required(true)
                .num_args(1..)
                .index(1),
        )
        .arg(
            Arg::new("len")
                .long("len")
                .short('l')
                .help("Filter reads less or equal to this length")
                .num_args(1)
                .default_value("60"),
        )
        .arg(
            Arg::new("qual")
                .long("qual")
                .short('q')
                .help("Quality score for 3' end")
                .num_args(1)
                .default_value("15"),
        )
        .arg(
            Arg::new("prefilter")
                .long("prefilter")
                .help("Prefilter=N (1 or 2) for tadpole and bbmerge")
                .num_args(1),
        )
        .arg(
            Arg::new("ecphase")
                .long("ecphase")
                .help("Error-correct phases. Phase 2 can be skipped")
                .num_args(1)
                .default_value("1 2 3"),
        )
        .arg(
            Arg::new("prefixm")
                .long("prefixm")
                .help("Prefix of merged reads")
                .num_args(1)
                .default_value("M"),
        )
        .arg(
            Arg::new("prefixu")
                .long("prefixu")
                .help("Prefix of unmerged reads")
                .num_args(1)
                .default_value("U"),
        )
        .arg(
            Arg::new("xmx")
                .long("xmx")
                .help("Set Java memory usage")
                .num_args(1),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .help("Number of threads")
                .num_args(1)
                .default_value("8"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .help("Output filename. [stdout] for screen")
                .num_args(1)
                .default_value("merge.sh"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    // context from args
    let mut opt = HashMap::new();
    opt.insert("len", args.get_one::<String>("len").unwrap());
    opt.insert("qual", args.get_one::<String>("qual").unwrap());
    opt.insert("ecphase", args.get_one::<String>("ecphase").unwrap());
    opt.insert("prefixm", args.get_one::<String>("prefixm").unwrap());
    opt.insert("prefixu", args.get_one::<String>("prefixu").unwrap());
    opt.insert("parallel", args.get_one::<String>("parallel").unwrap());

    let binding_0 = "0".to_string();
    opt.insert(
        "prefilter",
        if args.contains_id("prefilter") {
            args.get_one::<String>("prefilter").unwrap()
        } else {
            &binding_0
        },
    );

    opt.insert(
        "xmx",
        if args.contains_id("xmx") {
            args.get_one::<String>("xmx").unwrap()
        } else {
            &binding_0
        },
    );

    let infiles = args.get_many::<String>("infiles").unwrap().collect_vec();

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

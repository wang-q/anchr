use clap::*;
use itertools::Itertools;
use std::collections::HashMap;
use std::path::PathBuf;
use std::{env, fs};
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name("trim")
        .about("Trim Illumina PE/SE fastq files")
        .after_help(
            r#"
<file1> [file2]

Fastq files can be gzipped
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
                .help("Quality threshold")
                .long("qual")
                .short("q")
                .takes_value(true)
                .default_value("25")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("filter")
                .long("filter")
                .help("Adapter, artifact, or both")
                .takes_value(true)
                .default_value("adapter")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("trimq")
                .long("trimq")
                .help("Quality score for 3' end")
                .takes_value(true)
                .default_value("15")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("trimk")
                .long("trimk")
                .help("Kmer for 5' adapter trimming")
                .takes_value(true)
                .default_value("23")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("matchk")
                .long("matchk")
                .help("Kmer for decontamination")
                .takes_value(true)
                .default_value("27")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("cutk")
                .long("cutk")
                .help("Kmer for cutoff")
                .takes_value(true)
                .default_value("31")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("adapter")
                .long("adapter")
                .help("The adapter file")
                .takes_value(true)
                .empty_values(false),
        )
        .arg(
            Arg::with_name("artifact")
                .long("artifact")
                .help("The artifact file")
                .takes_value(true)
                .empty_values(false),
        )
        .arg(
            Arg::with_name("prefix")
                .long("prefix")
                .help("Prefix of trimmed reads")
                .takes_value(true)
                .default_value("R")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("dedupe")
                .long("dedupe")
                .help("Do the dedupe step"),
        )
        .arg(
            Arg::with_name("tile")
                .long("tile")
                .help("With normal Illumina names, do tile-based filtering"),
        )
        .arg(
            Arg::with_name("cutoff")
                .long("cutoff")
                .help("Min kmer depth cutoff")
                .takes_value(true)
                .empty_values(false),
        )
        .arg(
            Arg::with_name("sample")
                .long("sample")
                .help("The sampling step")
                .takes_value(true)
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
                .help("Number of threads")
                .long("parallel")
                .short("p")
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
                .default_value("stdout")
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
    opt.insert("filter", args.value_of("filter").unwrap());
    opt.insert("trimq", args.value_of("trimq").unwrap());
    opt.insert("trimk", args.value_of("trimk").unwrap());
    opt.insert("matchk", args.value_of("matchk").unwrap());
    opt.insert("cutk", args.value_of("cutk").unwrap());
    opt.insert("prefix", args.value_of("prefix").unwrap());
    opt.insert("dedupe", if args.is_present("dedupe") { "1" } else { "0" });
    opt.insert("tile", if args.is_present("tile") { "1" } else { "0" });
    opt.insert(
        "cutoff",
        if args.is_present("cutoff") {
            args.value_of("cutoff").unwrap()
        } else {
            "0"
        },
    );
    opt.insert(
        "sample",
        if args.is_present("sample") {
            args.value_of("sample").unwrap()
        } else {
            "0"
        },
    );
    opt.insert(
        "xmx",
        if args.is_present("xmx") {
            args.value_of("xmx").unwrap()
        } else {
            "0"
        },
    );
    opt.insert("parallel", args.value_of("parallel").unwrap());

    // Default adapter and artifact files
    let path = if args.is_present("adapter") {
        PathBuf::from(args.value_of("adapter").unwrap())
            .canonicalize()
            .unwrap()
    } else {
        // write default adapter file
        let file = "illumina_adapters.fa";
        fs::write(file, include_str!("../../templates/illumina_adapters.fa"))?;
        env::current_dir()?.join(file).canonicalize().unwrap()
    };
    opt.insert("adapter", path.to_str().unwrap());

    let path = if args.is_present("artifact") {
        PathBuf::from(args.value_of("artifact").unwrap())
            .canonicalize()
            .unwrap()
    } else {
        // write default adapter file
        let file = "sequencing_artifacts.fa";
        fs::write(
            file,
            include_str!("../../templates/sequencing_artifacts.fa"),
        )?;
        env::current_dir()?.join(file).canonicalize().unwrap()
    };
    opt.insert("artifact", path.to_str().unwrap());

    let infiles = args.values_of("infiles").unwrap().collect_vec();

    let mut context = Context::new();
    context.insert("opt", &opt);
    context.insert("args", &infiles);

    // eprintln!("{:#?}", context);

    // many templates
    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("trim", include_str!("../../templates/trim.tera.sh")),
    ])
    .unwrap();

    // eprintln!("{:#?}", tera);

    let rendered = tera.render("trim", &context).unwrap();

    writer.write_all(rendered.as_ref())?;

    Ok(())
}

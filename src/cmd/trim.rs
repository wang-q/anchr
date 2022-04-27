use clap::*;
use itertools::Itertools;
use std::collections::HashMap;
use std::path::PathBuf;
use std::{env, fs};
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand<'a>() -> Command<'a> {
    Command::new("trim")
        .about("Trim Illumina PE/SE fastq files")
        .after_help(
            r#"
<file1> [file2]

Fastq files can be gzipped
"#,
        )
        .arg(
            Arg::new("infiles")
                .help("Sets the input file to use")
                .required(true)
                .min_values(1)
                .index(1),
        )
        .arg(
            Arg::new("qual")
                .long("qual")
                .short('q')
                .help("Quality threshold")
                .takes_value(true)
                .default_value("25")
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("len")
                .long("len")
                .short('l')
                .help("Filter reads less or equal to this length")
                .takes_value(true)
                .default_value("60")
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("filter")
                .long("filter")
                .help("Adapter, artifact, or both")
                .takes_value(true)
                .default_value("adapter")
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("trimq")
                .long("trimq")
                .help("Quality score for 3' end")
                .takes_value(true)
                .default_value("15")
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("trimk")
                .long("trimk")
                .help("Kmer for 5' adapter trimming")
                .takes_value(true)
                .default_value("23")
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("matchk")
                .long("matchk")
                .help("Kmer for decontamination")
                .takes_value(true)
                .default_value("27")
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("cutk")
                .long("cutk")
                .help("Kmer for cutoff")
                .takes_value(true)
                .default_value("31")
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("adapter")
                .long("adapter")
                .help("The adapter file")
                .takes_value(true)
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("artifact")
                .long("artifact")
                .help("The artifact file")
                .takes_value(true)
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("prefix")
                .long("prefix")
                .help("Prefix of trimmed reads")
                .takes_value(true)
                .default_value("R")
                .forbid_empty_values(true),
        )
        .arg(Arg::new("dedupe").long("dedupe").help("Do the dedupe step"))
        .arg(
            Arg::new("tile")
                .long("tile")
                .help("With normal Illumina names, do tile-based filtering"),
        )
        .arg(
            Arg::new("cutoff")
                .long("cutoff")
                .help("Min kmer depth cutoff")
                .takes_value(true)
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("sample")
                .long("sample")
                .help("The sampling step")
                .takes_value(true)
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("xmx")
                .long("xmx")
                .help("Set Java memory usage")
                .takes_value(true)
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .help("Number of threads")
                .takes_value(true)
                .default_value("8")
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .help("Output filename. [stdout] for screen")
                .takes_value(true)
                .default_value("trim.sh")
                .forbid_empty_values(true),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> std::result::Result<(), std::io::Error> {
    let mut writer = intspan::writer(args.value_of("outfile").unwrap());

    // context from args
    let mut opt = HashMap::new();
    opt.insert("qual", args.value_of("qual").unwrap());
    opt.insert("len", args.value_of("len").unwrap());
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

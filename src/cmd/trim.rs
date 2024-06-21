use clap::*;
use itertools::Itertools;
use std::collections::HashMap;
use std::path::PathBuf;
use std::{env, fs};
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
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
                .num_args(1..)
                .index(1),
        )
        .arg(
            Arg::new("qual")
                .long("qual")
                .short('q')
                .help("Quality threshold")
                .num_args(1)
                .default_value("25"),
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
            Arg::new("filter")
                .long("filter")
                .help("Adapter, artifact, or both")
                .num_args(1)
                .default_value("adapter"),
        )
        .arg(
            Arg::new("trimq")
                .long("trimq")
                .help("Quality score for 3' end")
                .num_args(1)
                .default_value("15"),
        )
        .arg(
            Arg::new("trimk")
                .long("trimk")
                .help("Kmer for 5' adapter trimming")
                .num_args(1)
                .default_value("23"),
        )
        .arg(
            Arg::new("matchk")
                .long("matchk")
                .help("Kmer for decontamination")
                .num_args(1)
                .default_value("27"),
        )
        .arg(
            Arg::new("cutk")
                .long("cutk")
                .help("Kmer for cutoff")
                .num_args(1)
                .default_value("31"),
        )
        .arg(
            Arg::new("adapter")
                .long("adapter")
                .help("The adapter file")
                .num_args(1),
        )
        .arg(
            Arg::new("artifact")
                .long("artifact")
                .help("The artifact file")
                .num_args(1),
        )
        .arg(
            Arg::new("prefix")
                .long("prefix")
                .help("Prefix of trimmed reads")
                .num_args(1)
                .default_value("R"),
        )
        .arg(
            Arg::new("dedupe")
                .long("dedupe")
                .action(ArgAction::SetTrue)
                .help("Do the dedupe step"),
        )
        .arg(
            Arg::new("tile")
                .long("tile")
                .action(ArgAction::SetTrue)
                .help("With normal Illumina names, do tile-based filtering"),
        )
        .arg(
            Arg::new("cutoff")
                .long("cutoff")
                .help("Min kmer depth cutoff")
                .num_args(1),
        )
        .arg(
            Arg::new("sample")
                .long("sample")
                .help("The sampling step")
                .num_args(1),
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
                .default_value("trim.sh"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    // context from args
    let mut opt = HashMap::new();
    opt.insert("qual", args.get_one::<String>("qual").unwrap());
    opt.insert("len", args.get_one::<String>("len").unwrap());
    opt.insert("filter", args.get_one::<String>("filter").unwrap());
    opt.insert("trimq", args.get_one::<String>("trimq").unwrap());
    opt.insert("trimk", args.get_one::<String>("trimk").unwrap());
    opt.insert("matchk", args.get_one::<String>("matchk").unwrap());
    opt.insert("cutk", args.get_one::<String>("cutk").unwrap());
    opt.insert("prefix", args.get_one::<String>("prefix").unwrap());

    let binding_1 = "1".to_string();
    let binding_0 = "0".to_string();

    opt.insert(
        "dedupe",
        if args.get_flag("dedupe") {
            &binding_1
        } else {
            &binding_0
        },
    );
    opt.insert(
        "tile",
        if args.get_flag("tile") {
            &binding_1
        } else {
            &binding_0
        },
    );

    opt.insert(
        "cutoff",
        if args.contains_id("cutoff") {
            args.get_one::<String>("cutoff").unwrap()
        } else {
            &binding_0
        },
    );

    opt.insert(
        "sample",
        if args.contains_id("sample") {
            args.get_one::<String>("sample").unwrap()
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

    opt.insert("parallel", args.get_one::<String>("parallel").unwrap());

    // Default adapter and artifact files
    let path = if args.contains_id("adapter") {
        PathBuf::from(args.get_one::<String>("adapter").unwrap())
            .canonicalize()
            .unwrap()
    } else {
        // write default adapter file
        let file = "illumina_adapters.fa";
        fs::write(file, include_str!("../../templates/illumina_adapters.fa"))?;
        env::current_dir()?.join(file).canonicalize().unwrap()
    };
    let binding = path.to_str().unwrap().to_string();
    opt.insert("adapter", &binding);

    let path = if args.contains_id("artifact") {
        PathBuf::from(args.get_one::<String>("artifact").unwrap())
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
    let binding = path.to_str().unwrap().to_string();
    opt.insert("artifact", &binding);

    let infiles = args.get_many::<String>("infiles").unwrap().collect_vec();

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

use clap::*;
use itertools::Itertools;
use std::collections::HashMap;
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("unitigs")
        .about("Create unitigs from trimmed/merged reads")
        .after_help(
            r#"
<pe.cor.fa> <env.json>

Fasta files can't be gzipped
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
            Arg::new("unitigger")
                .long("unitigger")
                .short('u')
                .help("Which unitig constructor to use: bcalm, bifrost, superreads, or tadpole")
                .num_args(1)
                .default_value("superreads"),
        )
        .arg(
            Arg::new("estsize")
                .long("estsize")
                .help("Estimated genome size")
                .num_args(1)
                .default_value("auto"),
        )
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .help("K-mer size to be used")
                .num_args(1)
                .default_value("31"),
        )
        .arg(
            Arg::new("min")
                .long("min")
                .help("Minimal length of unitigs")
                .num_args(1)
                .default_value("1000"),
        )
        .arg(
            Arg::new("merge")
                .long("merge")
                .action(ArgAction::SetTrue)
                .help("Merge unitigs from all k-mers"),
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
                .default_value("unitigs.sh"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    // context from args
    let mut opt = HashMap::new();
    opt.insert("unitigger", args.get_one::<String>("unitigger").unwrap());
    opt.insert("estsize", args.get_one::<String>("estsize").unwrap());
    opt.insert("kmer", args.get_one::<String>("kmer").unwrap());
    opt.insert("min", args.get_one::<String>("min").unwrap());

    let binding = if args.get_flag("merge") {
        "1".to_string()
    } else {
        "0".to_string()
    };
    opt.insert("merge", &binding);
    opt.insert("parallel", args.get_one::<String>("parallel").unwrap());

    let infiles = args.get_many::<String>("infiles").unwrap().collect_vec();

    let mut context = Context::new();
    context.insert("opt", &opt);
    context.insert("args", &infiles);

    // eprintln!("{:#?}", context);

    // many templates
    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/unitigs.tera.sh")),
    ])
    .unwrap();

    // eprintln!("{:#?}", tera);

    let rendered = tera.render("t", &context).unwrap();

    writer.write_all(rendered.as_ref())?;

    Ok(())
}

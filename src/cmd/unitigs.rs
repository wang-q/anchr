use clap::*;
use itertools::Itertools;
use std::collections::HashMap;
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand<'a>() -> Command<'a> {
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
                .min_values(1)
                .index(1),
        )
        .arg(
            Arg::new("unitigger")
                .long("unitigger")
                .short('u')
                .help("Which unitig constructor to use: bcalm, bifrost, superreads, or tadpole")
                .takes_value(true)
                .default_value("superreads")
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("estsize")
                .long("estsize")
                .help("Estimated genome size")
                .takes_value(true)
                .default_value("auto")
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .help("K-mer size to be used")
                .takes_value(true)
                .default_value("31")
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("min")
                .long("min")
                .help("Minimal length of unitigs")
                .takes_value(true)
                .default_value("1000")
                .forbid_empty_values(true),
        )
        .arg(
            Arg::new("merge")
                .long("merge")
                .help("Merge unitigs from all k-mers")
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
                .default_value("unitigs.sh")
                .forbid_empty_values(true),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> std::result::Result<(), std::io::Error> {
    let mut writer = intspan::writer(args.value_of("outfile").unwrap());

    // context from args
    let mut opt = HashMap::new();
    opt.insert("unitigger", args.value_of("unitigger").unwrap());
    opt.insert("estsize", args.value_of("estsize").unwrap());
    opt.insert("kmer", args.value_of("kmer").unwrap());
    opt.insert("min", args.value_of("min").unwrap());

    opt.insert("merge", if args.is_present("merge") { "1" } else { "0" });
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
        ("t", include_str!("../../templates/unitigs.tera.sh")),
    ])
    .unwrap();

    // eprintln!("{:#?}", tera);

    let rendered = tera.render("t", &context).unwrap();

    writer.write_all(rendered.as_ref())?;

    Ok(())
}

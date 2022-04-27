use clap::*;
use itertools::Itertools;
use std::collections::HashMap;
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand<'a>() -> Command<'a> {
    Command::new("quorum")
        .about("Run quorum to discard bad reads")
        .after_help(
            r#"
<PE file1> [PE file2] [SE file]

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
            Arg::new("jf")
                .long("coverage")
                .help("Jellyfish hash size")
                .takes_value(true)
                .default_value("500000000")
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
            Arg::new("prefix")
                .long("prefix")
                .help("Prefix of .cor.fa.gz")
                .takes_value(true)
                .default_value("pe")
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
                .default_value("quorum.sh")
                .forbid_empty_values(true),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> std::result::Result<(), std::io::Error> {
    let mut writer = intspan::writer(args.value_of("outfile").unwrap());

    // context from args
    let mut opt = HashMap::new();
    opt.insert("jf", args.value_of("jf").unwrap());
    opt.insert("estsize", args.value_of("estsize").unwrap());
    opt.insert("prefix", args.value_of("prefix").unwrap());
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
        ("quorum", include_str!("../../templates/quorum.tera.sh")),
    ])
    .unwrap();

    // eprintln!("{:#?}", tera);

    let rendered = tera.render("quorum", &context).unwrap();

    writer.write_all(rendered.as_ref())?;

    Ok(())
}

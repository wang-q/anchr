use clap::*;
use itertools::Itertools;
use std::collections::HashMap;
use std::io::Write;
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
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
                .num_args(1..)
                .index(1),
        )
        .arg(
            Arg::new("prefix")
                .long("prefix")
                .help("Prefix of .cor.fa.gz")
                .num_args(1)
                .default_value("pe"),
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
                .default_value("quorum.sh"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = pgr::libs::io::writer(args.get_one::<String>("outfile").unwrap())?;

    // context from args
    let mut opt = HashMap::new();
    opt.insert("prefix", args.get_one::<String>("prefix").unwrap());
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
        ("quorum", include_str!("../../templates/quorum.tera.sh")),
    ])
    .unwrap();

    // eprintln!("{:#?}", tera);

    let rendered = tera.render("quorum", &context).unwrap();

    writer.write_all(rendered.as_ref())?;

    Ok(())
}

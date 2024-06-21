use clap::*;
use itertools::Itertools;
use std::collections::HashMap;
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("anchors")
        .about("Select anchors (proper covered regions) from contigs")
        .after_help(
            r#"
<contig.fasta> <pe.cor.fa> [more reads]

Fasta files can‘t be gzipped

To get single-copy regions, set --uscale to 1.5

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
            Arg::new("min")
                .long("min")
                .help("Minimal length of anchors")
                .num_args(1)
                .default_value("1000"),
        )
        .arg(
            Arg::new("mincov")
                .long("mincov")
                .help("Minimal coverage of reads")
                .num_args(1)
                .default_value("5"),
        )
        .arg(
            Arg::new("readl")
                .long("readl")
                .help("Length of reads")
                .num_args(1)
                .default_value("100"),
        )
        .arg(
            Arg::new("mscale")
                .long("mscale")
                .help("The scale factor for MAD, median +/- k * MAD")
                .num_args(1)
                .default_value("3"),
        )
        .arg(
            Arg::new("lscale")
                .long("lscale")
                .help("The scale factor for lower, (median - k * MAD) / l")
                .num_args(1)
                .default_value("3"),
        )
        .arg(
            Arg::new("uscale")
                .long("uscale")
                .help("The scale factor for upper, (median + k * MAD) * u")
                .num_args(1)
                .default_value("2"),
        )
        .arg(
            Arg::new("fill")
                .long("fill")
                .help("Fill holes short than or equal to this")
                .num_args(1)
                .default_value("1"),
        )
        .arg(
            Arg::new("ratio")
                .long("ratio")
                .help("Fill large holes (opt.fill * 10) when covered ratio larger than this")
                .num_args(1)
                .default_value("0.98"),
        )
        .arg(
            Arg::new("longest")
                .long("longest")
                .action(ArgAction::SetTrue)
                .help("Only keep the longest proper region"),
        )
        .arg(
            Arg::new("keepedge")
                .long("keepedge")
                .action(ArgAction::SetTrue)
                .help("Keep edges of anchors"),
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
                .default_value("anchors.sh"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    // context from args
    let mut opt = HashMap::new();
    opt.insert("min", args.get_one::<String>("min").unwrap());
    opt.insert("mincov", args.get_one::<String>("mincov").unwrap());
    opt.insert("readl", args.get_one::<String>("readl").unwrap());
    opt.insert("mscale", args.get_one::<String>("mscale").unwrap());
    opt.insert("lscale", args.get_one::<String>("lscale").unwrap());
    opt.insert("uscale", args.get_one::<String>("uscale").unwrap());
    opt.insert("ratio", args.get_one::<String>("ratio").unwrap());
    opt.insert("fill", args.get_one::<String>("fill").unwrap());

    let binding_1 = "1".to_string();
    let binding_0 = "0".to_string();

    opt.insert(
        "longest",
        if args.get_flag("longest") {
            &binding_1
        } else {
            &binding_0
        },
    );
    opt.insert(
        "keepedge",
        if args.get_flag("keepedge") {
            &binding_1
        } else {
            &binding_0
        },
    );

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
        ("t", include_str!("../../templates/anchors.tera.sh")),
    ])
    .unwrap();

    // eprintln!("{:#?}", tera);

    let rendered = tera.render("t", &context).unwrap();

    writer.write_all(rendered.as_ref())?;

    Ok(())
}

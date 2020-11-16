use clap::*;
use itertools::Itertools;
use std::collections::HashMap;
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name("anchors")
        .about("Select anchors (proper covered regions) from contigs")
        .after_help(
            r#"
<contig.fasta> <pe.cor.fa> [more reads]

Fasta files can‘t be gzipped

To get single-copy regions, set --uscale to 1.5

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
            Arg::with_name("min")
                .long("min")
                .help("Minimal length of anchors")
                .takes_value(true)
                .default_value("1000")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("mincov")
                .long("mincov")
                .help("Minimal coverage of reads")
                .takes_value(true)
                .default_value("3")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("readl")
                .long("readl")
                .help("Length of reads")
                .takes_value(true)
                .default_value("100")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("mscale")
                .long("mscale")
                .help("The scale factor for MAD, median +/- k * MAD")
                .takes_value(true)
                .default_value("3")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("lscale")
                .long("lscale")
                .help("The scale factor for lower, (median - k * MAD) / l")
                .takes_value(true)
                .default_value("3")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("uscale")
                .long("uscale")
                .help("The scale factor for upper, (median + k * MAD) * u")
                .takes_value(true)
                .default_value("3")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("fill")
                .long("fill")
                .help("Fill holes short than or equal to this")
                .takes_value(true)
                .default_value("1")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("ratio")
                .long("ratio")
                .help("Fill large holes (opt.fill * 10) when covered ratio larger than this")
                .takes_value(true)
                .default_value("0.98")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("longest")
                .long("longest")
                .help("Only keep the longest proper region"),
        )
        .arg(
            Arg::with_name("keepedge")
                .long("keepedge")
                .help("Keep edges of anchors"),
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
                .default_value("anchors.sh")
                .empty_values(false),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> std::result::Result<(), std::io::Error> {
    let mut writer = intspan::writer(args.value_of("outfile").unwrap());

    // context from args
    let mut opt = HashMap::new();
    opt.insert("min", args.value_of("min").unwrap());
    opt.insert("mincov", args.value_of("mincov").unwrap());
    opt.insert("readl", args.value_of("readl").unwrap());
    opt.insert("mscale", args.value_of("mscale").unwrap());
    opt.insert("lscale", args.value_of("lscale").unwrap());
    opt.insert("uscale", args.value_of("uscale").unwrap());
    opt.insert("ratio", args.value_of("ratio").unwrap());
    opt.insert("fill", args.value_of("fill").unwrap());

    opt.insert(
        "longest",
        if args.is_present("longest") { "1" } else { "0" },
    );
    opt.insert(
        "keepedge",
        if args.is_present("keepedge") {
            "1"
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
        ("t", include_str!("../../templates/anchors.tera.sh")),
    ])
    .unwrap();

    // eprintln!("{:#?}", tera);

    let rendered = tera.render("t", &context).unwrap();

    writer.write_all(rendered.as_ref())?;

    Ok(())
}

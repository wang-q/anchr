use clap::*;
use std::collections::HashMap;
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand<'a, 'b>() -> App<'a, 'b> {
    SubCommand::with_name("template")
        .about("Creates Bash scripts")
        .after_help(
            r#"
* Info

    * --genome
    * --se
    * --xmx
    * --parallel 8
    * --queue mpi

* Quality check

    * --fastqc
    * --kmergenie

* Trimming

    * --trim "--dedupe"
    * --qual "25 30"
    * --len "60"
    * --filter "adapter"

* Post-trimming

    * --quorum
    * --merge
    * --prefilter
    * --ecphase "1 2 3"

* Downsampling

    * --cov "40 80"
    * --splitp 10
    * --statp 2

"#,
        )
        // Info
        .arg(
            Arg::with_name("genome")
                .long("genome")
                .help("Your best guess of the haploid genome size")
                .takes_value(true)
                .empty_values(false),
        )
        .arg(Arg::with_name("se").long("se").help("Single end mode"))
        .arg(
            Arg::with_name("xmx")
                .long("xmx")
                .help("Set Java memory usage")
                .takes_value(true)
                .empty_values(false),
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
            Arg::with_name("queue")
                .long("queue")
                .help("Queue name of the LSF cluster")
                .takes_value(true)
                .default_value("mpi")
                .empty_values(false),
        )
        // Quality check
        .arg(Arg::with_name("fastqc").long("fastqc").help("Run FastQC"))
        .arg(
            Arg::with_name("kmergenie")
                .long("kmergenie")
                .help("Run KmerGenie"),
        )
        // Trimming
        .arg(
            Arg::with_name("trim")
                .long("trim")
                .help("Opts for trim")
                .takes_value(true)
                .default_value("--dedupe")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("sample")
                .long("sample")
                .help("Sampling coverage")
                .takes_value(true)
                .empty_values(false),
        )
        .arg(
            Arg::with_name("qual")
                .long("qual")
                .help("Quality threshold")
                .takes_value(true)
                .default_value("25 30")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("len")
                .long("len")
                .help("Filter reads less or equal to this length")
                .takes_value(true)
                .default_value("60")
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
        // Post-trimming
        .arg(Arg::with_name("quorum").long("quorum").help("Run quorum"))
        .arg(
            Arg::with_name("merge")
                .long("merge")
                .help("Run merge reads"),
        )
        .arg(
            Arg::with_name("prefilter")
                .long("prefilter")
                .help("Prefilter=N (1 or 2) for tadpole and bbmerge")
                .takes_value(true)
                .empty_values(false),
        )
        .arg(
            Arg::with_name("ecphase")
                .long("ecphase")
                .help("Error-correct phases. Phase 2 can be skipped")
                .takes_value(true)
                .default_value("1 2 3")
                .empty_values(false),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> std::result::Result<(), std::io::Error> {
    //----------------------------
    // context from args
    //----------------------------
    let mut opt = HashMap::new();

    opt.insert(
        "genome",
        if args.is_present("genome") {
            args.value_of("genome").unwrap()
        } else {
            "0"
        },
    );
    opt.insert("se", if args.is_present("se") { "1" } else { "0" });
    opt.insert(
        "xmx",
        if args.is_present("xmx") {
            args.value_of("xmx").unwrap()
        } else {
            "0"
        },
    );
    opt.insert("parallel", args.value_of("parallel").unwrap());
    opt.insert("queue", args.value_of("queue").unwrap());

    opt.insert("trim", args.value_of("trim").unwrap());
    opt.insert(
        "sample",
        if args.is_present("sample") {
            args.value_of("sample").unwrap()
        } else {
            "0"
        },
    );
    opt.insert("qual", args.value_of("qual").unwrap());
    opt.insert("len", args.value_of("len").unwrap());
    opt.insert("filter", args.value_of("filter").unwrap());

    opt.insert(
        "prefilter",
        if args.is_present("prefilter") {
            args.value_of("prefilter").unwrap()
        } else {
            "0"
        },
    );
    opt.insert("ecphase", args.value_of("ecphase").unwrap());

    let mut context = Context::new();
    context.insert("opt", &opt);

    //----------------------------
    // create scripts
    //----------------------------
    if args.is_present("fastqc") {
        gen_fastqc(&context)?;
    }

    if args.is_present("kmergenie") {
        gen_kmergenie(&context)?;
    }

    gen_trim(&context)?;

    if !args.is_present("se") && args.is_present("merge") {
        gen_merge(&context)?;
    }

    Ok(())
}

fn gen_fastqc(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "2_fastqc.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_fastqc.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_kmergenie(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "2_kmergenie.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_kmergenie.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_trim(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "2_trim.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_trim.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_merge(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "2_merge.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_merge.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

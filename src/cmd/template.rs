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
    * --queue mpi
    * --xmx
    * --parallel 8

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
}

// command implementation
pub fn execute(args: &ArgMatches) -> std::result::Result<(), std::io::Error> {
    // context from args
    let mut opt = HashMap::new();

    opt.insert("se", if args.is_present("se") { "1" } else { "0" });
    opt.insert("parallel", args.value_of("parallel").unwrap());

    opt.insert("fastqc", if args.is_present("fastqc") { "1" } else { "0" });
    opt.insert(
        "kmergenie",
        if args.is_present("kmergenie") {
            "1"
        } else {
            "0"
        },
    );

    opt.insert("trim", args.value_of("trim").unwrap());
    opt.insert("qual", args.value_of("qual").unwrap());
    opt.insert("len", args.value_of("len").unwrap());
    opt.insert("filter", args.value_of("filter").unwrap());

    let mut context = Context::new();
    context.insert("opt", &opt);

    // create scripts
    if args.is_present("fastqc") {
        gen_fastqc(&context)?;
    }

    if args.is_present("kmergenie") {
        gen_kmergenie(&context)?;
    }

    Ok(())
}

fn gen_fastqc(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "2_fastqc.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("template", include_str!("../../templates/2_fastqc.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("template", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_kmergenie(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "2_kmergenie.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("template", include_str!("../../templates/2_kmergenie.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("template", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

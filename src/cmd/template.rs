use clap::*;
use itertools::Itertools;
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
    * --kat
    * --fastk
    * --insertsize
    * --reads 1000000

* Trimming

    * --trim "--dedupe"
    * --sample "300"
    * --qual "25 30"
    * --len "60"
    * --filter "adapter"

* Post-trimming

    * --quorum
    * --merge
    * --prefilter
    * --ecphase "1 2 3"

* Mapping
    * --bwa
    * --gatk

* Down sampling, unitigs, and anchors

    * --cov "40 80"
    * --unitigger "bcalm"
    * --splitp 20
    * --statp 2
    * --readl 100
    * --uscale 2
    * --lscale 3
    * --redo

* Extend anchors

    * --extend
    * --gluemin 30
    * --fillmax 100

"#,
        )
        // Info
        .arg(
            Arg::with_name("genome")
                .long("genome")
                .help("Your best guess of the haploid genome size")
                .takes_value(true)
                .default_value("1000000")
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
                .empty_values(false),
        )
        // Quality check
        .arg(Arg::with_name("fastqc").long("fastqc").help("Run FastQC"))
        .arg(Arg::with_name("kat").long("kat").help("Run KAT"))
        .arg(Arg::with_name("fastk").long("fastk").help("Run FastK"))
        .arg(
            Arg::with_name("insertsize")
                .long("insertsize")
                .help("Calc insert sizes"),
        )
        .arg(
            Arg::with_name("reads")
                .long("reads")
                .help("How many reads to estimate insert sizes")
                .takes_value(true)
                .default_value("1000000")
                .empty_values(false),
        )
        // Trimming
        .arg(
            Arg::with_name("trim")
                .long("trim")
                .help("Opts for trim")
                .takes_value(true)
                .default_value("--dedupe")
                .empty_values(false)
                .allow_hyphen_values(true),
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
                .help("Prefilter=N (1 or 2) for tadpole and bbmerge, 1 use less memories")
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
        // Mapping
        .arg(
            Arg::with_name("bwa")
                .long("bwa")
                .help("Map trimmed reads to the genome")
                .takes_value(true)
                .empty_values(false),
        )
        .arg(
            Arg::with_name("gatk")
                .long("gatk")
                .help("Calling variants with GATK Mutect2"),
        )
        // Down sampling, unitigs, and anchors
        .arg(
            Arg::with_name("cov")
                .long("cov")
                .help("Down sampling coverages")
                .takes_value(true)
                .default_value("40 80")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("unitigger")
                .long("unitigger")
                .short("u")
                .help("Unitigger used: bcalm, bifrost, superreads, or tadpole")
                .takes_value(true)
                .default_value("bcalm")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("splitp")
                .long("splitp")
                .help("Parts of splitting")
                .takes_value(true)
                .default_value("20")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("statp")
                .long("statp")
                .help("Parts of stats")
                .takes_value(true)
                .default_value("2")
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
            Arg::with_name("uscale")
                .long("uscale")
                .help("The scale factor for upper, (median + k * MAD) * u")
                .takes_value(true)
                .default_value("2")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("lscale")
                .long("lscale")
                .help("The scale factor for upper, (median - k * MAD) / l")
                .takes_value(true)
                .default_value("3")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("redo")
                .long("redo")
                .help("Redo anchors when merging anchors"),
        )
        // Extend anchors
        .arg(
            Arg::with_name("extend")
                .long("extend")
                .help("Extend anchors with other contigs"),
        )
        .arg(
            Arg::with_name("gluemin")
                .long("gluemin")
                .help("Min length of overlaps to be glued")
                .takes_value(true)
                .default_value("30")
                .empty_values(false),
        )
        .arg(
            Arg::with_name("fillmax")
                .long("fillmax")
                .help("Max length of gaps")
                .takes_value(true)
                .default_value("100")
                .empty_values(false),
        )
        // Extend anchors
        .arg(Arg::with_name("busco").long("busco").help("Run busco"))
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
    opt.insert(
        "queue",
        if args.is_present("queue") {
            args.value_of("queue").unwrap()
        } else {
            "0"
        },
    );

    opt.insert("reads", args.value_of("reads").unwrap());

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

    opt.insert("merge", if args.is_present("merge") { "1" } else { "0" });
    opt.insert(
        "prefilter",
        if args.is_present("prefilter") {
            args.value_of("prefilter").unwrap()
        } else {
            "0"
        },
    );
    opt.insert("ecphase", args.value_of("ecphase").unwrap());

    opt.insert(
        "bwa",
        if args.is_present("bwa") {
            args.value_of("bwa").unwrap()
        } else {
            "0"
        },
    );
    opt.insert("gatk", if args.is_present("gatk") { "1" } else { "0" });

    opt.insert("cov", args.value_of("cov").unwrap());
    opt.insert("unitigger", args.value_of("unitigger").unwrap());
    opt.insert("splitp", args.value_of("splitp").unwrap());
    opt.insert("statp", args.value_of("statp").unwrap());
    opt.insert("readl", args.value_of("readl").unwrap());
    opt.insert("uscale", args.value_of("uscale").unwrap());
    opt.insert("lscale", args.value_of("lscale").unwrap());
    opt.insert("redo", if args.is_present("redo") { "1" } else { "0" });

    opt.insert("extend", if args.is_present("extend") { "1" } else { "0" });
    opt.insert("gluemin", args.value_of("gluemin").unwrap());
    opt.insert("fillmax", args.value_of("fillmax").unwrap());

    let mut context = Context::new();
    context.insert("opt", &opt);

    //----------------------------
    // create scripts
    //----------------------------
    if args.is_present("fastqc") {
        gen_fastqc(&context)?;
    }
    if args.is_present("insertsize") {
        gen_insert_size(&context)?;
    }
    if args.is_present("kat") {
        gen_kat(&context)?;
    }
    if args.is_present("fastk") {
        gen_fastk(&context)?;
    }

    gen_trim(&context)?;

    gen_stat_reads(&context)?;

    if args.is_present("bwa") {
        gen_bwa(&context)?;
    }
    if args.is_present("gatk") {
        gen_gatk(&context)?;
    }

    if args.is_present("quorum") {
        gen_quorum(&context)?;
    } else {
        gen_no_quorum(&context)?;
    }
    gen_down_sampling(&context)?;

    let unitiggers = args
        .value_of("unitigger")
        .unwrap()
        .split_ascii_whitespace()
        .collect_vec();

    for u in unitiggers.clone() {
        gen_unitigs(&context, u)?;
    }
    gen_anchors(&context)?;
    gen_stat_anchors(&context)?;

    if !args.is_present("se") && args.is_present("merge") {
        gen_merge(&context)?;
        gen_mr_down_sampling(&context)?;
        for u in unitiggers.clone() {
            gen_mr_unitigs(&context, u)?;
        }
        gen_mr_anchors(&context)?;
        gen_stat_mr_anchors(&context)?;
    }

    gen_merge_anchors(&context)?;
    gen_stat_merge_anchors(&context)?;

    gen_spades(&context)?;
    gen_megahit(&context)?;
    gen_platanus(&context)?;
    if !args.is_present("se") && args.is_present("merge") {
        gen_mr_spades(&context)?;
        gen_mr_megahit(&context)?;
    }
    gen_stat_other_anchors(&context)?;

    if args.is_present("extend") {
        gen_glue_anchors(&context)?;
        gen_fill_anchors(&context)?;
    }

    gen_quast(&context)?;
    gen_stat_final(&context)?;

    if args.is_present("busco") {
        gen_busco(&context)?;
    }

    gen_cleanup(&context)?;
    gen_real_clean(&context)?;
    gen_master(&context)?;
    if args.is_present("queue") {
        gen_bsub(&context)?;
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

fn gen_insert_size(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "2_insert_size.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_insert_size.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_kat(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "2_kat.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_kat.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_fastk(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "2_fastk.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_fastk.tera.sh")),
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

fn gen_stat_reads(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "9_stat_reads.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/9_stat_reads.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_quorum(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "2_quorum.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_quorum.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_no_quorum(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "2_quorum.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_no_quorum.tera.sh")),
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

fn gen_bwa(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "3_bwa.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/3_bwa.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_gatk(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "3_gatk.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/3_gatk.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_down_sampling(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "4_down_sampling.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/4_down_sampling.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_mr_down_sampling(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "6_down_sampling.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/6_down_sampling.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_unitigs(context: &Context, unitigger: &str) -> std::result::Result<(), std::io::Error> {
    let outname = format!("4_unitigs_{}.sh", unitigger);
    eprintln!("Create {}", outname);

    let mut con = Context::new();
    con.insert("outname", outname.as_str());
    con.insert("unitigger", unitigger);
    con.extend(context.clone());

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/4_unitigs.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &con).unwrap();
    intspan::write_lines(outname.as_str(), &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_mr_unitigs(context: &Context, unitigger: &str) -> std::result::Result<(), std::io::Error> {
    let outname = format!("6_unitigs_{}.sh", unitigger);
    eprintln!("Create {}", outname);

    let mut con = Context::new();
    con.insert("outname", outname.as_str());
    con.insert("unitigger", unitigger);
    con.extend(context.clone());

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/6_unitigs.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &con).unwrap();
    intspan::write_lines(outname.as_str(), &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_anchors(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "4_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/4_anchors.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_mr_anchors(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "6_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/6_anchors.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_stat_anchors(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "9_stat_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/9_stat_anchors.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_stat_mr_anchors(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "9_stat_mr_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        (
            "t",
            include_str!("../../templates/9_stat_mr_anchors.tera.sh"),
        ),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_merge_anchors(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "7_merge_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/7_merge_anchors.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_stat_merge_anchors(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "9_stat_merge_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        (
            "t",
            include_str!("../../templates/9_stat_merge_anchors.tera.sh"),
        ),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_spades(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "8_spades.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/8_spades.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_mr_spades(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "8_mr_spades.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/8_mr_spades.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_megahit(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "8_megahit.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/8_megahit.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_mr_megahit(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "8_mr_megahit.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/8_mr_megahit.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_platanus(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "8_platanus.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/8_platanus.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_stat_other_anchors(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "9_stat_other_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        (
            "t",
            include_str!("../../templates/9_stat_other_anchors.tera.sh"),
        ),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_glue_anchors(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "7_glue_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/7_glue_anchors.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_fill_anchors(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "7_fill_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/7_fill_anchors.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_quast(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "9_quast.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/9_quast.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_busco(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "9_busco.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/9_busco.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_stat_final(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "9_stat_final.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/9_stat_final.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_cleanup(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "0_cleanup.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/0_cleanup.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_real_clean(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "0_real_clean.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/0_real_clean.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_master(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "0_master.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/0_master.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_bsub(context: &Context) -> std::result::Result<(), std::io::Error> {
    let outname = "0_bsub.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/0_bsub.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

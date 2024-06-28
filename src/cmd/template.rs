use clap::*;
use itertools::Itertools;
use std::collections::HashMap;
use std::fs;
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("template")
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
            Arg::new("genome")
                .long("genome")
                .help("Your best guess of the haploid genome size")
                .num_args(1)
                .default_value("1000000"),
        )
        .arg(
            Arg::new("se")
                .long("se")
                .action(ArgAction::SetTrue)
                .help("Single end mode"),
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
            Arg::new("queue")
                .long("queue")
                .help("Queue name of the LSF cluster")
                .num_args(1),
        )
        // Quality check
        .arg(
            Arg::new("fastqc")
                .long("fastqc")
                .action(ArgAction::SetTrue)
                .help("Run FastQC"),
        )
        .arg(
            Arg::new("kat")
                .long("kat")
                .action(ArgAction::SetTrue)
                .help("Run KAT"),
        )
        .arg(
            Arg::new("fastk")
                .long("fastk")
                .action(ArgAction::SetTrue)
                .help("Run FastK"),
        )
        .arg(
            Arg::new("insertsize")
                .long("insertsize")
                .action(ArgAction::SetTrue)
                .help("Calc insert sizes"),
        )
        .arg(
            Arg::new("reads")
                .long("reads")
                .help("How many reads to estimate insert sizes")
                .num_args(1)
                .default_value("1000000"),
        )
        // Trimming
        .arg(
            Arg::new("trim")
                .long("trim")
                .help("Opts for trim")
                .num_args(1)
                .default_value("--dedupe")
                .allow_hyphen_values(true),
        )
        .arg(
            Arg::new("sample")
                .long("sample")
                .help("Sampling coverage")
                .num_args(1),
        )
        .arg(
            Arg::new("qual")
                .long("qual")
                .help("Quality threshold")
                .num_args(1)
                .default_value("25 30"),
        )
        .arg(
            Arg::new("len")
                .long("len")
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
        // Post-trimming
        .arg(
            Arg::new("quorum")
                .long("quorum")
                .action(ArgAction::SetTrue)
                .help("Run quorum"),
        )
        .arg(
            Arg::new("merge")
                .long("merge")
                .action(ArgAction::SetTrue)
                .help("Run merge reads"),
        )
        .arg(
            Arg::new("prefilter")
                .long("prefilter")
                .help("Prefilter=N (1 or 2) for tadpole and bbmerge, 1 use less memories")
                .num_args(1),
        )
        .arg(
            Arg::new("ecphase")
                .long("ecphase")
                .help("Error-correct phases. Phase 2 can be skipped")
                .num_args(1)
                .default_value("1 2 3"),
        )
        // Mapping
        .arg(
            Arg::new("bwa")
                .long("bwa")
                .help("Map trimmed reads to the genome")
                .num_args(1),
        )
        .arg(
            Arg::new("gatk")
                .long("gatk")
                .action(ArgAction::SetTrue)
                .help("Calling variants with GATK Mutect2"),
        )
        // Down sampling, unitigs, and anchors
        .arg(
            Arg::new("cov")
                .long("cov")
                .help("Down sampling coverages")
                .num_args(1)
                .default_value("40 80"),
        )
        .arg(
            Arg::new("unitigger")
                .long("unitigger")
                .short('u')
                .help("Unitigger used: bcalm, bifrost, superreads, or tadpole")
                .num_args(1)
                .default_value("bcalm"),
        )
        .arg(
            Arg::new("splitp")
                .long("splitp")
                .help("Parts of splitting")
                .num_args(1)
                .default_value("20"),
        )
        .arg(
            Arg::new("statp")
                .long("statp")
                .help("Parts of stats")
                .num_args(1)
                .default_value("2"),
        )
        .arg(
            Arg::new("readl")
                .long("readl")
                .help("Length of reads")
                .num_args(1)
                .default_value("100"),
        )
        .arg(
            Arg::new("uscale")
                .long("uscale")
                .help("The scale factor for upper, (median + k * MAD) * u")
                .num_args(1)
                .default_value("2"),
        )
        .arg(
            Arg::new("lscale")
                .long("lscale")
                .help("The scale factor for upper, (median - k * MAD) / l")
                .num_args(1)
                .default_value("3"),
        )
        .arg(
            Arg::new("redo")
                .long("redo")
                .action(ArgAction::SetTrue)
                .help("Redo anchors when merging anchors"),
        )
        // Extend anchors
        .arg(
            Arg::new("extend")
                .long("extend")
                .action(ArgAction::SetTrue)
                .help("Extend anchors with other contigs"),
        )
        .arg(
            Arg::new("gluemin")
                .long("gluemin")
                .help("Min length of overlaps to be glued")
                .num_args(1)
                .default_value("30"),
        )
        .arg(
            Arg::new("fillmax")
                .long("fillmax")
                .help("Max length of gaps")
                .num_args(1)
                .default_value("100"),
        )
        // Extend anchors
        .arg(
            Arg::new("busco")
                .long("busco")
                .action(ArgAction::SetTrue)
                .help("Run busco"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // context from args
    //----------------------------
    let mut opt = HashMap::new();

    let binding_1 = "1".to_string();
    let binding_0 = "0".to_string();

    opt.insert(
        "genome",
        if args.contains_id("genome") {
            args.get_one::<String>("genome").unwrap()
        } else {
            "0"
        },
    );
    opt.insert(
        "se",
        if args.get_flag("se") {
            &binding_1
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
    opt.insert(
        "queue",
        if args.contains_id("queue") {
            args.get_one::<String>("queue").unwrap()
        } else {
            &binding_0
        },
    );

    opt.insert("reads", args.get_one::<String>("reads").unwrap());

    opt.insert("trim", args.get_one::<String>("trim").unwrap());
    opt.insert(
        "sample",
        if args.contains_id("sample") {
            args.get_one::<String>("sample").unwrap()
        } else {
            &binding_0
        },
    );
    opt.insert("qual", args.get_one::<String>("qual").unwrap());
    opt.insert("len", args.get_one::<String>("len").unwrap());
    opt.insert("filter", args.get_one::<String>("filter").unwrap());

    opt.insert(
        "merge",
        if args.get_flag("merge") {
            &binding_1
        } else {
            &binding_0
        },
    );
    opt.insert(
        "prefilter",
        if args.contains_id("prefilter") {
            args.get_one::<String>("prefilter").unwrap()
        } else {
            "0"
        },
    );
    opt.insert("ecphase", args.get_one::<String>("ecphase").unwrap());

    opt.insert(
        "bwa",
        if args.contains_id("bwa") {
            args.get_one::<String>("bwa").unwrap()
        } else {
            &binding_0
        },
    );
    opt.insert(
        "gatk",
        if args.get_flag("gatk") {
            &binding_1
        } else {
            &binding_0
        },
    );

    opt.insert("cov", args.get_one::<String>("cov").unwrap());
    opt.insert("unitigger", args.get_one::<String>("unitigger").unwrap());
    opt.insert("splitp", args.get_one::<String>("splitp").unwrap());
    opt.insert("statp", args.get_one::<String>("statp").unwrap());
    opt.insert("readl", args.get_one::<String>("readl").unwrap());
    opt.insert("uscale", args.get_one::<String>("uscale").unwrap());
    opt.insert("lscale", args.get_one::<String>("lscale").unwrap());
    opt.insert(
        "redo",
        if args.get_flag("redo") {
            &binding_1
        } else {
            &binding_0
        },
    );
    opt.insert(
        "extend",
        if args.get_flag("extend") {
            &binding_1
        } else {
            &binding_0
        },
    );
    opt.insert("gluemin", args.get_one::<String>("gluemin").unwrap());
    opt.insert("fillmax", args.get_one::<String>("fillmax").unwrap());

    let mut context = Context::new();
    context.insert("opt", &opt);

    //----------------------------
    // create scripts
    //----------------------------
    fs::create_dir_all("9_markdown")?;

    if args.get_flag("fastqc") {
        gen_fastqc(&context)?;
    }
    if args.get_flag("insertsize") {
        gen_insert_size(&context)?;
    }
    if args.get_flag("kat") {
        gen_kat(&context)?;
    }
    if args.get_flag("fastk") {
        gen_fastk(&context)?;
        gen_genescopefk(&context)?;
    }

    gen_trim(&context)?;

    gen_stat_reads(&context)?;

    if args.contains_id("bwa") {
        gen_bwa(&context)?;
    }
    if args.get_flag("gatk") {
        gen_gatk(&context)?;
    }

    if args.get_flag("quorum") {
        gen_quorum(&context)?;
    } else {
        gen_no_quorum(&context)?;
    }
    gen_down_sampling(&context)?;

    let unitiggers = args
        .get_one::<String>("unitigger")
        .unwrap()
        .split_ascii_whitespace()
        .collect_vec();

    for u in unitiggers.clone() {
        gen_unitigs(&context, u)?;
    }
    gen_anchors(&context)?;
    gen_stat_anchors(&context)?;

    if !args.get_flag("se") && args.get_flag("merge") {
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
    if !args.get_flag("se") && args.get_flag("merge") {
        gen_mr_spades(&context)?;
        gen_mr_megahit(&context)?;
    }
    gen_stat_other_anchors(&context)?;

    if args.get_flag("extend") {
        gen_glue_anchors(&context)?;
        gen_fill_anchors(&context)?;
    }

    gen_quast(&context)?;
    gen_stat_final(&context)?;

    if args.get_flag("busco") {
        gen_busco(&context)?;
    }

    gen_cleanup(&context)?;
    gen_real_clean(&context)?;
    gen_master(&context)?;
    if args.contains_id("queue") {
        gen_bsub(&context)?;
    }

    Ok(())
}

fn gen_fastqc(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/2_fastqc.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_fastqc.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_insert_size(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/2_insert_size.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_insert_size.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_kat(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/2_kat.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_kat.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_fastk(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/2_fastk.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_fastk.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_genescopefk(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/genescopefk.R";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![("t", include_str!("../../templates/genescopefk.R"))])
        .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_trim(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/2_trim.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_trim.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_stat_reads(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/9_stat_reads.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/9_stat_reads.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_quorum(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/2_quorum.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_quorum.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_no_quorum(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/2_quorum.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_no_quorum.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_merge(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/2_merge.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/2_merge.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_bwa(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/3_bwa.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/3_bwa.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_gatk(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/3_gatk.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/3_gatk.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_down_sampling(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/4_down_sampling.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/4_down_sampling.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_mr_down_sampling(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/6_down_sampling.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/6_down_sampling.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_unitigs(context: &Context, unitigger: &str) -> anyhow::Result<()> {
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

fn gen_mr_unitigs(context: &Context, unitigger: &str) -> anyhow::Result<()> {
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

fn gen_anchors(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/4_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/4_anchors.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_mr_anchors(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/6_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/6_anchors.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_stat_anchors(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/9_stat_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/9_stat_anchors.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_stat_mr_anchors(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/9_stat_mr_anchors.sh";
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

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_merge_anchors(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/7_merge_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/7_merge_anchors.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_stat_merge_anchors(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/9_stat_merge_anchors.sh";
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

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_spades(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/8_spades.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/8_spades.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_mr_spades(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/8_mr_spades.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/8_mr_spades.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_megahit(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/8_megahit.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/8_megahit.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_mr_megahit(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/8_mr_megahit.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/8_mr_megahit.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_stat_other_anchors(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/9_stat_other_anchors.sh";
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

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_glue_anchors(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/7_glue_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/7_glue_anchors.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_fill_anchors(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/7_fill_anchors.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/7_fill_anchors.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_quast(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/9_quast.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/9_quast.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_busco(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/9_busco.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/9_busco.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_stat_final(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/9_stat_final.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/9_stat_final.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_cleanup(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/0_cleanup.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/0_cleanup.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_real_clean(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/0_real_clean.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/0_real_clean.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_master(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/0_master.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/0_master.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

fn gen_bsub(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/0_script/0_bsub.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/0_bsub.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    intspan::write_lines(outname, &vec![rendered.as_str()])?;

    Ok(())
}

use clap::*;
use std::collections::HashMap;
use std::fs;
use tera::{Context, Tera};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("template")
        .about("Creates scripts")
        .after_help(
            r#"
* Info
    * --genome
    * --se
    * --parallel 8
    * --queue mpi

* Resources
    * --repetitive

* Quality check
    * --fastqc
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
    * --merge

* Mapping
    * --bwa
    * --gatk

* Down sampling, unitigs, and anchors
    * --cov "40 80"
    * --splitp 20
    * --statp 2
    * --uscale 2
    * --lscale 3


* Validate assemblies
    * --busco

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
        // Resources
        .arg(
            Arg::new("repetitive")
                .long("repetitive")
                .action(ArgAction::SetTrue)
                .help("Find repetitive regions"),
        )
        // Quality check
        .arg(
            Arg::new("fastqc")
                .long("fastqc")
                .action(ArgAction::SetTrue)
                .help("Run FastQC"),
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
        .arg(
            Arg::new("merge")
                .long("merge")
                .action(ArgAction::SetTrue)
                .help("Run merge reads"),
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
        // Validate assemblies
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
    opt.insert("splitp", args.get_one::<String>("splitp").unwrap());
    opt.insert("statp", args.get_one::<String>("statp").unwrap());
    opt.insert("uscale", args.get_one::<String>("uscale").unwrap());
    opt.insert("lscale", args.get_one::<String>("lscale").unwrap());
    let mut context = Context::new();
    context.insert("opt", &opt);

    //----------------------------
    // create scripts
    //----------------------------
    fs::create_dir_all("0_script")?;
    fs::create_dir_all("9_markdown")?;

    if args.get_flag("repetitive") {
        gen_repetitive(&context)?;
    }

    if args.get_flag("fastqc") {
        gen_fastqc(&context)?;
    }
    if args.get_flag("insertsize") {
        gen_insert_size(&context)?;
    }
    if args.get_flag("fastk") {
        gen_fastk(&context)?;
    }

    gen_trim(&context)?;

    gen_stat_reads(&context)?;

    if args.contains_id("bwa") {
        gen_bwa(&context)?;
    }
    if args.get_flag("gatk") {
        gen_gatk(&context)?;
    }

    // s-filter replaced the external quorum (no_quorum fallback removed);
    // the filter step always runs.
    gen_quorum(&context)?;
    gen_down_sampling(&context)?;

    gen_unitigs(&context)?;
    gen_anchors(&context)?;
    gen_stat_anchors(&context)?;

    if !args.get_flag("se") && args.get_flag("merge") {
        gen_merge(&context)?;
        gen_mr_down_sampling(&context)?;
        gen_mr_unitigs(&context)?;
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

fn gen_repetitive(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/1_repetitive.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/1_repetitive.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

    Ok(())
}

fn gen_unitigs(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/4_unitigs_multik.sh";
    eprintln!("Create {}", outname);

    let mut con = Context::new();
    con.insert("outname", outname);
    con.extend(context.clone());

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/4_unitigs.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &con).unwrap();
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

    Ok(())
}

fn gen_mr_unitigs(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/6_unitigs_multik.sh";
    eprintln!("Create {}", outname);

    let mut con = Context::new();
    con.insert("outname", outname);
    con.extend(context.clone());

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/6_unitigs.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", &con).unwrap();
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

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
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

    Ok(())
}

fn gen_bsub(context: &Context) -> anyhow::Result<()> {
    let outname = "0_script/0_bsub.sh";
    eprintln!("Create {}", outname);

    let mut tera = Tera::default();
    tera.add_raw_templates(vec![
        ("header", include_str!("../../templates/header.tera.sh")),
        ("t", include_str!("../../templates/0_bsub.tera.sh")),
    ])
    .unwrap();

    let rendered = tera.render("t", context).unwrap();
    anchr::utils::write_lines(outname, &[rendered.as_str()])?;

    Ok(())
}

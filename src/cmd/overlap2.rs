use clap::*;
use cmd_lib::*;
use std::collections::HashMap;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("overlap2")
        .about("Detect overlaps between two (large) files by daligner")
        .after_help(
            r###"
* Since `daligner` cannot perform overlap comparisons between two databases,
  the current strategy is to place the two FASTA files into a single database
  and then split that database into multiple chunks. The chunk containing the
  first sequence will be compared against the other chunks. This approach can
  still significantly reduce the computational workload.

* The inputs can't be stdin
* All intermediate files (.fasta, .replace.tsv, .db, .las, .show.txt, .ovlp.tsv)
  are kept in the working directory.

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Sets the first input file"),
        )
        .arg(
            Arg::new("infile2")
                .required(true)
                .num_args(1)
                .index(2)
                .help("Sets the second input file"),
        )
        .arg(
            Arg::new("dir")
                .long("dir")
                .short('d')
                .num_args(1)
                .default_value(".")
                .help("Change working directory"),
        )
        .arg(
            Arg::new("p1")
                .long("p1")
                .num_args(1)
                .default_value("anchor")
                .help("Prefix of the first"),
        )
        .arg(
            Arg::new("p2")
                .long("p2")
                .num_args(1)
                .default_value("long")
                .help("Prefix of the second"),
        )
        .arg(
            Arg::new("pd")
                .long("pd")
                .num_args(1)
                .help("Prefix of the result files"),
        )
        .arg(
            Arg::new("block")
                .long("block")
                .short('b')
                .num_args(1)
                .default_value("20")
                .value_parser(value_parser!(usize))
                .help("Block size in Mbp"),
        )
        .arg(
            Arg::new("len")
                .long("len")
                .short('l')
                .num_args(1)
                .default_value("500")
                .value_parser(value_parser!(i32))
                .help("minimal length of overlaps"),
        )
        .arg(
            Arg::new("idt")
                .long("idt")
                .short('i')
                .num_args(1)
                .default_value("0.7")
                .value_parser(value_parser!(f32))
                .help("minimal identities of overlaps"),
        )
        .arg(
            Arg::new("all")
                .long("all")
                .action(ArgAction::SetTrue)
                .help("All overlaps instead of proper ones"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("8")
                .value_parser(value_parser!(usize))
                .help("Number of threads"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Args
    //----------------------------
    let opt_p1 = args.get_one::<String>("p1").unwrap();
    let opt_p2 = args.get_one::<String>("p2").unwrap();

    let opt_pd = if args.contains_id("pd") {
        args.get_one::<String>("pd").unwrap()
    } else {
        &format!("{}{}", opt_p1, anchr::ucfirst(opt_p2))
    };

    let opt_block = *args.get_one::<usize>("block").unwrap();
    let opt_min_len = *args.get_one::<i32>("len").unwrap();
    let opt_min_idt = *args.get_one::<f32>("idt").unwrap();

    let is_all = args.get_flag("all");

    let opt_parallel = *args.get_one::<usize>("parallel").unwrap();

    let curdir = std::env::current_dir()?;
    let anchr = std::env::current_exe().unwrap().display().to_string();

    let workdir = std::path::Path::new(args.get_one::<String>("dir").unwrap()).to_path_buf();
    let workdir_str = workdir.to_str().unwrap();
    std::fs::create_dir_all(workdir_str)?;

    run_cmd!(info "==> Paths")?;
    run_cmd!(info "    \"anchr\"   = ${anchr}")?;
    run_cmd!(info "    \"curdir\"  = ${curdir}")?;
    run_cmd!(info "    \"workdir\" = ${workdir_str}")?;

    //----------------------------
    // Paths
    //----------------------------
    run_cmd!(info "==> Absolute paths")?;
    let file1 = intspan::absolute_path(args.get_one::<String>("infile").unwrap())
        .unwrap()
        .display()
        .to_string();
    let file2 = intspan::absolute_path(args.get_one::<String>("infile2").unwrap())
        .unwrap()
        .display()
        .to_string();

    //----------------------------
    // Ops
    //----------------------------
    run_cmd!(info "==> Switch to workdir")?;
    std::env::set_current_dir(workdir_str)?;

    run_cmd!(info "==> Preprocess the first file")?;
    run_cmd!(
        ${anchr} dazzname ${file1} --prefix ${opt_p1} -o ${opt_p1}.fasta
    )?;
    let first_sum = run_fun!(
        hnsm n50 -H -N 0 -S ${opt_p1}.fasta
    )
    .unwrap()
    .parse::<usize>()
    .unwrap();
    let first_count = run_fun!(
        hnsm n50 -H -N 0 -C ${opt_p1}.fasta
    )
    .unwrap()
    .parse::<usize>()
    .unwrap();

    run_cmd!(info "==> Preprocess the second file")?;
    {
        let second_start = first_count + 1;
        run_cmd!(
            hnsm filter -u -a ${opt_min_len} ${file2} |
                ${anchr} dazzname stdin --prefix ${opt_p2} --start ${second_start} -o ${opt_p2}.fasta
        )?;
    }

    run_cmd!(info "==> Make the dazzler DB, each block is of size ${opt_block} MB")?;
    {
        if anchr::file_exists(&std::env::current_dir()?, &format!("{}.db", opt_pd))
            || anchr::file_exists(&std::env::current_dir()?, &format!(".{}.bps", opt_pd))
        {
            run_cmd!(info "    Remove existing DB")?;
            run_cmd!(
                DBrm ${opt_pd}
            )?;
        }

        run_cmd!(
            fasta2DB ${opt_pd} ${opt_p1}.fasta;
            fasta2DB ${opt_pd} ${opt_p2}.fasta;
            DBdust ${opt_pd};
            DBsplit -f -s${opt_block} ${opt_pd};
        )?;
    }

    run_cmd!(info "==> Run daligner")?;
    {
        let block_number = run_fun!(
            cat ${opt_pd}.db |
                perl -n -e r#"/^blocks\s+=\s+(\d+)\s*$/ and print $1"#
        )
        .unwrap()
        .parse::<usize>()
        .unwrap();
        let first_idx = (first_sum as f64 / 1_000_000.0 / opt_block as f64 + 1.0).floor() as usize;

        run_cmd!(info "    \"first_sum\"    = ${first_sum}")?;
        run_cmd!(info "    \"first_count\"  = ${first_count}")?;
        run_cmd!(info "    \"first_idx\"    = ${first_idx}")?;
        run_cmd!(info "    \"block_number\" = ${block_number}")?;
        run_cmd!(info "    \"block_size\"   = ${opt_block}")?;

        if anchr::file_exists(&std::env::current_dir()?, &format!("{}.las", opt_pd))
            || anchr::file_exists(&std::env::current_dir()?, &format!(".{}.1.las", opt_pd))
        {
            run_cmd!(info "    Remove existing alignments")?;
            run_cmd!(
                fd -g "${opt_pd}*.las" -x rm
            )?;
        }

        // Don't use HPC.daligner as we want to avoid all-vs-all comparisons.
        // HPC.daligner tries to give every sequences the same change to match with others.
        for i in 1..=first_idx {
            // Start from $i instead of $first_idx for conveniences of LAmerge
            for j in i..=block_number {
                run_cmd!(
                    daligner ${opt_pd}.${i}  ${opt_pd}.${j} -M16 -T${opt_parallel} -e${opt_min_idt} -l${opt_min_len} -s${opt_min_len} -mdust;
                    LAcheck -S ${opt_pd} ${opt_pd}.${i}.${opt_pd}.${j};
                    LAcheck -S ${opt_pd} ${opt_pd}.${j}.${opt_pd}.${i};
                )?;
            }
        }

        for i in 1..=first_idx {
            run_cmd!(
                LAmerge ${opt_pd}.$i ${opt_pd}.${i}.${opt_pd}.@;
                LAcheck -S ${opt_pd} ${opt_pd}.${i};
            )?;
        }
        run_cmd!(
            fd -g "${opt_pd}.*.${opt_pd}.*.las" -x rm
        )?;

        run_cmd!(
            LAcat ${opt_pd}.@.las > ${opt_pd}.las
        )?;
        run_cmd!(
            fd -g "${opt_pd}.*.las" -x rm
        )?;
    }

    run_cmd!(info "==> Outputs")?;
    if is_all {
        run_cmd!(
            LAshow ${opt_pd}.db ${opt_pd}.las > ${opt_pd}.show.txt
        )?;
    } else {
        run_cmd!(
            LAshow -o ${opt_pd}.db ${opt_pd}.las > ${opt_pd}.show.txt
        )?;
    }

    run_cmd!(
        cat ${opt_p1}.fasta.replace.tsv ${opt_p2}.fasta.replace.tsv |
            ${anchr} show2ovlp ${opt_pd}.show.txt stdin --orig -o ${opt_pd}.ovlp.tsv
    )?;
    // anchr::pause();

    //----------------------------
    // Done
    //----------------------------
    std::env::set_current_dir(&curdir)?;

    Ok(())
}

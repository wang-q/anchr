use crate::libs::asm::assemble::{assemble_unitigs, AssembleOptions};
use anyhow::Context;
use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for unitig.
pub fn make_subcommand() -> Command {
    Command::new("unitig")
        .about("Assembles reads into maximal unitigs (non-branching paths)")
        .after_help(
            r###"
This command assembles reads into maximal unitigs through the k-mer graph,
following the BCALM 2 compaction semantics (GATB `ograph.cpp` `graph3`):
every solid k-mer (count >= 3) extends in both directions only while it has
exactly one solid successor whose own predecessor is also unique, so the
assembly stops at branches, junctions, coverage gaps, and loops. Parallel
paths stay separate (no bubble popping), and the result is independent of
the k-mer scan order.

This is the strict graph-compression counterpart of `anchr asm contig`, whose
seeded contig mode keeps extending through weak branches (tadpole-compatible
behavior). Unitigs are best suited to high-coverage or error-corrected input,
such as the anchr `unitigs` step's `pe.cor.fa`.

Notes:
* Input is one or more FASTA/FASTQ files (plain or gzipped); pairing is
  irrelevant for assembly (BCALM semantics), and `--list-files` reads a
  one-path-per-line list
* Unitigs are written longest-first with a `unitig_<id>` FASTA header
  carrying length, coverage, GC, and dimer composition fields
* Processing is ordered and deterministic, independent of scan order
* Output sequences are wrapped at 70 columns, like BBTools FASTA output
* Supports both plain text and gzipped (.gz) files

Examples:
1. Assemble unitigs from corrected reads:
   anchr asm unitig pe.cor.fa -o unitigs_K31.fasta --kmer 31

2. Assemble from paired-end reads:
   anchr asm unitig R1.fq.gz R2.fq.gz -o unitigs.fasta

3. Raise the solid k-mer threshold (like bcalm `-abundance-min`):
   anchr asm unitig in.fq -o out.fasta --min-count-seed 5

4. Assemble from a list of files and emit every k-mer abundance:
   anchr asm unitig files.list -o out.fasta --list-files \
       --all-abundance-counts
"###,
        )
        .arg(crate::cmd::args::infiles_arg_with_numargs(
            "Input file(s): FASTA/FASTQ, plain or gzipped; use --list-files for a one-path-per-line list",
            1..,
        ))
        .arg(crate::cmd::args::outfile_arg())
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("31")
                .value_parser(clap::builder::RangedU64ValueParser::<usize>::new().range(1..))
                .help("K-mer length (1..=256)"),
        )
        .arg(
            Arg::new("min_contig_len")
                .long("min-contig-len")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Minimum unitig length (default: keep all unitigs, matching bcalm's lossless compaction)"),
        )
        .arg(
            Arg::new("min_count_seed")
                .long("min-count-seed")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Solid k-mer count threshold (default 3, like bcalm -abundance-min)"),
        )
        .arg(
            Arg::new("links")
                .long("links")
                .action(ArgAction::SetTrue)
                .help("Append L: links to unitig FASTA headers (bcalm format)"),
        )
        .arg(
            Arg::new("gfa")
                .long("gfa")
                .action(ArgAction::SetTrue)
                .help("Emit a GFA graph instead of FASTA"),
        )
        .arg(
            Arg::new("list_files")
                .long("list-files")
                .action(ArgAction::SetTrue)
                .help("Treat infiles as list files, one sequence file path per line"),
        )
        .arg(
            Arg::new("all_abundance_counts")
                .long("all-abundance-counts")
                .action(ArgAction::SetTrue)
                .help("Emit every k-mer abundance in FASTA headers (ab:Z:, like bcalm -all-abundance-counts)"),
        )
        .arg(
            Arg::new("dfa")
                .long("dfa")
                .action(ArgAction::SetTrue)
                .help("Use the DFA-state walk engine (default; kept for compatibility)"),
        )
        .arg(
            Arg::new("no_dfa")
                .long("no-dfa")
                .action(ArgAction::SetTrue)
                .help("Disable the DFA-state walk engine (use the bucket-scan walk)"),
        )
        .arg(
            Arg::new("supermer")
                .long("supermer")
                .action(ArgAction::SetTrue)
                .help("Use FastK-style super-mer two-stage counting (default for FASTA; kept for compatibility)"),
        )
        .arg(
            Arg::new("no_supermer")
                .long("no-supermer")
                .action(ArgAction::SetTrue)
                .help("Disable super-mer counting (use the direct streamed counter)"),
        )
        .arg(
            Arg::new("supermer_m")
                .long("supermer-m")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Minimizer length for --supermer (default: auto min(12, max(5, k/4)))"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("half")
                .help("Worker threads for counting + DFA classification (default: half of logical cores, capped at 8; auto = all cores)"),
        )
}

/// Execute the unitig command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let is_list = args.get_flag("list_files");
    let mut infiles: Vec<String> = Vec::new();
    for f in args.get_many::<String>("infiles").unwrap() {
        infiles.extend(pgr::libs::par::resolve_paths(f, is_list)?);
    }
    anyhow::ensure!(
        !infiles.is_empty(),
        "--list-files resolved to no input files"
    );
    let outfile = crate::cmd::args::get_outfile(args);
    // Reject `-o` that would overwrite an input file (the writer is opened
    // before the reads are consumed).
    crate::cmd::args::ensure_outfile_distinct(outfile, infiles.iter().map(|s| s.as_str()))?;
    let parallel = crate::cmd::args::parse_parallel(args.get_one::<String>("parallel").unwrap())?;
    let k = *args.get_one::<usize>("kmer").unwrap();
    let use_supermer = !args.get_flag("no_supermer");
    let opts = AssembleOptions {
        k,
        min_contig_len: args
            .get_one::<usize>("min_contig_len")
            .copied()
            .unwrap_or(0),
        min_count_seed: args
            .get_one::<usize>("min_count_seed")
            .copied()
            .unwrap_or(3),
        emit_links: args.get_flag("links"),
        emit_gfa: args.get_flag("gfa"),
        all_abundance_counts: args.get_flag("all_abundance_counts"),
        use_dfa: !args.get_flag("no_dfa"),
        use_supermer,
        supermer_m: args
            .get_one::<usize>("supermer_m")
            .copied()
            .or_else(|| use_supermer.then(|| (12).min((5).max(k / 4)))),
        parallel,
        ..AssembleOptions::default()
    };

    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    // One rayon pool for the whole pipeline (counting + DFA classification);
    // the streamed counting workers sort sequentially and never spawn rayon.
    let stats = if opts.parallel > 0 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(opts.parallel)
            .build()
            .with_context(|| "failed to build assembly thread pool")?;
        pool.install(|| assemble_unitigs(&infiles, &mut out, &opts))?
    } else {
        assemble_unitigs(&infiles, &mut out, &opts)?
    };
    out.flush()?;
    eprintln!(
        "Reads in: {}  Unitigs: {}  Bases: {}  Longest: {}",
        stats.reads_in, stats.contigs_built, stats.bases_built, stats.longest_contig
    );
    Ok(())
}

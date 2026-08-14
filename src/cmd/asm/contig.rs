use crate::libs::asm::assemble::{assemble, AssembleOptions};
use anyhow::Context;
use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for contig.
pub fn make_subcommand() -> Command {
    Command::new("contig")
        .about("Assembles reads into contigs via k-mer graph traversal (tadpole-compatible)")
        .after_help(
            r###"
This command assembles reads into contigs through the k-mer graph, reproducing
the BBTools `tadpole.sh` contig mode (the default mode when no `ecc`/`extend`
flag is set): k-mers are counted with a quality gate (`--min-prob`), contigs
are seeded from k-mers above a depth threshold and extended greedily in both
directions, stopping at branches and dead ends. This replaces the tadpole
assembly steps of the anchr `2_insert_size` and `unitigs` flows.

Notes:
* Input is one or more FASTA/FASTQ files (plain or gzipped); pairing is
  irrelevant for assembly, and `--list-files` reads a one-path-per-line list
* Contigs are written longest-first with a `contig_<id>` FASTA header carrying
  length, coverage, GC, and dimer composition fields (BBTools SHORT_NAMES)
* Processing is ordered and deterministic (equivalent to `threads=1`)
* Bubble-popping resolutions may differ slightly from BBTools on some
  overlapping structures (its expand order depends on a memory-dependent
  hash layout); the contig set and total bases match, and the output is
  reproducible across runs
* Bubble popping is on by default (tadpole `popbubbles=t`); pass
  `--no-bubbles` to keep parallel-path contigs separate (tadpole
  `popbubbles=f`)
* Output sequences are wrapped at 70 columns, like BBTools FASTA output
* Supports both plain text and gzipped (.gz) files

Examples:
1. Assemble contigs from corrected reads:
   anchr asm contig pe.cor.fa -o unitigs_K31.fasta --kmer 31

2. Assemble from paired-end reads (anchr 2_insert_size step):
   anchr asm contig R1.fq.gz R2.fq.gz -o contigs.fasta

3. Raise the minimum contig length:
   anchr asm contig in.fq -o out.fasta --min-contig-len 500

4. Raise the seeding depth threshold (tadpole `mincountseed`):
   anchr asm contig in.fq -o out.fasta --min-count-seed 5

5. Drop low-coverage contigs (tadpole `mincoverage`):
   anchr asm contig in.fq -o out.fasta --min-coverage 5

6. Assemble from a list of files:
   anchr asm contig files.list -o contigs.fasta --list-files
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
                .help("K-mer length (1..=128)"),
        )
        .arg(
            Arg::new("min_contig_len")
                .long("min-contig-len")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Minimum contig length (default: max(124, 2*k))"),
        )
        .arg(
            Arg::new("min_count_seed")
                .long("min-count-seed")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Minimum k-mer depth to seed a contig (tadpole mincountseed, default 3)"),
        )
        .arg(
            Arg::new("min_coverage")
                .long("min-coverage")
                .num_args(1)
                .value_parser(value_parser!(f32))
                .help(
                    "Minimum mean k-mer coverage for a contig (tadpole mincoverage, default 1.0)",
                ),
        )
        .arg(
            Arg::new("no_bubbles")
                .long("no-bubbles")
                .action(ArgAction::SetTrue)
                .help("Keep parallel-path contigs separate (disable bubble popping)"),
        )
        .arg(
            Arg::new("list_files")
                .long("list-files")
                .action(ArgAction::SetTrue)
                .help("Treat infiles as list files, one sequence file path per line"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("half")
                .help("Worker threads for k-mer counting (default: half of logical cores, capped at 8; auto = all cores); walk stays deterministic"),
        )
}

/// Execute the contig command.
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
    // Validate the thread-count value; processing stays deterministic
    // single-pass (see the design notes), so the result is not used.
    let parallel = crate::cmd::args::parse_parallel(args.get_one::<String>("parallel").unwrap())?;
    let k = *args.get_one::<usize>("kmer").unwrap();
    let opts = AssembleOptions {
        k,
        min_contig_len: args
            .get_one::<usize>("min_contig_len")
            .copied()
            // tadpole `mincontiglen` auto default (kept: contig mode is
            // tadpole-compatible; unitig mode deliberately has no filter)
            .unwrap_or_else(|| (124).max(2 * k)),
        min_count_seed: args
            .get_one::<usize>("min_count_seed")
            .copied()
            .unwrap_or(3),
        min_coverage: args.get_one::<f32>("min_coverage").copied().unwrap_or(1.0),
        pop_bubbles: !args.get_flag("no_bubbles"),
        parallel,
        ..AssembleOptions::default()
    };

    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    let stats = if opts.parallel > 0 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(opts.parallel)
            .build()
            .with_context(|| "failed to build assembly thread pool")?;
        pool.install(|| assemble(&infiles, &mut out, &opts))?
    } else {
        assemble(&infiles, &mut out, &opts)?
    };
    out.flush()?;
    eprintln!(
        "Reads in: {}  Contigs: {}  Bases: {}  Longest: {}",
        stats.reads_in, stats.contigs_built, stats.bases_built, stats.longest_contig
    );
    Ok(())
}

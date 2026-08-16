use crate::libs::asm::multik::{assemble_multik, MultikOptions};
use anyhow::Context;
use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for multik.
pub fn make_subcommand() -> Command {
    Command::new("multik")
        .about("Assembles reads into long unitigs by iterating over increasing k")
        .after_help(
            r###"
Iterates the unitig graph over increasing k-mer lengths (metaMDBG-style
cross-round validation): pass 0 assembles maximal unitigs at the first k,
then every later k validates the previous graph — the bridge k-mer covering
each unitig junction must be solid (count >= --min-count-extend) in the
reads plus the previous unitigs, and every internal k-mer of a long-enough
unitig must stay solid (chimeric-unitig cleanup). A progressive abundance
filter then drops the lowest-abundance unitigs (cutoff grows ~10% per round
from 1.1 up to the graph maximum) and recompacts surviving chains — only
branching junctions and isolated nodes are pruned by abundance, unique-chain
unitigs (the main path) are always kept, so coverage fluctuations within one
genome do not break the assembly. High-coverage strain divergence dissolves,
the main path compacts into longer unitigs, and dropped unitigs stay
independent output (low-abundance species are not lost). See
notes/design/asm-multik.md.

This is the iterative counterpart of `anchr asm olc` (parallel multi-k
pooling + heuristic layout): here each junction is accepted only when the
larger k's k-mer count supports it, which is how bubble branches are chosen
without heuristics.

Notes:
* Input is one or more FASTA/FASTQ files (plain or gzipped); pairing is
  irrelevant for assembly, and `--list-files` reads a one-path-per-line list
* k-mer lengths are sorted and deduplicated internally; the smallest k is
  the master k (pass 0 skeleton), each larger k is a slave k that validates
  the graph. Run one master per invocation (`-k 31,41,51,61,71,81` uses
  master 31; `-k 51,61,71,81` uses master 51) and let the template drive
  several masters in parallel, then merge their unitigs
* `--kmer auto` (default) derives the sequence from the read-length N50
  (~1/3 of the read length up to 51, e.g. 50/70/90/110 for 150 bp reads,
  51/81/111 for long reads)
* Output unitigs are written longest-first with `unitig_<id>` FASTA headers
  carrying length and coverage

Examples:
1. Iterate over three k values:
   anchr asm multik reads.fq.gz -o unitigs.fa --kmer 21,51,81
2. Raise the validation threshold:
   anchr asm multik reads.fa -o unitigs.fa --min-count-extend 3
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
                .default_value("auto")
                .help("Comma-separated increasing k-mer lengths (1..=256), or auto (default) to derive from the read-length N50"),
        )
        .arg(
            Arg::new("min_count_seed")
                .long("min-count-seed")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Solid k-mer count threshold for pass 0 (default 3)"),
        )
        .arg(
            Arg::new("min_count_extend")
                .long("min-count-extend")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Solid k-mer count threshold for cross-round validation (default 2)"),
        )
        .arg(
            Arg::new("merge_similar")
                .long("merge-similar")
                .num_args(1)
                .value_parser(value_parser!(f64))
                .help("Bubble merge: minimum similarity of collapsed alternative paths (default 0.95)"),
        )
        .arg(
            Arg::new("merge_len")
                .long("merge-len")
                .num_args(1)
                .value_parser(value_parser!(usize))
                .help("Bubble merge: maximum alternative-path length as k multiples (default 20)"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Worker threads for counting; 0 = all cores"),
        )
        .arg(
            Arg::new("guide_contigs")
                .long("guide-contigs")
                .num_args(1)
                .help(
                    "Previous-master unitigs for megahit-style guidance: each contig's full \
                     sequence feeds the master-k count as pseudo-reads (repeated to the solid \
                     threshold), so a low-k master's structure carries into higher-k rounds \
                     (e.g. K192 on 450 bp reads: unitig N50 37.6K -> 81.6K)",
                ),
        )
        .arg(
            Arg::new("print_ks")
                .long("print-ks")
                .action(ArgAction::SetTrue)
                .help(
                    "Print the auto-derived master-k sequence (from the read-length N50) \
                     and exit; lets templates drive per-master runs with an adaptive k \
                     list instead of hard-coded values",
                ),
        )
        .arg(
            Arg::new("list_files")
                .long("list-files")
                .action(ArgAction::SetTrue)
                .help("Treat infiles as list files, one sequence file path per line"),
        )
}

/// Execute the multik command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    if args.get_flag("print_ks") {
        let infiles: Vec<String> = args
            .get_many::<String>("infiles")
            .unwrap()
            .flat_map(|f| pgr::libs::par::resolve_paths(f, false).unwrap_or_default())
            .collect();
        let ks = crate::libs::asm::multik::auto_ks_for_reads(&infiles)
            .with_context(|| "failed to derive the k sequence from the reads")?;
        println!(
            "{}",
            ks.iter()
                .map(|k| k.to_string())
                .collect::<Vec<_>>()
                .join(" ")
        );
        return Ok(());
    }
    let is_list = args.get_flag("list_files");
    let mut infiles: Vec<String> = Vec::new();
    for f in args.get_many::<String>("infiles").unwrap() {
        infiles.extend(pgr::libs::par::resolve_paths(f, is_list)?);
    }
    let min_count_seed = args
        .get_one::<usize>("min_count_seed")
        .copied()
        .unwrap_or(3);
    // megahit `seq2sdbg --contig` guidance: the previous master's unitigs
    // seed the current master-k count as pseudo-reads (each full-length
    // contig repeated up to the solid threshold), so low-k structure
    // carries into higher-k rounds that would otherwise fragment from
    // insufficient read support.
    // The temp dir stays alive (as `_guide_keep`) until the assembly
    // finishes, so the pseudo-read file exists while it is read.
    let _guide_keep = if let Some(guide) = args.get_one::<String>("guide_contigs") {
        let (dir, guide_path) = write_guide_pseudo_reads(guide, min_count_seed)?;
        infiles.push(guide_path);
        Some(dir)
    } else {
        None
    };
    anyhow::ensure!(
        !infiles.is_empty(),
        "--list-files resolved to no input files"
    );
    let kmer_arg = args.get_one::<String>("kmer").unwrap();
    let ks: Vec<usize> = if kmer_arg.trim().eq_ignore_ascii_case("auto") {
        Vec::new() // auto-derive from the read-length N50 inside the lib
    } else {
        kmer_arg
            .split(',')
            .map(|s| {
                s.trim()
                    .parse::<usize>()
                    .with_context(|| format!("invalid --kmer value: {s}"))
            })
            .collect::<anyhow::Result<_>>()?
    };
    let opts = MultikOptions {
        ks,
        min_count_seed,
        min_count_extend: args
            .get_one::<usize>("min_count_extend")
            .copied()
            .unwrap_or(2),
        merge_similar: args
            .get_one::<f64>("merge_similar")
            .copied()
            .unwrap_or(0.95),
        merge_len: args.get_one::<usize>("merge_len").copied().unwrap_or(20),
        parallel: *args.get_one::<usize>("parallel").unwrap(),
    };
    let outfile = crate::cmd::args::get_outfile(args);
    crate::cmd::args::ensure_outfile_distinct(outfile, infiles.iter().map(|s| s.as_str()))?;
    let unitigs = assemble_multik(&infiles, &opts)?;
    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    for (i, u) in unitigs.iter().enumerate() {
        writeln!(
            out,
            ">unitig_{} len={},cov={}",
            i + 1,
            u.bases.len(),
            super::common::format_cov(u.coverage as f64)
        )?;
        for chunk in u.bases.chunks(70) {
            out.write_all(chunk)?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(())
}

/// Writes a previous master's unitigs as full-length pseudo-reads (each
/// contig repeated up to the solid threshold) and returns the temp dir
/// holding the file plus its path. Callers must keep the temp dir alive
/// until the assembly consumes the infiles.
fn write_guide_pseudo_reads(
    guide: &str,
    min_count_seed: usize,
) -> anyhow::Result<(tempfile::TempDir, String)> {
    use anyhow::Context;
    let contigs = crate::libs::map::read_fasta(std::slice::from_ref(&guide.to_string()))
        .with_context(|| format!("failed to read guide contigs {guide}"))?;
    let dir =
        tempfile::tempdir().with_context(|| "failed to create temp dir for guide pseudo-reads")?;
    let guide_file = dir.path().join("guide.fa");
    let mut out = pgr::libs::io::writer(guide_file.to_str().unwrap()).with_context(|| {
        format!(
            "failed to open guide pseudo-read file {}",
            guide_file.display()
        )
    })?;
    for (i, r) in contigs.iter().enumerate() {
        for rep in 0..min_count_seed {
            writeln!(out, ">g{i}_{rep}")?;
            for chunk in r.seq.chunks(70) {
                out.write_all(chunk)?;
                out.write_all(b"\n")?;
            }
        }
    }
    out.flush()?;
    Ok((dir, guide_file.to_string_lossy().into_owned()))
}

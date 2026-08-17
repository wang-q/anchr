use crate::libs::asm::assemble::{assemble_unitigs_buf, AssembleOptions};
use crate::libs::olc::consensus::consensus_with_ratio;
use crate::libs::olc::layout::build_layouts;
use crate::libs::olc::overlap::{
    drop_cross_chimeras, filter_contained, find_overlaps, CrossOptions, OverlapOptions, Unitig,
};
use anyhow::Context;
use clap::{value_parser, Arg, ArgAction, ArgMatches, Command};
use std::io::Write;
use std::path::Path;

/// Build the clap subcommand for olc.
pub fn make_subcommand() -> Command {
    Command::new("olc")
        .about("Assembles reads into contigs via multi-k unitig OLC")
        .after_help(
            r###"
Runs the full OLC pipeline in memory: for every k in --kmer the reads are
assembled into maximal unitigs (`anchr asm unitig` semantics), all unitigs are
pooled as pseudo-reads, exact overlaps are found (`anchr asm ovlp`), layouts
are built greedily (`anchr asm layout`), and each layout is stitched into a
consensus contig (`anchr asm cns`). With `--unitigs` the inputs are taken
as unitigs/contigs directly (no re-assembly), which is the intended path for
merging pre-assembled sets (e.g. per-coverage `asm multik` outputs); the
contained-unitig filter still runs, so overlapping sets are deduplicated.
See notes/design/asm-olc.md.

Unitigs are named `k<k>:unitig_<id>` so the per-k sets stay distinguishable
and reproducible. Overlaps are exact (error-free unitigs), layouts stop at
ambiguous junctions and non-reciprocal edges, and no bubble heuristics are
applied.

Notes:
* Input is one or more FASTA/FASTQ files (plain or gzipped); pairing is
  irrelevant for assembly, and `--list-files` reads a one-path-per-line list
* `--unitigs` treats the inputs as already-assembled sequences (one FASTA per
  file; names get the file stem as a tag, so separate files stay
  distinguishable — do not concatenate them into one file); `--kmer` and
  `--min-count-seed` are ignored in this mode
* --keep-dir writes the intermediate unitigs/overlap/layout files for
  debugging and inspection; the names there omit the `stem:` prefix that
  the standalone ovlp/layout/cns commands derive, so they are not directly
  re-runnable through those commands as-is
* Output contigs are written longest-first with `>contig_<id> len=...,cov=...`
  headers, 70-column wrapped

Examples:
1. Assemble a small metagenome with three k values:
   anchr asm olc reads.fq.gz -o contigs.fa --kmer 21,51,81
2. Keep the intermediates and raise the minimum contig length:
   anchr asm olc R1.fq.gz R2.fq.gz -o contigs.fa \
       --kmer 21,51,81 --min-contig-len 1000 --keep-dir stage/
3. Assemble from a list of files:
   anchr asm olc files.list -o contigs.fa --kmer 21,51,81 --list-files
4. Merge pre-assembled per-coverage unitigs (no re-assembly):
   anchr asm olc k40.fa k80.fa -o contigs.fa --unitigs
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
                .default_value("21,51,81")
                .help("Comma-separated k-mer lengths (1..=256) for the unitig sets"),
        )
        .arg(
            Arg::new("min_count_seed")
                .long("min-count-seed")
                .num_args(1)
                .default_value("3")
                .value_parser(value_parser!(usize))
                .help("Solid k-mer count threshold for unitig assembly"),
        )
        .arg(
            Arg::new("overlap_k")
                .long("overlap-k")
                .num_args(1)
                .default_value("17")
                .value_parser(value_parser!(usize))
                .help("Seed k-mer length for overlap detection (1..=256)"),
        )
        .arg(
            Arg::new("min_overlap")
                .long("min-overlap")
                .num_args(1)
                .default_value("34")
                .value_parser(value_parser!(usize))
                .help("Minimum accepted overlap length in bases"),
        )
        .arg(
            Arg::new("min_contig_len")
                .long("min-contig-len")
                .num_args(1)
                .default_value("500")
                .value_parser(value_parser!(usize))
                .help("Minimum output contig length in bases"),
        )
        .arg(
            Arg::new("unitigs")
                .long("unitigs")
                .action(ArgAction::SetTrue)
                .help("Inputs are already-assembled unitigs/contigs (skip the S0 unitig assembly; one FASTA per file, tagged by file stem)"),
        )
        .arg(
            Arg::new("cross_validate")
                .long("cross-validate")
                .action(ArgAction::SetTrue)
                .requires("unitigs")
                .help("With multiple input files: drop contigs whose two ends are each covered by >= 2 other files while no other file spans the middle junction (single-file chimeric joins)"),
        )
        .arg(
            Arg::new("keep_dir")
                .long("keep-dir")
                .num_args(1)
                .help("Directory for intermediate unitigs/ovlp/layout files (for inspection)"),
        )
        .arg(
            Arg::new("list_files")
                .long("list-files")
                .action(ArgAction::SetTrue)
                .help("Treat infiles as list files, one sequence file path per line"),
        )
}

/// Execute the olc command.
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
    let ks: Vec<usize> = args
        .get_one::<String>("kmer")
        .unwrap()
        .split(',')
        .map(|s| {
            s.trim()
                .parse::<usize>()
                .with_context(|| format!("invalid --kmer value: {s}"))
        })
        .collect::<anyhow::Result<_>>()?;
    anyhow::ensure!(!ks.is_empty(), "at least one k-mer length is required");
    let min_count_seed = *args.get_one::<usize>("min_count_seed").unwrap();
    let seed_k = *args.get_one::<usize>("overlap_k").unwrap();
    let min_overlap = *args.get_one::<usize>("min_overlap").unwrap();
    let min_contig_len = *args.get_one::<usize>("min_contig_len").unwrap();
    let keep_dir = args.get_one::<String>("keep_dir");
    let outfile = crate::cmd::args::get_outfile(args);
    // Reject `-o` that would overwrite an input read file.
    crate::cmd::args::ensure_outfile_distinct(outfile, infiles.iter().map(|s| s.as_str()))?;

    // S0: unitigs per k, or take the inputs directly in --unitigs mode.
    let mut unitigs = Vec::new();
    if args.get_flag("unitigs") {
        unitigs = super::common::read_unitigs(&infiles)?;
    } else {
        for &k in &ks {
            let opts = AssembleOptions {
                k,
                min_count_seed,
                ..AssembleOptions::default()
            };
            let bufs = assemble_unitigs_buf(&infiles, &opts)?;
            for (id, bases) in bufs {
                unitigs.push(Unitig {
                    name: format!("k{k}:unitig_{id}"),
                    seq: bases,
                });
            }
        }
    }
    // S1: exact overlaps.
    let overlaps = find_overlaps(
        &unitigs,
        &OverlapOptions {
            seed_k,
            min_overlap,
        },
    )?;

    // S1.25: cross-sample validation. Must run before `filter_contained`:
    // a single-file chimeric join prefix-contains the other files' correct
    // contigs, so the plain containment filter would drop the correct
    // versions and keep the chimera.
    let (unitigs, overlaps) = if args.get_flag("cross_validate") {
        drop_cross_chimeras(
            &unitigs,
            &overlaps,
            &CrossOptions {
                flank: min_overlap,
                span: min_overlap / 2,
                min_groups: 2,
            },
        )
    } else {
        (unitigs, overlaps)
    };

    // S1.5: drop unitigs fully contained in longer unitigs (multi-k
    // redundancy); layouts are unchanged, the graph shrinks.
    let (unitigs, overlaps) = filter_contained(&unitigs, &overlaps);
    if let Some(dir) = keep_dir {
        std::fs::create_dir_all(dir)
            .with_context(|| format!("failed to create --keep-dir {dir}"))?;
        dump_unitigs(dir, &unitigs)?;
        dump_paf(dir, &unitigs, &overlaps)?;
    }

    // S2: greedy layouts.
    let layouts = build_layouts(&unitigs, &overlaps)?;
    if let Some(dir) = keep_dir {
        dump_layouts(dir, &unitigs, &layouts)?;
    }

    // S3: consensus stitch.
    let dedup_ratio = if args.get_flag("unitigs") { 0.99 } else { 1.0 };
    let contigs = consensus_with_ratio(&unitigs, &layouts, min_contig_len, dedup_ratio)?;
    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    for (i, c) in contigs.iter().enumerate() {
        writeln!(
            out,
            ">contig_{} len={},cov={}",
            i + 1,
            c.seq.len(),
            super::common::format_cov(c.coverage)
        )?;
        for chunk in c.seq.chunks(70) {
            out.write_all(chunk)?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(())
}

/// Writes the pooled unitigs as one FASTA (k-tagged names).
fn dump_unitigs(dir: &str, unitigs: &[Unitig]) -> anyhow::Result<()> {
    let path = Path::new(dir).join("unitigs.fa");
    let path = path.to_str().unwrap();
    let mut out = pgr::libs::io::writer(path).with_context(|| format!("failed to open {path}"))?;
    for u in unitigs {
        writeln!(out, ">{}", u.name)?;
        for chunk in u.seq.chunks(70) {
            out.write_all(chunk)?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(())
}

/// Writes the overlap PAF.
fn dump_paf(
    dir: &str,
    unitigs: &[Unitig],
    overlaps: &[crate::libs::olc::overlap::Overlap],
) -> anyhow::Result<()> {
    let path = Path::new(dir).join("ovlp.paf");
    let path = path.to_str().unwrap();
    let mut out = pgr::libs::io::writer(path).with_context(|| format!("failed to open {path}"))?;
    for ov in overlaps {
        let rec = super::common::to_paf(ov, unitigs);
        pgr::libs::paf::record::write_paf_record(&mut out, &rec)?;
    }
    out.flush()?;
    Ok(())
}

/// Writes the layout TSV.
fn dump_layouts(
    dir: &str,
    unitigs: &[Unitig],
    layouts: &[crate::libs::olc::layout::Layout],
) -> anyhow::Result<()> {
    let path = Path::new(dir).join("layout.tsv");
    let path = path.to_str().unwrap();
    let mut out = pgr::libs::io::writer(path).with_context(|| format!("failed to open {path}"))?;
    super::common::write_layout_tsv(&mut out, unitigs, layouts)?;
    out.flush()?;
    Ok(())
}

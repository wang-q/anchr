use crate::libs::map::{map_files, read_fasta, MapOptions};
use crate::libs::olc::anchor::{
    anchor_regions, anchor_stats, coverage_from_alignments, extract_anchors, Alignment,
    AnchorOptions,
};
use anyhow::Context;
use clap::{value_parser, Arg, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for anchor.
pub fn make_subcommand() -> Command {
    Command::new("anchor")
        .about("Selects reliable anchors (properly covered regions) from unitigs")
        .after_help(
            r###"
Selects the reliable anchors of the unitigs by mapping the reads back with
perfect matches (`anchr asm map` semantics), computing the per-base depth,
and keeping positions inside `[lower, upper]`:
`lower = max(mincov, (median - mscale*MAD)/lscale)`,
`upper = (median + mscale*MAD)*uscale`.
Low-coverage positions are likely errors, high-coverage positions likely
repeats; both are excluded, and the remaining contiguous stretches are the
anchors (legacy `anchr anchors` flow, modern `asm` implementation).

Anchors are the reliable fragments to feed the OLC merge: assemble per
coverage set with `anchr asm multik`, select anchors per set, then merge
them all with `anchr asm olc --unitigs`.

Notes:
* The first input is the unitig/contig FASTA; the remaining inputs are the
  read files (FASTA/FASTQ, plain or gzipped) — use the SAME coverage subset
  that produced the unitigs (not the full read set), so the depth matches
* Output FASTA names carry the source unitig and interval
  (`<unitig>_<start>-<end>`), 70-column wrapped

Examples:
1. Anchors from one coverage set:
   anchr asm anchor k40.fa reads40.fq.gz -o anchors40.fa
2. Tune the coverage window:
   anchr asm anchor k40.fa reads40.fq.gz -o anchors40.fa \
       --mincov 10 --mscale 3 --lscale 3 --uscale 2 --min-anchor-len 1000
"###,
        )
        .arg(crate::cmd::args::infiles_arg_with_numargs(
            "Unitig/contig FASTA, then read file(s); use --list-files for a one-path-per-line list",
            2..,
        ))
        .arg(crate::cmd::args::outfile_arg())
        .arg(
            Arg::new("mincov")
                .long("mincov")
                .num_args(1)
                .default_value("5")
                .value_parser(value_parser!(u32))
                .help("Absolute floor for the lower coverage bound"),
        )
        .arg(
            Arg::new("mscale")
                .long("mscale")
                .num_args(1)
                .default_value("3")
                .value_parser(value_parser!(f64))
                .help("Median absolute deviation multiplier"),
        )
        .arg(
            Arg::new("lscale")
                .long("lscale")
                .num_args(1)
                .default_value("3")
                .value_parser(value_parser!(f64))
                .help("Lower-window divider"),
        )
        .arg(
            Arg::new("uscale")
                .long("uscale")
                .num_args(1)
                .default_value("2")
                .value_parser(value_parser!(f64))
                .help("Upper-window multiplier"),
        )
        .arg(
            Arg::new("min_anchor_len")
                .long("min-anchor-len")
                .num_args(1)
                .default_value("500")
                .value_parser(value_parser!(usize))
                .help("Minimum anchor length in bases"),
        )
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("31")
                .value_parser(value_parser!(usize))
                .help("Seed k-mer length for the perfect-match read mapping"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("8")
                .value_parser(value_parser!(usize))
                .help("Worker threads for the read mapping"),
        )
        .arg(
            Arg::new("stats")
                .long("stats")
                .num_args(1)
                .help("Output TSV with Mapped/median/MAD/lower/upper/SumOthers"),
        )
        .arg(
            Arg::new("list_files")
                .long("list-files")
                .action(clap::ArgAction::SetTrue)
                .help("Treat infiles as list files, one sequence file path per line"),
        )
}

/// Execute the anchor command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let is_list = args.get_flag("list_files");
    let mut infiles: Vec<String> = Vec::new();
    for f in args.get_many::<String>("infiles").unwrap() {
        infiles.extend(pgr::libs::par::resolve_paths(f, is_list)?);
    }
    anyhow::ensure!(
        infiles.len() >= 2,
        "need a unitig FASTA and at least one read file"
    );
    let unitig_file = infiles[0].clone();
    let read_files = infiles[1..].to_vec();
    let outfile = crate::cmd::args::get_outfile(args);
    crate::cmd::args::ensure_outfile_distinct(outfile, infiles.iter().map(|s| s.as_str()))?;

    let opts = AnchorOptions {
        mincov: *args.get_one::<u32>("mincov").unwrap(),
        mscale: *args.get_one::<f64>("mscale").unwrap(),
        lscale: *args.get_one::<f64>("lscale").unwrap(),
        uscale: *args.get_one::<f64>("uscale").unwrap(),
        min_len: *args.get_one::<usize>("min_anchor_len").unwrap(),
    };
    let k = *args.get_one::<usize>("kmer").unwrap();
    let parallel = *args.get_one::<usize>("parallel").unwrap();

    // Map the reads back to the unitigs (perfect matches only).
    let refs = read_fasta(&[unitig_file])?;
    let tempdir = tempfile::Builder::new()
        .prefix("anchr_asm_anchor_")
        .tempdir()
        .context("failed to create tempdir")?;
    let sam_path = tempdir.path().join("mapped.sam");
    let sam_str = sam_path.to_str().unwrap();
    let mapped = map_files(
        &refs,
        &read_files,
        &MapOptions {
            k,
            outm: Some(sam_str.to_string()),
            outu: None,
            paired: false,
            max_reads: None,
            parallel,
        },
    )?;
    let mapped_ratio = if mapped.reads_in > 0 {
        mapped.mapped as f64 / mapped.reads_in as f64
    } else {
        0.0
    };

    // Parse the mapped SAM into per-unitig alignments (perfect M CIGAR).
    let name_to_idx: std::collections::HashMap<String, usize> = refs
        .iter()
        .enumerate()
        .map(|(i, r)| (r.name.split(',').next().unwrap_or(&r.name).to_string(), i))
        .collect();
    let lens: Vec<usize> = refs.iter().map(|r| r.seq.len()).collect();
    let mut aligns: Vec<Alignment> = Vec::new();
    let mut reader =
        pgr::libs::io::reader(sam_str).with_context(|| format!("failed to open SAM {sam_str}"))?;
    let mut line = String::new();
    let mut line_no = 0usize;
    while reader.read_line(&mut line)? > 0 {
        line_no += 1;
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        if fields.len() < 6 || fields[0].starts_with('@') {
            line.clear();
            continue;
        }
        let Some(&ri) = name_to_idx.get(fields[2].split(',').next().unwrap_or("")) else {
            anyhow::bail!("SAM line {line_no}: unknown reference {}", fields[2]);
        };
        let pos: usize = fields[3]
            .parse()
            .with_context(|| format!("SAM line {line_no}: bad POS {}", fields[3]))?;
        let cigar = fields[5];
        // Perfect matches only: `<len>M`.
        let mlen = cigar
            .strip_suffix('M')
            .and_then(|s| s.parse::<usize>().ok())
            .with_context(|| format!("SAM line {line_no}: unexpected CIGAR {cigar}"))?;
        aligns.push((ri, pos, pos + mlen - 1));
        line.clear();
    }

    // Coverage window and anchor regions.
    let covs = coverage_from_alignments(&lens, &aligns);
    let stats = anchor_stats(&covs, &opts);
    let regions = anchor_regions(&covs, &opts, stats.lower, stats.upper);
    let seqs: Vec<Vec<u8>> = refs.iter().map(|r| r.seq.clone()).collect();
    let anchors = extract_anchors(&seqs, &regions);

    if let Some(stats_file) = args.get_one::<String>("stats") {
        let total_bases: usize = lens.iter().sum();
        let anchor_bases: usize = regions.iter().map(|&(_, a, b)| b - a + 1).sum();
        let mut out = pgr::libs::io::writer(stats_file)
            .with_context(|| format!("failed to open stats file {stats_file}"))?;
        writeln!(
            out,
            "{:.3}\t{}\t{}\t{:.1}\t{:.1}\t{}",
            mapped_ratio,
            stats.median,
            stats.mad,
            stats.lower,
            stats.upper,
            total_bases - anchor_bases
        )?;
    }

    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    for (i, (seq, &(ri, a, b))) in anchors.iter().zip(regions.iter()).enumerate() {
        let name = &refs[ri].name.split(',').next().unwrap_or(&refs[ri].name);
        writeln!(out, ">anchor_{}_{}_{}-{}", i + 1, name, a, b)?;
        for chunk in seq.chunks(70) {
            out.write_all(chunk)?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(())
}

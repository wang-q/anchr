use crate::libs::asm::assemble::read_records;
use crate::libs::asm::extend::{extend_contigs, ExtendOptions};
use crate::libs::map::read_fasta;
use anyhow::{ensure, Context};
use clap::{value_parser, Arg, ArgMatches, Command};
use std::io::Write;

/// Build the clap subcommand for extend.
pub fn make_subcommand() -> Command {
    Command::new("extend")
        .about("Extends contig ends along read-supported k-mer paths")
        .after_help(
            r###"
Walks each contig end base-by-base through the reads' k-mer graph: a base is
appended only when it has a strict majority of read support (>= --min-support
reads and >= 2x the runner-up), so junctions and repetitive contexts stop the
extension instead of joining distant loci. This closes small coverage gaps at
contig ends (megahit local-assembly goal) without reassembling the reads; the
output keeps the input contig order and leaves unsupported ends unchanged.

Notes:
* Input is one contigs FASTA plus one or more read FASTA/FASTQ files (plain
  or gzipped); pairing is irrelevant for the extension walk
* The extension is capped at --max-extend bases per end and discarded when
  both ends together extend less than --min-extend bases
* Output sequences are wrapped at 70 columns

Examples:
1. Extend unitigs with the reads they were assembled from:
   anchr asm extend unitigs.fasta pe.cor.fa.gz -o unitigs.ext.fasta
2. Raise the read-support bar:
   anchr asm extend unitigs.fasta reads.fq.gz --min-support 3 -o out.fasta
"###,
        )
        .arg(
            Arg::new("contigs")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Input contigs FASTA file"),
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(2)
                .help("Read file(s): FASTA/FASTQ, plain or gzipped"),
        )
        .arg(crate::cmd::args::outfile_arg())
        .arg(
            Arg::new("kmer")
                .long("kmer")
                .short('k')
                .num_args(1)
                .default_value("31")
                .value_parser(value_parser!(usize))
                .help("Seed k-mer length (default 31)"),
        )
        .arg(
            Arg::new("max_extend")
                .long("max-extend")
                .num_args(1)
                .default_value("500")
                .value_parser(value_parser!(usize))
                .help("Maximum extension in bases per contig end (default 500)"),
        )
        .arg(
            Arg::new("min_support")
                .long("min-support")
                .num_args(1)
                .default_value("2")
                .value_parser(value_parser!(u32))
                .help("Minimum read support for each appended base (default 2)"),
        )
        .arg(
            Arg::new("min_extend")
                .long("min-extend")
                .num_args(1)
                .default_value("0")
                .value_parser(value_parser!(usize))
                .help("Minimum total extension to keep (default 0)"),
        )
}

/// Execute the extend command.
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let contigs_file = args.get_one::<String>("contigs").unwrap();
    let infiles: Vec<String> = args
        .get_many::<String>("infiles")
        .unwrap()
        .cloned()
        .collect();
    let contigs: Vec<(String, Vec<u8>)> = read_fasta(std::slice::from_ref(contigs_file))?
        .into_iter()
        .map(|r| (r.name, r.seq))
        .collect();
    let reads = read_records(&infiles)?;
    ensure!(!reads.is_empty(), "no reads found in the input file(s)");
    let opts = ExtendOptions {
        k: *args.get_one::<usize>("kmer").unwrap(),
        max_extend: *args.get_one::<usize>("max_extend").unwrap(),
        min_support: *args.get_one::<u32>("min_support").unwrap(),
        min_extend: *args.get_one::<usize>("min_extend").unwrap(),
    };
    let extended = extend_contigs(&contigs, reads, &opts)
        .with_context(|| format!("failed to extend contigs in {contigs_file}"))?;
    let outfile = crate::cmd::args::get_outfile(args);
    crate::cmd::args::ensure_outfile_distinct(outfile, std::iter::once(contigs_file.as_str()))?;
    let mut out = pgr::libs::io::writer(outfile)
        .with_context(|| format!("failed to open output {outfile}"))?;
    for (i, (name, seq)) in extended.iter().enumerate() {
        writeln!(
            out,
            ">{}",
            if name.is_empty() {
                format!("extended_{i}")
            } else {
                name.clone()
            }
        )?;
        for chunk in seq.chunks(70) {
            out.write_all(chunk)?;
            out.write_all(b"\n")?;
        }
    }
    out.flush()?;
    Ok(())
}

//! Shared clap argument builders for subcommands.

use clap::{builder, Arg, ArgMatches};


/// Standard `-o/--outfile` argument defaulting to stdout.
pub fn outfile_arg() -> Arg {
    Arg::new("outfile")
        .long("outfile")
        .short('o')
        .num_args(1)
        .default_value("stdout")
        .help("Output filename. [stdout] for screen")
}

/// `-o/--outfile` with a custom default value.
/// `-o/--outfile` required (no default).
pub fn outfile_arg_required() -> Arg {
    Arg::new("outfile")
        .long("outfile")
        .short('o')
        .num_args(1)
        .required(true)
        .help("Output filename")
}

/// Standard `-o/--outdir` argument defaulting to stdout.
/// Required positional `infile` argument with a custom help text.
/// Index is auto-assigned by clap — do not add `.index(N)` to other positionals
/// unless this is the only positional or all positionals use explicit indices.
pub fn infile_arg_required_with_help(help: &'static str) -> Arg {
    Arg::new("infile").required(true).num_args(1).help(help)
}

/// Standard positional `infiles` argument (one or more, required) at index 1.
///
/// `label` is the format name used in the help text (e.g. `"FASTA"`,
/// `"block FA"`, `"2bit"`). Use inline definition with a different `.index(N)`
/// when another positional precedes `infiles`.
pub fn infiles_arg(label: &str) -> Arg {
    Arg::new("infiles")
        .required(true)
        .num_args(1..)
        .index(1)
        .help(format!("Input {label} file(s) to process"))
}

/// Positional `infiles` argument at a custom index (required, 1 or more files).
/// Use when another positional precedes `infiles`.
/// Positional `infiles` argument at index 1 with custom num_args and help.
/// Use when the default `1..` range doesn't fit (e.g., `2..`, `1..=2`, `1..=4`).
pub fn infiles_arg_with_numargs(
    help: &'static str,
    num_args: impl clap::builder::IntoResettable<builder::ValueRange>,
) -> Arg {
    Arg::new("infiles")
        .required(true)
        .num_args(num_args)
        .index(1)
        .help(help)
}

/// Positional `target` genome file argument (required, index 1).
/// Standard `-r/--rgfile` argument (file of regions, one per line).
pub fn rgfile_arg() -> Arg {
    Arg::new("rgfile")
        .long("rgfile")
        .short('r')
        .num_args(1)
        .help("File of regions, one per line")
}

/// Standard `-t/--t-sizes` argument (target chromosome sizes file).
/// `-p/--parallel` with a custom default value.
pub fn parallel_arg_with_default(default: &'static str) -> Arg {
    Arg::new("parallel")
        .long("parallel")
        .short('p')
        .num_args(1)
        .default_value(default)
        .value_parser(clap::builder::RangedU64ValueParser::<usize>::new().range(1..=1024))
        .help("Number of threads for parallel processing (1..=1024)")
}

/// Parse a `--parallel` value that may be `auto` (logical CPU count) or an
/// integer in `1..=1024`. Returns a friendly error for invalid or out-of-range
/// values before any thread pool is created.
pub fn parse_parallel_auto(s: &str) -> anyhow::Result<usize> {
    if s == "auto" {
        return Ok(pgr::libs::sys::logical_cpus());
    }
    let n: usize = s
        .parse()
        .map_err(|_| anyhow::anyhow!("invalid --parallel: {}", s))?;
    if !(1..=1024).contains(&n) {
        anyhow::bail!("--parallel must be in 1..=1024, got {}", n);
    }
    Ok(n)
}

/// `--no-ns` flag (output size without Ns).
/// `--name-prefix` argument with an optional default value.
pub fn name_prefix_arg(default: Option<&'static str>) -> Arg {
    let arg = Arg::new("name_prefix").long("name-prefix").num_args(1);
    match default {
        Some(d) => arg.default_value(d).help("Prefix of record names"),
        None => arg.help("Add prefix to sequence names"),
    }
}

/// Extract the `outfile` value from `args` as `&str`.
pub fn get_outfile(args: &ArgMatches) -> &str {
    args.get_one::<String>("outfile").unwrap()
}

/// Reject an `-o` path that would overwrite an input file before it has been
/// read (streaming commands open the output before consuming their inputs).
pub fn ensure_outfile_distinct<'a>(
    outfile: &str,
    inputs: impl IntoIterator<Item = &'a str>,
) -> anyhow::Result<()> {
    // `stdout` is the screen sentinel (`writer` never creates a file for
    // it), so it cannot overwrite an input even when a file literally named
    // `stdout` exists in the working directory.
    if outfile == "stdout" {
        return Ok(());
    }
    for input in inputs {
        // `stdin` is the stream sentinel (`reader` never opens a file for
        // it), so no input file exists that an output could overwrite.
        if input == "stdin" {
            continue;
        }
        if pgr::libs::io::same_path(outfile, input) {
            anyhow::bail!("output file {} is also an input file", outfile);
        }
    }
    Ok(())
}

/// Extract the `infile` value from `args` as `&str`.
/// Collect region strings from `ranges` (positional, optional) and `rgfile`
/// (`-r/--rgfile`) arguments. Returns the combined list.
pub fn collect_ranges(args: &ArgMatches) -> anyhow::Result<Vec<String>> {
    let mut ranges: Vec<String> = if args.contains_id("ranges") {
        args.get_many::<String>("ranges")
            .unwrap()
            .cloned()
            .collect()
    } else {
        Vec::new()
    };
    if args.contains_id("rgfile") {
        let mut rgs =
            pgr::libs::io::read_names::<Vec<String>>(args.get_one::<String>("rgfile").unwrap())?;
        ranges.append(&mut rgs);
    }
    Ok(ranges)
}

/// Add POA scoring arguments (`--match`, `--mismatch`, `--gap-open`,
/// `--gap-extend`) to `cmd`. When `with_shorts` is true, also registers the
/// `-m`/`-n`/`-g`/`-e` short flags (used by `fas consensus`; paf commands
/// pass false because `-m` collides with `--max-depth`).
/// `-k/--kmer` size argument with a custom default value.
pub fn kmer_arg_with_default(default: &'static str) -> Arg {
    Arg::new("kmer")
        .long("kmer")
        .short('k')
        .num_args(1)
        .default_value(default)
        .value_parser(clap::value_parser!(usize))
        .help("K-mer size")
}

/// `-w/--window` size argument (default: 1, for minimizers).
/// Positional `ranges` argument (optional, index 2).
pub fn ranges_arg() -> Arg {
    Arg::new("ranges")
        .required(false)
        .index(2)
        .num_args(0..)
        .help("Ranges of interest")
}

/// `--replace-tsv` argument (required) for replace commands.
/// Optional `-q/--qual-thresh` argument for quality-weighted commands
/// (migrated from `pgr cmd_pgr::kmer::qhist`).
pub fn qual_thresh_arg() -> Arg {
    Arg::new("qual_thresh")
        .long("qual-thresh")
        .short('q')
        .num_args(1)
        .value_parser(clap::value_parser!(u8))
        .help("Quality ASCII threshold (default: detected Phred offset + 5)")
}

/// Optional `-b/--bits` argument for quality-weighted commands
/// (migrated from `pgr cmd_pgr::kmer::qhist`).
pub fn bits_arg() -> Arg {
    Arg::new("bits")
        .long("bits")
        .short('b')
        .num_args(1)
        .default_value("7")
        .value_parser(clap::value_parser!(u8))
        .help("Count bits (quorum create_database -b; max count = 2^bits - 1)")
}

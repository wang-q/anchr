use clap::*;
use indexmap::IndexSet;
use pgr::libs::ds::IntSpan;
use pgr::libs::runlist::{depth_at_least, depth_by_level};
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io::{BufRead, Write};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("covered")
        .about("Covered regions from .ovlp.tsv files")
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input files to use"),
        )
        .arg(
            Arg::new("coverage")
                .long("coverage")
                .short('c')
                .num_args(1)
                .default_value("3")
                .value_parser(value_parser!(i32))
                .help("minimal coverage"),
        )
        .arg(
            Arg::new("len")
                .long("len")
                .short('l')
                .num_args(1)
                .default_value("1000")
                .value_parser(value_parser!(i32))
                .help("minimal length of overlaps"),
        )
        .arg(
            Arg::new("idt")
                .long("idt")
                .short('i')
                .num_args(1)
                .default_value("0.0")
                .value_parser(value_parser!(f32))
                .help("minimal identities of overlaps"),
        )
        .arg(
            Arg::new("paf")
                .long("paf")
                .action(ArgAction::SetTrue)
                .help("PAF as input format"),
        )
        .arg(
            Arg::new("longest")
                .long("longest")
                .action(ArgAction::SetTrue)
                .help("only keep the longest span"),
        )
        .arg(
            Arg::new("base")
                .long("base")
                .action(ArgAction::SetTrue)
                .help("per base coverage"),
        )
        .arg(
            Arg::new("mean")
                .long("mean")
                .action(ArgAction::SetTrue)
                .help("mean coverage"),
        )
        .arg(
            Arg::new("outfile")
                .long("outfile")
                .short('o')
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    //----------------------------
    // Loading
    //----------------------------
    let mut writer = pgr::libs::io::writer(args.get_one::<String>("outfile").unwrap())?;

    let coverage = *args.get_one::<i32>("coverage").unwrap();
    let min_len = *args.get_one::<i32>("len").unwrap();
    let min_idt = *args.get_one::<f32>("idt").unwrap();

    let is_paf = args.get_flag("paf");
    let is_longest = args.get_flag("longest");
    let is_base = args.get_flag("base");
    let is_mean = args.get_flag("mean");

    // seq_name => half-open intervals; seq_name => sequence length
    let mut ivs_of: HashMap<String, Vec<(u32, u32)>> = HashMap::new();
    let mut len_of: HashMap<String, u32> = HashMap::new();
    let mut index_of: IndexSet<String> = IndexSet::new();
    let mut seen: HashSet<(usize, usize)> = HashSet::new();

    for infile in args.get_many::<String>("infiles").unwrap() {
        let reader = pgr::libs::io::reader(infile)?;
        for line in reader.lines().map_while(Result::ok) {
            let ovlp = if is_paf {
                anchr::Overlap::from_paf(&line)
            } else {
                anchr::Overlap::new(&line)
            };
            let f_id = ovlp.f_id.to_string();
            let g_id = ovlp.g_id.to_string();

            // ignore self overlapping
            if f_id == g_id {
                continue;
            }

            // ignore poor overlaps
            if ovlp.len < min_len {
                continue;
            }
            if ovlp.idt < min_idt {
                continue;
            }

            // skip duplicated overlaps, i.e., f -> g and g -> f
            let (f_idx, _) = index_of.insert_full(f_id.clone());
            let (g_idx, _) = index_of.insert_full(g_id.clone());
            let tup = (f_idx.min(g_idx), f_idx.max(g_idx));
            // If the set did not have this value present, true is returned.
            let not_seen = seen.insert(tup);
            if !not_seen {
                continue;
            }

            // collect 1-based inclusive overlap intervals; pgr runlist takes
            // half-open [s, e) so the +1 conversion happens in the helpers
            ivs_of
                .entry(f_id.clone())
                .or_default()
                .push((ovlp.f_begin as u32, ovlp.f_end as u32));
            len_of.entry(f_id.clone()).or_insert(ovlp.f_len as u32);

            ivs_of
                .entry(g_id.clone())
                .or_default()
                .push((ovlp.g_begin as u32, ovlp.g_end as u32));
            len_of.entry(g_id.clone()).or_insert(ovlp.g_len as u32);
        }
    }

    //----------------------------
    // Output
    //----------------------------
    let mut keys = ivs_of.keys().map(|k| k.to_string()).collect::<Vec<String>>();
    keys.sort();

    for key in &keys {
        let mut _out_line = String::new();
        let ivs = &ivs_of[key];
        let len = *len_of.get(key).unwrap();

        if is_base {
            let tiers = depth_tiers(ivs, len, coverage);
            _out_line = base_lines(key, &tiers);
        } else if is_mean {
            let tiers = depth_tiers(ivs, len, coverage);
            _out_line = mean_line(key, &tiers);
        } else {
            let intspan = depth_at_least(ivs, coverage as u32);

            if !is_longest || intspan.span_size() <= 1 {
                _out_line = format!("{}:{}", key, intspan);
            } else {
                _out_line = longest_line(key, &intspan);
            }
        }

        if !_out_line.is_empty() {
            writer.write_all((_out_line + "\n").as_ref())?;
        }
    }

    Ok(())
}

/// Per-depth tiers over 1-based inclusive overlap intervals, using the pgr
/// sweep-line runlist. Depth is clamped to `max` and the `-1` (full-length)
/// / `0` (uncovered) tiers are added, matching the former vendored
/// `Coverage` semantics used by `--base` / `--mean`.
fn depth_tiers(ivs: &[(u32, u32)], len: u32, max: i32) -> BTreeMap<i32, IntSpan> {
    let half_open: Vec<(u32, u32)> = ivs.iter().map(|&(s, e)| (s, e + 1)).collect();
    let by_level = depth_by_level(&half_open, 1);

    let mut out: BTreeMap<i32, IntSpan> = BTreeMap::new();
    let mut covered = IntSpan::new();
    for (depth, is) in &by_level {
        let depth: i32 = depth.parse().unwrap();
        out.entry(depth.min(max)).or_default().merge(is);
        covered.merge(is);
    }

    let mut zero = IntSpan::from_pair(1, len as i32);
    zero.subtract(&covered);
    out.insert(-1, IntSpan::from_pair(1, len as i32));
    out.insert(0, zero);
    // the vendored `Coverage` pre-filled every tier 1..=max; consumers index
    // all of them (e.g. `mean_line` sums over `0..=max`)
    for i in 1..=max {
        out.entry(i).or_default();
    }
    out
}

fn base_lines(key: &str, tiers: &BTreeMap<i32, IntSpan>) -> String {
    let mut basecovs: HashMap<i32, i32> = HashMap::new();
    let max_tier = tiers.keys().max().unwrap();
    for i in 0..=*max_tier {
        for pos in tiers[&i].elements() {
            basecovs.insert(pos, i);
        }
    }

    let mut sorted: Vec<i32> = basecovs.keys().copied().collect();
    sorted.sort_unstable();

    let mut out_lines: Vec<String> = vec![];
    for pos in sorted {
        let line = format!("{}\t{}\t{}", key, pos - 1, basecovs[&pos]);
        out_lines.push(line);
    }

    out_lines.join("\n")
}

fn mean_line(key: &str, tiers: &BTreeMap<i32, IntSpan>) -> String {
    let total_len = tiers[&-1].cardinality();
    let max_tier = tiers.keys().max().unwrap();
    let mut sum = 0;
    for i in 0..=*max_tier {
        sum += i * tiers[&i].cardinality();
    }
    let mean_cov = sum as f32 / total_len as f32;

    format!("{}\t{}\t{:.1}", key, total_len, mean_cov)
}

fn longest_line(key: &str, intspan: &IntSpan) -> String {
    let ranges = intspan.ranges();

    let mut sizes: Vec<i32> = Vec::new();
    for i in 0..intspan.span_size() {
        let size = ranges[i * 2 + 1] - ranges[i * 2] + 1;
        sizes.push(size);
    }

    let mut max_i = 0;
    for i in 0..intspan.span_size() {
        let size = sizes[i];
        if size > sizes[max_i] {
            max_i = i;
        }
    }

    let mut longest = IntSpan::new();
    longest.add_pair(ranges[max_i * 2], ranges[max_i * 2 + 1]);

    format!("{}:{}", key, longest)
}

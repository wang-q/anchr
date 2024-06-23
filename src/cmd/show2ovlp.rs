use clap::*;
use regex::Regex;
use std::collections::BTreeMap;
use std::io::BufRead;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("show2ovlp")
        .about("Convert LAshow outputs to overlaps")
        .after_help(
            r#"
f_id and g_id are integers, --orig convert them to the original ones

"#,
        )
        .arg(
            Arg::new("show.txt")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("replace.tsv")
                .required(true)
                .num_args(1)
                .index(2)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("orig")
                .long("orig")
                .action(ArgAction::SetTrue)
                .help("Original names of sequences"),
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
    let is_orig = args.get_flag("orig");
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    lazy_static::lazy_static! {
        static ref RE_LEN: Regex = Regex::new(
            r"(?xi)
            \/(?<id>\d+)\/\d+_(?<len>\d+)
            ",
        )
        .unwrap();
    }

    let mut len_of = BTreeMap::new();
    let mut replace_of = BTreeMap::new();
    for line in &intspan::read_lines(args.get_one::<String>("replace.tsv").unwrap()) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 2 {
            continue;
        }

        let Some(caps) = RE_LEN.captures(fields[0]) else {
            continue;
        };

        let id = caps["id"].to_string();
        let len = caps["len"].parse::<i32>().unwrap();
        len_of.insert(id.to_string(), len);

        replace_of.insert(id.to_string(), fields[1].to_string());
    }
    // eprintln!("len_of = {:#?}", len_of);

    lazy_static::lazy_static! {
        static ref RE_SHOW: Regex = Regex::new(
            r"(?xi)
            ^\D*
            (?<f_id>\d+)
            \s+(?<g_id>\d+)
            \s+(?<g_orien>\w)
            \D+(?<f_B>\d+)
            \D+(?<f_E>\d+)
            \D+(?<g_B>\d+)
            \D+(?<g_E>\d+)
            \D+(?<idt>[\d.]+)
            .*$
            ",
        )
        .unwrap();
    }

    let reader = intspan::reader(args.get_one::<String>("show.txt").unwrap());
    for line in reader.lines().map_while(Result::ok) {
        let show = line.replace(',', "");
        let Some(caps) = RE_SHOW.captures(&show) else {
            continue;
        };

        let f_id = caps["f_id"].to_string();
        let g_id = caps["g_id"].to_string();

        if !len_of.contains_key(&f_id) {
            continue;
        }
        if !len_of.contains_key(&g_id) {
            continue;
        }

        let g_strand = if caps["g_orien"].eq("n") { 0 } else { 1 };
        let f_begin = caps["f_B"].parse::<i32>().unwrap();
        let f_end = caps["f_E"].parse::<i32>().unwrap();
        let g_begin = caps["g_B"].parse::<i32>().unwrap();
        let g_end = caps["g_E"].parse::<i32>().unwrap();
        let idt = caps["idt"].parse::<f32>().unwrap();
        let idt = (100.0 - idt) / 100.0;

        // relations
        let contained = if (len_of.get(&g_id).unwrap() < len_of.get(&f_id).unwrap())
            && (g_begin < 1)
            && (len_of.get(&g_id).unwrap() - g_end < 1)
        {
            "contains".to_string()
        } else if (len_of.get(&f_id).unwrap() < len_of.get(&g_id).unwrap())
            && (f_begin < 1)
            && (len_of.get(&f_id).unwrap() - f_end < 1)
        {
            "contained".to_string()
        } else {
            "overlap".to_string()
        };

        let ovlp = anchr::Overlap {
            f_id: if is_orig {
                replace_of.get(&f_id).unwrap().clone()
            } else {
                f_id.clone()
            },
            g_id: if is_orig {
                replace_of.get(&g_id).unwrap().clone()
            } else {
                g_id.clone()
            },
            len: f_end - f_begin,
            idt,
            f_strand: 0,
            f_begin,
            f_end,
            f_len: *len_of.get(&f_id).unwrap(),
            g_strand,
            g_begin,
            g_end,
            g_len: *len_of.get(&g_id).unwrap(),
            contained,
        };

        //----------------------------
        // Output
        //----------------------------
        writer.write_all((ovlp.to_string() + "\n").as_ref())?;
    }

    Ok(())
}

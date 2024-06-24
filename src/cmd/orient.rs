use clap::*;
use cmd_lib::*;
use petgraph::prelude::NodeIndex;
use petgraph::*;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, Write};
use std::path::Path;
use std::{env, fs};
use tempfile::{Builder, TempDir};

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("orient")
        .about("Orient overlapped sequences to the same strand")
        .after_help(
            r###"
* This command is for small files.
* All operations are running in a tempdir and no intermediate files are kept.

"###,
        )
        .arg(
            Arg::new("infiles")
                .required(true)
                .num_args(1..)
                .index(1)
                .help("Set the input files to use"),
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
                .default_value("0.85")
                .value_parser(value_parser!(f32))
                .help("minimal identities of overlaps"),
        )
        .arg(
            Arg::new("restrict")
                .long("restrict")
                .short('r')
                .num_args(1)
                .help("Restrict to known pairs"),
        )
        .arg(
            Arg::new("parallel")
                .long("parallel")
                .short('p')
                .num_args(1)
                .default_value("8")
                .value_parser(value_parser!(i32))
                .help("Number of threads"),
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
    // Args
    //----------------------------
    let min_len = *args.get_one::<i32>("len").unwrap();
    let min_idt = *args.get_one::<f32>("idt").unwrap();

    let parallel = *args.get_one::<i32>("parallel").unwrap();

    let outfile = args.get_one::<String>("outfile").unwrap();

    let curdir = env::current_dir()?;
    let anchr = env::current_exe().unwrap().display().to_string();
    let tempdir = Builder::new().prefix("anchr_orient_").tempdir()?;
    let tempdir_str = tempdir.path().to_str().unwrap();

    run_cmd!(info "==> Paths")?;
    run_cmd!(info "    \"anchr\"   = ${anchr}")?;
    run_cmd!(info "    \"curdir\"  = ${curdir}")?;
    run_cmd!(info "    \"tempdir\" = ${tempdir_str}")?;

    //----------------------------
    // Paths
    //----------------------------
    run_cmd!(info "==> Absolute paths")?;
    // basename => abs_path
    let mut abs_infiles = vec![];
    for infile in args.get_many::<String>("infiles").unwrap() {
        if infile == "stdin" {
            abs_infiles.push("stdin".to_string());
        } else {
            let absolute = intspan::absolute_path(infile)
                .unwrap()
                .display()
                .to_string();

            abs_infiles.push(absolute.to_string());
        }
    }

    let abs_outfile = if outfile == "stdout" {
        outfile.to_string()
    } else {
        intspan::absolute_path(outfile)
            .unwrap()
            .display()
            .to_string()
    };

    let mut abs_restrict = "".to_string();
    if args.contains_id("restrict") {
        let restrict = args.get_one::<String>("restrict").unwrap();
        abs_restrict = intspan::absolute_path(restrict)
            .unwrap()
            .display()
            .to_string();
    }

    //----------------------------
    // Operating
    //----------------------------
    run_cmd!(info "==> Switch to tempdir")?;
    env::set_current_dir(tempdir_str)?;

    run_cmd!(info "==> Sort sequences by lengths")?;
    let mut infiles = vec![];
    for (i, infile) in abs_infiles.iter().enumerate() {
        run_cmd!(
            faops size ${infile} |
                sort -n -r -k2,2 |
                cut -f 1 > infile.${i}.order.txt
        )?;
        run_cmd!(
            faops order ${infile} infile.${i}.order.txt infile.${i}.fasta
        )?;
        infiles.push(format!("infile.{}.fasta", i));
    }

    run_cmd!(info "==> Preprocess reads to format them for dazzler")?;
    run_cmd!(
        ${anchr} dazzname $[infiles] -o renamed.fasta
    )?;

    run_cmd!(info "==> `anchr overlap`")?;
    run_cmd!(
        ${anchr} overlap renamed.fasta --len ${min_len} --idt ${min_idt} --parallel ${parallel} -o renamed.ovlp.tsv
    )?;

    if !abs_restrict.is_empty() {
        run_cmd!(info "==> Filter overlaps")?;
        run_cmd!(
            rgr replace renamed.ovlp.tsv renamed.fasta.replace.tsv |
                ${anchr} restrict stdin ${abs_restrict} |
                rgr replace stdin renamed.fasta.replace.tsv -r -o restrict.ovlp.tsv
        )?;
    }

    run_cmd!(info "==> Build ovlp graph")?;
    // all nodes stored in one graph
    // Node weight - String - node names
    // Edge weight - i32 - g_strand
    let mut graph: Graph<String, i32, Undirected> = Graph::new_undirected();
    // cache node indices
    // petgraph use NodeIndex to store and identify nodes
    let mut idx_of_id: BTreeMap<String, NodeIndex> = BTreeMap::new();

    let lines = if !abs_restrict.is_empty() {
        intspan::read_lines("restrict.ovlp.tsv")
    } else {
        intspan::read_lines("renamed.ovlp.tsv")
    };
    for line in &lines {
        let ovlp = anchr::Overlap::new(line);
        if ovlp.is_empty() {
            continue;
        }

        // ignore self overlapping
        if ovlp.f_id() == ovlp.g_id() {
            continue;
        }

        // ignore poor overlaps
        if *ovlp.len() < min_len {
            continue;
        }
        if *ovlp.idt() < min_idt {
            continue;
        }

        let f_idx = graph.add_node(ovlp.f_id().to_string());
        idx_of_id.insert(ovlp.f_id().to_string(), f_idx);

        let g_idx = graph.add_node(ovlp.g_id().to_string());
        idx_of_id.insert(ovlp.g_id().to_string(), g_idx);

        if graph.find_edge(f_idx, g_idx).is_none() {
            graph.add_edge(f_idx, g_idx, *ovlp.g_strand());
        }
    }
    // eprintln!("graph = {:#?}", graph);
    run_cmd!(info "==> To positive strands in each SCC")?;
    // no edge weights
    let mut new_graph: Graph<String, (), Undirected> = Graph::new_undirected();
    let scc: Vec<Vec<NodeIndex>> = petgraph::algo::tarjan_scc(&graph);
    for cc_indices in &scc {
        let count = cc_indices.len();

        if count < 2 {
            continue;
        }

        // set first sequence to positive strand


        // assign strands to other nodes
        let mut assigned: BTreeSet<usize> = BTreeSet::new();
        assigned.insert(0);
        let mut unhandled: BTreeSet<usize> = (1..count).collect();

        let mut edges = vec![];
        while assigned.len() < count {
            for i in assigned.iter().cloned().collect::<Vec<usize>>() {
                for j in unhandled.iter().cloned().collect::<Vec<usize>>() {
                    // not connected in old graph
                    let edge = graph.find_edge(idx_of_part[&parts[i]], idx_of_part[&parts[j]]);
                    if edge.is_none() {
                        continue;
                    }
                }
            }
        } // end of assigning
    }


    pause();

    //----------------------------
    // Done
    //----------------------------
    env::set_current_dir(&curdir)?;

    Ok(())
}

use std::io::Read;

fn pause() {
    let mut stdin = std::io::stdin();
    let mut stdout = std::io::stdout();

    // We want the cursor to stay at the end of the line, so we print without a newline and flush manually.
    write!(stdout, "Press any key to continue...").unwrap();
    stdout.flush().unwrap();

    // Read a single byte and discard
    let _ = stdin.read(&mut [0u8]).unwrap();
}

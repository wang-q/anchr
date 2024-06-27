use clap::*;
use cmd_lib::*;
use petgraph::prelude::*;
use petgraph::*;
use std::collections::{BTreeMap, BTreeSet};
use std::env;
use bio::io::fasta;
use tempfile::Builder;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
    Command::new("merge")
        .about("Merge overlapped unitigs")
        .after_help(
            r###"
* All operations are running in a tempdir and no intermediate files are kept.

"###,
        )
        .arg(
            Arg::new("infile")
                .required(true)
                .num_args(1)
                .index(1)
                .help("Set the input file to use"),
        )
        .arg(
            Arg::new("len")
                .long("len")
                .short('l')
                .num_args(1)
                .default_value("500")
                .value_parser(value_parser!(i32))
                .help("Minimal length of overlaps"),
        )
        .arg(
            Arg::new("idt")
                .long("idt")
                .short('i')
                .num_args(1)
                .default_value("0.98")
                .value_parser(value_parser!(f32))
                .help("Minimal identities of overlaps"),
        )
        .arg(
            Arg::new("svg")
                .long("svg")
                .action(ArgAction::SetTrue)
                .help("Write a .svg file representing the merge graph"),
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

    let is_svg = args.get_flag("svg");
    let parallel = *args.get_one::<i32>("parallel").unwrap();

    let infile = args.get_one::<String>("infile").unwrap();
    let outfile = args.get_one::<String>("outfile").unwrap();

    let curdir = env::current_dir()?;
    let anchr = env::current_exe().unwrap().display().to_string();
    let tempdir = Builder::new().prefix("anchr_contained_").tempdir()?;
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
    let abs_infile = if infile == "stdin" {
        infile.to_string()
    } else {
        intspan::absolute_path(infile)
            .unwrap()
            .display()
            .to_string()
    };

    let abs_outfile = if outfile == "stdout" {
        outfile.to_string()
    } else {
        intspan::absolute_path(outfile)
            .unwrap()
            .display()
            .to_string()
    };

    //----------------------------
    // Operating
    //----------------------------
    run_cmd!(info "==> Switch to tempdir")?;
    env::set_current_dir(tempdir_str)?;

    run_cmd!(info "==> `anchr overlap`")?;
    let file = abs_infile.clone();
    run_cmd!(
        ${anchr} overlap ${file} --len ${min_len} --idt ${min_idt} --parallel ${parallel} -o merge.ovlp.tsv
    )?;

    run_cmd!(info "==> Build ovlp graph")?;
    let mut seen = BTreeSet::new();
    let mut ovlps = vec![]; // In order for these overlaps to live long enough
    for line in &intspan::read_lines("merge.ovlp.tsv") {
        let ovlp = anchr::Overlap::new(line);
        if ovlp.is_empty() {
            continue;
        }

        // ignore self overlapping
        if ovlp.f_id == ovlp.g_id {
            continue;
        }

        // ignore poor overlaps
        if ovlp.len < min_len {
            continue;
        }
        if ovlp.idt < min_idt {
            continue;
        }

        // we've orient overlapped sequences to the same strand
        if ovlp.g_strand == 1 {
            continue;
        }

        // skip duplicated overlaps, i.e., f -> g and g -> f
        let tup = (
            ovlp.f_id.to_string().min(ovlp.g_id.to_string()),
            ovlp.f_id.to_string().max(ovlp.g_id.to_string()),
        );
        // If the set did not have this value present, true is returned.
        let not_seen = seen.insert(tup);
        if !not_seen {
            continue;
        }

        ovlps.push(ovlp);
    }

    // Node weight - &str - node names
    // Edge weight - i32 - append length, back index of strings
    let mut graph = DiGraphMap::new();
    for ovlp in ovlps.iter() {
        // contained unitigs have been removed
        if ovlp.f_begin > 0 && ovlp.f_end == ovlp.f_len {
            //          f.B        f.E
            // f ========+---------->
            // g         -----------+=======>
            //          g.B        g.E
            graph.add_edge(
                ovlp.f_id.as_str(),
                ovlp.g_id.as_str(),
                ovlp.g_len - ovlp.g_end,
            );
        } else if ovlp.g_begin > 0 && ovlp.g_end == ovlp.g_len {
            //          f.B        f.E
            // f         -----------+=======>
            // g ========+---------->
            //          g.B        g.E
            graph.add_edge(
                ovlp.g_id.as_str(),
                ovlp.f_id.as_str(),
                ovlp.f_len - ovlp.f_end,
            );
        }
    }

    // Remove cyclic nodes
    while let Some(cy) = algo::toposort(&graph, None).err() {
        let v = cy.node_id();
        graph.remove_node(v);
    }
    // Branching nodes will stay

    intspan::write_lines("overlapped.list", &graph.nodes().collect::<Vec<_>>())?;
    if is_svg {
        run_cmd!(info "==> Write .svg file")?;
        g2gv(&graph, &abs_outfile)?;
    }

    run_cmd!(info "==> Write non-overlapped sequences")?;
    let file = abs_infile.clone();
    run_cmd!(
        faops some -i -l 0 ${file} overlapped.list non-overlapped.fasta
    )?;

    run_cmd!(info "==> Merge")?;
    let seq_of = anchr::read_fasta(&abs_infile);
    let topo_sorted = algo::toposort(&graph, None).unwrap();
    let mut merge_of = BTreeMap::new();
    let mut serial = 1;
    let wccs = rustworkx_core::connectivity::connected_components(&graph);
    for wcc in wccs.iter() {
        let pieces: Vec<&str> = wcc.iter().copied().collect();
        if pieces.len() < 2 {
            merge_of.insert(
                format!("merge_{}", serial),
                seq_of.get(&String::from(*pieces.first().unwrap())).unwrap().clone(),
            );
            serial += 1;
            continue;
        }

        let mut idx_of = BTreeMap::new();
        for piece in &pieces {
            let idx = topo_sorted.iter().position(|r| r == piece).unwrap();
            idx_of.insert(piece, idx);
        }
        let mut sorted = pieces.clone();
        sorted.sort_by_cached_key(|k| idx_of.get(k).unwrap());

        for i in 0..(sorted.len() - 1) {
            let p0 = sorted[i];
            let p1 = sorted[i + 1];

            if i == 0 {
                merge_of.insert(format!("merge_{}", serial), seq_of.get(p0).unwrap().clone());
            }

            if let Some(&weight) = graph.edge_weight(p0, p1) {
                let seq = seq_of.get(p1).unwrap();
                let len = seq.len();
                let ss: String = seq
                    .chars()
                    .skip(len - weight as usize)
                    .take(weight as usize)
                    .collect();
                merge_of
                    .get_mut(&format!("merge_{}", serial))
                    .unwrap()
                    .push_str(&ss);
            } else {
                serial += 1;
                merge_of.insert(format!("merge_{}", serial), seq_of.get(p1).unwrap().clone());
            }
        }
        serial += 1;
    }
    let mut fa_out = fasta::Writer::new(intspan::writer("merged.fasta"));
    for (k, v) in merge_of.iter() {
        let record = fasta::Record::with_attrs(k, None, v.as_ref());
        fa_out.write_record(&record)?;
    }

    run_cmd!(info "==> Outputs")?;
    run_cmd!(
        cat non-overlapped.fasta merged.fasta |
            faops filter -l 0 stdin ${abs_outfile}
    )?;

    //----------------------------
    // Done
    //----------------------------
    env::set_current_dir(&curdir)?;

    Ok(())
}

fn g2gv(g: &GraphMap<&str, i32, Directed>, file: &str) -> anyhow::Result<()> {
    let mut dot = "digraph {\n".to_string();

    for edge in g.all_edges() {
        dot += &format!(
            "    {} -> {} [label=\"{}\"];\n",
            edge.source(),
            edge.target(),
            edge.weight(),
        );
    }
    dot += "}\n";

    let mut parser = layout::gv::DotParser::new(&dot);
    let tree = parser.process().unwrap();
    let mut gb = layout::gv::GraphBuilder::new();
    gb.visit_graph(&tree);
    let mut vg = gb.get();
    let mut svg = layout::backends::svg::SVGWriter::new();
    vg.do_it(false, false, false, &mut svg);
    let content = svg.finalize();
    layout::core::utils::save_to_file(&format!("{}.svg", file), &content)?;

    Ok(())
}

// use std::io::{Read, Write};
//
// fn pause() {
//     let mut stdin = std::io::stdin();
//     let mut stdout = std::io::stdout();
//
//     // We want the cursor to stay at the end of the line, so we print without a newline and flush manually.
//     write!(stdout, "Press any key to continue...").unwrap();
//     stdout.flush().unwrap();
//
//     // Read a single byte and discard
//     let _ = stdin.read(&mut [0u8]).unwrap();
// }

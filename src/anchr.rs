extern crate clap;

use clap::*;

mod cmd;
mod libs;

fn main() -> anyhow::Result<()> {
    let app = Command::new("anchr")
        .version(crate_version!())
        .author(crate_authors!())
        .about("Anchr - the Assembler of N-free CHRomosomes")
        .propagate_version(true)
        .arg_required_else_help(true)
        .subcommand(cmd::anchors::make_subcommand())
        .subcommand(cmd::contained::make_subcommand())
        .subcommand(cmd::covered::make_subcommand())
        .subcommand(cmd::dazzname::make_subcommand())
        .subcommand(cmd::dep::make_subcommand())
        .subcommand(cmd::ena::make_subcommand())
        .subcommand(cmd::merge::make_subcommand())
        .subcommand(cmd::mergeread::make_subcommand())
        .subcommand(cmd::orient::make_subcommand())
        .subcommand(cmd::overlap::make_subcommand())
        .subcommand(cmd::paf2ovlp::make_subcommand())
        .subcommand(cmd::quorum::make_subcommand())
        .subcommand(cmd::restrict::make_subcommand())
        .subcommand(cmd::show2ovlp::make_subcommand())
        .subcommand(cmd::template::make_subcommand())
        .subcommand(cmd::trim::make_subcommand())
        .subcommand(cmd::unitigs::make_subcommand())
        .after_help(
            r###"
Subcommand groups:

* Dependence
    * dep check / dep install
* Download
    * ena info / ena prep
* Overlaps
    * Standalone
        * dazzname / show2ovlp / paf2ovlp / covered / restrict
    * Daligner pipelines
        * overlap / orient / contained / merge
        * overlap2
        * group
        * layout
* Assembling
    * trim / quorum / mergeread / unitigs / anchors
    * template

"###,
        );

    // Check which subcomamnd the user ran...
    match app.get_matches().subcommand() {
        // Dependence
        Some(("dep", sub_matches)) => cmd::dep::execute(sub_matches),
        // Download
        Some(("ena", sub_matches)) => cmd::ena::execute(sub_matches),
        // Overlaps - Standalone
        Some(("dazzname", sub_matches)) => cmd::dazzname::execute(sub_matches),
        Some(("show2ovlp", sub_matches)) => cmd::show2ovlp::execute(sub_matches),
        Some(("paf2ovlp", sub_matches)) => cmd::paf2ovlp::execute(sub_matches),
        Some(("covered", sub_matches)) => cmd::covered::execute(sub_matches),
        Some(("restrict", sub_matches)) => cmd::restrict::execute(sub_matches),
        // Overlaps - Daligner pipelines
        Some(("overlap", sub_matches)) => cmd::overlap::execute(sub_matches),
        Some(("orient", sub_matches)) => cmd::orient::execute(sub_matches),
        Some(("contained", sub_matches)) => cmd::contained::execute(sub_matches),
        Some(("merge", sub_matches)) => cmd::merge::execute(sub_matches),
        // Assembling
        Some(("trim", sub_matches)) => cmd::trim::execute(sub_matches),
        Some(("quorum", sub_matches)) => cmd::quorum::execute(sub_matches),
        Some(("mergeread", sub_matches)) => cmd::mergeread::execute(sub_matches),
        Some(("unitigs", sub_matches)) => cmd::unitigs::execute(sub_matches),
        Some(("anchors", sub_matches)) => cmd::anchors::execute(sub_matches),
        Some(("template", sub_matches)) => cmd::template::execute(sub_matches),
        _ => unreachable!(),
    }
    .unwrap();

    Ok(())
}

// TODO:
//  Replace `tsv-utils` with `rgr`

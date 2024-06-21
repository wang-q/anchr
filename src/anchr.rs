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
        .subcommand(cmd::covered::make_subcommand())
        .subcommand(cmd::dep::make_subcommand())
        .subcommand(cmd::ena::make_subcommand())
        .subcommand(cmd::merge::make_subcommand())
        .subcommand(cmd::paf2ovlp::make_subcommand())
        .subcommand(cmd::quorum::make_subcommand())
        .subcommand(cmd::restrict::make_subcommand())
        .subcommand(cmd::template::make_subcommand())
        .subcommand(cmd::trim::make_subcommand())
        .subcommand(cmd::unitigs::make_subcommand())
        .after_help(
            r###"
Subcommand groups:

* Overlaps
    * Standalone
        * paf2ovlp
        * show2ovlp
        * covered
        * restrict
        * dazzname
    * Daligner pipelines
        * overlap
        * contained
        * merge
        * orient
        * group
        * layout
        * overlap2

* Assembling
    * anchors
    * dep
    * ena
    * merge
    * quorum
    * trim
    * unitigs

"###,
        );

    // Check which subcomamnd the user ran...
    match app.get_matches().subcommand() {
        Some(("anchors", sub_matches)) => cmd::anchors::execute(sub_matches),
        Some(("covered", sub_matches)) => cmd::covered::execute(sub_matches),
        Some(("dep", sub_matches)) => cmd::dep::execute(sub_matches),
        Some(("ena", sub_matches)) => cmd::ena::execute(sub_matches),
        Some(("merge", sub_matches)) => cmd::merge::execute(sub_matches),
        Some(("paf2ovlp", sub_matches)) => cmd::paf2ovlp::execute(sub_matches),
        Some(("quorum", sub_matches)) => cmd::quorum::execute(sub_matches),
        Some(("restrict", sub_matches)) => cmd::restrict::execute(sub_matches),
        Some(("template", sub_matches)) => cmd::template::execute(sub_matches),
        Some(("trim", sub_matches)) => cmd::trim::execute(sub_matches),
        Some(("unitigs", sub_matches)) => cmd::unitigs::execute(sub_matches),
        _ => unreachable!(),
    }
    .unwrap();

    Ok(())
}

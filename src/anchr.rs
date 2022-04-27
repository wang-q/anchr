extern crate clap;
use clap::*;

mod cmd;

fn main() -> std::io::Result<()> {
    let app = Command::new("anchr")
        .version(crate_version!())
        .author(crate_authors!())
        .about("Anchr - the Assembler of N-free CHRomosomes")
        .propagate_version(true)
        .arg_required_else_help(true)
        .subcommand(cmd::anchors::make_subcommand())
        .subcommand(cmd::dep::make_subcommand())
        .subcommand(cmd::ena::make_subcommand())
        .subcommand(cmd::merge::make_subcommand())
        .subcommand(cmd::quorum::make_subcommand())
        .subcommand(cmd::template::make_subcommand())
        .subcommand(cmd::trim::make_subcommand())
        .subcommand(cmd::unitigs::make_subcommand());

    // Check which subcomamnd the user ran...
    match app.get_matches().subcommand() {
        Some(("anchors", sub_matches)) => cmd::anchors::execute(sub_matches),
        Some(("dep", sub_matches)) => cmd::dep::execute(sub_matches),
        Some(("ena", sub_matches)) => cmd::ena::execute(sub_matches),
        Some(("merge", sub_matches)) => cmd::merge::execute(sub_matches),
        Some(("quorum", sub_matches)) => cmd::quorum::execute(sub_matches),
        Some(("template", sub_matches)) => cmd::template::execute(sub_matches),
        Some(("trim", sub_matches)) => cmd::trim::execute(sub_matches),
        Some(("unitigs", sub_matches)) => cmd::unitigs::execute(sub_matches),
        _ => unreachable!(),
    }
    .unwrap();

    Ok(())
}

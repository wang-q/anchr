extern crate clap;

use clap::*;

mod cmd;
mod libs;

fn main() -> anyhow::Result<()> {
    // Default to `info` level so fq/asm progress messages remain visible by
    // default, matching pgr; users can override via RUST_LOG.
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let app = Command::new("anchr")
        .version(crate_version!())
        .author(crate_authors!())
        .about("Anchr - the Assembler of N-free CHRomosomes")
        .propagate_version(true)
        .arg_required_else_help(true)
        .subcommand(cmd::asm::make_subcommand())
        .subcommand(cmd::dep::make_subcommand())
        .subcommand(cmd::ena::make_subcommand())
        .subcommand(cmd::fq::make_subcommand())
        .subcommand(cmd::mergeread::make_subcommand())
        .subcommand(cmd::quorum::make_subcommand())
        .subcommand(cmd::sam::make_subcommand())
        .subcommand(cmd::template::make_subcommand())
        .subcommand(cmd::trim::make_subcommand())
        .after_help(
            r###"
Subcommand groups:

* Dependence
    * dep check / dep install
* Download
    * ena meta / ena manifest
* Assembling
    * trim / quorum / mergeread
    * template

"###,
        );

    // Check which subcomamnd the user ran...
    match app.get_matches().subcommand() {
        // Reads processing (migrated from pgr)
        Some(("fq", sub_matches)) => cmd::fq::execute(sub_matches),
        // Assembly (migrated from pgr)
        Some(("asm", sub_matches)) => cmd::asm::execute(sub_matches),
        // Dependence
        Some(("dep", sub_matches)) => cmd::dep::execute(sub_matches),
        // Download
        Some(("ena", sub_matches)) => cmd::ena::execute(sub_matches),
        Some(("sam", sub_matches)) => cmd::sam::execute(sub_matches),
        // Assembling
        Some(("trim", sub_matches)) => cmd::trim::execute(sub_matches),
        Some(("quorum", sub_matches)) => cmd::quorum::execute(sub_matches),
        Some(("mergeread", sub_matches)) => cmd::mergeread::execute(sub_matches),
        Some(("template", sub_matches)) => cmd::template::execute(sub_matches),
        _ => unreachable!(),
    }?;

    Ok(())
}

// TODO:
//  Replace `tsv-utils` with `rgr`

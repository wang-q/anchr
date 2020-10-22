extern crate clap;
use clap::*;

mod cmd;

fn main() -> std::io::Result<()> {
    let app = App::new("anchr")
        .version(crate_version!())
        .author(crate_authors!())
        .about("Anchr - Assembler of N-free CHRomosomes")
        .setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(cmd::template::make_subcommand())
        .subcommand(cmd::quorum::make_subcommand())
        .subcommand(cmd::trim::make_subcommand())
        .subcommand(cmd::dep::make_subcommand());

    // Check which subcomamnd the user ran...
    match app.get_matches().subcommand() {
        ("template", Some(sub_matches)) => cmd::template::execute(sub_matches),
        ("quorum", Some(sub_matches)) => cmd::quorum::execute(sub_matches),
        ("trim", Some(sub_matches)) => cmd::trim::execute(sub_matches),
        ("dep", Some(sub_matches)) => cmd::dep::execute(sub_matches),
        (_, _) => unreachable!(),
    }?;

    Ok(())
}

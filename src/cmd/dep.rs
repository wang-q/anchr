use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand<'a>() -> Command<'a> {
    Command::new("dep")
        .about("Dependencies")
        .after_help(
            r#"
* check   - check dependencies
* install - install dependencies
"#,
        )
        .arg(
            Arg::new("infile")
                .help("Sets the input file to use")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("outfile")
                .short('o')
                .long("outfile")
                .takes_value(true)
                .default_value("stdout")
                .forbid_empty_values(true)
                .help("Output filename. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> std::result::Result<(), std::io::Error> {
    let mut writer = intspan::writer(args.value_of("outfile").unwrap());

    let kb = match args.value_of("infile").unwrap() {
        "check" => include_str!("../../templates/check_dep.sh"),
        "install" => include_str!("../../templates/install_dep.sh"),
        _ => unreachable!(),
    };

    writer.write_all(kb.as_ref())?;

    Ok(())
}

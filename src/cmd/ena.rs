use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand<'a>() -> Command<'a> {
    Command::new("ena")
        .about("ENA scripts")
        .after_help(
            r#"
* info - Grab information from ENA
* prep - Create downloading scripts
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
pub fn execute(args: &ArgMatches) -> std::result::Result<(), Box<dyn std::error::Error>> {
    let mut writer = intspan::writer(args.value_of("outfile").unwrap());

    let kb = match args.value_of("infile").unwrap() {
        "info" => include_str!("../../templates/ena_info.pl"),
        "prep" => include_str!("../../templates/ena_prep.pl"),
        _ => unreachable!(),
    };

    writer.write_all(kb.as_ref())?;

    Ok(())
}

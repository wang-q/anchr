use clap::*;

// Create clap subcommand arguments
pub fn make_subcommand() -> Command {
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
                .num_args(1)
                .default_value("stdout")
                .help("Output filename. [stdout] for screen"),
        )
}

// command implementation
pub fn execute(args: &ArgMatches) -> anyhow::Result<()> {
    let mut writer = intspan::writer(args.get_one::<String>("outfile").unwrap());

    let kb = match args.get_one::<String>("infile").unwrap().as_ref() {
        "check" => include_str!("../../templates/check_dep.sh"),
        "install" => include_str!("../../templates/install_dep.sh"),
        _ => unreachable!(),
    };

    writer.write_all(kb.as_ref())?;

    Ok(())
}

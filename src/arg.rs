use ansi_term::Colour;
use clap_v3::{App, Arg, ArgMatches};

pub fn get_arg() -> ArgMatches {
    println!(
        "{} {}",
        Colour::Green.paint("Welcome to use kallisto2matrix!"),
        Colour::Red.paint("Author: Yu kaiquan <1962568272@qq.com>")
    );
    App::new("kallisto2matrix")
        .version("0.0.2")
        .author("Yu kaiquan <1962568272@qq.com>")
        .about("Convert kallisto/salmon output to Count and TPM matrix")
        .arg(
            Arg::with_name("input")
                .short('i')
                .long("input")
                .value_name("FILE")
                .help("kallisto output file")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("output")
                .short('o')
                .long("output")
                .value_name("FILE")
                .help("output file")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("type")
                .short('t')
                .long("type")
                .value_name("STRING")
                .help("input file type eg: kallisto or salmon (default: kallisto)")
                .takes_value(true)
                .required(false),
        )
        .arg(
            Arg::with_name("threads")
                .short('n')
                .long("threads")
                .value_name("STRING")
                .help("threads (default: 1)")
                .takes_value(true)
                .required(false),
        )
        .get_matches()
}

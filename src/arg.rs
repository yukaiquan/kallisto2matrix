use ansi_term::Colour;

use clap::{Arg, ArgAction, Command, ArgMatches};


pub fn get_arg() -> ArgMatches {
    println!(
        "{} {}",
        Colour::Green.paint("Welcome to use kallisto2matrix!"),
        Colour::Red.paint("Author: Yu kaiquan <1962568272@qq.com>")
    );


    Command::new("kallisto2matrix")
        .version("0.0.2")
        .author("Yu kaiquan <1962568272@qq.com>")
        .about("Convert kallisto/salmon output to Count and TPM matrix")
        // 核心参数：input（-i/--input），必填
        .arg(
            Arg::new("input")
                .short('i')
                .long("input")
                .value_name("FILE")
                .help("kallisto/salmon output file list (csv format, e.g., ./file.tsv,sample1)")
                .action(ArgAction::Set) 
                .required(true),      
        )
        // 核心参数：output（-o/--output），必填
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_name("PREFIX")
                .help("Output prefix for matrix files (e.g., ./result will generate result_count_matrix.txt)")
                .action(ArgAction::Set)
                .required(true),
        )
        // 可选参数：type（-t/--type），默认kallisto
        .arg(
            Arg::new("type")
                .short('t')
                .long("type")
                .value_name("STRING")
                .help("Input file type: kallisto or salmon (default: kallisto)")
                .action(ArgAction::Set)
                .default_value("kallisto"),  // 直接设置默认值，简化主程序逻辑
        )
        // 可选参数：threads（-n/--threads），默认1
        .arg(
            Arg::new("threads")
                .short('n')
                .long("threads")
                .value_name("UINT")
                .help("Number of threads to use (default: 1)")
                .action(ArgAction::Set)
                .default_value("1"),        // 直接设置默认值，避免主程序unwrap_or
        )
        .get_matches()
}
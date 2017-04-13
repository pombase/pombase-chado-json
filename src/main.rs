extern crate postgres;
extern crate regex;
extern crate getopts;
extern crate bit_set;
extern crate serde;
extern crate serde_json;
extern crate serde_yaml;

use postgres::{Connection, TlsMode};

use std::env;
use getopts::Options;
use std::process;
use std::fs;
use std::fs::File;
use std::io::BufReader;

extern crate pombase;

use pombase::db::*;
use pombase::web::data_build::*;

const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");
const VERSION: &'static str = env!("CARGO_PKG_VERSION");

fn make_subdirs(output_dir: &str) {
    let subdirs = vec!["gene", "genotype", "term", "reference"];

    for subdir in &subdirs {
        let dir = String::new() + output_dir + "/" + subdir;
        fs::create_dir_all(&dir).unwrap_or_else(|why| {
            println!("Creating output directory failed: {:?}", why.kind());
        });
    }
}

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} [options]", program);
    print!("{}", opts.usage(&brief));
}

fn read_config(config_file_name: &str) -> Config {
    let file = match File::open(config_file_name) {
        Ok(file) => file,
        Err(err) => {
            print!("Failed to read {}: {}\n", config_file_name, err);
            process::exit(1);
        }
    };
    let reader = BufReader::new(file);

    match serde_json::from_reader(reader) {
        Ok(config) => config,
        Err(err) => {
            print!("failed to parse {}: {}", config_file_name, err);
            process::exit(1);
        },
    }
}

fn main() {
    print!("{} v{}\n", PKG_NAME, VERSION);

    let args: Vec<String> = env::args().collect();
    let mut opts = Options::new();

    opts.optflag("h", "help", "print this help message");
    opts.optopt("c", "config-file", "Configuration file name", "CONFIG");
    opts.optopt("p", "postgresql-connection-string",
                "PostgresSQL connection string like: postgres://user:pass@host/db_name",
                "CONN_STR");
    opts.optopt("d", "output-directory",
                "Destination directory for JSON output", "DIR");
    opts.optflag("j", "store-json",
                 "optionally create a 'web_json' schema to store the generated JSON in the database");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!("Invalid options\n{}", f)
    };

    let program = args[0].clone();

    if matches.opt_present("help") {
        print_usage(&program, opts);
        process::exit(0);
    }

    if !matches.opt_present("config-file") {
        print!("no -c|--config-file option\n");
        print_usage(&program, opts);
        process::exit(1);
    }
    if !matches.opt_present("postgresql-connection-string") {
        print!("no -p|--postgresql-connection-string option\n");
        print_usage(&program, opts);
        process::exit(1);
    }
    if !matches.opt_present("output-directory") {
        print!("no -d|--output-directory option\n");
        print_usage(&program, opts);
        process::exit(1);
    }

    let config = read_config(&matches.opt_str("c").unwrap());
    let connection_string = matches.opt_str("p").unwrap();
    let output_dir = matches.opt_str("d").unwrap();

    make_subdirs(&output_dir);

    let conn = Connection::connect(connection_string.as_str(), TlsMode::None).unwrap();

    let raw = Raw::new(&conn);
    let mut web_data_build = WebDataBuild::new(&raw, &config);
    let web_data = web_data_build.get_web_data();

    web_data.write(&output_dir);

    if matches.opt_present("store-json") {
        conn.execute("DROP SCHEMA IF EXISTS web_json CASCADE", &[]).unwrap();
        conn.execute("CREATE SCHEMA web_json", &[]).unwrap();
        conn.execute("CREATE EXTENSION IF NOT EXISTS pg_trgm;", &[]).unwrap();
        conn.execute("CREATE TABLE web_json.gene (uniquename TEXT, data JSONB)", &[]).unwrap();
        conn.execute("CREATE TABLE web_json.term (termid TEXT, data JSONB)", &[]).unwrap();
        conn.execute("CREATE TABLE web_json.reference (uniquename TEXT, data JSONB)", &[]).unwrap();

        web_data.store_jsonb(&conn);

        print!("stored results as JSONB using {}\n", &connection_string);
    }
}

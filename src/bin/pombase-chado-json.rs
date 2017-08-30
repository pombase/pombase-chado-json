extern crate postgres;
extern crate getopts;

use postgres::{Connection, TlsMode};

use std::env;
use getopts::Options;
use std::process;

extern crate pombase;

use pombase::db::*;
use pombase::web::config::*;
use pombase::web::data_build::*;
use pombase::interpro::parse_interpro;

const PKG_NAME: &str = env!("CARGO_PKG_NAME");
const VERSION: &str = env!("CARGO_PKG_VERSION");

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} [options]", program);
    print!("{}", opts.usage(&brief));
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
    opts.optopt("i", "interpro-data-file",
                "The name of the file generated by 'pombase-interpro'", "FILE");
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

    let config = Config::read(&matches.opt_str("c").unwrap());
    let connection_string = matches.opt_str("p").unwrap();
    let interpro_json = matches.opt_str("i").unwrap();
    let output_dir = matches.opt_str("d").unwrap();

    let conn = Connection::connect(connection_string.as_str(), TlsMode::None).unwrap();

    let raw = Raw::new(&conn);
    let interpro_data = parse_interpro(&config, &interpro_json);
    let web_data_build = WebDataBuild::new(&raw, &interpro_data, &config);
    let web_data = web_data_build.get_web_data();

    web_data.write(&config, &output_dir);

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

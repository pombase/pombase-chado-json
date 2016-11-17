extern crate postgres;
extern crate regex;
extern crate getopts;
#[macro_use]
extern crate serde;

use postgres::{Connection, TlsMode};

use std::env;
use getopts::Options;
use std::fs;

extern crate pombase;

use pombase::db::*;
use pombase::web::data_build::*;

const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");
const VERSION: &'static str = env!("CARGO_PKG_VERSION");

fn make_subdirs(output_dir: &str) {
    let subdirs = vec!["gene", "term", "reference"];

    for subdir in &subdirs {
        let dir = String::new() + output_dir + "/" + subdir;
        fs::create_dir_all(&dir).unwrap_or_else(|why| {
            println!("Creating output directory failed: {:?}", why.kind());
        });
    }
}

fn main() {
    print!("{} v{}\n", PKG_NAME, VERSION);

    let args: Vec<String> = env::args().collect();
    let mut opts = Options::new();

    opts.reqopt("c", "connection-string",
                "PostgresSQL connection string like: postgres://user:pass@host/db_name",
                "CONN_STR");
    opts.reqopt("d", "output-directory",
                "Destination directory for JSON output", "DIR");
    opts.reqopt("O", "organism",
                "Only output genes from this organism (eg. 'Schizosaccharomyces_pombe')", "GENUS_SPECIES");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!("Invalid options\n{}", f)
    };

    let connection_string = matches.opt_str("c").unwrap();
    let output_dir = matches.opt_str("d").unwrap();
    let organism_genus_species = matches.opt_str("O").unwrap();

    make_subdirs(&output_dir);

    let conn = Connection::connect(connection_string.as_str(), TlsMode::None).unwrap();
    let raw = Raw::new(&conn);
    let mut web_data_build = WebDataBuild::new(&raw);
    let web_data = web_data_build.get_web_data();

    web_data.write(&output_dir, &organism_genus_species);
}

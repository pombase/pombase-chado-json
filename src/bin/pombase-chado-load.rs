extern crate pombase;

use std::io::Read;

use std::fs::File;
use std::io::BufReader;
use std::str::FromStr;

use std::error::Error;
use std::env;
use std::process;

use getopts::Options;

use deadpool_postgres::{Pool, Manager};

use getopts::ParsingStyle;
use pombase::data_types::AlleleShortMap;

use pombase::types::OrganismTaxonId;

use pombase::db::Raw;
use pombase::db::Processed;
use pombase::load::Loader;

const PKG_NAME: &str = env!("CARGO_PKG_NAME");
const VERSION: &str = env!("CARGO_PKG_VERSION");

fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} [options] [action] [file_name]", program);
    print!("{}", opts.usage(&brief));
}

// if mark_as_obsolete is true, is_obsolete will be set on new features
//
// later when the allele is read (in the Perl code) from any other source,
// the is_obsolete flag will be set to true
//
// See PhenotypeFeatureFinder.pm
//
fn read_allele_map(file_name: &str) -> AlleleShortMap {
    let file = match File::open(file_name) {
        Ok(file) => file,
        Err(err) => {
            eprintln!("Failed to read {}: {}", file_name, err);
            process::exit(1);
        }
    };

    let mut reader = BufReader::new(file);

    let mut decoded_json = String::new();
    reader.read_to_string(&mut decoded_json).unwrap();
    let serde_result = serde_json::from_str(&decoded_json);

    match serde_result {
        Ok(results) => results,
        Err(err) => {
            eprint!("failed to parse {}: {}", file_name, err);
            process::exit(1);
        },
    }
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn Error>> {
    println!("Loading using {} v{}", PKG_NAME, VERSION);

    let args: Vec<String> = env::args().collect();
    let mut opts = Options::new();
    let opts = opts.parsing_style(ParsingStyle::StopAtFirstFree);

    opts.optflag("h", "help", "print this help message");
    opts.optopt("p", "postgresql-connection-string",
                "PostgresSQL connection string like: postgres://user:pass@host/db_name",
                "CONN_STR");
    opts.optopt("t", "taxonid",
                "Taxon ID of the organism to load",
                "TAXONID");

    let program = args[0].clone();

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(e) => {
            print_usage(&program, opts);
            println!("\nerror: {}", e);
            process::exit(0);
        }
    };

    if matches.opt_present("help") {
        print_usage(&program, opts);
        process::exit(0);
    }

    if !matches.opt_present("postgresql-connection-string") {
        println!("no -p|--postgresql-connection-string option");
        print_usage(&program, opts);
        process::exit(1);
    }

    if !matches.opt_present("taxonid") {
        println!("no -t|--taxonid option");
        print_usage(&program, opts);
        process::exit(1);
    }

    let mut remaining_args = matches.free.clone();

    if remaining_args.len() < 2 {
        println!("needs [action] and [file_name] arguments");
        print_usage(&program, opts);
        process::exit(1);
    }

    let action = remaining_args.remove(0);

    if action != "allele-json" {
        println!("unknown action {}", action);
        print_usage(&program, opts);
        process::exit(1);
    }

    let mark_as_obsolete =
        if let Some(first_remaining) = remaining_args.first() {
            first_remaining == "--mark-as-obsolete"
        } else {
            false
        };
    if mark_as_obsolete {
        remaining_args.remove(0);
    }

    if remaining_args.is_empty() {
        println!("allele-json needs a [file_name] argument");
        print_usage(&program, opts);
        process::exit(1);
    }

    let allele_file_name = &remaining_args[0];

    let connection_string = matches.opt_str("p").unwrap();
    let taxonid_opt = matches.opt_str("t").unwrap();

    let taxonid_result = taxonid_opt.parse();

    let taxonid: OrganismTaxonId =
        taxonid_result.unwrap_or_else(|_| panic!("failed to parse taxon ID {}", taxonid_opt));

    let pg_config = tokio_postgres::Config::from_str(&connection_string)?;
    let manager = Manager::new(pg_config, tokio_postgres::NoTls);
    let pool = Pool::builder(manager).max_size(16).build().unwrap();

    let mut client = pool.get().await?;

    let raw = Raw::new(&mut client).await?;
    let pro = Processed::new(raw);

    let client = pool.get().await?;

    let mut loader = Loader::new(client, taxonid, pro);

    let allele_map = read_allele_map(allele_file_name);

    loader.load_alleles(&allele_map, mark_as_obsolete).await?;

    Ok(())
}

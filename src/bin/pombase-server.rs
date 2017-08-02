#![feature(plugin)]
#![plugin(rocket_codegen)]

extern crate getopts;
extern crate rocket;
extern crate serde_json;
#[macro_use] extern crate rocket_contrib;

extern crate pombase;

use std::sync::Mutex;
use std::process;
use std::env;

use getopts::Options;

use rocket_contrib::{Json, Value};

use pombase::api::query::Query;
use pombase::api::result::QueryAPIResult;
use pombase::api::query_exec::QueryExec;
use pombase::api::server_data::ServerData;

const PKG_NAME: &str = env!("CARGO_PKG_NAME");
const VERSION: &str = env!("CARGO_PKG_VERSION");

#[post("/query", data="<q>", format = "application/json")]
fn query_post(q: Json<Query>, state: rocket::State<Mutex<QueryExec>>)
              -> Option<Json<QueryAPIResult>>
{
    let query_exec = state.lock().expect("failed to lock");
    Some(Json(query_exec.exec(&q.into_inner())))
}

#[get ("/reload")]
fn reload(state: rocket::State<Mutex<QueryExec>>) {
    let mut query_exec = state.lock().expect("failed to lock");
    print!("reloading ...\n");
    query_exec.reload();
    print!("... done\n");
}

#[get ("/ping")]
fn ping() -> Option<String> {
    Some(String::from("OK") + " " + PKG_NAME + " " + VERSION)
}

#[error(404)]
fn not_found() -> Json<Value> {
    Json(json!({
        "status": "error",
        "reason": "Resource was not found."
    }))
}

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
    opts.optopt("m", "search-maps", "Search data", "MAPS_JSON_FILE");
    opts.optopt("s", "gene-subsets", "Gene subset data", "SUBSETS_JSON_FILE");

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
    if !matches.opt_present("search-maps") {
        print!("no --search-maps|-m option\n");
        print_usage(&program, opts);
        process::exit(1);
    }
    if !matches.opt_present("gene-subsets") {
        print!("no --gene-subsets|-s option\n");
        print_usage(&program, opts);
        process::exit(1);
    }

    let search_maps_filename = matches.opt_str("m").unwrap();
    let gene_subsets_filename = matches.opt_str("s").unwrap();
    println!("Reading data files ...");

    let config_file_name = matches.opt_str("c").unwrap();
    let server_data = ServerData::new(&config_file_name, &search_maps_filename,
                                      &gene_subsets_filename);
    let query_exec = QueryExec::new(server_data);

    println!("Starting server ...");
    rocket::ignite()
        .mount("/", routes![query_post, reload, ping])
        .catch(errors![not_found])
        .manage(Mutex::new(query_exec))
        .launch();
}

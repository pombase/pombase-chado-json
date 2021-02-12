#![feature(nll)]
#![feature(proc_macro_hygiene, decl_macro)]

extern crate getopts;

#[macro_use] extern crate rocket;
#[macro_use] extern crate rocket_contrib;
use rocket::response::content;

#[macro_use] extern crate serde_derive;

extern crate pombase;

use std::sync::Mutex;
use std::process;
use std::env;
use std::path::{Path, PathBuf};

use getopts::Options;
use rocket_contrib::json::{Json, JsonValue};
use rocket_contrib::serve::StaticFiles;
use rocket::response::NamedFile;

use pombase::api::query::Query as PomBaseQuery;
use pombase::api::result::QueryAPIResult;
use pombase::api::search::{Search, TermSearchMatch, RefSearchMatch, DocSearchMatch,
                           SolrSearchScope};
use pombase::api::query_exec::QueryExec;
use pombase::api_data::{api_maps_from_file, APIData};
use pombase::api::site_db::SiteDB;

use pombase::data_types::{SolrTermSummary, SolrReferenceSummary, GeneDetails, GenotypeDetails,
                          FeatureShort,
                         TermDetails, ReferenceDetails};
use pombase::web::simple_pages::{render_simple_gene_page, render_simple_reference_page,
                                 render_simple_term_page};
use pombase::web::config::Config;

const PKG_NAME: &str = env!("CARGO_PKG_NAME");
const VERSION: &str = env!("CARGO_PKG_VERSION");

struct StaticFileState {
    web_root_dir: String,
}


#[derive(Serialize, Debug)]
struct TermLookupResponse {
    status: String,
    summary: Option<SolrTermSummary>,
}



// If the path is a directory, return path+"/index.html".  Otherwise
// try the path, then try path + ".json", then default to loading the
// Angular app from /index.html
#[get("/<path..>", rank=3)]
fn get_misc(mut path: PathBuf, state: rocket::State<Mutex<StaticFileState>>,
            config: rocket::State<Config>)
            -> Option<NamedFile>
{
    let web_root_dir = &state.lock().expect("failed to lock").web_root_dir;
    let root_dir_path = Path::new("/").join(web_root_dir);
    let is_jbrowse_path = path.starts_with("jbrowse/");

    if let Some(path_str) = path.to_str() {
        if let Some(new_path_str) =
            config.doc_page_aliases.get(&format!("/{}", path_str)) {
                path = PathBuf::from(new_path_str);
            }
    }

    let full_path = root_dir_path.join(path);

    if full_path.is_dir() {
        let index_path = full_path.join("index.html");
        return NamedFile::open(index_path).ok();
    }

    if full_path.exists() {
        return NamedFile::open(full_path).ok();
    }

    let mut json_path_str = full_path.to_str().unwrap().to_owned();
    json_path_str += ".json";
    let json_path: PathBuf = json_path_str.into();

    if json_path.exists() {
        return NamedFile::open(json_path).ok();
    }

    // special case for missing JBrowse files - return 4044
    if is_jbrowse_path {
        return None;
    }

    NamedFile::open(root_dir_path.join("index.html")).ok()
}

#[get("/api/v1/dataset/latest/data/gene/<id>", rank=2)]
fn get_gene(id: String, state: rocket::State<Mutex<QueryExec>>) -> Option<Json<GeneDetails>> {
    let query_exec = state.lock().expect("failed to lock");
    if let Some(gene) = query_exec.get_api_data().get_full_gene_details(&id) {
        Some(Json(gene))
    } else {
        None
    }
}

#[get("/api/v1/dataset/latest/data/genotype/<id>", rank=2)]
fn get_genotype(id: String, state: rocket::State<Mutex<QueryExec>>) -> Option<Json<GenotypeDetails>> {
    let query_exec = state.lock().expect("failed to lock");
    if let Some(genotype) = query_exec.get_api_data().get_genotype_details(&id) {
        Some(Json(genotype))
    } else {
        None
    }
}

#[get("/api/v1/dataset/latest/data/term/<id>", rank=2)]
fn get_term(id: String, state: rocket::State<Mutex<QueryExec>>) -> Option<Json<TermDetails>> {
    let query_exec = state.lock().expect("failed to lock");
    if let Some(term) = query_exec.get_api_data().get_term_details(&id) {
        Some(Json(term))
    } else {
        None
    }
}

#[get("/api/v1/dataset/latest/summary/term/<id>", rank=2)]
fn get_term_summary_by_id(id: String, state: rocket::State<Mutex<Search>>)
                          -> Option<Json<TermLookupResponse>>
{
    let search = state.lock().expect("failed to lock");
    let res = search.term_summary_by_id(&id);

    let lookup_response =
        match res {
            Ok(summary) => {
                TermLookupResponse {
                    status: "Ok".to_owned(),
                    summary: summary,
                }
            },
            Err(err) => {
                println!("{:?}", err);
                TermLookupResponse {
                    status: "Error".to_owned(),
                    summary: None,
                }
            },
        };

    Some(Json(lookup_response))
}

#[get("/api/v1/dataset/latest/data/reference/<id>", rank=2)]
fn get_reference(id: String, state: rocket::State<Mutex<QueryExec>>) -> Option<Json<ReferenceDetails>> {
    let query_exec = state.lock().expect("failed to lock");
    if let Some(reference) = query_exec.get_api_data().get_reference_details(&id) {
        Some(Json(reference))
    } else {
        None
    }
}

#[get("/api/v1/dataset/latest/data/seq_feature_page_features", rank=2)]
fn seq_feature_page_features(state: rocket::State<Mutex<QueryExec>>) -> Json<Vec<FeatureShort>> {
    let query_exec = state.lock().expect("failed to lock");
    Json(query_exec.get_api_data().seq_feature_page_features())
}

#[get("/", rank=1)]
fn get_index(state: rocket::State<Mutex<StaticFileState>>) -> Option<NamedFile> {
    let web_root_dir = &state.lock().expect("failed to lock").web_root_dir;
    let root_dir_path = Path::new("/").join(web_root_dir);
    NamedFile::open(root_dir_path.join("index.html")).ok()
}

/*
Return a simple HTML version a gene page for search engines
*/
#[get("/simple/gene/<id>", rank=1)]
fn get_simple_gene(id: String, query_exec_state: rocket::State<Mutex<QueryExec>>,
                   config: rocket::State<Config>) -> Option<content::Html<String>> {
    let query_exec = query_exec_state.lock().expect("failed to lock QueryExec");
    if let Some(gene) = query_exec.get_api_data().get_full_gene_details(&id) {
        Some(content::Html(render_simple_gene_page(&config, &gene)))
    } else {
        None
    }
}

/*
Return a simple HTML version a reference page for search engines
*/
#[get("/simple/reference/<id>", rank=1)]
fn get_simple_reference(id: String, query_exec_state: rocket::State<Mutex<QueryExec>>,
                        config: rocket::State<Config>) -> Option<content::Html<String>> {
    let query_exec = query_exec_state.lock().expect("failed to lock QueryExec");
    if let Some(reference) = query_exec.get_api_data().get_reference_details(&id) {
        Some(content::Html(render_simple_reference_page(&config, &reference)))
    } else {
        None
    }
}
/*
Return a simple HTML version a term page for search engines
*/
#[get("/simple/term/<id>", rank=1)]
fn get_simple_term(id: String, query_exec_state: rocket::State<Mutex<QueryExec>>,
                   config: rocket::State<Config>) -> Option<content::Html<String>> {
    let query_exec = query_exec_state.lock().expect("failed to lock QueryExec");
    if let Some(term) = query_exec.get_api_data().get_term_details(&id) {
        Some(content::Html(render_simple_term_page(&config, &term)))
    } else {
        None
    }
}

#[post("/api/v1/dataset/latest/query", rank=1, data="<q>", format = "application/json")]
fn query_post(q: Json<PomBaseQuery>, state: rocket::State<Mutex<QueryExec>>)
              -> Option<Json<QueryAPIResult>>
{
    let query_exec = state.lock().expect("failed to lock");
    Some(Json(query_exec.exec(&q.into_inner())))
}

#[derive(Serialize, Debug)]
struct TermCompletionResponse {
    status: String,
    matches: Vec<SolrTermSummary>,
}

#[derive(Serialize, Debug)]
struct RefCompletionResponse {
    status: String,
    matches: Vec<SolrReferenceSummary>,
}

#[derive(Serialize, Debug)]
struct SolrSearchResponse  {
    status: String,
    term_matches: Vec<TermSearchMatch>,
    ref_matches: Vec<RefSearchMatch>,
    doc_matches: Vec<DocSearchMatch>,
}

#[get ("/api/v1/dataset/latest/complete/term/<cv_name>/<q>", rank=1)]
fn term_complete(cv_name: String, q: String, state: rocket::State<Mutex<Search>>)
              -> Option<Json<TermCompletionResponse>>
{
    let search = state.lock().expect("failed to lock");
    let res = search.term_complete(&cv_name, &q);

    let completion_response =
        match res {
            Ok(matches) => {
                TermCompletionResponse {
                    status: "Ok".to_owned(),
                    matches: matches,
                }
            },
            Err(err) => {
                println!("{:?}", err);
                TermCompletionResponse {
                    status: "Error".to_owned(),
                    matches: vec![],
                }
            },
        };

    Some(Json(completion_response))
}

#[get ("/api/v1/dataset/latest/complete/ref/<q>", rank=1)]
fn ref_complete(q: String, state: rocket::State<Mutex<Search>>)
                -> Option<Json<RefCompletionResponse>>
{
    let search = state.lock().expect("failed to lock");
    let res = search.ref_complete(&q);

    let completion_response =
        match res {
            Ok(matches) => {
                RefCompletionResponse {
                    status: "Ok".to_owned(),
                    matches: matches,
                }
            },
            Err(err) => {
                println!("{:?}", err);
                RefCompletionResponse {
                    status: "Error".to_owned(),
                    matches: vec![],
                }
            },
        };

    Some(Json(completion_response))
}

// search for terms, refs or docs that match the query
#[get ("/api/v1/dataset/latest/search/<scope>/<q>", rank=1)]
fn solr_search(scope: String, q: String, state: rocket::State<Mutex<Search>>)
    -> Option<Json<SolrSearchResponse>>
{
    let search = state.lock().expect("failed to lock");

    if let Some(parsed_scope) = SolrSearchScope::new_from_str(&scope) {
        let search_result = search.solr_search(&parsed_scope, &q);

        match search_result {
            Ok(search_all_result) => {
                Some(Json(SolrSearchResponse {
                    status: "Ok".to_owned(),
                    term_matches: search_all_result.term_matches,
                    ref_matches: search_all_result.ref_matches,
                    doc_matches: search_all_result.doc_matches,
                }))
            },
            Err(err) => {
                Some(Json(SolrSearchResponse {
                    status: err,
                    term_matches: vec![],
                    ref_matches: vec![],
                    doc_matches: vec![],
                }))
            },
        }
    } else {
        Some(Json(SolrSearchResponse {
            status: format!("no such search scope: {}", scope),
            term_matches: vec![],
            ref_matches: vec![],
            doc_matches: vec![],
        }))
    }
}

#[get ("/api/v1/dataset/latest/motif_search/<scope>/<q>", rank=1)]
fn motif_search(scope: String, q: String, state: rocket::State<Mutex<Search>>)
                -> Option<String>
{
    let search = state.lock().expect("failed to lock");
    let res = search.motif_search(&scope, &q);

    match res {
        Ok(search_result) => {
            Some(search_result)
        },
        Err(err) => {
            println!("Motif search error: {:?}", err);
            None
        },
    }
}

#[get ("/ping", rank=1)]
fn ping() -> Option<String> {
    Some(String::from("OK") + " " + PKG_NAME + " " + VERSION)
}

#[catch(404)]
fn not_found() -> JsonValue {
    json!({
        "status": "error",
        "reason": "Resource was not found."
    })
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
    opts.optopt("w", "web-root-dir", "Root web data directory", "WEB_ROOT_DIR");
    opts.optopt("s", "site-db", "Connection string for the site local databae", "SITE_DB");

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
    if !matches.opt_present("web-root-dir") {
        print!("no --web-root-dir|-w option\n");
        print_usage(&program, opts);
        process::exit(1);
    }

    let site_db_conn_string = matches.opt_str("s");

    let search_maps_filename = matches.opt_str("m").unwrap();
    println!("Reading data files ...");

    let config_file_name = matches.opt_str("c").unwrap();
    let config = Config::read(&config_file_name);
    let api_maps = api_maps_from_file(&search_maps_filename);
    let api_data = APIData::new(&config, &api_maps);
    let site_db = site_db_conn_string.map(|conn_str| SiteDB::new(&conn_str));
    let query_exec = QueryExec::new(api_data, site_db);
    let searcher = Search::new(&config);

    let web_root_dir = matches.opt_str("w").unwrap();
    let static_file_state = StaticFileState {
        web_root_dir: web_root_dir.clone(),
    };

    println!("Starting server ...");
    rocket::ignite()
        .mount("/", routes![get_index, get_misc, query_post,
                            get_gene, get_genotype, get_term, get_reference,
                            get_simple_gene, get_simple_reference, get_simple_term,
                            get_term_summary_by_id,
                            seq_feature_page_features,
                            term_complete, ref_complete,
                            solr_search, motif_search, ping])
        .mount("/jbrowse",
               StaticFiles::from(Path::new("/")
                                 .join(&web_root_dir)
                                 .join("jbrowse")))
        .register(catchers![not_found])
        .manage(Mutex::new(query_exec))
        .manage(Mutex::new(searcher))
        .manage(Mutex::new(static_file_state))
        .manage(config)
        .launch();
}

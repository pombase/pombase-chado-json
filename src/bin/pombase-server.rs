extern crate getopts;

#[macro_use] extern crate rocket;
use pombase::data_types::AlleleDetails;
use rocket::fs::NamedFile;
use rocket::response::content;
use rocket::response::content::RawHtml;

#[macro_use] extern crate serde_derive;

extern crate pombase;

use std::process;
use std::env;
use std::path::{Path, PathBuf};

use getopts::Options;

use rocket::serde::json::{Json, Value, json};
use rocket::fs::FileServer;

use pombase::api::query::Query as PomBaseQuery;
use pombase::api::result::QueryAPIResult;
use pombase::api::search::{Search, DocSearchMatch,
                           SolrSearchScope, PNGPlot};
use pombase::api::query_exec::QueryExec;
use pombase::api_data::{api_maps_from_file, APIData};
use pombase::api::site_db::SiteDB;

use pombase::data_types::{SolrTermSummary, SolrReferenceSummary, SolrAlleleSummary,
                          GeneDetails, GenotypeDetails, FeatureShort,
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
async fn get_misc(mut path: PathBuf, state: &rocket::State<StaticFileState>,
            config: &rocket::State<Config>)
            -> Option<NamedFile>
{
    let web_root_dir = &state.web_root_dir;
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
        return NamedFile::open(index_path).await.ok();
    }

    if full_path.exists() {
        return NamedFile::open(full_path).await.ok();
    }

    let mut json_path_str = full_path.to_str().unwrap().to_owned();
    json_path_str += ".json";
    let json_path: PathBuf = json_path_str.into();

    if json_path.exists() {
        return NamedFile::open(json_path).await.ok();
    }

    // special case for missing JBrowse files - return 4044
    if is_jbrowse_path {
        return None;
    }

    NamedFile::open(root_dir_path.join("index.html")).await.ok()
}

#[get("/api/v1/dataset/latest/data/gene/<id>", rank=2)]
async fn get_gene(id: String, query_exec: &rocket::State<QueryExec>) -> Option<Json<GeneDetails>> {
    query_exec.get_api_data().get_full_gene_details(&id).map(Json)
}

#[get("/api/v1/dataset/latest/data/genotype/<id>", rank=2)]
async fn get_genotype(id: String, query_exec: &rocket::State<QueryExec>) -> Option<Json<GenotypeDetails>> {
    query_exec.get_api_data().get_genotype_details(&id).map(Json)
}

#[get("/api/v1/dataset/latest/data/allele/<id>", rank=2)]
async fn get_allele(id: String, query_exec: &rocket::State<QueryExec>) -> Option<Json<AlleleDetails>> {
    query_exec.get_api_data().get_allele_details(&id).map(Json)
}

#[get("/api/v1/dataset/latest/data/term/<id>", rank=2)]
async fn get_term(id: String, query_exec: &rocket::State<QueryExec>) -> Option<Json<TermDetails>> {
    query_exec.get_api_data().get_term_details(&id).map(Json)
}

#[get("/api/v1/dataset/latest/summary/term/<id>", rank=2)]
async fn get_term_summary_by_id(id: String, search: &rocket::State<Search>)
                          -> Option<Json<TermLookupResponse>>
{
    let res = search.term_summary_by_id(&id);

    let lookup_response =
        match res {
            Ok(summary) => {
                TermLookupResponse {
                    status: "Ok".to_owned(),
                    summary,
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
async fn get_reference(id: String, query_exec: &rocket::State<QueryExec>) -> Option<Json<ReferenceDetails>> {
    query_exec.get_api_data().get_reference_details(&id).map(Json)
}

#[get("/api/v1/dataset/latest/data/seq_feature_page_features", rank=2)]
async fn seq_feature_page_features(query_exec: &rocket::State<QueryExec>) -> Json<Vec<FeatureShort>> {
    Json(query_exec.get_api_data().seq_feature_page_features())
}

#[get("/", rank=1)]
async fn get_index(state: &rocket::State<StaticFileState>) -> Option<NamedFile> {
    let web_root_dir = &state.web_root_dir;
    let root_dir_path = Path::new("/").join(web_root_dir);
    NamedFile::open(root_dir_path.join("index.html")).await.ok()
}


/*
Return a simple HTML version a gene page for search engines
*/
#[get("/structure_view/<protein_id>", rank=1)]
async fn structure_view(protein_id: String, config: &rocket::State<Config>)
                        -> Option<content::RawHtml<String>>
{
    let search_url = config.server.django_url.to_owned() + "/structure_view/";
    let params = [("protein_id", protein_id)];
    let client = reqwest::blocking::Client::new();
    let result = client.get(&search_url).query(&params).send();

    match result {
        Ok(res) => {
            match res.text() {
                Ok(text) => Some(RawHtml(text)),
                Err(err) => {
                    eprintln!("Error proxying to Django: {:?}", err);
                    None
                }

            }
        },
        Err(err) => {
            eprintln!("Error proxying to Django: {:?}", err);
            None
        }
    }

}

/*
Return a simple HTML version a gene page for search engines
*/
#[get("/simple/gene/<id>", rank=1)]
async fn get_simple_gene(id: String, query_exec: &rocket::State<QueryExec>,
                   config: &rocket::State<Config>) -> Option<content::RawHtml<String>> {
    query_exec.get_api_data().get_full_gene_details(&id)
        .map(|gene| content::RawHtml(render_simple_gene_page(config, &gene)))
}

/*
Return a simple HTML version a reference page for search engines
*/
#[get("/simple/reference/<id>", rank=1)]
async fn get_simple_reference(id: String, query_exec: &rocket::State<QueryExec>,
                        config: &rocket::State<Config>) -> Option<content::RawHtml<String>> {
    query_exec.get_api_data().get_reference_details(&id)
        .map(|reference| content::RawHtml(render_simple_reference_page(config, &reference)))
}
/*
Return a simple HTML version a term page for search engines
*/
#[get("/simple/term/<id>", rank=1)]
async fn get_simple_term(id: String, query_exec: &rocket::State<QueryExec>,
                   config: &rocket::State<Config>) -> Option<content::RawHtml<String>> {
    query_exec.get_api_data().get_term_details(&id)
        .map(|term| content::RawHtml(render_simple_term_page(config, &term)))
}

#[post("/api/v1/dataset/latest/query", rank=1, data="<q>", format = "application/json")]
async fn query_post(q: Json<PomBaseQuery>, query_exec: &rocket::State<QueryExec>)
              -> Option<Json<QueryAPIResult>>
{
    Some(Json(query_exec.exec(&q.into_inner()).await))
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
struct AlleleCompletionResponse {
    status: String,
    matches: Vec<SolrAlleleSummary>,
}

#[derive(Serialize, Debug)]
struct SolrSearchResponse  {
    status: String,
    term_matches: Vec<SolrTermSummary>,
    ref_matches: Vec<SolrReferenceSummary>,
    doc_matches: Vec<DocSearchMatch>,
}

#[get ("/api/v1/dataset/latest/complete/term/<cv_name>/<q>", rank=1)]
async fn term_complete(cv_name: String, q: String, search: &rocket::State<Search>)
              -> Option<Json<TermCompletionResponse>>
{
    let res = search.term_complete(&cv_name, &q);

    let completion_response =
        match res {
            Ok(matches) => {
                TermCompletionResponse {
                    status: "Ok".to_owned(),
                    matches,
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
async fn ref_complete(q: String, search: &rocket::State<Search>)
                -> Option<Json<RefCompletionResponse>>
{
    let res = search.ref_complete(&q);

    let completion_response =
        match res {
            Ok(matches) => {
                RefCompletionResponse {
                    status: "Ok".to_owned(),
                    matches,
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

#[get ("/api/v1/dataset/latest/complete/allele/<q>", rank=1)]
async fn allele_complete(q: String, search: &rocket::State<Search>)
                -> Option<Json<AlleleCompletionResponse>>
{
    let res = search.allele_complete(&q).await;

    let completion_response =
        match res {
            Ok(matches) => {
                AlleleCompletionResponse {
                    status: "Ok".to_owned(),
                    matches,
                }
            },
            Err(err) => {
                println!("{:?}", err);
                AlleleCompletionResponse {
                    status: "Error".to_owned(),
                    matches: vec![],
                }
            },
        };

    Some(Json(completion_response))
}

// search for terms, refs or docs that match the query
#[get ("/api/v1/dataset/latest/search/<scope>/<q>", rank=1)]
async fn solr_search(scope: String, q: String, search: &rocket::State<Search>)
    -> Option<Json<SolrSearchResponse>>
{
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
async fn motif_search(scope: String, q: String, search: &rocket::State<Search>)
                -> Option<String>
{
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


#[get ("/api/v1/dataset/latest/gene_ex_violin_plot/<plot_size>/<genes>", rank=1)]
async fn gene_ex_violin_plot(plot_size: String, genes: String,
                       search: &rocket::State<Search>)
                       -> Option<PNGPlot>
{
    let res = search.gene_ex_violin_plot(&plot_size, &genes);

    match res {
        Ok(png_plot) => {
            Some(png_plot)
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
fn not_found() -> Value {
    json!({
        "status": "error",
        "reason": "Resource was not found."
    })
}

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} [options]", program);
    print!("{}", opts.usage(&brief));
}

#[launch]
async fn rocket() -> _ {
    println!("{} v{}", PKG_NAME, VERSION);

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
        println!("no -c|--config-file option");
        print_usage(&program, opts);
        process::exit(1);
    }
    if !matches.opt_present("search-maps") {
        println!("no --search-maps|-m option");
        print_usage(&program, opts);
        process::exit(1);
    }
    if !matches.opt_present("web-root-dir") {
        println!("no --web-root-dir|-w option");
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
    let site_db =
        if let Some(conn_str) = site_db_conn_string {
            match SiteDB::new(&conn_str).await {
                Ok(site_db) => Some(site_db),
                Err(err) => panic!("failed to connect to site db: {}", err),
            }
        } else {
            None
        };


    let query_exec = QueryExec::new(api_data, site_db);
    let searcher = Search::new(&config);

    let web_root_dir = matches.opt_str("w").unwrap();
    let static_file_state = StaticFileState {
        web_root_dir: web_root_dir.clone(),
    };

    println!("Starting server ...");
    rocket::build()
        .mount("/", routes![get_index, get_misc, query_post,
                            get_gene, get_genotype, get_allele, get_term, get_reference,
                            get_simple_gene, get_simple_reference, get_simple_term,
                            get_term_summary_by_id,
                            seq_feature_page_features,
                            term_complete, ref_complete, allele_complete,
                            solr_search, motif_search, gene_ex_violin_plot,
                            structure_view,
                            ping])
        .mount("/jbrowse",
               FileServer::from(Path::new("/")
                                 .join(&web_root_dir)
                                 .join("jbrowse")))
        .register("/", catchers![not_found])
        .manage(query_exec)
        .manage(searcher)
        .manage(static_file_state)
        .manage(config)
}

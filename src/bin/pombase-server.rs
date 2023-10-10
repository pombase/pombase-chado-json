extern crate getopts;

use axum::{
    routing::{get, post},
    Json, Router, extract::{State, Path},
    http::{StatusCode, Response, Uri, Request, HeaderMap, header}, response::{Html, IntoResponse},
    body::{boxed, BoxBody, Body, Full},
    ServiceExt,
};

use percent_encoding::{utf8_percent_encode, AsciiSet, CONTROLS};

const FRAGMENT: &AsciiSet = &CONTROLS.add(b' ').add(b'"').add(b'<').add(b'>').add(b'`');

use bytes::Bytes;

use tower::layer::Layer;

use tower_http::services::ServeDir;
use tower_http::normalize_path::NormalizePathLayer;

use serde::Serialize;
use serde_json::{json, Value};

use pombase::data_types::ProteinViewType;

use rusqlite::Connection;

extern crate serde_derive;

extern crate pombase;

use std::{process, sync::Arc};
use std::env;

use getopts::Options;


use pombase::api::query::Query;

use pombase::api::search::{Search, DocSearchMatch, SolrSearchScope};
use pombase::api::query_exec::QueryExec;
use pombase::api_data::{api_maps_from_file, APIData};
use pombase::api::site_db::SiteDB;

use pombase::data_types::{SolrTermSummary, SolrReferenceSummary, SolrAlleleSummary};
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

async fn get_static_file(path: &str) -> Result<Response<BoxBody>, (StatusCode, String)> {
    let path = utf8_percent_encode(path, FRAGMENT).to_string();
    let uri: Uri = path.parse().unwrap();
    let req = Request::builder().uri(uri).body(Body::empty()).unwrap();

    use tower::ServiceExt;

    match ServeDir::new("/").oneshot(req).await {
        Ok(res) => Ok(res.map(boxed)),
        Err(err) => Err((
            StatusCode::INTERNAL_SERVER_ERROR,
            format!("Something went wrong: {}", err),
        )),
    }
}

struct AllState {
    query_exec: QueryExec,
    search: Search,
    static_file_state: StaticFileState,
    config: Config
}

// If the path is a directory, return path+"/index.html".  Otherwise
// try the path, then try path + ".json", then default to loading the
// Angular app from /index.html
async fn get_misc(Path(mut path): Path<String>,
                  State(all_state): State<Arc<AllState>>)
            -> Result<Response<BoxBody>, (StatusCode, String)>
{
    let static_file_state = &all_state.static_file_state;
    let config = &all_state.config;
    let web_root_dir = &static_file_state.web_root_dir;
    let is_jbrowse_path = path.starts_with("jbrowse/");

    let alias = config.doc_page_aliases.get(&format!("/{}", path));

    if let Some(new_path_str) = alias {
        path = new_path_str.clone();
    }

    let full_path = format!("{}/{}", web_root_dir, path);

    if std::path::Path::new(&full_path).is_dir() {
        let index_path = format!("{}/index.html", full_path);
        return get_static_file(&index_path).await;
    }

    if std::path::Path::new(&full_path).exists() {
        return get_static_file(&full_path).await;
    }

    let json_path = format!("{}.json", full_path);

    if std::path::Path::new(&json_path).exists() {
        return get_static_file(&json_path).await;
    }

    // special case for missing JBrowse files - return 404
    if is_jbrowse_path {
        return Err((StatusCode::NOT_FOUND, format!("File not found: {}", full_path)))
    }

    let file_name = format!("{}/index.html", web_root_dir);
    get_static_file(&file_name).await
}

fn option_json_to_result<T>(id: &str, opt: Option<Json<T>>) -> Result<(StatusCode, Json<T>), (StatusCode, String)> {
    if let Some(res) = opt {
        Ok((StatusCode::OK, res))
    } else {
        Err((StatusCode::NOT_FOUND, format!("no page for: {}", id)))
    }
}

async fn get_gene(Path(id): Path<String>, State(all_state): State<Arc<AllState>>) -> impl IntoResponse {
    let res = all_state.query_exec.get_api_data().get_full_gene_details(&id).map(Json);
    option_json_to_result(&id, res)
}

async fn get_genotype(Path(id): Path<String>, State(all_state): State<Arc<AllState>>) -> impl IntoResponse {
    let res = all_state.query_exec.get_api_data().get_genotype_details(&id).map(Json);
    option_json_to_result(&id, res)
}

async fn get_allele(Path(id): Path<String>, State(all_state): State<Arc<AllState>>) -> impl IntoResponse {
    let res = all_state.query_exec.get_api_data().get_allele_details(&id).map(Json);
    option_json_to_result(&id, res)
}

async fn get_term(Path(id): Path<String>, State(all_state): State<Arc<AllState>>) -> impl IntoResponse {
    let res = all_state.query_exec.get_api_data().get_term_details(&id).map(Json);
    option_json_to_result(&id, res)
}

async fn get_protein_features(Path((full_or_widget, gene_uniquename)): Path<(String, String)>,
                              State(all_state): State<Arc<AllState>>)
       -> impl IntoResponse
{
    let full_or_widget =
        match ProteinViewType::try_from(full_or_widget.as_ref()) {
            Ok(val) => val,
            Err(err) => {
                eprintln!("get_protein_features(): {}", err);
                return Err((StatusCode::INTERNAL_SERVER_ERROR, format!("Internal error: {}", err)));
            }
        };

    let res = all_state.query_exec.get_api_data().get_protein_features_of_gene(full_or_widget, &gene_uniquename)
              .map(|s| s.to_owned())
              .map(Json);
    option_json_to_result(&gene_uniquename, res)
}

async fn get_term_summary_by_id(Path(id): Path<String>, State(all_state): State<Arc<AllState>>)
   -> impl IntoResponse
{
    let res = all_state.search.term_summary_by_id(&id).await;

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

    Json(lookup_response)
}

async fn get_reference(Path(id): Path<String>, State(all_state): State<Arc<AllState>>) -> impl IntoResponse {
    let res = all_state.query_exec.get_api_data().get_reference_details(&id).map(Json);
    option_json_to_result(&id, res)
}

async fn seq_feature_page_features(State(all_state): State<Arc<AllState>>) -> impl IntoResponse {
    Json(all_state.query_exec.get_api_data().seq_feature_page_features())
}

async fn get_index(State(all_state): State<Arc<AllState>>) -> Result<Response<BoxBody>, (StatusCode, String)> {
    let web_root_dir = &all_state.static_file_state.web_root_dir;
    get_static_file(&format!("{}/index.html", web_root_dir)).await
}

async fn structure_view(Path((structure_type, id)): Path<(String, String)>,
                        State(all_state): State<Arc<AllState>>)
                        -> (StatusCode, Html<String>)
{
    let search_url = all_state.config.server.django_url.to_owned() + "/structure_view/";
    let params = [("structure_type", structure_type), ("id", id)];
    let client = reqwest::Client::new();
    let result = client.get(search_url).query(&params).send().await;

    match result {
        Ok(res) => {
            match res.text().await {
                Ok(text) => (StatusCode::OK, Html(text)),
                Err(err) => {
                    let err_mess = format!("Error proxying to Django: {:?}", err);
                    eprintln!("{}", err_mess);
                    (StatusCode::INTERNAL_SERVER_ERROR, Html(err_mess))
                }

            }
        },
        Err(err) => {
            eprintln!("Error proxying to Django: {:?}", err);
            (StatusCode::INTERNAL_SERVER_ERROR, Html(err.to_string()))
        }
    }

}

async fn protein_feature_view(Path((full_or_widget, gene_uniquename)): Path<(String, String)>,
                              State(all_state): State<Arc<AllState>>)
   -> (StatusCode, Html<String>)
{
    let search_url = all_state.config.server.django_url.to_owned() + "/protein_feature_view/";
    let params = [("full_or_widget", full_or_widget),
                  ("gene_uniquename", gene_uniquename)];
    let client = reqwest::Client::new();
    let result = client.get(search_url).query(&params).send().await;

    match result {
        Ok(res) => {
            match res.text().await {
                Ok(text) => (StatusCode::OK, Html(text)),
                Err(err) => {
                    let err_mess = format!("Error proxying to Django: {:?}", err);
                    eprintln!("{}", err_mess);
                    (StatusCode::INTERNAL_SERVER_ERROR, Html(err_mess))
                }
            }
        },
        Err(err) => {
            let err_mess = format!("Error proxying to Django: {:?}", err);
            eprintln!("{}", err_mess);
            (StatusCode::INTERNAL_SERVER_ERROR, Html(err_mess))
        }
    }
}

/*
Return a simple HTML version a gene page for search engines
*/
async fn get_simple_gene(Path(id): Path<String>,
                         State(all_state): State<Arc<AllState>>) -> (StatusCode, Html<String>) {
    if let Some(gene) = all_state.query_exec.get_api_data().get_full_gene_details(&id) {
        (StatusCode::OK, Html(render_simple_gene_page(&all_state.config, &gene)))
    } else {
        (StatusCode::NOT_FOUND, Html(format!("no page for: {}", id)))
    }
}

/*
Return a simple HTML version a reference page for search engines
*/
async fn get_simple_reference(Path(id): Path<String>,
                              State(all_state): State<Arc<AllState>>) -> (StatusCode, Html<String>) {
    if let Some(reference) = all_state.query_exec.get_api_data().get_reference_details(&id) {
        (StatusCode::OK, Html(render_simple_reference_page(&all_state.config, &reference)))
    } else {
        (StatusCode::NOT_FOUND, Html(format!("no page for: {}", id)))
    }
}
/*
Return a simple HTML version a term page for search engines
*/
async fn get_simple_term(Path(id): Path<String>,
                         State(all_state): State<Arc<AllState>>) -> (StatusCode, Html<String>) {
    if let Some(term) = all_state.query_exec.get_api_data().get_term_details(&id) {
        (StatusCode::OK, Html(render_simple_term_page(&all_state.config, &term)))
    } else {
        (StatusCode::NOT_FOUND, Html(format!("no page for: {}", id)))
    }
}

async fn query_post(State(all_state): State<Arc<AllState>>, Json(q): Json<Query>)
              -> impl IntoResponse
{
    Json(all_state.query_exec.exec(&q).await)
}

async fn query_get(State(all_state): State<Arc<AllState>>, Path(q): Path<String>)
              -> impl IntoResponse
{
    match serde_json::from_str::<Query>(&q) {
        Ok(q) => Ok(Json(all_state.query_exec.exec(&q).await)),
        Err(err) => Err((StatusCode::BAD_REQUEST, err.to_string()))
    }
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

async fn term_complete(Path((cv_name, q)): Path<(String, String)>,
                       State(all_state): State<Arc<AllState>>)
              -> Json<TermCompletionResponse>
{
    let res = all_state.search.term_complete(&cv_name, &q).await;

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

    Json(completion_response)
}

async fn ref_complete(Path(q): Path<String>, State(all_state): State<Arc<AllState>>)
                -> Json<RefCompletionResponse>
{
    let res = all_state.search.ref_complete(&q).await;

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

    Json(completion_response)
}

async fn allele_complete(Path(q): Path<String>, State(all_state): State<Arc<AllState>>)
                -> Json<AlleleCompletionResponse>
{
    let res = all_state.search.allele_complete(&q).await;

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

    Json(completion_response)
}

// search for terms, refs or docs that match the query
async fn solr_search(Path((scope, q)): Path<(String, String)>, State(all_state): State<Arc<AllState>>)
    -> Json<SolrSearchResponse>
{
    if let Some(parsed_scope) = SolrSearchScope::new_from_str(&scope) {
        let search_result = all_state.search.solr_search(&parsed_scope, &q).await;

        match search_result {
            Ok(search_all_result) => {
                Json(SolrSearchResponse {
                    status: "Ok".to_owned(),
                    term_matches: search_all_result.term_matches,
                    ref_matches: search_all_result.ref_matches,
                    doc_matches: search_all_result.doc_matches,
                })
            },
            Err(err) => {
                Json(SolrSearchResponse {
                    status: err.to_string(),
                    term_matches: vec![],
                    ref_matches: vec![],
                    doc_matches: vec![],
                })
            },
        }
    } else {
        Json(SolrSearchResponse {
            status: format!("no such search scope: {}", scope),
            term_matches: vec![],
            ref_matches: vec![],
            doc_matches: vec![],
        })
    }
}

async fn motif_search(Path((scope, q)): Path<(String, String)>, State(all_state): State<Arc<AllState>>)
                -> Result<(HeaderMap, String), StatusCode>
{
    let res = all_state.search.motif_search(&scope, &q).await;

    match res {
        Ok(search_result) => {
            let mut headers = HeaderMap::new();
            headers.insert(header::CONTENT_TYPE, "application/json".parse().unwrap());
            Ok((headers, search_result))

        },
        Err(err) => {
            println!("Motif search error: {:?}", err);
            Err(StatusCode::NOT_FOUND)
        },
    }
}


async fn gene_ex_violin_plot(Path((plot_size, genes)): Path<(String, String)>,
                             State(all_state): State<Arc<AllState>>)
             -> Result<(HeaderMap, Full<Bytes>), StatusCode>
{
    let res = all_state.search.gene_ex_violin_plot(&plot_size, &genes).await;

    match res {
        Ok(svg_plot) => {
            let mut headers = HeaderMap::new();
            headers.insert(header::CONTENT_TYPE, "image/svg+xml".parse().unwrap());
            Ok((headers, Full::new(svg_plot.bytes)))
        },
        Err(err) => {
            println!("Motif search error: {:?}", err);
            Err(StatusCode::INTERNAL_SERVER_ERROR)
        },
    }
}

async fn ping() -> String {
    String::from("OK") + " " + PKG_NAME + " " + VERSION
}

async fn not_found() -> Json<Value> {
    json!({
        "status": "error",
        "reason": "Resource was not found."
    }).into()
}

fn print_usage(program: &str, opts: Options) {
    let brief = format!("Usage: {} [options]", program);
    print!("{}", opts.usage(&brief));
}

#[tokio::main]
async fn main() {
    println!("{} v{}", PKG_NAME, VERSION);

    let args: Vec<String> = env::args().collect();
    let mut opts = Options::new();

    opts.optflag("h", "help", "print this help message");
    opts.optopt("c", "config-file", "Configuration file name", "CONFIG");
    opts.optopt("b", "bind-address-and-port", "The address:port to bind to", "BIND_ADDRESS_AND_PORT");
    opts.optopt("m", "search-maps", "Search data", "MAPS_JSON_FILE");
    opts.optopt("d", "api-maps-database", "SQLite3 database of API maps", "API_MAPS_DATABASE");
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
    if !matches.opt_present("api-maps-database") {
        println!("no --api-maps-database|-d option");
        print_usage(&program, opts);
        process::exit(1);
    }
    if !matches.opt_present("web-root-dir") {
        println!("no --web-root-dir|-w option");
        print_usage(&program, opts);
        process::exit(1);
    }

    let bind_address_and_port = matches.opt_str("bind-address-and-port");
    let socket_addr =
        if let Some(bind_address_and_port) = bind_address_and_port {
            match bind_address_and_port.parse() {
                Ok(sock) => sock,
                Err(err) => panic!("{}", err)
            }
        } else {
           "0.0.0.0:8500".parse().unwrap()
        };

    let site_db_conn_string = matches.opt_str("s");
    let api_maps_database_path = matches.opt_str("d").unwrap();

    let search_maps_filename = matches.opt_str("m").unwrap();
    println!("Reading data files ...");

    let site_db =
        if let Some(conn_str) = site_db_conn_string {
            match SiteDB::new(&conn_str).await {
                Ok(site_db) => Some(site_db),
                Err(err) => panic!("failed to connect to site db: {}", err),
            }
        } else {
            None
        };

    let config_file_name = matches.opt_str("c").unwrap();
    let config = Config::read(&config_file_name);
    let api_maps = api_maps_from_file(&search_maps_filename);
    let api_maps_database_conn = Connection::open(api_maps_database_path).unwrap();
    let api_data = APIData::new(&config, api_maps_database_conn, api_maps);

    let query_exec = QueryExec::new(api_data, site_db);
    let search = Search::new(&config);

    let web_root_dir = matches.opt_str("w").unwrap();
    let static_file_state = StaticFileState {
        web_root_dir: web_root_dir.clone(),
    };

    let all_state = AllState {
        query_exec,
        search,
        static_file_state,
        config,
    };

    println!("Starting server ...");
    let app = Router::new()
        .route("/*path", get(get_misc))
        .route("/", get(get_index))
        .route("/structure_view/:structure_type/:id", get(structure_view))
        .route("/protein_feature_view/:full_or_widget/:gene_uniquename", get(protein_feature_view))
        .route("/simple/gene/:id", get(get_simple_gene))
        .route("/simple/reference/:id", get(get_simple_reference))
        .route("/simple/term/:id", get(get_simple_term))
        .route("/api/v1/dataset/latest/complete/allele/*q", get(allele_complete))
        .route("/api/v1/dataset/latest/complete/ref/:q", get(ref_complete))
        .route("/api/v1/dataset/latest/complete/term/:cv_name/:q", get(term_complete))
        .route("/api/v1/dataset/latest/data/allele/:id", get(get_allele))
        .route("/api/v1/dataset/latest/data/gene/:id", get(get_gene))
        .route("/api/v1/dataset/latest/data/genotype/:id", get(get_genotype))
        .route("/api/v1/dataset/latest/data/reference/:id", get(get_reference))
        .route("/api/v1/dataset/latest/data/seq_feature_page_features", get(seq_feature_page_features))
        .route("/api/v1/dataset/latest/data/term/:id", get(get_term))
        .route("/api/v1/dataset/latest/gene_ex_violin_plot/:plot_size/:genes", get(gene_ex_violin_plot))
        .route("/api/v1/dataset/latest/motif_search/:scope/:q", get(motif_search))
        .route("/api/v1/dataset/latest/protein_features/:full_or_widget/:gene_uniquename", get(get_protein_features))
        .route("/api/v1/dataset/latest/query/:q", get(query_get))
        .route("/api/v1/dataset/latest/query", post(query_post))
        .route("/api/v1/dataset/latest/search/:scope/:q", get(solr_search))
        .route("/api/v1/dataset/latest/summary/term/:id", get(get_term_summary_by_id))
        .route("/ping", get(ping))
        .fallback(not_found)
        .with_state(Arc::new(all_state));

    let app = NormalizePathLayer::trim_trailing_slash().layer(app).into_make_service();

    axum::Server::bind(&socket_addr)
        .serve(app)
        .await
        .unwrap();
}

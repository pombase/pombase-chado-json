extern crate getopts;

use axum::{
    body::Body, extract::{Path, Request, State}, http::{header, HeaderMap, StatusCode}, response::{Html, IntoResponse, Response}, routing::{get, post}, Json, Router, ServiceExt
};

use axum_extra::{headers::Range, TypedHeader};
use axum_range::{KnownSize, Ranged};

use pombase_gocam::RemoveType;
use tracing_subscriber::EnvFilter;
use tokio::fs::{read, File};

use tower::{layer::Layer, timeout::TimeoutLayer};

use tower_http::normalize_path::NormalizePathLayer;
use tower_http::trace::TraceLayer;

use serde::Serialize;
use serde_json::{json, Value};

use rand::rng;
use rand::seq::IteratorRandom;

use pombase_gocam_process::{model_connections_to_cytoscope, model_to_cytoscape_simple, GoCamCytoscapeStyle};
use pombase::{bio::gocam_model_process::{read_connected_gocam_models, read_gocam_model,
                                         read_merged_gocam_model},
              data_types::{GoCamId, GoCamSummary, ProteinViewType}, web::config::{PanelConfig, Testimonial}};

use rusqlite::Connection;

extern crate serde_derive;

extern crate pombase;

use std::{collections::{HashMap, HashSet, BTreeSet}, process, sync::Arc, time::Duration};

use getopts::Options;


use pombase::api::query::Query;

use pombase::api::search::{Search, DocSearchMatch, SolrSearchScope};
use pombase::api::query_exec::QueryExec;
use pombase::api_data::{api_maps_from_file, APIData};
use pombase::api::site_db::SiteDB;
use pombase::api::stats_plot::StatsPlots;

use pombase::data_types::{SolrTermSummary, SolrReferenceSummary, SolrAlleleSummary};
use pombase::web::simple_pages::{render_simple_gene_page, render_simple_reference_page,
                                 render_simple_term_page, render_simple_genotype_page};
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


async fn get_file_range(range: Range, file_path: &str)
      -> Ranged<KnownSize<File>>
{
    let file = File::open(file_path).await.unwrap();
    let body = KnownSize::file(file).await.unwrap();

    Ranged::new(Some(range), body)
}

async fn get_static_file(path: &str) -> Response {
    let res = read(path).await;

    let content_type = mime_guess::from_path(path).first_raw().unwrap_or("text/plain");

    match res {
        Ok(bytes) => {
            (StatusCode::OK, [(header::CONTENT_TYPE, content_type.to_string())], bytes).into_response()
        },
        Err(_) => {
            (StatusCode::NOT_FOUND, [(header::CONTENT_TYPE, "text/plain".to_string())], "not found".to_string()).into_response()
        }
    }
}

struct AllState {
    query_exec: QueryExec,
    gocam_data: HashMap<GoCamId, GoCamSummary>,
    search: Search,
    stats_plots: StatsPlots,
    static_file_state: StaticFileState,

    // gene ID to name map for cytoscape JSON generation
    gene_name_map: HashMap<String, String>,

    // needed because this information isn't in the GO-CAM model JSON
    pro_term_to_gene_map: HashMap<String, String>,

    config: Config,

    front_page_testimonals: Vec<Testimonial>,
    full_testimonals: Vec<Testimonial>,

    front_page_spotlights: Vec<PanelConfig>,
    full_spotlights: Vec<PanelConfig>,

    front_page_explore: Vec<PanelConfig>,
    full_explore: Vec<PanelConfig>,
}

// If the path is a directory, return path+"/index.html".  Otherwise
// try the path, then try path + ".json", then default to loading the
// Angular app from /index.html
async fn get_misc(range: Option<TypedHeader<Range>>,
                  Path(mut path): Path<String>,
                  State(all_state): State<Arc<AllState>>)
            -> Response
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
        if let Some(TypedHeader(maybe_range_header)) = range {
           return get_file_range(maybe_range_header, &full_path).await.into_response();
        } else {
           return get_static_file(&full_path).await;
        }
    }

    let json_path = format!("{}.json", full_path);

    if std::path::Path::new(&json_path).exists() {
        return get_static_file(&json_path).await;
    }

    // special case for missing JBrowse files - return 404
    if is_jbrowse_path {
        return (StatusCode::NOT_FOUND, [(header::CONTENT_TYPE, "text/plain".to_string())],
                format!("File not found: {}", full_path)).into_response()
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

async fn get_protein_features(Path((scope, gene_uniquename)): Path<(String, String)>,
                              State(all_state): State<Arc<AllState>>)
       -> impl IntoResponse
{
    let scope =
        match ProteinViewType::try_from(scope.as_ref()) {
            Ok(val) => val,
            Err(err) => {
                eprintln!("get_protein_features(): {}", err);
                return Err((StatusCode::INTERNAL_SERVER_ERROR, format!("Internal error: {}", err)));
            }
        };

    let res = all_state.query_exec.get_api_data().get_protein_features_of_gene(scope, &gene_uniquename)
              .map(|s| s.to_owned())
              .map(Json);
    option_json_to_result(&gene_uniquename, res)
}

async fn get_gocam_data(Path((_full_or_widget, gene_uniquename)): Path<(String, String)>,
                        State(all_state): State<Arc<AllState>>)
        -> impl IntoResponse
{
   let res = all_state.query_exec.get_api_data().get_gocam_data_of_gene(&gene_uniquename)
              .map(|s| s.to_owned())
              .map(Json);
    option_json_to_result(&gene_uniquename, res)
}

async fn get_all_gocam_data(State(all_state): State<Arc<AllState>>)
        -> impl IntoResponse
{
    let res = all_state.query_exec.get_api_data().get_all_gocam_data()
        .values().cloned().collect();

    type DataResult = Result<(StatusCode, Json<Vec<GoCamSummary>>), (StatusCode, String)>;

    let res: DataResult = Ok((StatusCode::OK, Json(res)));

    res
}

async fn get_gocam_data_by_id(Path(gocam_ids): Path<String>,
                                  State(all_state): State<Arc<AllState>>)
        -> impl IntoResponse
{
    let details_list: Vec<_> =
        gocam_ids.split(",").map(|gocam_id| {
            all_state.query_exec.get_api_data().get_gocam_details_by_id(gocam_id)
        }).collect();

    if details_list.is_empty() {
        Err((StatusCode::NOT_FOUND, format!("no page for: {}", gocam_ids)))
    } else {
        Ok((StatusCode::OK, Json(details_list)))
    }
}

async fn get_gocam_overlaps(State(all_state): State<Arc<AllState>>)
     -> impl IntoResponse
{
    let overlaps = &all_state.query_exec.get_api_data().get_maps().gocam_overlaps;

    Json(overlaps.to_owned())
}

async fn get_gocam_holes(State(all_state): State<Arc<AllState>>)
     -> impl IntoResponse
{
    let holes = &all_state.query_exec.get_api_data().get_maps().gocam_holes;

    Json(holes.to_owned())
}

async fn get_cytoscape_gocam_by_id(Path(gocam_id_arg): Path<String>,
                                   State(all_state): State<Arc<AllState>>)
       -> impl IntoResponse
{
    get_cytoscape_gocam_by_id_retain_genes(Path((gocam_id_arg, "".to_owned())), State(all_state)).await
}

async fn get_cytoscape_gocam_by_id_retain_genes(Path((gocam_id_arg, gene_list)): Path<(String, String)>,
                                   State(all_state): State<Arc<AllState>>)
       -> impl IntoResponse
{
    let static_file_state = &all_state.static_file_state;
    let web_root_dir = &static_file_state.web_root_dir;

    let all_gocam_data = &all_state.gocam_data;
    let overlaps = &all_state.query_exec.get_api_data().get_maps().gocam_overlaps;

    let id_split: Vec<_> = gocam_id_arg.split(':').collect();
    let (gocam_id, flag_string) = (id_split[0], id_split.get(1));
    let flags: HashSet<_> =
        if let Some(flag_string) = flag_string {
            flag_string.split(",").map(String::from).collect()
        } else {
            HashSet::new()
        };
    let gene_set: BTreeSet<_> =
        gene_list.split(",").map(|g| g.to_owned()).collect();

    let mut read_res =
        if gocam_id.contains("+") {
            let filtered_data = gocam_id.split("+")
               .filter_map(|gocam_id| {
                    all_gocam_data.get(gocam_id).map(|model| (gocam_id.into(), model.clone()))
                }).collect();
            read_merged_gocam_model(web_root_dir, &filtered_data,
                                    &flags, &gene_set).await
        } else {
            match gocam_id {
                "ALL_MERGED" => {
                    read_merged_gocam_model(web_root_dir, all_gocam_data,
                                            &flags, &gene_set).await
                }
                "ALL_CONNECTED" => {
                    read_connected_gocam_models(web_root_dir, all_gocam_data,
                                                overlaps,
                                                &flags).await
                }
                _ => read_gocam_model(web_root_dir, gocam_id,
                                      &flags).await
            }
        };

    let style = if flags.contains("hide_models") {
        GoCamCytoscapeStyle::HideParents
    } else {
        GoCamCytoscapeStyle::IncludeParents
    };

    if flags.contains("retain_genes") {
        // we're hightlighting genes in the mega-models then hiding
        // models that don't have a hightlighted gene
        read_res = read_res.map(|model| {
            let mut remove_types = HashSet::new();
            remove_types.insert(RemoveType::Chemicals);
            remove_types.insert(RemoveType::Targets);
            model.remove_nodes(remove_types)
                .retain_enabling_genes(&gene_set)
        });
    }


    match read_res {
        anyhow::Result::Ok(mut model) => {
            let gene_name_map = &all_state.gene_name_map;
            model.add_gene_name_map(gene_name_map);
            let pro_term_to_gene_map = &all_state.pro_term_to_gene_map;
            model.add_pro_term_to_gene_map(pro_term_to_gene_map);
            let elements =  model_to_cytoscape_simple(&model, overlaps, style);

         (StatusCode::OK, Json(elements)).into_response()
        },
        anyhow::Result::Err(err) => {
            eprintln!("err: {}", err);
           (StatusCode::NOT_FOUND, [(header::CONTENT_TYPE, "text/plain".to_string())], "not found".to_string()).into_response()
        }
    }
}

async fn get_model_summary_for_cytoscape_all_no_flags(State(all_state): State<Arc<AllState>>)
       -> impl IntoResponse
{
    get_model_summary_for_cytoscape_all(Path(String::default()), State(all_state)).await
}

async fn get_model_summary_for_cytoscape_all(Path(flags): Path<String>,
                                             State(all_state): State<Arc<AllState>>)
       -> impl IntoResponse
{
    let api_maps = all_state.query_exec.get_api_data().get_maps();
    let overlaps =
       if flags.contains("merge_by_chemical") {
         &api_maps.gocam_overlaps_merge_by_chemical
       } else {
         &api_maps.gocam_overlaps
       };

    let all_gocam_data = &all_state.gocam_data;

    let ids_and_titles: Vec<(String, String)> =
        all_gocam_data.iter()
        .map(|(gocam_id, summary)| {
            (format!("gomodel:{}", gocam_id), summary.title.to_string())
        })
        .collect();

    let model_connections = model_connections_to_cytoscope(overlaps, &ids_and_titles);

    Json(model_connections.to_owned())
}
async fn get_model_summary_for_cytoscape_connected_no_flags(State(all_state): State<Arc<AllState>>)
       -> impl IntoResponse
{
    get_model_summary_for_cytoscape_connected(Path(String::default()), State(all_state)).await
}

async fn get_model_summary_for_cytoscape_connected(Path(flags): Path<String>,
                                                   State(all_state): State<Arc<AllState>>)
       -> impl IntoResponse
{
    let api_maps = all_state.query_exec.get_api_data().get_maps();
    let overlaps =
       if flags.contains("merge_by_chemical") {
         &api_maps.gocam_overlaps_merge_by_chemical
       } else {
         &api_maps.gocam_overlaps
       };

    let model_connections = model_connections_to_cytoscope(overlaps, &[]);

    Json(model_connections.to_owned())
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

async fn get_index(State(all_state): State<Arc<AllState>>) -> Response {
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

async fn rna_2d_structure(Path((gene_uniquename, urs_id)): Path<(String, String)>,
                        State(all_state): State<Arc<AllState>>)
                        -> (StatusCode, Html<String>)
{
    let search_url = all_state.config.server.django_url.to_owned() + "/rna_2d_structure/";
    let params = [("gene_uniquename", gene_uniquename), ("urs_id", urs_id)];
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

async fn protein_feature_view(Path((scope, gene_uniquename)): Path<(String, String)>,
                              State(all_state): State<Arc<AllState>>)
   -> (StatusCode, Html<String>)
{
    let search_url = all_state.config.server.django_url.to_owned() + "/protein_feature_view/";
    let params = [("scope", scope),
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

async fn gocam_viz_view(Path((viz_or_view, full_or_widget, gocam_id)):
                                       Path<(String, String, String)>,
                                  State(all_state): State<Arc<AllState>>)
   -> (StatusCode, Html<String>)
{
    gocam_viz_view_highlight(Path((viz_or_view, full_or_widget, gocam_id, "".to_owned())),
                             State(all_state)).await
}

async fn gocam_viz_view_highlight(Path((viz_or_view, full_or_widget, gocam_id, highlight_gene_ids)):
                                       Path<(String, String, String, String)>,
                                  State(all_state): State<Arc<AllState>>)
   -> (StatusCode, Html<String>)
{
    let search_url = format!("{}/gocam_{}/", all_state.config.server.django_url, viz_or_view);

    let params = [("full_or_widget", full_or_widget),
                  ("gocam_id", gocam_id),
                  ("model_view_path", "/gocam/pombase-view/docs".to_owned()),
                  ("api_path", "/api/v1/dataset/latest/data/go-cam-cytoscape".to_owned()),
                  ("highlight_gene_ids", highlight_gene_ids)];
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

async fn rhea_widget(Path((view_type, rhea_id)): Path<(String, String)>,
                     State(all_state): State<Arc<AllState>>)
   -> (StatusCode, Html<String>)
{
    let search_url = all_state.config.server.django_url.to_owned() + "/rhea_widget/";
    let params = [("view_type", view_type),
                  ("rhea_id", rhea_id)];
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

async fn get_gocam_summary(all_models: bool, flags_string: &str, all_state: Arc<AllState>)
   -> (StatusCode, Html<String>)
{
    let url = format!("{}/gocam_connections/", all_state.config.server.django_url);

    let summary_type =
        if all_models {
            "all_models"
        } else {
            "connected_only"
        };

    let api_path = "/api/v1/dataset/latest/data/gocam/model_summary";

    let params = [("summary_type", summary_type),
                  ("model_view_path", "/gocam/pombase-view/docs"),
                  ("flags", flags_string),
                  ("api_path", api_path)];

    let client = reqwest::Client::new();
    let result = client.get(url).query(&params).send().await;

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

async fn get_gocam_summary_connected_no_flags(State(all_state): State<Arc<AllState>>)
   -> (StatusCode, Html<String>)
{
    get_gocam_summary_connected(Path(String::default()), State(all_state)).await
}

async fn get_gocam_summary_connected(Path(flags): Path<String>,
                                     State(all_state): State<Arc<AllState>>)
   -> (StatusCode, Html<String>)
{
    get_gocam_summary(false, &flags, all_state).await
}

async fn get_gocam_summary_all_no_flags(State(all_state): State<Arc<AllState>>)
   -> (StatusCode, Html<String>)
{
    get_gocam_summary_all(Path(String::default()), State(all_state)).await
}

async fn get_gocam_summary_all(Path(flags): Path<String>,
                               State(all_state): State<Arc<AllState>>)
   -> (StatusCode, Html<String>)
{
    get_gocam_summary(true, &flags, all_state).await
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
Return a simple HTML version a genotype page for search engines
*/
async fn get_simple_genotype(Path(id): Path<String>,
                             State(all_state): State<Arc<AllState>>) -> (StatusCode, Html<String>) {
    if let Some(genotype) = all_state.query_exec.get_api_data().get_genotype_details(&id) {
        (StatusCode::OK, Html(render_simple_genotype_page(&all_state.config, &genotype)))
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

async fn get_config_testimonials(Path((all_or_random, location)): Path<(String, String)>,
                                 State(all_state): State<Arc<AllState>>)
    -> Json<Vec<Testimonial>>
{
    let testimonials = if &location == "front" {
        &all_state.front_page_testimonals
    } else {
        &all_state.full_testimonals
    };

    if &all_or_random == "all" {
        Json(testimonials.clone())
    } else {
        let mut rng = rng();
        if let Some(el) = testimonials.iter().choose(&mut rng) {
            Json(vec![el.to_owned()])
        } else {
            Json(vec![])
        }
    }
}

async fn get_panel_configs(Path((panel_type, location)): Path<(String, String)>,
                           State(all_state): State<Arc<AllState>>)
    -> Json<Vec<PanelConfig>>
{
    if &location == "front" {
        let mut rng = rng();
        let panel_configs =
            match panel_type.as_str() {
                "spotlight" => &all_state.front_page_spotlights,
                "explore" => &all_state.front_page_explore,
                _ => return Json(vec![])
            };

        if let Some(el) = panel_configs.iter().choose(&mut rng) {
            Json(vec![el.to_owned()])
        } else {
            Json(vec![])
        }
    } else {
        match panel_type.as_str() {
            "spotlight" => Json(all_state.full_spotlights.clone()),
            "explore" => Json(all_state.full_explore.clone()),
            _ => return Json(vec![])
        }
    }
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

async fn motif_search(Path((scope, q, max_gene_details)): Path<(String, String, String)>, State(all_state): State<Arc<AllState>>)
                -> Result<(HeaderMap, String), StatusCode>
{
    let res = all_state.search.motif_search(&scope, &q, &max_gene_details).await;

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
             -> impl IntoResponse
{
    let res = all_state.stats_plots.gene_ex_violin_plot(&plot_size, &genes).await;

    match res {
        Ok(svg_plot) => {
            (StatusCode::OK,
            [(header::CONTENT_TYPE, "image/svg+xml")],
            Body::from(svg_plot.bytes)).into_response()
        },
        Err(err) => {
            (StatusCode::INTERNAL_SERVER_ERROR,
             [(header::CONTENT_TYPE, "text/plain")],
             println!("Gene expression plot error: {:?}", err)).into_response()
        },
    }
}

async fn get_stats(Path(graph_type): Path<String>,
                   State(all_state): State<Arc<AllState>>)
          -> impl IntoResponse
{
    let res = match graph_type.as_ref() {
        "curated_by_year" |
        "curatable_by_year" |
        "cumulative_curated_by_year" |
        "ltp_genes_per_pub_per_year_range" |
        "ltp_annotations_per_pub_per_year_range" |
        "htp_annotations_per_pub_per_year_range" |
        "community_response_rates" |
        "cumulative_annotation_type_counts_by_year" |
        "cumulative_micropublications_by_year"
            => all_state.stats_plots.get_svg_graph(graph_type.as_ref()).await,
        _ => return StatusCode::NOT_FOUND.into_response(),
    };

    match res {
        Ok(svg_plot) => {
            (StatusCode::OK,
             [(header::CONTENT_TYPE, "image/svg+xml")],
              Body::from(svg_plot.bytes)).into_response()
        },
        Err(err) => {
            println!("Error getting stats graph: {:?}", err);
            StatusCode::INTERNAL_SERVER_ERROR.into_response()
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

    let args: Vec<String> = std::env::args().collect();
    let mut opts = Options::new();

    opts.optflag("h", "help", "print this help message");
    opts.optopt("c", "config-file", "Configuration file name", "CONFIG");
    opts.optopt("", "testimonials", "JSON file of testimonials", "TESTIMONIALS_FILE");
    opts.optopt("", "spotlight-config", "JSON file of spotlight item configuration", "SPOTLIGHT_CONFIG_FILE");
    opts.optopt("", "explore-config", "JSON file of explore item configuration", "EXPLORE_CONFIG_FILE");
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

    tracing_subscriber::fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env()
                .or_else(|_| EnvFilter::try_new("pombase-chado-json=warn,tower_http=warn"))
                .unwrap(),
        )
        .init();

    let bind_address_and_port = matches.opt_str("bind-address-and-port");
    let listener =
        if let Some(bind_address_and_port) = bind_address_and_port {
           tokio::net::TcpListener::bind(bind_address_and_port).await.unwrap()
        } else {
           tokio::net::TcpListener::bind("0.0.0.0:8500").await.unwrap()
        };

    let site_db_conn_string = matches.opt_str("s");
    let api_maps_database_path = matches.opt_str("d").unwrap();

    let search_maps_filename = matches.opt_str("m").unwrap();

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

    let gocam_data = api_data.get_all_gocam_data();

    let gene_name_map = api_data.get_maps().gene_summaries
        .iter()
        .filter_map(|(k, v)|
            v.name.as_ref().map(|name| (k.to_std_string(), name.to_std_string())))
        .collect();

    let pro_term_to_gene_map = api_data.get_maps().pro_term_to_gene_map.clone();

    let query_exec = QueryExec::new(api_data, site_db);
    let search = Search::new(&config);
    let stats_plots = StatsPlots::new(&config);

    let web_root_dir = matches.opt_str("w").unwrap();
    let static_file_state = StaticFileState {
        web_root_dir: web_root_dir.clone(),
    };

    let testimonals =
        if let Some(testimonials_filename) = matches.opt_str("testimonials") {
            Testimonial::read_testimonials(&testimonials_filename)
        } else {
            vec![]
        };

    let mut front_page_testimonals = vec![];
    let mut full_testimonals = vec![];

    for testimonal in testimonals.into_iter() {
        match testimonal.location.as_str() {
          "FRONT" => front_page_testimonals.push(testimonal),
          "FULL" => full_testimonals.push(testimonal),
          "BOTH" => {
            front_page_testimonals.push(testimonal.clone());
            full_testimonals.push(testimonal);
          },
          _ => panic!("unknown location in testimonial: {}", testimonal.location),
        }
    }

    let spotlights =
        if let Some(spotlights_filename) = matches.opt_str("spotlight-config") {
            PanelConfig::read_panel_configs(&spotlights_filename)
        } else {
            vec![]
        };

    let mut front_page_spotlights = vec![];
    let mut full_spotlights = vec![];

    for spotlight in spotlights.into_iter() {
        if spotlight.show_on_front_page {
            front_page_spotlights.push(spotlight.clone());
        }
        full_spotlights.push(spotlight)
    }

    let explore =
        if let Some(explore_filename) = matches.opt_str("explore-config") {
            PanelConfig::read_panel_configs(&explore_filename)
        } else {
            vec![]
        };

    let mut front_page_explore = vec![];
    let mut full_explore = vec![];

    for explore_conf in explore.into_iter() {
        if explore_conf.show_on_front_page {
            front_page_explore.push(explore_conf.clone());
        }
        full_explore.push(explore_conf)
    }

    let all_state = AllState {
        query_exec,
        gocam_data,
        search,
        stats_plots,
        static_file_state,
        gene_name_map,
        pro_term_to_gene_map,
        config,
        front_page_testimonals,
        full_testimonals,
        front_page_spotlights,
        full_spotlights,
        front_page_explore,
        full_explore,
    };

    println!("Starting server ...");
    let app = Router::new()
        .route("/{*path}", get(get_misc))
        .route("/", get(get_index))
        .route("/structure_view/{structure_type}/{id}", get(structure_view))
        .route("/rna_2d_structure/{gene_uniquename}/{urs_id}", get(rna_2d_structure))
        .route("/protein_feature_view/{scope}/{gene_uniquename}", get(protein_feature_view))
        .route("/gocam_{viz_or_view}/{full_or_widget}/{gocam_id}", get(gocam_viz_view))
        .route("/gocam_{viz_or_view}/{full_or_widget}/{gocam_id}/{highlight_gene_ids}", get(gocam_viz_view_highlight))
        .route("/rhea_widget/{view_type}/{rhea_id}", get(rhea_widget))
        .route("/gocam_summary/connected", get(get_gocam_summary_connected_no_flags))
        .route("/gocam_summary/connected:{flags}", get(get_gocam_summary_connected))
        .route("/gocam_summary/all", get(get_gocam_summary_all_no_flags))
        .route("/gocam_summary/all:{flags}", get(get_gocam_summary_all))
        .route("/simple/gene/{id}", get(get_simple_gene))
        .route("/simple/genotype/{id}", get(get_simple_genotype))
        .route("/simple/reference/{id}", get(get_simple_reference))
        .route("/simple/term/{id}", get(get_simple_term))
        .route("/api/v1/dataset/latest/complete/allele/{*q}", get(allele_complete))
        .route("/api/v1/dataset/latest/complete/ref/{q}", get(ref_complete))
        .route("/api/v1/dataset/latest/complete/term/{cv_name}/{q}", get(term_complete))
        .route("/api/v1/dataset/latest/data/allele/{id}", get(get_allele))
        .route("/api/v1/dataset/latest/data/gene/{id}", get(get_gene))
        .route("/api/v1/dataset/latest/data/genotype/{id}", get(get_genotype))
        .route("/api/v1/dataset/latest/data/reference/{id}", get(get_reference))
        .route("/api/v1/dataset/latest/data/seq_feature_page_features", get(seq_feature_page_features))
        .route("/api/v1/dataset/latest/data/term/{id}", get(get_term))
        .route("/api/v1/dataset/latest/data/gocam/{full_or_widget}/{gene_uniquename}", get(get_gocam_data))
        .route("/api/v1/dataset/latest/data/gocam/all", get(get_all_gocam_data))
        .route("/api/v1/dataset/latest/data/gocam/by_id/{gocam_id}", get(get_gocam_data_by_id))
        .route("/api/v1/dataset/latest/data/gocam/overlaps", get(get_gocam_overlaps))
        .route("/api/v1/dataset/latest/data/gocam/holes", get(get_gocam_holes))
        .route("/api/v1/dataset/latest/data/gocam/model_summary/all_models:{flags}", get(get_model_summary_for_cytoscape_all))
        .route("/api/v1/dataset/latest/data/gocam/model_summary/all_models", get(get_model_summary_for_cytoscape_all_no_flags))
        .route("/api/v1/dataset/latest/data/gocam/model_summary/connected_only:{flags}", get(get_model_summary_for_cytoscape_connected))
        .route("/api/v1/dataset/latest/data/gocam/model_summary/connected_only", get(get_model_summary_for_cytoscape_connected_no_flags))
        .route("/api/v1/dataset/latest/data/go-cam-cytoscape/{gocam_id}", get(get_cytoscape_gocam_by_id))
        .route("/api/v1/dataset/latest/data/go-cam-cytoscape/{gocam_id}/{retain_genes}", get(get_cytoscape_gocam_by_id_retain_genes))
        .route("/api/v1/dataset/latest/config/testimonials/{all_or_random}/{location}", get(get_config_testimonials))
        .route("/api/v1/dataset/latest/config/panels/{panel_type}/{location}", get(get_panel_configs))
        .route("/api/v1/dataset/latest/gene_ex_violin_plot/{plot_size}/{genes}", get(gene_ex_violin_plot))
        .route("/api/v1/dataset/latest/stats/{type}", get(get_stats))
        .route("/api/v1/dataset/latest/motif_search/{scope}/{q}/{max_gene_details}", get(motif_search))
        .route("/api/v1/dataset/latest/protein_features/{scope}/{gene_uniquename}", get(get_protein_features))
        .route("/api/v1/dataset/latest/query", post(query_post))
        .route("/api/v1/dataset/latest/query/{q}", get(query_get))
        .route("/api/v1/dataset/latest/search/{scope}/{q}", get(solr_search))
        .route("/api/v1/dataset/latest/summary/term/{id}", get(get_term_summary_by_id))
        .route("/ping", get(ping))
        .fallback(not_found)
        .with_state(Arc::new(all_state))
        .layer(TraceLayer::new_for_http());

    let app = NormalizePathLayer::trim_trailing_slash().layer(app);

    let _: () = axum::serve(listener, ServiceExt::<Request>::into_make_service(app))
    .await
    .unwrap();
    ()
    .layer(TimeoutLayer::new(Duration::from_secs(120)));

    tracing_subscriber::fmt()
        .with_max_level(tracing::Level::DEBUG)
        .init();
}

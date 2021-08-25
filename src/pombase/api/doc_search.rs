
use std::collections::HashMap;

use crate::web::config::ServerConfig;
use crate::api::search_types::*;
use crate::api::search_utils::{do_solr_request, clean_words};

#[derive(Serialize, Deserialize, Debug)]
pub struct DocSearchMatch {
    pub id: String,
    pub heading: String,
    pub hl: SolrMatchHighlight,
}

#[derive(Serialize, Deserialize, Debug)]
struct DocSearchRes {
    pub id: String,
    pub heading: String,
}

#[derive(Deserialize, Debug)]
struct SolrDocResponse {
    pub docs: Vec<DocSearchRes>,
}

#[derive(Deserialize, Debug)]
struct SolrDocSearchResponseContainer {
    pub highlighting: Option<HashMap<SolrMatchId, SolrMatchHighlight>>,
    pub response: SolrDocResponse,
}

fn make_docs_url(config: &ServerConfig, q: &str) -> Option<String> {
    if !q.is_empty() {
        let clean_words = clean_words(q);

        if clean_words.is_empty() {
            return None;
        }

        let clean_q = clean_words.join(" ");

        let prefix = format!("{}/docs/select?wt=json&q=", &config.solr_url);

        Some(format!("{}heading:({}*) OR content:({}*)^0.2", prefix,
                     clean_q, clean_q))
    } else {
        None
    }
}

pub fn search_docs(config: &ServerConfig, q: &str) -> Result<Vec<DocSearchMatch>, String> {
    if let Some(mut url) = make_docs_url(config, q) {
        url += "&hl=on&hl.fl=heading,content&fl=id,heading,authors_abbrev";
        let res = do_solr_request(&url)?;

        match serde_json::from_reader(res) {
            Ok(container) => {
                let response_container: SolrDocSearchResponseContainer = container;
                let mut hl_by_id =
                    response_container.highlighting.unwrap_or_else(HashMap::new);
                let matches: Vec<DocSearchMatch> = response_container.response.docs
                    .iter().map(|doc| DocSearchMatch {
                        id: String::from(&doc.id),
                        heading: String::from(&doc.heading),
                        hl: hl_by_id.remove(doc.id.as_str())
                            .unwrap_or_else(HashMap::new),
                    }).collect();
                Ok(matches)
            },
            Err(err) => {
                Err(format!("Error parsing response from Solr: {:?}", err))
            }
        }
    } else {
        Ok(vec![])
    }
}


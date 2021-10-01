use regex::Regex;

use crate::data_types::SolrTermSummary;

use std::collections::HashMap;

use crate::web::config::ServerConfig;
use crate::api::search_types::*;
use crate::api::search_utils::{do_solr_request, clean_words, get_query_part};

#[derive(Deserialize, Debug)]
struct SolrTermSearchResponse {
    pub docs: Vec<SolrTermSummary>,
}

#[derive(Deserialize, Debug)]
struct SolrTermSearchResponseContainer {
    pub highlighting: HashMap<SolrMatchId, SolrMatchHighlight>,
    pub response: SolrTermSearchResponse,
}

#[derive(Deserialize, Debug)]
struct SolrTermResponse {
    pub docs: Vec<SolrTermSummary>,
}

#[derive(Deserialize, Debug)]
struct SolrTermResponseContainer {
    pub response: SolrTermResponse,
}

pub fn search_terms(config: &ServerConfig, q: &str)
                    -> Result<Vec<SolrTermSummary>, String>
{
    let cv_name = &config.cv_name_for_terms_search;

    if let Some(url) = make_terms_url(config, cv_name, q) {
        let res = do_solr_request(&url)?;

        match serde_json::from_reader(res) {
            Ok(container) => {
                Ok(container_to_matches(container))
            },
            Err(err) => {
                Err(format!("Error parsing response from Solr: {:?}", err))
            }
        }
    } else {
        Ok(vec![])
    }
}

fn container_to_matches(container: SolrTermSearchResponseContainer) -> Vec<SolrTermSummary> {
    let mut response_container: SolrTermSearchResponseContainer = container;
    let mut hl_by_id = response_container.highlighting;
    let matches: Vec<SolrTermSummary> = response_container.response.docs
        .drain(0..).map(|mut doc: SolrTermSummary| {
            doc.highlighting = hl_by_id.remove(doc.id.as_str())
              .unwrap_or_else(HashMap::new);
            doc
        }).collect();
    matches
}

pub fn make_terms_url(config: &ServerConfig, cv_name: &str, q: &str) -> Option<String> {
    let mut terms_url =
        config.solr_url.to_owned() + "/terms/select?wt=json&hl=on&hl.fl=name,definition&q=";

    let termid_re_string = r"^(?P<prefix>[\w_]+):(?P<accession>\d+)$";
    let termid_re = Regex::new(termid_re_string).unwrap();

    let parent_re_string = r"^\[".to_owned() + termid_re_string + r"\]$";
    let parent_re = Regex::new(&parent_re_string).unwrap();

    let maybe_captures = parent_re.captures(cv_name);

    if let Some(captures) = maybe_captures {
        let prefix = captures.name("prefix").unwrap().as_str();
        let accession = captures.name("accession").unwrap().as_str();
        terms_url += &format!("(interesting_isa_parents:{}\\:{} OR id:{}\\:{})",
                              prefix, accession, prefix, accession);
    } else {
        terms_url += &format!("cv_name:+{}^=1", cv_name);
    }

    if let Some(captures) = termid_re.captures(q) {
        let prefix = captures.name("prefix").unwrap().as_str();
        let accession = captures.name("accession").unwrap().as_str();
        terms_url = format!(r"{} AND (id:{}\:{} OR secondary_identifiers:{}\:{})",
                            terms_url, prefix, accession, prefix, accession);
    } else {
        let clean_words = clean_words(q);

        if clean_words.is_empty() {
            return None;
        }

        terms_url += " AND annotation_count:[1 TO *] AND (name:(";

        let query_part = get_query_part(&clean_words);

        terms_url += &format!("{}) OR close_synonym_words:({})^{} OR distant_synonym_words:({})^{} OR definition:({})^{})",
                              query_part, query_part, config.close_synonym_boost,
                              query_part, config.distant_synonym_boost,
                              query_part, config.term_definition_boost);
    }
    Some(terms_url)
}

pub fn term_complete(config: &ServerConfig, cv_name: &str, q: &str)
                     -> Result<Vec<SolrTermSummary>, String>
{
    if let Some(terms_url) = make_terms_url(config, cv_name, q) {
        match reqwest::blocking::get(&terms_url) {
            Ok(res) => {
                if res.status().is_success() {
                    match serde_json::from_reader(res) {
                        Ok(container) => {
                            Ok(container_to_matches(container))
                        },
                        Err(err) => {
                            Err(format!("Error parsing response from Solr: {:?}", err))
                        }
                    }
                } else {
                    if let Some(reason) = res.status().canonical_reason() {
                        Err(format!("HTTP request to Solr failed: {} - {}", res.status(), reason))
                    } else {
                        Err(format!("HTTP request to Solr failed with status code: {}",
                                    res.status()))
                    }
                }
            },
            Err(err) => {
                Err(format!("Error from Reqwest: {:?}", err))
            }
        }
    } else {
        Ok(vec![])
    }
}

pub fn term_summary_by_id(config: &ServerConfig, termid: &str)
                          -> Result<Option<SolrTermSummary>, String>
{
    let term_url = format!("{}/terms/select?wt=json&q=id:{}",
                           config.solr_url.to_owned(), termid.replace(":", r"\:"));

    match reqwest::blocking::get(&term_url) {
        Ok(res) => {
            if res.status().is_success() {
                match serde_json::from_reader(res) {
                    Ok(container) => {
                        let solr_response_container: SolrTermResponseContainer = container;
                        let summaries = solr_response_container.response.docs;
                        if let Some(summary) = summaries.get(0) {
                            if summary.id == termid {
                                Ok(Some(summary.clone()))
                            } else {
                                Ok(None)
                            }
                        } else {
                            Ok(None)
                        }
                    },
                    Err(err) => {
                        Err(format!("Error parsing response from Solr: {:?}", err))
                    }
                }
            } else {
                if let Some(reason) = res.status().canonical_reason() {
                    Err(format!("HTTP request to Solr failed: {} - {}", res.status(), reason))
                } else {
                    Err(format!("HTTP request to Solr failed with status code: {}",
                                res.status()))
                }
            }
        },
        Err(err) => {
            Err(format!("Error from Reqwest: {:?}", err))
        }
    }
}

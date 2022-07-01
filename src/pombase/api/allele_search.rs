use regex::Regex;

use crate::data_types::SolrAlleleSummary;

use std::collections::HashMap;

use crate::web::config::ServerConfig;
use crate::api::search_types::*;
use crate::api::search_utils::do_solr_request;

#[derive(Deserialize, Debug)]
struct SolrAlleleResponse {
    pub docs: Vec<SolrAlleleSummary>,
}

#[derive(Deserialize, Debug)]
struct SolrAlleleResponseContainer {
    pub highlighting: HashMap<SolrMatchId, SolrMatchHighlight>,
    pub response: SolrAlleleResponse,
}

lazy_static! {
    static ref CLEAN_RE: Regex = Regex::new(r#"([\+\-\&\\!\(\)"\~\*\?:])"#).unwrap();
}

// remove non-alphanumeric chars and then split on spaces
pub fn clean_words(q: &str) -> Vec<String> {
    let substring = |s: &str, len: usize| s.chars().take(len).collect::<String>();
    let lower_q = substring(&q.to_lowercase(), 200);
    CLEAN_RE.replace_all(&lower_q, "\\$1")
       .split_whitespace().map(|s| s.to_owned()).collect()
}

fn make_allele_url(config: &ServerConfig, q: &str, query_field_names: &[&str])
                   -> Option<String>
{
    let joined_query_fields = query_field_names.join(",");
    let mut allele_url =
        format!("{}/alleles/select?wt=json&hl=on&hl.fl={}&q=",
                &config.solr_url,
                &joined_query_fields);

    let substring = |s: &str, len: usize| s.chars().take(len).collect::<String>();
    let lower_q = substring(&q.to_lowercase(), 200);

    let mut clean_words: Vec<String> = clean_words(&lower_q);

    if clean_words.is_empty() {
        return None
    }

    if clean_words.iter().any(|s| s == "Δ" || s == "δ") {
        clean_words.push("delta".to_owned());
    }
    if clean_words.iter().any(|s| s == "delete" || s =="delta") {
        clean_words.push("deletion".to_owned());
    }

    let url_parts =
        itertools::iproduct!(query_field_names, clean_words)
        .map(|(field_name, query_word)| {
            let weight = if *field_name == "name" || *field_name == "allele_type" {
                2.0
            } else {
                if *field_name == "synonyms" || *field_name == "gene_name" {
                    1.0
                } else {
                    0.2
                }
            };

            format!("{}:({})^{}", field_name, query_word, weight)
        })
        .collect::<Vec<String>>();

    allele_url += &url_parts.join(" OR ");

    Some(allele_url)
}

fn matches_from_container(container: SolrAlleleResponseContainer) -> Vec<SolrAlleleSummary> {
    let mut response_container: SolrAlleleResponseContainer = container;
    let mut hl_by_id = response_container.highlighting;
    let matches: Vec<SolrAlleleSummary> = response_container.response.docs
        .drain(0..)
        .map(|mut doc: SolrAlleleSummary| {
            doc.highlighting = hl_by_id.remove(doc.id.as_str())
                .unwrap_or_else(HashMap::new);
            doc
        }).collect();
    matches
}

pub fn search_alleles(config: &ServerConfig, q: &str) -> Result<Vec<SolrAlleleSummary>, String> {
    let query_field_names = ["name", "synonyms", "description", "allele_type",
                             "gene_name", "gene_uniquename"];
    let maybe_url = make_allele_url(config, q, &query_field_names);

    if let Some(url) = maybe_url {
        let res = do_solr_request(&url)?;

        match serde_json::from_reader(res) {
            Ok(container) => Ok(matches_from_container(container)),
            Err(err) => {
                Err(format!("Error parsing response from Solr: {:?}", err))
            }
        }
    } else {
        Ok(vec![])
    }
}

use bytes::Buf;
use regex::Regex;
use reqwest::Client;
use std::error::Error;

use crate::data_types::SolrAlleleSummary;

use std::collections::HashMap;

use crate::web::config::ServerConfig;
use crate::api::search_types::*;
use crate::api::search_utils::solr_request;

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
    static ref DELTA_RE: Regex = Regex::new(r#"[Δδ]"#).unwrap();
}

// tidy search string - don't split for name as there are allele synonyms that
// contain spaces
fn allele_words(q: &str) -> Vec<String> {
    let substring = |s: &str, len: usize| s.chars().take(len).collect::<String>();
    let lc_substr = substring(&q.to_lowercase(), 200);
    let q = DELTA_RE.replace_all(&lc_substr, "delta").replace("wild type", "wild_type");

//    q.split_whitespace().map(|s| s.to_owned()).collect()
    vec![q.to_owned()]
}

fn make_allele_url(q: &str) -> Option<String> {
    let substring = |s: &str, len: usize| s.chars().take(len).collect::<String>();
    let lower_q = substring(&q.to_lowercase(), 100);

    let mut words: Vec<String> = allele_words(&lower_q);

    if words.is_empty() {
        return None
    }

    if words.iter().any(|s| s == "delete" || s =="delta" || s == "deletion") {
        words.push("deletion".to_owned());
        words.push("delta".to_owned());
    }

    let mut url_parts = vec![];
    let mut clean_words = vec![];

    for word in &words {
        let clean_word = CLEAN_RE.replace_all(word, "\\$1");

        clean_words.push(clean_word.clone());

        url_parts.push(format!("name:({})^5.0", clean_word));
        url_parts.push(format!("synonyms:({})^4.0", clean_word));
        url_parts.push(format!("allele_type:({})^4.0", clean_word));
        url_parts.push(format!("gene_name:({})^3.0", clean_word));
        url_parts.push(format!("gene_uniquename:({})^3.0", clean_word));
    }

    for word in &clean_words {
        url_parts.push(format!("description:({})^1.0", word));
    }

    Some(url_parts.join(" OR "))
}

fn matches_from_container(container: SolrAlleleResponseContainer) -> Vec<SolrAlleleSummary> {
    let mut response_container: SolrAlleleResponseContainer = container;
    let mut hl_by_id = response_container.highlighting;
    let matches: Vec<SolrAlleleSummary> = response_container.response.docs
        .drain(0..)
        .map(|mut doc: SolrAlleleSummary| {
            doc.highlighting = hl_by_id.remove(doc.id.as_str()).unwrap_or_default();
            doc
        }).collect();
    matches
}

pub async fn search_alleles(config: &ServerConfig, reqwest_client: &Client, q: &str)
     -> Result<Vec<SolrAlleleSummary>, Box<dyn Error + Send + Sync>>
{
    let query_field_names = ["name", "synonyms", "description",
                             "gene_name", "gene_uniquename"];

    let maybe_url = make_allele_url(q);

    let base_url = format!("{}/alleles/select", &config.solr_url);

    if let Some(url) = maybe_url {
        let res = solr_request(reqwest_client, &base_url, &query_field_names, &url).await?;

        let bytes = res.bytes().await?;

        match serde_json::from_reader(bytes.reader()) {
            Ok(res) => {
                Ok(matches_from_container(res))
            },
            Err(err) => {
                let mess = format!("Error parsing response from Solr: {:?}", err);
                Err(mess.as_str().into())
            }
        }
    } else {
        Ok(vec![])
    }
}

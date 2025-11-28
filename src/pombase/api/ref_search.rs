use crate::data_types::SolrReferenceSummary;

use regex::Regex;

use std::collections::HashMap;

use once_cell::unsync::Lazy;

use anyhow::Result;

use crate::web::config::ServerConfig;
use crate::api::search_types::*;
use crate::api::search_utils::{do_solr_request, clean_words};

lazy_static! {
    static ref CLEAN_WORDS_RE: Regex = Regex::new(r"([\w\d\-]+)").unwrap();
}

#[derive(Deserialize, Debug)]
struct SolrRefResponse {
    pub docs: Vec<SolrReferenceSummary>,
}

#[derive(Deserialize, Debug)]
struct SolrRefResponseContainer {
    pub highlighting: HashMap<SolrMatchId, SolrMatchHighlight>,
    pub response: SolrRefResponse,
}

fn make_refs_url(config: &ServerConfig, q: &str, query_field_names: &[&str])
                 -> Option<String>
{
    let joined_query_fields = query_field_names.join(",");
    let mut refs_url =
        format!("{}/refs/select?wt=json&hl=on&hl.fl={},publication_year&q=",
                &config.solr_url,
                &joined_query_fields);

    let id_re_string = r"^(?:(?P<prefix>[\w_]+):\s*)?(?P<rest>\d\d\d\d\d+)$";
    let id_re = Regex::new(id_re_string).unwrap();

    if let Some(captures) = id_re.captures(q) {
        let prefix =
            if let Some(prefix_capture) = captures.name("prefix") {
                prefix_capture.as_str()
            } else {
                "PMID"
            };
        let accession = captures.name("rest").unwrap().as_str();
        refs_url += "id:";
        refs_url += prefix;
        refs_url += r"\:";
        refs_url += accession;
    } else {
        let substring = |s: &str, len: usize| s.chars().take(len).collect::<String>();
        let lower_q = substring(&q.to_lowercase(), 200);

        let gene_uniquename_re =
            Lazy::new(|| {
                let lower_gene_re = config.gene_uniquename_re.to_lowercase();
                Regex::new(&lower_gene_re).unwrap()
            });

        if gene_uniquename_re.is_match(&lower_q) {
            refs_url += & query_field_names
                .iter()
                .map(|field_name| {
                    format!("{}:(\"{}\")", field_name, lower_q)
                })
                .collect::<Vec<String>>()
                .join(" OR ");
        } else {
            let clean_words: Vec<String> = clean_words(&lower_q);

            if clean_words.is_empty() {
                return None
            }

            let clean_words_length = clean_words.len();

            let mut clean_words_for_url = String::new();

            for (i, word) in clean_words.iter().enumerate() {
                if i == clean_words_length - 1 {
                    clean_words_for_url += &format!("{} {}*", word, word);
                } else {
                    clean_words_for_url += &format!("{} ", word);
                }
            }

            let url_parts =
                query_field_names
                .iter()
                .map(|field_name| {
                    let weight = if *field_name == "title" || *field_name == "authors" {
                        2.0
                    } else {
                        0.5
                    };
                    format!("{}:({})^{}", field_name, clean_words_for_url, weight)
                })
                .collect::<Vec<String>>();

            refs_url += &url_parts.join(" OR ");

            for word in clean_words {
                if word.len() == 4 && (word.starts_with("19") || word.starts_with("20"))
                    && let Ok(num) = word.parse::<u32>() {
                        refs_url += &format!(" OR publication_year:{}^20", num);
                    }
            }
        }
    }

    Some(refs_url)
}

fn matches_from_container(container: SolrRefResponseContainer) -> Vec<SolrReferenceSummary> {
    let mut response_container: SolrRefResponseContainer = container;
    let mut hl_by_id = response_container.highlighting;
    let matches: Vec<SolrReferenceSummary> = response_container.response.docs
        .drain(0..)
        .map(|mut doc: SolrReferenceSummary| {
            doc.highlighting = hl_by_id.remove(doc.id.as_str()).unwrap_or_default();
            doc
        }).collect();
    matches
}

pub async fn search_refs(config: &ServerConfig, q: &str) -> Result<Vec<SolrReferenceSummary>> {
    let query_field_names =
      ["title", "citation", "authors", "pubmed_abstract", "authors_abbrev"];
    let maybe_url = make_refs_url(config, q, &query_field_names);
    if let Some(url) = maybe_url {
        let res = do_solr_request(&url).await?;

        let container = res.json::<SolrRefResponseContainer>().await?;
        Ok(matches_from_container(container))
    } else {
        Ok(vec![])
    }
}

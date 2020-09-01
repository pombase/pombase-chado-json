use crate::data_types::SolrReferenceSummary;

use regex::Regex;

use std::collections::HashMap;

use once_cell::unsync::Lazy;

use crate::web::config::ServerConfig;
use crate::api::search_types::*;
use crate::api::search_utils::do_solr_request;

lazy_static! {
    static ref CLEAN_WORDS_RE: Regex = Regex::new(r"([\w\d\-]+)").unwrap();
}

#[derive(Serialize, Deserialize, Debug)]
pub struct RefSearchMatch {
    pub id: String,
    pub authors_abbrev: Option<String>,
    pub title: Option<String>,
    pub citation: Option<String>,
    pub publication_year: Option<u32>,
    pub hl: SolrMatchHighlight,
}

#[derive(Deserialize, Debug)]
struct SolrRefResponse {
    pub docs: Vec<SolrReferenceSummary>,
}

#[derive(Deserialize, Debug)]
struct SolrRefResponseContainer {
    pub response: SolrRefResponse,
}

#[derive(Deserialize, Debug)]
struct SolrRefSearchResponse {
    pub docs: Vec<RefSearchRes>,
}

#[derive(Deserialize, Debug)]
struct SolrRefSearchResponseContainer {
    pub highlighting: HashMap<SolrMatchId, SolrMatchHighlight>,
    pub response: SolrRefSearchResponse,
}

#[derive(Serialize, Deserialize, Debug)]
struct RefSearchRes {
    pub id: String,
    pub authors_abbrev: Option<String>,
    pub citation: Option<String>,
    pub title: Option<String>,
    pub publication_year: Option<u32>,
}

pub fn ref_complete(config: &ServerConfig, q: &str)
                    -> Result<Vec<SolrReferenceSummary>, String>
{
    if let Some(refs_url) = make_refs_url(config, q, &["title", "citation", "authors"]) {
        let res = do_solr_request(&refs_url)?;

        match serde_json::from_reader(res) {
            Ok(container) => {
                let solr_response_container: SolrRefResponseContainer = container;
                Ok(solr_response_container.response.docs)
            },
            Err(err) => {
                Err(format!("Error parsing response from Solr: {:?}", err))
            }
        }
    } else {
        // query string is too short to do search
        Ok(vec![])
    }
}

fn make_refs_url(config: &ServerConfig, q: &str, query_field_names: &[&str])
                 -> Option<String>
{
    let mut refs_url =
        config.solr_url.to_owned() + "/refs/select?wt=json&q=";

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
            let clean_words: Vec<String> =
                CLEAN_WORDS_RE.captures_iter(&lower_q)
                .map(|cap| cap.get(1).unwrap().as_str().to_owned()).collect();

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
                    let weight = if *field_name == "title" {
                        2.0
                    } else {
                        0.5
                    };
                    format!("{}:({})^{}", field_name, clean_words_for_url, weight)
                })
                .collect::<Vec<String>>();

            refs_url += &url_parts.join(" OR ");

            for word in clean_words {
                if word.len() == 4 && (word.starts_with("19") || word.starts_with("20")) {
                    if let Ok(num) = word.parse::<u32>() {
                        refs_url += &format!(" OR publication_year:{}^20", num);
                    }
                }
            }
        }
    }

    print!("{:?}\n", refs_url);

    Some(refs_url)
}

pub fn search_refs(config: &ServerConfig, q: &str) -> Result<Vec<RefSearchMatch>, String> {
    let hl_field_names = ["title", "citation", "authors", "pubmed_abstract"];
    let maybe_url = make_refs_url(config, q, &hl_field_names);
    if let Some(mut url) = maybe_url {
        url += &format!("&hl=on&hl.fl={},publication_year&fl=id,authors_abbrev,title,publication_year",
                        hl_field_names.join(","));
        let res = do_solr_request(&url)?;

        match serde_json::from_reader(res) {
            Ok(container) => {
                let response_container: SolrRefSearchResponseContainer = container;
                let mut hl_by_id = response_container.highlighting;
                let str_from = |s| String::from(s);
                let matches: Vec<RefSearchMatch> = response_container.response.docs
                    .iter().map(|doc| RefSearchMatch {
                        id: String::from(&doc.id),
                        citation: doc.citation.as_ref().map(str_from),
                        authors_abbrev: doc.authors_abbrev.as_ref().map(str_from),
                        title: doc.title.as_ref().map(str_from),
                        publication_year: doc.publication_year,
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

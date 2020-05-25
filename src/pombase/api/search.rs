use regex::Regex;
use serde_json;

use std::collections::HashMap;

use reqwest::Response;

use crate::data_types::{SolrTermSummary, SolrReferenceSummary};
use crate::web::config::Config;

pub struct Search {
    solr_url: String,
    close_synonym_boost: f32,
    distant_synonym_boost: f32,
    django_url: String,
    cv_name_for_terms_search: String,
}

type SolrMatchId = String;
type SolrFieldName = String;
type SolrHighlightDetails = String;

type SolrMatchHighlight = HashMap<SolrFieldName, Vec<SolrHighlightDetails>>;

#[derive(Serialize, Deserialize, Debug)]
struct TermSearchRes {
    pub id: String,
    pub cv_name: String,
    pub name: String,
}

#[derive(Serialize, Deserialize, Debug)]
struct RefSearchRes {
    pub id: String,
    pub authors_abbrev: Option<String>,
    pub title: Option<String>,
    pub publication_year: Option<u32>,
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

#[derive(Deserialize, Debug)]
struct SolrTermSearchResponse {
    pub docs: Vec<TermSearchRes>,
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
pub struct TermSearchMatch {
    pub id: String,
    pub cv_name: String,
    pub name: String,
    pub hl: SolrMatchHighlight,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct RefSearchMatch {
    pub id: String,
    pub authors_abbrev: Option<String>,
    pub title: Option<String>,
    pub publication_year: Option<u32>,
    pub hl: SolrMatchHighlight,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct DocSearchMatch {
    pub id: String,
    pub heading: String,
    pub hl: SolrMatchHighlight,
}

#[derive(Serialize, Debug)]
pub struct SearchAllResult {
    pub term_matches: Vec<TermSearchMatch>,
    pub ref_matches: Vec<RefSearchMatch>,
    pub doc_matches: Vec<DocSearchMatch>,
}

impl Search {
    pub fn new(config: &Config) -> Search {
        Search {
            solr_url: config.server.solr_url.clone(),
            close_synonym_boost: config.server.close_synonym_boost,
            distant_synonym_boost: config.server.distant_synonym_boost,
            django_url: config.server.django_url.clone(),
            cv_name_for_terms_search: config.server.cv_name_for_terms_search.clone(),
        }
    }

    fn get_query_part(&self, words: &[String]) -> String {
        let mut ret = String::new();

        let words_length = words.len();

        for (i, word) in words.iter().enumerate() {
            if i == words_length - 1 {
                ret += &format!("{} {}*", word, word);
            } else {
                ret += &format!("{} ", word);
            }
        }

        ret
    }

    fn make_terms_url(&self, cv_name: &str, q: &str) -> Option<String> {
        let mut terms_url =
            self.solr_url.to_owned() + "/terms/select?wt=json&q=";

        let termid_re_string = r"^(?P<prefix>[\w_]+):(?P<accession>\d+)$";
        let termid_re = Regex::new(termid_re_string).unwrap();

        let parent_re_string = r"^\[".to_owned() + termid_re_string + r"\]$";
        let parent_re = Regex::new(&parent_re_string).unwrap();

        let maybe_captures = parent_re.captures(cv_name);

        if let Some(captures) = maybe_captures {
            let prefix = captures.name("prefix").unwrap().as_str();
            let accession = captures.name("accession").unwrap().as_str();
            terms_url += &format!("(interesting_parents:{}\\:{} OR id:{}\\:{})",
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
            let substring = |s: &str, len: usize| s.chars().take(len).collect::<String>();
            let lower_q = substring(&q.to_lowercase(), 200);

            terms_url += " AND (name:(";

            let clean_words: Vec<String> =
                Regex::new(r"([\w\d\-]+)").unwrap().captures_iter(&lower_q)
                .map(|cap| cap.get(1).unwrap().as_str().to_owned()).collect();

            if clean_words.is_empty() {
                return None;
            }

            let query_part = self.get_query_part(&clean_words);

            terms_url += &format!("{}) OR close_synonym_words:({})^{} OR distant_synonym_words:({})^{})",
                                  query_part, query_part, self.close_synonym_boost,
                                  query_part, self.distant_synonym_boost);
        }
        Some(terms_url)
    }

    pub fn term_complete(&self, cv_name: &str, q: &str)
                         -> Result<Vec<SolrTermSummary>, String>
    {
        if let Some(terms_url) = self.make_terms_url(cv_name, q) {
            match reqwest::get(&terms_url) {
                Ok(res) => {
                    if res.status().is_success() {
                        match serde_json::from_reader(res) {
                            Ok(container) => {
                                let solr_response_container: SolrTermResponseContainer = container;
                                Ok(solr_response_container.response.docs)
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

    pub fn ref_complete(&self, q: &str) -> Result<Vec<SolrReferenceSummary>, String> {
        if let Some(refs_url) = self.make_refs_url(q) {
            let res = self.do_solr_request(&refs_url)?;

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

    fn make_refs_url(&self, q: &str) -> Option<String> {
        let mut refs_url =
            self.solr_url.to_owned() + "/refs/select?wt=json&q=";

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

            let clean_words: Vec<String> =
                Regex::new(r"([\w\d\-]+)").unwrap().captures_iter(&lower_q)
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
                vec!["title", "citation", "authors"]
                .iter()
                .map(|field_name| format!("{}:({})", field_name, clean_words_for_url))
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

        Some(refs_url)
    }

    fn do_solr_request(&self, url: &str) -> Result<Response, String> {
        match reqwest::get(url) {
            Ok(res) => {
                if res.status().is_success() {
                    Ok(res)
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
                Err(format!("Error from Reqwest: {:?} for {}", err, url))
            }
        }
    }

    pub fn motif_search(&self, scope: &str, pattern: &str) -> Result<String, String> {
        let search_url = self.django_url.to_owned() + "/motifsearch/query/";
        let params = [("scope", scope), ("pattern", pattern)];
        let client = reqwest::Client::new();
        let request = client.get(&search_url).query(&params);
        let result = request.send();

        match result {
            Ok(mut res) => {
                match res.text() {
                    Ok(text) => Ok(text),
                    Err(err) => Err(format!("Error getting text from motif search: {:?}", err)),
                }
            },
            Err(err) => {
                Err(format!("Motif search error: {:?}", err))
            }
        }
    }

    fn search_terms(&self, q: &str) -> Result<Vec<TermSearchMatch>, String> {
        let cv_name = &self.cv_name_for_terms_search;

        if let Some(mut url) = self.make_terms_url(cv_name, q) {
            url += "&hl=on&hl.fl=name&fl=id,name,cv_name";
            let res = self.do_solr_request(&url)?;

            match serde_json::from_reader(res) {
                Ok(container) => {
                    let response_container: SolrTermSearchResponseContainer = container;
                    let mut hl_by_id = response_container.highlighting;
                    let matches: Vec<TermSearchMatch> = response_container.response.docs
                        .iter().map(|doc| TermSearchMatch {
                            id: String::from(&doc.id),
                            cv_name: String::from(doc.cv_name.as_str()),
                            name: String::from(doc.name.as_str()),
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

    fn search_refs(&self, q: &str) -> Result<Vec<RefSearchMatch>, String> {
        if let Some(mut url) = self.make_refs_url(q) {
            url += "&hl=on&hl.fl=title,authors,citation&fl=id,authors_abbrev,title,publication_year";
            let res = self.do_solr_request(&url)?;

            match serde_json::from_reader(res) {
                Ok(container) => {
                    let response_container: SolrRefSearchResponseContainer = container;
                    let mut hl_by_id = response_container.highlighting;
                    let str_from = |s| String::from(s);
                    let matches: Vec<RefSearchMatch> = response_container.response.docs
                        .iter().map(|doc| RefSearchMatch {
                            id: String::from(&doc.id),
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

    fn make_docs_url(&self, q: &str) -> Option<String> {
        if q.len() > 0 {
            let prefix =
                self.solr_url.to_owned() + "/docs/select?wt=json&q=";

            Some(format!("{}heading:({}*) OR content:({}*)^0.5", prefix, q, q))
        } else {
            None
        }
    }

    fn search_docs(&self, q: &str) -> Result<Vec<DocSearchMatch>, String> {
        if let Some(mut url) = self.make_docs_url(q) {
            url += "&hl=on&hl.fl=heading,content&fl=id,heading";
            let res = self.do_solr_request(&url)?;

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

    pub fn search_all(&self, q: &str) -> Result<SearchAllResult, String> {
        let trimmed_query = q.trim();

        let term_matches = self.search_terms(trimmed_query)?;
        let ref_matches = self.search_refs(trimmed_query)?;
        let doc_matches = self.search_docs(trimmed_query)?;

        Ok(SearchAllResult {
            term_matches,
            ref_matches,
            doc_matches,
        })
    }
}

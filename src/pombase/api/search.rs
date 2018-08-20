use regex::Regex;
use reqwest;
use serde_json;

use web::data::{SolrTermSummary, SolrReferenceSummary};
use web::config::Config;

pub struct Search {
    solr_url: String,
    close_synonym_boost: f32,
    distant_synonym_boost: f32,
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

impl Search {
    pub fn new(config: &Config) -> Search {
        Search {
            solr_url: config.server.solr_url.clone(),
            close_synonym_boost: config.server.close_synonym_boost,
            distant_synonym_boost: config.server.distant_synonym_boost,
        }
    }

    fn get_query_part(&self, words: &[String]) -> String {
        let mut ret = String::new();

        let words_length = words.len();

        for (i, word) in words.iter().enumerate() {
            if i == words_length - 1 {
                ret += &format!("{} {}~0.8 {}*", word, word, word);
            } else {
                ret += &format!("{} {}~0.8 ", word, word);
            }
        }

        ret
    }

    pub fn term_complete(&self, cv_name: &str, q: &str)
                         -> Result<Vec<SolrTermSummary>, String>
    {
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
            terms_url += "cv_name:";
            terms_url += cv_name;
        }

        if let Some(captures) = termid_re.captures(q) {
            let prefix = captures.name("prefix").unwrap().as_str();
            let accession = captures.name("accession").unwrap().as_str();
            terms_url += " AND id:";
            terms_url += prefix;
            terms_url += r"\:";
            terms_url += accession;
        } else {
            let substring = |s: &str, len: usize| s.chars().take(len).collect::<String>();
            let lower_q = substring(&q.to_lowercase(), 200);

            terms_url += " AND (name:(";

            let clean_words: Vec<String> =
                Regex::new(r"([\w\d\-]+)").unwrap().captures_iter(&lower_q)
                .map(|cap| cap.get(1).unwrap().as_str().to_owned()).collect();

            if clean_words.is_empty() {
                return Ok(vec![]);
            }

            let clean_words_length = clean_words.len();

            for (i, word) in clean_words.iter().enumerate() {
                if i == clean_words_length - 1 {
                    terms_url += &format!("{} {}*", word, word);
                } else {
                    terms_url += &format!("{} ", word);
                }
            }

            let query_part = self.get_query_part(&clean_words);

            terms_url += &format!(") OR close_synonym_words:({})^{} OR distant_synonym_words:({})^{})",
                                  query_part, self.close_synonym_boost,
                                  query_part, self.distant_synonym_boost);
        }
        print!("Solr URL: {:?}\n", terms_url);

        match reqwest::get(&terms_url) {
            Ok(res) => {
                println!("Status: {}", res.status());

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
    }

    pub fn ref_complete(&self, q: &str) -> Result<Vec<SolrReferenceSummary>, String> {
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
                return Ok(vec![]);
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
                vec!["title", "citation", "authors", "authors_abbrev"]
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

        print!("Solr refs URL: {:?}\n", refs_url);

        match reqwest::get(&refs_url) {
            Ok(res) => {
                println!("Status: {}", res.status());

                if res.status().is_success() {
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

}

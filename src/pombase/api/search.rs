use regex::Regex;
use reqwest;
use serde_json;

use web::data::SolrTermSummary;

pub struct Search {
    solr_url: String,
}

#[derive(Deserialize, Debug)]
struct SolrResponse {
    pub docs: Vec<SolrTermSummary>,
}

#[derive(Deserialize, Debug)]
struct SolrResponseContainer {
    pub response: SolrResponse,
}

impl Search {
    pub fn new(solr_url: String) -> Search {
        Search {
            solr_url: solr_url,
        }
    }

    fn get_query_part(&self, words: &Vec<String>) -> String {
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
                         -> Result<Vec<SolrTermSummary>, reqwest::Error>
    {
        let mut terms_url =
            self.solr_url.to_owned() + "/terms/select?wt=json&q=";

        let termid_re_string = r"(?P<prefix>[\w_]+):(?P<accession>\d+)";
        let termid_re = Regex::new(termid_re_string).unwrap();

        let parent_re_string = r"^\[".to_owned() + termid_re_string + r"\]$";
        let parent_re = Regex::new(&parent_re_string).unwrap();

        let maybe_captures = parent_re.captures(cv_name);

        if let Some(captures) = maybe_captures {
            let prefix = captures.name("prefix").unwrap();
            let accession = captures.name("accession").unwrap();
            terms_url += &format!("(interesting_parents:{}\\:{} OR id:{}\\:{})",
                                  prefix, accession, prefix, accession);
        } else {
            terms_url += "cv_name:";
            terms_url += cv_name;
        }

        if let Some(captures) = termid_re.captures(q) {
            let prefix = captures.name("prefix").unwrap();
            let accession = captures.name("accession").unwrap();
            terms_url += " AND id:";
            terms_url += prefix;
            terms_url += r"\:";
            terms_url += accession;
        } else {
            terms_url += " AND (name:(";

            let clean_words: Vec<String> =
                Regex::new(r"(\w+)").unwrap().captures_iter(q)
                .map(|cap| cap.at(1).unwrap().to_owned()).collect();

            if clean_words.len() == 0 {
                return Ok(vec![]);
            }

            let clean_words_length = clean_words.len();

            for (i, word) in clean_words.iter().enumerate() {
                if i == clean_words_length - 1 {
                    terms_url += &format!("{} {}~0.8 {}*", word, word, word);
                } else {
                    terms_url += &format!("{} {}~0.8 ", word, word);
                }
            }

            let query_part = self.get_query_part(&clean_words);

            terms_url += &format!(")^4 OR close_synonyms:({})^0.2 OR distant_synonyms:({})^0.1)",
                                 query_part, query_part);
        }
        print!("{:?}\n", terms_url);

        let res = reqwest::get(&terms_url)?;

        println!("Status: {}", res.status());
        println!("Headers:\n{}", res.headers());

        match serde_json::from_reader(res) {
            Ok(solr_response_container) => {
                let container: SolrResponseContainer = solr_response_container;
                Ok(container.response.docs)
            },
            Err(err) => {
                panic!(format!("{:?}", err));
            }
        }
    }
}

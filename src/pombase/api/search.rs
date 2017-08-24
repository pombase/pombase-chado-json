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

    fn get_query_part(&self, words: &Vec<String>, weight: f32) -> String {
        let mut ret = String::new();

        let words_length = words.len();

        for (i, word) in words.iter().enumerate() {
            if i == words_length - 1 {
                ret += &format!("{}^{} {}*", word, weight, word);
            } else {
                ret += &format!("{}^{} {}~0.8 ", word, weight, word);
            }
        }

        ret
    }

    pub fn term_complete(&self, cv_name: &str, q: &str)
                         -> Result<Vec<SolrTermSummary>, reqwest::Error>
    {
        let mut terms_url =
            self.solr_url.to_owned() + "/terms/select?wt=json&q=";

        if Regex::new(r"^\D+:\d+$").unwrap().is_match(cv_name) {
            terms_url += "interesting_parents:";
            terms_url += cv_name;
        } else {
            terms_url += "cv_name:";
            terms_url += cv_name;
        }

        terms_url += " AND name:(";

        let clean_words: Vec<String> =
            Regex::new(r"(\w+)").unwrap().captures_iter(q)
            .map(|cap| cap.at(1).unwrap().to_owned()).collect();

        if clean_words.len() == 0 {
            return Ok(vec![]);
        }

        let clean_words_length = clean_words.len();

        for (i, word) in clean_words.iter().enumerate() {
            if i == clean_words_length - 1 {
                terms_url += &format!("{}^4 {}*", word, word);
            } else {
                terms_url += &format!("{}^4 {}~0.8 ", word, word);
            }
        }

        terms_url += ") OR close_synonyms:(";
        terms_url += &self.get_query_part(&clean_words, 0.5);
        terms_url += ") OR distant_synonyms:(";
        terms_url += &self.get_query_part(&clean_words, 0.1);
        terms_url += ")";

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

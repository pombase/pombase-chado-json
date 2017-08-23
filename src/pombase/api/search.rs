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

    pub fn term_complete(&self, cv_name: &str, q: &str)
                         -> Result<Vec<SolrTermSummary>, reqwest::Error>
    {
        let mut terms_url =
            self.solr_url.to_owned() + "/terms/select?wt=json&q=cv_name:" +
            cv_name + " AND name:(";

        let clean_words: Vec<String> =
            Regex::new(r"(\w+)").unwrap().captures_iter(q)
            .map(|cap| cap.at(1).unwrap().to_owned()).collect();

        if clean_words.len() == 0 {
            return Ok(vec![]);
        }

        let clean_words_length = clean_words.len();

        for (i, word) in clean_words.into_iter().enumerate() {
            if i == clean_words_length - 1 {
                terms_url += &format!("{}^2 {}*", word, word);
            } else {
                terms_url += &format!("{}^2 {}~0.8 ", word, word);
            }
        }
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
                panic!(err);
            }
        }
    }
}

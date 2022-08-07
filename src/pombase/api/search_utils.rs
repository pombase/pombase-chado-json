use std::error::Error;

use reqwest::{Response, Client};
use regex::Regex;

pub fn get_query_part(words: &[String]) -> String {
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

lazy_static! {
    static ref WORD_RE: Regex = Regex::new(r"([\w\d\-]+)").unwrap();
    static ref SHORT_NUMBERS_RE: Regex = Regex::new(r"^\d\d?\d?$").unwrap();
}

// remove non-alphanumeric chars and then split on spaces
pub fn clean_words(q: &str) -> Vec<String> {
    let substring = |s: &str, len: usize| s.chars().take(len).collect::<String>();
    let lower_q = substring(&q.to_lowercase(), 200);

    WORD_RE.captures_iter(&lower_q)
        .map(|cap| cap.get(1).unwrap().as_str().to_owned())
        // fixes the "02" and "2" from gene IDs like "SPCC1620.02.2" matching
        // lots of titles and abstracts:
        .filter(|s| !SHORT_NUMBERS_RE.is_match(s))
        .collect()
}

pub fn do_solr_request(url: &str) -> Result<reqwest::blocking::Response, String> {
    println!("do_solr_request({:?})", url);
    match reqwest::blocking::get(url) {
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

pub async fn solr_request(reqwest_client: &Client, url_base: &str, highlight_fields: &[&str], query: &str)
   -> Result<Response, Box<dyn Error + Send + Sync>>
{
    eprintln!("querying with: {}", query);

    let highlighted_fields_str = highlight_fields.join(",");
    let params = vec![("wt", "json"), ("hl", "on"), ("hl.fl", &highlighted_fields_str), ("q", query)];

    let url = reqwest::Url::parse_with_params(url_base, &params)?;

    let req = reqwest_client
        .get(url)
        .header("Accepts", "application/json");

    match req.send().await {
        Ok(res) => {
            if res.status().is_success() {
                Ok(res)
            } else {
                if let Some(reason) = res.status().canonical_reason() {
                    let var_name = format!("HTTP request to Solr failed: {} - {}", res.status(), reason);
                    Err(var_name.as_str().into())
                } else {
                    let mess = format!("HTTP request to Solr failed with status code: {}", res.status());
                    Err(mess.as_str().into())
                 }
            }
        },
        Err(err) => {
            let mess = format!("Error from Reqwest: {:?} for {} q:{}", err, url_base, query);
            Err(mess.as_str().into())
        }
    }
}


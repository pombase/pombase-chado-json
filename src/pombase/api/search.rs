use bytes::Bytes;

use std::io::Cursor;

use rocket::response::{self, Response, Responder};
use rocket::request::Request;
use rocket::http::ContentType;

use crate::web::config::{Config, ServerConfig};

use crate::data_types::{SolrTermSummary, SolrReferenceSummary};

use crate::api::term_search::{search_terms, term_complete, term_summary_by_id};

use crate::api::ref_search::{search_refs};

use crate::api::doc_search::search_docs;
pub use crate::api::doc_search::DocSearchMatch;

pub struct Search {
    config: ServerConfig,
}

pub struct PNGPlot {
    pub bytes: Bytes
}

impl<'a> Responder<'a, 'a> for PNGPlot {
    fn respond_to(self, _: &Request) -> response::Result<'a> {
        Response::build()
            .sized_body(self.bytes.len(), Cursor::new(self.bytes))
            .header(ContentType::new("image", "png"))
            .ok()
    }
}


#[derive(Deserialize, Debug)]
pub enum SolrSearchScope {
#[serde(rename = "term")]
    Term,
#[serde(rename = "reference")]
    Reference,
#[serde(rename = "doc")]
    Documentation,
}

impl SolrSearchScope {
    pub fn new_from_str(scope_str: &str) -> Option<SolrSearchScope> {
        match scope_str {
            "term" => Some(SolrSearchScope::Term),
            "ref" => Some(SolrSearchScope::Reference),
            "doc" => Some(SolrSearchScope::Documentation),
            _ => None,
        }
    }

    pub fn is_term(&self) -> bool {
        matches!(self, SolrSearchScope::Term)
    }

    pub fn is_reference(&self) -> bool {
        matches!(self, SolrSearchScope::Reference)
    }

    pub fn is_documentation(&self) -> bool {
        matches!(self, SolrSearchScope::Documentation)
    }
}


#[derive(Serialize, Debug)]
pub struct SolrSearchResult {
    pub term_matches: Vec<SolrTermSummary>,
    pub ref_matches: Vec<SolrReferenceSummary>,
    pub doc_matches: Vec<DocSearchMatch>,
}

impl Search {
    pub fn new(config: &Config) -> Search {
        Search {
            config: config.server.clone(),
        }
    }

    pub fn motif_search(&self, scope: &str, pattern: &str) -> Result<String, String> {
        let search_url = self.config.django_url.to_owned() + "/motifsearch/query/";
        let params = [("scope", scope), ("pattern", pattern)];
        let client = reqwest::blocking::Client::new();
        let result = client.get(&search_url).query(&params).send();

        match result {
            Ok(res) => {
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

    pub fn gene_ex_violin_plot(&self, plot_size: &str, genes: &str)
                               -> Result<PNGPlot, String>
    {
        let plot_url = self.config.django_url.to_owned() + "/gene_ex/gene_ex_violin/";
        let params = [("plot_size", plot_size), ("genes", genes)];
        let client = reqwest::blocking::Client::new();
        let result = client.get(&plot_url).query(&params).send();

        match result {
            Ok(res) => {
                match res.bytes() {
                    Ok(bytes) => Ok(PNGPlot { bytes }),
                    Err(err) => Err(format!("Error getting violin plot image: {:?}", err)),
                }
            },
            Err(err) => {
                Err(format!("Gene expression violin plot error: {:?}", err))
            }
        }
    }

    pub fn term_complete(&self, cv_name: &str, q: &str)
                         -> Result<Vec<SolrTermSummary>, String>
    {
        term_complete(&self.config, cv_name, q)
    }

    pub fn term_summary_by_id(&self, termid: &str)
                             -> Result<Option<SolrTermSummary>, String>
    {
        term_summary_by_id(&self.config, termid)
    }

    pub fn ref_complete(&self, q: &str)
                        -> Result<Vec<SolrReferenceSummary>, String>
    {
        search_refs(&self.config, q)
    }

    pub fn solr_search(&self, scope: &SolrSearchScope, q: &str)
                       -> Result<SolrSearchResult, String>
    {
        let trimmed_query = q.trim();

        let term_matches =
            if scope.is_term() {
                search_terms(&self.config, trimmed_query)?
            } else {
                vec![]
            };
        let ref_matches =
            if scope.is_reference() {
                search_refs(&self.config, trimmed_query)?
            } else {
                vec![]
            };
        let doc_matches =
            if scope.is_documentation() {
                search_docs(&self.config, trimmed_query)?
            } else {
                vec![]
            };

        Ok(SolrSearchResult {
            term_matches,
            ref_matches,
            doc_matches,
        })
    }
}

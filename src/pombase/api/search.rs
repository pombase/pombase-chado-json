use bytes::Bytes;
use reqwest::Client;
use std::error::Error;

use anyhow::Result;

use crate::web::config::{Config, ServerConfig};

use crate::data_types::{SolrTermSummary, SolrReferenceSummary, SolrAlleleSummary};

use crate::api::term_search::{search_terms, term_complete, term_summary_by_id};

use crate::api::ref_search::search_refs;

use crate::api::allele_search::search_alleles;

use crate::api::doc_search::search_docs;
pub use crate::api::doc_search::DocSearchMatch;

pub struct Search {
    config: ServerConfig,
    reqwest_client: Client,
}


pub struct SVGPlot {
    pub bytes: Bytes
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
            reqwest_client: Client::new(),
        }
    }

    pub async fn motif_search(&self, scope: &str, pattern: &str) -> Result<String> {
        let search_url = self.config.django_url.to_owned() + "/motifsearch/query/";
        let params = [("scope", scope), ("pattern", pattern)];
        let client = reqwest::Client::new();
        let text = client.get(search_url).query(&params).send().await?.text().await?;
        Ok(text)
    }

    pub async fn term_complete(&self, cv_name: &str, q: &str)
                         -> Result<Vec<SolrTermSummary>>
    {
        Ok(term_complete(&self.config, cv_name, q).await?)
    }

    pub async fn term_summary_by_id(&self, termid: &str)
                             -> Result<Option<SolrTermSummary>>
    {
        Ok(term_summary_by_id(&self.config, termid).await?)
    }

    pub async fn ref_complete(&self, q: &str)
                        -> Result<Vec<SolrReferenceSummary>>
    {
        Ok(search_refs(&self.config, q).await?)
    }

    pub async fn allele_complete(&self, q: &str)
                           -> Result<Vec<SolrAlleleSummary>, Box<dyn Error + Send + Sync>>
    {
        search_alleles(&self.config, &self.reqwest_client, q).await
    }

    pub async fn solr_search(&self, scope: &SolrSearchScope, q: &str)
                       -> Result<SolrSearchResult>
    {
        let trimmed_query = q.trim();

        let term_matches =
            if scope.is_term() {
                search_terms(&self.config, trimmed_query).await?
            } else {
                vec![]
            };
        let ref_matches =
            if scope.is_reference() {
                search_refs(&self.config, trimmed_query).await?
            } else {
                vec![]
            };
        let doc_matches =
            if scope.is_documentation() {
                search_docs(&self.config, trimmed_query).await?
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

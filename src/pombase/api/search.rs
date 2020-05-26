use crate::web::config::{Config, ServerConfig};

use crate::data_types::{SolrTermSummary, SolrReferenceSummary};

use crate::api::term_search::{search_terms, term_complete};
pub use crate::api::term_search::TermSearchMatch;

use crate::api::ref_search::{search_refs, ref_complete};
pub use crate::api::ref_search::RefSearchMatch;

use crate::api::doc_search::search_docs;
pub use crate::api::doc_search::DocSearchMatch;

pub struct Search {
    config: ServerConfig,
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
            config: config.server.clone(),
        }
    }

    pub fn motif_search(&self, scope: &str, pattern: &str) -> Result<String, String> {
        let search_url = self.config.django_url.to_owned() + "/motifsearch/query/";
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

    pub fn term_complete(&self, cv_name: &str, q: &str)
                         -> Result<Vec<SolrTermSummary>, String>
    {
        term_complete(&self.config, cv_name, q)
    }

    pub fn ref_complete(&self, q: &str)
                        -> Result<Vec<SolrReferenceSummary>, String>
    {
        ref_complete(&self.config, q)
    }

    pub fn search_all(&self, q: &str) -> Result<SearchAllResult, String> {
        let trimmed_query = q.trim();

        let term_matches = search_terms(&self.config, trimmed_query)?;
        let ref_matches = search_refs(&self.config, trimmed_query)?;
        let doc_matches = search_docs(&self.config, trimmed_query)?;

        Ok(SearchAllResult {
            term_matches,
            ref_matches,
            doc_matches,
        })
    }
}

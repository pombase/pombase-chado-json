use anyhow::Result;

use crate::web::config::{ServerConfig, Config};


pub struct StatsPlots {
    config: ServerConfig,
}

use crate::api::search::SVGPlot;

impl StatsPlots {

    pub fn new(config: &Config) -> Self {
        Self {
            config: config.server.clone(),
        }
    }

    pub async fn gene_ex_violin_plot(&self, plot_size: &str, genes: &str)
                               -> Result<SVGPlot>
    {
        let plot_url = self.config.django_url.to_owned() + "/gene_ex/gene_ex_violin/";
        let params = [("plot_size", plot_size), ("genes", genes)];
        let client = reqwest::Client::new();
        let bytes = client.get(plot_url).query(&params).send().await?.bytes().await?;

        Ok(SVGPlot { bytes })
    }


    pub async fn get_svg_graph(&self, graph_type: &str) -> Result<SVGPlot>
    {
        let plot_url = self.config.django_url.to_owned() + "/stats/" + graph_type;
        let client = reqwest::Client::new();
        let params: &[(String, String)] = &[];
        let bytes = client.get(plot_url).query(params).send()
            .await?.bytes().await?;

        Ok(SVGPlot { bytes })
    }
}

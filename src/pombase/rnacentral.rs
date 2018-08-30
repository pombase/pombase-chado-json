use std::collections::HashSet;

use chrono::prelude::{Local, DateTime};

use web::config::Config;

use web::data::{UniquenameGeneMap, GeneDetails};

use pombase_rc_string::RcString;

pub type RNAcentralGeneSynonym = RcString;

#[derive(Serialize, Debug)]
pub struct RNAcentralGene {
#[serde(rename = "geneId")]
    pub gene_id: RcString,
#[serde(skip_serializing_if="Option::is_none")]
    pub symbol: Option<RcString>,
#[serde(skip_serializing_if="Option::is_none")]
    pub name: Option<RcString>,
    pub url: RcString,
    pub synonyms: Vec<RNAcentralGeneSynonym>,
}


#[derive(Serialize, Debug)]
pub struct RNAcentralGenomeLocation {
    pub exons: Vec<RNAcentralGenomeLocationExon>,
}

#[derive(Serialize, Debug)]
pub struct RNAcentralGenomeLocationExon {
    pub chromosome: RcString,
#[serde(rename = "startPosition")]
    pub start_position: u32,
#[serde(rename = "endPosition")]
    pub end_position: u32,
    pub strand: RcString,
}

#[derive(Serialize, Debug)]
pub struct RNAcentralNcRNA {
#[serde(rename = "primaryId")]
    pub primary_id: RcString,
#[serde(rename = "taxonId")]
    pub taxon_id: RcString,
#[serde(skip_serializing_if="Option::is_none")]
    pub symbol: Option<RcString>,
#[serde(rename = "symbolSynonyms")]
    pub symbol_synonyms: Vec<RcString>,
#[serde(rename = "soTermId")]
    pub so_term_id: RcString,
    pub sequence: Option<RcString>,
    pub url: RcString,
    pub gene: RNAcentralGene,
#[serde(skip_serializing_if="HashSet::is_empty", default)]
    pub publications: HashSet<RcString>,
}

#[derive(Serialize, Debug)]
pub struct RNAcentralMetadata {
#[serde(rename = "dateProduced")]
    pub date_produced: RcString,
#[serde(rename = "dataProvider")]
    pub data_provider: RcString,
#[serde(rename = "schemaVersion")]
    pub schema_version: RcString,
    pub publications: Vec<RcString>,
}

#[derive(Serialize, Debug)]
pub struct RNAcentral {
    data: Vec<RNAcentralNcRNA>,
#[serde(rename = "metaData")]
    metadata: RNAcentralMetadata,
}

fn gene_synonyms(gene_details: &GeneDetails) -> Vec<RcString> {
    gene_details.synonyms.iter().map(|syn| syn.name.clone()).collect::<Vec<_>>()
}

fn db_uniquename(config: &Config, gene_details: &GeneDetails) -> RcString {
    RcString::from(&format!("{}:{}", &config.database_name, gene_details.uniquename.clone()))
}

fn make_url(config: &Config, gene_details: &GeneDetails) -> RcString {
    RcString::from(&format!("{}/gene/{}", &config.base_url, &gene_details.uniquename))
}

fn make_gene_struct(config: &Config, gene_details: &GeneDetails) -> RNAcentralGene {
    RNAcentralGene {
        gene_id: db_uniquename(config, gene_details),
        symbol: gene_details.name.clone(),
        name: gene_details.product.clone(),
        url: make_url(config, gene_details),
        synonyms: gene_synonyms(gene_details),
    }
}

fn make_data(config: &Config, genes: &UniquenameGeneMap) -> Vec<RNAcentralNcRNA> {
    genes.values()
        .filter(|gene_details| {
            if let Some(ref rnacentral_config) = config.file_exports.rnacentral {
                rnacentral_config.export_so_ids.contains(&gene_details.transcript_so_termid)
            } else {
                panic!("no configuration for exporting RNAcentral data");
            }
        })
        .map(|gene_details: &GeneDetails| {
            let rnacentral_gene = make_gene_struct(config, &gene_details);
            let primary_id = db_uniquename(config, gene_details);
            let symbol_synonyms = gene_synonyms(gene_details);
            RNAcentralNcRNA {
                primary_id,
                taxon_id: RcString::from(&format!("NCBITaxon:{}", gene_details.taxonid)),
                symbol: gene_details.name.clone(),
                symbol_synonyms,
                so_term_id: gene_details.transcript_so_termid.clone(),
                sequence: gene_details.spliced_transcript_sequence()
                    .map(|s| RcString::from(&s.to_uppercase())),
                url: make_url(config, gene_details),
                gene: rnacentral_gene,
                publications: gene_details.feature_publications.clone(),
            }
        }).collect()
}

pub fn make_rnacentral_struct(config: &Config, genes: &UniquenameGeneMap) -> RNAcentral {
    let local: DateTime<Local> = Local::now();

    let database_name = config.database_name.clone();
    let data = make_data(&config, genes);

    RNAcentral {
        data,
        metadata: RNAcentralMetadata {
            date_produced: RcString::from(&local.to_rfc3339()),
            data_provider: database_name.clone(),
            schema_version: RcString::from("0.3.0"),
            publications: vec![config.database_citation.clone()],
        }
    }
}

use std::error::Error;
use std::collections::{HashSet, HashMap};
use std::fs::File;
use std::io::BufReader;

use chrono::prelude::{Local, DateTime};

use crate::web::config::Config;
use crate::data_types::{UniquenameGeneMap, GeneDetails,
                        TranscriptDetails, FeatureType,
                        RNAcentralAnnotations};

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
    pub start_position: usize,
#[serde(rename = "endPosition")]
    pub end_position: usize,
    pub strand: RcString,
}

#[derive(Serialize, Debug)]
pub struct RNAcentralNcRNALocationExon {
    pub chromosome: RcString,
#[serde(rename = "startPosition")]
    pub start_position: usize,
#[serde(rename = "endPosition")]
    pub end_position: usize,
    pub strand: RcString,
}

#[derive(Serialize, Debug)]
pub struct RNAcentralNcRNALocation {
    pub assembly: RcString,
    pub exons: Vec<RNAcentralNcRNALocationExon>,
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
    pub sequence: RcString,
    pub url: RcString,
    pub gene: RNAcentralGene,
#[serde(rename = "genomeLocations")]
    pub genome_locations: Vec<RNAcentralNcRNALocation>,
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

fn db_uniquename(config: &Config, uniquename: &RcString) -> RcString {
    RcString::from(&format!("{}:{}", &config.database_name, uniquename))
}

fn make_url(config: &Config, gene_details: &GeneDetails) -> RcString {
    RcString::from(&format!("{}/gene/{}", &config.base_url, &gene_details.uniquename))
}

fn make_gene_struct(config: &Config, gene_details: &GeneDetails) -> RNAcentralGene {
    RNAcentralGene {
        gene_id: db_uniquename(config, &gene_details.uniquename),
        symbol: gene_details.name.clone(),
        name: gene_details.product.clone(),
        url: make_url(config, gene_details),
        synonyms: gene_synonyms(gene_details),
    }
}

fn make_genome_location(config: &Config, gene_details: &GeneDetails,
                        transcript_details: &TranscriptDetails)
                        -> RNAcentralNcRNALocation
{
    let assembly =
        config.organisms.iter()
        .filter(|org| org.taxonid == gene_details.taxonid)
        .collect::<Vec<_>>().get(0)
        .unwrap_or_else(|| panic!("organism not found in configuration: {}", gene_details.taxonid))
        .assembly_version.clone()
        .unwrap_or_else(|| panic!("no assembly_version for: {}", gene_details.taxonid));

    let mut exons = vec![];
    for part in &transcript_details.parts {
        if part.feature_type == FeatureType::Exon {
            let mut start_position = part.location.start_pos - 1;
            let mut end_position = part.location.end_pos;
            if start_position > end_position {
                use std::mem;
                mem::swap(&mut start_position, &mut end_position);
            }
            let chromosome_name = &part.location.chromosome_name;
            let chromosome =
                config.find_chromosome_config(chromosome_name).export_id.clone();
            exons.push(RNAcentralNcRNALocationExon {
                chromosome,
                start_position,
                end_position,
                strand: RcString::from(part.location.strand.to_gff_str()),
            });
        }
    }

    RNAcentralNcRNALocation {
        assembly,
        exons,
    }
}

fn make_data(config: &Config, genes: &UniquenameGeneMap) -> Vec<RNAcentralNcRNA> {
    let rnacentral_config = config.file_exports.rnacentral.clone().unwrap();

    let mut ret = vec![];

    for gene_details in genes.values() {
        let so_term_id = &gene_details.transcript_so_termid;

        if !rnacentral_config.export_so_ids.contains(so_term_id) {
            continue;
        }

        for transcript_details in &gene_details.transcripts {
            let rnacentral_gene = make_gene_struct(config, gene_details);

            let primary_id = db_uniquename(config, &transcript_details.uniquename);
            let location = make_genome_location(config, gene_details, transcript_details);

            let uppercase_sequence =
                transcript_details.spliced_transcript_sequence().to_uppercase();

            ret.push(RNAcentralNcRNA {
                primary_id,
                taxon_id: RcString::from(&format!("NCBITaxon:{}", gene_details.taxonid)),
                symbol: None,
                symbol_synonyms: vec![],
                so_term_id: so_term_id.clone(),
                sequence: RcString::from(&uppercase_sequence),
                url: make_url(config, gene_details),
                gene: rnacentral_gene,
                genome_locations: vec![location],
                publications: gene_details.feature_publications.clone(),
            })
        }
    }

    ret
}

pub fn make_rnacentral_struct(config: &Config, genes: &UniquenameGeneMap) -> RNAcentral {
    let local: DateTime<Local> = Local::now();

    let data_provider = config.database_name.clone();
    let data = make_data(config, genes);

    RNAcentral {
        data,
        metadata: RNAcentralMetadata {
            date_produced: RcString::from(&local.to_rfc3339()),
            data_provider,
            schema_version: RcString::from("0.3.0"),
            publications: vec![config.database_citation.clone()],
        }
    }
}


#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RfamAnnotation {
#[serde(rename = "URS-Id")]
    pub urs_identifier: String,
#[serde(rename = "Rfam-Model-Id")]
    pub rfam_model_id: String,
#[serde(rename = "Score")]
    pub score: f32,
#[serde(rename = "E-value")]
    pub e_value: f32,
#[serde(rename = "Sequence-Start")]
    pub sequence_start: u32,
#[serde(rename = "Sequence-Stop")]
    pub sequence_stop: u32,
#[serde(rename = "Model-Start")]
    pub model_start: u32,
#[serde(rename = "Model-Stop")]
    pub model_stop: u32,
#[serde(rename = "Rfam-Model-Description")]
    pub rfam_model_description: String,
}

pub fn parse_annotation_json(file_name: &str)
    -> Result<RNAcentralAnnotations, Box<dyn Error>>
{
    let file = match File::open(file_name) {
        Ok(file) => file,
        Err(err) => {
            eprintln!("Failed to read {}: {}", file_name, err);
            return Err(Box::new(err));
        }
    };
    let reader = BufReader::new(file);

    let annotation: HashMap<RcString, Vec<RfamAnnotation>> =
        match serde_json::from_reader(reader) {
            Ok(annotation) => annotation,
            Err(err) => {
                eprint!("failed to parse {}: {}", file_name, err);
                return Err(Box::new(err))
            },
        };

    Ok(annotation)
}

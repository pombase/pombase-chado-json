extern crate tokio_postgres;
extern crate serde_json;

use std::cmp::min;
use std::fs::{File, create_dir_all};
use std::io::{Write, BufWriter};
use std::io;
use std::collections::{HashMap, HashSet};
use regex::Regex;

use flate2::Compression;
use flate2::write::GzEncoder;

use crate::bio::util::{format_fasta, format_gene_gff, format_misc_feature_gff};

use pombase_rc_string::RcString;

use crate::web::config::*;
use crate::rnacentral::*;

use crate::types::CvName;
use crate::data_types::*;
use crate::annotation_util::table_for_export;

use crate::bio::go_format_writer::*;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct WebData {
    pub metadata: Metadata,
    pub chromosomes: ChrNameDetailsMap,
    pub chromosome_summaries: Vec<ChromosomeShort>,
    pub recent_references: RecentReferences,
    pub all_community_curated: Vec<ReferenceShort>,
    pub all_admin_curated: Vec<ReferenceShort>,
    pub api_maps: APIMaps,
    pub solr_data: SolrData,
    pub search_gene_summaries: Vec<GeneSummary>,
    pub ont_annotations: Vec<OntAnnotation>,
    pub stats: Stats,
}

const FASTA_SEQ_COLUMNS: usize = 60;

fn write_as_fasta(writer: &mut dyn Write, id: &str, desc: Option<String>,
                  seq: &str)
{
    let fasta = format_fasta(id, desc, seq, FASTA_SEQ_COLUMNS);
    writer.write_all(fasta.as_bytes()).unwrap();
}

impl WebData {
    fn get_chromosomes(&self) -> &ChrNameDetailsMap {
        &self.chromosomes
    }

    fn create_dir(&self, output_dir: &str, dir_name: &str) -> String {
        let path = String::new() + output_dir + "/" + dir_name;
        create_dir_all(&path).unwrap_or_else(|why| {
            println!("Creating output directory failed: {:?}", why.kind());
        });
        path
    }

    fn write_chromosome_seq_chunks(&self, output_dir: &str, chunk_sizes: &[usize]) {
        for chunk_size in chunk_sizes {
            for (chromosome_uniquename, chromosome_details) in &self.chromosomes {
                let new_path_part = &format!("{}/sequence/{}", chromosome_uniquename, chunk_size);
                let chr_path = self.create_dir(output_dir, new_path_part);
                let mut index = 0;
                let max_index = chromosome_details.residues.len() / chunk_size;
                while index <= max_index {
                    let start_pos = index*chunk_size;
                    let end_pos = min(start_pos+chunk_size, chromosome_details.residues.len());
                    let chunk: String = chromosome_details.residues[start_pos..end_pos].into();
                    let file_name = format!("{}/chunk_{}", chr_path, index);
                    let f = File::create(file_name).expect("Unable to open file");
                    let mut writer = BufWriter::new(&f);
                    writer.write_all(chunk.as_bytes()).expect("Unable to write chromosome chunk");
                    index += 1;
                }
            }
        }
    }

    fn write_chromosome_json(&self, config: &Config, output_dir: &str) {
        let new_path = self.create_dir(output_dir, "chromosome");
        for (chromosome_uniquename, chromosome_details) in &self.chromosomes {
            let s = serde_json::to_string(&chromosome_details).unwrap();
            let file_name = format!("{}/{}.json", new_path, &chromosome_uniquename);
            let f = File::create(file_name).expect("Unable to open file");
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes()).expect("Unable to write chromosome JSON");
        }
        self.write_chromosome_seq_chunks(&new_path, &config.api_seq_chunk_sizes);
    }

    fn write_gene_summaries(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.search_gene_summaries).unwrap();
        let file_name = String::new() + output_dir + "/gene_summaries.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write gene_summaries.json");
    }

    fn write_metadata(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.metadata).unwrap();
        let file_name = String::new() + output_dir + "/metadata.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write metadata.json");
    }

    fn write_recent_references(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.recent_references).unwrap();
        let file_name = String::new() + output_dir + "/recent_references.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write recent references JSON");
    }

    fn write_all_community_curated(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.all_community_curated).unwrap();
        let file_name = String::new() + output_dir + "/community_curated_references.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write recent references JSON");
    }

    fn write_all_admin_curated(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.all_admin_curated).unwrap();
        let file_name = String::new() + output_dir + "/admin_curated_references.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write admin curated refs JSON");
    }

    fn write_api_maps(&self, output_dir: &str) {
        let file_name = String::new() + output_dir + "/api_maps.json.gz";
        let f = File::create(file_name).expect("Unable to open file");

        let mut compressor = GzEncoder::new(f, Compression::default());
        serde_json::ser::to_writer(&mut compressor, &self.api_maps).unwrap();

        compressor.finish().unwrap();
    }

    fn write_solr_data(&self, output_dir: &str) {
        let new_path = self.create_dir(output_dir, "solr_data/");

        let terms = self.solr_data.term_summaries.clone();

        let terms_json_text = serde_json::to_string(&terms).unwrap();
        let terms_file_name = format!("{}/terms.json.gz", new_path);
        let terms_file = File::create(terms_file_name).expect("Unable to open file");

        let mut terms_compressor = GzEncoder::new(terms_file, Compression::default());
        terms_compressor.write_all(terms_json_text.as_bytes()).expect("Unable to write terms as JSON");
        terms_compressor.finish().expect("Unable to write terms as JSON");

        let genes = self.solr_data.gene_summaries.clone();

        let genes_json_text = serde_json::to_string(&genes).unwrap();
        let genes_file_name = format!("{}/genes.json.gz", new_path);
        let genes_file = File::create(genes_file_name).expect("Unable to open file");

        let mut genes_compressor = GzEncoder::new(genes_file, Compression::default());
        genes_compressor.write_all(genes_json_text.as_bytes()).expect("Unable to write genes as JSON");
        genes_compressor.finish().expect("Unable to write genes as JSON");

        let references = self.solr_data.reference_summaries.clone();

        let references_json_text = serde_json::to_string(&references).unwrap();
        let references_file_name = format!("{}/references.json.gz", new_path);
        let references_file = File::create(references_file_name).expect("Unable to open file");

        let mut references_compressor = GzEncoder::new(references_file, Compression::default());
        references_compressor.write_all(references_json_text.as_bytes()).expect("Unable to write references as JSON");
        references_compressor.finish().expect("Unable to write references as JSON");
    }


    fn write_intermine_data(&self, config: &Config, output_dir: &str)
                            -> Result<(), io::Error>
    {
        let load_org_taxonid =
            if let Some(load_org_taxonid) = config.load_organism_taxonid {
                load_org_taxonid
            } else {
                return Ok(())
            };

        let intermine_genes: Vec<_> =
            self.api_maps.genes.values()
            .filter(|gene_details| gene_details.taxonid == load_org_taxonid)
            .map(|gene_details| {
                InterMineGeneDetails::from_gene_details(gene_details)
            }).collect();

        let genes_json_text = serde_json::to_string(&intermine_genes)?;
        let genes_file_name = format!("{}/pombemine_gene_details.gz", output_dir);
        let genes_file = File::create(genes_file_name)?;

        let mut genes_compressor = GzEncoder::new(genes_file, Compression::default());
        genes_compressor.write_all(genes_json_text.as_bytes())?;
        genes_compressor.finish()?;

        Ok(())
    }

    fn write_subsets(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.api_maps.term_subsets).unwrap();
        let file_name = String::new() + output_dir + "/term_subsets.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write");

        let s = serde_json::to_string(&self.api_maps.gene_subsets).unwrap();
        let file_name = String::new() + output_dir + "/gene_subsets.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write");
    }

    fn write_feature_sequences(&self, output_dir: &str) {
        let make_seq_writer = |name: &str| {
            let file_name = String::new() + output_dir + "/" + name;
            let file = File::create(file_name).expect("Unable to open file");
            BufWriter::new(file)
        };

        let mut cds_writer = make_seq_writer("cds.fa");
        let mut cds_introns_writer = make_seq_writer("cds+introns.fa");
        let mut cds_introns_utrs_writer = make_seq_writer("cds+introns+utrs.fa");
        let mut introns_writer = make_seq_writer("introns_within_cds.fa");
        let mut five_prime_utrs_writer = make_seq_writer("five_prime_utrs.fa");
        let mut three_prime_utrs_writer = make_seq_writer("three_prime_utrs.fa");
        let mut peptide_writer = make_seq_writer("peptide.fa");

        for (gene_uniquename, gene_details) in &self.api_maps.genes {
            if let Some(transcript) = gene_details.transcripts.get(0) {
                let mut cds_seq = String::new();
                let mut cds_introns_seq = String::new();
                let mut cds_introns_utrs_seq = String::new();
                let mut five_prime_utr_seq = String::new();
                let mut three_prime_utr_seq = String::new();
                for part in &transcript.parts {
                    if part.feature_type == FeatureType::Exon {
                        cds_seq += &part.residues;
                        cds_introns_seq += &part.residues;
                    }
                    if part.feature_type == FeatureType::CdsIntron {
                        cds_introns_seq += &part.residues;
                    }
                    if part.feature_type == FeatureType::CdsIntron {
                        write_as_fasta(&mut introns_writer, &part.uniquename,
                                       Some(String::from(gene_uniquename)),
                                       &part.residues);
                    }
                    cds_introns_utrs_seq += &part.residues;
                    if part.feature_type == FeatureType::FivePrimeUtr {
                        five_prime_utr_seq += &part.residues;
                    }
                    if part.feature_type == FeatureType::ThreePrimeUtr {
                        three_prime_utr_seq += &part.residues;
                    }
                }

                write_as_fasta(&mut cds_writer, gene_uniquename, None, &cds_seq);
                write_as_fasta(&mut cds_introns_writer, gene_uniquename, None, &cds_introns_seq);
                write_as_fasta(&mut cds_introns_utrs_writer,
                               gene_uniquename, None, &cds_introns_utrs_seq);
                if !five_prime_utr_seq.is_empty() {
                    write_as_fasta(&mut five_prime_utrs_writer,
                                   gene_uniquename, None, &five_prime_utr_seq);
                }
                if !three_prime_utr_seq.is_empty() {
                    write_as_fasta(&mut three_prime_utrs_writer,
                                   gene_uniquename, None, &three_prime_utr_seq);
                }
                if let Some(ref protein) = transcript.protein {
                    let name_and_product =
                        if gene_details.name.is_some() || gene_details.product.is_some() {
                            let mut buf = String::new();
                            if let Some(ref name) = gene_details.name {
                                buf.push_str(name);
                            }
                            buf.push('|');
                            if let Some(ref product) = gene_details.product {
                                buf.push_str(product);
                            }
                            Some(buf.to_owned())
                        } else {
                            None
                        };
                    write_as_fasta(&mut peptide_writer, &protein.uniquename,
                                   name_and_product, &protein.sequence);
                }
            }
        }

        cds_writer.flush().unwrap();
        cds_introns_writer.flush().unwrap();
        cds_introns_utrs_writer.flush().unwrap();
        introns_writer.flush().unwrap();
        peptide_writer.flush().unwrap();
        five_prime_utrs_writer.flush().unwrap();
        three_prime_utrs_writer.flush().unwrap();
    }

    pub fn write_chromosome_sequences(&self, config: &Config, output_dir: &str) {
        let make_seq_writer = |name: &str| {
            let file_name = String::new() + output_dir + "/" + name;
            let file = File::create(file_name).expect("Unable to open file");
            BufWriter::new(file)
        };

        if let Some(load_org) = config.load_organism() {
            let load_org_name = load_org.full_name();
            let chromosomes_file_name = load_org_name.clone() + "_all_chromosomes.fa";
            let mut chromosomes_writer = make_seq_writer(&chromosomes_file_name);

            for (uniquename, details) in &self.chromosomes {
                let chr_config = config.find_chromosome_config(uniquename);
                write_as_fasta(&mut chromosomes_writer, &chr_config.export_id,
                               Some(load_org_name.clone()), &details.residues);
                let this_chr_file_name =
                    load_org_name.clone() + "_" + &chr_config.export_file_id + ".fa";
                let mut this_chr_writer = make_seq_writer(&this_chr_file_name);
                write_as_fasta(&mut this_chr_writer, &chr_config.export_id,
                               Some(load_org_name.clone()), &details.residues);
                this_chr_writer.flush().unwrap();

            }

            chromosomes_writer.flush().unwrap();
        }
    }

    fn write_chromosome_summaries(&self, output_dir: &str) {
        let s = serde_json::to_string(&self.chromosome_summaries).unwrap();
        let file_name = String::new() + output_dir + "/chromosome_summaries.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write chromosome_summaries.json");
    }


    fn write_go_annotation_files(&self, config: &Config, go_eco_mappping: &GoEcoMapping,
                                 output_dir: &str)
                                 -> Result<(), io::Error>
    {
        let load_org_taxonid =
            if let Some(load_org_taxonid) = config.load_organism_taxonid {
                load_org_taxonid
            } else {
                return Ok(())
            };

        let database_name = &config.database_name;

        let gpi_file_name =
            format!("{}/gene_product_information_taxonid_{}.tsv", output_dir,
                    load_org_taxonid);
        let gpad_file_name =
            format!("{}/gene_product_annotation_data_taxonid_{}.tsv", output_dir,
                    load_org_taxonid);
        let pombase_gaf_file_name = format!("{}/pombase_style_gaf.tsv", output_dir);
        let standard_gaf_file_name = format!("{}/go_style_gaf.tsv", output_dir);

        let gpi_file = File::create(gpi_file_name).expect("Unable to open file");
        let gpad_file = File::create(gpad_file_name).expect("Unable to open file");
        let pombase_gaf_file =
            File::create(pombase_gaf_file_name).expect("Unable to open file");
        let standard_gaf_file =
            File::create(standard_gaf_file_name).expect("Unable to open file");
        let mut gpi_writer = BufWriter::new(&gpi_file);
        let mut gpad_writer = BufWriter::new(&gpad_file);
        let mut pombase_gaf_writer = BufWriter::new(&pombase_gaf_file);
        let mut standard_gaf_writer = BufWriter::new(&standard_gaf_file);

        let generated_by = format!("!generated-by: {}\n", database_name);
        let iso_date = self.metadata.db_creation_datetime.replace(" ", "T");
        let date_generated = format!("!date-generated: {}\n", &iso_date);
        let url_header = format!("!URL: {}\n", &config.base_url);
        let funding_header = format!("!funding: {}\n", &config.funder);

        gpi_writer.write_all("!gpi-version: 2.0\n".as_bytes())?;
        gpi_writer.write_all(format!("!namespace: {}\n", database_name).as_bytes())?;
        gpi_writer.write_all(generated_by.as_bytes())?;
        gpi_writer.write_all(date_generated.as_bytes())?;
        gpi_writer.write_all(url_header.as_bytes())?;
        gpi_writer.write_all(funding_header.as_bytes())?;

        gpad_writer.write_all("!gpa-version: 2.0\n".as_bytes())?;
        gpad_writer.write_all(generated_by.as_bytes())?;
        gpad_writer.write_all(date_generated.as_bytes())?;
        gpad_writer.write_all(url_header.as_bytes())?;
        gpad_writer.write_all(funding_header.as_bytes())?;

        standard_gaf_writer.write_all("!gaf-version: 2.2\n".as_bytes())?;
        standard_gaf_writer.write_all(generated_by.as_bytes())?;
        standard_gaf_writer.write_all(date_generated.as_bytes())?;
        standard_gaf_writer.write_all(url_header.as_bytes())?;
        let contact = format!("!contact: {}\n", &config.helpdesk_address);
        standard_gaf_writer.write_all(contact.as_bytes())?;

        for gene_details in self.api_maps.genes.values() {
            if gene_details.taxonid != load_org_taxonid {
                continue;
            }


            if gene_details.feature_type == "ncRNA gene" {
                let mut found = false;
                for aspect in GO_ASPECT_NAMES.iter() {
                    let term_annotations = gene_details.cv_annotations.get(&RcString::from(*aspect));
                    if term_annotations.is_some() {
                        found = true;
                        break;
                    }
                }
                if !found {
                    // special case for ncRNA genes: only write ND lines for an aspect
                    // if there are some annotations for the other aspects
                    continue;
                }
            }

            let db_object_id = format!("{}:{}", database_name, gene_details.uniquename);
            let db_object_symbol =
                gene_details.product.clone().unwrap_or_else(RcString::new);
            let db_object_name =
                gene_details.name.clone().unwrap_or_else(RcString::new);

            let db_object_synonyms =
                gene_details.synonyms.iter().filter(|synonym| {
                    synonym.synonym_type == "exact"
                })
                .map(|synonym| synonym.name.to_string())
                .collect::<Vec<String>>()
                .join("|");

            let db_object_type = config.file_exports.gpad_gpi.transcript_gene_so_term_map
                .get(gene_details.transcript_so_termid.as_str())
                .unwrap_or_else(|| {
                    panic!("failed for find configuration for {} in transcript_gene_so_term_map",
                           &gene_details.transcript_so_termid);
                });

            let db_object_taxon = format!("NCBITaxon:{}", load_org_taxonid);
            let db_xrefs =
                if let Some(ref uniprot_id) = gene_details.uniprot_identifier {
                    RcString::from(&format!("UniProtKB:{}", uniprot_id))
                } else {
                    RcString::from("")
                };

            let gpi_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t\t\t\t{}\tgo-annotation-summary={}\n",
                                   db_object_id,
                                   db_object_name,
                                   db_object_name,
                                   db_object_synonyms,
                                   db_object_type,
                                   db_object_taxon,
                                   db_xrefs,
                                   db_object_symbol);
            gpi_writer.write_all(gpi_line.as_bytes())?;

            write_gene_product_annotation(&mut gpad_writer, go_eco_mappping, config,
                                          &self.api_maps, gene_details)?;

            for aspect_name in &GO_ASPECT_NAMES {
                write_go_annotation_format(&mut pombase_gaf_writer, config,
                                           GafWriteMode::PomBase,
                                           &self.api_maps, gene_details,
                                           aspect_name)?;
                write_go_annotation_format(&mut standard_gaf_writer, config,
                                           GafWriteMode::Standard,
                                           &self.api_maps, gene_details,
                                           aspect_name)?;
            }
       }

        Ok(())
    }

    fn write_gene_id_table(&self, config: &Config, output_dir: &str) -> Result<(), io::Error> {
        let gene_file_name = output_dir.to_owned() + "/sysID2product.tsv";
        let rna_file_name = output_dir.to_owned() + "/sysID2product.rna.tsv";
        let pseudogenes_file_name = output_dir.to_owned() + "/pseudogeneIDs.tsv";
        let all_names_file_name = output_dir.to_owned() + "/gene_IDs_names.tsv";
        let all_ids_file_name = output_dir.to_owned() + "/gene_IDs_names_products.tsv";

        let gene_file = File::create(gene_file_name).expect("Unable to open file");
        let rna_file = File::create(rna_file_name).expect("Unable to open file");
        let pseudogenes_file = File::create(pseudogenes_file_name).expect("Unable to open file");
        let all_names_file = File::create(all_names_file_name).expect("Unable to open file");
        let all_ids_file = File::create(all_ids_file_name).expect("Unable to open file");

        let mut gene_writer = BufWriter::new(&gene_file);
        let mut rna_writer = BufWriter::new(&rna_file);
        let mut pseudogenes_writer = BufWriter::new(&pseudogenes_file);
        let mut all_names_writer = BufWriter::new(&all_names_file);
        let mut all_ids_writer = BufWriter::new(&all_ids_file);

        let db_version = format!("# Chado database date: {}\n", self.metadata.db_creation_datetime);
        gene_writer.write_all(db_version.as_bytes())?;
        rna_writer.write_all(db_version.as_bytes())?;
        pseudogenes_writer.write_all(db_version.as_bytes())?;
        all_names_writer.write_all(db_version.as_bytes())?;

        for gene_details in self.api_maps.genes.values() {
            if let Some(load_org_taxonid) = config.load_organism_taxonid {
                if gene_details.taxonid != load_org_taxonid {
                    continue;
                }
            }

            let synonyms =
                gene_details.synonyms.iter().filter(|synonym| {
                    synonym.synonym_type == "exact"
                })
                .map(|synonym| synonym.name.to_string())
                .collect::<Vec<String>>()
                .join(",");

            let line = format!("{}\t{}\t{}\n",
                               gene_details.uniquename,
                               gene_details.name.clone().unwrap_or_else(RcString::new),
                               synonyms);

            let gene_name = if let Some(ref gene_details_name) = gene_details.name {
                gene_details_name.clone()
            } else {
                RcString::new()
            };

            let gene_product = if let Some(ref gene_details_product) = gene_details.product {
                gene_details_product.clone()
            } else {
                RcString::new()
            };

            let line_with_product = format!("{}\t{}\t{}\t{}\n",
                                            gene_details.uniquename,
                                            gene_name,
                                            synonyms,
                                            gene_product);

            all_names_writer.write_all(line.as_bytes())?;

            if gene_details.feature_type == "pseudogene" {
                pseudogenes_writer.write_all(line.as_bytes())?;
            } else {
                if gene_details.feature_type == "mRNA gene" {
                    gene_writer.write_all(line_with_product.as_bytes())?;
                } else {
                    if gene_details.feature_type.contains("RNA") {
                        rna_writer.write_all(line_with_product.as_bytes())?;
                    }
                }
            }

            let uniprot_id =
                if let Some(ref gene_uniprot_id) = gene_details.uniprot_identifier {
                    gene_uniprot_id
                } else {
                    ""
                };

            let chromosome_name =
                if let Some(ref loc) = gene_details.location {
                    &loc.chromosome_name
                } else {
                    ""
                };

            let gene_type =
                if gene_details.feature_type == "mRNA gene" {
                    "protein coding gene"
                } else {
                    &gene_details.feature_type
                };

            let all_ids_line = format!("{}\t{}:{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                       gene_details.uniquename,
                                       config.database_name,
                                       gene_details.uniquename,
                                       gene_name,
                                       chromosome_name,
                                       gene_product,
                                       uniprot_id,
                                       gene_type,
                                       synonyms);
            all_ids_writer.write_all(all_ids_line.as_bytes())?;
        }

        gene_writer.flush()?;
        rna_writer.flush()?;
        pseudogenes_writer.flush()?;
        all_names_writer.flush()?;

        Ok(())
    }

    fn write_protein_features(&self, config: &Config, output_dir: &str)
                              -> Result<(), io::Error>
    {
        let peptide_stats_name = format!("{}/PeptideStats.tsv", output_dir);
        let peptide_stats_file = File::create(peptide_stats_name).expect("Unable to open file");
        let mut peptide_stats_writer = BufWriter::new(&peptide_stats_file);

        let peptide_stats_header = "Systematic_ID\tMass (kDa)\tpI\tCharge\tResidues\tCAI\n";
        peptide_stats_writer.write_all(peptide_stats_header.as_bytes())?;

        let protein_features_name = format!("{}/ProteinFeatures.tsv", output_dir);
        let protein_features_file = File::create(protein_features_name).expect("Unable to open file");
        let mut protein_features_writer = BufWriter::new(&protein_features_file);

        let disordered_regions_name = format!("{}/disordered_regions.tsv", output_dir);
        let disordered_regions_file = File::create(disordered_regions_name).expect("Unable to open file");
        let mut disordered_regions_writer = BufWriter::new(&disordered_regions_file);

        let aa_composition_name = format!("{}/aa_composition.tsv", output_dir);
        let aa_composition_file = File::create(aa_composition_name).expect("Unable to open file");
        let mut aa_composition_writer = BufWriter::new(&aa_composition_file);

        let protein_features_header =
            "systematic_id\tgene_name\tpeptide_id\tdomain_id\tdatabase\tseq_start\tseq_end\n";
        protein_features_writer.write_all(protein_features_header.as_bytes())?;

        let disordered_regions_header =
            "systematic_id\tgene_name\tpeptide_id\tseq_start\tseq_end\n";
        disordered_regions_writer.write_all(disordered_regions_header.as_bytes())?;

        let db_display_name = |db_alias: &str| {
            if let Some(name) = config.extra_database_aliases.get(&db_alias.to_lowercase()) {
                RcString::from(name)
            } else {
                RcString::from(db_alias)
            }
        };

        type AAComposition = HashMap<char, u32>;

        let mut total_composition: AAComposition = HashMap::new();

        let prot_composition =
            |total_composition: &mut AAComposition, protein: &ProteinDetails|
        {
            let mut composition = HashMap::new();
            for c in protein.sequence.chars() {
                let count = composition.entry(c).or_insert(0);
                *count += 1;
                let total_count = total_composition.entry(c).or_insert(0);
                *total_count += 1;
            }
            composition
        };

        let mut compositions_to_write = vec![];

        for (gene_uniquename, gene_details) in &self.api_maps.genes {
            if let Some(transcript) = gene_details.transcripts.get(0) {
                if let Some(ref protein) = transcript.protein {
                    let line = format!("{}\t{:.2}\t{}\t{}\t{}\t{}\n",
                                       gene_uniquename, protein.molecular_weight,
                                       protein.isoelectric_point,
                                       protein.charge_at_ph7,
                                       protein.sequence.len() - 1,
                                       protein.codon_adaptation_index);
                    peptide_stats_writer.write_all(line.as_bytes())?;

                    let gene_name = gene_details.name.clone().unwrap_or_else(RcString::new);
                    for interpro_match in &gene_details.interpro_matches {
                        let line_start = format!("{}\t{}\t{}\t{}\t{}",
                                                 gene_uniquename, gene_name,
                                                 protein.uniquename, interpro_match.id,
                                                 db_display_name(&interpro_match.dbname));
                        for location in &interpro_match.locations {
                            let line = format!("{}\t{}\t{}\n", line_start,
                                               location.start, location.end);
                            protein_features_writer.write_all(line.as_bytes())?;
                        }
                    }

                    for disordered_region in &gene_details.disordered_region_coords {
                        let (start_pos, end_pos) = disordered_region;
                        let line = format!("{}\t{}\t{}\t{}\t{}\n",
                                           gene_uniquename, gene_name,
                                           protein.uniquename,
                                           start_pos, end_pos);
                        disordered_regions_writer.write_all(line.as_bytes())?;
                    }

                    let composition = prot_composition(&mut total_composition, protein);

                    compositions_to_write.push((gene_uniquename.clone(), composition));
                }
            }

        }

        let mut all_composition_aa: Vec<char> = vec![];

        for ch in total_composition.keys() {
            if *ch != '*' {
                all_composition_aa.push(*ch);
            }
        }

        all_composition_aa.sort_unstable();

        let all_composition_string =
            all_composition_aa.iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>().join("\t");

        let composition_header = "Systematic_ID\t".to_owned() +
            &all_composition_string + "\n";
        aa_composition_writer.write_all(composition_header.as_bytes())?;

        let composition_line = |first_col_string: RcString, comp: &AAComposition| {
            let mut line = String::from(first_col_string);

            for ch in &all_composition_aa {
                line.push('\t');
                if let Some(count) = comp.get(ch) {
                    line.push_str(&count.to_string());
                } else {
                    line.push('0');
                }
            }
            line.push('\n');
            line
        };

        for (gene_uniquename, comp) in compositions_to_write.drain(0..) {
            let line = composition_line(gene_uniquename, &comp);
            aa_composition_writer.write_all(line.as_bytes())?;
        }

        let composition_total_line =
            composition_line(RcString::from("total"), &total_composition);
        aa_composition_writer.write_all(composition_total_line.as_bytes())?;

        peptide_stats_writer.flush()?;

        Ok(())
    }

    fn write_feature_coords(&self, config: &Config, output_dir: &str)
                            -> Result<(), io::Error>
    {
        let write_line =
            |uniquename: &str, location: &ChromosomeLocation,
             writer: &mut BufWriter<&File>| {
                let display_strand =
                    if location.strand == Strand::Forward {1} else {-1};
                let line = format!("{}\t{}\t{}\t{}\n",
                                   uniquename, location.start_pos,
                                   location.end_pos, display_strand);
                writer.write(line.as_bytes())
        };

        for (chr_uniquename, chr_details) in &self.chromosomes {
            if let Some(load_org_taxonid) = config.load_organism_taxonid {
                if chr_details.taxonid != load_org_taxonid {
                    continue;
                }
            }

            let gene_file_name = format!("{}/{}.gene.coords.tsv", output_dir, chr_uniquename);
            let cds_file_name = format!("{}/{}.cds.coords.tsv", output_dir, chr_uniquename);
            let exon_file_name = format!("{}/{}.exon.coords.tsv", output_dir, chr_uniquename);

            let gene_file = File::create(gene_file_name).expect("Unable to open file");
            let cds_file = File::create(cds_file_name).expect("Unable to open file");
            let exon_file = File::create(exon_file_name).expect("Unable to open file");

            let mut gene_writer = BufWriter::new(&gene_file);
            let mut cds_writer = BufWriter::new(&cds_file);
            let mut exon_writer = BufWriter::new(&exon_file);

            for gene_uniquename in &chr_details.gene_uniquenames {
                let gene = &self.api_maps.genes[gene_uniquename];
                if let Some(ref gene_location) = gene.location {
                    write_line(gene_uniquename, gene_location, &mut gene_writer)?;

                    for transcript in &gene.transcripts {
                        if let Some(ref cds_location) = transcript.cds_location {
                            write_line(gene_uniquename, cds_location, &mut cds_writer)?;
                        }

                        let is_forward =
                            transcript.parts[0].location.strand == Strand::Forward;

                        let parts: Vec<FeatureShort> = if is_forward {
                            transcript.parts.iter().cloned().filter(|part| {
                                part.feature_type == FeatureType::Exon ||
                                    part.feature_type == FeatureType::FivePrimeUtr ||
                                    part.feature_type == FeatureType::ThreePrimeUtr
                            }).collect()
                        } else {
                            transcript.parts.iter().cloned().rev().filter(|part| {
                                part.feature_type == FeatureType::Exon ||
                                    part.feature_type == FeatureType::FivePrimeUtr ||
                                    part.feature_type == FeatureType::ThreePrimeUtr
                            }).collect()
                        };

                        // merge Exon and (Three|Five)PrimeUtrs that abut
                        let mut merged_locs: Vec<ChromosomeLocation> = vec![];

                        for part in parts {
                            if let Some(prev) = merged_locs.pop() {
                                if prev.end_pos + 1 == part.location.start_pos {
                                    merged_locs.push(ChromosomeLocation {
                                        start_pos: prev.start_pos,
                                        end_pos: part.location.end_pos,
                                        chromosome_name: prev.chromosome_name,
                                        strand: prev.strand,
                                        phase: prev.phase,
                                    });
                                } else {
                                    merged_locs.push(prev);
                                    merged_locs.push(part.location);
                                }
                            } else {
                                merged_locs.push(part.location);
                            }
                        }

                        for loc in merged_locs {
                            write_line(gene_uniquename, &loc, &mut exon_writer)?;
                        }
                    }
                }
            }

            gene_writer.flush()?;
            cds_writer.flush()?;
            exon_writer.flush()?;
        }

        Ok(())
    }

    pub fn write_gff(&self, config: &Config, output_dir: &str)
                         -> Result<(), io::Error>
    {
        if let Some(load_org) = config.load_organism() {
            let load_org_name = load_org.full_name();

            let all_gff_name = format!("{}/{}_all_chromosomes.gff3", output_dir, load_org_name);
            let all_gff_file = File::create(all_gff_name).expect("Unable to open file");
            let mut all_gff_writer = BufWriter::new(&all_gff_file);

            let forward_features_gff_name =
                format!("{}/{}_all_chromosomes_forward_strand.gff3", output_dir, load_org_name);
            let forward_features_gff_file = File::create(forward_features_gff_name).expect("Unable to open file");
            let mut forward_features_gff_writer = BufWriter::new(&forward_features_gff_file);

            let reverse_features_gff_name =
                format!("{}/{}_all_chromosomes_reverse_strand.gff3", output_dir, load_org_name);
            let reverse_features_gff_file = File::create(reverse_features_gff_name).expect("Unable to open file");
            let mut reverse_features_gff_writer = BufWriter::new(&reverse_features_gff_file);

            let unstranded_features_gff_name =
                format!("{}/{}_all_chromosomes_unstranded.gff3", output_dir, load_org_name);
            let unstranded_features_gff_file = File::create(unstranded_features_gff_name).expect("Unable to open file");
            let mut unstranded_features_gff_writer = BufWriter::new(&unstranded_features_gff_file);

            all_gff_writer.write_all(b"##gff-version 3\n")?;
            forward_features_gff_writer.write_all(b"##gff-version 3\n")?;
            reverse_features_gff_writer.write_all(b"##gff-version 3\n")?;
            unstranded_features_gff_writer.write_all(b"##gff-version 3\n")?;

            let mut chr_writers = HashMap::new();

            let make_chr_gff_writer = |export_name: &str| {
                let file_name = String::new() +
                    output_dir + "/" + &load_org_name + "_" + export_name + ".gff3";
                let file = File::create(file_name).expect("Unable to open file");
                BufWriter::new(file)
            };

            for uniquename in self.chromosomes.keys() {
                let chr_config = config.find_chromosome_config(uniquename);
                chr_writers.insert(uniquename, make_chr_gff_writer(&chr_config.export_file_id));
            }

            for gene_details in self.api_maps.genes.values() {
                if let Some(ref gene_loc) = gene_details.location {
                    let chromosome_name = &gene_loc.chromosome_name;
                    let chromosome_export_id =
                        &config.find_chromosome_config(chromosome_name).export_id;
                    let gene_gff_lines =
                        format_gene_gff(chromosome_export_id, &config.database_name,
                                        gene_details);
                    for gff_line in gene_gff_lines {
                        all_gff_writer.write_all(gff_line.as_bytes())?;
                        all_gff_writer.write_all(b"\n")?;

                        match gene_loc.strand {
                            Strand::Forward => {
                                forward_features_gff_writer.write_all(gff_line.as_bytes())?;
                                forward_features_gff_writer.write_all(b"\n")?;
                            },
                            Strand::Reverse => {
                                reverse_features_gff_writer.write_all(gff_line.as_bytes())?;
                                reverse_features_gff_writer.write_all(b"\n")?;
                            }
                            Strand::Unstranded => {
                                unstranded_features_gff_writer.write_all(gff_line.as_bytes())?;
                                unstranded_features_gff_writer.write_all(b"\n")?;
                            }
                        }

                        if let Some(ref mut writer) = chr_writers.get_mut(chromosome_name) {
                            writer.write_all(gff_line.as_bytes())?;
                            writer.write_all(b"\n")?;
                        }
                    }
                }
            }

            for feature_short in self.api_maps.other_features.values() {
                let chromosome_name = &feature_short.location.chromosome_name;
                let chromosome_export_id =
                    &config.find_chromosome_config(chromosome_name).export_id;
                let gff_lines =
                    format_misc_feature_gff(chromosome_export_id, &config.database_name,
                                            feature_short);
                for gff_line in gff_lines {
                    all_gff_writer.write_all(gff_line.as_bytes())?;
                    all_gff_writer.write_all(b"\n")?;

                    match feature_short.location.strand {
                        Strand::Forward => {
                            forward_features_gff_writer.write_all(gff_line.as_bytes())?;
                            forward_features_gff_writer.write_all(b"\n")?;
                        },
                        Strand::Reverse => {
                            reverse_features_gff_writer.write_all(gff_line.as_bytes())?;
                            reverse_features_gff_writer.write_all(b"\n")?;
                        }
                        Strand::Unstranded => {
                            unstranded_features_gff_writer.write_all(gff_line.as_bytes())?;
                            unstranded_features_gff_writer.write_all(b"\n")?;
                        }
                    }

                    if let Some(ref mut writer) = chr_writers.get_mut(chromosome_name) {
                        writer.write_all(gff_line.as_bytes())?;
                        writer.write_all(b"\n")?;
                    }
                }
            }

            for writer in chr_writers.values_mut() {
                writer.flush().unwrap();
            }
        }

        Ok(())
    }

    pub fn write_macromolecular_complexes(&self, config: &Config, output_dir: &str)
                                          -> Result<(), io::Error>
    {
        let mut complex_data: HashMap<(TermShort, GeneShort, RcString), _> = HashMap::new();

        let no_evidence = RcString::from("NO_EVIDENCE");

        let make_key = |annotation: &OntAnnotation| {
            let evidence = annotation.evidence.clone().unwrap_or_else(|| no_evidence.clone());
            (annotation.term_short.clone(), annotation.genes.iter().next().unwrap().clone(),
             evidence)
        };

        if let Some(ref complexes_config) = config.file_exports.macromolecular_complexes {
            let check_parent_term = |el: &RcString| {
                *el == complexes_config.parent_complex_termid
            };
            'TERM: for annotation in &self.ont_annotations {
                let term_short = &annotation.term_short;
                let termid = &term_short.termid;

                if complexes_config.excluded_terms.contains(termid.as_str()) {
                    continue 'TERM;
                }
                if !term_short.interesting_parent_ids.iter().any(check_parent_term) {
                    continue 'TERM;
                }

                let key: (TermShort, GeneShort, RcString) = make_key(annotation);
                complex_data.entry(key)
                    .or_insert_with(Vec::new)
                    .push((annotation.reference_short.clone(), annotation.assigned_by.clone()));
            }
        }

        let complexes_file_name = format!("{}/Complex_annotation.tsv", output_dir);
        let complexes_file = File::create(complexes_file_name).expect("Unable to open file");
        let mut complexes_writer = BufWriter::new(&complexes_file);

        let header = "acc\tGO_name\tsystematic_id\tsymbol\tgene_product_description\tevidence_code\tsource\tassigned_by";

        complexes_writer.write_all(header.as_bytes())?;
        complexes_writer.write_all(b"\n")?;

        let mut lines = vec![];

        for (key, values) in complex_data.drain() {
            let (term_short, gene_short, evidence) = key;
            let mut refs = HashSet::new();
            let mut assigned_bys = HashSet::new();
            for (maybe_ref_short, maybe_assigned_by) in values {
                if let Some(ref_short) = maybe_ref_short {
                    refs.insert(ref_short.uniquename);
                }
                if let Some(assigned_by) = maybe_assigned_by {
                    assigned_bys.insert(assigned_by);
                }
            }

            let mut refs_vec = refs.into_iter().collect::<Vec<_>>();
            refs_vec.sort();
            let mut assigned_bys_vec = assigned_bys.into_iter().collect::<Vec<_>>();
            assigned_bys_vec.sort();

            let refs_string = refs_vec.join(",");
            let assigned_by_string = assigned_bys_vec.join(",");

            let line_bits = vec![term_short.termid.as_str(), term_short.name.as_str(),
                                 gene_short.uniquename.as_str(),
                                 gene_short.name.as_ref().map(RcString::as_str)
                                   .unwrap_or_else(|| gene_short.uniquename.as_str()),
                                 gene_short.product.as_ref().map(RcString::as_str).unwrap_or_else(|| ""),
                                 evidence.as_str(), refs_string.as_str(),
                                 assigned_by_string.as_str()];

            lines.push(line_bits.join("\t"));
        }

        lines.sort();

        for line in lines.drain(0..) {
            complexes_writer.write_all(line.as_bytes())?;
            complexes_writer.write_all(b"\n")?
        }

        Ok(())
    }

    fn write_rnacentral(&self, config: &Config, output_dir: &str) -> Result<(), io::Error> {
        if config.file_exports.rnacentral.is_some() {
            let rnacentral_file_name = format!("{}/rnacentral.json", output_dir);
            let rnacentral_file = File::create(rnacentral_file_name).expect("Unable to open file");
            let mut rnacentral_writer = BufWriter::new(&rnacentral_file);
            let rnacentral_struct = make_rnacentral_struct(config, &self.api_maps.genes);
            let s = serde_json::to_string(&rnacentral_struct).unwrap();

            rnacentral_writer.write_all(s.as_bytes())?;
            rnacentral_writer.write_all(b"\n")?;

            Ok(())
        } else {
            Ok(())
        }
    }

    pub fn write_deletion_viability(&self, config: &Config, output_dir: &str)
                                    -> Result<(), io::Error>
    {
        let deletion_viability_file_name = output_dir.to_owned() + "/FYPOviability.tsv";
        let deletion_viability_file =
            File::create(deletion_viability_file_name).expect("Unable to open file");
        let mut deletion_viability_writer = BufWriter::new(&deletion_viability_file);

        for gene_details in self.api_maps.genes.values() {
            if let Some(load_org_taxonid) = config.load_organism_taxonid {
                if gene_details.taxonid != load_org_taxonid {
                    continue;
                }
            }

            let line = format!("{}\t{}\n",
                               gene_details.uniquename,
                               match gene_details.deletion_viability {
                                   DeletionViability::Viable => "viable",
                                   DeletionViability::Inviable => "inviable",
                                   DeletionViability::DependsOnConditions =>
                                       "condition-dependent",
                                   DeletionViability::Unknown => "unknown",
                               });

            deletion_viability_writer.write_all(line.as_bytes())?;
        }

        deletion_viability_writer.flush()?;

        Ok(())
    }

    pub fn write_slim_ids_and_names(&self, config: &Config, output_dir: &str)
                                       -> Result<(), io::Error> {
        for (slim_name, slim_config) in &config.slims {
            let slim_file_name = format!("{}/{}_ids_and_names.tsv", output_dir, slim_name);

            let slim_file = File::create(slim_file_name).expect("Unable to open file");
            let mut slim_writer = BufWriter::new(&slim_file);

            for term_and_name in &slim_config.terms {
                let line = format!("{}\t{}\n", term_and_name.termid, term_and_name.name);

                slim_writer.write_all(line.as_bytes())?;
            }
        }

        Ok(())
    }

    pub fn write_transmembrane_domains(&self, config: &Config, output_dir: &str)
                                       -> Result<(), io::Error> {
        let tm_domain_file_name =
            output_dir.to_owned() + "/transmembrane_domain_coords_and_seqs.tsv";
        let tm_domain_file =
            File::create(tm_domain_file_name).expect("Unable to open file");
        let mut tm_domain_writer = BufWriter::new(&tm_domain_file);

        let coords_and_seqs = |coords: &[(usize, usize)], prot_seq: &str| {
            let mut coords_strings = vec![];
            let mut seqs = vec![];
            for (start, end) in coords {
                if prot_seq.len() >= *end {
                    coords_strings.push(format!("{}..{}", start, end));
                    let seq = &prot_seq[start-1..*end];
                    seqs.push(seq.to_owned());
                } else {
                    eprintln!("TM domain coord is outside of protein sequence, {}..{} for {}",
                             start, end, prot_seq);
                }
            }
            (coords_strings.join(","), seqs.join(","))
        };

        let star_re = Regex::new(r"\*$").unwrap();

        let format_one_gene = |gene_details: &GeneDetails, prot_seq: &str| {
            let prot_seq = star_re.replace_all(prot_seq, "");
            let (coords, seqs) =
                coords_and_seqs(&gene_details.tm_domain_coords, &prot_seq);
            format!("{}\t{}\t{}\t{}\t{}\n",
                    gene_details.uniquename,
                    gene_details.name.as_ref().map(|n| n.as_str()).unwrap_or(""),
                    prot_seq,
                    coords, seqs)
        };

        for gene_details in self.api_maps.genes.values() {
            if let Some(load_org_taxonid) = config.load_organism_taxonid {
                if gene_details.taxonid != load_org_taxonid {
                    continue;
                }
            }

            if gene_details.tm_domain_coords.is_empty() {
                continue;
            }

            if let Some(transcript) = gene_details.transcripts.get(0) {
                if let Some(ref protein) = transcript.protein {
                    let line = format_one_gene(gene_details, &protein.sequence);

                    tm_domain_writer.write_all(line.as_bytes())?;
                }
            }
        }

        tm_domain_writer.flush()?;

        Ok(())
    }

    fn write_gene_expression_table(&self, output_dir: &str) -> Result<(), io::Error> {
        let file_name = format!("{}/gene_expression_table.tsv", output_dir);
        let file = File::create(file_name).expect("Unable to open file for writing");
        let mut writer = BufWriter::new(&file);

        let header =
            "reference\tfirst_author\tgene\tscale\tterm_name\tterm_id\tduring_term_id\tduring_term_name\tcopies_per_cell\taverage_copies_per_cell\n";
        writer.write_all(header.as_bytes())?;

        for (gene_uniquename, datasets) in &self.api_maps.gene_expression_measurements {

            for measurement in datasets.values() {

                let ref_uniquename = &measurement.reference_uniquename;

                let ref_authors_abbrev =
                    if let Some(ref_details) = self.api_maps.references.get(ref_uniquename) {
                        if let Some(ref authors_abbrev) = ref_details.authors_abbrev {
                            authors_abbrev.split(char::is_whitespace).next().unwrap_or("NA")
                        } else {
                            "Author unknown"
                        }
                    } else {
                        "Author unknown"
                    };


                let termid = &measurement.level_type_termid;

                let term_name =
                    if let Some(term_details) = self.api_maps.terms.get(termid) {
                        term_details.name.as_str()
                    } else {
                        panic!("can't find term details for {}", termid);
                    };


                let during_termid = &measurement.during_termid;

                let during_term_name =
                    if let Some(term_details) = self.api_maps.terms.get(during_termid) {
                        term_details.name.as_str()
                    } else {
                        panic!("can't find term details for {}", during_termid);
                    };

                let scale = &measurement.scale;

                let copies_per_cell =
                    if let Some(ref copies_per_cell) = measurement.copies_per_cell {
                        copies_per_cell
                    } else {
                        "NA"
                    };

                let avg_copies_per_cell =
                    if let Some(ref avg_copies_per_cell) = measurement.avg_copies_per_cell {
                        avg_copies_per_cell
                    } else {
                        "NA"
                    };

                let line =
                    format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            ref_uniquename,
                            ref_authors_abbrev,
                            gene_uniquename,
                            scale,
                            term_name,
                            termid,
                            during_termid,
                            during_term_name,
                            copies_per_cell,
                            avg_copies_per_cell);
                writer.write_all(line.as_bytes())?
            }
        }

        Ok(())
    }

    fn write_site_map_txt(&self, config: &Config, doc_config: &DocConfig, output_dir: &str)
                          -> Result<(), io::Error>
    {
        let base_url = &config.base_url;

        let mut s = format!("{}\n", base_url);

        for page_name in doc_config.pages.keys() {
            if !page_name.ends_with("/index") {
                s += &format!("{}/{}\n", base_url, page_name);
            }
        }

        for gene_details in self.api_maps.genes.values() {
            if let Some(load_org_taxonid) = config.load_organism_taxonid {
                if gene_details.taxonid != load_org_taxonid {
                    continue;
                }
            }
            s += &format!("{}/gene/{}\n", base_url, gene_details.uniquename);
        }

        for term_details in self.api_maps.terms.values() {
            if term_details.gene_count == 0 && term_details.genotype_count == 0 {
                continue;
            }

            let mut found = false;

            for prefix in &config.file_exports.site_map_term_prefixes {
                if term_details.termid.starts_with(prefix.as_str()) {
                    found = true;
                    break;
                }
            }

            if !found {
                continue;
            }

            s += &format!("{}/term/{}\n", base_url, term_details.termid);
        }

        for ref_details in self.api_maps.references.values() {
            if ref_details.cv_annotations.is_empty() {
                continue;
            }

            let mut found = false;

            for prefix in &config.file_exports.site_map_reference_prefixes {
                if ref_details.uniquename.starts_with(prefix.as_str()) {
                    found = true;
                    break;
                }
            }

            if !found {
                continue;
            }


            s += &format!("{}/reference/{}\n", base_url, ref_details.uniquename);
        }

        let file_name = format!("{}/sitemap.txt", output_dir);
        let f = File::create(file_name).expect("Unable to open file");

        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write sitemap.xml");

        Ok(())
    }

    // write the subsets configured using "subset_export"
    fn write_annotation_subsets(&self, config: &Config, output_dir: &str)
                                -> Result<(), io::Error>
    {
        for subset_config in &config.file_exports.annotation_subsets {
            self.write_annotation_subset(&config.cv_config,
                                         subset_config, output_dir)?;
        }

        Ok(())
    }

    fn write_annotation_subset(&self,
                               cv_config_map: &HashMap<CvName, CvConfig>,
                               subset_config: &AnnotationSubsetConfig,
                               output_dir: &str)
                               -> Result<(), io::Error>
    {
        let file_name = format!("{}/{}", output_dir, subset_config.file_name);

        let file = File::create(file_name).expect("Unable to open file for writing");

        let mut writer = BufWriter::new(&file);

        let table = table_for_export(&self.api_maps, cv_config_map, subset_config);

        for row in table {
            let line = row.join("\t") + "\n";
            writer.write_all(line.as_bytes())?;
        }

        Ok(())
    }


    pub fn write_stats(&self, output_dir: &str) -> Result<(), io::Error> {
        let s = serde_json::to_string(&self.stats).unwrap();
        let file_name = String::new() + output_dir + "/stats.json";
        let f = File::create(file_name).expect("Unable to open file");
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes()).expect("Unable to write stats.json");

        Ok(())
    }

    pub fn write(&self, config: &Config, go_eco_mappping: &GoEcoMapping,
                 doc_config: &DocConfig, output_dir: &str)
                 -> Result<(), io::Error>
    {
        let web_json_path = self.create_dir(output_dir, "web-json");

        self.write_chromosome_json(config, &web_json_path);
        println!("wrote {} chromosomes", self.get_chromosomes().len());
        self.write_gene_summaries(&web_json_path);
        self.write_chromosome_summaries(&web_json_path);
        println!("wrote summaries");
        self.write_metadata(&web_json_path);
        println!("wrote metadata");
        self.write_recent_references(&web_json_path);
        self.write_all_community_curated(&web_json_path);
        self.write_all_admin_curated(&web_json_path);
        println!("wrote references");
        self.write_api_maps(&web_json_path);
        self.write_solr_data(&web_json_path);
        println!("wrote search data");

        let intermine_data_path = self.create_dir(output_dir, "intermine_data");
        self.write_intermine_data(config, &intermine_data_path)?;
        println!("wrote intermine data");

        self.write_subsets(&web_json_path);
        println!("wrote subsets");

        let fasta_path = self.create_dir(output_dir, "fasta");
        let feature_sequences_path = self.create_dir(&fasta_path, "feature_sequences");
        self.write_feature_sequences(&feature_sequences_path);
        let chromosomes_path = self.create_dir(&fasta_path, "chromosomes");
        self.write_chromosome_sequences(config, &chromosomes_path);
        println!("wrote fasta");

        let misc_path = self.create_dir(output_dir, "misc");

        self.write_go_annotation_files(config, go_eco_mappping, &misc_path)?;

        self.write_gene_id_table(config, &misc_path)?;
        self.write_protein_features(config, &misc_path)?;
        self.write_feature_coords(config, &misc_path)?;
        self.write_macromolecular_complexes(config, &misc_path)?;
        self.write_rnacentral(config, &misc_path)?;
        self.write_deletion_viability(config, &misc_path)?;
        self.write_slim_ids_and_names(config, &misc_path)?;
        self.write_transmembrane_domains(config, &misc_path)?;
        self.write_gene_expression_table(&misc_path)?;
        self.write_site_map_txt(config, doc_config, &misc_path)?;

        self.write_annotation_subsets(config, &misc_path)?;

        self.write_stats(&web_json_path)?;

        let gff_path = self.create_dir(output_dir, "gff");
        self.write_gff(config, &gff_path)?;

        Ok(())
    }
}


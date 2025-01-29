extern crate tokio_postgres;
extern crate serde_json;

use std::cmp::min;
use std::fs::{File, create_dir_all};
use std::io::{Write, BufWriter};
use std::io;
use std::collections::HashMap;
use std::sync::{Arc, RwLock};
use regex::Regex;

use flate2::Compression;
use flate2::write::GzEncoder;
use zstd::stream::Encoder;

use rusqlite::Connection;

use flexstr::{SharedStr as FlexStr, shared_str as flex_str, ToSharedStr, shared_fmt as flex_fmt};

use crate::bio::util::{format_fasta, format_gene_gff, format_misc_feature_gff,
                       process_modification_ext};

use crate::constants::*;

use crate::web::config::*;
use crate::rnacentral::*;

use crate::types::{CvName, TermId, GenotypeDisplayUniquename, GeneUniquename, AlleleUniquename,
                   ReferenceUniquename};
use crate::data_types::*;
use crate::annotation_util::table_for_export;

use crate::bio::go_format_writer::write_go_annotation_files;
use crate::bio::phenotype_format_writer::write_phenotype_annotation_files;
use crate::bio::macromolecular_complexes::write_macromolecular_complexes;
use crate::bio::gene_expression_writer::write_gene_expression_row;

use crate::utils::{join, make_maps_database_tables, store_maps_into_database};

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct WebData {
    pub metadata: Metadata,
    pub chromosomes: ChrNameDetailsMap,
    pub chromosome_summaries: Vec<ChromosomeShort>,
    pub recent_references: RecentReferences,
    pub all_community_curated: Vec<ReferenceShort>,
    pub all_admin_curated: Vec<ReferenceShort>,
    pub api_maps: APIMaps,
    pub genes: UniquenameGeneMap,
    pub alleles: UniquenameAlleleMap,
    pub genotypes: IdGenotypeMap,
    pub terms: TermIdDetailsMap,
    pub references: UniquenameReferenceMap,
    pub annotation_details: IdOntAnnotationDetailMap,
    pub termid_genotype_annotation: HashMap<TermId, Vec<APIGenotypeAnnotation>>,

    pub solr_data: SolrData,
    pub search_gene_summaries: Vec<GeneSummary>,
    pub ont_annotations: Vec<OntAnnotation>,
    pub stats: Stats,
    pub detailed_stats: DetailedStats,

    pub physical_interaction_annotations: Vec<InteractionAnnotation>,
    pub genetic_interaction_annotations: Vec<InteractionAnnotation>,

    pub arc_terms: Arc<RwLock<HashMap<TermId, Arc<TermDetails>>>>,
    pub arc_genes: Arc<RwLock<HashMap<GeneUniquename, Arc<GeneDetails>>>>,
    pub arc_alleles: Arc<RwLock<HashMap<AlleleUniquename, Arc<AlleleDetails>>>>,
    pub arc_references: Arc<RwLock<HashMap<ReferenceUniquename, Arc<ReferenceDetails>>>>,
    pub arc_genotypes: Arc<RwLock<HashMap<GenotypeDisplayUniquename, Arc<GenotypeDetails>>>>,
    pub arc_annotation_details: Arc<RwLock<HashMap<OntAnnotationId, Arc<OntAnnotationDetail>>>>,
}

impl DataLookup for WebData {
    fn get_term(&self, termid: &TermId) -> Option<Arc<TermDetails>> {
        let mut arc_terms = self.arc_terms.write().unwrap();
        if let Some(arc_term_details) = arc_terms.get(termid) {
            Some(arc_term_details.to_owned())
        } else {
           if let Some(term_details) = self.terms.get(termid) {
               let arc_term_details = Arc::new(term_details.to_owned());
               arc_terms.insert(termid.to_owned(), arc_term_details.clone());
               Some(arc_term_details)
           } else {
               None
           }
        }
    }

    fn get_gene(&self, gene_uniquename: &GeneUniquename) -> Option<Arc<GeneDetails>> {
        let mut arc_genes = self.arc_genes.write().unwrap();
        if let Some(arc_gene_details) = arc_genes.get(gene_uniquename) {
            Some(arc_gene_details.to_owned())
        } else {
           if let Some(gene_details) = self.genes.get(gene_uniquename) {
               let arc_gene_details = Arc::new(gene_details.to_owned());
               arc_genes.insert(gene_uniquename.to_owned(), arc_gene_details.clone());
               Some(arc_gene_details)
           } else {
               None
           }
        }
    }

    fn get_allele(&self, allele_uniquename: &AlleleUniquename) -> Option<Arc<AlleleDetails>> {
        let mut arc_alleles = self.arc_alleles.write().unwrap();
        if let Some(arc_allele_details) = arc_alleles.get(allele_uniquename) {
            Some(arc_allele_details.to_owned())
        } else {
           if let Some(allele_details) = self.alleles.get(allele_uniquename) {
               let arc_allele_details = Arc::new(allele_details.to_owned());
               arc_alleles.insert(allele_uniquename.to_owned(), arc_allele_details.clone());
               Some(arc_allele_details)
           } else {
               None
           }
        }
    }

    fn get_reference(&self, reference_uniquename: &ReferenceUniquename) -> Option<Arc<ReferenceDetails>> {
        let mut arc_references = self.arc_references.write().unwrap();
        if let Some(arc_reference_details) = arc_references.get(reference_uniquename) {
            Some(arc_reference_details.to_owned())
        } else {
           if let Some(reference_details) = self.references.get(reference_uniquename) {
               let arc_reference_details = Arc::new(reference_details.to_owned());
               arc_references.insert(reference_uniquename.to_owned(), arc_reference_details.clone());
               Some(arc_reference_details)
           } else {
               None
           }
        }
    }

    fn get_genotype(&self, genotype_display_uniquename: &GenotypeDisplayUniquename)
           -> Option<Arc<GenotypeDetails>>
    {
        let mut arc_genotypes = self.arc_genotypes.write().unwrap();
        if let Some(arc_genotype_details) = arc_genotypes.get(genotype_display_uniquename) {
            Some(arc_genotype_details.to_owned())
        } else {
           if let Some(genotype_details) = self.genotypes.get(genotype_display_uniquename) {
               let arc_genotype_details = Arc::new(genotype_details.to_owned());
               arc_genotypes.insert(genotype_display_uniquename.to_owned(),
                                         arc_genotype_details.clone());
               Some(arc_genotype_details)
           } else {
               None
           }
        }
    }

    fn get_annotation_detail(&self, annotation_id: OntAnnotationId)
           -> Option<Arc<OntAnnotationDetail>>
    {
        let mut arc_annotation_details = self.arc_annotation_details.write().unwrap();
        if let Some(arc_annotation_detail_details) = arc_annotation_details.get(&annotation_id) {
            Some(arc_annotation_detail_details.to_owned())
        } else {
           if let Some(annotation_detail_details) = self.annotation_details.get(&annotation_id) {
               let arc_annotation_detail_details = Arc::new(annotation_detail_details.to_owned());
               arc_annotation_details.insert(annotation_id,
                                             arc_annotation_detail_details.clone());
               Some(arc_annotation_detail_details)
           } else {
               None
           }
        }
    }
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

    fn create_dir(&self, output_dir: &str, dir_name: &str) -> Result<String, io::Error> {
        let path = String::new() + output_dir + "/" + dir_name;
        create_dir_all(&path)?;
        Ok(path)
    }

    fn write_chromosome_seq_chunks(&self, output_dir: &str, chunk_sizes: &[usize]) -> Result<(), io::Error> {
        for chunk_size in chunk_sizes {
            for (chromosome_uniquename, chromosome_details) in &self.chromosomes {
                let new_path_part = &format!("{}/sequence/{}", chromosome_uniquename, chunk_size);
                let chr_path = self.create_dir(output_dir, new_path_part)?;
                let mut index = 0;
                let max_index = chromosome_details.residues.len() / chunk_size;
                while index <= max_index {
                    let start_pos = index*chunk_size;
                    let end_pos = min(start_pos+chunk_size, chromosome_details.residues.len());
                    let chunk: String = chromosome_details.residues[start_pos..end_pos].into();
                    let file_name = format!("{}/chunk_{}", chr_path, index);
                    let f = File::create(file_name)?;
                    let mut writer = BufWriter::new(&f);
                    writer.write_all(chunk.as_bytes())?;
                    index += 1;
                }
            }
        }
        Ok(())
    }

    fn write_chromosome_json(&self, config: &Config, output_dir: &str) -> Result<(), io::Error> {
        let new_path = self.create_dir(output_dir, "chromosome")?;
        for (chromosome_uniquename, chromosome_details) in &self.chromosomes {
            let s = serde_json::to_string(&chromosome_details).unwrap();
            let file_name = format!("{}/{}.json", new_path, &chromosome_uniquename);
            let f = File::create(file_name)?;
            let mut writer = BufWriter::new(&f);
            writer.write_all(s.as_bytes())?
        }
        self.write_chromosome_seq_chunks(&new_path, &config.api_seq_chunk_sizes)?;
        Ok(())
    }

    fn write_gene_summaries(&self, output_dir: &str) -> Result<(), io::Error> {
        let s = serde_json::to_string(&self.search_gene_summaries).unwrap();
        let file_name = String::new() + output_dir + "/gene_summaries.json";
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes())?;
        Ok(())
    }

    fn write_metadata(&self, output_dir: &str) -> Result<(), io::Error> {
        let s = serde_json::to_string(&self.metadata).unwrap();
        let file_name = String::new() + output_dir + "/metadata.json";
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes())?;
        Ok(())
    }

    fn write_recent_references(&self, output_dir: &str) -> Result<(), io::Error> {
        let s = serde_json::to_string(&self.recent_references).unwrap();
        let file_name = String::new() + output_dir + "/recent_references.json";
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes())?;
        Ok(())
    }

    fn write_all_community_curated(&self, output_dir: &str) -> Result<(), io::Error> {
        let s = serde_json::to_string(&self.all_community_curated).unwrap();
        let file_name = String::new() + output_dir + "/community_curated_references.json";
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes())?;
        Ok(())
    }

    fn write_all_admin_curated(&self, output_dir: &str) -> Result<(), io::Error> {
        let s = serde_json::to_string(&self.all_admin_curated).unwrap();
        let file_name = String::new() + output_dir + "/admin_curated_references.json";
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes())?;
        Ok(())
    }

    fn write_api_maps(&self, output_dir: &str) -> Result<(), io::Error> {
        let file_name = String::new() + output_dir + "/api_maps.json.zst";
        let f = File::create(file_name)?;

        let mut compressor = Encoder::new(f, 12)?;
        compressor.multithread(8).unwrap();
        compressor.long_distance_matching(true)?;
        serde_json::ser::to_writer(&mut compressor, &self.api_maps).unwrap();

        compressor.finish()?;

        Ok(())
    }

    fn write_solr_data(&self, output_dir: &str) -> Result<(), io::Error> {
        let new_path = self.create_dir(output_dir, "solr_data/")?;

        let terms = self.solr_data.term_summaries.clone();

        let terms_json_text = serde_json::to_string(&terms).unwrap();
        let terms_file_name = format!("{}/terms.json.gz", new_path);
        let terms_file = File::create(terms_file_name)?;

        let mut terms_compressor = GzEncoder::new(terms_file, Compression::default());
        terms_compressor.write_all(terms_json_text.as_bytes())?;
        terms_compressor.finish()?;

        let genes = self.solr_data.gene_summaries.clone();

        let genes_json_text = serde_json::to_string(&genes).unwrap();
        let genes_file_name = format!("{}/genes.json.gz", new_path);
        let genes_file = File::create(genes_file_name)?;

        let mut genes_compressor = GzEncoder::new(genes_file, Compression::default());
        genes_compressor.write_all(genes_json_text.as_bytes())?;
        genes_compressor.finish()?;


        let alleles = self.solr_data.allele_summaries.clone();

        let alleles_json_text = serde_json::to_string(&alleles).unwrap();
        let alleles_file_name = format!("{}/alleles.json.gz", new_path);
        let alleles_file = File::create(alleles_file_name)?;

        let mut alleles_compressor = GzEncoder::new(alleles_file, Compression::default());
        alleles_compressor.write_all(alleles_json_text.as_bytes())?;
        alleles_compressor.finish()?;

        let references = self.solr_data.reference_summaries.clone();

        let references_json_text = serde_json::to_string(&references).unwrap();
        let references_file_name = format!("{}/references.json.gz", new_path);
        let references_file = File::create(references_file_name)?;

        let mut references_compressor = GzEncoder::new(references_file, Compression::default());
        references_compressor.write_all(references_json_text.as_bytes())?;
        references_compressor.finish()?;
        Ok(())
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
            self.genes.values()
            .filter(|gene_details| gene_details.taxonid == load_org_taxonid)
            .map(|gene_details| {
                InterMineGeneDetails::from_gene_details(&self.api_maps.transcripts,
                                                        gene_details)
            }).collect();

        let genes_json_text = serde_json::to_string(&intermine_genes)?;
        let genes_file_name = format!("{}/pombemine_gene_details.gz", output_dir);
        let genes_file = File::create(genes_file_name)?;

        let mut genes_compressor = GzEncoder::new(genes_file, Compression::default());
        genes_compressor.write_all(genes_json_text.as_bytes())?;
        genes_compressor.finish()?;

        Ok(())
    }

    fn write_subsets(&self, output_dir: &str) -> Result<(), io::Error> {
        let s = serde_json::to_string(&self.api_maps.term_subsets).unwrap();
        let file_name = String::new() + output_dir + "/term_subsets.json";
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes())?;

        let s = serde_json::to_string(&self.api_maps.gene_subsets).unwrap();
        let file_name = String::new() + output_dir + "/gene_subsets.json";
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes())?;
        Ok(())
    }

    fn write_feature_sequences(&self, output_dir: &str) -> Result<(), io::Error> {
        let make_seq_writer = |name: &str| -> Result<BufWriter<File>, io::Error> {
            let file_name = String::new() + output_dir + "/" + name;
            let file = File::create(file_name)?;
            Ok(BufWriter::new(file))
        };

        let mut cds_writer = make_seq_writer("cds.fa")?;
        let mut cds_introns_writer = make_seq_writer("cds+introns.fa")?;
        let mut cds_introns_utrs_writer = make_seq_writer("cds+introns+utrs.fa")?;
        let mut introns_writer = make_seq_writer("introns_within_cds.fa")?;
        let mut five_prime_utrs_writer = make_seq_writer("five_prime_utrs.fa")?;
        let mut three_prime_utrs_writer = make_seq_writer("three_prime_utrs.fa")?;
        let mut peptide_writer = make_seq_writer("peptide.fa")?;

        for (gene_uniquename, gene_details) in &self.genes {
            if let Some(transcript_uniquename) =
                gene_details.transcripts.get(0)
            {
                let transcript = self.api_maps.transcripts.get(transcript_uniquename)
                    .expect(&format!("internal error, can't find transcript details for {}",
                                     transcript_uniquename));
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
                                       Some(gene_uniquename.to_string()),
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

        cds_writer.flush()?;
        cds_introns_writer.flush()?;
        cds_introns_utrs_writer.flush()?;
        introns_writer.flush()?;
        peptide_writer.flush()?;
        five_prime_utrs_writer.flush()?;
        three_prime_utrs_writer.flush()?;
        Ok(())
    }

    pub fn write_chromosome_sequences(&self, config: &Config, output_dir: &str) -> Result<(), io::Error> {
        let make_seq_writer = |name: &str| -> Result<BufWriter<File>, io::Error> {
            let file_name = String::new() + output_dir + "/" + name;
            let file = File::create(file_name)?;
            Ok(BufWriter::new(file))
        };

        if let Some(load_org) = config.load_organism() {
            let load_org_name = load_org.full_name();
            let chromosomes_file_name = load_org_name.clone() + "_all_chromosomes.fa";
            let mut chromosomes_writer = make_seq_writer(&chromosomes_file_name)?;

            for (uniquename, details) in &self.chromosomes {
                let chr_config = config.find_chromosome_config(uniquename);
                write_as_fasta(&mut chromosomes_writer, &chr_config.export_id,
                               Some(load_org_name.clone()), &details.residues);
                let this_chr_file_name =
                    load_org_name.clone() + "_" + &chr_config.export_file_id + ".fa";
                let mut this_chr_writer = make_seq_writer(&this_chr_file_name)?;
                write_as_fasta(&mut this_chr_writer, &chr_config.export_id,
                               Some(load_org_name.clone()), &details.residues);
                this_chr_writer.flush().unwrap();

            }

            chromosomes_writer.flush().unwrap();
        }
        Ok(())
    }

    fn write_chromosome_summaries(&self, output_dir: &str) -> Result<(), io::Error> {
        let s = serde_json::to_string(&self.chromosome_summaries).unwrap();
        let file_name = String::new() + output_dir + "/chromosome_summaries.json";
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes())?;
        Ok(())
    }

    fn write_alleles_json(&self, output_dir: &str) -> Result<(), io::Error> {
        let allele_summaries: AlleleShortMap =
            self.alleles.iter().map(|(uniquename, details)| {
                (uniquename.clone(), details.into())
            }).collect();

        let s = serde_json::to_string(&allele_summaries).unwrap();
        let file_name = String::new() + output_dir + "/allele_summaries.json";
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(&f);

        writer.write_all(s.as_bytes())?;

        Ok(())
    }

    fn write_gene_id_table(&self, config: &Config, output_dir: &str) -> Result<(), io::Error> {
        let gene_file_name = output_dir.to_owned() + "/sysID2product.tsv";
        let rna_file_name = output_dir.to_owned() + "/sysID2product.rna.tsv";
        let pseudogenes_file_name = output_dir.to_owned() + "/pseudogeneIDs.tsv";
        let all_names_file_name = output_dir.to_owned() + "/gene_IDs_names.tsv";
        let all_ids_file_name = output_dir.to_owned() + "/gene_IDs_names_products.tsv";
        let uniprot_ids_file_name = output_dir.to_owned() + "/uniprot_id_mapping.tsv";

        let gene_file = File::create(gene_file_name)?;
        let rna_file = File::create(rna_file_name)?;
        let pseudogenes_file = File::create(pseudogenes_file_name)?;
        let all_names_file = File::create(all_names_file_name)?;
        let all_ids_file = File::create(all_ids_file_name)?;
        let uniprot_ids_file = File::create(uniprot_ids_file_name)?;

        let mut gene_writer = BufWriter::new(&gene_file);
        let mut rna_writer = BufWriter::new(&rna_file);
        let mut pseudogenes_writer = BufWriter::new(&pseudogenes_file);
        let mut all_names_writer = BufWriter::new(&all_names_file);
        let mut all_ids_writer = BufWriter::new(&all_ids_file);
        let mut uniprot_ids_writer = BufWriter::new(&uniprot_ids_file);

        let db_version = format!("# Chado database date: {}\n", self.metadata.db_creation_datetime);
        gene_writer.write_all(db_version.as_bytes())?;
        rna_writer.write_all(db_version.as_bytes())?;
        pseudogenes_writer.write_all(db_version.as_bytes())?;
        all_names_writer.write_all(db_version.as_ref())?;

        let small_header =
            "gene_systematic_id\tgene_name\tsynonyms\n";
        write!(pseudogenes_writer, "{}", small_header)?;
        write!(all_names_writer, "{}", small_header)?;

        let header_with_product =
            "gene_systematic_id\tgene_name\tsynonyms\tgene_product\n";
        write!(gene_writer, "{}", header_with_product)?;
        write!(rna_writer, "{}", header_with_product)?;

        let big_header =
            "gene_systematic_id\tgene_systematic_id_with_prefix\tgene_name\tchromosome_id\tgene_product\tuniprot_id\tgene_type\tsynonyms\n";
        write!(all_ids_writer, "{}", big_header)?;

        for gene_details in self.genes.values() {
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
                               gene_details.name.clone().unwrap_or(flex_str!("")),
                               synonyms);

            let gene_name = if let Some(ref gene_details_name) = gene_details.name {
                gene_details_name.clone()
            } else {
                flex_str!("")
            };

            let gene_product = if let Some(ref gene_details_product) = gene_details.product {
                gene_details_product.clone()
            } else {
                flex_str!("")
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
                if let Some(gene_uniprot_id) = gene_details.uniprot_identifier.as_ref() {
                    gene_uniprot_id
                } else {
                    ""
                };


            if uniprot_id.len() > 0 && gene_details.feature_type == "mRNA gene" {
                let gene_name_or_dash =
                    if gene_name.len() == 0 {
                        "-"
                    } else {
                        gene_name.as_str()
                    };
                let uniprot_ids_line = format!("{}\t{}\t{}\n",
                                               uniprot_id,
                                               gene_details.uniquename,
                                               gene_name_or_dash);
                uniprot_ids_writer.write_all(uniprot_ids_line.as_bytes())?;
            }

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
        all_ids_writer.flush()?;
        uniprot_ids_writer.flush()?;

        Ok(())
    }

    fn write_protein_features(&self, config: &Config, output_dir: &str)
                              -> Result<(), io::Error>
    {
        let peptide_stats_name = format!("{}/PeptideStats.tsv", output_dir);
        let peptide_stats_file = File::create(peptide_stats_name)?;
        let mut peptide_stats_writer = BufWriter::new(&peptide_stats_file);

        let peptide_stats_header = "Systematic_ID\tMass (kDa)\tpI\tCharge\tResidues\tCAI\n";
        peptide_stats_writer.write_all(peptide_stats_header.as_bytes())?;

        let protein_features_name = format!("{}/ProteinFeatures.tsv", output_dir);
        let protein_features_file = File::create(protein_features_name)?;
        let mut protein_features_writer = BufWriter::new(&protein_features_file);

        let disordered_regions_name = format!("{}/disordered_regions.tsv", output_dir);
        let disordered_regions_file = File::create(disordered_regions_name)?;
        let mut disordered_regions_writer = BufWriter::new(&disordered_regions_file);

        let aa_composition_name = format!("{}/aa_composition.tsv", output_dir);
        let aa_composition_file = File::create(aa_composition_name)?;
        let mut aa_composition_writer = BufWriter::new(&aa_composition_file);

        let protein_features_header =
            "systematic_id\tgene_name\tpeptide_id\tdomain_id\tdatabase\tseq_start\tseq_end\n";
        protein_features_writer.write_all(protein_features_header.as_bytes())?;

        let disordered_regions_header =
            "systematic_id\tgene_name\tpeptide_id\tseq_start\tseq_end\n";
        disordered_regions_writer.write_all(disordered_regions_header.as_bytes())?;

        let db_display_name = |db_alias: &FlexStr| {
            let lower_db_alias = db_alias.to_lowercase().to_shared_str();
            if let Some(name) = config.extra_database_aliases.get(&lower_db_alias) {
                name.clone()
            } else {
                db_alias.clone()
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

        for (gene_uniquename, gene_details) in &self.genes {
            if let Some(transcript_uniquename) =
                gene_details.transcripts.get(0)
            {
                let transcript = self.api_maps.transcripts.get(transcript_uniquename)
                    .expect(&format!("internal error, can't find transcript details for {}",
                                     transcript_uniquename));

                if let Some(ref protein) = transcript.protein {
                    let line = format!("{}\t{:.2}\t{}\t{}\t{}\t{}\n",
                                       gene_uniquename, protein.molecular_weight,
                                       protein.isoelectric_point,
                                       protein.charge_at_ph7,
                                       protein.sequence.len() - 1,
                                       protein.codon_adaptation_index);
                    peptide_stats_writer.write_all(line.as_bytes())?;

                    let gene_name = gene_details.name.clone().unwrap_or_else(|| flex_str!(""));
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
                        let range = &disordered_region.range;
                        let start_pos = range.start;
                        let end_pos = range.end;
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

        let composition_line = |first_col_string: FlexStr, comp: &AAComposition| {
            let mut line = first_col_string.to_string();

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
            composition_line(flex_str!("total"), &total_composition);
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

            let gene_file = File::create(gene_file_name)?;
            let cds_file = File::create(cds_file_name)?;
            let exon_file = File::create(exon_file_name)?;

            let mut gene_writer = BufWriter::new(&gene_file);
            let mut cds_writer = BufWriter::new(&cds_file);
            let mut exon_writer = BufWriter::new(&exon_file);

            for gene_uniquename in &chr_details.gene_uniquenames {
                let gene = &self.genes[gene_uniquename];
                if let Some(ref gene_location) = gene.location {
                    write_line(gene_uniquename, gene_location, &mut gene_writer)?;

                    for transcript_uniquename in &gene.transcripts {
                        let transcript = self.api_maps.transcripts
                            .get(transcript_uniquename)
                            .expect(&format!("internal error, can't find transcript details for {}",
                                             transcript_uniquename));

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
            let all_gff_file = File::create(all_gff_name)?;
            let mut all_gff_writer = BufWriter::new(&all_gff_file);

            let forward_features_gff_name =
                format!("{}/{}_all_chromosomes_forward_strand.gff3", output_dir, load_org_name);
            let forward_features_gff_file = File::create(forward_features_gff_name)?;
            let mut forward_features_gff_writer = BufWriter::new(&forward_features_gff_file);

            let reverse_features_gff_name =
                format!("{}/{}_all_chromosomes_reverse_strand.gff3", output_dir, load_org_name);
            let reverse_features_gff_file = File::create(reverse_features_gff_name)?;
            let mut reverse_features_gff_writer = BufWriter::new(&reverse_features_gff_file);

            let unstranded_features_gff_name =
                format!("{}/{}_all_chromosomes_unstranded.gff3", output_dir, load_org_name);
            let unstranded_features_gff_file = File::create(unstranded_features_gff_name)?;
            let mut unstranded_features_gff_writer = BufWriter::new(&unstranded_features_gff_file);

            all_gff_writer.write_all(b"##gff-version 3\n")?;
            forward_features_gff_writer.write_all(b"##gff-version 3\n")?;
            reverse_features_gff_writer.write_all(b"##gff-version 3\n")?;
            unstranded_features_gff_writer.write_all(b"##gff-version 3\n")?;

            let mut chr_writers = HashMap::new();

            let make_chr_gff_writer = |export_name: &str| -> Result<BufWriter<File>, io::Error> {
                let file_name = String::new() +
                    output_dir + "/" + &load_org_name + "_" + export_name + ".gff3";
                let file = File::create(file_name)?;
                Ok(BufWriter::new(file))
            };

            for uniquename in self.chromosomes.keys() {
                let chr_config = config.find_chromosome_config(uniquename);
                chr_writers.insert(uniquename, make_chr_gff_writer(&chr_config.export_file_id)?);
            }

            for gene_details in self.genes.values() {
                if let Some(ref gene_loc) = gene_details.location {
                    let chromosome_name = &gene_loc.chromosome_name;
                    let chromosome_export_id =
                        &config.find_chromosome_config(chromosome_name).export_id;
                    let gene_gff_lines =
                        format_gene_gff(chromosome_export_id, &config.database_name,
                                        &self.api_maps.transcripts, gene_details);
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

    fn write_rnacentral(&self, config: &Config, output_dir: &str) -> Result<(), io::Error> {
        if config.file_exports.rnacentral.is_some() {
            let rnacentral_file_name = format!("{}/rnacentral.json", output_dir);
            let rnacentral_file = File::create(rnacentral_file_name)?;
            let mut rnacentral_writer = BufWriter::new(&rnacentral_file);
            let rnacentral_struct = make_rnacentral_struct(config, &self.api_maps.transcripts,
                                                           &self.genes);
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
            File::create(deletion_viability_file_name)?;
        let mut deletion_viability_writer = BufWriter::new(&deletion_viability_file);

        for gene_details in self.genes.values() {
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

            let slim_file = File::create(slim_file_name)?;
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
            File::create(tm_domain_file_name)?;
        let mut tm_domain_writer = BufWriter::new(&tm_domain_file);

        let coords_and_seqs = |assigned_ranges: &[AssignedByPeptideRange], prot_seq: &str| {
            let mut coords_strings = vec![];
            let mut seqs = vec![];
            for assigned_range in assigned_ranges {
                let range = &assigned_range.range;
                let start = range.start;
                let end = range.end;
                if prot_seq.len() >= end {
                    coords_strings.push(format!("{}..{}", start, end));
                    let seq = &prot_seq[start-1..end];
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

        for gene_details in self.genes.values() {
            if let Some(load_org_taxonid) = config.load_organism_taxonid {
                if gene_details.taxonid != load_org_taxonid {
                    continue;
                }
            }

            if gene_details.tm_domain_coords.is_empty() {
                continue;
            }

            if let Some(transcript_uniquename) =
                gene_details.transcripts.get(0)
            {
                let transcript = self.api_maps.transcripts.get(transcript_uniquename)
                    .expect(&format!("internal error, can't find transcript details for {}",
                                     transcript_uniquename));

                if let Some(ref protein) = transcript.protein {
                    let line = format_one_gene(gene_details, &protein.sequence);

                    tm_domain_writer.write_all(line.as_bytes())?;
                }
            }
        }

        tm_domain_writer.flush()?;

        Ok(())
    }

    fn write_htp_gene_expression_table(&self, output_dir: &str) -> Result<(), io::Error> {
        let file_name = format!("{}/htp_gene_expression_table.tsv", output_dir);
        let file = File::create(file_name)?;
        let mut writer = BufWriter::new(&file);

        let header =
            "reference\tfirst_author\tgene\tscale\tterm_name\tduring_term_id\tduring_term_name\tcopies_per_cell\taverage_copies_per_cell\n";
        writer.write_all(header.as_bytes())?;

        for (gene_uniquename, datasets) in &self.api_maps.gene_expression_measurements {

            for measurement in datasets.values() {

                let ref_uniquename = &measurement.reference_uniquename;

                let ref_authors_abbrev =
                    if let Some(ref_details) = self.references.get(ref_uniquename) {
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
                    if let Some(term_details) = self.terms.get(termid) {
                        term_details.name.as_str()
                    } else {
                        panic!("can't find term details for {}", termid);
                    };


                let during_termid = &measurement.during_termid;

                let during_term_name =
                    if let Some(term_details) = self.terms.get(during_termid) {
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
                    format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                            ref_uniquename,
                            ref_authors_abbrev,
                            gene_uniquename,
                            scale,
                            term_name,
                            during_termid,
                            during_term_name,
                            copies_per_cell,
                            avg_copies_per_cell);
                writer.write_all(line.as_bytes())?
            }
        }

        Ok(())
    }

    fn write_full_gene_expression_table(&self, output_dir: &str)
           -> Result<(), io::Error>
    {
        let file_name = format!("{}/full_gene_expression_table.tsv", output_dir);
        let file = File::create(file_name)?;
        let mut writer = BufWriter::new(&file);

        let header = "gene_systematic_id\tgene_name\ttype\textension\tcopies_per_cell\trange\tevidence\tscale\tcondition\treference\ttaxon\tdate\n";
        writer.write_all(header.as_bytes())?;

        for gene_details in self.genes.values() {
            let Some(term_annotations) = gene_details.cv_annotations.get("quantitative_gene_expression")
            else {
                continue;
            };

            for term_annotation in term_annotations {
                for annotation_id in &term_annotation.annotations {
                    let annotation_detail = self.get_annotation_detail(*annotation_id)
                        .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));

                    write_gene_expression_row(&mut writer,
                                              &self.terms,
                                              gene_details,
                                              &term_annotation.term,
                                              &annotation_detail)?
                }
            }
        }

        Ok(())
    }

    fn write_site_map_txt(&self, config: &Config, doc_config: &DocConfig,
                          references: &UniquenameReferenceMap, output_dir: &str)
                          -> Result<(), io::Error>
    {
        let base_url = &config.base_url;

        let mut s = format!("{}\n", base_url);

        for page_name in doc_config.pages.keys() {
            if !page_name.ends_with("/index") {
                s += &format!("{}/{}\n", base_url, page_name);
            }
        }

        for gene_details in self.genes.values() {
            if let Some(load_org_taxonid) = config.load_organism_taxonid {
                if gene_details.taxonid != load_org_taxonid {
                    continue;
                }
            }
            s += &format!("{}/gene/{}\n", base_url, gene_details.uniquename);
        }

        for term_details in self.terms.values() {
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

        for ref_details in references.values() {
            if ref_details.cv_annotations.is_empty() &&
               ref_details.physical_interactions.is_empty() &&
               ref_details.genetic_interactions.is_empty() &&
               ref_details.ortholog_annotations.is_empty() &&
               ref_details.paralog_annotations.is_empty() &&
               ref_details.pdb_entries.is_empty() &&
               ref_details.canto_triage_status != Some(flex_str!("Curatable")) {
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
        let f = File::create(file_name)?;

        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes())?;

        Ok(())
    }

    fn write_allele_tsv(&self, output_dir: &str) -> Result<(), io::Error> {
        let file_name = format!("{}/all_alleles.tsv", output_dir);
        let file = File::create(file_name)?;
        let mut writer = BufWriter::new(&file);

        let header = "#gene_systematic_id\tgene_name\tcurrent_internal_id\tallele_name\tallele_type\tallele_description\tsynonyms\n";
        writer.write_all(header.as_bytes())?;

        let empty_string = flex_str!("");

        for allele_short in self.alleles.values() {
            let gene = &allele_short.gene;
            let synonyms = allele_short.synonyms
                .iter().map(|s| s.name.as_str()).collect::<Vec<_>>().join("|");

            let line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                      gene.uniquename,
                                      gene.name.as_ref().unwrap_or(&empty_string),
                                      allele_short.uniquename,
                                      allele_short.name.as_ref().unwrap_or(&empty_string),
                                      allele_short.allele_type,
                                      allele_short.description.as_ref().unwrap_or(&empty_string),
                                      synonyms);

            writer.write_all(line.as_bytes())?;
        }

        Ok(())
    }

    fn write_disease_association(&self, config: &Config, output_dir: &str) -> Result<(), io::Error> {
        let file_name = format!("{}/disease_association.tsv", output_dir);
        let file = File::create(file_name)?;
        let mut writer = BufWriter::new(&file);

        let load_org_taxonid =
            if let Some(load_org_taxonid) = config.load_organism_taxonid {
                load_org_taxonid
            } else {
                return Ok(())
            };

        let empty_string = flex_str!("");

        let header = "#gene_systematic_id\tgene_name\tmondo_id\treference\tdate\n";
        writer.write_all(header.as_bytes())?;

        for gene_details in self.genes.values() {
            if gene_details.taxonid != load_org_taxonid {
                continue;
            }

            if let Some(term_annotations) = gene_details.cv_annotations.get(&flex_str!("mondo")) {
                for term_annotation in term_annotations {

                    for annotation_id in &term_annotation.annotations {
                        let annotation_detail = self.get_annotation_detail(*annotation_id)
                            .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));

                        let gene_name = gene_details.name.as_ref().unwrap_or(&empty_string);
                        let reference =
                            annotation_detail.reference.as_ref().unwrap_or(&empty_string);
                        let date = annotation_detail.date.as_ref().unwrap_or(&empty_string);
                        let line = format!("{}\t{}\t{}\t{}\t{}\n",
                                           gene_details.uniquename,
                                           gene_name,
                                           term_annotation.term,
                                           reference,
                                           date);

                       writer.write_all(line.as_bytes())?;
                    }
                }
            }
        }

        Ok(())
    }

    fn write_modifications(&self, config: &Config, output_dir: &str)
        -> Result<(), io::Error>
    {
        let load_org_taxonid =
            if let Some(load_org_taxonid) = config.load_organism_taxonid {
                load_org_taxonid
            } else {
                return Ok(())
            };

        let file_name = format!("{}/modifications.tsv", output_dir);
        let file = File::create(file_name)?;
        let mut writer = BufWriter::new(&file);

        let header = "#gene_systematic_id\tgene_name\tmodification_term_id\tevidence\tmodification\textension\treference\ttaxon_id\tdate\tassigned_by\n";
        writer.write_all(header.as_bytes())?;

        let empty_string = flex_str!("");

        for gene_details in self.genes.values() {
            if gene_details.taxonid != load_org_taxonid {
                continue;
            }

            if let Some(term_annotations) = gene_details.cv_annotations.get(&flex_str!("PSI-MOD")) {
                for term_annotation in term_annotations {
                    if term_annotation.is_not {
                        continue;
                    }
                    for annotation_id in &term_annotation.annotations {
                        let annotation_detail = self.get_annotation_detail(*annotation_id)
                            .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));

                        let gene_name = gene_details.name.as_ref().unwrap_or(&empty_string);
                        let mut maybe_evidence = annotation_detail.evidence.clone();
                        if let Some(ref evidence) = maybe_evidence {
                            if let Some(ev_config) = config.evidence_types.get(evidence) {
                                maybe_evidence = Some(ev_config.long.to_shared_str());
                            }
                        }
                        let (modification, extension) =
                            process_modification_ext(config, self, &gene_details.uniquename,
                                                     &annotation_detail.extension);

                        let reference =
                            annotation_detail.reference.as_ref().unwrap_or(&empty_string);
                        let date = annotation_detail.date.as_ref().unwrap_or(&empty_string);
                        let assigned_by = annotation_detail.assigned_by.as_ref().unwrap_or(&empty_string);
                        let line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                           gene_details.uniquename,
                                           gene_name,
                                           term_annotation.term,
                                           maybe_evidence.unwrap_or_else(|| empty_string.clone()),
                                           modification,
                                           extension,
                                           reference,
                                           load_org_taxonid,
                                           date,
                                           assigned_by);

                       writer.write_all(line.as_bytes())?;
                    }
                }
            }
        }

        Ok(())
    }

    fn write_interactions(&self, config: &Config, output_dir: &str)
                          -> Result<(), io::Error>
    {
        let load_org_taxonid =
            if let Some(load_org_taxonid) = config.load_organism_taxonid {
                load_org_taxonid
            } else {
                return Ok(())
            };

        let file_name = format!("{}/interactions.tsv", output_dir);
        let file = File::create(file_name)?;
        let mut writer = BufWriter::new(&file);

        let header = "#InteractorA\tInteractorB\tInteractorA_taxon\tInteractorB_taxon\tEvidence_Code\tPubmedID\tScore\tModification\tPhenotypes\tComment\n";
        writer.write_all(header.as_bytes())?;

        let mut write_interaction = |int_type: &str, int: &InteractionAnnotation|
             -> Result<(), io::Error>
        {
            let mut comment = format!("interaction_type: {}", int_type);

            if let Some(ref source_database) = int.source_database {
                comment.push_str("|source: ");
                comment.push_str(source_database.as_str());
            }

            if let Some(ref annotation_date) = int.annotation_date {
                comment.push_str("|annotation_date: ");
                comment.push_str(annotation_date.as_str());
            }

            let line = format!("{}\t{}\t{}\t{}\t{}\t{}\t\t\t\t{}\n",
                               int.gene_uniquename,
                               int.interactor_uniquename,
                               load_org_taxonid, load_org_taxonid,
                               int.evidence,
                               int.reference_uniquename.as_deref().unwrap_or_default(),
                               comment);
            writer.write_all(line.as_bytes())?;
            Ok(())
        };

        for phys_interaction in &self.physical_interaction_annotations {
            write_interaction("physical", phys_interaction)?;
        }

        for genetic_interaction in &self.genetic_interaction_annotations {
            write_interaction("genetic", genetic_interaction)?;
        }

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

        let file = File::create(file_name)?;

        let mut writer = BufWriter::new(&file);

        let table = table_for_export(self, cv_config_map, subset_config);

        for row in table {
            let line = join(&row, "\t") + "\n";
            writer.write_all(line.as_bytes())?;
        }

        Ok(())
    }

    pub fn write_apicuron_files(&self, config: &Config,
                                references: &UniquenameReferenceMap,
                                output_dir: &str)
             -> Result<(), io::Error>
    {
        let mut curation_reports = vec![];

        for ref_details in references.values() {
            if ref_details.cv_annotations.is_empty() {
                continue;
            }

            let Some(ref canto_approved_date) = ref_details.canto_approved_date
            else {
                continue
            };

            let mut approved_parts = canto_approved_date.split(" ");
            let Some(approved_date) = approved_parts.next()
            else {
                continue;
            };

            let approved_time = approved_parts.next().unwrap_or_else(|| "00:00:00");

            let timestamp = flex_fmt!("{}T{}.000Z", approved_date, approved_time);

            let mut add_record = |orcid: &FlexStr, activity_term: FlexStr| {
                curation_reports.push(ApicuronReport {
                    activity_term,
                    curator_orcid: orcid.clone(),
                    timestamp: timestamp.clone(),
                    entity_uri: flex_fmt!("{}/reference/{}", config.base_url, ref_details.uniquename),
                });
            };

            let mut maybe_add_record =
                |annotation_curator: &AnnotationCurator, activity_term: FlexStr| {
                    if let Some(ref curator_orcid) = annotation_curator.orcid {
                        add_record(curator_orcid, activity_term);
                    }
                };

            for annotation_curator in &ref_details.annotation_curators {
                maybe_add_record(&annotation_curator, flex_str!("publication_curated"));
            }
            for annotation_curator in &ref_details.annotation_file_curators {
                maybe_add_record(&annotation_curator, flex_str!("provided_dataset"));
            }
            if let Some(ref canto_approver_orcid) = ref_details.canto_approver_orcid {
                add_record(canto_approver_orcid, flex_str!("approved_publication"));
            }
        }

        let apicuron_data = ApicuronData {
            resource_id: config.apicuron.resource_id.clone(),
            reports: curation_reports,
        };

        let s = serde_json::to_string(&apicuron_data).unwrap();
        let file_name = String::new() + output_dir + "/apicuron_data.json";
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes())?;

        Ok(())
    }

    pub fn write_stats(&self, output_dir: &str) -> Result<(), io::Error> {
        let s = serde_json::to_string(&self.stats).unwrap();
        let file_name = String::new() + output_dir + "/stats.json";
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes())?;

        Ok(())
    }

    pub fn write_detailed_stats(&self, output_dir: &str) -> Result<(), io::Error> {
        let s = serde_json::to_string(&self.detailed_stats).unwrap();
        let file_name = String::new() + output_dir + "/detailed_stats.json";
        let f = File::create(file_name)?;
        let mut writer = BufWriter::new(&f);
        writer.write_all(s.as_bytes())?;

        Ok(())
    }

    fn write_sqlite_db(&self, output_dir: &str) -> anyhow::Result<()> {
        let db_path = String::new() + output_dir + "/" + API_MAPS_SQLITE3_FILE_NAME;
        let mut conn = Connection::open(db_path)?;

        make_maps_database_tables(&mut conn)?;

        store_maps_into_database(&mut conn, &self.terms, &self.genes, &self.alleles,
                                 &self.references, &self.genotypes,
                                 &self.annotation_details,
                                 &self.termid_genotype_annotation)?;

        Ok(())
    }

    pub fn write(self, config: &Config, go_eco_mappping: &GoEcoMapping,
                 doc_config: &DocConfig, output_dir: &str)
                 -> Result<(), io::Error>
    {
        let web_json_path = self.create_dir(output_dir, "web-json")?;

        self.write_chromosome_json(config, &web_json_path)?;
        println!("wrote {} chromosomes", self.get_chromosomes().len());
        self.write_gene_summaries(&web_json_path)?;
        self.write_chromosome_summaries(&web_json_path)?;
        println!("wrote summaries");
        self.write_metadata(&web_json_path)?;
        println!("wrote metadata");
        self.write_recent_references(&web_json_path)?;
        self.write_all_community_curated(&web_json_path)?;
        self.write_all_admin_curated(&web_json_path)?;
        println!("wrote references");

        self.write_sqlite_db(&output_dir).unwrap();
        println!("wrote SQLite DB");

        self.write_api_maps(&web_json_path)?;
        println!("wrote API maps");
        self.write_solr_data(&web_json_path)?;
        println!("wrote search data");

        let intermine_data_path = self.create_dir(output_dir, "intermine_data")?;
        self.write_intermine_data(config, &intermine_data_path)?;
        println!("wrote intermine data");

        self.write_subsets(&web_json_path)?;
        println!("wrote subsets");

        let fasta_path = self.create_dir(output_dir, "fasta")?;
        let feature_sequences_path = self.create_dir(&fasta_path, "feature_sequences")?;
        self.write_feature_sequences(&feature_sequences_path)?;
        let chromosomes_path = self.create_dir(&fasta_path, "chromosomes")?;
        self.write_chromosome_sequences(config, &chromosomes_path)?;
        println!("wrote fasta");

        let misc_path = self.create_dir(output_dir, "misc")?;

        write_go_annotation_files(&self.api_maps, config, &self,
                                  &self.metadata.db_creation_datetime,
                                  go_eco_mappping, &self.genes,
                                  &self.api_maps.transcripts,
                                  &self.api_maps.protein_complexes,
                                  &misc_path)?;

        write_phenotype_annotation_files(&self, &self.genotypes, config, false, &misc_path)?;
        write_phenotype_annotation_files(&self, &self.genotypes, config, true, &misc_path)?;
        println!("wrote GAF and PHAF files");

        self.write_alleles_json(&misc_path)?;
        println!("wrote allele data to JSON");

        self.write_gene_id_table(config, &misc_path)?;
        self.write_protein_features(config, &misc_path)?;
        self.write_feature_coords(config, &misc_path)?;
        write_macromolecular_complexes(&self.api_maps.protein_complex_data, &misc_path)?;
        self.write_rnacentral(config, &misc_path)?;
        self.write_deletion_viability(config, &misc_path)?;
        self.write_slim_ids_and_names(config, &misc_path)?;
        self.write_transmembrane_domains(config, &misc_path)?;
        self.write_htp_gene_expression_table(&misc_path)?;
        self.write_full_gene_expression_table(&misc_path)?;
        self.write_site_map_txt(config, doc_config, &self.references, &misc_path)?;
        self.write_allele_tsv(&misc_path)?;
        self.write_disease_association(config, &misc_path)?;
        self.write_modifications(config, &misc_path)?;
        self.write_interactions(config, &misc_path)?;

        self.write_annotation_subsets(config, &misc_path)?;

        self.write_apicuron_files(config, &self.references, &misc_path)?;

        println!("wrote miscellaneous files");

        self.write_stats(&web_json_path)?;
        self.write_detailed_stats(&web_json_path)?;

        println!("wrote stats files");

        let gff_path = self.create_dir(output_dir, "gff")?;
        self.write_gff(config, &gff_path)?;

        println!("wrote GFF files");

        Ok(())
    }
}


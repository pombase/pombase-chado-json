use std::io::{self, BufWriter};



use std::io::Write;
use std::fs::File;

use flexstr::{shared_fmt as flex_fmt, ToSharedStr};

use crate::web::config::*;
use crate::data_types::*;


use super::go_format_writer::GpadGafWriteMode;
use super::util::make_extension_string;


pub fn write_phenotype_annotation_files(api_maps: &APIMaps,
                                        config: &Config,
                                        use_eco_evidence: bool,
                                        output_dir: &str)
  -> Result<(), io::Error>
{
    let load_org_taxonid =
        if let Some(load_org_taxonid) = config.load_organism_taxonid {
            load_org_taxonid
        } else {
            return Ok(())
        };

    let phaf_cv_name = config.file_exports.phaf_cv_name.to_shared_str();

    let phaf_parental_strain = config.file_exports.phaf_parental_strain.get(&load_org_taxonid)
        .unwrap_or_else(|| panic!("no phaf_parental_strain configured for {}", load_org_taxonid));

    let database_name = &config.database_name;

    let eco_ev_bit =
        if use_eco_evidence {
            "_eco_evidence"
        } else {
            ""
        };
    let phaf_file_name =
        format!("{}/single_locus_phenotype_annotations_taxon_{}{}.phaf", output_dir,
                load_org_taxonid, eco_ev_bit);
    let phaf_file = File::create(phaf_file_name).expect("Unable to open file");
    let mut phaf_writer = BufWriter::new(&phaf_file);

    let header = "#Database name\tGene systematic ID\tFYPO ID\tAllele description\tExpression\tParental strain\tStrain name (background)\tGenotype description\tGene name\tAllele name\tAllele synonym\tAllele type\tEvidence\tCondition\tPenetrance\tSeverity\tExtension\tReference\tTaxon\tDate\tPloidy\n";

    phaf_writer.write_all(header.as_bytes())?;

    'GENOTYPES:
    for genotype_details in api_maps.genotypes.values() {
        if genotype_details.taxonid != load_org_taxonid {
            continue 'GENOTYPES;
        }

        // only export single locus genotypes that are haploid or homozygous diploid
        let expressed_allele =
            if genotype_details.loci.len() == 1 {
                let locus = &genotype_details.loci[0];

                let expressed_alleles = &locus.expressed_alleles;

                for test_expressed_allele in &expressed_alleles[1..] {
                    if expressed_alleles[0] != *test_expressed_allele {
                        continue 'GENOTYPES;
                    }
                }

                &expressed_alleles[0]
            } else {
                continue;
            };

        let is_homozygous_diploid = genotype_details.loci[0].expressed_alleles.len() > 1;

        let locus_allele =
            api_maps.alleles.get(&expressed_allele.allele_uniquename)
            .unwrap_or_else(|| panic!("no allele found for {}", expressed_allele.allele_uniquename));

        let locus_allele_synonyms =
            locus_allele.synonyms.iter().map(|s| s.name.clone()).collect::<Vec<_>>().join("|");

        let locus_gene = &locus_allele.gene;

        let locus_gene_name_or_uniquename =
            if let Some(ref name) = locus_gene.name {
                name
            } else {
                &locus_gene.uniquename
            };

        let expression =
            if let Some(ref expression) = expressed_allele.expression {
                if expression == "Null" && locus_allele.allele_type == "deletion" {
                    flex_fmt!("")
                } else {
                    expression.clone()
                }
            } else {
                flex_fmt!("")
            };

        if let Some(term_annotations) = genotype_details.cv_annotations.get(&phaf_cv_name) {
            for term_annotation in term_annotations {
                let term =
                    api_maps.terms.get(&term_annotation.term)
                    .unwrap_or_else(|| panic!("failed to find term summary for {}",
                                              term_annotation.term));

                for annotation_id in &term_annotation.annotations {
                    let annotation_detail = api_maps.annotation_details
                        .get(annotation_id)
                        .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));

                    let evidence =
                        if use_eco_evidence {
                            annotation_detail.eco_evidence.clone().unwrap_or_else(|| flex_fmt!(""))
                        } else {
                            annotation_detail.evidence.clone().unwrap_or_else(|| flex_fmt!(""))
                        };

                    let conditions =
                        annotation_detail.conditions.iter()
                        .map(|fs| fs.as_str())
                        .collect::<Vec<_>>()
                        .join(",");

                    let penetrance =
                        annotation_detail.extension.iter()
                        .filter(|bit| {
                            bit.rel_type_name == "has_penetrance"
                        })
                        .map(|bit| if let ExtRange::Term(ref termid) = bit.ext_range {
                            let bit_term = api_maps.terms.get(termid)
                                .unwrap_or_else(|| panic!("can't find term for {}", termid));
                            bit_term.name.to_string()
                        } else {
                            bit.ext_range.to_string()
                        })
                        .collect::<Vec<_>>()
                        .join(",");

                    let severity =
                        annotation_detail.extension.iter()
                        .filter(|bit| {
                            bit.rel_type_name == "has_severity"
                        })
                        .map(|bit| if let ExtRange::Term(ref termid) = bit.ext_range {
                            let bit_term = api_maps.terms.get(termid)
                                .unwrap_or_else(|| panic!("can't find term for {}", termid));
                            bit_term.name.to_string()
                        } else {
                            bit.ext_range.to_string()
                        })
                        .collect::<Vec<_>>()
                        .join(",");

                    let extension_bits =
                        annotation_detail.extension.iter()
                        .filter(|bit| {
                            bit.rel_type_name != "has_penetrance" && bit.rel_type_name != "has_severity"
                        })
                        .map(|bit| bit.clone())
                        .collect::<Vec<_>>();

                    let extension =
                        make_extension_string(config, &api_maps.terms,
                                              &GpadGafWriteMode::PomBaseGaf,
                                              &extension_bits);

                    let date =
                        annotation_detail.date.clone().unwrap_or_else(|| flex_fmt!("NO_DATE"));

                    let reference_uniquename =
                        if let Some(ref reference) = annotation_detail.reference {
                            reference.clone()
                        } else {
                            flex_fmt!("")
                        };

                    let line =
                        format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                database_name,
                                locus_gene.uniquename,
                                term.termid,
                                locus_allele.description.clone().unwrap_or_else(|| flex_fmt!("")),
                                expression,
                                phaf_parental_strain,
                                "", //annotation_detail.genotype_background.clone().unwrap_or_else(|| flex_fmt!("")),
                                locus_gene_name_or_uniquename,
                                locus_allele.name.clone().unwrap_or_else(|| flex_fmt!("")),
                                locus_allele_synonyms,
                                locus_allele.allele_type,
                                evidence,
                                conditions,
                                penetrance,
                                severity,
                                extension,
                                reference_uniquename,
                                load_org_taxonid,
                                date,
                                if is_homozygous_diploid { "homozygous diploid" } else { "haploid" },
                        );
                    phaf_writer.write_all(line.as_bytes())?;
                }
            }

        }
    }

  Ok(())
}

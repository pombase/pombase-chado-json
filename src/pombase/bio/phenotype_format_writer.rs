use std::io::{self, BufWriter};

use std::io::Write;
use std::fs::File;

use flexstr::{shared_fmt as flex_fmt, FlexStr, ToSharedStr};

use crate::web::config::*;
use crate::data_types::*;

use itertools::Itertools;

use super::go_format_writer::GpadGafWriteMode;
use crate::bio::util::{make_extension_string, ANNOTATION_COMMENT_NESTED_BRACKETS_RE, ANNOTATION_COMMENT_RE};

#[derive(Debug, PartialEq, Eq)]
pub enum FypoEvidenceType {
    PomBase,
    Eco,
}

#[derive(Debug, PartialEq, Eq)]
pub enum FypoExportComments {
    Export,
    NoExport,
}

pub fn write_phenotype_annotation_files(data_lookup: &dyn DataLookup,
                                        genotypes_map: &IdGenotypeMap,
                                        config: &Config,
                                        evidence_type: FypoEvidenceType,
                                        export_comments: FypoExportComments,
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


    let phaf_file_name =
        if export_comments == FypoExportComments::Export {
            format!("{}/canto_fypo_annotations_with_comments.tsv", output_dir)
        } else {
            let eco_ev_bit =
                if evidence_type == FypoEvidenceType::Eco {
                    "_eco_evidence"
                } else {
                    ""
                };
            format!("{}/single_locus_phenotype_annotations_taxon_{}{}.phaf", output_dir,
                    load_org_taxonid, eco_ev_bit)
        };

    let phaf_file = File::create(phaf_file_name).expect("Unable to open file");
    let mut phaf_writer = BufWriter::new(&phaf_file);

    let comment_header =
        if export_comments == FypoExportComments::Export {
            "\tAnnotation comment"
        } else {
            ""
        };

    let header = format!("#Database name\tGene systematic ID\tFYPO ID\tAllele description\tExpression\tParental strain\tStrain name (background)\tGenotype description\tGene symbol\tAllele name\tAllele synonym\tAllele type\tEvidence\tCondition\tPenetrance\tSeverity\tExtension\tReference\tTaxon\tDate\tPloidy{}\n",
                         comment_header);

    phaf_writer.write_all(header.as_bytes())?;

  'GENOTYPES:
    for genotype_details in genotypes_map.values() {
        if genotype_details.taxonid != load_org_taxonid {
            continue 'GENOTYPES;
        }

        // only export single locus genotypes that are haploid or homozygous diploid
        let Some(expressed_allele) = get_expressed_allele(genotype_details)
        else {
            continue 'GENOTYPES;
        };

        let is_homozygous_diploid = genotype_details.loci[0].expressed_alleles.len() > 1;

        let locus_allele =
            data_lookup.get_allele(&expressed_allele.allele_uniquename)
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
                    data_lookup.get_term(&term_annotation.term)
                    .unwrap_or_else(|| panic!("failed to find term summary for {}",
                                              term_annotation.term));

                for annotation_id in &term_annotation.annotations {
                    let annotation_detail = data_lookup.get_annotation_detail(*annotation_id)
                        .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));

                    let reference_uniquename =
                        if let Some(ref reference) = annotation_detail.reference {
                            reference.clone()
                        } else {
                            flex_fmt!("")
                        };

                    if export_comments == FypoExportComments::Export {
                        if let Some(reference_details) = data_lookup.get_reference(&reference_uniquename) {
                            if !reference_details.is_canto_curated() {
                                continue;
                            }
                        } else {
                            continue;
                        }
                    }

                    let comment_field =
                        if export_comments == FypoExportComments::Export {
                            if let Some(submitter_comment) = get_submitter_comment(annotation_detail.as_ref()) {
                                submitter_comment
                            } else {
                                format!("\t")
                            }
                        } else {
                            String::default()
                        };

                    let evidence =
                        if evidence_type == FypoEvidenceType::Eco {
                            annotation_detail.eco_evidence.clone().unwrap_or_else(|| flex_fmt!(""))
                        } else {
                            annotation_detail.evidence.clone().unwrap_or_else(|| flex_fmt!(""))
                        };

                    let conditions =
                        annotation_detail.conditions.iter()
                        .map(|fs| fs.as_str())
                        .sorted()
                        .join(",");

                    let penetrance = annotation_penetrance(data_lookup, &annotation_detail);

                    let severity = annotation_severity(data_lookup, &annotation_detail);

                    let extension_bits =
                        annotation_detail.extension.iter()
                        .filter(|bit| {
                            bit.rel_type_name != "has_penetrance" && bit.rel_type_name != "has_severity"
                        }).cloned()
                        .collect::<Vec<_>>();

                    let extension =
                        make_extension_string(config, data_lookup,
                                              &GpadGafWriteMode::PomBaseGaf,
                                              &extension_bits);

                    let date =
                        annotation_detail.date.clone().unwrap_or_else(|| flex_fmt!("NO_DATE"));

                    let line =
                        format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}{}\n",
                                database_name,
                                locus_gene.uniquename,
                                term.termid,
                                locus_allele.description.clone().unwrap_or_else(|| flex_fmt!("")),
                                expression,
                                phaf_parental_strain,
                                "",
                                annotation_detail.genotype_background.clone().unwrap_or_default(),
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
                                comment_field,
                        );
                    phaf_writer.write_all(line.as_bytes())?;
                }
            }

        }
    }

  Ok(())
}

pub fn write_heterozygous_diploid_annotations(data_lookup: &dyn DataLookup,
                                              genotypes_map: &IdGenotypeMap,
                                              config: &Config,
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

    let database_name = &config.database_name;

    let dominant_file_name = format!("{}/dominant_diploid_annotations.tsv", output_dir);
    let dominant_file = File::create(dominant_file_name).expect("Unable to open file");
    let mut dominant_writer = BufWriter::new(&dominant_file);

    let all_heterozygous_file_name = format!("{}/single_locus_heterozygous_diploid_annotations.tsv", output_dir);
    let all_heterozygous_file = File::create(all_heterozygous_file_name).expect("Unable to open file");
    let mut all_heterozygous_writer = BufWriter::new(&all_heterozygous_file);

    let header = "#database_name\tgene_systematic_id\tgene_name\tfypo_term_id\tfypo_term_name\tallele_1_name\tallele_1_description\tallele_1_type\tallele_1_expression\tallele_2_name\tallele_2_description\tallele_2_type\tallele_2_expression\tevidence\tconditions\tpenetrance\tseverity\textension\treference\ttaxon_id\tdate";
    writeln!(dominant_writer, "{}", header)?;
    writeln!(all_heterozygous_writer, "{}", header)?;

  'GENOTYPES:
    for genotype_details in genotypes_map.values() {
        if genotype_details.taxonid != load_org_taxonid {
            continue 'GENOTYPES;
        }

        // only export single locus genotypes that are heterozygous diploids
        let Some((expressed_allele_1, expressed_allele_2)) = get_heterozygous_alleles(genotype_details)
        else {
            continue 'GENOTYPES;
        };

        let allele_1 =
            data_lookup.get_allele(&expressed_allele_1.allele_uniquename)
            .unwrap_or_else(|| panic!("no allele found for {}", expressed_allele_1.allele_uniquename));

        let allele_2 =
            data_lookup.get_allele(&expressed_allele_2.allele_uniquename)
            .unwrap_or_else(|| panic!("no allele found for {}", expressed_allele_2.allele_uniquename));

        let mut is_dominant_diploid = true;

        // it's not a dominant allele if the other allele is a deletion
        if allele_1.allele_type == "deletion" || allele_2.allele_type == "deletion" {
            is_dominant_diploid = false;
        }

        // it's not a dominant allele unless the other allele is WT
        if allele_1.allele_type != "wild_type" && allele_2.allele_type != "wild_type" {
            is_dominant_diploid = false;
        }

        if allele_1.gene != allele_2.gene {
            panic!("differing genes for {} {}", expressed_allele_1.allele_uniquename,
                   expressed_allele_2.allele_uniquename);
        }

        let gene = &allele_1.gene;

        let gene_name_or_uniquename =
            if let Some(ref name) = gene.name {
                name
            } else {
                &gene.uniquename
            };

        let expression_1 =
            if let Some(ref expression) = expressed_allele_1.expression {
                if expression == "Null"  {
                    is_dominant_diploid = false;
                }
                expression.clone()
            } else {
                flex_fmt!("")
            };

        let expression_2 =
            if let Some(ref expression) = expressed_allele_2.expression {
                if expression == "Null"  {
                    is_dominant_diploid = false;
                }
                expression.clone()
            } else {
                flex_fmt!("")
            };

        if let Some(term_annotations) = genotype_details.cv_annotations.get(&phaf_cv_name) {
            for term_annotation in term_annotations {
                let term =
                    data_lookup.get_term(&term_annotation.term)
                    .unwrap_or_else(|| panic!("failed to find term summary for {}",
                                              term_annotation.term));

                for annotation_id in &term_annotation.annotations {
                    let annotation_detail = data_lookup.get_annotation_detail(*annotation_id)
                        .unwrap_or_else(|| panic!("can't find annotation {}", annotation_id));

                    let evidence =
                        annotation_detail.evidence.clone().unwrap_or_else(FlexStr::default);

                    let conditions = annotation_detail.conditions.iter()
                        .map(|fs| fs.as_str())
                        .sorted()
                        .join(",");

                    let penetrance = annotation_penetrance(data_lookup, &annotation_detail);

                    let severity = annotation_severity(data_lookup, &annotation_detail);

                    let extension_bits =
                        annotation_detail.extension.iter()
                        .filter(|bit| {
                            bit.rel_type_name != "has_penetrance" && bit.rel_type_name != "has_severity"
                        }).cloned()
                        .collect::<Vec<_>>();

                    let extension =
                        make_extension_string(config, data_lookup,
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
                        format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                                database_name,
                                gene.uniquename,
                                gene_name_or_uniquename,
                                term.termid,
                                term.name,
                                allele_1.name.clone().unwrap_or_else(FlexStr::default),
                                allele_1.description.clone().unwrap_or_else(FlexStr::default),
                                allele_1.allele_type,
                                expression_1,
                                allele_2.name.clone().unwrap_or_else(FlexStr::default),
                                allele_2.description.clone().unwrap_or_else(FlexStr::default),
                                allele_2.allele_type,
                                expression_2,
                                evidence,
                                conditions,
                                penetrance,
                                severity,
                                extension,
                                reference_uniquename,
                                load_org_taxonid,
                                date,
                        );

                    if is_dominant_diploid {
                        dominant_writer.write_all(line.as_bytes())?;
                    }
                    all_heterozygous_writer.write_all(line.as_bytes())?;
                }
            }

        }
    }

  Ok(())
}

fn annotation_penetrance(data_lookup: &dyn DataLookup, annotation_detail: &OntAnnotationDetail)
    -> String
{
    annotation_detail.extension.iter()
        .filter(|bit| {
            bit.rel_type_name == "has_penetrance"
        })
        .map(|bit| if let ExtRange::Term(ref termid) = bit.ext_range {
            let bit_term = data_lookup.get_term(termid)
                .unwrap_or_else(|| panic!("can't find term for {}", termid));
            bit_term.name.to_string()
        } else {
            bit.ext_range.to_string()
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn annotation_severity(data_lookup: &dyn DataLookup, annotation_detail: &OntAnnotationDetail)
    -> String
{
    annotation_detail.extension.iter()
        .filter(|bit| {
            bit.rel_type_name == "has_severity"
        })
        .map(|bit| if let ExtRange::Term(ref termid) = bit.ext_range {
            let bit_term = data_lookup.get_term(termid)
                .unwrap_or_else(|| panic!("can't find term for {}", termid));
            bit_term.name.to_string()
        } else {
            bit.ext_range.to_string()
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn get_submitter_comment(annotation_detail: &OntAnnotationDetail)
   -> Option<String>
{
    let Some(ref submitter_comment) = annotation_detail.submitter_comment
    else {
        return None;
    };

    // remove comments while handling nested brackets (in a hacky way)
    let submitter_comment =
        ANNOTATION_COMMENT_NESTED_BRACKETS_RE.replace_all(submitter_comment, "");
    let submitter_comment =
        ANNOTATION_COMMENT_RE.replace_all(&submitter_comment, "");
    let submitter_comment = submitter_comment.trim();

    if submitter_comment.is_empty() {
        None
    } else {
        Some(format!("\t{}", submitter_comment))
    }
}

fn get_expressed_allele(genotype_details: &GenotypeDetails) -> Option<&ExpressedAllele> {
    if genotype_details.loci.len() == 1 {
        let locus = &genotype_details.loci[0];

        let expressed_alleles = &locus.expressed_alleles;

        for test_expressed_allele in &expressed_alleles[1..] {
            if expressed_alleles[0] != *test_expressed_allele {
                return None;
            }
        }

        Some(&expressed_alleles[0])
    } else {
        None
    }
}

fn get_heterozygous_alleles(genotype_details: &GenotypeDetails)
    -> Option<(&ExpressedAllele, &ExpressedAllele)>
{
    if genotype_details.loci.len() == 1 {
        let locus = &genotype_details.loci[0];

        let expressed_alleles = &locus.expressed_alleles;

        if expressed_alleles.len() != 2 {
            return None;
        }

        let allele_1 = &expressed_alleles[0];
        let allele_2 = &expressed_alleles[1];

        if allele_1.allele_uniquename == allele_2.allele_uniquename {
            return None;
        }

        Some((allele_1, allele_2))
    } else {
        None
    }
}

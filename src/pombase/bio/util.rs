
use crate::web::config::RelationOrder;
use crate::types::GeneUniquename;

use flexstr::{SharedStr as FlexStr, shared_fmt as flex_fmt, ToSharedStr};

use crate::web::config::Config;
use crate::data_types::*;

use super::go_format_writer::GpadGafWriteMode;


fn complement_char(base: char) -> char {
    match base {
        'a' => 't',
        'A' => 'T',
        't' => 'a',
        'T' => 'A',
        'g' => 'c',
        'G' => 'C',
        'c' => 'g',
        'C' => 'G',
        _ => 'n',
    }
}

pub fn rev_comp(residues: &str) -> FlexStr {
    let residues: String = residues.chars()
        .rev().map(complement_char)
        .collect();
    residues.to_shared_str()
}


pub fn format_fasta(id: &str, maybe_desc: Option<String>,
                    seq: &str, width: usize) -> String {
    let mut ret = ">".to_owned() + id;

    if let Some(desc) = maybe_desc {
        ret.push(' ');
        ret.push_str(&desc);
    }

    ret.push('\n');

    if seq.is_empty() {
        ret.push('\n');
    } else {
        let mut count = 0;
        for c in seq.chars() {
            ret.push(c);
            count += 1;
            if count % width == 0 {
                ret.push('\n');
            }
        }

        if count % width != 0 {
            ret.push('\n');
        }
    }

    ret
}

fn to_gff(chromosome_export_id: &str,
          source: &str, feat_id: &str, maybe_name: Option<&str>, feat_type: &str,
          location: &ChromosomeLocation, maybe_parent: Option<&str>) -> String {
    let phase_char =
        if let Some(ref phase) = location.phase {
            phase.to_gff_str()
        } else if feat_type == "CDS" {
            // phase is required for CDS features
            "0"
        } else {
            "."
        };
    let start_pos_str = format!("{}", location.start_pos);
    let end_pos_str = format!("{}", location.end_pos);
    let mut ret_val = String::new();
    let bits: Vec<&str> =
        vec![chromosome_export_id, source, feat_type,
             &start_pos_str, &end_pos_str,
             ".", location.strand.to_gff_str(), phase_char];
    for s in bits {
        ret_val.push_str(s);
        ret_val.push('\t');
    }

    ret_val.push_str("ID=");
    ret_val.push_str(feat_id);

    if let Some(parent_id) = maybe_parent {
        ret_val.push_str(";Parent=");
        ret_val.push_str(parent_id);
    };

    if let Some(name) = maybe_name {
        ret_val.push_str(";Name=");
        ret_val.push_str(name);
    };

    ret_val
}

fn get_gff_gene_type(gene: &GeneDetails) -> &'static str {
    match gene.feature_type.as_ref() {
        "gene" => "gene",
        "lncRNA gene" => "lncRNA_gene",
        "mRNA gene" => "protein_coding_gene",
        "ncRNA gene" => "ncRNA_gene",
        "pseudogene" => "pseudogene",
        "rRNA gene" => "rRNA_gene",
        "snRNA gene" => "snRNA_gene",
        "sncRNA gene" => "sncRNA_gene",
        "snoRNA gene" => "snoRNA_gene",
        "tRNA gene" => "tRNA_gene",
        _ => panic!("unknown gene feature_type: {}", gene.feature_type),
    }
}

pub fn format_gene_gff(chromosome_export_id: &str,
                       source: &str, transcripts: &UniquenameTranscriptMap,
                       gene: &GeneDetails) -> Vec<String> {
    let mut ret_val = vec![];
    if let Some (ref gene_loc) = gene.location {
        let maybe_gene_name =
            gene.name.as_ref().map(|gene_name| gene_name as &str);

        let gene_type = get_gff_gene_type(gene);

        let gene_gff_line =
            to_gff(chromosome_export_id, source, &gene.uniquename, maybe_gene_name,
                   "gene", gene_loc, None);

        ret_val.push(format!("{};so_term_name={}", gene_gff_line, gene_type));

        for transcript_uniquename in &gene.transcripts {
            let transcript_details = transcripts
                .get(transcript_uniquename)
                .unwrap_or_else(|| panic!("internal error, failed to find transcript: {}",
                                 transcript_uniquename));

            ret_val.push(to_gff(chromosome_export_id,
                                source, transcript_uniquename, None,
                                &transcript_details.transcript_type,
                                &transcript_details.location,
                                Some(&gene.uniquename)));
            for part in &transcript_details.parts {
                let gff_feat_type =
                    match part.feature_type {
                        FeatureType::Exon => "CDS".to_owned(),
                        FeatureType::CdsIntron | FeatureType::FivePrimeUtrIntron |
                        FeatureType::ThreePrimeUtrIntron => "intron".to_owned(),
                        FeatureType::FivePrimeUtr => "five_prime_UTR".to_owned(),
                        FeatureType::ThreePrimeUtr => "three_prime_UTR".to_owned(),
                        _ => format!("{}", part.feature_type),
                    };
                ret_val.push(to_gff(chromosome_export_id,
                                    source, &part.uniquename, None, &gff_feat_type,
                                    &part.location,
                                    Some(transcript_uniquename)));
            }
        }
    }
    ret_val
}

pub fn format_misc_feature_gff(chromosome_export_id: &str,
                               source: &str, feature_short: &FeatureShort) -> Vec<String> {
    let mut ret_val = vec![];
    let feature_type_name = format!("{}", feature_short.feature_type);
    ret_val.push(to_gff(chromosome_export_id,
                        source, &feature_short.uniquename, None,
                        &feature_type_name, &feature_short.location, None));
    ret_val
}

pub fn process_modification_ext(config: &Config, data_lookup: &dyn DataLookup,
                                gene_uniquename: &GeneUniquename, extension: &[ExtPart])
    -> (FlexStr, FlexStr)
{
    let mut remaining = vec![];
    let mut mod_res: Vec<&str> = vec![];

    for ext_part in extension {
        if ext_part.rel_type_name == "residue" ||
           ext_part.rel_type_name == "modified residue" {
            let ExtRange::ModifiedResidues(ref modified_residues) = ext_part.ext_range
            else {
                panic!("unknown ext_range for: {:?}", ext_part);
            };

            let Some(ref cv_conf) = config.cv_config.get("PSI-MOD")
            else {
                panic!("no modification CV configured - can't find \"PSI-MOD\"");
            };

            for res in modified_residues {
                let Some(mod_abbrev_config) = cv_conf.modification_abbreviations.get(gene_uniquename)
                else {
                    mod_res.push(res);
                    continue;
                };

                let Some(replacement_residues) = mod_abbrev_config.get(res)
                else {
                    mod_res.push(res);
                    continue;
                };

                for replacement_res in replacement_residues.split(",") {
                    mod_res.push(replacement_res);
                }
            }
        } else {
            remaining.push(ext_part.clone());
        }
    }

    let compare_ext_part_func =
      |e1: &ExtPart, e2: &ExtPart| {
        e1.rel_type_name.cmp(&e2.rel_type_name)
      };

    remaining.sort_by(compare_ext_part_func);

    let remaining_str = make_extension_string(config, data_lookup, &GpadGafWriteMode::StandardGaf,
                                              &remaining);

    (mod_res.join(",").to_shared_str(), remaining_str)
}

pub fn make_extension_string(config: &Config, data_lookup: &dyn DataLookup,
                         write_mode: &GpadGafWriteMode, extension: &[ExtPart])
                         -> FlexStr
{
    let rel_mapping = &config.file_exports.gpad_gpi.extension_relation_mappings;
    let get_rel_term = |ext_part: &ExtPart| {
        if *write_mode == GpadGafWriteMode::PomBaseGaf ||
           *write_mode == GpadGafWriteMode::ExtendedPomBaseGaf {
            ext_part.rel_type_name.clone()
        } else {
            let rel_term_id =
                if let Some(map_termid) = rel_mapping.get(&ext_part.rel_type_name) {
                    map_termid.clone().unwrap_or_else(|| panic!("internal error, no mapping for {}",
                                                                &ext_part.rel_type_name))
                } else {
                    ext_part.rel_type_id.clone().unwrap()
                };

            if *write_mode == GpadGafWriteMode::Gpad {
                rel_term_id
            } else {
                data_lookup.get_term(&rel_term_id)
                    .unwrap_or_else(|| panic!("internal error, can't find term {}",
                                              rel_term_id))
                    .name.clone()
            }
        }
    };

    let get_range = |ext_part: &ExtPart| {
        match ext_part.ext_range {
            ExtRange::Gene(ref gene_uniquename) |
            ExtRange::Promoter(ref gene_uniquename) => {
                if !gene_uniquename.contains(':') {
                    let new_uniquename =
                        flex_fmt!("{}:{}", config.database_name, gene_uniquename);
                    ExtRange::Gene(new_uniquename)
                } else {
                    ext_part.ext_range.clone()
                }
            },
            ExtRange::GeneAndGeneProduct(GeneAndGeneProduct { gene_uniquename: _, ref product }) => {
                // avoid "SPCC622.08c (PR:000027564)" in extensions
                ExtRange::GeneProduct(product.clone())
            },
            _ => ext_part.ext_range.clone(),
        }
    };

    let mut parts = extension.iter()
        .filter(|ext_part| {
            if ext_part.rel_type_name == "binding site" {
                return false;
            }
            if *write_mode == GpadGafWriteMode::PomBaseGaf ||
               *write_mode == GpadGafWriteMode::ExtendedPomBaseGaf {
                true
            } else {
                if ext_part.rel_type_name == "residue" ||
                   ext_part.rel_type_name == "modified_residue" {
                    return false;
                }
                if let Some(map_termid) = rel_mapping.get(&ext_part.rel_type_name) {
                    map_termid.is_some()
                } else {
                    ext_part.rel_type_id.is_some()
                }
            }
        })
        .map(|ext_part| format!("{}({})", get_rel_term(ext_part),
                                get_range(ext_part)))
        .collect::<Vec<_>>();

    parts.sort();

    parts.join(",").to_shared_str()
}

pub fn compare_ext_part_with_config(extension_relation_order: &RelationOrder,
                                    ep1: &ExtPart, ep2: &ExtPart) -> Ordering {
    let rel_order_conf = extension_relation_order;
    let order_conf = &rel_order_conf.relation_order;
    let always_last_conf = &rel_order_conf.always_last;

    let maybe_ep1_index = order_conf.iter().position(|r| *r == ep1.rel_type_name);
    let maybe_ep2_index = order_conf.iter().position(|r| *r == ep2.rel_type_name);

    if let Some(ep1_index) = maybe_ep1_index {
        if let Some(ep2_index) = maybe_ep2_index {
            ep1_index.cmp(&ep2_index)
        } else {
            Ordering::Less
        }
    } else {
        if maybe_ep2_index.is_some() {
            Ordering::Greater
        } else {
            let maybe_ep1_last_index = always_last_conf.iter().position(|r| *r == ep1.rel_type_name);
            let maybe_ep2_last_index = always_last_conf.iter().position(|r| *r == ep2.rel_type_name);

            if let Some(ep1_last_index) = maybe_ep1_last_index {
                if let Some(ep2_last_index) = maybe_ep2_last_index {
                    ep1_last_index.cmp(&ep2_last_index)
                } else {
                    Ordering::Greater
                }
            } else {
                if maybe_ep2_last_index.is_some() {
                    Ordering::Less
                } else {
                    let name_cmp = ep1.rel_type_name.cmp(&ep2.rel_type_name);

                    if name_cmp == Ordering::Equal {
                        if ep1.ext_range.is_gene() && !ep2.ext_range.is_gene() {
                            Ordering::Less
                        } else {
                            if !ep1.ext_range.is_gene() && ep2.ext_range.is_gene() {
                                Ordering::Greater
                            } else {
                                Ordering::Equal
                            }
                        }
                    } else {
                        name_cmp
                    }
                }
            }
        }
    }
}

#[test]
fn test_format_fasta() {
    assert!(format_fasta("id1", None, "", 8) == ">id1\n\n");
    assert!(format_fasta("id2", Some("desc1".to_owned()), "atgc", 4) ==
            ">id2 desc1\natgc\n");
    assert!(format_fasta("id3", None, "atgc", 4) == ">id3\natgc\n");
    let input1 = "acgtacgattattaccggttacgcatccgtgtaaca";
    let expected1 = ">id4 desc4\nacgtacga\nttattacc\nggttacgc\natccgtgt\naaca\n";
    assert!(format_fasta("id4", Some("desc4".to_owned()), input1, 8) == expected1);
    let input2 = "acgtacgattattaccggttacgcatccgtgt";
    let expected2 = ">id5\nacgtacga\nttattacc\nggttacgc\natccgtgt\n";
    assert!(format_fasta("id5", None, input2, 8) == expected2);
}

#[test]
fn test_format_gff() {
    let gene = make_test_gene();
    let mut transcripts = std::collections::HashMap::new();
    transcripts.insert( flex_str!("SPCC18B5.06.1"),
                       gene.transcripts_by_uniquename[&flex_str!("SPCC18B5.06.1")].clone().unwrap());
    let gene_gff_lines = format_gene_gff("chromosome_3", "PomBase",
                         &transcripts, &gene);

    assert_eq!(gene_gff_lines.len(), 13);
    assert_eq!(gene_gff_lines[0],
               "chromosome_3\tPomBase\tgene\t729054\t730829\t.\t+\t.\tID=SPCC18B5.06;Name=dom34;so_term_name=protein_coding_gene");
    assert_eq!(gene_gff_lines[1],
               "chromosome_3\tPomBase\tmRNA\t729054\t730829\t.\t+\t.\tID=SPCC18B5.06.1;Parent=SPCC18B5.06");
    assert_eq!(gene_gff_lines[2],
               "chromosome_3\tPomBase\tfive_prime_UTR\t729054\t729132\t.\t+\t.\tID=SPCC18B5.06.1:five_prime_UTR:1;Parent=SPCC18B5.06.1");
    assert_eq!(gene_gff_lines[3],
               "chromosome_3\tPomBase\tCDS\t729133\t729212\t.\t+\t0\tID=SPCC18B5.06.1:exon:1;Parent=SPCC18B5.06.1");
    assert_eq!(gene_gff_lines[4],
               "chromosome_3\tPomBase\tintron\t729213\t729265\t.\t+\t.\tID=SPCC18B5.06.1:intron:1;Parent=SPCC18B5.06.1");
    assert_eq!(gene_gff_lines[5],
               "chromosome_3\tPomBase\tCDS\t729266\t729319\t.\t+\t0\tID=SPCC18B5.06.1:exon:2;Parent=SPCC18B5.06.1")
}


use std::cmp::Ordering;
#[cfg(test)]
use std::collections::HashSet;
#[cfg(test)]
use std::num::NonZeroUsize;
#[cfg(test)]
use flexstr::shared_str as flex_str;
#[cfg(test)]
fn make_test_gene() -> GeneDetails {
    GeneDetails {
        uniquename: flex_str!("SPCC18B5.06"),
        name: Some(flex_str!("dom34")),
        taxonid: 4896,
        product: Some(flex_str!("Dom34-Hbs1 translation release factor complex subunit, peloto ortholog Dom34")),
        deletion_viability: DeletionViability::Viable,
        uniprot_identifier: Some(flex_str!("Q9USL5")),
        secondary_identifier: None,
        biogrid_interactor_id: Some(275937),
        rnacentral_urs_identifier: None,
        pdb_entries: vec![],
        interpro_matches: vec![],
        tm_domain_coords: vec![],
        disordered_region_coords: vec![],
        low_complexity_region_coords: vec![],
        coiled_coil_coords: vec![],
        has_protein_features: false,
        rfam_annotations: vec![],
        orfeome_identifier: None,
        pombephosphoproteomics_unige_ch_gene: None,
        name_descriptions: vec![],
        synonyms: vec![],
        dbxrefs: HashSet::new(),
        gocam_ids: HashSet::new(),

        flags: HashSet::new(),

        feature_type: flex_str!("mRNA gene"),
        feature_so_termid: flex_str!("SO:0000704"),
        transcript_so_termid: Some(flex_str!("SO:0001217")),
        characterisation_status: Some(flex_str!("published")),
        taxonomic_distribution: Some(flex_str!("dubious")),
        location: Some(ChromosomeLocation {
            chromosome_name: flex_str!("chromosome_3"),
            start_pos: 729_054,
            end_pos: 730_829,
            strand: Strand::Forward,
            phase: None,
        }),
        gene_neighbourhood: vec![],
        transcripts: vec![
            flex_str!("SPCC18B5.06.1"),
        ],
        transcripts_by_uniquename: std::collections::HashMap::from([
            (flex_str!("SPCC18B5.06.1"),
             Some(TranscriptDetails {
                uniquename: flex_str!("SPCC18B5.06.1"),
                name: Some(flex_str!("dom34.1")),
                location: ChromosomeLocation {
                    chromosome_name: flex_str!("chromosome_3"),
                    start_pos: 729_054,
                    end_pos: 730_829,
                    strand: Strand::Forward,
                    phase: None,
                },
                parts: vec![
                    FeatureShort {
                        feature_type: FeatureType::FivePrimeUtr,
                        uniquename: flex_str!("SPCC18B5.06.1:five_prime_UTR:1"),
                        name: None,
                        location: ChromosomeLocation {
                            chromosome_name: flex_str!("chromosome_3"),
                            start_pos: 729_054,
                            end_pos: 729_132,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: flex_str!("AAAATTTCGAGATATTATGTCAGTCAAAATAGCCTAAAAAATTCCTGTTCACTTAAATTCTTCGTCAACCACATTCAAT"),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: flex_str!("SPCC18B5.06.1:exon:1"),
                        name: None,
                        location: ChromosomeLocation {
                            chromosome_name: flex_str!("chromosome_3"),
                            start_pos: 729_133,
                            end_pos: 729_212,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: flex_str!("ATGAAGTTGATTCAAAAAAACATCGAAAAAAATGGCTCCGGATGGATAACCATGTGCCCTGAAGAGCCAGAAGATATGTG"),
                    },
                    FeatureShort {
                        feature_type: FeatureType::CdsIntron,
                        uniquename: flex_str!("SPCC18B5.06.1:intron:1"),
                        name: None,
                        location: ChromosomeLocation {
                            chromosome_name: flex_str!("chromosome_3"),
                            start_pos: 729_213,
                            end_pos: 729_265,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: flex_str!("GTATGCCTGGTTGCTGTATATCTGCAATCAATGCACAAACTAACTTAGTTTAG"),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: flex_str!("SPCC18B5.06.1:exon:2"),
                        name: None,
                        location: ChromosomeLocation {
                            chromosome_name: flex_str!("chromosome_3"),
                            start_pos: 729_266,
                            end_pos: 729_319,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: flex_str!("GCATTTGTATAATATTCTTCAAGTTGGAGATCAGCTGAAGGCTTCTACAGTTCG"),
                    },
                    FeatureShort {
                        feature_type: FeatureType::CdsIntron,
                        uniquename: flex_str!("SPCC18B5.06.1:intron:2"),
                        name: None,
                        location: ChromosomeLocation {
                            chromosome_name: flex_str!("chromosome_3"),
                            start_pos: 729_320,
                            end_pos: 729_379,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: flex_str!("GTATGAAATAAAGTGTATCTTCAATGTTCGATAATAACACCTGTTTTTTCTGTTGTTTAG"),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: flex_str!("SPCC18B5.06.1:exon:3"),
                        name: None,
                        location: ChromosomeLocation {
                            chromosome_name: flex_str!("chromosome_3"),
                            start_pos: 729_380,
                            end_pos: 729_678,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: flex_str!("TCGTGTAGTGAAAGTTGGCGCTACAGGAAGTACGTCAGGTTCAAGAGTTGTGATGAAACTACGTATTTTAGTTGAGAATATGGACTTTGATACAAAGGCTGCTCAATTGCACATCAAAGGACGGACAACTGAATACCATCCTGAAGTTAAGATGGGATCCTACCATACCTTGGACTTAGAACTACATCGCAATTTTACTCTATATAAAAATGAATGGGATGCATTTGCATTGGACCGTGTAGATGCTGCTTGTAATCCTTCAAGAAATGCTGAAATAGGTGCTGTGGTTTTAGATGAAG"),
                    },
                    FeatureShort {
                        feature_type: FeatureType::CdsIntron,
                        uniquename: flex_str!("SPCC18B5.06.1:intron:3"),
                        name: None,
                        location: ChromosomeLocation {
                            chromosome_name: flex_str!("chromosome_3"),
                            start_pos: 729_679,
                            end_pos: 729_734,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: flex_str!("GTAAATTTCTTGGCTTACATCGCTTTATTAGTTCGTGCATGTTTACTGACTTGCAG"),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: flex_str!("SPCC18B5.06.1:exon:4"),
                        name: None,
                        location: ChromosomeLocation {
                            chromosome_name: flex_str!("chromosome_3"),
                            start_pos: 729_735,
                            end_pos: 729_841,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: flex_str!("GTCTTGCCAATATTTGTCTCATTACAGATTATATGACCATCCTGCGTCAAAGGATTGATCAAGTGATTCCAAGGAAACGGAGAGGGGACAGCAGCGCTTACCAAAAG"),
                    },
                    FeatureShort {
                        feature_type: FeatureType::CdsIntron,
                        uniquename: flex_str!("SPCC18B5.06.1:intron:4"),
                        name: None,
                        location: ChromosomeLocation {
                            chromosome_name: flex_str!("chromosome_3"),
                            start_pos: 729_842,
                            end_pos: 729_889,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: flex_str!("GTAAATTTTTAGACTTTGATTTTTCGTCCGTACTGATCAATTTTTTAG"),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: flex_str!("SPCC18B5.06.1:exon:5"),
                        name: None,
                        location: ChromosomeLocation {
                            chromosome_name: flex_str!("chromosome_3"),
                            start_pos: 729_890,
                            end_pos: 730_522,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: flex_str!("GGCCTTGATAAATTTTATGACTCTGTTTTTCAATCCATTAACTCAGAATTTGATTTTGATAAATTGAAAGTCGTTATTCTTGCTTCACCAGGATTTGTGGCTCGAGGCTTGTATGACTACATATTCAGCATGGCCGTGAAGTTAGACTTGAAACAAATTGTTAAATCAAAGAATAAATTTGTCATCCTTCATTCTAGCACTGGTCATATTCATTCCCTTAATGAAATTTTGAAGGACCCTGCTGTTGAATCAAAACTAGCCGACACAAAATACGTACAAGAAATTCGCGTTCTGAATAAATTTTACGATGTCATGAATGAAGATGATAGAAAGGCATGGTATGGTCCAAATCATGTTTTGAAGGCTTTTGAACTTGGCGCGATCGGAGAACTTCTGATTAGCGATTCTCTGTTCAGGAGTTCTGACATTGCTACTAGAAAAAAATGGGTTTCATTAGTAGAAGGTGTTAAGGAGATTAACTGTCCTGTTTATATTTTCAGTAGTTTGCATGAGTCAGGGAAGCAGCTGGATCTGTTGTCAGGTATTGCCGCCATTCTCACTTACCCAGTCGATGAAGAGGATATATCAGAAGATGAAGAGGATGAGGAATCCCAAAATTTTGAACATAGTTAA"),
                    },
                    FeatureShort {
                        feature_type: FeatureType::ThreePrimeUtr,
                        uniquename: flex_str!("SPCC18B5.06.1:three_prime_UTR:1"),
                        name: None,
                        location: ChromosomeLocation {
                            chromosome_name: flex_str!("chromosome_3"),
                            start_pos: 730_523,
                            end_pos: 730_829,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: flex_str!("AGTTCATCAGTATCCGAATTGTCATGAATCTAATTATTGCTAAGCCAATATTTCATACTTTAAGCTCGGTTAGAACAATTTGTTTCATTCTCTTAAAAAATTTATTTATGGGCTCGTTTTGGTAGTCATTATTTATGCTTTTACTTGGATGTTTTAGGGTATTTACTATGATAAACATGCAAAAAATTAGGTGTTAGAATGGTCAAAAATTGATACCCTAAAATGATTTATACTTATCGATTATAACTCTTAACTTGTAAAATTAAGCTGTTAATTATAGCCGGTCCAATACAGTATTCAATTACAG"),
                    }
                ],
                transcript_type: flex_str!("mRNA"),
                protein: Some(ProteinDetails {
                    uniquename: flex_str!("SPCC18B5.06.1:pep"),
                    sequence: flex_str!("MKLIQKNIEKNGSGWITMCPEEPEDMWHLYNILQVGDQLKASTVRRVVKVGATGSTSGSRVVMKLRILVENMDFDTKAAQLHIKGRTTEYHPEVKMGSYHTLDLELHRNFTLYKNEWDAFALDRVDAACNPSRNAEIGAVVLDEGLANICLITDYMTILRQRIDQVIPRKRRGDSSAYQKGLDKFYDSVFQSINSEFDFDKLKVVILASPGFVARGLYDYIFSMAVKLDLKQIVKSKNKFVILHSSTGHIHSLNEILKDPAVESKLADTKYVQEIRVLNKFYDVMNEDDRKAWYGPNHVLKAFELGAIGELLISDSLFRSSDIATRKKWVSLVEGVKEINCPVYIFSSLHESGKQLDLLSGIAAILTYPVDEEDISEDEEDEESQNFEHS*"),
                    number_of_residues: 390,
                    product: None,
                    molecular_weight: 44.3361,
                    average_residue_weight: 0.11339,
                    charge_at_ph7: -8.45,
                    isoelectric_point: 5.64,
                    codon_adaptation_index: 0.59
                }),
                cds_location: Some(ChromosomeLocation {
                    chromosome_name: flex_str!("chromosome_3"),
                    start_pos: 729_133,
                    end_pos: 730_522,
                    strand: Strand::Forward,
                    phase: None,
                }),
                gene_uniquename: flex_str!("SPCC18B5.06"),
                rna_seq_length_spliced: NonZeroUsize::new(1170),
                rna_seq_length_unspliced: NonZeroUsize::new(730_829 - 729_054),
            })
        )]),
        cv_annotations: std::collections::HashMap::new(),
        physical_interactions: vec![],
        genetic_interactions: std::collections::HashMap::new(),
        ortholog_annotations: vec![],
        paralog_annotations: vec![],
        target_of_annotations: vec![],
        references_by_uniquename: std::collections::HashMap::new(),
        genes_by_uniquename: std::collections::HashMap::new(),
        genotypes_by_uniquename: std::collections::HashMap::new(),
        alleles_by_uniquename: std::collections::HashMap::new(),
        terms_by_termid: std::collections::HashMap::new(),
        annotation_details: std::collections::HashMap::new(),
        feature_publications: HashSet::new(),
        subset_termids: HashSet::new(),
        gene_history: vec![],
    }
}

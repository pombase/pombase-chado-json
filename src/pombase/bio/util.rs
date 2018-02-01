use std::collections::{HashSet, HashMap};
use web::data::{GeneDetails, ChromosomeLocation, DeletionViability, Strand,
                ChromosomeShort, TranscriptDetails, FeatureShort,
                ProteinDetails, FeatureType};

pub fn format_fasta(id: &str, maybe_desc: Option<String>,
                    seq: &str, width: usize) -> String {
    let mut ret = ">".to_owned() + id;

    if let Some(desc) = maybe_desc {
        ret.push(' ');
        ret.push_str(&desc);
    }

    ret.push('\n');

    if seq.len() == 0 {
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
        } else {
            "."
        };
    let start_pos_str = format!("{}", location.start_pos);
    let end_pos_str = format!("{}", location.end_pos);
    let mut ret_val = String::new();
    let bits: Vec<&str> =
        vec![&chromosome_export_id, source, feat_type,
             &start_pos_str, &end_pos_str,
             ".", location.strand.to_gff_str(), phase_char];
    for s in bits {
        ret_val.push_str(s);
        ret_val.push_str("\t");
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


pub fn format_gene_gff(chromosome_export_id: &str,
                       source: &str, gene: &GeneDetails) -> Vec<String> {
    let mut ret_val = vec![];
    if let Some (ref gene_loc) = gene.location {
        let maybe_gene_name =
            if let Some(gene_name) = gene.name.as_ref() {
                Some(gene_name as &str)
            } else {
                None
            };

        ret_val.push(to_gff(chromosome_export_id, source, &gene.uniquename, maybe_gene_name,
                            "gene", &gene_loc, None));
        for tr in &gene.transcripts {
            ret_val.push(to_gff(chromosome_export_id,
                                source, &tr.uniquename, None, &tr.transcript_type,
                                &tr.location, Some(&gene.uniquename)));
            for part in &tr.parts {
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
                                    &part.location, Some(&tr.uniquename)));
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
    let gene_gff_lines = format_gene_gff("chromosome_3", "PomBase", &gene);

    assert_eq!(gene_gff_lines.len(), 13);
    assert_eq!(gene_gff_lines[0],
               "chromosome_3\tPomBase\tgene\t729054\t730829\t.\t+\t.\tID=SPCC18B5.06;Name=dom34");
    assert_eq!(gene_gff_lines[1],
               "chromosome_3\tPomBase\tmRNA\t729054\t730829\t.\t+\t.\tID=SPCC18B5.06.1;Parent=SPCC18B5.06");
    assert_eq!(gene_gff_lines[2],
               "chromosome_3\tPomBase\tfive_prime_UTR\t729054\t729132\t.\t+\t.\tID=SPCC18B5.06.1:five_prime_UTR:1;Parent=SPCC18B5.06.1");
    assert_eq!(gene_gff_lines[3],
               "chromosome_3\tPomBase\tCDS\t729133\t729212\t.\t+\t.\tID=SPCC18B5.06.1:exon:1;Parent=SPCC18B5.06.1");
    assert_eq!(gene_gff_lines[4],
               "chromosome_3\tPomBase\tintron\t729213\t729265\t.\t+\t.\tID=SPCC18B5.06.1:intron:1;Parent=SPCC18B5.06.1");
    assert_eq!(gene_gff_lines[5],
               "chromosome_3\tPomBase\tCDS\t729266\t729319\t.\t+\t.\tID=SPCC18B5.06.1:exon:2;Parent=SPCC18B5.06.1")
}

#[allow(dead_code)]
fn make_test_gene() -> GeneDetails {
    GeneDetails {
        uniquename: "SPCC18B5.06".to_owned(),
        name: Some("dom34".to_owned()),
        taxonid: 4896,
        product: Some("Dom34-Hbs1 translation release factor complex subunit, peloto ortholog Dom34".to_owned()),
        deletion_viability: DeletionViability::Viable,
        uniprot_identifier: Some("Q9USL5".to_owned()),
        interpro_matches: vec![],
        tm_domain_coords: vec![],
        orfeome_identifier: None,
        name_descriptions: vec![],
        synonyms: vec![],
        dbxrefs: HashSet::new(),
        feature_type: "mRNA gene".to_owned(),
        characterisation_status: Some("published".to_owned()),
        location: Some(ChromosomeLocation {
            chromosome: ChromosomeShort {
                name: "chromosome_3".to_owned(),
                length: 2452883,
                ena_identifier: "CU329672.1".to_owned(),
            },
            start_pos: 729054,
            end_pos: 730829,
            strand: Strand::Forward,
            phase: None,
        }),
        gene_neighbourhood: vec![],
        transcripts: vec![
            TranscriptDetails {
                uniquename:"SPCC18B5.06.1".to_owned(),
                location: ChromosomeLocation {
                    chromosome: ChromosomeShort {
                        name: "chromosome_3".to_owned(),
                        length: 2452883,
                        ena_identifier: "CU329672.1".to_owned(),
                    },
                    start_pos: 729054,
                    end_pos: 730829,
                    strand: Strand::Forward,
                    phase: None,
                },
                parts: vec![
                    FeatureShort {
                        feature_type: FeatureType::FivePrimeUtr,
                        uniquename: "SPCC18B5.06.1:five_prime_UTR:1".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729054,
                            end_pos: 729132,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "AAAATTTCGAGATATTATGTCAGTCAAAATAGCCTAAAAAATTCCTGTTCACTTAAATTCTTCGTCAACCACATTCAAT".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: "SPCC18B5.06.1:exon:1".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729133,
                            end_pos: 729212,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "ATGAAGTTGATTCAAAAAAACATCGAAAAAAATGGCTCCGGATGGATAACCATGTGCCCTGAAGAGCCAGAAGATATGTG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::CdsIntron,
                        uniquename: "SPCC18B5.06.1:intron:1".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729213,
                            end_pos: 729265,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GTATGCCTGGTTGCTGTATATCTGCAATCAATGCACAAACTAACTTAGTTTAG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: "SPCC18B5.06.1:exon:2".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729266,
                            end_pos: 729319,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GCATTTGTATAATATTCTTCAAGTTGGAGATCAGCTGAAGGCTTCTACAGTTCG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::CdsIntron,
                        uniquename: "SPCC18B5.06.1:intron:2".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729320,
                            end_pos: 729379,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GTATGAAATAAAGTGTATCTTCAATGTTCGATAATAACACCTGTTTTTTCTGTTGTTTAG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: "SPCC18B5.06.1:exon:3".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729380,
                            end_pos: 729678,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "TCGTGTAGTGAAAGTTGGCGCTACAGGAAGTACGTCAGGTTCAAGAGTTGTGATGAAACTACGTATTTTAGTTGAGAATATGGACTTTGATACAAAGGCTGCTCAATTGCACATCAAAGGACGGACAACTGAATACCATCCTGAAGTTAAGATGGGATCCTACCATACCTTGGACTTAGAACTACATCGCAATTTTACTCTATATAAAAATGAATGGGATGCATTTGCATTGGACCGTGTAGATGCTGCTTGTAATCCTTCAAGAAATGCTGAAATAGGTGCTGTGGTTTTAGATGAAG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::CdsIntron,
                        uniquename: "SPCC18B5.06.1:intron:3".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729679,
                            end_pos: 729734,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GTAAATTTCTTGGCTTACATCGCTTTATTAGTTCGTGCATGTTTACTGACTTGCAG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: "SPCC18B5.06.1:exon:4".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729735,
                            end_pos: 729841,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GTCTTGCCAATATTTGTCTCATTACAGATTATATGACCATCCTGCGTCAAAGGATTGATCAAGTGATTCCAAGGAAACGGAGAGGGGACAGCAGCGCTTACCAAAAG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::CdsIntron,
                        uniquename: "SPCC18B5.06.1:intron:4".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729842,
                            end_pos: 729889,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GTAAATTTTTAGACTTTGATTTTTCGTCCGTACTGATCAATTTTTTAG".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::Exon,
                        uniquename: "SPCC18B5.06.1:exon:5".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 729890,
                            end_pos: 730522,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "GGCCTTGATAAATTTTATGACTCTGTTTTTCAATCCATTAACTCAGAATTTGATTTTGATAAATTGAAAGTCGTTATTCTTGCTTCACCAGGATTTGTGGCTCGAGGCTTGTATGACTACATATTCAGCATGGCCGTGAAGTTAGACTTGAAACAAATTGTTAAATCAAAGAATAAATTTGTCATCCTTCATTCTAGCACTGGTCATATTCATTCCCTTAATGAAATTTTGAAGGACCCTGCTGTTGAATCAAAACTAGCCGACACAAAATACGTACAAGAAATTCGCGTTCTGAATAAATTTTACGATGTCATGAATGAAGATGATAGAAAGGCATGGTATGGTCCAAATCATGTTTTGAAGGCTTTTGAACTTGGCGCGATCGGAGAACTTCTGATTAGCGATTCTCTGTTCAGGAGTTCTGACATTGCTACTAGAAAAAAATGGGTTTCATTAGTAGAAGGTGTTAAGGAGATTAACTGTCCTGTTTATATTTTCAGTAGTTTGCATGAGTCAGGGAAGCAGCTGGATCTGTTGTCAGGTATTGCCGCCATTCTCACTTACCCAGTCGATGAAGAGGATATATCAGAAGATGAAGAGGATGAGGAATCCCAAAATTTTGAACATAGTTAA".to_owned(),
                    },
                    FeatureShort {
                        feature_type: FeatureType::ThreePrimeUtr,
                        uniquename: "SPCC18B5.06.1:three_prime_UTR:1".to_owned(),
                        location: ChromosomeLocation {
                            chromosome: ChromosomeShort {
                                name: "chromosome_3".to_owned(),
                                length: 2452883,
                                ena_identifier: "CU329672.1".to_owned(),
                            },
                            start_pos: 730523,
                            end_pos: 730829,
                            strand: Strand::Forward,
                            phase: None,
                        },
                        residues: "AGTTCATCAGTATCCGAATTGTCATGAATCTAATTATTGCTAAGCCAATATTTCATACTTTAAGCTCGGTTAGAACAATTTGTTTCATTCTCTTAAAAAATTTATTTATGGGCTCGTTTTGGTAGTCATTATTTATGCTTTTACTTGGATGTTTTAGGGTATTTACTATGATAAACATGCAAAAAATTAGGTGTTAGAATGGTCAAAAATTGATACCCTAAAATGATTTATACTTATCGATTATAACTCTTAACTTGTAAAATTAAGCTGTTAATTATAGCCGGTCCAATACAGTATTCAATTACAG".to_owned(),
                    }
                ],
                transcript_type: "mRNA".to_owned(),
                protein: Some(ProteinDetails {
                    uniquename: "SPCC18B5.06.1:pep".to_owned(),
                    sequence: "MKLIQKNIEKNGSGWITMCPEEPEDMWHLYNILQVGDQLKASTVRRVVKVGATGSTSGSRVVMKLRILVENMDFDTKAAQLHIKGRTTEYHPEVKMGSYHTLDLELHRNFTLYKNEWDAFALDRVDAACNPSRNAEIGAVVLDEGLANICLITDYMTILRQRIDQVIPRKRRGDSSAYQKGLDKFYDSVFQSINSEFDFDKLKVVILASPGFVARGLYDYIFSMAVKLDLKQIVKSKNKFVILHSSTGHIHSLNEILKDPAVESKLADTKYVQEIRVLNKFYDVMNEDDRKAWYGPNHVLKAFELGAIGELLISDSLFRSSDIATRKKWVSLVEGVKEINCPVYIFSSLHESGKQLDLLSGIAAILTYPVDEEDISEDEEDEESQNFEHS*".to_owned(),
                    molecular_weight: 44.3361,
                    average_residue_weight: 0.11339,
                    charge_at_ph7: -8.45,
                    isoelectric_point: 5.64,
                    codon_adaptation_index: 0.59
                }),
                cds_location: Some(ChromosomeLocation {
                    chromosome: ChromosomeShort {
                        name: "chromosome_3".to_owned(),
                        length: 2452883,
                        ena_identifier: "CU329672.1".to_owned(),
                    },
                    start_pos: 729133,
                    end_pos: 730522,
                    strand: Strand::Forward,
                    phase: None,
                }),
            },
        ],
        cv_annotations: HashMap::new(),
        physical_interactions: vec![],
        genetic_interactions: vec![],
        ortholog_annotations: vec![],
        paralog_annotations: vec![],
        target_of_annotations: vec![],
        references_by_uniquename: HashMap::new(),
        genes_by_uniquename: HashMap::new(),
        genotypes_by_uniquename: HashMap::new(),
        alleles_by_uniquename: HashMap::new(),
        terms_by_termid: HashMap::new(),
    }
}

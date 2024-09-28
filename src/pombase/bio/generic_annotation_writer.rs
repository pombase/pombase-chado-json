use std::io::{self, BufWriter, Write};

use crate::uniprot::UniProtDataMap;

use chrono::prelude::{Local, DateTime};

pub struct UniProtTermidMap {
  pub glycosylation_site_termid: String,
  pub disulphide_bond_termid: String,
}

pub fn write_from_uniprot_map(uniprot_data_map: &UniProtDataMap,
                              uniprot_pmid: &str,
                              termid_map: &UniProtTermidMap,
                              assigned_by: &str,
                              out: &mut dyn Write)
   -> Result<(), io::Error>
{
  let mut writer = BufWriter::new(out);

  let local: DateTime<Local> = Local::now();
  let date = local.format("%F").to_string();

  for uniprot_data in uniprot_data_map.values() {
    for site in &uniprot_data.glycosylation_sites {
      let evidence = site.evidence.as_deref().unwrap_or_default();
      let reference = site.reference.as_deref().unwrap_or(uniprot_pmid);
      let termid = &termid_map.glycosylation_site_termid;
      let residue_extension = format!("residue(N{})", site.range);
      write_generic_annotation(&mut writer,
                               &uniprot_data.gene_uniquename,
                               "",
                               termid,
                               evidence,
                               reference,
                               &date,
                               "",
                               &residue_extension,
                               assigned_by)?;
    }
    for site in &uniprot_data.disulfide_bonds {
      let evidence = site.evidence.as_deref().unwrap_or_default();
      let reference = site.reference.as_deref().unwrap_or(uniprot_pmid);
      let termid = &termid_map.disulphide_bond_termid;
      let residue_extension = format!("residue({})", site.range);
      write_generic_annotation(&mut writer,
                               &uniprot_data.gene_uniquename,
                               "",
                               termid,
                               evidence,
                               reference,
                               &date,
                               "",
                               &residue_extension,
                               assigned_by)?;
    }
    for site in &uniprot_data.lipidation_sites {
      let evidence = site.evidence.as_deref().unwrap_or_default();
      let reference = site.reference.as_deref().unwrap_or(uniprot_pmid);
      let termid = &site.termid;
      let residue_extension = format!("residue({})", site.range);
      write_generic_annotation(&mut writer,
                               &uniprot_data.gene_uniquename,
                               "",
                               termid,
                               evidence,
                               reference,
                               &date,
                               "",
                               &residue_extension,
                               assigned_by)?;
    }
    for site in &uniprot_data.modified_residues {
      let evidence = site.evidence.as_deref().unwrap_or_default();
      let reference = site.reference.as_deref().unwrap_or(uniprot_pmid);
      let termid = &site.termid;
      let residue_abbrev =
          match termid.as_str() {
              "MOD:00042"
                  => "D",
              "MOD:00047"
                  => "T",
              "MOD:00048" | "MOD:00156"
                  => "Y",
              "MOD:00050"
                  => "A",
              "MOD:00046" | "MOD:00060" | "MOD:00159" | "MOD:00171"
                  => "S",
              "MOD:00080"
                  => "Q",
              "MOD:00111" | "MOD:00114" | "MOD:00115" | "MOD:00234" |
              "MOD:00257" | "MOD:00441" | "MOD:00689"
                  => "C",
              "MOD:00064" | "MOD:00083" | "MOD:00084" | "MOD:00085" |
              "MOD:00123" | "MOD:00125" | "MOD:00126" | "MOD:00127" |
              "MOD:00128"
                  => "K",
              "MOD:00153" | "MOD:00226" | "MOD:00890"
                  => "H",
              "MOD:00304"
                  => "L",
              "MOD:00310"
                  => "R",
              "MOD:00068" | "MOD:00351" | "MOD:01625" | "MOD:01982"
                  => "G",
              _ => "",
          };
      let residue_extension = format!("residue({}{})", residue_abbrev, site.range);
      write_generic_annotation(&mut writer,
                               &uniprot_data.gene_uniquename,
                               "",
                               termid,
                               evidence,
                               reference,
                               &date,
                               "",
                               &residue_extension,
                               assigned_by)?;
    }
  }

  Ok(())
}

fn write_generic_annotation(out: &mut BufWriter<&mut dyn Write>,
                            uniquename: &str, name: &str, termid: &str,
                            evidence: &str,
                            reference_uniquename: &str,
                            date: &str,
                            qualifiers: &str,
                            extension: &str,
                            assigned_by: &str)
         -> Result<(), io::Error>
{
  writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
           uniquename, name, termid, evidence,
           reference_uniquename, date, qualifiers, extension, assigned_by)?;
  Ok(())
}

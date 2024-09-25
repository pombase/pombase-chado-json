use std::io::{self, BufWriter, Write};

use crate::uniprot::UniProtDataMap;

use chrono::prelude::{Local, DateTime};

pub struct UniProtTermidMap {
  pub glycosylation_site_termid: String,
  pub disulphide_bond_termid: String,
}

pub fn write_from_uniprot_map(uniprot_data_map: &UniProtDataMap,
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
      let reference = site.reference.as_deref().unwrap_or_default();
      let termid = &termid_map.glycosylation_site_termid;
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
    for site in &uniprot_data.disulfide_bonds {
      let evidence = site.evidence.as_deref().unwrap_or_default();
      let reference = site.reference.as_deref().unwrap_or_default();
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
      let reference = site.reference.as_deref().unwrap_or_default();
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
      let reference = site.reference.as_deref().unwrap_or_default();
      let termid = &site.termid;
      let residue_abbrev =
          if termid == "MOD:00047" {
            "T"
          } else {
            if termid == "MOD:00046" {
              "S"
            } else {
              ""
            }
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

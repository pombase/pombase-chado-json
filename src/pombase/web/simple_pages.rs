use pombase_rc_string::RcString;

use crate::web::config::Config;
use crate::web::data::{GeneDetails, OntAnnotationId};


fn make_title(config: &Config, gene_details: &GeneDetails) -> String {
    let name_and_uniquename =
        if let Some(ref name) = gene_details.name {
            format!("{} ({})", name, gene_details.uniquename)
        } else {
            format!("{}", gene_details.uniquename)
        };

    let feature_type =
        if gene_details.feature_type.ends_with(" gene") {
            String::from(gene_details.feature_type.as_str())
        } else {
            format!("{} gene", gene_details.feature_type)
        };

    if let Some(ref product) = gene_details.product {
        format!("{} - {} {} - {}", config.database_name, feature_type,
                name_and_uniquename, product)
    } else {
        format!("{} - {} {}", config.database_name, feature_type, name_and_uniquename)
    }
}

fn gene_header(config: &Config, title: &str) -> String  {
    format!(r##"
  <meta charset="utf-8">
  <title>{}</title>
  <base href="/">
  <meta http-equiv="X-UA-Compatible" content="IE=edge">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="title" content="{}">
  <meta name="author" content="{}">
  <meta name="description" content="{}">
  <meta name="google-site-verification" content="H8vhCdG7XsSbh-8f20VSipvu2PnZ22YSGZm81jE0-Pk" />
  <link rel="shortcut icon" href="/assets/favicon.ico">
  <meta name="apple-mobile-web-app-title" content="PomBase">
  <meta name="application-name" content="PomBase">
"##,
            title, title, config.database_name, title
    )
}

fn gene_summary(config: &Config, gene_details: &GeneDetails) -> String {
    let mut summ = String::new();

    summ += "<dl>\n";
    if let Some(ref name) = gene_details.name {
        summ += &format!("  <dt>Gene standard name</dt> <dd>{}</dd>\n",  name);
    }
    summ += &format!("  <dt>Systematic ID</dt> <dd>{}</dd>\n", gene_details.uniquename);

    if gene_details.synonyms.len() > 0 {
        let synonyms: Vec<RcString> =
            gene_details.synonyms.iter().map(|s| s.name.clone()).collect();
        summ += &format!("  <dt>Synonyms</dt> <dd>{}</dd>\n", synonyms.join(", "));
    }

    if let Some(ref uniprot_identifier) = gene_details.uniprot_identifier {
        summ += &format!("  <dt>UniProt ID</dt> <dd>{}</dd>\n", uniprot_identifier);
    }

    if let Some(ref orfeome_identifier) = gene_details.orfeome_identifier {
        summ += &format!("  <dt>ORFeome ID</dt> <dd>{}</dd>\n", orfeome_identifier);
    }

    if let Some(gene_organism) = config.organism_by_taxonid(gene_details.taxonid) {
        summ += &format!("  <dt>Organism</dt> <dd>{}\n", gene_organism.scientific_name());

        if gene_organism.alternative_names.len() > 0 {
            summ += &format!(" ({})", gene_organism.alternative_names.join(", "));
        }
        summ += "</dd>\n";
    }

    if let Some(ref characterisation_status) = gene_details.characterisation_status {
        summ += &format!("  <dt>Characterisation status</dt> <dd>{}\n",
                         characterisation_status);
    }

    summ += "</dl>\n";

    summ
}

fn get_annotations(ont_annotation_ids: &Vec<OntAnnotationId>,
                   gene_details: &GeneDetails) -> String {
    let mut ret = String::new();
    let mut references = vec![];
    for ont_annotation_id in ont_annotation_ids {
        if let Some(annotation_details) = gene_details.annotation_details.get(ont_annotation_id) {
            if let Some(ref reference) = annotation_details.reference {
                references.push(reference);
            }
        }
    }

    if references.len() > 0 {
        references.sort();
        references.dedup();
        ret += "<p>References:</p>\n<ul>";
        for reference in references {
            ret += &format!("<li><a href='/reference/{}'>{}</a></li>\n", reference, reference);
        }
        ret += "</ul>\n";
    }

    ret
}

fn gene_annotation(config: &Config, gene_details: &GeneDetails) -> String {
    let mut annotation_html = String::new();

    let mut cv_names: Vec<RcString> = vec![];

    for cv_name in gene_details.cv_annotations.keys() {
        cv_names.push(cv_name.clone());
    }

    let cmp_display_names =
        |cv_name_a: &RcString, cv_name_b: &RcString| {
            let cv_config_a = config.cv_config_by_name(&cv_name_a);
            let cv_display_name_a = cv_config_a.display_name;
            let cv_config_b = config.cv_config_by_name(&cv_name_b);
            let cv_display_name_b = cv_config_b.display_name;

            cv_display_name_a.cmp(&cv_display_name_b)
        };

    cv_names.sort_by(cmp_display_names);

    for cv_name in cv_names {
        let term_annotations = gene_details.cv_annotations.get(&cv_name).unwrap();
        let cv_config = config.cv_config_by_name(&cv_name);
        let cv_display_name = cv_config.display_name;
        annotation_html += &format!("<sect>\n<h3>{}</h3>\n", cv_display_name);
        for term_annotation in term_annotations {
            let term_name =
                if let Some(term_short_opt) =
                    gene_details.terms_by_termid.get(&term_annotation.term) {
                        if let Some(term_short) = term_short_opt {
                            String::from(term_short.name.as_str())
                        } else {
                            String::from("(unknown term name)")
                        }
                    } else {
                        String::from("(unknown term name)")
                    };
            annotation_html += &format!("<sect><h4><a href='/term/{}'>{}</a> - \
<a href='/term/{}'>{}</a></h4>\n{}</sect>\n",
                                        term_annotation.term,
                                        term_annotation.term,
                                        term_annotation.term,
                                        term_name,
                                        get_annotations(&term_annotation.annotations,
                                                        gene_details));
        }
        annotation_html += "</sect>\n";
    }

    annotation_html
}

fn orthologs(config: &Config, gene_details: &GeneDetails) -> String {
    let mut orth_html = String::new();

    if gene_details.ortholog_annotations.is_empty() {
        return orth_html;
    }

    orth_html += "<sect><h2>Orthologs</h2>\n<ul>\n";

    for orth_annotation in &gene_details.ortholog_annotations {
        let maybe_orth_org = config.organism_by_taxonid(orth_annotation.ortholog_taxonid);
        let orth_org_scientific_name =
            if let Some(orth_org) = maybe_orth_org {
                orth_org.scientific_name()
            } else {
                String::from("Unknown organism")
            };
        if let Some(maybe_orth_gene_short) =
            gene_details.genes_by_uniquename.get(&orth_annotation.ortholog_uniquename) {
                if let Some(orth_gene_short) = maybe_orth_gene_short {
                    orth_html += &format!("<li>{} {}</li>\n", orth_gene_short.display_name(),
                                          orth_org_scientific_name);
                }
            }
    }

    orth_html += "</ul></sect>\n";

    orth_html
}

fn gene_body(config: &Config, title: &str, gene_details: &GeneDetails) -> String {
    let mut body = String::new();

    body += &format!("<h1>{}</h1>\n", title);

    body += &format!("<sect><h2>Welcome to PomBase</h2>\n<p>{}</p></sect>\n",
                     config.site_description);

    body += &format!("<sect><h2>Gene summary</h2>\n{}</sect>\n",
                     gene_summary(config, gene_details));

    body += &format!("<sect><h2>Annotation</h2>\n{}</sect>\n",
                     gene_annotation(config, gene_details));

    body += &orthologs(config, gene_details);

    body
}

pub fn render_simple_gene_page(config: &Config, gene_details: &GeneDetails) -> String  {
    let title = make_title(config, gene_details);

    format! ("
<!DOCTYPE html>
<html>
  <head>
    {}
  </head>
  <body>
    {}
  </body>
</html>",
          gene_header(config, &title),
          gene_body(config, &title, gene_details))
}

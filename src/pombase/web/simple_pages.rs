use pombase_rc_string::RcString;

use crate::web::config::Config;
use crate::web::data::{GeneDetails, ReferenceDetails, TermDetails,
                       ContainerType, GenotypeShort,
                       OntAnnotationId, AnnotationContainer, OrthologAnnotationContainer};


fn format_page(header: &str, body: &str) -> String {
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
             header, body)
}

fn make_gene_title(config: &Config, gene_details: &GeneDetails) -> String {
    let name_and_uniquename =
        if let Some(ref name) = gene_details.name {
            format!("{} ({})", name, gene_details.uniquename)
        } else {
            format!("{}", gene_details.uniquename)
        };

    let feature_type =
        if gene_details.feature_type.ends_with("gene") {
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

fn header(config: &Config, title: &str) -> String  {
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
  <meta name="apple-mobile-web-app-title" content="{}">
  <meta name="application-name" content="{}">
"##,
            title, title, config.database_name, title, config.database_name,
            config.database_name
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

fn genotype_display_name(genotype_short: &GenotypeShort) -> String {
    let mut bits: Vec<_> = genotype_short.display_uniquename.split("-").collect();
    bits.dedup();
    bits.join(" ")
}

fn get_annotations(ont_annotation_ids: &Vec<OntAnnotationId>,
                   container: &dyn AnnotationContainer) -> String {
    let mut ret = String::new();
    let mut references = vec![];
    let mut genes = vec![];
    let mut genotypes = vec![];

    for ont_annotation_id in ont_annotation_ids {
        if let Some(annotation_details) = container.annotation_details().get(ont_annotation_id) {
            if let Some(ref reference) = annotation_details.reference {
                references.push(reference);
            }

            if let Some(ref genotype) = annotation_details.genotype {
                genotypes.push(genotype);
            }

            for gene_uniquename in &annotation_details.genes {
                if let Some(maybe_gene_short) =
                    container.genes_by_uniquename().get(gene_uniquename.as_str()) {
                        if let Some(gene_short) = maybe_gene_short {
                            genes.push(gene_short)
                        }
                    }
            }
        }
    }

    if container.container_type() != ContainerType::Reference && references.len() > 0 {
        references.sort();
        references.dedup();
        ret += "<p>References:</p>\n<ul>";
        for reference in references {
            ret += &format!("<li><a href='/reference/{}'>{}</a></li>\n", reference, reference);
        }
        ret += "</ul>\n";
    }

    if container.container_type() != ContainerType::Gene && genes.len() > 0 {
        genes.sort();
        genes.dedup();
        ret += "<p>Genes:</p>\n<ul>";
        for gene in genes {
            ret += &format!("<li><a href='/gene/{}'>{}</a></li>\n", gene.uniquename,
                            gene.display_name());
        }
        ret += "</ul>\n";
    }

    if container.container_type() != ContainerType::Genotype && genotypes.len() > 0 {
        genotypes.sort();
        genotypes.dedup();
        if let Some(genotypes_by_uniquename) = container.genotypes_by_uniquename() {
            ret += "<p>Genotypes:</p>\n<ul>";
            for uniquename in genotypes {
                if let Some(genotype_short) = genotypes_by_uniquename.get(uniquename.as_str()) {
                    ret += &format!("<li><a href='/genotype/{}'>{}</a></li>\n",
                                    uniquename,
                                    genotype_display_name(genotype_short));
                }
            }
            ret += "</ul>\n";
        }
    }

    ret
}

fn annotation_section(config: &Config, container: &dyn AnnotationContainer) -> String {
    let mut annotation_html = String::new();

    let mut cv_names: Vec<RcString> = vec![];

    for cv_name in container.cv_annotations().keys() {
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
        let term_annotations = container.cv_annotations().get(&cv_name).unwrap();
        let cv_config = config.cv_config_by_name(&cv_name);
        let cv_display_name = cv_config.display_name;
        annotation_html += &format!("<sect>\n<h3>{}</h3>\n", cv_display_name);
        for term_annotation in term_annotations {
            let term_name =
                if let Some(term_short_opt) =
                    container.terms_by_termid().get(&term_annotation.term) {
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
                                                        container));
        }
        annotation_html += "</sect>\n";
    }

    annotation_html
}

fn orthologs(config: &Config, container: &dyn OrthologAnnotationContainer) -> String {
    let mut orth_html = String::new();

    if container.ortholog_annotations().is_empty() {
        return orth_html;
    }

    orth_html += "<sect><h2>Orthologs</h2>\n<ul>\n";

    for orth_annotation in container.ortholog_annotations().iter() {
        let maybe_orth_org = config.organism_by_taxonid(orth_annotation.ortholog_taxonid);
        let orth_org_scientific_name =
            if let Some(orth_org) = maybe_orth_org {
                orth_org.scientific_name()
            } else {
                String::from("Unknown organism")
            };
        if let Some(maybe_orth_gene_short) =
            container.genes_by_uniquename().get(&orth_annotation.ortholog_uniquename) {
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

    body += &format!("<sect><h2>Welcome to <a href='/'>{}</a></h2>\n<p>{}</p></sect>\n",
                     config.database_name, config.site_description);

    body += &format!("<sect><h2>Gene summary</h2>\n{}</sect>\n",
                     gene_summary(config, gene_details));

    body += &format!("<sect><h2>Annotation</h2>\n{}</sect>\n",
                     annotation_section(config, gene_details));

    body += &orthologs(config, gene_details);

    body
}

pub fn render_simple_gene_page(config: &Config, gene_details: &GeneDetails) -> String  {
    let title = make_gene_title(config, gene_details);

    format_page(&header(config, &title), &gene_body(config, &title, gene_details))
}

fn make_reference_title(config: &Config, reference_details: &ReferenceDetails) -> String {
    if let Some(ref title) = reference_details.title {
        format!("{} - {} - {}", config.database_name, reference_details.uniquename, title)
    } else {
        format!("{} - {}", config.database_name, reference_details.uniquename)
    }
}

fn reference_summary(reference_details: &ReferenceDetails) -> String {
    let mut summ = String::new();

    summ += "<dl>\n";

    summ += &format!("<dt>PubMed ID</dt> <dd>{}</dd>\n", reference_details.uniquename);

    if let Some(ref title) = reference_details.title {
        summ += &format!("<dt>Title</dt> <dd>{}</dd>\n", title);
    }

    if let Some(ref authors) = reference_details.authors {
        summ += &format!("<dt>Authors</dt> <dd>{}</dd>\n", authors);
    }

    if let Some(ref citation) = reference_details.citation {
        summ += &format!("<dt>Citation</dt> <dd>{}</dd>\n", citation);
    }

    if let Some(ref publication_year) = reference_details.publication_year {
        summ += &format!("<dt>Publication year</dt> <dd>{}</dd>\n", publication_year);
    }

    summ += "</dl>\n";

    summ
}

fn reference_body(config: &Config, title: &str, reference_details: &ReferenceDetails) -> String {
    let mut body = String::new();

    body += &format!("<h1>{}</h1>\n", title);

    body += &format!("<sect><h2>Welcome to <a href='/'>{}</a></h2>\n<p>{}</p></sect>\n",
                     config.database_name, config.site_description);

    body += &format!("<sect><h2>Reference summary</h2>\n{}</sect>\n",
                     reference_summary(reference_details));

    body += &format!("<sect><h2>Annotation</h2>\n{}</sect>\n",
                     annotation_section(config, reference_details));

    body += &orthologs(config, reference_details);

    body
}

pub fn render_simple_reference_page(config: &Config, reference_details: &ReferenceDetails) -> String  {
    let title = make_reference_title(config, reference_details);

    format_page(&header(config, &title), &reference_body(config, &title, reference_details))
}

fn make_term_title(config: &Config, term_details: &TermDetails) -> String {
    let cv_config = config.cv_config_by_name(&term_details.cv_name);
    let cv_display_name = cv_config.display_name;

    format!("{} - {} - {} - {}", config.database_name, term_details.termid,
            term_details.name, cv_display_name)
}

fn term_summary(config: &Config, term_details: &TermDetails) -> String {
    let mut summ = String::new();

    summ += "<dl>\n";

    summ += &format!("<dt>Term ID</dt> <dd>{}</dd>\n", term_details.termid);
    summ += &format!("<dt>Term name</dt> <dd>{}</dd>\n", term_details.name);

    let cv_config = config.cv_config_by_name(&term_details.cv_name);
    let cv_display_name = cv_config.display_name;

    summ += &format!("<dt>CV name</dt> <dd>{}</dd>\n", cv_display_name);

    summ += "</dl>\n";

    summ
}

fn term_body(config: &Config, title: &str, term_details: &TermDetails) -> String {
    let mut body = String::new();

    body += &format!("<h1>{}</h1>\n", title);

    body += &format!("<sect><h2>Welcome to <a href='/'>{}</a></h2>\n<p>{}</p></sect>\n",
                     config.database_name, config.site_description);

    body += &format!("<sect><h2>Term summary</h2>\n{}</sect>\n",
                     term_summary(config, term_details));

    body += &format!("<sect><h2>Annotation</h2>\n{}</sect>\n",
                     annotation_section(config, term_details));

    body
}

pub fn render_simple_term_page(config: &Config, term_details: &TermDetails) -> String  {
    let title = make_term_title(config, term_details);

    format_page(&header(config, &title), &term_body(config, &title, term_details))
}

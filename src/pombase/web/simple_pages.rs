use flexstr::{SharedStr as FlexStr, shared_str as flex_str};

use crate::web::config::Config;
use crate::data_types::{GeneDetails, ReferenceDetails, TermDetails,
                        ContainerType, OntAnnotationId, AnnotationContainer,
                        OrthologAnnotationContainer,
                        Strand, GenotypeDetails, Ploidiness};

fn format_page(header: &str, body: &str) -> String {
    format! ("
<!DOCTYPE html>
<html>
  <head>
    {}
  </head>
  <body style='font-size: 140%;'>
    {}
  </body>
</html>",
             header, body)
}

fn make_gene_title(gene_details: &GeneDetails) -> String {
    let name_and_uniquename =
        if let Some(ref name) = gene_details.name {
            format!("{} ({})", name, gene_details.uniquename)
        } else {
            format!("{}", gene_details.uniquename)
        };

    let feature_type =
        if gene_details.feature_type == "mRNA gene" {
            String::from("protein coding gene")
        } else {
            if gene_details.feature_type.ends_with("gene") {
                String::from(gene_details.feature_type.as_str())
            } else {
                format!("{} gene", gene_details.feature_type)
            }
        };

    if let Some(ref product) = gene_details.product {
        format!("{} - {} - {}", feature_type,
                name_and_uniquename, product)
    } else {
        format!("{} - {}", feature_type, name_and_uniquename)
    }
}

fn make_genotype_title(genotype_details: &GenotypeDetails) -> String {
    let genotype_and_ploidiness =
        if genotype_details.ploidiness == Ploidiness::Haploid {
            "Haploid genotype"
        } else {
            "Diploid genotype"
        };

    format!("{} - {}", genotype_and_ploidiness, genotype_details.display_name)
}

fn head(config: &Config, title: &str) -> String  {
    let title = format!("{} - {} - {}", config.database_name, title, config.database_long_name);

    let style = "
body {
  padding: 1em 5em;
}
dd {
  margin-left: 2em;
}
tbody:nth-child(2n+1) {
  background-color: #f3f3f3;
}
th, td {
  border: 1px solid lightgrey;
  padding: 0.2em 0.4em;
}
section {
  margin-top: 1em;
}
h1,h2,h3,h4,h5 {
  margin: 0.5em 0;
}
.summary {
  display: none;
}
";

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
  <link href="https://cdn.jsdelivr.net/npm/fastbootstrap@1.1.2/dist/css/fastbootstrap.min.css" rel="stylesheet" integrity="sha256-xLGBU65wCDv2/qEdq3ZYw2Qdiia/wxxeGepRyZmpQdY=" crossorigin="anonymous">
  <style>
{}
  </style>
"##,
            title, title, config.database_name, title, config.database_name,
            config.database_name, style
    )
}

fn header(config: &Config) -> String {
    format!("
 <div>
   <a href='/'><img src='/assets/{}' title='{}' alt='{} home'></a>
 </div>
 ", config.logo_file_name, config.database_long_name, config.site_name)
}

fn gene_summary(config: &Config, gene_details: &GeneDetails) -> String {
    let mut summ = String::new();

    summ += "<dl>\n";
    if let Some(ref name) = gene_details.name {
        summ += &format!("  <dt>Standard name</dt> <dd>{}</dd>\n",  name);
    }
    summ += &format!("  <dt>Systematic ID</dt> <dd>{}</dd>\n", gene_details.uniquename);

    if let Some(ref product) = gene_details.product {
        summ += &format!("  <dt>Product</dt> <dd>{}</dd>\n", product);
    }

    if let Some(gene_organism) = config.organism_by_taxonid(gene_details.taxonid) {
        summ += &format!("  <dt>Organism</dt> <dd>{}\n", gene_organism.scientific_name());

        if !gene_organism.alternative_names.is_empty() {
            summ += &format!(" ({})", gene_organism.alternative_names.iter()
                         .map(FlexStr::to_string).collect::<Vec<_>>().join(", "));
        }
        summ += "</dd>\n";
    }

    if !gene_details.synonyms.is_empty() {
        let synonyms: Vec<FlexStr> =
            gene_details.synonyms.iter().map(|s| s.name.clone()).collect();
        summ += &format!("  <dt>Synonyms</dt> <dd>{}</dd>\n",
                         synonyms.iter().map(FlexStr::to_string).collect::<Vec<_>>().join(", "));
    }

    if let Some(ref uniprot_identifier) = gene_details.uniprot_identifier {
        summ += &format!("  <dt>UniProt ID</dt> <dd>{}</dd>\n", uniprot_identifier);
    }

    if let Some(ref orfeome_identifier) = gene_details.orfeome_identifier {
        summ += &format!("  <dt>ORFeome ID</dt> <dd>{}</dd>\n", orfeome_identifier);
    }

    if let Some(ref characterisation_status) = gene_details.characterisation_status {
        summ += &format!("  <dt>Characterisation status</dt> <dd>{}</dd>\n",
                         characterisation_status);
    }

    summ += &format!("  <dt>Feature type</dt> <dd>{}</dd>\n", gene_details.feature_type);

    if let Some(ref location) = gene_details.location {
        let strand_str = match location.strand {
            Strand::Forward => "forward strand",
            Strand::Reverse => "reverse strand",
            Strand::Unstranded => "",
        };
        let chr_config = config.find_chromosome_config(&location.chromosome_name);

        summ += &format!("  <dt>Genomic location</dt> <dd>chromosome {}: {}..{} {}</dd>\n",
                         chr_config.short_display_name, location.start_pos,
                         location.end_pos, strand_str);
    }

    summ += "</dl>\n";

    summ
}

fn genotype_alleles_summary(genotype_details: &GenotypeDetails) -> String {
    let mut summ = String::new();

    let empty_string = flex_str!("");

    let locus_type_name =
        if genotype_details.ploidiness == Ploidiness::Haploid {
            "gene"
        } else {
            "locus"
        };

    summ += &format!("<table><thead><th>{locus_type_name}</th> <th>Product</th> <th>Allele name</th> <th>Description</th> <th>Type</th> <th>Expression</th></thead>");

    for locus in &genotype_details.loci {
        summ += "<tbody>";

        for (idx, expressed_allele) in locus.expressed_alleles.iter().enumerate() {
            let Some(allele) = genotype_details.get_allele(&expressed_allele.allele_uniquename)
            else {
                continue;
            };
            let Some(gene) = genotype_details.get_gene_short(&allele.gene_uniquename)
            else {
                continue;
            };

            summ += "<tr>";

            let gene_product = gene.product.as_ref().unwrap_or_else(|| &empty_string);

            let allele_expression =
                expressed_allele.expression.as_ref().unwrap_or_else(|| &empty_string);
            let allele_name = allele.name.as_ref().unwrap_or_else(|| &empty_string);
            let allele_description = allele.description.as_ref().unwrap_or_else(|| &empty_string);

            if idx == 0 {
                summ += &format!("<td><a href='/gene/{}'>{}</a></td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td>\n",
                                 gene.uniquename, gene.display_name(), gene_product, allele_name,
                                 allele_description,
                                 allele.allele_type, allele_expression);
            } else {
                summ += &format!("<td></td><td></td><td>{}</td><td>{}</td><td>{}</td><td>{}</td>\n",
                                 allele_name, allele_description, allele.allele_type,
                                 allele_expression);
            }

            summ += "</tr>";
        }

        summ += "</tbody>";
    }

    summ += "</table>";

    summ
}

fn genotype_summary(config: &Config, genotype_details: &GenotypeDetails) -> String {
    let mut summ = String::new();

    summ += "<dl>\n";

    summ += &format!("  <dt>Description</dt> <dd>{}</dd>\n",  genotype_details.display_name);

    if let Some(genotype_organism) = config.organism_by_taxonid(genotype_details.taxonid) {
        summ += &format!("  <dt>Organism</dt> <dd>{}\n", genotype_organism.scientific_name());

        if !genotype_organism.alternative_names.is_empty() {
            summ += &format!(" ({})", genotype_organism.alternative_names.iter()
                         .map(FlexStr::to_string).collect::<Vec<_>>().join(", "));
        }
        summ += "</dd>\n";
    }

    summ += "</dl>\n";

    summ
}

fn get_annotations(ont_annotation_ids: &[OntAnnotationId],
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
                if let Some(Some(gene_short)) =
                    container.genes_by_uniquename().get(gene_uniquename) {
                        genes.push(gene_short)
                    }
            }
        }
    }

    if container.container_type() != ContainerType::Reference && !references.is_empty() {
        references.sort();
        references.dedup();
        ret += "<p>References:</p>\n<ul>";
        for reference in references {
            ret += &format!("<li><a href='/reference/{}'>{}</a></li>\n", reference, reference);
        }
        ret += "</ul>\n";
    }

    if container.container_type() != ContainerType::Gene && !genes.is_empty() {
        genes.sort();
        genes.dedup();
        ret += "<p>Genes:</p>\n<ul>";
        for gene in genes {
            ret += &format!("<li><a href='/gene/{}'>{}</a></li>\n", gene.uniquename,
                            gene.display_name());
        }
        ret += "</ul>\n";
    }

    if container.container_type() != ContainerType::Genotype && !genotypes.is_empty() {
        genotypes.sort();
        genotypes.dedup();
        let genotypes_by_uniquename = container.genotypes_by_uniquename();

        ret += "<p>Genotypes:</p>\n<ul>";
        for uniquename in genotypes {
            if let Some(genotype_short) = genotypes_by_uniquename.get(uniquename) {
                ret += &format!("<li><a href='/genotype/{}'>{}</a></li>\n",
                                uniquename, genotype_short.display_name);
            }
        }
        ret += "</ul>\n";
    }

    ret
}

fn annotation_section(config: &Config, container: &dyn AnnotationContainer) -> String {
    let mut annotation_html = String::new();

    let mut cv_names: Vec<FlexStr> = vec![];

    for cv_name in container.cv_annotations().keys() {
        cv_names.push(cv_name.clone());
    }

    let cmp_display_names =
        |cv_name_a: &FlexStr, cv_name_b: &FlexStr| {
            let cv_config_a = config.cv_config_by_name_with_default(cv_name_a);
            let cv_display_name_a = cv_config_a.display_name;
            let cv_config_b = config.cv_config_by_name_with_default(cv_name_b);
            let cv_display_name_b = cv_config_b.display_name;

            cv_display_name_a.cmp(&cv_display_name_b)
        };

    cv_names.sort_by(cmp_display_names);

    for cv_name in cv_names {
        let term_annotations = container.cv_annotations().get(&cv_name).unwrap();
        let cv_config = config.cv_config_by_name_with_default(&cv_name);

        if let Some(cv_display_name) = cv_config.display_name {
            annotation_html += &format!("<section>\n<h3>{}</h3>\n", cv_display_name);
            for term_annotation in term_annotations {
                let term_name =
                    if let Some(Some(term_short)) =
                    container.terms_by_termid().get(&term_annotation.term) {
                        String::from(term_short.name.as_str())
                    } else {
                        String::from("(unknown term name)")
                    };
                annotation_html += &format!("<section><h4><a href='/term/{}'>{}</a> - \
                                             <a href='/term/{}'>{}</a></h4>\n{}</section>\n",
                                            term_annotation.term,
                                            term_annotation.term,
                                            term_annotation.term,
                                            term_name,
                                            get_annotations(&term_annotation.annotations,
                                                            container));
            }
            annotation_html += "</section>\n";
        }
    }

    annotation_html
}

fn protein_features(gene_details: &GeneDetails) -> String {
    let mut protein_feature_html = String::new();

    protein_feature_html += "<table>\n<thead>\n<tr>\n";
    protein_feature_html += "<th>ID</th><th>Name</th><th>InterPro name</th><th>DB name</th>\n";
    protein_feature_html += "</tr>\n</thead>\n";
    protein_feature_html += "<tbody>\n";

    let empty_str = &"".into();

    for interpro_match in &gene_details.interpro_matches {
        protein_feature_html +=
            &format!("<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>\n",
                     interpro_match.id,
                     interpro_match.name.as_ref().unwrap_or(empty_str),
                     interpro_match.interpro_name.as_ref().unwrap_or(empty_str),
                     interpro_match.dbname);
    }

    protein_feature_html += "</tbody>\n</table>\n";

    protein_feature_html
}

fn orthologs(config: &Config, container: &dyn OrthologAnnotationContainer) -> String {
    let mut orth_html = String::new();

    if container.ortholog_annotations().is_empty() {
        return orth_html;
    }

    orth_html += "<section><h2>Orthologs</h2>\n<ul>\n";

    for orth_annotation in container.ortholog_annotations().iter() {
        let maybe_orth_org = config.organism_by_taxonid(orth_annotation.ortholog_taxonid);
        let orth_org_scientific_name =
            if let Some(orth_org) = maybe_orth_org {
                orth_org.scientific_name()
            } else {
                String::from("Unknown organism")
            };
        if let Some(Some(orth_gene_short)) =
            container.genes_by_uniquename().get(&orth_annotation.ortholog_uniquename) {
                orth_html += &format!("<li>{} {}</li>\n", orth_gene_short.display_name(),
                                      orth_org_scientific_name);
            }
    }

    orth_html += "</ul></section>\n";

    orth_html
}

fn gene_references(gene_details: &GeneDetails) -> String {
    let mut refs_html = String::new();

    refs_html += "<dl>\n";

    let unknown_title = flex_str!("Unknown title");
    let empty = flex_str!("");

    for ref_uniquename in gene_details.references_by_uniquename.keys() {
        let ref_by_uniquename = &gene_details.references_by_uniquename;
        if let Some(ref_short_opt) = ref_by_uniquename.get(ref_uniquename) {
            if let Some(ref_short) = ref_short_opt {
                refs_html += &format!("<dt><a href='/reference/{}'>{}</a> - <a href='/reference/{}'>{}</a></dt> <dd>{} {}</dd>\n",
                                      ref_short.uniquename, ref_short.uniquename,
                                      ref_short.uniquename,
                                      ref_short.title.as_ref().unwrap_or_else(|| &unknown_title),
                                      ref_short.authors_abbrev.as_ref().unwrap_or_else(|| &empty),
                                      ref_short.citation.as_ref().unwrap_or_else(|| &empty));
            }
        }
    }

    refs_html += "</dl>\n";

    refs_html
}

fn gene_body(config: &Config, title: &str, gene_details: &GeneDetails) -> String {
    let mut body = String::new();

    body += &header(config);

    body += &format!("<h1>{}</h1>\n", title);

    body += &format!("<section><h2 class='summary'>Gene summary</h2>\n{}</section>\n",
                     gene_summary(config, gene_details));

    body += &format!("<section><h2>Annotation</h2>\n{}</section>\n",
                     annotation_section(config, gene_details));

    if gene_details.interpro_matches.len() > 0 {
        body += &format!("<section><h2>Protein features</h2>\n{}</section>\n",
                         protein_features(gene_details));
    }

    body += &orthologs(config, gene_details);

    body += &format!("<section><h2>References / Literature</h2>\n{}</section>\n",
                     gene_references(gene_details));

    body
}

fn genotype_body(config: &Config, title: &str, genotype_details: &GenotypeDetails) -> String {
    let mut body = String::new();

    body += &header(config);

    body += &format!("<h1>{}</h1>\n", title);

    body += &format!("<section><h2 class='summary'>Genotype summary</h2>\n{}</section>\n",
                     genotype_summary(config, genotype_details));

    body += &format!("<section><h2>Allele summary</h2>\n{}</section>\n",
                     genotype_alleles_summary(genotype_details));

    body += &format!("<section><h2>Annotation</h2>\n{}</section>\n",
                     annotation_section(config, genotype_details));

    body
}

pub fn render_simple_gene_page(config: &Config, gene_details: &GeneDetails) -> String  {
    let title = make_gene_title(gene_details);

    format_page(&head(config, &title), &gene_body(config, &title, gene_details))
}

pub fn render_simple_genotype_page(config: &Config, genotype_details: &GenotypeDetails) -> String  {
    let title = make_genotype_title(genotype_details);

    format_page(&head(config, &title), &genotype_body(config, &title, genotype_details))
}

fn make_reference_title(reference_details: &ReferenceDetails) -> String {
    if let Some(ref title) = reference_details.title {
        format!("Reference - {} - {}", reference_details.uniquename,
                title)
    } else {
        format!("Reference - {}", reference_details.uniquename)
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

    if let Some(ref pubmed_abstract) = reference_details.pubmed_abstract {
        summ += &format!("<dt>Abstract</dt> <dd>{}</dd>\n", pubmed_abstract);
    }

    summ += "</dl>\n";

    summ
}

fn reference_body(config: &Config, title: &str, reference_details: &ReferenceDetails) -> String {
    let mut body = String::new();

    body += &header(config);

    body += &format!("<h1>{}</h1>\n", title);

    body += &format!("<section><h2 class='summary'>Reference summary</h2>\n{}</section>\n",
                     reference_summary(reference_details));

    body += &format!("<section><h2>Annotation</h2>\n{}</section>\n",
                     annotation_section(config, reference_details));

    body += &orthologs(config, reference_details);

    body
}

pub fn render_simple_reference_page(config: &Config, reference_details: &ReferenceDetails) -> String  {
    let title = make_reference_title(reference_details);

    format_page(&head(config, &title), &reference_body(config, &title, reference_details))
}

fn make_term_title(config: &Config, term_details: &TermDetails) -> String {
    let cv_config = config.cv_config_by_name_with_default(&term_details.cv_name);
    let cv_display_name = cv_config.display_name.unwrap_or_else(|| flex_str!("DEFAULT"));

    format!("{} ontology term - {} - {}", cv_display_name, term_details.termid,
            term_details.name)
}

fn term_summary(config: &Config, term_details: &TermDetails) -> String {
    let mut summ = String::new();

    summ += "<dl>\n";

    summ += &format!("<dt>ID</dt> <dd>{}</dd>\n", term_details.termid);
    summ += &format!("<dt>Name</dt> <dd>{}</dd>\n", term_details.name);

    let cv_config = config.cv_config_by_name_with_default(&term_details.cv_name);
    let cv_display_name = cv_config.display_name.unwrap_or_else(|| flex_str!("DEFAULT"));

    summ += &format!("<dt>Ontology or CV name</dt> <dd>{}</dd>\n", cv_display_name);

    if let Some(ref definition) = term_details.definition {
        summ += &format!("<dt>Definition</dt> <dd>{}</dd>\n", definition);
    }

    summ += "</dl>\n";

    summ
}

fn term_parents(term_details: &TermDetails) -> String {
    let mut parent_string: String = "<ul>".into();

    for parent in &term_details.direct_ancestors {
        parent_string += &format!("<li><b>{}</b> <a href='/term/{}'>{}</a></li>",
                                  parent.relation_name, parent.termid,
                                  parent.term_name);
    }

    parent_string += "</ul>";

    parent_string
}

fn term_body(config: &Config, title: &str, term_details: &TermDetails) -> String {
    let mut body = String::new();

    body += &header(config);

    body += &format!("<h1>{}</h1>\n", title);

    body += &format!("<section><h2 class='summary'>Term summary</h2>\n{}</section>\n",
                     term_summary(config, term_details));

    body += &format!("<section><h2>Parents</h2>\n{}</section>\n",
                     term_parents(term_details));

    body += &format!("<section><h2>Annotation</h2>\n{}</section>\n",
                     annotation_section(config, term_details));

    body
}

pub fn render_simple_term_page(config: &Config, term_details: &TermDetails) -> String  {
    let title = make_term_title(config, term_details);

    format_page(&head(config, &title), &term_body(config, &title, term_details))
}

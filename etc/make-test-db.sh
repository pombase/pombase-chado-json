#!/bin/sh -

# create a small Chado database from a PomBase Chado DB

slim_file_name=/tmp/slim_ids_$$.txt

jq -r '.go_slim_terms | .[] | .name' /var/pomcur/sources/pombe-embl/website/pombase_v2_config.json > $slim_file_name

ids="'SPAC19G12.04','SPAC27E2.05','SPBC1D7.02c','SPBC11G11.01','SPBC18H10.02','SPBC18H10.20c',\
'SPBC9B6.04c','SPAC2F3.09','SPAC11G7.02','SPAC27E2.06c','SPAC19G12.03', 'SPACUNK4.17-antisense-1',\
'SPAC2E1P3.03c','SPRRNA.26','snR94','SPATRNALYS.01','SPAC2E12.05'"

delete_leaf_children="
DELETE FROM cvterm
  WHERE is_relationshiptype = 0
    AND cvterm_id NOT IN (SELECT cvterm_id FROM feature_cvterm)
    AND cvterm_id NOT IN (SELECT object_id FROM cvterm_relationship)
    AND cvterm_id NOT IN (SELECT type_id FROM feature)
    AND name NOT IN (SELECT name FROM go_slim_term_names)
    AND cv_id IN (SELECT cv_id FROM cv
          WHERE name in ('PomBase gene products', 'biological_process', 'cellular_component',
                         'molecular_function', 'fission_yeast_phenotype', 'sequence',
                         'PSI-MOD', 'protein'));
"


psql -e $1 -c "
CREATE TEMP TABLE go_slim_term_names(name text);
COPY go_slim_term_names FROM '$slim_file_name';

DELETE FROM feature
 WHERE type_id IN (SELECT cvterm_id FROM cvterm WHERE name in ('gene','pseudogene'))
   AND uniquename NOT IN ($ids)
   AND feature_id NOT IN
     (SELECT object_id FROM feature_relationship
       WHERE subject_id IN (SELECT feature_id FROM feature WHERE uniquename IN ($ids)) limit 2)
   AND feature_id NOT IN
     (SELECT subject_id FROM feature_relationship
       WHERE object_id IN (SELECT feature_id FROM feature WHERE uniquename IN ($ids)) limit 2);

DELETE FROM feature_cvterm WHERE
   feature_cvterm_id IN
     (SELECT feature_cvterm_id FROM feature_cvtermprop prop
        JOIN cvterm type ON prop.type_id = type.cvterm_id
       WHERE value NOT IN ($ids) AND type.name = 'with');

DELETE FROM cvterm WHERE cvterm_id IN
   (SELECT cvterm_id FROM cvtermprop
     WHERE type_id IN (SELECT cvterm_id FROM cvterm WHERE name LIKE 'annotation_extension_relation-%')
       AND value NOT IN ($ids) AND NOT value ~ '^[0-9]');

DELETE FROM cvterm WHERE cv_id IN (SELECT cv_id FROM cv WHERE name = 'quality');

DELETE FROM feature
 WHERE type_id IN (SELECT cvterm_id FROM cvterm WHERE name = 'allele')
   AND feature_id NOT IN
     (SELECT subject_id FROM feature_relationship WHERE object_id IN
         (SELECT feature_id FROM feature WHERE type_id IN
              (SELECT cvterm_id FROM cvterm WHERE name in ('gene','pseudogene'))));

DELETE FROM feature
 WHERE type_id IN (SELECT cvterm_id FROM cvterm WHERE name = 'genotype')
   AND feature_id NOT IN
    (SELECT object_id FROM feature_relationship
      WHERE subject_id IN
         (SELECT feature_id FROM feature WHERE type_id IN
              (SELECT cvterm_id FROM cvterm WHERE name = 'allele')));

DELETE FROM feature
 WHERE type_id IN (SELECT cvterm_id FROM cvterm WHERE name like '%RNA' OR name = 'pseudogenic_transcript')
    AND feature_id NOT IN
      (SELECT subject_id FROM feature_relationship
        WHERE type_id IN
            (SELECT cvterm_id FROM cvterm WHERE name = 'part_of')
          AND object_id IN
            (SELECT feature_id FROM feature
              WHERE type_id IN (SELECT cvterm_id FROM cvterm WHERE name in ('gene','pseudogene'))));

DELETE FROM feature
 WHERE type_id IN
    (SELECT cvterm_id FROM cvterm
      WHERE name in ('exon', 'intron', 'three_prime_UTR', 'five_prime_UTR'))
    AND feature_id NOT IN
      (SELECT subject_id FROM feature_relationship
        WHERE type_id IN
            (SELECT cvterm_id FROM cvterm WHERE name = 'part_of')
          AND object_id IN
            (SELECT feature_id FROM feature
              WHERE type_id IN (SELECT cvterm_id FROM cvterm
                                 WHERE name like '%RNA' OR name = 'pseudogenic_transcript')));

DELETE FROM feature
 WHERE type_id IN
    (SELECT cvterm_id FROM cvterm WHERE name = 'polypeptide')
    AND feature_id NOT IN
      (SELECT subject_id FROM feature_relationship
        WHERE type_id IN
            (SELECT cvterm_id FROM cvterm WHERE name = 'derives_from')
          AND object_id IN
            (SELECT feature_id FROM feature
              WHERE type_id IN (SELECT cvterm_id FROM cvterm
                                 WHERE name = 'mRNA')));

DELETE FROM cvterm
   WHERE cv_id = (SELECT cv_id FROM cv WHERE name = 'PomBase annotation extension terms')
     AND cvterm_id NOT IN (SELECT cvterm_id FROM feature_cvterm);

$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children
$delete_leaf_children

DELETE FROM cvterm
 WHERE cv_id in (SELECT cv_id FROM cv
                  WHERE name = 'cvterm_property_type')
                    AND cvterm_id NOT IN (SELECT type_id FROM cvtermprop);

DELETE FROM pub
   WHERE pub_id NOT IN (SELECT pub_id FROM feature_cvterm)
     AND pub_id NOT IN (SELECT pub_id FROM feature_pub)
     AND pub_id NOT IN (SELECT pub_id FROM feature_synonym)
     AND pub_id NOT IN (SELECT pub_id FROM feature_relationship_pub);
" < /dev/null;

rm $slim_file_name

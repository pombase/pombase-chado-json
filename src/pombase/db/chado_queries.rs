extern crate tokio_postgres;

use std::collections::{BTreeSet, HashMap};

use crate::data_types::StatsIntegerTable;

use self::tokio_postgres::Client;

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct CommunityResponseRate {
  pub year: i32,
  pub submitted: i64,
  pub sent_sessions: i64,
  pub response_rate: f32,
}

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct ChadoQueries {
  pub community_response_rates: Vec<CommunityResponseRate>,
  pub annotation_type_counts_by_year: StatsIntegerTable,
}

const RESPONSE_RATE_SQL: &str = r#"
WITH
  session_years AS
     (SELECT canto_session, canto_curator_role,
          coalesce(canto_first_sent_to_curator_year, least(canto_session_accepted_year, canto_session_submitted_year, canto_first_approved_year)) AS sent_or_accepted_year,
          coalesce(canto_session_submitted_year, canto_first_approved_year) AS submitted_year
      FROM pombase_publication_curation_summary
      WHERE canto_curator_role = 'community'),
  year_series AS (SELECT * FROM generate_series(2013, (SELECT extract(YEAR FROM CURRENT_DATE))::integer) YEAR),
  cumulative_years AS
    (SELECT year_series.year,
        (SELECT count(*) FROM session_years WHERE session_years.submitted_year IS NOT NULL
     AND session_years.submitted_year <= year_series.year) as submitted_count,
        (SELECT count(*) FROM session_years WHERE session_years.sent_or_accepted_year IS NOT NULL
      AND session_years.sent_or_accepted_year <= year_series.year) as sent_or_accepted_count
   FROM year_series)
SELECT year_series.year,
       submitted_count, sent_or_accepted_count,
       CASE WHEN sent_or_accepted_count > 0 THEN trunc(100.0* submitted_count / sent_or_accepted_count, 1)::real ELSE 0 END AS response_rate
  FROM year_series join cumulative_years on cumulative_years.year = year_series.year;
"#;

async fn get_community_response_rates(conn: &mut Client)
  -> Result<Vec<CommunityResponseRate>, tokio_postgres::Error>
{
  let mut ret = vec![];

  let result = conn.query(RESPONSE_RATE_SQL, &[]).await?;

  for row in &result {
    let el = CommunityResponseRate {
      year: row.get(0),
      submitted: row.get(1),
      sent_sessions: row.get(2),
      response_rate: row.get(3),
    };

    ret.push(el);
  }

  Ok(ret)
}

const ANNOTATION_TYPE_COUNT_SQL: &str = r#"
SELECT annotation_year, annotation_type, count(DISTINCT id)
  FROM pombase_genes_annotations_dates
 WHERE annotation_year IS NOT NULL
   AND (evidence_code IS NULL OR evidence_code <> 'Inferred from Electronic Annotation')
  AND (annotation_type NOT IN ('cat_act', 'subunit_composition', 'external_link', 'pathway'))
  AND (annotation_source IS NULL OR annotation_source <> 'BIOGRID')
 GROUP BY annotation_year, annotation_type
 ORDER BY annotation_year, annotation_type
"#;

async fn get_annotation_type_counts(conn: &mut Client)
    -> Result<StatsIntegerTable, tokio_postgres::Error>
{
  let mut res = HashMap::new();

  let mut first_year: i32 = 9999;
  let mut last_year: i32 = 0;
  let mut annotation_types = BTreeSet::new();

  let result = conn.query(ANNOTATION_TYPE_COUNT_SQL, &[]).await?;

  for row in &result {

    let year: i32 = row.get(0);
    let annotation_type: String = row.get(1);
    let count: i64 = row.get(2);

    if year > last_year {
      last_year = year;
    }
    if year < first_year {
      first_year = year;
    }

    annotation_types.insert(annotation_type.clone());

    res.entry(year)
       .or_insert_with(HashMap::new)
       .insert(annotation_type, count);
  }

  let header = annotation_types.iter().cloned().collect::<Vec<_>>();
  let mut data = vec![];

  for year in first_year..=last_year {
    let mut data_row = vec![];

    if let Some(year_data) = res.get(&year) {
      for annotation_type in annotation_types.iter() {
        if let Some(count) = year_data.get(annotation_type) {
          data_row.push(*count as usize);
        } else {
          data_row.push(0usize);
        }
      }
    } else {
      data_row = vec![0; annotation_types.len()]
    }

    data.push((year.to_string(), data_row));
  }

  Ok(StatsIntegerTable {
    header,
    data,
  })
}

impl ChadoQueries {
  pub async fn new(conn: &mut Client) -> Result<ChadoQueries, tokio_postgres::Error> {
    let community_response_rates = get_community_response_rates(conn).await?;
    let annotation_type_counts_by_year = get_annotation_type_counts(conn).await?;
    Ok(ChadoQueries {
      community_response_rates,
      annotation_type_counts_by_year,
    })
  }
}

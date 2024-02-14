extern crate tokio_postgres;

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

impl ChadoQueries {
  pub async fn new(conn: &mut Client) -> Result<ChadoQueries, tokio_postgres::Error> {
    let community_response_rates = get_community_response_rates(conn).await?;
    Ok(ChadoQueries {
      community_response_rates,
    })
  }
}

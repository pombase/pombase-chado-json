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

const RESPONSE_RATE_SQL: &str = "
WITH counts AS
  (SELECT YEAR,
     (SELECT COUNT (*)
      FROM pombase_publication_curation_summary
      WHERE canto_curator_role = 'community'
        AND (canto_annotation_status = 'NEEDS_APPROVAL'
             OR canto_annotation_status = 'APPROVAL_IN_PROGRESS'
             OR canto_annotation_status = 'APPROVED')
        AND (canto_session_submitted_date IS NOT NULL
             AND canto_session_submitted_date <= (YEAR || '-12-30')::date)) AS submitted,
     (SELECT COUNT (*)
      FROM pombase_publication_curation_summary
      WHERE canto_curator_role = 'community'
        AND ((canto_first_sent_to_curator_year IS NOT NULL
              AND canto_first_sent_to_curator_year <= YEAR)
             OR (canto_session_accepted_year IS NOT NULL
                 AND canto_session_accepted_year <= YEAR))) AS sent_sessions
   FROM generate_series(2013, (SELECT extract(YEAR FROM CURRENT_DATE))::integer) AS YEAR)
SELECT YEAR, submitted, sent_sessions, trunc(100.0*submitted/sent_sessions, 1)::real AS response_rate
FROM counts;
";

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

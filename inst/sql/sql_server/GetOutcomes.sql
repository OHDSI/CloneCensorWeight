-- Outcome events for subjects in the temp cohort table (#ccw_cohorts).
--
-- Parameters:
--   @outcome_database_schema   Schema containing the outcome cohort table
--   @outcome_table             Outcome cohort table name
--   @outcome_ids               Comma-separated outcome cohort definition IDs

SELECT
  o.subject_id,
  o.cohort_start_date AS outcome_date,
  o.cohort_definition_id AS outcome_id
FROM @outcome_database_schema.@outcome_table o
INNER JOIN #ccw_cohorts c
  ON o.subject_id = c.subject_id
WHERE o.cohort_definition_id IN (@outcome_ids)
;

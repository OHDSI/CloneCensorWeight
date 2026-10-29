SELECT
  outcome.subject_id,
  outcome.cohort_start_date AS outcome_date,
  outcome.cohort_definition_id AS outcome_id
FROM @outcome_database_schema.@outcome_table AS outcome
WHERE outcome.cohort_definition_id IN (@outcome_ids);
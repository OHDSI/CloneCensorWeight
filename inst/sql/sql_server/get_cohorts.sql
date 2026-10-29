SELECT
  target.subject_id,
  target.cohort_start_date,
  target.cohort_end_date,
  exposure.exposure_start_date AS initiation_date
FROM @target_database_schema.@target_table AS target
LEFT JOIN @exposure_database_schema.@exposure_table AS exposure
  ON target.subject_id = exposure.subject_id
  AND exposure.cohort_definition_id = @exposure_id
WHERE target.cohort_definition_id = @target_id;
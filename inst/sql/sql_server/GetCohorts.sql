-- Eligible target population with treatment initiation dates.
--
-- Parameters:
--   @cdm_database_schema        CDM schema
--   @target_database_schema     Schema containing the target cohort table
--   @target_table               Target cohort table name
--   @target_id                  Target cohort definition ID
--   @exposure_database_schema   Schema containing the exposure cohort table
--   @exposure_table             Exposure cohort table name
--   @exposure_id                Exposure cohort definition ID
--   @washout_period             Required days of prior observation

SELECT
  t.subject_id,
  t.cohort_start_date,
  t.cohort_end_date,
  e.cohort_start_date AS initiation_date
FROM @target_database_schema.@target_table t
INNER JOIN @cdm_database_schema.observation_period op
  ON  t.subject_id = op.person_id
  AND t.cohort_start_date >= DATEADD(DAY, @washout_period, op.observation_period_start_date)
  AND t.cohort_start_date <= op.observation_period_end_date
LEFT JOIN @exposure_database_schema.@exposure_table e
  ON  t.subject_id            = e.subject_id
  AND e.cohort_definition_id  = @exposure_id
  AND e.cohort_start_date    >= t.cohort_start_date
  AND e.cohort_start_date    <= t.cohort_end_date
WHERE t.cohort_definition_id = @target_id
;

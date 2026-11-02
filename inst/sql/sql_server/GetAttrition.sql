-- Return per-step subject counts for the attrition table.
-- Additional rows for prior-outcome exclusion and min-days-at-risk filter
-- are added in Phase 5 when createStudyPopulation() writes #ccw_excluded.
--
-- Requires: #ccw_cohorts, #ccw_clones

SELECT
  step_order,
  description,
  n_subjects,
  n_excluded
FROM (

  SELECT
    1 AS step_order,
    'Target cohort after washout period' AS description,
    COUNT(DISTINCT subject_id) AS n_subjects,
    0 AS n_excluded
  FROM #ccw_cohorts

  UNION ALL

  SELECT
    2 AS step_order,
    'Clones created (both arms)' AS description,
    COUNT(*) AS n_subjects,
    0 AS n_excluded
  FROM #ccw_clones

  UNION ALL

  SELECT
    3 AS step_order,
    'Artificially censored clones' AS description,
    SUM(artificially_censored) AS n_subjects,
    SUM(artificially_censored) AS n_excluded
  FROM #ccw_clones

  UNION ALL

  SELECT
    4 AS step_order,
    'Clones retained after artificial censoring' AS description,
    SUM(1 - artificially_censored) AS n_subjects,
    0 AS n_excluded
  FROM #ccw_clones

) attrition
ORDER BY step_order
;

-- Duplicate every subject across two treatment-strategy arms and apply
-- artificial censoring. The result is written to #ccw_clones.
--
-- For 'initiateWithin': clone is artificially censored at day @grace_period if
-- no initiation has occurred by then.
-- For 'deferBeyond': clone is artificially censored at the day of initiation
-- if initiation occurs within day @grace_period.
--
-- Parameters:
-- @arm1_id Arm 1 identifier string
-- @arm1_strategy 'initiateWithin' or 'deferBeyond'
-- @arm1_grace_period Grace period in days for arm 1
-- @arm2_id Arm 2 identifier string
-- @arm2_strategy 'initiateWithin' or 'deferBeyond'
-- @arm2_grace_period Grace period in days for arm 2
--
-- Requires: #ccw_cohorts (subject_id, cohort_start_date, cohort_end_date,
-- initiation_date)

SELECT
  subject_id,
  arm_id,
  cohort_start_date,
  cohort_end_date,
  artificially_censored
INTO #ccw_clones
FROM (

  -- ── Arm 1 clones ──────────────────────────────────────────────────────────
  SELECT
    c.subject_id,
    CAST('@arm1_id' AS VARCHAR(255)) AS arm_id,
    c.cohort_start_date,

{@arm1_strategy == 'initiateWithin'} ? {
    CASE
      WHEN c.initiation_date IS NULL
        OR DATEDIFF(DAY, c.cohort_start_date, c.initiation_date) > @arm1_grace_period
      THEN
        CASE
          WHEN DATEADD(DAY, @arm1_grace_period, c.cohort_start_date) < c.cohort_end_date
            THEN DATEADD(DAY, @arm1_grace_period, c.cohort_start_date)
          ELSE c.cohort_end_date
        END
      ELSE c.cohort_end_date
    END AS cohort_end_date,
    CASE
      WHEN c.initiation_date IS NULL
        OR DATEDIFF(DAY, c.cohort_start_date, c.initiation_date) > @arm1_grace_period
      THEN 1
      ELSE 0
    END AS artificially_censored
} : {
    CASE
      WHEN c.initiation_date IS NOT NULL
        AND DATEDIFF(DAY, c.cohort_start_date, c.initiation_date)
              BETWEEN 0 AND @arm1_grace_period
      THEN c.initiation_date
      ELSE c.cohort_end_date
    END AS cohort_end_date,
    CASE
      WHEN c.initiation_date IS NOT NULL
        AND DATEDIFF(DAY, c.cohort_start_date, c.initiation_date)
              BETWEEN 0 AND @arm1_grace_period
      THEN 1
      ELSE 0
    END AS artificially_censored
}

  FROM #ccw_cohorts c

  UNION ALL

  -- ── Arm 2 clones ──────────────────────────────────────────────────────────
  SELECT
    c.subject_id,
    CAST('@arm2_id' AS VARCHAR(255)) AS arm_id,
    c.cohort_start_date,

{@arm2_strategy == 'initiateWithin'} ? {
    CASE
      WHEN c.initiation_date IS NULL
        OR DATEDIFF(DAY, c.cohort_start_date, c.initiation_date) > @arm2_grace_period
      THEN
        CASE
          WHEN DATEADD(DAY, @arm2_grace_period, c.cohort_start_date) < c.cohort_end_date
            THEN DATEADD(DAY, @arm2_grace_period, c.cohort_start_date)
          ELSE c.cohort_end_date
        END
      ELSE c.cohort_end_date
    END AS cohort_end_date,
    CASE
      WHEN c.initiation_date IS NULL
        OR DATEDIFF(DAY, c.cohort_start_date, c.initiation_date) > @arm2_grace_period
      THEN 1
      ELSE 0
    END AS artificially_censored
} : {
    CASE
      WHEN c.initiation_date IS NOT NULL
        AND DATEDIFF(DAY, c.cohort_start_date, c.initiation_date)
              BETWEEN 0 AND @arm2_grace_period
      THEN c.initiation_date
      ELSE c.cohort_end_date
    END AS cohort_end_date,
    CASE
      WHEN c.initiation_date IS NOT NULL
        AND DATEDIFF(DAY, c.cohort_start_date, c.initiation_date)
              BETWEEN 0 AND @arm2_grace_period
      THEN 1
      ELSE 0
    END AS artificially_censored
}

  FROM #ccw_cohorts c

) all_clones
;

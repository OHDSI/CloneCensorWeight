-- Expand clone follow-up into fixed-length time intervals and join outcomes.
-- The result is written to #ccw_person_time.
--
-- Intervals are numbered 0, 1, 2, ... from cohort_start_date.
-- The final interval is truncated at cohort_end_date.
-- event = 1 when the first qualifying outcome falls within the interval.
-- artificially_censored is propagated only to the final interval of each
-- clone; earlier intervals carry 0 (the clone was still under observation).
--
-- Parameters:
--   @interval_length    Length of each interval in days (> 0)
--   @outcome_id         Outcome cohort_definition_id to flag as events
--
-- Requires:
--   #ccw_clones         (subject_id, arm_id, cohort_start_date,
--                         cohort_end_date, artificially_censored)
--   #ccw_outcomes       (subject_id, outcome_date, outcome_id)
--   #ccw_intervals      (n INT) — integer sequence 0, 1, ..., max_intervals
--                         uploaded by R before this SQL runs

WITH first_outcomes AS (

  SELECT subject_id,
         MIN(outcome_date) AS first_outcome_date
  FROM   #ccw_outcomes
  WHERE  outcome_id = @outcome_id
  GROUP BY subject_id

)

SELECT
  c.subject_id,
  c.arm_id,
  i.n                                                      AS interval_num,
  DATEADD(DAY,  i.n       * @interval_length,
          c.cohort_start_date)                             AS interval_start_date,
  CASE
    WHEN DATEADD(DAY, (i.n + 1) * @interval_length,
                 c.cohort_start_date) < c.cohort_end_date
      THEN DATEADD(DAY, (i.n + 1) * @interval_length,
                   c.cohort_start_date)
    ELSE c.cohort_end_date
  END                                                      AS interval_end_date,
  CASE
    WHEN fo.first_outcome_date IS NOT NULL
      AND fo.first_outcome_date >= DATEADD(DAY, i.n * @interval_length,
                                            c.cohort_start_date)
      AND fo.first_outcome_date <  CASE
            WHEN DATEADD(DAY, (i.n + 1) * @interval_length,
                         c.cohort_start_date) < c.cohort_end_date
              THEN DATEADD(DAY, (i.n + 1) * @interval_length,
                           c.cohort_start_date)
            ELSE c.cohort_end_date
          END
    THEN 1
    ELSE 0
  END                                                      AS event,
  -- Carry artificial-censoring flag only to the last interval of this clone
  CASE
    WHEN DATEADD(DAY, (i.n + 1) * @interval_length,
                 c.cohort_start_date) < c.cohort_end_date
      THEN 0
    ELSE c.artificially_censored
  END                                                      AS artificially_censored
INTO #ccw_person_time
FROM   #ccw_clones c
CROSS JOIN #ccw_intervals i
-- Include only intervals that start before the clone's end of follow-up
WHERE  DATEADD(DAY, i.n * @interval_length,
               c.cohort_start_date) < c.cohort_end_date
LEFT JOIN first_outcomes fo
  ON fo.subject_id = c.subject_id
;

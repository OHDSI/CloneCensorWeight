/************************************************************************
@file CreateOrCountCohorts.sql

Copyright 2026 Observational Health Data Sciences and Informatics

This file is part of CloneCensorWeight

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Modelled on CohortMethod's CreateOrCountCohorts.sql. Builds the eligible
cohort (action = 'CREATE', into #ccw_cohort) or counts it at every step
(action = 'COUNT', for the attrition table).

Prior use of a treatment is an era of one of the treatment cohorts that
overlaps days -L to -1 relative to the index date (L = prior use lookback),
so starts on day 0 are allowed. The R code only asks for the removal of
prior users when L > 0.
************************************************************************/

{DEFAULT @cdm_database_schema = 'CDM_SIM'}
{DEFAULT @target_database_schema = 'CDM_SIM'}
{DEFAULT @target_table = 'cohort'}
{DEFAULT @exposure_database_schema = 'CDM_SIM'}
{DEFAULT @exposure_table = 'cohort'}
{DEFAULT @target_id = ''}
{DEFAULT @treatment_ids = ''}
{DEFAULT @study_start_date = ''}
{DEFAULT @study_end_date = ''}
{DEFAULT @first_only = FALSE}
{DEFAULT @washout_period = 0}
{DEFAULT @remove_subjects_with_prior_treatment = FALSE}
{DEFAULT @prior_use_lookback = 365}
{DEFAULT @action = 'CREATE'}

{@action == 'CREATE'} ? {
DROP TABLE IF EXISTS #ccw_cohort;
}

-- raw_cohorts
WITH raw_cohorts AS (
	SELECT target_table.subject_id,
		target_table.cohort_definition_id,
		target_table.cohort_start_date,
		target_table.cohort_end_date
	FROM @target_database_schema.@target_table target_table
	WHERE target_table.cohort_definition_id = @target_id
)

{@first_only} ? {
-- first_only
, first_only AS (
	SELECT subject_id,
		cohort_definition_id,
		cohort_start_date,
		cohort_end_date
	FROM (
		SELECT subject_id,
			cohort_definition_id,
			cohort_start_date,
			cohort_end_date,
			ROW_NUMBER() OVER (PARTITION BY subject_id ORDER BY cohort_start_date) AS cohort_number
		FROM raw_cohorts
	) tmp
	WHERE cohort_number = 1
)
}

-- washout_period: index date in an observation period, with enough prior observation
, washout_period AS (
	SELECT subject_id,
		cohort_definition_id,
		cohort_start_date,
		cohort_end_date,
		observation_period_start_date,
		observation_period_end_date
	FROM {@first_only} ? {first_only} : {raw_cohorts} cohorts
	INNER JOIN @cdm_database_schema.observation_period
		ON cohorts.subject_id = observation_period.person_id
	WHERE cohort_start_date >= observation_period_start_date
		AND cohort_start_date <= observation_period_end_date
		AND DATEDIFF(DAY, observation_period_start_date, cohort_start_date) >= @washout_period
)

{@study_start_date != '' | @study_end_date != ''} ? {
-- restrict_to_study_period
, restrict_to_study_period AS (
	SELECT subject_id,
		cohort_definition_id,
		cohort_start_date,
		cohort_end_date,
		observation_period_start_date,
		observation_period_end_date
	FROM washout_period
{@study_start_date != ''} ? {	WHERE cohort_start_date >= CAST('@study_start_date' AS DATE)}
{@study_end_date != ''} ? {
{@study_start_date != ''} ? {		AND} : {	WHERE} cohort_start_date <= CAST('@study_end_date' AS DATE)
}
)
}

{@remove_subjects_with_prior_treatment} ? {
-- no_prior_treatment
, no_prior_treatment AS (
	SELECT subject_id,
		cohort_definition_id,
		cohort_start_date,
		cohort_end_date,
		observation_period_start_date,
		observation_period_end_date
	FROM {@study_start_date != '' | @study_end_date != ''} ? {restrict_to_study_period} : {washout_period} cohorts
	WHERE NOT EXISTS (
		SELECT 1
		FROM @exposure_database_schema.@exposure_table exposure_table
		WHERE exposure_table.subject_id = cohorts.subject_id
			AND exposure_table.cohort_definition_id IN (@treatment_ids)
			AND exposure_table.cohort_start_date < cohorts.cohort_start_date
			AND exposure_table.cohort_end_date >= DATEADD(DAY, -@prior_use_lookback, cohorts.cohort_start_date)
	)
)
}

{@action == 'CREATE'} ? {
-- Generating a unique row_id that can be used instead of the compound key
-- (subject_id, cohort_start_date) for easier joining etc.
SELECT ROW_NUMBER() OVER (ORDER BY subject_id, cohort_start_date) AS row_id,
	subject_id,
	cohort_definition_id,
	cohort_start_date,
	DATEDIFF(DAY, observation_period_start_date, cohort_start_date) AS days_from_obs_start,
{@study_end_date != ''} ? {
	CASE
		WHEN cohort_end_date <= CAST('@study_end_date' AS DATE)
			THEN DATEDIFF(DAY, cohort_start_date, cohort_end_date)
		ELSE
			DATEDIFF(DAY, cohort_start_date, CAST('@study_end_date' AS DATE))
	END
} : {
	DATEDIFF(DAY, cohort_start_date, cohort_end_date)
} AS days_to_cohort_end,
{@study_end_date != ''} ? {
	CASE
		WHEN observation_period_end_date <= CAST('@study_end_date' AS DATE)
			THEN DATEDIFF(DAY, cohort_start_date, observation_period_end_date)
		ELSE
			DATEDIFF(DAY, cohort_start_date, CAST('@study_end_date' AS DATE))
	END
} : {
	DATEDIFF(DAY, cohort_start_date, observation_period_end_date)
} AS days_to_obs_end,
	cohort_end_date
INTO #ccw_cohort
FROM {@remove_subjects_with_prior_treatment} ? {no_prior_treatment} : {{@study_start_date != '' | @study_end_date != ''} ? {restrict_to_study_period} : {washout_period}} cohorts;
} : {
-- Count the cohort at the various stages
SELECT COUNT(*) AS exposures,
	COUNT(DISTINCT subject_id) AS persons,
	'Original cohorts' AS description,
	CAST(1 AS FLOAT) AS seq_id
FROM raw_cohorts

{@first_only} ? {
UNION ALL

SELECT COUNT(*) AS exposures,
	COUNT(DISTINCT subject_id) AS persons,
	'First cohort entry only' AS description,
	CAST(2 AS FLOAT) AS seq_id
FROM first_only
}

UNION ALL

SELECT COUNT(*) AS exposures,
	COUNT(DISTINCT subject_id) AS persons,
	'@washout_period days of prior observation' AS description,
	CAST(3 AS FLOAT) AS seq_id
FROM washout_period

{@study_start_date != '' | @study_end_date != ''} ? {
UNION ALL

SELECT COUNT(*) AS exposures,
	COUNT(DISTINCT subject_id) AS persons,
	'Restrict to study period' AS description,
	CAST(4 AS FLOAT) AS seq_id
FROM restrict_to_study_period
}

{@remove_subjects_with_prior_treatment} ? {
UNION ALL

SELECT COUNT(*) AS exposures,
	COUNT(DISTINCT subject_id) AS persons,
	'No prior treatment use in the @prior_use_lookback days before index' AS description,
	CAST(5 AS FLOAT) AS seq_id
FROM no_prior_treatment
}
;
}

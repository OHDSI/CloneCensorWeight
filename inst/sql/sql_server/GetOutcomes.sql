/************************************************************************
@file GetOutcomes.sql

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

Modelled on CohortMethod's GetOutcomes.sql (cohort-table branch): every
outcome cohort start in the index observation period, relative to the index
date.
************************************************************************/

{DEFAULT @outcome_database_schema = 'CDM_SIM'}
{DEFAULT @outcome_table = 'cohort'}
{DEFAULT @outcome_ids = ''}
{DEFAULT @cohort_table = '#ccw_cohort'}

SELECT DISTINCT cohort.row_id,
	outcome.cohort_definition_id AS outcome_id,
	DATEDIFF(DAY, cohort.cohort_start_date, outcome.cohort_start_date) AS days_to_event
FROM @cohort_table cohort
INNER JOIN @outcome_database_schema.@outcome_table outcome
	ON cohort.subject_id = outcome.subject_id
WHERE DATEDIFF(DAY, outcome.cohort_start_date, cohort.cohort_start_date) <= cohort.days_from_obs_start
	AND DATEDIFF(DAY, cohort.cohort_start_date, outcome.cohort_start_date) <= cohort.days_to_obs_end
	AND outcome.cohort_definition_id IN (@outcome_ids);

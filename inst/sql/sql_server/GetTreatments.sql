/************************************************************************
@file GetTreatments.sql

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

Every era of the treatment cohorts that overlaps days -L to days_to_obs_end
relative to the index date (L = prior use lookback).
************************************************************************/

{DEFAULT @exposure_database_schema = 'CDM_SIM'}
{DEFAULT @exposure_table = 'cohort'}
{DEFAULT @treatment_ids = ''}
{DEFAULT @prior_use_lookback = 365}
{DEFAULT @cohort_table = '#ccw_cohort'}

SELECT DISTINCT cohort.row_id,
	exposure_table.cohort_definition_id AS treatment_id,
	DATEDIFF(DAY, cohort.cohort_start_date, exposure_table.cohort_start_date) AS days_to_start,
	DATEDIFF(DAY, cohort.cohort_start_date, exposure_table.cohort_end_date) AS days_to_end
FROM @cohort_table cohort
INNER JOIN @exposure_database_schema.@exposure_table exposure_table
	ON cohort.subject_id = exposure_table.subject_id
WHERE exposure_table.cohort_definition_id IN (@treatment_ids)
	AND DATEDIFF(DAY, cohort.cohort_start_date, exposure_table.cohort_start_date) <= cohort.days_to_obs_end
	AND DATEDIFF(DAY, cohort.cohort_start_date, exposure_table.cohort_end_date) >= -@prior_use_lookback;

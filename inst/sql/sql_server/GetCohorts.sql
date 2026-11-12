/************************************************************************
@file GetCohorts.sql

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

Modelled on CohortMethod's GetCohorts.sql. Adds days_to_death, the days
from the index date to the person's earliest death date (NULL if none).
************************************************************************/

{DEFAULT @cdm_database_schema = 'CDM_SIM'}
{DEFAULT @cohort_table = '#ccw_cohort'}

-- Storing person_id as a VARCHAR because R doesn't support 64-bit integers. Generating
-- a proper integer (person_seq_id) on the fly for easy joining, etc.
SELECT cohort.row_id,
	unique_ids.person_seq_id,
	CAST(cohort.subject_id AS VARCHAR(30)) AS person_id,
	CAST(1 AS INT) AS treatment,
	cohort.cohort_start_date,
	cohort.days_from_obs_start,
	cohort.days_to_cohort_end,
	cohort.days_to_obs_end,
	DATEDIFF(DAY, cohort.cohort_start_date, death.death_date) AS days_to_death
FROM @cohort_table cohort
INNER JOIN (
	SELECT subject_id,
		ROW_NUMBER() OVER (ORDER BY subject_id) AS person_seq_id
	FROM (
		SELECT DISTINCT subject_id
		FROM @cohort_table
	) tmp
) unique_ids
	ON cohort.subject_id = unique_ids.subject_id
LEFT JOIN (
	SELECT person_id,
		MIN(death_date) AS death_date
	FROM @cdm_database_schema.death
	GROUP BY person_id
) death
	ON cohort.subject_id = death.person_id
ORDER BY cohort.row_id;

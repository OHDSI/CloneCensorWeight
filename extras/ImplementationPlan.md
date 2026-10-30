# CloneCensorWeight — Implementation Plan

Prepared for agent execution. Each phase is self-contained; implement phases in order.
All code must pass `R CMD check --as-cran` with 0 errors, 0 warnings, and at most 1 NOTE
(the standard "New submission" CRAN note).

---

## Confirmed design decisions

| Decision | Choice |
|---|---|
| Arm types | Generic: `"initiateWithin"`, `"neverInitiate"`, `"sustainedUse"` |
| Arms per analysis | Exactly 2 |
| Classical censoring | Second separate IPCW model; combined weight = product |
| Clustered SEs | Both robust sandwich (Lin-Wei) and bootstrap (subject-level) |
| Baseline confounding | IPCW only — no PS trimming/matching/stratification |
| Empirical calibration | Integrated from v1 via `EmpiricalCalibration` package |
| Person-time expansion | SQL/server-side via `SqlRender` + `DatabaseConnector` |
| Strategus module | Standalone R6 class; no `Strategus` package import |
| Results schema | CohortMethod pattern; `ResultModelManager`-compatible CSVs |
| `sustainedUse` strategy | Deferred to v1.1; not in v1 scope |
| Style | HADES / CohortMethod conventions throughout |
| CRAN target | Yes — v0.1.0 |

---

## Package metadata (`DESCRIPTION`)

Update the following fields:

```
Package: CloneCensorWeight
Version: 0.1.0
Title: Clone-Censor-Weight Target Trial Emulation
Maintainer: William J. O'Brien <wobrien8@its.jnj.com>
Authors@R:
    c(person("William J.", "O'Brien", email = "wobrien8@its.jnj.com",
             role = c("aut", "cre")),
      person("Patrick", "Ryan", email = "ryan@ohdsi.org",
             role = "ctb"),
      person("Observational Health Data Sciences and Informatics",
             role = "cph"))
Description: Implements clone-censor-weight target trial emulation for
    observational studies in the OMOP Common Data Model. The package
    follows the OHDSI HADES framework and provides data extraction,
    cloning, artificial and classical inverse probability of censoring
    weighting, clustered-variance outcome models, covariate balance
    diagnostics, empirical calibration, and a Strategus-compatible
    module interface.
License: Apache License (>= 2)
Encoding: UTF-8
LazyData: false
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.3.2
Depends:
    R (>= 4.1.0),
    Andromeda (>= 0.6.3)
Imports:
    checkmate,
    Cyclops (>= 3.6.0),
    DatabaseConnector (>= 6.0.0),
    dplyr,
    EmpiricalCalibration (>= 3.1.0),
    FeatureExtraction (>= 3.5.0),
    jsonlite,
    methods,
    ParallelLogger (>= 3.4.2),
    R6,
    readr,
    rlang,
    SqlRender (>= 1.18.0),
    survival,
    utils
Suggests:
    boot,
    Eunomia,
    ggplot2,
    knitr,
    ResultModelManager,
    rmarkdown,
    testthat (>= 3.0.0),
    withr
VignetteBuilder: knitr
Config/testthat/edition: 3
URL: https://ohdsi.github.io/CloneCensorWeight
BugReports: https://github.com/ohdsi/CloneCensorWeight/issues
```

---

## File map

### Keep and refactor
| Existing file | Action |
|---|---|
| `R/CloneCensorWeightData.R` | Major refactor — full FeatureExtraction integration |
| `R/Cloning.R` | Refactor — delegate to SQL temp table |
| `R/Weighting.R` | Refactor — dual IPCW, SQL-expanded intervals |
| `R/OutcomeModel.R` | Refactor — robust sandwich + bootstrap SEs |
| `R/Analysis.R` | Replace with `R/RunAnalyses.R` (CohortMethod naming) |
| `R/zzz.R` | Keep, add `utils::globalVariables()` declarations |
| `inst/sql/sql_server/get_cohorts.sql` | Rename to `GetCohorts.sql`, fix column names |
| `inst/sql/sql_server/get_outcomes.sql` | Rename to `GetOutcomes.sql` |

### New files to create
```
R/HelperFunctions.R
R/SettingsObjects.R
R/StudyPopulation.R
R/Balance.R
R/KaplanMeier.R
R/Power.R
R/EmpiricalCalibration.R
R/ResultsDataModel.R
R/Export.R
R/RunAnalyses.R
R/Simulation.R
R/Strategus.R
inst/sql/sql_server/CreateClones.sql
inst/sql/sql_server/DiscretizePersonTime.sql
inst/sql/sql_server/GetAttrition.sql
inst/sql/sql_server/CreateResultsDataModel.sql
inst/csv/resultsDataModelSpecification.csv
man-roxygen/param_connection_details.R
man-roxygen/param_cdm_database_schema.R
man-roxygen/param_cohort_table.R
NEWS.md
vignettes/SingleStudy.Rmd         (replace existing stub)
vignettes/MultipleAnalyses.Rmd    (replace existing stub)
vignettes/StrategusModule.Rmd     (replace existing stub)
```

### Delete
```
R/Analysis.R                      (replaced by RunAnalyses.R)
vignettes/MethodChoice.Rmd        (content merged into SingleStudy)
vignettes/ResultsSchema.Rmd       (content merged into Export.R docs)
vignettes/clone-censor-weight-workflow.Rmd  (replaced by SingleStudy)
docs/hades-comparison.md          (superseded by vignettes)
```

---

## Phase 0 — Package housekeeping

1. Apply all `DESCRIPTION` changes above.
2. Create `NEWS.md`:
   ```markdown
   # CloneCensorWeight 0.1.0
   
   * Initial CRAN release.
   ```
3. Add `.Rbuildignore` entries for `extras/` and `docs/`.
4. In `R/zzz.R`, declare all dplyr column-name globals with
   `utils::globalVariables(c(...))`. Include every bare name used
   in `dplyr` verbs across all R files.
5. Rename SQL files to PascalCase as listed above.

---

## Phase 1 — `HelperFunctions.R`

Internal utility functions. Do **not** export any of these.

```r
# Null-coalescing operator (already in Cloning.R — move here, delete from Cloning.R)
`%||%` <- function(left, right) if (!is.null(left)) left else right

# Coerce a Date column in a data frame to Date class
.asDate <- function(x) as.Date(x)

# Low-level wrapper: render + execute a package SQL file
.renderTranslateExecuteSql <- function(connection, sqlFile, ...) { ... }

# Low-level wrapper: render + query a package SQL file, returns data frame
.renderTranslateQuerySql <- function(connection, sqlFile, snakeCaseToCamelCase = TRUE, ...) { ... }

# Resolve the SQL folder (handles installed vs development layout)
.sqlFolder <- function() {
  folder <- system.file("sql", "sql_server", package = "CloneCensorWeight")
  if (!nzchar(folder)) folder <- file.path("inst", "sql", "sql_server")
  folder
}

# Insert a temporary cohort table into the database, return table name "#ccw_cohorts"
.uploadTempCohortTable <- function(connection, cohorts) { ... }
```

---

## Phase 2 — `SettingsObjects.R`

Implement the full settings object hierarchy following CohortMethod's
`AbstractSerializableSettings` pattern.

### 2a — `AbstractSerializableSettings` base class

```r
AbstractSerializableSettings <- R6::R6Class(
  "AbstractSerializableSettings",
  public = list(
    initialize = function(typedList = NULL, untypedList = NULL, json = NULL) { ... },
    toList    = function() { ... },
    toJson    = function() { ... },
    validate  = function() invisible(self)
  )
)
```

`toJson()` wraps `jsonlite::toJSON(self$toList(), auto_unbox = TRUE, null = "null")`.
`initialize()` accepts exactly one of `typedList`, `untypedList`, or `json` and
populates public fields from it.

### 2b — `ArmDefinition`

```r
# Factory function (exported)
createArmDefinition <- function(
  armId,           # character(1), unique within an analysis
  strategyType,    # "initiateWithin" | "deferBeyond"
  gracePeriod   = 0,    # days from cohort start; required for both strategy types
  label         = NULL  # optional display label
)
```

Strategy type semantics (v1):
- `"initiateWithin"`: clone is censored if treatment has NOT occurred by day `gracePeriod`.
- `"deferBeyond"`: clone is censored if treatment occurs WITHIN day `gracePeriod`.
  Treatment after `gracePeriod` days is permitted and does not trigger censoring.
  This is the complement of `"initiateWithin"`; together they form a grace period
  comparison (early vs. deferred/no treatment).

`"sustainedUse"` is deferred to v1.1 (requires `exposure_end_date` in the data model).

R6 class `ArmDefinition` extends `AbstractSerializableSettings`.
`validate()` asserts: `armId` is non-empty character; `strategyType` is one of the
two permitted values; `gracePeriod >= 0`.

### 2c — `GetDbCloneCensorWeightDataArgs`

```r
createGetDbCloneCensorWeightDataArgs <- function(
  washoutPeriod          = 365,
  covariateSettings,           # FeatureExtraction covariateSettings object
  maxCohortSize          = 0   # 0 = no limit; sample if exceeded
)
```

Special `toList()` override preserves `covariateSettings` via
`.convertAttrToMember()` (copy the helper from CohortMethod).

### 2d — `CreateStudyPopulationArgs`

```r
createCreateStudyPopulationArgs <- function(
  removeSubjectsWithPriorOutcome = TRUE,
  priorOutcomeLookback           = 99999,
  minDaysAtRisk                  = 1,
  riskWindowStart                = 0,
  startAnchor                    = "cohort start",
  riskWindowEnd                  = 0,
  endAnchor                      = "cohort end"
)
```

### 2e — `CreateClonesAndCensorArgs`

```r
createCreateClonesAndCensorArgs <- function(
  arms,          # list of exactly 2 ArmDefinition objects
  intervalLength = 30   # person-time interval in days
)
```

`validate()` asserts exactly 2 arms with distinct `armId`.

### 2f — `CreateCensoringWeightsArgs`

```r
createCreateCensoringWeightsArgs <- function(
  artificialCensoringFormula = NULL,  # NULL → auto-build from covariates
  classicalCensoringFormula  = NULL,  # NULL → auto-build from covariates
  trim                       = c(0.01, 0.99),
  maxWeight                  = 10
)
```

### 2g — `FitOutcomeModelArgs`

```r
createFitOutcomeModelArgs <- function(
  modelType          = "cox",         # "cox" | "logistic"
  varianceEstimator  = "robust",      # "robust" | "bootstrap" | "conventional"
  nBootstraps        = 500,
  seed               = 1
)
```

### 2h — `CcwAnalysis` (the per-analysis specification)

```r
createCcwAnalysis <- function(
  analysisId,
  description                   = "",
  getDbCloneCensorWeightDataArgs,
  createStudyPopulationArgs,
  createClonesAndCensorArgs,
  createCensoringWeightsArgs,
  fitOutcomeModelArgs
)
```

R6 class `CcwAnalysis` extends `AbstractSerializableSettings`.
Save/load helpers:
```r
saveCcwAnalysisList(ccwAnalysisList, file)
loadCcwAnalysisList(file)
```

### 2i — `TargetArmOutcome`

```r
createTargetArmOutcome <- function(
  targetId,          # eligible population cohort id
  arm1CohortId,      # exposure cohort id for arm 1
  arm2CohortId,      # exposure cohort id for arm 2
  outcomeIds,
  excludedCovariateConceptIds = c(),
  includedCovariateConceptIds = c()
)
```

Save/load:
```r
saveTargetArmOutcomeList(targetArmOutcomeList, file)
loadTargetArmOutcomeList(file)
```

### 2j — `CcwAnalysesSpecifications` (top-level Strategus-compatible spec)

```r
CcwAnalysesSpecifications <- R6::R6Class(
  "CcwAnalysesSpecifications",
  inherit = AbstractSerializableSettings,
  public = list(
    ccwAnalysisList           = NULL,
    targetArmOutcomeList      = NULL,
    initialize                = function(...) { ... },
    validate                  = function() { ... }
  )
)

createCcwAnalysesSpecifications <- function(ccwAnalysisList, targetArmOutcomeList)
```

---

## Phase 3 — `CloneCensorWeightData.R` (refactor)

### 3a — Keep existing class definition and save/load unchanged.

### 3b — Full `getDbCloneCensorWeightData()` implementation

```r
getDbCloneCensorWeightData <- function(
  connectionDetails,
  cdmDatabaseSchema,
  tempEmulationSchema   = getOption("sqlRenderTempEmulationSchema"),
  targetDatabaseSchema,
  targetTable,
  targetId,
  outcomeDatabaseSchema,
  outcomeTable,
  outcomeIds,
  exposureDatabaseSchema,
  exposureTable,
  exposureId,
  covariateSettings,
  washoutPeriod         = 365,
  maxCohortSize         = 0,
  cdmVersion            = "5"
)
```

Implementation steps:
1. Connect, upload nothing yet.
2. Render + execute `GetCohorts.sql` into temp table `#ccw_cohorts`.
   - SQL joins `@target_table` (by `target_id`) with `@exposure_table` (by `exposure_id`)
     LEFT JOIN on `subject_id` where exposure `cohort_start_date` is within
     the washout-adjusted follow-up window.
   - Returns: `subject_id`, `cohort_start_date`, `cohort_end_date`, `initiation_date`
     (NULL if never initiated).
3. Render + query `GetOutcomes.sql`.
4. Call `FeatureExtraction::getDbCovariateData()` with `cohortTable = "#ccw_cohorts"`,
   `cohortTableIsTemp = TRUE`, `rowIdField = "subject_id"`.
5. If `maxCohortSize > 0` and n > `maxCohortSize`, sample subjects.
6. Assemble into `CloneCensorWeightData` via `createCloneCensorWeightData()`.
7. Drop `#ccw_cohorts`.
8. Disconnect.

### 3c — `createSampleCloneCensorWeightData()` — keep but add `covariateRef` and `analysisRef` from FeatureExtraction-compatible structure.

---

## Phase 4 — SQL files

### `GetCohorts.sql`

```sql
-- @target_database_schema, @target_table, @target_id
-- @exposure_database_schema, @exposure_table, @exposure_id
-- @washout_period (days)

SELECT
  t.subject_id,
  t.cohort_start_date,
  t.cohort_end_date,
  e.cohort_start_date AS initiation_date
FROM @target_database_schema.@target_table t
LEFT JOIN @exposure_database_schema.@exposure_table e
  ON  t.subject_id = e.subject_id
  AND e.cohort_definition_id = @exposure_id
  AND e.cohort_start_date >= t.cohort_start_date
  AND e.cohort_start_date <= t.cohort_end_date
WHERE t.cohort_definition_id = @target_id
;
```

### `GetOutcomes.sql`

```sql
-- @outcome_database_schema, @outcome_table, @outcome_ids (comma list)
SELECT
  subject_id,
  cohort_start_date AS outcome_date,
  cohort_definition_id AS outcome_id
FROM @outcome_database_schema.@outcome_table
WHERE cohort_definition_id IN (@outcome_ids)
;
```

### `GetAttrition.sql`

Returns counts at each attrition step.
Parameters: `@temp_table_name`, `@step_id`, `@description`.

### `CreateClones.sql`

Executed once per analysis after `#ccw_cohorts` exists.
Parameters: `@arm1_id`, `@arm1_strategy` (`'initiateWithin'|'neverInitiate'|'sustainedUse'`),
`@arm1_grace_period`, `@arm1_min_duration`,
`@arm2_id`, `@arm2_strategy`, `@arm2_grace_period`, `@arm2_min_duration`.

Logic (per arm):
- For `initiateWithin`: censor day = `gracePeriod` if `initiation_date IS NULL`
  or `DATEDIFF(day, cohort_start_date, initiation_date) > gracePeriod`;
  else follow-up end = natural end.
- For `deferBeyond`: censor day = `DATEDIFF(day, cohort_start_date, initiation_date)`
  if `initiation_date IS NOT NULL` and that day <= `gracePeriod`; else natural end.
  Initiation after `gracePeriod` days does NOT trigger censoring.
- `"sustainedUse"` is not implemented in v1.

Output temp table `#ccw_clones`:
`subject_id`, `clone_id` (varchar), `arm_id` (varchar), `cohort_start_date`,
`cohort_end_date` (natural or censored), `artificially_censored` (0/1).

### `DiscretizePersonTime.sql`

Uses a recursive CTE or tally CTE to generate interval numbers.
Parameters: `@interval_length`, `@max_follow_up_days` (computed from data).

Joins `#ccw_clones` with interval table.
Output temp table `#ccw_person_time`:
`subject_id`, `clone_id`, `arm_id`, `interval_number`,
`t_start`, `t_stop`, `person_time_days`,
`censor_in_interval` (1 if artificially censored in this interval),
`classical_censor_in_interval` (1 if observation end in this interval — computed from
original `cohort_end_date` vs natural maximum).

---

## Phase 5 — `StudyPopulation.R`

```r
createStudyPopulation <- function(
  cloneCensorWeightData,
  outcomeId,
  createStudyPopulationArgs
)
```

Applies:
- `removeSubjectsWithPriorOutcome`: exclude subjects with outcome before `cohortStartDate`
  within `priorOutcomeLookback` days.
- `minDaysAtRisk`: exclude subjects with follow-up < threshold.
- Risk window trimming: `riskWindowStart`/`riskWindowEnd` with anchor.

Returns `cloneCensorWeightData` with `cohorts` filtered and an attrition attribute.

```r
getAttritionTable <- function(population)
```

Returns a data frame with `sequenceNumber`, `description`, `subjects`.

---

## Phase 6 — `Cloning.R` (refactor)

`createClonesAndCensor()` now:
1. Validates that `cloneCensorWeightData` has `cohorts` with required columns.
2. Validates `arms` as a list of exactly 2 `ArmDefinition` objects.
3. Uploads `cohorts` to `#ccw_cohorts` (if not already a DB analysis — for
   pure in-R mode, keep existing logic as fallback for `createSampleCloneCensorWeightData`
   workflows).
4. Executes `CreateClones.sql` with arm parameters.
5. Reads `#ccw_clones` back into R and stores in `cloneCensorWeightData[["cohorts"]]`.
6. Drops temp table.

For the **in-memory/sample-data path** (no DB connection), retain the existing
pure-R logic as an internal function `_.createClonesR()` called when
`connection` is `NULL`.

---

## Phase 7 — `Weighting.R` (refactor)

### 7a — `discretizePersonTime()`

When a `connection` is available: executes `DiscretizePersonTime.sql`.
When `connection` is NULL (in-memory): keeps current R-side loop as fallback.
Returns a `tbl_df` / data frame with columns:
`subjectId`, `cloneId`, `armId`, `intervalNumber`, `tStart`, `tStop`,
`personTimeDays`, `censorInInterval`, `classicalCensorInInterval`, `event`.

### 7b — `createCensoringWeights()`

Two separate pooled logistic models:

**Artificial censoring weight:**
- Outcome: `censorInInterval == 0` (uncensored at this interval)
- Numerator formula: `~ factor(armId) + ns(intervalNumber, df = 3)`
- Denominator formula: adds all baseline covariate columns

**Classical censoring weight:**
- Outcome: `classicalCensorInInterval == 0` (in follow-up at this interval)  
- Same numerator/denominator structure

Combined weight per interval = `artificialWeight * classicalWeight`.
Cumulative weight = `cumprod()` over intervals within `cloneId`.

Both models use Cyclops with a mild Laplace prior (variance = 0.1).
GLM fallback if Cyclops fails.

Stabilized form:
`w_artificial(t) = P(uncensored | arm, t) / P(uncensored | arm, t, covariates)`

`trimWeights()` clips at specified percentiles AFTER computing cumulative weights.

Function signature:
```r
createCensoringWeights <- function(
  personTime,
  cloneCensorWeightData,
  artificialCensoringFormula = NULL,
  classicalCensoringFormula  = NULL,
  trim = c(0.01, 0.99),
  maxWeight = 10
)
```

Returns `personTime` data frame with added columns:
`artificialWeight`, `classicalWeight`, `weight` (combined, cumulative, trimmed).
Attaches model objects and engine names as attributes.

### 7c — Existing diagnostic functions

Keep `computeWeightDistribution()`, `plotWeightDistribution()`, `trimWeights()`
with minor signature updates to match new column names.

---

## Phase 8 — `OutcomeModel.R` (refactor)

```r
fitOutcomeModel <- function(
  personTime,
  modelType         = "cox",       # "cox" | "logistic"
  varianceEstimator = "robust",    # "robust" | "bootstrap" | "conventional"
  nBootstraps       = 500,
  seed              = 1
)
```

### Point estimate

Unchanged from current implementation: Cyclops Cox or logistic with
`weights = personTime$weight`.

### Variance: `"robust"` (Lin-Wei sandwich)

For Cox: use `survival::coxph(Surv(tStop, event) ~ arm, data, weights, cluster = subjectId)`.
Extract `coef()` and `vcov()` for CI.

For logistic: use `sandwich::sandwich()` on a `glm` fit with `cluster = subjectId`.
Note: `sandwich` is a suggested package; fall back to conventional SE if not available
and warn.

### Variance: `"bootstrap"` (subject-level)

```
for (b in 1:nBootstraps) {
  subjects_b <- sample(unique(personTime$subjectId), replace = TRUE)
  personTime_b <- personTime[personTime$subjectId %in% subjects_b, ]
  # re-run createCensoringWeights() + fitOutcomeModel(varianceEstimator = "conventional")
  # store log(estimate_b)
}
CI = exp(quantile(log_estimates, c(0.025, 0.975)))
```

### Return value

`CloneCensorWeightOutcomeModel` S3 list (unchanged class name) with added fields:
`varianceEstimator`, `nBootstraps`, `bootstrapDistribution` (NULL unless bootstrap).

---

## Phase 9 — `Balance.R`

```r
computeCovariateBalance <- function(
  cloneCensorWeightData,
  personTime = NULL,    # if supplied, compute weighted balance
  targetArmId,
  comparatorArmId
)
```

Returns data frame with columns matching `ccw_covariate_balance` schema:
`covariateId`, `covariateName`, `targetMeanBefore`, `comparatorMeanBefore`,
`stdDiffBefore`, `targetMeanAfter`, `comparatorMeanAfter`, `stdDiffAfter`.

Standardized mean difference formula (Austin 2009):
`SMD = (mean_target - mean_comparator) / sqrt((sd_target^2 + sd_comparator^2) / 2)`

```r
plotCovariateBalance <- function(balance, beforeLabel = "Before weighting",
                                 afterLabel = "After weighting")
```

Love plot style (scatter of SMD before vs after). Requires `ggplot2`.

---

## Phase 10 — `KaplanMeier.R`

```r
computeWeightedKaplanMeier <- function(personTime, targetArmId, comparatorArmId)
```

Uses `survival::survfit()` with `weights = personTime$weight` and
`cluster = subjectId` for confidence intervals.
Returns data frame with columns matching `ccw_kaplan_meier_dist` schema.

```r
plotKaplanMeier <- function(kmDist, targetLabel = NULL, comparatorLabel = NULL)
```

Step function plot with shaded CI bands. Requires `ggplot2`.

---

## Phase 11 — `Power.R`

```r
computeMdrr <- function(
  personTime,
  targetArmId,
  comparatorArmId,
  alpha     = 0.05,
  power     = 0.8,
  twoSided  = TRUE,
  modelType = "cox"
)
```

Follow CohortMethod's `computeMdrr()` logic exactly but adapted to CCW person-time
(effective sample size uses the sum of weights, not raw counts).

---

## Phase 12 — `EmpiricalCalibration.R`

```r
calibrateCcwResults <- function(
  outcomeModels,        # named list of CloneCensorWeightOutcomeModel
  negativeControlIds,   # outcomeIds designated as negative controls
  positiveControlIds    = NULL,
  alpha                 = 0.05
)
```

Implementation:
1. Extract `log_rr` and `se_log_rr` for negative control outcomes.
2. Call `EmpiricalCalibration::fitNull()` to estimate null distribution.
3. Call `EmpiricalCalibration::calibrateP()` for calibrated p-values.
4. Call `EmpiricalCalibration::calibrateConfidenceInterval()` for calibrated CIs.
5. Add calibrated columns to each `outcomeModel$summary`.

Returns the same list with `summary` extended by:
`calibratedRr`, `calibratedCi95Lb`, `calibratedCi95Ub`, `calibratedP`,
`calibratedLogRr`, `calibratedSeLogRr`.

---

## Phase 13 — `ResultsDataModel.R` and schema files

### `inst/csv/resultsDataModelSpecification.csv`

Create with the following tables (snake_case column names):

| Table | Key columns |
|---|---|
| `ccw_analysis` | `analysis_id` (PK), `description`, `definition` |
| `ccw_arm` | `target_arm_id` (PK), `target_id`, `arm1_cohort_id`, `arm2_cohort_id`, `nesting_cohort_id` |
| `ccw_arm_outcome` | `target_arm_id` (PK), `outcome_id` (PK), `outcome_of_interest`, `true_effect_size` |
| `ccw_attrition` | `sequence_number` (PK), `description`, `subjects` (min_cell), `analysis_id` (PK), `target_arm_id` (PK), `outcome_id` (PK), `database_id` (PK) |
| `ccw_follow_up_dist` | `target_arm_id` (PK), `outcome_id` (PK), `analysis_id` (PK), `database_id` (PK), `arm1_min_days`, `arm1_p10_days`, `arm1_median_days`, `arm1_p90_days`, `arm1_max_days`, `arm2_*` (same) |
| `ccw_result` | `analysis_id` (PK), `target_arm_id` (PK), `outcome_id` (PK), `database_id` (PK), `rr`, `ci_95_lb`, `ci_95_ub`, `p`, `one_sided_p`, `arm1_subjects` (min_cell), `arm2_subjects` (min_cell), `arm1_days`, `arm2_days`, `arm1_outcomes` (min_cell), `arm2_outcomes` (min_cell), `log_rr`, `se_log_rr`, `llr`, `calibrated_rr`, `calibrated_ci_95_lb`, `calibrated_ci_95_ub`, `calibrated_p`, `calibrated_log_rr`, `calibrated_se_log_rr`, `target_estimator` |
| `ccw_weight_distribution` | `analysis_id` (PK), `target_arm_id` (PK), `outcome_id` (PK), `database_id` (PK), `arm_id` (PK), `n`, `mean_weight`, `sd_weight`, `min_weight`, `p01_weight`, `median_weight`, `p99_weight`, `max_weight`, `effective_sample_size` |
| `ccw_covariate_analysis` | `covariate_analysis_id` (PK), `covariate_analysis_name`, `analysis_id` (PK), `database_id` (PK) |
| `ccw_covariate` | `covariate_id` (PK), `covariate_name`, `analysis_id` (PK), `covariate_analysis_id`, `database_id` (PK) |
| `ccw_covariate_balance` | `database_id` (PK), `target_arm_id` (PK), `outcome_id` (PK), `analysis_id` (PK), `covariate_id` (PK), `target_mean_before` (min_cell), `comparator_mean_before` (min_cell), `std_diff_before`, `target_mean_after` (min_cell), `comparator_mean_after` (min_cell), `std_diff_after` |
| `ccw_kaplan_meier_dist` | `time_day` (PK), `arm1_survival`, `arm1_survival_lb`, `arm1_survival_ub`, `arm2_survival`, `arm2_survival_lb`, `arm2_survival_ub`, `arm1_at_risk` (min_cell), `arm2_at_risk` (min_cell), `target_arm_id` (PK), `outcome_id` (PK), `analysis_id` (PK), `database_id` (PK) |
| `ccw_likelihood_profile` | `log_rr` (PK), `log_likelihood`, `gradient`, `target_arm_id` (PK), `outcome_id` (PK), `analysis_id` (PK), `database_id` (PK) |
| `ccw_diagnostics_summary` | `analysis_id` (PK), `target_arm_id` (PK), `outcome_id` (PK), `database_id` (PK), `max_sdm`, `mdrr`, `ease`, `balance_diagnostic`, `mdrr_diagnostic`, `ease_diagnostic`, `unblind`, `unblind_for_evidence_synthesis` |

### `inst/sql/sql_server/CreateResultsDataModel.sql`

Standard `CREATE TABLE IF NOT EXISTS` DDL for each table above.
Use `SqlRender`-compatible dialect (no database-specific syntax).

### `R/ResultsDataModel.R`

```r
getResultsDataModelSpecifications <- function() {
  readr::read_csv(system.file("csv", "resultsDataModelSpecification.csv",
                              package = "CloneCensorWeight"),
                  show_col_types = FALSE)
}

createResultsDataModel <- function(connectionDetails, databaseSchema,
                                   tablePrefix = "") { ... }
```

---

## Phase 14 — `Export.R`

```r
exportToCsv <- function(
  outputFolder,
  databaseId,
  databaseName         = databaseId,
  databaseDescription  = databaseId,
  minCellCount         = 5,
  maxCores             = 1
)
```

Reads result RDS/CSV files from `outputFolder` (written by `runCcwAnalyses()`).
Applies `minCellCount` censoring to all columns flagged `min_cell_count = Yes`
in the spec CSV.
Writes one CSV per table into `outputFolder/export/`.
Zips all CSVs to `outputFolder/export/Results_{databaseId}.zip`.

---

## Phase 15 — `RunAnalyses.R`

### Multi-threading settings

```r
createMultiThreadingSettings <- function(
  getDbCcwDataThreads        = 1,
  createStudyPopThreads      = 1,
  createClonesThreads        = 1,
  createWeightsThreads       = 1,
  fitOutcomeModelThreads     = 1,
  calibrationThreads         = 1
)

createDefaultMultiThreadingSettings <- function(maxCores)
```

### Main runner

```r
runCcwAnalyses <- function(
  connectionDetails,
  cdmDatabaseSchema,
  tempEmulationSchema     = getOption("sqlRenderTempEmulationSchema"),
  targetDatabaseSchema,
  targetTable,
  outcomeDatabaseSchema,
  outcomeTable,
  outputFolder,
  ccwAnalysesSpecifications,
  multiThreadingSettings  = createMultiThreadingSettings()
)
```

Execution pattern follows CohortMethod exactly:
1. Loop over unique (targetId, arm1CohortId, arm2CohortId, covariateSettings) tuples.
   Cache `CloneCensorWeightData` objects by hash; skip if already on disk.
2. Loop over (cachedData, analysisId, outcomeId) tuples in parallel.
   - `createStudyPopulation()`
   - `createClonesAndCensor()`
   - `discretizePersonTime()`
   - `createCensoringWeights()`
   - `fitOutcomeModel()`
   - `computeCovariateBalance()`
   - `computeWeightedKaplanMeier()`
   - `computeMdrr()`
3. Run `calibrateCcwResults()` once per (targetArmId) after all outcomes complete.
4. Write results to RDS files under `outputFolder/{analysisId}_{targetArmId}_{outcomeId}/`.
5. Return an `CcwAnalysisReference` data frame mapping all (analysisId, targetArmId, outcomeId)
   to their output folders and convergence status.

```r
summarizeCcwAnalyses <- function(referenceTable, outputFolder)
```

Reads all result files and returns a single summary data frame — one row per
(analysisId, targetArmId, outcomeId) with point estimate, CI, p-value, and calibrated equivalents.

```r
getCcwAnalysisResults <- function(referenceTable, analysisId, targetArmId, outcomeId, outputFolder)
```

Returns the full `CloneCensorWeightOutcomeModel` for one specific combination.

---

## Phase 16 — `Simulation.R`

```r
createCloneCensorWeightDataSimulationProfile <- function(
  cloneCensorWeightData,
  outcomeId,
  minCellCount = 5
)
```

Mirrors CohortMethod's simulation profile approach:
- Covariate prevalences (binary covariates only, suppressed below `minCellCount`)
- Temporal rates: `cohortEndRate`, `initiationRate` modeled as exponential
- Outcome model coefficients from a fitted Cyclops Poisson model
- Pre-index outcome rates

```r
simulateCloneCensorWeightData <- function(
  simulationProfile,
  n = 1000,
  seed = 1
)
```

Returns a `CloneCensorWeightData` object with realistic characteristics.
Store a built-in `ccwDataSimulationProfile` data object in `data/` for use
in examples and tests.

---

## Phase 17 — `Strategus.R`

Standalone R6 class — do **not** import or depend on the `Strategus` package.
The class mirrors the interface Strategus expects (`toList()`, `toJson()`,
`$moduleType`) so a Strategus runner can deserialize it, but `Strategus` is
not listed in `Imports` or `Suggests`.

```r
CcwModuleSpecifications <- R6::R6Class(
  "CcwModuleSpecifications",
  public = list(
    ccwAnalysesSpecifications = NULL,
    initialize = function(ccwAnalysesSpecifications) { ... },
    toList     = function() { ... },
    toJson     = function() { ... }
  )
)

createCcwModuleSpecifications <- function(ccwAnalysesSpecifications)
```

Module execution function (called by Strategus runner):
```r
runCcwModule <- function(
  jobContext    # Strategus job context with connectionDetails, analysisSpecifications, etc.
)
```

This function:
1. Deserializes `ccwAnalysesSpecifications` from JSON in `jobContext`.
2. Calls `runCcwAnalyses()` with connection info from `jobContext`.
3. Calls `exportToCsv()` to write standardized outputs.
4. Does NOT extract data (that is a pre-step in the Strategus pipeline).

---

## Phase 18 — CRAN compliance

Work through this checklist in order:

1. **`utils::globalVariables()`** in `R/zzz.R`: declare every bare name used in
   `dplyr` verbs and `rlang::.data` references across all files.

2. **Documentation completeness**: every exported function must have
   `@param`, `@return`, and `@examples`. Database-dependent examples use
   `\dontrun{}`. Pure-R examples (e.g., `createArmDefinition`, `createSampleCloneCensorWeightData`)
   must be runnable without wrapping.

3. **`man-roxygen/` templates**: create shared `@param` files for
   `param_connection_details.R`, `param_cdm_database_schema.R`,
   `param_cohort_table.R`, `param_temp_emulation_schema.R`.
   Reference them with `@template`.

4. **No `T`/`F`**: search and replace all `T` and `F` used as logical
   constants with `TRUE` and `FALSE`.

5. **No `require()` inside functions**: replace any `require()` with
   `requireNamespace(..., quietly = TRUE)`.

6. **No global state mutations**: no `options()`, `par()`, or `setwd()` calls
   without restoring in an `on.exit()`.

7. **No writing outside `tempdir()`** in tests: all test files that write to disk
   must use `withr::local_tempdir()` or `tempfile()`.

8. **`LazyData: false`** in DESCRIPTION (since we have `data/` objects, they should
   be lazy-loaded; set `LazyData: true`). Actually: if using `data/` objects,
   set `LazyData: true`.

9. **`NAMESPACE`**: regenerated by `devtools::document()`. Do not edit manually.

10. **Encoding**: all files must be UTF-8. Author string uses Unicode for special
    characters if needed (none in this case).

11. **`inst/CITATION`**: add a CITATION file referencing the PMC paper
    (Hernan et al. approach) and the package itself.

12. **`R CMD check --as-cran`**: must pass with 0 errors, 0 warnings.
    Acceptable: one NOTE for "New submission".

---

## Test strategy

### File: `tests/testthat/test-settings-objects.R`
- Round-trip JSON serialization for each `createXxx` factory function.
- `validate()` catches invalid inputs.
- `ArmDefinition` rejects `sustainedUse` without `minDuration`.

### File: `tests/testthat/test-sample-pipeline.R` (update existing)
- Full pipeline using `createSampleCloneCensorWeightData()`:
  `createStudyPopulation` → `createClonesAndCensor` → `discretizePersonTime` →
  `createCensoringWeights` → `fitOutcomeModel` (both `"robust"` and `"bootstrap"`)
- Weight diagnostics produce finite positive weights.
- Both arm strategies (`"initiateWithin"` and `"neverInitiate"`) create valid clones.
- Bootstrap SE is wider than or equal to conventional SE.

### File: `tests/testthat/test-balance.R`
- SMD before weighting is non-zero when covariates differ.
- SMD after weighting is closer to zero.
- `plotCovariateBalance()` returns a `ggplot` object.

### File: `tests/testthat/test-export.R`
- `exportToCsv()` creates a ZIP with the expected table files.
- `minCellCount` censoring applies to flagged columns.

### File: `tests/testthat/test-public-api.R` (update existing)
- `runCcwAnalyses()` end-to-end with sample data writes expected files.
- `summarizeCcwAnalyses()` returns one row per analysis-outcome.

### General rules
- No internet access in any test.
- All disk writes use `withr::local_tempdir()`.
- No `skip_on_cran()` for pure-R tests.
- Use `skip_if_not_installed("Eunomia")` / `skip_if_not_installed("ggplot2")` for
  optional-dependency tests.

---

## Vignette outline

### `vignettes/SingleStudy.Rmd`
1. Background: what CCW answers, the three steps
2. Setting up a `TargetArmOutcome`
3. Defining arm strategies with `createArmDefinition()`
4. Running `getDbCloneCensorWeightData()` (code in `\dontrun{}`)
5. `createStudyPopulation()`, `createClonesAndCensor()`, `discretizePersonTime()`
6. `createCensoringWeights()` — inspecting weight distributions
7. `fitOutcomeModel()` — reading the result
8. `computeCovariateBalance()`, `plotCovariateBalance()`
9. `plotKaplanMeier()`
10. Empirical calibration

### `vignettes/MultipleAnalyses.Rmd`
1. Creating a list of `CcwAnalysis` objects
2. Running `runCcwAnalyses()`
3. Summarizing results with `summarizeCcwAnalyses()`
4. Exporting with `exportToCsv()`

### `vignettes/StrategusModule.Rmd`
1. Creating `CcwModuleSpecifications`
2. Serializing to JSON
3. Integrating into a Strategus study package

---

## Notes for the implementing agent

- Never use `T` or `F` as logical literals.
- `dplyr` column references must always use `.data$columnName` or be declared in
  `utils::globalVariables()`.
- All SQL files must be `SqlRender`-parameterized (no string pasting of schema names).
- Cyclops is the primary fitting engine; GLM is only a fallback. Always record which
  engine was used in the result metadata.
- The `connection` argument to internal functions should always be paired with
  `on.exit(DatabaseConnector::disconnect(connection), add = TRUE)` in the calling
  function.
- When a function both modifies `cloneCensorWeightData` and returns it, it should
  return the modified object invisibly if the primary output is a side effect
  (e.g., adding a table). When it returns a new derived object (e.g., `personTime`
  data frame), return visibly.
- Robust SE for Cox models: use `survival::coxph(..., cluster = subjectId)` with
  `subjectId` as the clustering variable. This gives the Lin-Wei-Ying sandwich
  variance that accounts for within-subject clone correlation.
- Bootstrap: always resample at the **subject** level, not the clone level, to
  preserve the dependence structure.
- `EmpiricalCalibration::fitNull()` expects a vector of `logRr` and `seLogRr`
  from negative controls. The function returns a null distribution object.
  Apply `EmpiricalCalibration::calibrateP()` per outcome after fitting.

# CloneCensorWeight

> **Experimental — under active development.** APIs may change without notice
> and the package has not yet been validated for production use. Feedback and
> issue reports are welcome.

<!-- badges: start -->
[![R-CMD-check](https://github.com/ohdsi/CloneCensorWeight/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ohdsi/CloneCensorWeight/actions/workflows/R-CMD-check.yaml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
<!-- badges: end -->

**CloneCensorWeight** is an OHDSI R package for
**clone-censor-weight (CCW) target trial emulation** in the
[OMOP Common Data Model](https://ohdsi.github.io/CommonDataModel/).

CCW eliminates immortal time bias by cloning every subject into both treatment
arms at baseline, then using inverse probability of censoring weighting (IPCW)
to correct for the informative nature of the artificial censoring that follows.
The result is an unbiased per-protocol estimate of a treatment effect from
observational data.

---

## The CCW Pipeline

```
            ┌──────────────────────────────┐
            │   getDbCloneCensorWeightData  │  Extract from OMOP CDM
            └──────────────┬───────────────┘
                           │
            ┌──────────────▼───────────────┐
            │     createStudyPopulation    │  Washout · eligibility · risk window
            └──────────────┬───────────────┘
                           │
            ┌──────────────▼───────────────┐
            │     createClonesAndCensor    │  Clone × arm; artificial censoring
            └──────────────┬───────────────┘
                           │
            ┌──────────────▼───────────────┐
            │     discretizePersonTime     │  Interval-format person-time rows
            └──────────────┬───────────────┘
                           │
            ┌──────────────▼───────────────┐
            │    createCensoringWeights    │  Fit IPCW models per arm
            └──────────────┬───────────────┘
                           │
            ┌──────────────▼───────────────┐
            │       fitOutcomeModel        │  Cox / logistic with robust SE
            └──────────────┬───────────────┘
                           │
            ┌──────────────▼───────────────┐
            │     calibrateEstimates       │  Empirical calibration (optional)
            └──────────────────────────────┘
```

---

## Installation

```r
install.packages("remotes")
remotes::install_github("ohdsi/CloneCensorWeight")
```

---

## Quick Start

The example below emulates a target trial of recombinant zoster vaccination
(RZV) in new JAK inhibitor users, using Capr to define cohorts and
CohortGenerator to instantiate them.  See `vignette("SingleStudy")` for the
full walkthrough.

### 1. Define cohorts with Capr

```r
library(CloneCensorWeight)
library(Capr)
library(CohortGenerator)

jakInhibitors <- cs(descendants(42873017L, 1510448L, 2213473L), name = "JAK inhibitors")
rzv           <- cs(descendants(1160056L, 1160057L),            name = "RZV (Shingrix)")
herpesZoster  <- cs(descendants(4067760L),                      name = "Herpes zoster")

jakNewUsers <- cohort(
  entry     = entry(drugExposure(jakInhibitors, firstOccurrence()),
                    observationWindow = continuousObservation(priorDays = 365L)),
  attrition = attrition(
    "No prior JAK" = withAll(exactly(0, drugExposure(jakInhibitors),
                                     duringInterval(eventStarts(-365L, -1L)))),
    "No prior RZV" = withAll(exactly(0, drugExposure(rzv),
                                     duringInterval(eventStarts(-1825L, -1L))))
  ),
  exit = exit(endStrategy = observationExit())
)

rzvCohort          <- cohort(entry = entry(drugExposure(rzv, firstOccurrence())),
                             exit  = exit(endStrategy = fixedExit("startDate", 0L)))
herpesZosterCohort <- cohort(entry = entry(conditionOccurrence(herpesZoster, firstOccurrence())),
                             exit  = exit(endStrategy = fixedExit("startDate", 0L)))

# Cohort IDs: 1 = JAK new users, 2 = RZV, 3 = herpes zoster
allCohorts <- makeCohortSet(jakNewUsers, rzvCohort, herpesZosterCohort)
```

### 2. Instantiate cohorts on the CDM

```r
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql", server = "localhost/ohdsi",
  user = "joe", password = "secret"
)
cdmDatabaseSchema    <- "cdm_v54"
cohortDatabaseSchema <- "scratch"
cohortTable          <- "ccw_jak_rzv"

cohortTableNames <- CohortGenerator::getCohortTableNames(cohortTable)
CohortGenerator::createCohortTables(connectionDetails, cohortDatabaseSchema, cohortTableNames)
CohortGenerator::generateCohortSet(connectionDetails, cdmDatabaseSchema,
                                   cohortDatabaseSchema, cohortTableNames, allCohorts)
```

### 3. Specify the analysis

```r
# Two arms: vaccinate within 90 days vs. defer
armVaccinate <- createArmDefinition("rzv_within_90",    "initiateWithin", gracePeriod = 90L)
armDeferred  <- createArmDefinition("no_rzv_within_90", "deferBeyond",    gracePeriod = 90L)

analysis <- createCcwAnalysis(
  analysisId                     = "main_90day",
  description                    = "RZV within 90 days vs. deferred",
  getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(
    washoutPeriod     = 365L,
    covariateSettings = FeatureExtraction::createDefaultCovariateSettings()
  ),
  createStudyPopulationArgs  = createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = TRUE,
    riskWindowEnd = 730L, endAnchor = "cohort start"
  ),
  createClonesAndCensorArgs  = createCreateClonesAndCensorArgs(
    arms = list(armVaccinate, armDeferred), intervalLength = 30L
  ),
  createCensoringWeightsArgs = createCreateCensoringWeightsArgs(trim = c(0.01, 0.99)),
  fitOutcomeModelArgs        = createFitOutcomeModelArgs("cox", "robust")
)

# Link cohort IDs to the analysis
# indexCohortId = eligible population; targetCohortId / comparatorCohortId = exposure cohort
# (both arms share the same exposure cohort — the arm strategy distinguishes them)
tao <- createTargetComparatorOutcome(
  indexCohortId      = 1L,  # JAK new users
  targetCohortId     = 2L,  # RZV
  comparatorCohortId = 2L,  # RZV (deferBeyond arm)
  outcomeIds         = 3L   # herpes zoster
)

specs <- createCcwAnalysesSpecifications(list(analysis), list(tao))
```

### 4. Run

```r
runCcwAnalyses(
  connectionDetails         = connectionDetails,
  cdmDatabaseSchema         = cdmDatabaseSchema,
  targetDatabaseSchema      = cohortDatabaseSchema,
  targetTable               = cohortTableNames$cohortTable,
  outcomeDatabaseSchema     = cohortDatabaseSchema,
  outcomeTable              = cohortTableNames$cohortTable,
  exposureDatabaseSchema    = cohortDatabaseSchema,
  exposureTable             = cohortTableNames$cohortTable,
  ccwAnalysesSpecifications = specs,
  outputFolder              = "jak_rzv_results",
  databaseId                = "synpuf"
)

results <- importFromCsv("jak_rzv_results")
results$ccw_outcome_model[, c("analysis_id", "estimate", "lower_95", "upper_95", "p_value")]
```

---

## Vignettes

| Vignette | Description |
|---|---|
| [Single Study](vignettes/SingleStudy.Rmd) | End-to-end walkthrough: JAK inhibitor + recombinant zoster vaccine |
| [Multiple Analyses](vignettes/MultipleAnalyses.Rmd) | Sensitivity analyses with `runCcwAnalysesOnData()` |
| [Simulation Studies](vignettes/Simulation.Rmd) | Power analysis and pipeline validation with `simulateCcwData()` |
| [Developer Guide](vignettes/DeveloperGuide.Rmd) | Internal architecture and data structures at each pipeline step |

---

## License

Apache License 2.0 — see [LICENSE](LICENSE).

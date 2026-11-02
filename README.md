# CloneCensorWeight

> **Experimental — under active development.** APIs may change without notice
> and the package has not yet been validated for production use. Feedback and
> issue reports are welcome.

<!-- badges: start -->
[![R-CMD-check](https://github.com/ohdsi/CloneCensorWeight/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ohdsi/CloneCensorWeight/actions/workflows/R-CMD-check.yaml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
<!-- badges: end -->

**CloneCensorWeight** is an OHDSI [HADES](https://ohdsi.github.io/Hades/) R package for
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

## Features

| Capability | Functions |
|---|---|
| OMOP CDM data extraction | `getDbCloneCensorWeightData()` |
| Flexible arm definitions | `createArmDefinition()` — `initiateWithin` / `deferBeyond` strategies |
| Inverse probability of censoring weighting | `createCensoringWeights()`, `trimWeights()` |
| Outcome models | Cox proportional hazards · logistic regression; robust or bootstrap SE |
| Covariate balance diagnostics | `computeStandardizedDifferences()`, `plotLovePlot()` |
| Weighted Kaplan–Meier curves | `computeKaplanMeier()`, `plotKaplanMeier()` |
| Minimum detectable relative risk | `computeMdrr()`, `plotMdrr()` |
| Empirical calibration | `calibrateEstimates()` with negative control outcomes |
| Multi-analysis runner with resume | `runCcwAnalyses()`, `runCcwAnalysesOnData()` |
| Simulation & power analysis | `simulateCcwData()`, `runSimulation()`, `summarizeSimulation()` |
| Reproducible JSON specifications | `saveCcwAnalysisList()`, `loadCcwAnalysisList()` |
| Results database export | `exportToCsv()`, `uploadResults()` |

---

## Installation

```r
install.packages("remotes")
remotes::install_github("ohdsi/CloneCensorWeight")
```

---

## Quick Start

### With sample data (no database required)

```r
library(CloneCensorWeight)

# 1. Load built-in synthetic data
ccwData <- createSampleCloneCensorWeightData(n = 500, seed = 42)

# 2. Define the two treatment arms
arms <- list(
  createArmDefinition("vaccinate", "initiateWithin", gracePeriod = 90),
  createArmDefinition("defer",     "deferBeyond",    gracePeriod = 90)
)

# 3. Build step-level settings objects
popArgs     <- createCreateStudyPopulationArgs(riskWindowEnd = 365, endAnchor = "cohort start")
cloneArgs   <- createCreateClonesAndCensorArgs(arms = arms, intervalLength = 30)
weightsArgs <- createCreateCensoringWeightsArgs()
modelArgs   <- createFitOutcomeModelArgs("cox", "robust")

# 4. Run the pipeline
createStudyPopulation(ccwData, outcomeId = 1L, popArgs)
cloned   <- createClonesAndCensor(ccwData, cloneArgs)
pt       <- discretizePersonTime(cloned, intervalLength = 30)
weighted <- createCensoringWeights(pt, cloned, weightsArgs)
model    <- fitOutcomeModel(weighted, fitOutcomeModelArgs = modelArgs)

print(model$summary)
#>   estimate lower95 upper95 pValue
#>      0.812   0.651   1.013  0.065
```

### With a real CDM

Cohort IDs are specified when extracting data from the database.
`indexCohortId` is the eligible population (e.g. new users of a drug of
interest); `targetCohortId` / `comparatorCohortId` identify the treatment
exposure cohort whose `cohort_start_date` marks treatment initiation.

```r
ccwData <- getDbCloneCensorWeightData(
  connectionDetails      = connectionDetails,
  cdmDatabaseSchema      = "cdm",
  targetDatabaseSchema   = "results",
  targetTable            = "cohort",
  targetId               = 101L,   # indexCohortId  — eligible population
  exposureDatabaseSchema = "results",
  exposureTable          = "cohort",
  exposureId             = 102L,   # targetCohortId — treatment initiation events
  outcomeDatabaseSchema  = "results",
  outcomeTable           = "cohort",
  outcomeIds             = c(201L, 301L:320L),
  getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(
    washoutPeriod     = 365,
    covariateSettings = FeatureExtraction::createDefaultCovariateSettings()
  )
)
```

The rest of the pipeline is identical to the sample-data example above.

---

## Running Multiple Analyses

`runCcwAnalysesOnData()` iterates over every combination of target cohort ×
analysis specification × outcome, appending results to CSV files with resume
support.

```r
analysisList <- list(
  createCcwAnalysis("gp30", "Grace period 30 d",
    createStudyPopulationArgs  = popArgs,
    createClonesAndCensorArgs  = createCreateClonesAndCensorArgs(arms = arms, intervalLength = 30),
    createCensoringWeightsArgs = weightsArgs,
    fitOutcomeModelArgs        = modelArgs
  ),
  createCcwAnalysis("gp90", "Grace period 90 d",
    createStudyPopulationArgs  = popArgs,
    createClonesAndCensorArgs  = createCreateClonesAndCensorArgs(arms = arms, intervalLength = 30),
    createCensoringWeightsArgs = weightsArgs,
    fitOutcomeModelArgs        = modelArgs
  )
)

tao <- createTargetComparatorOutcome(
  indexCohortId      = 1L,   # eligible population cohort (sample data uses 1)
  targetCohortId     = 1L,   # treatment-exposure cohort
  comparatorCohortId = 1L,   # comparator cohort (same ID; arm strategy distinguishes them)
  outcomeIds         = 1L
)
specs <- createCcwAnalysesSpecifications(analysisList, list(tao))

runCcwAnalysesOnData(ccwData, specs, outputFolder = "results", databaseId = "synpuf")

results <- importFromCsv("results")
results$ccw_outcome_model[, c("analysis_id", "estimate", "lower_95", "upper_95")]
```

---

## Simulation & Power Analysis

```r
# Check operating characteristics before running on real data
sim <- runSimulation(
  nSims           = 200,
  n               = 500,
  trueHazardRatio = 0.75,
  baselineHazard  = 0.02,
  confoundingStrength = 0.5,
  seed = 1L
)

summarizeSimulation(sim, trueHazardRatio = 0.75)
#> $power
#> [1] 0.81
#> $coverage
#> [1] 0.94
#> $bias
#> [1] 0.003
```

---

## Diagnostics

```r
# Covariate balance before and after IPCW
smd_before <- computeStandardizedDifferences(ccwData)
smd_after  <- computeStandardizedDifferences(ccwData, weights = weighted)
plotLovePlot(smd_before, smd_after)

# Weighted Kaplan–Meier curves
plotKaplanMeier(weighted)

# IPCW weight distribution
plotWeightDistribution(weighted)
```

---

## Vignettes

| Vignette | Description |
|---|---|
| [Single Study](vignettes/SingleStudy.Rmd) | End-to-end walkthrough: JAK inhibitor + recombinant zoster vaccine |
| [Multiple Analyses](vignettes/MultipleAnalyses.Rmd) | Sensitivity analyses with `runCcwAnalysesOnData()` |
| [Simulation Studies](vignettes/Simulation.Rmd) | Power analysis and pipeline validation with `simulateCcwData()` |

---

## HADES

CloneCensorWeight is part of [HADES](https://ohdsi.github.io/Hades/) — the
OHDSI Health Analytics Data-to-Evidence Suite.

<img src="https://ohdsi.github.io/Hades/images/hadesBadge.png" width="200" alt="HADES"/>

---

## License

Apache License 2.0 — see [LICENSE](LICENSE).

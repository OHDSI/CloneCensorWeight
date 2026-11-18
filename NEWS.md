# CloneCensorWeight 0.1.0

* Initial development

Changes that alter results of existing code, or of settings saved before them (which load with the
new defaults):

* `createFitCensoringModelsArgs(maxCohortSizeForFitting = 250000)` now also applies to the
  initiation hazard models (pooled logistic and Cox), not only to the binary models: a hazard model
  whose risk set holds more persons is fit on a uniform sample of that many persons (all their rows),
  drawn with the seed of `control` without changing the caller's random number stream, and then
  predicts the linear predictor of every person. The numerator and denominator models, the
  cross-validation of the prior and the high-correlation check use the sample; a uniform sample of
  persons needs no correction of the intercept. A sample with no initiations is not used (the model
  is fit on the whole risk set). Each bootstrap replicate draws its own sample. The
  diagnostics of every censoring model record `sampled`, `sampleSize` and `sampleEvents`, and the
  fitting message says when a hazard model is sampled. `maxCohortSizeForFitting = 0` fits the hazard
  models to their whole risk set, as before; risk sets of at most 250,000 persons are not affected.
  `computeCcwDiagnostics()` and `computeCcwSensitivity()` refit the hazard models given to them when
  `maxCohortSizeForFitting` differs.
* The Cox and pooled logistic outcome models are fitted with unstabilized weights by default
  (`createFitOutcomeModelArgs(msmWeights = "unstabilized")`), whatever `createComputeWeightsArgs()`
  `stabilized` is. `msmWeights = "stabilized"` reproduces the earlier fits.
* The truncation rule also caps the weights on the other scale (the unstabilized ones with the default
  `stabilized = TRUE`), on their own scale (`createComputeWeightsArgs(truncateMsmWeights = TRUE)`), so
  the weights of the Cox and pooled logistic models are capped whatever `msmWeights` is.
  `truncateMsmWeights = FALSE` forms them from the truncated weights, as before. The segments gain
  `otherScaleWeight`, and `fitOutcomeModel()` reports the cap, largest weight and share truncated of
  the weights its MSM used (`metaData$msmWeightSummary`).
* The robust variance of the Cox and pooled logistic outcome models comes from Cyclops'
  `vcov(fit, type = "robust")` by default (`createFitOutcomeModelArgs(robustVarianceEngine =
  "cyclops")`): the Cox model is fitted with the persons as Cyclops subjects, and the pooled logistic
  models are refitted, stacked, to rows per person. `robustVarianceEngine = "ccw"` gives this
  package's own sandwich, as before. Both are the HC0 sandwich clustered on person and agree to
  rounding when the outcome model is converged, as with the default `control` (about 1e-12 in the
  tests), and with `survival::coxph(cluster = )` and `sandwich::vcovCL(type = "HC0", cadjust =
  FALSE)` to 1e-6. With a loose convergence criterion the refit moves the estimates, and the two
  differ accordingly. For the pooled logistic models `"cyclops"` takes more time and memory than
  `"ccw"`, because of the refit. The standard errors of the Kaplan-Meier and Aalen-Johansen estimates
  stay this package's own.
* The default robust variance needs a Cyclops whose `vcov()` has `type = "robust"`. No released
  Cyclops has it; it comes with Cyclops 4.0.2 (on Cyclops develop), so the DESCRIPTION requires
  `Cyclops (>= 4.0.2)`, with `Remotes: ohdsi/Cyclops@develop` (plan decision 9 item 2). Every Cyclops
  4.0.2 on develop has it, so `fitOutcomeModel()` and `runCcwAnalyses()` no longer check the loaded
  Cyclops for it.
* The minimum detectable relative risk of logistic outcome models (`computeMdrr(modelType =
  "logistic")`, the `mdrr` of `runCcwAnalyses()` with `modelType = "logistic"`) comes from CohortMethod,
  whose left-censored branch now includes the treatment allocation in the standard error of the log
  odds ratio. It is larger than before (it had about 30% power instead of 80% with equal allocation),
  so the MDRR diagnostic can fail where it passed.
* The empirical calibration is seeded: `createCcwAnalysesSpecifications(calibrationSeed = 1)` fits the
  null distribution of each calibration group (`EmpiricalCalibration::fitMcmcNull()`, an MCMC) with
  this seed and R's default random number kinds, whatever the caller's `RNGkind()`, and then restores
  the caller's random number stream, so the calibrated p-values and confidence intervals and the EASE
  are the same from run to run and for any number of calibration threads; unseeded, they moved between
  runs. `calibrationSeed = NULL` runs unseeded, as before, and specifications saved without it load
  with `NULL`. `runCcwAnalyses()` treats a change of `calibrationSeed` (including an output folder
  whose specifications were saved without it) as a change of the summary only: it recomputes the
  results and diagnostics summaries and keeps the other outputs.
* The initiation hazard model has covariate-by-day-0 interactions by default
  (`createFitCensoringModelsArgs(covariateTimeInteraction = "day0")`; `"none"` with
  `hazardModelType = "cox"`), also with longer intervals, where they multiply the first interval.
  `covariateTimeInteraction = "none"` reproduces the earlier models.
* After the grace period, the weights of `"notInitiateWithin"` strategies are `1 / S(G | X)` from the
  initiation hazard model by default (`postGraceModel = "hazard"`). `postGraceModel = "binary"` uses the
  binary model at the end of the grace period, which is biased when starting the treatment changes the
  risk during the grace period.
* Settings saved without a field that has since been added (for example `truncateMsmWeights`,
  `postGraceModel`, `covariateTimeInteraction`, `msmWeights` or `robustVarianceEngine`) load with its
  current default, as listed above, and loading them prints a message naming the fields filled in this
  way. To reproduce an earlier analysis, set those fields to the earlier behaviour.
* The grace-period selection check (`computeGracePeriodSelection()`) compares initiated and
  not-initiated person-time on the same day of the grace period by default (`stratifyByDay = TRUE`), so
  a baseline risk that changes over the grace period does not raise the flag. `stratifyByDay = FALSE`
  compares them within groups of the predicted probability of starting only.
* A fixed truncation cap is applied as the same number on both weight scales with
  `truncateMsmWeights = TRUE`; `createWeightTruncation(otherScaleUpper = )` sets the cap on the other
  scale, and `computeCensoringWeights()` warns when the same cap truncates a much larger share on the
  other scale.
* Percentile truncation is taken over the clones whose weight varies
  (`createWeightTruncation(basis = "varying")`). `basis = "allClones"` reproduces the earlier caps.
  `metaData$weights$caps` reports the cap and the share of clones truncated per strategy.
* The covariate analyses of `simulateCloneCensorWeightData()` have windows that end on day -1
  (`endDay` in `analysisRef`), not day 0.
* The cross-validation of the pooled logistic initiation hazard model (`hazardModelType = "logistic"`)
  keeps all interval rows of a person in one fold. It split them over folds, so held-out rows were
  scored by a model fitted to other rows of the same person, which favours weaker penalties. Its
  selected prior variance, coefficients and the weights change. In bootstrap replicates with
  `fixHyperparameters = FALSE`, every censoring model (binary, Cox and pooled logistic) also keeps the
  copies of a person in one fold; the binary models pass them to `CohortMethod::createPs()` as
  `population$subjectId`. Cyclops' minimum for cross-validation (`createControl(minCVData = 100)`) then
  counts persons, not rows, so small replicates can fail where they did not. This needs Cyclops 4.0.1,
  which keeps the rows of a subject together for every model type.
* New arguments were inserted before existing ones: `covariateTimeInteraction`, `postGraceModel`,
  `initiateWithinModel` and `graceRegime` after `timeBasis` in `createFitCensoringModelsArgs()`, and
  `msmWeights` after `armByTimeInteraction` and `robustVarianceEngine` after `varianceEstimator` in
  `createFitOutcomeModelArgs()`. Calls that pass the later arguments by position must name them.
* With the indicator time basis (the default), the numerator of the pooled logistic initiation hazard
  model is computed in closed form, as the observed hazard of each interval, instead of fitted by
  Cyclops. It is the same model, so its coefficients and the stabilized weights change only within
  Cyclops' convergence tolerance, and its fit can no longer fail.
* The default control of `createFitCensoringModelsArgs()` cross-validates the prior variance with 3
  repetitions of the folds (`cvRepetitions = 3`) instead of 10, so cross-validation takes about a third
  of the time. It applies to the binary models (through `CohortMethod::createPs()`) and to the hazard
  models. The selected prior variances, the coefficients and the weights change. A control with
  `cvRepetitions = 10` reproduces the earlier fits.
* `createTreatmentStrategy()` stops when `graceStart` is not 0 for a `"notInitiateWithin"` or
  `"neverInitiate"` strategy (it applies to `"initiateWithin"` only).
* `computeCensoringModelCalibration()` and `computePositivity()` gain a `graceRegime` column, and the
  `variant` labels of `computeCcwDiagnostics()` name the `graceRegime` when a strategy is
  `"initiateWithin"`. The diagnostics of analyses with an `"initiateWithin"` strategy compare twice as
  many choices (both regimes), and the default variants of `computeCcwSensitivity()` add
  `graceRegime`.
* The pooled logistic outcome model with `armByTimeInteraction = TRUE` reports the first-interval odds
  ratio (`"or"`) as `NA` when either strategy has no outcomes, or only outcomes, in the first interval.
  Those cells have a fixed hazard and are left out of the fit, and the `target` coefficient reported as
  the odds ratio was then the odds ratio of a later interval, chosen by the column pivoting, with a
  finite standard error. The risks are unchanged. In `runCcwAnalyses()` such outcomes no longer enter
  the calibration of `or`.
* The pooled logistic outcome model warns when a strategy censors clones at a start on days that end no
  interval of its grid (`intervalLength > 1` with a `"notInitiateWithin"` or `"neverInitiate"` strategy,
  or an `"initiateWithin"` strategy with `graceStart > 0`). A clone censored inside an interval counted
  as at risk for the whole interval, while the clones that remain carry the weights that make up for
  it, which biases the risks and the first-interval odds ratio. Daily intervals reproduce the weighted
  Kaplan-Meier risks; the warning about a grid that ends no interval at the end of a grace period no
  longer says that any grid with such an interval does.
* `bootstrapCcwAnalysis()` analyzes each replicate with frequency weights by default
  (`createBootstrapArgs(method = "reweight")`, plan section 4.4) instead of copying the drawn persons
  into new data: the study population and clones of the full data are built once, the censoring models
  are refit with each person's count as its case weight (persons drawn 0 times are left out), on the
  full data's Cyclops data where the fit's structure allows (`Cyclops::updateWeights()`), and the
  outcome model is fit with each clone's weight times its person's count. Every step counts a person as
  often as its copies would (covariate tidying, the samples of `maxCohortSizeForFitting`, the percentile
  caps of the weights), so with `fixHyperparameters = TRUE` and cold starts the replicates equal those
  of the copy bootstrap to the convergence tolerance of the fits (1e-6 or closer in the tests, for
  pooled logistic and Cox censoring models, binary models, and Kaplan-Meier, Cox and pooled logistic
  outcome models). By default the refits start from the full-data coefficients
  (`createBootstrapArgs(warmStart = TRUE)`), which moves the replicates within that tolerance.
  `method = "copy"` gives the earlier bootstrap. With `fixHyperparameters = FALSE` the two methods
  differ, because cross-validation on weighted persons draws other folds than on copies. In the
  reweight bootstrap the binary models' sampling (`maxCohortSizeForFitting`) is seeded per replicate.

New features:

* `bootstrapCcwAnalysis(threads = )` runs the replicates on parallel workers (and
  `createCcwMultiThreadingSettings(bootstrapThreads = )` in `runCcwAnalyses()`, default 1). Each
  replicate depends only on its counts, so the results do not depend on the number of workers (up to
  the last bits of DuckDB's sums, which also vary from run to run). `metaData$bootstrap$replicateStatus`
  gives each replicate's status, the reason it failed, and the statuses of its censoring model fits
  that did not succeed (`censoringStatus`); failed replicates are counted and reported, as before.

* `createFitCensoringModelsArgs(pooledLogisticEngine = "plr")` fits the pooled logistic initiation
  hazard model with Cyclops's pooled logistic regression (`modelType = "plr"`), which keeps the
  covariates one row per person instead of expanding them to one row per person and interval. It needs
  a Cyclops with plr (not yet released) and stops with a clear error otherwise. With a fixed prior
  variance, or no prior when the maximum likelihood estimate exists, it gives the coefficients and
  weights of the default, `"lr"` (to about 1e-6); with cross-validation, a plr Cyclops before the
  separation and starting-variance fixes (`2b211c6a`, now on Cyclops develop) can choose a different
  variance. `"lr"` stays the default while plr is not in a released Cyclops.
  The diagnostics of each pooled logistic model record the `engine`.
* `createFitCensoringModelsArgs(covariateTimeInteraction = "day0" or "interval")` adds
  covariate-by-time interactions to the pooled logistic initiation hazard model.
* `computeCensoringModelCalibration()` and `plotCensoringModelCalibration()` compare the predicted
  and observed probabilities of starting by the end of the grace period.
* `computeWeightDistribution()` reports the exact weight of every day it reports.
* Diagnostics to choose these options for the data at hand: `computeCcwDiagnostics()` (all of them,
  for each `postGraceModel` and `covariateTimeInteraction`), `computeGracePeriodSelection()` (does
  starting change the grace-period risk?), `computePositivity()`, `computeWeightSummary()` (both scales,
  before and after truncation), and the opt-in `computeCcwSensitivity()` (estimates under alternative
  choices), with print methods and `plotGracePeriodSelection()`, `plotPositivity()`,
  `plotWeightSummary()` and `plotCcwSensitivity()`. `computeCensoringModelCalibration()` adds the
  calibration in the large, covariate strata and the selection check.
* Warnings for pooled logistic outcome grids that do not end an interval at the end of a grace period
  or use a smooth time basis with a grace period, for covariate windows that end on or after day 0,
  and for bootstrap replicates that sit predominantly on one side of the full-data estimate.
* `createFitCensoringModelsArgs(initiateWithinModel = "hazard")` gives the weights of
  `"initiateWithin(G)"` clones after the grace period as `1 / (1 - S(G | X))` from the initiation
  hazard model (shared with a `"notInitiateWithin(G)"` strategy of the same treatment), stabilized by
  `1 - S-bar(G)`. The default, `"binary"`, keeps the binary model at `G`. The calibration, positivity
  and grace-period selection diagnostics, `computeCcwDiagnostics()` and the default variants of
  `computeCcwSensitivity()` include both.
* `computeCcwDiagnostics(censoringModels = )` reuses fitted censoring models, and the diagnostics and
  `computeCcwSensitivity()` fit each censoring model once and share it between the choices that use it
  (before, the specified models were refit, and the binary model at `G` once per interaction setting).
  Choices that give the same censoring models are compared once.
* `fitCensoringModels()` says when the columns of `covariateTimeInteraction` are left out because they
  take a single value over the intervals in the fit (for example `"day0"` with 7-day intervals and a
  grace period ending on day 7 when nobody starts on day 7). For `"day0"` it says why: nobody (or
  everybody) starts in the first interval, so there is no first-interval effect to model, or the first
  interval is the only one in the fit, which shorter intervals can split.
* `runCcwAnalyses()` runs every analysis of a `createCcwAnalysesSpecifications()` against every
  target-treatment-outcome combination, in stages modelled on CohortMethod's `runCmAnalyses()`: data,
  study populations, censoring models, covariate balance, outcome models and a summary. Files are
  shared between analyses with the same settings for a stage, a censoring model shared by strategies
  of different analyses is fit once, and a run resumes from the files that exist. With
  `refitCensoringModelsForEveryOutcome = FALSE` (the default) the censoring models are fit once, on
  the study population without an outcome, and shared by every outcome; `TRUE` fits them on each
  outcome's clones. When the specifications in the output folder have changed, an interactive session
  asks whether to delete the old outputs, and other sessions delete them with a message.
  `createCcwMultiThreadingSettings()` and `createDefaultCcwMultiThreadingSettings()` set the threads
  of each stage. `createCcwAnalysesSpecifications()` stops when target-treatment-outcomes with the
  same cohorts have different covariate concept IDs or repeat an outcome, and `runCcwAnalyses()` when
  `analysesToExclude` leaves nothing to run.
* `getResultsSummary()` (one row per analysis, target-treatment-outcome and estimand, with the ratio
  estimands (`hr`, the first-interval `or`, and `rr_<h>`) empirically calibrated with the negative controls, and the confidence intervals of the
  hazard ratio also with at least 5 positive controls, whose true effect sizes are on the hazard
  scale), `getDiagnosticsSummary()` (MDRR, covariate balance, EASE, effective sample size and maximum
  weight (for every outcome, so they also decide which negative controls calibrate) and positivity against
  `createCcwDiagnosticThresholds()`, with `unblind` and `unblindForEvidenceSynthesis` as in
  CohortMethod) and `getFileReference()` read the outputs of `runCcwAnalyses()`.
* `computeCensoringWeights()` and the diagnostics accept censoring models fit on clones without an
  outcome (`metaData$outcomeId` `NA`), as `runCcwAnalyses()` shares them across outcomes.
  `bootstrapCcwAnalysis()` given such models refits them without an outcome in every replicate, so
  the replicates repeat the estimator of the full-data estimate.
* `computeCensoringWeights()` stops when the clones are followed past the horizon of a
  `"neverInitiate"` censoring model (models fit on a shorter risk window), whose weights would
  otherwise stay at their value on the horizon while clones are still censored at initiation.
* The results data model, modelled on CohortMethod's: `exportToCsv()` writes the results of
  `runCcwAnalyses()` to CSV files and `Results_<databaseId>.zip`, with counts below `minCellCount`
  written as `-minCellCount`; `getResultsDataModelSpecifications()` describes the tables, all
  prefixed `ccw_`, with `database_id` in the primary key of every site-specific table and strategy
  IDs in the keys of the tables per strategy (`strategy_id`) and per contrast (`target_strategy_id`,
  `comparator_strategy_id`); `createResultsDataModel()`, `uploadResults()`, `migrateDataModel()` and
  `getDataMigrator()` create, fill and migrate the tables with ResultModelManager. The DDL of this
  first version is frozen; later changes go through migrations. `inst/Results_Eunomia.zip` holds an
  export of a Eunomia run for tests.
* `findResultsZipFile()` returns the one `Results_<databaseId>.zip` in an export folder, and stops,
  naming them, when there is more than one (as after an earlier run with another `databaseId`). The
  Strategus module uploads the zip it returns, rather than the first zip file of its results folder.
* `convertUntypedListToCcwAnalysesSpecifications()` documents the two forms in which a Strategus
  module gets the specifications back, and stops with a clear message when given the
  `CcwAnalysesSpecifications` object instead of its `$toList()`, or a list without
  `ccwAnalysisList` or `targetTreatmentOutcomesList`.
* `getResultsSummary()` gains the weighted days of follow-up and weighted outcomes of each strategy
  (`targetWeightedDays`, `targetWeightedOutcomes` and the comparator's), and the outcome model files
  of `runCcwAnalyses()` keep the attrition and the follow-up distribution of the clones of outcomes
  of interest.
* `"initiateWithin"` strategies with `graceStart > 0` (plan section 2.2): a start before day
  `graceStart` is a deviation. `createClones()`, `fitCensoringModels()`, `computeCensoringWeights()` and
  `computeTrueRisks()` support them. The weights of the days before `graceStart` come from the
  initiation hazard model through `graceStart - 1` (with `initiateWithinModel = "binary"`; shared with a
  `"notInitiateWithin(graceStart - 1)"` strategy) or through `G` (with `"hazard"`), and the binary
  model at `G` is fit on the persons who have not started before `graceStart`.
* `createFitCensoringModelsArgs(graceRegime = "naturalGrace")` gives `"initiateWithin"` strategies the
  natural-grace (last-interval) weights of Cain et al. (2010) and Wanis et al. (2024): 1 for clones that
  started before `G`, and the inverse probability of starting on day `G` of the persons who reach it
  untreated for those that start on day `G` (binary model fit on them, or the day-`G` initiation hazard
  with `initiateWithinModel = "hazard"`). The default, `"binaryAtG"`, keeps the binary-at-`G` weights
  (plan decision 9 item 7). `computeTrueRisks(graceRegime = )` gives the truth of either regime. The
  calibration and positivity diagnostics, `computeCcwDiagnostics()` and the default variants of
  `computeCcwSensitivity()` include both regimes (leaving out a regime or `initiateWithinModel` that
  the hazard model's `intervalLength` does not allow). `fitCensoringModels()` warns when, under
  `"naturalGrace"`, persons reach day `G` untreated but nobody starts on day `G`, so that the regime is
  not identified.

Dependencies:

* Requires Andromeda 1.2.2 (`Remotes: ohdsi/Andromeda@develop`), whose `flushAndromeda()` no longer
  creates a temporary database on every flush. The test-only workaround for Andromeda 1.2.0 and 1.2.1
  (`tests/testthat/setup-andromeda.R`) is removed.
* Requires CohortMethod 6.0.3 from its left-censored branch, and uses its exports instead of copies:
  the settings objects inherit from `CohortMethod::AbstractSerializableSettings` (keeping this
  package's message about fields filled with defaults and its 15-digit JSON), the covariate settings
  are serialized with `CohortMethod::.convertAttrToMember()`, `.convertMemberToAttr()` and
  `.convertEmptyListToLogical()`, and `exportToCsv()` censors small counts with
  `CohortMethod::enforceMinCellValue()`. The exported results are unchanged.

Documentation:

* New vignettes `SingleStudy` (one target trial on Synthea: first AMI, metoprolol within 7 days vs not;
  every step of the pipeline with its diagnostics, and a comparison with naive and landmark analyses),
  `MultipleAnalyses` (two grace periods x three outcome models x the outcome and 18 negative controls
  through `runCcwAnalyses()`, empirical calibration, `exportToCsv()`, and `createResultsDataModel()` and
  `uploadResults()` into SQLite) and `ResultsSchema` (the results data model, generated from
  `inst/csv/resultsDataModelSpecification.csv`). The database-dependent code is not evaluated when the
  vignettes are built; the outputs shown are precomputed from runs on Synthea and stored in
  `vignettes/data/`, with counts below 11 suppressed. They replace the earlier `SingleStudy` and
  `MultipleAnalyses` vignettes, which no longer matched the package.
* The `PackageComparison` vignette (plan sections 6.4 and 7.3) shows precomputed end-to-end comparisons:
  the What If G = 0 reduction (Program 17.3 on NHEFS; gated at 1e-6 on coefficients and risks), the
  simulated DGP (a) against `computeTrueRisks()` (gated; the g-formula packages were not installed, so
  their half is skipped), `ccwr` on its lung cancer data (outcome stage gated, full pipeline reported),
  and Synthea AMI against `ccwr` (reported). Comparisons whose package or data were missing are listed
  as skipped. The scripts that make the results are kept with the developers' validation materials.
* The cross-check `test-crossChecks-whatif.R` repeats the What If G = 0 reduction in the cross-checks,
  which run on demand (`CCW_RUN_CROSSCHECKS=true`, as the developers' `extras/nightly/runNightly.sh`
  sets; nothing schedules them), when the environment variable `CCW_NHEFS_PATH` names a local copy of
  the NHEFS data, which is not distributed with the package; it skips otherwise.

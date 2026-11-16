# CloneCensorWeight 0.1.0

* Initial development

Changes that alter results of existing code, or of settings saved before them (which load with the
new defaults):

* The Cox and pooled logistic outcome models are fitted with unstabilized weights by default
  (`createFitOutcomeModelArgs(msmWeights = "unstabilized")`), whatever `createComputeWeightsArgs()`
  `stabilized` is. `msmWeights = "stabilized"` reproduces the earlier fits.
* The truncation rule also caps the weights on the other scale (the unstabilized ones with the default
  `stabilized = TRUE`), on their own scale (`createComputeWeightsArgs(truncateMsmWeights = TRUE)`), so
  the weights of the Cox and pooled logistic models are capped whatever `msmWeights` is.
  `truncateMsmWeights = FALSE` forms them from the truncated weights, as before. The segments gain
  `otherScaleWeight`, and `fitOutcomeModel()` reports the cap, largest weight and share truncated of
  the weights its MSM used (`metaData$msmWeightSummary`).
* The initiation hazard model has covariate-by-day-0 interactions by default
  (`createFitCensoringModelsArgs(covariateTimeInteraction = "day0")`; `"none"` with
  `hazardModelType = "cox"`), also with longer intervals, where they multiply the first interval.
  `covariateTimeInteraction = "none"` reproduces the earlier models.
* After the grace period, the weights of `"notInitiateWithin"` strategies are `1 / S(G | X)` from the
  initiation hazard model by default (`postGraceModel = "hazard"`). `postGraceModel = "binary"` uses the
  binary model at the end of the grace period, which is biased when starting the treatment changes the
  risk during the grace period.
* Settings saved without a field that has since been added (for example `truncateMsmWeights`,
  `postGraceModel`, `covariateTimeInteraction` or `msmWeights`) load with its current default, as listed
  above, and loading them prints a message naming the fields filled in this way. To reproduce an
  earlier analysis, set those fields to the earlier behaviour.
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
* New arguments were inserted before existing ones: `covariateTimeInteraction` and `postGraceModel`
  after `timeBasis` in `createFitCensoringModelsArgs()`, and `msmWeights` after
  `armByTimeInteraction` in `createFitOutcomeModelArgs()`. Calls that pass the later arguments by
  position must name them.

New features:

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

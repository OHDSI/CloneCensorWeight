# CloneCensorWeight 0.1.0

* Initial development

Changes that alter results of existing code, or of settings saved before them (which load with the
new defaults):

* The Cox and pooled logistic outcome models are fitted with unstabilized weights by default
  (`createFitOutcomeModelArgs(msmWeights = "unstabilized")`), whatever `createComputeWeightsArgs()`
  `stabilized` is. `msmWeights = "stabilized"` reproduces the earlier fits. The truncation cap applies
  to the weights of `computeCensoringWeights()`, so with the default `stabilized = TRUE` it does not
  bound the unstabilized weights of these models.
* After the grace period, the weights of `"notInitiateWithin"` strategies come from the binary model at
  the end of the grace period (`createFitCensoringModelsArgs(postGraceModel = "binary")`).
  `postGraceModel = "hazard"` reproduces the earlier weights.
* Percentile truncation is taken over the clones whose weight varies
  (`createWeightTruncation(basis = "varying")`). `basis = "allClones"` reproduces the earlier caps.
  `metaData$weights$caps` reports the cap and the share of clones truncated per strategy.
* The covariate analyses of `simulateCloneCensorWeightData()` have windows that end on day -1
  (`endDay` in `analysisRef`), not day 0.
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
* Warnings for pooled logistic outcome grids that do not end an interval at the end of a grace period
  or use a smooth time basis with a grace period, for covariate windows that end on or after day 0,
  and for bootstrap replicates that sit predominantly on one side of the full-data estimate.

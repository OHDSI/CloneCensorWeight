# The censoring models (without the baseline covariates, to keep the example fast; by default they
# use all covariates) and the weighted clones, also split on the daily grid of the pooled logistic
# outcome model:
censoringModels <- fitCensoringModels(
  clones = clones,
  cloneCensorWeightData = data,
  fitCensoringModelsArgs = createFitCensoringModelsArgs(useCovariates = FALSE)
)
weightedClones <- computeCensoringWeights(
  clones = clones,
  censoringModels = censoringModels,
  computeWeightsArgs = createComputeWeightsArgs(),
  outcomeTimeBasis = createTimeBasis(type = "naturalSpline", intervalLength = 1, df = 4)
)

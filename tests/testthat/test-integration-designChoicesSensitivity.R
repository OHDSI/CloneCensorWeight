library(testthat)

# computeCcwSensitivity() (contract section 12.5), split from test-integration-designChoices.R so that the
# files run in parallel. The helpers are in helper-designChoices.R.

test_that("The sensitivity analysis refits every default variant and custom ones", {
  data <- simulateCloneCensorWeightData(n = 300, maxDays = 60, graceEnd = 7, seed = 7)
  clones <- designClones(data)
  on.exit({
    close(clones)
    close(data)
  })
  fitCensoringModelsArgs <- designArgs()
  fitOutcomeModelArgs <- createFitOutcomeModelArgs(horizons = 30)
  sensitivity <- suppressMessages(computeCcwSensitivity(clones, data, fitCensoringModelsArgs, fitOutcomeModelArgs = fitOutcomeModelArgs))
  expect_s3_class(sensitivity, "CcwSensitivity")
  variants <- attr(sensitivity, "variants")
  expect_equal(variants$variant, c(
    "specified", "postGraceModel = binary", "initiateWithinModel = hazard", "graceRegime = naturalGrace", "covariateTimeInteraction = none",
    "covariateTimeInteraction = interval", "truncateMsmWeights = FALSE", "basis = allClones", "msmWeights = stabilized"
  ))
  expect_equal(variants$postGraceModel, c("hazard", "binary", rep("hazard", 7)))
  expect_equal(variants$initiateWithinModel, c("binary", "binary", "hazard", rep("binary", 6)))
  expect_equal(variants$graceRegime, c(rep("binaryAtG", 3), "naturalGrace", rep("binaryAtG", 5)))
  expect_equal(variants$covariateTimeInteraction, c("day0", "day0", "day0", "day0", "none", "interval", "day0", "day0", "day0"))
  expect_equal(variants$truncateMsmWeights, c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE))
  expect_equal(variants$truncationBasis, c(rep("varying", 7), "allClones", "varying"))
  expect_equal(variants$msmWeights, c(rep("unstabilized", 8), "stabilized"))
  # Each variant equals the pipeline run directly with its settings:
  direct <- function(fitCensoringModelsArgs, computeWeightsArgs, fitOutcomeModelArgs) {
    models <- suppressMessages(fitCensoringModels(clones, data, fitCensoringModelsArgs))
    weighted <- suppressMessages(computeCensoringWeights(clones, models, computeWeightsArgs))
    on.exit(close(weighted))
    return(suppressMessages(fitOutcomeModel(weighted, fitOutcomeModelArgs))$estimates)
  }
  pick <- function(variant) {
    estimates <- as.data.frame(sensitivity[sensitivity$variant == variant, colnames(sensitivity)[-(1:2)]])
    attr(estimates, "variants") <- NULL
    return(estimates)
  }
  directSettings <- list(
    "specified" = list(fitCensoringModelsArgs, createComputeWeightsArgs(), fitOutcomeModelArgs),
    "postGraceModel = binary" = list(designArgs(postGraceModel = "binary"), createComputeWeightsArgs(), fitOutcomeModelArgs),
    "initiateWithinModel = hazard" = list(designArgs(initiateWithinModel = "hazard"), createComputeWeightsArgs(), fitOutcomeModelArgs),
    "graceRegime = naturalGrace" = list(designArgs(graceRegime = "naturalGrace"), createComputeWeightsArgs(), fitOutcomeModelArgs),
    "covariateTimeInteraction = none" = list(designArgs(covariateTimeInteraction = "none"), createComputeWeightsArgs(), fitOutcomeModelArgs),
    "truncateMsmWeights = FALSE" = list(fitCensoringModelsArgs, createComputeWeightsArgs(truncateMsmWeights = FALSE), fitOutcomeModelArgs),
    "msmWeights = stabilized" = list(fitCensoringModelsArgs, createComputeWeightsArgs(), createFitOutcomeModelArgs(horizons = 30, msmWeights = "stabilized"))
  )
  # The default tier compares one variant of each kind of settings (censoring models, weights, outcome
  # model) besides the specified one:
  compared <- testConfigurations(names(directSettings), c("specified", "initiateWithinModel = hazard", "truncateMsmWeights = FALSE", "msmWeights = stabilized"))
  for (variant in compared) {
    expect_equal(pick(variant), as.data.frame(do.call(direct, directSettings[[variant]])), info = variant)
  }
  # The KM risks do not depend on msmWeights, the hazard ratio does:
  specified <- pick("specified")
  stabilized <- pick("msmWeights = stabilized")
  expect_equal(specified$estimate[-1], stabilized$estimate[-1])
  expect_output(print(sensitivity), "specified")

  custom <- suppressMessages(computeCcwSensitivity(
    clones, data, fitCensoringModelsArgs,
    fitOutcomeModelArgs = createFitOutcomeModelArgs(modelType = "km", horizons = 30),
    variants = list("binary, none, unstabilized" = list(
      fitCensoringModelsArgs = list(postGraceModel = "binary", covariateTimeInteraction = "none"),
      computeWeightsArgs = list(stabilized = FALSE)
    ))
  ))
  expect_equal(unique(custom$variant), c("specified", "binary, none, unstabilized"))
  expect_equal(unique(custom$changes), c("", "postGraceModel = binary; covariateTimeInteraction = none; stabilized = FALSE"))
  # Invalid variants stop. computeCcwSensitivity() runs the specified analysis before it applies the
  # variant, so the default tier checks the function that applies it, and the full tier the whole call:
  if (isFullTestRun()) {
    expect_error(computeCcwSensitivity(clones, data, fitCensoringModelsArgs, variants = list(a = list(fitCensoringModelsArgs = list(nonsense = 1)))), "no setting 'nonsense'")
    expect_error(computeCcwSensitivity(clones, data, fitCensoringModelsArgs, variants = list(a = list(other = list()))), "unknown settings")
    expect_error(suppressMessages(computeCcwSensitivity(clones, data, fitCensoringModelsArgs, variants = list(a = list(fitCensoringModelsArgs = list(postGraceModel = "x"))))), "postGraceModel")
  } else {
    base <- list(fitCensoringModelsArgs = fitCensoringModelsArgs, computeWeightsArgs = createComputeWeightsArgs(), fitOutcomeModelArgs = createFitOutcomeModelArgs())
    expect_error(.applyVariant(base, list(fitCensoringModelsArgs = list(nonsense = 1)), "a"), "no setting 'nonsense'")
    expect_error(.applyVariant(base, list(other = list()), "a"), "unknown settings")
    expect_error(.applyVariant(base, list(fitCensoringModelsArgs = list(postGraceModel = "x")), "a"), "postGraceModel")
  }
  skip_if_not_installed("ggplot2")
  expect_s3_class(plotCcwSensitivity(sensitivity), "ggplot")
  expect_s3_class(plotCcwSensitivity(sensitivity, estimands = "hr"), "ggplot")
})

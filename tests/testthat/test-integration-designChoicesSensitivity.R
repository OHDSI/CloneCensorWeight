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

test_that("Interaction settings that leave out every interaction column share one hazard model fit", {
  data <- simulateCloneCensorWeightData(n = 300, maxDays = 60, graceEnd = 7, seed = 7)
  clones <- designClones(data)
  on.exit({
    close(clones)
    close(data)
  })
  variants <- list(
    "covariateTimeInteraction = none" = list(fitCensoringModelsArgs = list(covariateTimeInteraction = "none")),
    "covariateTimeInteraction = interval" = list(fitCensoringModelsArgs = list(covariateTimeInteraction = "interval"))
  )
  fitOutcomeModelArgs <- createFitOutcomeModelArgs(modelType = "km", horizons = 30)
  run <- function(intervalLength, variants) {
    messages <- capture_messages(sensitivity <- computeCcwSensitivity(
      clones, data, designArgs(timeBasis = createTimeBasis(intervalLength = intervalLength)),
      fitOutcomeModelArgs = fitOutcomeModelArgs, variants = variants
    ))
    return(list(fits = sum(grepl("^Fitting logistic initiation hazard model", messages)), messages = messages, sensitivity = sensitivity))
  }
  estimates <- function(sensitivity, variant) {
    rows <- as.data.frame(sensitivity[sensitivity$variant == variant, colnames(sensitivity)[-(1:2)]])
    attr(rows, "variants") <- NULL
    return(rows)
  }
  # With 8-day intervals the grace period (days 0-7) is one interval, so "day0" and "interval" leave out
  # their interaction columns: the hazard model is the model without interactions, fit once (for the
  # specified "day0") and reused by both variants, which give the specified estimates:
  weekly <- run(8, variants)
  expect_equal(weekly$fits, 1)
  expect_true(any(grepl("with \"interval\" it has no interaction columns either", weekly$messages, fixed = TRUE)))
  # ... and share the specified variant's weights:
  expect_equal(sum(grepl("so it shares its weights", weekly$messages, fixed = TRUE)), 2)
  specified <- estimates(weekly$sensitivity, "specified")
  expect_identical(estimates(weekly$sensitivity, "covariateTimeInteraction = none"), specified)
  expect_identical(estimates(weekly$sensitivity, "covariateTimeInteraction = interval"), specified)
  # ... which are those of a separate fit without interactions:
  models <- suppressMessages(fitCensoringModels(clones, data, designArgs(timeBasis = createTimeBasis(intervalLength = 8), covariateTimeInteraction = "none")))
  weighted <- suppressMessages(computeCensoringWeights(clones, models))
  on.exit(close(weighted), add = TRUE)
  expect_equal(specified, as.data.frame(suppressMessages(fitOutcomeModel(weighted, fitOutcomeModelArgs))$estimates))
  # With 4-day intervals "day0" and "interval" keep interaction columns (of the first and the second
  # interval), so each has its own fit, also when the model without interactions was fit before:
  expect_equal(run(4, variants)$fits, 3)
})

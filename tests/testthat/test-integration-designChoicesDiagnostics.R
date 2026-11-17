library(testthat)

# computeCcwDiagnostics() (contract section 12.5), split from test-integration-designChoices.R so that the
# files run in parallel. The helpers are in helper-designChoices.R.

# The data, clones and diagnostics (with the messages of computeCcwDiagnostics()) that the tests below
# share, computed once: each computeCcwDiagnostics() call weights the clones and computes the balance of
# 16 variants. The tests do not change or close them.
sharedDiagnostics <- local({
  shared <- NULL
  function() {
    if (is.null(shared)) {
      data <- simulateCloneCensorWeightData(n = 300, maxDays = 60, graceEnd = 7, seed = 8)
      clones <- designClones(data)
      messages <- capture_messages(diagnostics <- computeCcwDiagnostics(clones, data, designArgs(), numberOfGroups = 5, strataCovariateIds = 1001))
      shared <<- list(data = data, clones = clones, messages = messages, diagnostics = diagnostics)
    }
    return(shared)
  }
})

# Evaluates `expr` until it emits a message that contains `text`, and returns that message; NULL when
# `expr` finishes without one.
firstMessageWith <- function(expr, text) {
  return(tryCatch(
    withCallingHandlers(
      {
        expr
        NULL
      },
      message = function(m) {
        if (grepl(text, conditionMessage(m), fixed = TRUE)) {
          signalCondition(structure(class = c("firstMessageWith", "condition"), list(message = conditionMessage(m), call = NULL)))
        }
        invokeRestart("muffleMessage")
      }
    ),
    firstMessageWith = function(condition) conditionMessage(condition)
  ))
}

test_that("computeCcwDiagnostics compares postGraceModel, covariateTimeInteraction, initiateWithinModel and graceRegime choices", {
  shared <- sharedDiagnostics()
  data <- shared$data
  clones <- shared$clones
  # Four models are fit: the binary model at G, that of the persons untreated at G (graceRegime =
  # "naturalGrace"), and the hazard model through G with each interaction setting, shared by both
  # strategies with initiateWithinModel = "hazard" (under both regimes):
  messages <- shared$messages
  diagnostics <- shared$diagnostics
  expect_equal(sum(grepl("^Fitting", messages)), 4)
  expect_s3_class(diagnostics, "CcwDiagnostics")
  expect_equal(diagnostics$variants$postGraceModel, rep(c("hazard", "binary"), 8))
  expect_equal(diagnostics$variants$covariateTimeInteraction, rep(c("day0", "day0", "none", "none"), 4))
  expect_equal(diagnostics$variants$initiateWithinModel, rep(rep(c("binary", "hazard"), each = 4), 2))
  expect_equal(diagnostics$variants$graceRegime, rep(c("binaryAtG", "naturalGrace"), each = 8))
  expect_equal(diagnostics$variants$specified, c(TRUE, rep(FALSE, 15)))
  expect_equal(diagnostics$variants$variant[1], "postGraceModel = hazard, covariateTimeInteraction = day0, initiateWithinModel = binary, graceRegime = binaryAtG")
  expect_equal(diagnostics$balance$variant, diagnostics$variants$variant)
  expect_equal(colnames(diagnostics$balance), c("variant", "covariates", "maxAbsStdDiffBefore", "maxAbsStdDiffAfter", "aboveThresholdBefore", "aboveThresholdAfter"))
  expect_equal(sort(unique(diagnostics$weightSummary$variant)), sort(diagnostics$variants$variant))
  # The hazard model with each interaction setting and the binary model once per strategy:
  overall <- diagnostics$calibration[diagnostics$calibration$grouping == "overall", ]
  expect_equal(
    paste(overall$strategyId, overall$model, overall$covariateTimeInteraction, overall$graceRegime),
    c(
      "1 binary NA binaryAtG", "2 hazard day0 NA", "2 binary NA NA", "2 hazard none NA", "1 hazard day0 binaryAtG", "1 hazard none binaryAtG",
      "1 binary NA naturalGrace", "1 hazard day0 naturalGrace", "1 hazard none naturalGrace"
    )
  )
  expect_true(any(diagnostics$calibration$grouping == "covariate"))
  expect_s3_class(diagnostics$gracePeriodSelection, "CcwGracePeriodSelection")
  expect_s3_class(diagnostics$positivity, "CcwPositivity")
  # Positivity of both models of each strategy, for the specified interaction setting:
  expect_equal(
    paste(diagnostics$positivity$strategyId, diagnostics$positivity$model, diagnostics$positivity$graceRegime),
    c("1 binary binaryAtG", "1 hazard binaryAtG", "1 binary naturalGrace", "1 hazard naturalGrace", "2 hazard NA", "2 binary NA")
  )
  predictions <- attr(diagnostics$positivity, "predictions")
  expect_equal(nrow(unique(predictions[, c("strategyId", "model", "graceRegime")])), 6)
  expect_false(anyDuplicated(predictions[, c("strategyId", "model", "graceRegime", "rowId")]) > 0)
  expect_null(diagnostics$sensitivity)
  expect_output(print(diagnostics), "Clone-censor-weight diagnostics")

  # Each variant's models equal a direct fit with its settings:
  hazard <- suppressMessages(fitCensoringModels(clones, data, designArgs()))
  fits <- .censoringFitsFromModels(suppressMessages(fitCensoringModels(clones, data, designArgs(postGraceModel = "binary"))))
  messages <- capture_messages(fromFits <- .fitCensoringModelsWithFits(clones, data, designArgs(), fits = fits)$models)
  expect_false(any(grepl("^Fitting", messages)))
  # (The order of the interaction coefficients can differ between fits of the same model.)
  sortCoefficients <- function(models) {
    lapply(models, function(model) {
      model$coefficients <- model$coefficients[order(names(model$coefficients))]
      model
    })
  }
  expect_equal(sortCoefficients(fromFits), sortCoefficients(hazard))
  expect_equal(attr(fromFits, "metaData"), attr(hazard, "metaData"))
  # The specified variant's balance equals computeCovariateBalance() on its weights:
  weighted <- suppressMessages(computeCensoringWeights(clones, hazard))
  balance <- suppressMessages(computeCovariateBalance(weighted, data))
  close(weighted)
  expect_equal(diagnostics$balance$maxAbsStdDiffAfter[1], max(abs(balance$afterMatchingStdDiff), na.rm = TRUE))
  expect_equal(diagnostics$balance$aboveThresholdAfter[1], sum(abs(balance$afterMatchingStdDiff) > 0.1, na.rm = TRUE))
})

test_that("computeCcwDiagnostics with the Cox hazard model compares fewer choices, and adds the sensitivity analysis", {
  # Full tier: another computeCcwDiagnostics() call, with the sensitivity analysis. The default tier
  # runs sensitivity = TRUE on the smaller design of the last test, and the Cox grid with the Cox models
  # given below.
  skipUnlessFullTests()
  shared <- sharedDiagnostics()
  data <- shared$data
  clones <- shared$clones
  withSensitivity <- suppressMessages(computeCcwDiagnostics(
    clones, data, designArgs(hazardModelType = "cox"),
    numberOfGroups = 5, sensitivity = TRUE,
    fitOutcomeModelArgs = createFitOutcomeModelArgs(modelType = "km", horizons = 30)
  ))
  # The Cox model has no interactions, so only postGraceModel, initiateWithinModel and graceRegime vary:
  expect_equal(withSensitivity$variants$covariateTimeInteraction, rep("none", 8))
  expect_equal(withSensitivity$variants$initiateWithinModel, rep(c("binary", "binary", "hazard", "hazard"), 2))
  expect_s3_class(withSensitivity$sensitivity, "CcwSensitivity")
  expect_output(print(withSensitivity), "Estimates under the specified settings")
})

test_that("computeCcwDiagnostics reuses given censoring models and fits only the ones it lacks", {
  shared <- sharedDiagnostics()
  data <- shared$data
  clones <- shared$clones
  args <- designArgs()
  models <- suppressMessages(fitCensoringModels(clones, data, args))
  messages <- capture_messages(reused <- computeCcwDiagnostics(clones, data, args, numberOfGroups = 5, strataCovariateIds = 1001, censoringModels = models))
  # The specified variant's models (the binary model at G and the hazard model with "day0") are not
  # refit; the binary model at G also serves postGraceModel = "binary", and the hazard model through G
  # initiateWithinModel = "hazard" (under both regimes). Only the hazard model without interactions and
  # the binary model of graceRegime = "naturalGrace" are fit:
  fitting <- grep("^Fitting", messages, value = TRUE)
  expect_length(fitting, 2)
  expect_match(fitting[1], "logistic initiation hazard model through day 7 for strategy 2")
  expect_match(fitting[2], "binary censoring model at day 7 of the persons who have not started before day 7 for strategy 1")
  expect_true(any(grepl("^Reusing the fitted binary censoring model at day 7 for strategy 1", messages)))
  # The same diagnostics as without the models:
  fresh <- shared$diagnostics
  expect_equal(reused$variants, fresh$variants)
  expect_equal(reused$balance, fresh$balance)
  expect_equal(reused$calibration, fresh$calibration)
  expect_equal(reused$weightSummary, fresh$weightSummary)
  expect_equal(reused$positivity, fresh$positivity)
  expect_equal(reused$gracePeriodSelection, fresh$gracePeriodSelection)
  # Models fit with other settings are refit, with a message (which comes before any model is fit; the
  # full tier runs the call to the end):
  expect_match(
    firstMessageWith(computeCcwDiagnostics(clones, data, designArgs(useCovariates = FALSE), numberOfGroups = 5, censoringModels = models), "fit with other settings"),
    "fit with other settings than fitCensoringModelsArgs"
  )
  expect_error(computeCcwDiagnostics(clones, data, args, censoringModels = list()), "CcwCensoringModels")
})

test_that("computeCcwDiagnostics takes the specified settings from Cox models given, and compares fewer choices with them", {
  # The light form of the next test (and of the Cox grid of the test with the sensitivity analysis):
  # Cox hazard models give 8 variants rather than 16.
  shared <- sharedDiagnostics()
  data <- shared$data
  clones <- shared$clones
  cox <- suppressMessages(fitCensoringModels(clones, data, designArgs(hazardModelType = "cox", initiateWithinModel = "hazard")))
  # Without fitCensoringModelsArgs, the specified settings are those of the models:
  messages <- capture_messages(fromModels <- computeCcwDiagnostics(clones, data, numberOfGroups = 5, censoringModels = cox))
  variants <- fromModels$variants
  specified <- variants[variants$specified, ]
  expect_equal(nrow(specified), 1)
  expect_equal(c(specified$postGraceModel, specified$covariateTimeInteraction, specified$initiateWithinModel, specified$graceRegime), c("hazard", "none", "hazard", "binaryAtG"))
  expect_equal(variants$variant[1], specified$variant)
  # The Cox model has no interactions, so only postGraceModel, initiateWithinModel and graceRegime vary:
  expect_equal(variants$covariateTimeInteraction, rep("none", 8))
  expect_equal(nrow(unique(variants[, c("postGraceModel", "initiateWithinModel", "graceRegime")])), 8)
  # The hazard models are reused; only the binary models at G of strategy 1 (of both regimes) are fit:
  fitting <- grep("^Fitting", messages, value = TRUE)
  expect_length(fitting, 2)
  expect_match(fitting[1], "binary censoring model at day 7 for strategy 1")
  expect_match(fitting[2], "binary censoring model at day 7 of the persons who have not started before day 7 for strategy 1")
  expect_true(any(grepl("^Reusing the fitted cox initiation hazard model through day 7 for strategy 1", messages)))
})

test_that("computeCcwDiagnostics takes the specified settings from the models given, and refits only the models whose settings differ", {
  # Full tier: three more computeCcwDiagnostics() calls of 16 variants each (the test above is the
  # light form of the first).
  skipUnlessFullTests()
  shared <- sharedDiagnostics()
  data <- shared$data
  clones <- shared$clones
  fresh <- shared$diagnostics
  args <- designArgs()
  models <- suppressMessages(fitCensoringModels(clones, data, args))
  # Without fitCensoringModelsArgs, the specified settings are those of the models:
  other <- suppressMessages(fitCensoringModels(clones, data, designArgs(covariateTimeInteraction = "none", initiateWithinModel = "hazard")))
  messages <- capture_messages(fromModels <- computeCcwDiagnostics(clones, data, numberOfGroups = 5, censoringModels = other))
  specified <- fromModels$variants[fromModels$variants$specified, ]
  expect_equal(c(specified$postGraceModel, specified$covariateTimeInteraction, specified$initiateWithinModel), c("hazard", "none", "hazard"))
  expect_equal(fromModels$variants$variant[1], specified$variant)
  # The hazard model without interactions is reused; the binary models at G (of both regimes) and the
  # hazard model with "day0" are fit:
  expect_equal(sum(grepl("^Fitting", messages)), 3)
  # Models fit with other settings are refit, with a message:
  expect_message(
    computeCcwDiagnostics(clones, data, designArgs(useCovariates = FALSE), numberOfGroups = 5, censoringModels = models),
    "fit with other settings than fitCensoringModelsArgs"
  )
  # The binary model at G does not depend on the hazard model's settings, so with only timeBasis
  # different it is reused and the two hazard models are fit again:
  weekly <- suppressMessages(fitCensoringModels(clones, data, designArgs(timeBasis = createTimeBasis(intervalLength = 7))))
  messages <- capture_messages(fromWeekly <- computeCcwDiagnostics(clones, data, args, numberOfGroups = 5, censoringModels = weekly))
  expect_true(any(grepl("fit with other settings than fitCensoringModelsArgs", messages)))
  expect_true(any(grepl("^Reusing the fitted binary censoring model at day 7 for strategy 1", messages)))
  expect_equal(sum(grepl("^Fitting", messages)), 3)
  expect_false(any(grepl("^Fitting binary censoring model at day 7 for strategy", messages)))
  expect_equal(fromWeekly$weightSummary, fresh$weightSummary)
})

test_that("computeCcwDiagnostics compares choices that give the same models once", {
  data <- simulateCloneCensorWeightData(n = 300, maxDays = 60, graceEnd = 7, seed = 10)
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowStart = 0, riskWindowEnd = 60, endAnchor = "cohort start")
  ))
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
    createTreatmentStrategy(2, "Initiate within 14 days", "initiateWithin", graceEnd = 14)
  )
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = strategies)))
  on.exit({
    close(clones)
    close(data)
  })
  diagnostics <- suppressMessages(computeCcwDiagnostics(
    clones, data, designArgs(),
    numberOfGroups = 5, sensitivity = TRUE,
    fitOutcomeModelArgs = createFitOutcomeModelArgs(modelType = "km", horizons = 30)
  ))
  expect_s3_class(diagnostics$sensitivity, "CcwSensitivity")
  expect_output(print(diagnostics), "Estimates under the specified settings")
  # With initiateWithinModel = "binary" no hazard model is used, so the interaction settings give the
  # same models:
  expect_equal(diagnostics$variants$initiateWithinModel, rep(c("binary", "hazard", "hazard"), 2))
  expect_equal(diagnostics$variants$covariateTimeInteraction, rep(c("day0", "day0", "none"), 2))
  expect_equal(diagnostics$variants$graceRegime, rep(c("binaryAtG", "naturalGrace"), each = 3))
  expect_equal(diagnostics$variants$specified, c(TRUE, rep(FALSE, 5)))
  expect_equal(diagnostics$balance$variant, diagnostics$variants$variant)
  overall <- diagnostics$calibration[diagnostics$calibration$grouping == "overall", ]
  expect_equal(
    paste(overall$strategyId, overall$model, overall$covariateTimeInteraction, overall$graceRegime),
    paste(c("1 binary NA", "2 binary NA", "1 hazard day0", "2 hazard day0", "1 hazard none", "2 hazard none"), rep(c("binaryAtG", "naturalGrace"), each = 6))
  )
})

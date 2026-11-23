library(testthat)

# Regression tests for issues found by the Phase 1 validation runs (simulated data only).

test_that("A variance that is 0 up to rounding gives a standard error of 0 and no p-value", {
  vcov <- diag(c(0.25, 1))
  expect_equal(.standardError(-3e-13, vcov), 0)
  expect_equal(.standardError(0.04, vcov), 0.2)
  expect_equal(.standardError(-1e-12, vcov, gradient = c(1, -1)), 0)
  expect_true(is.na(.standardError(NA_real_, vcov)))
  # A variance that is clearly negative is not hidden:
  expect_warning(se <- .standardError(-0.01, vcov), "NaN")
  expect_true(is.nan(se))

  row <- .ratioRow("or", -8.5e-13, 0)
  expect_equal(row$estimate, 1)
  expect_equal(row$lb95, row$estimate)
  expect_equal(row$ub95, row$estimate)
  expect_true(is.na(row$p))
  row <- .differenceRow("rd_0", 0, 0)
  expect_equal(c(row$lb95, row$ub95), c(0, 0))
  expect_true(is.na(row$p))
  expect_equal(.ratioRow("hr", log(2), 0.5)$p, 2 * pnorm(-log(2) / 0.5))
})

test_that("The saturated pooled logistic model reports the null first-interval odds ratio without NaN", {
  # initiateWithin(7) vs notInitiateWithin(7) with "eventsFirst": on day 0 both strategies have the same
  # rows, weights and outcomes, so the odds ratio of the first interval is 1 and its robust variance 0
  # (-3e-13 after rounding on these data, whose square root was NaN with a warning).
  data <- simulateCloneCensorWeightData(n = 1000, nCovariates = 1, graceEnd = 7, maxDays = 30, initiationHazard = 0.02, outcomeHazard = 0.004, seed = 1)
  population <- suppressMessages(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 30, endAnchor = "cohort start")))
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = list(
    createTreatmentStrategy(1, "Initiate within 7", "initiateWithin", graceEnd = 7),
    createTreatmentStrategy(2, "Not within 7", "notInitiateWithin", graceEnd = 7)
  ))))
  # The settings this test was written with (the defaults before contract section 12): with other
  # censoring models the rounding noise of the null variance can be positive (about 1e-12), which
  # .standardError() does not set to 0.
  models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(prior = Cyclops::createPrior("none"), postGraceModel = "binary", covariateTimeInteraction = "none")))
  basis <- createTimeBasis(type = "indicator", intervalLength = 1)
  weighted <- suppressMessages(computeCensoringWeights(clones, models, outcomeTimeBasis = basis))
  close(clones)
  close(data)
  warnings <- character(0)
  fit <- withCallingHandlers(
    suppressMessages(fitOutcomeModel(weighted, createFitOutcomeModelArgs(modelType = "logistic", horizons = c(0, 30), timeBasis = basis, armByTimeInteraction = TRUE))),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  close(weighted)
  expect_false(any(grepl("NaN", warnings)))
  estimates <- as.data.frame(fit$estimates)
  rownames(estimates) <- estimates$estimand
  expect_lt(abs(fit$vcov["target", "target"]), 1e-10)
  expect_equal(estimates["or", "estimate"], 1, tolerance = 1e-8)
  expect_equal(estimates["or", "seLogEstimate"], 0)
  expect_equal(estimates["or", "lb95"], estimates["or", "estimate"])
  expect_equal(estimates["or", "ub95"], estimates["or", "estimate"])
  expect_true(is.na(estimates["or", "p"]))
  # The day-0 risks are equal, so their difference and ratio have no variance either:
  expect_equal(estimates["risk_target_0", "estimate"], estimates["risk_comparator_0", "estimate"], tolerance = 1e-8)
  expect_lt(estimates["rd_0", "se"], 1e-8)
  expect_lt(estimates["rr_0", "seLogEstimate"], 1e-6)
  expect_false(any(is.nan(unlist(estimates[, c("seLogEstimate", "se", "lb95", "ub95", "p")]))))
  expect_gt(estimates["rd_30", "se"], 0)
})

test_that("The first-interval odds ratio is NA when the first interval has no outcomes", {
  # With seed 6 no clone has the outcome on day 0. Both first-interval cells then have a fixed hazard of
  # 0 and are left out of the fit, and the target coefficient was the log odds ratio of a later interval
  # chosen by the column pivoting (0.80, SE 1.23), reported as the first-interval odds ratio.
  data <- simulateCloneCensorWeightData(n = 1000, nCovariates = 1, graceEnd = 7, maxDays = 30, initiationHazard = 0.05, outcomeHazard = 0.002, seed = 6)
  population <- suppressMessages(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 30, endAnchor = "cohort start")))
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = list(
    createTreatmentStrategy(1, "Initiate within 7", "initiateWithin", graceEnd = 7),
    createTreatmentStrategy(2, "Not within 7", "notInitiateWithin", graceEnd = 7)
  ))))
  on.exit({
    close(clones)
    close(data)
  })
  expect_equal(clones$clones |> dplyr::filter(outcomeCount > 0, survivalTime <= 1) |> dplyr::count() |> dplyr::pull(), 0)
  models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(useCovariates = FALSE)))
  basis <- createTimeBasis(type = "indicator", intervalLength = 1)
  weighted <- suppressMessages(computeCensoringWeights(clones, models, outcomeTimeBasis = basis))
  on.exit(close(weighted), add = TRUE)
  logistic <- suppressMessages(fitOutcomeModel(weighted, createFitOutcomeModelArgs(modelType = "logistic", horizons = 30, timeBasis = basis, armByTimeInteraction = TRUE)))
  km <- suppressMessages(fitOutcomeModel(weighted, createFitOutcomeModelArgs(modelType = "km", horizons = 30)))
  estimates <- as.data.frame(logistic$estimates)
  rownames(estimates) <- estimates$estimand
  expect_true(all(is.na(unlist(estimates["or", c("estimate", "logEstimate", "seLogEstimate", "lb95", "ub95", "p")]))))
  # The target coefficient is another interval's log odds ratio, so it has no likelihood profile either:
  expect_null(logistic$logLikelihoodProfile)
  # The risks are not affected, and equal the weighted Kaplan-Meier risks:
  kmEstimates <- as.data.frame(km$estimates)
  rownames(kmEstimates) <- kmEstimates$estimand
  for (estimand in c("risk_target_30", "risk_comparator_30")) {
    expect_equal(estimates[estimand, "estimate"], kmEstimates[estimand, "estimate"], tolerance = 1e-6, info = estimand)
    expect_gt(estimates[estimand, "se"], 0)
  }
  # Without the interaction, the odds ratio is the common odds ratio of all intervals, which the other
  # intervals identify:
  common <- suppressMessages(fitOutcomeModel(weighted, createFitOutcomeModelArgs(modelType = "logistic", horizons = 30, timeBasis = basis, armByTimeInteraction = FALSE)))
  expect_true(is.finite(common$estimates$logEstimate[common$estimates$estimand == "or"]))
  expect_gt(common$estimates$seLogEstimate[common$estimates$estimand == "or"], 0)
})

test_that("A censoring model that cannot be fitted is named in the error", {
  fixture <- createM3Fixture(150, riskStart = 0L, seed = 12)
  # A covariate that is the indicator of starting by day 7:
  initiated <- fixture$persons$rowId[!is.na(fixture$persons$treatmentDay) & fixture$persons$treatmentDay <= 7]
  fixture$covariates <- rbind(fixture$covariates[fixture$covariates$covariateId != 1001, ], data.frame(rowId = initiated, covariateId = 1001, covariateValue = 1))
  data <- createM3Data(fixture)
  args <- m3Args("logistic")
  args$errorOnHighCorrelation <- TRUE
  # The binary model is fitted by CohortMethod::createPs(), whose error speaks of the "treatment":
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = 7), "eventsFirst")
  expect_error(
    suppressMessages(fitCensoringModels(clones, data, args)),
    "^Fitting the binary censoring model at day 7 for strategy 1 \\(initiateWithin 1\\) failed: .*High correlation"
  )
  close(clones)
  clones <- buildM3Clones(fixture$persons, m3Strategies("notInitiateWithin", "neverInitiate", graceEnd = 7), "eventsFirst")
  expect_error(
    suppressMessages(fitCensoringModels(clones, data, args)),
    "^Fitting the logistic initiation hazard model through day 7 for strategy 1 \\(notInitiateWithin 1\\) failed: High correlation"
  )
  close(clones)
  close(data)
})

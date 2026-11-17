library(testthat)

# The robust variance of the outcome MSMs from Cyclops' vcov(type = "robust") (plan section 5, C2; Phase 2
# item 1). Exit test of Phase 2 (plan section 8, item 8, first bullet): on clone data made by the
# package, the robust standard errors match survival::coxph(cluster = , weights = , ties = "breslow") and
# sandwich::vcovCL(glm(), cluster = , type = "HC0", cadjust = FALSE) to 1e-6, for both outcome models,
# with and without strata, and match this package's own sandwich (robustVarianceEngine = "ccw").

# Simulated data with two covariates, initiateWithin(5) vs notInitiateWithin(5), censoring models without
# a prior, and unstabilized weights (the weights of the MSMs, msmWeights = "unstabilized") on daily
# segments, which are the clone x day rows of the pooled logistic model. Every person has two clones, one
# per strategy, so the clusters span the strategies. The caller closes the result.
robustWeightedClones <- function(timeBasis) {
  data <- simulateCloneCensorWeightData(
    n = 400, nCovariates = 2, graceEnd = 5, maxDays = 30, initiationHazard = 0.08, outcomeHazard = 0.012,
    deathHazard = 0.005, seed = 11
  )
  on.exit(close(data))
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 30, endAnchor = "cohort start")
  ))
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within 5 days", "initiateWithin", graceEnd = 5),
    createTreatmentStrategy(2, "Do not initiate within 5 days", "notInitiateWithin", graceEnd = 5)
  )
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = strategies)))
  on.exit(close(clones), add = TRUE)
  models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(prior = Cyclops::createPrior("none"))))
  return(suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = FALSE), outcomeTimeBasis = timeBasis)))
}

robustSegments <- function(weightedClones) {
  segments <- as.data.frame(dplyr::collect(weightedClones$segments))
  segments <- segments[segments$weight > 0, ]
  segments <- segments[order(segments$cloneRowId, segments$startTime), ]
  rownames(segments) <- NULL
  segments$x <- as.integer(segments$strategyId == 1)
  # Day of follow-up (0-based) of a daily segment, and a stratum per person:
  segments$k <- segments$stopTime - 1
  segments$stratum <- segments$personSeqId %% 3 + 1
  return(segments)
}

robustControl <- Cyclops::createControl(algorithm = "newton", convergenceType = "gradient", tolerance = 1e-12, noiseLevel = "silent")

test_that("Cox MSM: Cyclops robust SE matches coxph(cluster) and the package's own sandwich, with and without strata", {
  skip_if_not_installed("survival")
  weightedClones <- robustWeightedClones(createTimeBasis("linear", intervalLength = 1))
  on.exit(close(weightedClones))
  segments <- robustSegments(weightedClones)
  # Every person has clones in both strategies:
  expect_true(all(table(unique(segments[, c("personSeqId", "strategyId")])$personSeqId) == 2))
  expect_gt(sum(segments$outcome), 20)
  expect_gt(sum(segments$competing), 5)

  fitArgs <- function(...) createFitOutcomeModelArgs(modelType = "cox", horizons = 20, ...)
  cyclopsModel <- suppressMessages(fitOutcomeModel(weightedClones, fitArgs()))
  ccwModel <- suppressMessages(fitOutcomeModel(weightedClones, fitArgs(robustVarianceEngine = "ccw")))
  control <- survival::coxph.control(eps = 1e-12, toler.chol = 1e-14, iter.max = 100)
  reference <- survival::coxph(
    survival::Surv(startTime, stopTime, outcome) ~ x,
    data = segments, weights = weight, cluster = personSeqId, ties = "breslow", control = control
  )
  expect_equal(cyclopsModel$coefficients[["target"]], unname(coef(reference)), tolerance = 1e-8)
  expect_equal(unname(cyclopsModel$vcov), unname(vcov(reference)), tolerance = 1e-6)
  expect_equal(getM4Estimate(cyclopsModel, "hr", "seLogEstimate"), sqrt(vcov(reference)[1, 1]), tolerance = 1e-6)
  # The robust variance differs from the model-based one, so the test can tell them apart:
  expect_gt(abs(vcov(reference)[1, 1] / reference$naive.var[1, 1] - 1), 0.01)
  expect_equal(cyclopsModel$vcov, ccwModel$vcov, tolerance = 1e-8)
  expect_equal(cyclopsModel$estimates, ccwModel$estimates, tolerance = 1e-8)

  # With strata: Cyclops on the same clone segments, with the persons as subjects.
  outcomes <- data.frame(
    rowId = seq_len(nrow(segments)), stratumId = segments$stratum, subjectId = segments$personSeqId,
    y = segments$outcome, time = segments$stopTime, startTime = segments$startTime, weights = segments$weight
  )
  covariates <- data.frame(rowId = which(segments$x == 1), covariateId = 1, covariateValue = 1)
  cyclopsData <- Cyclops::convertToCyclopsData(outcomes, covariates, modelType = "cox_counting", addIntercept = FALSE, quiet = TRUE)
  fit <- Cyclops::fitCyclopsModel(cyclopsData, control = robustControl)
  stratifiedReference <- survival::coxph(
    survival::Surv(startTime, stopTime, outcome) ~ x + strata(stratum),
    data = segments, weights = weight, cluster = personSeqId, ties = "breslow", control = control
  )
  expect_equal(unname(coef(fit)), unname(coef(stratifiedReference)), tolerance = 1e-8)
  expect_equal(unname(vcov(fit, type = "robust")), unname(vcov(stratifiedReference)), tolerance = 1e-6)
  expect_gt(abs(vcov(stratifiedReference)[1, 1] / vcov(reference)[1, 1] - 1), 1e-3)
})

test_that("Pooled logistic MSM: Cyclops robust covariance matches vcovCL and the package's own sandwich, with and without strata", {
  skip_if_not_installed("sandwich")
  timeBasis <- createTimeBasis("linear", intervalLength = 1)
  weightedClones <- robustWeightedClones(timeBasis)
  on.exit(close(weightedClones))
  rows <- robustSegments(weightedClones)
  # Daily segments are the clone x day rows:
  expect_true(all(rows$stopTime - rows$startTime == 1))

  fitArgs <- function(...) {
    createFitOutcomeModelArgs(modelType = "logistic", horizons = c(5, 20), timeBasis = timeBasis, armByTimeInteraction = TRUE, ...)
  }
  cyclopsModel <- withoutGracePeriodWarnings(suppressMessages(fitOutcomeModel(weightedClones, fitArgs())))
  ccwModel <- withoutGracePeriodWarnings(suppressMessages(fitOutcomeModel(weightedClones, fitArgs(robustVarianceEngine = "ccw"))))

  # Young et al. (2020): the outcome model on the rows without death and the death model on all rows,
  # whose joint covariance is vcovCL of the stacked glm.
  design <- cbind(1, rows$x, rows$k, rows$x * rows$k)
  noDeath <- rows$competing == 0
  stackedDesign <- rbind(cbind(design[noDeath, ], 0 * design[noDeath, ]), cbind(0 * design, design))
  stackedY <- c(rows$outcome[noDeath], rows$competing)
  stackedWeight <- c(rows$weight[noDeath], rows$weight)
  stackedCluster <- c(rows$personSeqId[noDeath], rows$personSeqId)
  glmControl <- glm.control(epsilon = 1e-14, maxit = 100)
  stackedFit <- glm(stackedY ~ stackedDesign - 1, family = quasibinomial(), weights = stackedWeight, control = glmControl)
  reference <- sandwich::vcovCL(stackedFit, cluster = stackedCluster, type = "HC0", cadjust = FALSE)
  expect_equal(unname(cyclopsModel$coefficients), unname(coef(stackedFit)), tolerance = 1e-8)
  expect_equal(unname(cyclopsModel$vcov), unname(reference), tolerance = 1e-6)
  expect_equal(getM4Estimate(cyclopsModel, "or", "seLogEstimate"), sqrt(reference[2, 2]), tolerance = 1e-6)
  # Clustering on person matters here (the clones of a person share it):
  rowLevel <- sandwich::vcovCL(stackedFit, cluster = seq_along(stackedY), type = "HC0", cadjust = FALSE)
  expect_gt(max(abs(diag(rowLevel) / diag(reference) - 1)), 0.01)
  expect_equal(cyclopsModel$vcov, ccwModel$vcov, tolerance = 1e-8)
  expect_equal(cyclopsModel$estimates, ccwModel$estimates, tolerance = 1e-8)

  # With strata (an intercept per stratum): Cyclops on the same rows, with the persons as subjects.
  strataDesign <- cbind(design[, -1], rows$stratum == 2, rows$stratum == 3)
  outcomes <- data.frame(rowId = seq_len(nrow(rows)), y = rows$outcome, weights = rows$weight, subjectId = rows$personSeqId)[noDeath, ]
  nonZero <- which(strataDesign[noDeath, ] != 0, arr.ind = TRUE)
  covariates <- data.frame(
    rowId = outcomes$rowId[nonZero[, 1]], covariateId = nonZero[, 2], covariateValue = strataDesign[noDeath, ][nonZero]
  )
  cyclopsData <- Cyclops::convertToCyclopsData(outcomes, covariates, modelType = "lr", addIntercept = TRUE, quiet = TRUE)
  fit <- Cyclops::fitCyclopsModel(cyclopsData, control = robustControl)
  stratifiedFit <- glm(rows$outcome[noDeath] ~ strataDesign[noDeath, ], family = quasibinomial(), weights = rows$weight[noDeath], control = glmControl)
  stratifiedReference <- sandwich::vcovCL(stratifiedFit, cluster = rows$personSeqId[noDeath], type = "HC0", cadjust = FALSE)
  expect_equal(unname(coef(fit)), unname(coef(stratifiedFit)), tolerance = 1e-8)
  expect_equal(unname(vcov(fit, type = "robust")), unname(stratifiedReference), tolerance = 1e-6)
})

test_that("Without competing events the pooled logistic robust covariance is that of the outcome model alone", {
  skip_if_not_installed("sandwich")
  fixture <- createM4Fixture(n = 200, riskStart = 1, competingEvents = "none", seed = 3)
  timeBasis <- createTimeBasis("polynomial", intervalLength = 2, degree = 2)
  cyclopsModel <- fitM4(fixture, modelType = "logistic", horizons = c(6, 12), timeBasis = timeBasis, armByTimeInteraction = TRUE)
  ccwModel <- fitM4(fixture, modelType = "logistic", horizons = c(6, 12), timeBasis = timeBasis, armByTimeInteraction = TRUE, robustVarianceEngine = "ccw")
  expect_equal(dim(cyclopsModel$vcov), c(6, 6))
  expect_equal(cyclopsModel$vcov, ccwModel$vcov, tolerance = 1e-8)
  expect_equal(cyclopsModel$estimates, ccwModel$estimates, tolerance = 1e-8)
})

test_that("Cells with a fixed hazard and coefficients that are not identified are left out of Cyclops' robust covariance", {
  # A saturated daily model (indicator basis with the arm-by-time interaction) with competing death, where
  # some strategy x day cells have no outcome or no death (a fixed hazard of 0):
  fixture <- createM4Fixture(n = 200, riskStart = 0, competingEvents = "death", seed = 13, maxDays = 10)
  rows <- expandM4Rows(fixture$segments, 1)
  cells <- aggregate(cbind(y, death) ~ x + k, data = rows, FUN = sum)
  expect_true(any(cells$y == 0) && any(cells$death == 0))
  timeBasis <- createTimeBasis("indicator", intervalLength = 1)
  fit <- function(fixture, timeBasis, engine) {
    suppressWarnings(fitM4(fixture, modelType = "logistic", horizons = c(4, 9), timeBasis = timeBasis, armByTimeInteraction = TRUE, robustVarianceEngine = engine))
  }
  cyclopsModel <- fit(fixture, timeBasis, "cyclops")
  ccwModel <- fit(fixture, timeBasis, "ccw")
  expect_identical(is.na(cyclopsModel$vcov), is.na(ccwModel$vcov))
  expect_equal(cyclopsModel$vcov, ccwModel$vcov, tolerance = 1e-8)
  # The odds ratio of the first interval, where both strategies have the same rows, is null by design: its
  # variance is 0 up to rounding (here about 1e-13 with either engine, of either sign), so its standard
  # error is 0 or of the order of 1e-7.
  expect_equal(cyclopsModel$estimates$estimate, ccwModel$estimates$estimate, tolerance = 1e-8)
  expect_lt(abs(cyclopsModel$vcov["target", "target"]), 1e-12)
  expect_lt(getM4Estimate(cyclopsModel, "or", "seLogEstimate"), 1e-6)
  expect_equal(cyclopsModel$estimates[-1, ], ccwModel$estimates[-1, ], tolerance = 1e-8)

  # Target clones that end by t = 10 while comparator clones run to t = 15: the target's hazard in the
  # last 5-day interval is not identified (test-m4-outcomes.R).
  fixture <- createM4Fixture(n = 200, riskStart = 0, competingEvents = "none", seed = 3, maxDays = 14, outcomeRate = 0.7)
  drop <- fixture$clones$cloneRowId[fixture$clones$strategyId == 1 & fixture$clones$survivalTime > 10]
  fixture$segments <- fixture$segments[!fixture$segments$cloneRowId %in% drop, ]
  fixture$clones <- fixture$clones[!fixture$clones$cloneRowId %in% drop, ]
  oldIds <- fixture$clones$cloneRowId
  fixture$clones$cloneRowId <- seq_along(oldIds)
  fixture$segments$cloneRowId <- match(fixture$segments$cloneRowId, oldIds)
  fixture$segments$segmentId <- seq_len(nrow(fixture$segments))
  fixture$personWeights <- m4PersonWeights(fixture$segments)
  timeBasis <- createTimeBasis("indicator", intervalLength = 5)
  cyclopsModel <- fit(fixture, timeBasis, "cyclops")
  ccwModel <- fit(fixture, timeBasis, "ccw")
  expect_true(all(is.na(cyclopsModel$vcov["target:time2", ])))
  expect_identical(is.na(cyclopsModel$vcov), is.na(ccwModel$vcov))
  expect_equal(cyclopsModel$vcov, ccwModel$vcov, tolerance = 1e-8)
  expect_equal(cyclopsModel$estimates, ccwModel$estimates, tolerance = 1e-8)
})

test_that("The default robust variance stops with a clear error when the loaded Cyclops has no vcov(type = \"robust\")", {
  # The loaded (development) Cyclops has it:
  expect_silent(.checkCyclopsRobustVcov())
  # vcov.cyclopsFit() of Cyclops releases up to 4.0.1, which ignores `type`:
  oldVcov <- function(object, control, overrideNoRegularization = FALSE, ...) NULL
  expect_error(.checkCyclopsRobustVcov(oldVcov), "robustVarianceEngine = \"ccw\"")
  expect_error(.checkCyclopsRobustVcov(oldVcov), "does not have")

  fixture <- createM4Fixture(n = 60, riskStart = 0, competingEvents = "none")
  realCheck <- .checkCyclopsRobustVcov
  local_mocked_bindings(.checkCyclopsRobustVcov = function() realCheck(oldVcov))
  expect_error(fitM4(fixture, modelType = "cox", horizons = 10), "robustVarianceEngine = \"ccw\"")
  expect_error(fitM4(fixture, modelType = "logistic", horizons = 10, timeBasis = createTimeBasis("linear", intervalLength = 1)), "robustVarianceEngine = \"ccw\"")
  # Not needed by the package's own sandwich, the model-based variance or the Kaplan-Meier:
  expect_false(is.na(getM4Estimate(fitM4(fixture, modelType = "cox", horizons = 10, robustVarianceEngine = "ccw"), "hr", "seLogEstimate")))
  expect_false(is.na(getM4Estimate(fitM4(fixture, modelType = "cox", horizons = 10, varianceEstimator = "modelBased"), "hr", "seLogEstimate")))
  expect_false(is.na(getM4Estimate(fitM4(fixture, modelType = "km", horizons = 10), "risk_target_10", "se")))
})

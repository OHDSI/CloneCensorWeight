# The stacked sandwich variance (varianceEstimator = "stackedSandwich"; Phase4Designs.md,
# stacked-sandwich, verification tier 1): the binary model at G and the daily logistic hazard model,
# unpenalized, against the hand-written per-entry reference refStackedSandwich().

test_that("The stacked sandwich equals the hand-written stacked influence, also for a person with two entries, and its Kaplan-Meier SEs do not depend on stabilization", {
  graceEnd <- 3
  horizon <- 10
  fixture <- createM3Fixture(150, seed = 1, deathFraction = 0)
  # Entries 1 and 2 (other covariates and days) are two index dates of one person (firstOnly = FALSE):
  fixture$persons$personSeqId[2] <- 1L
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd), competingEvents = "none")
  on.exit({
    close(clones)
    close(data)
  })
  args <- m3Args("logistic", timeBasis = createTimeBasis(type = "indicator"), covariateTimeInteraction = "day0")
  models <- suppressMessages(fitCensoringModels(clones, data, args))
  fitStacked <- function(stabilized) {
    weightedClones <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(
      stabilized = stabilized,
      truncation = createWeightTruncation(type = "none")
    )))
    on.exit(close(weightedClones))
    fitOutcomeModelArgs <- createFitOutcomeModelArgs(horizons = horizon, varianceEstimator = "stackedSandwich")
    return(suppressMessages(fitOutcomeModel(weightedClones, fitOutcomeModelArgs, cloneCensorWeightData = data)))
  }
  model <- fitStacked(stabilized = TRUE)
  reference <- refStackedSandwich(fixture, graceEnd, horizon)

  estimate <- function(estimand, column) model$estimates[[column]][model$estimates$estimand == estimand]
  risks <- sprintf(c("risk_target_%s", "risk_comparator_%s"), horizon)
  expect_equal(c(estimate("hr", "logEstimate"), estimate(risks[1], "estimate"), estimate(risks[2], "estimate")), reference$estimates, tolerance = 1e-6)
  se <- function(influence) sqrt(sum(influence^2))
  # The reference's influences are per entry (rowId); the person is the cluster:
  person <- fixture$persons$personSeqId[order(fixture$persons$rowId)]
  influence <- rowsum(reference$influence, person)
  expect_equal(estimate("hr", "seLogEstimate"), se(influence[, 1]), tolerance = 1e-6)
  expect_equal(estimate(risks[1], "se"), se(influence[, 2]), tolerance = 1e-6)
  expect_equal(estimate(risks[2], "se"), se(influence[, 3]), tolerance = 1e-6)
  expect_equal(estimate(sprintf("rd_%s", horizon), "se"), se(influence[, 2] - influence[, 3]), tolerance = 1e-6)
  expect_lt(max(abs(colSums(influence))), 1e-10)
  # The diagnostics: the robust SEs (of the IJ alone) and the ratio:
  diagnostics <- model$varianceDiagnostics$estimands
  expect_equal(diagnostics$estimand, model$estimates$estimand)
  robust <- diagnostics$seRobust[match(c("hr", risks), diagnostics$estimand)]
  expect_equal(robust, apply(rowsum(reference$robustInfluence, person), 2, se), tolerance = 1e-6)
  expect_equal(diagnostics$ratio, diagnostics$seStacked / diagnostics$seRobust)
  expect_message(utils::capture.output(print(model)), "Stacked sandwich SE / robust SE")

  # Untruncated Kaplan-Meier risks do not depend on the numerators, nor do their stacked SEs (the weights
  # once more, so in the full tier):
  if (isFullTestRun()) {
    unstabilized <- fitStacked(stabilized = FALSE)
    kmRows <- !model$estimates$estimand %in% "hr"
    expect_equal(unstabilized$estimates$se[kmRows], model$estimates$se[kmRows], tolerance = 1e-10)
  }
})

test_that("The generalized inverse drops rounding-level pivots of an exactly singular information", {
  # A degree-5 polynomial time basis on 7 intervals with rows in only the first 4: 6 columns of rank 4.
  # LAPACK's default tolerance keeps a rounding-level pivot here, and the inverse is then not a
  # generalized inverse.
  basis <- cbind(1, stats::poly(0:6, degree = 5))[1:4, ]
  information <- crossprod(basis * sqrt(c(10, 5, 2, 1)))
  inverse <- .generalizedInverse(information)
  expect_equal(inverse$rank, 4L)
  expect_equal(information %*% inverse$inverse %*% information, information, tolerance = 1e-10)
})

test_that("A Laplace prior's support treats the interactions of an excluded covariate as penalized, as the fit does, and its label does not depend on the support", {
  # 1003 is excluded from the prior (kept at 0), but its interaction is penalized (dropped at 0):
  coefficients <- c("(Intercept)" = -2, "1003:day0" = 0, "1002:day0" = 0.1, "1001" = 0, "1003" = 0)
  prior <- list(type = "laplace", exclude = c("0", "1003"))
  expect_equal(.stackedSupport(coefficients, prior, "(Intercept)"), c("(Intercept)", "1002:day0", "1003"))
  # A Laplace model whose support has no penalized parameter (every covariate at 0) is still conditional on it:
  block <- .stackedBlock(1L, FALSE, "binary model at G", 1:4, 1:4, matrix(1, 4, 1), matrix(1, 1, 1), 1L, 1L,
                         matrix(0.5, 4, 1), matrix(1, 4, 1), matrix(c(0, 1, 0, 1)), list(),
                         penalty = .stackedPenalty(prior, 0.1, NA))
  expect_equal(block$diagnostics[c("prior", "conditionalOnSupport", "maxKktGap")], list(prior = "laplace", conditionalOnSupport = TRUE, maxKktGap = NA_real_))
})

test_that("The stacked sandwich needs the covariates, and stops for what it does not support", {
  expect_error(createFitOutcomeModelArgs(modelType = "fineGray", varianceEstimator = "stackedSandwich"), "fineGray")
  expect_error(
    createCcwAnalysis(
      analysisId = 1,
      getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(),
      createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(),
      createClonesArgs = createCreateClonesArgs(strategies = m3Strategies()),
      fitCensoringModelsArgs = createFitCensoringModelsArgs(hazardModelType = "cox"),
      computeWeightsArgs = createComputeWeightsArgs(),
      fitOutcomeModelArgs = createFitOutcomeModelArgs(varianceEstimator = "stackedSandwich")
    ),
    "Cox censoring model"
  )
  weightedClones <- createM4WeightedClones(createM4Fixture())
  on.exit(close(weightedClones))
  expect_error(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(varianceEstimator = "stackedSandwich")), "cloneCensorWeightData")
})

# The strategy and weight-scale cases of the weight derivatives (design sections 2, 4 and 5), each against
# numerical derivatives of the weights of computeCensoringWeights() (refStackedWeightDerivatives()). Each
# recomputes the weights three times (about 5 s), so they run in the full tier.

test_that("The stacked sandwich's weight derivatives are those of the weights: graceStart, naturalGrace, postGraceModel = \"binary\", truncateMsmWeights = FALSE", {
  skipUnlessFullTests()
  fixture <- createM3Fixture(150, seed = 1, deathFraction = 0)
  data <- createM3Data(fixture)
  strategies <- list(
    createTreatmentStrategy(strategyId = 1, label = "start on days 2-5", type = "initiateWithin", graceStart = 2, graceEnd = 5),
    createTreatmentStrategy(strategyId = 2, label = "no start by day 5", type = "notInitiateWithin", graceEnd = 5)
  )
  clones <- buildM3Clones(fixture$persons, strategies, competingEvents = "none")
  on.exit({
    close(clones)
    close(data)
  })
  args <- m3Args("logistic", timeBasis = createTimeBasis(type = "indicator"), covariateTimeInteraction = "day0", postGraceModel = "binary")
  args$graceRegime <- "naturalGrace"
  models <- suppressWarnings(suppressMessages(fitCensoringModels(clones, data, args)))
  # A fixed cap that binds on the unstabilized weights, with the stabilized MSM weights formed from them:
  weightsArgs <- createComputeWeightsArgs(stabilized = FALSE, truncation = createWeightTruncation(type = "fixed", upper = 3), truncateMsmWeights = FALSE)
  set.seed(1)
  # The numerator of the naturalGrace binary model is the share of starters of the clones' persons, not
  # stored in the model:
  derivatives <- refStackedWeightDerivatives(clones, models, weightsArgs, data, skip = "1 binary model at G numerator")
  expect_lt(derivatives$weight$difference, 1e-5)
  expect_lt(derivatives$otherScaleWeight$difference, 1e-5)
  expect_gt(derivatives$weight$truncated, 0)
  diagnostics <- derivatives$diagnostics
  expect_setequal(diagnostics$component[diagnostics$strategyId == 1], c("binary model at G", "binary model at G numerator", "hazard model before graceStart", "hazard model before graceStart numerator"))
  expect_true("binary model at G (postGraceModel)" %in% diagnostics$component[diagnostics$strategyId == 2])
  caps <- derivatives$caps
  expect_equal(diagnostics$truncatedShare, caps$truncatedShare[diagnostics$strategyId])
})

test_that("The stacked sandwich's weight derivatives are those of the weights: initiateWithinModel = \"hazard\", deviationFirst, a Laplace prior that excludes a covariate, a sampled fit, truncateMsmWeights = TRUE", {
  skipUnlessFullTests()
  fixture <- createM3Fixture(150, seed = 1, deathFraction = 0)
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "neverInitiate", 3), sameDayOrdering = "deviationFirst", competingEvents = "none")
  on.exit({
    close(clones)
    close(data)
  })
  # Covariate 1003 is excluded from the prior, but not its day-0 interaction, which the fit penalizes (and
  # sets to a nonzero value at this variance):
  prior <- Cyclops::createPrior("laplace", variance = 0.5, exclude = c(0, 1003))
  args <- m3Args("logistic", timeBasis = createTimeBasis(type = "polynomial", degree = 2), covariateTimeInteraction = "day0", prior = prior)
  args$initiateWithinModel <- "hazard"
  args$maxCohortSizeForFitting <- 100
  models <- suppressMessages(fitCensoringModels(clones, data, args))
  # A fixed cap that binds on the unstabilized MSM weights (their own scale):
  weightsArgs <- createComputeWeightsArgs(stabilized = TRUE, truncation = createWeightTruncation(type = "fixed", upper = 3))
  set.seed(1)
  derivatives <- suppressWarnings(refStackedWeightDerivatives(clones, models, weightsArgs, data))
  expect_lt(derivatives$weight$difference, 1e-5)
  expect_lt(derivatives$otherScaleWeight$difference, 1e-5)
  expect_gt(derivatives$otherScaleWeight$truncated, 0)
  # The covariate models are conditional on their support, at the optimum (KKT gap 0, which includes the
  # interaction of the excluded covariate), and fit on a sample:
  expect_true(all(vapply(models, function(model) model$coefficients[["1003:day0"]] != 0, logical(1))))
  diagnostics <- derivatives$diagnostics
  covariateModels <- diagnostics[diagnostics$prior == "laplace", ]
  expect_equal(nrow(covariateModels), 2)
  expect_true(all(covariateModels$conditionalOnSupport))
  expect_lt(max(covariateModels$maxKktGap), 1e-4)
  expect_true(all(diagnostics$sampled))
  expect_true(all(diagnostics$rank == diagnostics$supportSize))
})

test_that("The stacked sandwich supports hazard models without covariate-by-time interactions, and stops for sampled binary models, plug-in intervals and censoring models of an earlier version", {
  fixture <- createM3Fixture(150, seed = 1, deathFraction = 0)
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", 3), competingEvents = "none")
  # covariateTimeInteraction = "none" leaves the hazard model an interaction basis without columns:
  args <- m3Args("logistic", useCovariates = FALSE, timeBasis = createTimeBasis(type = "indicator"), covariateTimeInteraction = "none")
  models <- suppressMessages(fitCensoringModels(clones, data, args))
  weightedClones <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(truncation = createWeightTruncation(type = "none"))))
  on.exit({
    close(weightedClones)
    close(clones)
    close(data)
  })
  outcomeData <- .getOutcomeData(weightedClones)
  blocks <- function(censoringModels) {
    outcomeData$weightsMetaData$censoringModels <- censoringModels
    return(.stackedCensoringBlocks(outcomeData, data, weightedClones))
  }
  expect_silent(blocks(models))
  model <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(horizons = 10, varianceEstimator = "stackedSandwich"), cloneCensorWeightData = data))
  expect_true(all(is.finite(ifelse(is.na(model$estimates$seLogEstimate), model$estimates$se, model$estimates$seLogEstimate))))
  sampled <- models
  sampled[[1]]$diagnostics$sampled <- TRUE
  expect_error(blocks(sampled), "maxCohortSizeForFitting = 0, or use varianceEstimator = \"bootstrap\"")
  # A sample without the initiations of an indicator interval (its hazard is then a plug-in):
  plugIn <- models
  predictions <- plugIn[[2]]$predictions
  plugIn[[2]]$predictions$inFit <- !(predictions$y == 1 & predictions$lastDay == 1)
  expect_error(blocks(plugIn), "plug-ins")
  earlier <- models
  earlier[[2]]$predictions$inFit <- NULL
  expect_error(blocks(earlier), "refit the censoring models with this version")
})

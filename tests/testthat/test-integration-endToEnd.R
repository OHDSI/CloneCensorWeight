library(testthat)

# End-to-end test of the single-study pipeline on simulated data with confounded initiation (DGP "a",
# riskStart = 0): study population, clones, censoring models, weights, covariate balance, the three
# outcome models, cumulative incidence and MDRR, for initiateWithin(G) vs notInitiateWithin(G) and
# initiateWithin(G) vs neverInitiate, under both same-day orderings. The weighted Kaplan-Meier
# (Aalen-Johansen) risks, of "km" and of the default "cox" report, must be within 4 robust standard
# errors of the true risks, and the unweighted analysis must not be. Initiation and the outcome both
# depend on the covariates, so the clones that follow each strategy are selected on prognosis: starters
# by G have a higher outcome risk, and persons who do not start a lower one.

test_that("The weighted pipeline recovers the true risks that the unweighted analysis misses", {
  skip_on_cran()
  # Full tier: a statistical check at n = 3000 with 5e5 persons for the truth. The default tier runs
  # every step of this pipeline in the other tests.
  skipUnlessFullTests()
  graceEnd <- 5
  maxDays <- 20
  horizons <- c(3, 13, 20)
  data <- simulateCloneCensorWeightData(
    n = 3000,
    nCovariates = 5,
    dgp = "a",
    graceEnd = graceEnd,
    maxDays = maxDays,
    initiationHazard = 0.03,
    outcomeHazard = 0.004,
    deathHazard = 0.001,
    observationEndHazard = 0.002,
    initiationCoefficients = rep(0.5, 5),
    outcomeCoefficients = rep(0.8, 5),
    deathCoefficients = rep(0.2, 5),
    treatmentEffect = log(0.5),
    seed = 3
  )
  on.exit(close(data))
  settings <- attr(data, "metaData")$simulationSettings
  initiate <- createTreatmentStrategy(1, "Initiate within 5 days", "initiateWithin", graceEnd = graceEnd)
  comparators <- list(
    notInitiateWithin = createTreatmentStrategy(2, "Do not initiate within 5 days", "notInitiateWithin", graceEnd = graceEnd),
    neverInitiate = createTreatmentStrategy(2, "Never initiate", "neverInitiate")
  )
  truth <- computeTrueRisks(settings, c(list(initiate), comparators), horizons = horizons, nTruth = 5e5, seed = 1, riskStart = 0)
  # Under DGP "a" the treatment acts only after G through "started by G", so not starting by G and never
  # starting have the same truth:
  expect_equal(truth$risk[truth$strategyId == 2], truth$risk[truth$strategyId == 3])

  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
      riskWindowStart = 0,
      riskWindowEnd = maxDays,
      endAnchor = "cohort start"
    )
  ))
  expect_equal(nrow(population), 3000)
  expect_equal(unique(population$riskStart), 0)
  fitCensoringModelsArgs <- createFitCensoringModelsArgs(prior = Cyclops::createPrior("none"))
  outcomeArgs <- list(
    km = createFitOutcomeModelArgs(modelType = "km", horizons = horizons),
    # The default report:
    cox = createFitOutcomeModelArgs(horizons = horizons),
    logistic = createFitOutcomeModelArgs(modelType = "logistic", horizons = horizons)
  )
  expect_equal(outcomeArgs$cox$modelType, "cox")

  for (comparator in names(comparators)) {
    strategies <- list(initiate, comparators[[comparator]])
    truthOf <- function(strategyId) {
      truthId <- if (strategyId == 1) 1 else if (comparator == "notInitiateWithin") 2 else 3
      return(truth$risk[truth$strategyId == truthId])
    }
    for (sameDayOrdering in c("eventsFirst", "deviationFirst")) {
      label <- paste(comparator, sameDayOrdering)
      clones <- suppressMessages(createClones(
        population = population,
        cloneCensorWeightData = data,
        createClonesArgs = createCreateClonesArgs(strategies = strategies, sameDayOrdering = sameDayOrdering)
      ))
      cloneTable <- dplyr::collect(clones$clones)
      expect_setequal(unique(cloneTable$strategyId), 1:2)
      expect_true(all(cloneTable$endReason %in% c("outcome", "death", "deviation", "observationEnd", "riskEnd")))
      expect_true(all(cloneTable$survivalTime >= 1 & cloneTable$endDay <= maxDays))
      expect_equal(cloneTable$outcomeCount, as.integer(cloneTable$endReason == "outcome"))
      if (sameDayOrdering == "eventsFirst") {
        expect_equal(nrow(cloneTable), 6000)
      } else {
        # Clones that deviate on day 0 have no follow-up:
        expect_lt(nrow(cloneTable), 6000)
      }
      attrition <- attr(clones, "metaData")$attrition
      expect_equal(utils::tail(attrition$description, 2), c("Clones created", "Clones with follow-up"))

      models <- suppressMessages(fitCensoringModels(clones, data, fitCensoringModelsArgs))
      expect_equal(vapply(models, function(model) model$modelType, character(1)), c("binary", "logistic"))
      expect_equal(vapply(models, function(model) model$diagnostics$status, character(1)), c("OK", "OK"))
      expect_gt(models[[1]]$diagnostics$auc, 0.6)
      # The hazard models' covariates are all associated with faster initiation:
      expect_true(all(models[[2]]$coefficients[c("1001", "2002", "3001", "4002", "5001")] > 0))

      weighted <- suppressMessages(computeCensoringWeights(
        clones, models, createComputeWeightsArgs(),
        outcomeTimeBasis = outcomeArgs$logistic$timeBasis
      ))
      weightsMetaData <- attr(weighted, "metaData")$weights
      expect_equal(weightsMetaData$cutIntervalLength, 1)
      expect_true(weightsMetaData$stabilized)
      segmentWeights <- dplyr::pull(weighted$segments, "weight")
      expect_true(all(is.finite(segmentWeights) & segmentWeights > 0))
      expect_true(all(segmentWeights <= max(weightsMetaData$caps$cap)))

      # Weighting balances the covariates of the clones still following each strategy after G:
      balance <- suppressMessages(computeCovariateBalance(weighted, data))
      expect_equal(nrow(balance), 5)
      expect_lt(max(abs(balance$afterMatchingStdDiff)), max(abs(balance$beforeMatchingStdDiff)) / 2)
      expect_lt(max(abs(balance$afterMatchingStdDiff)), 0.1)

      outcomeModels <- lapply(outcomeArgs, function(args) withoutGracePeriodWarnings(suppressMessages(fitOutcomeModel(weighted, args))))
      for (modelType in names(outcomeModels)) {
        estimates <- outcomeModels[[modelType]]$estimates
        for (h in horizons) {
          for (role in c("target", "comparator")) {
            row <- estimates[estimates$estimand == sprintf("risk_%s_%d", role, h), ]
            expect_true(row$estimate > 0 && row$estimate < 1, info = paste(label, modelType, h, role))
            expect_true(row$se > 0 && row$lb95 < row$estimate && row$estimate < row$ub95, info = paste(label, modelType, h, role))
          }
        }
        expect_false(anyNA(estimates$estimate), info = paste(label, modelType))
      }
      km <- outcomeModels$km$estimates
      expect_equal(outcomeModels$km$metaData$competingEvents, "death")
      # The effect: starting by G halves the odds of the outcome after G.
      expect_lt(km$estimate[km$estimand == "rd_20"], 0)
      expect_lt(outcomeModels$cox$estimates$estimate[1], 1)
      expect_lt(outcomeModels$logistic$estimates$estimate[1], 1)

      # The weighted Aalen-Johansen risks, of km and of the default cox report, are within 4 robust
      # standard errors of the truth:
      for (modelType in c("km", "cox")) {
        estimates <- outcomeModels[[modelType]]$estimates
        for (strategyId in 1:2) {
          role <- c("target", "comparator")[strategyId]
          rows <- match(sprintf("risk_%s_%d", role, horizons), estimates$estimand)
          z <- (estimates$estimate[rows] - truthOf(strategyId)) / estimates$se[rows]
          expect_true(all(abs(z) < 4), info = paste(label, modelType, role, paste(round(z, 2), collapse = ", ")))
        }
      }

      curves <- computeCumulativeIncidence(weighted, horizons = horizons)
      expect_equal(sort(unique(curves$strategyId)), 1:2)
      for (strategyId in 1:2) {
        curve <- curves[curves$strategyId == strategyId, ]
        expect_true(all(diff(curve$risk) >= 0))
        role <- c("target", "comparator")[strategyId]
        atHorizons <- curve$risk[match(horizons, curve$day)]
        expect_equal(atHorizons, km$estimate[match(sprintf("risk_%s_%d", role, horizons), km$estimand)], tolerance = 1e-10)
      }
      mdrr <- computeMdrr(weighted)
      expect_gt(mdrr$mdrr, 1)
      expect_equal(mdrr$targetExposures + mdrr$comparatorExposures, nrow(cloneTable))

      # The unweighted analysis: stabilized weights of models without covariates are all 1.
      unweightedModels <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(useCovariates = FALSE)))
      unweighted <- suppressMessages(computeCensoringWeights(clones, unweightedModels, createComputeWeightsArgs()))
      expect_equal(dplyr::pull(unweighted$segments, "weight"), rep(1, dplyr::pull(dplyr::count(unweighted$segments))), tolerance = 1e-12)
      unweightedKm <- suppressMessages(fitOutcomeModel(unweighted, outcomeArgs$km))$estimates
      zUnweighted <- c()
      for (strategyId in 1:2) {
        rows <- match(sprintf("risk_%s_%d", c("target", "comparator")[strategyId], horizons), unweightedKm$estimand)
        zUnweighted <- c(zUnweighted, (unweightedKm$estimate[rows] - truthOf(strategyId)) / unweightedKm$se[rows])
      }
      expect_gt(max(abs(zUnweighted)), 4, label = paste(label, paste(round(zUnweighted, 2), collapse = ", ")))
      close(clones)
    }
  }
})

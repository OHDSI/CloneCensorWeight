library(testthat)

# Contract section 11.5: createWeightTruncation(type = "percentile", basis = "varying") takes the
# percentile of the maximum weights over the clones under follow-up on day tau + 1 (whose weight varies);
# basis = "allClones" is the original rule. metaData$weights$caps reports the cap and the share of
# clones whose weight is truncated per strategy.

test_that("Percentile truncation over the clones whose weight varies truncates about the nominal share", {
  # The setting of the validation repro truncationPercentileBasis.R: few persons start by G = 1, so most
  # clones of initiateWithin(1) keep a weight of 1, and a percentile over all of them truncates many
  # more than 1% of the weights that vary.
  data <- simulateCloneCensorWeightData(
    n = 5000, nCovariates = 5, dgp = "a", graceEnd = 1, maxDays = 30, initiationHazard = 0.03,
    outcomeHazard = 0.005, initiationCoefficients = rep(0.6, 5), outcomeCoefficients = rep(0.4, 5), seed = 1
  )
  on.exit(close(data))
  population <- suppressMessages(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
    riskWindowStart = 0, riskWindowEnd = 30, endAnchor = "cohort start"
  )))
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = list(
    createTreatmentStrategy(1, "Initiate by day 1", "initiateWithin", graceEnd = 1),
    createTreatmentStrategy(2, "Not by day 1", "notInitiateWithin", graceEnd = 1)
  ))))
  on.exit(close(clones), add = TRUE)
  models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(prior = Cyclops::createPrior("none"))))

  untruncated <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(truncation = createWeightTruncation("none"))))
  raw <- collectSegments(untruncated)
  personWeights <- as.data.frame(dplyr::collect(untruncated$personWeights))
  cloneTable <- collectClones(untruncated)
  expect_equal(attr(untruncated, "metaData")$weights$caps$truncatedShare, c(0, 0))
  close(untruncated)
  personWeights$varying <- cloneTable$endDay[match(personWeights$cloneRowId, cloneTable$cloneRowId)] >= 2

  truncate <- function(basis) {
    weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(truncation = createWeightTruncation("percentile", 0.99, basis = basis))))
    on.exit(close(weighted))
    return(list(segments = collectSegments(weighted), caps = attr(weighted, "metaData")$weights$caps))
  }
  varying <- truncate("varying")
  allClones <- truncate("allClones")
  for (strategyId in 1:2) {
    strategy <- personWeights$strategyId == strategyId
    # "allClones" reproduces the original caps exactly:
    expect_equal(allClones$caps$cap[strategyId], unname(quantile(personWeights$maxWeight[strategy], 0.99, type = 7)))
    # "varying" takes the percentile over the clones under follow-up on day tau + 1 = 2:
    weights <- personWeights$maxWeight[strategy & personWeights$varying]
    cap <- varying$caps$cap[strategyId]
    expect_equal(cap, unname(quantile(weights, 0.99, type = 7)))
    # (quantile type 7 interpolates between order statistics, so up to one more clone can be above it)
    expect_lte(mean(weights > cap), 0.01 + 1 / length(weights))
    expect_gte(mean(weights > cap), 0.01 - 1 / length(weights))
    for (result in list(varying, allClones)) {
      cap <- result$caps$cap[strategyId]
      s <- raw$strategyId == strategyId
      expect_equal(result$segments$weight[s], pmin(raw$weight[s], cap))
      expect_equal(result$caps$truncatedShare[strategyId], mean(personWeights$maxWeight[strategy] > cap))
    }
  }
  # Over all clones, the cap of initiateWithin(1) truncates far more than 1% of the weights that vary:
  weights <- personWeights$maxWeight[personWeights$strategyId == 1 & personWeights$varying]
  expect_gt(mean(weights > allClones$caps$cap[1]), 0.05)
  expect_gt(varying$caps$cap[1], allClones$caps$cap[1])
  expect_gt(varying$caps$truncatedShare[1], 0)
  expect_lt(varying$caps$truncatedShare[1], allClones$caps$truncatedShare[1])
  # The fixed rule reports its share too (the same cap on the other scale, set explicitly so it does not warn):
  weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(truncation = createWeightTruncation("fixed", 2, otherScaleUpper = 2))))
  caps <- attr(weighted, "metaData")$weights$caps
  close(weighted)
  expect_equal(caps$cap, c(2, 2))
  expect_equal(caps$truncatedShare, as.vector(tapply(personWeights$maxWeight > 2, personWeights$strategyId, mean)))
})

test_that("The truncated share does not depend on how the clones are cut", {
  # Review of contract 11.5: a share of segments changed with the cuts of outcomeTimeBasis (0.11 vs
  # 0.21 for the same weights and cap). The share of clones does not.
  data <- simulateCloneCensorWeightData(
    n = 1500, nCovariates = 1, dgp = "a", graceEnd = 7, maxDays = 40, initiationHazard = 0.02,
    outcomeHazard = 0.004, initiationCoefficients = 1, outcomeCoefficients = 1, seed = 21
  )
  on.exit(close(data))
  population <- suppressMessages(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
    riskWindowStart = 0, riskWindowEnd = 40, endAnchor = "cohort start"
  )))
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = list(
    createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
    createTreatmentStrategy(2, "Do not initiate within 7 days", "notInitiateWithin", graceEnd = 7)
  ))))
  on.exit(close(clones), add = TRUE)
  models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(prior = Cyclops::createPrior("none"))))
  for (stabilized in c(FALSE, TRUE)) {
    for (truncation in list(createWeightTruncation("fixed", 1.5, otherScaleUpper = 1.5), createWeightTruncation("percentile", 0.9))) {
      # (The default tier runs the fixed cap unstabilized and the percentile stabilized.)
      if (!isFullTestRun() && stabilized != (truncation$type == "percentile")) next
      args <- createComputeWeightsArgs(stabilized = stabilized, truncation = truncation)
      coarse <- suppressMessages(computeCensoringWeights(clones, models, args))
      daily <- suppressMessages(computeCensoringWeights(clones, models, args, outcomeTimeBasis = createTimeBasis(intervalLength = 1)))
      label <- sprintf("stabilized = %s, %s", stabilized, truncation$type)
      coarseCaps <- attr(coarse, "metaData")$weights$caps
      dailyCaps <- attr(daily, "metaData")$weights$caps
      expect_gt(nrow(collectSegments(daily)), nrow(collectSegments(coarse)))
      expect_equal(dailyCaps$cap, coarseCaps$cap, label = label)
      expect_equal(dailyCaps$truncatedShare, coarseCaps$truncatedShare, label = label)
      personWeights <- as.data.frame(dplyr::collect(coarse$personWeights))
      expect_equal(coarseCaps$truncatedShare, as.vector(tapply(personWeights$maxWeight > coarseCaps$cap[personWeights$strategyId], personWeights$strategyId, mean)), label = label)
      close(coarse)
      close(daily)
    }
  }
  # The fixed cap truncates some clones of both strategies (unstabilized weights):
  weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = FALSE, truncation = createWeightTruncation("fixed", 1.5))))
  expect_true(all(attr(weighted, "metaData")$weights$caps$truncatedShare > 0))
  close(weighted)
})

test_that("Without clones under follow-up on day tau + 1, the percentile is taken over all clones", {
  fixture <- createM3Fixture(150, riskStart = 0L, seed = 23)
  # Follow-up ends by day 5, before tau + 1 = 8:
  fixture$persons$riskEnd <- pmin(fixture$persons$riskEnd, 5L)
  fixture$persons$outcomeDay[!is.na(fixture$persons$outcomeDay) & fixture$persons$outcomeDay > 5] <- NA
  fixture$persons$deathDay[!is.na(fixture$persons$deathDay) & fixture$persons$deathDay > 5] <- NA
  data <- createM3Data(fixture)
  on.exit(close(data))
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "neverInitiate", graceEnd = 7), "eventsFirst")
  on.exit(close(clones), add = TRUE)
  models <- suppressMessages(fitCensoringModels(clones, data, m3Args("cox")))
  untruncated <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = FALSE, truncation = createWeightTruncation("none"))))
  personWeights <- as.data.frame(dplyr::collect(untruncated$personWeights))
  close(untruncated)
  weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = FALSE, truncation = createWeightTruncation("percentile", 0.9))))
  on.exit(close(weighted), add = TRUE)
  expected <- as.vector(tapply(personWeights$maxWeight, personWeights$strategyId, quantile, probs = 0.9, type = 7))
  expect_equal(attr(weighted, "metaData")$weights$caps$cap, expected)
})

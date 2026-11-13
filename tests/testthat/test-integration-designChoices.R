library(testthat)

# Contract section 12: the defaults postGraceModel = "hazard" and covariateTimeInteraction = "day0", the
# truncation of the MSM-scale weights (truncateMsmWeights), every alternative still reachable, and the
# diagnostics that help choose between them (section 12.5).

designStrategies <- function(second = "notInitiateWithin", graceEnd = 7) {
  list(
    createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = graceEnd),
    createTreatmentStrategy(2, "Comparator", second, graceEnd = if (second == "neverInitiate") NULL else graceEnd)
  )
}

designClones <- function(data, second = "notInitiateWithin", sameDayOrdering = "eventsFirst") {
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowStart = 0, riskWindowEnd = 60, endAnchor = "cohort start")
  ))
  return(suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = designStrategies(second), sameDayOrdering = sameDayOrdering))))
}

# A fixed prior: unpenalized fits of the covariate-by-day-0 interactions can separate in small samples.
designArgs <- function(..., prior = Cyclops::createPrior("laplace", exclude = c(0), variance = 1)) {
  createFitCensoringModelsArgs(prior = prior, ...)
}

# The grace-period DGP of the statistical review in which starting the treatment cuts the outcome hazard
# by 80% from the next day ("bStrong"), and its null version.
selectionData <- function(treatmentEffect, n = 5000, seed = 1) {
  return(simulateCloneCensorWeightData(
    n = n, seed = seed, nCovariates = 1, dgp = "b", graceEnd = 7, maxDays = 60, binaryPrevalence = 0.4,
    initiationHazard = 0.08, outcomeHazard = 0.03, deathHazard = 0.002, observationEndHazard = 0.003,
    initiationCoefficients = 1.5, outcomeCoefficients = 1.5, deathCoefficients = 0.5,
    treatmentEffect = treatmentEffect, effectRampDays = 1
  ))
}

test_that("The new defaults fit day-0 interactions and hazard weights after G; with longer intervals day0 is the first interval", {
  data <- simulateCloneCensorWeightData(n = 300, maxDays = 60, graceEnd = 7, seed = 2)
  clones <- designClones(data)
  on.exit({
    close(clones)
    close(data)
  })
  models <- suppressMessages(fitCensoringModels(clones, data, designArgs()))
  expect_null(models[[2]]$postGrace)
  expect_equal(models[[2]]$interactions$type, "day0")
  expect_true(any(grepl(":day0$", names(models[[2]]$coefficients))))
  # Every alternative stays reachable:
  binary <- suppressMessages(fitCensoringModels(clones, data, designArgs(postGraceModel = "binary", covariateTimeInteraction = "none")))
  expect_false(is.null(binary[[2]]$postGrace))
  expect_false(any(grepl(":", names(binary[[2]]$coefficients))))
  interval <- suppressMessages(fitCensoringModels(clones, data, designArgs(covariateTimeInteraction = "interval")))
  expect_true(any(grepl(":timeBasis", names(interval[[2]]$coefficients))))
  # Weekly intervals: "day0" multiplies the indicator of the first interval (days 0-6):
  weekly <- suppressMessages(fitCensoringModels(clones, data, designArgs(timeBasis = createTimeBasis(intervalLength = 7))))
  expect_equal(unname(weekly[[2]]$interactions$basis[, "day0"]), c(1, 0))
  expect_true(any(grepl(":day0$", names(weekly[[2]]$coefficients))))
  weights <- suppressMessages(computeCensoringWeights(clones, weekly, createComputeWeightsArgs()))
  expect_true(all(is.finite(collectSegments(weights)$weight)))
  close(weights)
})

test_that("truncateMsmWeights = TRUE caps the other-scale weights on their own scale; FALSE reproduces section 11.1", {
  data <- simulateCloneCensorWeightData(n = 500, maxDays = 60, graceEnd = 7, seed = 3)
  clones <- designClones(data, second = "neverInitiate")
  on.exit({
    close(clones)
    close(data)
  })
  models <- suppressMessages(fitCensoringModels(clones, data, designArgs()))
  cloneTable <- collectClones(clones)
  for (stabilized in c(TRUE, FALSE)) {
    label <- sprintf("stabilized = %s", stabilized)
    untruncated <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = stabilized, truncation = createWeightTruncation("none"))))
    u <- collectSegments(untruncated)
    close(untruncated)
    # Without truncation the other-scale weight is the weight divided (or multiplied) by the numerator:
    other <- if (stabilized) u$weight / u$stabilizationNumerator else u$weight * u$stabilizationNumerator
    expect_equal(u$otherScaleWeight, other, tolerance = 1e-10, info = label)

    on <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = stabilized)))
    off <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = stabilized, truncateMsmWeights = FALSE)))
    capsOn <- attr(on, "metaData")$weights$caps
    capsOff <- attr(off, "metaData")$weights$caps
    expect_true(attr(on, "metaData")$weights$truncateMsmWeights, info = label)
    expect_false(attr(off, "metaData")$weights$truncateMsmWeights, info = label)
    sOn <- collectSegments(on)
    sOff <- collectSegments(off)
    # The weights of computeCensoringWeights() do not depend on truncateMsmWeights:
    expect_equal(sOn$weight, sOff$weight, info = label)
    expect_equal(capsOn$cap, capsOff$cap, info = label)
    # TRUE: the 99th percentile of the clones' largest other-scale weights over the clones under follow-up
    # on day tau + 1 = 8 (survivalTime >= 9), and each other-scale weight capped there:
    maxOther <- tapply(u$otherScaleWeight, u$cloneRowId, max)
    survivalTime <- cloneTable$survivalTime[match(as.integer(names(maxOther)), cloneTable$cloneRowId)]
    strategyId <- cloneTable$strategyId[match(as.integer(names(maxOther)), cloneTable$cloneRowId)]
    for (s in 1:2) {
      varying <- strategyId == s & survivalTime >= 9
      expectedCap <- unname(stats::quantile(maxOther[varying], 0.99, type = 7))
      expect_equal(capsOn$otherScaleCap[capsOn$strategyId == s], expectedCap, tolerance = 1e-10, info = label)
      inStrategy <- u$strategyId == s
      expect_equal(sOn$otherScaleWeight[inStrategy], pmin(u$otherScaleWeight[inStrategy], expectedCap), tolerance = 1e-10, info = label)
      expect_equal(capsOn$otherScaleTruncatedShare[capsOn$strategyId == s], mean(maxOther[strategyId == s] > expectedCap), info = label)
    }
    # FALSE: formed from the truncated weights, with no cap of its own:
    expectedOff <- if (stabilized) sOff$weight / sOff$stabilizationNumerator else sOff$weight * sOff$stabilizationNumerator
    expect_equal(sOff$otherScaleWeight, expectedOff, tolerance = 1e-10, info = label)
    expect_true(all(is.na(capsOff$otherScaleCap)), info = label)
    expect_equal(capsOff$otherScaleTruncatedShare, capsOff$truncatedShare, info = label)

    # The outcome model uses and reports the weights of msmWeights:
    for (msmWeights in c("unstabilized", "stabilized")) {
      info <- paste(label, msmWeights)
      args <- createFitOutcomeModelArgs(horizons = 30, msmWeights = msmWeights)
      fitOn <- suppressMessages(fitOutcomeModel(on, args))
      fitOff <- suppressMessages(fitOutcomeModel(off, args))
      sameScale <- (msmWeights == "stabilized") == stabilized
      summaryOn <- fitOn$metaData$msmWeightSummary
      summaryOff <- fitOff$metaData$msmWeightSummary
      expect_equal(summaryOn$scale, rep(msmWeights, 2), info = info)
      usedOn <- if (sameScale) sOn$weight else sOn$otherScaleWeight
      usedOff <- if (sameScale) sOff$weight else sOff$otherScaleWeight
      expect_equal(summaryOn$maxWeight, as.vector(tapply(usedOn, sOn$strategyId, max)), info = info)
      expect_equal(summaryOff$maxWeight, as.vector(tapply(usedOff, sOff$strategyId, max)), info = info)
      expect_equal(summaryOn$cap, if (sameScale) capsOn$cap else capsOn$otherScaleCap, info = info)
      expect_equal(summaryOn$truncatedShare, if (sameScale) capsOn$truncatedShare else capsOn$otherScaleTruncatedShare, info = info)
      expect_equal(summaryOff$cap, if (sameScale) capsOff$cap else rep(NA_real_, 2), info = info)
      expect_true(all(summaryOn$maxWeight <= summaryOn$cap + 1e-12), info = info)
      if (sameScale) {
        expect_equal(fitOn$estimates, fitOff$estimates, info = info)
      }
      # The Cox MSM is fitted with exactly those weights:
      outcomeData <- .withMsmWeights(.getOutcomeData(on), msmWeights)
      expect_equal(outcomeData$segments$weight, usedOn, info = info)
      # Frequency weights multiply the other-scale weights too:
      counts <- data.frame(personSeqId = unique(sOn$personSeqId), count = rep(c(0, 2), length.out = length(unique(sOn$personSeqId))))
      counted <- .withMsmWeights(.applyPersonCounts(.getOutcomeData(on), counts), msmWeights)
      expect_equal(counted$segments$weight, usedOn * counts$count[match(counted$segments$personSeqId, counts$personSeqId)], info = info)
    }
    expect_null(suppressMessages(fitOutcomeModel(on, createFitOutcomeModelArgs(modelType = "km", horizons = 30)))$metaData$msmWeightSummary)
    close(on)
    close(off)
  }
})

test_that("Calibration reports groups, calibration in the large, covariate strata and the selection check", {
  data <- simulateCloneCensorWeightData(n = 400, maxDays = 60, graceEnd = 7, seed = 4)
  clones <- designClones(data)
  on.exit({
    close(clones)
    close(data)
  })
  models <- suppressMessages(fitCensoringModels(clones, data, designArgs(postGraceModel = "binary")))
  calibration <- computeCensoringModelCalibration(clones, models, numberOfGroups = 5, cloneCensorWeightData = data, strataCovariateIds = c(1001, 2002))
  expect_s3_class(calibration, "CcwCensoringModelCalibration")
  expect_equal(colnames(calibration), c("strategyId", "model", "covariateTimeInteraction", "grouping", "covariateId", "group", "persons", "predicted", "observed"))
  expect_equal(unique(calibration$covariateTimeInteraction[calibration$model == "hazard"]), "day0")
  expect_true(all(is.na(calibration$covariateTimeInteraction[calibration$model == "binary"])))

  persons <- dplyr::collect(clones$persons)
  riskSet <- .getBinaryRiskSet(persons, persons$treatmentDay, 7L, "eventsFirst")
  hazard <- 1 - .computeHazardSurvival(models[[2]], riskSet$rowId, 7L)
  rows <- calibration[calibration$strategyId == 2 & calibration$model == "hazard", ]
  overall <- rows[rows$grouping == "overall", ]
  expect_equal(overall$persons, nrow(riskSet))
  expect_equal(overall$observed, mean(riskSet$y))
  expect_equal(overall$predicted, mean(hazard))
  covariates <- dplyr::collect(data$covariates)
  has1001 <- riskSet$rowId %in% covariates$rowId[covariates$covariateId == 1001 & covariates$covariateValue != 0]
  strata <- rows[rows$grouping == "covariate" & rows$covariateId == 1001, ]
  expect_equal(strata$group, c(0L, 1L))
  expect_equal(strata$persons, c(sum(!has1001), sum(has1001)))
  expect_equal(strata$observed, c(mean(riskSet$y[!has1001]), mean(riskSet$y[has1001])))
  expect_equal(strata$predicted, c(mean(hazard[!has1001]), mean(hazard[has1001])))
  groups <- rows[rows$grouping == "prediction", ]
  expect_equal(sum(groups$persons), nrow(riskSet))
  expect_equal(groups$predicted, as.vector(tapply(hazard, .groupByRank(hazard, 5), mean)))
  # The binary model reproduces the observed share in the large:
  binaryOverall <- calibration[calibration$model == "binary" & calibration$grouping == "overall", ]
  expect_equal(binaryOverall$predicted, binaryOverall$observed, tolerance = 1e-6)
  expect_s3_class(attr(calibration, "gracePeriodSelection"), "CcwGracePeriodSelection")
  expect_output(print(calibration), "Calibration in the large")
  expect_error(computeCensoringModelCalibration(clones, models, strataCovariateIds = 1001), "cloneCensorWeightData")
  skip_if_not_installed("ggplot2")
  expect_s3_class(plotCensoringModelCalibration(calibration), "ggplot")
})

test_that("Calibration skips a strategy with nobody under follow-up at G (empty risk set guard)", {
  # Everybody's follow-up ends before G = 12, so the binary risk set at G is empty:
  persons <- data.frame(
    rowId = 1:6,
    personSeqId = 1:6,
    riskStart = 0L,
    riskEnd = c(10L, 9L, 10L, 6L, 10L, 8L),
    nominalRiskEnd = 10L,
    outcomeDay = c(NA, NA, 1L, NA, NA, NA),
    deathDay = NA_integer_,
    treatmentDay = c(0L, 1L, NA, 2L, NA, 4L),
    comparatorTreatmentDay = NA_integer_
  )
  data <- createM3Data(list(persons = persons, covariates = data.frame(rowId = integer(0), covariateId = numeric(0), covariateValue = numeric(0))))
  clones <- buildM3Clones(persons, m3Strategies("neverInitiate", "notInitiateWithin", graceEnd = 12))
  on.exit({
    close(clones)
    close(data)
  })
  models <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic", useCovariates = FALSE)))
  calibration <- computeCensoringModelCalibration(clones, models)
  expect_equal(nrow(calibration), 0)
  expect_equal(nrow(computeGracePeriodSelection(clones, models)[0, ]), 0)
})

test_that("The grace-period selection check counts initiated and not-initiated person-time by hand", {
  # G = 2. Person 8 starts on day 0 and has the outcome on day 2; person 9 starts on day 1 and dies on
  # day 1. With eventsFirst, follow-up starts on day 1 and initiated time on the day after initiation:
  # initiated days 2 (person 1) + 1 (2) + 2 (8) = 5; not-initiated days 0 + 1 + 1 (3: outcome day 1) + 2
  # + 2 + 2 + 2 (7: outcome day 2) + 0 + 1 (9: death day 1, before its initiation counts) = 11.
  # With deviationFirst, from day 0 and the day of initiation: initiated 3 + 2 + 1 (4) + 3 + 1 (9) = 10,
  # not initiated 0 + 1 + 2 + 2 + 3 + 3 + 3 + 0 + 1 = 15, and the death of person 9 is initiated.
  persons <- data.frame(
    rowId = 1:9,
    personSeqId = 1:9,
    riskStart = 0L,
    riskEnd = c(10L, 10L, 10L, 6L, 10L, 10L, 10L, 10L, 1L),
    nominalRiskEnd = 10L,
    outcomeDay = c(NA, NA, 1L, NA, NA, NA, 2L, 2L, NA),
    deathDay = c(NA, NA, NA, NA, NA, 8L, NA, NA, 1L),
    treatmentDay = c(0L, 1L, NA, 2L, NA, 4L, NA, 0L, 1L),
    comparatorTreatmentDay = NA_integer_
  )
  expected <- list(
    eventsFirst = c(persons = 9, initiators = 5, initiatedDays = 5, notInitiatedDays = 11, initiatedOutcomes = 1, notInitiatedOutcomes = 2, initiatedDeaths = 0, notInitiatedDeaths = 1),
    deviationFirst = c(persons = 9, initiators = 5, initiatedDays = 10, notInitiatedDays = 15, initiatedOutcomes = 1, notInitiatedOutcomes = 2, initiatedDeaths = 1, notInitiatedDeaths = 0)
  )
  for (sameDayOrdering in names(expected)) {
    clones <- buildM3Clones(persons, m3Strategies(graceEnd = 2), sameDayOrdering)
    selection <- computeGracePeriodSelection(clones)
    expect_s3_class(selection, "CcwGracePeriodSelection")
    expect_equal(selection$grouping, c("overall", "overall"))
    row <- unlist(selection[1, names(expected[[sameDayOrdering]])])
    expect_equal(row, expected[[sameDayOrdering]], info = sameDayOrdering)
    e <- expected[[sameDayOrdering]]
    rateRatio <- ((e[["initiatedOutcomes"]] + e[["initiatedDeaths"]]) / e[["initiatedDays"]]) / ((e[["notInitiatedOutcomes"]] + e[["notInitiatedDeaths"]]) / e[["notInitiatedDays"]])
    expect_equal(selection$rateRatio[1], rateRatio, info = sameDayOrdering)
    summary <- attr(selection, "summary")
    expect_false(any(summary$stratified))
    expect_true(all(summary$byDay))
    # Without the day strata the summary is the crude rate ratio:
    crude <- attr(computeGracePeriodSelection(clones, stratifyByDay = FALSE), "summary")
    expect_false(any(crude$byDay))
    expect_equal(crude$rateRatio[crude$strategyId == 1 & crude$event == "outcomeOrDeath"], rateRatio, info = sameDayOrdering)
    close(clones)
    if (sameDayOrdering == "eventsFirst") {
      # By day: day 1 has initiated days 2 (persons 1, 8) and not-initiated days 7 with 2 events (3, 9);
      # day 2 has initiated days 3 (1, 2, 8) with 1 event (8) and not-initiated days 4 with 1 event (7).
      # Mantel-Haenszel: (0 + 1 * 4 / 7) / (2 * 2 / 9 + 1 * 3 / 7) = 36 / 55.
      expect_equal(summary$rateRatio[summary$strategyId == 1 & summary$event == "outcomeOrDeath"], 36 / 55)
    }
  }
  # The Mantel-Haenszel rate ratio of two strata with equal rates within strata is 1, while the crude
  # one is not:
  mh <- .mantelHaenszelRateRatio(a = c(10, 40), b = c(100, 40), t1 = c(100, 100), t0 = c(1000, 100))
  expect_equal(mh$rateRatio, 1)
  expect_gt(mh$p, 0.5)
  expect_output(print(selection), "Grace-period selection check")
})

test_that("bStrong: the selection check flags, binary is biased and hazard is not; the null DGP is not flagged", {
  strategies <- designStrategies()
  results <- list()
  for (effect in c("bStrong", "null")) {
    data <- selectionData(if (effect == "bStrong") log(0.2) else 0)
    settings <- attr(data, "metaData")$simulationSettings
    clones <- designClones(data)
    risks <- c()
    for (postGraceModel in c("hazard", "binary")) {
      models <- suppressMessages(fitCensoringModels(clones, data, designArgs(postGraceModel = postGraceModel, prior = Cyclops::createPrior("none"))))
      if (postGraceModel == "hazard") {
        summary <- attr(computeGracePeriodSelection(clones, models), "summary")
        results[[effect]] <- summary[summary$strategyId == 2 & summary$event == "outcomeOrDeath", ]
      }
      weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(truncation = createWeightTruncation("none"))))
      fit <- suppressMessages(fitOutcomeModel(weighted, createFitOutcomeModelArgs(modelType = "km", horizons = 60)))
      close(weighted)
      row <- fit$estimates[fit$estimates$estimand == "risk_comparator_60", ]
      risks[postGraceModel] <- row$estimate
      if (postGraceModel == "hazard") {
        se <- row$se
      }
    }
    close(clones)
    close(data)
    if (effect == "bStrong") {
      truth <- computeTrueRisks(settings, strategies, horizons = 60, nTruth = 2e5, seed = 1)
      truth <- truth$risk[truth$strategyId == 2]
      # Binary over-weights the selected survivors (review: bias +0.020); hazard is unbiased:
      expect_lt(abs(risks[["hazard"]] - truth), 2 * se)
      expect_gt(risks[["binary"]] - risks[["hazard"]], 0.01)
      expect_gt(abs(risks[["binary"]] - truth), abs(risks[["hazard"]] - truth))
    }
  }
  expect_true(results$bStrong$flag)
  expect_true(results$bStrong$stratified)
  expect_lt(results$bStrong$ub95, 0.5)
  expect_false(results$null$flag)
  expect_true(results$null$lb95 < 1 && results$null$ub95 > 1)
})

test_that("A baseline risk that changes over the grace period does not raise the selection flag with day strata", {
  # The null DGP (starting has no effect) with a falling baseline outcome risk: 75% of the outcomes after
  # day 3 removed at random, whoever started. Initiated person-time falls later in the grace period, so
  # without the day strata it has the lower risk (review of section 12: rate ratio 0.63, p = 1e-109 at
  # n = 50,000).
  data <- selectionData(0)
  outcomes <- data$outcomes |>
    collect()
  set.seed(1)
  data$outcomes <- outcomes[!(outcomes$daysToEvent > 3 & runif(nrow(outcomes)) < 0.75), ]
  clones <- designClones(data)
  on.exit({
    close(clones)
    close(data)
  })
  models <- suppressMessages(fitCensoringModels(clones, data, designArgs(prior = Cyclops::createPrior("none"))))
  pick <- function(selection) {
    summary <- attr(selection, "summary")
    return(summary[summary$strategyId == 2 & summary$event == "outcomeOrDeath", ])
  }
  byDay <- pick(computeGracePeriodSelection(clones, models))
  expect_true(byDay$byDay)
  expect_false(byDay$flag)
  expect_true(byDay$lb95 < 1 && byDay$ub95 > 1)
  byGroup <- pick(computeGracePeriodSelection(clones, models, stratifyByDay = FALSE))
  expect_true(byGroup$flag)
  expect_lt(byGroup$ub95, 0.8)
})

test_that("The selection flag advises on postGraceModel only for notInitiateWithin strategies", {
  selection <- dplyr::tibble(strategyId = 1:2)
  attr(selection, "summary") <- dplyr::tibble(
    strategyId = 1:2, event = "outcomeOrDeath", stratified = TRUE, byDay = TRUE,
    rateRatio = 0.3, lb95 = 0.25, ub95 = 0.35, p = 1e-20, flag = TRUE
  )
  attr(selection, "strategyTypes") <- c("initiateWithin", "notInitiateWithin")
  output <- capture.output(.printSelectionFlags(selection))
  expect_match(output[2], "Strategy 1: stratified by group and day 0.30.*FLAG")
  expect_no_match(output[2], "postGraceModel")
  expect_match(output[3], "Strategy 2: .*FLAG.*postGraceModel = \"binary\" is likely biased")
})

test_that("Positivity summarizes the predicted probabilities of every model", {
  data <- simulateCloneCensorWeightData(n = 300, maxDays = 60, graceEnd = 7, seed = 5)
  clones <- designClones(data)
  on.exit({
    close(clones)
    close(data)
  })
  models <- suppressMessages(fitCensoringModels(clones, data, designArgs(postGraceModel = "binary")))
  positivity <- computePositivity(clones, models, threshold = 0.4)
  expect_s3_class(positivity, "CcwPositivity")
  expect_equal(positivity$model, c("binary", "hazard", "binary"))
  expect_equal(positivity$day, c(7L, 7L, 7L))
  pi <- models[[1]]$predictions$propensityScore
  expect_equal(positivity$median[1], stats::median(pi))
  expect_equal(positivity$p99[1], unname(stats::quantile(pi, 0.99)))
  expect_equal(positivity$shareExtreme[1], mean(pi < 0.4))
  hazard <- 1 - .computeHazardSurvival(models[[2]], models[[2]]$predictions$rowId, 7L)
  expect_equal(positivity$mean[2], mean(hazard))
  expect_equal(positivity$persons[2], length(hazard))
  expect_equal(positivity$shareExtreme[2], mean(1 - hazard < 0.4))
  expect_equal(positivity$shareAbove[2], mean(hazard > 0.6))
  # The binary model of the weights after G of strategy 2 is that of strategy 1:
  columns <- c("day", "persons", "mean", "min", "median", "max", "shareBelow", "shareAbove")
  expect_equal(positivity[3, columns], positivity[1, columns])
  expect_equal(positivity$shareExtreme[3], mean(1 - pi < 0.4))
  expect_equal(nrow(attr(positivity, "predictions")), length(pi) * 2 + length(hazard))
  expect_output(print(positivity), "Predicted probability")
  skip_if_not_installed("ggplot2")
  expect_s3_class(plotPositivity(positivity), "ggplot")
})

test_that("The weight summary matches the daily weights on both scales, before and after truncation", {
  data <- simulateCloneCensorWeightData(n = 400, maxDays = 60, graceEnd = 7, seed = 6)
  clones <- designClones(data, second = "neverInitiate")
  on.exit({
    close(clones)
    close(data)
  })
  # Without interactions, so refDayWeights() gives the exact daily weights:
  models <- suppressMessages(fitCensoringModels(clones, data, designArgs(covariateTimeInteraction = "none")))
  days <- c(3, 8, 20, 40)
  reference <- refDayWeights(models, collectClones(clones), days, "eventsFirst")
  for (truncateMsmWeights in c(TRUE, FALSE)) {
    weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(truncation = createWeightTruncation("percentile", 0.9), truncateMsmWeights = truncateMsmWeights)))
    caps <- attr(weighted, "metaData")$weights$caps
    summary <- computeWeightSummary(weighted, days = days, dayBandStarts = c(0, 8, 9))
    expect_s3_class(summary, "CcwWeightSummary")
    expect_equal(sort(unique(summary$dayBand)), c("0-7", "8", "9-40"))
    for (s in 1:2) {
      cap <- caps$cap[caps$strategyId == s]
      otherCap <- caps$otherScaleCap[caps$strategyId == s]
      for (band in list(c(0, 7), c(8, 8), c(9, 40))) {
        r <- reference[reference$strategyId == s & reference$day >= band[1] & reference$day <= band[2], ]
        stabilizedBefore <- r$weight * r$numerator
        expected <- list(
          stabilized.before = stabilizedBefore,
          stabilized.after = pmin(stabilizedBefore, cap),
          unstabilized.before = r$weight,
          unstabilized.after = if (truncateMsmWeights) pmin(r$weight, otherCap) else pmin(stabilizedBefore, cap) / r$numerator
        )
        for (key in names(expected)) {
          parts <- strsplit(key, ".", fixed = TRUE)[[1]]
          row <- summary[summary$strategyId == s & summary$scale == parts[1] & summary$truncation == parts[2] & summary$bandStart == band[1], ]
          w <- expected[[key]]
          info <- sprintf("truncateMsmWeights = %s, strategy %d, band %d, %s", truncateMsmWeights, s, band[1], key)
          expect_equal(row$cloneDays, length(w), info = info)
          expect_equal(row$meanWeight, mean(w), tolerance = 1e-6, info = info)
          expect_equal(row$maxWeight, max(w), tolerance = 1e-6, info = info)
          expect_equal(row$p99Weight, unname(stats::quantile(w, 0.99)), tolerance = 1e-6, info = info)
          ess <- vapply(split(w, r$day), function(x) sum(x)^2 / sum(x^2), numeric(1))
          expect_equal(row$effectiveSampleSize, mean(ess), tolerance = 1e-6, info = info)
          if (parts[2] == "after") {
            before <- expected[[paste0(parts[1], ".before")]]
            expect_equal(row$shareTruncated, mean(before > w * (1 + 1e-9)), info = info)
            expectedCap <- if (parts[1] == "stabilized") cap else if (truncateMsmWeights) otherCap else NA_real_
            expect_equal(row$cap, expectedCap, info = info)
          } else {
            expect_true(is.na(row$shareTruncated) && is.na(row$cap), info = info)
          }
        }
      }
    }
    # The weights of computeCensoringWeights() after truncation agree with computeWeightDistribution():
    distribution <- computeWeightDistribution(weighted, days = 20)
    row <- summary[summary$scale == "stabilized" & summary$truncation == "after" & summary$bandStart == 9, ]
    single <- computeWeightSummary(weighted, days = 20, dayBandStarts = 9)
    single <- single[single$scale == "stabilized" & single$truncation == "after", ]
    expect_equal(single$meanWeight, distribution$meanWeight, tolerance = 1e-10)
    expect_equal(single$effectiveSampleSize, distribution$effectiveSampleSize, tolerance = 1e-10)
    close(weighted)
  }
  expect_output(print(summary), "Weights per strategy")
  # The default bands: the grace period, tau + 1, and later days:
  weighted <- suppressMessages(computeCensoringWeights(clones, models))
  defaults <- computeWeightSummary(weighted)
  close(weighted)
  expect_equal(unique(defaults$bandStart), c(0L, 8L, 9L))
  expect_equal(unique(defaults$dayBand)[1:2], c("0-7", "8"))
  skip_if_not_installed("ggplot2")
  expect_s3_class(plotWeightSummary(defaults), "ggplot")
})

test_that("A fixed cap on the other scale is the same number unless otherScaleUpper is set, with a warning", {
  data <- simulateCloneCensorWeightData(n = 2000, maxDays = 60, graceEnd = 7, seed = 3, initiationHazard = 0.02)
  clones <- designClones(data)
  on.exit({
    close(clones)
    close(data)
  })
  models <- suppressMessages(fitCensoringModels(clones, data, designArgs()))
  untruncated <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(truncation = createWeightTruncation("none"))))
  u <- collectSegments(untruncated)
  close(untruncated)
  # The same cap 5 on both scales truncates many more unstabilized weights, and warns:
  expect_warning(
    same <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(truncation = createWeightTruncation("fixed", 5)))),
    "otherScaleUpper"
  )
  caps <- attr(same, "metaData")$weights$caps
  expect_equal(caps$otherScaleCap, c(5, 5))
  expect_gt(max(caps$otherScaleTruncatedShare - caps$truncatedShare), 0.01)
  s <- collectSegments(same)
  close(same)
  expect_equal(s$otherScaleWeight, pmin(u$otherScaleWeight, 5))
  expect_equal(s$weight, pmin(u$weight, 5))
  # otherScaleUpper sets the cap on the other scale only:
  expect_no_warning(
    own <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(truncation = createWeightTruncation("fixed", 5, otherScaleUpper = 50))))
  )
  caps <- attr(own, "metaData")$weights$caps
  expect_equal(caps$cap, c(5, 5))
  expect_equal(caps$otherScaleCap, c(50, 50))
  s <- collectSegments(own)
  close(own)
  expect_equal(s$otherScaleWeight, pmin(u$otherScaleWeight, 50))
  expect_equal(s$weight, pmin(u$weight, 5))
})

test_that("The default sensitivity variants leave out truncateMsmWeights for Kaplan-Meier", {
  variants <- function(modelType) {
    names(.defaultSensitivityVariants(designStrategies(), designArgs(), createComputeWeightsArgs(), createFitOutcomeModelArgs(modelType = modelType)))
  }
  expect_false("truncateMsmWeights = FALSE" %in% variants("km"))
  expect_true("truncateMsmWeights = FALSE" %in% variants("cox"))
})

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
    "specified", "postGraceModel = binary", "covariateTimeInteraction = none", "covariateTimeInteraction = interval",
    "truncateMsmWeights = FALSE", "basis = allClones", "msmWeights = stabilized"
  ))
  expect_equal(variants$postGraceModel, c("hazard", "binary", rep("hazard", 5)))
  expect_equal(variants$covariateTimeInteraction, c("day0", "day0", "none", "interval", "day0", "day0", "day0"))
  expect_equal(variants$truncateMsmWeights, c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE))
  expect_equal(variants$truncationBasis, c(rep("varying", 5), "allClones", "varying"))
  expect_equal(variants$msmWeights, c(rep("unstabilized", 6), "stabilized"))
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
  expect_equal(pick("specified"), as.data.frame(direct(fitCensoringModelsArgs, createComputeWeightsArgs(), fitOutcomeModelArgs)))
  expect_equal(pick("postGraceModel = binary"), as.data.frame(direct(designArgs(postGraceModel = "binary"), createComputeWeightsArgs(), fitOutcomeModelArgs)))
  expect_equal(pick("truncateMsmWeights = FALSE"), as.data.frame(direct(fitCensoringModelsArgs, createComputeWeightsArgs(truncateMsmWeights = FALSE), fitOutcomeModelArgs)))
  expect_equal(pick("msmWeights = stabilized"), as.data.frame(direct(fitCensoringModelsArgs, createComputeWeightsArgs(), createFitOutcomeModelArgs(horizons = 30, msmWeights = "stabilized"))))
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
  expect_error(computeCcwSensitivity(clones, data, fitCensoringModelsArgs, variants = list(a = list(fitCensoringModelsArgs = list(nonsense = 1)))), "no setting 'nonsense'")
  expect_error(computeCcwSensitivity(clones, data, fitCensoringModelsArgs, variants = list(a = list(other = list()))), "unknown settings")
  expect_error(suppressMessages(computeCcwSensitivity(clones, data, fitCensoringModelsArgs, variants = list(a = list(fitCensoringModelsArgs = list(postGraceModel = "x"))))), "postGraceModel")
  skip_if_not_installed("ggplot2")
  expect_s3_class(plotCcwSensitivity(sensitivity), "ggplot")
  expect_s3_class(plotCcwSensitivity(sensitivity, estimands = "hr"), "ggplot")
})

test_that("computeCcwDiagnostics compares postGraceModel and covariateTimeInteraction choices", {
  data <- simulateCloneCensorWeightData(n = 300, maxDays = 60, graceEnd = 7, seed = 8)
  clones <- designClones(data)
  on.exit({
    close(clones)
    close(data)
  })
  diagnostics <- suppressMessages(computeCcwDiagnostics(clones, data, designArgs(), numberOfGroups = 5, strataCovariateIds = 1001))
  expect_s3_class(diagnostics, "CcwDiagnostics")
  expect_equal(diagnostics$variants$postGraceModel, c("hazard", "binary", "hazard", "binary"))
  expect_equal(diagnostics$variants$covariateTimeInteraction, c("day0", "day0", "none", "none"))
  expect_equal(diagnostics$variants$specified, c(TRUE, FALSE, FALSE, FALSE))
  expect_equal(diagnostics$balance$variant, diagnostics$variants$variant)
  expect_equal(colnames(diagnostics$balance), c("variant", "covariates", "maxAbsStdDiffBefore", "maxAbsStdDiffAfter", "aboveThresholdBefore", "aboveThresholdAfter"))
  expect_equal(sort(unique(diagnostics$weightSummary$variant)), sort(diagnostics$variants$variant))
  # The hazard model with each interaction setting and the binary model once:
  overall <- diagnostics$calibration[diagnostics$calibration$grouping == "overall", ]
  expect_equal(paste(overall$strategyId, overall$model, overall$covariateTimeInteraction), c("1 binary NA", "2 hazard day0", "2 binary NA", "2 hazard none"))
  expect_true(any(diagnostics$calibration$grouping == "covariate"))
  expect_s3_class(diagnostics$gracePeriodSelection, "CcwGracePeriodSelection")
  expect_s3_class(diagnostics$positivity, "CcwPositivity")
  expect_null(diagnostics$sensitivity)
  expect_output(print(diagnostics), "Clone-censor-weight diagnostics")

  # A fit with "binary" without its binary models at G gives the weights of a fit with "hazard":
  binary <- suppressMessages(fitCensoringModels(clones, data, designArgs(postGraceModel = "binary")))
  hazard <- suppressMessages(fitCensoringModels(clones, data, designArgs()))
  fromBinary <- suppressMessages(computeCensoringWeights(clones, .withPostGraceModel(binary, "hazard")))
  direct <- suppressMessages(computeCensoringWeights(clones, hazard))
  expect_equal(collectSegments(fromBinary), collectSegments(direct))
  close(fromBinary)
  close(direct)
  # The specified variant's balance equals computeCovariateBalance() on its weights:
  weighted <- suppressMessages(computeCensoringWeights(clones, hazard))
  balance <- suppressMessages(computeCovariateBalance(weighted, data))
  close(weighted)
  expect_equal(diagnostics$balance$maxAbsStdDiffAfter[1], max(abs(balance$afterMatchingStdDiff), na.rm = TRUE))
  expect_equal(diagnostics$balance$aboveThresholdAfter[1], sum(abs(balance$afterMatchingStdDiff) > 0.1, na.rm = TRUE))

  withSensitivity <- suppressMessages(computeCcwDiagnostics(
    clones, data, designArgs(hazardModelType = "cox"),
    numberOfGroups = 5, sensitivity = TRUE,
    fitOutcomeModelArgs = createFitOutcomeModelArgs(modelType = "km", horizons = 30)
  ))
  # The Cox model has no interactions, so only postGraceModel varies:
  expect_equal(withSensitivity$variants$covariateTimeInteraction, c("none", "none"))
  expect_s3_class(withSensitivity$sensitivity, "CcwSensitivity")
  expect_output(print(withSensitivity), "Estimates under the specified settings")
})

test_that("Diagnostic plots accept minimal input", {
  skip_if_not_installed("ggplot2")
  selection <- data.frame(
    strategyId = 2, grouping = c("overall", "prediction", "prediction"), group = c(1, 1, 2),
    initiatedDays = c(400, 100, 300), notInitiatedDays = c(2000, 1500, 500),
    initiatedOutcomes = c(1, 0, 1), notInitiatedOutcomes = c(20, 12, 8), initiatedDeaths = 0, notInitiatedDeaths = 1
  )
  expect_s3_class(plotGracePeriodSelection(selection), "ggplot")
  expect_s3_class(plotGracePeriodSelection(selection[1, ], labels = c("A", "B")), "ggplot")
  expect_error(plotGracePeriodSelection(selection[, -1]), "strategyId")
})

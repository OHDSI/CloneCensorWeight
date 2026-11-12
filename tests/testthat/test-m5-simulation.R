library(testthat)

# simulateCloneCensorWeightData ------------------------------------------------------------------

test_that("simulateCloneCensorWeightData returns valid CloneCensorWeightData with settings", {
  data <- simulateCloneCensorWeightData(n = 200, seed = 1)
  on.exit(close(data))

  expect_true(isCloneCensorWeightData(data))
  expect_true(.checkCloneCensorWeightData(data))
  metaData <- attr(data, "metaData")
  expect_equal(metaData$targetId, 1)
  expect_equal(metaData$treatmentIds, c(treatment = 2))
  expect_equal(metaData$outcomeIds, 3)
  expect_equal(metaData$populationSize, 200)
  settings <- metaData$simulationSettings
  expect_true(all(.simulationSettingsFields %in% names(settings)))
  expect_equal(settings$dgp, "a")
  expect_equal(settings$initiationCoefficients, rep(0.4, 5))

  cohorts <- dplyr::collect(data$cohorts)
  expect_equal(cohorts$rowId, 1:200)
  expect_true(all(cohorts$daysToObsEnd >= 0 & cohorts$daysToObsEnd <= settings$maxDays))
  # Death ends observation:
  dead <- !is.na(cohorts$daysToDeath)
  expect_equal(cohorts$daysToDeath[dead], cohorts$daysToObsEnd[dead])

  outcomes <- dplyr::collect(data$outcomes)
  expect_equal(anyDuplicated(outcomes$rowId), 0)
  expect_true(all(outcomes$daysToEvent >= 0))
  expect_true(all(outcomes$daysToEvent <= cohorts$daysToObsEnd[outcomes$rowId]))

  treatments <- dplyr::collect(data$treatments)
  expect_equal(anyDuplicated(treatments$rowId), 0)
  expect_true(all(treatments$treatmentId == 2))
  expect_true(all(treatments$daysToStart >= 0))
  expect_equal(treatments$daysToEnd, cohorts$daysToObsEnd[treatments$rowId])
  expect_true(all(treatments$daysToStart <= treatments$daysToEnd))
})

test_that("simulateCloneCensorWeightData is reproducible with a seed", {
  tables <- function(seed, ...) {
    data <- simulateCloneCensorWeightData(n = 100, seed = seed, ...)
    on.exit(close(data))
    result <- lapply(.cloneCensorWeightDataTables, function(table) as.data.frame(dplyr::collect(data[[table]])))
    names(result) <- .cloneCensorWeightDataTables
    result$metaData <- attr(data, "metaData")
    return(result)
  }
  expect_identical(tables(123), tables(123))
  expect_identical(tables(123, dgp = "b"), tables(123, dgp = "b"))
  expect_false(identical(tables(123)$covariates, tables(124)$covariates))
  expect_false(identical(tables(123)$cohorts, tables(124)$cohorts))
})

test_that("covariates are in FeatureExtraction format", {
  data <- simulateCloneCensorWeightData(n = 200, nCovariates = 4, seed = 1)
  on.exit(close(data))

  covariateRef <- dplyr::collect(data$covariateRef)
  expect_equal(covariateRef$covariateId, c(1001, 2002, 3001, 4002))
  expect_equal(covariateRef$analysisId, c(1, 2, 1, 2))
  expect_equal(covariateRef$conceptId, 1:4)
  analysisRef <- dplyr::collect(data$analysisRef)
  expect_equal(analysisRef$analysisId, c(1, 2))
  expect_equal(analysisRef$isBinary, c("Y", "N"))
  expect_true(FeatureExtraction::isCovariateData(data))

  covariates <- dplyr::collect(data$covariates)
  binary <- covariates[covariates$covariateId %in% c(1001, 3001), ]
  expect_true(all(binary$covariateValue == 1))
  normal <- covariates[covariates$covariateId %in% c(2002, 4002), ]
  expect_equal(nrow(normal), 400)

  aggregated <- FeatureExtraction::aggregateCovariates(data)
  on.exit(close(aggregated), add = TRUE)
  expect_equal(sort(dplyr::collect(aggregated$covariates)$covariateId), c(1001, 3001))
  expect_equal(sort(dplyr::collect(aggregated$covariatesContinuous)$covariateId), c(2002, 4002))

  # No covariates:
  empty <- simulateCloneCensorWeightData(n = 20, nCovariates = 0, seed = 1)
  on.exit(close(empty), add = TRUE)
  expect_equal(nrow(dplyr::collect(empty$covariates)), 0)
  expect_equal(nrow(dplyr::collect(empty$covariateRef)), 0)
})

test_that("CohortMethod::createStudyPopulation runs on simulated data", {
  data <- simulateCloneCensorWeightData(n = 200, seed = 1)
  on.exit(close(data))
  cohorts <- dplyr::collect(data$cohorts)
  outcomes <- dplyr::collect(data$outcomes)
  for (riskStart in 0:1) {
    population <- CohortMethod::createStudyPopulation(
      cohortMethodData = data,
      outcomeId = 3,
      createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
        riskWindowStart = riskStart,
        riskWindowEnd = 365,
        endAnchor = "cohort start"
      )
    )
    expect_true(all(c("daysToDeath", "riskStart", "riskEnd", "daysToEvent", "outcomeCount") %in% colnames(population)))
    # Removed: outcome before riskStart, or no day at risk.
    priorOutcome <- outcomes$rowId[outcomes$daysToEvent < riskStart]
    expected <- setdiff(cohorts$rowId[cohorts$daysToObsEnd >= riskStart], priorOutcome)
    expect_equal(sort(population$rowId), expected)
    expect_equal(population$riskEnd, pmin(population$daysToObsEnd, 365))
    withOutcome <- population$rowId %in% outcomes$rowId
    expect_equal(population$outcomeCount, as.numeric(withOutcome))
    expect_equal(population$daysToEvent[withOutcome], outcomes$daysToEvent[match(population$rowId[withOutcome], outcomes$rowId)])
  }
})

test_that("event rates match the DGP hazards", {
  n <- 20000
  # No death and no loss to follow-up, so every latent day up to maxDays is observed:
  data <- simulateCloneCensorWeightData(
    n = n, deathHazard = 0, observationEndHazard = 0, treatmentEffect = 0, maxDays = 60, seed = 1
  )
  on.exit(close(data))
  settings <- attr(data, "metaData")$simulationSettings
  x <- m5CovariateMatrix(data)
  hazard <- function(baseline, coefficients) stats::plogis(stats::qlogis(baseline) + as.vector(x %*% coefficients))
  expect_equal(dplyr::collect(data$cohorts)$daysToObsEnd, rep(60L, n))

  treatments <- dplyr::collect(data$treatments)
  outcomes <- dplyr::collect(data$outcomes)
  hA <- hazard(settings$initiationHazard, settings$initiationCoefficients)
  hY <- hazard(settings$outcomeHazard, settings$outcomeCoefficients)
  for (day in c(0, 10, 30, 60)) {
    expected <- mean(1 - (1 - hA)^(day + 1))
    observed <- sum(treatments$daysToStart <= day) / n
    expect_lt(abs(observed - expected), 4 * m5MonteCarloSe(expected, n))
    expected <- mean(1 - (1 - hY)^(day + 1))
    observed <- sum(outcomes$daysToEvent <= day) / n
    expect_lt(abs(observed - expected), 4 * m5MonteCarloSe(expected, n))
  }

  # Death and observation end, with no outcome:
  data2 <- simulateCloneCensorWeightData(
    n = n, outcomeHazard = 0, deathHazard = 0.002, observationEndHazard = 0.003, maxDays = 60, seed = 2
  )
  on.exit(close(data2), add = TRUE)
  x <- m5CovariateMatrix(data2)
  hD <- hazard(0.002, settings$deathCoefficients)
  cohorts <- dplyr::collect(data2$cohorts)
  expect_equal(nrow(dplyr::collect(data2$outcomes)), 0)
  for (day in c(0, 30, 59)) {
    # Death is observed when it comes on or before the observation end:
    expected <- mean(rowSums(sapply(0:day, function(d) (1 - hD)^d * hD * (1 - 0.003)^d)))
    observed <- mean(!is.na(cohorts$daysToDeath) & cohorts$daysToDeath <= day)
    expect_lt(abs(observed - expected), 4 * m5MonteCarloSe(expected, n))
    # Observation ends (by death or loss to follow-up):
    expected <- 1 - mean(((1 - hD) * (1 - 0.003))^(day + 1))
    observed <- mean(cohorts$daysToObsEnd <= day)
    expect_lt(abs(observed - expected), 4 * m5MonteCarloSe(expected, n))
  }
})

test_that("initiation and the outcome are confounded in the specified direction", {
  confounding <- function(sign) {
    data <- simulateCloneCensorWeightData(
      n = 5000,
      initiationCoefficients = sign * rep(0.5, 5),
      outcomeCoefficients = rep(0.5, 5),
      seed = 1
    )
    on.exit(close(data))
    x <- m5CovariateMatrix(data)
    treatments <- dplyr::collect(data$treatments)
    outcomes <- dplyr::collect(data$outcomes)
    startedBy30 <- seq_len(nrow(x)) %in% treatments$rowId[treatments$daysToStart <= 30]
    withOutcome <- seq_len(nrow(x)) %in% outcomes$rowId
    return(list(
      initiation = colMeans(x[startedBy30, ]) - colMeans(x[!startedBy30, ]),
      outcome = colMeans(x[withOutcome, ]) - colMeans(x[!withOutcome, ])
    ))
  }
  positive <- confounding(1)
  expect_true(all(positive$initiation > 0))
  expect_true(all(positive$outcome > 0))
  negative <- confounding(-1)
  expect_true(all(negative$initiation < 0))
  expect_true(all(negative$outcome > 0))
})

test_that("the treatment changes the outcome hazard only when the DGP says so", {
  # A near-complete protective effect leaves (almost) no outcome where the treatment acts.
  for (dgp in c("a", "b")) {
    data <- simulateCloneCensorWeightData(
      n = 5000, dgp = dgp, graceEnd = 10, outcomeHazard = 0.01, treatmentEffect = -30, effectRampDays = 1, seed = 1
    )
    cohorts <- dplyr::collect(data$cohorts)
    outcomes <- dplyr::collect(data$outcomes)
    treatments <- dplyr::collect(data$treatments)
    close(data)
    startDay <- treatments$daysToStart[match(outcomes$rowId, treatments$rowId)]
    if (dgp == "a") {
      # Only on days after G, and only for persons who started by G:
      acted <- !is.na(startDay) & startDay <= 10 & outcomes$daysToEvent > 10
      expect_true(any(!is.na(startDay) & startDay <= 10 & outcomes$daysToEvent > startDay))
      expect_true(any(!is.na(startDay) & startDay > 10 & outcomes$daysToEvent > startDay))
    } else {
      # From the day after initiation:
      acted <- !is.na(startDay) & outcomes$daysToEvent > startDay
      expect_true(any(!is.na(startDay) & outcomes$daysToEvent == startDay))
    }
    expect_equal(sum(acted), 0)
  }
})

test_that("simulateCloneCensorWeightData checks its arguments", {
  expect_error(simulateCloneCensorWeightData(n = 0), "n")
  expect_error(simulateCloneCensorWeightData(n = 10, dgp = "c"), "dgp")
  expect_error(simulateCloneCensorWeightData(n = 10, initiationHazard = 0), "initiationHazard")
  expect_error(simulateCloneCensorWeightData(n = 10, deathHazard = 1), "deathHazard")
  expect_error(simulateCloneCensorWeightData(n = 10, outcomeCoefficients = c(1, 2)), "outcomeCoefficients")
  expect_error(simulateCloneCensorWeightData(n = 10, effectRampDays = 0), "effectRampDays")
})

# computeTrueRisks ------------------------------------------------------------------------------

test_that("computeTrueRisks equals the closed form for neverInitiate vs initiateWithin(0)", {
  nTruth <- 2e5
  strategies <- list(
    createTreatmentStrategy(1, "Initiate on day 0", "initiateWithin", graceEnd = 0),
    createTreatmentStrategy(2, "Never initiate", "neverInitiate")
  )
  for (dgp in c("a", "b")) {
    settings <- m5TinySettings(dgp, graceEnd = 0)
    for (graceRegime in c("binaryAtG", "naturalGrace")) {
      truth <- computeTrueRisks(settings, strategies, horizons = 0:8, graceRegime = graceRegime, nTruth = nTruth, seed = 1)
      expect_equal(colnames(truth), c("strategyId", "day", "risk"))
      expect_type(truth$strategyId, "integer")
      expect_type(truth$day, "integer")
      expect_equal(truth$strategyId, rep(1:2, each = 9))
      expect_equal(truth$day, rep(0:8, 2))
      for (i in 1:2) {
        exact <- m5ExactRisks(settings, strategies[[i]], graceRegime)
        estimate <- truth$risk[truth$strategyId == i]
        expect_true(all(abs(estimate - exact) < 4 * m5MonteCarloSe(exact, nTruth)), label = sprintf("dgp %s, %s, strategy %d", dgp, graceRegime, i))
      }
    }
  }
  # With everybody starting on day 0, a protective effect lowers the risk after day 0:
  truth <- computeTrueRisks(m5TinySettings("b", 0), strategies, horizons = c(0, 8), nTruth = nTruth, seed = 1)
  expect_equal(truth$risk[1], truth$risk[3])
  expect_lt(truth$risk[2], truth$risk[4])
})

test_that("computeTrueRisks equals a brute-force calculation for grace regimes with G = 2", {
  nTruth <- 2e5
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within 2 days", "initiateWithin", graceEnd = 2),
    createTreatmentStrategy(2, "Do not initiate within 2 days", "notInitiateWithin", graceEnd = 2),
    createTreatmentStrategy(3, "Never initiate", "neverInitiate"),
    createTreatmentStrategy(4, "Initiate within 4 days", "initiateWithin", graceEnd = 4),
    createTreatmentStrategy(5, "Initiate within 1 day", "initiateWithin", graceEnd = 1)
  )
  cases <- list(
    list(dgp = "a", treatmentEffect = log(0.3)),
    list(dgp = "b", treatmentEffect = log(0.3)),
    list(dgp = "b", treatmentEffect = log(5))
  )
  for (case in cases) {
    settings <- m5TinySettings(case$dgp, graceEnd = 2, treatmentEffect = case$treatmentEffect)
    for (graceRegime in c("binaryAtG", "naturalGrace")) {
      truth <- computeTrueRisks(settings, strategies, horizons = 0:8, graceRegime = graceRegime, nTruth = nTruth, seed = 2)
      for (i in seq_along(strategies)) {
        exact <- m5ExactRisks(settings, strategies[[i]], graceRegime)
        estimate <- truth$risk[truth$strategyId == i]
        label <- sprintf("dgp %s, effect %.1f, %s, strategy %d", case$dgp, exp(case$treatmentEffect), graceRegime, i)
        expect_true(all(abs(estimate - exact) < 4 * m5MonteCarloSe(exact, nTruth)), label = label)
      }
    }
  }
  # The check can tell the two regimes apart:
  binaryAtG <- m5ExactRisks(settings, strategies[[1]], "binaryAtG")
  naturalGrace <- m5ExactRisks(settings, strategies[[1]], "naturalGrace")
  expect_gt(max(abs(binaryAtG - naturalGrace) / m5MonteCarloSe(binaryAtG, nTruth)), 12)
})

test_that("under DGP (a) the two grace regimes give the same truth", {
  settings <- m5Settings(dgp = "a", graceEnd = 30)
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within 30 days", "initiateWithin", graceEnd = 30),
    createTreatmentStrategy(2, "Do not initiate within 30 days", "notInitiateWithin", graceEnd = 30),
    createTreatmentStrategy(3, "Never initiate", "neverInitiate")
  )
  horizons <- c(0, 15, 30, 31, 60, 180, 365)
  binaryAtG <- computeTrueRisks(settings, strategies, horizons, graceRegime = "binaryAtG", nTruth = 1e5, seed = 1)
  naturalGrace <- computeTrueRisks(settings, strategies, horizons, graceRegime = "naturalGrace", nTruth = 1e5, seed = 1)
  # The same random numbers, and the same outcome hazards after G, so the truths are identical:
  expect_equal(binaryAtG, naturalGrace, tolerance = 0)
  # The treatment acts only after G:
  expect_equal(binaryAtG$risk[binaryAtG$strategyId == 1 & binaryAtG$day <= 30], binaryAtG$risk[binaryAtG$strategyId == 2 & binaryAtG$day <= 30])
  # Starting after G has no effect under DGP (a):
  expect_equal(binaryAtG$risk[binaryAtG$strategyId == 2], binaryAtG$risk[binaryAtG$strategyId == 3])
  # The protective effect:
  expect_lt(binaryAtG$risk[binaryAtG$strategyId == 1 & binaryAtG$day == 365], binaryAtG$risk[binaryAtG$strategyId == 2 & binaryAtG$day == 365])
})

test_that("under DGP (b) the grace regimes differ in the expected direction", {
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within 30 days", "initiateWithin", graceEnd = 30),
    createTreatmentStrategy(2, "Do not initiate within 30 days", "notInitiateWithin", graceEnd = 30)
  )
  horizons <- c(0, 15, 30, 31, 45, 60, 180, 365)
  for (treatmentEffect in c(log(0.5), log(2))) {
    settings <- m5Settings(dgp = "b", treatmentEffect = treatmentEffect, effectRampDays = 30)
    binaryAtG <- computeTrueRisks(settings, strategies, horizons, graceRegime = "binaryAtG", nTruth = 1e5, seed = 1)
    naturalGrace <- computeTrueRisks(settings, strategies, horizons, graceRegime = "naturalGrace", nTruth = 1e5, seed = 1)
    # Both follow the natural course through G, and other strategy types ignore the regime:
    expect_equal(binaryAtG$risk[binaryAtG$day <= 30], naturalGrace$risk[naturalGrace$day <= 30], tolerance = 0)
    expect_equal(binaryAtG[binaryAtG$strategyId == 2, ], naturalGrace[naturalGrace$strategyId == 2, ], tolerance = 0)
    # Binary-at-G moves the start of late starters earlier, where the effect has built up longer:
    later <- binaryAtG$strategyId == 1 & binaryAtG$day > 30
    difference <- binaryAtG$risk[later] - naturalGrace$risk[later]
    if (treatmentEffect < 0) {
      expect_true(all(difference < 0))
    } else {
      expect_true(all(difference > 0))
    }
  }
})

test_that("computeTrueRisks is reproducible and within Monte Carlo error across seeds", {
  settings <- m5Settings(dgp = "b")
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within 30 days", "initiateWithin", graceEnd = 30),
    createTreatmentStrategy(2, "Never initiate", "neverInitiate")
  )
  nTruth <- 5e4
  run <- function(seed) computeTrueRisks(settings, strategies, c(30, 90, 365), nTruth = nTruth, seed = seed)
  expect_identical(run(1), run(1))
  first <- run(1)
  second <- run(2)
  expect_false(identical(first$risk, second$risk))
  expect_true(all(abs(first$risk - second$risk) < 4 * sqrt(2) * m5MonteCarloSe(first$risk, nTruth)))
})

test_that("computeTrueRisks runs with nTruth = 1e6 in well under a minute", {
  skip_on_cran()
  settings <- m5Settings(dgp = "b")
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within 30 days", "initiateWithin", graceEnd = 30),
    createTreatmentStrategy(2, "Do not initiate within 30 days", "notInitiateWithin", graceEnd = 30)
  )
  timing <- system.time(truth <- computeTrueRisks(settings, strategies, horizons = c(90, 365), nTruth = 1e6, seed = 1))
  expect_lt(timing[["elapsed"]], 30)
  expect_equal(nrow(truth), 4)
})

test_that("computeTrueRisks checks its arguments", {
  settings <- m5Settings()
  initiate <- createTreatmentStrategy(1, "Initiate within 30 days", "initiateWithin", graceEnd = 30)
  expect_error(computeTrueRisks(settings, list(initiate), horizons = 366, nTruth = 10), "horizons")
  expect_error(computeTrueRisks(settings, list(initiate), horizons = 10, graceRegime = "other", nTruth = 10), "graceRegime")
  expect_error(computeTrueRisks(settings[-1], list(initiate), horizons = 10, nTruth = 10), "settings")
  expect_error(computeTrueRisks(settings, list(initiate$toList()), horizons = 10, nTruth = 10), "strategies\\[\\[1\\]\\]")
  lateStart <- createTreatmentStrategy(1, "Initiate on days 5 to 30", "initiateWithin", graceStart = 5, graceEnd = 30)
  expect_error(computeTrueRisks(settings, list(lateStart), horizons = 10, nTruth = 10), "Not yet supported in this version: graceStart")
  prohibited <- createTreatmentStrategy(1, "Initiate within 30 days", "initiateWithin", graceEnd = 30, prohibitedTreatments = "comparatorTreatment")
  expect_error(computeTrueRisks(settings, list(prohibited), horizons = 10, nTruth = 10), "Not yet supported in this version: prohibitedTreatments")
  comparator <- createTreatmentStrategy(1, "Initiate within 30 days", "initiateWithin", treatment = "comparatorTreatment", graceEnd = 30)
  expect_error(computeTrueRisks(settings, list(comparator), horizons = 10, nTruth = 10), "comparatorTreatment")
})

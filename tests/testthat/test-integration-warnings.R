library(testthat)

# Tests of the warnings of Phase1Contracts section 11.6: each fires when it should and is silent
# otherwise.
# - fitOutcomeModel() with modelType = "logistic": a grid with no interval ending at the end of the
#   grace period, and a smooth time basis with a grace period;
# - fitCensoringModels(): censoring models of initiation on day 0 (every strategy type) with covariate
#   windows that end on day 0 or later;
# - bootstrapCcwAnalysis(): replicates that sit predominantly on one side of the full-data estimate.

gridWarning <- "No interval of the pooled logistic outcome model's grid"
smoothWarning <- "smooth time basis"
day0Warning <- "models? initiation on day 0"
shiftWarning <- "sit predominantly on one side"

# Pooled logistic outcome model --------------------------------------------------------------------

# Simulated data, clones (without competing events) and censoring models (without covariates) of
# initiateWithin(7) vs notInitiateWithin(7). The caller closes `data` and `clones`.
pooledAnalysis <- function(riskStart = 0, sameDayOrdering = "eventsFirst") {
  data <- simulateCloneCensorWeightData(n = 300, maxDays = 30, graceEnd = 7, seed = 1)
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowStart = riskStart, riskWindowEnd = 30, endAnchor = "cohort start")
  ))
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
    createTreatmentStrategy(2, "Do not initiate within 7 days", "notInitiateWithin", graceEnd = 7)
  )
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = strategies, sameDayOrdering = sameDayOrdering, competingEvents = "none")))
  models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(useCovariates = FALSE)))
  return(list(data = data, clones = clones, models = models, riskStart = riskStart))
}

closeAnalysis <- function(analysis) {
  close(analysis$clones)
  close(analysis$data)
}

# The warnings of the pooled logistic outcome model with `timeBasis`, at a horizon that ends an interval
# (so that the model does not warn about the horizon).
pooledWarnings <- function(analysis, timeBasis) {
  weighted <- suppressMessages(computeCensoringWeights(analysis$clones, analysis$models, createComputeWeightsArgs(), outcomeTimeBasis = timeBasis))
  on.exit(close(weighted))
  intervalLength <- timeBasis$intervalLength
  horizon <- analysis$riskStart - 1 + intervalLength * (20 %/% intervalLength)
  args <- createFitOutcomeModelArgs(modelType = "logistic", horizons = horizon, timeBasis = timeBasis, armByTimeInteraction = TRUE)
  return(capture_warnings(suppressMessages(fitOutcomeModel(weighted, args))))
}

test_that("The pooled logistic outcome model warns when no interval ends at the end of the grace period", {
  analysis <- pooledAnalysis()
  on.exit(closeAnalysis(analysis))
  # The grace period of both strategies ends at the end of day 7, t = 8:
  warnings <- pooledWarnings(analysis, createTimeBasis("indicator", intervalLength = 3))
  expect_length(warnings, 1)
  expect_match(warnings, gridWarning)
  expect_match(warnings, "intervals of 3 days counted from riskStart = 0", fixed = TRUE)
  expect_match(warnings, "strategy 1 (Initiate within 7 days) at the end of day 7 (t = 8); strategy 2 (Do not initiate within 7 days) at the end of day 7 (t = 8)", fixed = TRUE)
  expect_match(warnings, "indicator time basis on such a grid, with armByTimeInteraction = TRUE, reproduces the weighted Kaplan-Meier risks", fixed = TRUE)
  # (Grids of 2 and 8 days in the default tier; 1, 2, 4 and 8 in the full tier.)
  for (intervalLength in testConfigurations(c(1, 2, 4, 8), c(2, 8))) {
    expect_length(pooledWarnings(analysis, createTimeBasis("indicator", intervalLength = intervalLength)), 0)
  }

  # The end of day 7 is at t = 6 with riskStart = 2:
  later <- pooledAnalysis(riskStart = 2)
  on.exit(closeAnalysis(later), add = TRUE)
  expect_length(pooledWarnings(later, createTimeBasis("indicator", intervalLength = 3)), 0)
  warnings <- pooledWarnings(later, createTimeBasis("indicator", intervalLength = 4))
  expect_length(warnings, 1)
  expect_match(warnings, "at the end of day 7 (t = 6)", fixed = TRUE)

  # With deviationFirst, the clones that deviate on day 7 end on day 6, t = 7:
  deviationFirst <- pooledAnalysis(sameDayOrdering = "deviationFirst")
  on.exit(closeAnalysis(deviationFirst), add = TRUE)
  expect_length(pooledWarnings(deviationFirst, createTimeBasis("indicator", intervalLength = 7)), 0)
  warnings <- pooledWarnings(deviationFirst, createTimeBasis("indicator", intervalLength = 8))
  expect_length(warnings, 1)
  expect_match(warnings, "at the end of day 6 (t = 7)", fixed = TRUE)

  # A grace period that ends before riskStart does not censor clones during follow-up:
  afterGrace <- pooledAnalysis(riskStart = 9)
  on.exit(closeAnalysis(afterGrace), add = TRUE)
  expect_length(pooledWarnings(afterGrace, createTimeBasis("indicator", intervalLength = 3)), 0)
  # (The odds ratio of this fit is not identified, no target outcome falling in its reference interval, and
  # its variance can come out negative, with an "NaNs produced" warning, depending on the weights.)
  warnings <- pooledWarnings(afterGrace, createTimeBasis("naturalSpline", intervalLength = 1, df = 3))
  expect_false(any(grepl(gridWarning, warnings) | grepl(smoothWarning, warnings)))
})

test_that("The pooled logistic outcome model warns about a smooth time basis with a grace period", {
  analysis <- pooledAnalysis()
  on.exit(closeAnalysis(analysis))
  smoothBases <- list(
    naturalSpline = createTimeBasis("naturalSpline", intervalLength = 1, df = 4),
    polynomial = createTimeBasis("polynomial", intervalLength = 1, degree = 2),
    linear = createTimeBasis("linear", intervalLength = 1)
  )
  for (type in names(smoothBases)) {
    warnings <- pooledWarnings(analysis, smoothBases[[type]])
    expect_length(warnings, 1)
    expect_match(warnings, sprintf("smooth time basis (%s) cannot follow the change in the hazards at the end of the grace period: strategy 1 (Initiate within 7 days) at the end of day 7 (t = 8)", type), fixed = TRUE, info = type)
    expect_match(warnings, "Use an indicator time basis (createTimeBasis(type = \"indicator\")) on a grid with an interval ending there, with armByTimeInteraction = TRUE", fixed = TRUE, info = type)
  }
  # A smooth basis on a grid with no interval ending at the end of the grace period gets both warnings:
  warnings <- pooledWarnings(analysis, createTimeBasis("naturalSpline", intervalLength = 3, df = 3))
  expect_length(warnings, 2)
  expect_match(warnings[1], gridWarning)
  expect_match(warnings[2], smoothWarning)
  # The indicator basis on a daily grid gets neither:
  expect_length(pooledWarnings(analysis, createTimeBasis("indicator", intervalLength = 1)), 0)
})

# Covariate windows that end on day 0 ----------------------------------------------------------------

# Simulated data (whose covariate windows end on day -1) with the covariate windows of analysis 1
# (binary covariates) and 2 (normal covariates) ending on `endDays`. The caller closes it.
windowData <- function(endDays = c(-1, -1)) {
  data <- simulateCloneCensorWeightData(n = 200, maxDays = 20, graceEnd = 5, seed = 1)
  analysisRef <- dplyr::collect(data$analysisRef)
  analysisRef$endDay <- as.numeric(endDays[match(analysisRef$analysisId, 1:2)])
  data$analysisRef <- analysisRef
  return(data)
}

windowClones <- function(data, comparatorType = "notInitiateWithin", firstType = "initiateWithin") {
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 20, endAnchor = "cohort start")
  ))
  graceEnd <- function(type) if (type == "neverInitiate") NULL else 5
  strategies <- list(
    createTreatmentStrategy(1, "First", firstType, graceEnd = graceEnd(firstType)),
    createTreatmentStrategy(2, "Second", comparatorType, graceEnd = graceEnd(comparatorType))
  )
  return(suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = strategies))))
}

windowWarnings <- function(clones, data, ...) {
  args <- createFitCensoringModelsArgs(prior = Cyclops::createPrior("none"), errorOnHighCorrelation = FALSE, ...)
  return(capture_warnings(suppressMessages(fitCensoringModels(clones, data, args))))
}

test_that("fitCensoringModels warns once when initiation on day 0 meets covariate windows that end on day 0", {
  # The simulated covariates are baseline covariates, with windows that end on day -1:
  data <- simulateCloneCensorWeightData(n = 200, maxDays = 20, graceEnd = 5, seed = 1)
  on.exit(close(data))
  expect_equal(dplyr::pull(data$analysisRef, "endDay"), c(-1, -1))
  clones <- windowClones(data)
  on.exit(close(clones), add = TRUE)
  expect_length(windowWarnings(clones, data), 0)

  # Windows that end on day 0, as with FeatureExtraction's default endDays = 0. The models of both
  # strategies predict initiation on day 0; one warning, although three models are fitted:
  day0 <- windowData(c(0, 0))
  on.exit(close(day0), add = TRUE)
  warnings <- windowWarnings(clones, day0)
  expect_length(warnings, 1)
  expect_match(warnings, "The censoring models of strategies 1 (First) and 2 (Second) model initiation on day 0, but the covariate windows of analyses SimulatedBinary (endDay = 0), SimulatedNormal (endDay = 0) end on day 0 or later", fixed = TRUE)
  expect_match(warnings, "endDays = -1", fixed = TRUE)
  # Silent without covariates:
  expect_length(windowWarnings(clones, day0, useCovariates = FALSE), 0)

  # Only the covariates the models use count: analysis 1 (covariates 1001, 3001, 5001) ends on day 0,
  # analysis 2 on day -1:
  mixed <- windowData(c(0, -1))
  on.exit(close(mixed), add = TRUE)
  warnings <- windowWarnings(clones, mixed)
  expect_length(warnings, 1)
  expect_match(warnings, "analyses SimulatedBinary (endDay = 0) end", fixed = TRUE)
  expect_length(windowWarnings(clones, mixed, excludeCovariateIds = c(1001, 3001, 5001)), 0)
  expect_length(windowWarnings(clones, mixed, includeCovariateIds = c(2002, 4002)), 0)
  expect_length(windowWarnings(clones, mixed, includeCovariateIds = c(1001, 2002)), 1)
  # Analyses without a window (such as demographics) are not a concern:
  noWindow <- windowData(c(NA, NA))
  on.exit(close(noWindow), add = TRUE)
  expect_length(windowWarnings(clones, noWindow), 0)

  # The initiation hazard models of notInitiateWithin and neverInitiate strategies also start on day 0,
  # so a design without an initiateWithin strategy warns too (review of contract 11.6):
  noInitiation <- windowClones(day0, comparatorType = "neverInitiate", firstType = "notInitiateWithin")
  on.exit(close(noInitiation), add = TRUE)
  warnings <- windowWarnings(noInitiation, day0)
  expect_length(warnings, 1)
  expect_match(warnings, "The censoring models of strategies 1 (First) and 2 (Second) model initiation on day 0", fixed = TRUE)
  expect_length(windowWarnings(noInitiation, windowData()), 0)

  # The hazard models of the days before graceStart > 0 of initiateWithin strategies start on day 0
  # too, also with initiateWithinModel = "binary":
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = day0,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 20, endAnchor = "cohort start")
  ))
  lateStarts <- suppressMessages(createClones(population, day0, createCreateClonesArgs(strategies = list(
    createTreatmentStrategy(1, "First", "initiateWithin", graceStart = 2, graceEnd = 5),
    createTreatmentStrategy(2, "Second", "initiateWithin", graceStart = 3, graceEnd = 5)
  ))))
  on.exit(close(lateStarts), add = TRUE)
  warnings <- windowWarnings(lateStarts, day0)
  expect_length(warnings, 1)
  expect_match(warnings, "The censoring models of strategies 1 (First) and 2 (Second) model initiation on day 0", fixed = TRUE)
})

test_that("bootstrapCcwAnalysis warns about day-0 covariate windows once, for the full data", {
  data <- windowData(c(0, 0))
  on.exit(close(data))
  warnings <- capture_warnings(suppressMessages(bootstrapCcwAnalysis(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 20, endAnchor = "cohort start"),
    createClonesArgs = createCreateClonesArgs(strategies = list(
      createTreatmentStrategy(1, "First", "initiateWithin", graceEnd = 5),
      createTreatmentStrategy(2, "Second", "notInitiateWithin", graceEnd = 5)
    )),
    fitCensoringModelsArgs = createFitCensoringModelsArgs(prior = Cyclops::createPrior("none")),
    computeWeightsArgs = createComputeWeightsArgs(),
    # (One replicate in the default tier, which also shows that the replicate does not warn again.)
    fitOutcomeModelArgs = createFitOutcomeModelArgs(modelType = "km", horizons = 20, bootstrapArgs = createBootstrapArgs(replicates = testConfigurations(2, 1), seed = 1))
  )))
  expect_length(grep(day0Warning, warnings), 1)
})

# Bootstrap replicates on one side of the estimate -----------------------------------------------------

test_that("The bootstrap shift criterion: the median replicate, with a sign test, with at least 20 replicates", {
  t0 <- dplyr::tibble(
    estimand = c("hr", "risk_target_30", "rd_30", "rr_30", "rd_60"),
    estimate = c(0.8, 0.1, 0.02, 1.25, 0.03),
    logEstimate = c(log(0.8), NA, NA, log(1.25), NA)
  )
  medianOnly <- c(rep(-3, 7), rep(1, 13))
  replicates <- cbind(
    # All 20 below the estimate (on the log scale):
    hr = exp(log(0.8) - seq(0.01, 0.2, length.out = 20)),
    risk_target_30 = rep(0.5, 20),
    # Symmetric around the estimate:
    rd_30 = 0.02 + c(-10:-1, 1:10) / 1000,
    rr_30 = rep(5, 20),
    # 13 above and 7 below (sign test p = 0.26), but the median 1 / sd = 0.51 SD above:
    rd_60 = 0.03 + medianOnly / 100
  )
  shift <- .bootstrapReplicateShift(t0, replicates)
  # Only the hazard ratio and the risk differences:
  expect_equal(shift$estimand, c("hr", "rd_30", "rd_60"))
  expect_equal(shift$replicates, c(20, 20, 20))
  expect_equal(shift$above, c(0, 10, 13))
  expect_equal(shift$below, c(20, 10, 7))
  expect_equal(shift$signTestP, c(binom.test(0, 20)$p.value, 1, binom.test(13, 20)$p.value))
  logHr <- log(replicates[, "hr"])
  expect_equal(shift$medianShift, c((median(logHr) - log(0.8)) / sd(logHr), 0, 1 / sd(medianOnly)))
  expect_gt(shift$signTestP[3], 0.01)
  expect_gt(shift$medianShift[3], 0.5)
  expect_equal(shift$shifted, c(TRUE, FALSE, TRUE))

  # Non-finite replicates are left out, and with fewer than 20 finite replicates nothing is assessed:
  replicates[1, "hr"] <- NA
  replicates[2, "rd_30"] <- Inf
  shift <- .bootstrapReplicateShift(t0, replicates)
  expect_equal(shift$replicates, c(19, 19, 20))
  expect_equal(shift$below[1], 19)
  expect_equal(shift$shifted, c(NA, NA, TRUE))

  # With many replicates, the sign test is significant for a negligible shift, which is not flagged
  # without a median shift of more than 0.25 SD (review of contract 11.6):
  small <- .bootstrapReplicateShift(t0[3, ], cbind(rd_30 = 0.02 + (stats::qnorm(stats::ppoints(1000)) + 0.15) / 100))
  expect_lt(small$signTestP, 0.001)
  expect_equal(small$medianShift, 0.15, tolerance = 0.01)
  expect_false(small$shifted)
  large <- .bootstrapReplicateShift(t0[3, ], cbind(rd_30 = 0.02 + (stats::qnorm(stats::ppoints(1000)) + 0.3) / 100))
  expect_lt(large$signTestP, 0.001)
  expect_true(large$shifted)

  # Replicates equal to the estimate are not shifted:
  same <- .bootstrapReplicateShift(t0[3, ], cbind(rd_30 = rep(0.02, 25)))
  expect_equal(same$above + same$below, 0)
  expect_true(is.na(same$signTestP))
  expect_true(is.na(same$medianShift))
  expect_false(same$shifted)
})

test_that("bootstrapCcwAnalysis warns when the replicates sit on one side of the estimate, and records the statistics", {
  # 20 replicates take about a minute. With fewer than 20 replicates nothing is assessed, so the other
  # bootstrap tests (test-integration-bootstrap.R) do not warn.
  skip_on_cran()
  # Full tier: the criterion itself (.bootstrapReplicateShift()) is tested by hand in the default tier
  # (the test above).
  skipUnlessFullTests()
  data <- simulateCloneCensorWeightData(
    n = 120, maxDays = 15, graceEnd = 3, initiationHazard = 0.1, outcomeHazard = 0.02, deathHazard = 0.005,
    initiationCoefficients = rep(0.4, 5), seed = 1
  )
  on.exit(close(data))
  cohorts <- dplyr::collect(data$cohorts)
  treatments <- dplyr::collect(data$treatments)
  outcomes <- dplyr::collect(data$outcomes)
  personSeqIds <- sort(unique(cohorts$personSeqId))
  # The persons who start by day 3 and have the outcome after day 3 count only in the initiateWithin(3)
  # strategy after day 3: every replicate draws them twice (and two other persons not at all and two
  # others twice), so every replicate's target risk at day 15, and so its risk difference, is above the
  # full-data estimate.
  starters <- treatments$rowId[treatments$daysToStart <= 3]
  lateOutcomes <- outcomes$rowId[outcomes$daysToEvent > 3]
  shifters <- cohorts$personSeqId[cohorts$rowId %in% intersect(starters, lateOutcomes)]
  expect_gt(length(shifters), 2)
  others <- setdiff(personSeqIds, shifters)
  set.seed(1)
  counts <- matrix(1L, length(personSeqIds), 20)
  for (replicate in 1:20) {
    counts[personSeqIds %in% shifters, replicate] <- 2L
    changed <- sample(others, 4)
    counts[match(changed, personSeqIds), replicate] <- c(0L, 0L, 2L, 2L)
  }
  settings <- list(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 15, endAnchor = "cohort start"),
    createClonesArgs = createCreateClonesArgs(
      strategies = list(
        createTreatmentStrategy(1, "Initiate within 3 days", "initiateWithin", graceEnd = 3),
        createTreatmentStrategy(2, "Do not initiate within 3 days", "notInitiateWithin", graceEnd = 3)
      ),
      competingEvents = "none"
    ),
    fitCensoringModelsArgs = createFitCensoringModelsArgs(useCovariates = FALSE),
    computeWeightsArgs = createComputeWeightsArgs(),
    fitOutcomeModelArgs = createFitOutcomeModelArgs(modelType = "cox", horizons = 15)
  )
  warnings <- capture_warnings(result <- suppressMessages(do.call(bootstrapCcwAnalysis, c(settings, list(counts = counts)))))
  shift <- result$metaData$bootstrap$replicateShift
  expect_equal(shift$estimand, c("hr", "rd_15"))
  expect_equal(shift$above[2], 20)
  expect_true(shift$shifted[2])
  # The statistics, by hand from the replicates:
  for (i in 1:2) {
    t <- result$replicates[, shift$estimand[i]]
    estimate <- result$estimates$estimate[result$estimates$estimand == shift$estimand[i]]
    if (i == 1) {
      t <- log(t)
      estimate <- log(estimate)
    }
    expect_equal(shift$replicates[i], sum(is.finite(t)))
    expect_equal(shift$above[i], sum(t > estimate))
    expect_equal(shift$below[i], sum(t < estimate))
    expect_equal(shift$signTestP[i], binom.test(sum(t > estimate), sum(t != estimate))$p.value)
    expect_equal(shift$medianShift[i], (median(t) - estimate) / sd(t))
  }
  shiftWarnings <- grep(shiftWarning, warnings, value = TRUE)
  expect_length(shiftWarnings, 1)
  expect_match(shiftWarnings, "rd_15 (20 replicates above it and 0 below, sign test p = 1.9e-06; median replicate", fixed = TRUE)
  expect_match(shiftWarnings, "metaData$bootstrap$replicateShift", fixed = TRUE)
  expect_length(setdiff(warnings, shiftWarnings), 0)
})

test_that("bootstrapCcwAnalysis words the shift warning from the statistics it records", {
  # The light form of the test above: 2 replicates, with the criterion (tested by hand above) mocked to
  # flag the risk difference.
  data <- simulateCloneCensorWeightData(
    n = 120, maxDays = 15, graceEnd = 3, initiationHazard = 0.1, outcomeHazard = 0.02, deathHazard = 0.005,
    initiationCoefficients = rep(0.4, 5), seed = 1
  )
  on.exit(close(data))
  realShift <- .bootstrapReplicateShift
  local_mocked_bindings(.bootstrapReplicateShift = function(t0, replicates) {
    shift <- realShift(t0, replicates)
    shift$shifted[shift$estimand == "rd_15"] <- TRUE
    return(shift)
  })
  warnings <- capture_warnings(result <- suppressMessages(bootstrapCcwAnalysis(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 15, endAnchor = "cohort start"),
    createClonesArgs = createCreateClonesArgs(
      strategies = list(
        createTreatmentStrategy(1, "Initiate within 3 days", "initiateWithin", graceEnd = 3),
        createTreatmentStrategy(2, "Do not initiate within 3 days", "notInitiateWithin", graceEnd = 3)
      ),
      competingEvents = "none"
    ),
    fitCensoringModelsArgs = createFitCensoringModelsArgs(useCovariates = FALSE),
    computeWeightsArgs = createComputeWeightsArgs(),
    fitOutcomeModelArgs = createFitOutcomeModelArgs(modelType = "cox", horizons = 15, bootstrapArgs = createBootstrapArgs(replicates = 2, seed = 1))
  )))
  shift <- result$metaData$bootstrap$replicateShift
  expect_equal(shift, realShift(result$estimates, result$replicates) |> dplyr::mutate(shifted = c(NA, TRUE)))
  expect_equal(shift$estimand, c("hr", "rd_15"))
  shiftWarnings <- grep(shiftWarning, warnings, value = TRUE)
  expect_length(shiftWarnings, 1)
  expect_match(shiftWarnings, sprintf(
    "of rd_15 (%d replicates above it and %d below, sign test p = %.2g; median replicate %.2f bootstrap SD %s it)",
    shift$above[2], shift$below[2], shift$signTestP[2], abs(shift$medianShift[2]), if (shift$medianShift[2] >= 0) "above" else "below"
  ), fixed = TRUE)
  expect_match(shiftWarnings, "metaData$bootstrap$replicateShift", fixed = TRUE)
  expect_false(grepl("of hr (", shiftWarnings, fixed = TRUE))
})

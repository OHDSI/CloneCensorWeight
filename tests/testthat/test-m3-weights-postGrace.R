library(testthat)

# Contract section 11.2: with postGraceModel = "binary", a notInitiateWithin(G) clone's weight after G is
# 1 / (1 - pi(X)) from the binary model at G (the model of initiateWithin(G)); the initiation hazard model
# gives the weights on days up to G. With "deviationFirst", the binary weight applies from day G on.

postGraceNoTruncation <- function(stabilized = FALSE) {
  createComputeWeightsArgs(stabilized = stabilized, truncation = createWeightTruncation("none"))
}

# Seven persons under initiateWithin(2) vs notInitiateWithin(2), riskStart = 0, nominal risk end day 10.
# Person 7 has the outcome on day G = 2, so the weight of day G is a segment weight.
postGracePersons <- function() {
  data.frame(
    rowId = 1:7,
    personSeqId = 1:7,
    riskStart = 0L,
    riskEnd = c(10L, 10L, 10L, 6L, 10L, 10L, 10L),
    nominalRiskEnd = 10L,
    outcomeDay = c(NA, NA, 1L, NA, NA, NA, 2L),
    deathDay = c(NA, NA, NA, NA, NA, 8L, NA),
    treatmentDay = c(0L, 1L, NA, 2L, NA, 4L, NA),
    comparatorTreatmentDay = NA_integer_
  )
}

postGraceData <- function(persons = postGracePersons()) {
  createM3Data(list(
    persons = persons,
    covariates = data.frame(rowId = integer(0), covariateId = numeric(0), covariateValue = numeric(0))
  ))
}

test_that("notInitiateWithin weights after G equal hand-computed 1 / (1 - pi) under both orderings", {
  data <- postGraceData()
  on.exit(close(data))
  columns <- c("cloneRowId", "startTime", "stopTime", "weight", "outcome", "competing")

  # eventsFirst. Initiation hazards (daily indicators): day 0: 1 of 7 (persons 1-7); day 1: 1 of 5 (2, 4, 5,
  # 6, 7; 3 has the outcome on day 1); day 2: 1 of 3 (4, 5, 6; 7 has the outcome on day 2). So S(0) = 6/7
  # and S(1) = 24/35. Binary risk set at day 2: persons 1, 2, 4, 5, 6 (no event on or before day 2), of
  # whom 1, 2 and 4 started: pi = 3/5. The defer clones (8-14 for persons 1-7) are cut at the event
  # times t = 2 (day 1) and t = 3 (day 2) inside the window, and at the end of day tau + 1 = 3 (t = 4).
  # W(0) = 1, W(1) = 1 / S(0) = 7/6, W(2) = 1 / S(1) = 35/24, and W(d > 2) = 1 / (1 - 3/5) = 5/2.
  clones <- buildM3Clones(postGracePersons(), m3Strategies(graceEnd = 2), "eventsFirst")
  defer <- data.frame(
    cloneRowId = as.integer(c(8, 9, 10, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14)),
    startTime = c(0, 0, 0, 0, 2, 0, 2, 3, 4, 0, 2, 3, 4, 0, 2),
    stopTime = c(1, 2, 2, 2, 3, 2, 3, 4, 11, 2, 3, 4, 9, 2, 3),
    weight = c(1, 7 / 6, 7 / 6, 7 / 6, 35 / 24, 7 / 6, 35 / 24, 5 / 2, 5 / 2, 7 / 6, 35 / 24, 5 / 2, 5 / 2, 7 / 6, 35 / 24),
    outcome = as.integer(c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)),
    competing = as.integer(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))
  )
  binaryModels <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic", useCovariates = FALSE, postGraceModel = "binary")))
  hazardModels <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic", useCovariates = FALSE, postGraceModel = "hazard")))
  expect_null(hazardModels[[2]]$postGrace)
  expect_equal(binaryModels[[2]]$postGrace$predictions$propensityScore, rep(3 / 5, 5))
  binary <- collectSegments(suppressMessages(computeCensoringWeights(clones, binaryModels, postGraceNoTruncation())))
  hazard <- collectSegments(suppressMessages(computeCensoringWeights(clones, hazardModels, postGraceNoTruncation())))
  # Cyclops converges to about 1e-7 relative on the saturated numerator model:
  expect_equal(binary[binary$strategyId == 2, columns], defer, tolerance = 1e-6, ignore_attr = TRUE)
  # With "hazard", the weight after G is 1 / S(2) = 35/16; everything else is the same:
  expected <- defer
  expected$weight[expected$weight == 5 / 2] <- 35 / 16
  expect_equal(hazard[hazard$strategyId == 2, columns], expected, tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(binary[binary$strategyId == 1, ], hazard[hazard$strategyId == 1, ])
  # The Cox hazard model gives other weights up to G, and the same binary weight after G:
  coxModels <- suppressMessages(fitCensoringModels(clones, data, m3Args("cox", useCovariates = FALSE, postGraceModel = "binary")))
  cox <- collectSegments(suppressMessages(computeCensoringWeights(clones, coxModels, postGraceNoTruncation())))
  after <- cox$strategyId == 2 & cox$stopTime >= 4
  expect_equal(cox$weight[after], rep(5 / 2, sum(after)), tolerance = 1e-12)
  # Stabilized weights of models without covariates are exactly 1 (numerator 1 - pi-bar after G):
  stabilized <- collectSegments(suppressMessages(computeCensoringWeights(clones, binaryModels, postGraceNoTruncation(stabilized = TRUE))))
  expect_equal(stabilized$weight, rep(1, nrow(stabilized)), tolerance = 1e-6)
  close(clones)

  # deviationFirst. Initiation hazards: day 0: 1 of 7; day 1: 1 of 6 (3 is at risk on its outcome day);
  # day 2: 1 of 4 (4, 5, 6, 7). S(0) = 6/7, S(1) = 5/7, S(2) = 15/28. Binary risk set at day 2: persons
  # 1, 2, 4, 5, 6, 7 (no event before day 2), of whom 1, 2 and 4 started: pi = 1/2. A deviation on day d
  # comes before an outcome that day, so clone 8 (starts on day 0) has no segments, and 9 and 11 end
  # the day before they start. W(0) = 1 / S(0) = 7/6, W(1) = 1 / S(1) = 7/5, and from day G = 2 on
  # W = 1 / (1 - 1/2) = 2: also for person 7's outcome on day 2.
  clones <- buildM3Clones(postGracePersons(), m3Strategies(graceEnd = 2), "deviationFirst")
  defer <- data.frame(
    cloneRowId = as.integer(c(9, 10, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14)),
    startTime = c(0, 0, 0, 0, 2, 3, 4, 0, 2, 3, 4, 0, 2),
    stopTime = c(1, 2, 2, 2, 3, 4, 11, 2, 3, 4, 9, 2, 3),
    weight = c(7 / 6, 7 / 5, 7 / 5, 7 / 5, 2, 2, 2, 7 / 5, 2, 2, 2, 7 / 5, 2),
    outcome = as.integer(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)),
    competing = as.integer(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0))
  )
  binaryModels <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic", useCovariates = FALSE, postGraceModel = "binary")))
  hazardModels <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic", useCovariates = FALSE, postGraceModel = "hazard")))
  expect_equal(sort(binaryModels[[2]]$postGrace$predictions$rowId), c(1L, 2L, 4L, 5L, 6L, 7L))
  binary <- collectSegments(suppressMessages(computeCensoringWeights(clones, binaryModels, postGraceNoTruncation())))
  hazard <- collectSegments(suppressMessages(computeCensoringWeights(clones, hazardModels, postGraceNoTruncation())))
  expect_equal(binary[binary$strategyId == 2, columns], defer, tolerance = 1e-6, ignore_attr = TRUE)
  expected <- defer
  expected$weight[expected$weight == 2] <- 28 / 15
  expect_equal(hazard[hazard$strategyId == 2, columns], expected, tolerance = 1e-6, ignore_attr = TRUE)
  # The defer clones at risk on day G weigh as much as the binary risk set (6 persons):
  expect_equal(sum(binary$weight[binary$strategyId == 2 & binary$startTime < 3 & binary$stopTime >= 3]), 6, tolerance = 1e-6)
  close(clones)
})

test_that("The post-grace binary model is shared with initiateWithin(G), and fitted without it", {
  fixture <- createM3Fixture(150, riskStart = 0L, seed = 31)
  data <- createM3Data(fixture)
  on.exit(close(data))
  for (sameDayOrdering in c("eventsFirst", "deviationFirst")) {
    reference <- refBinaryRiskSet(fixture$persons, 7, sameDayOrdering)
    x <- m3CovariateMatrix(fixture, reference$rowId)
    fitted <- as.vector(fitted(glm(reference$y ~ x, family = binomial(), control = glm.control(epsilon = 1e-12))))
    propensityScore <- function(model) model$predictions$propensityScore[match(reference$rowId, model$predictions$rowId)]

    # initiateWithin(7) first: the notInitiateWithin strategy shares its binary model.
    clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = 7), sameDayOrdering)
    messages <- capture_messages(models <- fitCensoringModels(clones, data, m3Args("logistic", postGraceModel = "binary")))
    expect_equal(sum(grepl("^Fitting binary", messages)), 1)
    expect_true(any(grepl("Strategy 2 shares the binary model at day 7 of an earlier strategy", messages)))
    expect_identical(models[[2]]$postGrace$predictions, models[[1]]$predictions)
    expect_equal(models[[2]]$postGrace$coefficients, models[[1]]$coefficients)
    expect_setequal(models[[2]]$postGrace$predictions$rowId, reference$rowId)
    expect_equal(propensityScore(models[[2]]$postGrace), fitted, tolerance = 1e-6)
    expect_equal(models[[2]]$postGrace$numerator$propensityScore, mean(reference$y))
    close(clones)
    if (sameDayOrdering == "deviationFirst") {
      next
    }

    # notInitiateWithin first: initiateWithin shares the binary model fitted for its post-grace weights.
    clones <- buildM3Clones(fixture$persons, m3Strategies("notInitiateWithin", "initiateWithin", graceEnd = 7), sameDayOrdering)
    messages <- capture_messages(models <- fitCensoringModels(clones, data, m3Args("logistic", postGraceModel = "binary")))
    expect_true(any(grepl("binary censoring model at day 7 for the weights after the grace period of strategy 1", messages)))
    expect_true(any(grepl("Strategy 2 shares the censoring model of an earlier strategy", messages)))
    expect_identical(models[[1]]$postGrace$predictions, models[[2]]$predictions)
    close(clones)

    # Without an initiateWithin strategy, the binary model is fitted for the notInitiateWithin strategy:
    clones <- buildM3Clones(fixture$persons, m3Strategies("notInitiateWithin", "neverInitiate", graceEnd = 7), sameDayOrdering)
    models <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic", postGraceModel = "binary")))
    expect_equal(models[[1]]$postGrace$modelType, "binary")
    expect_equal(models[[1]]$postGrace$horizon, 7L)
    expect_equal(propensityScore(models[[1]]$postGrace), fitted, tolerance = 1e-6)
    # neverInitiate is unchanged:
    expect_null(models[[2]]$postGrace)
    close(clones)
  }
})

test_that("Post-grace binary weights match glm, and weights up to G match the hazard model", {
  skip_if_not_installed("survival")
  configurations <- list(
    list(riskStart = 0L, sameDayOrdering = "eventsFirst", hazardModelType = "logistic"),
    list(riskStart = 1L, sameDayOrdering = "deviationFirst", hazardModelType = "cox")
  )
  for (configuration in configurations) {
    sameDayOrdering <- configuration$sameDayOrdering
    offset <- if (sameDayOrdering == "eventsFirst") 1 else 0
    fixture <- createM3Fixture(150, riskStart = configuration$riskStart, seed = 32)
    data <- createM3Data(fixture)
    clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = 7), sameDayOrdering)
    reference <- refBinaryRiskSet(fixture$persons, 7, sameDayOrdering)
    x <- m3CovariateMatrix(fixture, reference$rowId)
    ps <- fitted(glm(reference$y ~ x, family = binomial(), control = glm.control(epsilon = 1e-12)))
    binaryModels <- suppressMessages(fitCensoringModels(clones, data, m3Args(configuration$hazardModelType, postGraceModel = "binary")))
    hazardModels <- suppressMessages(fitCensoringModels(clones, data, m3Args(configuration$hazardModelType, postGraceModel = "hazard")))
    for (stabilized in c(FALSE, TRUE)) {
      binary <- collectSegments(suppressMessages(computeCensoringWeights(clones, binaryModels, postGraceNoTruncation(stabilized))))
      hazard <- collectSegments(suppressMessages(computeCensoringWeights(clones, hazardModels, postGraceNoTruncation(stabilized))))
      expect_equal(binary[, c("cloneRowId", "startTime", "stopTime", "outcome", "competing")], hazard[, c("cloneRowId", "startTime", "stopTime", "outcome", "competing")])
      stopDay <- binary$stopTime + configuration$riskStart - 1
      postGrace <- binary$strategyId == 2 & stopDay - offset >= 7
      expect_true(sum(postGrace) > 50)
      expect_equal(binary$weight[!postGrace], hazard$weight[!postGrace], tolerance = 1e-10)
      expected <- 1 / (1 - ps[match(binary$rowId[postGrace], reference$rowId)])
      if (stabilized) {
        expected <- expected * (1 - mean(reference$y))
      }
      expect_equal(binary$weight[postGrace], as.vector(expected), tolerance = 1e-6)
      # The post-grace weights differ from 1 / S(G | X):
      expect_gt(max(abs(binary$weight[postGrace] - hazard$weight[postGrace])), 0.01)
    }
    close(clones)
    close(data)
  }
})

# The setting of the validation repro repro_hazardModelBalance.R: one binary X that changes when
# persons start (X = 1 starts on day 0 or never; X = 0 starts on days 1-3), and no
# deaths or outcomes before day 20, so everyone is under follow-up at G = 7.
createDay0Data <- function(n = 4000, seed = 1) {
  set.seed(seed)
  x <- rbinom(n, 1, 0.4)
  startDay <- ifelse(x == 1,
    ifelse(runif(n) < 0.63, 0L, NA_integer_),
    sample(c(1L, 2L, 3L, NA_integer_), n, replace = TRUE, prob = c(0.2, 0.2, 0.14, 0.46))
  )
  deathDay <- ifelse(runif(n) < 0.1, sample(20:300, n, replace = TRUE), NA_integer_)
  cohorts <- data.frame(
    rowId = seq_len(n), personSeqId = seq_len(n), personId = as.character(seq_len(n)),
    cohortStartDate = as.Date("2020-01-01"), daysFromObsStart = 400L, daysToCohortEnd = 0L,
    daysToObsEnd = 365L, daysToDeath = as.integer(deathDay)
  )
  started <- !is.na(startDay)
  data <- createCloneCensorWeightData(
    cohorts = cohorts,
    outcomes = data.frame(rowId = which(!is.na(deathDay)), outcomeId = 3, daysToEvent = deathDay[!is.na(deathDay)]),
    treatments = data.frame(rowId = which(started), treatmentId = 2, daysToStart = startDay[started], daysToEnd = startDay[started] + 30L),
    covariates = data.frame(rowId = which(x == 1), covariateId = 1001, covariateValue = 1),
    covariateRef = data.frame(covariateId = 1001, covariateName = "X", analysisId = 1, conceptId = 0),
    analysisRef = data.frame(
      analysisId = 1, analysisName = "X", domainId = "Demographics", startDay = 0, endDay = 0,
      isBinary = "Y", missingMeansZero = "Y"
    ),
    metaData = list(targetId = 1, treatmentIds = c(treatment = 2), outcomeIds = 3)
  )
  return(list(data = data, x = x))
}

test_that("With postGraceModel = 'binary', the weighted comparator at tau + 1 reproduces the population (regression)", {
  day0 <- createDay0Data()
  on.exit(close(day0$data))
  population <- suppressMessages(createStudyPopulation(day0$data, outcomeId = 3, createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
    riskWindowStart = 0, riskWindowEnd = 365, endAnchor = "cohort start"
  )))
  for (sameDayOrdering in c("eventsFirst", "deviationFirst")) {
    clones <- suppressMessages(createClones(population, day0$data, createCreateClonesArgs(
      strategies = list(
        createTreatmentStrategy(1, "Start within 7 days", "initiateWithin", graceEnd = 7),
        createTreatmentStrategy(2, "No start within 7 days", "notInitiateWithin", graceEnd = 7)
      ),
      competingEvents = "none",
      sameDayOrdering = sameDayOrdering
    )))
    comparatorMean <- list()
    for (postGraceModel in c("binary", "hazard")) {
      models <- suppressMessages(fitCensoringModels(clones, day0$data, createFitCensoringModelsArgs(prior = Cyclops::createPrior("none"), postGraceModel = postGraceModel)))
      weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(truncation = createWeightTruncation("none"))))
      balance <- suppressMessages(computeCovariateBalance(weighted, day0$data))
      comparatorMean[[postGraceModel]] <- balance$afterMatchingMeanComparator
      expect_equal(balance$afterMatchingMeanTarget, mean(day0$x), tolerance = 1e-6)
      close(weighted)
    }
    # The binary model is saturated in X, so 1 / (1 - pi(X)) reweights the non-starters to everyone under
    # follow-up at G, here the whole population:
    expect_equal(comparatorMean$binary, mean(day0$x), tolerance = 1e-6)
    # 1 / S(7 | X) of the hazard model, which assumes the same effect of X on every day, does not:
    expect_gt(abs(comparatorMean$hazard - mean(day0$x)), 0.05)
    close(clones)
  }
})

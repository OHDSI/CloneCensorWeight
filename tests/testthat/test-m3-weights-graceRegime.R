library(testthat)

# Plan section 2.2 and decision 9 item 7 (Phase 2 item 6): initiateWithin strategies with graceStart > 0
# (a start before graceStart is a deviation, weighted by the initiation hazard of the days before
# graceStart), and the natural-grace (last-interval) weight, graceRegime = "naturalGrace".

gNoTruncation <- function(stabilized = FALSE) {
  createComputeWeightsArgs(stabilized = stabilized, truncation = createWeightTruncation("none"))
}

gArgs <- function(initiateWithinModel = "binary", graceRegime = "binaryAtG", hazardModelType = "logistic", ...) {
  args <- m3Args(hazardModelType, ...)
  args$initiateWithinModel <- initiateWithinModel
  args$graceRegime <- graceRegime
  return(args)
}

# The seven persons of test-m3-weights-initiateWithinModel.R (riskStart = 0, nominal risk end day 10, no
# covariates): starts on days 0 (1), 1 (2), 2 (4) and 4 (6); outcomes on days 1 (3) and 2 (7); death on
# day 8 (6); observation ends on day 6 (4).
gPersons <- function() {
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

gData <- function() {
  createM3Data(list(
    persons = gPersons(),
    covariates = data.frame(rowId = integer(0), covariateId = numeric(0), covariateValue = numeric(0))
  ))
}

gStrategies <- function(graceStart = 0, graceEnd = 2, second = createTreatmentStrategy(2, "Not within 2", "notInitiateWithin", graceEnd = 2)) {
  list(createTreatmentStrategy(1, "Initiate", "initiateWithin", graceStart = graceStart, graceEnd = graceEnd), second)
}

gColumns <- c("cloneRowId", "startTime", "stopTime", "weight", "outcome")

# The weights of the initiate clones of gPersons() with graceStart = 1, G = 2, eventsFirst, by clone,
# from the weight before the grace period (`early`, W on days 1..2) and after it (`after1`, `after2` of
# persons 2 and 4). Person 1 starts on day 0 and deviates; 3 and 7 have the outcome on days 1 and 2; 5
# and 6 deviate on day 2. Cuts: end of day G = 2 (t = 3), day tau + 1 = 3 (t = 4), and the ends.
gSegments <- function(early, after2, after4) {
  data.frame(
    cloneRowId = as.integer(c(1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 7)),
    startTime = c(0, 0, 3, 4, 0, 0, 3, 4, 0, 0, 0),
    stopTime = c(1, 3, 4, 11, 2, 3, 4, 7, 3, 3, 3),
    weight = c(1, early, after2, after2, early, early, after4, after4, early, early, early),
    outcome = as.integer(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1))
  )
}

test_that("graceStart > 0: weights before the grace period from the hazard of the days before it, and after G", {
  data <- gData()
  clones <- buildM3Clones(gPersons(), gStrategies(graceStart = 1), "eventsFirst")
  on.exit({
    close(clones)
    close(data)
  })
  # The early hazard model through day 0: person 1 of 7 starts on day 0, S(0) = 6/7, so W = 7/6 on days
  # 1..2. graceRegime = "binaryAtG": the binary model at G = 2 is fit on persons 2, 4, 5 and 6 (under
  # follow-up at G, not started before day 1), of whom 2 and 4 started by G: pi = 1/2, W = 7/6 * 2 = 7/3
  # after G. The numerator is the same product (no covariates), so the stabilized weights are 1.
  models <- suppressMessages(fitCensoringModels(clones, data, gArgs()))
  expect_equal(models[[1]]$kind, "binary")
  expect_equal(models[[1]]$firstDay, 1L)
  expect_equal(models[[1]]$graceStart, 1L)
  expect_equal(models[[1]]$graceRegime, "binaryAtG")
  expect_equal(models[[1]]$diagnostics$riskSetSize, 4)
  expect_equal(models[[1]]$early$kind, "hazard")
  expect_equal(models[[1]]$early$horizon, 0L)
  segments <- collectSegments(suppressMessages(computeCensoringWeights(clones, models, gNoTruncation())))
  expect_equal(segments[segments$strategyId == 1, gColumns], gSegments(7 / 6, 7 / 3, 7 / 3), tolerance = 1e-6, ignore_attr = TRUE)
  stabilized <- collectSegments(suppressMessages(computeCensoringWeights(clones, models, gNoTruncation(stabilized = TRUE))))
  expect_equal(stabilized$weight[stabilized$strategyId == 1], rep(1, 11), tolerance = 1e-6)

  # graceRegime = "naturalGrace": the binary model is fit on persons 4, 5 and 6 (untreated at G), of whom
  # 4 starts on day G: p = 1/3. Person 2 started before G: W = 7/6 after G; person 4: 7/6 * 3 = 7/2.
  # The numerator after G is S-bar(0) times the share who started by G of the persons of the
  # "binaryAtG" risk set, 6/7 * 1/2 = 3/7, the same for every clone.
  models <- suppressMessages(fitCensoringModels(clones, data, gArgs(graceRegime = "naturalGrace")))
  expect_equal(models[[1]]$firstDay, 2L)
  expect_equal(models[[1]]$diagnostics$riskSetSize, 3)
  expect_equal(models[[1]]$diagnostics$events, 1)
  segments <- collectSegments(suppressMessages(computeCensoringWeights(clones, models, gNoTruncation())))
  expect_equal(segments[segments$strategyId == 1, gColumns], gSegments(7 / 6, 7 / 6, 7 / 2), tolerance = 1e-6, ignore_attr = TRUE)
  stabilized <- collectSegments(suppressMessages(computeCensoringWeights(clones, models, gNoTruncation(stabilized = TRUE))))
  expected <- gSegments(1, 3 / 7 * 7 / 6, 3 / 7 * 7 / 2)
  expect_equal(stabilized$weight[stabilized$strategyId == 1], expected$weight, tolerance = 1e-6)
  expect_equal(stabilized$otherScaleWeight[stabilized$strategyId == 1], gSegments(7 / 6, 7 / 6, 7 / 2)$weight, tolerance = 1e-6)

  # initiateWithinModel = "hazard": the hazard model through G gives everything. Daily hazards 1/7 (day
  # 0), 1/5 (day 1) and 1/3 (day 2), so S(0) = 6/7, S(1) = 24/35 and S(2) = 16/35. "binaryAtG": after G,
  # W = 1 / (S(0) - S(2)) = 35/14 = 5/2. "naturalGrace": p = 1 - S(2) / S(1) = 1/3, as the binary model.
  for (graceRegime in c("binaryAtG", "naturalGrace")) {
    models <- suppressMessages(fitCensoringModels(clones, data, gArgs(initiateWithinModel = "hazard", graceRegime = graceRegime)))
    expect_equal(models[[1]]$kind, "hazard")
    expect_null(models[[1]]$early)
    segments <- collectSegments(suppressMessages(computeCensoringWeights(clones, models, gNoTruncation())))
    expected <- if (graceRegime == "binaryAtG") gSegments(7 / 6, 5 / 2, 5 / 2) else gSegments(7 / 6, 7 / 6, 7 / 2)
    expect_equal(segments[segments$strategyId == 1, gColumns], expected, tolerance = 1e-6, ignore_attr = TRUE)
    # The comparator (notInitiateWithin(2)) shares this hazard model:
    expect_equal(models[[2]][c("coefficients", "intervals")], models[[1]][c("coefficients", "intervals")])
  }

  # computeWeightDistribution() evaluates the same weights on any day (day 1 and day 5):
  models <- suppressMessages(fitCensoringModels(clones, data, gArgs(graceRegime = "naturalGrace")))
  weighted <- suppressMessages(computeCensoringWeights(clones, models, gNoTruncation()))
  distribution <- computeWeightDistribution(weighted, days = c(0, 1, 5))
  close(weighted)
  initiate <- distribution[distribution$strategyId == 1, ]
  expect_equal(initiate$clones, c(7L, 6L, 2L))
  expect_equal(initiate$maxWeight, c(1, 7 / 6, 7 / 2), tolerance = 1e-6)
  expect_equal(initiate$meanWeight[3], (7 / 6 + 7 / 2) / 2, tolerance = 1e-6)
})

test_that("graceStart > 0 with deviationFirst: the switch at G applies from day G", {
  data <- gData()
  clones <- buildM3Clones(gPersons(), gStrategies(graceStart = 1), "deviationFirst")
  on.exit({
    close(clones)
    close(data)
  })
  # deviationFirst: person 1 deviates on day 0 before any follow-up (no clone). The early model through
  # day 0 has all 7 at risk: S(0) = 6/7, and W(d) = 1 / S(min(d, 0)) = 7/6 from day 0 on. The binary model
  # at G = 2 is fit on the persons with no event before day 2 and not started before day 1: 2, 4, 5, 6
  # and 7 (outcome on day 2), of whom 2 and 4 started: pi = 2/5, so W = 7/6 * 5/2 = 35/12 from day G on.
  # Clones 5, 6 and 7 (persons 5, 6, 7) deviate on day 2 and end on day 1.
  models <- suppressMessages(fitCensoringModels(clones, data, gArgs()))
  segments <- collectSegments(suppressMessages(computeCensoringWeights(clones, models, gNoTruncation())))
  initiate <- segments[segments$strategyId == 1, ]
  expect_equal(initiate$weight[initiate$stopTime <= 2], rep(7 / 6, sum(initiate$stopTime <= 2)), tolerance = 1e-6)
  expect_equal(unique(round(initiate$weight[initiate$startTime >= 2], 6)), round(35 / 12, 6))
})

test_that("Natural-grace weights with graceStart = 0 by hand, and KM risks invariant to stabilization", {
  data <- gData()
  clones <- buildM3Clones(gPersons(), gStrategies(), "eventsFirst")
  on.exit({
    close(clones)
    close(data)
  })
  # G = 2. Persons 4, 5 and 6 reach day 2 untreated and under follow-up; 4 starts on day 2: p = 1/3
  # (binary model, and the day-2 hazard of the hazard model, 1 of 3). Clones 1 and 2 started before G:
  # W = 1; clone 4: W = 3 after G. Numerators: the share who started by G of the persons under follow-up
  # at G (1, 2, 4, 5, 6), 3/5, with the binary model; 1 - S-bar(2) = 19/35 with the hazard model.
  expectedAfter <- c(`1` = 1, `2` = 1, `4` = 3)
  for (initiateWithinModel in c("binary", "hazard")) {
    models <- suppressMessages(fitCensoringModels(clones, data, gArgs(initiateWithinModel = initiateWithinModel, graceRegime = "naturalGrace")))
    weighted <- suppressMessages(computeCensoringWeights(clones, models, gNoTruncation()))
    segments <- collectSegments(weighted)
    after <- segments[segments$strategyId == 1 & segments$startTime >= 3, ]
    expect_equal(after$weight, unname(expectedAfter[as.character(after$cloneRowId)]), tolerance = 1e-6)
    expect_equal(segments$weight[segments$strategyId == 1 & segments$stopTime <= 3], rep(1, 7))
    numerator <- if (initiateWithinModel == "binary") 3 / 5 else 19 / 35
    expect_equal(unique(round(after$stabilizationNumerator, 6)), round(numerator, 6))
    stabilizedWeighted <- suppressMessages(computeCensoringWeights(clones, models, gNoTruncation(stabilized = TRUE)))
    km <- function(w) fitOutcomeModel(w, createFitOutcomeModelArgs(modelType = "km", horizons = c(1, 5, 9)))$estimates
    expect_equal(km(stabilizedWeighted)$estimate, km(weighted)$estimate, tolerance = 1e-8)
    close(weighted)
    close(stabilizedWeighted)
  }
})

test_that("Models of the days before graceStart are shared with a notInitiateWithin(graceStart - 1) strategy and reused", {
  data <- gData()
  strategies <- gStrategies(graceStart = 2, graceEnd = 4, second = createTreatmentStrategy(2, "Not within 1", "notInitiateWithin", graceEnd = 1))
  clones <- buildM3Clones(gPersons(), strategies, "eventsFirst")
  on.exit({
    close(clones)
    close(data)
  })
  expect_message(models <- fitCensoringModels(clones, data, gArgs()), "Strategy 2 shares the censoring model of an earlier strategy")
  expect_equal(models[[1]]$early$horizon, 1L)
  expect_equal(models[[1]]$early[c("coefficients", "intervals")], models[[2]][c("coefficients", "intervals")])
  # The fits of these models are reused, and "naturalGrace" refits only its binary model:
  fits <- .censoringFitsFromModels(models)
  expect_true(any(startsWith(names(fits), "binaryFrom2 ")))
  messages <- testthat::capture_messages(refit <- .fitCensoringModelsWithFits(clones, data, gArgs(graceRegime = "naturalGrace"), fits = fits)$models)
  expect_true(any(grepl("Fitting binary censoring model at day 4 of the persons who have not started before day 4", messages, fixed = TRUE)))
  expect_false(any(grepl("Fitting logistic initiation hazard model", messages, fixed = TRUE)))
  expect_equal(refit[[1]]$early, models[[1]]$early)
  # The bootstrap's fixed prior variance of the early model is accepted:
  expect_s3_class(suppressMessages(fitCensoringModels(clones, data, gArgs(), .fixedPriorVariances = list("1.early" = 1))), "CcwCensoringModels")
})

test_that("The natural-grace binary model is not shared with the binary model at G of other strategies", {
  data <- gData()
  clones <- buildM3Clones(gPersons(), gStrategies(), "eventsFirst")
  on.exit({
    close(clones)
    close(data)
  })
  # postGraceModel = "binary" fits the binary model at G = 2 for the notInitiateWithin strategy (persons
  # 1, 2, 4, 5 and 6 under follow-up at G); the natural-grace model of the initiateWithin strategy is fit
  # on the persons untreated at G (4, 5 and 6):
  messages <- capture_messages(models <- fitCensoringModels(clones, data, gArgs(graceRegime = "naturalGrace", postGraceModel = "binary")))
  expect_equal(sum(grepl("^Fitting binary", messages)), 2)
  expect_equal(models[[1]]$diagnostics$riskSetSize, 3)
  expect_equal(models[[2]]$postGrace$diagnostics$riskSetSize, 5)
  # Fits of "binaryAtG" are not reused for "naturalGrace", and the other way round:
  fits <- .censoringFitsFromModels(suppressMessages(fitCensoringModels(clones, data, gArgs())))
  messages <- capture_messages(natural <- .fitCensoringModelsWithFits(clones, data, gArgs(graceRegime = "naturalGrace"), fits = fits)$models)
  expect_true(any(grepl("^Fitting binary censoring model at day 2 of the persons who have not started before day 2", messages)))
  expect_equal(natural[[1]]$diagnostics$riskSetSize, 3)
  messages <- capture_messages(atG <- .fitCensoringModelsWithFits(clones, data, gArgs(), fits = .censoringFitsFromModels(natural))$models)
  expect_true(any(grepl("^Fitting binary censoring model at day 2 for strategy 1", messages)))
  expect_equal(atG[[1]]$diagnostics$riskSetSize, 5)
})

test_that("The hazard model's interval grid must start an interval at graceStart, and at G with naturalGrace", {
  data <- gData()
  on.exit(close(data))
  weekly <- createTimeBasis(type = "indicator", intervalLength = 2)
  clones <- buildM3Clones(gPersons(), gStrategies(graceStart = 1, graceEnd = 2), "eventsFirst")
  on.exit(close(clones), add = TRUE)
  expect_error(fitCensoringModels(clones, data, gArgs(initiateWithinModel = "hazard", timeBasis = weekly)), "graceStart of strategy 'Initiate' \\(1\\) must be a multiple")
  clones3 <- buildM3Clones(gPersons(), gStrategies(graceStart = 0, graceEnd = 3), "eventsFirst")
  on.exit(close(clones3), add = TRUE)
  expect_error(fitCensoringModels(clones3, data, gArgs(initiateWithinModel = "hazard", graceRegime = "naturalGrace", timeBasis = weekly)), "graceEnd of strategy 'Initiate' \\(3\\) must be a multiple")
  # The binary model at G has no grid, and "binaryAtG" with graceStart = 0 needs none (nobody starts on
  # day 3, so "naturalGrace" warns):
  expect_warning(models <- suppressMessages(fitCensoringModels(clones3, data, gArgs(graceRegime = "naturalGrace", timeBasis = weekly))), "none of the 2 persons")
  expect_s3_class(models, "CcwCensoringModels")
  expect_s3_class(suppressMessages(fitCensoringModels(clones3, data, gArgs(initiateWithinModel = "hazard", timeBasis = weekly))), "CcwCensoringModels")
})

test_that("Models fit for another graceStart do not match the clones", {
  data <- gData()
  clones0 <- buildM3Clones(gPersons(), gStrategies(graceStart = 0), "eventsFirst")
  clones1 <- buildM3Clones(gPersons(), gStrategies(graceStart = 1), "eventsFirst")
  on.exit({
    close(clones0)
    close(clones1)
    close(data)
  })
  models <- suppressMessages(fitCensoringModels(clones0, data, gArgs()))
  expect_error(computeCensoringWeights(clones1, models, gNoTruncation()), "does not match the strategy")
})

test_that("Calibration and positivity report the models of graceStart > 0 and naturalGrace", {
  data <- gData()
  clones <- buildM3Clones(gPersons(), gStrategies(graceStart = 1), "eventsFirst")
  on.exit({
    close(clones)
    close(data)
  })
  for (graceRegime in c("binaryAtG", "naturalGrace")) {
    for (initiateWithinModel in c("binary", "hazard")) {
      # (The default tier runs binaryAtG with the hazard model and naturalGrace with the binary model.)
      if (!isFullTestRun() && (graceRegime == "naturalGrace") == (initiateWithinModel == "hazard")) next
      models <- suppressMessages(fitCensoringModels(clones, data, gArgs(initiateWithinModel = initiateWithinModel, graceRegime = graceRegime)))
      calibration <- suppressMessages(computeCensoringModelCalibration(clones, models))
      overall <- calibration[calibration$grouping == "overall" & calibration$strategyId == 1, ]
      expect_equal(overall$model, initiateWithinModel)
      expect_equal(overall$graceRegime, graceRegime)
      # The persons not started before day 1 (binaryAtG: 2, 4, 5, 6, of whom 2 started by G) or before
      # day G (naturalGrace: 4, 5, 6, of whom 4 started on G):
      expect_equal(overall$persons, if (graceRegime == "binaryAtG") 4L else 3L)
      expect_equal(overall$observed, if (graceRegime == "binaryAtG") 1 / 2 else 1 / 3)
      expectedPredicted <- if (graceRegime == "naturalGrace") 1 / 3 else if (initiateWithinModel == "binary") 1 / 2 else (6 / 7 - 16 / 35) / (6 / 7)
      expect_equal(overall$predicted, expectedPredicted, tolerance = 1e-6)
      expect_true(all(is.na(calibration$graceRegime[calibration$strategyId == 2])))
      positivity <- computePositivity(clones, models)
      initiate <- positivity[positivity$strategyId == 1, ]
      # The model at G and the probability of starting before graceStart (day 0: 1/7):
      expect_equal(initiate$day, c(2L, 0L))
      expect_equal(initiate$graceRegime, rep(graceRegime, 2))
      expect_equal(initiate$mean[2], 1 / 7, tolerance = 1e-6)
      skip_if_not_installed("ggplot2")
      expect_s3_class(plotPositivity(positivity), "ggplot")
      expect_s3_class(plotCensoringModelCalibration(calibration), "ggplot")
    }
  }
})

test_that("graceStart > 0 weights by hand differ by covariate stratum under both regimes", {
  # gPersons() (X = 0, rowIds 1-7) plus seven persons with X = 1 (covariate 1001, rowIds 8-14): starts
  # on days 0 (8, 9), 1 (10) and 2 (11); 12, 13 and 14 do not start, and 14 has the outcome on day 5.
  # With one binary covariate and no prior, the binary models and the early hazard model through day 0
  # (one interval) are saturated, so they give the empirical shares and hazards of each stratum.
  persons <- rbind(gPersons(), data.frame(
    rowId = 8:14,
    personSeqId = 8:14,
    riskStart = 0L,
    riskEnd = 10L,
    nominalRiskEnd = 10L,
    outcomeDay = c(NA, NA, NA, NA, NA, NA, 5L),
    deathDay = NA_integer_,
    treatmentDay = c(0L, 0L, 1L, 2L, NA, NA, NA),
    comparatorTreatmentDay = NA_integer_
  ))
  data <- createM3Data(list(
    persons = persons,
    covariates = data.frame(rowId = 8:14, covariateId = 1001, covariateValue = 1)
  ))
  clones <- buildM3Clones(persons, gStrategies(graceStart = 1), "eventsFirst")
  on.exit({
    close(clones)
    close(data)
  })
  # Early hazard model (day 0): X = 0: 1 of 7 starts, W = 7/6 on days 1..2; X = 1: 2 of 7, W = 7/5.
  # "binaryAtG", persons under follow-up at G = 2 not started before day 1: X = 0: 2, 4, 5, 6 (2 and 4
  # started), pi = 1/2, W = 7/3 after G; X = 1: 10-14 (10 and 11 started), pi = 2/5, W = 7/2.
  # "naturalGrace", persons untreated at G: X = 0: 4, 5, 6 (4 starts on day 2), p = 1/3; X = 1: 11-14
  # (11 starts on day 2), p = 1/4. After G, W = 7/6 (clone 2) and 7/2 (clone 4); 7/5 (clone 10) and
  # 28/5 (clone 11). The clones that start on day 0 (1, 8, 9) deviate on day 0.
  early <- c(1, rep(7 / 6, 6), 1, 1, rep(7 / 5, 5))
  after <- list(
    binaryAtG = c(`2` = 7 / 3, `4` = 7 / 3, `10` = 7 / 2, `11` = 7 / 2),
    naturalGrace = c(`2` = 7 / 6, `4` = 7 / 2, `10` = 7 / 5, `11` = 28 / 5)
  )
  for (graceRegime in names(after)) {
    models <- suppressMessages(fitCensoringModels(clones, data, gArgs(graceRegime = graceRegime)))
    segments <- collectSegments(suppressMessages(computeCensoringWeights(clones, models, gNoTruncation())))
    initiate <- segments[segments$strategyId == 1, ]
    before <- initiate[initiate$startTime < 3, ]
    expect_equal(before$weight, early[before$rowId], tolerance = 1e-6)
    afterG <- initiate[initiate$startTime >= 3, ]
    expect_setequal(as.character(afterG$rowId), names(after[[graceRegime]]))
    expect_equal(afterG$weight, unname(after[[graceRegime]][as.character(afterG$rowId)]), tolerance = 1e-6)
  }
})

test_that("graceStart > 0 and naturalGrace with the Cox hazard model: Breslow survival by hand", {
  data <- gData()
  clones <- buildM3Clones(gPersons(), gStrategies(graceStart = 1), "eventsFirst")
  on.exit({
    close(clones)
    close(data)
  })
  # Breslow cumulative hazards without covariates: 1/7 (day 0), + 1/5 (day 1), + 1/3 (day 2), and
  # S(d) = exp(-Lambda(d)). The early model through day 0 is that of a notInitiateWithin(0) strategy:
  # W = exp(1/7) on days 1..2.
  notWithin0 <- buildM3Clones(gPersons(), gStrategies(graceStart = 1, second = createTreatmentStrategy(2, "Not within 0", "notInitiateWithin", graceEnd = 0)), "eventsFirst")
  on.exit(close(notWithin0), add = TRUE)
  early <- exp(1 / 7)
  for (initiateWithinModel in c("binary", "hazard")) {
    for (graceRegime in c("binaryAtG", "naturalGrace")) {
      # (The default tier runs the binary model with binaryAtG and the hazard model with naturalGrace.)
      if (!isFullTestRun() && (initiateWithinModel == "hazard") != (graceRegime == "naturalGrace")) next
      args <- gArgs(initiateWithinModel = initiateWithinModel, graceRegime = graceRegime, hazardModelType = "cox")
      models <- suppressMessages(fitCensoringModels(clones, data, args))
      segments <- collectSegments(suppressMessages(computeCensoringWeights(clones, models, gNoTruncation())))
      initiate <- segments[segments$strategyId == 1, ]
      # After G: "binaryAtG" with the binary model, W = early / pi = early * 2; with the hazard model,
      # 1 / (S(0) - S(2)). "naturalGrace": clone 2 (started on day 1) W = early, clone 4 (day 2) W =
      # early / p with p = 1/3 (binary model) or 1 - S(2) / S(1) = 1 - exp(-1/3) (hazard model).
      after <- if (graceRegime == "binaryAtG") {
        rep(if (initiateWithinModel == "binary") 2 * early else 1 / (exp(-1 / 7) - exp(-(1 / 7 + 1 / 5 + 1 / 3))), 2)
      } else {
        c(early, early / (if (initiateWithinModel == "binary") 1 / 3 else 1 - exp(-1 / 3)))
      }
      expected <- gSegments(early, after[1], after[2])
      expect_equal(initiate[, gColumns], expected, tolerance = 1e-6, ignore_attr = TRUE)
    }
  }
  # The early weight equals the weight of a notInitiateWithin(0) strategy after day 0:
  models <- suppressMessages(fitCensoringModels(notWithin0, data, gArgs(hazardModelType = "cox")))
  segments <- collectSegments(suppressMessages(computeCensoringWeights(notWithin0, models, gNoTruncation())))
  comparator <- segments[segments$strategyId == 2 & segments$startTime >= 1, ]
  expect_equal(unique(round(comparator$weight, 8)), round(early, 8))
})

test_that("naturalGrace warns when nobody starts on day G", {
  data <- gData()
  # G = 3: persons 5 and 6 reach day 3 untreated and under follow-up, and nobody starts on day 3 (6
  # starts on day 4), so nobody represents them after G.
  clones <- buildM3Clones(gPersons(), gStrategies(graceEnd = 3, second = createTreatmentStrategy(2, "Not within 3", "notInitiateWithin", graceEnd = 3)), "eventsFirst")
  on.exit({
    close(clones)
    close(data)
  })
  for (initiateWithinModel in c("binary", "hazard")) {
    warnings <- capture_warnings(suppressMessages(fitCensoringModels(clones, data, gArgs(initiateWithinModel = initiateWithinModel, graceRegime = "naturalGrace"))))
    expect_length(warnings, 1)
    expect_match(warnings, "none of the 2 persons who reach day 3 (graceEnd) of strategy 1 (Initiate) untreated and under follow-up starts on day 3", fixed = TRUE)
    expect_match(warnings, "computePositivity()", fixed = TRUE)
    expect_length(capture_warnings(suppressMessages(fitCensoringModels(clones, data, gArgs(initiateWithinModel = initiateWithinModel)))), 0)
  }
  # Callers that fit several variants warn once:
  expect_length(capture_warnings(suppressMessages(.fitCensoringModelsWithFits(clones, data, gArgs(graceRegime = "naturalGrace"), warnNaturalGrace = FALSE))), 0)
})

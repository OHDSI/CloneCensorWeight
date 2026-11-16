library(testthat)

# Contract section 13.1: with initiateWithinModel = "hazard", an initiateWithin(G) clone's weight after G
# is 1 / (1 - S(G | X)) from the initiation hazard model through G, and its stabilization numerator is
# 1 - S-bar(G). With "deviationFirst", the weight applies from day G on.

iwNoTruncation <- function(stabilized = FALSE) {
  createComputeWeightsArgs(stabilized = stabilized, truncation = createWeightTruncation("none"))
}

iwArgs <- function(hazardModelType = "logistic", initiateWithinModel = "hazard", ...) {
  args <- m3Args(hazardModelType, ...)
  args$initiateWithinModel <- initiateWithinModel
  return(args)
}

# The seven persons of test-m3-weights-postGrace.R: initiateWithin(2) vs notInitiateWithin(2),
# riskStart = 0, nominal risk end day 10, no covariates.
iwPersons <- function() {
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

iwData <- function() {
  createM3Data(list(
    persons = iwPersons(),
    covariates = data.frame(rowId = integer(0), covariateId = numeric(0), covariateValue = numeric(0))
  ))
}

test_that("initiateWithin weights after G equal hand-computed 1 / (1 - S(G)) under both orderings", {
  data <- iwData()
  on.exit(close(data))
  columns <- c("cloneRowId", "startTime", "stopTime", "weight", "stabilizationNumerator", "outcome", "competing")

  # eventsFirst. Initiation hazards (daily indicators): day 0: 1 of 7; day 1: 1 of 5 (3 has the outcome
  # on day 1); day 2: 1 of 3 (4, 5, 6; 7 has the outcome on day 2). S(2) = 6/7 * 4/5 * 2/3 = 16/35, so
  # 1 - S(2) = 19/35, against pi = 3/5 of the binary model at G (persons 1, 2, 4, 5, 6 under follow-up
  # at G, of whom 1, 2 and 4 started). The initiate clones (1-7) are cut at the end of day G = 2 (t = 3)
  # and of day tau + 1 = 3 (t = 4). Clones 1, 2 and 4 started by G and are followed after it with
  # W = 35/19; 3 and 7 end with the outcome, and 5 and 6 deviate on day G.
  clones <- buildM3Clones(iwPersons(), m3Strategies(graceEnd = 2), "eventsFirst")
  initiate <- data.frame(
    cloneRowId = as.integer(c(1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 7)),
    startTime = c(0, 3, 4, 0, 3, 4, 0, 0, 3, 4, 0, 0, 0),
    stopTime = c(3, 4, 11, 3, 4, 11, 2, 3, 4, 7, 3, 3, 3),
    weight = c(1, 35 / 19, 35 / 19, 1, 35 / 19, 35 / 19, 1, 1, 35 / 19, 35 / 19, 1, 1, 1),
    stabilizationNumerator = c(1, 19 / 35, 19 / 35, 1, 19 / 35, 19 / 35, 1, 1, 19 / 35, 19 / 35, 1, 1, 1),
    outcome = as.integer(c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1)),
    competing = 0L
  )
  hazardModels <- suppressMessages(fitCensoringModels(clones, data, iwArgs()))
  binaryModels <- suppressMessages(fitCensoringModels(clones, data, iwArgs(initiateWithinModel = "binary")))
  expect_equal(hazardModels[[1]]$kind, "hazard")
  expect_equal(hazardModels[[1]]$type, "initiateWithin")
  expect_equal(hazardModels[[1]]$horizon, 2L)
  hazard <- collectSegments(suppressMessages(computeCensoringWeights(clones, hazardModels, iwNoTruncation())))
  binary <- collectSegments(suppressMessages(computeCensoringWeights(clones, binaryModels, iwNoTruncation())))
  expect_equal(hazard[hazard$strategyId == 1, columns], initiate, tolerance = 1e-6, ignore_attr = TRUE)
  # The binary model gives 1 / pi = 5/3 after G; the segments are the same:
  expected <- initiate
  expected$weight[expected$weight > 1] <- 5 / 3
  expected$stabilizationNumerator[expected$stabilizationNumerator < 1] <- 3 / 5
  expect_equal(binary[binary$strategyId == 1, columns], expected, tolerance = 1e-6, ignore_attr = TRUE)
  # The notInitiateWithin strategy is unchanged:
  expect_equal(hazard[hazard$strategyId == 2, ], binary[binary$strategyId == 2, ])
  # Stabilized weights of models without covariates are exactly 1 (numerator 1 - S-bar(G) = 1 - S(G)):
  stabilized <- collectSegments(suppressMessages(computeCensoringWeights(clones, hazardModels, iwNoTruncation(stabilized = TRUE))))
  expect_equal(stabilized$weight, rep(1, nrow(stabilized)), tolerance = 1e-6)
  expect_equal(stabilized$otherScaleWeight[stabilized$strategyId == 1], initiate$weight, tolerance = 1e-6)
  # computeWeightDistribution() reports the same weight on a day inside a segment (day 5 of clone 1):
  weighted <- suppressMessages(computeCensoringWeights(clones, hazardModels, iwNoTruncation()))
  distribution <- computeWeightDistribution(weighted, days = c(1, 5))
  expect_equal(distribution$meanWeight[distribution$strategyId == 1], c(1, 35 / 19), tolerance = 1e-6)
  close(weighted)
  close(clones)

  # deviationFirst. Initiation hazards: day 0: 1 of 7; day 1: 1 of 6 (3 is at risk on its outcome day);
  # day 2: 1 of 4 (4, 5, 6, 7). S(2) = 15/28, so 1 - S(2) = 13/28 (pi = 1/2 of the binary model). The
  # clones are cut at the end of day G - 1 = 1 (t = 2), and the weight 28/13 applies from day G on.
  # Clones 5, 6 and 7 deviate on day 2 before any outcome that day, so they end on day 1.
  clones <- buildM3Clones(iwPersons(), m3Strategies(graceEnd = 2), "deviationFirst")
  initiate <- data.frame(
    cloneRowId = as.integer(c(1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 7)),
    startTime = c(0, 2, 4, 0, 2, 4, 0, 0, 2, 4, 0, 0, 0),
    stopTime = c(2, 4, 11, 2, 4, 11, 2, 2, 4, 7, 2, 2, 2),
    weight = c(1, 28 / 13, 28 / 13, 1, 28 / 13, 28 / 13, 1, 1, 28 / 13, 28 / 13, 1, 1, 1),
    stabilizationNumerator = c(1, 13 / 28, 13 / 28, 1, 13 / 28, 13 / 28, 1, 1, 13 / 28, 13 / 28, 1, 1, 1),
    outcome = as.integer(c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)),
    competing = 0L
  )
  hazardModels <- suppressMessages(fitCensoringModels(clones, data, iwArgs()))
  hazard <- collectSegments(suppressMessages(computeCensoringWeights(clones, hazardModels, iwNoTruncation())))
  expect_equal(hazard[hazard$strategyId == 1, columns], initiate, tolerance = 1e-6, ignore_attr = TRUE)
  stabilized <- collectSegments(suppressMessages(computeCensoringWeights(clones, hazardModels, iwNoTruncation(stabilized = TRUE))))
  expect_equal(stabilized$weight, rep(1, nrow(stabilized)), tolerance = 1e-6)
  close(clones)
})

test_that("Truncation caps the initiateWithin hazard weights on both scales", {
  data <- iwData()
  clones <- buildM3Clones(iwPersons(), m3Strategies(graceEnd = 2), "eventsFirst")
  on.exit({
    close(clones)
    close(data)
  })
  models <- suppressMessages(fitCensoringModels(clones, data, iwArgs()))
  after <- function(segments) segments$strategyId == 1 & segments$startTime >= 3
  # A fixed cap of 1.8 below 35/19, unstabilized weights; the stabilized weights (1) stay below it:
  for (truncateMsmWeights in c(TRUE, FALSE)) {
    weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(
      stabilized = FALSE, truncation = createWeightTruncation("fixed", upper = 1.8), truncateMsmWeights = truncateMsmWeights
    )))
    segments <- collectSegments(weighted)
    caps <- attr(weighted, "metaData")$weights$caps
    close(weighted)
    expect_equal(segments$weight[after(segments)], rep(1.8, 6))
    # Clones 1, 2 and 4 of the 7 initiate clones are truncated:
    expect_equal(caps$truncatedShare[1], 3 / 7)
    otherScale <- if (truncateMsmWeights) 1 else 1.8 * 19 / 35
    expect_equal(segments$otherScaleWeight[after(segments)], rep(otherScale, 6), tolerance = 1e-6)
  }
  # Stabilized weights with a fixed cap on the unstabilized scale (otherScaleUpper):
  weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(
    stabilized = TRUE, truncation = createWeightTruncation("fixed", upper = 5, otherScaleUpper = 1.5)
  )))
  segments <- collectSegments(weighted)
  close(weighted)
  expect_equal(segments$weight[after(segments)], rep(1, 6), tolerance = 1e-6)
  expect_equal(segments$otherScaleWeight[after(segments)], rep(1.5, 6))
  # A percentile cap over the clones whose weight varies (basis = "varying"), all at 35/19 here:
  weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = FALSE)))
  expect_equal(attr(weighted, "metaData")$weights$caps$cap[1], 35 / 19, tolerance = 1e-6)
  close(weighted)
})

# 1 - S(G | X) per person (rows of `rowIds`) from a pooled logistic glm with daily indicators (or a
# Cox model with Breslow ties) on hand-built rows through G, and 1 - S-bar(G) from the model without
# covariates (contract section 5.3).
referenceStartByG <- function(fixture, rowIds, g, sameDayOrdering, hazardModelType) {
  if (hazardModelType == "cox") {
    rows <- refCoxRows(fixture$persons, g, sameDayOrdering)
    x <- m3CovariateMatrix(fixture, rows$rowId)
    fit <- survival::coxph(survival::Surv(rows$time, rows$initiated) ~ x, ties = "breslow", control = m3CoxControl())
    newX <- m3CovariateMatrix(fixture, rowIds)
    survival <- refCoxSurvival(fit, data.frame(x = I(newX)), g)[, 1]
    marginal <- survival::survfit(survival::Surv(rows$time, rows$initiated) ~ 1, stype = 2, ctype = 1)
    numerator <- summary(marginal, times = g + 1, extend = TRUE)$surv
    return(list(p = 1 - survival, pBar = 1 - numerator))
  }
  rows <- refIntervalRows(fixture$persons, g, sameDayOrdering, 1)
  initiations <- tapply(rows$initiated, rows$k, sum)
  atRisk <- tapply(rows$initiated, rows$k, length)
  intervals <- as.numeric(names(initiations))
  fitIntervals <- intervals[initiations > 0 & initiations < atRisk]
  fitRows <- rows[rows$k %in% fitIntervals, ]
  x <- m3CovariateMatrix(fixture, fitRows$rowId)
  k <- factor(fitRows$k)
  baselineOf <- function(coefs) {
    baseline <- rep(-Inf, g + 1)
    baseline[fitIntervals + 1] <- coefs[1] + c(0, coefs[grepl("^k", names(coefs))])
    baseline[intervals[initiations > 0 & initiations == atRisk] + 1] <- Inf
    return(baseline)
  }
  fit <- glm(fitRows$initiated ~ k + x, family = quasibinomial(), control = glm.control(epsilon = 1e-12))
  baseline <- baselineOf(coef(fit))
  linearPredictor <- as.vector(m3CovariateMatrix(fixture, rowIds) %*% coef(fit)[c("xx1", "xx2", "xx3")])
  survival <- vapply(linearPredictor, function(lp) prod(1 - plogis(baseline + lp)), numeric(1))
  numeratorBaseline <- baselineOf(coef(glm(fitRows$initiated ~ k, family = quasibinomial(), control = glm.control(epsilon = 1e-12))))
  return(list(p = 1 - survival, pBar = 1 - prod(1 - plogis(numeratorBaseline))))
}

test_that("initiateWithin hazard weights match glm and coxph references, stabilized and not", {
  skip_if_not_installed("survival")
  configurations <- list(
    list(riskStart = 0L, sameDayOrdering = "eventsFirst", hazardModelType = "logistic"),
    list(riskStart = 1L, sameDayOrdering = "deviationFirst", hazardModelType = "logistic"),
    list(riskStart = 0L, sameDayOrdering = "deviationFirst", hazardModelType = "cox"),
    list(riskStart = 1L, sameDayOrdering = "eventsFirst", hazardModelType = "cox")
  )
  for (configuration in configurations) {
    info <- paste(unlist(configuration), collapse = ", ")
    sameDayOrdering <- configuration$sameDayOrdering
    fixture <- createM3Fixture(150, riskStart = configuration$riskStart, seed = 41)
    data <- createM3Data(fixture)
    clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "neverInitiate", graceEnd = 7), sameDayOrdering)
    models <- suppressMessages(fitCensoringModels(clones, data, iwArgs(configuration$hazardModelType)))
    binaryModels <- suppressMessages(fitCensoringModels(clones, data, iwArgs(configuration$hazardModelType, initiateWithinModel = "binary")))
    # initiateWithin(7) and neverInitiate have different horizons, so they do not share a fit:
    expect_equal(models[[1]]$horizon, 7L, info = info)
    expect_gt(models[[2]]$horizon, 7L)
    cutDay <- if (sameDayOrdering == "eventsFirst") 7 else 6
    for (stabilized in c(FALSE, TRUE)) {
      segments <- collectSegments(suppressMessages(computeCensoringWeights(clones, models, iwNoTruncation(stabilized))))
      binary <- collectSegments(suppressMessages(computeCensoringWeights(clones, binaryModels, iwNoTruncation(stabilized))))
      expect_equal(segments[, c("cloneRowId", "startTime", "stopTime", "outcome", "competing")], binary[, c("cloneRowId", "startTime", "stopTime", "outcome", "competing")], info = info)
      expect_equal(segments[segments$strategyId == 2, ], binary[binary$strategyId == 2, ], info = info)
      initiate <- segments[segments$strategyId == 1, ]
      stopDay <- initiate$stopTime + configuration$riskStart - 1
      after <- stopDay > cutDay
      expect_true(sum(after) > 20, info = info)
      expect_equal(initiate$weight[!after], rep(1, sum(!after)), info = info)
      reference <- referenceStartByG(fixture, initiate$rowId[after], 7, sameDayOrdering, configuration$hazardModelType)
      expected <- 1 / reference$p
      if (stabilized) {
        expected <- expected * reference$pBar
      }
      expect_equal(initiate$weight[after], expected, tolerance = 1e-6, info = info)
      expect_equal(initiate$stabilizationNumerator[after], rep(reference$pBar, sum(after)), tolerance = 1e-6, info = info)
      # The hazard model gives other weights than the binary model at G:
      expect_gt(max(abs(initiate$weight[after] - binary$weight[binary$strategyId == 1][after])), 0.01)
    }
    close(clones)
    close(data)
  }
})

test_that("The initiateWithin hazard model is shared with notInitiateWithin(G) only", {
  fixture <- createM3Fixture(150, riskStart = 0L, seed = 42)
  data <- createM3Data(fixture)
  on.exit(close(data))
  for (order in list(c("initiateWithin", "notInitiateWithin"), c("notInitiateWithin", "initiateWithin"))) {
    clones <- buildM3Clones(fixture$persons, m3Strategies(order[1], order[2], graceEnd = 7), "eventsFirst")
    messages <- capture_messages(models <- fitCensoringModels(clones, data, iwArgs()))
    expect_equal(sum(grepl("^Fitting", messages)), 1)
    expect_true(any(grepl("Strategy 2 shares the censoring model of an earlier strategy", messages)))
    expect_equal(models[[1]][-(1:2)], models[[2]][-(1:2)])
    expect_equal(c(models[[1]]$type, models[[2]]$type), order)
    expect_equal(c(models[[1]]$kind, models[[2]]$kind), c("hazard", "hazard"))
    # With postGraceModel = "binary", the binary model at G is fit for the notInitiateWithin strategy:
    messages <- capture_messages(binary <- fitCensoringModels(clones, data, iwArgs(postGraceModel = "binary")))
    expect_equal(sum(grepl("^Fitting", messages)), 2)
    notInitiate <- which(order == "notInitiateWithin")
    expect_equal(binary[[notInitiate]]$postGrace$kind, "binary")
    expect_null(binary[[3 - notInitiate]]$postGrace)
    close(clones)
  }
  # A fixed prior variance (the bootstrap) applies to the hazard model of the initiateWithin strategy:
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "neverInitiate", graceEnd = 7), "eventsFirst")
  args <- iwArgs(prior = Cyclops::createPrior("laplace", exclude = c(0), useCrossValidation = TRUE))
  fixed <- suppressMessages(fitCensoringModels(clones, data, args, .fixedPriorVariances = list("1" = 0.37, "2" = 0.21)))
  expect_equal(fixed[[1]]$priorVariance, 0.37)
  expect_equal(fixed[[2]]$priorVariance, 0.21)
  close(clones)
})

test_that("Calibration and positivity report the initiateWithin hazard model", {
  fixture <- createM3Fixture(150, riskStart = 0L, seed = 43)
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = 7), "eventsFirst")
  on.exit({
    close(clones)
    close(data)
  })
  models <- suppressMessages(fitCensoringModels(clones, data, iwArgs(covariateTimeInteraction = "day0")))
  calibration <- computeCensoringModelCalibration(clones, models, numberOfGroups = 4)
  expect_equal(unique(calibration$model), "hazard")
  expect_equal(unique(calibration$covariateTimeInteraction), "day0")
  persons <- dplyr::collect(clones$persons)
  riskSet <- .getBinaryRiskSet(persons, persons$treatmentDay, 7L, "eventsFirst")
  overall <- calibration[calibration$grouping == "overall", ]
  expect_equal(overall$strategyId, 1:2)
  expect_equal(overall$predicted[1], mean(1 - .computeHazardSurvival(models[[1]], riskSet$rowId, 7L)))
  expect_equal(overall$observed[1], mean(riskSet$y))
  # The shared model and the same persons give strategy 1 the rows of strategy 2:
  columns <- c("grouping", "group", "persons", "predicted", "observed")
  expect_equal(calibration[calibration$strategyId == 1, columns], calibration[calibration$strategyId == 2, columns])

  positivity <- computePositivity(clones, models, threshold = 0.3)
  expect_equal(positivity$model, c("hazard", "hazard"))
  p <- 1 - .computeHazardSurvival(models[[1]], models[[1]]$predictions$rowId, 7L)
  # Following initiateWithin is starting by G (p), and notInitiateWithin not starting (1 - p):
  expect_equal(positivity$shareExtreme, c(mean(p < 0.3), mean(1 - p < 0.3)))
  expect_equal(positivity$mean[1], mean(p))
})

test_that("A message says when covariate-by-time interactions are dropped as constant", {
  # Nobody starts after day 6, so with 7-day intervals and G = 7 the second interval (day 7) has no
  # initiations and only the first one (days 0-6) is fit: its "day0" column is constant over the fit.
  fixture <- createM3Fixture(150, riskStart = 0L, seed = 44, maxTreatmentDay = 6L)
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = 7), "eventsFirst")
  on.exit({
    close(clones)
    close(data)
  })
  weekly <- createTimeBasis(intervalLength = 7)
  messages <- capture_messages(models <- fitCensoringModels(clones, data, m3Args("logistic", timeBasis = weekly, covariateTimeInteraction = "day0")))
  dropped <- grep("covariate-by-time interactions", messages, value = TRUE)
  expect_length(dropped, 1)
  expect_match(dropped, "logistic initiation hazard model through day 7 for strategy 2 .* has no covariate-by-time interactions although covariateTimeInteraction = \"day0\"")
  expect_match(dropped, "interval 0 of 2, each 7 days")
  expect_equal(ncol(models[[2]]$interactions$basis), 0)
  expect_false(any(grepl(":", names(models[[2]]$coefficients))))
  # "interval" has no column either with a single interval in the fit:
  messages <- capture_messages(fitCensoringModels(clones, data, m3Args("logistic", timeBasis = weekly, covariateTimeInteraction = "interval")))
  expect_length(grep("covariateTimeInteraction = \"interval\"", messages), 1)
  # Neither "none" nor daily intervals, where day 0 has its own interval, give the message:
  messages <- capture_messages(fitCensoringModels(clones, data, m3Args("logistic", timeBasis = weekly, covariateTimeInteraction = "none")))
  expect_false(any(grepl("covariate-by-time interactions", messages)))
  messages <- capture_messages(models <- fitCensoringModels(clones, data, m3Args("logistic", covariateTimeInteraction = "day0")))
  expect_false(any(grepl("covariate-by-time interactions", messages)))
  expect_equal(colnames(models[[2]]$interactions$basis), "day0")
  # Without covariates there are no interactions to drop:
  messages <- capture_messages(fitCensoringModels(clones, data, m3Args("logistic", useCovariates = FALSE, timeBasis = weekly, covariateTimeInteraction = "day0")))
  expect_false(any(grepl("covariate-by-time interactions", messages)))
})

library(testthat)

# Tier 1 references for computeCensoringWeights() (plan section 6.4): hand-computed W on a small
# fixture; W from glm (binary and pooled logistic) and from survfit (Cox) on larger fixtures;
# segments vs survival::survSplit with t-scale cuts, and survcheck.

noTruncation <- function(stabilized = FALSE) {
  createComputeWeightsArgs(stabilized = stabilized, truncation = createWeightTruncation("none"))
}

# Six persons under initiateWithin(2) vs notInitiateWithin(2), riskStart = 0, nominal risk end day 10.
handPersons <- function() {
  data.frame(
    rowId = 1:6,
    personSeqId = 1:6,
    riskStart = 0L,
    riskEnd = c(10L, 10L, 10L, 6L, 10L, 10L),
    nominalRiskEnd = 10L,
    outcomeDay = c(NA, 5L, 1L, NA, NA, NA),
    deathDay = c(NA, NA, NA, NA, NA, 8L),
    treatmentDay = c(0L, 1L, NA, 2L, NA, 4L),
    comparatorTreatmentDay = NA_integer_
  )
}

handData <- function() {
  createM3Data(list(
    persons = handPersons(),
    covariates = data.frame(rowId = integer(0), covariateId = numeric(0), covariateValue = numeric(0))
  ))
}

handSegments <- function(cloneRowId, startTime, stopTime, weight, outcome, competing) {
  data.frame(
    cloneRowId = as.integer(cloneRowId),
    startTime = as.numeric(startTime),
    stopTime = as.numeric(stopTime),
    weight = weight,
    outcome = as.integer(outcome),
    competing = as.integer(competing)
  )
}

test_that("Weights equal hand-computed values on a six-person fixture", {
  data <- handData()
  on.exit(close(data))
  columns <- c("cloneRowId", "startTime", "stopTime", "weight", "outcome", "competing")

  # eventsFirst. Binary risk set at day 2: persons 1, 2, 4, 5, 6 (3 has the outcome on day 1), of whom
  # 1, 2, 4 started: pi = 3/5. Initiation risk set (time = last day at risk + 1): 1 (t = 1, starts),
  # 2 (t = 2, starts), 3 (t = 1, outcome on day 1), 4 (t = 3, starts), 5 and 6 (t = 3).
  # Breslow: Lambda(1) = 1/6, Lambda(2) = 1/6 + 1/4, Lambda(3) = 5/12 + 1/3.
  # Pooled logistic, daily indicators: h = 1/6, 1/4, 1/3.
  clones <- buildM3Clones(handPersons(), m3Strategies(graceEnd = 2), "eventsFirst")
  initiate <- handSegments(
    cloneRowId = c(1, 1, 2, 2, 3, 4, 4, 5, 6),
    startTime = c(0, 3, 0, 3, 0, 0, 3, 0, 0),
    stopTime = c(3, 11, 3, 6, 2, 3, 7, 3, 3),
    weight = c(1, 5 / 3, 1, 5 / 3, 1, 1, 5 / 3, 1, 1),
    outcome = c(0, 0, 0, 1, 1, 0, 0, 0, 0),
    competing = 0
  )
  deferCox <- handSegments(
    cloneRowId = c(7, 8, 9, 10, 10, 11, 11, 12, 12),
    startTime = c(0, 0, 0, 0, 2, 0, 2, 0, 2),
    stopTime = c(1, 2, 2, 2, 3, 2, 11, 2, 9),
    weight = exp(c(0, 1 / 6, 1 / 6, 1 / 6, 5 / 12, 1 / 6, 3 / 4, 1 / 6, 3 / 4)),
    outcome = c(0, 0, 1, 0, 0, 0, 0, 0, 0),
    competing = c(0, 0, 0, 0, 0, 0, 0, 0, 1)
  )
  models <- suppressMessages(fitCensoringModels(clones, data, m3Args("cox", useCovariates = FALSE)))
  weighted <- suppressMessages(computeCensoringWeights(clones, models, noTruncation()))
  expect_equal(collectSegments(weighted)[, columns], rbind(initiate, deferCox), tolerance = 1e-12)

  deferLogistic <- deferCox
  deferLogistic$weight <- c(1, 6 / 5, 6 / 5, 6 / 5, 8 / 5, 6 / 5, 12 / 5, 6 / 5, 12 / 5)
  models <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic", useCovariates = FALSE)))
  weighted <- suppressMessages(computeCensoringWeights(clones, models, noTruncation()))
  # Cyclops converges to about 1e-7 relative on the saturated numerator model:
  expect_equal(collectSegments(weighted)[, columns], rbind(initiate, deferLogistic), tolerance = 1e-6)
  close(clones)

  # deviationFirst. Same binary risk set. Initiation risk set: 3 is at risk on day 1 (t = 2).
  # Breslow: Lambda(1) = 1/6, Lambda(2) = 1/6 + 1/5, Lambda(3) = 11/30 + 1/3; W(d) = 1/S(d).
  # A clone that deviates on day d does so before any outcome on day d, so day d is not in its
  # segments: non-starters 5 and 6 end on day 1 (t = 2) in the initiate arm, and in the defer arm
  # clone 7 (initiates on day 0) has no segments, 8 ends on day 0 and 10 on day 1.
  clones <- buildM3Clones(handPersons(), m3Strategies(graceEnd = 2), "deviationFirst")
  initiate <- handSegments(
    cloneRowId = c(1, 1, 2, 2, 3, 4, 4, 5, 6),
    startTime = c(0, 2, 0, 2, 0, 0, 2, 0, 0),
    stopTime = c(2, 11, 2, 6, 2, 2, 7, 2, 2),
    weight = c(1, 5 / 3, 1, 5 / 3, 1, 1, 5 / 3, 1, 1),
    outcome = c(0, 0, 0, 1, 1, 0, 0, 0, 0),
    competing = 0
  )
  defer <- handSegments(
    cloneRowId = c(8, 9, 10, 11, 11, 12, 12),
    startTime = c(0, 0, 0, 0, 2, 0, 2),
    stopTime = c(1, 2, 2, 2, 11, 2, 9),
    weight = exp(c(1 / 6, 11 / 30, 11 / 30, 11 / 30, 7 / 10, 11 / 30, 7 / 10)),
    outcome = c(0, 1, 0, 0, 0, 0, 0),
    competing = c(0, 0, 0, 0, 0, 0, 1)
  )
  models <- suppressMessages(fitCensoringModels(clones, data, m3Args("cox", useCovariates = FALSE)))
  weighted <- suppressMessages(computeCensoringWeights(clones, models, noTruncation()))
  expect_equal(collectSegments(weighted)[, columns], rbind(initiate, defer), tolerance = 1e-12)
  expect_equal(sort(dplyr::pull(weighted$personWeights, "cloneRowId")), c(1:6, 8:12))

  # Stabilized weights of models without covariates are exactly 1:
  weighted <- suppressMessages(computeCensoringWeights(clones, models, noTruncation(stabilized = TRUE)))
  expect_equal(collectSegments(weighted)$weight, rep(1, nrow(initiate) + nrow(defer)), tolerance = 1e-12)
  close(clones)
})

# Reference W for every segment, from glm (binary) and coxph/survfit (Cox), per contract section 5.4.
# Strategy 1 is initiateWithin(graceEnd); strategy 2 a hazard strategy with the given horizon.
referenceWeights <- function(segments, clones, fixture, sameDayOrdering, graceEnd, stabilized, horizon = graceEnd) {
  persons <- fixture$persons
  offset <- if (sameDayOrdering == "eventsFirst") 1 else 0
  riskStart <- clones$riskStart[match(segments$cloneRowId, clones$cloneRowId)]
  stopDay <- segments$stopTime + riskStart - 1
  weight <- rep(NA_real_, nrow(segments))

  riskSet <- refBinaryRiskSet(persons, graceEnd, sameDayOrdering)
  x <- m3CovariateMatrix(fixture, riskSet$rowId)
  binary <- glm(riskSet$y ~ x, family = binomial(), control = glm.control(epsilon = 1e-12))
  ps <- fitted(binary)[match(segments$rowId, riskSet$rowId)]
  s1 <- segments$strategyId == 1
  after <- stopDay > graceEnd - (1 - offset)
  weight[s1] <- ifelse(after[s1], (if (stabilized) mean(riskSet$y) else 1) / ps[s1], 1)

  rows <- refCoxRows(persons, horizon, sameDayOrdering)
  xc <- m3CovariateMatrix(fixture, rows$rowId)
  cox <- survival::coxph(survival::Surv(rows$time, rows$initiated) ~ xc, ties = "breslow", control = m3CoxControl())
  days <- -1:horizon
  survival <- refCoxSurvival(cox, data.frame(xc = I(xc)), days)
  cox0 <- survival::coxph(survival::Surv(rows$time, rows$initiated) ~ 1, ties = "breslow")
  survival0 <- c(1, summary(survival::survfit(cox0, stype = 2, ctype = 1), times = days[-1] + 1, extend = TRUE)$surv)
  s2 <- which(segments$strategyId == 2)
  u <- pmin(stopDay[s2] - offset, horizon)
  person <- match(segments$rowId[s2], rows$rowId)
  weight[s2] <- 1 / survival[cbind(person, u + 2)]
  if (stabilized) {
    weight[s2] <- weight[s2] * survival0[u + 2]
  }
  return(weight)
}

test_that("Binary and Cox weights match glm and survfit references, with segments as survSplit", {
  skip_if_not_installed("survival")
  configurations <- list(
    list(riskStart = 0L, sameDayOrdering = "eventsFirst", type = "notInitiateWithin"),
    list(riskStart = 0L, sameDayOrdering = "deviationFirst", type = "notInitiateWithin"),
    list(riskStart = 1L, sameDayOrdering = "eventsFirst", type = "notInitiateWithin"),
    list(riskStart = 1L, sameDayOrdering = "deviationFirst", type = "notInitiateWithin"),
    list(riskStart = 0L, sameDayOrdering = "deviationFirst", type = "neverInitiate"),
    list(riskStart = 1L, sameDayOrdering = "eventsFirst", type = "neverInitiate")
  )
  for (configuration in configurations) {
    sameDayOrdering <- configuration$sameDayOrdering
    fixture <- createM3Fixture(150, riskStart = configuration$riskStart, seed = 21)
    data <- createM3Data(fixture)
    strategies <- m3Strategies("initiateWithin", configuration$type, graceEnd = 7)
    horizon <- if (configuration$type == "neverInitiate") max(fixture$persons$riskEnd) else 7L
    clones <- buildM3Clones(fixture$persons, strategies, sameDayOrdering)
    cloneTable <- collectClones(clones)
    followUp <- refOutcomeFollowUp(cloneTable, sameDayOrdering)
    models <- suppressMessages(fitCensoringModels(clones, data, m3Args("cox")))
    for (stabilized in c(FALSE, TRUE)) {
      weighted <- suppressMessages(computeCensoringWeights(clones, models, noTruncation(stabilized)))
      segments <- collectSegments(weighted)

      expected <- refSegments(followUp, refCuts(followUp, strategies, sameDayOrdering, list(7, horizon)))
      expect_equal(segments[, names(expected)], expected)
      expect_equal(segments$segmentId, seq_len(nrow(segments)))
      expect_equal(segments$rowId, cloneTable$rowId[segments$cloneRowId])
      expect_equal(segments$personSeqId, cloneTable$personSeqId[segments$cloneRowId])
      expect_equal(segments$strategyId, cloneTable$strategyId[segments$cloneRowId])

      expectedWeight <- referenceWeights(segments, cloneTable, fixture, sameDayOrdering, 7, stabilized, horizon)
      expect_equal(segments$weight, expectedWeight, tolerance = 1e-6)
    }
    check <- survival::survcheck(survival::Surv(startTime, stopTime, outcome) ~ 1, id = cloneRowId, data = segments)
    expect_true(all(check$flag == 0))
    expect_equal(sum(segments$outcome), sum(cloneTable$endReason == "outcome"))
    expect_equal(sum(segments$competing), sum(cloneTable$endReason == "death"))
    close(clones)
    close(data)
  }
})

test_that("With deviationFirst, weighted risk sets equal the event-free persons when initiation is unconfounded", {
  # Saturated daily pooled logistic models without covariates: on day 0 the defer arm's clones at risk
  # (no initiation on day 0) weigh 1 / (1 - h0), so they add up to everyone at risk on day 0. In the
  # initiate arm, the starters at risk on day G weigh 1 / pi, so they add up to the binary risk set.
  # Deviation-ended clones must not be in the risk set on their deviation day.
  for (riskStart in 0:1) {
    fixture <- createM3Fixture(200, riskStart = riskStart, seed = 26)
    persons <- fixture$persons
    # An outcome on the first day of follow-up, so that W on that day is a segment weight:
    persons$outcomeDay[1] <- riskStart
    persons$deathDay[1] <- NA
    persons$treatmentDay[1] <- NA
    data <- createM3Data(list(persons = persons, covariates = fixture$covariates))
    clones <- buildM3Clones(persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = 5), "deviationFirst")
    models <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic", useCovariates = FALSE)))
    segments <- collectSegments(suppressMessages(computeCensoringWeights(clones, models, noTruncation())))
    weightedAtRisk <- function(strategyId, day) {
      t <- day - riskStart + 1
      sum(segments$weight[segments$strategyId == strategyId & segments$startTime < t & t <= segments$stopTime])
    }
    riskSet <- refBinaryRiskSet(persons, 5, "deviationFirst")
    expect_equal(weightedAtRisk(1, 5), nrow(riskSet), tolerance = 1e-10)
    if (riskStart == 0) {
      expect_equal(weightedAtRisk(2, 0), nrow(persons), tolerance = 1e-6)
    }
    close(clones)
    close(data)
  }
})

# Reference W of the strategy 2 segments from pooled logistic glm fits on hand-built K_g rows (indicator
# basis): hazard 0 in intervals without initiations and 1 in intervals where everyone at risk
# initiates; W(d) = 1 / prod of (1 - h) over the intervals ending on or before d - 1 (d with
# deviationFirst) and on or before the horizon. Stabilized: times the same product from glm(y ~ k).
referenceLogisticWeights <- function(segments, fixture, sameDayOrdering, horizon, intervalLength, stabilized) {
  offset <- if (sameDayOrdering == "eventsFirst") 1 else 0
  rows <- refIntervalRows(fixture$persons, horizon, sameDayOrdering, intervalLength)
  initiations <- tapply(rows$initiated, rows$k, sum)
  atRisk <- tapply(rows$initiated, rows$k, length)
  intervals <- as.numeric(names(initiations))
  fitIntervals <- intervals[initiations > 0 & initiations < atRisk]
  fullIntervals <- intervals[initiations > 0 & initiations == atRisk]
  fitRows <- rows[rows$k %in% fitIntervals, ]
  x <- m3CovariateMatrix(fixture, fitRows$rowId)
  k <- factor(fitRows$k)
  nIntervals <- horizon %/% intervalLength + 1
  intervalEnd <- pmin((0:(nIntervals - 1) + 1) * intervalLength - 1, horizon)
  getBaseline <- function(coefs) {
    baseline <- rep(-Inf, nIntervals)
    baseline[fitIntervals + 1] <- coefs[1] + c(0, coefs[grepl("^k", names(coefs))])
    baseline[fullIntervals + 1] <- Inf
    return(baseline)
  }
  reference <- glm(fitRows$initiated ~ k + x, family = quasibinomial(), control = glm.control(epsilon = 1e-12))
  baseline <- getBaseline(coef(reference))
  numeratorBaseline <- getBaseline(coef(glm(fitRows$initiated ~ k, family = quasibinomial(), control = glm.control(epsilon = 1e-12))))
  s2 <- which(segments$strategyId == 2)
  personX <- m3CovariateMatrix(fixture, segments$rowId[s2])
  linearPredictor <- as.vector(personX %*% coef(reference)[c("xx1", "xx2", "xx3")])
  riskStart <- fixture$persons$riskStart[match(segments$rowId[s2], fixture$persons$rowId)]
  u <- pmin(segments$stopTime[s2] + riskStart - 1 - offset, horizon)
  return(vapply(seq_along(s2), function(i) {
    complete <- intervalEnd <= u[i]
    weight <- 1 / prod(1 - plogis(baseline[complete] + linearPredictor[i]))
    if (stabilized) {
      weight <- weight * prod(1 - plogis(numeratorBaseline[complete]))
    }
    weight
  }, numeric(1)))
}

test_that("Pooled logistic weights match products of glm hazards", {
  skip_if_not_installed("survival")
  intervalLength <- 3
  # notInitiateWithin(7) with 3-day intervals has a partial last interval (days 6-7), and its weight
  # after G uses all of it.
  configurations <- list(
    list(riskStart = 0L, sameDayOrdering = "eventsFirst", type = "neverInitiate"),
    list(riskStart = 1L, sameDayOrdering = "deviationFirst", type = "neverInitiate"),
    list(riskStart = 0L, sameDayOrdering = "deviationFirst", type = "notInitiateWithin"),
    list(riskStart = 1L, sameDayOrdering = "eventsFirst", type = "notInitiateWithin")
  )
  for (configuration in configurations) {
    sameDayOrdering <- configuration$sameDayOrdering
    fixture <- createM3Fixture(150, riskStart = configuration$riskStart, seed = 22)
    data <- createM3Data(fixture)
    strategies <- m3Strategies("initiateWithin", configuration$type, graceEnd = 7)
    horizon <- if (configuration$type == "neverInitiate") max(fixture$persons$riskEnd) else 7L
    clones <- buildM3Clones(fixture$persons, strategies, sameDayOrdering)
    followUp <- refOutcomeFollowUp(collectClones(clones), sameDayOrdering)
    models <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic", timeBasis = createTimeBasis(intervalLength = intervalLength))))
    for (stabilized in c(FALSE, TRUE)) {
      weighted <- suppressMessages(computeCensoringWeights(clones, models, noTruncation(stabilized)))
      segments <- collectSegments(weighted)
      expected <- refSegments(followUp, refCuts(followUp, strategies, sameDayOrdering, list(7, horizon)))
      expect_equal(segments[, names(expected)], expected)
      expectedWeight <- referenceLogisticWeights(segments, fixture, sameDayOrdering, horizon, intervalLength, stabilized)
      expect_equal(segments$weight[segments$strategyId == 2], expectedWeight, tolerance = 1e-6)
    }
    if (configuration$type == "notInitiateWithin") {
      # Some clones are followed after G, with the weight of the full partial interval:
      expect_true(any(segments$strategyId == 2 & segments$stopTime + configuration$riskStart - 1 > 7))
    }
    close(clones)
    close(data)
  }
})

test_that("Truncation caps weights per strategy at a percentile of the clones' maximum weight", {
  fixture <- createM3Fixture(150, riskStart = 0L, seed = 23)
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "neverInitiate", graceEnd = 7), "eventsFirst")
  models <- suppressMessages(fitCensoringModels(clones, data, m3Args("cox")))

  untruncated <- suppressMessages(computeCensoringWeights(clones, models, noTruncation()))
  raw <- collectSegments(untruncated)
  personWeights <- as.data.frame(dplyr::collect(untruncated$personWeights))
  personWeights <- personWeights[order(personWeights$cloneRowId), ]
  expect_equal(personWeights$cloneRowId, seq_len(nrow(collectClones(clones))))
  expect_equal(personWeights$maxWeight, as.vector(tapply(raw$weight, raw$cloneRowId, max)))
  lastSegment <- raw[!duplicated(raw$cloneRowId, fromLast = TRUE), ]
  expect_equal(personWeights$weightAtEnd, lastSegment$weight)
  metaData <- attr(untruncated, "metaData")$weights
  expect_false(metaData$stabilized)
  expect_true(all(is.na(metaData$caps$cap)))
  expect_equal(metaData$untruncatedSummary$maxMaxWeight, as.vector(tapply(personWeights$maxWeight, personWeights$strategyId, max)))
  expect_true(checkmate::testTRUE(.checkCcwClones(untruncated, weighted = TRUE)))

  percentile <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = FALSE, truncation = createWeightTruncation("percentile", 0.9))))
  caps <- tapply(personWeights$maxWeight, personWeights$strategyId, quantile, probs = 0.9, type = 7)
  expect_equal(attr(percentile, "metaData")$weights$caps$cap, as.vector(caps))
  truncated <- collectSegments(percentile)
  expect_equal(truncated$weight, pmin(raw$weight, caps[raw$strategyId]))
  # personWeights keep the untruncated values:
  expect_equal(dplyr::collect(percentile$personWeights) |> dplyr::arrange(.data$cloneRowId) |> dplyr::pull("maxWeight"), personWeights$maxWeight)

  fixed <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = FALSE, truncation = createWeightTruncation("fixed", 1.5))))
  expect_equal(collectSegments(fixed)$weight, pmin(raw$weight, 1.5))
  expect_equal(attr(fixed, "metaData")$weights$caps$cap, c(1.5, 1.5))
  close(clones)
  close(data)
})

test_that("Weights are 1 when stabilized models have no covariates, and inputs are checked", {
  fixture <- createM3Fixture(100, riskStart = 1L, seed = 24)
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = 5), "deviationFirst")
  for (hazardModelType in c("cox", "logistic")) {
    models <- suppressMessages(fitCensoringModels(clones, data, m3Args(hazardModelType, useCovariates = FALSE, timeBasis = createTimeBasis(intervalLength = 2))))
    weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = TRUE)))
    expect_equal(collectSegments(weighted)$weight, rep(1, nrow(collectSegments(weighted))), tolerance = 1e-12)
  }
  # No temporary tables are left behind:
  expect_false(any(startsWith(names(weighted), "ccw")))

  otherClones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = 5), "eventsFirst")
  expect_error(computeCensoringWeights(otherClones, models), "sameDayOrdering")
  otherClones2 <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = 6), "deviationFirst")
  expect_error(computeCensoringWeights(otherClones2, models), "does not match")
  expect_error(computeCensoringWeights(clones, list()), "CcwCensoringModels")
  close(otherClones)
  close(otherClones2)
  close(clones)
  close(data)
})

test_that("Censoring models fit on other clones are rejected", {
  fixture <- createM3Fixture(100, riskStart = 0L, seed = 28)
  personsA <- fixture$persons
  # With outcome A, persons 1-5 have the outcome on day 0, so they are outside the eventsFirst
  # initiation risk set. Outcome B (same outcome ID here) leaves them under follow-up.
  personsA$outcomeDay[1:5] <- 0L
  personsA$deathDay[1:5] <- NA
  personsA$treatmentDay[1:5] <- NA
  personsB <- personsA
  personsB$outcomeDay[1:5] <- NA
  data <- createM3Data(list(persons = personsA, covariates = fixture$covariates))
  strategies <- m3Strategies("initiateWithin", "neverInitiate", graceEnd = 7)
  clonesA <- buildM3Clones(personsA, strategies, "eventsFirst")
  clonesB <- buildM3Clones(personsB, strategies, "eventsFirst")
  for (hazardModelType in c("cox", "logistic")) {
    models <- suppressMessages(fitCensoringModels(clonesA, data, m3Args(hazardModelType)))
    expect_false(any(1:5 %in% models[[2]]$predictions$rowId))
    expect_error(suppressMessages(computeCensoringWeights(clonesB, models)), "outside its risk set")
  }
  attr(clonesB, "metaData")$outcomeId <- 4
  expect_error(computeCensoringWeights(clonesB, models), "different outcome")
  close(clonesA)
  close(clonesB)
  close(data)
})

test_that("Edge cases: nobody starts by G, everyone starts on day 0, an interval where everyone initiates", {
  fixture <- createM3Fixture(150, riskStart = 0L, seed = 27)
  # Fits the models, and returns them with the segments (unstabilized and stabilized weights).
  computeSegments <- function(persons, strategies, sameDayOrdering, hazardModelType) {
    data <- createM3Data(list(persons = persons, covariates = fixture$covariates))
    on.exit(close(data))
    clones <- buildM3Clones(persons, strategies, sameDayOrdering)
    on.exit(close(clones), add = TRUE)
    models <- suppressMessages(fitCensoringModels(clones, data, m3Args(hazardModelType)))
    segments <- lapply(c(FALSE, TRUE), function(stabilized) {
      collectSegments(suppressMessages(computeCensoringWeights(clones, models, noTruncation(stabilized))))
    })
    return(list(models = models, segments = segments[[1]], stabilizedSegments = segments[[2]], followUp = refOutcomeFollowUp(collectClones(clones), sameDayOrdering)))
  }
  strategies <- m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = 7)

  # Nobody starts by G (some start later): no clone is followed after G in the initiate arm, and
  # nobody initiates in the defer arm's risk set, so every weight is 1.
  persons <- fixture$persons
  persons$treatmentDay[!is.na(persons$treatmentDay) & persons$treatmentDay <= 7] <- NA
  for (configuration in list(c("eventsFirst", "cox"), c("deviationFirst", "logistic"))) {
    result <- computeSegments(persons, strategies, configuration[1], configuration[2])
    expect_equal(result$models[[1]]$numerator$propensityScore, 0)
    expect_equal(result$models[[2]]$diagnostics$events, 0)
    expect_equal(result$segments$weight, rep(1, nrow(result$segments)))
    lastDay <- if (configuration[1] == "eventsFirst") 7 else 6
    expect_true(all(result$segments$stopTime[result$segments$strategyId == 1] - 1 <= lastDay))
    expect_setequal(unique(result$segments$cloneRowId), result$followUp$cloneRowId)
  }

  # Everyone starts on day 0: pi = 1, and the hazard of day 0 is 1. Defer clones end on day 0
  # (eventsFirst, W = 1) or have no segments (deviationFirst).
  persons <- fixture$persons
  persons$treatmentDay <- 0L
  for (configuration in list(c("eventsFirst", "logistic"), c("deviationFirst", "logistic"), c("deviationFirst", "cox"))) {
    result <- computeSegments(persons, strategies, configuration[1], configuration[2])
    expect_equal(result$models[[1]]$numerator$propensityScore, 1)
    if (configuration[2] == "logistic") {
      expect_equal(result$models[[2]]$timeBasis$fullIntervals, 0)
      expect_equal(result$models[[2]]$intervals$baseline, c(Inf, rep(-Inf, 7)))
    }
    expect_equal(result$segments$weight, rep(1, nrow(result$segments)))
    defer <- result$segments[result$segments$strategyId == 2, ]
    if (configuration[1] == "eventsFirst") {
      expect_equal(defer$stopTime, rep(1, nrow(defer)))
    } else {
      expect_equal(nrow(defer), 0)
    }
  }

  # An interval where everyone at risk initiates (person 1, alone at risk on day 60, initiates then):
  persons <- fixture$persons
  persons$riskEnd[1] <- 60L
  persons$nominalRiskEnd <- 60L
  persons$outcomeDay[1] <- NA
  persons$deathDay[1] <- NA
  persons$treatmentDay[1] <- 60L
  neverStrategies <- m3Strategies("initiateWithin", "neverInitiate", graceEnd = 7)
  for (sameDayOrdering in c("eventsFirst", "deviationFirst")) {
    result <- computeSegments(persons, neverStrategies, sameDayOrdering, "logistic")
    expect_equal(result$models[[2]]$timeBasis$fullIntervals, 60)
    for (stabilized in c(FALSE, TRUE)) {
      segments <- if (stabilized) result$stabilizedSegments else result$segments
      expectedWeight <- referenceLogisticWeights(segments, list(persons = persons, x = fixture$x), sameDayOrdering, 60, 1, stabilized)
      expect_equal(segments$weight[segments$strategyId == 2], expectedWeight, tolerance = 1e-6)
    }
  }
})

test_that("computeWeightDistribution summarizes the weights of clones under follow-up", {
  fixture <- createM3Fixture(120, riskStart = 1L, seed = 25)
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies(graceEnd = 7), "eventsFirst")
  models <- suppressMessages(fitCensoringModels(clones, data, m3Args("cox")))
  weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs()))
  segments <- collectSegments(weighted)
  segments$riskStart <- 1L

  distribution <- computeWeightDistribution(weighted, days = c(1, 5, 8, 20))
  expect_equal(colnames(distribution), c(
    "strategyId", "day", "clones", "sumWeights", "meanWeight", "minWeight", "p25Weight",
    "medianWeight", "p75Weight", "maxWeight", "effectiveSampleSize"
  ))
  for (i in seq_len(nrow(distribution))) {
    t <- distribution$day[i] - 1 + 1
    w <- segments$weight[segments$strategyId == distribution$strategyId[i] & segments$startTime < t & t <= segments$stopTime]
    expect_equal(distribution$clones[i], length(w))
    expect_equal(distribution$sumWeights[i], sum(w))
    expect_equal(distribution$medianWeight[i], median(w))
    expect_equal(distribution$p25Weight[i], unname(quantile(w, 0.25)))
    expect_equal(distribution$maxWeight[i], max(w))
    expect_equal(distribution$effectiveSampleSize[i], sum(w)^2 / sum(w^2))
  }
  expect_equal(nrow(distribution), 8)

  defaultDays <- computeWeightDistribution(weighted)
  expect_equal(sort(unique(defaultDays$day)), sort(unique(as.integer(segments$stopTime))))
  expect_false(any(startsWith(names(weighted), "ccw")))

  skip_if_not_installed("ggplot2")
  fileName <- tempfile(fileext = ".png")
  plot <- plotWeightDistribution(defaultDays, strategyLabels = c("Initiate", "Defer"), fileName = fileName)
  expect_true(file.exists(fileName))
  unlink(fileName)
  expect_s3_class(plot, "ggplot")
  close(clones)
  close(data)
})

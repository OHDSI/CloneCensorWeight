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
  clones <- buildM3Clones(handPersons(), m3Strategies(graceEnd = 2), "deviationFirst")
  initiate <- handSegments(
    cloneRowId = c(1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 6),
    startTime = c(0, 2, 0, 2, 0, 0, 2, 0, 2, 0, 2),
    stopTime = c(2, 11, 2, 6, 2, 2, 7, 2, 3, 2, 3),
    weight = c(1, 5 / 3, 1, 5 / 3, 1, 1, 5 / 3, 1, 5 / 3, 1, 5 / 3),
    outcome = c(0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0),
    competing = 0
  )
  defer <- handSegments(
    cloneRowId = c(7, 8, 9, 10, 10, 11, 11, 12, 12),
    startTime = c(0, 0, 0, 0, 2, 0, 2, 0, 2),
    stopTime = c(1, 2, 2, 2, 3, 2, 11, 2, 9),
    weight = exp(c(1 / 6, 11 / 30, 11 / 30, 11 / 30, 7 / 10, 11 / 30, 7 / 10, 11 / 30, 7 / 10)),
    outcome = c(0, 0, 1, 0, 0, 0, 0, 0, 0),
    competing = c(0, 0, 0, 0, 0, 0, 0, 0, 1)
  )
  models <- suppressMessages(fitCensoringModels(clones, data, m3Args("cox", useCovariates = FALSE)))
  weighted <- suppressMessages(computeCensoringWeights(clones, models, noTruncation()))
  expect_equal(collectSegments(weighted)[, columns], rbind(initiate, defer), tolerance = 1e-12)

  # Stabilized weights of models without covariates are exactly 1:
  weighted <- suppressMessages(computeCensoringWeights(clones, models, noTruncation(stabilized = TRUE)))
  expect_equal(collectSegments(weighted)$weight, rep(1, nrow(initiate) + nrow(defer)), tolerance = 1e-12)
  close(clones)
})

# Reference W for every segment, from glm (binary) and coxph/survfit (Cox), per contract section 5.4.
referenceWeights <- function(segments, clones, fixture, sameDayOrdering, graceEnd, stabilized) {
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

  rows <- refCoxRows(persons, graceEnd, sameDayOrdering)
  xc <- m3CovariateMatrix(fixture, rows$rowId)
  cox <- survival::coxph(survival::Surv(rows$time, rows$initiated) ~ xc, ties = "breslow", control = m3CoxControl())
  days <- -1:graceEnd
  survival <- refCoxSurvival(cox, data.frame(xc = I(xc)), days)
  cox0 <- survival::coxph(survival::Surv(rows$time, rows$initiated) ~ 1, ties = "breslow")
  survival0 <- c(1, summary(survival::survfit(cox0, stype = 2, ctype = 1), times = days[-1] + 1, extend = TRUE)$surv)
  s2 <- which(segments$strategyId == 2)
  u <- pmin(stopDay[s2] - offset, graceEnd)
  person <- match(segments$rowId[s2], rows$rowId)
  weight[s2] <- 1 / survival[cbind(person, u + 2)]
  if (stabilized) {
    weight[s2] <- weight[s2] * survival0[u + 2]
  }
  return(weight)
}

test_that("Binary and Cox weights match glm and survfit references, with segments as survSplit", {
  skip_if_not_installed("survival")
  for (riskStart in 0:1) {
    for (sameDayOrdering in c("eventsFirst", "deviationFirst")) {
      fixture <- createM3Fixture(150, riskStart = riskStart, seed = 21)
      data <- createM3Data(fixture)
      strategies <- m3Strategies(graceEnd = 7)
      clones <- buildM3Clones(fixture$persons, strategies, sameDayOrdering)
      cloneTable <- collectClones(clones)
      models <- suppressMessages(fitCensoringModels(clones, data, m3Args("cox")))
      for (stabilized in c(FALSE, TRUE)) {
        weighted <- suppressMessages(computeCensoringWeights(clones, models, noTruncation(stabilized)))
        segments <- collectSegments(weighted)

        expected <- refSegments(cloneTable, refCuts(cloneTable, strategies, sameDayOrdering, list(7, 7)))
        expect_equal(segments[, names(expected)], expected)
        expect_equal(segments$segmentId, seq_len(nrow(segments)))
        expect_equal(segments$rowId, cloneTable$rowId[segments$cloneRowId])
        expect_equal(segments$personSeqId, cloneTable$personSeqId[segments$cloneRowId])
        expect_equal(segments$strategyId, cloneTable$strategyId[segments$cloneRowId])

        expectedWeight <- referenceWeights(segments, cloneTable, fixture, sameDayOrdering, 7, stabilized)
        expect_equal(segments$weight, expectedWeight, tolerance = 1e-6)
      }
      check <- survival::survcheck(survival::Surv(startTime, stopTime, outcome) ~ 1, id = cloneRowId, data = segments)
      expect_true(all(check$flag == 0))
      expect_equal(sum(segments$outcome), sum(cloneTable$endReason == "outcome"))
      expect_equal(sum(segments$competing), sum(cloneTable$endReason == "death"))
      close(clones)
      close(data)
    }
  }
})

test_that("Pooled logistic weights for neverInitiate match products of glm hazards", {
  skip_if_not_installed("survival")
  intervalLength <- 3
  for (riskStart in 0:1) {
    sameDayOrdering <- if (riskStart == 0) "eventsFirst" else "deviationFirst"
    offset <- if (sameDayOrdering == "eventsFirst") 1 else 0
    fixture <- createM3Fixture(150, riskStart = riskStart, seed = 22)
    data <- createM3Data(fixture)
    strategies <- m3Strategies("initiateWithin", "neverInitiate", graceEnd = 5)
    clones <- buildM3Clones(fixture$persons, strategies, sameDayOrdering)
    cloneTable <- collectClones(clones)
    models <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic", timeBasis = createTimeBasis(intervalLength = intervalLength))))
    horizon <- max(fixture$persons$riskEnd)
    weighted <- suppressMessages(computeCensoringWeights(clones, models, noTruncation()))
    segments <- collectSegments(weighted)
    expected <- refSegments(cloneTable, refCuts(cloneTable, strategies, sameDayOrdering, list(5, horizon)))
    expect_equal(segments[, names(expected)], expected)

    # Hazards of every person and interval from glm on the hand-built rows (0 in intervals without
    # initiations), then W(d) = 1 / prod of (1 - h) over the intervals ending on or before d - 1.
    rows <- refIntervalRows(fixture$persons, horizon, sameDayOrdering, intervalLength)
    eventIntervals <- sort(unique(rows$k[rows$initiated == 1]))
    fitRows <- rows[rows$k %in% eventIntervals, ]
    x <- m3CovariateMatrix(fixture, fitRows$rowId)
    k <- factor(fitRows$k)
    reference <- glm(fitRows$initiated ~ k + x, family = quasibinomial(), control = glm.control(epsilon = 1e-12))
    nIntervals <- horizon %/% intervalLength + 1
    intervalEnd <- pmin((0:(nIntervals - 1) + 1) * intervalLength - 1, horizon)
    coefs <- coef(reference)
    gamma <- c(0, coefs[grepl("^k", names(coefs))])
    baseline <- rep(-Inf, nIntervals)
    baseline[eventIntervals + 1] <- coefs[1] + gamma
    s2 <- which(segments$strategyId == 2)
    personX <- m3CovariateMatrix(fixture, segments$rowId[s2])
    linearPredictor <- as.vector(personX %*% coefs[c("xx1", "xx2", "xx3")])
    stopDay <- segments$stopTime[s2] + riskStart - 1
    u <- pmin(stopDay - offset, horizon)
    expectedWeight <- vapply(seq_along(s2), function(i) {
      complete <- intervalEnd <= u[i]
      1 / prod(1 - plogis(baseline[complete] + linearPredictor[i]))
    }, numeric(1))
    expect_equal(segments$weight[s2], expectedWeight, tolerance = 1e-6)
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

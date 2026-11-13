library(testthat)

# Integration tests of the Phase1Contracts section 9 amendments (1-8 and 11). Amendments 1-3 (the clone
# rules) are tested with hand-written tables in test-m2-clones.R; here only their consequence downstream.
# The bootstrap driver (amendment 9) is tested in test-integration-bootstrap.R.

# A small simulated analysis: the study population, clones and censoring models (with covariates and
# no prior) of initiateWithin(graceEnd) vs `comparatorType`. The caller closes `data` and `clones`.
integrationAnalysis <- function(comparatorType = "neverInitiate",
                                graceEnd = 5,
                                riskStart = 0,
                                sameDayOrdering = "eventsFirst",
                                n = 300,
                                seed = 1) {
  data <- simulateCloneCensorWeightData(
    n = n, maxDays = 20, graceEnd = max(graceEnd, 1), initiationHazard = 0.05, outcomeHazard = 0.01,
    deathHazard = 0.005, seed = seed
  )
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
      riskWindowStart = riskStart,
      riskWindowEnd = 20,
      endAnchor = "cohort start"
    )
  ))
  strategies <- list(
    createTreatmentStrategy(1, "Initiate", "initiateWithin", graceEnd = graceEnd),
    createTreatmentStrategy(2, "Comparator", comparatorType, graceEnd = if (comparatorType == "neverInitiate") NULL else graceEnd)
  )
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = strategies, sameDayOrdering = sameDayOrdering)))
  models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(prior = Cyclops::createPrior("none"))))
  return(list(data = data, clones = clones, models = models))
}

unweightedArgs <- function() createComputeWeightsArgs(stabilized = FALSE, truncation = createWeightTruncation("none"))

# W(d) of a clone (unstabilized, "eventsFirst"; Phase1Contracts section 5.4) from the fitted censoring
# model: 1 up to day G and 1 / pi after it (binary), or 1 / prod_{k < d} (1 - h_k) (pooled logistic).
handWeight <- function(model, rowId, day) {
  day <- rep_len(day, length(rowId))
  if (model$kind == "binary") {
    pi <- model$predictions$propensityScore[match(rowId, model$predictions$rowId)]
    return(ifelse(day > model$horizon, 1 / pi, 1))
  }
  linearPredictor <- model$predictions$linearPredictor[match(rowId, model$predictions$rowId)]
  return(vapply(seq_along(day), function(i) {
    k <- seq_len(min(day[i], model$horizon + 1)) - 1
    1 / prod(1 - plogis(model$intervals$baseline[k + 1] + linearPredictor[i]))
  }, numeric(1)))
}

collectTable <- function(table, orderBy) {
  table <- as.data.frame(dplyr::collect(table))
  table <- table[do.call(order, table[orderBy]), ]
  rownames(table) <- NULL
  return(table)
}

# Amendment 1 ------------------------------------------------------------------------------------------

test_that("With deviationFirst, segments cover each clone's follow-up exactly", {
  # M2 ends a deviating clone on the day before its deviation (amendment 1), so M3's drop of the
  # deviation day is a no-op: every clone's segments run from 0 to its survivalTime.
  for (comparatorType in c("notInitiateWithin", "neverInitiate")) {
    analysis <- integrationAnalysis(comparatorType, sameDayOrdering = "deviationFirst")
    weighted <- suppressMessages(computeCensoringWeights(analysis$clones, analysis$models, unweightedArgs()))
    clones <- collectTable(weighted$clones, "cloneRowId")
    segments <- collectTable(weighted$segments, c("cloneRowId", "startTime"))
    expect_true(any(clones$endReason == "deviation"))
    expect_true(all(is.na(clones$deviationDay) | clones$deviationDay > clones$endDay))
    expect_equal(as.vector(tapply(segments$stopTime, segments$cloneRowId, max)), as.numeric(clones$survivalTime))
    expect_equal(as.vector(tapply(segments$startTime, segments$cloneRowId, min)), rep(0, nrow(clones)))
    close(analysis$clones)
    close(analysis$data)
  }
})

# Amendment 4 ------------------------------------------------------------------------------------------

test_that("The default censoring hazard model is the exact daily pooled logistic model", {
  args <- createFitCensoringModelsArgs()
  expect_equal(args$hazardModelType, "logistic")
  expect_equal(args$timeBasis$type, "indicator")
  expect_equal(args$timeBasis$intervalLength, 1)

  # Without covariates, the daily indicator model's S(d) is the product-limit survival of
  # non-initiation, here with many initiations on day 0, which exp(-Breslow) of the Cox model is not.
  persons <- data.frame(
    rowId = 1:10,
    personSeqId = 1:10,
    personId = as.character(1:10),
    cohortStartDate = as.Date("2020-01-01"),
    daysFromObsStart = 365L,
    daysToCohortEnd = 0L,
    daysToObsEnd = 10L,
    daysToDeath = NA_integer_
  )
  data <- createCloneCensorWeightData(
    cohorts = persons,
    outcomes = data.frame(rowId = 10L, outcomeId = 3, daysToEvent = 5L),
    treatments = data.frame(rowId = 1:6, treatmentId = 2, daysToStart = c(0L, 0L, 0L, 0L, 1L, 3L), daysToEnd = 10L),
    metaData = list(targetId = 1, treatmentIds = c(treatment = 2))
  )
  on.exit(close(data))
  population <- suppressMessages(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 10, endAnchor = "cohort start")))
  strategies <- list(
    createTreatmentStrategy(1, "Initiate", "initiateWithin", graceEnd = 3),
    createTreatmentStrategy(2, "Never", "neverInitiate")
  )
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = strategies)))
  on.exit(close(clones), add = TRUE)
  logistic <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(useCovariates = FALSE)))[[2]]
  expect_equal(logistic$modelType, "logistic")
  survival <- cumprod(1 - plogis(logistic$numerator$intervals$baseline))
  # At risk on days 0, 1, 2, 3: 10, 6, 5, 5 persons; initiations 4, 1, 0, 1.
  expect_equal(survival[1:4], cumprod(1 - c(4 / 10, 1 / 6, 0, 1 / 5)), tolerance = 1e-6)
  cox <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(hazardModelType = "cox", useCovariates = FALSE)))[[2]]
  expect_equal(exp(-.evaluateCumulativeHazard(cox$numerator$baselineHazard, 1)), exp(-4 / 10))
  expect_gt(exp(-4 / 10) - survival[1], 0.07)
})

# Amendment 5 ------------------------------------------------------------------------------------------

test_that("Clones are split at the end of day tau + 1, so balance reads W(tau + 1) exactly", {
  skip_if_not_installed("survival")
  configurations <- list(
    list(graceEnd = 5, riskStart = 0),
    list(graceEnd = 5, riskStart = 1),
    # tau + 1 = 1 is before riskStart = 2: the split is at the end of day riskStart (t = 1).
    list(graceEnd = 0, riskStart = 2)
  )
  for (configuration in configurations) {
    analysis <- integrationAnalysis("neverInitiate", graceEnd = configuration$graceEnd, riskStart = configuration$riskStart)
    weighted <- suppressMessages(computeCensoringWeights(analysis$clones, analysis$models, unweightedArgs()))
    clones <- collectTable(weighted$clones, "cloneRowId")
    segments <- collectTable(weighted$segments, c("cloneRowId", "startTime"))
    balanceDay <- max(configuration$graceEnd + 1, configuration$riskStart)
    cutTime <- balanceDay - configuration$riskStart + 1
    followed <- clones$cloneRowId[clones$survivalTime >= cutTime]
    expect_gt(length(followed), 100)
    expect_setequal(segments$cloneRowId[segments$stopTime == cutTime], followed)

    # The weight of the segment that contains day tau + 1 is W(tau + 1), also for neverInitiate:
    compliant <- clones[clones$endDay >= balanceDay, ]
    containing <- segments[segments$startTime < cutTime & segments$stopTime >= cutTime, ]
    compliant$weight <- containing$weight[match(compliant$cloneRowId, containing$cloneRowId)]
    expected <- numeric(nrow(compliant))
    for (strategyId in 1:2) {
      inStrategy <- compliant$strategyId == strategyId
      expected[inStrategy] <- handWeight(analysis$models[[strategyId]], compliant$rowId[inStrategy], balanceDay)
    }
    expect_equal(compliant$weight, expected, tolerance = 1e-10)
    expect_gt(length(unique(round(compliant$weight[compliant$strategyId == 2], 8))), 50)

    # computeCovariateBalance's weighted means use those weights:
    balance <- suppressMessages(computeCovariateBalance(weighted, analysis$data))
    covariates <- dplyr::collect(dplyr::filter(analysis$data$covariates, .data$covariateId == 1001))
    x <- covariates$covariateValue[match(compliant$rowId, covariates$rowId)]
    x[is.na(x)] <- 0
    row <- balance[balance$covariateId == 1001, ]
    for (strategyId in 1:2) {
      inStrategy <- compliant$strategyId == strategyId
      weightedMean <- sum(compliant$weight[inStrategy] * x[inStrategy]) / sum(compliant$weight[inStrategy])
      column <- if (strategyId == 1) "afterMatchingMeanTarget" else "afterMatchingMeanComparator"
      expect_equal(row[[column]], weightedMean, tolerance = 1e-10)
    }
    close(analysis$clones)
    close(analysis$data)
  }
})

test_that("With outcomeTimeBasis, clones are split on its grid, which leaves KM and Cox unchanged", {
  analysis <- integrationAnalysis("neverInitiate", riskStart = 1)
  on.exit({
    close(analysis$clones)
    close(analysis$data)
  })
  weighted <- suppressMessages(computeCensoringWeights(analysis$clones, analysis$models, unweightedArgs()))
  expect_true(is.na(attr(weighted, "metaData")$weights$cutIntervalLength))
  args <- lapply(c("km", "cox"), function(modelType) createFitOutcomeModelArgs(modelType = modelType, horizons = c(8, 20)))
  reference <- lapply(args, function(arg) suppressMessages(fitOutcomeModel(weighted, arg)))
  clones <- collectTable(weighted$clones, "cloneRowId")
  noGrid <- collectTable(weighted$segments, c("cloneRowId", "startTime"))

  for (intervalLength in c(3, 1)) {
    gridWeighted <- suppressMessages(computeCensoringWeights(
      analysis$clones, analysis$models, unweightedArgs(),
      outcomeTimeBasis = createTimeBasis(intervalLength = intervalLength)
    ))
    expect_equal(attr(gridWeighted, "metaData")$weights$cutIntervalLength, intervalLength)
    segments <- collectTable(gridWeighted$segments, c("cloneRowId", "startTime"))
    # Every segment lies within one interval, and every interval boundary inside a clone's follow-up
    # ends a segment:
    expect_equal(floor(segments$startTime / intervalLength), ceiling(segments$stopTime / intervalLength) - 1)
    boundaries <- do.call(rbind, lapply(seq_len(nrow(clones)), function(i) {
      t <- seq_len((clones$survivalTime[i] - 1) %/% intervalLength) * intervalLength
      if (length(t) == 0) NULL else data.frame(cloneRowId = clones$cloneRowId[i], stopTime = t)
    }))
    expect_true(all(paste(boundaries$cloneRowId, boundaries$stopTime) %in% paste(segments$cloneRowId, segments$stopTime)))
    # Each segment's weight is W on its last day (riskStart = 1):
    segments$day <- segments$stopTime + 1 - 1
    for (strategyId in 1:2) {
      inStrategy <- segments$strategyId == strategyId
      expected <- handWeight(analysis$models[[strategyId]], segments$rowId[inStrategy], segments$day[inStrategy])
      expect_equal(segments$weight[inStrategy], expected, tolerance = 1e-10)
    }
    if (intervalLength == 1) {
      expect_true(all(segments$stopTime - segments$startTime == 1))
      expect_equal(nrow(segments), sum(clones$survivalTime))
    }
    expect_equal(sum(segments$outcome), sum(noGrid$outcome))
    expect_equal(sum(segments$competing), sum(noGrid$competing))
    for (i in seq_along(args)) {
      gridModel <- suppressMessages(fitOutcomeModel(gridWeighted, args[[i]]))
      expect_equal(gridModel$estimates, reference[[i]]$estimates, tolerance = 1e-12)
    }
  }
})

# Amendment 6 ------------------------------------------------------------------------------------------

# Weighted clones with one segment per day and a weight that changes every day, w(t) = exp(a + b t) per
# clone, built from an M4 fixture's clones.
dailyWeightFixture <- function(riskStart, competingEvents, seed, n = 200, maxDays = 30, outcomeRate = 0.6, deathRate = 0.3) {
  fixture <- createM4Fixture(n = n, riskStart = riskStart, competingEvents = competingEvents, seed = seed, maxDays = maxDays, outcomeRate = outcomeRate, deathRate = deathRate)
  clones <- fixture$clones
  set.seed(seed)
  a <- stats::runif(nrow(clones), -0.2, 0.5)
  b <- stats::runif(nrow(clones), -0.02, 0.05)
  fixture$weightAt <- function(i, t) exp(a[i] + b[i] * t)
  segments <- lapply(seq_len(nrow(clones)), function(i) {
    t <- seq_len(clones$survivalTime[i])
    last <- t == clones$survivalTime[i]
    data.frame(
      cloneRowId = clones$cloneRowId[i],
      rowId = clones$rowId[i],
      personSeqId = clones$personSeqId[i],
      strategyId = clones$strategyId[i],
      startTime = as.numeric(t - 1),
      stopTime = as.numeric(t),
      weight = fixture$weightAt(i, t),
      outcome = as.integer(last & clones$endReason[i] == "outcome"),
      competing = as.integer(last & clones$endReason[i] == "death")
    )
  })
  fixture$segments <- cbind(segmentId = seq_len(sum(clones$survivalTime)), do.call(rbind, segments))
  fixture$personWeights <- m4PersonWeights(fixture$segments)
  return(fixture)
}

# Clone x interval rows built by hand from the fixture's weight function: the weight of a clone's last
# day in the interval.
handIntervalRows <- function(fixture, intervalLength) {
  clones <- fixture$clones
  rows <- lapply(seq_len(nrow(clones)), function(i) {
    endTime <- clones$survivalTime[i]
    k <- seq_len(ceiling(endTime / intervalLength)) - 1
    last <- k == max(k)
    data.frame(
      personSeqId = clones$personSeqId[i],
      x = as.integer(clones$strategyId[i] == 1),
      k = k,
      weight = fixture$weightAt(i, pmin((k + 1) * intervalLength, endTime)),
      y = as.integer(last & clones$endReason[i] == "outcome"),
      death = as.integer(last & clones$endReason[i] == "death")
    )
  })
  return(do.call(rbind, rows))
}

test_that("The pooled logistic model on daily weighted clones equals glm on hand-built daily rows", {
  control <- glm.control(epsilon = 1e-14, maxit = 100)
  for (riskStart in 0:1) {
    fixture <- dailyWeightFixture(riskStart, "death", seed = 11)
    weightedClones <- createM4WeightedClones(fixture, 1)
    horizons <- riskStart + c(9, 24)
    # The default time basis (natural spline with 4 df of daily intervals) and the default control:
    model <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(modelType = "logistic", horizons = horizons)))
    rows <- handIntervalRows(fixture, 1)
    basis <- splines::ns(rows$k, df = 4)
    noDeath <- rows$death == 0
    outcomeFit <- glm(rows$y[noDeath] ~ rows$x[noDeath] + basis[noDeath, ], family = quasibinomial(), weights = rows$weight[noDeath], control = control)
    deathFit <- glm(rows$death ~ rows$x + basis, family = quasibinomial(), weights = rows$weight, control = control)
    expect_equal(unname(model$coefficients), unname(c(coef(outcomeFit), coef(deathFit))), tolerance = 1e-8)
    for (h in horizons) {
      k <- seq_len(h - riskStart + 1) - 1
      newBasis <- predict(basis, k)
      for (x in 1:0) {
        hazardY <- plogis(as.vector(cbind(1, x, newBasis) %*% coef(outcomeFit)))
        hazardD <- plogis(as.vector(cbind(1, x, newBasis) %*% coef(deathFit)))
        estimand <- sprintf("risk_%s_%s", if (x == 1) "target" else "comparator", h)
        expect_equal(getM4Estimate(model, estimand), utils::tail(youngCif(hazardY, hazardD), 1), tolerance = 1e-8)
      }
    }

    close(weightedClones)

    # Intervals of 3 days: each row takes the weight of the clone's last day in the interval.
    weightedClones <- createM4WeightedClones(fixture, 3)
    model <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(
      modelType = "logistic", horizons = riskStart - 1 + c(9, 24), timeBasis = createTimeBasis("linear", intervalLength = 3), armByTimeInteraction = TRUE
    )))
    rows <- handIntervalRows(fixture, 3)
    noDeath <- rows$death == 0
    design <- cbind(1, rows$x, rows$k, rows$x * rows$k)
    outcomeFit <- glm(rows$y[noDeath] ~ design[noDeath, ] - 1, family = quasibinomial(), weights = rows$weight[noDeath], control = control)
    deathFit <- glm(rows$death ~ design - 1, family = quasibinomial(), weights = rows$weight, control = control)
    expect_equal(unname(model$coefficients), unname(c(coef(outcomeFit), coef(deathFit))), tolerance = 1e-8)
    close(weightedClones)
  }
})

test_that("A saturated daily pooled logistic model gives the weighted Aalen-Johansen risks", {
  # An indicator per day and strategy, on a fixture where every day of both strategies has deaths, and
  # outcomes but not only outcomes among the rows without death (Young et al. 2020).
  for (riskStart in 0:1) {
    fixture <- dailyWeightFixture(riskStart, "death", seed = 13, n = 1000, maxDays = 4, outcomeRate = 0.9, deathRate = 0.6)
    rows <- handIntervalRows(fixture, 1)
    cells <- aggregate(cbind(y, death, rows = 1) ~ x + k, data = rows, FUN = sum)
    outcomeCells <- aggregate(cbind(y, rows = 1) ~ x + k, data = rows[rows$death == 0, ], FUN = sum)
    expect_true(all(cells$death > 0 & cells$death < cells$rows) && all(outcomeCells$y > 0 & outcomeCells$y < outcomeCells$rows))
    weightedClones <- createM4WeightedClones(fixture, 1)
    horizons <- riskStart + c(2, 4)
    saturated <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(
      modelType = "logistic", horizons = horizons, timeBasis = createTimeBasis("indicator", intervalLength = 1), armByTimeInteraction = TRUE
    )))
    km <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(modelType = "km", horizons = horizons)))
    close(weightedClones)
    risks <- grepl("^risk_", saturated$estimates$estimand)
    expect_equal(sum(risks), 4)
    expect_equal(
      saturated$estimates$estimate[risks],
      km$estimates$estimate[match(saturated$estimates$estimand[risks], km$estimates$estimand)],
      tolerance = 1e-8
    )
  }
})

test_that("The pooled logistic model stops on weighted clones not split on its grid", {
  fixture <- dailyWeightFixture(0, "none", seed = 12)
  args <- createFitOutcomeModelArgs(modelType = "logistic", horizons = 10)
  for (cutIntervalLength in c(NA, 3)) {
    weightedClones <- createM4WeightedClones(fixture, cutIntervalLength)
    expect_error(suppressMessages(fitOutcomeModel(weightedClones, args)), "computeCensoringWeights\\(outcomeTimeBasis")
    expect_error(.fitOutcomeModelPoint(weightedClones, args), "intervals of 1 days")
    close(weightedClones)
  }
  # KM and Cox do not need the grid:
  weightedClones <- createM4WeightedClones(fixture)
  on.exit(close(weightedClones))
  expect_s3_class(suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(modelType = "km", horizons = 10))), "CcwOutcomeModel")

  defaults <- createFitOutcomeModelArgs()
  expect_equal(defaults$timeBasis$type, "naturalSpline")
  expect_equal(defaults$timeBasis$intervalLength, 1)
  expect_equal(defaults$timeBasis$df, 4)
})

# Amendment 7 ------------------------------------------------------------------------------------------

test_that("The default outcome-model control agrees with the references to 1e-8 for every model type", {
  skip_if_not_installed("survival")
  coxControl <- survival::coxph.control(eps = 1e-14, toler.chol = 1e-15, iter.max = 100)
  glmControl <- glm.control(epsilon = 1e-15, maxit = 100)
  for (riskStart in 0:1) {
    fixture <- createM4Fixture(n = 200, riskStart = riskStart, competingEvents = "death", seed = 3, maxDays = 30, outcomeRate = 0.7)
    data <- toSurvivalData(fixture$segments)

    # Cox: the cause-specific outcome model vs coxph.
    weightedClones <- createM4WeightedClones(fixture)
    cox <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(modelType = "cox", horizons = 20)))
    causeSpecific <- survival::coxph(
      survival::Surv(startTime, stopTime, outcome) ~ x,
      data = data, weights = weight, ties = "breslow", control = coxControl
    )
    expect_equal(unname(cox$coefficients), unname(coef(causeSpecific)), tolerance = 1e-8)

    # KM / Aalen-Johansen (no model fit) vs survfit:
    km <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(modelType = "km", horizons = riskStart + c(9, 24))))
    for (strategyId in 1:2) {
      fit <- survival::survfit(survival::Surv(startTime, stopTime, state) ~ 1, data = data[data$strategyId == strategyId, ], weights = weight, id = cloneRowId)
      role <- c("target", "comparator")[strategyId]
      expect_equal(
        c(getM4Estimate(km, sprintf("risk_%s_%d", role, riskStart + 9)), getM4Estimate(km, sprintf("risk_%s_%d", role, riskStart + 24))),
        summary(fit, times = c(10, 25))$pstate[, 2],
        tolerance = 1e-8
      )
    }
    close(weightedClones)

    # Pooled logistic, the default basis and an ill-scaled raw cubic with the interaction:
    rows <- expandM4Rows(fixture$segments, 1)
    noDeath <- rows$death == 0
    for (timeBasis in list(createTimeBasis("naturalSpline", intervalLength = 1, df = 4), createTimeBasis("polynomial", intervalLength = 1, degree = 3, raw = TRUE))) {
      interaction <- timeBasis$type == "polynomial"
      weightedClones <- createM4WeightedClones(fixture, 1)
      model <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(modelType = "logistic", horizons = 20, timeBasis = timeBasis, armByTimeInteraction = interaction)))
      close(weightedClones)
      basis <- if (interaction) cbind(rows$k, rows$k^2, rows$k^3) else unclass(splines::ns(rows$k, df = 4))
      design <- cbind(1, rows$x, basis)
      if (interaction) {
        design <- cbind(design, rows$x * basis)
      }
      outcomeFit <- glm(rows$y[noDeath] ~ design[noDeath, ] - 1, family = quasibinomial(), weights = rows$weight[noDeath], control = glmControl)
      deathFit <- glm(rows$death ~ design - 1, family = quasibinomial(), weights = rows$weight, control = glmControl)
      expect_equal(unname(model$coefficients), unname(c(coef(outcomeFit), coef(deathFit))), tolerance = 1e-8)
    }
  }
  control <- createFitOutcomeModelArgs()$control
  expect_equal(control$algorithm, "newton")
  expect_equal(control$convergenceType, "gradient")
  expect_equal(control$tolerance, 1e-10)
})

# Amendment 8 ------------------------------------------------------------------------------------------

test_that("computeTrueRisks with riskStart conditions on no outcome or death before riskStart", {
  data <- simulateCloneCensorWeightData(n = 10, dgp = "a", graceEnd = 5, maxDays = 30, outcomeHazard = 0.01, deathHazard = 0.005, seed = 1)
  settings <- attr(data, "metaData")$simulationSettings
  close(data)
  strategies <- list(
    createTreatmentStrategy(1, "Initiate", "initiateWithin", graceEnd = 5),
    createTreatmentStrategy(2, "Never", "neverInitiate")
  )
  nTruth <- 2e5
  fromZero <- computeTrueRisks(settings, strategies, horizons = c(0, 1, 2, 10, 30), nTruth = nTruth, seed = 2)
  expect_identical(computeTrueRisks(settings, strategies, horizons = c(0, 1, 2, 10, 30), nTruth = nTruth, seed = 2, riskStart = 0), fromZero)

  # The conversion (risk(h) - risk(r - 1)) / (1 - p), with p from the untreated hazards of the same X:
  set.seed(2)
  x <- .simulateCovariates(nTruth, settings)
  eta <- .linearPredictors(x, settings)
  noEventDay <- (1 - plogis(eta$outcome)) * (1 - plogis(eta$death))
  for (riskStart in 1:2) {
    p <- mean(1 - noEventDay^riskStart)
    conditional <- computeTrueRisks(settings, strategies, horizons = c(10, 30), nTruth = nTruth, seed = 2, riskStart = riskStart)
    for (strategyId in 1:2) {
      risk <- fromZero$risk[fromZero$strategyId == strategyId]
      expected <- (risk[4:5] - risk[riskStart]) / (1 - p)
      expect_equal(conditional$risk[conditional$strategyId == strategyId], expected, tolerance = 1e-12)
    }
  }

  # Against a direct simulation of the never-treated outcome process, among persons with no outcome or
  # death before day 2 (within 4 Monte Carlo standard errors):
  set.seed(3)
  x <- .simulateCovariates(nTruth, settings)
  eta <- .linearPredictors(x, settings)
  outcomeDay <- stats::rgeom(nTruth, plogis(eta$outcome))
  deathDay <- stats::rgeom(nTruth, plogis(eta$death))
  atRisk <- pmin(outcomeDay, deathDay) >= 2
  for (h in c(10, 30)) {
    direct <- mean(outcomeDay[atRisk] <= h & outcomeDay[atRisk] <= deathDay[atRisk])
    truth <- computeTrueRisks(settings, strategies[2], horizons = h, nTruth = nTruth, seed = 4, riskStart = 2)$risk
    expect_lt(abs(direct - truth), 4 * sqrt(2 * truth * (1 - truth) / sum(atRisk)))
  }

  expect_error(computeTrueRisks(settings, strategies, horizons = 10, nTruth = 100, riskStart = 7), "riskStart = 7 is not supported")
  settingsB <- settings
  settingsB$dgp <- "b"
  expect_error(computeTrueRisks(settingsB, strategies, horizons = 10, nTruth = 100, riskStart = 2), "riskStart = 2 is not supported")
  expect_error(computeTrueRisks(settings, strategies, horizons = c(1, 10), nTruth = 100, riskStart = 2), "Horizons must be on or after riskStart")
})

# Amendment 11 -----------------------------------------------------------------------------------------

test_that("The censoring and outcome models build the same time basis for the same grid and counts", {
  for (nIntervals in c(1, 2, 3, 8, 30)) {
    k <- seq_len(nIntervals) - 1
    for (timeBasis in list(
      createTimeBasis("indicator"),
      createTimeBasis("linear"),
      createTimeBasis("polynomial", degree = 3, raw = TRUE),
      createTimeBasis("polynomial", degree = 2, raw = FALSE),
      createTimeBasis("naturalSpline", df = 4)
    )) {
      # The censoring models' basis on the grid (each interval once) is the shared basis of rows with
      # one row per interval, the outcome model's construction:
      censoring <- .createTimeBasisMatrix(timeBasis, nIntervals)
      outcome <- .timeBasisColumns(timeBasis, k)$build(k)
      expect_equal(unname(censoring), unname(outcome))
      expect_equal(ncol(censoring), if (nIntervals == 1) 0 else switch(timeBasis$type,
        indicator = nIntervals - 1,
        linear = 1,
        polynomial = min(timeBasis$degree, nIntervals - 1),
        naturalSpline = min(timeBasis$df, nIntervals - 1)
      ))
      if (nIntervals > 1 && timeBasis$type != "indicator") {
        reference <- switch(timeBasis$type,
          linear = matrix(k),
          polynomial = stats::poly(k, degree = min(timeBasis$degree, nIntervals - 1), raw = timeBasis$raw),
          naturalSpline = splines::ns(k, df = min(timeBasis$df, nIntervals - 1))
        )
        expect_equal(unname(censoring), unname(matrix(as.numeric(reference), nrow = nIntervals)))
      }
    }
  }

  # Counts: rows counted `multiplicity` times give the basis of the rows copied, whose knots (natural
  # spline) and centring (orthogonal polynomial) follow the distribution of the rows' intervals:
  set.seed(1)
  interval <- sample(0:19, 200, replace = TRUE, prob = 20:1)
  multiplicity <- sample(0:3, 200, replace = TRUE)
  copied <- rep(interval, multiplicity)
  for (timeBasis in list(createTimeBasis("naturalSpline", df = 4), createTimeBasis("polynomial", degree = 2, raw = FALSE))) {
    shared <- .timeBasisColumns(timeBasis, interval, multiplicity)$build(0:19)
    reference <- if (timeBasis$type == "naturalSpline") splines::ns(copied, df = 4) else stats::poly(copied, degree = 2)
    expect_equal(unname(shared), unname(unclass(predict(reference, 0:19))[, , drop = FALSE]), tolerance = 1e-12)
    expect_equal(unname(.timeBasisColumns(timeBasis, copied)$build(0:19)), unname(shared), tolerance = 1e-12)
  }
  # The indicator basis of the censoring models keeps only the given intervals:
  censoring <- .createTimeBasisMatrix(createTimeBasis("indicator"), 6, c(0, 2, 5))
  expect_equal(colnames(censoring), c("timeBasis2", "timeBasis5"))
  expect_equal(unname(censoring), unname(.timeBasisColumns(createTimeBasis("indicator"), c(0, 2, 5))$build(0:5)))
})

# Review of the integrated package --------------------------------------------------------------------

# A dailyWeightFixture whose deaths (`reason = "death"`) or outcomes (`"outcome"`) are kept only on the
# given days (since index); the other clones with that end reason end follow-up as the end of
# observation instead.
keepEventsOn <- function(fixture, reason, days) {
  drop <- fixture$clones$endReason == reason & !(fixture$clones$endDay %in% days)
  fixture$clones$endReason[drop] <- "observationEnd"
  fixture$clones$outcomeCount[drop] <- 0L
  column <- if (reason == "outcome") "outcome" else "competing"
  fixture$segments[[column]][fixture$segments$cloneRowId %in% fixture$clones$cloneRowId[drop]] <- 0L
  return(fixture)
}

test_that("The default control fits a weighted Cox MSM on a million segments as coxph does", {
  skip_on_cran()
  skip_if_not_installed("survival")
  # The clones of 20,000 simulated persons split on the daily grid (about 1e6 segments). Newton steps
  # stopped by the standardized-gradient ("kkt") criterion reached maxIterations at the optimum here.
  data <- simulateCloneCensorWeightData(
    n = 20000, graceEnd = 10, maxDays = 60, initiationHazard = 0.03, outcomeHazard = 0.003, deathHazard = 0.001,
    observationEndHazard = 0.003, initiationCoefficients = rep(0.5, 5), outcomeCoefficients = rep(0.5, 5),
    deathCoefficients = rep(0.3, 5), treatmentEffect = 0, seed = 1
  )
  on.exit(close(data))
  population <- suppressMessages(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
    riskWindowEnd = 60, endAnchor = "cohort start"
  )))
  strategies <- list(
    createTreatmentStrategy(1, "Initiate", "initiateWithin", graceEnd = 10),
    createTreatmentStrategy(2, "Defer", "notInitiateWithin", graceEnd = 10)
  )
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = strategies, sameDayOrdering = "deviationFirst")))
  on.exit(close(clones), add = TRUE)
  models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(hazardModelType = "cox", prior = Cyclops::createPrior("none"))))
  weightedClones <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(truncation = createWeightTruncation("none")), outcomeTimeBasis = createTimeBasis("indicator", 1)))
  on.exit(close(weightedClones), add = TRUE)
  segments <- as.data.frame(dplyr::collect(weightedClones$segments))
  expect_gt(nrow(segments), 1e6)

  expect_silent(point <- .fitOutcomeModelPoint(weightedClones, createFitOutcomeModelArgs(modelType = "cox", horizons = 60)))
  segments$x <- as.integer(segments$strategyId == 1)
  reference <- survival::coxph(
    survival::Surv(startTime, stopTime, outcome) ~ x,
    data = segments[segments$weight > 0, ], weights = weight, ties = "breslow",
    control = survival::coxph.control(eps = 1e-12, toler.chol = 1e-14, iter.max = 50)
  )
  expect_equal(point$estimates$logEstimate[point$estimates$estimand == "hr"], unname(coef(reference)), tolerance = 1e-8)
  expect_false(anyNA(point$estimates$estimate))
})

test_that("The default control fits a quasi-separated death model, whose risks equal glm", {
  # Deaths on days 1, 2 and 10 only: the natural spline death model is quasi-separated (its hazard tends
  # to 0 away from those days). A standardized-gradient criterion was not met at its optimum.
  glmControl <- glm.control(epsilon = 1e-14, maxit = 100)
  for (riskStart in 0:1) {
    fixture <- keepEventsOn(dailyWeightFixture(riskStart, "death", seed = 3, n = 150, maxDays = 15, outcomeRate = 0.5, deathRate = 0.3), "death", riskStart + c(1, 2, 10))
    rows <- handIntervalRows(fixture, 1)
    expect_equal(sort(unique(rows$k[rows$death == 1])), c(1, 2, 10))
    weightedClones <- createM4WeightedClones(fixture, 1)
    horizons <- riskStart + c(7, 14)
    expect_silent(model <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(modelType = "logistic", horizons = horizons))))
    close(weightedClones)
    basis <- splines::ns(rows$k, df = 4)
    noDeath <- rows$death == 0
    outcomeFit <- glm(rows$y[noDeath] ~ rows$x[noDeath] + basis[noDeath, ], family = quasibinomial(), weights = rows$weight[noDeath], control = glmControl)
    deathFit <- glm(rows$death ~ rows$x + basis, family = quasibinomial(), weights = rows$weight, control = glmControl)
    expect_true(deathFit$converged)
    # The outcome model is not separated, so its coefficients agree; the death model's agree only where
    # its likelihood is not flat, which its hazards and so the risks show:
    expect_equal(unname(model$coefficients[1:6]), unname(coef(outcomeFit)), tolerance = 1e-8)
    for (h in horizons) {
      newBasis <- predict(basis, seq_len(h - riskStart + 1) - 1)
      for (x in 1:0) {
        hazardY <- plogis(as.vector(cbind(1, x, newBasis) %*% coef(outcomeFit)))
        hazardD <- plogis(as.vector(cbind(1, x, newBasis) %*% coef(deathFit)))
        estimand <- sprintf("risk_%s_%s", if (x == 1) "target" else "comparator", h)
        expect_equal(getM4Estimate(model, estimand), utils::tail(youngCif(hazardY, hazardD), 1), tolerance = 1e-8)
      }
    }
    expect_false(anyNA(model$estimates$estimate))
    expect_false(anyNA(model$estimates$se[grepl("^(risk|rd)_", model$estimates$estimand)]))
  }
})

test_that("Indicator bases give intervals or cells without events a hazard of 0", {
  skip_if_not_installed("survival")
  for (riskStart in 0:1) {
    # Saturated (indicator x strategy): some strategy x day cells have no deaths. Risks and their robust
    # standard errors equal the weighted Aalen-Johansen ones.
    fixture <- dailyWeightFixture(riskStart, "death", seed = 13, n = 300, maxDays = 8, outcomeRate = 0.6, deathRate = 0.25)
    rows <- handIntervalRows(fixture, 1)
    cells <- aggregate(cbind(death, rows = 1) ~ x + k, data = rows, FUN = sum)
    expect_true(any(cells$death == 0))
    weightedClones <- createM4WeightedClones(fixture, 1)
    horizons <- riskStart + c(3, 8)
    expect_silent(saturated <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(
      modelType = "logistic", horizons = horizons, timeBasis = createTimeBasis("indicator", intervalLength = 1), armByTimeInteraction = TRUE
    ))))
    km <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(modelType = "km", horizons = horizons)))
    close(weightedClones)
    contrasts <- grepl("^(risk|rd)_", saturated$estimates$estimand)
    expect_equal(sum(contrasts), 6)
    reference <- km$estimates[match(saturated$estimates$estimand[contrasts], km$estimates$estimand), ]
    expect_equal(saturated$estimates$estimate[contrasts], reference$estimate, tolerance = 1e-8)
    expect_equal(saturated$estimates$se[contrasts], reference$se, tolerance = 1e-6)

    # Without the interaction and with outcomes on only a few days: the other days have a hazard of 0,
    # and the model is the glm of the days with outcomes.
    fixture <- keepEventsOn(dailyWeightFixture(riskStart, "none", seed = 14, n = 200, maxDays = 20, outcomeRate = 0.5), "outcome", riskStart + c(2, 5, 6, 11, 17))
    rows <- handIntervalRows(fixture, 1)
    eventIntervals <- sort(unique(rows$k[rows$y == 1]))
    expect_lt(length(eventIntervals), 10)
    weightedClones <- createM4WeightedClones(fixture, 1)
    horizons <- riskStart + c(9, 19)
    expect_silent(model <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(
      modelType = "logistic", horizons = horizons, timeBasis = createTimeBasis("indicator", intervalLength = 1)
    ))))
    close(weightedClones)
    fitRows <- rows[rows$k %in% eventIntervals, ]
    referenceFit <- glm(y ~ x + factor(k), family = quasibinomial(), weights = weight, data = fitRows, control = glm.control(epsilon = 1e-14, maxit = 100))
    expect_equal(model$coefficients[["target"]], unname(coef(referenceFit)[["x"]]), tolerance = 1e-8)
    for (h in horizons) {
      k <- seq_len(h - riskStart + 1) - 1
      for (x in 1:0) {
        hazard <- numeric(length(k))
        hazard[k %in% eventIntervals] <- predict(referenceFit, newdata = data.frame(x = x, k = k[k %in% eventIntervals]), type = "response")
        estimand <- sprintf("risk_%s_%s", if (x == 1) "target" else "comparator", h)
        expect_equal(getM4Estimate(model, estimand), 1 - prod(1 - hazard), tolerance = 1e-8)
      }
    }
  }
})

test_that("A fit that fails with the given control is refitted with Cyclops' default control", {
  skip_if_not_installed("survival")
  fixture <- createM4Fixture(n = 200, riskStart = 0, competingEvents = "none", seed = 3)
  weightedClones <- createM4WeightedClones(fixture)
  on.exit(close(weightedClones))
  control <- Cyclops::createControl(algorithm = "newton", convergenceType = "gradient", tolerance = 1e-10, maxIterations = 1, noiseLevel = "silent")
  expect_warning(
    model <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(modelType = "cox", horizons = 10, control = control))),
    "Fitting the outcome Cox model failed \\(MAX_ITERATIONS\\), so it was fitted with Cyclops' default"
  )
  data <- toSurvivalData(fixture$segments)
  reference <- survival::coxph(survival::Surv(startTime, stopTime, outcome) ~ x, data = data, weights = weight, ties = "breslow")
  expect_equal(model$coefficients[["target"]], unname(coef(reference)), tolerance = 1e-6)
})

test_that("When only the pooled logistic death model fails, the odds ratio is kept and the risks are NA", {
  # (The Cox model fits no death model: its risks are those of the weighted Aalen-Johansen.)
  fixture <- createM4Fixture(n = 200, riskStart = 0, competingEvents = "death", seed = 3, maxDays = 15, outcomeRate = 0.7)
  full <- list(
    logistic = fitM4(fixture, modelType = "logistic", horizons = c(5, 10), timeBasis = createTimeBasis("linear", intervalLength = 1))
  )
  realPooled <- .fitPooledModel
  local_mocked_bindings(
    .fitPooledModel = function(rows, y, design, prior, control, label, variance, groups = NULL) {
      if (label == "death") NULL else realPooled(rows, y, design, prior, control, label, variance, groups)
    }
  )
  for (modelType in names(full)) {
    model <- fitM4(fixture, modelType = modelType, horizons = c(5, 10), timeBasis = createTimeBasis("linear", intervalLength = 1))
    ratio <- model$estimates[1, ]
    expect_equal(ratio, full[[modelType]]$estimates[1, ], tolerance = 1e-10)
    expect_false(is.na(ratio$estimate) || is.na(ratio$seLogEstimate))
    expect_true(all(is.na(model$estimates$estimate[-1])))
    expect_equal(model$estimates$estimand, full[[modelType]]$estimates$estimand)
    death <- startsWith(names(model$coefficients), "death:")
    expect_equal(names(model$coefficients), names(full[[modelType]]$coefficients))
    expect_true(all(is.na(model$coefficients[death])))
    expect_equal(model$coefficients[!death], full[[modelType]]$coefficients[!death], tolerance = 1e-10)
    expect_true(all(is.na(model$vcov[death, ])))
    expect_equal(model$vcov[!death, !death], full[[modelType]]$vcov[!death, !death], tolerance = 1e-8)
  }
})

test_that("The target and time-basis coefficients of the outcome models are not regularized", {
  fixture <- createM4Fixture(n = 200, riskStart = 0, competingEvents = "death", seed = 3, maxDays = 15, outcomeRate = 0.7)
  prior <- Cyclops::createPrior("laplace", variance = 0.01)
  for (modelType in c("cox", "logistic")) {
    timeBasis <- createTimeBasis("naturalSpline", intervalLength = 1, df = 3)
    unpenalized <- fitM4(fixture, modelType = modelType, horizons = c(5, 10), timeBasis = timeBasis)
    penalized <- fitM4(fixture, modelType = modelType, horizons = c(5, 10), timeBasis = timeBasis, prior = prior)
    expect_equal(penalized$coefficients, unpenalized$coefficients, tolerance = 1e-8)
    expect_equal(penalized$estimates, unpenalized$estimates, tolerance = 1e-6)
  }
})

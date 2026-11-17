library(testthat)

# Contract section 11.3: covariateTimeInteraction = "day0" adds covariate x I(day 0) columns (daily
# intervals only), and "interval" adds covariate x time-basis columns, to the pooled logistic initiation
# hazard model; they are penalized like the covariates, and the time-basis columns are not.

# Hand-built design of the pooled logistic hazard model on K_g rows (indicator basis): intervals without
# initiations, or where everyone at risk initiates, are left out; the first remaining interval is the
# reference. Returns the rows, the design matrix (named as the model's coefficients) and the intervals.
interactionDesign <- function(fixture, horizon, sameDayOrdering, intervalLength, covariateTimeInteraction) {
  rows <- refIntervalRows(fixture$persons, horizon, sameDayOrdering, intervalLength)
  initiations <- tapply(rows$initiated, rows$k, sum)
  atRisk <- tapply(rows$initiated, rows$k, length)
  intervals <- as.numeric(names(initiations))
  fitIntervals <- intervals[initiations > 0 & initiations < atRisk]
  fullIntervals <- intervals[initiations > 0 & initiations == atRisk]
  rows <- rows[rows$k %in% fitIntervals, ]
  basis <- sapply(fitIntervals[-1], function(k) as.numeric(rows$k == k))
  basis <- matrix(basis, nrow = nrow(rows), dimnames = list(NULL, sprintf("timeBasis%d", fitIntervals[-1])))
  x <- m3CovariateMatrix(fixture, rows$rowId)
  colnames(x) <- c("1001", "1002", "1003")
  design <- cbind(basis, x)
  if (covariateTimeInteraction == "day0") {
    interactions <- x * as.numeric(rows$k == 0)
    colnames(interactions) <- paste0(colnames(x), ":day0")
    design <- cbind(design, interactions)
  } else if (covariateTimeInteraction == "interval") {
    for (column in colnames(basis)) {
      interactions <- x * basis[, column]
      colnames(interactions) <- paste0(colnames(x), ":", column)
      design <- cbind(design, interactions)
    }
  }
  return(list(rows = rows, design = design, fitIntervals = fitIntervals, fullIntervals = fullIntervals))
}

# W of the strategy 2 segments from the glm coefficients: 1 / prod over the intervals ending on or before
# min(day - 1, horizon) (day with deviationFirst) of 1 - expit(eta_k), with eta_k the linear predictor
# of the person in interval k (hazard 0 in left-out intervals without initiations, 1 in full ones).
interactionReferenceWeights <- function(segments, fixture, design, coefficients, horizon, intervalLength, sameDayOrdering) {
  offset <- if (sameDayOrdering == "eventsFirst") 1 else 0
  nIntervals <- horizon %/% intervalLength + 1
  intervalEnd <- pmin((seq_len(nIntervals)) * intervalLength - 1, horizon)
  s2 <- which(segments$strategyId == 2)
  riskStart <- fixture$persons$riskStart[match(segments$rowId[s2], fixture$persons$rowId)]
  u <- pmin(segments$stopTime[s2] + riskStart - 1 - offset, horizon)
  x <- m3CovariateMatrix(fixture, segments$rowId[s2])
  return(vapply(seq_along(s2), function(i) {
    logSurvival <- 0
    for (k in seq_len(nIntervals) - 1) {
      if (intervalEnd[k + 1] > u[i]) {
        break
      }
      if (k %in% design$fullIntervals) {
        return(Inf)
      }
      if (!k %in% design$fitIntervals) {
        next
      }
      row <- c(as.numeric(design$fitIntervals[-1] == k), x[i, ])
      names(row) <- c(sprintf("timeBasis%d", design$fitIntervals[-1]), "1001", "1002", "1003")
      eta <- coefficients[["(Intercept)"]] + sum(row * coefficients[names(row)])
      day0 <- grepl(":day0$", names(coefficients))
      if (any(day0)) {
        eta <- eta + (k == 0) * sum(x[i, ] * coefficients[day0])
      }
      for (column in names(row)[startsWith(names(row), "timeBasis")]) {
        interaction <- paste0(c("1001", "1002", "1003"), ":", column)
        if (all(interaction %in% names(coefficients))) {
          eta <- eta + row[[column]] * sum(x[i, ] * coefficients[interaction])
        }
      }
      logSurvival <- logSurvival + log(1 - plogis(eta))
    }
    exp(-logSurvival)
  }, numeric(1)))
}

test_that("Interaction coefficients and weights match glm on the expanded rows without a penalty", {
  configurations <- list(
    list(covariateTimeInteraction = "day0", sameDayOrdering = "eventsFirst", intervalLength = 1),
    list(covariateTimeInteraction = "day0", sameDayOrdering = "deviationFirst", intervalLength = 1),
    list(covariateTimeInteraction = "interval", sameDayOrdering = "eventsFirst", intervalLength = 2),
    list(covariateTimeInteraction = "interval", sameDayOrdering = "deviationFirst", intervalLength = 1)
  )
  fixture <- createM3Fixture(1500, riskStart = 0L, seed = 41)
  data <- createM3Data(fixture)
  on.exit(close(data))
  horizon <- 3
  # The default tier runs two configurations that cover each interaction type, ordering and interval
  # length:
  for (configuration in testConfigurations(configurations, configurations[2:3])) {
    label <- paste(configuration, collapse = ", ")
    sameDayOrdering <- configuration$sameDayOrdering
    clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = horizon), sameDayOrdering)
    args <- m3Args("logistic", timeBasis = createTimeBasis(intervalLength = configuration$intervalLength))
    args$covariateTimeInteraction <- configuration$covariateTimeInteraction
    models <- suppressMessages(fitCensoringModels(clones, data, args))
    model <- models[[2]]
    design <- interactionDesign(fixture, horizon, sameDayOrdering, configuration$intervalLength, configuration$covariateTimeInteraction)
    reference <- glm(design$rows$initiated ~ design$design, family = binomial(), control = glm.control(epsilon = 1e-14, maxit = 100))
    expected <- coef(reference)
    names(expected) <- c("(Intercept)", colnames(design$design))
    expect_setequal(names(model$coefficients), names(expected))
    expect_equal(model$coefficients[names(expected)], expected, tolerance = 1e-6, label = label)
    expect_equal(model$interactions$type, configuration$covariateTimeInteraction)
    expect_equal(ncol(model$interactions$basis), if (configuration$covariateTimeInteraction == "day0") 1 else length(design$fitIntervals) - 1)

    # The weights, also after G (postGraceModel = "hazard"), from the same linear predictors:
    weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = FALSE, truncation = createWeightTruncation("none"))))
    segments <- collectSegments(weighted)
    expectedWeights <- interactionReferenceWeights(segments, fixture, design, expected, horizon, configuration$intervalLength, sameDayOrdering)
    expect_equal(segments$weight[segments$strategyId == 2], expectedWeights, tolerance = 1e-6, label = label)
    # 1 - S(G | X) for the calibration diagnostic uses the same linear predictors:
    afterGrace <- segments$strategyId == 2 & segments$stopTime + 0 - 1 - (sameDayOrdering == "eventsFirst") >= horizon
    survival <- .computeHazardSurvival(model, segments$rowId[afterGrace], horizon)
    expect_equal(1 / survival, segments$weight[afterGrace], tolerance = 1e-10)
    close(weighted)
    close(clones)
  }
  # Without interactions, the model and its interaction element are as before:
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = horizon), "eventsFirst")
  model <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic")))[[2]]
  expect_equal(model$interactions$type, "none")
  expect_equal(ncol(model$interactions$basis), 0)
  expect_equal(nrow(model$interactions$predictions), 0)
  expect_false(any(grepl(":", names(model$coefficients))))
  close(clones)
})

test_that("Interaction columns are penalized and time-basis columns are not", {
  fixture <- createM3Fixture(400, riskStart = 0L, seed = 42)
  data <- createM3Data(fixture)
  on.exit(close(data))
  # Hazard strategies only (no binary model, whose all-zero fit would warn):
  clones <- buildM3Clones(fixture$persons, m3Strategies("notInitiateWithin", "neverInitiate", graceEnd = 3), "eventsFirst")
  on.exit(close(clones), add = TRUE)
  # A very strong Laplace prior sets every penalized coefficient to 0, so the fit is the time-only model:
  args <- m3Args("logistic", prior = Cyclops::createPrior("laplace", variance = 1e-8, exclude = c(0)))
  args$covariateTimeInteraction <- "interval"
  model <- suppressMessages(fitCensoringModels(clones, data, args))[[1]]
  interactions <- grepl(":", names(model$coefficients))
  covariates <- names(model$coefficients) %in% c("1001", "1002", "1003")
  expect_true(sum(interactions) > 0)
  expect_equal(unname(model$coefficients[interactions | covariates]), rep(0, sum(interactions | covariates)))
  # The time basis is not penalized: the fit equals the numerator model.
  basis <- startsWith(names(model$coefficients), "timeBasis") | names(model$coefficients) == "(Intercept)"
  expect_equal(model$coefficients[basis], model$numerator$coefficients[names(model$coefficients)[basis]], tolerance = 1e-5)
  expect_true(all(model$coefficients[startsWith(names(model$coefficients), "timeBasis")] != 0))
})

test_that("covariateTimeInteraction = 'day0' with longer intervals multiplies the first interval", {
  # Contract 12.1 (superseding the review of 11.3, which stopped here): "day0" is the default, so it is
  # accepted with longer intervals, where it is the indicator of the first interval (days 0..L - 1).
  for (intervalLength in c(2, 7)) {
    args <- createFitCensoringModelsArgs(timeBasis = createTimeBasis(intervalLength = intervalLength), covariateTimeInteraction = "day0")
    expect_equal(args$covariateTimeInteraction, "day0")
    basis <- .createInteractionBasis("day0", matrix(0, nrow = 4, ncol = 3), 0:3, rowIntervals = 0:3)
    expect_equal(unname(basis[, "day0"]), c(1, 0, 0, 0))
  }
  expect_equal(createFitCensoringModelsArgs(covariateTimeInteraction = "day0")$covariateTimeInteraction, "day0")
  expect_equal(createFitCensoringModelsArgs(timeBasis = createTimeBasis(intervalLength = 7), covariateTimeInteraction = "interval")$covariateTimeInteraction, "interval")
})

test_that("covariateTimeInteraction is not available for the Cox model", {
  for (covariateTimeInteraction in c("day0", "interval")) {
    error <- tryCatch(createFitCensoringModelsArgs(hazardModelType = "cox", covariateTimeInteraction = covariateTimeInteraction), error = function(e) conditionMessage(e))
    error <- gsub("\\s*\\n\\s*\\*\\s*", " ", error)
    expect_match(error, "covariateTimeInteraction must be \"none\" with hazardModelType = \"cox\"", fixed = TRUE)
    expect_match(error, "available only in the pooled logistic initiation hazard model", fixed = TRUE)
  }
})

test_that("With day-0 initiation effects, the 'day0' model's 1 - S(G | X) is calibrated by X stratum", {
  # The covariate makes persons start on day 0 (log odds ratio 3) and less often later (-1.5). Under
  # DGP (a), starting during the grace period does not change the outcome or death hazards, so the
  # share who started by G among persons under follow-up at G estimates 1 - S(G | X).
  data <- simulateCloneCensorWeightData(
    n = 10000, nCovariates = 1, dgp = "a", graceEnd = 7, maxDays = 20, binaryPrevalence = 0.4,
    initiationHazard = 0.06, initiationCoefficients = -1.5, initiationHazardDay0 = 0.1, initiationCoefficientsDay0 = 3,
    outcomeHazard = 0.005, deathHazard = 0.002, outcomeCoefficients = 0.5, deathCoefficients = 0.3, seed = 1
  )
  on.exit(close(data))
  population <- suppressMessages(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
    riskWindowEnd = 20, endAnchor = "cohort start"
  )))
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = list(
    createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
    createTreatmentStrategy(2, "Do not initiate within 7 days", "notInitiateWithin", graceEnd = 7)
  ))))
  on.exit(close(clones), add = TRUE)
  # Observed by X stratum, by hand: persons with no outcome or death on or before day 7 and riskEnd >= 7.
  persons <- as.data.frame(dplyr::collect(clones$persons))
  x <- as.numeric(persons$rowId %in% dplyr::pull(dplyr::filter(data$covariates, .data$covariateId == 1001), "rowId"))
  eventDay <- pmin(persons$outcomeDay, persons$deathDay, na.rm = TRUE)
  atG <- (is.na(eventDay) | eventDay > 7) & persons$riskEnd >= 7
  started <- !is.na(persons$treatmentDay) & persons$treatmentDay <= 7
  observed <- tapply(started[atG], x[atG], mean)
  se <- sqrt(observed * (1 - observed) / as.vector(table(x[atG])))

  # 1 - S(7 | X) by hand from the coefficients: daily intervals 0..7, X effect beta on every day plus
  # delta on day 0.
  predicted <- function(model) {
    coefficients <- model$coefficients
    delta <- if ("1001:day0" %in% names(coefficients)) coefficients[["1001:day0"]] else 0
    sapply(0:1, function(xValue) {
      eta <- model$intervals$baseline + coefficients[["1001"]] * xValue + delta * xValue * (model$intervals$k == 0)
      1 - prod(1 - plogis(eta))
    })
  }
  z <- list()
  for (covariateTimeInteraction in c("none", "day0")) {
    models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(
      prior = Cyclops::createPrior("none"),
      covariateTimeInteraction = covariateTimeInteraction
    )))
    z[[covariateTimeInteraction]] <- (observed - predicted(models[[2]])) / se
    # The diagnostic reports the same predictions and observed shares (its groups are the X strata):
    calibration <- computeCensoringModelCalibration(clones, models)
    hazard <- calibration[calibration$model == "hazard" & calibration$grouping == "prediction", ]
    expect_equal(hazard$observed, as.vector(observed), tolerance = 1e-12)
    expect_equal(hazard$predicted, predicted(models[[2]]), tolerance = 1e-10)
  }
  expect_true(all(abs(z$day0) < 3), label = paste(round(z$day0, 2), collapse = ", "))
  # Without the interaction the hazard model assumes the same effect of X on every day, and is
  # miscalibrated:
  expect_gt(max(abs(z$none)), 5)
})

library(testthat)

# Tier 1 references for fitCensoringModels() (plan section 6.4): glm for the binary and pooled logistic
# models on hand-built risk sets, coxph(ties = "breslow") and survfit(stype = 2, ctype = 1) for the Cox
# initiation hazard, and glmnet at a fixed prior variance.

glmControl <- glm.control(epsilon = 1e-12, maxit = 100)

test_that("Binary model at G matches glm on the hand-built risk set", {
  for (riskStart in 0:1) {
    for (sameDayOrdering in c("eventsFirst", "deviationFirst")) {
      fixture <- createM3Fixture(150, riskStart = riskStart, seed = 2)
      data <- createM3Data(fixture)
      clones <- buildM3Clones(fixture$persons, m3Strategies(), sameDayOrdering)
      models <- suppressMessages(fitCensoringModels(clones, data, m3Args("cox")))
      model <- models[[1]]
      expect_equal(model$kind, "binary")
      expect_equal(model$horizon, 7L)

      riskSet <- refBinaryRiskSet(fixture$persons, 7, sameDayOrdering)
      expect_equal(sort(model$predictions$rowId), riskSet$rowId)
      x <- m3CovariateMatrix(fixture, riskSet$rowId)
      reference <- glm(riskSet$y ~ x, family = binomial(), control = glmControl)
      expect_equal(unname(model$coefficients), unname(coef(reference)), tolerance = 1e-6)
      expect_equal(names(model$coefficients), c("(Intercept)", "1001", "1002", "1003"))
      ps <- model$predictions$propensityScore[match(riskSet$rowId, model$predictions$rowId)]
      expect_equal(ps, unname(fitted(reference)), tolerance = 1e-6)
      expect_equal(model$numerator$propensityScore, mean(riskSet$y))
      expect_true(is.na(model$priorVariance))
      expect_equal(model$diagnostics$riskSetSize, nrow(riskSet))
      expect_equal(model$diagnostics$events, sum(riskSet$y))
      # AUC as the Mann-Whitney statistic:
      p1 <- ps[riskSet$y == 1]
      p0 <- ps[riskSet$y == 0]
      auc <- mean(outer(p1, p0, ">") + 0.5 * outer(p1, p0, "=="))
      expect_equal(model$diagnostics$auc, auc, tolerance = 1e-8)
      close(clones)
      close(data)
    }
  }
})

test_that("Cox initiation hazard matches coxph and survfit, per strategy horizon", {
  skip_if_not_installed("survival")
  for (riskStart in 0:1) {
    for (sameDayOrdering in c("eventsFirst", "deviationFirst")) {
      fixture <- createM3Fixture(150, riskStart = riskStart, seed = 3)
      data <- createM3Data(fixture)
      clones <- buildM3Clones(fixture$persons, m3Strategies("notInitiateWithin", "neverInitiate", graceEnd = 10), sameDayOrdering)
      models <- suppressMessages(fitCensoringModels(clones, data, m3Args("cox")))
      horizons <- c(10L, max(fixture$persons$riskEnd))
      for (strategyId in 1:2) {
        model <- models[[strategyId]]
        expect_equal(model$kind, "hazard")
        expect_equal(model$modelType, "cox")
        expect_equal(model$horizon, horizons[strategyId])

        rows <- refCoxRows(fixture$persons, horizons[strategyId], sameDayOrdering)
        expect_equal(sort(model$predictions$rowId), sort(rows$rowId))
        expect_equal(model$diagnostics$events, sum(rows$initiated))
        x <- m3CovariateMatrix(fixture, rows$rowId)
        reference <- survival::coxph(survival::Surv(rows$time, rows$initiated) ~ x, ties = "breslow", control = m3CoxControl())
        expect_equal(unname(model$coefficients), unname(coef(reference)), tolerance = 1e-6)
        expect_equal(names(model$coefficients), c("1001", "1002", "1003"))

        # S(d | X) = exp(-Lambda0(d + 1) exp(X beta)) vs survfit(newdata, stype = 2, ctype = 1) at t = d + 1:
        days <- 0:horizons[strategyId]
        newdata <- data.frame(x = I(x))
        expected <- refCoxSurvival(reference, newdata, days)
        linearPredictor <- model$predictions$linearPredictor[match(rows$rowId, model$predictions$rowId)]
        cumHazard <- c(0, model$baselineHazard$cumHazard)[findInterval(days + 1, model$baselineHazard$time) + 1]
        survival <- exp(-outer(exp(linearPredictor), cumHazard))
        expect_equal(survival, expected, tolerance = 1e-6)

        # Numerator: Breslow of the model without covariates (exp(-Breslow), not the product-limit KM):
        reference0 <- survival::coxph(survival::Surv(rows$time, rows$initiated) ~ 1, ties = "breslow")
        sf0 <- survival::survfit(reference0, stype = 2, ctype = 1)
        numeratorCumHazard <- c(0, model$numerator$baselineHazard$cumHazard)[findInterval(days + 1, model$numerator$baselineHazard$time) + 1]
        expect_equal(exp(-numeratorCumHazard), summary(sf0, times = days + 1, extend = TRUE)$surv, tolerance = 1e-10)
      }
      close(clones)
      close(data)
    }
  }
})

test_that("Weighted Breslow baseline hazard matches survfit on a weighted coxph", {
  skip_if_not_installed("survival")
  set.seed(5)
  n <- 120
  time <- sample(1:15, n, replace = TRUE)
  y <- rbinom(n, 1, 0.6)
  x <- rbinom(n, 1, 0.5)
  w <- runif(n, 0.2, 3)
  fit <- survival::coxph(survival::Surv(time, y) ~ x, weights = w, ties = "breslow", control = m3CoxControl())
  hazard <- .computeBreslowHazard(time, y, coef(fit) * x, w)
  sf <- survival::survfit(fit, newdata = data.frame(x = 0), stype = 2, ctype = 1)
  expected <- summary(sf, times = hazard$time)$cumhaz
  expect_equal(hazard$cumHazard, expected, tolerance = 1e-10)
  expect_equal(nrow(.computeBreslowHazard(time, rep(0, n))), 0)
})

test_that("Pooled logistic hazard matches glm on hand-built K_g rows", {
  configurations <- list(
    list(riskStart = 0L, sameDayOrdering = "eventsFirst", intervalLength = 1, types = c("notInitiateWithin", "notInitiateWithin")),
    list(riskStart = 1L, sameDayOrdering = "deviationFirst", intervalLength = 3, types = c("notInitiateWithin", "neverInitiate")),
    list(riskStart = 0L, sameDayOrdering = "deviationFirst", intervalLength = 1, types = c("neverInitiate", "initiateWithin"))
  )
  for (configuration in configurations) {
    fixture <- createM3Fixture(150, riskStart = configuration$riskStart, seed = 4)
    data <- createM3Data(fixture)
    strategies <- m3Strategies(configuration$types[1], configuration$types[2], graceEnd = 8)
    clones <- buildM3Clones(fixture$persons, strategies, configuration$sameDayOrdering)
    args <- m3Args("logistic", timeBasis = createTimeBasis(intervalLength = configuration$intervalLength))
    models <- suppressMessages(fitCensoringModels(clones, data, args))
    for (strategyId in which(configuration$types != "initiateWithin")) {
      model <- models[[strategyId]]
      horizon <- if (configuration$types[strategyId] == "neverInitiate") max(fixture$persons$riskEnd) else 8L
      expect_equal(model$horizon, horizon)
      nIntervals <- horizon %/% configuration$intervalLength + 1
      expect_equal(model$intervals$k, 0:(nIntervals - 1))
      expect_equal(max(model$intervals$endDay), horizon)

      rows <- refIntervalRows(fixture$persons, horizon, configuration$sameDayOrdering, configuration$intervalLength)
      # Intervals without initiations have hazard 0 and are left out of the fit:
      eventIntervals <- sort(unique(rows$k[rows$initiated == 1]))
      expect_equal(model$timeBasis$emptyIntervals, setdiff(0:(nIntervals - 1), eventIntervals))
      expect_true(all(model$intervals$baseline[model$timeBasis$emptyIntervals + 1] == -Inf))
      rows <- rows[rows$k %in% eventIntervals, ]
      expect_equal(model$diagnostics$rows, nrow(rows))
      x <- m3CovariateMatrix(fixture, rows$rowId)
      k <- factor(rows$k)
      reference <- glm(rows$initiated ~ k + x, family = quasibinomial(link = "logit"), control = glmControl)
      expect_equal(unname(model$coefficients), unname(coef(reference)), tolerance = 1e-6)
      expect_equal(names(model$coefficients), c("(Intercept)", paste0("timeBasis", eventIntervals[-1]), "1001", "1002", "1003"))
      # Linear predictor of every row from the stored baseline and per-person part:
      eta <- model$intervals$baseline[rows$k + 1] + model$predictions$linearPredictor[match(rows$rowId, model$predictions$rowId)]
      expect_equal(eta, unname(predict(reference, type = "link")), tolerance = 1e-6)

      reference0 <- glm(rows$initiated ~ k, family = quasibinomial(link = "logit"), control = glmControl)
      expect_equal(unname(model$numerator$coefficients), unname(coef(reference0)), tolerance = 1e-6)
      # The saturated numerator is the observed hazard of each interval:
      hazard <- tapply(rows$initiated, rows$k, mean)
      expect_equal(plogis(model$numerator$intervals$baseline[eventIntervals + 1]), unname(as.vector(hazard)), tolerance = 1e-6)
    }
    close(clones)
    close(data)
  }
})

test_that("Pooled logistic hazard with smooth time bases matches glm", {
  fixture <- createM3Fixture(150, riskStart = 0L, seed = 6)
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies("notInitiateWithin", "neverInitiate", graceEnd = 12), "eventsFirst")
  timeBases <- list(
    createTimeBasis("naturalSpline", intervalLength = 1, df = 3),
    createTimeBasis("polynomial", intervalLength = 2, degree = 2, raw = TRUE),
    createTimeBasis("polynomial", intervalLength = 1, degree = 3, raw = FALSE),
    createTimeBasis("linear", intervalLength = 1)
  )
  for (timeBasis in timeBases) {
    models <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic", timeBasis = timeBasis)))
    for (strategyId in 1:2) {
      model <- models[[strategyId]]
      horizon <- model$horizon
      nIntervals <- horizon %/% timeBasis$intervalLength + 1
      k <- 0:(nIntervals - 1)
      expectedBasis <- switch(timeBasis$type,
        naturalSpline = splines::ns(k, df = timeBasis$df),
        polynomial = stats::poly(k, degree = timeBasis$degree, raw = timeBasis$raw),
        linear = matrix(k, ncol = 1)
      )
      expect_equal(unname(model$timeBasis$basis), unname(matrix(as.numeric(expectedBasis), nrow = nIntervals)))
      expect_length(model$timeBasis$emptyIntervals, 0)

      rows <- refIntervalRows(fixture$persons, horizon, "eventsFirst", timeBasis$intervalLength)
      basis <- model$timeBasis$basis[rows$k + 1, , drop = FALSE]
      x <- m3CovariateMatrix(fixture, rows$rowId)
      reference <- glm(rows$initiated ~ basis + x, family = quasibinomial(link = "logit"), control = glmControl)
      expect_equal(unname(model$coefficients), unname(coef(reference)), tolerance = 1e-6)
      reference0 <- glm(rows$initiated ~ basis, family = quasibinomial(link = "logit"), control = glmControl)
      expect_equal(unname(model$numerator$coefficients), unname(coef(reference0)), tolerance = 1e-6)
      expect_equal(model$numerator$intervals$baseline, unname(as.vector(cbind(1, model$timeBasis$basis) %*% coef(reference0))), tolerance = 1e-6)
    }
  }
  close(clones)
  close(data)
})

test_that("A single interval (notInitiateWithin(0)) gives an intercept-only time model", {
  fixture <- createM3Fixture(150, riskStart = 1L, seed = 7)
  fixture$persons$treatmentDay[1:40] <- 0L
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = 0), "eventsFirst")
  models <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic")))
  model <- models[[2]]
  expect_equal(nrow(model$intervals), 1)
  expect_equal(ncol(model$timeBasis$basis), 0)
  rows <- refIntervalRows(fixture$persons, 0, "eventsFirst", 1)
  x <- m3CovariateMatrix(fixture, rows$rowId)
  reference <- glm(rows$initiated ~ x, family = binomial(), control = glmControl)
  expect_equal(unname(model$coefficients), unname(coef(reference)), tolerance = 1e-6)
  expect_equal(unname(model$numerator$coefficients), qlogis(mean(rows$initiated)), tolerance = 1e-10)
  # At K = {0} the binary model at G = 0 is the complementary model on the same persons:
  expect_equal(sort(models[[1]]$predictions$rowId), sort(rows$rowId))
  ps <- models[[1]]$predictions$propensityScore[match(rows$rowId, models[[1]]$predictions$rowId)]
  expect_equal(ps, unname(fitted(reference)), tolerance = 1e-6)

  # With riskStart = 1 all follow-up is after G = 0, so each clone has one segment with W = 1 / pi in
  # the initiate arm and W = 1 / (1 - pi) in the defer arm:
  weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = FALSE, truncation = createWeightTruncation("none"))))
  segments <- collectSegments(weighted)
  expect_false(any(duplicated(segments$cloneRowId)))
  pi <- unname(fitted(reference))[match(segments$rowId, rows$rowId)]
  expected <- ifelse(segments$strategyId == 1, 1 / pi, 1 / (1 - pi))
  expect_equal(segments$weight, expected, tolerance = 1e-6)
  close(clones)
  close(data)
})

test_that("Covariates can be included or excluded, and highly correlated ones are caught", {
  fixture <- createM3Fixture(150, riskStart = 0L, seed = 12)
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies("notInitiateWithin", "neverInitiate", graceEnd = 7), "eventsFirst")
  args <- m3Args("cox")
  args$excludeCovariateIds <- 1003
  models <- suppressMessages(fitCensoringModels(clones, data, args))
  rows <- refCoxRows(fixture$persons, 7, "eventsFirst")
  x <- m3CovariateMatrix(fixture, rows$rowId)[, 1:2]
  if (requireNamespace("survival", quietly = TRUE)) {
    reference <- survival::coxph(survival::Surv(rows$time, rows$initiated) ~ x, ties = "breslow", control = m3CoxControl())
    expect_equal(unname(models[[1]]$coefficients), unname(coef(reference)), tolerance = 1e-6)
  }
  expect_equal(names(models[[1]]$coefficients), c("1001", "1002"))

  args <- m3Args("logistic", timeBasis = createTimeBasis("linear"))
  args$includeCovariateIds <- 1002
  models <- suppressMessages(fitCensoringModels(clones, data, args))
  expect_equal(names(models[[1]]$coefficients), c("(Intercept)", "timeBasis1", "1002"))
  close(data)

  # A covariate that is the initiation indicator:
  initiated <- fixture$persons$rowId[!is.na(fixture$persons$treatmentDay) & fixture$persons$treatmentDay <= 7]
  fixture$covariates <- rbind(fixture$covariates[fixture$covariates$covariateId != 1001, ], data.frame(rowId = initiated, covariateId = 1001, covariateValue = 1))
  data <- createM3Data(fixture)
  args <- m3Args("cox")
  args$errorOnHighCorrelation <- TRUE
  expect_error(suppressMessages(fitCensoringModels(clones, data, args)), "High correlation")
  args$stopOnError <- FALSE
  # (The neverInitiate model passes the check but separates on this covariate; Cyclops warns.)
  models <- suppressWarnings(suppressMessages(fitCensoringModels(clones, data, args)))
  expect_match(models[[1]]$diagnostics$status, "High correlation")
  expect_length(models[[1]]$coefficients, 0)
  expect_equal(models[[1]]$baselineHazard, models[[1]]$numerator$baselineHazard)
  close(clones)
  close(data)
})

test_that("useCovariates = FALSE fits intercept and time only", {
  fixture <- createM3Fixture(120, riskStart = 0L, seed = 8)
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies(), "eventsFirst")
  models <- suppressMessages(fitCensoringModels(clones, data, m3Args("cox", useCovariates = FALSE)))
  riskSet <- refBinaryRiskSet(fixture$persons, 7, "eventsFirst")
  expect_equal(models[[1]]$coefficients, c("(Intercept)" = qlogis(mean(riskSet$y))))
  expect_true(all(models[[1]]$predictions$propensityScore == mean(riskSet$y)))
  expect_length(models[[2]]$coefficients, 0)
  expect_true(all(models[[2]]$predictions$linearPredictor == 0))
  expect_equal(models[[2]]$baselineHazard, models[[2]]$numerator$baselineHazard)

  modelsLogistic <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic", useCovariates = FALSE)))
  expect_equal(modelsLogistic[[2]]$coefficients, modelsLogistic[[2]]$numerator$coefficients)
  expect_equal(modelsLogistic[[2]]$intervals, modelsLogistic[[2]]$numerator$intervals)
  close(clones)
  close(data)
})

test_that("Fits are shared only between strategies with the same type, treatment and horizon", {
  fixture <- createM3Fixture(100, riskStart = 0L, seed = 9)
  data <- createM3Data(fixture)
  maxRiskEnd <- max(fixture$persons$riskEnd)
  same <- list(
    createTreatmentStrategy(1, "Defer A", "notInitiateWithin", graceEnd = 5),
    createTreatmentStrategy(2, "Defer B", "notInitiateWithin", graceEnd = 5)
  )
  clones <- buildM3Clones(fixture$persons, same, "eventsFirst")
  messages <- capture_messages(models <- fitCensoringModels(clones, data, m3Args("cox", useCovariates = FALSE)))
  expect_true(any(grepl("Strategy 2 shares the censoring model", messages)))
  expect_equal(models[[1]][-1], models[[2]][-1])
  close(clones)

  different <- list(
    createTreatmentStrategy(1, "Defer", "notInitiateWithin", graceEnd = maxRiskEnd),
    createTreatmentStrategy(2, "Never", "neverInitiate")
  )
  clones <- buildM3Clones(fixture$persons, different, "eventsFirst")
  messages <- capture_messages(models <- fitCensoringModels(clones, data, m3Args("cox", useCovariates = FALSE)))
  expect_false(any(grepl("shares", messages)))
  expect_equal(models[[1]]$horizon, models[[2]]$horizon)
  close(clones)
  close(data)
})

test_that("Prior variance: fixed variances skip cross-validation, cross-validation stores its choice", {
  fixture <- createM3Fixture(150, riskStart = 0L, seed = 10)
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies(), "eventsFirst")
  args <- m3Args("logistic", prior = Cyclops::createPrior("laplace", exclude = c(0), useCrossValidation = TRUE))
  args$control <- Cyclops::createControl(noiseLevel = "silent", cvType = "auto", fold = 5, cvRepetitions = 1, seed = 1, tolerance = 1e-8)
  models <- suppressMessages(fitCensoringModels(clones, data, args))
  for (model in models) {
    expect_true(is.finite(model$priorVariance) && model$priorVariance > 0)
  }
  fixed <- suppressMessages(fitCensoringModels(clones, data, args, .fixedPriorVariances = list("1" = 0.37, "2" = 0.21)))
  expect_equal(fixed[[1]]$priorVariance, 0.37)
  expect_equal(fixed[[2]]$priorVariance, 0.21)
  # Replaying the cross-validated variances reproduces the fit:
  replay <- suppressMessages(fitCensoringModels(clones, data, args, .fixedPriorVariances = list("1" = models[[1]]$priorVariance, "2" = models[[2]]$priorVariance)))
  expect_equal(replay[[1]]$coefficients, models[[1]]$coefficients, tolerance = 1e-5)
  expect_equal(replay[[2]]$coefficients, models[[2]]$coefficients, tolerance = 1e-5)

  expect_error(fitCensoringModels(clones, data, args, .fixedPriorVariances = list("3" = 1)), "names")
  expect_error(fitCensoringModels(data, data, args), "CcwClones")
  close(clones)
  close(data)
})

test_that("Laplace and normal priors at a fixed variance match glmnet", {
  skip_if_not_installed("glmnet", minimum_version = "5.0")
  skip_if_not_installed("survival")
  expect_equal(eval(formals(glmnet::glmnet)$cox.ties)[1], "breslow")
  fixture <- createM3Fixture(200, riskStart = 0L, seed = 11)
  data <- createM3Data(fixture)
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = 6), "eventsFirst")
  variance <- 0.5
  # glmnet minimizes -loglik / N + lambda * penalty; Cyclops -loglik + sqrt(2 / v) |b| (Laplace) or
  # b^2 / (2 v) (normal). glmnet rescales penalty factors to sum to the number of columns.
  glmnetCoefficients <- function(x, y, family, priorType, penaltyFactor = rep(1, ncol(x))) {
    nobs <- if (family == "cox") nrow(y) else length(y)
    target <- if (priorType == "laplace") sqrt(2 / variance) / nobs else 1 / (variance * nobs)
    target <- target * sum(penaltyFactor) / length(penaltyFactor)
    path <- target * exp(seq(log(20), 0, length.out = 30))
    fit <- glmnet::glmnet(
      x = x, y = y, family = family, alpha = if (priorType == "laplace") 1 else 0,
      lambda = path, standardize = FALSE, penalty.factor = penaltyFactor, cox.ties = "breslow",
      control = list(thresh = 1e-14, maxit = 1e6)
    )
    return(as.vector(as.matrix(stats::coef(fit, s = target, exact = FALSE))))
  }
  for (priorType in c("laplace", "normal")) {
    for (hazardModelType in c("logistic", "cox")) {
      prior <- Cyclops::createPrior(priorType, exclude = c(0), useCrossValidation = TRUE)
      args <- m3Args(hazardModelType, prior = prior, timeBasis = createTimeBasis("linear"))
      args$control <- Cyclops::createControl(noiseLevel = "silent", tolerance = 1e-12, convergenceType = "gradient", maxIterations = 10000)
      models <- suppressMessages(fitCensoringModels(clones, data, args, .fixedPriorVariances = list("1" = variance, "2" = variance)))

      # Binary model through CohortMethod::createPs:
      riskSet <- refBinaryRiskSet(fixture$persons, 6, "eventsFirst")
      x <- m3CovariateMatrix(fixture, riskSet$rowId)
      expected <- glmnetCoefficients(x, riskSet$y, "binomial", priorType)
      expect_equal(unname(models[[1]]$coefficients), expected, tolerance = 1e-5)
      expect_equal(models[[1]]$priorVariance, variance)

      if (hazardModelType == "logistic") {
        rows <- refIntervalRows(fixture$persons, 6, "eventsFirst", 1)
        x <- cbind(k = rows$k, m3CovariateMatrix(fixture, rows$rowId))
        expected <- glmnetCoefficients(x, rows$initiated, "binomial", priorType, penaltyFactor = c(0, 1, 1, 1))
      } else {
        rows <- refCoxRows(fixture$persons, 6, "eventsFirst")
        x <- m3CovariateMatrix(fixture, rows$rowId)
        expected <- glmnetCoefficients(x, survival::Surv(rows$time, rows$initiated), "cox", priorType)
      }
      expect_equal(unname(models[[2]]$coefficients), expected, tolerance = 1e-5)
      expect_equal(models[[2]]$priorVariance, variance)
      if (priorType == "laplace") {
        expect_equal(abs(unname(models[[2]]$coefficients)) > 1e-6, abs(expected) > 1e-6)
      }
    }
  }
  close(clones)
  close(data)
})

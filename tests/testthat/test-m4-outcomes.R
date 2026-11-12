# Tier 1 reference tests of the outcome models (Phase1Contracts section 5.5; plan sections 2.6, 2.7 and
# 6.4). Fixtures run with riskStart = 0 and 1; horizons h are read at t_h = h - riskStart + 1.

library(testthat)

coxControl <- function() survival::coxph.control(eps = 1e-12, toler.chol = 1e-14, iter.max = 100)

# Cox ---------------------------------------------------------------------------------------------

test_that("Cox MSM matches coxph: log HR, robust and model-based variance", {
  skip_if_not_installed("survival")
  for (riskStart in 0:1) {
    for (competingEvents in c("none", "death")) {
      fixture <- createM4Fixture(riskStart = riskStart, competingEvents = competingEvents)
      data <- toSurvivalData(fixture$segments)
      reference <- survival::coxph(
        survival::Surv(startTime, stopTime, outcome) ~ x,
        data = data, weights = weight, cluster = personSeqId, ties = "breslow", control = coxControl()
      )
      robust <- fitM4(fixture, modelType = "cox", horizons = 10)
      modelBased <- fitM4(fixture, modelType = "cox", horizons = 10, varianceEstimator = "modelBased")

      expect_equal(getM4Estimate(robust, "hr", "logEstimate"), unname(coef(reference)), tolerance = 1e-6)
      expect_equal(getM4Estimate(robust, "hr", "seLogEstimate"), sqrt(vcov(reference)[1, 1]), tolerance = 1e-6)
      expect_equal(getM4Estimate(modelBased, "hr", "seLogEstimate"), sqrt(reference$naive.var[1, 1]), tolerance = 1e-6)
      expect_equal(getM4Estimate(robust, "hr"), exp(unname(coef(reference))), tolerance = 1e-6)
      expect_equal(getM4Estimate(robust, "hr", "lb95"), exp(coef(reference) - qnorm(0.975) * sqrt(vcov(reference)[1, 1]))[[1]], tolerance = 1e-6)
      expect_equal(getM4Estimate(robust, "hr", "p"), 2 * pnorm(-abs(coef(reference) / sqrt(vcov(reference)[1, 1])))[[1]], tolerance = 1e-6)

      if (competingEvents == "death") {
        # The two cause-specific Cox MSMs and their joint robust covariance equal the multi-state coxph.
        multiState <- survival::coxph(
          survival::Surv(startTime, stopTime, state) ~ x,
          data = data, weights = weight, id = cloneRowId, cluster = personSeqId, ties = "breslow", control = coxControl()
        )
        expect_equal(unname(robust$coefficients), unname(coef(multiState)), tolerance = 1e-6)
        expect_equal(unname(robust$vcov), unname(vcov(multiState)), tolerance = 1e-6)
        expect_equal(unname(diag(modelBased$vcov)), unname(diag(multiState$naive.var)), tolerance = 1e-6)
      } else {
        expect_equal(names(robust$coefficients), "target")
      }
    }
  }
})

test_that("Cox MSM risks without competing events match survfit.coxph (ctype = 1, stype = 2)", {
  skip_if_not_installed("survival")
  horizonsByRiskStart <- list(c(0, 4, 12, 30), c(1, 5, 13, 31))
  for (riskStart in 0:1) {
    horizons <- horizonsByRiskStart[[riskStart + 1]]
    fixture <- createM4Fixture(riskStart = riskStart, competingEvents = "none")
    data <- toSurvivalData(fixture$segments)
    model <- fitM4(fixture, modelType = "cox", horizons = horizons)
    # Evaluate the reference curves at our coefficient ("given identical coefficients"; the
    # coefficient itself is checked against coxph above):
    reference <- survival::coxph(
      survival::Surv(startTime, stopTime, outcome) ~ x,
      data = data, weights = weight, cluster = personSeqId, ties = "breslow",
      init = model$coefficients[["target"]], control = survival::coxph.control(iter.max = 0)
    )
    curves <- survival::survfit(reference, newdata = data.frame(x = 1:0), ctype = 1, stype = 2)
    summaryAtHorizons <- summary(curves, times = horizons - riskStart + 1, extend = TRUE)
    for (i in seq_along(horizons)) {
      risks <- unname(1 - summaryAtHorizons$surv[i, ])
      expect_equal(getM4Estimate(model, sprintf("risk_target_%s", horizons[i])), risks[1], tolerance = 1e-8)
      expect_equal(getM4Estimate(model, sprintf("risk_comparator_%s", horizons[i])), risks[2], tolerance = 1e-8)
      expect_equal(getM4Estimate(model, sprintf("rd_%s", horizons[i])), risks[1] - risks[2], tolerance = 1e-8)
      expect_equal(getM4Estimate(model, sprintf("rr_%s", horizons[i])), risks[1] / risks[2], tolerance = 1e-8)
    }
  }
})

test_that("Cox MSM risks with competing death match the multi-state coxph by Aalen-Johansen", {
  skip_if_not_installed("survival")
  for (riskStart in 0:1) {
    horizons <- riskStart + c(0, 5, 15, 40)
    fixture <- createM4Fixture(riskStart = riskStart, competingEvents = "death")
    data <- toSurvivalData(fixture$segments)
    model <- fitM4(fixture, modelType = "cox", horizons = horizons)
    multiState <- survival::coxph(
      survival::Surv(startTime, stopTime, state) ~ x,
      data = data, weights = weight, id = cloneRowId, cluster = personSeqId, ties = "breslow",
      init = unname(model$coefficients), control = survival::coxph.control(iter.max = 0)
    )
    # stype = 1 gives the Aalen-Johansen product of Phase1Contracts section 5.5; the survfit default
    # (stype = 2) uses the matrix exponential of the hazard increments instead.
    curves <- survival::survfit(multiState, newdata = data.frame(x = 1:0), stype = 1)
    summaryAtHorizons <- summary(curves, times = horizons - riskStart + 1, extend = TRUE)
    for (i in seq_along(horizons)) {
      risks <- unname(summaryAtHorizons$pstate[i, , 2])
      expect_equal(getM4Estimate(model, sprintf("risk_target_%s", horizons[i])), risks[1], tolerance = 1e-8)
      expect_equal(getM4Estimate(model, sprintf("risk_comparator_%s", horizons[i])), risks[2], tolerance = 1e-8)
    }
  }
})

test_that("Cox MSM risk standard errors equal the numerical infinitesimal jackknife", {
  skip_if_not_installed("survival")
  for (competingEvents in c("none", "death")) {
    riskStart <- if (competingEvents == "none") 0 else 1
    horizons <- riskStart + c(4, 12)
    times <- horizons - riskStart + 1
    fixture <- createM4Fixture(n = 50, riskStart = riskStart, competingEvents = competingEvents, sharedPersons = 6, seed = 2)
    model <- fitM4(fixture, modelType = "cox", horizons = horizons)
    statistic <- function(segments) {
      data <- toSurvivalData(segments)
      control <- survival::coxph.control(eps = 1e-14, toler.chol = 1e-15, iter.max = 100)
      if (competingEvents == "none") {
        fit <- survival::coxph(survival::Surv(startTime, stopTime, outcome) ~ x, data = data, weights = weight, ties = "breslow", control = control)
        curves <- survival::survfit(fit, newdata = data.frame(x = 1:0), ctype = 1, stype = 2)
        risks <- 1 - summary(curves, times = times, extend = TRUE)$surv
      } else {
        fit <- survival::coxph(survival::Surv(startTime, stopTime, state) ~ x, data = data, weights = weight, id = cloneRowId, ties = "breslow", control = control)
        curves <- survival::survfit(fit, newdata = data.frame(x = 1:0), stype = 1)
        risks <- summary(curves, times = times, extend = TRUE)$pstate[, , 2]
      }
      return(as.vector(t(risks)))
    }
    risks <- statistic(fixture$segments)
    influence <- numericalPersonInfluence(fixture$segments, statistic, epsilon = 1e-4)
    for (i in seq_along(horizons)) {
      target <- influence[, 2 * i - 1]
      comparator <- influence[, 2 * i]
      expect_equal(getM4Estimate(model, sprintf("risk_target_%s", horizons[i]), "se"), sqrt(sum(target^2)), tolerance = 1e-6)
      expect_equal(getM4Estimate(model, sprintf("risk_comparator_%s", horizons[i]), "se"), sqrt(sum(comparator^2)), tolerance = 1e-6)
      expect_equal(getM4Estimate(model, sprintf("rd_%s", horizons[i]), "se"), sqrt(sum((target - comparator)^2)), tolerance = 1e-6)
      logRrInfluence <- target / risks[2 * i - 1] - comparator / risks[2 * i]
      expect_equal(getM4Estimate(model, sprintf("rr_%s", horizons[i]), "seLogEstimate"), sqrt(sum(logRrInfluence^2)), tolerance = 1e-6)
    }
  }
})

# Kaplan-Meier and Aalen-Johansen --------------------------------------------------------------------

test_that("Weighted Kaplan-Meier matches survfit: risks, IJ SEs, CIs, RD, RR and RMST", {
  skip_if_not_installed("survival")
  for (riskStart in 0:1) {
    horizons <- riskStart + c(0, 3, 9, 20, 40)
    times <- horizons - riskStart + 1
    fixture <- createM4Fixture(riskStart = riskStart, competingEvents = "none")
    data <- toSurvivalData(fixture$segments)
    model <- fitM4(fixture, modelType = "km", horizons = horizons)
    allPersons <- sort(unique(data$personSeqId))
    fits <- list()
    for (strategyId in 1:2) {
      strategyData <- data[data$strategyId == strategyId, ]
      fit <- survival::survfit(
        survival::Surv(startTime, stopTime, outcome) ~ 1,
        data = strategyData, weights = weight, id = cloneRowId, cluster = personSeqId, robust = TRUE, influence = TRUE, model = TRUE
      )
      riskInfluence <- -residuals(fit, times = times, type = "pstate", collapse = FALSE, weighted = TRUE)
      rmstInfluence <- sapply(times, function(t) residuals(fit, times = t, type = "rmst", collapse = FALSE, weighted = TRUE))
      fits[[strategyId]] <- list(
        summary = summary(fit, times = times, extend = TRUE),
        rmst = sapply(times, function(t) summary(fit, rmean = t)$table[["rmean"]]),
        riskInfluence = apply(riskInfluence, 2, personInfluence, personSeqId = strategyData$personSeqId, allPersons = allPersons),
        rmstInfluence = apply(rmstInfluence, 2, personInfluence, personSeqId = strategyData$personSeqId, allPersons = allPersons)
      )
    }
    for (i in seq_along(horizons)) {
      h <- horizons[i]
      for (strategyId in 1:2) {
        role <- c("target", "comparator")[strategyId]
        s <- fits[[strategyId]]$summary
        expect_equal(getM4Estimate(model, sprintf("risk_%s_%s", role, h)), 1 - s$surv[i], tolerance = 1e-10)
        expect_equal(getM4Estimate(model, sprintf("risk_%s_%s", role, h), "se"), unname(s$std.err[i]), tolerance = 1e-8)
        expect_equal(getM4Estimate(model, sprintf("risk_%s_%s", role, h), "lb95"), 1 - unname(s$upper[i]), tolerance = 1e-8)
        expect_equal(getM4Estimate(model, sprintf("risk_%s_%s", role, h), "ub95"), 1 - unname(s$lower[i]), tolerance = 1e-8)
        expect_equal(getM4Estimate(model, sprintf("rmst_%s_%s", role, h)), fits[[strategyId]]$rmst[i], tolerance = 1e-10)
        expect_equal(getM4Estimate(model, sprintf("rmst_%s_%s", role, h), "se"), sqrt(sum(fits[[strategyId]]$rmstInfluence[, i]^2)), tolerance = 1e-8)
        expect_equal(getM4Estimate(model, sprintf("rmtl_%s_%s", role, h)), times[i] - fits[[strategyId]]$rmst[i], tolerance = 1e-10)
      }
      risks <- 1 - c(fits[[1]]$summary$surv[i], fits[[2]]$summary$surv[i])
      target <- fits[[1]]$riskInfluence[, i]
      comparator <- fits[[2]]$riskInfluence[, i]
      expect_equal(getM4Estimate(model, sprintf("rd_%s", h)), risks[1] - risks[2], tolerance = 1e-10)
      expect_equal(getM4Estimate(model, sprintf("rd_%s", h), "se"), sqrt(sum((target - comparator)^2)), tolerance = 1e-8)
      if (all(risks > 0)) {
        expect_equal(getM4Estimate(model, sprintf("rr_%s", h)), risks[1] / risks[2], tolerance = 1e-10)
        expect_equal(
          getM4Estimate(model, sprintf("rr_%s", h), "seLogEstimate"),
          sqrt(sum((target / risks[1] - comparator / risks[2])^2)),
          tolerance = 1e-8
        )
      } else {
        expect_true(is.na(getM4Estimate(model, sprintf("rr_%s", h))))
      }
      rmstDifference <- fits[[1]]$rmst[i] - fits[[2]]$rmst[i]
      expect_equal(getM4Estimate(model, sprintf("rmst_%s", h)), rmstDifference, tolerance = 1e-10)
      expect_equal(
        getM4Estimate(model, sprintf("rmst_%s", h), "se"),
        sqrt(sum((fits[[1]]$rmstInfluence[, i] - fits[[2]]$rmstInfluence[, i])^2)),
        tolerance = 1e-8
      )
      expect_equal(getM4Estimate(model, sprintf("rmtl_%s", h)), -rmstDifference, tolerance = 1e-10)
    }
  }
})

test_that("Weighted Aalen-Johansen matches the multi-state survfit: CIF, IJ SEs, CIs, RD, RR and RMTL", {
  skip_if_not_installed("survival")
  for (riskStart in 0:1) {
    horizons <- riskStart + c(0, 3, 9, 20, 40)
    times <- horizons - riskStart + 1
    fixture <- createM4Fixture(riskStart = riskStart, competingEvents = "death")
    data <- toSurvivalData(fixture$segments)
    model <- fitM4(fixture, modelType = "km", horizons = horizons)
    expect_false(any(grepl("^rmst", model$estimates$estimand)))
    allPersons <- sort(unique(data$personSeqId))
    fits <- list()
    for (strategyId in 1:2) {
      strategyData <- data[data$strategyId == strategyId, ]
      fit <- survival::survfit(
        survival::Surv(startTime, stopTime, state) ~ 1,
        data = strategyData, weights = weight, id = cloneRowId, cluster = personSeqId, robust = TRUE, influence = TRUE, model = TRUE
      )
      riskInfluence <- residuals(fit, times = times, type = "pstate", collapse = FALSE, weighted = TRUE)[, "outcome", ]
      aucInfluence <- sapply(times, function(t) residuals(fit, times = t, type = "auc", collapse = FALSE, weighted = TRUE)[, "outcome"])
      fits[[strategyId]] <- list(
        summary = summary(fit, times = times, extend = TRUE),
        rmtl = sapply(times, function(t) summary(fit, rmean = t)$table["outcome", "rmean"]),
        riskInfluence = apply(riskInfluence, 2, personInfluence, personSeqId = strategyData$personSeqId, allPersons = allPersons),
        aucInfluence = apply(aucInfluence, 2, personInfluence, personSeqId = strategyData$personSeqId, allPersons = allPersons)
      )
    }
    for (i in seq_along(horizons)) {
      h <- horizons[i]
      for (strategyId in 1:2) {
        role <- c("target", "comparator")[strategyId]
        s <- fits[[strategyId]]$summary
        expect_equal(getM4Estimate(model, sprintf("risk_%s_%s", role, h)), unname(s$pstate[i, 2]), tolerance = 1e-10)
        expect_equal(getM4Estimate(model, sprintf("risk_%s_%s", role, h), "se"), unname(s$std.err[i, 2]), tolerance = 1e-8)
        expect_equal(getM4Estimate(model, sprintf("risk_%s_%s", role, h), "lb95"), unname(s$lower[i, 2]), tolerance = 1e-8)
        expect_equal(getM4Estimate(model, sprintf("risk_%s_%s", role, h), "ub95"), unname(s$upper[i, 2]), tolerance = 1e-8)
        expect_equal(getM4Estimate(model, sprintf("rmtl_%s_%s", role, h)), fits[[strategyId]]$rmtl[i], tolerance = 1e-10)
        expect_equal(getM4Estimate(model, sprintf("rmtl_%s_%s", role, h), "se"), sqrt(sum(fits[[strategyId]]$aucInfluence[, i]^2)), tolerance = 1e-8)
      }
      risks <- unname(c(fits[[1]]$summary$pstate[i, 2], fits[[2]]$summary$pstate[i, 2]))
      target <- fits[[1]]$riskInfluence[, i]
      comparator <- fits[[2]]$riskInfluence[, i]
      expect_equal(getM4Estimate(model, sprintf("rd_%s", h)), risks[1] - risks[2], tolerance = 1e-10)
      expect_equal(getM4Estimate(model, sprintf("rd_%s", h), "se"), sqrt(sum((target - comparator)^2)), tolerance = 1e-8)
      if (all(risks > 0)) {
        expect_equal(
          getM4Estimate(model, sprintf("rr_%s", h), "seLogEstimate"),
          sqrt(sum((target / risks[1] - comparator / risks[2])^2)),
          tolerance = 1e-8
        )
      }
      expect_equal(getM4Estimate(model, sprintf("rmtl_%s", h)), fits[[1]]$rmtl[i] - fits[[2]]$rmtl[i], tolerance = 1e-10)
      expect_equal(
        getM4Estimate(model, sprintf("rmtl_%s", h), "se"),
        sqrt(sum((fits[[1]]$aucInfluence[, i] - fits[[2]]$aucInfluence[, i])^2)),
        tolerance = 1e-8
      )
    }
  }
})

test_that("computeCumulativeIncidence matches survfit curves", {
  skip_if_not_installed("survival")
  for (competingEvents in c("none", "death")) {
    for (riskStart in 0:1) {
      fixture <- createM4Fixture(riskStart = riskStart, competingEvents = competingEvents)
      data <- toSurvivalData(fixture$segments)
      weightedClones <- createM4WeightedClones(fixture)
      curves <- computeCumulativeIncidence(weightedClones)
      horizons <- riskStart + c(2, 10, 33)
      atHorizons <- computeCumulativeIncidence(weightedClones, horizons = horizons)
      close(weightedClones)
      expect_equal(attr(curves, "strategyLabels"), c("Initiate within 4 days", "Do not initiate within 4 days"))
      for (strategyId in 1:2) {
        strategyData <- data[data$strategyId == strategyId, ]
        if (competingEvents == "none") {
          fit <- survival::survfit(
            survival::Surv(startTime, stopTime, outcome) ~ 1,
            data = strategyData, weights = weight, id = cloneRowId, cluster = personSeqId, robust = TRUE
          )
        } else {
          fit <- survival::survfit(
            survival::Surv(startTime, stopTime, state) ~ 1,
            data = strategyData, weights = weight, id = cloneRowId, cluster = personSeqId, robust = TRUE
          )
        }
        curve <- curves[curves$strategyId == strategyId, ]
        expect_equal(curve$day[1], riskStart - 1)
        expect_equal(curve$risk[1], 0)
        expect_equal(curve$atRiskCount[1], sum(fixture$clones$strategyId == strategyId))
        curve <- curve[-1, ]
        eventTimes <- sort(unique(strategyData$stopTime[strategyData$outcome == 1 | strategyData$competing == 1]))
        expect_equal(curve$day, eventTimes + riskStart - 1)
        for (rows in list(list(curve = curve, times = eventTimes), list(curve = atHorizons[atHorizons$strategyId == strategyId, ], times = horizons - riskStart + 1))) {
          s <- summary(fit, times = rows$times, extend = TRUE)
          if (competingEvents == "none") {
            expect_equal(rows$curve$risk, 1 - s$surv, tolerance = 1e-10)
            expect_equal(rows$curve$lb95, 1 - s$upper, tolerance = 1e-8)
            expect_equal(rows$curve$ub95, 1 - s$lower, tolerance = 1e-8)
            expect_equal(rows$curve$atRiskWeighted, s$n.risk, tolerance = 1e-10)
            expect_equal(rows$curve$events, s$n.event, tolerance = 1e-10)
          } else {
            expect_equal(rows$curve$risk, unname(s$pstate[, 2]), tolerance = 1e-10)
            expect_equal(rows$curve$lb95, unname(s$lower[, 2]), tolerance = 1e-8)
            expect_equal(rows$curve$ub95, unname(s$upper[, 2]), tolerance = 1e-8)
            expect_equal(rows$curve$atRiskWeighted, unname(s$n.risk[, 1]), tolerance = 1e-10)
            expect_equal(rows$curve$events, unname(s$n.event[, 2]), tolerance = 1e-10)
          }
          atRiskCount <- sapply(rows$times, function(t) sum(strategyData$startTime < t & strategyData$stopTime >= t))
          expect_equal(rows$curve$atRiskCount, atRiskCount)
        }
      }
    }
  }
})

# Pooled logistic ---------------------------------------------------------------------------------

# The time basis built independently of the package from the 0-based interval index.
referenceBasis <- function(type, k, kNew, df = 3, degree = 2, raw = TRUE) {
  if (type == "indicator") {
    levels <- sort(unique(k))[-1]
    build <- function(x) sapply(levels, function(level) as.numeric(x == level))
  } else if (type == "linear") {
    build <- function(x) cbind(as.numeric(x))
  } else if (type == "polynomial") {
    basis <- stats::poly(k, degree = degree, raw = raw)
    build <- function(x) unclass(predict(basis, x))
  } else {
    basis <- splines::ns(k, df = df)
    build <- function(x) unclass(predict(basis, x))
  }
  return(list(columns = matrix(build(k), nrow = length(k)), build = function(x) matrix(build(x), nrow = length(x))))
}

test_that("Pooled logistic MSM matches glm and vcovCL (HC0, clustered on person) without competing events", {
  skip_if_not_installed("sandwich")
  settings <- list(
    # Indicator bases on fixtures with outcomes in every strategy and interval:
    list(type = "indicator", intervalLength = 5, interaction = FALSE, riskStart = 0, maxDays = 14, horizons = c(4, 9, 14)),
    list(type = "indicator", intervalLength = 7, interaction = TRUE, riskStart = 1, maxDays = 20, horizons = c(6, 13, 20)),
    list(type = "linear", intervalLength = 2, interaction = TRUE, riskStart = 0),
    list(type = "polynomial", intervalLength = 3, interaction = FALSE, riskStart = 1, raw = TRUE, horizons = c(5, 11, 20)),
    list(type = "polynomial", intervalLength = 1, interaction = TRUE, riskStart = 0, raw = FALSE),
    list(type = "naturalSpline", intervalLength = 2, interaction = FALSE, riskStart = 1)
  )
  for (setting in settings) {
    raw <- if (is.null(setting$raw)) TRUE else setting$raw
    riskStart <- setting$riskStart
    intervalLength <- setting$intervalLength
    horizons <- riskStart + if (is.null(setting$horizons)) c(5, 13, 21) else setting$horizons
    if (is.null(setting$maxDays)) {
      fixture <- createM4Fixture(n = 200, riskStart = riskStart, competingEvents = "none", seed = 3)
    } else {
      fixture <- createM4Fixture(n = 200, riskStart = riskStart, competingEvents = "none", seed = 3, maxDays = setting$maxDays, outcomeRate = 0.7)
    }
    model <- fitM4(
      fixture,
      modelType = "logistic",
      horizons = horizons,
      timeBasis = createTimeBasis(type = setting$type, intervalLength = intervalLength, df = 3, degree = 2, raw = raw),
      armByTimeInteraction = setting$interaction
    )
    rows <- expandM4Rows(fixture$segments, intervalLength)
    basis <- referenceBasis(setting$type, rows$k, df = 3, degree = 2, raw = raw)
    design <- cbind(1, rows$x, basis$columns)
    if (setting$interaction) {
      design <- cbind(design, rows$x * basis$columns)
    }
    referenceFit <- glm(rows$y ~ design - 1, family = quasibinomial(), weights = rows$weight, control = glm.control(epsilon = 1e-14, maxit = 100))
    expect_equal(unname(model$coefficients), unname(coef(referenceFit)), tolerance = 1e-7)
    vcovCl <- sandwich::vcovCL(referenceFit, cluster = rows$personSeqId, type = "HC0", cadjust = FALSE)
    expect_equal(unname(model$vcov), unname(vcovCl), tolerance = 1e-6)
    expect_equal(getM4Estimate(model, "or", "logEstimate"), unname(coef(referenceFit)[2]), tolerance = 1e-7)
    expect_equal(getM4Estimate(model, "or", "seLogEstimate"), sqrt(vcovCl[2, 2]), tolerance = 1e-6)

    # Risks from our coefficients through the independently built design, with the analytic gradient
    # (1 - risk) sum_k h_k X_k for the delta method.
    for (h in horizons) {
      nIntervals <- floor((h - riskStart + 1) / intervalLength)
      k <- seq_len(nIntervals) - 1
      newBasis <- basis$build(k)
      risks <- numeric(2)
      gradients <- matrix(0, 2, length(model$coefficients))
      for (x in 1:0) {
        newDesign <- cbind(1, x, newBasis)
        if (setting$interaction) {
          newDesign <- cbind(newDesign, x * newBasis)
        }
        hazard <- plogis(as.vector(newDesign %*% model$coefficients))
        risk <- 1 - prod(1 - hazard)
        risks[2 - x] <- risk
        gradients[2 - x, ] <- (1 - risk) * colSums(hazard * newDesign)
      }
      expect_equal(getM4Estimate(model, sprintf("risk_target_%s", h)), risks[1], tolerance = 1e-8)
      expect_equal(getM4Estimate(model, sprintf("risk_comparator_%s", h)), risks[2], tolerance = 1e-8)
      se <- sqrt(c(
        t(gradients[1, ]) %*% vcovCl %*% gradients[1, ],
        t(gradients[2, ]) %*% vcovCl %*% gradients[2, ],
        t(gradients[1, ] - gradients[2, ]) %*% vcovCl %*% (gradients[1, ] - gradients[2, ])
      ))
      logGradient <- gradients[1, ] / risks[1] - gradients[2, ] / risks[2]
      expect_equal(getM4Estimate(model, sprintf("risk_target_%s", h), "se"), se[1], tolerance = 1e-6)
      expect_equal(getM4Estimate(model, sprintf("risk_comparator_%s", h), "se"), se[2], tolerance = 1e-6)
      expect_equal(getM4Estimate(model, sprintf("rd_%s", h), "se"), se[3], tolerance = 1e-6)
      expect_equal(getM4Estimate(model, sprintf("rr_%s", h), "seLogEstimate"), sqrt(as.numeric(t(logGradient) %*% vcovCl %*% logGradient)), tolerance = 1e-6)
    }
  }
})

test_that("Young et al. nonparametric CIF equals the survfit Aalen-Johansen at intervalLength = 1", {
  skip_if_not_installed("survival")
  for (riskStart in 0:1) {
    fixture <- createM4Fixture(riskStart = riskStart, competingEvents = "death")
    data <- toSurvivalData(fixture$segments)
    rows <- expandM4Rows(fixture$segments, 1)
    times <- c(1, 6, 15, 22)
    for (strategyId in 1:2) {
      fit <- survival::survfit(
        survival::Surv(startTime, stopTime, state) ~ 1,
        data = data[data$strategyId == strategyId, ], weights = weight, id = cloneRowId
      )
      cif <- youngNonparametricCif(rows[rows$x == as.integer(strategyId == 1), ], max(times))
      expect_equal(cif[times], summary(fit, times = times)$pstate[, 2], tolerance = 1e-10)
    }
  }
})

test_that("Pooled logistic MSM with competing death follows Young et al. (2020)", {
  skip_if_not_installed("sandwich")
  # Saturated per-strategy hazards (indicator basis with interaction) reproduce the nonparametric CIF.
  fixture <- createM4Fixture(n = 200, riskStart = 1, competingEvents = "death", seed = 3, maxDays = 14, outcomeRate = 0.6, deathRate = 0.4)
  intervalLength <- 5
  rows <- expandM4Rows(fixture$segments, intervalLength)
  cells <- aggregate(cbind(y, death) ~ x + k, data = rows, FUN = sum)
  expect_true(all(cells$y > 0) && all(cells$death > 0))
  horizons <- 1 + c(4, 9, 14)
  model <- fitM4(
    fixture,
    modelType = "logistic", horizons = horizons, armByTimeInteraction = TRUE,
    timeBasis = createTimeBasis(type = "indicator", intervalLength = intervalLength)
  )
  for (strategyId in 1:2) {
    cif <- youngNonparametricCif(rows[rows$x == as.integer(strategyId == 1), ], 3)
    role <- c("target", "comparator")[strategyId]
    expect_equal(sapply(horizons, function(h) getM4Estimate(model, sprintf("risk_%s_%s", role, h))), cif, tolerance = 1e-6)
  }

  # A parametric basis: our CIF equals the Young CIF of the glm hazards, and the joint robust covariance
  # of the two models equals vcovCL of the stacked glm.
  for (riskStart in 0:1) {
    fixture <- createM4Fixture(n = 200, riskStart = riskStart, competingEvents = "death", seed = 5)
    intervalLength <- 2
    horizons <- riskStart + c(3, 11, 19)
    model <- fitM4(
      fixture,
      modelType = "logistic", horizons = horizons, armByTimeInteraction = TRUE,
      timeBasis = createTimeBasis(type = "polynomial", intervalLength = intervalLength, degree = 2)
    )
    rows <- expandM4Rows(fixture$segments, intervalLength)
    design <- cbind(1, rows$x, rows$k, rows$k^2, rows$x * rows$k, rows$x * rows$k^2)
    noDeath <- rows$death == 0
    control <- glm.control(epsilon = 1e-14, maxit = 100)
    outcomeFit <- glm(rows$y[noDeath] ~ design[noDeath, ] - 1, family = quasibinomial(), weights = rows$weight[noDeath], control = control)
    deathFit <- glm(rows$death ~ design - 1, family = quasibinomial(), weights = rows$weight, control = control)
    expect_equal(unname(model$coefficients), unname(c(coef(outcomeFit), coef(deathFit))), tolerance = 1e-7)
    expect_equal(names(model$coefficients)[7:8], c("death:(Intercept)", "death:target"))

    # Stacked glm: the outcome rows (no death) and the death rows, each with its own coefficients.
    stackedDesign <- rbind(
      cbind(design[noDeath, ], 0 * design[noDeath, ]),
      cbind(0 * design, design)
    )
    stackedY <- c(rows$y[noDeath], rows$death)
    stackedWeight <- c(rows$weight[noDeath], rows$weight)
    stackedCluster <- c(rows$personSeqId[noDeath], rows$personSeqId)
    stackedFit <- glm(stackedY ~ stackedDesign - 1, family = quasibinomial(), weights = stackedWeight, control = control)
    vcovCl <- sandwich::vcovCL(stackedFit, cluster = stackedCluster, type = "HC0", cadjust = FALSE)
    expect_equal(unname(model$vcov), unname(vcovCl), tolerance = 1e-6)

    riskFunction <- function(theta) {
      result <- c()
      for (h in horizons) {
        nIntervals <- floor((h - riskStart + 1) / intervalLength)
        k <- seq_len(nIntervals) - 1
        for (x in 1:0) {
          newDesign <- cbind(1, x, k, k^2, x * k, x * k^2)
          hazardY <- plogis(as.vector(newDesign %*% theta[1:6]))
          hazardD <- plogis(as.vector(newDesign %*% theta[7:12]))
          result <- c(result, utils::tail(youngCif(hazardY, hazardD), 1))
        }
      }
      return(result)
    }
    risks <- riskFunction(model$coefficients)
    # Reference Jacobian by Richardson-extrapolated central differences:
    jacobian <- sapply(seq_along(model$coefficients), function(j) {
      difference <- function(step) {
        up <- model$coefficients
        down <- model$coefficients
        up[j] <- up[j] + step
        down[j] <- down[j] - step
        return((riskFunction(up) - riskFunction(down)) / (2 * step))
      }
      return((4 * difference(1e-4) - difference(2e-4)) / 3)
    })
    for (i in seq_along(horizons)) {
      h <- horizons[i]
      expect_equal(getM4Estimate(model, sprintf("risk_target_%s", h)), risks[2 * i - 1], tolerance = 1e-8)
      expect_equal(getM4Estimate(model, sprintf("risk_comparator_%s", h)), risks[2 * i], tolerance = 1e-8)
      gradientTarget <- jacobian[2 * i - 1, ]
      gradientComparator <- jacobian[2 * i, ]
      expect_equal(getM4Estimate(model, sprintf("risk_target_%s", h), "se"), sqrt(as.numeric(t(gradientTarget) %*% vcovCl %*% gradientTarget)), tolerance = 1e-5)
      gradientRd <- gradientTarget - gradientComparator
      expect_equal(getM4Estimate(model, sprintf("rd_%s", h), "se"), sqrt(as.numeric(t(gradientRd) %*% vcovCl %*% gradientRd)), tolerance = 1e-5)
      # CIs on the cumulative incidence, log scale:
      se <- getM4Estimate(model, sprintf("risk_target_%s", h), "se")
      expect_equal(getM4Estimate(model, sprintf("risk_target_%s", h), "lb95"), risks[2 * i - 1] * exp(-qnorm(0.975) * se / risks[2 * i - 1]), tolerance = 1e-8)
    }
  }
})

test_that("Pooled logistic risks beyond the follow-up are NA with a warning", {
  fixture <- createM4Fixture(n = 60, riskStart = 0, competingEvents = "none", maxDays = 12)
  expect_warning(
    model <- fitM4(fixture, modelType = "logistic", horizons = c(5, 40), timeBasis = createTimeBasis("linear", intervalLength = 2)),
    "No follow-up"
  )
  expect_false(is.na(getM4Estimate(model, "risk_target_5")))
  expect_true(is.na(getM4Estimate(model, "risk_target_40")))
})

test_that("Pooled logistic risks are read at the end of the last complete interval by the horizon", {
  for (riskStart in 0:1) {
    fixture <- createM4Fixture(n = 200, riskStart = riskStart, competingEvents = "none", seed = 3, maxDays = 40, outcomeRate = 0.7)
    # t_h = 5 (no complete interval), 10 and 15 (one), 20 and 29 (two):
    horizons <- riskStart - 1 + c(5, 10, 15, 20, 29)
    warnings <- capture_warnings(
      model <- fitM4(fixture, modelType = "logistic", horizons = horizons, timeBasis = createTimeBasis("linear", intervalLength = 10))
    )
    expect_match(warnings, sprintf("No complete interval of 10 days by horizon\\(s\\) %d:", horizons[1]), all = FALSE)
    expect_match(warnings, sprintf("Horizon\\(s\\) %d, %d do not end .* \\(day\\(s\\) %d, %d\\)", horizons[3], horizons[5], horizons[2], horizons[4]), all = FALSE)
    expect_length(warnings, 2)
    estimates <- model$estimates
    for (estimand in sprintf(c("risk_target_%d", "risk_comparator_%d", "rd_%d", "rr_%d"), horizons[1])) {
      expect_true(all(is.na(unlist(estimates[estimates$estimand == estimand, -1]))))
    }
    for (pair in list(horizons[2:3], horizons[4:5])) {
      for (format in c("risk_target_%d", "risk_comparator_%d", "rd_%d", "rr_%d")) {
        early <- estimates[estimates$estimand == sprintf(format, pair[1]), -1]
        expect_false(is.na(early$estimate))
        expect_equal(early, estimates[estimates$estimand == sprintf(format, pair[2]), -1])
      }
    }
    expect_false(isTRUE(all.equal(getM4Estimate(model, sprintf("risk_target_%d", horizons[2])), getM4Estimate(model, sprintf("risk_target_%d", horizons[4])))))
  }
})

test_that("Pooled logistic hazards of intervals one strategy never enters are not identified with arm x time", {
  skip_if_not_installed("sandwich")
  # Target clones end by t = 10 and comparator clones run to t = 15: with an indicator basis and the
  # interaction, the target's hazard in interval 2 (t in (10, 15]) is not identified.
  fixture <- createM4Fixture(n = 200, riskStart = 0, competingEvents = "none", seed = 3, maxDays = 14, outcomeRate = 0.7)
  drop <- fixture$clones$cloneRowId[fixture$clones$strategyId == 1 & fixture$clones$survivalTime > 10]
  fixture$segments <- fixture$segments[!fixture$segments$cloneRowId %in% drop, ]
  fixture$clones <- fixture$clones[!fixture$clones$cloneRowId %in% drop, ]
  oldIds <- fixture$clones$cloneRowId
  fixture$clones$cloneRowId <- seq_along(oldIds)
  fixture$segments$cloneRowId <- match(fixture$segments$cloneRowId, oldIds)
  fixture$segments$segmentId <- seq_len(nrow(fixture$segments))
  fixture$personWeights <- m4PersonWeights(fixture$segments)
  expect_equal(max(fixture$segments$stopTime[fixture$segments$strategyId == 1]), 10)
  expect_equal(max(fixture$segments$stopTime[fixture$segments$strategyId == 2]), 15)
  horizons <- c(4, 9, 14)
  timeBasis <- createTimeBasis("indicator", intervalLength = 5)

  warnings <- capture_warnings(
    model <- fitM4(fixture, modelType = "logistic", horizons = horizons, timeBasis = timeBasis, armByTimeInteraction = TRUE)
  )
  expect_match(warnings, "target:time2 of the outcome pooled logistic model are not identified", all = FALSE)
  expect_match(warnings, "target strategy are not identified in some intervals: its risks at horizon\\(s\\) 14 are NA", all = FALSE)
  expect_true(is.na(model$coefficients[["target:time2"]]))
  expect_true(all(is.na(model$vcov["target:time2", ])))
  for (estimand in c("risk_target_14", "rd_14", "rr_14")) {
    expect_true(is.na(getM4Estimate(model, estimand)))
  }
  expect_true(is.na(getM4Estimate(model, "risk_target_14", "se")))

  # The identified coefficients and their robust covariance equal glm and vcovCL without the column,
  # and the identified risks follow from them.
  rows <- expandM4Rows(fixture$segments, 5)
  design <- cbind(1, rows$x, rows$k == 1, rows$k == 2, rows$x * (rows$k == 1))
  referenceFit <- glm(rows$y ~ design - 1, family = quasibinomial(), weights = rows$weight, control = glm.control(epsilon = 1e-14, maxit = 100))
  vcovCl <- sandwich::vcovCL(referenceFit, cluster = rows$personSeqId, type = "HC0", cadjust = FALSE)
  identified <- c("(Intercept)", "target", "time1", "time2", "target:time1")
  expect_equal(unname(model$coefficients[identified]), unname(coef(referenceFit)), tolerance = 1e-7)
  expect_equal(unname(model$vcov[identified, identified]), unname(vcovCl), tolerance = 1e-6)
  beta <- coef(referenceFit)
  hazardTarget <- plogis(beta[1] + beta[2] + c(0, beta[3] + beta[5]))
  hazardComparator <- plogis(beta[1] + c(0, beta[3], beta[4]))
  expect_equal(getM4Estimate(model, "risk_target_9"), 1 - prod(1 - hazardTarget), tolerance = 1e-7)
  expect_equal(getM4Estimate(model, "risk_comparator_14"), 1 - prod(1 - hazardComparator), tolerance = 1e-7)
  expect_false(is.na(getM4Estimate(model, "risk_comparator_14", "se")))
  expect_false(is.na(getM4Estimate(model, "rd_9", "se")))

  # The point estimates for the bootstrap driver are the same:
  weightedClones <- createM4WeightedClones(fixture)
  on.exit(close(weightedClones))
  args <- createFitOutcomeModelArgs(modelType = "logistic", horizons = horizons, timeBasis = timeBasis, armByTimeInteraction = TRUE, control = m4Control)
  point <- suppressWarnings(.fitOutcomeModelPoint(weightedClones, args))
  expect_equal(point$estimates, model$estimates[, c("estimand", "estimate", "logEstimate")])

  # Without the interaction, the proportional-odds model identifies every hazard:
  expect_silent(withoutInteraction <- fitM4(fixture, modelType = "logistic", horizons = horizons, timeBasis = timeBasis))
  expect_false(anyNA(withoutInteraction$estimates$estimate))
})

test_that("A singular information matrix gives NA variance with a warning", {
  expect_warning(inverse <- .invertInformation(matrix(1, 2, 2), "test"), "singular")
  expect_equal(dim(inverse), c(2, 2))
  expect_true(all(is.na(inverse)))
  expect_equal(.invertInformation(diag(c(2, 4)), "test"), diag(c(0.5, 0.25)))
})

test_that("Pooled logistic model-based variance equals the inverse weighted information of glm", {
  for (riskStart in 0:1) {
    for (competingEvents in c("none", "death")) {
      fixture <- createM4Fixture(n = 200, riskStart = riskStart, competingEvents = competingEvents, seed = 5)
      intervalLength <- 2
      horizons <- riskStart + c(5, 13)
      model <- fitM4(
        fixture,
        modelType = "logistic", horizons = horizons, varianceEstimator = "modelBased", armByTimeInteraction = TRUE,
        timeBasis = createTimeBasis("linear", intervalLength = intervalLength)
      )
      rows <- expandM4Rows(fixture$segments, intervalLength)
      design <- cbind(1, rows$x, rows$k, rows$x * rows$k)
      control <- glm.control(epsilon = 1e-14, maxit = 100)
      outcomeRows <- if (competingEvents == "death") rows$death == 0 else rep(TRUE, nrow(rows))
      outcomeFit <- glm(rows$y[outcomeRows] ~ design[outcomeRows, ] - 1, family = quasibinomial(), weights = rows$weight[outcomeRows], control = control)
      # cov.unscaled is the binomial covariance (dispersion 1): the inverse weighted information.
      expected <- summary(outcomeFit)$cov.unscaled
      coefficients <- coef(outcomeFit)
      if (competingEvents == "death") {
        deathFit <- glm(rows$death ~ design - 1, family = quasibinomial(), weights = rows$weight, control = control)
        expected <- rbind(cbind(expected, 0 * expected), cbind(0 * expected, summary(deathFit)$cov.unscaled))
        coefficients <- c(coefficients, coef(deathFit))
      }
      expect_equal(unname(model$coefficients), unname(coefficients), tolerance = 1e-7)
      expect_equal(unname(model$vcov), unname(expected), tolerance = 1e-6)
      expect_equal(getM4Estimate(model, "or", "seLogEstimate"), sqrt(expected[2, 2]), tolerance = 1e-6)

      # Delta-method risk SEs with that covariance (reference Jacobian by Richardson extrapolation):
      riskFunction <- function(theta) {
        result <- c()
        for (h in horizons) {
          k <- seq_len((h - riskStart + 1) / intervalLength) - 1
          for (x in 1:0) {
            newDesign <- cbind(1, x, k, x * k)
            hazardY <- plogis(as.vector(newDesign %*% theta[1:4]))
            hazardD <- if (competingEvents == "death") plogis(as.vector(newDesign %*% theta[5:8])) else 0 * hazardY
            result <- c(result, utils::tail(youngCif(hazardY, hazardD), 1))
          }
        }
        return(result)
      }
      jacobian <- sapply(seq_along(coefficients), function(j) {
        difference <- function(step) {
          up <- coefficients
          down <- coefficients
          up[j] <- up[j] + step
          down[j] <- down[j] - step
          return((riskFunction(up) - riskFunction(down)) / (2 * step))
        }
        return((4 * difference(1e-4) - difference(2e-4)) / 3)
      })
      risks <- riskFunction(coefficients)
      for (i in seq_along(horizons)) {
        h <- horizons[i]
        gradientTarget <- jacobian[2 * i - 1, ]
        gradientComparator <- jacobian[2 * i, ]
        expect_equal(getM4Estimate(model, sprintf("risk_target_%s", h)), risks[2 * i - 1], tolerance = 1e-7)
        expect_equal(getM4Estimate(model, sprintf("risk_target_%s", h), "se"), sqrt(as.numeric(t(gradientTarget) %*% expected %*% gradientTarget)), tolerance = 1e-5)
        expect_equal(getM4Estimate(model, sprintf("risk_comparator_%s", h), "se"), sqrt(as.numeric(t(gradientComparator) %*% expected %*% gradientComparator)), tolerance = 1e-5)
        gradientRd <- gradientTarget - gradientComparator
        expect_equal(getM4Estimate(model, sprintf("rd_%s", h), "se"), sqrt(as.numeric(t(gradientRd) %*% expected %*% gradientRd)), tolerance = 1e-5)
      }
    }
  }
})

test_that("The pooled logistic sandwich does not depend on the person block size", {
  fixture <- createM4Fixture(n = 200, riskStart = 1, competingEvents = "death", seed = 5)
  weightedClones <- createM4WeightedClones(fixture)
  on.exit(close(weightedClones))
  data <- .getOutcomeData(weightedClones)
  rows <- .expandToIntervals(data$segments, 2)
  basis <- .timeBasisColumns(createTimeBasis("polynomial", intervalLength = 2, degree = 2), rows$interval)
  nIntervals <- max(rows$interval) + 1
  design <- .pooledDesign(rep(0:1, each = nIntervals), basis$build(rep(seq_len(nIntervals) - 1, 2)), basis$names, TRUE)
  rows$cell <- as.integer(rows$x * nIntervals + rows$interval + 1)
  noDeath <- rows$death == 0
  fits <- list(
    .fitPooledModel(rows[noDeath, ], rows$y[noDeath], design, Cyclops::createPrior("none"), m4Control, "outcome", "robust"),
    .fitPooledModel(rows, rows$death, design, Cyclops::createPrior("none"), m4Control, "death", "robust")
  )
  nPersons <- length(data$personSeqIds)
  oneBlock <- .jointSandwich(fits, nPersons)
  expect_equal(.jointSandwich(fits, nPersons, maxEntries = 100), oneBlock, tolerance = 1e-12)
  expect_equal(.jointSandwich(fits, nPersons, maxEntries = 1), oneBlock, tolerance = 1e-12)
})

# Properties ------------------------------------------------------------------------------------------

test_that("Robust SEs are unchanged when all weights are multiplied by a constant", {
  fixture <- createM4Fixture(riskStart = 1, competingEvents = "death")
  scaled <- fixture
  scaled$segments$weight <- scaled$segments$weight * 3.7
  for (modelType in c("cox", "logistic", "km")) {
    timeBasis <- createTimeBasis("linear", intervalLength = 3)
    model <- fitM4(fixture, modelType = modelType, horizons = c(6, 15), timeBasis = timeBasis)
    scaledModel <- fitM4(scaled, modelType = modelType, horizons = c(6, 15), timeBasis = timeBasis)
    expect_equal(scaledModel$estimates$estimate, model$estimates$estimate, tolerance = 1e-6)
    expect_equal(scaledModel$estimates$se, model$estimates$se, tolerance = 1e-6)
    expect_equal(scaledModel$estimates$seLogEstimate, model$estimates$seLogEstimate, tolerance = 1e-6)
  }
  # The model-based variance does change:
  modelBased <- fitM4(fixture, modelType = "cox", horizons = 6, varianceEstimator = "modelBased")
  scaledModelBased <- fitM4(scaled, modelType = "cox", horizons = 6, varianceEstimator = "modelBased")
  expect_equal(scaledModelBased$vcov, modelBased$vcov / 3.7, tolerance = 1e-6)
})

test_that("Horizons are read at t = horizon - riskStart + 1", {
  # The same data on the survival time scale with riskStart = 0 and riskStart = 1:
  fixture0 <- createM4Fixture(riskStart = 0, competingEvents = "death")
  fixture1 <- fixture0
  fixture1$persons[, c("riskStart", "riskEnd", "nominalRiskEnd", "outcomeDay", "deathDay")] <-
    fixture1$persons[, c("riskStart", "riskEnd", "nominalRiskEnd", "outcomeDay", "deathDay")] + 1L
  fixture1$clones$riskStart <- 1L
  fixture1$clones$endDay <- fixture1$clones$endDay + 1L
  for (modelType in c("cox", "logistic", "km")) {
    timeBasis <- createTimeBasis("polynomial", intervalLength = 2)
    model0 <- fitM4(fixture0, modelType = modelType, horizons = c(3, 9, 13), timeBasis = timeBasis)
    model1 <- fitM4(fixture1, modelType = modelType, horizons = c(4, 10, 14), timeBasis = timeBasis)
    expect_equal(model1$estimates$estimate, model0$estimates$estimate, tolerance = 1e-10)
    expect_equal(model1$estimates$se, model0$estimates$se, tolerance = 1e-10)
    expect_equal(model1$metaData$riskStart, 1L)
  }
  weightedClones <- createM4WeightedClones(fixture1)
  on.exit(close(weightedClones))
  expect_error(suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(modelType = "km", horizons = 0))), "riskStart")
})

test_that("Splitting a segment without changing its weight leaves all estimates unchanged", {
  fixture <- createM4Fixture(riskStart = 0, competingEvents = "death", extraSplits = FALSE)
  split <- fixture
  segments <- split$segments
  long <- segments$stopTime - segments$startTime >= 2
  middle <- floor((segments$startTime[long] + segments$stopTime[long]) / 2)
  first <- segments[long, ]
  second <- segments[long, ]
  first$stopTime <- middle
  first$outcome <- 0L
  first$competing <- 0L
  second$startTime <- middle
  segments <- rbind(segments[!long, ], first, second)
  segments <- segments[order(segments$cloneRowId, segments$startTime), ]
  segments$segmentId <- seq_len(nrow(segments))
  split$segments <- segments
  for (modelType in c("cox", "logistic", "km")) {
    timeBasis <- createTimeBasis("naturalSpline", intervalLength = 2, df = 3)
    model <- fitM4(fixture, modelType = modelType, horizons = c(5, 11), timeBasis = timeBasis)
    splitModel <- fitM4(split, modelType = modelType, horizons = c(5, 11), timeBasis = timeBasis)
    expect_equal(splitModel$estimates, model$estimates, tolerance = 1e-8)
  }
})

test_that("Frequency weights in .fitOutcomeModelPoint equal copying persons", {
  # riskStart = 1, so horizons 6 and 15 end pooled logistic intervals of 3 days (t = 6 and 15).
  fixture <- createM4Fixture(n = 80, riskStart = 1, competingEvents = "death", seed = 6)
  persons <- sort(unique(fixture$segments$personSeqId))
  set.seed(7)
  personCounts <- data.frame(personSeqId = persons, count = sample(0:3, length(persons), replace = TRUE))
  # Copy-based data: each person copied `count` times with new person and clone IDs.
  copies <- list()
  for (copy in 1:3) {
    copiedPersons <- personCounts$personSeqId[personCounts$count >= copy]
    segments <- fixture$segments[fixture$segments$personSeqId %in% copiedPersons, ]
    segments$personSeqId <- segments$personSeqId + copy * 10000L
    segments$cloneRowId <- segments$cloneRowId + copy * 10000L
    copies[[copy]] <- segments
  }
  copySegments <- do.call(rbind, copies)
  copyClones <- fixture$clones[match(copySegments$cloneRowId %% 10000L, fixture$clones$cloneRowId), ]
  copyClones$cloneRowId <- copySegments$cloneRowId
  copyClones$personSeqId <- copySegments$personSeqId
  copyClones$rowId <- copySegments$rowId + (copySegments$cloneRowId %/% 10000L) * 10000L
  copyClones <- copyClones[!duplicated(copyClones$cloneRowId), ]
  # Renumber clones 1..n as the contract requires:
  newCloneRowId <- match(copySegments$cloneRowId, copyClones$cloneRowId)
  copyClones$cloneRowId <- seq_len(nrow(copyClones))
  copySegments$cloneRowId <- newCloneRowId
  copySegments$segmentId <- seq_len(nrow(copySegments))
  copyPersons <- fixture$persons[match(copyClones$rowId %% 10000L, fixture$persons$rowId), ]
  copyPersons$rowId <- copyClones$rowId
  copyPersons$personSeqId <- copyClones$personSeqId
  copyPersons <- copyPersons[!duplicated(copyPersons$rowId), ]
  copyFixture <- fixture
  copyFixture$persons <- copyPersons
  copyFixture$clones <- copyClones
  copyFixture$segments <- copySegments
  copyFixture$personWeights <- m4PersonWeights(copySegments)

  weightedClones <- createM4WeightedClones(fixture)
  copyWeightedClones <- createM4WeightedClones(copyFixture)
  on.exit({
    close(weightedClones)
    close(copyWeightedClones)
  })
  outcomeData <- .getOutcomeData(weightedClones)
  settings <- list(
    list(modelType = "cox"),
    list(modelType = "km"),
    list(modelType = "logistic", type = "linear", interaction = FALSE),
    # Bases whose definition depends on the distribution of the rows' intervals (knots at quantiles;
    # orthogonal polynomials), with and without the arm x time interaction:
    list(modelType = "logistic", type = "naturalSpline", interaction = FALSE),
    list(modelType = "logistic", type = "naturalSpline", interaction = TRUE),
    list(modelType = "logistic", type = "polynomial", interaction = FALSE),
    list(modelType = "logistic", type = "polynomial", interaction = TRUE)
  )
  for (setting in settings) {
    type <- if (is.null(setting$type)) "linear" else setting$type
    args <- createFitOutcomeModelArgs(
      modelType = setting$modelType, horizons = c(6, 15), varianceEstimator = "bootstrap",
      timeBasis = createTimeBasis(type, intervalLength = 3, df = 3, degree = 2, raw = FALSE),
      armByTimeInteraction = isTRUE(setting$interaction), control = m4Control
    )
    weighted <- .fitOutcomeModelPoint(weightedClones, args, personCounts = personCounts)
    copied <- .fitOutcomeModelPoint(copyWeightedClones, args)
    expect_false(anyNA(weighted$estimates$estimate[!grepl("^rr", weighted$estimates$estimand)]))
    expect_equal(weighted$estimates, copied$estimates, tolerance = 1e-7)
    expect_equal(weighted$coefficients, copied$coefficients, tolerance = 1e-7)
    expect_equal(colnames(weighted$estimates), c("estimand", "estimate", "logEstimate"))
    # The same from segments collected once:
    fromData <- .fitOutcomeModelPoint(NULL, args, personCounts = personCounts, outcomeData = outcomeData)
    expect_equal(fromData, weighted)
  }
  expect_error(.fitOutcomeModelPoint(weightedClones, args, personCounts = data.frame(personSeqId = 1, count = 0.5)), "count")
})

test_that(".fitOutcomeModelPoint equals the point estimates of fitOutcomeModel", {
  fixture <- createM4Fixture(riskStart = 0, competingEvents = "none")
  weightedClones <- createM4WeightedClones(fixture)
  on.exit(close(weightedClones))
  for (modelType in c("cox", "logistic", "km")) {
    args <- createFitOutcomeModelArgs(modelType = modelType, horizons = c(5, 14), timeBasis = createTimeBasis("linear", intervalLength = 3), control = m4Control)
    model <- suppressMessages(fitOutcomeModel(weightedClones, args))
    point <- .fitOutcomeModelPoint(weightedClones, args)
    expect_equal(point$estimates, model$estimates[, c("estimand", "estimate", "logEstimate")])
    expect_equal(point$coefficients, model$coefficients)
  }
})

test_that("Zero-weight segments are dropped", {
  fixture <- createM4Fixture(riskStart = 0, competingEvents = "death")
  zeroed <- fixture
  zeroClones <- fixture$clones$cloneRowId[fixture$clones$rowId %in% c(5, 17, 33)]
  zeroed$segments$weight[zeroed$segments$cloneRowId %in% zeroClones] <- 0
  for (modelType in c("cox", "logistic", "km")) {
    timeBasis <- createTimeBasis("linear", intervalLength = 3)
    model <- fitM4(zeroed, modelType = modelType, horizons = c(5, 11), timeBasis = timeBasis)
    # The same data without those clones: zero-weight persons contribute nothing to the robust SEs.
    removed <- fixture
    removed$segments <- removed$segments[!removed$segments$cloneRowId %in% zeroClones, ]
    removed$clones <- removed$clones[!removed$clones$cloneRowId %in% zeroClones, ]
    oldIds <- removed$clones$cloneRowId
    removed$clones$cloneRowId <- seq_along(oldIds)
    removed$segments$cloneRowId <- match(removed$segments$cloneRowId, oldIds)
    removed$personWeights <- m4PersonWeights(removed$segments)
    removedModel <- fitM4(removed, modelType = modelType, horizons = c(5, 11), timeBasis = timeBasis)
    expect_equal(model$estimates, removedModel$estimates, tolerance = 1e-8)
  }
})

# Interface ---------------------------------------------------------------------------------------

test_that("fitOutcomeModel returns the contract's object and counts", {
  fixture <- createM4Fixture(riskStart = 1, competingEvents = "death")
  model <- fitM4(fixture, modelType = "km", horizons = c(10, 20))
  expect_s3_class(model, "CcwOutcomeModel")
  expect_named(model, c("modelType", "varianceEstimator", "estimates", "coefficients", "vcov", "counts", "metaData"))
  expect_equal(colnames(model$estimates), c("estimand", "estimate", "logEstimate", "seLogEstimate", "se", "lb95", "ub95", "p"))
  expect_equal(
    model$estimates$estimand,
    c(
      "risk_target_10", "risk_comparator_10", "rd_10", "rr_10", "rmtl_target_10", "rmtl_comparator_10", "rmtl_10",
      "risk_target_20", "risk_comparator_20", "rd_20", "rr_20", "rmtl_target_20", "rmtl_comparator_20", "rmtl_20"
    )
  )
  expect_null(model$coefficients)
  expect_equal(model$metaData$riskStart, 1L)
  expect_equal(model$metaData$competingEvents, "death")
  expect_equal(model$metaData$outcomeId, 3)
  expect_equal(model$metaData$args$modelType, "km")

  counts <- model$counts
  clones <- fixture$clones
  segments <- fixture$segments
  for (strategyId in 1:2) {
    inStrategy <- clones$strategyId == strategyId
    segmentsInStrategy <- segments$strategyId == strategyId
    row <- counts[counts$strategyId == strategyId, ]
    expect_equal(row$clones, sum(inStrategy))
    expect_equal(row$persons, length(unique(clones$personSeqId[inStrategy])))
    expect_equal(row$outcomes, sum(clones$outcomeCount[inStrategy]))
    expect_equal(row$competingEvents, sum(clones$endReason[inStrategy] == "death"))
    expect_equal(row$personDays, sum(clones$survivalTime[inStrategy]))
    expect_equal(row$weightedPersonDays, sum((segments$weight * (segments$stopTime - segments$startTime))[segmentsInStrategy]))
    expect_equal(row$weightedOutcomes, sum((segments$weight * segments$outcome)[segmentsInStrategy]))
    expect_equal(row$weightedCompetingEvents, sum((segments$weight * segments$competing)[segmentsInStrategy]))
  }

  cox <- fitM4(fixture, modelType = "cox", horizons = 10)
  expect_equal(cox$estimates$estimand, c("hr", "risk_target_10", "risk_comparator_10", "rd_10", "rr_10"))
  expect_equal(names(cox$coefficients), c("target", "death:target"))
  expect_equal(dimnames(cox$vcov), list(names(cox$coefficients), names(cox$coefficients)))
  logistic <- fitM4(fixture, modelType = "logistic", horizons = 10, timeBasis = createTimeBasis("linear", intervalLength = 5), armByTimeInteraction = TRUE)
  expect_equal(logistic$estimates$estimand[1], "or")
  expect_equal(names(logistic$coefficients)[1:4], c("(Intercept)", "target", "time", "target:time"))
})

test_that("Cox risks have no SEs with the model-based variance", {
  fixture <- createM4Fixture(riskStart = 0, competingEvents = "none")
  model <- fitM4(fixture, modelType = "cox", horizons = 10, varianceEstimator = "modelBased")
  expect_false(is.na(getM4Estimate(model, "hr", "seLogEstimate")))
  expect_true(is.na(getM4Estimate(model, "risk_target_10", "se")))
  expect_false(is.na(getM4Estimate(model, "risk_target_10")))
})

test_that("fitOutcomeModel rejects the bootstrap and unsupported settings", {
  fixture <- createM4Fixture(riskStart = 0, competingEvents = "none")
  weightedClones <- createM4WeightedClones(fixture)
  on.exit(close(weightedClones))
  expect_error(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(varianceEstimator = "bootstrap")), "Use the bootstrap driver")
  expect_error(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(modelType = "km", varianceEstimator = "modelBased")), "not available")
  crossValidated <- createFitOutcomeModelArgs(modelType = "logistic", prior = Cyclops::createPrior("laplace", useCrossValidation = TRUE))
  expect_error(suppressMessages(fitOutcomeModel(weightedClones, crossValidated)), "Cross-validation is not supported")
  expect_error(fitOutcomeModel(fixture$clones, createFitOutcomeModelArgs()), "CcwClones")
  ccwClones <- .createCcwClones(fixture$persons, fixture$clones, fixture$metaData)
  on.exit(close(ccwClones), add = TRUE)
  expect_error(fitOutcomeModel(ccwClones, createFitOutcomeModelArgs()), "segments")
})

test_that("A model without outcome events gives NA estimates with a warning", {
  fixture <- createM4Fixture(n = 40, riskStart = 0, competingEvents = "none")
  fixture$segments$outcome <- 0L
  expect_warning(model <- fitM4(fixture, modelType = "cox", horizons = 10), "No outcome events")
  expect_true(all(is.na(model$estimates$estimate)))
  km <- fitM4(fixture, modelType = "km", horizons = 10)
  expect_equal(getM4Estimate(km, "risk_target_10"), 0)
  expect_equal(getM4Estimate(km, "rd_10"), 0)
})

test_that("print, coef and confint methods", {
  fixture <- createM4Fixture(riskStart = 0, competingEvents = "none")
  model <- fitM4(fixture, modelType = "cox", horizons = c(5, 10))
  messages <- capture_messages(output <- capture.output(print(model)))
  expect_true(any(grepl("Model type: cox", messages)))
  expect_true(any(grepl("Strategy 1 \\(target\\): Initiate within 4 days", messages)))
  expect_true(any(grepl("risk_target_5", output)))
  estimates <- coef(model)
  expect_equal(names(estimates), model$estimates$estimand)
  expect_equal(unname(estimates), model$estimates$estimate)
  intervals <- confint(model)
  expect_equal(dim(intervals), c(nrow(model$estimates), 2))
  expect_equal(unname(intervals[, 1]), model$estimates$lb95)
  expect_equal(unname(confint(model, parm = "hr")[1, ]), c(getM4Estimate(model, "hr", "lb95"), getM4Estimate(model, "hr", "ub95")))
  expect_error(confint(model, level = 0.9), "95%")
})

test_that("plotCumulativeIncidence returns a ggplot", {
  skip_if_not_installed("ggplot2")
  fixture <- createM4Fixture(riskStart = 0, competingEvents = "death")
  weightedClones <- createM4WeightedClones(fixture)
  on.exit(close(weightedClones))
  curves <- computeCumulativeIncidence(weightedClones)
  plot <- plotCumulativeIncidence(curves)
  expect_s3_class(plot, "ggplot")
  fileName <- tempfile(fileext = ".png")
  plot <- plotCumulativeIncidence(curves, labels = c("A", "B"), fileName = fileName)
  expect_true(file.exists(fileName))
  unlink(fileName)
  expect_error(plotCumulativeIncidence(curves, labels = "A"), "labels")
})

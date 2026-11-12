# Tier 2 cross-checks of the outcome models against independent implementations (plan section 6.4).
# Run nightly with CCW_RUN_CROSSCHECKS=true. Fixtures use riskStart = 0.

library(testthat)

skipCrossChecks <- function(package, minimumVersion) {
  skip_on_cran()
  skip_if(Sys.getenv("CCW_RUN_CROSSCHECKS") != "true", "Set CCW_RUN_CROSSCHECKS=true to run the cross-checks")
  skip_if_not_installed(package, minimum_version = minimumVersion)
}

test_that("Weighted KM and Aalen-Johansen match TTE::curve_tte (point estimates)", {
  skipCrossChecks("TTE", "1.1.2")
  horizons <- c(2, 7, 15, 24)
  for (competingEvents in c("none", "death")) {
    fixture <- createM4Fixture(riskStart = 0, competingEvents = competingEvents)
    model <- fitM4(fixture, modelType = "km", horizons = horizons)
    # One row per clone-day (k = t - 1) with its segment's weight, and the event code on the last row.
    rows <- expandM4Rows(fixture$segments, 1)
    rows$event <- ifelse(rows$y == 1, 1L, ifelse(rows$death == 1, 2L, 0L))
    type <- if (competingEvents == "none") "km" else "aj"
    curves <- TTE::curve_tte(rows, time = k, event = event, treatment = x, weights = "weight", type = type, cause = 1)
    for (h in horizons) {
      reference <- summary(curves, time = h + 1)
      expect_equal(getM4Estimate(model, sprintf("risk_target_%s", h)), reference$risk[reference$treatment == "1"], tolerance = 1e-10)
      expect_equal(getM4Estimate(model, sprintf("risk_comparator_%s", h)), reference$risk[reference$treatment == "0"], tolerance = 1e-10)
    }
  }
})

test_that("Weighted KM matches trtswitch::kmest (point estimates)", {
  skipCrossChecks("trtswitch", "0.2.8")
  horizons <- c(2, 7, 15, 24)
  fixture <- createM4Fixture(riskStart = 0, competingEvents = "none")
  model <- fitM4(fixture, modelType = "km", horizons = horizons)
  curves <- trtswitch::kmest(fixture$segments, stratum = "strategyId", time = "startTime", time2 = "stopTime", event = "outcome", weight = "weight")
  for (h in horizons) {
    for (strategyId in 1:2) {
      curve <- curves[curves$strategyId == strategyId & curves$time <= h + 1, ]
      survival <- if (nrow(curve) == 0) 1 else curve$surv[nrow(curve)]
      role <- c("target", "comparator")[strategyId]
      expect_equal(getM4Estimate(model, sprintf("risk_%s_%s", role, h)), 1 - survival, tolerance = 1e-10)
    }
  }
})

test_that("Cox MSM log HR and robust SE match trtswitch::phregr", {
  skipCrossChecks("trtswitch", "0.2.8")
  for (competingEvents in c("none", "death")) {
    fixture <- createM4Fixture(riskStart = 0, competingEvents = competingEvents)
    model <- fitM4(fixture, modelType = "cox", horizons = 10)
    data <- toSurvivalData(fixture$segments)
    reference <- trtswitch::phregr(
      data,
      time = "startTime", time2 = "stopTime", event = "outcome", covariates = "x", weight = "weight",
      id = "personSeqId", ties = "breslow", robust = TRUE
    )
    expect_equal(getM4Estimate(model, "hr", "logEstimate"), reference$parest$beta, tolerance = 1e-6)
    expect_equal(getM4Estimate(model, "hr", "seLogEstimate"), reference$parest$sebeta, tolerance = 1e-6)
  }
})

test_that("Pooled logistic robust covariance matches geepack::geeglm (independence, sandwich)", {
  skipCrossChecks("geepack", "1.3.13")
  fixture <- createM4Fixture(n = 200, riskStart = 0, competingEvents = "none", seed = 3)
  intervalLength <- 2
  model <- fitM4(
    fixture,
    modelType = "logistic", horizons = 9, armByTimeInteraction = TRUE,
    timeBasis = createTimeBasis("linear", intervalLength = intervalLength)
  )
  rows <- expandM4Rows(fixture$segments, intervalLength)
  # geeglm needs each person's rows to be contiguous:
  rows <- rows[order(rows$personSeqId), ]
  rows$xk <- rows$x * rows$k
  reference <- suppressWarnings(geepack::geeglm(
    y ~ x + k + xk,
    family = binomial(), data = rows, weights = weight, id = personSeqId,
    corstr = "independence", std.err = "san.se", scale.fix = TRUE,
    control = geepack::geese.control(epsilon = 1e-12, maxit = 100)
  ))
  expect_equal(unname(model$coefficients), unname(coef(reference)), tolerance = 1e-6)
  expect_equal(unname(model$vcov), unname(stats::vcov(reference)), tolerance = 1e-6)
})

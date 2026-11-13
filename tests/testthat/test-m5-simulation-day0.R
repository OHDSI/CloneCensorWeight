library(testthat)

# Day-0 initiation parameters of simulateCloneCensorWeightData() (initiationHazardDay0 and
# initiationCoefficientsDay0; NULL means the same as on the other days), and computeTrueRisks() with them.

# The tiny DGP of helper-m5-simulation.R with day-0 initiation parameters.
day0TinySettings <- function(dgp, graceEnd = 2, treatmentEffect = log(0.3), initiationHazardDay0 = 0.6, initiationCoefficientsDay0 = c(-1, 1)) {
  settings <- m5TinySettings(dgp, graceEnd = graceEnd, treatmentEffect = treatmentEffect)
  settings$initiationHazardDay0 <- initiationHazardDay0
  settings$initiationCoefficientsDay0 <- initiationCoefficientsDay0
  return(settings)
}

# The CIF under a strategy for one person (scalar linear predictors), by enumerating start days with
# their probabilities, as m5StrategyCif() but with the initiation hazard hA0 on day 0:
# P(A = 0) = hA0 and P(A = a) = (1 - hA0) (1 - hA)^(a - 1) hA for a >= 1.
day0StrategyCif <- function(etaInitiationDay0, etaInitiation, etaOutcome, etaDeath, strategy, graceRegime, settings) {
  hA0 <- stats::plogis(etaInitiationDay0)
  hA <- stats::plogis(etaInitiation)
  hD <- stats::plogis(etaDeath)
  maxDays <- settings$maxDays
  cif <- function(a, untreatedThrough = a) {
    return(m5Cif(m5OutcomeHazards(etaOutcome, a, settings, untreatedThrough), hD))
  }
  pStart <- function(a) if (a == 0) hA0 else (1 - hA0) * (1 - hA)^(a - 1) * hA
  pNotBy <- function(g) (1 - hA0) * (1 - hA)^g
  if (strategy$type == "neverInitiate") {
    return(cif(Inf))
  }
  g <- strategy$graceEnd
  if (strategy$type == "notInitiateWithin") {
    # Starts after g >= 0 follow the hazard of the days after day 0:
    later <- if (g < maxDays) (g + 1):maxDays else integer(0)
    result <- (1 - hA)^length(later) * cif(Inf)
    for (a in later) {
      result <- result + (1 - hA)^(a - g - 1) * hA * cif(a)
    }
    return(result)
  }
  if (graceRegime == "naturalGrace") {
    result <- (if (g == 0) 1 else pNotBy(g - 1)) * cif(g)
    for (a in seq_len(g) - 1) {
      result <- result + pStart(a) * cif(a)
    }
    return(result)
  }
  early <- 0:g
  result <- 0
  for (a in early) {
    result <- result + pStart(a) * cif(a)
  }
  q <- sapply(early, function(a) pStart(a) * prod(1 - m5OutcomeHazards(etaOutcome, a, settings)[1:(g + 1)]))
  q <- q / sum(q)
  for (i in seq_along(early)) {
    result <- result + pNotBy(g) * q[i] * cif(early[i], untreatedThrough = g)
  }
  return(result)
}

# The population risk on each day 0..maxDays (one binary and one normal covariate), as m5ExactRisks().
day0ExactRisks <- function(settings, strategy, graceRegime) {
  eta <- function(hazard, coefficients, x1, x2) stats::qlogis(hazard) + coefficients[1] * x1 + coefficients[2] * x2
  quadrature <- m5GaussHermite()
  risks <- numeric(settings$maxDays + 1)
  for (x1 in 0:1) {
    px1 <- ifelse(x1 == 1, settings$binaryPrevalence, 1 - settings$binaryPrevalence)
    for (k in seq_along(quadrature$nodes)) {
      x2 <- quadrature$nodes[k]
      cif <- day0StrategyCif(
        etaInitiationDay0 = eta(settings$initiationHazardDay0, settings$initiationCoefficientsDay0, x1, x2),
        etaInitiation = eta(settings$initiationHazard, settings$initiationCoefficients, x1, x2),
        etaOutcome = eta(settings$outcomeHazard, settings$outcomeCoefficients, x1, x2),
        etaDeath = eta(settings$deathHazard, settings$deathCoefficients, x1, x2),
        strategy = strategy,
        graceRegime = graceRegime,
        settings = settings
      )
      risks <- risks + px1 * quadrature$weights[k] * cif
    }
  }
  return(risks)
}

test_that("Day-0 initiation parameters are stored, checked, and NULL means the same as other days", {
  data <- simulateCloneCensorWeightData(n = 100, seed = 1)
  settings <- attr(data, "metaData")$simulationSettings
  close(data)
  expect_true(all(c("initiationHazardDay0", "initiationCoefficientsDay0") %in% names(settings)))
  expect_null(settings$initiationHazardDay0)
  expect_null(settings$initiationCoefficientsDay0)

  data <- simulateCloneCensorWeightData(n = 100, nCovariates = 2, initiationHazardDay0 = 0.3, initiationCoefficientsDay0 = c(1, -1), seed = 1)
  settings <- attr(data, "metaData")$simulationSettings
  close(data)
  expect_equal(settings$initiationHazardDay0, 0.3)
  expect_equal(settings$initiationCoefficientsDay0, c(1, -1))

  # Day-0 parameters equal to those of the other days give the same data as NULL:
  treatments <- function(...) {
    data <- simulateCloneCensorWeightData(n = 2000, nCovariates = 2, initiationHazard = 0.05, initiationCoefficients = c(0.5, -0.3), seed = 3, ...)
    on.exit(close(data))
    return(dplyr::collect(data$treatments) |> dplyr::arrange(.data$rowId))
  }
  expect_equal(treatments(initiationHazardDay0 = 0.05, initiationCoefficientsDay0 = c(0.5, -0.3)), treatments())
  expect_equal(treatments(initiationHazardDay0 = 0.05), treatments())

  expect_error(simulateCloneCensorWeightData(n = 10, initiationHazardDay0 = 1), "initiationHazardDay0")
  expect_error(simulateCloneCensorWeightData(n = 10, initiationHazardDay0 = c(0.1, 0.2)), "initiationHazardDay0")
  expect_error(simulateCloneCensorWeightData(n = 10, nCovariates = 2, initiationCoefficientsDay0 = 1), "initiationCoefficientsDay0")
  # Settings created before the day-0 parameters existed are valid for computeTrueRisks():
  settings <- m5TinySettings("a", graceEnd = 2)
  settings$initiationHazardDay0 <- NULL
  settings$initiationCoefficientsDay0 <- NULL
  expect_false("initiationHazardDay0" %in% names(settings))
  strategies <- list(createTreatmentStrategy(1, "Initiate within 2 days", "initiateWithin", graceEnd = 2))
  expect_equal(
    computeTrueRisks(settings, strategies, horizons = 0:8, nTruth = 1000, seed = 1),
    computeTrueRisks(m5TinySettings("a", graceEnd = 2), strategies, horizons = 0:8, nTruth = 1000, seed = 1)
  )
})

test_that("Start days follow the day-0 hazard on day 0 and the daily hazard after it", {
  n <- 1e5
  data <- simulateCloneCensorWeightData(
    n = n, nCovariates = 1, maxDays = 30, binaryPrevalence = 0.4, initiationHazard = 0.05, initiationCoefficients = -1,
    initiationHazardDay0 = 0.2, initiationCoefficientsDay0 = 1.5, seed = 4
  )
  on.exit(close(data))
  cohorts <- dplyr::collect(data$cohorts)
  treatments <- dplyr::collect(data$treatments)
  x <- rep(0, n)
  x[dplyr::pull(data$covariates, "rowId")] <- 1
  startDay <- rep(Inf, n)
  startDay[treatments$rowId] <- treatments$daysToStart
  for (xValue in 0:1) {
    stratum <- x == xValue
    # Day 0 (always observed):
    expected <- plogis(qlogis(0.2) + 1.5 * xValue)
    observed <- mean(startDay[stratum] == 0)
    expect_lt(abs(observed - expected), 4 * sqrt(expected * (1 - expected) / sum(stratum)))
    # Days 1 to 3, among persons not started and still observed:
    expected <- plogis(qlogis(0.05) - xValue)
    for (day in 1:3) {
      atRisk <- stratum & startDay >= day & cohorts$daysToObsEnd >= day
      observed <- mean(startDay[atRisk] == day)
      expect_lt(abs(observed - expected), 4 * sqrt(expected * (1 - expected) / sum(atRisk)), label = sprintf("x = %d, day %d", xValue, day))
    }
  }
})

test_that("computeTrueRisks with day-0 initiation parameters equals a brute-force calculation", {
  nTruth <- 2e5
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within 2 days", "initiateWithin", graceEnd = 2),
    createTreatmentStrategy(2, "Do not initiate within 2 days", "notInitiateWithin", graceEnd = 2),
    createTreatmentStrategy(3, "Never initiate", "neverInitiate"),
    createTreatmentStrategy(4, "Initiate within 4 days", "initiateWithin", graceEnd = 4),
    createTreatmentStrategy(5, "Initiate on day 0", "initiateWithin", graceEnd = 0)
  )
  for (dgp in c("a", "b")) {
    settings <- day0TinySettings(dgp)
    for (graceRegime in c("binaryAtG", "naturalGrace")) {
      truth <- computeTrueRisks(settings, strategies, horizons = 0:8, graceRegime = graceRegime, nTruth = nTruth, seed = 5)
      for (i in seq_along(strategies)) {
        exact <- day0ExactRisks(settings, strategies[[i]], graceRegime)
        estimate <- truth$risk[truth$strategyId == i]
        label <- sprintf("dgp %s, %s, strategy %d", dgp, graceRegime, i)
        expect_true(all(abs(estimate - exact) < 4 * m5MonteCarloSe(exact, nTruth)), label = label)
      }
    }
  }
  # The brute force can tell the day-0 parameters apart under DGP (b):
  withoutDay0 <- m5ExactRisks(m5TinySettings("b", graceEnd = 2), strategies[[1]], "binaryAtG")
  withDay0 <- day0ExactRisks(day0TinySettings("b"), strategies[[1]], "binaryAtG")
  expect_gt(max(abs(withDay0 - withoutDay0) / m5MonteCarloSe(withDay0, nTruth)), 8)
})

test_that("Under DGP (a) the truths of strategies with G = graceEnd do not depend on the initiation process", {
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within 2 days", "initiateWithin", graceEnd = 2),
    createTreatmentStrategy(2, "Do not initiate within 2 days", "notInitiateWithin", graceEnd = 2),
    createTreatmentStrategy(3, "Never initiate", "neverInitiate")
  )
  truths <- function(settings, graceRegime) {
    computeTrueRisks(settings, strategies, horizons = 0:8, graceRegime = graceRegime, nTruth = 1e5, seed = 6)
  }
  for (graceRegime in c("binaryAtG", "naturalGrace")) {
    reference <- truths(m5TinySettings("a", graceEnd = 2), graceRegime)
    expect_equal(truths(day0TinySettings("a"), graceRegime), reference, tolerance = 0)
    other <- day0TinySettings("a", initiationHazardDay0 = 0.05, initiationCoefficientsDay0 = c(2, 0))
    other$initiationHazard <- 0.1
    other$initiationCoefficients <- c(-1, 0.5)
    expect_equal(truths(other, graceRegime), reference, tolerance = 0)
    # Under DGP (b) the initiateWithin truths change (and neverInitiate's does not):
    b <- truths(m5TinySettings("b", graceEnd = 2), graceRegime)
    bDay0 <- truths(day0TinySettings("b"), graceRegime)
    expect_false(isTRUE(all.equal(b$risk[b$strategyId == 1], bDay0$risk[bDay0$strategyId == 1])))
    expect_equal(b[b$strategyId == 3, ], bDay0[bDay0$strategyId == 3, ], tolerance = 0)
  }
})

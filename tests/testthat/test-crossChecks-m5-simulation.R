library(testthat)

# Tier 2 (plan section 6.4): computeTrueRisks() and the natural course of
# simulateCloneCensorWeightData() against simcausal, on the tiny DGP (days 0..8, G = 2). The DAG has,
# for each day t, the outcome Y_t, then death D_t (both ending follow-up), then the absorbing treatment
# A_t, so the outcome comes before death on the same day and the treatment acts from day t + 1 on.

# The simcausal DAG of a tiny DGP. Numbers are inlined into the node formulas with bquote().
crossCheckDag <- function(settings) {
  maxDays <- settings$maxDays
  g <- settings$graceEnd
  eta <- function(hazard, coefficients) {
    return(bquote(.(stats::qlogis(hazard)) + .(coefficients[1]) * X1 + .(coefficients[2]) * X2))
  }
  etaInitiation <- eta(settings$initiationHazard, settings$initiationCoefficients)
  # The day-0 initiation parameters (NULL: as on the other days):
  etaInitiationDay0 <- eta(
    if (is.null(settings$initiationHazardDay0)) settings$initiationHazard else settings$initiationHazardDay0,
    if (is.null(settings$initiationCoefficientsDay0)) settings$initiationCoefficients else settings$initiationCoefficientsDay0
  )
  etaOutcome <- eta(settings$outcomeHazard, settings$outcomeCoefficients)
  etaDeath <- eta(settings$deathHazard, settings$deathCoefficients)
  tau <- settings$treatmentEffect
  rampDays <- settings$effectRampDays
  addNode <- function(dag, ...) {
    return(dag + do.call(simcausal::node, list(...)))
  }

  dag <- simcausal::DAG.empty()
  dag <- addNode(dag, "X1", distr = "rbern", prob = settings$binaryPrevalence)
  dag <- addNode(dag, "X2", distr = "rnorm", mean = 0, sd = 1)
  for (t in 0:maxDays) {
    if (t == 0) {
      outcomeProb <- bquote(plogis(.(etaOutcome)))
      initiationProb <- bquote(plogis(.(etaInitiationDay0)))
    } else {
      if (settings$dgp == "a") {
        effect <- if (t > g) bquote(.(tau) * A[.(g)]) else 0
      } else {
        # Days since initiation: the number of earlier days on treatment.
        effect <- bquote(.(tau) * pmin(sum(A[0:(t - 1)]), .(rampDays)) / .(rampDays))
      }
      outcomeProb <- bquote(plogis(.(etaOutcome) + .(effect)))
      initiationProb <- bquote(ifelse(A[t - 1] == 1, 1, plogis(.(etaInitiation))))
    }
    dag <- addNode(dag, "Y", t = t, distr = "rbern", prob = outcomeProb, EFU = TRUE)
    dag <- addNode(dag, "D", t = t, distr = "rbern", prob = bquote(plogis(.(etaDeath))), EFU = TRUE)
    dag <- addNode(dag, "A", t = t, distr = "rbern", prob = initiationProb)
  }
  # set.DAG() prints the list of vectorized functions:
  utils::capture.output(dag <- suppressMessages(simcausal::set.DAG(dag, vecfun = "pmin", verbose = FALSE)))
  setA <- function(dag, name, days, value) {
    return(dag + simcausal::action(name, nodes = simcausal::node("A", t = days, distr = "rbern", prob = value)))
  }
  dag <- setA(dag, "initiateWithinNaturalGrace", g, 1)
  dag <- setA(dag, "notInitiateWithin", 0:g, 0)
  dag <- setA(dag, "neverInitiate", 0:maxDays, 0)
  return(dag)
}

# Risk of the outcome on each day 0..maxDays in simcausal data.
crossCheckRisks <- function(data, maxDays) {
  outcome <- as.matrix(data[, paste0("Y_", 0:maxDays)]) == 1
  outcome[is.na(outcome)] <- FALSE
  return(colMeans(t(apply(outcome, 1, cumsum)) > 0))
}

test_that("computeTrueRisks agrees with simcausal", {
  skip_on_cran()
  skip_if(Sys.getenv("CCW_RUN_CROSSCHECKS") != "true")
  skip_if_not_installed("simcausal", minimum_version = "0.5.7")

  nTruth <- 2e5
  nSimcausal <- 1e5
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within 2 days", "initiateWithin", graceEnd = 2),
    createTreatmentStrategy(2, "Do not initiate within 2 days", "notInitiateWithin", graceEnd = 2),
    createTreatmentStrategy(3, "Never initiate", "neverInitiate")
  )
  actions <- c("initiateWithinNaturalGrace", "notInitiateWithin", "neverInitiate")
  for (dgp in c("a", "b")) {
    settings <- m5TinySettings(dgp, graceEnd = 2)
    dag <- crossCheckDag(settings)
    simulated <- suppressMessages(simcausal::sim(dag, actions = actions, n = nSimcausal, rndseed = 1, verbose = FALSE))
    # Binary-at-G equals the natural grace period only under DGP (a); simcausal cannot express it.
    graceRegimes <- if (dgp == "a") c("naturalGrace", "binaryAtG") else "naturalGrace"
    for (graceRegime in graceRegimes) {
      truth <- computeTrueRisks(settings, strategies, horizons = 0:settings$maxDays, graceRegime = graceRegime, nTruth = nTruth, seed = 1)
      for (i in seq_along(strategies)) {
        reference <- crossCheckRisks(simulated[[actions[i]]], settings$maxDays)
        estimate <- truth$risk[truth$strategyId == i]
        se <- sqrt(m5MonteCarloSe(reference, nSimcausal)^2 + m5MonteCarloSe(reference, nTruth)^2)
        expect_true(all(abs(estimate - reference) < 4 * se), label = sprintf("dgp %s, %s, strategy %d", dgp, graceRegime, i))
      }
    }
  }
})

test_that("the natural course of simulateCloneCensorWeightData agrees with simcausal", {
  skip_on_cran()
  skip_if(Sys.getenv("CCW_RUN_CROSSCHECKS") != "true")
  skip_if_not_installed("simcausal", minimum_version = "0.5.7")

  n <- 1e5
  for (dgp in c("a", "b")) {
    settings <- m5TinySettings(dgp, graceEnd = 2)
    args <- settings[setdiff(names(settings), "observationEndHazard")]
    data <- do.call(simulateCloneCensorWeightData, c(list(n = n, observationEndHazard = 0, seed = 1), args))
    cohorts <- dplyr::collect(data$cohorts)
    outcomes <- dplyr::collect(data$outcomes)
    treatments <- dplyr::collect(data$treatments)
    close(data)
    simulated <- suppressMessages(simcausal::sim(crossCheckDag(settings), n = n, rndseed = 1, verbose = FALSE))
    se <- function(risk) sqrt(2) * m5MonteCarloSe(risk, n)

    # The outcome (recorded when on or before death):
    risk <- cumsum(tabulate(outcomes$daysToEvent + 1, nbins = settings$maxDays + 1)) / n
    reference <- crossCheckRisks(simulated, settings$maxDays)
    expect_true(all(abs(risk - reference) < 4 * se(reference)), label = sprintf("outcome, dgp %s", dgp))

    # Started on or before day t and still followed after day t (no outcome or death on or before
    # day t). In simcausal, A_t is 1 exactly then (it is NA once follow-up has ended).
    outcomeDay <- rep(Inf, n)
    outcomeDay[outcomes$rowId] <- outcomes$daysToEvent
    deathDay <- ifelse(is.na(cohorts$daysToDeath), Inf, cohorts$daysToDeath)
    startDay <- rep(Inf, n)
    startDay[treatments$rowId] <- treatments$daysToStart
    for (t in 0:settings$maxDays) {
      observed <- mean(startDay <= t & outcomeDay > t & deathDay > t)
      reference <- mean(simulated[[paste0("A_", t)]] %in% 1)
      expect_lt(abs(observed - reference), 4 * se(reference), label = sprintf("started by day %d, dgp %s", t, dgp))
    }
  }
})

test_that("computeTrueRisks and the natural course agree with simcausal with day-0 initiation parameters", {
  skip_on_cran()
  skip_if(Sys.getenv("CCW_RUN_CROSSCHECKS") != "true")
  skip_if_not_installed("simcausal", minimum_version = "0.5.7")

  n <- 1e5
  # DGP (b), where the start day within the grace period changes later risk:
  settings <- m5TinySettings("b", graceEnd = 2)
  settings$initiationHazardDay0 <- 0.6
  settings$initiationCoefficientsDay0 <- c(-1, 1)
  dag <- crossCheckDag(settings)
  actions <- c("initiateWithinNaturalGrace", "notInitiateWithin", "neverInitiate")
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within 2 days", "initiateWithin", graceEnd = 2),
    createTreatmentStrategy(2, "Do not initiate within 2 days", "notInitiateWithin", graceEnd = 2),
    createTreatmentStrategy(3, "Never initiate", "neverInitiate")
  )
  simulated <- suppressMessages(simcausal::sim(dag, actions = actions, n = n, rndseed = 1, verbose = FALSE))
  truth <- computeTrueRisks(settings, strategies, horizons = 0:settings$maxDays, graceRegime = "naturalGrace", nTruth = 2e5, seed = 1)
  for (i in seq_along(strategies)) {
    reference <- crossCheckRisks(simulated[[actions[i]]], settings$maxDays)
    se <- sqrt(m5MonteCarloSe(reference, n)^2 + m5MonteCarloSe(reference, 2e5)^2)
    expect_true(all(abs(truth$risk[truth$strategyId == i] - reference) < 4 * se), label = sprintf("strategy %d", i))
  }

  # The natural course: started on or before day t and still followed after day t.
  args <- settings[setdiff(names(settings), "observationEndHazard")]
  data <- do.call(simulateCloneCensorWeightData, c(list(n = n, observationEndHazard = 0, seed = 1), args))
  cohorts <- dplyr::collect(data$cohorts)
  outcomes <- dplyr::collect(data$outcomes)
  treatments <- dplyr::collect(data$treatments)
  close(data)
  natural <- suppressMessages(simcausal::sim(dag, n = n, rndseed = 2, verbose = FALSE))
  outcomeDay <- rep(Inf, n)
  outcomeDay[outcomes$rowId] <- outcomes$daysToEvent
  deathDay <- ifelse(is.na(cohorts$daysToDeath), Inf, cohorts$daysToDeath)
  startDay <- rep(Inf, n)
  startDay[treatments$rowId] <- treatments$daysToStart
  for (t in 0:settings$maxDays) {
    observed <- mean(startDay <= t & outcomeDay > t & deathDay > t)
    reference <- mean(natural[[paste0("A_", t)]] %in% 1)
    expect_lt(abs(observed - reference), 4 * sqrt(2) * m5MonteCarloSe(reference, n), label = sprintf("started by day %d", t))
  }
  risk <- cumsum(tabulate(outcomes$daysToEvent + 1, nbins = settings$maxDays + 1)) / n
  reference <- crossCheckRisks(natural, settings$maxDays)
  expect_true(all(abs(risk - reference) < 4 * sqrt(2) * m5MonteCarloSe(reference, n)))
})

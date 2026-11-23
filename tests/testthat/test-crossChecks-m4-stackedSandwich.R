# Tier 2 cross-checks of varianceEstimator = "stackedSandwich" (Phase4Designs.md, stacked-sandwich,
# verification tier 2) against independent implementations. Full tier: set CCW_FULL_TESTS=true.
#
# (a) G = 0: initiateWithin(0) against notInitiateWithin(0) is a point treatment with a propensity score,
#     so the stacked SEs equal the corrected sandwich of ipwCoxCSV (HR) and WeightIt's M-estimation
#     (RD and RR at a horizon with no censoring before it).
# (b) geex: an M-estimator written from the person-level data (binary model at G, daily logistic hazard
#     model, their numerators, weighted Kaplan-Meier increments per arm, Cox beta and Breslow increments),
#     with the weights recomputed from the censoring-model parameters.

library(testthat)

skipStackedCrossChecks <- function(package, minimumVersion) {
  skip_on_cran()
  skipUnlessFullTests()
  skip_if_not_installed(package, minimum_version = minimumVersion)
}

stackedEstimate <- function(model, estimand, column) {
  return(model$estimates[[column]][model$estimates$estimand == estimand])
}

test_that("Stacked sandwich at G = 0 matches ipwCoxCSV::ipwCoxInd (HR) and WeightIt::glm_weightit (RD, RR)", {
  skipStackedCrossChecks("ipwCoxCSV", "1.2")
  skipStackedCrossChecks("WeightIt", "2.1.0")
  # Starting on day 0 is the treatment; outcomes from day 1 on; follow-up of at least 15 days, so no
  # censoring before the horizon:
  n <- 200
  horizon <- 10
  set.seed(1)
  x <- cbind(rbinom(n, 1, 0.4), rbinom(n, 1, 0.5), rbinom(n, 1, 0.3))
  start <- runif(n) < plogis(-0.3 + 0.9 * x[, 1] - 0.7 * x[, 2] + 0.5 * x[, 3])
  riskEnd <- sample(15:30, n, replace = TRUE)
  outcomeDay <- 1L + as.integer(floor(runif(n) * riskEnd))
  outcomeDay[runif(n) > plogis(-0.8 + 0.6 * start + 0.5 * x[, 1])] <- NA
  persons <- data.frame(
    rowId = seq_len(n), personSeqId = seq_len(n), riskStart = 0L, riskEnd = as.integer(riskEnd), nominalRiskEnd = 30L,
    outcomeDay = outcomeDay, deathDay = NA_integer_, treatmentDay = ifelse(start, 0L, NA_integer_), comparatorTreatmentDay = NA_integer_
  )
  covariates <- data.frame(rowId = rep(seq_len(n), 3), covariateId = rep(c(1001, 1002, 1003), each = n), covariateValue = as.vector(x))
  data <- createM3Data(list(persons = persons, covariates = covariates[covariates$covariateValue != 0, ]))
  clones <- buildM3Clones(persons, m3Strategies("initiateWithin", "notInitiateWithin", 0), competingEvents = "none")
  on.exit({
    close(clones)
    close(data)
  })
  args <- m3Args("logistic", timeBasis = createTimeBasis(type = "indicator"))
  args$includeCovariateIds <- c(1001, 1002, 1003)
  args$maxCohortSizeForFitting <- 0
  models <- suppressMessages(fitCensoringModels(clones, data, args))
  weightedClones <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(truncation = createWeightTruncation(type = "none"))))
  on.exit(close(weightedClones), add = TRUE)
  fit <- function(msmWeights) {
    fitOutcomeModelArgs <- createFitOutcomeModelArgs(horizons = horizon, varianceEstimator = "stackedSandwich", msmWeights = msmWeights)
    return(suppressMessages(fitOutcomeModel(weightedClones, fitOutcomeModelArgs, cloneCensorWeightData = data)))
  }
  unstabilized <- fit("unstabilized")
  stabilized <- fit("stabilized")

  people <- data.frame(
    A = as.integer(start), X1 = x[, 1], X2 = x[, 2], X3 = x[, 3],
    time = pmin(ifelse(is.na(outcomeDay), Inf, outcomeDay), riskEnd) + 1, status = as.integer(!is.na(outcomeDay))
  )
  # Row 1: unstabilized ATE weights; row 2: stabilized:
  reference <- ipwCoxCSV::ipwCoxInd(people, indA = "A", indX = c("X1", "X2", "X3"), indStatus = "status", indTime = "time", weightType = "ATE", ties = "breslow")
  expect_equal(stackedEstimate(unstabilized, "hr", "logEstimate"), reference[1, 1], tolerance = 1e-6)
  expect_equal(stackedEstimate(unstabilized, "hr", "seLogEstimate"), reference[1, 2], tolerance = 1e-6)
  expect_equal(stackedEstimate(stabilized, "hr", "seLogEstimate"), reference[2, 2], tolerance = 1e-6)

  people$Y <- as.integer(!is.na(outcomeDay) & outcomeDay <= horizon)
  weightIt <- WeightIt::weightit(A ~ X1 + X2 + X3, data = people, method = "glm", estimand = "ATE")
  rd <- WeightIt::glm_weightit(Y ~ A, data = people, weightit = weightIt, family = gaussian, vcov = "asympt")
  rr <- WeightIt::glm_weightit(Y ~ A, data = people, weightit = weightIt, family = quasipoisson, vcov = "asympt")
  expect_equal(stackedEstimate(unstabilized, sprintf("rd_%s", horizon), "estimate"), unname(stats::coef(rd)[2]), tolerance = 1e-6)
  expect_equal(stackedEstimate(unstabilized, sprintf("rd_%s", horizon), "se"), sqrt(stats::vcov(rd)[2, 2]), tolerance = 1e-6)
  expect_equal(stackedEstimate(unstabilized, sprintf("rr_%s", horizon), "logEstimate"), unname(stats::coef(rr)[2]), tolerance = 1e-6)
  expect_equal(stackedEstimate(unstabilized, sprintf("rr_%s", horizon), "seLogEstimate"), sqrt(stats::vcov(rr)[2, 2]), tolerance = 1e-6)
})

# geex reference ---------------------------------------------------------------------------------------
# Strategy 1 = initiateWithin(G) with the binary model at G; strategy 2 = notInitiateWithin(G) with a
# daily logistic hazard (one effect per day 0..G, the covariates, and covariate x day-0 interactions);
# riskStart 0, "eventsFirst", no deaths. Weights recomputed from the parameters on complete days:
# strategy 1, W = 1 / pi(x) after G; strategy 2, W(d) = 1 / S(min(d - 1, G) | x); the numerators pi-bar
# and S-bar (the hazard numerator has one effect per day). The weight at outcome time t is that of day
# t - 1.

# Persons, with two binary covariates and one continuous one (so the normalization s_j is not 1).
simulateStackedPersons <- function(n, seed) {
  set.seed(seed)
  x <- cbind(rbinom(n, 1, 0.4), rbinom(n, 1, 0.5), round(runif(n, 0.5, 3), 1))
  treatmentDay <- rgeom(n, plogis(-2.7 + 0.9 * x[, 1] - 0.6 * x[, 2] + 0.35 * x[, 3]))
  treatmentDay[treatmentDay > 12] <- NA
  riskEnd <- pmin(15L, sample(5:30, n, replace = TRUE))
  outcomeDay <- as.integer(floor(runif(n) * (riskEnd + 1)))
  outcomeDay[runif(n) > plogis(-1.2 + 0.8 * x[, 1] + 0.3 * x[, 3])] <- NA
  persons <- data.frame(
    rowId = seq_len(n), personSeqId = seq_len(n), riskStart = 0L, riskEnd = as.integer(riskEnd), nominalRiskEnd = 15L,
    outcomeDay = outcomeDay, deathDay = NA_integer_, treatmentDay = as.integer(treatmentDay), comparatorTreatmentDay = NA_integer_
  )
  covariates <- data.frame(rowId = rep(seq_len(n), 3), covariateId = rep(c(1001, 1002, 1003), each = n), covariateValue = as.vector(x))
  covariates <- covariates[covariates$covariateValue != 0, ]
  cohorts <- data.frame(
    rowId = persons$rowId, personSeqId = persons$personSeqId, personId = as.character(1000 + persons$rowId),
    cohortStartDate = as.Date("2020-01-01"), daysFromObsStart = 400L, daysToCohortEnd = 0L,
    daysToObsEnd = persons$riskEnd, daysToDeath = persons$deathDay
  )
  data <- createCloneCensorWeightData(
    cohorts = cohorts, outcomes = NULL, treatments = NULL, covariates = covariates,
    covariateRef = data.frame(covariateId = c(1001, 1002, 1003), covariateName = paste("Covariate", 1:3), analysisId = 1:3, conceptId = 0),
    analysisRef = data.frame(
      analysisId = 1:3, analysisName = paste("Analysis", 1:3), domainId = "Condition",
      startDay = -365, endDay = -1, isBinary = c("Y", "Y", "N"), missingMeansZero = "Y"
    ),
    metaData = list(targetId = 1, treatmentIds = c(treatment = 21))
  )
  return(list(persons = persons, x = x, data = data))
}

# The person-level stack: risk sets, clones, the support of each censoring model (a Laplace prior's
# zero coefficients are held at 0), its roots mapped from CCW's stored predictions (exact linear maps),
# and the parameter layout.
geexStackedSetup <- function(sim, models, G, tau) {
  x <- sim$x
  p <- ncol(x)
  inf <- function(v) ifelse(is.na(v), Inf, v)
  a <- inf(sim$persons$treatmentDay)
  o <- inf(sim$persons$outcomeDay)
  e <- sim$persons$riskEnd
  nK <- G + 1
  inBinary <- !(o <= G) & e >= G
  lastDay <- pmin(a, o - 1, e, G)
  dev1 <- ifelse(a <= G, Inf, G)
  dev2 <- ifelse(a <= G, a, Inf)
  survivalTime <- cbind(pmin(o, dev1, e), pmin(o, dev2, e)) + 1
  event <- cbind(o <= pmin(dev1, e), o <= pmin(dev2, e))
  coxTimes <- sort(unique(survivalTime[event]))

  strategyModel <- function(id) models[[which(vapply(models, function(m) m$strategyId, numeric(1)) == id)]]
  binary <- strategyModel(1)
  hazard <- strategyModel(2)
  stopifnot(all(is.finite(hazard$intervals$baseline[seq_len(nK)])))
  supportOf <- function(coefficients, names) names[coefficients[names] != 0]
  ids <- as.character(c(1001, 1002, 1003))
  binarySupport <- match(supportOf(binary$coefficients, ids), ids)
  hazardSupport <- match(supportOf(hazard$coefficients, ids), ids)
  interactionNames <- grep(":", names(hazard$coefficients), fixed = TRUE, value = TRUE)
  interactionSupport <- match(sub(":.*$", "", supportOf(hazard$coefficients, interactionNames)), ids)

  xb <- cbind(1, x[binary$predictions$personSeqId, binarySupport, drop = FALSE])
  binaryFit <- stats::lm.fit(xb, stats::qlogis(binary$predictions$propensityScore))
  xh <- x[hazard$predictions$personSeqId, , drop = FALSE]
  hazardFit <- stats::lm.fit(cbind(1, xh[, hazardSupport, drop = FALSE]), hazard$predictions$linearPredictor)
  # The day-0 interaction (column 1 of the interaction basis):
  stopifnot(all(hazard$interactions$basis[seq_len(nK), 1] == (seq_len(nK) == 1)))
  interactions <- hazard$interactions$predictions
  value <- numeric(nrow(hazard$predictions))
  value[match(interactions$rowId, hazard$predictions$rowId)] <- interactions$value
  interactionFit <- stats::lm.fit(cbind(1, xh[, interactionSupport, drop = FALSE]), value)
  shift <- c(interactionFit$coefficients[1], rep(0, nK - 1))
  roots <- unname(c(
    binaryFit$coefficients, stats::qlogis(binary$numerator$propensityScore),
    hazard$intervals$baseline[seq_len(nK)] + hazardFit$coefficients[1] + shift, hazardFit$coefficients[-1], interactionFit$coefficients[-1],
    hazard$numerator$intervals$baseline[seq_len(nK)]
  ))
  mapResidual <- max(abs(c(binaryFit$residuals, hazardFit$residuals, interactionFit$residuals)))

  sizes <- c(
    alphaB = 1 + length(binarySupport), aB = 1, alphaH = nK + length(hazardSupport) + length(interactionSupport), aH = nK,
    lambda1 = tau, lambda2 = tau, beta = 1, lambda0 = length(coxTimes)
  )
  ends <- cumsum(sizes)
  layout <- lapply(seq_along(sizes), function(j) (ends[j] - sizes[j] + 1):ends[j])
  names(layout) <- names(sizes)
  layout$total <- sum(sizes)
  return(list(
    x = x, n = nrow(x), G = G, nK = nK, tau = tau, inBinary = inBinary, yBinary = as.integer(a <= G),
    lastDay = lastDay, yHazard = as.integer(a == lastDay), survivalTime = survivalTime, event = event, coxTimes = coxTimes,
    binarySupport = binarySupport, hazardSupport = hazardSupport, interactionSupport = interactionSupport,
    layout = layout, censoringRoots = roots, mapResidual = mapResidual
  ))
}

# The hazard design of covariates x on days `ks`: day effects, the covariates and their day-0 interactions.
geexHazardDesign <- function(s, x, ks) {
  return(cbind(outer(ks, seq_len(s$nK) - 1, "==") * 1, outer(rep(1, length(ks)), x[s$hazardSupport]), outer(ks == 0, x[s$interactionSupport]) * 1))
}

# The weights of the Kaplan-Meier rows (stabilized, truncated at `cap` when not NULL) and of the Cox MSM
# (unstabilized: with a cap, the truncated stabilized weight divided by the numerator, as
# truncateMsmWeights = FALSE does) of person i at outcome times 1..maxT.
geexWeights <- function(s, theta, i, maxT, cap) {
  L <- s$layout
  d <- seq_len(maxT) - 1
  x <- s$x[i, ]
  pi <- stats::plogis(sum(c(1, x[s$binarySupport]) * theta[L$alphaB]))
  piBar <- stats::plogis(theta[L$aB])
  h <- stats::plogis(as.vector(geexHazardDesign(s, x, seq_len(s$nK) - 1) %*% theta[L$alphaH]))
  hBar <- stats::plogis(theta[L$aH])
  complete <- outer(pmin(d - 1, s$G), seq_len(s$nK) - 1, ">=")
  w <- list(ifelse(d > s$G, 1 / pi, 1), exp(-as.vector(complete %*% log(1 - h))))
  numerator <- list(ifelse(d > s$G, piBar, 1), exp(as.vector(complete %*% log(1 - hBar))))
  km <- lapply(1:2, function(arm) if (is.null(cap)) w[[arm]] * numerator[[arm]] else pmin(w[[arm]] * numerator[[arm]], cap))
  cox <- lapply(1:2, function(arm) km[[arm]] / numerator[[arm]])
  return(list(km = km, cox = cox))
}

# The estimating function of person i. A Laplace prior (lambda) adds -lambda s_j sign(alpha_j) / N_m to
# the score of each penalized coefficient in the support (s_j the covariate's largest value over model
# m's risk set of N_m persons): psi = U_S - lambda s_S / N on the normalized scale (x / s_j, coefficient
# alpha_j s_j), written on the raw scale (multiplied by s_j, which leaves the SEs unchanged).
geexEstFun <- function(s, cap, lambda) {
  L <- s$layout
  maxT <- max(c(s$tau, s$coxTimes))
  scaleOf <- function(persons, columns) vapply(columns, function(j) max(s$x[persons, j]), numeric(1))
  inHazard <- s$lastDay >= 0
  penaltyB <- if (lambda > 0) lambda * c(0, scaleOf(s$inBinary, s$binarySupport)) / sum(s$inBinary) else 0
  penaltyH <- if (lambda > 0) lambda * c(rep(0, s$nK), scaleOf(inHazard, s$hazardSupport), scaleOf(inHazard, s$interactionSupport)) / sum(inHazard) else 0
  function(data) {
    i <- data$personSeqId[1]
    xb <- c(1, s$x[i, s$binarySupport])
    ks <- if (inHazard[i]) 0:s$lastDay[i] else integer(0)
    Z <- geexHazardDesign(s, s$x[i, ], ks)
    E <- outer(ks, seq_len(s$nK) - 1, "==") * 1
    yRows <- as.integer(ks == s$lastDay[i] & s$yHazard[i] == 1)
    function(theta) {
      out <- numeric(L$total)
      if (s$inBinary[i]) {
        out[L$alphaB] <- xb * (s$yBinary[i] - stats::plogis(sum(xb * theta[L$alphaB]))) - penaltyB * sign(theta[L$alphaB])
        out[L$aB] <- s$yBinary[i] - stats::plogis(theta[L$aB])
      }
      if (inHazard[i]) {
        out[L$alphaH] <- colSums(Z * (yRows - stats::plogis(as.vector(Z %*% theta[L$alphaH])))) - penaltyH * sign(theta[L$alphaH])
        out[L$aH] <- colSums(E * (yRows - stats::plogis(as.vector(E %*% theta[L$aH]))))
      }
      w <- geexWeights(s, theta, i, maxT, cap)
      for (arm in 1:2) {
        Ts <- s$survivalTime[i, arm]
        lambdaIndex <- L[[c("lambda1", "lambda2")[arm]]]
        t <- seq_len(min(s$tau, Ts))
        dN <- as.numeric(s$event[i, arm] & t == Ts)
        out[lambdaIndex[t]] <- out[lambdaIndex[t]] + w$km[[arm]][t] * (dN - theta[lambdaIndex[t]])
        atRisk <- which(s$coxTimes <= Ts)
        tc <- s$coxTimes[atRisk]
        xArm <- as.numeric(arm == 1)
        r <- w$cox[[arm]][tc] * (as.numeric(s$event[i, arm] & tc == Ts) - exp(theta[L$beta] * xArm) * theta[L$lambda0[atRisk]])
        out[L$beta] <- out[L$beta] + xArm * sum(r)
        out[L$lambda0[atRisk]] <- out[L$lambda0[atRisk]] + r
      }
      return(out)
    }
  }
}

# The roots: the censoring roots, CCW's log HR, and the Kaplan-Meier and Breslow increments solved at
# those (closed form).
geexRoots <- function(s, logHr, cap) {
  L <- s$layout
  theta <- numeric(L$total)
  theta[seq_along(s$censoringRoots)] <- s$censoringRoots
  maxT <- max(c(s$tau, s$coxTimes))
  kmNumerator <- kmDenominator <- matrix(0, 2, s$tau)
  coxNumerator <- coxDenominator <- numeric(length(s$coxTimes))
  for (i in seq_len(s$n)) {
    w <- geexWeights(s, theta, i, maxT, cap)
    for (arm in 1:2) {
      Ts <- s$survivalTime[i, arm]
      t <- seq_len(min(s$tau, Ts))
      kmDenominator[arm, t] <- kmDenominator[arm, t] + w$km[[arm]][t]
      if (s$event[i, arm] && Ts <= s$tau) kmNumerator[arm, Ts] <- kmNumerator[arm, Ts] + w$km[[arm]][Ts]
      atRisk <- which(s$coxTimes <= Ts)
      coxDenominator[atRisk] <- coxDenominator[atRisk] + w$cox[[arm]][s$coxTimes[atRisk]] * exp(logHr * (arm == 1))
      if (s$event[i, arm]) {
        j <- match(Ts, s$coxTimes)
        coxNumerator[j] <- coxNumerator[j] + w$cox[[arm]][Ts]
      }
    }
  }
  theta[L$lambda1] <- kmNumerator[1, ] / kmDenominator[1, ]
  theta[L$lambda2] <- kmNumerator[2, ] / kmDenominator[2, ]
  theta[L$beta] <- logHr
  theta[L$lambda0] <- coxNumerator / coxDenominator
  return(theta)
}

# The geex SEs of the log HR, the two risks at the horizon (1 - prod(1 - lambda_t)) and their difference
# (delta method), and the sum of the estimating functions at the roots.
geexStackedSEs <- function(s, cap = NULL, lambda = 0, logHr) {
  L <- s$layout
  theta <- geexRoots(s, logHr, cap)
  estFun <- geexEstFun(s, cap, lambda)
  # Two Richardson steps (default 4) halve the time and keep the bread accurate to about 1e-10:
  m <- geex::m_estimate(
    estFUN = estFun, data = data.frame(personSeqId = seq_len(s$n)), units = "personSeqId", roots = theta, compute_roots = FALSE,
    deriv_control = geex::setup_deriv_control(method.args = list(r = 2))
  )
  psiSum <- Reduce(`+`, lapply(seq_len(s$n), function(i) estFun(data.frame(personSeqId = i))(theta)))
  jacobian <- matrix(0, 4, L$total)
  jacobian[1, L$beta] <- 1
  for (arm in 1:2) {
    lambdaArm <- theta[L[[c("lambda1", "lambda2")[arm]]]]
    jacobian[1 + arm, L[[c("lambda1", "lambda2")[arm]]]] <- prod(1 - lambdaArm) / (1 - lambdaArm)
  }
  jacobian[4, ] <- jacobian[2, ] - jacobian[3, ]
  risk <- 1 - c(prod(1 - theta[L$lambda1]), prod(1 - theta[L$lambda2]))
  return(list(
    estimates = c(logHr, risk, risk[1] - risk[2]),
    se = sqrt(diag(jacobian %*% geex::vcov(m) %*% t(jacobian))),
    psiSum = psiSum
  ))
}

test_that("Stacked sandwich matches a geex stack: a binding fixed cap with truncateMsmWeights = FALSE, and a Laplace prior conditional on its support", {
  skipStackedCrossChecks("geex", "1.1.1")
  G <- 3
  horizon <- 6
  sim <- simulateStackedPersons(150, seed = 2026)
  clones <- buildM3Clones(sim$persons, m3Strategies("initiateWithin", "notInitiateWithin", G), competingEvents = "none")
  on.exit({
    close(clones)
    close(sim$data)
  })
  estimands <- c("hr", sprintf(c("risk_target_%s", "risk_comparator_%s", "rd_%s"), horizon))
  stackedModel <- function(models, truncation) {
    weightsArgs <- createComputeWeightsArgs(stabilized = TRUE, truncation = truncation, truncateMsmWeights = FALSE)
    weightedClones <- suppressMessages(computeCensoringWeights(clones, models, weightsArgs))
    on.exit(close(weightedClones))
    fitOutcomeModelArgs <- createFitOutcomeModelArgs(horizons = horizon, varianceEstimator = "stackedSandwich", msmWeights = "unstabilized")
    model <- suppressMessages(fitOutcomeModel(weightedClones, fitOutcomeModelArgs, cloneCensorWeightData = sim$data))
    value <- ifelse(is.na(model$estimates$logEstimate), model$estimates$estimate, model$estimates$logEstimate)
    se <- ifelse(is.na(model$estimates$seLogEstimate), model$estimates$se, model$estimates$seLogEstimate)
    index <- match(estimands, model$estimates$estimand)
    return(list(estimates = value[index], se = se[index], caps = attr(weightedClones, "metaData")$weights$caps, diagnostics = model$varianceDiagnostics$models))
  }

  # Configuration 1: prior "none", stabilized, a fixed cap on `weight` that binds (the 80th percentile of
  # the weights above 1); the unstabilized Cox weights are the truncated ones divided by the numerator,
  # so a truncated segment keeps only its numerator term.
  args <- m3Args("logistic", timeBasis = createTimeBasis(type = "indicator"), covariateTimeInteraction = "day0")
  models <- suppressMessages(fitCensoringModels(clones, sim$data, args))
  untruncated <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(truncation = createWeightTruncation(type = "none"))))
  weights <- dplyr::pull(untruncated$segments, "weight")
  close(untruncated)
  cap <- round(as.numeric(stats::quantile(weights[weights > 1], 0.8)), 2)
  ccw <- stackedModel(models, createWeightTruncation(type = "fixed", upper = cap))
  expect_gt(min(ccw$caps$truncatedShare), 0)
  setup <- geexStackedSetup(sim, models, G, horizon + 1)
  expect_lt(setup$mapResidual, 1e-8)
  reference <- geexStackedSEs(setup, cap = cap, logHr = ccw$estimates[1])
  expect_equal(ccw$estimates, reference$estimates, tolerance = 1e-6)
  expect_equal(ccw$se, reference$se, tolerance = 1e-6)
  expect_lt(max(abs(reference$psiSum)), 1e-6)

  # Configuration 2: a Laplace prior at a fixed variance on both censoring models, conditional on its
  # support, untruncated:
  variance <- 0.1
  args$prior <- Cyclops::createPrior("laplace", variance = variance, exclude = 0)
  models <- suppressMessages(fitCensoringModels(clones, sim$data, args))
  ccw <- stackedModel(models, createWeightTruncation(type = "none"))
  setup <- geexStackedSetup(sim, models, G, horizon + 1)
  # Some covariate coefficients are 0 and some are not:
  supportSize <- c(length(setup$binarySupport), length(setup$hazardSupport) + length(setup$interactionSupport))
  expect_true(all(supportSize > 0) && sum(supportSize) < 9)
  expect_true(all(ccw$diagnostics$conditionalOnSupport[ccw$diagnostics$prior == "laplace"]))
  expect_lt(setup$mapResidual, 1e-8)
  reference <- geexStackedSEs(setup, lambda = sqrt(2 / variance), logHr = ccw$estimates[1])
  expect_equal(ccw$estimates, reference$estimates, tolerance = 1e-6)
  expect_equal(ccw$se, reference$se, tolerance = 1e-6)
  # The penalized estimating equations hold at CCW's roots (support and s_j as Cyclops used them):
  expect_lt(max(abs(reference$psiSum)), 1e-4)
})

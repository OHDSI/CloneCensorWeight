# Hand-written references that cut across modules.

# The stacked sandwich (Phase4Designs.md, stacked-sandwich, verification tier 1), written per person,
# for the clones of the M3 fixture `fixture` (no deaths, riskStart 0, "eventsFirst") under
# initiateWithin(G) (strategy 1, binary model at G) and notInitiateWithin(G) (strategy 2, daily logistic
# hazard model with an indicator basis and day-0 covariate interactions), both unpenalized:
# - the censoring models are glm fits on the hand-written risk sets, and the unstabilized weights are
#   computed in closed form from their coefficients alpha: W(d) = 1 / pi(X) after G (strategy 1), and
#   W(d) = 1 / S(min(d - 1, G) | X) (strategy 2), on daily segments;
# - the estimates are the log hazard ratio of a weighted Cox model (coxph, Breslow ties) and the weighted
#   Kaplan-Meier risks of both strategies at `horizon`;
# - IJ_i is their infinitesimal-jackknife influence (Cox: weighted score residuals times the inverse
#   information; Kaplan-Meier: written out), U_i the glm scores, A the glm information (block diagonal),
#   and D the numerical Jacobian of the refitted estimates as functions of alpha.
# Returns the `estimates` (log HR, target risk, comparator risk), the stacked influence of each entry
# (rowId) IF_i = IJ_i + U_i' A^-1 D (`influence`, entries x estimates; sum them by personSeqId when a
# person has several entries) and `robustInfluence` (IJ).
refStackedSandwich <- function(fixture, graceEnd, horizon) {
  persons <- fixture$persons
  x <- fixture$x[persons$rowId, , drop = FALSE]
  n <- nrow(persons)
  control <- list(epsilon = 1e-14, maxit = 100)

  # Censoring models: the binary model at G, and the daily hazard model through G.
  binary <- refBinaryRiskSet(persons, graceEnd, "eventsFirst")
  binaryX <- cbind(1, x[binary$rowId, , drop = FALSE])
  binaryFit <- stats::glm.fit(binaryX, binary$y, family = stats::binomial(), control = control)
  hazardZ <- function(day, rowId) {
    xi <- x[rowId, , drop = FALSE]
    return(cbind(1, outer(day, seq_len(graceEnd), "=="), xi, xi * (day == 0)))
  }
  days <- refHazardDays(persons, graceEnd, "eventsFirst")
  hazardX <- hazardZ(days$day, days$rowId)
  hazardFit <- stats::glm.fit(hazardX, days$initiated, family = stats::binomial(), control = control)
  nHazard <- ncol(hazardX)
  alpha <- c(hazardFit$coefficients, binaryFit$coefficients)

  # Daily segments of the clones (t = day + 1), with the outcome on the last day of a clone that ends
  # with it:
  clones <- buildM3CloneTable(persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd))
  segments <- do.call(rbind, lapply(seq_len(nrow(clones)), function(c) {
    t <- seq_len(clones$survivalTime[c])
    data.frame(
      person = clones$rowId[c],
      strategyId = clones$strategyId[c],
      day = t - 1,
      event = as.integer(t == clones$survivalTime[c] & clones$endReason[c] == "outcome")
    )
  }))
  weightsOf <- function(alpha) {
    hazard <- matrix(stats::plogis(hazardZ(rep(0:graceEnd, each = n), rep(seq_len(n), graceEnd + 1)) %*% alpha[seq_len(nHazard)]), n)
    logSurvival <- t(apply(log(1 - hazard), 1, cumsum))
    pi <- as.vector(stats::plogis(cbind(1, x) %*% alpha[-seq_len(nHazard)]))
    last <- pmin(segments$day - 1, graceEnd)
    hazardWeight <- exp(-ifelse(last < 0, 0, logSurvival[cbind(segments$person, pmax(last, 0) + 1)]))
    binaryWeight <- ifelse(segments$day > graceEnd, 1 / pi[segments$person], 1)
    return(ifelse(segments$strategyId == 1, binaryWeight, hazardWeight))
  }
  tau <- horizon + 1
  coxFit <- function(w) {
    survival::coxph(
      survival::Surv(day, day + 1, event) ~ I(strategyId == 1),
      data = segments, weights = w, ties = "breslow", robust = FALSE,
      control = survival::coxph.control(eps = 1e-14, iter.max = 100, toler.chol = 1e-15)
    )
  }
  # The Kaplan-Meier risk at tau of each strategy, and the influence of each segment on it:
  kaplanMeier <- function(w) {
    risk <- numeric(2)
    influence <- matrix(0, nrow(segments), 2)
    for (s in 1:2) {
      atRisk <- segments$strategyId == s & segments$day + 1 <= tau
      time <- segments$day[atRisk] + 1
      y <- tapply(w[atRisk], time, sum)
      d <- tapply(w[atRisk] * segments$event[atRisk], time, sum)
      lambda <- d / y
      survival <- prod(1 - lambda)
      k <- match(time, as.numeric(names(y)))
      influence[atRisk, s] <- survival * w[atRisk] * (segments$event[atRisk] - lambda[k]) / (y[k] * (1 - lambda[k]))
      risk[s] <- 1 - survival
    }
    return(list(risk = risk, influence = influence))
  }
  estimatesOf <- function(alpha) {
    w <- weightsOf(alpha)
    return(c(stats::coef(coxFit(w))[[1]], kaplanMeier(w)$risk))
  }

  w <- weightsOf(alpha)
  fit <- coxFit(w)
  km <- kaplanMeier(w)
  segmentInfluence <- cbind(stats::residuals(fit, type = "score", weighted = TRUE) * fit$var[1, 1], km$influence)
  ij <- matrix(0, n, 3)
  sums <- rowsum(segmentInfluence, segments$person)
  ij[as.integer(rownames(sums)), ] <- sums

  # Scores (persons x parameters, 0 outside a risk set), information and the derivatives:
  scores <- matrix(0, n, length(alpha))
  hazardScores <- rowsum(hazardX * (days$initiated - hazardFit$fitted.values), days$rowId)
  scores[as.integer(rownames(hazardScores)), seq_len(nHazard)] <- hazardScores
  scores[binary$rowId, -seq_len(nHazard)] <- binaryX * (binary$y - binaryFit$fitted.values)
  information <- matrix(0, length(alpha), length(alpha))
  information[seq_len(nHazard), seq_len(nHazard)] <- crossprod(hazardX * sqrt(hazardFit$weights))
  information[-seq_len(nHazard), -seq_len(nHazard)] <- crossprod(binaryX * sqrt(binaryFit$weights))
  derivative <- t(.numericalJacobian(estimatesOf, alpha))
  return(list(
    estimates = estimatesOf(alpha),
    influence = ij + scores %*% solve(information, derivative),
    robustInfluence = ij
  ))
}

# The derivatives of the log weights of the stacked sandwich (Phase4Designs.md, stacked-sandwich, design
# section 2) against numerical derivatives of the weights that computeCensoringWeights() computes. Every
# censoring model and numerator (the stacked sandwich's blocks) is moved along a random direction `v` of
# its parameters other than the covariate-by-time interactions, by changing its stored predictions and
# baselines (`skip`: the "<strategyId> <component>" of blocks left out, such as numerators that are not
# stored in the models). With random segment influences IJ (segments x 2), returns per weight scale
# ("weight", "otherScaleWeight") the largest relative difference between v' D from
# .stackedBlockDerivative() (times the scale's multipliers, .stackedScaleMultipliers()) and
# IJ' dlog w / dt by central differences, and the number of truncated segments per scale; and the
# `diagnostics` of the censoring models (.stackedModelDiagnostics()) with the weights' `caps`.
refStackedWeightDerivatives <- function(clones, models, computeWeightsArgs, cloneCensorWeightData, skip = character(0), step = 1e-6) {
  weighted <- suppressMessages(computeCensoringWeights(clones, models, computeWeightsArgs))
  on.exit(close(weighted))
  data <- .getOutcomeData(weighted)
  stacked <- .stackedCensoringBlocks(data, cloneCensorWeightData, weighted)
  directions <- lapply(stacked$blocks, function(block) {
    v <- stats::rnorm(length(block$parameterX))
    v[block$parameterX > 1 & block$parameterTime > 1] <- 0
    if (paste(block$strategyId, block$component) %in% skip) {
      v[] <- 0
    }
    return(v)
  })
  # The stored pieces of a block's model: the strategy's model, or its `early` or `postGrace` model.
  pieces <- c("hazard model before graceStart" = "early", "binary model at G (postGraceModel)" = "postGrace")
  move <- function(model, block, v) {
    personPart <- block$parameterTime == 1 & block$parameterX > 1
    timePart <- block$parameterX == 1
    a <- as.vector(block$x[, block$parameterX[personPart], drop = FALSE] %*% v[personPart])
    b <- as.vector(block$timeFunctions[, block$parameterTime[timePart], drop = FALSE] %*% v[timePart])
    if (model$modelType == "binary" && block$numerator) {
      model$numerator$propensityScore <- stats::plogis(stats::qlogis(model$numerator$propensityScore) + b)
    } else if (model$modelType == "binary") {
      model$predictions$propensityScore <- stats::plogis(stats::qlogis(model$predictions$propensityScore) + a + b)
    } else if (block$numerator) {
      model$numerator$intervals$baseline <- model$numerator$intervals$baseline + b
    } else {
      model$predictions$linearPredictor <- model$predictions$linearPredictor + a
      model$intervals$baseline <- model$intervals$baseline + b
    }
    return(model)
  }
  segmentsAt <- function(t) {
    moved <- models
    for (i in seq_along(stacked$blocks)) {
      block <- stacked$blocks[[i]]
      piece <- pieces[sub(" numerator$", "", block$component)]
      if (is.na(piece)) {
        moved[[block$strategyId]] <- move(moved[[block$strategyId]], block, t * directions[[i]])
      } else {
        moved[[block$strategyId]][[piece]] <- move(moved[[block$strategyId]][[piece]], block, t * directions[[i]])
      }
    }
    movedWeights <- suppressMessages(computeCensoringWeights(clones, moved, computeWeightsArgs))
    on.exit(close(movedWeights))
    return(.getOutcomeData(movedWeights)$segments)
  }
  plus <- segmentsAt(step)
  minus <- segmentsAt(-step)
  stopifnot(identical(plus$cloneRowId, data$segments$cloneRowId), identical(minus$startTime, data$segments$startTime))
  segments <- data.frame(rowId = data$segments$rowId, personSeqId = data$segments$personSeqId, strategyId = data$segments$strategyId, lastDay = data$segments$stopTime + data$riskStart - 1)
  influence <- matrix(stats::rnorm(2 * nrow(segments)), ncol = 2)
  result <- list(diagnostics = stacked$diagnostics, caps = data$weightsMetaData$caps)
  for (scale in c("weight", "otherScaleWeight")) {
    stabilized <- if (scale == "weight") data$stabilized else !data$stabilized
    multipliers <- .stackedScaleMultipliers(stacked, segments, stabilized)
    analytic <- 0
    for (i in seq_along(stacked$blocks)) {
      block <- stacked$blocks[[i]]
      multiplier <- if (block$numerator) multipliers$numerator else multipliers$denominator
      analytic <- analytic + as.vector(crossprod(directions[[i]], .stackedBlockDerivative(block, segments, influence * multiplier)))
    }
    numerical <- as.vector(crossprod(influence, (log(plus[[scale]]) - log(minus[[scale]])) / (2 * step)))
    result[[scale]] <- list(difference = max(abs(analytic - numerical) / abs(numerical)), truncated = sum(multipliers$denominator == 0))
  }
  return(result)
}

# Copyright 2026 Observational Health Data Sciences and Informatics
#
# This file is part of CloneCensorWeight
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Robust (person-clustered) variance for the outcome models.
#
# All robust variances are HC0 sandwiches clustered on personSeqId with no M / (M - 1) factor, so all
# clones of one person share a cluster. Influence functions are infinitesimal-jackknife (IJ)
# derivatives with respect to multiplying all of a person's segment weights by (1 + epsilon), so the
# variance of an estimate is the sum over persons of its squared influence. This is the robust variance
# of survival::coxph(cluster = personSeqId) and of survival::survfit(id = cloneRowId, cluster =
# personSeqId, influence = TRUE). The weights are treated as fixed.

# Grouped sums ------------------------------------------------------------------------------------

# Sums `values` by the integer `index` (in 1..n) and returns a vector of length n.
.sumBy <- function(values, index, n) {
  result <- numeric(n)
  if (length(values) > 0) {
    sums <- rowsum(values, index)
    result[as.integer(rownames(sums))] <- sums[, 1]
  }
  return(result)
}

# Event-time grid ----------------------------------------------------------------------------------

# Adds, for each segment, the range of grid indices (`idxStart`..`idxStop`) whose risk set it belongs
# to (grid times in (startTime, stopTime]), and the grid index of its stop time (`idxEvent`, NA when
# the stop time is not on the grid). `times` must be sorted and unique.
.indexSegments <- function(segments, times) {
  segments$idxStart <- findInterval(segments$startTime, times) + 1L
  segments$idxStop <- findInterval(segments$stopTime, times)
  segments$idxEvent <- match(segments$stopTime, times)
  return(segments)
}

# Sum of `values` over the risk set at each of the `nTimes` grid times: segment k contributes
# values[k] at grid indices idxStart[k]..idxStop[k].
.riskSetSums <- function(values, idxStart, idxStop, nTimes) {
  delta <- .sumBy(values, idxStart, nTimes + 1) - .sumBy(values, idxStop + 1, nTimes + 1)
  return(cumsum(delta)[seq_len(nTimes)])
}

# Weighted hazard increments ------------------------------------------------------------------------

# Weighted Breslow quantities on the grid for a Cox model with the single binary covariate `x`
# (target indicator) and coefficient `beta`, with events where `isEvent` is TRUE. With beta = 0 and
# the segments of one strategy, `dLambda0` is the weighted Nelson-Aalen increment of that strategy.
.breslow <- function(segments, isEvent, beta, nTimes) {
  relRisk <- exp(beta * segments$x)
  s0 <- .riskSetSums(segments$weight * relRisk, segments$idxStart, segments$idxStop, nTimes)
  s1 <- .riskSetSums(segments$weight * relRisk * segments$x, segments$idxStart, segments$idxStop, nTimes)
  events <- .sumBy(segments$weight[isEvent], segments$idxEvent[isEvent], nTimes)
  positive <- s0 > 0
  dLambda0 <- numeric(nTimes)
  dLambda0[positive] <- events[positive] / s0[positive]
  xBar <- numeric(nTimes)
  xBar[positive] <- s1[positive] / s0[positive]
  return(list(
    beta = beta,
    isEvent = isEvent,
    relRisk = relRisk,
    s0 = s0,
    xBar = xBar,
    events = events,
    dLambda0 = dLambda0,
    # Breslow information for a binary covariate (x^2 = x):
    information = sum(events * (xBar - xBar^2))
  ))
}

# Cox score residuals and robust variance -----------------------------------------------------------

# Weighted score residuals of each segment, including the risk-set compensator:
# w_k [delta_k (x_k - xBar(stop_k)) - exp(beta x_k) sum_{s in (start_k, stop_k]} (x_k - xBar(s)) dLambda0(s)].
.coxScoreResiduals <- function(segments, breslow) {
  cumHazard <- c(0, cumsum(breslow$dLambda0))
  cumXBarHazard <- c(0, cumsum(breslow$xBar * breslow$dLambda0))
  hazardIn <- cumHazard[segments$idxStop + 1] - cumHazard[segments$idxStart]
  xBarHazardIn <- cumXBarHazard[segments$idxStop + 1] - cumXBarHazard[segments$idxStart]
  residuals <- -breslow$relRisk * (segments$x * hazardIn - xBarHazardIn)
  isEvent <- breslow$isEvent
  residuals[isEvent] <- residuals[isEvent] + segments$x[isEvent] - breslow$xBar[segments$idxEvent[isEvent]]
  return(segments$weight * residuals)
}

# Person-level influence of the Cox coefficient: I^-1 times the score residuals summed by person.
.coxCoefficientInfluence <- function(segments, breslow, nPersons) {
  scores <- .sumBy(.coxScoreResiduals(segments, breslow), segments$personIndex, nPersons)
  return(scores / breslow$information)
}

# Influence of functionals of hazard increments ----------------------------------------------------

# Person-level influence of sum_j a[j] * dLambda(s_j), where dLambda(s_j) = scale * dLambda0(s_j) is
# the weighted Breslow (or Nelson-Aalen) increment of `breslow`, holding its coefficient fixed:
# for each segment, w_k a(s) scale / S0(s) for its event, minus
# w_k exp(beta x_k) sum_{s in (start_k, stop_k]} a(s) scale dLambda0(s) / S0(s).
.hazardInfluence <- function(segments, breslow, a, scale, nPersons) {
  positive <- breslow$s0 > 0
  eventCoefficient <- numeric(length(a))
  eventCoefficient[positive] <- a[positive] * scale / breslow$s0[positive]
  riskCoefficient <- eventCoefficient * breslow$dLambda0
  cumRisk <- c(0, cumsum(riskCoefficient))
  contribution <- -segments$weight * breslow$relRisk * (cumRisk[segments$idxStop + 1] - cumRisk[segments$idxStart])
  isEvent <- breslow$isEvent
  contribution[isEvent] <- contribution[isEvent] + segments$weight[isEvent] * eventCoefficient[segments$idxEvent[isEvent]]
  return(.sumBy(contribution, segments$personIndex, nPersons))
}

# Derivative of sum_j a[j] * dLambda(s_j | x) with respect to the Cox coefficient, where
# dLambda(s | x) = exp(beta x) dLambda0(s; beta): d dLambda(s | x) / d beta = dLambda(s | x) (x - xBar(s)).
.hazardBetaDerivative <- function(breslow, a, x) {
  return(sum(a * exp(breslow$beta * x) * breslow$dLambda0 * (x - breslow$xBar)))
}

# Pooled logistic sandwich -------------------------------------------------------------------------

# Bread ((X'WVX)^-1) and per-person scores (sum of w (y - p) X) of a weighted logistic regression.
.logisticScores <- function(design, y, weights, coefficients, personIndex, nPersons) {
  p <- as.vector(1 / (1 + exp(-(design %*% coefficients))))
  information <- crossprod(design, design * (weights * p * (1 - p)))
  scores <- rowsum(design * (weights * (y - p)), personIndex, reorder = TRUE)
  personScores <- matrix(0, nPersons, ncol(design))
  personScores[as.integer(rownames(scores)), ] <- scores
  return(list(bread = solve(information), scores = personScores))
}

# Joint HC0 sandwich clustered on person for one or more models fitted to rows of the same persons
# (each element of `fits` is the output of .logisticScores()): V = B (sum_i U_i U_i') B, with B block
# diagonal. Equal, block by block, to sandwich::vcovCL(glm, cluster = personSeqId, type = "HC0",
# cadjust = FALSE).
.jointSandwich <- function(fits) {
  bread <- .blockDiagonal(lapply(fits, function(fit) fit$bread))
  scores <- do.call(cbind, lapply(fits, function(fit) fit$scores))
  return(bread %*% crossprod(scores) %*% bread)
}

.blockDiagonal <- function(matrices) {
  sizes <- vapply(matrices, nrow, integer(1))
  result <- matrix(0, sum(sizes), sum(sizes))
  offset <- 0
  for (m in matrices) {
    idx <- offset + seq_len(nrow(m))
    result[idx, idx] <- m
    offset <- offset + nrow(m)
  }
  return(result)
}

# Delta method --------------------------------------------------------------------------------------

# Central-difference Jacobian of the vector-valued function f at theta.
.numericalJacobian <- function(f, theta, relativeStep = 1e-6) {
  f0 <- f(theta)
  jacobian <- matrix(0, length(f0), length(theta))
  for (i in seq_along(theta)) {
    step <- relativeStep * max(1, abs(theta[i]))
    up <- theta
    down <- theta
    up[i] <- up[i] + step
    down[i] <- down[i] - step
    jacobian[, i] <- (f(up) - f(down)) / (2 * step)
  }
  return(jacobian)
}

# Inference rows ------------------------------------------------------------------------------------

.z95 <- stats::qnorm(0.975)

.estimateRow <- function(estimand, estimate, logEstimate = NA_real_, seLogEstimate = NA_real_, se = NA_real_,
                         lb95 = NA_real_, ub95 = NA_real_, p = NA_real_) {
  return(tibble(
    estimand = estimand,
    estimate = estimate,
    logEstimate = logEstimate,
    seLogEstimate = seLogEstimate,
    se = se,
    lb95 = lb95,
    ub95 = ub95,
    p = p
  ))
}

# A ratio (HR, OR, RR) with a Wald interval on the log scale.
.ratioRow <- function(estimand, logEstimate, seLogEstimate) {
  if (!is.finite(logEstimate)) {
    logEstimate <- NA_real_
  }
  return(.estimateRow(
    estimand = estimand,
    estimate = exp(logEstimate),
    logEstimate = logEstimate,
    seLogEstimate = seLogEstimate,
    lb95 = exp(logEstimate - .z95 * seLogEstimate),
    ub95 = exp(logEstimate + .z95 * seLogEstimate),
    p = 2 * stats::pnorm(-abs(logEstimate / seLogEstimate))
  ))
}

# A difference (RD, RMST or RMTL difference) with a Wald interval and a test of no difference.
.differenceRow <- function(estimand, estimate, se) {
  return(.estimateRow(
    estimand = estimand,
    estimate = estimate,
    se = se,
    lb95 = estimate - .z95 * se,
    ub95 = estimate + .z95 * se,
    p = 2 * stats::pnorm(-abs(estimate / se))
  ))
}

# An amount (RMST or RMTL of one strategy) with a Wald interval.
.amountRow <- function(estimand, estimate, se) {
  return(.estimateRow(
    estimand = estimand,
    estimate = estimate,
    se = se,
    lb95 = estimate - .z95 * se,
    ub95 = estimate + .z95 * se
  ))
}

# Limits of a probability p with standard error se on the log scale, as survfit's conf.type = "log":
# p exp(-+ z se / p), with the upper limit capped at 1; both limits equal p when se = 0 (below 1e-12,
# to absorb rounding), and NA when p = 0 < se.
.logLimits <- function(p, se) {
  lower <- rep(NA_real_, length(p))
  upper <- rep(NA_real_, length(p))
  known <- !is.na(p) & !is.na(se)
  noVariance <- known & se < 1e-12
  lower[noVariance] <- p[noVariance]
  upper[noVariance] <- p[noVariance]
  positive <- known & !noVariance & p > 0
  lower[positive] <- p[positive] * exp(-.z95 * se[positive] / p[positive])
  upper[positive] <- pmin(1, p[positive] * exp(.z95 * se[positive] / p[positive]))
  return(list(lower = lower, upper = upper))
}

# Risk limits following survfit: conf.type = "log" on the survival S = 1 - risk without competing
# events (risk limits are 1 minus the S limits), and on the cumulative incidence with competing events.
.riskLimits <- function(risk, se, competing) {
  if (competing) {
    return(.logLimits(risk, se))
  }
  limits <- .logLimits(1 - risk, se)
  return(list(lower = 1 - limits$upper, upper = 1 - limits$lower))
}

# A risk with its log-scale interval (see .riskLimits()).
.riskRow <- function(estimand, risk, se, competing) {
  limits <- .riskLimits(risk, se, competing)
  return(.estimateRow(estimand = estimand, estimate = risk, se = se, lb95 = limits$lower, ub95 = limits$upper))
}

# Risk, RD and RR rows for one horizon from the risks and person-level influences of both strategies
# (influence vectors may be NULL for point estimates only).
.riskContrastRows <- function(horizon, risks, influences, competing) {
  se <- function(influence) if (is.null(influence)) NA_real_ else sqrt(sum(influence^2))
  target <- influences[[1]]
  comparator <- influences[[2]]
  hasInfluence <- !is.null(target) && !is.null(comparator)
  logRr <- log(risks[1]) - log(risks[2])
  seLogRr <- if (hasInfluence && is.finite(logRr)) se(target / risks[1] - comparator / risks[2]) else NA_real_
  return(dplyr::bind_rows(
    .riskRow(sprintf("risk_target_%s", horizon), risks[1], se(target), competing),
    .riskRow(sprintf("risk_comparator_%s", horizon), risks[2], se(comparator), competing),
    .differenceRow(sprintf("rd_%s", horizon), risks[1] - risks[2], if (hasInfluence) se(target - comparator) else NA_real_),
    .ratioRow(sprintf("rr_%s", horizon), logRr, seLogRr)
  ))
}

# Risk, RD and RR rows for one horizon by the delta method from a covariance matrix of the risks:
# `gradients` holds the Jacobians (rows: target risk, comparator risk) with respect to the parameters
# whose covariance is `vcov` (NULL for point estimates only).
.deltaContrastRows <- function(horizon, risks, gradients, vcov, competing) {
  if (is.null(vcov) || anyNA(gradients)) {
    return(.riskContrastRows(horizon, risks, list(NULL, NULL), competing))
  }
  variance <- function(gradient) as.numeric(t(gradient) %*% vcov %*% gradient)
  gradientTarget <- gradients[1, ]
  gradientComparator <- gradients[2, ]
  logRr <- log(risks[1]) - log(risks[2])
  seLogRr <- if (is.finite(logRr)) sqrt(variance(gradientTarget / risks[1] - gradientComparator / risks[2])) else NA_real_
  return(dplyr::bind_rows(
    .riskRow(sprintf("risk_target_%s", horizon), risks[1], sqrt(variance(gradientTarget)), competing),
    .riskRow(sprintf("risk_comparator_%s", horizon), risks[2], sqrt(variance(gradientComparator)), competing),
    .differenceRow(sprintf("rd_%s", horizon), risks[1] - risks[2], sqrt(variance(gradientTarget - gradientComparator))),
    .ratioRow(sprintf("rr_%s", horizon), logRr, seLogRr)
  ))
}

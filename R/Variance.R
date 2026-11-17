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
#
# By default (createFitOutcomeModelArgs(robustVarianceEngine = "cyclops")) the robust variances of the Cox
# and pooled logistic coefficients come from Cyclops' vcov(type = "robust"); the sandwiches here are used
# with robustVarianceEngine = "ccw", for the Kaplan-Meier and Aalen-Johansen estimates, and as references
# in the tests.

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

# Person-level influence of sum_j a[j] * dLambda0(s_j), where dLambda0 is the weighted Breslow (or
# Nelson-Aalen) increment of `breslow`, holding its coefficient fixed: for each segment,
# w_k a(s) / S0(s) for its event, minus
# w_k exp(beta x_k) sum_{s in (start_k, stop_k]} a(s) dLambda0(s) / S0(s).
# `a` may also be a matrix with one column per functional; the result is then a persons x functionals
# matrix.
.hazardInfluence <- function(segments, breslow, a, nPersons) {
  isMatrix <- is.matrix(a)
  a <- as.matrix(a)
  positive <- breslow$s0 > 0
  eventCoefficient <- matrix(0, nrow(a), ncol(a))
  eventCoefficient[positive, ] <- a[positive, , drop = FALSE] / breslow$s0[positive]
  cumRisk <- apply(rbind(0, eventCoefficient * breslow$dLambda0), 2, cumsum)
  dim(cumRisk) <- c(nrow(a) + 1, ncol(a))
  contribution <- -(segments$weight * breslow$relRisk) *
    (cumRisk[segments$idxStop + 1, , drop = FALSE] - cumRisk[segments$idxStart, , drop = FALSE])
  isEvent <- which(breslow$isEvent)
  contribution[isEvent, ] <- contribution[isEvent, , drop = FALSE] +
    segments$weight[isEvent] * eventCoefficient[segments$idxEvent[isEvent], , drop = FALSE]
  influence <- matrix(0, nPersons, ncol(a))
  if (nrow(contribution) > 0) {
    sums <- rowsum(contribution, segments$personIndex)
    influence[as.integer(rownames(sums)), ] <- sums
  }
  return(if (isMatrix) influence else influence[, 1])
}

# Pooled logistic sandwich -------------------------------------------------------------------------

# Inverse of an information matrix, or a matrix of NA (with a warning) when it is singular.
.invertInformation <- function(information, label) {
  inverse <- tryCatch(solve(information), error = function(e) NULL)
  if (is.null(inverse)) {
    warning(sprintf("The information matrix of the %s model is singular: its variance is NA", label))
    inverse <- matrix(NA_real_, nrow(information), ncol(information))
  }
  return(inverse)
}

# Sums `values` by person (1..nPersons) and cell (1..nCells): a data frame with columns person, cell and
# value, one row per pair present, sorted by person and cell.
.sumByPersonAndCell <- function(values, personIndex, cell, nCells) {
  key <- (as.numeric(personIndex) - 1) * nCells + (cell - 1)
  keys <- sort(unique(key))
  sums <- rowsum(values, match(key, keys), reorder = TRUE)
  return(data.frame(person = keys %/% nCells + 1, cell = keys %% nCells + 1, value = sums[, 1]))
}

# Joint HC0 sandwich clustered on person for pooled models fitted to rows of the same persons (each
# element of `fits` is a .fitPooledModel() result): V = B (sum_i U_i U_i') B, with B block diagonal and
# U_i = sum_c r_ic X_c the person's score, where r_ic is the sum of the person's score residuals
# w (y - p) in cell c and X_c the cell's design row. Equal, block by block, to
# sandwich::vcovCL(glm, cluster = personSeqId, type = "HC0", cadjust = FALSE). The scores are formed for
# blocks of persons, so memory stays O(block size x coefficients) and no row-level design is built.
.jointSandwich <- function(fits, nPersons, maxEntries = 1e7) {
  bread <- .blockDiagonal(lapply(fits, function(fit) fit$bread))
  blockSize <- max(1, floor(maxEntries / max(nrow(bread), nrow(fits[[1]]$design))))
  meat <- matrix(0, nrow(bread), nrow(bread))
  for (first in seq(0, nPersons - 1, by = blockSize)) {
    size <- min(blockSize, nPersons - first)
    scores <- lapply(fits, function(fit) {
      range <- findInterval(c(first, first + size), fit$residuals$person)
      inBlock <- seq_len(range[2] - range[1]) + range[1]
      residuals <- matrix(0, size, nrow(fit$design))
      residuals[cbind(fit$residuals$person[inBlock] - first, fit$residuals$cell[inBlock])] <- fit$residuals$value[inBlock]
      return(residuals %*% fit$design)
    })
    meat <- meat + crossprod(do.call(cbind, scores))
  }
  return(bread %*% meat %*% bread)
}

# The same joint sandwich from Cyclops' vcov(type = "robust"), for the pooled models `models` (as in
# .jointSandwich()) with identified coefficients `theta` (those of every model, in order), fitted to the
# clone x interval `rows` given by `modelRows` (per model, a list of its logical `subset` of the rows and
# the name of its `response` column). The models are stacked into one lr model
# with a block of columns (each model's identified design columns, intercept included) per model and the
# persons as subjects, on the rows of each model's fitted cells (not those with a fixed hazard of 0 or
# 1) summed by person, cell and response. Its log likelihood is the sum of the models' log likelihoods,
# so its estimates are theirs (it is refitted from `theta`) and its robust variance is the sandwich with
# a block diagonal bread and each person's scores of all models: the joint covariance clustered on
# person. Returns a matrix of NA (with a warning) when the stacked model cannot be fitted.
.cyclopsJointSandwich <- function(models, rows, modelRows, theta, control) {
  outcomes <- list()
  covariates <- list()
  columnOffset <- 0
  rowOffset <- 0
  for (m in seq_along(models)) {
    model <- models[[m]]
    modelRow <- modelRows[[m]]$subset & is.na(model$fixedHazard[rows$cell])
    y <- rows[[modelRows[[m]]$response]][modelRow]
    nCells <- nrow(model$design)
    sums <- lapply(c(1, 0), function(response) {
      sums <- .sumByPersonAndCell(rows$weight[modelRow] * (y == response), rows$personIndex[modelRow], rows$cell[modelRow], nCells)
      sums$y <- response
      return(sums[sums$value > 0, ])
    })
    sums <- do.call(rbind, sums)
    rowIds <- rowOffset + seq_len(nrow(sums))
    outcomes[[m]] <- data.frame(rowId = rowIds, y = sums$y, weights = sums$value, subjectId = sums$person)
    cellDesign <- model$design[sums$cell, , drop = FALSE]
    nonZero <- which(cellDesign != 0, arr.ind = TRUE)
    covariates[[m]] <- data.frame(
      rowId = rowIds[nonZero[, 1]],
      covariateId = columnOffset + nonZero[, 2],
      covariateValue = cellDesign[nonZero]
    )
    columnOffset <- columnOffset + ncol(model$design)
    rowOffset <- rowOffset + nrow(sums)
  }
  outcomes <- do.call(rbind, outcomes)
  covariates <- do.call(rbind, covariates)
  cyclopsData <- Cyclops::convertToCyclopsData(outcomes, covariates, modelType = "lr", addIntercept = FALSE, quiet = TRUE)
  ids <- as.numeric(Cyclops::getCovariateIds(cyclopsData))
  fit <- .fitOutcomeCyclopsModel(
    cyclopsData, Cyclops::createPrior("none"), control,
    unpenalizedIds = ids, description = "stacked pooled logistic model of the robust variance",
    startingCoefficients = theta[ids]
  )
  if (is.null(fit)) {
    return(matrix(NA_real_, length(theta), length(theta)))
  }
  vcov <- matrix(NA_real_, length(theta), length(theta))
  robust <- stats::vcov(fit, type = "robust")
  vcov[ids, ids] <- robust[as.character(ids), as.character(ids)]
  return(vcov)
}

# Stops unless the loaded Cyclops has vcov(fit, type = "robust") (Cyclops C2, plan section 5), that is,
# unless its `vcovMethod` has a `type` argument. An older vcov.cyclopsFit() would silently ignore `type`
# (it takes `...`) and return the model-based variance. Its version cannot tell: the development
# versions with and without it are both 4.0.1.
.checkCyclopsRobustVcov <- function(vcovMethod = utils::getS3method("vcov", "cyclopsFit")) {
  if (!"type" %in% names(formals(vcovMethod))) {
    stop(sprintf(
      "The robust variance with robustVarianceEngine = \"cyclops\" (the default of createFitOutcomeModelArgs()) needs a Cyclops whose vcov() has type = \"robust\" (not yet in a released Cyclops), which the loaded Cyclops %s (%s) does not have. Install a Cyclops that has it, or use createFitOutcomeModelArgs(robustVarianceEngine = \"ccw\"), this package's own sandwich, which gives the same variance",
      utils::packageVersion("Cyclops"), find.package("Cyclops")
    ), call. = FALSE)
  }
  invisible(NULL)
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

# Standard error from `variance`, a diagonal element (`gradient = 1`) or a quadratic form
# t(gradient) %*% vcov %*% gradient of a covariance matrix `vcov` that is positive semi-definite by
# construction (a sandwich or an inverse information), so non-negative in exact arithmetic. When the
# true variance is 0, for example that of the odds ratio of an interval in which both strategies have
# the same rows, rounding can make it slightly negative: a negative value within rounding of 0, relative
# to the largest variance in `vcov` (times the squared size of the gradient), is taken as 0.
.standardError <- function(variance, vcov, gradient = 1) {
  if (!is.na(variance) && variance < 0 && -variance <= 1e-8 * sum(abs(gradient))^2 * max(abs(diag(vcov)), na.rm = TRUE)) {
    return(0)
  }
  return(sqrt(variance))
}

# Two-sided p-value of a Wald test of no effect; NA when the standard error is 0 (an effect that is 0 by
# design, whose estimate is 0 up to rounding).
.waldP <- function(estimate, se) {
  if (isTRUE(se == 0)) {
    return(NA_real_)
  }
  return(2 * stats::pnorm(-abs(estimate / se)))
}

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
    p = .waldP(logEstimate, seLogEstimate)
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
    p = .waldP(estimate, se)
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
# whose covariance is `vcov` (NULL for point estimates only). NA risks, gradients or covariances give NA
# standard errors.
.deltaContrastRows <- function(horizon, risks, gradients, vcov, competing) {
  if (is.null(vcov)) {
    return(.riskContrastRows(horizon, risks, list(NULL, NULL), competing))
  }
  se <- function(gradient) .standardError(as.numeric(t(gradient) %*% vcov %*% gradient), vcov, gradient)
  gradientTarget <- gradients[1, ]
  gradientComparator <- gradients[2, ]
  logRr <- log(risks[1]) - log(risks[2])
  seLogRr <- if (is.finite(logRr)) se(gradientTarget / risks[1] - gradientComparator / risks[2]) else NA_real_
  return(dplyr::bind_rows(
    .riskRow(sprintf("risk_target_%s", horizon), risks[1], se(gradientTarget), competing),
    .riskRow(sprintf("risk_comparator_%s", horizon), risks[2], se(gradientComparator), competing),
    .differenceRow(sprintf("rd_%s", horizon), risks[1] - risks[2], se(gradientTarget - gradientComparator)),
    .ratioRow(sprintf("rr_%s", horizon), logRr, seLogRr)
  ))
}

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

# Product-limit and Aalen-Johansen functionals ------------------------------------------------------

# Point estimates, and gradients with respect to the hazard increments, of the cumulative incidence
# CIF(tau) and of its integral over (0, tau], given the outcome (qY) and competing-event (qD) hazard
# increments on the sorted grid `times`:
#   S(s_j) = prod_{l <= j} (1 - qY_l - qD_l),  CIF(s_j) = sum_{l <= j} S(s_{l-1}) qY_l.
# With qD = 0 this is the Kaplan-Meier: CIF = 1 - S, and the integral of S is tau minus that of CIF.
# Gradients are zero at grid times after tau.
.ajFunctionals <- function(times, qY, qD, tau) {
  nTimes <- length(times)
  result <- list(
    risk = 0,
    survival = 1,
    auc = 0,
    gradRiskY = numeric(nTimes),
    gradRiskD = numeric(nTimes),
    gradAucY = numeric(nTimes),
    gradAucD = numeric(nTimes)
  )
  n <- sum(times <= tau)
  if (n == 0) {
    return(result)
  }
  idx <- seq_len(n)
  qy <- qY[idx]
  survival <- cumprod(1 - qy - qD[idx])
  survivalBefore <- c(1, survival[-n])
  cif <- cumsum(survivalBefore * qy)
  width <- diff(c(times[idx], tau))
  cifWidth <- cif * width
  auc <- sum(cifWidth)
  risk <- cif[n]
  remaining <- tau - times[idx]
  # S(s_{j-1}) / S(s_j); when S(s_j) = 0 nothing can happen after s_j, so the later terms are 0.
  ratio <- ifelse(survival != 0, survivalBefore / survival, 0)
  laterRisk <- ratio * (risk - cif)
  laterAuc <- ratio * (auc - cumsum(cifWidth) - cif * (remaining - width))
  result$risk <- risk
  result$survival <- survival[n]
  result$auc <- auc
  result$gradRiskY[idx] <- survivalBefore - laterRisk
  result$gradRiskD[idx] <- -laterRisk
  result$gradAucY[idx] <- survivalBefore * remaining - laterAuc
  result$gradAucD[idx] <- -laterAuc
  return(result)
}

# Nonparametric curves per strategy ----------------------------------------------------------------

# Weighted Nelson-Aalen increments of one strategy on the pooled grid, for the outcome and (with
# competing events) the competing event.
.strategyHazards <- function(segments, strategyId, nTimes, competing) {
  strategySegments <- segments[segments$strategyId == strategyId, ]
  return(list(
    segments = strategySegments,
    outcome = .breslow(strategySegments, strategySegments$outcome == 1, 0, nTimes),
    death = if (competing) .breslow(strategySegments, strategySegments$competing == 1, 0, nTimes) else NULL
  ))
}

# Kaplan-Meier (no competing events) or Aalen-Johansen risk, survival and area under the cumulative
# incidence at tau for one strategy, with their gradients (see .ajFunctionals()).
.strategyFunctionals <- function(hazards, times, tau) {
  qD <- if (is.null(hazards$death)) numeric(length(times)) else hazards$death$dLambda0
  return(.ajFunctionals(times, hazards$outcome$dLambda0, qD, tau))
}

# Person-level influence of a functional of one strategy's hazard increments with the given gradients
# (or of several functionals, with one gradient column each; see .hazardInfluence()).
.strategyInfluence <- function(hazards, gradientY, gradientD, nPersons) {
  influence <- .hazardInfluence(hazards$segments, hazards$outcome, gradientY, nPersons)
  if (!is.null(hazards$death)) {
    influence <- influence + .hazardInfluence(hazards$segments, hazards$death, gradientD, nPersons)
  }
  return(influence)
}

# Weighted sums over segments at risk at each time in `tau` (startTime < tau <= stopTime).
.atRiskAt <- function(startTime, stopTime, values, tau) {
  sumBelow <- function(time) {
    o <- order(time)
    cumulative <- c(0, cumsum(values[o]))
    return(cumulative[findInterval(tau, time[o], left.open = TRUE) + 1])
  }
  return(sumBelow(startTime) - sumBelow(stopTime))
}

#' Compute cumulative incidence curves
#'
#' @description
#' Computes the weighted Kaplan-Meier risk (`competingEvents = "none"`) or the weighted Aalen-Johansen
#' cumulative incidence of the outcome with death as a competing event (`competingEvents = "death"`)
#' for each strategy from the weighted clone segments, with confidence limits based on the
#' infinitesimal-jackknife standard error clustered on person.
#'
#' @details
#' Confidence limits follow the `survival::survfit()` defaults: `conf.type = "log"` on the survival
#' for the Kaplan-Meier (the risk limits are 1 minus the survival limits), and on the cumulative
#' incidence for the Aalen-Johansen. The standard error is the infinitesimal-jackknife (robust)
#' standard error of `survfit(..., id = cloneRowId, cluster = personSeqId, robust = TRUE)`, treating
#' the weights as fixed.
#'
#' Days are days since index. The risk on day `d` includes events on day `d`, and is read on the
#' survival time scale at `t = d - riskStart + 1`.
#'
#' @param weightedClones A `CcwClones` object with weighted segments, as created by
#'                       `computeCensoringWeights()`.
#' @param horizons       Days since index at which to report the curves. If `NULL`, the curves are
#'                       reported at time zero (day `riskStart - 1`) and at every day with an outcome or
#'                       competing event in the strategy.
#'
#' @return
#' A tibble with one row per strategy and day, with columns `strategyId`, `day`, `atRiskWeighted` (the
#' weighted number of clones at risk), `atRiskCount` (the number of clones at risk), `events` (the
#' weighted number of outcomes since the previous row), `risk`, `lb95` and `ub95`. The attribute
#' `strategyLabels` holds the strategy labels.
#'
#' `atRiskWeighted` sums the weights of the segments at risk, as the `n.risk` of
#' `survival::survfit()` on the segments does: each segment's weight is the weight of its last day
#' (see [computeCensoringWeights()]). That is the clone's weight on the row's day when an outcome or
#' competing event happens on that day, but not in general on other days. On the row at time zero, in
#' particular, it is the sum of the weights of the clones' first segments, which can differ from the
#' sum of their weights on day `riskStart`.
#'
#' Curves computed at `horizons` hold only those days. To plot a curve, compute it with
#' `horizons = NULL`.
#'
#' @example inst/examples/clones.R
#' @example inst/examples/weightedClones.R
#' @examples
#' curves <- computeCumulativeIncidence(weightedClones, horizons = c(7, 30, 60))
#' curves
#' close(weightedClones)
#' close(clones)
#' close(data)
#'
#' @export
computeCumulativeIncidence <- function(weightedClones, horizons = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(weightedClones, weighted = TRUE, add = errorMessages)
  checkmate::assertIntegerish(horizons, lower = 0, min.len = 1, any.missing = FALSE, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  data <- .getOutcomeData(weightedClones)
  competing <- data$competingEvents == "death"
  segments <- data$segments[data$segments$weight > 0, ]
  times <- .eventTimes(segments)
  segments <- .indexSegments(segments, times)
  nPersons <- length(data$personSeqIds)
  if (!is.null(horizons)) {
    .checkHorizons(horizons, data$riskStart)
  }

  curves <- list()
  for (strategyId in 1:2) {
    hazards <- .strategyHazards(segments, strategyId, length(times), competing)
    strategySegments <- hazards$segments
    if (is.null(horizons)) {
      hasEvent <- hazards$outcome$events > 0
      if (competing) {
        hasEvent <- hasEvent | hazards$death$events > 0
      }
      taus <- times[hasEvent]
    } else {
      taus <- horizons - data$riskStart + 1
    }
    risk <- numeric(length(taus))
    se <- numeric(length(taus))
    # The influences of a block of time points are computed in one pass over the segments:
    blockSize <- max(1, floor(1e7 / nrow(strategySegments)))
    for (block in split(seq_along(taus), (seq_along(taus) - 1) %/% blockSize)) {
      functionals <- lapply(taus[block], .strategyFunctionals, hazards = hazards, times = times)
      risk[block] <- vapply(functionals, function(f) f$risk, numeric(1))
      gradient <- function(name) matrix(unlist(lapply(functionals, function(f) f[[name]])), nrow = length(times), ncol = length(block))
      influence <- .strategyInfluence(hazards, gradient("gradRiskY"), gradient("gradRiskD"), nPersons)
      se[block] <- sqrt(colSums(influence^2))
    }
    limits <- .riskLimits(risk, se, competing)
    isOutcome <- strategySegments$outcome == 1
    eventTime <- strategySegments$stopTime[isOutcome]
    eventWeight <- strategySegments$weight[isOutcome]
    o <- order(eventTime)
    cumulativeEvents <- c(0, cumsum(eventWeight[o]))[findInterval(taus, eventTime[o]) + 1]
    curve <- tibble(
      strategyId = strategyId,
      day = taus + data$riskStart - 1,
      atRiskWeighted = .atRiskAt(strategySegments$startTime, strategySegments$stopTime, strategySegments$weight, taus),
      atRiskCount = .atRiskAt(strategySegments$startTime, strategySegments$stopTime, rep(1, nrow(strategySegments)), taus),
      events = diff(c(0, cumulativeEvents)),
      risk = risk,
      lb95 = limits$lower,
      ub95 = limits$upper
    )
    if (is.null(horizons)) {
      firstSegment <- strategySegments$startTime == 0
      origin <- tibble(
        strategyId = strategyId,
        day = data$riskStart - 1,
        atRiskWeighted = sum(strategySegments$weight[firstSegment]),
        atRiskCount = sum(firstSegment),
        events = 0,
        risk = 0,
        lb95 = 0,
        ub95 = 0
      )
      curve <- dplyr::bind_rows(origin, curve)
    }
    curves[[strategyId]] <- curve
  }
  curves <- dplyr::bind_rows(curves)
  attr(curves, "strategyLabels") <- vapply(data$strategies, function(strategy) as.character(strategy$label), character(1))
  return(curves)
}

#' Plot cumulative incidence curves
#'
#' @description
#' Plots the cumulative incidence (risk) curves of both strategies with their 95 percent confidence
#' bands, as computed by [computeCumulativeIncidence()]. Requires the `ggplot2` package.
#'
#' @param curves   A tibble as created by [computeCumulativeIncidence()] with `horizons = NULL`. Curves
#'                 computed at `horizons` hold only those days, and are drawn as steps between them.
#' @param labels   The labels of strategies 1 (target) and 2 (comparator). If `NULL`, the
#'                 `strategyLabels` attribute of `curves` is used.
#' @param fileName Name of the file where the plot should be saved, for example 'plot.png'. See
#'                 [ggplot2::ggsave()] for supported file formats.
#'
#' @return
#' A ggplot object. Use the [ggplot2::ggsave()] function to save to file in a different format.
#'
#' @example inst/examples/clones.R
#' @example inst/examples/weightedClones.R
#' @examples
#' curves <- computeCumulativeIncidence(weightedClones)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plotCumulativeIncidence(curves, labels = c("Initiate within 7 days", "Do not initiate"))
#' }
#' close(weightedClones)
#' close(clones)
#' close(data)
#'
#' @export
plotCumulativeIncidence <- function(curves, labels = NULL, fileName = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(curves, min.rows = 1, add = errorMessages)
  checkmate::assertNames(colnames(curves), must.include = c("strategyId", "day", "risk", "lb95", "ub95"), add = errorMessages)
  checkmate::assertCharacter(labels, len = 2, any.missing = FALSE, null.ok = TRUE, add = errorMessages)
  checkmate::assertCharacter(fileName, len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required to plot cumulative incidence curves. Please install it.")
  }
  if (is.null(labels)) {
    labels <- attr(curves, "strategyLabels")
    if (is.null(labels)) {
      labels <- c("Target", "Comparator")
    }
  }

  # Step functions: each value holds from its day until the next day.
  stepData <- list()
  for (strategyId in 1:2) {
    curve <- curves[curves$strategyId == strategyId, ]
    curve <- curve[order(curve$day), ]
    n <- nrow(curve)
    if (n == 0) {
      next
    }
    dayIdx <- c(1, rep(seq_len(n)[-1], each = 2))
    valueIdx <- c(rep(seq_len(n - 1), each = 2), n)
    stepData[[strategyId]] <- data.frame(
      strategy = labels[strategyId],
      day = curve$day[dayIdx],
      risk = curve$risk[valueIdx],
      lb95 = curve$lb95[valueIdx],
      ub95 = curve$ub95[valueIdx]
    )
  }
  stepData <- do.call(rbind, stepData)
  stepData$strategy <- factor(stepData$strategy, levels = labels)

  plot <- ggplot2::ggplot(stepData, ggplot2::aes(x = .data$day, y = .data$risk, color = .data$strategy, fill = .data$strategy)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lb95, ymax = .data$ub95), alpha = 0.2, color = NA) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::scale_y_continuous("Cumulative incidence", limits = c(0, NA)) +
    ggplot2::scale_x_continuous("Days since index") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 7, height = 5, dpi = 400)
  }
  return(plot)
}

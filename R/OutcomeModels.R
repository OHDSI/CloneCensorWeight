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

#' Fit the outcome model
#'
#' @description
#' Fits the weighted (marginal structural) outcome model to the weighted clones, and estimates risks,
#' risk differences and risk ratios at each horizon, plus the hazard ratio (`"cox"`), odds ratio
#' (`"logistic"`), or restricted mean survival and lost time (`"km"`).
#'
#' @details
#' The model is fitted to the clone `segments` (time scale `t = day - riskStart + 1`, intervals
#' `(startTime, stopTime]`, with the segment weight applying to risk-set membership and to events at
#' `stopTime`). The covariate is the target indicator (strategy 1). Effects are target versus
#' comparator.
#'
#' - `"cox"`: weighted Cox regression (Cyclops `cox_counting`, Breslow ties), with competing events
#'   treated as censoring (cause-specific hazard ratio). Risks are `1 - exp(-Lambda(t | x))` from the
#'   weighted Breslow estimator when `competingEvents = "none"`. When `competingEvents = "death"` a Cox
#'   model of death is also fitted, and risks combine the two cause-specific Breslow hazards by
#'   Aalen-Johansen.
#' - `"logistic"`: weighted pooled logistic regression (Cyclops `lr`) on one row per clone and interval
#'   of `timeBasis$intervalLength` days that the clone enters, with the weight that applies on the
#'   clone's first day in the interval, and covariates intercept, target indicator, time basis, and (if
#'   `armByTimeInteraction`) their interaction. Risks are `1 - prod(1 - h)` without competing events;
#'   with competing death, the outcome model is fitted to rows with no death in the interval, a death
#'   model to all rows, and the risk is that of Young et al. (2020). The model is fitted to the
#'   weighted numbers of rows with and without an event per strategy and interval, which have the
#'   same likelihood as the rows. The risk at a horizon is read at the end of the last complete
#'   interval on or before it, with a warning when the horizon does not end an interval; it is `NA`
#'   (with a warning) when no interval is complete by the horizon or the horizon is beyond the
#'   follow-up. Coefficients the data do not identify (for example the arm-by-time interaction of an
#'   interval only one strategy enters) are `NA`, and so are the risks that need them.
#' - `"km"`: weighted Kaplan-Meier (no competing events) or Aalen-Johansen (competing death) per
#'   strategy.
#'
#' Horizons are days since index and are read at `t = horizon - riskStart + 1`, so the risk includes
#' events on the horizon day.
#'
#' The robust variance (`varianceEstimator = "robust"`) is the HC0 sandwich clustered on person, so all
#' clones of a person share a cluster, treating the weights as fixed. For `"cox"` and `"km"` the risk
#' standard errors are infinitesimal-jackknife standard errors (for `"cox"`, including the estimation
#' of the hazard ratio and the baseline hazards); for `"logistic"` they use the delta method with a
#' numerical Jacobian. With `varianceEstimator = "modelBased"` the hazard ratio or odds ratio uses the
#' inverse of the weighted information; the `"cox"` risks then have no standard errors, and `"km"` is
#' not supported. `varianceEstimator = "bootstrap"` is computed by the bootstrap driver, not by this
#' function.
#'
#' @param weightedClones      A `CcwClones` object with weighted segments, as created by
#'                            `computeCensoringWeights()`.
#' @param fitOutcomeModelArgs An object of type `FitOutcomeModelArgs` as created using
#'                            [createFitOutcomeModelArgs()].
#'
#' @return
#' An object of class `CcwOutcomeModel`: a list with
#' - `modelType` and `varianceEstimator`;
#' - `estimates`: a tibble with one row per estimand (`"hr"` or `"or"`, and per horizon `h`
#'   `"risk_target_<h>"`, `"risk_comparator_<h>"`, `"rd_<h>"`, `"rr_<h>"`, and for `"km"`
#'   `"rmst_<h>"` (without competing events) and `"rmtl_<h>"`, the target minus comparator
#'   differences, with the per-strategy values in `"rmst_target_<h>"`, `"rmst_comparator_<h>"`,
#'   `"rmtl_target_<h>"`, and `"rmtl_comparator_<h>"`), and columns `estimand`, `estimate`,
#'   `logEstimate` and `seLogEstimate` (ratios), `se` (other estimands), `lb95`, `ub95`, and `p` (ratios
#'   and differences);
#' - `coefficients` and `vcov`: the model coefficients and their (robust or model-based) covariance
#'   (`NULL` for `"km"`); death-model coefficients are prefixed with `"death:"`;
#' - `counts`: per strategy, the number of clones, persons, outcomes and competing events, the
#'   person-days, and the weighted person-days, outcomes and competing events;
#' - `metaData`: the strategies, outcome ID, target ID, competing events, risk start, and the
#'   arguments.
#'
#' @references
#' Young JG, Stensrud MJ, Tchetgen Tchetgen EJ, Hernan MA. A causal framework for classical statistical
#' estimands in failure-time settings with competing events. Statistics in Medicine. 2020;39:1199-1236.
#'
#' @examples
#' \dontrun{
#' outcomeModel <- fitOutcomeModel(
#'   weightedClones = weightedClones,
#'   fitOutcomeModelArgs = createFitOutcomeModelArgs(modelType = "km", horizons = c(30, 365))
#' )
#' outcomeModel
#' confint(outcomeModel)
#' }
#'
#' @export
fitOutcomeModel <- function(weightedClones, fitOutcomeModelArgs = createFitOutcomeModelArgs()) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(weightedClones, weighted = TRUE, add = errorMessages)
  checkmate::assertR6(fitOutcomeModelArgs, "FitOutcomeModelArgs", add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (fitOutcomeModelArgs$varianceEstimator == "bootstrap") {
    stop("varianceEstimator = \"bootstrap\" requires refitting the censoring models in every replicate. Use the bootstrap driver")
  }
  if (fitOutcomeModelArgs$varianceEstimator == "modelBased" && fitOutcomeModelArgs$modelType == "km") {
    stop("varianceEstimator = \"modelBased\" is not available for modelType = \"km\"")
  }
  start <- Sys.time()

  data <- .getOutcomeData(weightedClones)
  model <- .fitOutcomeModelCore(data, fitOutcomeModelArgs, computeVariance = TRUE)
  metaData <- attr(weightedClones, "metaData")
  result <- list(
    modelType = fitOutcomeModelArgs$modelType,
    varianceEstimator = fitOutcomeModelArgs$varianceEstimator,
    estimates = model$estimates,
    coefficients = model$coefficients,
    vcov = model$vcov,
    counts = .getOutcomeCounts(weightedClones),
    metaData = list(
      strategies = metaData$strategies,
      outcomeId = metaData$outcomeId,
      targetId = metaData$targetId,
      competingEvents = data$competingEvents,
      riskStart = data$riskStart,
      args = fitOutcomeModelArgs$toList()
    )
  )
  class(result) <- "CcwOutcomeModel"
  delta <- Sys.time() - start
  message(paste("Fitting outcome model took", signif(delta, 3), attr(delta, "units")))
  return(result)
}

# Point estimates only, for the bootstrap driver (R/Bootstrap.R). `personCounts` (optional) is a data
# frame with columns personSeqId and count (non-negative integers): each person's segment weights are
# multiplied by their count (frequency weights; persons not listed get count 0), with the same result
# as copying each person `count` times. `outcomeData` (optional) is the output of
# .getOutcomeData(weightedClones), so that a driver refitting the same segments many times collects them
# only once (`weightedClones` is then not used). The variance estimator in `fitOutcomeModelArgs` is
# ignored. Returns a list with `estimates` (a tibble with estimand, estimate and logEstimate, in the
# same order as fitOutcomeModel()) and `coefficients`.
.fitOutcomeModelPoint <- function(weightedClones, fitOutcomeModelArgs, personCounts = NULL, outcomeData = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  if (is.null(outcomeData)) {
    .assertCcwClones(weightedClones, weighted = TRUE, add = errorMessages)
  }
  checkmate::assertR6(fitOutcomeModelArgs, "FitOutcomeModelArgs", add = errorMessages)
  checkmate::assertDataFrame(personCounts, null.ok = TRUE, add = errorMessages)
  if (!is.null(personCounts)) {
    checkmate::assertNames(colnames(personCounts), must.include = c("personSeqId", "count"), add = errorMessages)
    checkmate::assertIntegerish(personCounts$count, lower = 0, any.missing = FALSE, add = errorMessages)
  }
  checkmate::assertList(outcomeData, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  data <- if (is.null(outcomeData)) .getOutcomeData(weightedClones) else outcomeData
  if (!is.null(personCounts)) {
    data <- .applyPersonCounts(data, personCounts)
  }
  model <- .fitOutcomeModelCore(data, fitOutcomeModelArgs, computeVariance = FALSE)
  return(list(
    estimates = model$estimates[, c("estimand", "estimate", "logEstimate")],
    coefficients = model$coefficients
  ))
}

# Data --------------------------------------------------------------------------------------------

# Collects the segments (all of them, including zero-weight ones) and the settings the outcome models
# need. Each segment has `count`, its person's frequency weight (see .applyPersonCounts()).
.getOutcomeData <- function(weightedClones) {
  metaData <- attr(weightedClones, "metaData")
  riskStart <- weightedClones$clones |>
    dplyr::distinct(.data$riskStart) |>
    pull()
  if (length(riskStart) != 1) {
    stop("All clones must have the same riskStart")
  }
  segments <- weightedClones$segments |>
    select("cloneRowId", "personSeqId", "strategyId", "startTime", "stopTime", "weight", "outcome", "competing") |>
    collect()
  segments <- as.data.frame(segments[order(segments$cloneRowId, segments$startTime), ])
  segments$count <- rep(1, nrow(segments))
  segments$x <- as.integer(segments$strategyId == 1L)
  personSeqIds <- sort(unique(segments$personSeqId))
  segments$personIndex <- match(segments$personSeqId, personSeqIds)
  return(list(
    segments = segments,
    personSeqIds = personSeqIds,
    riskStart = riskStart,
    competingEvents = metaData$competingEvents,
    strategies = metaData$strategies
  ))
}

# Frequency weights: multiplies each segment's weight by its person's count (data frame personSeqId,
# count; persons not listed get count 0) and keeps the count in `count`. The count is kept apart from
# the weight because quantities computed from the rows without weights (the time basis of the pooled
# logistic model) must count each person `count` times, as copies of the person would be counted.
.applyPersonCounts <- function(data, personCounts) {
  count <- personCounts$count[match(data$segments$personSeqId, personCounts$personSeqId)]
  count[is.na(count)] <- 0
  data$segments$count <- count
  data$segments$weight <- data$segments$weight * count
  return(data)
}

# Sorted distinct times of outcomes and competing events with positive weight.
.eventTimes <- function(segments) {
  isEvent <- (segments$outcome == 1 | segments$competing == 1) & segments$weight > 0
  return(sort(unique(segments$stopTime[isEvent])))
}

.checkHorizons <- function(horizons, riskStart) {
  if (any(horizons < riskStart)) {
    stop(sprintf("Horizons must be on or after riskStart (%s)", riskStart))
  }
}

# Per-strategy counts, aggregated in the database.
.getOutcomeCounts <- function(weightedClones) {
  counts <- weightedClones$clones |>
    dplyr::group_by(.data$strategyId) |>
    dplyr::summarise(
      clones = dplyr::n(),
      persons = dplyr::n_distinct(.data$personSeqId),
      outcomes = sum(.data$outcomeCount, na.rm = TRUE),
      competingEvents = sum(dplyr::if_else(.data$endReason == "death", 1, 0), na.rm = TRUE),
      personDays = sum(as.double(.data$survivalTime), na.rm = TRUE),
      .groups = "drop"
    ) |>
    collect()
  weighted <- weightedClones$segments |>
    dplyr::group_by(.data$strategyId) |>
    dplyr::summarise(
      weightedPersonDays = sum(.data$weight * (.data$stopTime - .data$startTime), na.rm = TRUE),
      weightedOutcomes = sum(.data$weight * .data$outcome, na.rm = TRUE),
      weightedCompetingEvents = sum(.data$weight * .data$competing, na.rm = TRUE),
      .groups = "drop"
    ) |>
    collect()
  counts <- dplyr::left_join(counts, weighted, by = "strategyId")
  counts <- counts[order(counts$strategyId), ]
  integerColumns <- c("strategyId", "clones", "persons", "outcomes", "competingEvents")
  counts[integerColumns] <- lapply(counts[integerColumns], as.integer)
  numericColumns <- c("personDays", "weightedPersonDays", "weightedOutcomes", "weightedCompetingEvents")
  counts[numericColumns] <- lapply(counts[numericColumns], as.numeric)
  return(counts)
}

# Models ------------------------------------------------------------------------------------------

.fitOutcomeModelCore <- function(data, fitOutcomeModelArgs, computeVariance) {
  segments <- data$segments
  if (!all(1:2 %in% segments$strategyId[segments$weight > 0])) {
    stop("Both strategies must have clones with positive weights")
  }
  horizons <- fitOutcomeModelArgs$horizons
  .checkHorizons(horizons, data$riskStart)
  competing <- data$competingEvents == "death"
  if (fitOutcomeModelArgs$modelType == "km") {
    return(.fitKm(data, horizons, competing, computeVariance))
  }
  variance <- if (computeVariance) fitOutcomeModelArgs$varianceEstimator else "none"
  if (fitOutcomeModelArgs$modelType == "cox") {
    return(.fitCoxMsm(data, fitOutcomeModelArgs, competing, variance))
  }
  return(.fitLogisticMsm(data, fitOutcomeModelArgs, competing, variance))
}

# Kaplan-Meier / Aalen-Johansen ------------------------------------------------------------------

.fitKm <- function(data, horizons, competing, computeVariance) {
  segments <- data$segments[data$segments$weight > 0, ]
  times <- .eventTimes(segments)
  segments <- .indexSegments(segments, times)
  nPersons <- length(data$personSeqIds)
  hazards <- lapply(1:2, function(strategyId) .strategyHazards(segments, strategyId, length(times), competing))

  se <- function(influence) if (is.null(influence)) NA_real_ else sqrt(sum(influence^2))
  rows <- list()
  for (horizon in horizons) {
    tau <- horizon - data$riskStart + 1
    functionals <- lapply(hazards, .strategyFunctionals, times = times, tau = tau)
    riskInfluences <- list(NULL, NULL)
    aucInfluences <- list(NULL, NULL)
    if (computeVariance) {
      for (i in 1:2) {
        riskInfluences[[i]] <- .strategyInfluence(hazards[[i]], functionals[[i]]$gradRiskY, functionals[[i]]$gradRiskD, nPersons)
        aucInfluences[[i]] <- .strategyInfluence(hazards[[i]], functionals[[i]]$gradAucY, functionals[[i]]$gradAucD, nPersons)
      }
    }
    risks <- vapply(functionals, function(f) f$risk, numeric(1))
    rows[[length(rows) + 1]] <- .riskContrastRows(horizon, risks, riskInfluences, competing)

    auc <- vapply(functionals, function(f) f$auc, numeric(1))
    if (!competing) {
      # RMST = tau - (area under the risk), so its influence is minus that of the area.
      rows[[length(rows) + 1]] <- dplyr::bind_rows(
        .amountRow(sprintf("rmst_target_%s", horizon), tau - auc[1], se(aucInfluences[[1]])),
        .amountRow(sprintf("rmst_comparator_%s", horizon), tau - auc[2], se(aucInfluences[[2]])),
        .differenceRow(sprintf("rmst_%s", horizon), auc[2] - auc[1], .seDifference(aucInfluences))
      )
    }
    rows[[length(rows) + 1]] <- dplyr::bind_rows(
      .amountRow(sprintf("rmtl_target_%s", horizon), auc[1], se(aucInfluences[[1]])),
      .amountRow(sprintf("rmtl_comparator_%s", horizon), auc[2], se(aucInfluences[[2]])),
      .differenceRow(sprintf("rmtl_%s", horizon), auc[1] - auc[2], .seDifference(aucInfluences))
    )
  }
  return(list(estimates = dplyr::bind_rows(rows), coefficients = NULL, vcov = NULL))
}

# Standard error of the target minus comparator difference from person-level influences merged by
# person (persons absent from a strategy have influence 0 there).
.seDifference <- function(influences) {
  if (is.null(influences[[1]]) || is.null(influences[[2]])) {
    return(NA_real_)
  }
  return(sqrt(sum((influences[[1]] - influences[[2]])^2)))
}

# Cox ------------------------------------------------------------------------------------------------

# Fits the weighted counting-process Cox model with the target indicator as the only covariate, using
# Cyclops cox_counting (Breslow ties). Returns the coefficient, or NA (with a warning) when there are
# no events or the fit does not converge.
.fitCyclopsCox <- function(segments, y, prior, control, label) {
  if (!any(y == 1)) {
    warning(sprintf("No %s events: cannot fit the %s Cox model", label, label))
    return(NA_real_)
  }
  outcomes <- data.frame(
    rowId = seq_len(nrow(segments)),
    y = as.numeric(y),
    time = segments$stopTime,
    startTime = segments$startTime,
    weights = segments$weight
  )
  covariates <- data.frame(rowId = which(segments$x == 1), covariateId = 1, covariateValue = 1)
  cyclopsData <- Cyclops::convertToCyclopsData(outcomes, covariates, modelType = "cox_counting", addIntercept = FALSE, quiet = TRUE)
  fit <- tryCatch(
    Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control),
    error = function(e) e$message
  )
  if (is.character(fit) || fit$return_flag != "SUCCESS") {
    warning(sprintf("Fitting the %s Cox model failed: %s", label, if (is.character(fit)) fit else fit$return_flag))
    return(NA_real_)
  }
  return(as.numeric(stats::coef(fit)[["1"]]))
}

.fitCoxMsm <- function(data, fitOutcomeModelArgs, competing, variance) {
  segments <- data$segments[data$segments$weight > 0, ]
  times <- .eventTimes(segments)
  segments <- .indexSegments(segments, times)
  nTimes <- length(times)
  nPersons <- length(data$personSeqIds)
  horizons <- fitOutcomeModelArgs$horizons

  beta <- .fitCyclopsCox(segments, segments$outcome, fitOutcomeModelArgs$prior, fitOutcomeModelArgs$control, "outcome")
  if (is.na(beta)) {
    return(.naEstimates("hr", horizons, coefficientNames = "target"))
  }
  outcome <- .breslow(segments, segments$outcome == 1, beta, nTimes)
  coefficients <- c(target = beta)
  death <- NULL
  if (competing && any(segments$competing == 1)) {
    deathBeta <- .fitCyclopsCox(segments, segments$competing, fitOutcomeModelArgs$prior, fitOutcomeModelArgs$control, "death")
    if (is.na(deathBeta)) {
      return(.naEstimates("hr", horizons, coefficientNames = c("target", "death:target")))
    }
    death <- .breslow(segments, segments$competing == 1, deathBeta, nTimes)
    coefficients <- c(coefficients, `death:target` = deathBeta)
  }

  coefficientInfluence <- NULL
  if (variance == "robust") {
    coefficientInfluence <- cbind(.coxCoefficientInfluence(segments, outcome, nPersons))
    if (!is.null(death)) {
      coefficientInfluence <- cbind(coefficientInfluence, .coxCoefficientInfluence(segments, death, nPersons))
    }
    vcov <- crossprod(coefficientInfluence)
  } else if (variance == "modelBased") {
    information <- c(outcome$information, if (!is.null(death)) death$information)
    vcov <- diag(1 / information, nrow = length(information))
  } else {
    vcov <- NULL
  }
  if (!is.null(vcov)) {
    dimnames(vcov) <- list(names(coefficients), names(coefficients))
  }

  rows <- list(.ratioRow("hr", beta, if (is.null(vcov)) NA_real_ else sqrt(vcov[1, 1])))
  for (horizon in horizons) {
    tau <- horizon - data$riskStart + 1
    risks <- numeric(2)
    influences <- list(NULL, NULL)
    for (strategyId in 1:2) {
      x <- as.integer(strategyId == 1)
      qY <- exp(beta * x) * outcome$dLambda0
      if (competing) {
        qD <- if (is.null(death)) numeric(nTimes) else exp(death$beta * x) * death$dLambda0
        functionals <- .ajFunctionals(times, qY, qD, tau)
        risks[strategyId] <- functionals$risk
        gradientY <- functionals$gradRiskY
        gradientD <- functionals$gradRiskD
      } else {
        inWindow <- times <= tau
        survival <- exp(-sum(qY[inWindow]))
        risks[strategyId] <- 1 - survival
        gradientY <- survival * inWindow
      }
      if (variance == "robust") {
        influence <- .hazardInfluence(segments, outcome, gradientY, exp(beta * x), nPersons) +
          .hazardBetaDerivative(outcome, gradientY, x) * coefficientInfluence[, 1]
        if (!is.null(death)) {
          influence <- influence +
            .hazardInfluence(segments, death, gradientD, exp(death$beta * x), nPersons) +
            .hazardBetaDerivative(death, gradientD, x) * coefficientInfluence[, 2]
        }
        influences[[strategyId]] <- influence
      }
    }
    rows[[length(rows) + 1]] <- .riskContrastRows(horizon, risks, influences, competing)
  }
  return(list(estimates = dplyr::bind_rows(rows), coefficients = coefficients, vcov = vcov))
}

# All estimates NA, when the model cannot be fitted.
.naEstimates <- function(ratio, horizons, coefficientNames) {
  estimands <- c(ratio, unlist(lapply(horizons, function(horizon) {
    sprintf(c("risk_target_%s", "risk_comparator_%s", "rd_%s", "rr_%s"), horizon)
  })))
  estimates <- .estimateRow(estimands, NA_real_)
  coefficients <- stats::setNames(rep(NA_real_, length(coefficientNames)), coefficientNames)
  return(list(estimates = estimates, coefficients = coefficients, vcov = NULL))
}

# Pooled logistic ----------------------------------------------------------------------------------

# One row per clone and interval of `intervalLength` days the clone enters (interval k covers t in
# (k L, (k + 1) L], counted from riskStart), with the weight that applies on the clone's first day in
# the interval (t = k L + 1), the person's frequency weight `count`, and indicators of the outcome and
# of a competing event in the interval. `segments` must hold all segments of each clone, sorted by
# cloneRowId and startTime.
.expandToIntervals <- function(segments, intervalLength) {
  cloneIndex <- match(segments$cloneRowId, unique(segments$cloneRowId))
  last <- !duplicated(cloneIndex, fromLast = TRUE)
  lastSegments <- segments[last, ]
  nIntervals <- ceiling(lastSegments$stopTime / intervalLength)
  rowClone <- rep(seq_len(nrow(lastSegments)), nIntervals)
  interval <- sequence(nIntervals) - 1L
  firstDay <- interval * intervalLength + 1
  # Find the segment with startTime < firstDay <= stopTime of each row's clone:
  span <- max(lastSegments$stopTime) + 1
  segmentKey <- cloneIndex * span + segments$startTime
  rowKey <- rowClone * span + firstDay
  segment <- findInterval(rowKey, segmentKey, left.open = TRUE)
  isLast <- interval == nIntervals[rowClone] - 1
  return(data.frame(
    personIndex = lastSegments$personIndex[rowClone],
    x = lastSegments$x[rowClone],
    interval = interval,
    weight = segments$weight[segment],
    count = lastSegments$count[rowClone],
    y = as.integer(isLast & lastSegments$outcome[rowClone] == 1),
    death = as.integer(isLast & lastSegments$competing[rowClone] == 1)
  ))
}

# Time-basis columns (Phase1Contracts section 4) of the 0-based interval index: `build(k)` computes them
# for the interval indices k, and `names` names them. Columns that are zero at every interval in
# `interval` are dropped. The natural-spline knots (quantiles) and the orthogonal polynomial (centring
# and scaling) depend on the distribution of the rows' interval indices `interval`; each row counts
# `multiplicity` times (its person's frequency weight), so frequency weights give the basis of the data
# with persons copied.
.timeBasisColumns <- function(timeBasis, interval, multiplicity = rep(1, length(interval))) {
  intervals <- sort(unique(interval))
  maxInterval <- max(intervals)
  if (timeBasis$type == "indicator") {
    build <- function(k) {
      columns <- matrix(0, length(k), maxInterval)
      later <- k >= 1 & k <= maxInterval
      columns[cbind(which(later), k[later])] <- 1
      return(columns)
    }
    names <- sprintf("time%d", seq_len(maxInterval))
  } else if (timeBasis$type == "linear") {
    build <- function(k) matrix(as.numeric(k), ncol = 1)
    names <- "time"
  } else {
    sample <- rep(intervals, .sumBy(multiplicity, match(interval, intervals), length(intervals)))
    if (timeBasis$type == "polynomial") {
      basis <- stats::poly(sample, degree = timeBasis$degree, raw = timeBasis$raw)
      names <- sprintf("time%d", seq_len(timeBasis$degree))
    } else {
      basis <- splines::ns(sample, df = timeBasis$df)
      names <- sprintf("time%d", seq_len(ncol(basis)))
    }
    build <- function(k) unclass(stats::predict(basis, k))[, , drop = FALSE]
  }
  keep <- colSums(build(intervals) != 0) > 0
  return(list(
    build = function(k) build(k)[, keep, drop = FALSE],
    names = names[keep]
  ))
}

# Design matrix: intercept, target indicator, time basis and, optionally, target x time basis.
.pooledDesign <- function(x, basisColumns, basisNames, interaction) {
  design <- cbind(1, x, basisColumns)
  names <- c("(Intercept)", "target", basisNames)
  if (interaction && ncol(basisColumns) > 0) {
    design <- cbind(design, x * basisColumns)
    names <- c(names, paste0("target:", basisNames))
  }
  colnames(design) <- names
  return(design)
}

# Columns of `design` identified by the cells with rows (`observed`): a maximal set of linearly
# independent columns, keeping earlier columns first (intercept, target, time basis, interaction), as
# glm does. A cell's linear predictor is identified (`estimable`) when its design row is a linear
# combination of the rows of the observed cells.
.identifiedColumns <- function(design, observed) {
  observedDesign <- design[observed, , drop = FALSE]
  decomposition <- qr(observedDesign)
  kept <- sort(decomposition$pivot[seq_len(decomposition$rank)])
  estimable <- rep(TRUE, nrow(design))
  if (length(kept) < ncol(design)) {
    residual <- qr.resid(qr(t(observedDesign)), t(design))
    estimable <- observed | apply(abs(residual), 2, max) <= 1e-7 * apply(abs(design), 1, max)
  }
  return(list(kept = kept, estimable = estimable))
}

# Fits a weighted pooled logistic model (Cyclops lr) with response `y` to the clone x interval `rows`.
# The design depends on a row only through its cell (strategy and interval; `design` has one row per
# cell), so the model is fitted to the weighted numbers of rows with and without an event in each cell,
# which have the same likelihood as the rows (plan section 4.3). Columns that the cells with rows do not
# identify are dropped and get coefficient NA (see .identifiedColumns()). Unless `variance` is "none",
# also returns the bread (inverse weighted information) and the score residuals summed by person and
# cell. Returns NULL (with a warning) if the model cannot be fitted.
.fitPooledModel <- function(rows, y, design, prior, control, label, variance) {
  if (!any(y == 1)) {
    warning(sprintf("No %s events: cannot fit the %s pooled logistic model", label, label))
    return(NULL)
  }
  nCells <- nrow(design)
  events <- .sumBy(rows$weight * y, rows$cell, nCells)
  nonEvents <- .sumBy(rows$weight * (1 - y), rows$cell, nCells)
  identification <- .identifiedColumns(design, events + nonEvents > 0)
  kept <- identification$kept

  cells <- which(events + nonEvents > 0)
  outcomes <- data.frame(cell = c(cells, cells), y = rep(c(1, 0), each = length(cells)), weights = c(events[cells], nonEvents[cells]))
  outcomes <- outcomes[outcomes$weights > 0, ]
  outcomes$rowId <- seq_len(nrow(outcomes))
  # kept[1] is the intercept:
  covariateColumns <- design[outcomes$cell, kept[-1], drop = FALSE]
  nonZero <- which(covariateColumns != 0, arr.ind = TRUE)
  covariates <- data.frame(
    rowId = nonZero[, 1],
    covariateId = kept[-1][nonZero[, 2]],
    covariateValue = covariateColumns[nonZero]
  )
  cyclopsData <- Cyclops::convertToCyclopsData(outcomes[, c("rowId", "y", "weights")], covariates, modelType = "lr", addIntercept = TRUE, quiet = TRUE)
  fit <- tryCatch(
    Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control),
    error = function(e) e$message
  )
  if (is.character(fit) || fit$return_flag != "SUCCESS") {
    warning(sprintf("Fitting the %s pooled logistic model failed: %s", label, if (is.character(fit)) fit else fit$return_flag))
    return(NULL)
  }
  fitted <- stats::coef(fit)
  coefficients <- stats::setNames(rep(NA_real_, ncol(design)), colnames(design))
  coefficients[1] <- fitted[["(Intercept)"]]
  ids <- setdiff(names(fitted), "(Intercept)")
  coefficients[as.integer(ids)] <- fitted[ids]
  if (length(kept) < ncol(design)) {
    warning(sprintf(
      "Coefficient(s) %s of the %s pooled logistic model are not identified by the data and are NA",
      paste(colnames(design)[-kept], collapse = ", "), label
    ))
  }
  result <- list(coefficients = coefficients, kept = kept, estimable = identification$estimable, design = design[, kept, drop = FALSE])
  if (variance != "none") {
    probability <- stats::plogis(as.vector(result$design %*% coefficients[kept]))
    information <- crossprod(result$design, result$design * ((events + nonEvents) * probability * (1 - probability)))
    result$bread <- .invertInformation(information, sprintf("%s pooled logistic", label))
    result$residuals <- .sumByPersonAndCell(rows$weight * (y - probability[rows$cell]), rows$personIndex, rows$cell, nCells)
  }
  return(result)
}

# Number of complete intervals of `intervalLength` days by each horizon, read at t_h = h - riskStart +
# 1. NA, with a warning, when there is none or when it exceeds the follow-up (`nIntervals`). Warns when
# t_h does not end an interval: the risk is then read at the end of the last complete interval.
.pooledHorizonIntervals <- function(horizons, riskStart, intervalLength, nIntervals) {
  tHorizons <- horizons - riskStart + 1
  intervals <- floor(tHorizons / intervalLength)
  tooEarly <- intervals == 0
  beyond <- intervals > nIntervals
  offGrid <- !tooEarly & !beyond & tHorizons %% intervalLength != 0
  if (any(tooEarly)) {
    warning(sprintf(
      "No complete interval of %d days by horizon(s) %s: their pooled logistic risks are NA",
      intervalLength, paste(horizons[tooEarly], collapse = ", ")
    ))
  }
  if (any(beyond)) {
    warning(sprintf("No follow-up to estimate pooled logistic risks at horizon(s) %s", paste(horizons[beyond], collapse = ", ")))
  }
  if (any(offGrid)) {
    warning(sprintf(
      "Horizon(s) %s do not end an interval of %d days: their pooled logistic risks are read at the end of the last complete interval (day(s) %s)",
      paste(horizons[offGrid], collapse = ", "), intervalLength,
      paste(intervals[offGrid] * intervalLength + riskStart - 1, collapse = ", ")
    ))
  }
  intervals[tooEarly | beyond] <- NA
  return(intervals)
}

# Risks at each horizon (`horizonIntervals` complete intervals, or NA), ordered (target, comparator)
# per horizon, as a function of `theta`, the identified coefficients of the outcome model and then of
# the death model (`layout$index`). A risk is NA when it needs the hazard of a cell that is not
# identified.
.pooledRisks <- function(theta, layout, horizonIntervals) {
  hazard <- function(i) {
    result <- stats::plogis(as.vector(layout$designs[[i]] %*% theta[layout$index[[i]]]))
    result[!layout$estimable[[i]]] <- NA
    return(result)
  }
  nIntervals <- nrow(layout$designs[[1]]) / 2
  hazardY <- hazard(1)
  hazardD <- if (length(layout$designs) > 1) hazard(2) else NULL
  risks <- matrix(NA_real_, 2, length(horizonIntervals))
  for (strategyId in 1:2) {
    cells <- (strategyId == 1) * nIntervals + seq_len(nIntervals)
    if (is.null(hazardD)) {
      cif <- 1 - cumprod(1 - hazardY[cells])
    } else {
      survival <- cumprod((1 - hazardY[cells]) * (1 - hazardD[cells]))
      cif <- cumsum(hazardY[cells] * (1 - hazardD[cells]) * c(1, survival[-nIntervals]))
    }
    risks[strategyId, ] <- cif[horizonIntervals]
  }
  return(as.vector(risks))
}

.fitLogisticMsm <- function(data, fitOutcomeModelArgs, competing, variance) {
  timeBasis <- fitOutcomeModelArgs$timeBasis
  intervalLength <- timeBasis$intervalLength
  horizons <- fitOutcomeModelArgs$horizons
  prior <- fitOutcomeModelArgs$prior
  control <- fitOutcomeModelArgs$control
  if (isTRUE(prior$useCrossValidation)) {
    stop("Cross-validation is not supported for the pooled logistic outcome model, which is fitted to rows aggregated by strategy and interval. Use a prior with a fixed variance")
  }
  rows <- .expandToIntervals(data$segments, intervalLength)
  rows <- rows[rows$weight > 0, ]
  basis <- .timeBasisColumns(timeBasis, rows$interval, rows$count)
  # The design depends on a row only through its cell: the comparator's intervals 0, ..., nIntervals - 1
  # are cells 1, ..., nIntervals, and the target's intervals the next nIntervals cells.
  nIntervals <- max(rows$interval) + 1
  cellInterval <- rep(seq_len(nIntervals) - 1, 2)
  design <- .pooledDesign(rep(0:1, each = nIntervals), basis$build(cellInterval), basis$names, fitOutcomeModelArgs$armByTimeInteraction)
  rows$cell <- as.integer(rows$x * nIntervals + rows$interval + 1)

  outcomeRows <- if (competing) rows$death == 0 else rep(TRUE, nrow(rows))
  models <- list(.fitPooledModel(rows[outcomeRows, ], rows$y[outcomeRows], design, prior, control, "outcome", variance))
  if (is.null(models[[1]])) {
    return(.naEstimates("or", horizons, coefficientNames = colnames(design)))
  }
  if (competing && any(rows$death == 1)) {
    deathModel <- .fitPooledModel(rows, rows$death, design, prior, control, "death", variance)
    if (is.null(deathModel)) {
      return(.naEstimates("or", horizons, coefficientNames = c(colnames(design), paste0("death:", colnames(design)))))
    }
    names(deathModel$coefficients) <- paste0("death:", names(deathModel$coefficients))
    models[[2]] <- deathModel
  }
  coefficients <- unlist(lapply(models, function(model) model$coefficients))
  # The identified coefficients, which the variance and the risks are functions of:
  theta <- unlist(lapply(models, function(model) model$coefficients[model$kept]))
  sizes <- vapply(models, function(model) length(model$kept), integer(1))
  layout <- list(
    designs = lapply(models, function(model) model$design),
    estimable = lapply(models, function(model) model$estimable),
    index = split(seq_along(theta), rep(seq_along(models), sizes))
  )

  thetaVcov <- NULL
  vcov <- NULL
  if (variance != "none") {
    if (variance == "robust") {
      thetaVcov <- .jointSandwich(models, length(data$personSeqIds))
    } else {
      thetaVcov <- .blockDiagonal(lapply(models, function(model) model$bread))
    }
    vcov <- matrix(NA_real_, length(coefficients), length(coefficients), dimnames = list(names(coefficients), names(coefficients)))
    vcov[names(theta), names(theta)] <- thetaVcov
  }

  estimates <- list(.ratioRow("or", coefficients[["target"]], if (is.null(vcov)) NA_real_ else sqrt(vcov["target", "target"])))
  horizonIntervals <- .pooledHorizonIntervals(horizons, data$riskStart, intervalLength, nIntervals)
  riskFunction <- function(theta) .pooledRisks(theta, layout, horizonIntervals)
  risks <- matrix(riskFunction(theta), nrow = 2)
  notIdentified <- is.na(risks) & rep(!is.na(horizonIntervals), each = 2)
  for (strategyId in which(rowSums(notIdentified) > 0)) {
    warning(sprintf(
      "The pooled logistic hazards of the %s strategy are not identified in some intervals: its risks at horizon(s) %s are NA",
      c("target", "comparator")[strategyId], paste(horizons[notIdentified[strategyId, ]], collapse = ", ")
    ))
  }
  jacobian <- if (is.null(thetaVcov)) NULL else .numericalJacobian(riskFunction, theta)
  for (i in seq_along(horizons)) {
    gradients <- if (is.null(jacobian)) NULL else jacobian[2 * i - c(1, 0), , drop = FALSE]
    estimates[[length(estimates) + 1]] <- .deltaContrastRows(horizons[i], risks[, i], gradients, thetaVcov, competing)
  }
  return(list(estimates = dplyr::bind_rows(estimates), coefficients = coefficients, vcov = vcov))
}

# Methods ------------------------------------------------------------------------------------------

#' @export
print.CcwOutcomeModel <- function(x, ...) {
  message(paste("Model type:", x$modelType))
  message(paste("Variance estimator:", x$varianceEstimator))
  message(paste("Competing events:", x$metaData$competingEvents))
  for (i in seq_along(x$metaData$strategies)) {
    message(sprintf("Strategy %d (%s): %s", i, c("target", "comparator")[i], x$metaData$strategies[[i]]$label))
  }
  message("")
  d <- x$estimates
  se <- ifelse(is.na(d$seLogEstimate), d$se, d$seLogEstimate)
  output <- data.frame(d$estimate, d$lb95, d$ub95, se, d$p)
  colnames(output) <- c("Estimate", "lower .95", "upper .95", "SE", "p")
  rownames(output) <- d$estimand
  printCoefmat(output)
  invisible(x)
}

#' @export
coef.CcwOutcomeModel <- function(object, ...) {
  return(stats::setNames(object$estimates$estimate, object$estimates$estimand))
}

#' @export
confint.CcwOutcomeModel <- function(object, parm, level = 0.95, ...) {
  if (level != 0.95) {
    stop("Only supporting 95% confidence interval")
  }
  estimates <- object$estimates
  if (!missing(parm)) {
    estimates <- estimates[estimates$estimand %in% parm, ]
  }
  result <- cbind(estimates$lb95, estimates$ub95)
  dimnames(result) <- list(estimates$estimand, c("2.5 %", "97.5 %"))
  return(result)
}

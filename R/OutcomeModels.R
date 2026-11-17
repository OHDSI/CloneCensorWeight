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
#' risk differences and risk ratios at each horizon, plus the hazard ratio (`"cox"`, the default), the
#' odds ratio (`"logistic"`), and the restricted mean survival and lost time (`"km"` and `"cox"`).
#'
#' @details
#' The model is fitted to the clone `segments` (time scale `t = day - riskStart + 1`, intervals
#' `(startTime, stopTime]`, with the segment weight applying to risk-set membership and to events at
#' `stopTime`). The covariate is the target indicator (strategy 1). Effects are target versus
#' comparator.
#'
#' - `"cox"` (the default): the hazard ratio of a weighted Cox regression (Cyclops `cox_counting`,
#'   Breslow ties), with competing events treated as censoring (cause-specific hazard ratio), and, at
#'   the horizons, exactly the rows of `"km"` on the same weighted clones: the weighted Kaplan-Meier
#'   (or, with competing death, Aalen-Johansen) risks, risk differences, risk ratios and restricted
#'   means. Risks from the Cox model itself (its Breslow baseline hazard) are not reported: they assume
#'   proportional hazards, which a grace period violates by design (see below), so they are biased.
#'   The hazard ratio is reported because meta-analysis and empirical calibration use it. It has
#'   caveats of its own: when the hazards are not proportional it averages the effect over follow-up
#'   in a way that depends on the follow-up and the censoring, and hazards among survivors have a
#'   built-in selection bias (Hernan 2010; Stensrud and Hernan 2020).
#' - `"logistic"`: weighted pooled logistic regression (Cyclops `lr`) on one row per clone and interval
#'   of `timeBasis$intervalLength` days that the clone enters (intervals counted from `riskStart`),
#'   with the weight that applies on the clone's last day in the interval (the weight of its last
#'   segment in the interval), and covariates intercept, target indicator, time basis, and (if
#'   `armByTimeInteraction`) their interaction. The weighted clones must be split on the same interval
#'   grid, by `computeCensoringWeights(outcomeTimeBasis = timeBasis)`. With daily intervals (the
#'   default) every row carries the exact weight of its day, so a saturated time basis per strategy
#'   (an indicator basis with `armByTimeInteraction = TRUE`) gives the weighted Kaplan-Meier or
#'   Aalen-Johansen risks. With an indicator basis, an interval (or, with the interaction, a strategy
#'   and interval) with rows but no events has a hazard of 0, and one in which every row has the event
#'   a hazard of 1, the limits of the maximum likelihood fit; their rows are left out of the fit, and the
#'   coefficients that only they inform (whose estimates are infinite) are `NA`. Risks are
#'   `1 - prod(1 - h)` without competing events;
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
#' The Cox and pooled logistic models, and their robust and model-based variances, use the weights
#' chosen by `fitOutcomeModelArgs$msmWeights` (unstabilized by default; see
#' [createFitOutcomeModelArgs()] for why): the segments' `weight` when it is on that scale, and their
#' `otherScaleWeight` otherwise (see [computeCensoringWeights()]). The Kaplan-Meier and Aalen-Johansen
#' estimates (of `"km"`, and the rows at the horizons of `"cox"`) use `weight`; they are the same with
#' either weight.
#'
#' Horizons are days since index and are read at `t = horizon - riskStart + 1`, so the risk includes
#' events on the horizon day.
#'
#' The `"logistic"` risks are model-based. Without `armByTimeInteraction` they assume an effect of the
#' strategy that is constant over follow-up (proportional odds). A grace period violates this
#' assumption by design: the clones of both strategies follow the same course during the grace period.
#' Their risks, risk differences and risk ratios are then biased, also at horizons within the grace
#' period, so with a grace period use `armByTimeInteraction = TRUE`. With a smooth time basis, such as
#' the default natural spline, the risks can remain biased even then, because each strategy's hazard
#' can change abruptly after the grace period, and with intervals longer than a day, clones that end
#' inside an interval count as at risk for the whole interval (see [createFitOutcomeModelArgs()]). The
#' `"km"` risks, and so the `"cox"` risks, do not make these assumptions. With a strategy whose grace
#' period ends on or after `riskStart`, `"logistic"` warns when the time basis is smooth (natural
#' spline, polynomial or linear), and when no interval ends at the end of the grace period, the end of
#' day `G` (`t = G - riskStart + 1`; day `G - 1` with `sameDayOrdering = "deviationFirst"`), where the
#' clones of an `"initiateWithin(G)"` strategy that have not started are censored. An indicator basis on
#' a grid with an interval ending there, with `armByTimeInteraction = TRUE`, avoids both.
#'
#' The models are fitted with Cyclops using `fitOutcomeModelArgs$control`; a fit that fails is refitted
#' with Cyclops' default algorithm and convergence criterion, with a warning. When the Cox model cannot
#' be fitted (for example because there are no outcomes), the hazard ratio is `NA` and the `"km"` rows
#' are still reported. When the death model of the pooled logistic model with competing death cannot be
#' fitted, the odds ratio is still reported, and the risks are `NA`.
#'
#' The robust variance (`varianceEstimator = "robust"`) is the HC0 sandwich clustered on person, so all
#' clones of a person share a cluster, treating the weights as fixed. For the hazard ratio and the
#' pooled logistic coefficients it is Cyclops' `vcov(fit, type = "robust")` by default, or this
#' package's own sandwich with `robustVarianceEngine = "ccw"` (see [createFitOutcomeModelArgs()]). For
#' `"cox"` and `"km"` the standard errors of the risks, their differences and ratios, and the restricted
#' means are the infinitesimal-jackknife standard errors of the weighted Kaplan-Meier or Aalen-Johansen
#' estimates; for `"logistic"` they use the delta method with a numerical Jacobian. With
#' `varianceEstimator = "modelBased"` the hazard ratio or odds ratio uses the inverse of the weighted
#' information; the `"cox"` rows at the horizons then have no standard errors or limits (the
#' Kaplan-Meier has no model-based variance), and `"km"` is not supported. `varianceEstimator =
#' "bootstrap"` is computed by [bootstrapCcwAnalysis()], not by this function.
#'
#' @param weightedClones      A `CcwClones` object with weighted segments, as created by
#'                            `computeCensoringWeights()`.
#' @param fitOutcomeModelArgs An object of type `FitOutcomeModelArgs` as created using
#'                            [createFitOutcomeModelArgs()].
#'
#' @return
#' An object of class `CcwOutcomeModel`: a list with
#' - `modelType` and `varianceEstimator`;
#' - `estimates`: a tibble with one row per estimand (`"hr"` or `"or"` (with `armByTimeInteraction`,
#'   the odds ratio where the time basis is 0; see [createFitOutcomeModelArgs()]), and per horizon `h`
#'   `"risk_target_<h>"`, `"risk_comparator_<h>"`, `"rd_<h>"`, `"rr_<h>"`, and for `"km"` and `"cox"`
#'   `"rmst_<h>"` (without competing events) and `"rmtl_<h>"`, the target minus comparator
#'   differences, with the per-strategy values in `"rmst_target_<h>"`, `"rmst_comparator_<h>"`,
#'   `"rmtl_target_<h>"`, and `"rmtl_comparator_<h>"`), and columns `estimand`, `estimate`,
#'   `logEstimate` and `seLogEstimate` (ratios), `se` (other estimands), `lb95`, `ub95`, and `p` (ratios
#'   and differences). A contrast that is null by design has a standard error of 0 and `p = NA`: for
#'   example the `"or"` row with `armByTimeInteraction = TRUE` (the odds ratio of the first interval)
#'   when both strategies have the same rows in that interval, as an `"initiateWithin"` strategy and a
#'   strategy of not initiating do on day 0 with `"eventsFirst"`;
#' - `coefficients` and `vcov`: the model coefficients and their (robust or model-based) covariance
#'   (`NULL` for `"km"`): for `"cox"` the target coefficient (log hazard ratio) of the outcome model;
#'   for `"logistic"` those of the outcome model and, with competing death, of the death model,
#'   prefixed with `"death:"`;
#' - `counts`: per strategy, the number of clones, persons, outcomes and competing events, the
#'   person-days, and the weighted person-days, outcomes and competing events (weighted with the
#'   segments' `weight` from [computeCensoringWeights()], not with the weights of `msmWeights`);
#' - `metaData`: the strategies, outcome ID, target ID, competing events, risk start, the arguments,
#'   and, for `"cox"` and `"logistic"`, `msmWeightSummary`: per strategy, the `scale` of the weights the
#'   model used (`msmWeights`), their `cap` on that scale (`NA` when they were not truncated on their own
#'   scale, see `truncateMsmWeights` in [createComputeWeightsArgs()]), `maxWeight`, the largest weight
#'   used (the effective maximum), and `truncatedShare`, the share of the strategy's clones whose weight
#'   was truncated.
#'
#' @references
#' Hernan MA. The hazards of hazard ratios. Epidemiology. 2010;21:13-15.
#'
#' Stensrud MJ, Hernan MA. Why test for proportional hazards? JAMA. 2020;323:1401-1402.
#'
#' Young JG, Stensrud MJ, Tchetgen Tchetgen EJ, Hernan MA. A causal framework for classical statistical
#' estimands in failure-time settings with competing events. Statistics in Medicine. 2020;39:1199-1236.
#'
#' @example inst/examples/clones.R
#' @example inst/examples/weightedClones.R
#' @examples
#' outcomeModel <- fitOutcomeModel(
#'   weightedClones = weightedClones,
#'   fitOutcomeModelArgs = createFitOutcomeModelArgs(modelType = "km", horizons = c(30, 60))
#' )
#' outcomeModel
#' confint(outcomeModel)
#' # The pooled logistic model with a daily indicator basis and strategy-by-time interactions, which
#' # follow the change in the hazards at the end of the grace period:
#' logisticModel <- fitOutcomeModel(
#'   weightedClones = weightedClones,
#'   fitOutcomeModelArgs = createFitOutcomeModelArgs(
#'     modelType = "logistic",
#'     horizons = c(30, 60),
#'     timeBasis = createTimeBasis(type = "indicator", intervalLength = 1),
#'     armByTimeInteraction = TRUE
#'   )
#' )
#' coef(logisticModel)
#' close(weightedClones)
#' close(clones)
#' close(data)
#'
#' @export
fitOutcomeModel <- function(weightedClones, fitOutcomeModelArgs = createFitOutcomeModelArgs()) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(weightedClones, weighted = TRUE, add = errorMessages)
  checkmate::assertR6(fitOutcomeModelArgs, "FitOutcomeModelArgs", add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (fitOutcomeModelArgs$varianceEstimator == "bootstrap") {
    stop("varianceEstimator = \"bootstrap\" requires refitting the censoring models in every replicate. Use the bootstrap driver, bootstrapCcwAnalysis()")
  }
  if (fitOutcomeModelArgs$varianceEstimator == "modelBased" && fitOutcomeModelArgs$modelType == "km") {
    stop("varianceEstimator = \"modelBased\" is not available for modelType = \"km\"")
  }
  if (.needsCyclopsRobustVcov(fitOutcomeModelArgs)) {
    .checkCyclopsRobustVcov()
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
  if (fitOutcomeModelArgs$modelType != "km") {
    result$metaData$msmWeightSummary <- .summarizeMsmWeights(data, fitOutcomeModelArgs$msmWeights)
  }
  class(result) <- "CcwOutcomeModel"
  delta <- Sys.time() - start
  message(paste("Fitting outcome model took", signif(delta, 3), attr(delta, "units")))
  return(result)
}

# Point estimates only. The bootstrap driver (R/Bootstrap.R) resamples by copying persons and calls this
# with the weighted clones of each replicate. `personCounts` and `outcomeData` serve the frequency-weight
# bootstrap of plan section 4.4, which is not built yet; the test that frequency weights equal copies
# (test-m4-outcomes.R) uses them. `personCounts` (optional) is a data frame with columns personSeqId and
# count (non-negative integers): each person's segment weights are multiplied by their count (persons
# not listed get count 0), with the same result as copying each person `count` times. `outcomeData`
# (optional) is the output of .getOutcomeData(weightedClones), so that segments refitted many times are
# collected only once (`weightedClones` is then not used). The variance estimator in
# `fitOutcomeModelArgs` is ignored; its `msmWeights` is used, as in fitOutcomeModel(). Returns a list
# with `estimates` (a tibble with estimand, estimate and logEstimate, in the same order as
# fitOutcomeModel()) and `coefficients`.
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
  if (!"stabilizationNumerator" %in% colnames(weightedClones$segments)) {
    stop("The weighted clones have no stabilizationNumerator: they were created by an earlier version of computeCensoringWeights(). Rerun computeCensoringWeights()")
  }
  # otherScaleWeight is written by computeCensoringWeights() (Phase1Contracts section 12.2); segments
  # made otherwise have only the numerator:
  columns <- c("cloneRowId", "personSeqId", "strategyId", "startTime", "stopTime", "weight", "stabilizationNumerator", "otherScaleWeight", "outcome", "competing")
  segments <- weightedClones$segments |>
    select(dplyr::any_of(columns)) |>
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
    strategies = metaData$strategies,
    sameDayOrdering = metaData$sameDayOrdering,
    cutIntervalLength = metaData$weights$cutIntervalLength,
    stabilized = metaData$weights$stabilized,
    weightsMetaData = metaData$weights
  ))
}

# The data with the segment weights of the Cox and pooled logistic MSMs (`msmWeights`, Phase1Contracts
# sections 11.1 and 12.2): the weights of computeCensoringWeights() when they are already of that kind,
# and otherwise the segments' otherScaleWeight (truncated on its own scale, or formed from the truncated
# weights, as computeWeightsArgs$truncateMsmWeights set). Segments without otherScaleWeight (not made by
# computeCensoringWeights()) get the weights divided (to unstabilize) or multiplied (to stabilize) by
# their stabilization numerators.
.withMsmWeights <- function(data, msmWeights) {
  stabilize <- msmWeights == "stabilized"
  if (stabilize == data$stabilized) {
    return(data)
  }
  if (!is.null(data$segments$otherScaleWeight)) {
    data$segments$weight <- data$segments$otherScaleWeight * data$segments$count
  } else {
    numerator <- data$segments$stabilizationNumerator
    data$segments$weight <- if (stabilize) data$segments$weight * numerator else data$segments$weight / numerator
  }
  data$stabilized <- stabilize
  return(data)
}

# The weights of the Cox or pooled logistic model (msmWeights) per strategy: their scale, the cap on that
# scale (NA when they are not truncated on their own scale), the largest weight used, and the share of
# the strategy's clones whose weight was truncated (Phase1Contracts section 12.2).
.summarizeMsmWeights <- function(data, msmWeights) {
  msmData <- .withMsmWeights(data, msmWeights)
  caps <- data$weightsMetaData$caps
  strategyIds <- sort(unique(data$segments$strategyId))
  sameScale <- (msmWeights == "stabilized") == isTRUE(data$stabilized)
  pick <- function(column) {
    if (is.null(caps[[column]])) {
      return(rep(NA_real_, length(strategyIds)))
    }
    return(as.numeric(caps[[column]][match(strategyIds, caps$strategyId)]))
  }
  return(tibble(
    strategyId = as.integer(strategyIds),
    scale = msmWeights,
    cap = if (sameScale) pick("cap") else pick("otherScaleCap"),
    maxWeight = as.numeric(tapply(msmData$segments$weight, factor(msmData$segments$strategyId, levels = strategyIds), max)),
    truncatedShare = if (sameScale) pick("truncatedShare") else pick("otherScaleTruncatedShare")
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
  # The MSMs use the weights chosen by msmWeights; the Kaplan-Meier and Aalen-Johansen estimates, which
  # are the same with either weight, use the weights of computeCensoringWeights():
  msmData <- .withMsmWeights(data, fitOutcomeModelArgs$msmWeights)
  if (fitOutcomeModelArgs$modelType == "cox") {
    # The hazard ratio of the Cox MSM, and the rows of "km" for the horizons (Phase1Contracts section
    # 10.5), with their IJ standard errors only under the robust variance:
    model <- .fitCoxMsm(msmData, fitOutcomeModelArgs, variance)
    km <- .fitKm(data, horizons, competing, computeVariance = variance == "robust")
    model$estimates <- dplyr::bind_rows(model$estimates, km$estimates)
    return(model)
  }
  return(.fitLogisticMsm(msmData, fitOutcomeModelArgs, competing, variance))
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

# Fits a Cyclops model of an outcome model with `control`. When that fit fails (an error, or a status
# other than SUCCESS, such as MAX_ITERATIONS), it is refitted with Cyclops' default algorithm and
# convergence criterion (the other settings of `control` kept), with a warning. Returns the fit, or NULL
# with a warning when both fail. The covariates `unpenalizedIds` (the target indicator, the time basis
# and their interactions) are excluded from the prior, as CohortMethod excludes the treatment variable.
# `startingCoefficients` (optional) are passed to Cyclops::fitCyclopsModel().
.fitOutcomeCyclopsModel <- function(cyclopsData, prior, control, unpenalizedIds, description, startingCoefficients = NULL) {
  if (prior$priorType != "none") {
    prior$exclude <- unique(c(prior$exclude, unpenalizedIds))
  }
  attempt <- .tryFitCyclopsModel(cyclopsData, prior, control, startingCoefficients)
  if (!is.null(attempt$failure)) {
    defaultControl <- control
    defaultControl[c("algorithm", "convergenceType", "tolerance", "maxIterations")] <- list(NULL)
    retry <- .tryFitCyclopsModel(cyclopsData, prior, defaultControl, startingCoefficients)
    if (is.null(retry$failure)) {
      warning(sprintf("Fitting the %s failed (%s), so it was fitted with Cyclops' default algorithm and convergence criterion", description, attempt$failure))
      attempt <- retry
    }
  }
  for (message in attempt$warnings) {
    warning(message, call. = FALSE)
  }
  if (!is.null(attempt$failure)) {
    warning(sprintf("Fitting the %s failed: %s", description, attempt$failure))
    return(NULL)
  }
  return(attempt$fit)
}

# One Cyclops fit, with its warnings collected (not signalled) and `failure`: NULL, or the error message
# or return flag of a failed fit.
.tryFitCyclopsModel <- function(cyclopsData, prior, control, startingCoefficients = NULL) {
  warnings <- character(0)
  fit <- withCallingHandlers(
    tryCatch(
      Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control, startingCoefficients = startingCoefficients),
      error = function(e) conditionMessage(e)
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  failure <- if (is.character(fit)) fit else if (fit$return_flag != "SUCCESS") fit$return_flag else NULL
  return(list(fit = fit, failure = failure, warnings = warnings))
}

# Fits the weighted counting-process Cox model with the target indicator (covariate 1) as the only
# covariate, using Cyclops cox_counting (Breslow ties), with `subjectId` (one per segment, for
# vcov(type = "robust")) as the Cyclops subjects when given. Returns the fit, or NULL (with a warning)
# when there are no events or the fit fails.
.fitCyclopsCox <- function(segments, y, prior, control, label, subjectId = NULL) {
  if (!any(y == 1)) {
    warning(sprintf("No %s events: cannot fit the %s Cox model", label, label))
    return(NULL)
  }
  outcomes <- data.frame(
    rowId = seq_len(nrow(segments)),
    y = as.numeric(y),
    time = segments$stopTime,
    startTime = segments$startTime,
    weights = segments$weight
  )
  if (!is.null(subjectId)) {
    outcomes$subjectId <- subjectId
  }
  covariates <- data.frame(rowId = which(segments$x == 1), covariateId = 1, covariateValue = 1)
  cyclopsData <- Cyclops::convertToCyclopsData(outcomes, covariates, modelType = "cox_counting", addIntercept = FALSE, quiet = TRUE)
  return(.fitOutcomeCyclopsModel(cyclopsData, prior, control, unpenalizedIds = 1, description = sprintf("%s Cox model", label)))
}

# The hazard ratio row of the weighted Cox MSM of the outcome, with competing events treated as
# censoring, its coefficient, and its (robust or model-based) variance. Risks from this model are not
# reported (see fitOutcomeModel()).
.fitCoxMsm <- function(data, fitOutcomeModelArgs, variance) {
  segments <- data$segments[data$segments$weight > 0, ]
  times <- .eventTimes(segments)
  segments <- .indexSegments(segments, times)

  # Cyclops' robust variance clusters on the subjects of the data, so the persons are the subjects:
  cyclopsRobust <- variance == "robust" && fitOutcomeModelArgs$robustVarianceEngine == "cyclops"
  fit <- .fitCyclopsCox(
    segments, segments$outcome, fitOutcomeModelArgs$prior, fitOutcomeModelArgs$control, "outcome",
    subjectId = if (cyclopsRobust) segments$personIndex else NULL
  )
  beta <- if (is.null(fit)) NA_real_ else as.numeric(stats::coef(fit)[["1"]])
  coefficients <- c(target = beta)
  if (is.na(beta)) {
    return(list(estimates = .ratioRow("hr", NA_real_, NA_real_), coefficients = coefficients, vcov = NULL))
  }
  outcome <- .breslow(segments, segments$outcome == 1, beta, length(times))
  if (cyclopsRobust) {
    vcov <- unname(stats::vcov(fit, type = "robust"))
  } else if (variance == "robust") {
    vcov <- crossprod(cbind(.coxCoefficientInfluence(segments, outcome, length(data$personSeqIds))))
  } else if (variance == "modelBased") {
    vcov <- matrix(1 / outcome$information)
  } else {
    vcov <- NULL
  }
  if (!is.null(vcov)) {
    dimnames(vcov) <- list(names(coefficients), names(coefficients))
  }
  hr <- .ratioRow("hr", beta, if (is.null(vcov)) NA_real_ else sqrt(vcov[1, 1]))
  return(list(estimates = hr, coefficients = coefficients, vcov = vcov))
}

# The result when the outcome model was fitted but the death model failed: the ratio row(s) `rows`, NA
# risks, RD and RR at every horizon, the coefficients with the death-model ones (`deathNames`) NA, and
# the covariance with NA rows and columns for them.
.withoutRisks <- function(rows, horizons, coefficients, vcov, deathNames) {
  risks <- .naEstimates("ratio", horizons, coefficientNames = character(0))$estimates[-1, ]
  allNames <- c(names(coefficients), deathNames)
  if (!is.null(vcov)) {
    fullVcov <- matrix(NA_real_, length(allNames), length(allNames), dimnames = list(allNames, allNames))
    fullVcov[names(coefficients), names(coefficients)] <- vcov[names(coefficients), names(coefficients)]
    vcov <- fullVcov
  }
  coefficients <- c(coefficients, stats::setNames(rep(NA_real_, length(deathNames)), deathNames))
  return(list(estimates = dplyr::bind_rows(c(rows, list(risks))), coefficients = coefficients, vcov = vcov))
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
# (k L, (k + 1) L], counted from riskStart), with the weight that applies on the clone's last day in
# the interval (t = min((k + 1) L, end of follow-up)), the person's frequency weight `count`, and
# indicators of the outcome and of a competing event in the interval. `segments` must hold all segments
# of each clone, sorted by cloneRowId and startTime.
.expandToIntervals <- function(segments, intervalLength) {
  cloneIndex <- match(segments$cloneRowId, unique(segments$cloneRowId))
  last <- !duplicated(cloneIndex, fromLast = TRUE)
  lastSegments <- segments[last, ]
  nIntervals <- ceiling(lastSegments$stopTime / intervalLength)
  rowClone <- rep(seq_len(nrow(lastSegments)), nIntervals)
  interval <- sequence(nIntervals) - 1L
  lastDay <- pmin((interval + 1) * intervalLength, lastSegments$stopTime[rowClone])
  # Find the segment with startTime < lastDay <= stopTime of each row's clone:
  span <- max(lastSegments$stopTime) + 1
  segmentKey <- cloneIndex * span + segments$startTime
  rowKey <- rowClone * span + lastDay
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
# identify are dropped and get coefficient NA (see .identifiedColumns()).
#
# With an indicator time basis, `groups` (one per cell) gives the cells that share a free hazard: the
# interval (without the arm-by-time interaction) or the cell itself (with it). A group with rows but no
# events has a hazard of 0, and one in which every row has the event a hazard of 1: the limits of the
# maximum likelihood fit, whose coefficient would be -Inf or Inf. Their rows carry no information on the
# other coefficients, so they are left out of the fit, and their cells get the fixed hazard
# (`fixedHazard`, NA for the other cells), as the censoring models do. The coefficients that only these
# cells inform are NA (their estimates are infinite), without a warning. Unless `variance` is "none", also
# returns the bread (inverse weighted information) and the score residuals summed by person and cell.
# Returns NULL (with a warning) if the model cannot be fitted.
.fitPooledModel <- function(rows, y, design, prior, control, label, variance, groups = NULL) {
  if (!any(y == 1)) {
    warning(sprintf("No %s events: cannot fit the %s pooled logistic model", label, label))
    return(NULL)
  }
  nCells <- nrow(design)
  events <- .sumBy(rows$weight * y, rows$cell, nCells)
  nonEvents <- .sumBy(rows$weight * (1 - y), rows$cell, nCells)
  fixedHazard <- rep(NA_real_, nCells)
  if (!is.null(groups)) {
    groupEvents <- .sumBy(events, groups, max(groups))[groups]
    groupNonEvents <- .sumBy(nonEvents, groups, max(groups))[groups]
    fixedHazard[groupEvents + groupNonEvents > 0 & groupEvents == 0] <- 0
    fixedHazard[groupEvents > 0 & groupNonEvents == 0] <- 1
  }
  fitted <- events + nonEvents > 0 & is.na(fixedHazard)
  if (!any(fitted)) {
    warning(sprintf("No cell of the %s pooled logistic model has both rows with and without the event: it cannot be fitted", label))
    return(NULL)
  }
  identification <- .identifiedColumns(design, fitted)
  kept <- identification$kept
  # Columns that no cell with rows identifies (columns that only cells with a fixed hazard identify have
  # infinite estimates instead; both are NA):
  notIdentified <- setdiff(seq_len(ncol(design)), .identifiedColumns(design, events + nonEvents > 0)$kept)

  cells <- which(fitted)
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
  # 0 is the intercept:
  fit <- .fitOutcomeCyclopsModel(cyclopsData, prior, control, unpenalizedIds = c(0, kept[-1]), description = sprintf("%s pooled logistic model", label))
  if (is.null(fit)) {
    return(NULL)
  }
  estimates <- stats::coef(fit)
  coefficients <- stats::setNames(rep(NA_real_, ncol(design)), colnames(design))
  coefficients[1] <- estimates[["(Intercept)"]]
  ids <- setdiff(names(estimates), "(Intercept)")
  coefficients[as.integer(ids)] <- estimates[ids]
  if (length(notIdentified) > 0) {
    warning(sprintf(
      "Coefficient(s) %s of the %s pooled logistic model are not identified by the data and are NA",
      paste(colnames(design)[notIdentified], collapse = ", "), label
    ))
  }
  result <- list(
    coefficients = coefficients,
    kept = kept,
    estimable = identification$estimable,
    fixedHazard = fixedHazard,
    design = design[, kept, drop = FALSE]
  )
  if (variance != "none") {
    probability <- stats::plogis(as.vector(result$design %*% coefficients[kept]))
    # The cells with a fixed hazard add nothing to the information or the scores (p (1 - p) = 0, and
    # y - p = 0 on their rows):
    probability[!is.na(fixedHazard)] <- fixedHazard[!is.na(fixedHazard)]
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
# the death model (`layout$index`). Cells with a fixed hazard of 0 or 1 (`layout$fixedHazards`) get it.
# A risk is NA when it needs the hazard of a cell that is not identified.
.pooledRisks <- function(theta, layout, horizonIntervals) {
  hazard <- function(i) {
    result <- stats::plogis(as.vector(layout$designs[[i]] %*% theta[layout$index[[i]]]))
    result[!layout$estimable[[i]]] <- NA
    fixed <- !is.na(layout$fixedHazards[[i]])
    result[fixed] <- layout$fixedHazards[[i]][fixed]
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

# Warnings for a pooled logistic outcome model with a grace period (Phase1Contracts section 11.6). At the
# end of the grace period of a strategy (the end of day G, or of day G - 1 with "deviationFirst", the
# last day of follow-up of the clones that deviate on day G), the clones of an "initiateWithin(G)"
# strategy that have not started are censored, and the hazards of both strategies can change abruptly.
# A row counts its clone as at risk for the whole interval, so on a grid with no interval ending there
# the censored clones dilute the hazard of the rest of that interval; and a smooth time basis cannot
# follow the change. Grace periods that end before riskStart are left out.
.warnPooledGracePeriod <- function(data, timeBasis) {
  ends <- list()
  for (strategyId in seq_along(data$strategies)) {
    strategy <- data$strategies[[strategyId]]
    if (is.null(strategy$graceEnd)) {
      next
    }
    endDay <- strategy$graceEnd - if (identical(data$sameDayOrdering, "deviationFirst")) 1 else 0
    t <- endDay - data$riskStart + 1
    if (t >= 1) {
      ends[[length(ends) + 1]] <- list(description = sprintf("strategy %d (%s) at the end of day %d (t = %d)", strategyId, strategy$label, endDay, t), t = t)
    }
  }
  intervalLength <- timeBasis$intervalLength
  offGrid <- Filter(function(end) end$t %% intervalLength != 0, ends)
  if (length(offGrid) > 0) {
    warning(sprintf(
      "No interval of the pooled logistic outcome model's grid (intervals of %d days counted from riskStart = %d) ends where the grace period ends (where clones are censored or their weights change): %s. A clone counts as at risk for the whole interval it is censored in, which dilutes that interval's hazard and biases the risks. Use an intervalLength that divides t (for example 1): an indicator time basis on such a grid, with armByTimeInteraction = TRUE, reproduces the weighted Kaplan-Meier risks (modelType = \"km\")",
      intervalLength, data$riskStart, paste(vapply(offGrid, function(end) end$description, character(1)), collapse = "; ")
    ), call. = FALSE)
  }
  if (length(ends) > 0 && timeBasis$type != "indicator") {
    warning(sprintf(
      "The pooled logistic outcome model's smooth time basis (%s) cannot follow the change in the hazards at the end of the grace period: %s, so its risks can be biased. Use an indicator time basis (createTimeBasis(type = \"indicator\")) on a grid with an interval ending there, with armByTimeInteraction = TRUE, which reproduces the weighted Kaplan-Meier risks (modelType = \"km\")",
      timeBasis$type, paste(vapply(ends, function(end) end$description, character(1)), collapse = "; ")
    ), call. = FALSE)
  }
  invisible(NULL)
}

.fitLogisticMsm <- function(data, fitOutcomeModelArgs, competing, variance) {
  timeBasis <- fitOutcomeModelArgs$timeBasis
  intervalLength <- timeBasis$intervalLength
  horizons <- fitOutcomeModelArgs$horizons
  prior <- fitOutcomeModelArgs$prior
  control <- fitOutcomeModelArgs$control
  if (isTRUE(prior$useCrossValidation)) {
    stop("Cross-validation is not supported for the pooled logistic outcome model, which is fitted to rows aggregated by strategy and interval. Use a prior without cross-validation")
  }
  # Each row takes the weight of the clone's last segment in its interval, which is the weight of its
  # last day in the interval only when the segments are split on the interval grid:
  if (!isTRUE(data$cutIntervalLength == intervalLength)) {
    stop(sprintf(
      "The pooled logistic outcome model with intervals of %d days needs weighted clones split on that grid, but they are split on %s. Use computeCensoringWeights(outcomeTimeBasis = ...) with the outcome model's timeBasis",
      intervalLength,
      if (is.null(data$cutIntervalLength) || is.na(data$cutIntervalLength)) "no interval grid" else sprintf("intervals of %d days", data$cutIntervalLength)
    ))
  }
  .warnPooledGracePeriod(data, timeBasis)
  rows <- .expandToIntervals(data$segments, intervalLength)
  rows <- rows[rows$weight > 0, ]
  basis <- .timeBasisColumns(timeBasis, rows$interval, rows$count)
  basisNames <- if (timeBasis$type == "linear") "time" else sprintf("time%d", basis$labels)
  # The design depends on a row only through its cell: the comparator's intervals 0, ..., nIntervals - 1
  # are cells 1, ..., nIntervals, and the target's intervals the next nIntervals cells.
  nIntervals <- max(rows$interval) + 1
  cellInterval <- rep(seq_len(nIntervals) - 1, 2)
  design <- .pooledDesign(rep(0:1, each = nIntervals), basis$build(cellInterval), basisNames, fitOutcomeModelArgs$armByTimeInteraction)
  rows$cell <- as.integer(rows$x * nIntervals + rows$interval + 1)
  # With an indicator basis, the cells that share a free hazard (see .fitPooledModel()): each cell with
  # the arm-by-time interaction, each interval without it.
  groups <- NULL
  if (timeBasis$type == "indicator") {
    groups <- if (fitOutcomeModelArgs$armByTimeInteraction) seq_len(2 * nIntervals) else cellInterval + 1
  }

  outcomeRows <- if (competing) rows$death == 0 else rep(TRUE, nrow(rows))
  models <- list(.fitPooledModel(rows[outcomeRows, ], rows$y[outcomeRows], design, prior, control, "outcome", variance, groups))
  # The rows and response column of each model, for Cyclops' robust variance:
  modelRows <- list(list(subset = outcomeRows, response = "y"))
  if (is.null(models[[1]])) {
    return(.naEstimates("or", horizons, coefficientNames = colnames(design)))
  }
  deathFailed <- FALSE
  if (competing && any(rows$death == 1)) {
    deathModel <- .fitPooledModel(rows, rows$death, design, prior, control, "death", variance, groups)
    if (is.null(deathModel)) {
      # The odds ratio comes from the outcome model alone; the risks need the death model too.
      deathFailed <- TRUE
    } else {
      names(deathModel$coefficients) <- paste0("death:", names(deathModel$coefficients))
      models[[2]] <- deathModel
      modelRows[[2]] <- list(subset = rep(TRUE, nrow(rows)), response = "death")
    }
  }
  coefficients <- unlist(lapply(models, function(model) model$coefficients))
  # The identified coefficients, which the variance and the risks are functions of:
  theta <- unlist(lapply(models, function(model) model$coefficients[model$kept]))
  sizes <- vapply(models, function(model) length(model$kept), integer(1))
  layout <- list(
    designs = lapply(models, function(model) model$design),
    estimable = lapply(models, function(model) model$estimable),
    fixedHazards = lapply(models, function(model) model$fixedHazard),
    index = split(seq_along(theta), rep(seq_along(models), sizes))
  )

  thetaVcov <- NULL
  vcov <- NULL
  if (variance != "none") {
    if (variance == "robust" && fitOutcomeModelArgs$robustVarianceEngine == "cyclops") {
      thetaVcov <- .cyclopsJointSandwich(models, rows, modelRows, theta, control)
    } else if (variance == "robust") {
      thetaVcov <- .jointSandwich(models, length(data$personSeqIds))
    } else {
      thetaVcov <- .blockDiagonal(lapply(models, function(model) model$bread))
    }
    vcov <- matrix(NA_real_, length(coefficients), length(coefficients), dimnames = list(names(coefficients), names(coefficients)))
    vcov[names(theta), names(theta)] <- thetaVcov
  }

  estimates <- list(.ratioRow("or", coefficients[["target"]], if (is.null(vcov)) NA_real_ else .standardError(vcov["target", "target"], thetaVcov)))
  if (deathFailed) {
    return(.withoutRisks(estimates, horizons, coefficients, vcov, paste0("death:", colnames(design))))
  }
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
  if (x$modelType == "cox") {
    message(sprintf("Risks: weighted %s", if (identical(x$metaData$competingEvents, "death")) "Aalen-Johansen" else "Kaplan-Meier"))
  }
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

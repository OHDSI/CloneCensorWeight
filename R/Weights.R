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

#' Compute the censoring weights
#'
#' @description
#' Computes the inverse probability of censoring weights of every clone from the censoring models, and
#' splits each clone's follow-up into segments on which the weight is constant at the times that matter
#' for the outcome models.
#'
#' Let `W(d)` be the weight on day `d` (days since index). With `sameDayOrdering = "eventsFirst"`:
#'
#' - `"initiateWithin"` (grace period ending on day `G`): `W(d) = 1` for `d <= G` and `1 / pi(X)` for
#'   `d > G`, where `pi(X)` is the probability of starting by `G` from the binary model.
#' - `"notInitiateWithin"`: `W(d) = 1 / S(d - 1 | X)` for `d <= G` (with `S(-1) = 1`), where
#'   `S(d | X)` is the probability of not initiating on days `0..d` from the initiation hazard model.
#'   For `d > G`, `W(d) = 1 / (1 - pi(X))`, with `pi(X)` from the binary model at `G`, when the models
#'   were fit with `postGraceModel = "binary"` (the default), and `1 / S(G | X)` with
#'   `postGraceModel = "hazard"` (see [createFitCensoringModelsArgs()] for the assumption that makes
#'   the two agree).
#' - `"neverInitiate"`: `W(d) = 1 / S(d - 1 | X)`.
#'
#' With `"deviationFirst"`, `S(d | X)` is used in place of `S(d - 1 | X)`, and the `"initiateWithin"`
#' weight is `1 / pi(X)` from day `G` on (and the `"notInitiateWithin"` weight `1 / (1 - pi(X))` from
#' day `G` on, with `postGraceModel = "binary"`). A clone that deviates on day `d` then does so before
#' any outcome on that day, so it is not at risk of the outcome on day `d`: its segments end on day
#' `d - 1`, and a clone that deviates on its first day of follow-up has no segments. Stabilized
#' weights multiply by the same quantity from the numerator model without covariates (the marginal
#' `pi-bar` of the binary model in place of `pi(X)`, also in `1 - pi(X)`).
#'
#' Each clone's follow-up is split, on the survival time scale `t = day - riskStart + 1`, at the end of
#' day `G` (day `G - 1` with `"deviationFirst"`) for `"initiateWithin"`, and, for the hazard
#' strategies, at every distinct outcome or competing-event time of any clone (of either strategy)
#' that falls inside the clone's window where the weight varies (days `riskStart..G` for
#' `"notInitiateWithin"`, all follow-up for `"neverInitiate"`). A segment `(startTime, stopTime]` gets
#' the weight of its last day, which is exact for estimators that use the weights only at event times
#' (Cox, Kaplan-Meier and Aalen-Johansen).
#'
#' Every clone is also split at the end of day `tau + 1` (or of day `riskStart`, if later), where `tau`
#' is the largest `graceEnd` of the two strategies, so the weight of the segment that contains that day
#' is exactly `W(tau + 1)`, the weight [computeCovariateBalance()] uses. With `outcomeTimeBasis`, every
#' clone is also split at the end of every interval of `outcomeTimeBasis$intervalLength` days (counted
#' from `riskStart`), so the weight of a clone's last segment in an interval is its weight on its last
#' day in that interval, as the pooled logistic outcome model requires. With daily intervals every
#' segment is one day long. These extra splits do not change any weight at an event time, so the Cox,
#' Kaplan-Meier and Aalen-Johansen estimates are the same with and without them.
#'
#' Weights are truncated per strategy: each segment weight above the cap is set to the cap. The cap
#' is the fixed value, or the `upper` percentile (`quantile` type 7) of the clones' maximum weights,
#' taken over the clones under follow-up on day `tau + 1` (or `riskStart`, if later), whose weight
#' varies, with `basis = "varying"`, or over all clones of the strategy with `basis = "allClones"` (see
#' [createWeightTruncation()]).
#'
#' Every segment also stores its stabilization numerator (the quantity from the numerator model that
#' stabilizes its weight), whether or not the weights are stabilized, so that the outcome models can
#' use either weight (see `msmWeights` in [createFitOutcomeModelArgs()]): the other weight is `weight /
#' stabilizationNumerator` (stabilized weights) or `weight * stabilizationNumerator` (unstabilized
#' weights). Truncation applies to `weight`, and the other weight is formed from the truncated
#' `weight`, so the same segments are truncated under both weights: the other weight's cap on a
#' segment is the cap divided (or multiplied) by the segment's numerator. Within a strategy, the
#' clones at risk at an event time share the numerator, so the weighted Kaplan-Meier and
#' Aalen-Johansen estimates are the same with either weight.
#'
#' @param clones             An object of type [CcwClones] as created by `createClones()`.
#' @param censoringModels    An object of class `CcwCensoringModels` as created by
#'                           [fitCensoringModels()] on the same clones (same outcome, same-day
#'                           ordering and competing events).
#' @param computeWeightsArgs An object of type `ComputeWeightsArgs` as created using
#'                           [createComputeWeightsArgs()].
#' @param outcomeTimeBasis   `NULL`, or the object of type `TimeBasis` (as created using
#'                           [createTimeBasis()]) of the pooled logistic outcome model to be fitted to
#'                           these clones: the clones are then also split on its interval grid.
#'
#' @return
#' A new [CcwClones] object, in its own Andromeda (close it when done), with copies of the tables and
#' metaData of `clones` and two more tables. `clones` itself is not changed, so the same clones can be
#' weighted several times (for example with and without stabilization).
#'
#' - `segments`: `segmentId`, `cloneRowId`, `rowId`, `personSeqId`, `strategyId`, `startTime`
#'   (exclusive), `stopTime` (inclusive), `weight` (truncated; stabilized when
#'   `computeWeightsArgs$stabilized`), `stabilizationNumerator` (the stabilized weight is the
#'   unstabilized one times it, see above), `outcome` and `competing` (1 only on the clone's last
#'   segment when it ends with the outcome or death).
#' - `personWeights`: `cloneRowId`, `strategyId`, `maxWeight` and `weightAtEnd`, the largest and the
#'   last segment weight of each clone with segments, before truncation.
#'
#' `metaData$weights` holds whether the weights are stabilized, the truncation rule, `caps` (a tibble
#' with, per strategy, the `cap` used, `NA` without truncation, and `truncatedShare`, the share of the
#' strategy's segments whose weight was truncated), a summary of the untruncated weights,
#' `cutIntervalLength`, the interval length of the `outcomeTimeBasis` grid the clones are split on
#' (`NA` without one), and `censoringModels`, from which [computeWeightDistribution()] evaluates the
#' weight of any day.
#'
#' @example inst/examples/clones.R
#' @examples
#' censoringModels <- fitCensoringModels(
#'   clones = clones,
#'   cloneCensorWeightData = data,
#'   fitCensoringModelsArgs = createFitCensoringModelsArgs(useCovariates = FALSE)
#' )
#' weightedClones <- computeCensoringWeights(
#'   clones = clones,
#'   censoringModels = censoringModels,
#'   computeWeightsArgs = createComputeWeightsArgs(stabilized = FALSE)
#' )
#' attr(weightedClones, "metaData")$weights$untruncatedSummary
#' close(weightedClones)
#' close(clones)
#' close(data)
#'
#' @export
computeCensoringWeights <- function(clones,
                                    censoringModels,
                                    computeWeightsArgs = createComputeWeightsArgs(),
                                    outcomeTimeBasis = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(clones, add = errorMessages)
  checkmate::assertClass(censoringModels, "CcwCensoringModels", add = errorMessages)
  checkmate::assertR6(computeWeightsArgs, "ComputeWeightsArgs", add = errorMessages)
  checkmate::assertR6(outcomeTimeBasis, "TimeBasis", null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  metaData <- attr(clones, "metaData")
  .checkModelsMatchClones(censoringModels, metaData)
  graceEnds <- unlist(lapply(metaData$strategies, function(strategy) strategy$graceEnd))
  # The day after the grace periods, whose weight computeCovariateBalance() reads (NA: no grace period):
  balanceDay <- if (length(graceEnds) == 0) NA_integer_ else as.integer(max(graceEnds) + 1)
  cutIntervalLength <- if (is.null(outcomeTimeBasis)) NA_integer_ else as.integer(outcomeTimeBasis$intervalLength)

  start <- Sys.time()
  # The weights are computed on a copy of the clones, so the input is not changed:
  weightedClones <- .copyCcwClones(clones)
  completed <- FALSE
  on.exit(if (completed) .dropTables(weightedClones, "ccwW") else close(weightedClones))
  .writeModelTables(weightedClones, censoringModels)
  deviationFirst <- metaData$sameDayOrdering == "deviationFirst"
  offset <- if (deviationFirst) 0L else 1L

  # The follow-up of each clone in the outcome risk set (t scale). With deviationFirst, a clone that
  # deviates on its last day does so before any outcome that day, so that day is left out:
  .andromedaExecute(weightedClones, sprintf("
    CREATE TABLE ccwWClones AS
    SELECT cloneRowId,
      rowId,
      personSeqId,
      strategyId,
      riskStart,
      endReason,
      survivalTime - CASE WHEN %s AND deviationDay = endDay THEN 1 ELSE 0 END AS survivalTime
    FROM clones;", if (deviationFirst) "TRUE" else "FALSE"))

  # Cut times (t scale) of each clone, including its end:
  .andromedaExecute(weightedClones, "
    CREATE TABLE ccwWEventTimes AS
    SELECT DISTINCT survivalTime AS eventTime
    FROM ccwWClones
    WHERE endReason IN ('outcome', 'death');")
  .andromedaExecute(weightedClones, "
    CREATE TABLE ccwWCuts AS
    SELECT ccwWClones.cloneRowId,
      ccwWEventTimes.eventTime AS cutTime
    FROM ccwWClones
    INNER JOIN ccwWParams
      ON ccwWClones.strategyId = ccwWParams.strategyId
    INNER JOIN ccwWEventTimes
      ON ccwWEventTimes.eventTime < ccwWClones.survivalTime
        AND ccwWEventTimes.eventTime + ccwWClones.riskStart - 1 <= ccwWParams.windowEnd
    WHERE ccwWParams.kind = 'hazard'
    UNION
    SELECT ccwWClones.cloneRowId,
      ccwWParams.cutDay - ccwWClones.riskStart + 1 AS cutTime
    FROM ccwWClones
    INNER JOIN ccwWParams
      ON ccwWClones.strategyId = ccwWParams.strategyId
    WHERE ccwWParams.kind = 'binary'
      AND ccwWParams.cutDay - ccwWClones.riskStart + 1 > 0
      AND ccwWParams.cutDay - ccwWClones.riskStart + 1 < ccwWClones.survivalTime
    UNION
    SELECT cloneRowId,
      survivalTime AS cutTime
    FROM ccwWClones
    WHERE survivalTime > 0;")
  # The end of the day after the grace periods (or of riskStart, if later):
  if (!is.na(balanceDay)) {
    .andromedaExecute(weightedClones, sprintf("
      INSERT INTO ccwWCuts
      SELECT cloneRowId,
        GREATEST(%1$d, riskStart) - riskStart + 1 AS cutTime
      FROM ccwWClones
      WHERE GREATEST(%1$d, riskStart) - riskStart + 1 < survivalTime
      EXCEPT
      SELECT cloneRowId, cutTime FROM ccwWCuts;", balanceDay))
  }
  # The ends of the intervals of the outcome model's grid:
  if (!is.na(cutIntervalLength)) {
    .andromedaExecute(weightedClones, sprintf("
      INSERT INTO ccwWCuts
      SELECT cloneRowId,
        cutTime
      FROM (
        SELECT cloneRowId,
          CAST(UNNEST(range(%1$d, survivalTime, %1$d)) AS INTEGER) AS cutTime
        FROM ccwWClones
        WHERE survivalTime > %1$d
      )
      EXCEPT
      SELECT cloneRowId, cutTime FROM ccwWCuts;", cutIntervalLength))
  }

  # Segments, with the day whose survival gives the weight of their last day:
  .andromedaExecute(weightedClones, sprintf("
    CREATE TABLE ccwWSegments AS
    SELECT ccwWCuts.cloneRowId,
      ccwWClones.rowId,
      ccwWClones.personSeqId,
      ccwWClones.strategyId,
      COALESCE(LAG(ccwWCuts.cutTime) OVER (PARTITION BY ccwWCuts.cloneRowId ORDER BY ccwWCuts.cutTime), 0) AS startTime,
      ccwWCuts.cutTime AS stopTime,
      CASE WHEN ccwWCuts.cutTime = ccwWClones.survivalTime AND ccwWClones.endReason = 'outcome' THEN 1 ELSE 0 END AS outcome,
      CASE WHEN ccwWCuts.cutTime = ccwWClones.survivalTime AND ccwWClones.endReason = 'death' THEN 1 ELSE 0 END AS competing,
      %s AS evalDay
    FROM ccwWCuts
    INNER JOIN ccwWClones
      ON ccwWCuts.cloneRowId = ccwWClones.cloneRowId
    INNER JOIN ccwWParams
      ON ccwWClones.strategyId = ccwWParams.strategyId;", .evalDaySql("ccwWCuts.cutTime + ccwWClones.riskStart - 1", offset)))

  .checkPredictionsExist(weightedClones, "ccwWSegments")
  .computeSurvival(weightedClones, "ccwWSegments")

  # The numerator is stored whether or not the weights are stabilized, so the outcome models can use
  # either weight (Phase1Contracts section 11.1):
  .andromedaExecute(weightedClones, sprintf("
    CREATE TABLE ccwWWeighted AS
    SELECT ccwWSegments.*,
      %s AS weight,
      ccwWSurvival.numeratorSurvival AS stabilizationNumerator
    FROM ccwWSegments
    INNER JOIN ccwWSurvival
      ON ccwWSegments.strategyId = ccwWSurvival.strategyId
        AND ccwWSegments.rowId = ccwWSurvival.rowId
        AND ccwWSegments.evalDay = ccwWSurvival.evalDay;", .weightSql(computeWeightsArgs$stabilized)))
  badWeights <- .andromedaSql(weightedClones, "
    SELECT COUNT(*) AS n
    FROM ccwWWeighted
    WHERE weight IS NULL OR isinf(weight) OR isnan(weight);")$n
  if (badWeights > 0) {
    stop(sprintf("%d segment(s) have an infinite or missing weight: a clone remains uncensored where the censoring model gives it no chance of doing so (positivity violation)", badWeights))
  }

  .andromedaExecute(weightedClones, "
    CREATE TABLE ccwWPersonWeights AS
    SELECT cloneRowId,
      strategyId,
      MAX(weight) AS maxWeight,
      arg_max(weight, stopTime) AS weightAtEnd
    FROM ccwWWeighted
    GROUP BY cloneRowId,
      strategyId;")

  caps <- .computeCaps(weightedClones, computeWeightsArgs$truncation, balanceDay)
  weightedClones$ccwWCaps <- caps |>
    filter(!is.na(.data$cap))
  truncated <- .andromedaSql(weightedClones, "
    SELECT ccwWWeighted.strategyId,
      AVG(CASE WHEN ccwWCaps.cap IS NOT NULL AND ccwWWeighted.weight > ccwWCaps.cap THEN 1.0 ELSE 0.0 END) AS truncatedShare
    FROM ccwWWeighted
    LEFT JOIN ccwWCaps
      ON ccwWWeighted.strategyId = ccwWCaps.strategyId
    GROUP BY ccwWWeighted.strategyId;")
  caps$truncatedShare <- as.numeric(truncated$truncatedShare[match(caps$strategyId, truncated$strategyId)])
  .andromedaExecute(weightedClones, sprintf("
    CREATE TABLE ccwWFinalSegments AS
    SELECT CAST(ROW_NUMBER() OVER (ORDER BY ccwWWeighted.cloneRowId, ccwWWeighted.startTime) AS INTEGER) AS segmentId,
      ccwWWeighted.cloneRowId,
      ccwWWeighted.rowId,
      ccwWWeighted.personSeqId,
      ccwWWeighted.strategyId,
      CAST(ccwWWeighted.startTime AS DOUBLE) AS startTime,
      CAST(ccwWWeighted.stopTime AS DOUBLE) AS stopTime,
      %s AS weight,
      ccwWWeighted.stabilizationNumerator,
      CAST(ccwWWeighted.outcome AS INTEGER) AS outcome,
      CAST(ccwWWeighted.competing AS INTEGER) AS competing
    FROM ccwWWeighted
    LEFT JOIN ccwWCaps
      ON ccwWWeighted.strategyId = ccwWCaps.strategyId;", .truncatedWeightSql("ccwWWeighted.weight")))

  weightsMetaData <- list(
    stabilized = computeWeightsArgs$stabilized,
    truncation = computeWeightsArgs$truncation$toList(),
    caps = caps,
    untruncatedSummary = .summarizePersonWeights(weightedClones),
    cutIntervalLength = cutIntervalLength,
    # So computeWeightDistribution() can evaluate the weight of any day (Phase1Contracts section 11.7):
    censoringModels = censoringModels
  )
  weightedClones <- .addSegments(weightedClones, weightedClones$ccwWFinalSegments, weightedClones$ccwWPersonWeights, weightsMetaData)
  completed <- TRUE
  delta <- Sys.time() - start
  message("Computing censoring weights took ", signif(delta, 3), " ", attr(delta, "units"))
  return(weightedClones)
}

.checkModelsMatchClones <- function(censoringModels, metaData) {
  modelsMetaData <- attr(censoringModels, "metaData")
  if (!identical(modelsMetaData$sameDayOrdering, metaData$sameDayOrdering)) {
    stop("The censoring models were fit with a different sameDayOrdering than the clones use")
  }
  # Death ends the risk sets differently under competingEvents = "death" and "none" (Phase1Contracts
  # sections 5.3 and 9, amendment 2):
  if (!identical(modelsMetaData$competingEvents, metaData$competingEvents)) {
    stop("The censoring models were fit on clones with a different competingEvents than the clones use")
  }
  # The risk sets depend on the outcome (contract section 5.3), so models are not shared across outcomes:
  if (!identical(as.numeric(modelsMetaData$outcomeId), as.numeric(metaData$outcomeId))) {
    stop("The censoring models were fit on the clones of a different outcome")
  }
  strategies <- metaData$strategies
  if (length(censoringModels) != length(strategies)) {
    stop("The censoring models must have one element per strategy of the clones")
  }
  for (strategyId in seq_along(strategies)) {
    model <- censoringModels[[strategyId]]
    strategy <- strategies[[strategyId]]
    if (model$strategyId != strategyId || model$type != strategy$type || model$treatment != strategy$treatment ||
      (strategy$type != "neverInitiate" && model$horizon != strategy$graceEnd)) {
      stop(sprintf("The censoring model of strategy %d does not match the strategy of the clones", strategyId))
    }
  }
}

.dropTables <- function(andromeda, prefix) {
  if (!Andromeda::isValidAndromeda(andromeda)) {
    return(invisible(NULL))
  }
  tables <- DBI::dbListTables(andromeda)
  for (table in tables[startsWith(tables, prefix)]) {
    DBI::dbRemoveTable(andromeda, table)
  }
  invisible(NULL)
}

# Step function of a cumulative hazard (time, cumHazard) evaluated at times `at` (0 before the first
# event time).
.evaluateCumulativeHazard <- function(baselineHazard, at) {
  idx <- findInterval(at, baselineHazard$time)
  return(c(0, baselineHazard$cumHazard)[idx + 1])
}

# Writes the per-strategy model parameters, per-person predictions, and baseline tables into the
# clones Andromeda.
.writeModelTables <- function(clones, censoringModels) {
  params <- list()
  values <- list()
  coxBaseline <- list()
  intervals <- list()
  postGraceValues <- list()
  interactionValues <- list()
  interactionBasis <- list()
  for (model in censoringModels) {
    strategyId <- as.integer(model$strategyId)
    isLogistic <- model$modelType == "logistic"
    hasPostGrace <- !is.null(model$postGrace)
    params[[strategyId]] <- tibble(
      strategyId = strategyId,
      kind = model$kind,
      modelType = model$modelType,
      horizon = as.integer(model$horizon),
      windowEnd = if (model$kind == "hazard") as.integer(model$horizon) else NA_integer_,
      cutDay = NA_integer_,
      numerator = if (model$kind == "binary") as.numeric(model$numerator$propensityScore) else NA_real_,
      intervalLength = if (isLogistic) as.integer(model$timeBasis$intervalLength) else NA_integer_,
      nIntervals = if (isLogistic) as.integer(model$timeBasis$nIntervals) else NA_integer_,
      postGrace = hasPostGrace,
      postGraceNumerator = if (hasPostGrace) as.numeric(model$postGrace$numerator$propensityScore) else NA_real_
    )
    if (model$kind == "binary") {
      values[[strategyId]] <- tibble(strategyId = strategyId, rowId = model$predictions$rowId, value = model$predictions$propensityScore)
    } else {
      values[[strategyId]] <- tibble(strategyId = strategyId, rowId = model$predictions$rowId, value = model$predictions$linearPredictor)
    }
    if (hasPostGrace) {
      postGraceValues[[strategyId]] <- tibble(strategyId = strategyId, rowId = model$postGrace$predictions$rowId, value = model$postGrace$predictions$propensityScore)
    }
    interactions <- model$interactions
    if (isLogistic && !is.null(interactions) && nrow(interactions$predictions) > 0) {
      interactionValues[[strategyId]] <- tibble(
        strategyId = strategyId,
        rowId = interactions$predictions$rowId,
        columnIndex = interactions$predictions$columnIndex,
        value = interactions$predictions$value
      )
      interactionBasis[[strategyId]] <- tibble(
        strategyId = strategyId,
        k = rep(as.integer(model$intervals$k), ncol(interactions$basis)),
        columnIndex = rep(seq_len(ncol(interactions$basis)), each = nrow(interactions$basis)),
        zValue = as.vector(interactions$basis)
      ) |>
        filter(.data$zValue != 0)
    }
    if (model$modelType == "cox") {
      days <- seq(-1L, model$horizon)
      coxBaseline[[strategyId]] <- tibble(
        strategyId = strategyId,
        day = days,
        cumHazard = .evaluateCumulativeHazard(model$baselineHazard, days + 1),
        numeratorCumHazard = .evaluateCumulativeHazard(model$numerator$baselineHazard, days + 1)
      )
    }
    if (isLogistic) {
      intervals[[strategyId]] <- tibble(
        strategyId = strategyId,
        k = as.integer(model$intervals$k),
        baseline = model$intervals$baseline,
        numeratorBaseline = model$numerator$intervals$baseline
      )
    }
  }
  params <- dplyr::bind_rows(params)
  sameDayOrdering <- attr(censoringModels, "metaData")$sameDayOrdering
  binary <- params$kind == "binary"
  # The last day with W = 1 for initiateWithin:
  params$cutDay[binary] <- params$horizon[binary] - if (sameDayOrdering == "eventsFirst") 0L else 1L
  clones$ccwWParams <- params
  clones$ccwWValues <- dplyr::bind_rows(values)
  clones$ccwWCoxBaseline <- dplyr::bind_rows(c(
    list(tibble(strategyId = integer(0), day = integer(0), cumHazard = numeric(0), numeratorCumHazard = numeric(0))),
    coxBaseline
  ))
  clones$ccwWIntervals <- dplyr::bind_rows(c(
    list(tibble(strategyId = integer(0), k = integer(0), baseline = numeric(0), numeratorBaseline = numeric(0))),
    intervals
  ))
  clones$ccwWPostGrace <- dplyr::bind_rows(c(
    list(tibble(strategyId = integer(0), rowId = integer(0), value = numeric(0))),
    postGraceValues
  ))
  clones$ccwWInteractionValues <- dplyr::bind_rows(c(
    list(tibble(strategyId = integer(0), rowId = integer(0), columnIndex = integer(0), value = numeric(0))),
    interactionValues
  ))
  clones$ccwWInteractionBasis <- dplyr::bind_rows(c(
    list(tibble(strategyId = integer(0), k = integer(0), columnIndex = integer(0), zValue = numeric(0))),
    interactionBasis
  ))
  invisible(NULL)
}

# The day whose survival gives the weight W(d) of the day `day` (an SQL expression of days since index),
# from the strategy's row of ccwWParams: for binary models 1 after the switch and 0 before (see
# .computeSurvival()), and for hazard models d - `offset` (d - 1 with "eventsFirst", d with
# "deviationFirst"), capped at the end of the window in which the weight varies.
.evalDaySql <- function(day, offset) {
  return(sprintf("
      CASE
        WHEN ccwWParams.kind = 'binary' THEN
          CASE WHEN %1$s > ccwWParams.cutDay THEN 1 ELSE 0 END
        ELSE LEAST(%1$s - %2$d, ccwWParams.windowEnd)
      END", day, offset))
}

# The weight from the survival in ccwWSurvival (see .computeSurvival()), stabilized or not.
.weightSql <- function(stabilized) {
  if (stabilized) {
    return("ccwWSurvival.numeratorSurvival / ccwWSurvival.survival")
  }
  return("1.0 / ccwWSurvival.survival")
}

# The weight `weight` (an SQL expression) truncated at the strategy's cap in ccwWCaps (no row: no cap).
.truncatedWeightSql <- function(weight) {
  return(sprintf("CASE WHEN ccwWCaps.cap IS NOT NULL AND %1$s > ccwWCaps.cap THEN ccwWCaps.cap ELSE %1$s END", weight))
}

# Stops when a weight needs the prediction of a person outside the model's risk set (models fit on
# other clones). `table` holds the strategyId, rowId and evalDay of the weights. Persons outside the
# risk set need none: their segments use S(-1) = 1 or, for the binary model, the days up to the switch.
# So the COALESCE in .computeSurvival never changes a weight.
.checkPredictionsExist <- function(clones, table) {
  missing <- .andromedaSql(clones, sprintf("
    SELECT COUNT(*) AS n
    FROM %1$s
    INNER JOIN ccwWParams
      ON %1$s.strategyId = ccwWParams.strategyId
    LEFT JOIN ccwWValues
      ON %1$s.strategyId = ccwWValues.strategyId
        AND %1$s.rowId = ccwWValues.rowId
    WHERE ccwWValues.rowId IS NULL
      AND CASE WHEN ccwWParams.kind = 'binary' THEN %1$s.evalDay = 1 ELSE %1$s.evalDay >= 0 END;", table))$n
  # The weights after G of notInitiateWithin with postGraceModel = "binary" (evalDay = G) need the binary
  # model at G:
  missing <- missing + .andromedaSql(clones, sprintf("
    SELECT COUNT(*) AS n
    FROM %1$s
    INNER JOIN ccwWParams
      ON %1$s.strategyId = ccwWParams.strategyId
    LEFT JOIN ccwWPostGrace
      ON %1$s.strategyId = ccwWPostGrace.strategyId
        AND %1$s.rowId = ccwWPostGrace.rowId
    WHERE ccwWPostGrace.rowId IS NULL
      AND ccwWParams.postGrace
      AND %1$s.evalDay = ccwWParams.windowEnd;", table))$n
  if (missing > 0) {
    stop(sprintf("%d segment(s) need the censoring model's prediction for a person outside its risk set. Were the censoring models fit on these clones?", missing))
  }
  invisible(NULL)
}

# log(1 - expit(eta)) = -log(1 + exp(eta)), computed without overflow.
.logOneMinusExpitSql <- function(eta) {
  return(sprintf("-(CASE WHEN %1$s > 30 THEN %1$s ELSE LN(1 + EXP(%1$s)) END)", eta))
}

# Fills ccwWSurvival (strategyId, rowId, evalDay, survival, numeratorSurvival) for the distinct
# strategyId, rowId and evalDay of `table`: for binary models the probability of starting by G after the
# switch (1 before), and for hazard models S(evalDay | X).
.computeSurvival <- function(clones, table) {
  .andromedaExecute(clones, sprintf("
    CREATE TABLE ccwWNeeded AS
    SELECT DISTINCT %1$s.strategyId,
      %1$s.rowId,
      %1$s.evalDay,
      ccwWParams.modelType,
      CASE
        WHEN ccwWParams.modelType <> 'logistic' THEN NULL
        WHEN %1$s.evalDay >= ccwWParams.horizon THEN ccwWParams.nIntervals - 1
        ELSE CAST(FLOOR((%1$s.evalDay + 1) / ccwWParams.intervalLength) AS INTEGER) - 1
      END AS lastInterval
    FROM %1$s
    INNER JOIN ccwWParams
      ON %1$s.strategyId = ccwWParams.strategyId;", table))
  .andromedaExecute(clones, "
    CREATE TABLE ccwWSurvival (
      strategyId INTEGER,
      rowId INTEGER,
      evalDay INTEGER,
      survival DOUBLE,
      numeratorSurvival DOUBLE
    );")
  .andromedaExecute(clones, "
    INSERT INTO ccwWSurvival
    SELECT ccwWNeeded.strategyId,
      ccwWNeeded.rowId,
      ccwWNeeded.evalDay,
      CASE WHEN ccwWNeeded.evalDay = 1 THEN ccwWValues.value ELSE 1 END,
      CASE WHEN ccwWNeeded.evalDay = 1 THEN ccwWParams.numerator ELSE 1 END
    FROM ccwWNeeded
    INNER JOIN ccwWParams
      ON ccwWNeeded.strategyId = ccwWParams.strategyId
    LEFT JOIN ccwWValues
      ON ccwWNeeded.strategyId = ccwWValues.strategyId
        AND ccwWNeeded.rowId = ccwWValues.rowId
    WHERE ccwWNeeded.modelType = 'binary';")
  .andromedaExecute(clones, "
    INSERT INTO ccwWSurvival
    SELECT ccwWNeeded.strategyId,
      ccwWNeeded.rowId,
      ccwWNeeded.evalDay,
      EXP(-ccwWCoxBaseline.cumHazard * EXP(COALESCE(ccwWValues.value, 0))),
      EXP(-ccwWCoxBaseline.numeratorCumHazard)
    FROM ccwWNeeded
    INNER JOIN ccwWCoxBaseline
      ON ccwWNeeded.strategyId = ccwWCoxBaseline.strategyId
        AND ccwWNeeded.evalDay = ccwWCoxBaseline.day
    LEFT JOIN ccwWValues
      ON ccwWNeeded.strategyId = ccwWValues.strategyId
        AND ccwWNeeded.rowId = ccwWValues.rowId
    WHERE ccwWNeeded.modelType = 'cox';")
  # Logistic: the covariate-by-time interaction part of the linear predictor per person and interval
  # (covariateTimeInteraction), and the cumulative log survival over the intervals, per person, up to
  # the last interval needed.
  .andromedaExecute(clones, "
    CREATE TABLE ccwWInteractionTerms AS
    SELECT ccwWInteractionValues.strategyId,
      ccwWInteractionValues.rowId,
      ccwWInteractionBasis.k,
      SUM(ccwWInteractionValues.value * ccwWInteractionBasis.zValue) AS value
    FROM ccwWInteractionValues
    INNER JOIN ccwWInteractionBasis
      ON ccwWInteractionValues.strategyId = ccwWInteractionBasis.strategyId
        AND ccwWInteractionValues.columnIndex = ccwWInteractionBasis.columnIndex
    GROUP BY ccwWInteractionValues.strategyId,
      ccwWInteractionValues.rowId,
      ccwWInteractionBasis.k;")
  .andromedaExecute(clones, sprintf("
    CREATE TABLE ccwWCumulative AS
    SELECT persons.strategyId,
      persons.rowId,
      ccwWIntervals.k,
      SUM(%s) OVER (PARTITION BY persons.strategyId, persons.rowId ORDER BY ccwWIntervals.k) AS logSurvival
    FROM (
      SELECT strategyId,
        rowId,
        MAX(lastInterval) AS maxInterval
      FROM ccwWNeeded
      WHERE modelType = 'logistic'
      GROUP BY strategyId,
        rowId
    ) persons
    INNER JOIN ccwWIntervals
      ON persons.strategyId = ccwWIntervals.strategyId
        AND ccwWIntervals.k <= persons.maxInterval
    LEFT JOIN ccwWValues
      ON persons.strategyId = ccwWValues.strategyId
        AND persons.rowId = ccwWValues.rowId
    LEFT JOIN ccwWInteractionTerms
      ON persons.strategyId = ccwWInteractionTerms.strategyId
        AND persons.rowId = ccwWInteractionTerms.rowId
        AND ccwWIntervals.k = ccwWInteractionTerms.k;", .logOneMinusExpitSql("(ccwWIntervals.baseline + COALESCE(ccwWValues.value, 0) + COALESCE(ccwWInteractionTerms.value, 0))")))
  .andromedaExecute(clones, sprintf("
    CREATE TABLE ccwWNumeratorCumulative AS
    SELECT strategyId,
      k,
      SUM(%s) OVER (PARTITION BY strategyId ORDER BY k) AS logSurvival
    FROM ccwWIntervals;", .logOneMinusExpitSql("numeratorBaseline")))
  .andromedaExecute(clones, "
    INSERT INTO ccwWSurvival
    SELECT ccwWNeeded.strategyId,
      ccwWNeeded.rowId,
      ccwWNeeded.evalDay,
      CASE WHEN ccwWNeeded.lastInterval < 0 THEN 1 ELSE EXP(ccwWCumulative.logSurvival) END,
      CASE WHEN ccwWNeeded.lastInterval < 0 THEN 1 ELSE EXP(ccwWNumeratorCumulative.logSurvival) END
    FROM ccwWNeeded
    LEFT JOIN ccwWCumulative
      ON ccwWNeeded.strategyId = ccwWCumulative.strategyId
        AND ccwWNeeded.rowId = ccwWCumulative.rowId
        AND ccwWNeeded.lastInterval = ccwWCumulative.k
    LEFT JOIN ccwWNumeratorCumulative
      ON ccwWNeeded.strategyId = ccwWNumeratorCumulative.strategyId
        AND ccwWNeeded.lastInterval = ccwWNumeratorCumulative.k
    WHERE ccwWNeeded.modelType = 'logistic';")
  # notInitiateWithin with postGraceModel = "binary": from the day the weight no longer changes
  # (evalDay = G: day G + 1 with eventsFirst, day G with deviationFirst), 1 - pi(X) from the binary model
  # at G replaces S(G | X), and 1 - pi-bar replaces the numerator's S-bar(G):
  .andromedaExecute(clones, "
    UPDATE ccwWSurvival
    SET survival = 1 - ccwWPostGrace.value,
      numeratorSurvival = 1 - ccwWParams.postGraceNumerator
    FROM ccwWParams, ccwWPostGrace
    WHERE ccwWSurvival.strategyId = ccwWParams.strategyId
      AND ccwWParams.postGrace
      AND ccwWSurvival.evalDay = ccwWParams.windowEnd
      AND ccwWSurvival.strategyId = ccwWPostGrace.strategyId
      AND ccwWSurvival.rowId = ccwWPostGrace.rowId;")
  invisible(NULL)
}

# The truncation cap per strategy (NA: no truncation). A percentile with basis = "varying" is taken
# over the clones under follow-up on day `balanceDay` (or riskStart, if later), falling back to all
# clones without a grace period (balanceDay NA) or for a strategy with no such clone; with
# basis = "allClones", over all clones.
.computeCaps <- function(clones, truncation, balanceDay) {
  strategyIds <- .andromedaSql(clones, "SELECT DISTINCT strategyId FROM ccwWPersonWeights ORDER BY strategyId;")$strategyId
  if (truncation$type == "none") {
    return(tibble(strategyId = as.integer(strategyIds), cap = NA_real_))
  }
  if (truncation$type == "fixed") {
    return(tibble(strategyId = as.integer(strategyIds), cap = as.numeric(truncation$upper)))
  }
  caps <- .andromedaSql(clones, sprintf("
    SELECT strategyId,
      quantile_cont(maxWeight, %.17g) AS cap
    FROM ccwWPersonWeights
    GROUP BY strategyId
    ORDER BY strategyId;", truncation$upper))
  caps <- tibble(strategyId = as.integer(caps$strategyId), cap = as.numeric(caps$cap))
  if (truncation$basis == "varying" && !is.na(balanceDay)) {
    varyingCaps <- .andromedaSql(clones, sprintf("
      SELECT ccwWPersonWeights.strategyId,
        quantile_cont(ccwWPersonWeights.maxWeight, %.17g) AS cap
      FROM ccwWPersonWeights
      INNER JOIN ccwWClones
        ON ccwWPersonWeights.cloneRowId = ccwWClones.cloneRowId
      WHERE GREATEST(%d, ccwWClones.riskStart) - ccwWClones.riskStart + 1 <= ccwWClones.survivalTime
      GROUP BY ccwWPersonWeights.strategyId;", truncation$upper, balanceDay))
    idx <- match(caps$strategyId, varyingCaps$strategyId)
    caps$cap[!is.na(idx)] <- as.numeric(varyingCaps$cap[idx[!is.na(idx)]])
  }
  return(caps)
}

.summarizePersonWeights <- function(clones) {
  summary <- .andromedaSql(clones, "
    SELECT strategyId,
      COUNT(*) AS clones,
      AVG(maxWeight) AS meanMaxWeight,
      MIN(maxWeight) AS minMaxWeight,
      quantile_cont(maxWeight, 0.25) AS p25MaxWeight,
      quantile_cont(maxWeight, 0.5) AS medianMaxWeight,
      quantile_cont(maxWeight, 0.75) AS p75MaxWeight,
      quantile_cont(maxWeight, 0.99) AS p99MaxWeight,
      MAX(maxWeight) AS maxMaxWeight,
      AVG(weightAtEnd) AS meanWeightAtEnd
    FROM ccwWPersonWeights
    GROUP BY strategyId
    ORDER BY strategyId;")
  summary$strategyId <- as.integer(summary$strategyId)
  summary$clones <- as.integer(summary$clones)
  return(summary)
}

# Weight diagnostics -------------------------------------------------------------------------------

#' Compute the distribution of the censoring weights over time
#'
#' @description
#' For each strategy and day (since index), summarizes the weights of the clones under follow-up on
#' that day (those with a segment that contains `t = d - riskStart + 1`). The weight of a clone on day
#' `d` is its weight `W(d)` (see [computeCensoringWeights()]), stabilized or not as the weights of the
#' clones are, and truncated at the strategy's cap. It is evaluated from the censoring models stored in
#' the weighted clones, with the same code that computes the segment weights, so it is exact on every
#' day, also on the days inside a segment, whose weight is that of its last day (for example `W(0) =
#' 1` on day 0, inside the first segment of a `"notInitiateWithin"` clone).
#'
#' @param weightedClones An object of type [CcwClones] as returned by [computeCensoringWeights()].
#' @param days           The days (since index) at which to summarize the weights. If `NULL`, the
#'                       last days of all segments (the event days, the days on which follow-up ends,
#'                       and the days after which [computeCensoringWeights()] splits the clones).
#'
#' @return
#' A tibble with columns `strategyId`, `day`, `clones` (number under follow-up), `sumWeights`,
#' `meanWeight`, `minWeight`, `p25Weight`, `medianWeight`, `p75Weight`, `maxWeight`, and
#' `effectiveSampleSize` (`sum(w)^2 / sum(w^2)`). The attribute `strategyLabels` holds the strategy
#' labels.
#'
#' @example inst/examples/clones.R
#' @example inst/examples/weightedClones.R
#' @examples
#' weightDistribution <- computeWeightDistribution(weightedClones, days = c(0, 7, 8, 30))
#' weightDistribution
#' close(weightedClones)
#' close(clones)
#' close(data)
#'
#' @export
computeWeightDistribution <- function(weightedClones, days = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(weightedClones, weighted = TRUE, add = errorMessages)
  checkmate::assertIntegerish(days, any.missing = FALSE, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  metaData <- attr(weightedClones, "metaData")
  censoringModels <- metaData$weights$censoringModels
  if (is.null(censoringModels)) {
    stop("The weighted clones do not hold the censoring models of their weights: use the weighted clones returned by computeCensoringWeights()")
  }

  on.exit({
    .dropTables(weightedClones, "ccwD")
    .dropTables(weightedClones, "ccwW")
  })
  if (is.null(days)) {
    .andromedaExecute(weightedClones, "
      CREATE TABLE ccwDDays AS
      SELECT DISTINCT CAST(segments.stopTime + clones.riskStart - 1 AS INTEGER) AS day
      FROM segments
      INNER JOIN clones
        ON segments.cloneRowId = clones.cloneRowId;")
  } else {
    weightedClones$ccwDDays <- tibble(day = as.integer(unique(days)))
  }
  # The clones under follow-up on each day, with the day whose survival gives their weight, evaluated
  # as in computeCensoringWeights():
  .writeModelTables(weightedClones, censoringModels)
  offset <- if (metaData$sameDayOrdering == "deviationFirst") 0L else 1L
  .andromedaExecute(weightedClones, sprintf("
    CREATE TABLE ccwDCloneDays AS
    SELECT segments.strategyId,
      segments.rowId,
      ccwDDays.day,
      %s AS evalDay
    FROM segments
    INNER JOIN clones
      ON segments.cloneRowId = clones.cloneRowId
    INNER JOIN ccwDDays
      ON segments.startTime < ccwDDays.day - clones.riskStart + 1
        AND ccwDDays.day - clones.riskStart + 1 <= segments.stopTime
    INNER JOIN ccwWParams
      ON segments.strategyId = ccwWParams.strategyId;", .evalDaySql("ccwDDays.day", offset)))
  .checkPredictionsExist(weightedClones, "ccwDCloneDays")
  .computeSurvival(weightedClones, "ccwDCloneDays")
  weightedClones$ccwWCaps <- metaData$weights$caps |>
    filter(!is.na(.data$cap)) |>
    select("strategyId", "cap")
  distribution <- .andromedaSql(weightedClones, sprintf("
    SELECT strategyId,
      day,
      COUNT(*) AS clones,
      SUM(weight) AS sumWeights,
      AVG(weight) AS meanWeight,
      MIN(weight) AS minWeight,
      quantile_cont(weight, 0.25) AS p25Weight,
      quantile_cont(weight, 0.5) AS medianWeight,
      quantile_cont(weight, 0.75) AS p75Weight,
      MAX(weight) AS maxWeight,
      POWER(SUM(weight), 2) / SUM(POWER(weight, 2)) AS effectiveSampleSize
    FROM (
      SELECT ccwDCloneDays.strategyId,
        ccwDCloneDays.day,
        %s AS weight
      FROM ccwDCloneDays
      INNER JOIN ccwWSurvival
        ON ccwDCloneDays.strategyId = ccwWSurvival.strategyId
          AND ccwDCloneDays.rowId = ccwWSurvival.rowId
          AND ccwDCloneDays.evalDay = ccwWSurvival.evalDay
      LEFT JOIN ccwWCaps
        ON ccwDCloneDays.strategyId = ccwWCaps.strategyId
    ) cloneDayWeights
    GROUP BY strategyId,
      day
    ORDER BY strategyId,
      day;", .truncatedWeightSql(.weightSql(metaData$weights$stabilized))))
  distribution$strategyId <- as.integer(distribution$strategyId)
  distribution$day <- as.integer(distribution$day)
  distribution$clones <- as.integer(distribution$clones)
  strategies <- metaData$strategies
  attr(distribution, "strategyLabels") <- vapply(strategies, function(strategy) as.character(strategy$label), character(1))
  return(distribution)
}

#' Plot the distribution of the censoring weights over time
#'
#' @description
#' Plots, per strategy, the median weight (line), the interquartile range (band), and the maximum
#' weight (dashed line) of the clones under follow-up on each day. Requires the `ggplot2` package.
#'
#' @param weightDistribution A tibble as created by [computeWeightDistribution()].
#' @param labels             The labels of strategies 1 (target) and 2 (comparator). If `NULL`, the
#'                           `strategyLabels` attribute of `weightDistribution` is used.
#' @param fileName           Name of the file where the plot should be saved, for example
#'                           'plot.png'. See [ggplot2::ggsave()] for supported file formats. If
#'                           `NULL`, the plot is not saved.
#'
#' @return
#' A ggplot object. Use the [ggplot2::ggsave()] function to save to file in a different format.
#'
#' @examples
#' weightDistribution <- data.frame(
#'   strategyId = rep(1:2, each = 3),
#'   day = rep(c(0, 7, 30), 2),
#'   p25Weight = c(1, 1, 1.2, 1, 1.1, 1.3),
#'   medianWeight = c(1, 1, 1.5, 1, 1.2, 1.6),
#'   p75Weight = c(1, 1, 2, 1, 1.4, 2.1),
#'   maxWeight = c(1, 1, 6, 1, 3, 7)
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plotWeightDistribution(weightDistribution)
#' }
#'
#' @export
plotWeightDistribution <- function(weightDistribution,
                                   labels = NULL,
                                   fileName = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(weightDistribution, add = errorMessages)
  checkmate::assertNames(colnames(weightDistribution), must.include = c("strategyId", "day", "p25Weight", "medianWeight", "p75Weight", "maxWeight"), add = errorMessages)
  checkmate::assertCharacter(labels, len = 2, any.missing = FALSE, null.ok = TRUE, add = errorMessages)
  checkmate::assertCharacter(fileName, len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required to plot the weight distribution. Please install it.")
  }
  if (is.null(labels)) {
    labels <- attr(weightDistribution, "strategyLabels")
    if (is.null(labels)) {
      labels <- c("Target", "Comparator")
    }
  }

  data <- weightDistribution
  data$strategy <- factor(labels[data$strategyId], levels = labels)
  colors <- c("#2a78d6", "#eb6834")
  names(colors) <- labels
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$day, group = .data$strategy)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$p25Weight, ymax = .data$p75Weight, fill = .data$strategy), alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = .data$medianWeight, color = .data$strategy), linewidth = 0.7) +
    ggplot2::geom_line(ggplot2::aes(y = .data$maxWeight, color = .data$strategy), linewidth = 0.5, linetype = "dashed") +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      x = "Days since index",
      y = "Weight (median, interquartile range, maximum dashed)",
      color = "Strategy",
      fill = "Strategy"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 7, height = 4.5, dpi = 400)
  }
  return(plot)
}

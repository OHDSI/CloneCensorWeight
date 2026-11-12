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
#' - `"notInitiateWithin"`: `W(d) = 1 / S(d - 1 | X)` for `d <= G` (with `S(-1) = 1`) and
#'   `1 / S(G | X)` for `d > G`, where `S(d | X)` is the probability of not initiating on days `0..d`.
#' - `"neverInitiate"`: `W(d) = 1 / S(d - 1 | X)`.
#'
#' With `"deviationFirst"`, `S(d | X)` is used in place of `S(d - 1 | X)`, and the `"initiateWithin"`
#' weight is `1 / pi(X)` from day `G` on. Stabilized weights multiply by the same quantity from the
#' numerator model without covariates.
#'
#' Each clone's follow-up is split, on the survival time scale `t = day - riskStart + 1`, at the end of
#' day `G` (day `G - 1` with `"deviationFirst"`) for `"initiateWithin"`, and, for the hazard
#' strategies, at every distinct outcome or competing-event time of any clone (of either strategy)
#' that falls inside the clone's window where the weight varies (days `riskStart..G` for
#' `"notInitiateWithin"`, all follow-up for `"neverInitiate"`). A segment `(startTime, stopTime]` gets
#' the weight of its last day, which is exact for estimators that use the weights only at event times
#' (Cox, Kaplan-Meier and Aalen-Johansen).
#'
#' Weights are truncated per strategy: the cap is the `upper` percentile (`quantile` type 7) of the
#' clones' maximum weight, or the fixed value.
#'
#' @param clones             An object of type [CcwClones] as created by `createClones()`.
#' @param censoringModels    An object of class `CcwCensoringModels` as created by
#'                           [fitCensoringModels()] on the same clones.
#' @param computeWeightsArgs An object of type `ComputeWeightsArgs` as created using
#'                           [createComputeWeightsArgs()].
#'
#' @return
#' The `clones` object, with two more tables:
#'
#' - `segments`: `segmentId`, `cloneRowId`, `rowId`, `personSeqId`, `strategyId`, `startTime`
#'   (exclusive), `stopTime` (inclusive), `weight` (truncated), `outcome` and `competing` (1 only on the
#'   clone's last segment when it ends with the outcome or death).
#' - `personWeights`: `cloneRowId`, `strategyId`, `maxWeight` and `weightAtEnd`, the largest and the
#'   last segment weight of each clone before truncation.
#'
#' `metaData$weights` holds whether the weights are stabilized, the truncation rule, the cap used per
#' strategy, and a summary of the untruncated weights. The tables are added in place, but the returned
#' object must be used, because only it carries `metaData$weights`.
#'
#' @examples
#' \dontrun{
#' weightedClones <- computeCensoringWeights(
#'   clones = clones,
#'   censoringModels = censoringModels,
#'   computeWeightsArgs = createComputeWeightsArgs()
#' )
#' }
#'
#' @export
computeCensoringWeights <- function(clones,
                                    censoringModels,
                                    computeWeightsArgs = createComputeWeightsArgs()) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(clones, add = errorMessages)
  checkmate::assertClass(censoringModels, "CcwCensoringModels", add = errorMessages)
  checkmate::assertR6(computeWeightsArgs, "ComputeWeightsArgs", add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  metaData <- attr(clones, "metaData")
  .checkModelsMatchClones(censoringModels, metaData)

  start <- Sys.time()
  on.exit(.dropTables(clones, "ccwW"))
  .writeModelTables(clones, censoringModels)
  offset <- if (metaData$sameDayOrdering == "eventsFirst") 1L else 0L

  # Cut times (t scale) of each clone, including its end:
  .andromedaExecute(clones, "
    CREATE TABLE ccwWEventTimes AS
    SELECT DISTINCT survivalTime AS eventTime
    FROM clones
    WHERE endReason IN ('outcome', 'death');")
  .andromedaExecute(clones, "
    CREATE TABLE ccwWCuts AS
    SELECT clones.cloneRowId,
      ccwWEventTimes.eventTime AS cutTime
    FROM clones
    INNER JOIN ccwWParams
      ON clones.strategyId = ccwWParams.strategyId
    INNER JOIN ccwWEventTimes
      ON ccwWEventTimes.eventTime < clones.survivalTime
        AND ccwWEventTimes.eventTime + clones.riskStart - 1 <= ccwWParams.windowEnd
    WHERE ccwWParams.kind = 'hazard'
    UNION
    SELECT clones.cloneRowId,
      ccwWParams.cutDay - clones.riskStart + 1 AS cutTime
    FROM clones
    INNER JOIN ccwWParams
      ON clones.strategyId = ccwWParams.strategyId
    WHERE ccwWParams.kind = 'binary'
      AND ccwWParams.cutDay - clones.riskStart + 1 > 0
      AND ccwWParams.cutDay - clones.riskStart + 1 < clones.survivalTime
    UNION
    SELECT cloneRowId,
      survivalTime AS cutTime
    FROM clones;")

  # Segments, with the day whose survival gives the weight (binary: 1 after the switch, 0 before):
  .andromedaExecute(clones, sprintf("
    CREATE TABLE ccwWSegments AS
    SELECT ccwWCuts.cloneRowId,
      clones.rowId,
      clones.personSeqId,
      clones.strategyId,
      COALESCE(LAG(ccwWCuts.cutTime) OVER (PARTITION BY ccwWCuts.cloneRowId ORDER BY ccwWCuts.cutTime), 0) AS startTime,
      ccwWCuts.cutTime AS stopTime,
      CASE WHEN ccwWCuts.cutTime = clones.survivalTime AND clones.endReason = 'outcome' THEN 1 ELSE 0 END AS outcome,
      CASE WHEN ccwWCuts.cutTime = clones.survivalTime AND clones.endReason = 'death' THEN 1 ELSE 0 END AS competing,
      CASE
        WHEN ccwWParams.kind = 'binary' THEN
          CASE WHEN ccwWCuts.cutTime + clones.riskStart - 1 > ccwWParams.cutDay THEN 1 ELSE 0 END
        ELSE LEAST(ccwWCuts.cutTime + clones.riskStart - 1 - %d, ccwWParams.windowEnd)
      END AS evalDay
    FROM ccwWCuts
    INNER JOIN clones
      ON ccwWCuts.cloneRowId = clones.cloneRowId
    INNER JOIN ccwWParams
      ON clones.strategyId = ccwWParams.strategyId;", offset))

  .computeSurvival(clones)

  weightExpression <- if (computeWeightsArgs$stabilized) {
    "ccwWSurvival.numeratorSurvival / ccwWSurvival.survival"
  } else {
    "1.0 / ccwWSurvival.survival"
  }
  .andromedaExecute(clones, sprintf("
    CREATE TABLE ccwWWeighted AS
    SELECT ccwWSegments.*,
      %s AS weight
    FROM ccwWSegments
    INNER JOIN ccwWSurvival
      ON ccwWSegments.strategyId = ccwWSurvival.strategyId
        AND ccwWSegments.rowId = ccwWSurvival.rowId
        AND ccwWSegments.evalDay = ccwWSurvival.evalDay;", weightExpression))
  badWeights <- .andromedaSql(clones, "
    SELECT COUNT(*) AS n
    FROM ccwWWeighted
    WHERE weight IS NULL OR isinf(weight) OR isnan(weight);")$n
  if (badWeights > 0) {
    stop(sprintf("%d segment(s) have an infinite or missing weight: a clone remains uncensored where the censoring model gives it no chance of doing so (positivity violation)", badWeights))
  }

  .andromedaExecute(clones, "
    CREATE TABLE ccwWPersonWeights AS
    SELECT cloneRowId,
      strategyId,
      MAX(weight) AS maxWeight,
      arg_max(weight, stopTime) AS weightAtEnd
    FROM ccwWWeighted
    GROUP BY cloneRowId,
      strategyId;")

  caps <- .computeCaps(clones, computeWeightsArgs$truncation)
  clones$ccwWCaps <- caps |>
    filter(!is.na(.data$cap))
  .andromedaExecute(clones, "
    CREATE TABLE ccwWFinalSegments AS
    SELECT CAST(ROW_NUMBER() OVER (ORDER BY ccwWWeighted.cloneRowId, ccwWWeighted.startTime) AS INTEGER) AS segmentId,
      ccwWWeighted.cloneRowId,
      ccwWWeighted.rowId,
      ccwWWeighted.personSeqId,
      ccwWWeighted.strategyId,
      CAST(ccwWWeighted.startTime AS DOUBLE) AS startTime,
      CAST(ccwWWeighted.stopTime AS DOUBLE) AS stopTime,
      CASE
        WHEN ccwWCaps.cap IS NOT NULL AND ccwWWeighted.weight > ccwWCaps.cap THEN ccwWCaps.cap
        ELSE ccwWWeighted.weight
      END AS weight,
      CAST(ccwWWeighted.outcome AS INTEGER) AS outcome,
      CAST(ccwWWeighted.competing AS INTEGER) AS competing
    FROM ccwWWeighted
    LEFT JOIN ccwWCaps
      ON ccwWWeighted.strategyId = ccwWCaps.strategyId;")

  weightsMetaData <- list(
    stabilized = computeWeightsArgs$stabilized,
    truncation = computeWeightsArgs$truncation$toList(),
    caps = caps,
    untruncatedSummary = .summarizePersonWeights(clones)
  )
  clones <- .addSegments(clones, clones$ccwWFinalSegments, clones$ccwWPersonWeights, weightsMetaData)
  delta <- Sys.time() - start
  message("Computing censoring weights took ", signif(delta, 3), " ", attr(delta, "units"))
  return(clones)
}

.checkModelsMatchClones <- function(censoringModels, metaData) {
  modelsMetaData <- attr(censoringModels, "metaData")
  if (!identical(modelsMetaData$sameDayOrdering, metaData$sameDayOrdering)) {
    stop("The censoring models were fit with a different sameDayOrdering than the clones use")
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
  for (model in censoringModels) {
    strategyId <- as.integer(model$strategyId)
    isLogistic <- model$modelType == "logistic"
    params[[strategyId]] <- tibble(
      strategyId = strategyId,
      kind = model$kind,
      modelType = model$modelType,
      horizon = as.integer(model$horizon),
      windowEnd = if (model$kind == "hazard") as.integer(model$horizon) else NA_integer_,
      cutDay = NA_integer_,
      numerator = if (model$kind == "binary") as.numeric(model$numerator$propensityScore) else NA_real_,
      intervalLength = if (isLogistic) as.integer(model$timeBasis$intervalLength) else NA_integer_,
      nIntervals = if (isLogistic) as.integer(model$timeBasis$nIntervals) else NA_integer_
    )
    if (model$kind == "binary") {
      values[[strategyId]] <- tibble(strategyId = strategyId, rowId = model$predictions$rowId, value = model$predictions$propensityScore)
    } else {
      values[[strategyId]] <- tibble(strategyId = strategyId, rowId = model$predictions$rowId, value = model$predictions$linearPredictor)
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
  invisible(NULL)
}

# log(1 - expit(eta)) = -log(1 + exp(eta)), computed without overflow.
.logOneMinusExpitSql <- function(eta) {
  return(sprintf("-(CASE WHEN %1$s > 30 THEN %1$s ELSE LN(1 + EXP(%1$s)) END)", eta))
}

# Fills ccwWSurvival (strategyId, rowId, evalDay, survival, numeratorSurvival): for binary models the
# probability of starting by G after the switch (1 before), and for hazard models S(evalDay | X).
.computeSurvival <- function(clones) {
  .andromedaExecute(clones, "
    CREATE TABLE ccwWNeeded AS
    SELECT DISTINCT ccwWSegments.strategyId,
      ccwWSegments.rowId,
      ccwWSegments.evalDay,
      ccwWParams.modelType,
      CASE
        WHEN ccwWParams.modelType <> 'logistic' THEN NULL
        WHEN ccwWSegments.evalDay >= ccwWParams.horizon THEN ccwWParams.nIntervals - 1
        ELSE CAST(FLOOR((ccwWSegments.evalDay + 1) / ccwWParams.intervalLength) AS INTEGER) - 1
      END AS lastInterval
    FROM ccwWSegments
    INNER JOIN ccwWParams
      ON ccwWSegments.strategyId = ccwWParams.strategyId;")
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
  # Logistic: cumulative log survival over the intervals, per person, up to the last interval needed.
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
        AND persons.rowId = ccwWValues.rowId;", .logOneMinusExpitSql("(ccwWIntervals.baseline + COALESCE(ccwWValues.value, 0))")))
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
  invisible(NULL)
}

# The truncation cap per strategy (NA: no truncation).
.computeCaps <- function(clones, truncation) {
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
  return(tibble(strategyId = as.integer(caps$strategyId), cap = as.numeric(caps$cap)))
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
#' that day: the weight of a clone on day `d` is the weight of its segment that contains
#' `t = d - riskStart + 1`. Weights are those used by the outcome models (after truncation).
#'
#' @param weightedClones An object of type [CcwClones] as returned by [computeCensoringWeights()].
#' @param days           The days (since index) at which to summarize the weights. If `NULL`, the
#'                       last days of all segments (the days on which a weight can change or follow-up
#'                       ends).
#'
#' @return
#' A tibble with columns `strategyId`, `day`, `clones` (number under follow-up), `sumWeights`,
#' `meanWeight`, `minWeight`, `p25Weight`, `medianWeight`, `p75Weight`, `maxWeight`, and
#' `effectiveSampleSize` (`sum(w)^2 / sum(w^2)`).
#'
#' @examples
#' \dontrun{
#' weightDistribution <- computeWeightDistribution(weightedClones)
#' }
#'
#' @export
computeWeightDistribution <- function(weightedClones, days = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(weightedClones, weighted = TRUE, add = errorMessages)
  checkmate::assertIntegerish(days, any.missing = FALSE, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  on.exit(.dropTables(weightedClones, "ccwD"))
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
  distribution <- .andromedaSql(weightedClones, "
    SELECT segments.strategyId,
      ccwDDays.day,
      COUNT(*) AS clones,
      SUM(segments.weight) AS sumWeights,
      AVG(segments.weight) AS meanWeight,
      MIN(segments.weight) AS minWeight,
      quantile_cont(segments.weight, 0.25) AS p25Weight,
      quantile_cont(segments.weight, 0.5) AS medianWeight,
      quantile_cont(segments.weight, 0.75) AS p75Weight,
      MAX(segments.weight) AS maxWeight,
      POWER(SUM(segments.weight), 2) / SUM(POWER(segments.weight, 2)) AS effectiveSampleSize
    FROM segments
    INNER JOIN clones
      ON segments.cloneRowId = clones.cloneRowId
    INNER JOIN ccwDDays
      ON segments.startTime < ccwDDays.day - clones.riskStart + 1
        AND ccwDDays.day - clones.riskStart + 1 <= segments.stopTime
    GROUP BY segments.strategyId,
      ccwDDays.day
    ORDER BY segments.strategyId,
      ccwDDays.day;")
  distribution$strategyId <- as.integer(distribution$strategyId)
  distribution$day <- as.integer(distribution$day)
  distribution$clones <- as.integer(distribution$clones)
  return(distribution)
}

#' Plot the distribution of the censoring weights over time
#'
#' @description
#' Plots, per strategy, the median weight (line), the interquartile range (band), and the maximum
#' weight (dashed line) of the clones under follow-up on each day. Requires the `ggplot2` package.
#'
#' @param weightDistribution A tibble as created by [computeWeightDistribution()].
#' @param strategyLabels     Labels of the two strategies (target first).
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
                                   strategyLabels = c("Target", "Comparator"),
                                   fileName = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(weightDistribution, add = errorMessages)
  checkmate::assertNames(colnames(weightDistribution), must.include = c("strategyId", "day", "p25Weight", "medianWeight", "p75Weight", "maxWeight"), add = errorMessages)
  checkmate::assertCharacter(strategyLabels, len = 2, any.missing = FALSE, add = errorMessages)
  checkmate::assertCharacter(fileName, len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required to plot the weight distribution. Please install it.")
  }

  data <- weightDistribution
  data$strategy <- factor(strategyLabels[data$strategyId], levels = strategyLabels)
  colors <- c("#2a78d6", "#eb6834")
  names(colors) <- strategyLabels
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

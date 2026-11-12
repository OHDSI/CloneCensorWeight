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

#' Fit the censoring models
#'
#' @description
#' Fits, for each strategy, the model of artificial censoring (deviation from the strategy) used to
#' compute the inverse probability of censoring weights. Days are relative to the index date (day 0).
#'
#' - `"initiateWithin"` (grace period ending on day `G`): a binary model at `G`. The risk set holds the
#'   persons with no outcome or death on or before day `G` (before day `G` with
#'   `sameDayOrdering = "deviationFirst"`) and `riskEnd >= G`. The response is "started the treatment
#'   on or before day `G`". The model is fit by [CohortMethod::createPs()] on a `CohortMethodData` view
#'   of the risk set, so it gets CohortMethod's covariate tidying, high-correlation check, sampling
#'   and regularization, and the AUC is computed with [CohortMethod::computePsAuc()].
#' - `"notInitiateWithin"` and `"neverInitiate"`: a model of the hazard of initiation. Each person is
#'   at risk of initiating from day 0 until the first of: initiation, the day before an outcome or
#'   death (that day with `"deviationFirst"`), `riskEnd`, and the strategy's horizon `H` (`G`, or the
#'   largest `riskEnd` for `"neverInitiate"`). Days outside this window are never part of the risk
#'   set (structural zeros). With `hazardModelType = "cox"` there is one row per person with
#'   time = last day at risk + 1, fit with a Cyclops Cox model and a Breslow baseline hazard, so the
#'   probability of not initiating on days `0..d` is `S(d | X) = exp(-Lambda0(d + 1) exp(X beta))`.
#'   With `hazardModelType = "logistic"` there is one row per person and interval of
#'   `timeBasis$intervalLength` days (intervals counted from day 0, the last one ending at `H`), fit by
#'   pooled logistic regression with an intercept, the time basis (never penalized) and the
#'   covariates; `S(d | X)` is the product of `1 - h` over the intervals ending on or before day `d`.
#'   With the indicator basis, an interval without initiations gets a hazard of 0, and an interval in
#'   which everyone at risk initiates a hazard of 1 (the limits of the maximum likelihood fit). The rows
#'   of these intervals are left out of the fit, and the first remaining interval is the reference.
#'
#' Hazard models are fit per strategy horizon, and a fit is shared only between strategies with the
#' same type, treatment and horizon. Each model also has a numerator model without covariates, fit on
#' the same rows, used to stabilize the weights: the marginal probability of starting by `G` (binary),
#' the Breslow cumulative hazard of a model without covariates (Cox), or the pooled logistic model
#' with the time basis only. If the pooled logistic numerator model fails and `stopOnError = FALSE`,
#' an intercept-only numerator is used and the error is recorded in the diagnostics.
#'
#' The Cox model treats initiation days as continuous times with Breslow ties, so
#' `exp(-Lambda0(d + 1) exp(X beta))` is not the probability of no initiation on days `0..d` when
#' the daily probability of initiating is large (for example `exp(-0.3) = 0.74` against `1 - 0.3 =
#' 0.7` for one day). The resulting weights do not reproduce the population, and weighted risks can
#' be biased when initiation is confounded. Prefer `hazardModelType = "logistic"` for short grace
#' periods with many initiations per day.
#'
#' With `useCovariates = FALSE` the models have only an intercept and time, so they equal their
#' numerator models.
#'
#' @param clones                 An object of type [CcwClones] as created by `createClones()`.
#' @param cloneCensorWeightData  An object of type [CloneCensorWeightData] holding the baseline
#'                               covariates of the persons in `clones`.
#' @param fitCensoringModelsArgs An object of type `FitCensoringModelsArgs` as created using
#'                               [createFitCensoringModelsArgs()].
#' @param .fixedPriorVariances   For internal use (the bootstrap): `NULL`, or a list named by
#'                               `strategyId` (`"1"`, `"2"`) of prior variances. For a strategy in the
#'                               list, its covariate model is fit with that prior variance and without
#'                               cross-validation.
#'
#' @return
#' An object of class `CcwCensoringModels`: a list with one element per strategy (in strategy order),
#' each a list with
#'
#' - `strategyId`, `type`, `treatment`, `kind` (`"binary"` or `"hazard"`), `modelType` (`"binary"`,
#'   `"cox"` or `"logistic"`), and `horizon` (`G`, or the largest `riskEnd` for `"neverInitiate"`).
#' - `coefficients` (named by covariate ID; `"(Intercept)"` and the time-basis columns for the logistic
#'   model) and `priorVariance` (`NA` without a prior).
#' - `predictions`: a tibble with `rowId` and `propensityScore` (binary) or `linearPredictor` (the
#'   covariate part of the linear predictor, hazard models) for every person in the risk set.
#' - For Cox: `baselineHazard`, a tibble `time, cumHazard` (Breslow, at the initiation times, on the
#'   scale time = day + 1). For logistic: `timeBasis` (the settings, interval length, number of
#'   intervals, the basis matrix, and the intervals with a hazard of 0 or 1) and `intervals`, a tibble
#'   `k, startDay, endDay, baseline` with the intercept plus time-basis part of the linear predictor of
#'   each interval (`-Inf` or `Inf` for a hazard of 0 or 1).
#' - `numerator`: the model without covariates (`propensityScore`; `baselineHazard`; or
#'   `coefficients` and `intervals`).
#' - `diagnostics`: risk-set size, number of rows, events (initiations or starts by `G`), the AUC
#'   (binary models), and the fitting status.
#'
#' The `metaData` attribute holds the outcome ID, same-day ordering, strategies and the settings.
#'
#' @examples
#' \dontrun{
#' censoringModels <- fitCensoringModels(
#'   clones = clones,
#'   cloneCensorWeightData = cloneCensorWeightData,
#'   fitCensoringModelsArgs = createFitCensoringModelsArgs()
#' )
#' }
#'
#' @export
fitCensoringModels <- function(clones,
                               cloneCensorWeightData,
                               fitCensoringModelsArgs = createFitCensoringModelsArgs(),
                               .fixedPriorVariances = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(clones, add = errorMessages)
  .assertCloneCensorWeightData(cloneCensorWeightData, add = errorMessages)
  checkmate::assertR6(fitCensoringModelsArgs, "FitCensoringModelsArgs", add = errorMessages)
  checkmate::assertList(.fixedPriorVariances, types = "numeric", null.ok = TRUE, add = errorMessages)
  if (!is.null(.fixedPriorVariances)) {
    checkmate::assertSubset(names(.fixedPriorVariances), c("1", "2"), empty.ok = FALSE, add = errorMessages)
  }
  checkmate::reportAssertions(collection = errorMessages)

  metaData <- attr(clones, "metaData")
  strategies <- metaData$strategies
  for (strategy in strategies) {
    if (strategy$graceStart > 0 || length(strategy$prohibitedTreatments) > 0) {
      stop("Not yet supported in this version: censoring models for strategies with graceStart > 0 or prohibitedTreatments")
    }
  }
  sameDayOrdering <- metaData$sameDayOrdering
  persons <- clones$persons |>
    collect()
  maxRiskEnd <- if (nrow(persons) == 0) 0L else max(persons$riskEnd)

  models <- list()
  fitted <- list()
  for (strategyId in seq_along(strategies)) {
    strategy <- strategies[[strategyId]]
    if (strategy$type == "initiateWithin") {
      kind <- "binary"
      horizon <- as.integer(strategy$graceEnd)
    } else {
      kind <- "hazard"
      horizon <- if (strategy$type == "neverInitiate") as.integer(maxRiskEnd) else as.integer(strategy$graceEnd)
    }
    fixedVariance <- .fixedPriorVariances[[as.character(strategyId)]]
    key <- paste(strategy$type, strategy$treatment, horizon, if (is.null(fixedVariance)) "cv" else fixedVariance)
    if (key %in% names(fitted)) {
      message(sprintf("Strategy %d shares the censoring model of an earlier strategy", strategyId))
      model <- fitted[[key]]
    } else {
      treatmentDay <- if (strategy$treatment == "treatment") persons$treatmentDay else persons$comparatorTreatmentDay
      prior <- .getCensoringPrior(fitCensoringModelsArgs$prior, fixedVariance)
      if (kind == "binary") {
        message(sprintf("Fitting binary censoring model at day %d for strategy %d (%s)", horizon, strategyId, strategy$label))
        riskSet <- .getBinaryRiskSet(persons, treatmentDay, horizon, sameDayOrdering)
        model <- .fitBinaryCensoringModel(riskSet, cloneCensorWeightData, fitCensoringModelsArgs, prior)
      } else {
        message(sprintf(
          "Fitting %s initiation hazard model through day %d for strategy %d (%s)",
          fitCensoringModelsArgs$hazardModelType, horizon, strategyId, strategy$label
        ))
        riskSet <- .getHazardRiskSet(persons, treatmentDay, horizon, sameDayOrdering)
        if (fitCensoringModelsArgs$hazardModelType == "cox") {
          model <- .fitCoxCensoringModel(riskSet, cloneCensorWeightData, fitCensoringModelsArgs, prior)
        } else {
          model <- .fitLogisticCensoringModel(riskSet, horizon, cloneCensorWeightData, fitCensoringModelsArgs, prior)
        }
      }
      model <- c(list(kind = kind, horizon = horizon), model)
      fitted[[key]] <- model
    }
    models[[strategyId]] <- c(
      list(strategyId = strategyId, type = strategy$type, treatment = strategy$treatment),
      model
    )
  }
  attr(models, "metaData") <- list(
    outcomeId = metaData$outcomeId,
    sameDayOrdering = sameDayOrdering,
    strategies = strategies,
    fitCensoringModelsArgs = fitCensoringModelsArgs$toList()
  )
  class(models) <- "CcwCensoringModels"
  return(models)
}

# Risk sets ----------------------------------------------------------------------------------------

# Persons under follow-up at day G, with y = started on or before day G.
.getBinaryRiskSet <- function(persons, treatmentDay, graceEnd, sameDayOrdering) {
  eventDay <- pmin(persons$outcomeDay, persons$deathDay, na.rm = TRUE)
  if (sameDayOrdering == "eventsFirst") {
    noEvent <- is.na(eventDay) | eventDay > graceEnd
  } else {
    noEvent <- is.na(eventDay) | eventDay >= graceEnd
  }
  inRiskSet <- noEvent & persons$riskEnd >= graceEnd
  return(tibble(
    rowId = persons$rowId[inRiskSet],
    personSeqId = persons$personSeqId[inRiskSet],
    y = as.integer(!is.na(treatmentDay[inRiskSet]) & treatmentDay[inRiskSet] <= graceEnd)
  ))
}

# Persons at risk of initiating on day 0, with the last day at risk (through the horizon) and whether
# they initiated on that day.
.getHazardRiskSet <- function(persons, treatmentDay, horizon, sameDayOrdering) {
  eventDay <- pmin(persons$outcomeDay, persons$deathDay, na.rm = TRUE)
  if (sameDayOrdering == "eventsFirst") {
    eventDay <- eventDay - 1L
  }
  lastDay <- pmin(treatmentDay, eventDay, persons$riskEnd, horizon, na.rm = TRUE)
  initiated <- !is.na(treatmentDay) & treatmentDay == lastDay
  inRiskSet <- lastDay >= 0
  return(tibble(
    rowId = persons$rowId[inRiskSet],
    personSeqId = persons$personSeqId[inRiskSet],
    lastDay = as.integer(lastDay[inRiskSet]),
    y = as.integer(initiated[inRiskSet])
  ))
}

# Priors -------------------------------------------------------------------------------------------

.getCensoringPrior <- function(prior, fixedVariance) {
  if (!is.null(fixedVariance) && prior$priorType != "none") {
    prior$variance <- fixedVariance
    prior$useCrossValidation <- FALSE
  }
  return(prior)
}

# Restricts the prior's excluded covariates to the ones in the model (Cyclops stops on unknown IDs),
# and adds the always-unpenalized ones (intercept and time basis).
.restrictPrior <- function(prior, covariateIds, unpenalizedIds) {
  if (prior$priorType == "none") {
    return(prior)
  }
  prior$exclude <- unique(c(intersect(prior$exclude, covariateIds), unpenalizedIds))
  return(prior)
}

.getPriorVariance <- function(prior, cyclopsFit) {
  if (prior$priorType == "none" || is.null(cyclopsFit)) {
    return(NA_real_)
  }
  return(cyclopsFit$variance[1])
}

# Covariates ---------------------------------------------------------------------------------------

# Returns a CovariateData (Andromeda) with the tidied (normalized, without redundant or infrequent
# covariates) baseline covariates of the risk set, as CohortMethod::createPs does, or NULL when
# covariates are not used. The caller closes it.
.getRiskSetCovariateData <- function(riskSet, cloneCensorWeightData, fitCensoringModelsArgs) {
  if (!fitCensoringModelsArgs$useCovariates || nrow(riskSet) == 0) {
    return(NULL)
  }
  view <- .createCmDataView(
    cloneCensorWeightData,
    cohorts = tibble(rowId = riskSet$rowId, personSeqId = riskSet$personSeqId, treatment = riskSet$y)
  )
  on.exit(close(view))
  covariates <- view$covariates
  if (length(fitCensoringModelsArgs$includeCovariateIds) != 0) {
    covariates <- covariates |>
      filter(.data$covariateId %in% !!fitCensoringModelsArgs$includeCovariateIds)
  }
  if (length(fitCensoringModelsArgs$excludeCovariateIds) != 0) {
    covariates <- covariates |>
      filter(!.data$covariateId %in% !!fitCensoringModelsArgs$excludeCovariateIds)
  }
  covariateData <- Andromeda::andromeda(
    covariates = covariates,
    covariateRef = view$covariateRef,
    analysisRef = view$analysisRef
  )
  attr(covariateData, "metaData") <- list(populationSize = nrow(riskSet))
  class(covariateData) <- "CovariateData"
  attr(class(covariateData), "package") <- "FeatureExtraction"
  tidyCovariateData <- FeatureExtraction::tidyCovariateData(covariateData)
  close(covariateData)
  return(tidyCovariateData)
}

# Stops (or returns an error message, with stopOnError = FALSE) when a covariate is highly correlated
# with initiation at the person level, as CohortMethod::createPs does for the treatment.
.checkHighCorrelation <- function(covariateData, riskSet, fitCensoringModelsArgs) {
  if (!fitCensoringModelsArgs$errorOnHighCorrelation || length(unique(riskSet$y)) < 2) {
    return(NULL)
  }
  covariateData$correlationOutcomes <- riskSet |>
    select("rowId", "y")
  cyclopsData <- Cyclops::convertToCyclopsData(
    covariateData$correlationOutcomes,
    covariateData$covariates,
    modelType = "lr",
    quiet = TRUE
  )
  covariateData$correlationOutcomes <- NULL
  suspect <- Cyclops::getUnivariableCorrelation(cyclopsData, threshold = 0.5)
  suspect <- suspect[!is.na(suspect)]
  if (length(suspect) == 0) {
    return(NULL)
  }
  message("High correlation between covariate(s) and initiation detected: ", paste(names(suspect), collapse = ", "))
  error <- "High correlation between covariate(s) and initiation detected. Perhaps you forgot to exclude part of the treatment definition from the covariates?"
  if (fitCensoringModelsArgs$stopOnError) {
    stop(error)
  }
  return(error)
}

# Fits a Cyclops model, returning the fit, or an error message when stopOnError = FALSE.
.fitCyclops <- function(cyclopsData, prior, control, stopOnError) {
  cyclopsFit <- tryCatch(
    Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control),
    error = function(e) conditionMessage(e)
  )
  if (!is.character(cyclopsFit) && cyclopsFit$return_flag != "SUCCESS") {
    cyclopsFit <- cyclopsFit$return_flag
  }
  if (is.character(cyclopsFit) && stopOnError) {
    stop(cyclopsFit)
  }
  return(cyclopsFit)
}

# The covariate part of the linear predictor for every row of the risk set (0 without covariates).
.computeLinearPredictor <- function(covariateData, riskSet, coefficients) {
  linearPredictor <- rep(0, nrow(riskSet))
  coefficients <- coefficients[!names(coefficients) %in% "(Intercept)"]
  coefficients <- coefficients[coefficients != 0]
  if (is.null(covariateData) || length(coefficients) == 0) {
    return(linearPredictor)
  }
  covariateData$coefficients <- tibble(covariateId = as.numeric(names(coefficients)), beta = as.vector(coefficients))
  sums <- .andromedaSql(covariateData, "
    SELECT covariates.rowId,
      SUM(covariates.covariateValue * coefficients.beta) AS linearPredictor
    FROM covariates
    INNER JOIN coefficients
      ON covariates.covariateId = coefficients.covariateId
    GROUP BY covariates.rowId;")
  covariateData$coefficients <- NULL
  idx <- match(sums$rowId, riskSet$rowId)
  linearPredictor[idx[!is.na(idx)]] <- sums$linearPredictor[!is.na(idx)]
  return(linearPredictor)
}

# Binary model at G --------------------------------------------------------------------------------

.fitBinaryCensoringModel <- function(riskSet, cloneCensorWeightData, fitCensoringModelsArgs, prior) {
  marginal <- if (nrow(riskSet) == 0) NA_real_ else mean(riskSet$y)
  coefficients <- c("(Intercept)" = stats::qlogis(marginal))
  priorVariance <- NA_real_
  error <- "OK"
  auc <- NA_real_
  propensityScore <- rep(marginal, nrow(riskSet))
  hasCovariates <- FALSE
  if (fitCensoringModelsArgs$useCovariates && length(unique(riskSet$y)) == 2) {
    view <- .createCmDataView(
      cloneCensorWeightData,
      cohorts = tibble(rowId = riskSet$rowId, personSeqId = riskSet$personSeqId, treatment = riskSet$y)
    )
    on.exit(close(view))
    hasCovariates <- (view$covariates |> count() |> pull()) > 0
  }
  if (hasCovariates) {
    createPsArgs <- CohortMethod::createCreatePsArgs(
      excludeCovariateIds = fitCensoringModelsArgs$excludeCovariateIds,
      includeCovariateIds = fitCensoringModelsArgs$includeCovariateIds,
      maxCohortSizeForFitting = fitCensoringModelsArgs$maxCohortSizeForFitting,
      errorOnHighCorrelation = fitCensoringModelsArgs$errorOnHighCorrelation,
      stopOnError = fitCensoringModelsArgs$stopOnError,
      prior = prior,
      control = fitCensoringModelsArgs$control
    )
    population <- riskSet |>
      dplyr::rename(treatment = "y")
    ps <- CohortMethod::createPs(cohortMethodData = view, population = population, createPsArgs = createPsArgs)
    psMetaData <- attr(ps, "metaData")
    error <- psMetaData$psError
    if (error == "OK") {
      coefficients <- psMetaData$psModelCoef
      priorVariance <- if (prior$priorType == "none") NA_real_ else psMetaData$psModelPriorVariance
      propensityScore <- ps$propensityScore[match(riskSet$rowId, ps$rowId)]
    } else {
      message("Using the censoring model without covariates: ", error)
    }
  }
  if (length(unique(riskSet$y)) == 2) {
    # computePsAuc is O(starters x non-starters); keep its default down-sampling to 100,000 rows.
    auc <- as.numeric(CohortMethod::computePsAuc(
      tibble(treatment = riskSet$y, propensityScore = propensityScore)
    ))
  }
  return(list(
    modelType = "binary",
    coefficients = coefficients,
    priorVariance = priorVariance,
    predictions = tibble(rowId = riskSet$rowId, propensityScore = propensityScore),
    numerator = list(propensityScore = marginal),
    diagnostics = list(
      riskSetSize = nrow(riskSet),
      rows = nrow(riskSet),
      events = sum(riskSet$y),
      auc = auc,
      status = error
    )
  ))
}

# Cox initiation hazard ----------------------------------------------------------------------------

# Breslow cumulative baseline hazard at the event times, with row weights:
# Lambda0(t) = sum over event times s <= t of (weighted events at s) / (sum over time >= s of w exp(lp)).
.computeBreslowHazard <- function(time, y, linearPredictor = rep(0, length(time)), weights = rep(1, length(time))) {
  if (!any(y == 1 & weights > 0)) {
    return(tibble(time = numeric(0), cumHazard = numeric(0)))
  }
  times <- sort(unique(time))
  group <- match(time, times)
  # Weighted events and risk weights per distinct time; the risk set at s sums the times >= s.
  events <- as.vector(rowsum(weights * (y == 1), group))
  atRisk <- rev(cumsum(rev(as.vector(rowsum(weights * exp(linearPredictor), group)))))
  isEventTime <- events > 0
  return(tibble(time = as.numeric(times[isEventTime]), cumHazard = cumsum(events[isEventTime] / atRisk[isEventTime])))
}

.fitCoxCensoringModel <- function(riskSet, cloneCensorWeightData, fitCensoringModelsArgs, prior) {
  riskSet$time <- riskSet$lastDay + 1
  numeratorHazard <- .computeBreslowHazard(riskSet$time, riskSet$y)
  coefficients <- numeric(0)
  priorVariance <- NA_real_
  error <- "OK"
  linearPredictor <- rep(0, nrow(riskSet))
  covariateData <- NULL
  if (sum(riskSet$y) > 0) {
    covariateData <- .getRiskSetCovariateData(riskSet, cloneCensorWeightData, fitCensoringModelsArgs)
  }
  if (!is.null(covariateData)) {
    on.exit(close(covariateData))
    covariateIds <- covariateData$covariates |>
      dplyr::distinct(.data$covariateId) |>
      pull()
    if (length(covariateIds) > 0) {
      error <- .checkHighCorrelation(covariateData, riskSet, fitCensoringModelsArgs)
      if (is.null(error)) {
        error <- "OK"
        covariateData$outcomes <- riskSet |>
          select("rowId", "y", "time")
        cyclopsData <- Cyclops::convertToCyclopsData(
          covariateData$outcomes,
          covariateData$covariates,
          modelType = "cox",
          quiet = TRUE
        )
        prior <- .restrictPrior(prior, covariateIds, unpenalizedIds = c())
        cyclopsFit <- .fitCyclops(cyclopsData, prior, fitCensoringModelsArgs$control, fitCensoringModelsArgs$stopOnError)
        if (is.character(cyclopsFit)) {
          error <- cyclopsFit
        } else {
          coefficients <- stats::coef(cyclopsFit)
          priorVariance <- .getPriorVariance(prior, cyclopsFit)
          linearPredictor <- .computeLinearPredictor(covariateData, riskSet, coefficients)
        }
      }
    }
  }
  if (error != "OK") {
    message("Using the censoring model without covariates: ", error)
  }
  baselineHazard <- if (length(coefficients) == 0) {
    numeratorHazard
  } else {
    .computeBreslowHazard(riskSet$time, riskSet$y, linearPredictor)
  }
  return(list(
    modelType = "cox",
    coefficients = coefficients,
    priorVariance = priorVariance,
    predictions = tibble(rowId = riskSet$rowId, linearPredictor = linearPredictor),
    baselineHazard = baselineHazard,
    numerator = list(baselineHazard = numeratorHazard),
    diagnostics = list(
      riskSetSize = nrow(riskSet),
      rows = nrow(riskSet),
      events = sum(riskSet$y),
      auc = NA_real_,
      status = error
    )
  ))
}

# Pooled logistic initiation hazard ----------------------------------------------------------------

# The time-basis columns (without intercept) for intervals k = 0..nIntervals - 1. The basis is built
# once on this grid, so its knots do not depend on which rows are in a model. The indicator basis has
# a column for each of `indicatorIntervals` except the first (the reference), named by its interval.
.createTimeBasisMatrix <- function(timeBasis, nIntervals, indicatorIntervals = seq_len(nIntervals) - 1) {
  k <- seq_len(nIntervals) - 1
  if (nIntervals <= 1 || (timeBasis$type == "indicator" && length(indicatorIntervals) <= 1)) {
    return(matrix(numeric(0), nrow = nIntervals, ncol = 0))
  }
  basis <- switch(timeBasis$type,
    indicator = outer(k, indicatorIntervals[-1], "==") * 1,
    linear = matrix(k, ncol = 1),
    polynomial = stats::poly(k, degree = min(timeBasis$degree, nIntervals - 1), raw = timeBasis$raw),
    naturalSpline = splines::ns(k, df = min(timeBasis$df, nIntervals - 1))
  )
  basis <- matrix(as.numeric(basis), nrow = nIntervals)
  if (timeBasis$type == "indicator") {
    colnames(basis) <- paste0("timeBasis", indicatorIntervals[-1])
  } else {
    colnames(basis) <- paste0("timeBasis", seq_len(ncol(basis)))
  }
  return(basis)
}

# Intercept-only logistic MLE on the hazard rows.
.fitInterceptOnly <- function(covariateData) {
  meanY <- .andromedaSql(covariateData, "SELECT AVG(y) AS meanY FROM hazardOutcomes;")$meanY
  return(c("(Intercept)" = stats::qlogis(meanY)))
}

.fitLogisticCensoringModel <- function(riskSet, horizon, cloneCensorWeightData, fitCensoringModelsArgs, prior) {
  intervalLength <- fitCensoringModelsArgs$timeBasis$intervalLength
  nIntervals <- horizon %/% intervalLength + 1L
  k <- seq_len(nIntervals) - 1L
  intervals <- tibble(
    k = k,
    startDay = k * intervalLength,
    endDay = pmin((k + 1L) * intervalLength - 1L, horizon)
  )
  riskSet$lastInterval <- riskSet$lastDay %/% intervalLength
  # With the indicator basis, the hazard of an interval without initiations is 0, and that of an
  # interval in which everyone at risk initiates is 1 (its coefficient is -Inf or Inf at the maximum
  # likelihood). The rows of these intervals carry no information on the other coefficients. So they
  # are left out of the fit, and the first remaining interval is the reference. This is the exact
  # limit of the fit, which would otherwise not converge.
  fitIntervals <- k
  fullIntervals <- integer(0)
  if (fitCensoringModelsArgs$timeBasis$type == "indicator" && sum(riskSet$y) > 0) {
    initiations <- tabulate(riskSet$lastInterval[riskSet$y == 1] + 1L, nbins = nIntervals)
    atRisk <- rev(cumsum(rev(tabulate(riskSet$lastInterval + 1L, nbins = nIntervals))))
    fullIntervals <- k[initiations > 0 & initiations == atRisk]
    fitIntervals <- k[initiations > 0 & initiations < atRisk]
  }
  emptyIntervals <- setdiff(k, c(fitIntervals, fullIntervals))
  basis <- .createTimeBasisMatrix(fitCensoringModelsArgs$timeBasis, nIntervals, fitIntervals)
  basisIds <- -seq_len(ncol(basis))

  covariateData <- .getRiskSetCovariateData(riskSet, cloneCensorWeightData, fitCensoringModelsArgs)
  if (is.null(covariateData)) {
    covariateData <- Andromeda::andromeda(covariates = tibble(rowId = integer(0), covariateId = numeric(0), covariateValue = numeric(0)))
  }
  on.exit(close(covariateData))
  covariateIds <- covariateData$covariates |>
    dplyr::distinct(.data$covariateId) |>
    pull()
  covariateData$riskSet <- riskSet |>
    select("rowId", "lastInterval", "y")
  covariateData$intervals <- tibble(k = fitIntervals)
  covariateData$basis <- tibble(
    k = rep(k, ncol(basis)),
    covariateId = rep(as.numeric(basisIds), each = nIntervals),
    covariateValue = as.vector(basis)
  ) |>
    filter(.data$covariateValue != 0)
  .andromedaExecute(covariateData, "
    CREATE TABLE hazardRows AS
    SELECT CAST(ROW_NUMBER() OVER (ORDER BY riskSet.rowId, intervals.k) AS DOUBLE) AS rowId,
      riskSet.rowId AS personRowId,
      intervals.k,
      CASE WHEN riskSet.y = 1 AND intervals.k = riskSet.lastInterval THEN 1 ELSE 0 END AS y
    FROM riskSet
    INNER JOIN intervals
      ON intervals.k <= riskSet.lastInterval;")
  .andromedaExecute(covariateData, "
    CREATE TABLE hazardOutcomes AS
    SELECT rowId, y
    FROM hazardRows;")
  .andromedaExecute(covariateData, "
    CREATE TABLE hazardBasis AS
    SELECT hazardRows.rowId,
      basis.covariateId,
      basis.covariateValue
    FROM hazardRows
    INNER JOIN basis
      ON hazardRows.k = basis.k;")
  .andromedaExecute(covariateData, "
    CREATE TABLE hazardCovariates AS
    SELECT hazardRows.rowId,
      covariates.covariateId,
      covariates.covariateValue
    FROM hazardRows
    INNER JOIN covariates
      ON hazardRows.personRowId = covariates.rowId
    UNION ALL
    SELECT rowId, covariateId, covariateValue
    FROM hazardBasis;")
  nRows <- covariateData$hazardRows |>
    count() |>
    pull()
  events <- sum(riskSet$y)
  # Basis columns with no rows (intervals nobody reaches) are not in the model:
  presentBasisIds <- covariateData$hazardBasis |>
    dplyr::distinct(.data$covariateId) |>
    pull()

  fitLr <- function(covariates, prior) {
    cyclopsData <- Cyclops::convertToCyclopsData(covariateData$hazardOutcomes, covariates, modelType = "lr", quiet = TRUE)
    return(.fitCyclops(cyclopsData, prior, fitCensoringModelsArgs$control, fitCensoringModelsArgs$stopOnError))
  }

  # Numerator: intercept and time basis.
  numeratorError <- NULL
  if (nRows == 0) {
    # Every interval has a hazard of 0 or 1, so there is nothing to fit. The intercept is the first
    # interval's log-odds.
    numeratorCoefficients <- c("(Intercept)" = if (0L %in% fullIntervals) Inf else -Inf)
  } else if (length(presentBasisIds) == 0 || events == 0) {
    numeratorCoefficients <- .fitInterceptOnly(covariateData)
  } else {
    numeratorFit <- fitLr(covariateData$hazardBasis, Cyclops::createPrior("none"))
    if (is.character(numeratorFit)) {
      numeratorError <- paste("Fitting the numerator model failed:", numeratorFit)
      message(numeratorError, ". Using an intercept-only numerator model")
      numeratorCoefficients <- .fitInterceptOnly(covariateData)
    } else {
      numeratorCoefficients <- stats::coef(numeratorFit)
    }
  }

  # Denominator: intercept, time basis and covariates.
  coefficients <- numeratorCoefficients
  priorVariance <- NA_real_
  error <- "OK"
  linearPredictor <- rep(0, nrow(riskSet))
  if (length(covariateIds) > 0 && events > 0 && nRows > 0) {
    error <- .checkHighCorrelation(covariateData, riskSet, fitCensoringModelsArgs)
    if (is.null(error)) {
      error <- "OK"
      prior <- .restrictPrior(prior, covariateIds, unpenalizedIds = c(0, presentBasisIds))
      cyclopsFit <- fitLr(covariateData$hazardCovariates, prior)
      if (is.character(cyclopsFit)) {
        error <- cyclopsFit
      } else {
        coefficients <- stats::coef(cyclopsFit)
        priorVariance <- .getPriorVariance(prior, cyclopsFit)
        covariateCoefficients <- coefficients[!names(coefficients) %in% c("(Intercept)", as.character(basisIds))]
        linearPredictor <- .computeLinearPredictor(covariateData, riskSet, covariateCoefficients)
      }
    }
    if (error != "OK") {
      message("Using the censoring model without covariates: ", error)
    }
  }
  return(list(
    modelType = "logistic",
    coefficients = .renameBasisCoefficients(coefficients, basisIds, colnames(basis)),
    priorVariance = priorVariance,
    predictions = tibble(rowId = riskSet$rowId, linearPredictor = linearPredictor),
    timeBasis = list(
      settings = fitCensoringModelsArgs$timeBasis$toList(),
      intervalLength = intervalLength,
      nIntervals = nIntervals,
      basis = basis,
      emptyIntervals = emptyIntervals,
      fullIntervals = fullIntervals
    ),
    intervals = dplyr::mutate(intervals, baseline = .computeBaseline(coefficients, basis, basisIds, emptyIntervals, fullIntervals)),
    numerator = list(
      coefficients = .renameBasisCoefficients(numeratorCoefficients, basisIds, colnames(basis)),
      intervals = dplyr::mutate(intervals, baseline = .computeBaseline(numeratorCoefficients, basis, basisIds, emptyIntervals, fullIntervals))
    ),
    diagnostics = list(
      riskSetSize = nrow(riskSet),
      rows = as.integer(nRows),
      events = events,
      auc = NA_real_,
      status = paste(c(numeratorError, if (error != "OK" || is.null(numeratorError)) error), collapse = "; ")
    )
  ))
}

# The intercept plus time-basis part of the linear predictor for each interval (-Inf or Inf, a hazard
# of 0 or 1, for intervals left out of an indicator fit).
.computeBaseline <- function(coefficients, basis, basisIds, emptyIntervals, fullIntervals) {
  baseline <- rep(coefficients[["(Intercept)"]], nrow(basis))
  if (ncol(basis) > 0) {
    gamma <- coefficients[as.character(basisIds)]
    gamma[is.na(gamma)] <- 0
    baseline <- baseline + as.vector(basis %*% gamma)
  }
  baseline[emptyIntervals + 1] <- -Inf
  baseline[fullIntervals + 1] <- Inf
  return(baseline)
}

# Names the time-basis coefficients and orders them: intercept, time basis, covariates.
.renameBasisCoefficients <- function(coefficients, basisIds, basisNames) {
  idx <- match(names(coefficients), as.character(basisIds))
  names(coefficients)[!is.na(idx)] <- basisNames[idx[!is.na(idx)]]
  order <- c("(Intercept)", basisNames[basisNames %in% names(coefficients)])
  return(coefficients[c(order, setdiff(names(coefficients), order))])
}

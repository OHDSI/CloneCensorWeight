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
#'   With `covariateTimeInteraction = "day0"` or `"interval"`, the pooled logistic model also has the
#'   products of the covariates with an indicator of the first interval, or with each time-basis
#'   column, penalized like the covariates (see [createFitCensoringModelsArgs()]).
#' - `"notInitiateWithin"` with `postGraceModel = "binary"` (the default): also the binary model at `G`
#'   of an `"initiateWithin(G)"` strategy of the same treatment (shared with such a strategy when the
#'   analysis has one), which gives the weights after `G` (see [createFitCensoringModelsArgs()]).
#'
#' Hazard models are fit per strategy horizon, and a fit is shared only between strategies with the
#' same type, treatment and horizon. Each model also has a numerator model without covariates, fit on
#' the same rows, used to stabilize the weights: the marginal probability of starting by `G` (binary),
#' the Breslow cumulative hazard of a model without covariates (Cox), or the pooled logistic model
#' with the time basis only. If the pooled logistic numerator model fails and `stopOnError = FALSE`,
#' an intercept-only numerator is used and the error is recorded in the diagnostics.
#'
#' The default, `hazardModelType = "logistic"` with a daily indicator basis, is the exact
#' discrete-time hazard of initiation. The Cox model treats initiation days as continuous times with
#' Breslow ties, so `exp(-Lambda0(d + 1) exp(X beta))` is not the probability of no initiation on days
#' `0..d` when the daily probability of initiating is large (for example `exp(-0.3) = 0.74` against
#' `1 - 0.3 = 0.7` for one day, as when many persons start on day 0). The error builds up over the
#' days of the horizon, so the bias of the Cox weights grows with the horizon and with each person's
#' daily probability of initiating, and is largest for `"neverInitiate"` strategies. The resulting
#' weights do not reproduce the population, and weighted risks can be biased when initiation is
#' confounded. For long horizons, prefer the pooled logistic model, with a longer `intervalLength` if
#' daily intervals are too costly.
#'
#' `maxCohortSizeForFitting` applies only to the binary models (through [CohortMethod::createPs()]);
#' in this version the hazard models are always fitted to their whole risk set.
#'
#' With `useCovariates = FALSE` the models have only an intercept and time, so they equal their
#' numerator models.
#'
#' Records of day 0 can follow the treatment decision of that day (for example the diagnoses of the
#' visit at which a drug is prescribed). So when a strategy allows initiation on day 0 (an
#' `"initiateWithin"` strategy with `graceStart = 0`), the function warns (once) if the models use
#' covariates whose windows end on day 0 or later: those of analyses with `endDay >= 0` in
#' `cloneCensorWeightData$analysisRef`, as with FeatureExtraction's default `endDays = 0`. Extract the
#' covariates with `endDays = -1` instead.
#'
#' @param clones                 An object of type [CcwClones] as created by `createClones()`.
#' @param cloneCensorWeightData  An object of type [CloneCensorWeightData] holding the baseline
#'                               covariates of the persons in `clones`.
#' @param fitCensoringModelsArgs An object of type `FitCensoringModelsArgs` as created using
#'                               [createFitCensoringModelsArgs()].
#' @param .fixedPriorVariances   For internal use (the bootstrap): `NULL`, or a list named by
#'                               `strategyId` (`"1"`, `"2"`) of prior variances, and by
#'                               `"<strategyId>.postGrace"` for the binary model of a
#'                               `"notInitiateWithin"` strategy's weights after the grace period. For a
#'                               model in the list, its covariate model is fit with that prior variance
#'                               and without cross-validation.
#'
#' @return
#' An object of class `CcwCensoringModels`: a list with one element per strategy (in strategy order),
#' each a list with
#'
#' - `strategyId`, `type`, `treatment`, `kind` (`"binary"` or `"hazard"`), `modelType` (`"binary"`,
#'   `"cox"` or `"logistic"`), and `horizon` (`G`, or the largest `riskEnd` for `"neverInitiate"`).
#' - `coefficients` (named by covariate ID; `"(Intercept)"`, the time-basis columns and the
#'   interactions, named `"<covariateId>:<column>"`, for the logistic model) and `priorVariance` (`NA`
#'   without a prior).
#' - `predictions`: a tibble with `rowId` and `propensityScore` (binary) or `linearPredictor` (the
#'   covariate part of the linear predictor, hazard models) for every person in the risk set.
#' - For Cox: `baselineHazard`, a tibble `time, cumHazard` (Breslow, at the initiation times, on the
#'   scale time = day + 1). For logistic: `timeBasis` (the settings, interval length, number of
#'   intervals, the basis matrix, and the intervals with a hazard of 0 or 1), `intervals`, a tibble
#'   `k, startDay, endDay, baseline` with the intercept plus time-basis part of the linear predictor of
#'   each interval (`-Inf` or `Inf` for a hazard of 0 or 1), and `interactions`: the
#'   `covariateTimeInteraction` setting (`type`), the interaction `basis` (a matrix with one row per
#'   interval and one column per interaction column, which multiplies every covariate), and
#'   `predictions`, a tibble `rowId, columnIndex, value` with the covariate part of each interaction
#'   column's coefficients per person (non-zero values only). The linear predictor of a person in
#'   interval `k` is `baseline + linearPredictor + sum(basis[k, ] * value)`.
#' - `numerator`: the model without covariates (`propensityScore`; `baselineHazard`; or
#'   `coefficients` and `intervals`).
#' - `diagnostics`: risk-set size, number of rows, events (initiations or starts by `G`), the AUC
#'   (binary models), and the fitting status.
#' - `postGrace` (`"notInitiateWithin"` strategies with `postGraceModel = "binary"` only): the binary
#'   model at `G` whose `1 - propensityScore` gives the weights after `G`, a list with the elements of
#'   a binary model above (`kind`, `horizon`, `modelType`, `coefficients`, `priorVariance`,
#'   `predictions`, `numerator`, `diagnostics`).
#'
#' The `metaData` attribute holds the outcome ID, same-day ordering, competing events, strategies and
#' the settings.
#'
#' @example inst/examples/clones.R
#' @examples
#' # Without the baseline covariates, to keep the example fast (by default the models use all
#' # covariates, with a cross-validated prior):
#' censoringModels <- fitCensoringModels(
#'   clones = clones,
#'   cloneCensorWeightData = data,
#'   fitCensoringModelsArgs = createFitCensoringModelsArgs(useCovariates = FALSE)
#' )
#' censoringModels[[2]]$diagnostics
#' close(clones)
#' close(data)
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
    checkmate::assertSubset(names(.fixedPriorVariances), c("1", "2", "1.postGrace", "2.postGrace"), empty.ok = FALSE, add = errorMessages)
  }
  checkmate::reportAssertions(collection = errorMessages)

  metaData <- attr(clones, "metaData")
  strategies <- metaData$strategies
  for (strategy in strategies) {
    if (strategy$graceStart > 0 || length(strategy$prohibitedTreatments) > 0) {
      stop("Not yet supported in this version: censoring models for strategies with graceStart > 0 or prohibitedTreatments")
    }
  }
  .warnDay0CovariateWindows(strategies, cloneCensorWeightData, fitCensoringModelsArgs)
  sameDayOrdering <- metaData$sameDayOrdering
  persons <- clones$persons |>
    collect()
  maxRiskEnd <- if (nrow(persons) == 0) 0L else max(persons$riskEnd)

  models <- list()
  fitted <- list()
  # Fits a model, or reuses the fit of an earlier strategy with the same key. The binary model at G has
  # the key of an initiateWithin(G) strategy, so the binary model that gives a notInitiateWithin(G)
  # strategy's weights after G is shared with an initiateWithin(G) strategy of the same treatment: their
  # risk sets are the same (both use the analysis' sameDayOrdering).
  fitModel <- function(kind, treatment, horizon, fixedVariance, description) {
    treatmentDay <- if (treatment == "treatment") persons$treatmentDay else persons$comparatorTreatmentDay
    prior <- .getCensoringPrior(fitCensoringModelsArgs$prior, fixedVariance)
    message("Fitting ", description)
    # Errors (for example CohortMethod::createPs()'s high-correlation check, which speaks of the
    # "treatment") name the model that failed:
    model <- tryCatch(
      {
        if (kind == "binary") {
          riskSet <- .getBinaryRiskSet(persons, treatmentDay, horizon, sameDayOrdering)
          .fitBinaryCensoringModel(riskSet, cloneCensorWeightData, fitCensoringModelsArgs, prior)
        } else {
          riskSet <- .getHazardRiskSet(persons, treatmentDay, horizon, sameDayOrdering)
          if (fitCensoringModelsArgs$hazardModelType == "cox") {
            .fitCoxCensoringModel(riskSet, cloneCensorWeightData, fitCensoringModelsArgs, prior)
          } else {
            .fitLogisticCensoringModel(riskSet, horizon, cloneCensorWeightData, fitCensoringModelsArgs, prior)
          }
        }
      },
      error = function(e) stop(sprintf("Fitting the %s failed: %s", description, conditionMessage(e)), call. = FALSE)
    )
    return(c(list(kind = kind, horizon = horizon), model))
  }
  varianceKey <- function(fixedVariance) if (is.null(fixedVariance)) "cv" else fixedVariance
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
    key <- paste(strategy$type, strategy$treatment, horizon, varianceKey(fixedVariance))
    if (kind == "binary") {
      description <- sprintf("binary censoring model at day %d for strategy %d (%s)", horizon, strategyId, strategy$label)
    } else {
      description <- sprintf(
        "%s initiation hazard model through day %d for strategy %d (%s)",
        fitCensoringModelsArgs$hazardModelType, horizon, strategyId, strategy$label
      )
    }
    if (key %in% names(fitted)) {
      message(sprintf("Strategy %d shares the censoring model of an earlier strategy", strategyId))
    } else {
      fitted[[key]] <- fitModel(kind, strategy$treatment, horizon, fixedVariance, description)
    }
    models[[strategyId]] <- c(
      list(strategyId = strategyId, type = strategy$type, treatment = strategy$treatment),
      fitted[[key]]
    )
    if (strategy$type == "notInitiateWithin" && fitCensoringModelsArgs$postGraceModel == "binary") {
      fixedVariance <- .fixedPriorVariances[[paste0(strategyId, ".postGrace")]]
      key <- paste("initiateWithin", strategy$treatment, horizon, varianceKey(fixedVariance))
      if (key %in% names(fitted)) {
        message(sprintf("Strategy %d shares the binary model at day %d of an earlier strategy for its weights after the grace period", strategyId, horizon))
      } else {
        description <- sprintf(
          "binary censoring model at day %d for the weights after the grace period of strategy %d (%s)",
          horizon, strategyId, strategy$label
        )
        fitted[[key]] <- fitModel("binary", strategy$treatment, horizon, fixedVariance, description)
      }
      models[[strategyId]]$postGrace <- fitted[[key]]
    }
  }
  attr(models, "metaData") <- list(
    outcomeId = metaData$outcomeId,
    sameDayOrdering = sameDayOrdering,
    competingEvents = metaData$competingEvents,
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

# Warns, once, when a strategy allows initiation on day 0 ("initiateWithin" with graceStart = 0) while
# the censoring models use covariates whose window ends on day 0 or later (`endDay >= 0` in analysisRef,
# as with FeatureExtraction's default endDays = 0): records of day 0 may follow the treatment decision
# of that day (Phase1Contracts section 11.6). Analyses without a window (endDay NA, such as demographics)
# are not a concern.
.warnDay0CovariateWindows <- function(strategies, cloneCensorWeightData, fitCensoringModelsArgs) {
  day0 <- vapply(strategies, function(strategy) strategy$type == "initiateWithin" && strategy$graceStart == 0, logical(1))
  if (!fitCensoringModelsArgs$useCovariates || !any(day0)) {
    return(invisible(NULL))
  }
  covariateRef <- cloneCensorWeightData$covariateRef
  if (length(fitCensoringModelsArgs$includeCovariateIds) != 0) {
    covariateRef <- covariateRef |>
      filter(.data$covariateId %in% !!fitCensoringModelsArgs$includeCovariateIds)
  }
  if (length(fitCensoringModelsArgs$excludeCovariateIds) != 0) {
    covariateRef <- covariateRef |>
      filter(!.data$covariateId %in% !!fitCensoringModelsArgs$excludeCovariateIds)
  }
  analysisIds <- covariateRef |>
    dplyr::distinct(.data$analysisId) |>
    pull()
  if (length(analysisIds) == 0) {
    return(invisible(NULL))
  }
  windows <- cloneCensorWeightData$analysisRef |>
    filter(.data$analysisId %in% !!analysisIds, !is.na(.data$endDay), .data$endDay >= 0) |>
    dplyr::arrange(.data$analysisId) |>
    collect()
  if (nrow(windows) == 0) {
    return(invisible(NULL))
  }
  ids <- which(day0)
  labels <- vapply(strategies[ids], function(strategy) strategy$label, character(1))
  analyses <- sprintf("%s (endDay = %d)", windows$analysisName, as.integer(windows$endDay))
  if (length(analyses) > 5) {
    analyses <- c(analyses[1:5], sprintf("%d more", length(analyses) - 5))
  }
  warning(sprintf(
    "%s %s allow%s initiation on day 0, but the covariate windows of analyses %s end on day 0 or later, so the censoring models can use records of day 0 that may follow the treatment decision. Extract the covariates with windows that end before day 0, for example with endDays = -1 in FeatureExtraction::createCovariateSettings()",
    if (length(ids) == 1) "Strategy" else "Strategies",
    paste(sprintf("%d (%s)", ids, labels), collapse = " and "),
    if (length(ids) == 1) "s" else "",
    paste(analyses, collapse = ", ")
  ), call. = FALSE)
  invisible(NULL)
}

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

# The time-basis columns (without intercept) for intervals k = 0..nIntervals - 1, from the shared
# .timeBasisColumns(). The basis is built once on this grid (each interval counted once), so its knots
# do not depend on which rows are in a model. The indicator basis has a column for each of
# `indicatorIntervals` except the first (the reference), named by its interval.
.createTimeBasisMatrix <- function(timeBasis, nIntervals, indicatorIntervals = seq_len(nIntervals) - 1) {
  k <- seq_len(nIntervals) - 1
  basis <- .timeBasisColumns(timeBasis, if (timeBasis$type == "indicator") indicatorIntervals else k)
  columns <- basis$build(k)
  colnames(columns) <- sprintf("timeBasis%d", basis$labels)
  return(columns)
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
  # Covariate-by-time interactions: column c of `interactionBasis` multiplies every covariate.
  interactionBasis <- .createInteractionBasis(
    fitCensoringModelsArgs$covariateTimeInteraction,
    basis,
    k,
    rowIntervals = fitIntervals[fitIntervals <= max(c(-1L, riskSet$lastInterval))]
  )
  interactionRef <- tibble(
    covariateId = rep(as.numeric(covariateIds), ncol(interactionBasis)),
    columnIndex = rep(seq_len(ncol(interactionBasis)), each = length(covariateIds))
  )
  interactionRef$interactionId <- -(ncol(basis) + seq_len(nrow(interactionRef)))
  covariateData$interactionRef <- interactionRef
  covariateData$interactionBasis <- tibble(
    k = rep(k, ncol(interactionBasis)),
    columnIndex = rep(seq_len(ncol(interactionBasis)), each = nIntervals),
    zValue = as.vector(interactionBasis)
  ) |>
    filter(.data$zValue != 0)
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
    FROM hazardBasis
    UNION ALL
    SELECT hazardRows.rowId,
      interactionRef.interactionId AS covariateId,
      covariates.covariateValue * interactionBasis.zValue AS covariateValue
    FROM hazardRows
    INNER JOIN interactionBasis
      ON hazardRows.k = interactionBasis.k
    INNER JOIN covariates
      ON hazardRows.personRowId = covariates.rowId
    INNER JOIN interactionRef
      ON covariates.covariateId = interactionRef.covariateId
        AND interactionBasis.columnIndex = interactionRef.columnIndex;")
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
  interactionPredictions <- tibble(rowId = integer(0), columnIndex = integer(0), value = numeric(0))
  if (length(covariateIds) > 0 && events > 0 && nRows > 0) {
    error <- .checkHighCorrelation(covariateData, riskSet, fitCensoringModelsArgs)
    if (is.null(error)) {
      error <- "OK"
      # The interaction columns are not excluded, so they are penalized like the covariates:
      prior <- .restrictPrior(prior, covariateIds, unpenalizedIds = c(0, presentBasisIds))
      cyclopsFit <- fitLr(covariateData$hazardCovariates, prior)
      if (is.character(cyclopsFit)) {
        error <- cyclopsFit
      } else {
        coefficients <- stats::coef(cyclopsFit)
        priorVariance <- .getPriorVariance(prior, cyclopsFit)
        isInteraction <- names(coefficients) %in% as.character(interactionRef$interactionId)
        covariateCoefficients <- coefficients[!names(coefficients) %in% c("(Intercept)", as.character(basisIds)) & !isInteraction]
        linearPredictor <- .computeLinearPredictor(covariateData, riskSet, covariateCoefficients)
        interactionPredictions <- .computeInteractionPredictions(covariateData, interactionRef, coefficients[isInteraction])
      }
    }
    if (error != "OK") {
      message("Using the censoring model without covariates: ", error)
    }
  }
  # Interaction coefficients are named covariateId:column, for example "1001:day0":
  idx <- match(names(coefficients), as.character(interactionRef$interactionId))
  names(coefficients)[!is.na(idx)] <- sprintf(
    "%s:%s",
    as.character(interactionRef$covariateId[idx[!is.na(idx)]]),
    colnames(interactionBasis)[interactionRef$columnIndex[idx[!is.na(idx)]]]
  )
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
    interactions = list(
      type = fitCensoringModelsArgs$covariateTimeInteraction,
      basis = interactionBasis,
      predictions = interactionPredictions
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

# The interaction basis (intervals x columns) whose columns multiply every covariate: none, an indicator
# of the first interval ("day0"), or the time-basis columns ("interval"). A column that takes a single
# value over `rowIntervals`, the intervals with rows in the fit, would duplicate the covariates, so it is
# left out (for example "day0" when the first interval is not fit, or nobody reaches a later one).
.createInteractionBasis <- function(covariateTimeInteraction, basis, k, rowIntervals) {
  if (covariateTimeInteraction == "day0") {
    interactionBasis <- matrix(as.numeric(k == 0), ncol = 1, dimnames = list(NULL, "day0"))
  } else if (covariateTimeInteraction == "interval") {
    interactionBasis <- basis
  } else {
    interactionBasis <- matrix(0, nrow = length(k), ncol = 0)
  }
  varies <- vapply(seq_len(ncol(interactionBasis)), function(column) {
    length(unique(interactionBasis[rowIntervals + 1L, column])) > 1
  }, logical(1))
  return(interactionBasis[, varies, drop = FALSE])
}

# The covariate part of the interaction coefficients, per person and interaction column with a non-zero
# value (rowId, columnIndex, value): in interval k, the linear predictor of a person gains the sum over
# columns c of interactionBasis[k, c] * value_c.
.computeInteractionPredictions <- function(covariateData, interactionRef, coefficients) {
  coefficients <- coefficients[coefficients != 0]
  if (length(coefficients) == 0) {
    return(tibble(rowId = integer(0), columnIndex = integer(0), value = numeric(0)))
  }
  idx <- match(names(coefficients), as.character(interactionRef$interactionId))
  covariateData$interactionCoefficients <- tibble(
    covariateId = interactionRef$covariateId[idx],
    columnIndex = interactionRef$columnIndex[idx],
    beta = as.vector(coefficients)
  )
  predictions <- .andromedaSql(covariateData, "
    SELECT covariates.rowId,
      interactionCoefficients.columnIndex,
      SUM(covariates.covariateValue * interactionCoefficients.beta) AS value
    FROM covariates
    INNER JOIN interactionCoefficients
      ON covariates.covariateId = interactionCoefficients.covariateId
    GROUP BY covariates.rowId,
      interactionCoefficients.columnIndex
    ORDER BY covariates.rowId,
      interactionCoefficients.columnIndex;")
  covariateData$interactionCoefficients <- NULL
  return(tibble(rowId = as.integer(predictions$rowId), columnIndex = as.integer(predictions$columnIndex), value = predictions$value))
}

# S(day | X) of a hazard model, the probability of no initiation on days 0..day (day <= horizon), for
# the persons `rowIds` of its risk set.
.computeHazardSurvival <- function(model, rowIds, day) {
  person <- match(rowIds, model$predictions$rowId)
  linearPredictor <- model$predictions$linearPredictor[person]
  if (model$modelType == "cox") {
    return(exp(-.evaluateCumulativeHazard(model$baselineHazard, day + 1) * exp(linearPredictor)))
  }
  complete <- which(model$intervals$endDay <= day)
  eta <- outer(linearPredictor, model$intervals$baseline[complete], "+")
  interactions <- model$interactions
  if (!is.null(interactions) && nrow(interactions$predictions) > 0) {
    values <- matrix(0, nrow = nrow(model$predictions), ncol = ncol(interactions$basis))
    values[cbind(match(interactions$predictions$rowId, model$predictions$rowId), interactions$predictions$columnIndex)] <- interactions$predictions$value
    eta <- eta + values[person, , drop = FALSE] %*% t(interactions$basis[complete, , drop = FALSE])
  }
  # log(1 - expit(eta)) = -log(1 + exp(eta)):
  logSurvival <- -ifelse(eta > 30, eta, log1p(exp(eta)))
  return(exp(rowSums(matrix(logSurvival, nrow = length(rowIds)))))
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

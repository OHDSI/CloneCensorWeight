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
#' - `"initiateWithin"` (grace period ending on day `G`), with `initiateWithinModel = "binary"` (the
#'   default): a binary model at `G`. The risk set holds the
#'   persons with no outcome or death on or before day `G` (before day `G` with
#'   `sameDayOrdering = "deviationFirst"`) and `riskEnd >= G`. The response is "started the treatment
#'   on or before day `G`". The model is fit by [CohortMethod::createPs()] on a `CohortMethodData` view
#'   of the risk set, so it gets CohortMethod's covariate tidying, high-correlation check, sampling
#'   and regularization, and the AUC is computed with [CohortMethod::computePsAuc()].
#' - `"initiateWithin"` with `initiateWithinModel = "hazard"`: the model of the hazard of initiation
#'   through `G` described next, the same model as that of a `"notInitiateWithin(G)"` strategy of the
#'   same treatment (shared with such a strategy when the analysis has one), whose `1 - S(G | X)`
#'   gives the weights after `G` (see [createFitCensoringModelsArgs()]).
#' - `"initiateWithin"` with `graceStart = S > 0`: a start before day `S` is a deviation. The binary
#'   model at `G` is fit on the persons who have not started before day `S`, with the response
#'   "started on days `S..G`", and the weights of days up to `S` come from the initiation hazard model
#'   through day `S - 1` (element `early`; the model of a `"notInitiateWithin(S - 1)"` strategy, shared
#'   with one in the analysis). With `initiateWithinModel = "hazard"` the hazard model through `G`
#'   gives both.
#' - `"initiateWithin"` with `graceRegime = "naturalGrace"`: the binary model at `G` is fit on the
#'   persons who reach day `G` untreated (no start before day `G`), with the response "started on day
#'   `G`", the probability of the last-interval weight (with `initiateWithinModel = "hazard"`, the
#'   hazard model through `G` gives it). Only the clones that start on day `G` represent the persons
#'   who reach day `G` untreated after `G`, so when none of these persons starts on day `G` the
#'   regime is not identified: those persons are all censored at `G`, and the risks after `G` stand
#'   for the persons who started before `G` only. `fitCensoringModels()` then warns. A covariate
#'   stratum without starts on day `G` shows in [computePositivity()].
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
#'   With `covariateTimeInteraction = "day0"` (the default) or `"interval"`, the pooled logistic model
#'   also has the products of the covariates with an indicator of day 0 (of the first interval with
#'   longer intervals), or with each time-basis column, penalized like the covariates (see
#'   [createFitCensoringModelsArgs()]). An interaction column that takes a single value over the
#'   intervals in the fit would duplicate the covariates, so it is left out, with a message: for
#'   example `"day0"` with 7-day intervals and `G = 7` when nobody starts on day 7, so the first
#'   interval is the only one fit.
#' - `"notInitiateWithin"` with `postGraceModel = "binary"` (not the default, `"hazard"`): also the binary model at `G`
#'   of an `"initiateWithin(G)"` strategy of the same treatment (shared with such a strategy when the
#'   analysis has one), which gives the weights after `G` (see [createFitCensoringModelsArgs()]).
#'
#' Hazard models are fit per strategy horizon, and a fit is shared only between strategies with the
#' same treatment and horizon whose risk sets are the same: `"notInitiateWithin"` strategies, and
#' `"initiateWithin"` strategies with `initiateWithinModel = "hazard"`, with the same `G`, or
#' `"neverInitiate"` strategies. Each model also has a numerator model without covariates, fit on
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
#' `maxCohortSizeForFitting` (`0`: no sampling) applies to every model. The binary models sample
#' their starters and non-starters to at most that many each (through [CohortMethod::createPs()]). A
#' hazard model whose risk set holds more persons than that is fit on a uniform sample of that many
#' persons, with all their rows, drawn with the seed of `control` (1 when it has none) without
#' changing the caller's random number stream: its numerator and denominator models, the
#' cross-validation of the prior and the high-correlation check use the sample, and the linear
#' predictor is then computed for every person of the risk set. A uniform sample of persons needs no
#' correction of the intercept. With the indicator time basis, an interval in which no sampled person,
#' or every one, at risk initiates, but not so in the whole risk set, gets the observed hazard of the
#' whole risk set (shifted by the mean linear predictor of its persons in the model with covariates),
#' so that persons outside the sample keep finite weights. A sample with no initiations, from a risk
#' set with some, is not used: the model is then fit on the whole risk set. The diagnostics of each
#' model record whether it was sampled.
#'
#' With `useCovariates = FALSE` the models have only an intercept and time, so they equal their
#' numerator models.
#'
#' Records of day 0 can follow the treatment decision of that day (for example the diagnoses of the
#' visit at which a drug is prescribed). Every censoring model predicts initiation on day 0: the binary
#' model of an `"initiateWithin"` strategy (with `graceStart = 0`) the probability of starting on days
#' `0..G`, and the hazard models of `"notInitiateWithin"` and `"neverInitiate"` strategies the hazard
#' of starting from day 0 on. So the function warns (once, naming the strategies) if the models use
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
#'                               `"notInitiateWithin"` strategy's weights after the grace period, and
#'                               `"<strategyId>.early"` for the hazard model of the days before the
#'                               `graceStart` of an `"initiateWithin"` strategy. For a model in the list, its covariate model is fit with that prior variance
#'                               and without cross-validation.
#' @param .sourcePersonSeqIds    For internal use (the bootstrap): `NULL`, or a data frame with the
#'                               `personSeqId` of each person and the `sourcePersonSeqId` of the person
#'                               it is a copy of. Copies of one person stay in one cross-validation fold
#'                               of every censoring model.
#' @param .subsampleSeed         For internal use (the bootstrap): `NULL`, or the seed of the samples
#'                               of the hazard models' risk sets (see `maxCohortSizeForFitting`), in
#'                               place of the seed of `control`.
#'
#' @return
#' An object of class `CcwCensoringModels`: a list with one element per strategy (in strategy order),
#' each a list with
#'
#' - `strategyId`, `type`, `treatment`, `kind` (`"binary"` or `"hazard"`; `"hazard"` for an
#'   `"initiateWithin"` strategy with `initiateWithinModel = "hazard"`), `modelType` (`"binary"`,
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
#' - `diagnostics`: risk-set size, number of rows in the fit, events (initiations or starts by `G`),
#'   the AUC (binary models), the fitting status, and the sample the model was fit on (see
#'   `maxCohortSizeForFitting`): `sampled` (`TRUE` when it was fit on a sample of the risk set),
#'   `sampleSize` (persons) and `sampleEvents`. Pooled logistic models also record the `engine` that fit
#'   them (`pooledLogisticEngine`).
#' - `postGrace` (`"notInitiateWithin"` strategies with `postGraceModel = "binary"` only): the binary
#'   model at `G` whose `1 - propensityScore` gives the weights after `G`, a list with the elements of
#'   a binary model above (`kind`, `horizon`, `modelType`, `coefficients`, `priorVariance`,
#'   `predictions`, `numerator`, `diagnostics`).
#' - `"initiateWithin"` strategies only: `graceStart` and `graceRegime`; `firstDay` (binary models fit
#'   on the persons who have not started before a day `firstDay > 0`: `graceStart`, or `G` with
#'   `graceRegime = "naturalGrace"`); and `early` (with `graceStart > 0` and
#'   `initiateWithinModel = "binary"`), the initiation hazard model through `graceStart - 1`, a list
#'   with the elements of a hazard model above.
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
                               .fixedPriorVariances = NULL,
                               .sourcePersonSeqIds = NULL,
                               .subsampleSeed = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(clones, add = errorMessages)
  .assertCloneCensorWeightData(cloneCensorWeightData, add = errorMessages)
  checkmate::assertR6(fitCensoringModelsArgs, "FitCensoringModelsArgs", add = errorMessages)
  checkmate::assertList(.fixedPriorVariances, types = "numeric", null.ok = TRUE, add = errorMessages)
  if (!is.null(.fixedPriorVariances)) {
    checkmate::assertSubset(names(.fixedPriorVariances), c("1", "2", "1.postGrace", "2.postGrace", "1.early", "2.early"), empty.ok = FALSE, add = errorMessages)
  }
  checkmate::assertDataFrame(.sourcePersonSeqIds, null.ok = TRUE, add = errorMessages)
  checkmate::assertInt(.subsampleSeed, null.ok = TRUE, add = errorMessages)
  if (!is.null(.sourcePersonSeqIds)) {
    checkmate::assertNames(names(.sourcePersonSeqIds), must.include = c("personSeqId", "sourcePersonSeqId"), add = errorMessages)
  }
  checkmate::reportAssertions(collection = errorMessages)

  return(.fitCensoringModelsWithFits(
    clones = clones,
    cloneCensorWeightData = cloneCensorWeightData,
    fitCensoringModelsArgs = fitCensoringModelsArgs,
    fixedPriorVariances = .fixedPriorVariances,
    sourcePersonSeqIds = .sourcePersonSeqIds,
    subsampleSeed = .subsampleSeed
  )$models)
}

# fitCensoringModels() with a cache of fits: `fits` is a list of fitted models named by
# .censoringFitKey() (from earlier calls, possibly with other design choices, or from
# .censoringFitsFromModels()), and a model whose key is in it is reused, not refit. Returns a list with
# the CcwCensoringModels `models` and the `fits`, with the new ones added. Callers that fit several
# variants of the same analysis warn about day-0 covariate windows, and about natural-grace regimes
# without starts on day G, once.
.fitCensoringModelsWithFits <- function(clones,
                                        cloneCensorWeightData,
                                        fitCensoringModelsArgs,
                                        fixedPriorVariances = NULL,
                                        sourcePersonSeqIds = NULL,
                                        subsampleSeed = NULL,
                                        fits = list(),
                                        warnDay0CovariateWindows = TRUE,
                                        warnNaturalGrace = TRUE) {
  metaData <- attr(clones, "metaData")
  strategies <- metaData$strategies
  for (strategy in strategies) {
    if (length(strategy$prohibitedTreatments) > 0) {
      stop("Not yet supported in this version: censoring models for strategies with prohibitedTreatments")
    }
  }
  .checkInitiateWithinHazardGrid(strategies, fitCensoringModelsArgs)
  # Before any fit, so that a long run does not stop late:
  if (fitCensoringModelsArgs$hazardModelType == "logistic" && fitCensoringModelsArgs$pooledLogisticEngine == "plr" && !.cyclopsHasPooledLogistic()) {
    stop(sprintf("pooledLogisticEngine = \"plr\" needs a Cyclops with pooled logistic regression (modelType = \"plr\"), which the loaded Cyclops %s does not have. Install one that has it, or use pooledLogisticEngine = \"lr\" (the default), which gives the same fit with more memory", utils::packageVersion("Cyclops")))
  }
  if (warnDay0CovariateWindows) {
    .warnDay0CovariateWindows(strategies, cloneCensorWeightData, fitCensoringModelsArgs)
  }
  sameDayOrdering <- metaData$sameDayOrdering
  persons <- clones$persons |>
    collect()
  # The unit of the censoring models' cross-validation folds: the person, or in a bootstrap replicate
  # the person a copy was drawn from, so that the rows of one person never score each other.
  persons$subjectId <- persons$personSeqId
  if (!is.null(sourcePersonSeqIds)) {
    persons$subjectId <- sourcePersonSeqIds$sourcePersonSeqId[match(persons$personSeqId, sourcePersonSeqIds$personSeqId)]
  }
  if (warnNaturalGrace && fitCensoringModelsArgs$graceRegime == "naturalGrace") {
    .warnNoNaturalGraceStarts(strategies, persons, sameDayOrdering)
  }
  maxRiskEnd <- if (nrow(persons) == 0) 0L else max(persons$riskEnd)
  argsList <- fitCensoringModelsArgs$toList()
  if (is.null(subsampleSeed)) {
    subsampleSeed <- .defaultSubsampleSeed(fitCensoringModelsArgs)
  }

  models <- list()
  # The keys of the fits used by this call: a key used by an earlier strategy is shared with it.
  used <- character(0)
  # Fits a model, or reuses the fit of an earlier strategy, or from `fits`, with the same key. The
  # binary model at G that gives a notInitiateWithin(G) strategy's weights after G is shared with an
  # initiateWithin(G) strategy of the same treatment: their risk sets are the same (both use the
  # analysis' sameDayOrdering). So is the hazard model through G of a notInitiateWithin(G) strategy and
  # of an initiateWithin(G) strategy with initiateWithinModel = "hazard".
  getFit <- function(component, kind, treatment, horizon, fixedVariance, description, shareMessage, firstDay = 0L) {
    key <- .censoringFitKey(component, treatment, horizon, fixedVariance, argsList)
    if (key %in% used) {
      message(shareMessage)
    } else if (key %in% names(fits)) {
      message("Reusing the fitted ", description)
    } else {
      fits[[key]] <<- fitModel(kind, treatment, horizon, fixedVariance, description, firstDay)
    }
    used <<- c(used, key)
    return(fits[[key]])
  }
  fitModel <- function(kind, treatment, horizon, fixedVariance, description, firstDay) {
    treatmentDay <- if (treatment == "treatment") persons$treatmentDay else persons$comparatorTreatmentDay
    prior <- .getCensoringPrior(fitCensoringModelsArgs$prior, fixedVariance)
    message("Fitting ", description)
    # Errors (for example CohortMethod::createPs()'s high-correlation check, which speaks of the
    # "treatment") name the model that failed:
    model <- tryCatch(
      {
        if (kind == "binary") {
          riskSet <- .getBinaryRiskSet(persons, treatmentDay, horizon, sameDayOrdering, firstDay)
          .fitBinaryCensoringModel(riskSet, cloneCensorWeightData, fitCensoringModelsArgs, prior)
        } else {
          riskSet <- .getHazardRiskSet(persons, treatmentDay, horizon, sameDayOrdering)
          fitSet <- .sampleRiskSet(riskSet, fitCensoringModelsArgs$maxCohortSizeForFitting, subsampleSeed, description)
          if (fitCensoringModelsArgs$hazardModelType == "cox") {
            .fitCoxCensoringModel(riskSet, cloneCensorWeightData, fitCensoringModelsArgs, prior, fitSet)
          } else {
            .fitLogisticCensoringModel(riskSet, horizon, cloneCensorWeightData, fitCensoringModelsArgs, prior, description, fitSet)
          }
        }
      },
      error = function(e) stop(sprintf("Fitting the %s failed: %s", description, conditionMessage(e)), call. = FALSE)
    )
    # The binary model of the persons who have not started before day firstDay > 0 (graceStart, or G
    # with graceRegime = "naturalGrace") records that day:
    return(c(list(kind = kind, horizon = horizon), if (firstDay > 0) list(firstDay = as.integer(firstDay)), model))
  }
  for (strategyId in seq_along(strategies)) {
    strategy <- strategies[[strategyId]]
    horizon <- if (strategy$type == "neverInitiate") as.integer(maxRiskEnd) else as.integer(strategy$graceEnd)
    kind <- if (strategy$type == "initiateWithin" && fitCensoringModelsArgs$initiateWithinModel == "binary") "binary" else "hazard"
    firstDay <- if (kind == "binary") .binaryModelFirstDay(strategy, fitCensoringModelsArgs$graceRegime) else 0L
    if (kind == "binary") {
      description <- sprintf(
        "binary censoring model at day %d%s for strategy %d (%s)",
        horizon, if (firstDay > 0) sprintf(" of the persons who have not started before day %d", firstDay) else "", strategyId, strategy$label
      )
    } else {
      description <- sprintf(
        "%s initiation hazard model through day %d for strategy %d (%s)",
        fitCensoringModelsArgs$hazardModelType, horizon, strategyId, strategy$label
      )
    }
    models[[strategyId]] <- c(
      list(strategyId = strategyId, type = strategy$type, treatment = strategy$treatment),
      getFit(
        component = .censoringFitComponent(kind, strategy$type, firstDay),
        kind = kind,
        treatment = strategy$treatment,
        horizon = horizon,
        fixedVariance = fixedPriorVariances[[as.character(strategyId)]],
        description = description,
        shareMessage = sprintf("Strategy %d shares the censoring model of an earlier strategy", strategyId),
        firstDay = firstDay
      )
    )
    if (strategy$type == "initiateWithin") {
      models[[strategyId]]$graceStart <- as.integer(strategy$graceStart)
      models[[strategyId]]$graceRegime <- fitCensoringModelsArgs$graceRegime
      # With graceStart > 0 and the binary model at G, the weights of the days before graceStart come from
      # the initiation hazard model through graceStart - 1 (that of a notInitiateWithin(graceStart - 1)
      # strategy). With initiateWithinModel = "hazard", the hazard model through G gives them.
      if (strategy$graceStart > 0 && kind == "binary") {
        earlyHorizon <- as.integer(strategy$graceStart - 1L)
        models[[strategyId]]$early <- getFit(
          component = "hazard",
          kind = "hazard",
          treatment = strategy$treatment,
          horizon = earlyHorizon,
          fixedVariance = fixedPriorVariances[[paste0(strategyId, ".early")]],
          description = sprintf(
            "%s initiation hazard model through day %d for the days before the grace period of strategy %d (%s)",
            fitCensoringModelsArgs$hazardModelType, earlyHorizon, strategyId, strategy$label
          ),
          shareMessage = sprintf("Strategy %d shares the initiation hazard model through day %d of an earlier strategy for the days before its grace period", strategyId, earlyHorizon)
        )
      }
    }
    if (strategy$type == "notInitiateWithin" && fitCensoringModelsArgs$postGraceModel == "binary") {
      models[[strategyId]]$postGrace <- getFit(
        component = "binary",
        kind = "binary",
        treatment = strategy$treatment,
        horizon = horizon,
        fixedVariance = fixedPriorVariances[[paste0(strategyId, ".postGrace")]],
        description = sprintf(
          "binary censoring model at day %d for the weights after the grace period of strategy %d (%s)",
          horizon, strategyId, strategy$label
        ),
        shareMessage = sprintf("Strategy %d shares the binary model at day %d of an earlier strategy for its weights after the grace period", strategyId, horizon)
      )
    }
  }
  attr(models, "metaData") <- list(
    outcomeId = metaData$outcomeId,
    sameDayOrdering = sameDayOrdering,
    competingEvents = metaData$competingEvents,
    strategies = strategies,
    fitCensoringModelsArgs = argsList
  )
  class(models) <- "CcwCensoringModels"
  return(list(models = models, fits = fits))
}

# Fits shared between strategies and calls -----------------------------------------------------------

# The fitted model a strategy's censoring model is: the binary model at G ("binary"), the initiation
# hazard model through G of notInitiateWithin(G) and of initiateWithin(G) with initiateWithinModel =
# "hazard" ("hazard"), or the initiation hazard model of neverInitiate ("neverInitiate"), which never
# shares a fit with the others (Phase1Contracts section 5.3).
.censoringFitComponent <- function(kind, type, firstDay = 0L) {
  if (kind == "binary") {
    # The binary model of the persons who have not started before day firstDay > 0 has another risk set:
    return(if (firstDay > 0) sprintf("binaryFrom%d", as.integer(firstDay)) else "binary")
  }
  return(if (type == "neverInitiate") "neverInitiate" else "hazard")
}

# The key of a fitted model: its component (.censoringFitComponent()), treatment, horizon, fixed prior
# variance ("cv": none) and the settings (`argsList`, from FitCensoringModelsArgs$toList()) the fit
# depends on. The design choices postGraceModel and initiateWithinModel choose between fits, so they
# are not part of it, and the binary model at G does not depend on the hazard model's settings. Every
# model depends on maxCohortSizeForFitting (and the seed of control, which draws the hazard models'
# samples).
.censoringFitKey <- function(component, treatment, horizon, fixedVariance, argsList) {
  argsList$postGraceModel <- NULL
  argsList$initiateWithinModel <- NULL
  argsList$graceRegime <- NULL
  # The number of threads (set by runCcwAnalyses() from its multi-threading settings) does not change
  # the fit:
  argsList$control$threads <- NULL
  if (startsWith(component, "binary")) {
    argsList$hazardModelType <- NULL
    argsList$timeBasis <- NULL
    argsList$covariateTimeInteraction <- NULL
    argsList$pooledLogisticEngine <- NULL
  }
  # The default engine is left out, so that the keys of models saved before the setting existed match:
  if (identical(argsList$pooledLogisticEngine, "lr")) {
    argsList$pooledLogisticEngine <- NULL
  }
  return(paste(
    component,
    treatment,
    as.integer(horizon),
    if (is.null(fixedVariance)) "cv" else sprintf("%.17g", fixedVariance),
    .settingsJson(argsList)
  ))
}

# A settings list (from toList()) as JSON with its named elements sorted by name at every level, because
# the order of the fields of settings objects (R6) depends on how they were made (for example clone()).
.settingsJson <- function(settingsList) {
  sortNames <- function(x) {
    if (!is.list(x)) {
      return(x)
    }
    if (!is.null(names(x))) {
      x <- x[order(names(x))]
    }
    return(lapply(x, sortNames))
  }
  return(as.character(jsonlite::toJSON(sortNames(settingsList), auto_unbox = TRUE, digits = NA, null = "null")))
}

# The fits (named by .censoringFitKey()) of CcwCensoringModels, so they can be reused by
# .fitCensoringModelsWithFits(). Models fit with .fixedPriorVariances cannot be told from
# cross-validated ones, so pass only models fit without them.
.censoringFitsFromModels <- function(censoringModels) {
  argsList <- attr(censoringModels, "metaData")$fitCensoringModelsArgs
  fits <- list()
  for (model in censoringModels) {
    fit <- model[setdiff(names(model), c("strategyId", "type", "treatment", "postGrace", "graceStart", "graceRegime", "early"))]
    firstDay <- if (is.null(model$firstDay)) 0L else model$firstDay
    fits[[.censoringFitKey(.censoringFitComponent(model$kind, model$type, firstDay), model$treatment, model$horizon, NULL, argsList)]] <- fit
    if (!is.null(model$postGrace)) {
      fits[[.censoringFitKey("binary", model$treatment, model$postGrace$horizon, NULL, argsList)]] <- model$postGrace
    }
    if (!is.null(model$early)) {
      fits[[.censoringFitKey("hazard", model$treatment, model$early$horizon, NULL, argsList)]] <- model$early
    }
  }
  return(fits)
}

# Risk sets ----------------------------------------------------------------------------------------

# Persons under follow-up at day G who have not started before day `firstDay` (graceStart, or G with
# graceRegime = "naturalGrace"), with their subject (the unit of cross-validation folds, when `persons`
# has one) and y = started on or before day G.
.getBinaryRiskSet <- function(persons, treatmentDay, graceEnd, sameDayOrdering, firstDay = 0L) {
  eventDay <- pmin(persons$outcomeDay, persons$deathDay, na.rm = TRUE)
  if (sameDayOrdering == "eventsFirst") {
    noEvent <- is.na(eventDay) | eventDay > graceEnd
  } else {
    noEvent <- is.na(eventDay) | eventDay >= graceEnd
  }
  inRiskSet <- noEvent & persons$riskEnd >= graceEnd & (is.na(treatmentDay) | treatmentDay >= firstDay)
  return(tibble(
    rowId = persons$rowId[inRiskSet],
    personSeqId = persons$personSeqId[inRiskSet],
    subjectId = persons[["subjectId"]][inRiskSet], # NULL (no column) for persons without subjects
    y = as.integer(!is.na(treatmentDay[inRiskSet]) & treatmentDay[inRiskSet] <= graceEnd)
  ))
}

# Persons at risk of initiating on day 0, with their subject (the unit of cross-validation folds), the
# last day at risk (through the horizon) and whether they initiated on that day.
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
    subjectId = persons$subjectId[inRiskSet],
    lastDay = as.integer(lastDay[inRiskSet]),
    y = as.integer(initiated[inRiskSet])
  ))
}

# The seed of the samples of the hazard models' risk sets: that of control (modulo
# .Machine$integer.max), or 1 when it has none.
.defaultSubsampleSeed <- function(fitCensoringModelsArgs) {
  seed <- fitCensoringModelsArgs$control$seed
  return(if (is.null(seed)) 1L else as.integer(as.numeric(seed) %% .Machine$integer.max))
}

# The persons a hazard model is fit on: its whole risk set or, when that holds more than
# maxCohortSizeForFitting (> 0) persons, a uniform sample of that many persons (with all their rows),
# drawn with `seed` and R's default random number kinds, leaving the caller's stream as it was. A
# sample without initiations, from a risk set with some, would give a hazard of 0 to everyone, so the
# whole risk set is used instead.
.sampleRiskSet <- function(riskSet, maxCohortSizeForFitting, seed, description) {
  if (maxCohortSizeForFitting == 0 || nrow(riskSet) <= maxCohortSizeForFitting) {
    return(riskSet)
  }
  index <- withr::with_seed(
    seed,
    sample.int(nrow(riskSet), maxCohortSizeForFitting),
    .rng_kind = "Mersenne-Twister",
    .rng_normal_kind = "Inversion",
    .rng_sample_kind = "Rejection"
  )
  fitSet <- riskSet[sort(index), , drop = FALSE]
  if (sum(fitSet$y) == 0 && sum(riskSet$y) > 0) {
    message(sprintf(
      "Fitting the %s on all its %d persons: the sample of %d (maxCohortSizeForFitting) has no initiations",
      description, nrow(riskSet), as.integer(maxCohortSizeForFitting)
    ))
    return(riskSet)
  }
  message(sprintf(
    "Fitting the %s on a uniform sample of %d of its %d persons (maxCohortSizeForFitting)",
    description, as.integer(maxCohortSizeForFitting), nrow(riskSet)
  ))
  return(fitSet)
}

# The sample diagnostics of a model fit on `fitSet`, the risk set or a sample of it.
.sampleDiagnostics <- function(riskSet, fitSet) {
  return(list(sampled = nrow(fitSet) < nrow(riskSet), sampleSize = nrow(fitSet), sampleEvents = sum(fitSet$y)))
}

# The first day of the binary model at G of an initiateWithin strategy: its risk set holds the persons
# who have not started before it. That is graceStart (graceRegime = "binaryAtG"), or G
# ("naturalGrace": the persons who reach day G untreated, whose probability of starting on day G gives the
# last-interval weight).
.binaryModelFirstDay <- function(strategy, graceRegime) {
  if (identical(graceRegime, "naturalGrace")) {
    return(as.integer(strategy$graceEnd))
  }
  return(as.integer(strategy$graceStart))
}

# Warns when, with graceRegime = "naturalGrace", persons reach day G of an initiateWithin strategy
# untreated and under follow-up but none of them starts on day G. Only the clones that start on day G
# represent these persons after G, so the regime is not identified: they are all censored at G, the fit
# of the day-G model succeeds without events, no weight is infinite, and the risks after G stand for the
# persons who started before G only.
.warnNoNaturalGraceStarts <- function(strategies, persons, sameDayOrdering) {
  for (strategyId in seq_along(strategies)) {
    strategy <- strategies[[strategyId]]
    if (strategy$type != "initiateWithin") {
      next
    }
    treatmentDay <- if (strategy$treatment == "treatment") persons$treatmentDay else persons$comparatorTreatmentDay
    graceEnd <- as.integer(strategy$graceEnd)
    riskSet <- .getBinaryRiskSet(persons, treatmentDay, graceEnd, sameDayOrdering, graceEnd)
    if (nrow(riskSet) > 0 && sum(riskSet$y) == 0) {
      warning(sprintf(
        "With graceRegime = \"naturalGrace\", none of the %d persons who reach day %d (graceEnd) of strategy %d (%s) untreated and under follow-up starts on day %d, so no clone represents them after day %d: the regime is not identified, and the risks after day %d stand for the persons who started before day %d only. See computePositivity(), or use graceRegime = \"binaryAtG\"",
        nrow(riskSet), graceEnd, strategyId, strategy$label, graceEnd, graceEnd, graceEnd, graceEnd
      ), call. = FALSE)
    }
  }
  invisible(NULL)
}

# With initiateWithinModel = "hazard", the weights of an initiateWithin strategy use S(graceStart - 1 | X)
# and, with graceRegime = "naturalGrace", S(G - 1 | X) of the pooled logistic hazard model, which it
# evaluates at the ends of its intervals (counted from day 0). So graceStart, and G for "naturalGrace",
# must start an interval: otherwise the weights would use the survival through an earlier interval end.
.checkInitiateWithinHazardGrid <- function(strategies, fitCensoringModelsArgs) {
  if (fitCensoringModelsArgs$initiateWithinModel != "hazard" || fitCensoringModelsArgs$hazardModelType != "logistic") {
    return(invisible(NULL))
  }
  intervalLength <- fitCensoringModelsArgs$timeBasis$intervalLength
  for (strategy in strategies) {
    if (strategy$type != "initiateWithin") {
      next
    }
    if (strategy$graceStart %% intervalLength != 0) {
      stop(sprintf(
        "With initiateWithinModel = \"hazard\", the graceStart of strategy '%s' (%d) must be a multiple of the hazard model's intervalLength (%d), so that the days before it are whole intervals. Use daily intervals or initiateWithinModel = \"binary\"",
        strategy$label, as.integer(strategy$graceStart), as.integer(intervalLength)
      ))
    }
    if (identical(fitCensoringModelsArgs$graceRegime, "naturalGrace") && strategy$graceEnd %% intervalLength != 0) {
      stop(sprintf(
        "With initiateWithinModel = \"hazard\" and graceRegime = \"naturalGrace\", the graceEnd of strategy '%s' (%d) must be a multiple of the hazard model's intervalLength (%d), so that the last interval is day %d alone. Use daily intervals or initiateWithinModel = \"binary\"",
        strategy$label, as.integer(strategy$graceEnd), as.integer(intervalLength), as.integer(strategy$graceEnd)
      ))
    }
  }
  invisible(NULL)
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

# Warns, once, when a strategy's censoring model predicts initiation on day 0 while the censoring models
# use covariates whose window ends on day 0 or later (`endDay >= 0` in analysisRef, as with
# FeatureExtraction's default endDays = 0): records of day 0 may follow the treatment decision of that
# day (Phase1Contracts section 11.6). The hazard models of "notInitiateWithin" and "neverInitiate"
# strategies have a risk set that starts on day 0, and so have those of "initiateWithin" strategies
# (with initiateWithinModel = "hazard", or for the days before graceStart > 0). The binary model at G of
# an "initiateWithin" strategy predicts initiation on days firstDay..G, so it predicts day 0 unless it
# is fit from a later day (graceRegime = "naturalGrace" with G > 0). Analyses without a window (endDay
# NA, such as demographics) are not a concern.
.warnDay0CovariateWindows <- function(strategies, cloneCensorWeightData, fitCensoringModelsArgs) {
  day0 <- vapply(strategies, function(strategy) {
    strategy$type != "initiateWithin" || fitCensoringModelsArgs$initiateWithinModel == "hazard" ||
      strategy$graceStart > 0 || .binaryModelFirstDay(strategy, fitCensoringModelsArgs$graceRegime) == 0
  }, logical(1))
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
    "The censoring model%s of %s %s model%s initiation on day 0, but the covariate windows of analyses %s end on day 0 or later, so the censoring models can use records of day 0 that may follow the treatment decision. Extract the covariates with windows that end before day 0, for example with endDays = -1 in FeatureExtraction::createCovariateSettings()",
    if (length(ids) == 1) "" else "s",
    if (length(ids) == 1) "strategy" else "strategies",
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
.fitCyclops <- function(cyclopsData, prior, control, stopOnError, fixedCoefficients = NULL) {
  cyclopsFit <- tryCatch(
    Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control, fixedCoefficients = fixedCoefficients),
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
  sample <- .sampleDiagnostics(riskSet, riskSet)
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
    # One row per person, so subjects matter only when persons share one (copies in a bootstrap
    # replicate). createPs() passes subjectId on to Cyclops as the unit of its folds.
    if (anyDuplicated(population$subjectId) == 0) {
      population$subjectId <- NULL
    }
    ps <- CohortMethod::createPs(cohortMethodData = view, population = population, createPsArgs = createPsArgs)
    # createPs() samples the starters and the non-starters to at most maxCohortSizeForFitting each
    # (only when it reaches the fit, not on its "..., so cannot fit model" paths):
    maxSize <- fitCensoringModelsArgs$maxCohortSizeForFitting
    if (maxSize > 0 && !grepl("so cannot fit model$", attr(ps, "metaData")$psError)) {
      sampleEvents <- as.integer(min(sum(riskSet$y), maxSize))
      sample <- list(
        sampled = max(sum(riskSet$y), sum(1 - riskSet$y)) > maxSize,
        sampleSize = sampleEvents + as.integer(min(sum(1 - riskSet$y), maxSize)),
        sampleEvents = sampleEvents
      )
    }
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
    diagnostics = c(
      list(
        riskSetSize = nrow(riskSet),
        rows = nrow(riskSet),
        events = sum(riskSet$y),
        auc = auc,
        status = error
      ),
      sample
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

# Fit on `fitSet` (the risk set, or a sample of its persons), with the linear predictor of every person
# of the risk set.
.fitCoxCensoringModel <- function(riskSet, cloneCensorWeightData, fitCensoringModelsArgs, prior, fitSet = riskSet) {
  riskSet$time <- riskSet$lastDay + 1
  fitSet$time <- fitSet$lastDay + 1
  numeratorHazard <- .computeBreslowHazard(fitSet$time, fitSet$y)
  coefficients <- numeric(0)
  priorVariance <- NA_real_
  error <- "OK"
  linearPredictor <- rep(0, nrow(riskSet))
  covariateData <- NULL
  if (sum(fitSet$y) > 0) {
    covariateData <- .getRiskSetCovariateData(riskSet, cloneCensorWeightData, fitCensoringModelsArgs)
  }
  if (!is.null(covariateData)) {
    on.exit(close(covariateData))
    covariateIds <- covariateData$covariates |>
      dplyr::distinct(.data$covariateId) |>
      pull()
    if (length(covariateIds) > 0) {
      error <- .checkHighCorrelation(covariateData, fitSet, fitCensoringModelsArgs)
      if (is.null(error)) {
        error <- "OK"
        # One row per person, so subjects matter only when persons share one (copies in a bootstrap
        # replicate); without them each row is its own fold unit. Cyclops leaves out the covariates of
        # persons outside the sample.
        outcomeColumns <- c("rowId", "y", "time", if (anyDuplicated(fitSet$subjectId) > 0) "subjectId")
        covariateData$outcomes <- fitSet |>
          select(dplyr::all_of(outcomeColumns))
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
    .computeBreslowHazard(fitSet$time, fitSet$y, linearPredictor[match(fitSet$rowId, riskSet$rowId)])
  }
  return(list(
    modelType = "cox",
    coefficients = coefficients,
    priorVariance = priorVariance,
    predictions = tibble(rowId = riskSet$rowId, linearPredictor = linearPredictor),
    baselineHazard = baselineHazard,
    numerator = list(baselineHazard = numeratorHazard),
    diagnostics = c(
      list(
        riskSetSize = nrow(riskSet),
        rows = nrow(fitSet),
        events = sum(riskSet$y),
        auc = NA_real_,
        status = error
      ),
      .sampleDiagnostics(riskSet, fitSet)
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

# TRUE when the loaded Cyclops has pooled logistic regression (modelType = "plr"). Not every Cyclops
# 4.0.2 has it: the first commits of 4.0.2 on Cyclops develop came before plr (45957c11).
.cyclopsHasPooledLogistic <- function() {
  return("convertToLongOutcome" %in% getNamespaceExports("Cyclops"))
}

# The outcomes of the pooled logistic regression with Cyclops's modelType = "plr" (interval k is plr
# interval k + 1): a person's rows are the intervals of `fitIntervals` up to their last interval. When
# the fit intervals are consecutive, a person is one row (implicit rows: entry, exit, y); otherwise a
# person has one row per interval (explicit rows: interval, y), still without covariates. Returns the
# `outcomes`, the `persons` with rows (rowId, lastInterval) and the number of person-interval rows.
.createPooledRows <- function(fitSet, fitIntervals) {
  if (length(fitIntervals) > 0 && all(diff(fitIntervals) == 1L)) {
    first <- fitIntervals[1]
    last <- fitIntervals[length(fitIntervals)]
    persons <- fitSet[fitSet$lastInterval >= first, , drop = FALSE]
    outcomes <- tibble(
      rowId = persons$rowId,
      subjectId = persons$subjectId,
      entry = first + 1L,
      exit = pmin(persons$lastInterval, last) + 1L,
      y = as.numeric(persons$y == 1 & persons$lastInterval <= last)
    )
    nRows <- sum(outcomes$exit - outcomes$entry + 1)
  } else {
    n <- findInterval(fitSet$lastInterval, fitIntervals)
    persons <- fitSet[n > 0, , drop = FALSE]
    person <- rep(seq_len(nrow(fitSet)), n)
    k <- fitIntervals[sequence(n)]
    outcomes <- tibble(
      rowId = fitSet$rowId[person],
      subjectId = fitSet$subjectId[person],
      interval = k + 1L,
      y = as.numeric(fitSet$y[person] == 1 & k == fitSet$lastInterval[person])
    )
    nRows <- nrow(outcomes)
  }
  return(list(outcomes = outcomes, persons = persons[, c("rowId", "lastInterval")], nRows = nRows))
}

# Fits the pooled logistic regression of the "plr" engine on `pooled` (.createPooledRows()): the
# covariates (an Andromeda table, or NULL) stay one row per person, and the time basis (its columns with
# rows, `presentBasisIds`) and the covariate-by-time interactions (`interactionRef`: column c of
# `interactionBasis` multiplies every covariate) are given per interval. An interaction column that is
# not a time-basis column (such as "day0") enters the basis with its coefficient fixed at 0, so that the
# model is that of the "lr" engine, and an interaction without rows (no person with the covariate reaches
# an interval where its column is not 0) is left out, as the "lr" engine has no column for it. Returns
# the fit and its coefficients (the intercept, then by id, as Cyclops orders those of "lr"), or the error
# message.
.fitPooledRows <- function(pooled, covariates, basis, basisIds, presentBasisIds, interactionRef, interactionBasis, rowIntervals, prior, fitCensoringModelsArgs) {
  outcomes <- pooled$outcomes
  if (is.null(covariates)) {
    covariates <- tibble(rowId = numeric(0), covariateId = numeric(0), covariateValue = numeric(0))
  } else {
    covariates <- collect(covariates)
  }
  columnIds <- basisIds[match(colnames(interactionBasis), colnames(basis))]
  extra <- is.na(columnIds)
  columnIds[extra] <- min(c(0, basisIds, interactionRef$interactionId)) - seq_len(sum(extra))
  values <- cbind(basis[, match(presentBasisIds, basisIds), drop = FALSE], interactionBasis[, extra, drop = FALSE])
  timeBasis <- tibble(
    interval = rep(seq_len(nrow(values)), ncol(values)),
    covariateId = rep(c(presentBasisIds, columnIds[extra]), each = nrow(values)),
    covariateValue = as.vector(values)
  ) |>
    filter(.data$covariateValue != 0)
  firstRow <- vapply(seq_len(ncol(interactionBasis)), function(column) {
    min(c(Inf, rowIntervals[interactionBasis[rowIntervals + 1L, column] != 0]))
  }, numeric(1))
  lastInterval <- tapply(pooled$persons$lastInterval[match(covariates$rowId, pooled$persons$rowId)], covariates$covariateId, max)
  covariateLast <- lastInterval[as.character(interactionRef$covariateId)]
  keep <- !is.na(covariateLast) & firstRow[interactionRef$columnIndex] <= covariateLast
  timeEffectMap <- tibble(
    covariateId = interactionRef$covariateId,
    timeCovariateId = columnIds[interactionRef$columnIndex],
    interactionId = interactionRef$interactionId
  )[keep, ]
  cyclopsData <- Cyclops::convertToCyclopsData(
    outcomes,
    covariates,
    modelType = "plr",
    timeBasis = if (nrow(timeBasis) > 0) timeBasis,
    timeEffectMap = if (nrow(timeEffectMap) > 0) timeEffectMap,
    quiet = TRUE
  )
  fixed <- as.numeric(Cyclops::getCovariateIds(cyclopsData)) %in% columnIds[extra]
  cyclopsFit <- .fitCyclops(cyclopsData, prior, fitCensoringModelsArgs$control, fitCensoringModelsArgs$stopOnError, fixedCoefficients = if (any(fixed)) fixed)
  if (is.character(cyclopsFit)) {
    return(cyclopsFit)
  }
  coefficients <- stats::coef(cyclopsFit)
  coefficients <- coefficients[!names(coefficients) %in% as.character(columnIds[extra])]
  others <- setdiff(names(coefficients), "(Intercept)")
  coefficients <- coefficients[c("(Intercept)", others[order(as.numeric(others))])]
  return(list(fit = cyclopsFit, coefficients = coefficients))
}

# Fit on `fitSet` (the risk set, or a sample of its persons), with the linear predictor of every person
# of the risk set.
.fitLogisticCensoringModel <- function(riskSet, horizon, cloneCensorWeightData, fitCensoringModelsArgs, prior, description = "initiation hazard model", fitSet = riskSet) {
  intervalLength <- fitCensoringModelsArgs$timeBasis$intervalLength
  nIntervals <- horizon %/% intervalLength + 1L
  k <- seq_len(nIntervals) - 1L
  intervals <- tibble(
    k = k,
    startDay = k * intervalLength,
    endDay = pmin((k + 1L) * intervalLength - 1L, horizon)
  )
  riskSet$lastInterval <- riskSet$lastDay %/% intervalLength
  fitSet$lastInterval <- fitSet$lastDay %/% intervalLength
  # With the indicator basis, the hazard of an interval without initiations is 0, and that of an
  # interval in which everyone at risk initiates is 1 (its coefficient is -Inf or Inf at the maximum
  # likelihood). The rows of these intervals carry no information on the other coefficients. So they
  # are left out of the fit, and the first remaining interval is the reference. This is the exact
  # limit of the fit, which would otherwise not converge.
  fitIntervals <- k
  fullIntervals <- integer(0)
  if (fitCensoringModelsArgs$timeBasis$type == "indicator" && sum(fitSet$y) > 0) {
    initiations <- tabulate(fitSet$lastInterval[fitSet$y == 1] + 1L, nbins = nIntervals)
    atRisk <- rev(cumsum(rev(tabulate(fitSet$lastInterval + 1L, nbins = nIntervals))))
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
  rowIntervals <- fitIntervals[fitIntervals <= max(c(-1L, fitSet$lastInterval))]
  interactionBasis <- .createInteractionBasis(fitCensoringModelsArgs$covariateTimeInteraction, basis, k, rowIntervals)
  interactionRef <- tibble(
    covariateId = rep(as.numeric(covariateIds), ncol(interactionBasis)),
    columnIndex = rep(seq_len(ncol(interactionBasis)), each = length(covariateIds))
  )
  interactionRef$interactionId <- -(ncol(basis) + seq_len(nrow(interactionRef)))
  # The engine: "lr" expands the rows (person x interval) with their covariates in Andromeda; "plr" keeps
  # the covariates one row per person and gives Cyclops the intervals (.createPooledRows()).
  engine <- fitCensoringModelsArgs$pooledLogisticEngine
  if (engine == "lr") {
    covariateData$interactionRef <- interactionRef
    covariateData$interactionBasis <- tibble(
      k = rep(k, ncol(interactionBasis)),
      columnIndex = rep(seq_len(ncol(interactionBasis)), each = nIntervals),
      zValue = as.vector(interactionBasis)
    ) |>
      filter(.data$zValue != 0)
    covariateData$riskSet <- fitSet |>
      select("rowId", "subjectId", "lastInterval", "y")
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
        riskSet.subjectId,
        intervals.k,
        CASE WHEN riskSet.y = 1 AND intervals.k = riskSet.lastInterval THEN 1 ELSE 0 END AS y
      FROM riskSet
      INNER JOIN intervals
        ON intervals.k <= riskSet.lastInterval;")
    # subjectId keeps the rows of a person in one cross-validation fold (Cyclops 4.0.1 or later)
    .andromedaExecute(covariateData, "
      CREATE TABLE hazardOutcomes AS
      SELECT rowId, y, subjectId
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
    # Basis columns with no rows (intervals nobody reaches) are not in the model:
    presentBasisIds <- covariateData$hazardBasis |>
      dplyr::distinct(.data$covariateId) |>
      pull()

    fitLr <- function(covariates, prior) {
      cyclopsData <- Cyclops::convertToCyclopsData(covariateData$hazardOutcomes, covariates, modelType = "lr", quiet = TRUE)
      cyclopsFit <- .fitCyclops(cyclopsData, prior, fitCensoringModelsArgs$control, fitCensoringModelsArgs$stopOnError)
      if (is.character(cyclopsFit)) {
        return(cyclopsFit)
      }
      return(list(fit = cyclopsFit, coefficients = stats::coef(cyclopsFit)))
    }
    fitNumerator <- function() fitLr(covariateData$hazardBasis, Cyclops::createPrior("none"))
    fitDenominator <- function(prior) fitLr(covariateData$hazardCovariates, prior)
    interceptOnly <- function() .fitInterceptOnly(covariateData)
  } else {
    pooled <- .createPooledRows(fitSet, fitIntervals)
    nRows <- pooled$nRows
    presentBasisIds <- basisIds[colSums(basis[rowIntervals + 1L, , drop = FALSE] != 0) > 0]
    interceptOnly <- function() c("(Intercept)" = stats::qlogis(sum(pooled$outcomes$y) / nRows))
    fitNumerator <- function() {
      .fitPooledRows(pooled, NULL, basis, basisIds, presentBasisIds, interactionRef[0, ], interactionBasis[, 0, drop = FALSE], rowIntervals, Cyclops::createPrior("none"), fitCensoringModelsArgs)
    }
    # Only the covariates of the persons with rows (the fit sample) leave Andromeda:
    covariateData$pooledPersons <- pooled$persons |>
      select("rowId")
    fitCovariates <- covariateData$covariates |>
      dplyr::inner_join(covariateData$pooledPersons, by = "rowId")
    fitDenominator <- function(prior) {
      .fitPooledRows(pooled, fitCovariates, basis, basisIds, presentBasisIds, interactionRef, interactionBasis, rowIntervals, prior, fitCensoringModelsArgs)
    }
  }
  events <- sum(fitSet$y)
  if (length(covariateIds) > 0 && events > 0 && nRows > 0) {
    requested <- if (fitCensoringModelsArgs$covariateTimeInteraction == "day0") "day0" else colnames(basis)
    .messageDroppedInteractions(fitCensoringModelsArgs$covariateTimeInteraction, requested, colnames(interactionBasis), rowIntervals, nIntervals, intervalLength, description)
  }

  # Numerator: intercept and time basis.
  numeratorError <- NULL
  if (nRows == 0) {
    # Every interval has a hazard of 0 or 1, so there is nothing to fit. The intercept is the first
    # interval's log-odds.
    numeratorCoefficients <- c("(Intercept)" = if (0L %in% fullIntervals) Inf else -Inf)
  } else if (length(presentBasisIds) == 0 || events == 0) {
    numeratorCoefficients <- interceptOnly()
  } else if (fitCensoringModelsArgs$timeBasis$type == "indicator") {
    # The indicator basis is saturated in the fitted intervals, so the maximum likelihood fit is the
    # observed hazard of each: the intercept is the log-odds of the reference (first) interval, and the
    # coefficient of each later interval's column (in interval order) its log-odds minus the reference's.
    logOdds <- stats::qlogis(initiations[fitIntervals + 1L] / atRisk[fitIntervals + 1L])
    numeratorCoefficients <- c(logOdds[1], logOdds[-1] - logOdds[1])
    names(numeratorCoefficients) <- c("(Intercept)", as.character(basisIds))
  } else {
    numeratorFit <- fitNumerator()
    if (is.character(numeratorFit)) {
      numeratorError <- paste("Fitting the numerator model failed:", numeratorFit)
      message(numeratorError, ". Using an intercept-only numerator model")
      numeratorCoefficients <- interceptOnly()
    } else {
      numeratorCoefficients <- numeratorFit$coefficients
    }
  }

  # Denominator: intercept, time basis and covariates.
  coefficients <- numeratorCoefficients
  priorVariance <- NA_real_
  error <- "OK"
  linearPredictor <- rep(0, nrow(riskSet))
  interactionPredictions <- tibble(rowId = integer(0), columnIndex = integer(0), value = numeric(0))
  if (length(covariateIds) > 0 && events > 0 && nRows > 0) {
    error <- .checkHighCorrelation(covariateData, fitSet, fitCensoringModelsArgs)
    if (is.null(error)) {
      error <- "OK"
      # The interaction columns are not excluded, so they are penalized like the covariates:
      prior <- .restrictPrior(prior, covariateIds, unpenalizedIds = c(0, presentBasisIds))
      cyclopsFit <- fitDenominator(prior)
      if (is.character(cyclopsFit)) {
        error <- cyclopsFit
      } else {
        coefficients <- cyclopsFit$coefficients
        priorVariance <- .getPriorVariance(prior, cyclopsFit$fit)
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
  baseline <- .computeBaseline(coefficients, basis, basisIds, emptyIntervals, fullIntervals)
  numeratorBaseline <- .computeBaseline(numeratorCoefficients, basis, basisIds, emptyIntervals, fullIntervals)
  if (nrow(fitSet) < nrow(riskSet) && fitCensoringModelsArgs$timeBasis$type == "indicator") {
    # An interval with a hazard of 0 or 1 in the sample but not in the whole risk set gets the observed
    # hazard of the whole risk set (shifted by the mean linear predictor of its persons in the model
    # with covariates), so that persons outside the sample keep finite weights:
    allInitiations <- tabulate(riskSet$lastInterval[riskSet$y == 1] + 1L, nbins = nIntervals)
    allAtRisk <- rev(cumsum(rev(tabulate(riskSet$lastInterval + 1L, nbins = nIntervals))))
    refit <- k[k %in% c(emptyIntervals, fullIntervals) & allInitiations > 0 & allInitiations < allAtRisk]
    for (interval in refit) {
      logOdds <- stats::qlogis(allInitiations[interval + 1L] / allAtRisk[interval + 1L])
      numeratorBaseline[interval + 1L] <- logOdds
      baseline[interval + 1L] <- logOdds - mean(linearPredictor[riskSet$lastInterval >= interval])
    }
    emptyIntervals <- setdiff(emptyIntervals, refit)
    fullIntervals <- setdiff(fullIntervals, refit)
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
    intervals = dplyr::mutate(intervals, baseline = baseline),
    numerator = list(
      coefficients = .renameBasisCoefficients(numeratorCoefficients, basisIds, colnames(basis)),
      intervals = dplyr::mutate(intervals, baseline = numeratorBaseline)
    ),
    diagnostics = c(
      list(
        riskSetSize = nrow(riskSet),
        rows = as.integer(nRows),
        events = sum(riskSet$y),
        auc = NA_real_,
        status = paste(c(numeratorError, if (error != "OK" || is.null(numeratorError)) error), collapse = "; "),
        engine = engine
      ),
      .sampleDiagnostics(riskSet, fitSet)
    )
  ))
}

# The interaction basis (intervals x columns) whose columns multiply every covariate: none, an indicator
# of the first interval ("day0": day 0 with daily intervals, the first interval otherwise), or the
# time-basis columns ("interval"). A column that takes a single value over `rowIntervals`, the intervals
# with rows in the fit, would duplicate the covariates, so it is left out (for example "day0" when the
# first interval is not fit, or nobody reaches a later one).
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

# Says when covariateTimeInteraction asks for interactions but .createInteractionBasis() leaves out some
# or all of their columns (`requested`; `keptColumns` are left), because they take a single value over
# the intervals in the fit (`rowIntervals`), so the model has fewer (or no) covariate-by-time
# interactions than asked for. For "day0", it also says why: nobody (or everybody) initiates in the first
# interval, so there is no first-interval effect to model, or the first interval is the only one fit,
# which shorter intervals can split.
.messageDroppedInteractions <- function(covariateTimeInteraction, requested, keptColumns, rowIntervals, nIntervals, intervalLength, description) {
  if (covariateTimeInteraction == "none") {
    return(invisible(NULL))
  }
  dropped <- setdiff(requested, keptColumns)
  if (length(dropped) == 0 && length(keptColumns) > 0) {
    return(invisible(NULL))
  }
  firstInterval <- if (intervalLength == 1) "day 0" else sprintf("the first interval (days 0-%d)", intervalLength - 1)
  reason <- ""
  if (covariateTimeInteraction == "day0") {
    if (!0L %in% rowIntervals) {
      reason <- sprintf(" Nobody (or everybody) at risk initiates on %s, so there is no effect of the covariates there to model.", firstInterval)
    } else if (intervalLength == 1) {
      reason <- " Day 0 is the only day in the fit, so the covariates' main effects are already their effects on day 0."
    } else {
      reason <- sprintf(" The first interval (days 0-%d) is the only one in the fit, for example when nobody starts in the days of the grace period after it; use shorter intervals to model the effects of the covariates on day 0.", intervalLength - 1)
    }
  }
  message(sprintf(
    "The %s %s covariateTimeInteraction = \"%s\": %s a single value over the intervals in the fit (%s of %d%s; with the indicator basis, intervals in which nobody or everybody at risk initiates are not fit), so they would duplicate the covariates and are left out.%s",
    description,
    if (length(keptColumns) == 0) "has no covariate-by-time interactions although" else sprintf("leaves out the interaction columns %s of", paste(dropped, collapse = ", ")),
    covariateTimeInteraction,
    if (length(keptColumns) == 0) "its interaction columns would take" else "they take",
    if (length(rowIntervals) == 0) "no intervals" else sprintf("interval%s %s", if (length(rowIntervals) == 1) "" else "s", paste(rowIntervals, collapse = ", ")),
    nIntervals,
    if (intervalLength == 1) "" else sprintf(", each %d days", intervalLength),
    reason
  ))
  invisible(NULL)
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

# S-bar(day) of a hazard model's numerator model (no covariates), the probability of no initiation on
# days 0..day (day <= horizon).
.computeNumeratorSurvival <- function(model, day) {
  if (model$modelType == "cox") {
    return(exp(-.evaluateCumulativeHazard(model$numerator$baselineHazard, day + 1)))
  }
  eta <- model$numerator$intervals$baseline[model$numerator$intervals$endDay <= day]
  return(exp(-sum(ifelse(eta > 30, eta, log1p(exp(eta))))))
}

# The probability of starting on days firstDay..G among the persons who have not started before day
# firstDay, from an initiation hazard model through G: 1 - S(G | X) / S(firstDay - 1 | X), with S(-1) = 1
# (1 - S(G | X) for firstDay = 0), for the persons `rowIds` of its risk set; and the same from its
# numerator model.
.hazardStartProbability <- function(model, rowIds, firstDay, graceEnd) {
  return(1 - .computeHazardSurvival(model, rowIds, graceEnd) / .computeHazardSurvival(model, rowIds, firstDay - 1L))
}

.hazardNumeratorStartProbability <- function(model, firstDay, graceEnd) {
  return(1 - .computeNumeratorSurvival(model, graceEnd) / .computeNumeratorSurvival(model, firstDay - 1L))
}

# The initiation hazard model of an initiateWithin(G) strategy with initiateWithinModel = "hazard" as the
# binary model at G whose weights after G it gives: pi(X) = 1 - S(G | X) for the persons of its risk
# set, and pi-bar = 1 - S-bar(G).
.hazardAsBinaryModel <- function(model) {
  return(list(
    strategyId = model$strategyId,
    type = model$type,
    treatment = model$treatment,
    kind = "binary",
    modelType = "binary",
    horizon = model$horizon,
    predictions = tibble(
      rowId = model$predictions$rowId,
      propensityScore = 1 - .computeHazardSurvival(model, model$predictions$rowId, model$horizon)
    ),
    numerator = list(propensityScore = 1 - .computeNumeratorSurvival(model, model$horizon))
  ))
}

# The model that gives the weights of an initiateWithin strategy, for .writeModelTables(). With
# graceStart = 0 and graceRegime = "binaryAtG": the binary model at G, or the hazard model as one
# (.hazardAsBinaryModel()). Otherwise the probability that a clone still followed after G followed the
# strategy through G is computed per person (the clones that started on days graceStart..G and are under
# follow-up at G): P = E(X) * p, where E(X) = S(graceStart - 1 | X) is the probability of not starting
# before graceStart (1 with graceStart = 0), from the hazard model of the days before graceStart (the
# `early` model with the binary model at G, else the hazard model through G), and p the probability of
# then following the strategy at G:
# - "binaryAtG": p = pi(X), the probability of starting by G of the persons under follow-up at G who
#   have not started before graceStart (binary model, or 1 - S(G | X) / S(graceStart - 1 | X));
# - "naturalGrace": p = 1 for a clone that started before G, and for a clone that started on day G the
#   probability of starting on day G of the persons who reach day G untreated (the last-interval weight:
#   binary model fit on them, or 1 - S(G | X) / S(G - 1 | X)).
# The numerator is the same product from the numerator models: S-bar(graceStart - 1) times the marginal
# probability of starting by G of the persons under follow-up at G who have not started before graceStart
# (binary model: their observed share; hazard model: 1 - S-bar(G) / S-bar(graceStart - 1)), under both
# regimes, so that it is the same for every clone of the strategy after G and the weighted Kaplan-Meier
# risks do not depend on stabilization. With graceStart = 0 the result is a binary model whose
# propensityScore is P; with graceStart > 0 it is the hazard model of the days before graceStart with
# windowEnd = graceStart - 1, switchDay = G, and the survival after the switch `postSurvival` (rowId,
# survival = P) with its numerator `postNumerator`.
.initiateWithinWeightModel <- function(model, persons, sameDayOrdering) {
  graceStart <- if (is.null(model$graceStart)) 0L else as.integer(model$graceStart)
  graceRegime <- if (is.null(model$graceRegime)) "binaryAtG" else model$graceRegime
  if (graceStart == 0 && graceRegime == "binaryAtG") {
    return(if (model$kind == "hazard") .hazardAsBinaryModel(model) else model)
  }
  graceEnd <- as.integer(model$horizon)
  treatmentDay <- if (model$treatment == "treatment") persons$treatmentDay else persons$comparatorTreatmentDay
  riskSet <- .getBinaryRiskSet(persons, treatmentDay, graceEnd, sameDayOrdering, graceStart)
  followers <- riskSet[riskSet$y == 1, ]
  startDay <- treatmentDay[match(followers$rowId, persons$rowId)]
  earlyModel <- if (model$kind == "hazard") model else model$early
  early <- rep(1, nrow(followers))
  earlyNumerator <- 1
  if (graceStart > 0) {
    early <- .computeHazardSurvival(earlyModel, followers$rowId, graceStart - 1L)
    earlyNumerator <- .computeNumeratorSurvival(earlyModel, graceStart - 1L)
  }
  firstDay <- .binaryModelFirstDay(list(graceStart = graceStart, graceEnd = graceEnd), graceRegime)
  # The clones at risk of deviating at G (all of them with "binaryAtG"; those that start on day G with
  # "naturalGrace"):
  atRisk <- startDay >= firstDay
  following <- rep(1, nrow(followers))
  if (model$kind == "binary") {
    following[atRisk] <- model$predictions$propensityScore[match(followers$rowId[atRisk], model$predictions$rowId)]
    numeratorFollowing <- if (graceRegime == "naturalGrace") mean(riskSet$y) else model$numerator$propensityScore
  } else {
    following[atRisk] <- .hazardStartProbability(model, followers$rowId[atRisk], firstDay, graceEnd)
    numeratorFollowing <- .hazardNumeratorStartProbability(model, graceStart, graceEnd)
  }
  if (graceStart == 0) {
    return(list(
      strategyId = model$strategyId,
      type = model$type,
      treatment = model$treatment,
      kind = "binary",
      modelType = "binary",
      horizon = graceEnd,
      predictions = tibble(rowId = followers$rowId, propensityScore = following),
      numerator = list(propensityScore = numeratorFollowing)
    ))
  }
  earlyModel$strategyId <- model$strategyId
  earlyModel$type <- model$type
  earlyModel$treatment <- model$treatment
  earlyModel$windowEnd <- graceStart - 1L
  earlyModel$switchDay <- graceEnd
  earlyModel$postSurvival <- tibble(rowId = followers$rowId, survival = early * following)
  earlyModel$postNumerator <- earlyNumerator * numeratorFollowing
  return(earlyModel)
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

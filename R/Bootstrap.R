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

#' Bootstrap a clone-censor-weight analysis
#'
#' @description
#' Estimates the outcome model's estimates with a person-level nonparametric bootstrap that repeats the
#' analysis in every replicate: the censoring models, the weights and the outcome model. The bootstrap
#' standard errors and confidence intervals therefore include the estimation of the censoring weights,
#' which the robust variance of [fitOutcomeModel()] ignores.
#'
#' @details
#' The persons (`personSeqId`) of the `cohorts` table of `cloneCensorWeightData` are resampled with
#' replacement: each replicate draws multinomial counts (`n` draws among the `n` persons). All counts
#' are drawn with the seed of `fitOutcomeModelArgs$bootstrapArgs` before any model is fitted.
#'
#' `bootstrapArgs$method` sets how a replicate is analyzed:
#'
#' - `"reweight"` (the default): each person's count is its frequency weight, and persons with count 0
#'   are left out. The study population and clones of the full data are built once. The censoring
#'   models are refit with the counts as case weights, on the Cyclops data of the full data
#'   (built once per worker, with [Cyclops::updateWeights()]), the weights of the drawn persons' clones
#'   are computed from them, and the outcome model is fit with each clone's weight times its person's
#'   count. With `bootstrapArgs$warmStart = TRUE` (the default) the refits start from the coefficients
#'   of the full-data censoring models.
#' - `"copy"`: a new [CloneCensorWeightData] holds each drawn person as many times as drawn, every copy
#'   with its own row IDs and `personSeqId`, so copies are distinct persons in every step of the
#'   analysis, which is repeated from the study population on.
#'
#' Every step of the reweight bootstrap counts a person as many times as its copies would count, so the
#' two give the same replicates (to the convergence tolerance of the fits) with the same counts,
#' `fixHyperparameters = TRUE` and `warmStart = FALSE`. In particular the covariates of the censoring
#' models are tidied as FeatureExtraction tidies those of the copies (a covariate the copies would
#' drop as redundant or infrequent is fixed at 0, and covariates are normalized by their largest value
#' among the drawn persons, through scaled prior variances); a hazard model whose risk set holds more
#' copies than `maxCohortSizeForFitting` is fit on the persons of the sample of copies the copy
#' bootstrap would draw, each weighted by its number of sampled copies; and the percentile caps of the
#' weights count each clone as many times as its person was drawn. A replicate refits the full data's
#' Cyclops data only when its fit has the same structure (with the indicator time basis, the same
#' intervals with initiations, and the same interaction columns) and, under sampling, persons the
#' full-data fit has; otherwise it builds Cyclops data of its own. The binary models are fit by
#' [CohortMethod::createPs()], which takes no case weights, on copies of the drawn persons in both
#' methods (their sampling is seeded per replicate in the reweight bootstrap). The two methods differ
#' (beyond the tolerance) when the hyperparameters are not fixed: with `fixHyperparameters = FALSE`
#' the reweight bootstrap cross-validates on the weighted persons, whose folds differ from those of the
#' copies.
#'
#' With `bootstrapArgs$fixHyperparameters = TRUE`, the covariate models of the censoring
#' models use the prior variances of the full-data fit (chosen by cross-validation, for example) in
#' every replicate, without cross-validation. With `FALSE`, the cross-validation of every censoring
#' model (binary, Cox and pooled logistic) keeps the copies of a person in one fold, so that no copy is
#' scored by a model fitted to another copy. The
#' pooled logistic outcome model gets clones split on its interval grid
#' (`computeCensoringWeights(outcomeTimeBasis = )`).
#'
#' In the replicates, the censoring models are fitted with `errorOnHighCorrelation = FALSE`: the check
#' for covariates highly correlated with initiation is a design diagnostic, which the full-data fit has
#' already passed, and duplicated persons can make it fail by chance in small data. A hazard model
#' whose risk set holds more than `maxCohortSizeForFitting` persons (see [fitCensoringModels()]) is fit
#' in each replicate on a new uniform sample of the replicate's persons, drawn with the seed of the
#' censoring models' `control` plus the replicate number, so the replicates stay reproducible. The warnings of the
#' replicates are not shown; those of a failed replicate are reported with its failure, except those
#' that the analysis of the full data has already given (such as warnings about the design). The
#' progress messages of every replicate (including those of CohortMethod and FeatureExtraction) are
#' shown; wrap the call in [suppressMessages()] to hide them. With `threads > 1`, the workers' messages
#' are not shown.
#'
#' A replicate fails when the analysis stops with an error (for example because a censoring model
#' cannot be fitted) or gives no finite estimate (for example because the pooled logistic outcome
#' model cannot be fitted). Failed replicates get `NA` estimates and are counted in
#' `metaData$bootstrap$failedReplicates`, with a warning. Replicates with a missing or infinite
#' estimate (for example the risks of a replicate whose pooled logistic death model cannot be fitted,
#' or the hazard ratio of a replicate whose Cox model cannot be fitted, whose Kaplan-Meier risks are
#' still estimated) are left out of that estimate's standard error and limits, as in
#' [boot::boot.ci()]; `metaData$bootstrap$finiteReplicates` gives the number of replicates used for
#' each estimate, with a warning when it is smaller than the number of replicates that did not fail.
#' Ratios (hazard, odds and risk
#' ratios) are summarized on the log scale. The standard error is the standard deviation of the
#' replicates. The 95% limits follow [boot::boot.ci()]: with `bootstrapArgs$ci = "percentile"` they are
#' its `"perc"` limits, the (R + 1) 0.025 and (R + 1) 0.975 order statistics of the R replicates,
#' interpolated on the normal quantile scale (the smallest or largest replicate when R < 39); with
#' `"normal"` its `"norm"` limits, the bias-corrected `2 t0 - mean(t) -+ 1.96 sd(t)`. The p-value tests
#' no effect with the estimate divided by its bootstrap standard error.
#'
#' The driver warns when the replicates sit predominantly on one side of the full-data estimate of the
#' hazard ratio (on the log scale) or of a risk difference: then the replicates do not vary around the
#' estimate as its sampling distribution would, and the limits may not cover it. One known cause is
#' censoring models with many covariates and a penalty, which the duplicated persons of each replicate
#' let overfit at the full-data prior variance, moving the replicates toward the unweighted estimate.
#' The criterion, assessed for an estimate with at least 20 finite replicates, is a median replicate
#' more than 0.5 bootstrap standard deviations from the estimate, or more than 0.25 standard deviations
#' with a two-sided sign test of the numbers of replicates above and below the estimate (ties left out)
#' with p < 0.01. The sign test alone would flag ever smaller shifts as the number of replicates grows
#' (with 1000 replicates, about half of the estimates shifted by 0.1 standard deviation, which does not
#' matter for the limits), hence the floor of 0.25. Without a shift, the criterion still flags about 7%
#' of estimates with 20 replicates, 1% with 50, and fewer than 1% with 100 or more; it flags about half
#' of the estimates shifted by 0.5 standard deviation with 20 replicates, nearly 90% with 100, and
#' nearly all with 200 or more.
#' `metaData$bootstrap$replicateShift` holds the statistics for these estimands.
#'
#' @param cloneCensorWeightData     An object of type [CloneCensorWeightData].
#' @param outcomeId                 The ID of the outcome.
#' @param createStudyPopulationArgs An object of type `CreateStudyPopulationArgs` as created using
#'                                  [CohortMethod::createCreateStudyPopulationArgs()].
#' @param createClonesArgs          An object of type `CreateClonesArgs` as created using
#'                                  [createCreateClonesArgs()].
#' @param fitCensoringModelsArgs    An object of type `FitCensoringModelsArgs` as created using
#'                                  [createFitCensoringModelsArgs()].
#' @param computeWeightsArgs        An object of type `ComputeWeightsArgs` as created using
#'                                  [createComputeWeightsArgs()].
#' @param fitOutcomeModelArgs       An object of type `FitOutcomeModelArgs` as created using
#'                                  [createFitOutcomeModelArgs()]. Its `bootstrapArgs` set the number of
#'                                  replicates, the seed, the type of interval, and whether the
#'                                  hyperparameters are fixed. Its `varianceEstimator` is not used.
#' @param censoringModels           The censoring models of the full data, as created by
#'                                  [fitCensoringModels()] with `fitCensoringModelsArgs` on the clones of
#'                                  this analysis, or `NULL` to fit them. Models fit on clones without an
#'                                  outcome (`metaData$outcomeId` `NA`), which [runCcwAnalyses()] shares
#'                                  across outcomes with `refitCensoringModelsForEveryOutcome = FALSE`, are
#'                                  refit the same way in every replicate: on the clones of the
#'                                  replicate's study population without an outcome, so the replicates
#'                                  repeat the estimator of the full-data estimate.
#' @param counts                    For testing: `NULL`, or a matrix of non-negative integer counts with
#'                                  one row per person (in increasing order of `personSeqId`) and one
#'                                  column per replicate, used instead of drawing counts. It sets the
#'                                  number of replicates.
#' @param returnCounts              Return the counts of the replicates?
#' @param threads                   The number of parallel workers for the replicates. Each worker loads
#'                                  its own copy of `cloneCensorWeightData` (and, with the reweight
#'                                  bootstrap, builds the Cyclops data of the censoring models once), so
#'                                  memory grows with the number of workers. With the reweight bootstrap,
#'                                  each replicate depends only on its counts, so the results do not
#'                                  depend on the number of workers; with `method = "copy"`, a binary
#'                                  model sampled by `maxCohortSizeForFitting` draws from the random
#'                                  number generator of its worker, so they can. The Cyclops threads of each fit are set by the `control` of
#'                                  the censoring and outcome model arguments.
#'
#' @return
#' An object of class `CcwOutcomeModel` (see [fitOutcomeModel()]) with `varianceEstimator =
#' "bootstrap"`: the full-data estimates with their bootstrap standard errors, limits and p-values,
#' the full-data coefficients with their bootstrap covariance (`vcov`), the counts of the full-data
#' clones, and `metaData$bootstrap` (the settings, the number of replicates, the number of failed
#' replicates, the number of finite replicates of each estimate, the fixed prior variances,
#' `replicateStatus`, a tibble with, per `replicate`, its `status` (`"OK"` or `"failed"`), the `reason`
#' it failed, and `censoringStatus`, the statuses of its censoring model fits other than `"OK"` (such
#' as Cyclops's return flag of a fit that did not converge, after which the model without covariates
#' was used) or `"OK"`, and
#' `replicateShift`, a tibble with, for the hazard ratio and each risk difference, `estimand`, the number
#' of finite `replicates`, the numbers `above` and `below` the estimate, the sign test p-value
#' `signTestP`, `medianShift` (the median replicate minus the estimate, in bootstrap standard deviations,
#' on the log scale for the hazard ratio), and `shifted` (whether the criterion above is met; `NA` with
#' fewer than 20 finite replicates)). For the
#' pooled logistic model with a natural spline or orthogonal polynomial time basis, whose knots or
#' polynomial coefficients are recomputed from the data of each replicate, the time-basis coefficients
#' (and the intercepts, and with an orthogonal polynomial and the interaction also the target
#' coefficients) mean different things in different replicates: their rows and columns of `vcov` are
#' `NA`. It also holds `replicates`, a matrix with one row per
#' replicate and one column per estimand of the estimates, and, with `returnCounts = TRUE`,
#' `bootstrapCounts`, the counts (persons by replicates, with the `personSeqId` as row names).
#'
#' @references
#' Davison AC, Hinkley DV. Bootstrap Methods and Their Application. Cambridge University Press; 1997.
#'
#' @examples
#' \dontrun{
#' data <- simulateCloneCensorWeightData(n = 1000, seed = 1)
#' result <- bootstrapCcwAnalysis(
#'   cloneCensorWeightData = data,
#'   outcomeId = 3,
#'   createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
#'     riskWindowEnd = 365,
#'     endAnchor = "cohort start"
#'   ),
#'   createClonesArgs = createCreateClonesArgs(
#'     strategies = list(
#'       createTreatmentStrategy(1, "Initiate within 30 days", "initiateWithin", graceEnd = 30),
#'       createTreatmentStrategy(2, "Not within 30 days", "notInitiateWithin", graceEnd = 30)
#'     )
#'   ),
#'   fitCensoringModelsArgs = createFitCensoringModelsArgs(),
#'   computeWeightsArgs = createComputeWeightsArgs(),
#'   fitOutcomeModelArgs = createFitOutcomeModelArgs(
#'     modelType = "km",
#'     horizons = 365,
#'     bootstrapArgs = createBootstrapArgs(replicates = 200)
#'   )
#' )
#' result
#' }
#'
#' @export
bootstrapCcwAnalysis <- function(cloneCensorWeightData,
                                 outcomeId,
                                 createStudyPopulationArgs,
                                 createClonesArgs,
                                 fitCensoringModelsArgs,
                                 computeWeightsArgs,
                                 fitOutcomeModelArgs,
                                 censoringModels = NULL,
                                 counts = NULL,
                                 returnCounts = FALSE,
                                 threads = 1) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCloneCensorWeightData(cloneCensorWeightData, add = errorMessages)
  checkmate::assertNumber(outcomeId, add = errorMessages)
  checkmate::assertR6(createStudyPopulationArgs, "CreateStudyPopulationArgs", add = errorMessages)
  checkmate::assertR6(createClonesArgs, "CreateClonesArgs", add = errorMessages)
  checkmate::assertR6(fitCensoringModelsArgs, "FitCensoringModelsArgs", add = errorMessages)
  checkmate::assertR6(computeWeightsArgs, "ComputeWeightsArgs", add = errorMessages)
  # CohortMethod has an R6 class with the same name; require members only this package's class has.
  checkmate::assertR6(fitOutcomeModelArgs, "FitOutcomeModelArgs", public = c("horizons", "bootstrapArgs"), add = errorMessages)
  checkmate::assertClass(censoringModels, "CcwCensoringModels", null.ok = TRUE, add = errorMessages)
  checkmate::assertMatrix(counts, min.cols = 1, null.ok = TRUE, add = errorMessages)
  if (!is.null(counts)) {
    checkmate::assertIntegerish(counts, lower = 0, any.missing = FALSE, add = errorMessages)
  }
  checkmate::assertFlag(returnCounts, add = errorMessages)
  checkmate::assertInt(threads, lower = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  bootstrapArgs <- fitOutcomeModelArgs$bootstrapArgs

  personSeqIds <- cloneCensorWeightData$cohorts |>
    dplyr::distinct(.data$personSeqId) |>
    dplyr::arrange(.data$personSeqId) |>
    pull()
  nPersons <- length(personSeqIds)
  if (is.null(counts)) {
    # Drawn before any model is fitted: model fitting can use the random number generator.
    set.seed(bootstrapArgs$seed)
    counts <- stats::rmultinom(bootstrapArgs$replicates, size = nPersons, prob = rep(1, nPersons))
  } else {
    if (nrow(counts) != nPersons) {
      stop(sprintf("counts must have one row per person (%d), but has %d", nPersons, nrow(counts)))
    }
    if (!is.null(rownames(counts)) && !identical(rownames(counts), as.character(personSeqIds))) {
      stop("The row names of counts must be the personSeqIds in increasing order")
    }
  }
  storage.mode(counts) <- "integer"
  dimnames(counts) <- list(as.character(personSeqIds), NULL)
  nReplicates <- ncol(counts)

  start <- Sys.time()
  message("Fitting the analysis on the full data")
  # Its warnings (for example about the design) are shown, and left out of the failure reasons of
  # replicates, which repeat them:
  fullWarnings <- character(0)
  full <- withCallingHandlers(.runBootstrapAnalysis(
    cloneCensorWeightData = cloneCensorWeightData,
    outcomeId = outcomeId,
    createStudyPopulationArgs = createStudyPopulationArgs,
    createClonesArgs = createClonesArgs,
    fitCensoringModelsArgs = fitCensoringModelsArgs,
    computeWeightsArgs = computeWeightsArgs,
    fitOutcomeModelArgs = fitOutcomeModelArgs,
    censoringModels = censoringModels
  ), warning = function(w) {
    fullWarnings <<- c(fullWarnings, conditionMessage(w))
  })
  fixedPriorVariances <- NULL
  if (bootstrapArgs$fixHyperparameters) {
    variances <- lapply(full$censoringModels, function(model) model$priorVariance)
    names(variances) <- as.character(seq_along(variances))
    # The binary models of notInitiateWithin strategies' weights after the grace period:
    postGraceVariances <- lapply(full$censoringModels, function(model) model$postGrace$priorVariance)
    names(postGraceVariances) <- paste0(seq_along(postGraceVariances), ".postGrace")
    # The hazard models of the days before graceStart of initiateWithin strategies:
    earlyVariances <- lapply(full$censoringModels, function(model) model$early$priorVariance)
    names(earlyVariances) <- paste0(seq_along(earlyVariances), ".early")
    variances <- c(variances, postGraceVariances, earlyVariances)
    variances <- variances[!vapply(variances, function(variance) is.null(variance) || is.na(variance), logical(1))]
    if (length(variances) > 0) {
      fixedPriorVariances <- variances
    }
  }

  estimands <- full$estimates$estimand
  replicates <- matrix(NA_real_, nReplicates, length(estimands), dimnames = list(NULL, estimands))
  coefficientReplicates <- NULL
  if (!is.null(full$coefficients)) {
    coefficientReplicates <- matrix(NA_real_, nReplicates, length(full$coefficients), dimnames = list(NULL, names(full$coefficients)))
  }
  # The full-data fit has passed the high-correlation diagnostic; in replicates, duplicated persons can
  # fail it by chance, which would drop those replicates selectively.
  replicateCensoringModelsArgs <- fitCensoringModelsArgs$clone()
  replicateCensoringModelsArgs$errorOnHighCorrelation <- FALSE
  # Censoring models fit without an outcome (shared across outcomes by runCcwAnalyses()) are refit
  # without an outcome in every replicate, as the full-data models were:
  censoringWithoutOutcome <- !is.null(censoringModels) && isTRUE(is.na(attr(censoringModels, "metaData")$outcomeId))
  setup <- list(
    method = bootstrapArgs$method,
    warmStart = bootstrapArgs$warmStart,
    outcomeId = outcomeId,
    createStudyPopulationArgs = createStudyPopulationArgs,
    createClonesArgs = createClonesArgs,
    fitCensoringModelsArgs = replicateCensoringModelsArgs,
    computeWeightsArgs = computeWeightsArgs,
    fitOutcomeModelArgs = fitOutcomeModelArgs,
    fixedPriorVariances = fixedPriorVariances,
    censoringWithoutOutcome = censoringWithoutOutcome,
    personSeqIds = personSeqIds,
    nReplicates = nReplicates,
    estimands = estimands,
    fullWarnings = fullWarnings,
    # Warm starts of the reweight bootstrap's refits:
    startingFits = .startingFitsFromModels(full$censoringModels),
    # Each replicate draws its own samples of the hazard models' risk sets (maxCohortSizeForFitting):
    subsampleSeed = .defaultSubsampleSeed(fitCensoringModelsArgs)
  )
  replicateResults <- .runBootstrapReplicates(cloneCensorWeightData, setup, counts, threads)
  failures <- character(0)
  for (replicate in seq_len(nReplicates)) {
    replicateResult <- replicateResults[[replicate]]
    if (!is.null(replicateResult$failure)) {
      failures <- c(failures, sprintf("replicate %d: %s", replicate, replicateResult$failure))
      next
    }
    replicates[replicate, ] <- replicateResult$estimates
    if (!is.null(coefficientReplicates)) {
      coefficientReplicates[replicate, ] <- replicateResult$coefficients[colnames(coefficientReplicates)]
    }
  }
  replicateStatus <- tibble(
    replicate = seq_len(nReplicates),
    status = vapply(replicateResults, function(result) if (is.null(result$failure)) "OK" else "failed", character(1)),
    reason = vapply(replicateResults, function(result) if (is.null(result$failure)) NA_character_ else result$failure, character(1)),
    censoringStatus = vapply(replicateResults, function(result) result$censoringStatus, character(1))
  )
  if (length(failures) > 0) {
    warning(sprintf(
      "The analysis failed in %d of %d bootstrap replicates, which are left out: %s",
      length(failures), nReplicates, paste(utils::head(failures, 3), collapse = "; ")
    ))
  }
  finiteReplicates <- colSums(is.finite(replicates))
  incomplete <- finiteReplicates < nReplicates - length(failures)
  if (any(incomplete)) {
    warning(sprintf(
      "Some estimates are missing or infinite in bootstrap replicates that did not fail, and are left out of their standard errors and limits (finite replicates of %d): %s",
      nReplicates - length(failures),
      paste(sprintf("%s (%d)", estimands[incomplete], finiteReplicates[incomplete]), collapse = ", ")
    ))
  }

  replicateShift <- .bootstrapReplicateShift(full$estimates, replicates)
  shifted <- replicateShift[replicateShift$shifted %in% TRUE, ]
  if (nrow(shifted) > 0) {
    medianShift <- ifelse(
      is.na(shifted$medianShift),
      "",
      sprintf("; median replicate %.2f bootstrap SD %s it", abs(shifted$medianShift), ifelse(shifted$medianShift >= 0, "above", "below"))
    )
    warning(sprintf(
      "The bootstrap replicates sit predominantly on one side of the full-data estimate of %s, so the bootstrap limits may not cover the estimate. A known cause is censoring models with many covariates and a penalty, which the duplicated persons of each replicate let overfit at the full-data prior variance. Check the replicates (replicates and metaData$bootstrap$replicateShift) before using the bootstrap limits",
      paste(sprintf(
        "%s (%d replicates above it and %d below, sign test p = %.2g%s)",
        shifted$estimand, shifted$above, shifted$below, shifted$signTestP, medianShift
      ), collapse = "; ")
    ), call. = FALSE)
  }

  vcov <- NULL
  if (!is.null(coefficientReplicates)) {
    vcov <- stats::cov(coefficientReplicates, use = "pairwise.complete.obs")
    variable <- !(colnames(vcov) %in% .replicateInvariantCoefficients(colnames(vcov), fitOutcomeModelArgs))
    vcov[variable, ] <- NA
    vcov[, variable] <- NA
  }
  metaData <- full$metaData
  result <- list(
    modelType = fitOutcomeModelArgs$modelType,
    varianceEstimator = "bootstrap",
    estimates = .bootstrapEstimates(full$estimates, replicates, bootstrapArgs$ci),
    coefficients = full$coefficients,
    vcov = vcov,
    counts = full$counts,
    metaData = list(
      strategies = metaData$strategies,
      outcomeId = metaData$outcomeId,
      targetId = metaData$targetId,
      competingEvents = metaData$competingEvents,
      riskStart = full$riskStart,
      args = fitOutcomeModelArgs$toList(),
      bootstrap = list(
        replicates = nReplicates,
        failedReplicates = length(failures),
        finiteReplicates = finiteReplicates,
        seed = bootstrapArgs$seed,
        ci = bootstrapArgs$ci,
        fixHyperparameters = bootstrapArgs$fixHyperparameters,
        method = bootstrapArgs$method,
        warmStart = bootstrapArgs$warmStart,
        fixedPriorVariances = fixedPriorVariances,
        replicateShift = replicateShift,
        replicateStatus = replicateStatus
      )
    ),
    replicates = replicates
  )
  if (returnCounts) {
    result$bootstrapCounts <- counts
  }
  class(result) <- "CcwOutcomeModel"
  delta <- Sys.time() - start
  message(paste("Bootstrapping the analysis took", signif(delta, 3), attr(delta, "units")))
  return(result)
}

# The analysis of one (full or resampled) data set, up to the point estimates of the outcome model.
.runBootstrapAnalysis <- function(cloneCensorWeightData,
                                  outcomeId,
                                  createStudyPopulationArgs,
                                  createClonesArgs,
                                  fitCensoringModelsArgs,
                                  computeWeightsArgs,
                                  fitOutcomeModelArgs,
                                  censoringModels = NULL,
                                  fixedPriorVariances = NULL,
                                  sourcePersonSeqIds = NULL,
                                  censoringWithoutOutcome = FALSE,
                                  subsampleSeed = NULL) {
  population <- createStudyPopulation(
    cloneCensorWeightData = cloneCensorWeightData,
    outcomeId = outcomeId,
    createStudyPopulationArgs = createStudyPopulationArgs
  )
  clones <- createClones(population, cloneCensorWeightData, createClonesArgs)
  on.exit(close(clones))
  if (is.null(censoringModels)) {
    censoringClones <- clones
    if (censoringWithoutOutcome) {
      censoringClones <- createClones(
        .createStudyPopulationWithoutOutcome(cloneCensorWeightData, createStudyPopulationArgs),
        cloneCensorWeightData,
        createClonesArgs
      )
      on.exit(close(censoringClones), add = TRUE)
    }
    censoringModels <- fitCensoringModels(censoringClones, cloneCensorWeightData, fitCensoringModelsArgs,
                                          .fixedPriorVariances = fixedPriorVariances,
                                          .sourcePersonSeqIds = sourcePersonSeqIds,
                                          .subsampleSeed = subsampleSeed)
  }
  outcomeTimeBasis <- if (fitOutcomeModelArgs$modelType == "logistic") fitOutcomeModelArgs$timeBasis else NULL
  weightedClones <- computeCensoringWeights(clones, censoringModels, computeWeightsArgs, outcomeTimeBasis = outcomeTimeBasis)
  on.exit(close(weightedClones), add = TRUE)
  point <- .fitOutcomeModelPoint(weightedClones, fitOutcomeModelArgs)
  return(list(
    estimates = point$estimates,
    coefficients = point$coefficients,
    censoringModels = censoringModels,
    counts = .getOutcomeCounts(weightedClones),
    metaData = attr(weightedClones, "metaData"),
    riskStart = createStudyPopulationArgs$riskWindowStart
  ))
}

# Replicates ---------------------------------------------------------------------------------------

# Runs the replicates (`counts`: one column per replicate) of `setup` (see bootstrapCcwAnalysis()) on
# `threads` parallel workers, each with its own copy of the data, and returns one result per replicate
# (.bootstrapReplicate()). Every reweight replicate is a function of its counts only (the warm starts
# are the full-data coefficients in `setup`), so its results do not depend on the number of workers.
.runBootstrapReplicates <- function(cloneCensorWeightData, setup, counts, threads) {
  replicateIds <- seq_len(ncol(counts))
  threads <- min(threads, length(replicateIds))
  if (threads == 1) {
    return(.bootstrapReplicates(cloneCensorWeightData, setup, counts, replicateIds))
  }
  dataFile <- tempfile(fileext = ".zip")
  on.exit(unlink(dataFile))
  Andromeda::saveAndromeda(cloneCensorWeightData, dataFile, maintainConnection = TRUE)
  # Contiguous blocks of replicates, one per worker:
  blocks <- split(replicateIds, cut(replicateIds, threads, labels = FALSE))
  tasks <- lapply(blocks, function(ids) list(dataFile = dataFile, setup = setup, counts = counts[, ids, drop = FALSE], replicateIds = ids))
  cluster <- ParallelLogger::makeCluster(threads)
  on.exit(ParallelLogger::stopCluster(cluster), add = TRUE)
  .requirePackageOnCluster(cluster)
  results <- ParallelLogger::clusterApply(cluster, tasks, .bootstrapReplicatesTask)
  return(do.call(c, unname(results)))
}

# A worker's task of .runBootstrapReplicates().
.bootstrapReplicatesTask <- function(task) {
  data <- loadCloneCensorWeightData(task$dataFile)
  on.exit(close(data))
  counts <- matrix(0L, nrow(task$counts), max(task$replicateIds))
  counts[, task$replicateIds] <- task$counts
  return(.bootstrapReplicates(data, task$setup, counts, task$replicateIds))
}

# Loads this package on the workers of a cluster: the installed package, or, when this process loaded
# it from its source folder (pkgload::load_all(), as in the tests), that folder.
.requirePackageOnCluster <- function(cluster) {
  path <- getNamespaceInfo("CloneCensorWeight", "path")
  if (isNamespaceLoaded("pkgload") && pkgload::is_dev_package("CloneCensorWeight")) {
    loadSource <- function(i, path) {
      suppressMessages(pkgload::load_all(path, quiet = TRUE, export_all = FALSE))
      NULL
    }
    # Sent without this namespace as its environment, which the workers do not have yet:
    environment(loadSource) <- globalenv()
    ParallelLogger::clusterApply(cluster, seq_along(cluster), loadSource, path = path)
  } else {
    ParallelLogger::clusterRequire(cluster, "CloneCensorWeight")
  }
  invisible(NULL)
}

# The replicates `replicateIds` (columns of `counts`) in this process. The reweight bootstrap builds the
# study population and clones of the full data once, and keeps the Cyclops data of the censoring models
# for the replicates after the first.
.bootstrapReplicates <- function(cloneCensorWeightData, setup, counts, replicateIds) {
  context <- NULL
  if (setup$method == "reweight") {
    context <- suppressMessages(.reweightContext(cloneCensorWeightData, setup))
    on.exit(.closeReweightContext(context))
  }
  return(lapply(replicateIds, function(replicate) {
    .bootstrapReplicate(cloneCensorWeightData, setup, context, counts[, replicate], replicate)
  }))
}

# One replicate: a list with `estimates` (one per estimand of the full data), `coefficients`, `failure`
# (NULL, or why the replicate failed) and `censoringStatus` (the statuses of the censoring models other
# than "OK", or "OK"). The replicate's warnings are not shown; those of a failed replicate are its
# failure, except those that the full-data analysis has already given.
.bootstrapReplicate <- function(cloneCensorWeightData, setup, context, counts, replicate) {
  message(sprintf("Bootstrap replicate %d of %d", replicate, setup$nReplicates))
  subsampleSeed <- as.integer((as.numeric(setup$subsampleSeed) + replicate) %% .Machine$integer.max)
  replicateWarnings <- character(0)
  result <- withCallingHandlers(
    tryCatch(
      suppressMessages(if (setup$method == "reweight") {
        .reweightReplicate(cloneCensorWeightData, setup, context, counts, subsampleSeed)
      } else {
        .copyReplicate(cloneCensorWeightData, setup, counts, subsampleSeed)
      }),
      error = function(e) conditionMessage(e)
    ),
    warning = function(w) {
      if (!conditionMessage(w) %in% setup$fullWarnings) {
        replicateWarnings <<- c(replicateWarnings, conditionMessage(w))
      }
      invokeRestart("muffleWarning")
    }
  )
  if (is.character(result)) {
    return(list(failure = result, censoringStatus = NA_character_))
  }
  censoringStatus <- .censoringModelStatus(result$censoringModels)
  estimates <- result$estimates$estimate[match(setup$estimands, result$estimates$estimand)]
  if (!any(is.finite(estimates))) {
    # The outcome model could not be fitted (it warns why):
    return(list(
      failure = sprintf(
        "no finite estimate%s",
        if (length(replicateWarnings) == 0) "" else sprintf(" (%s)", paste(unique(replicateWarnings), collapse = "; "))
      ),
      censoringStatus = censoringStatus
    ))
  }
  return(list(estimates = estimates, coefficients = result$coefficients, censoringStatus = censoringStatus))
}

# A replicate of the copy bootstrap: the analysis of new data with each drawn person copied as many
# times as drawn.
.copyReplicate <- function(cloneCensorWeightData, setup, counts, subsampleSeed) {
  resampled <- .resampleCloneCensorWeightData(cloneCensorWeightData, setup$personSeqIds, counts)
  on.exit(close(resampled))
  return(.runBootstrapAnalysis(
    cloneCensorWeightData = resampled,
    outcomeId = setup$outcomeId,
    createStudyPopulationArgs = setup$createStudyPopulationArgs,
    createClonesArgs = setup$createClonesArgs,
    fitCensoringModelsArgs = setup$fitCensoringModelsArgs,
    computeWeightsArgs = setup$computeWeightsArgs,
    fitOutcomeModelArgs = setup$fitOutcomeModelArgs,
    fixedPriorVariances = setup$fixedPriorVariances,
    sourcePersonSeqIds = .copySourcePersonSeqIds(setup$personSeqIds, counts),
    censoringWithoutOutcome = setup$censoringWithoutOutcome,
    subsampleSeed = subsampleSeed
  ))
}

# The statuses of the censoring models (and of their postGrace and early models) other than "OK", or
# "OK": a model that could not be fit with covariates has its error or Cyclops's return flag.
.censoringModelStatus <- function(censoringModels) {
  statuses <- unlist(lapply(censoringModels, function(model) {
    c(model$diagnostics$status, model$postGrace$diagnostics$status, model$early$diagnostics$status)
  }))
  statuses <- unique(statuses[!is.na(statuses) & statuses != "OK"])
  return(if (length(statuses) == 0) "OK" else paste(statuses, collapse = "; "))
}

# The study population and clones of the full data (and the clones without an outcome the censoring
# models are fit on, with `censoringWithoutOutcome`), and the cache of the censoring models' Cyclops
# data, for the replicates of the reweight bootstrap.
.reweightContext <- function(cloneCensorWeightData, setup) {
  population <- createStudyPopulation(
    cloneCensorWeightData = cloneCensorWeightData,
    outcomeId = setup$outcomeId,
    createStudyPopulationArgs = setup$createStudyPopulationArgs
  )
  context <- new.env(parent = emptyenv())
  context$clones <- createClones(population, cloneCensorWeightData, setup$createClonesArgs)
  context$censoringClones <- context$clones
  context$separateCensoringClones <- setup$censoringWithoutOutcome
  if (setup$censoringWithoutOutcome) {
    context$censoringClones <- createClones(
      .createStudyPopulationWithoutOutcome(cloneCensorWeightData, setup$createStudyPopulationArgs),
      cloneCensorWeightData,
      setup$createClonesArgs
    )
  }
  context$cache <- new.env(parent = emptyenv())
  return(context)
}

.closeReweightContext <- function(context) {
  close(context$clones)
  if (context$separateCensoringClones) {
    close(context$censoringClones)
  }
  invisible(NULL)
}

# A replicate of the reweight bootstrap: the censoring models refit with each person's count as its
# frequency weight (persons with count 0 left out), the weights of the clones of the drawn persons, and
# the outcome model with each clone's weight times its person's count.
.reweightReplicate <- function(cloneCensorWeightData, setup, context, counts, subsampleSeed) {
  personCounts <- tibble(personSeqId = setup$personSeqIds, count = as.integer(counts))
  censoringModels <- .fitCensoringModelsWithFits(
    clones = context$censoringClones,
    cloneCensorWeightData = cloneCensorWeightData,
    fitCensoringModelsArgs = setup$fitCensoringModelsArgs,
    fixedPriorVariances = setup$fixedPriorVariances,
    subsampleSeed = subsampleSeed,
    refit = list(counts = personCounts, cache = context$cache, warmStart = setup$warmStart, startingFits = setup$startingFits)
  )$models
  fitOutcomeModelArgs <- setup$fitOutcomeModelArgs
  outcomeTimeBasis <- if (fitOutcomeModelArgs$modelType == "logistic") fitOutcomeModelArgs$timeBasis else NULL
  weightedClones <- computeCensoringWeights(context$clones, censoringModels, setup$computeWeightsArgs, outcomeTimeBasis = outcomeTimeBasis, .personCounts = personCounts)
  on.exit(close(weightedClones))
  point <- .fitOutcomeModelPoint(weightedClones, fitOutcomeModelArgs, personCounts = personCounts)
  return(list(estimates = point$estimates, coefficients = point$coefficients, censoringModels = censoringModels))
}

# The coefficients (of `coefficientNames`) that mean the same in every replicate. The time basis of the
# pooled logistic model is rebuilt from the data of each replicate for a natural spline (knots at
# quantiles) or an orthogonal polynomial (centring and scaling), so then only the target coefficients
# are comparable: the odds ratios without the interaction and, with it, where the basis is 0 (the first
# interval for the natural spline; a data-dependent point for the orthogonal polynomial, so none).
.replicateInvariantCoefficients <- function(coefficientNames, fitOutcomeModelArgs) {
  timeBasis <- fitOutcomeModelArgs$timeBasis
  dataDependent <- timeBasis$type == "naturalSpline" || (timeBasis$type == "polynomial" && !timeBasis$raw)
  if (fitOutcomeModelArgs$modelType != "logistic" || !dataDependent) {
    return(coefficientNames)
  }
  if (timeBasis$type == "polynomial" && fitOutcomeModelArgs$armByTimeInteraction) {
    return(character(0))
  }
  return(intersect(coefficientNames, c("target", "death:target")))
}

# The source person of each copy made by .resampleCloneCensorWeightData(personSeqIds, counts), which
# numbers the copies by source personSeqId and copy: a tibble `personSeqId, sourcePersonSeqId`.
.copySourcePersonSeqIds <- function(personSeqIds, counts) {
  order <- order(personSeqIds)
  sources <- rep(as.integer(personSeqIds[order]), counts[order])
  return(tibble(personSeqId = seq_along(sources), sourcePersonSeqId = sources))
}

# A new CloneCensorWeightData in which the person with personSeqIds[i] appears counts[i] times, each
# copy with new personSeqIds (numbered by person and copy) and row IDs (numbered by copy and source
# row), and the outcomes, treatments and covariates of its rows. The caller closes it.
.resampleCloneCensorWeightData <- function(cloneCensorWeightData, personSeqIds, counts) {
  drawn <- counts > 0
  resampled <- Andromeda::andromeda(
    copies = tibble(personSeqId = as.integer(personSeqIds[drawn]), copies = as.integer(counts[drawn]))
  )
  for (table in c("cohorts", "outcomes", "treatments", "covariates")) {
    resampled[[paste0("source", table)]] <- cloneCensorWeightData[[table]]
  }
  .andromedaExecute(resampled, "
    CREATE TABLE personCopies AS
    SELECT personSeqId,
      copy,
      CAST(ROW_NUMBER() OVER (ORDER BY personSeqId, copy) AS INTEGER) AS newPersonSeqId
    FROM (
      SELECT personSeqId,
        UNNEST(range(copies)) AS copy
      FROM copies
    );")
  .andromedaExecute(resampled, "
    CREATE TABLE rowIdMap AS
    SELECT sourceCohorts.rowId,
      personCopies.newPersonSeqId,
      CAST(ROW_NUMBER() OVER (ORDER BY personCopies.newPersonSeqId, sourceCohorts.rowId) AS INTEGER) AS newRowId
    FROM sourceCohorts
    INNER JOIN personCopies
      ON sourceCohorts.personSeqId = personCopies.personSeqId;")
  .andromedaExecute(resampled, "
    CREATE TABLE cohorts AS
    SELECT rowIdMap.newRowId AS rowId,
      rowIdMap.newPersonSeqId AS personSeqId,
      sourceCohorts.* EXCLUDE (rowId, personSeqId)
    FROM sourceCohorts
    INNER JOIN rowIdMap
      ON sourceCohorts.rowId = rowIdMap.rowId
    ORDER BY rowIdMap.newRowId;")
  for (table in c("outcomes", "treatments", "covariates")) {
    .andromedaExecute(resampled, sprintf("
      CREATE TABLE %1$s AS
      SELECT rowIdMap.newRowId AS rowId,
        source%1$s.* EXCLUDE (rowId)
      FROM source%1$s
      INNER JOIN rowIdMap
        ON source%1$s.rowId = rowIdMap.rowId
      ORDER BY rowIdMap.newRowId;", table))
  }
  for (table in c("copies", "personCopies", "rowIdMap", "sourcecohorts", "sourceoutcomes", "sourcetreatments", "sourcecovariates")) {
    .andromedaExecute(resampled, sprintf("DROP TABLE %s;", table))
  }
  resampled$covariateRef <- cloneCensorWeightData$covariateRef
  resampled$analysisRef <- cloneCensorWeightData$analysisRef

  metaData <- attr(cloneCensorWeightData, "metaData")
  metaData$populationSize <- .andromedaSql(resampled, "SELECT COUNT(*) AS n FROM cohorts;")$n
  attr(resampled, "metaData") <- metaData
  class(resampled) <- "CloneCensorWeightData"
  attr(class(resampled), "package") <- "CloneCensorWeight"
  return(resampled)
}

# The estimates table of the bootstrap: the full-data estimates `t0` (estimand, estimate, logEstimate)
# with the standard errors, limits and p-values from the replicates (one column per estimand). Ratios
# are summarized on the log scale.
.bootstrapEstimates <- function(t0, replicates, ci) {
  rows <- list()
  for (i in seq_len(nrow(t0))) {
    estimand <- t0$estimand[i]
    estimate <- t0$estimate[i]
    t <- replicates[, estimand]
    if (grepl("^(hr|or|rr_.*)$", estimand)) {
      logEstimate <- log(estimate)
      logT <- log(t)
      seLogEstimate <- .bootstrapSe(logT)
      limits <- exp(.bootstrapLimits(logEstimate, logT, ci))
      rows[[i]] <- .estimateRow(
        estimand = estimand,
        estimate = estimate,
        logEstimate = logEstimate,
        seLogEstimate = seLogEstimate,
        lb95 = limits[1],
        ub95 = limits[2],
        p = 2 * stats::pnorm(-abs(logEstimate / seLogEstimate))
      )
    } else {
      se <- .bootstrapSe(t)
      limits <- .bootstrapLimits(estimate, t, ci)
      isDifference <- grepl("^(rd|rmst|rmtl)_[0-9]+$", estimand)
      rows[[i]] <- .estimateRow(
        estimand = estimand,
        estimate = estimate,
        se = se,
        lb95 = limits[1],
        ub95 = limits[2],
        p = if (isDifference) 2 * stats::pnorm(-abs(estimate / se)) else NA_real_
      )
    }
  }
  return(dplyr::bind_rows(rows))
}

# Whether the replicates sit predominantly on one side of the full-data estimate (Phase1Contracts section
# 11.6), for the primary estimands of `t0` (estimand, estimate): the hazard ratio (on the log scale) and
# the risk differences at the horizons. For each, from its finite replicates: the numbers above and below
# the estimate (ties left out), the two-sided sign test p-value of those counts, the median replicate
# minus the estimate in bootstrap standard deviations (NA when the replicates do not vary), and
# `shifted`, TRUE when the median is more than 0.5 standard deviations from the estimate, or more than
# 0.25 with a sign test p < 0.01 (the sign test alone flags negligible shifts with many replicates),
# assessed only with at least 20 finite replicates (NA otherwise).
.bootstrapReplicateShift <- function(t0, replicates) {
  estimands <- t0$estimand[grepl("^(hr|shr|rd_[0-9]+)$", t0$estimand)]
  rows <- list()
  for (estimand in estimands) {
    estimate <- t0$estimate[t0$estimand == estimand]
    t <- replicates[, estimand]
    if (estimand %in% c("hr", "shr")) {
      estimate <- log(estimate)
      t <- log(t)
    }
    t <- t[is.finite(t)]
    above <- sum(t > estimate)
    below <- sum(t < estimate)
    signTestP <- if (is.finite(estimate) && above + below > 0) stats::binom.test(above, above + below)$p.value else NA_real_
    sd <- if (length(t) >= 2) stats::sd(t) else NA_real_
    medianShift <- if (is.finite(estimate) && !is.na(sd) && sd > 0) (stats::median(t) - estimate) / sd else NA_real_
    shifted <- NA
    if (length(t) >= 20 && is.finite(estimate)) {
      shifted <- isTRUE(abs(medianShift) > 0.5) || (isTRUE(signTestP < 0.01) && isTRUE(abs(medianShift) > 0.25))
    }
    rows[[length(rows) + 1]] <- tibble(
      estimand = estimand,
      replicates = length(t),
      above = above,
      below = below,
      signTestP = signTestP,
      medianShift = medianShift,
      shifted = shifted
    )
  }
  if (length(rows) == 0) {
    return(tibble(
      estimand = character(0), replicates = integer(0), above = integer(0), below = integer(0),
      signTestP = numeric(0), medianShift = numeric(0), shifted = logical(0)
    ))
  }
  return(dplyr::bind_rows(rows))
}

.bootstrapSe <- function(t) {
  t <- t[is.finite(t)]
  if (length(t) < 2) {
    return(NA_real_)
  }
  return(stats::sd(t))
}

# The 95% limits of boot::boot.ci() from the full-data estimate t0 and the replicates t (non-finite
# replicates are left out). "percentile" (type "perc"): the (R + 1) alpha order statistics, alpha =
# 0.025 and 0.975, interpolated on the normal quantile scale when (R + 1) alpha is not an integer
# (Davison and Hinkley 1997, equation 5.8), and the smallest or largest replicate when (R + 1) alpha is
# below 1 or above R. "normal" (type "norm"): t0 - bias -+ z sd(t), with bias = mean(t) - t0.
.bootstrapLimits <- function(t0, t, ci, conf = 0.95) {
  t <- t[is.finite(t)]
  nReplicates <- length(t)
  if (nReplicates < 2 || !is.finite(t0)) {
    return(c(NA_real_, NA_real_))
  }
  if (ci == "normal") {
    bias <- mean(t) - t0
    margin <- stats::sd(t) * stats::qnorm((1 + conf) / 2)
    return(c(t0 - bias - margin, t0 - bias + margin))
  }
  alpha <- (1 + c(-conf, conf)) / 2
  rank <- (nReplicates + 1) * alpha
  k <- trunc(rank)
  sorted <- sort(t)
  limits <- numeric(2)
  for (i in 1:2) {
    if (k[i] == 0) {
      limits[i] <- sorted[1]
    } else if (k[i] >= nReplicates) {
      limits[i] <- sorted[nReplicates]
    } else if (k[i] == rank[i]) {
      limits[i] <- sorted[k[i]]
    } else {
      # Interpolation between the k-th and (k + 1)-th order statistics on the normal quantile scale:
      z <- stats::qnorm(c(alpha[i], k[i] / (nReplicates + 1), (k[i] + 1) / (nReplicates + 1)))
      limits[i] <- sorted[k[i]] + (z[1] - z[2]) / (z[3] - z[2]) * (sorted[k[i] + 1] - sorted[k[i]])
    }
  }
  return(limits)
}

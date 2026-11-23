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

# Stacked sandwich variance (varianceEstimator = "stackedSandwich").
#
# The outcome estimates and the censoring models (with their stabilization numerators) are one joint
# M-estimator, clustered on person. For estimates e with infinitesimal-jackknife influence IJ (Variance.R),
# person i's influence is
#   IF_i = IJ_i + sum_m psi_im' c_m,   c_m = A_m^- D_m,   D_m = sum_k IJ_k dlog w_k / d alpha_m,
# where IJ_k = d est / d log w_k is the influence of segment k, psi_im the centered score of model m (its
# rows' U_r - mean(U_m), summed over the person's rows: a person with several index dates, firstOnly =
# FALSE, has a row per entry; 0 for persons outside the persons it was fit on) and A_m its information (plus
# diag(s_j^2 / sigma^2) for a normal prior). The variance is assembled as the quadratic form
#   V = IJ'IJ + C'G + G'C + C'BC,   G_m = sum_i psi_im IJ_i',   B = sum_i psi_i psi_i',
# which equals sum_i IF_i IF_i'. A, B and the scores depend only on the censoring models and come from
# their stored predictions and the covariates (the Cyclops fits are not kept); each set of estimates
# needs only D and G.
#
# Each model is a "block" on a row x interval grid (one interval for a binary model): the rows (rowId, an
# entry of a person) of its risk set, their fitted probabilities, rows at risk in the fit and responses
# give the scores and the information; segments are matched to the block's rows by rowId. A parameter
# is a column of the covariates (`x`, the first column the intercept) times a function of the interval (`timeFunctions`, the first column 1): its value on row (i, l) is
# x[i, a] * timeFunctions[l, t]. A block's `terms` give its share of the log weights: per term, `kappa`
# (the derivative of the log weight per interval, on the scale of the linear predictor),
# `cumulativeLog` (the log weight up to each interval) and `lastInterval` (per segment day, the last
# interval the term sums over, -1 for none). So segment k's derivative is sum over the terms of
# sum_{l <= lastInterval(d_k)} kappa_il z_il.
#
# The weight scales (Phase4Designs.md, stacked-sandwich, design section 2): with g the unstabilized log
# weight, log n the log numerator and sigma = 1 when computeWeightsArgs$stabilized,
# - `weight` (Kaplan-Meier, Aalen-Johansen, and the MSMs when msmWeights is of its kind):
#   dlog w = 1[not truncated] (dg + sigma dlog n);
# - `otherScaleWeight` truncated on its own scale (truncateMsmWeights = TRUE):
#   dlog w = 1[not truncated on its cap] (dg + (1 - sigma) dlog n);
# - `otherScaleWeight` formed from the truncated `weight` (truncateMsmWeights = FALSE):
#   dlog w = 1[weight not truncated] (dg + sigma dlog n) + (1 - 2 sigma) dlog n.
# The caps of metaData$weights$caps are held fixed (exact for fixed caps, conditional for percentiles): a
# segment is truncated when its weight is the cap, that is when its untruncated weight, which the blocks
# reproduce, exceeds it.
#
# Priors: none and a normal prior at its (fixed or cross-validated) variance are exact. A Laplace prior is
# conditional on its support (the nonzero coefficients, plus the intercept and the time basis) and on
# lambda: the centering absorbs lambda sign(beta) s / N and an inexact optimum (diagnostic kktGap).
# Covariates enter on their raw scale: the influence does not depend on the rescaling of the covariates by
# FeatureExtraction's tidying (a linear reparameterization), except for the prior's curvature, which
# uses s_j, the covariate's largest value over the model's risk set (the normalization of the tidying).
#
# Censoring models fit on a uniform sample of their risk set (hazard models, maxCohortSizeForFitting) have
# scores only for the persons in the sample, and their information is that of the sample. Binary models
# that CohortMethod::createPs() fit on a sample, and indicator intervals that a sample leaves without rows
# (which get a plug-in hazard of the whole risk set), are not supported.

# Generalized inverse ------------------------------------------------------------------------------

# Generalized inverse of a positive semi-definite information matrix, as Cyclops' .invertInformation():
# a pivoted Cholesky decomposition of its correlation form gives the inverse on a maximal set of
# linearly independent coordinates (in pivot order), and the other coordinates get 0. Scores and
# derivatives lie in the row space, so U' A^- D does not depend on which inverse is used. That needs
# the rank to be right: an exactly singular information (a time basis collinear with the intercept over
# the intervals with rows) leaves rounding-level pivots (about 1e-15) that LAPACK's default tolerance
# keeps, so pivots below 1e-10 (on the correlation scale) count as 0. Returns the `inverse` and
# the `rank`.
.generalizedInverse <- function(information) {
  n <- nrow(information)
  inverse <- matrix(0, n, n)
  informative <- which(diag(information) > 0)
  if (length(informative) == 0) {
    return(list(inverse = inverse, rank = 0L))
  }
  scale <- 1 / sqrt(diag(information)[informative])
  factor <- suppressWarnings(chol(information[informative, informative, drop = FALSE] * outer(scale, scale), pivot = TRUE, tol = 1e-10))
  rank <- attr(factor, "rank")
  kept <- attr(factor, "pivot")[seq_len(rank)]
  inverse[informative[kept], informative[kept]] <- chol2inv(factor[seq_len(rank), seq_len(rank), drop = FALSE]) * outer(scale[kept], scale[kept])
  return(list(inverse = inverse, rank = rank))
}

# The condition number (largest over smallest eigenvalue) of an information matrix on its correlation
# scale (so it does not depend on the scale of the covariates); Inf when it is singular.
.conditionNumber <- function(information) {
  informative <- diag(information) > 0
  if (!any(informative)) {
    return(NA_real_)
  }
  scale <- 1 / sqrt(diag(information)[informative])
  values <- eigen(information[informative, informative, drop = FALSE] * outer(scale, scale), symmetric = TRUE, only.values = TRUE)$values
  return(if (min(values) <= 0) Inf else max(values) / min(values))
}

# Censoring blocks ------------------------------------------------------------------------------------

# The blocks of the censoring models of the outcome data `data` (.getOutcomeData()), with the covariates
# of `cloneCensorWeightData`: per strategy, its censoring models and their numerator models. The persons
# of `weightedClones` are collected only for the binary model at G with graceRegime = "naturalGrace",
# whose numerator is the share of starters of these persons. Returns the `blocks`, the
# `persons` (personSeqId) of all their populations, their centered scores `psi` (persons x parameters,
# the blocks' columns in order), the block diagonal generalized inverse `inverse` of their information,
# `B` = psi' psi, the `weights` settings (stabilized, truncateMsmWeights, caps), and the `diagnostics` of
# every model. Stops for the settings this version does not support, and when the weights these blocks
# give are not those of the segments.
.stackedCensoringBlocks <- function(data, cloneCensorWeightData, weightedClones) {
  censoringModels <- data$weightsMetaData$censoringModels
  if (is.null(censoringModels)) {
    stop("varianceEstimator = \"stackedSandwich\" needs the censoring models of the weights: rerun computeCensoringWeights()")
  }
  context <- list(
    sameDayOrdering = data$sameDayOrdering,
    offset = if (data$sameDayOrdering == "eventsFirst") 1L else 0L,
    prior = .stackedPrior(censoringModels),
    getPersons = function() collect(weightedClones$persons)
  )
  blocks <- list()
  for (model in censoringModels) {
    blocks <- c(blocks, .modelStackedBlocks(model, cloneCensorWeightData, context))
  }
  persons <- sort(unique(unlist(lapply(blocks, function(block) block$personSeqId))))
  sizes <- vapply(blocks, function(block) ncol(block$psi), integer(1))
  psi <- matrix(0, length(persons), sum(sizes))
  offset <- 0
  for (block in blocks) {
    # The scores of a person's entries (firstOnly = FALSE) are summed:
    summed <- rowsum(block$psi, block$personSeqId)
    psi[match(as.numeric(rownames(summed)), persons), offset + seq_len(ncol(block$psi))] <- summed
    offset <- offset + ncol(block$psi)
  }
  caps <- data$weightsMetaData$caps
  result <- list(
    blocks = lapply(blocks, function(block) block[!names(block) %in% c("psi", "diagnostics")]),
    persons = persons,
    psi = psi,
    inverse = .blockDiagonal(lapply(blocks, function(block) block$inverse)),
    B = crossprod(psi),
    weights = list(
      stabilized = isTRUE(data$weightsMetaData$stabilized),
      truncateMsmWeights = !isFALSE(data$weightsMetaData$truncateMsmWeights),
      caps = caps
    )
  )
  result$diagnostics <- .stackedModelDiagnostics(blocks, caps)
  .checkStackedWeights(result, data)
  return(result)
}

# The prior of the censoring models (attr(censoringModels, "metaData")): its type and the covariates it
# excludes. Models saved without these settings have no prior when none has a prior variance.
.stackedPrior <- function(censoringModels) {
  prior <- attr(censoringModels, "metaData")$fitCensoringModelsArgs$prior
  if (is.null(prior$priorType)) {
    if (all(vapply(censoringModels, function(model) is.na(model$priorVariance), logical(1)))) {
      return(list(type = "none", exclude = NULL))
    }
    stop("The prior of the censoring models is not known: refit the censoring models with this version to use varianceEstimator = \"stackedSandwich\"")
  }
  if (!prior$priorType %in% c("none", "normal", "laplace")) {
    stop(sprintf("varianceEstimator = \"stackedSandwich\" supports censoring models with no prior, a normal or a Laplace prior, not a \"%s\" prior", prior$priorType))
  }
  return(list(type = prior$priorType, exclude = as.character(prior$exclude)))
}

# The blocks of one strategy's censoring model (and of the models its weights also use: `early`,
# `postGrace`), by strategy type (see computeCensoringWeights() for the weights):
# - "initiateWithin": after the switch at G (cutDay), -log P with P = E(X) p: E(X) = S(graceStart - 1 | X)
#   from the prefix of the hazard model of the days before graceStart (its weights before the switch),
#   and p from the binary model at G (-log pi) or from the hazard model (-log(1 - S(G) / S(firstDay - 1)),
#   for the persons at risk of deviating at G);
# - "notInitiateWithin": the prefix of the hazard model through G, and with postGraceModel = "binary",
#   -log(1 - pi(X)) from the binary model at G in its place from evalDay G on;
# - "neverInitiate": the prefix of the hazard model.
.modelStackedBlocks <- function(model, cloneCensorWeightData, context) {
  strategyId <- as.integer(model$strategyId)
  .checkStackedFit(model, strategyId, "censoring model")
  offset <- context$offset
  if (identical(model$type, "initiateWithin")) {
    graceStart <- if (is.null(model$graceStart)) 0L else as.integer(model$graceStart)
    graceRegime <- if (is.null(model$graceRegime)) "binaryAtG" else model$graceRegime
    graceEnd <- as.integer(model$horizon)
    cutDay <- graceEnd - if (context$sameDayOrdering == "eventsFirst") 0L else 1L
    firstDay <- .binaryModelFirstDay(list(graceStart = graceStart, graceEnd = graceEnd), graceRegime)
    if (model$kind == "binary") {
      numeratorPopulation <- NULL
      if (graceRegime == "naturalGrace") {
        # The numerator is the share of starters by G of the persons (of these clones) under follow-up at G
        # who have not started before graceStart (.initiateWithinWeightModel()):
        persons <- context$getPersons()
        treatmentDay <- if (model$treatment == "treatment") persons$treatmentDay else persons$comparatorTreatmentDay
        numeratorPopulation <- .getBinaryRiskSet(persons, treatmentDay, graceEnd, context$sameDayOrdering, graceStart)
      }
      blocks <- .binaryStackedBlocks(model, strategyId, "binary model at G", cloneCensorWeightData, context,
                                     function(day) ifelse(day > cutDay, 0L, -1L), complement = FALSE, numeratorPopulation)
      if (graceStart > 0) {
        .checkStackedFit(model$early, strategyId, "hazard model of the days before graceStart")
        blocks <- c(blocks, .logisticStackedBlocks(model$early, strategyId, "hazard model before graceStart", cloneCensorWeightData, context,
                                                   prefix = list(windowEnd = graceStart - 1L)))
      }
      return(blocks)
    }
    return(.logisticStackedBlocks(model, strategyId, "hazard model", cloneCensorWeightData, context,
                                  prefix = if (graceStart > 0) list(windowEnd = graceStart - 1L),
                                  following = list(firstDay = firstDay, numeratorFirstDay = graceStart, cutDay = cutDay)))
  }
  if (!is.null(model$postGrace)) {
    .checkStackedFit(model$postGrace, strategyId, "binary model at G of postGraceModel = \"binary\"")
    graceEnd <- as.integer(model$horizon)
    return(c(
      .logisticStackedBlocks(model, strategyId, "hazard model", cloneCensorWeightData, context,
                             prefix = list(windowEnd = graceEnd, stopAt = graceEnd)),
      .binaryStackedBlocks(model$postGrace, strategyId, "binary model at G (postGraceModel)", cloneCensorWeightData, context,
                           function(day) ifelse(day - offset >= graceEnd, 0L, -1L), complement = TRUE)
    ))
  }
  return(.logisticStackedBlocks(model, strategyId, "hazard model", cloneCensorWeightData, context, prefix = list(windowEnd = as.integer(model$horizon))))
}

# Stops for a censoring model (or a model its weights also use, `label`) that the stacked sandwich does
# not support.
.checkStackedFit <- function(model, strategyId, label) {
  if (model$modelType == "cox") {
    stop(sprintf("varianceEstimator = \"stackedSandwich\" does not yet support the Cox censoring model (hazardModelType = \"cox\", strategy %d)", strategyId))
  }
  columns <- c("personSeqId", "y", "inFit", if (model$modelType == "logistic") "lastDay")
  if (!all(columns %in% colnames(model$predictions))) {
    stop(sprintf("The %s of strategy %d was fit by an earlier version: refit the censoring models with this version to use varianceEstimator = \"stackedSandwich\"", label, strategyId))
  }
  if (model$modelType == "binary" && (isTRUE(model$diagnostics$sampled) || anyNA(model$predictions$inFit))) {
    stop(sprintf(
      "The %s of strategy %d was fit by CohortMethod::createPs() on a sample (maxCohortSizeForFitting), which varianceEstimator = \"stackedSandwich\" does not support: set maxCohortSizeForFitting = 0, or use varianceEstimator = \"bootstrap\"",
      label, strategyId
    ))
  }
  invisible(NULL)
}

# The covariates `covariateIds` (raw values) of the persons `rowIds`: a dense persons x covariates matrix.
.stackedCovariateMatrix <- function(cloneCensorWeightData, rowIds, covariateIds) {
  x <- matrix(0, length(rowIds), length(covariateIds))
  if (length(covariateIds) == 0) {
    return(x)
  }
  ids <- as.numeric(covariateIds)
  covariates <- cloneCensorWeightData$covariates |>
    filter(.data$covariateId %in% ids) |>
    collect()
  person <- match(covariates$rowId, rowIds)
  keep <- !is.na(person)
  x[cbind(person[keep], match(covariates$covariateId[keep], ids))] <- covariates$covariateValue[keep]
  return(x)
}

# The coefficients of a model that are its parameters: all of them, except, with a Laplace prior, the
# penalized ones at 0 (the support is held fixed). As in the fit (.restrictPrior()), the prior excludes
# covariates, not their interactions with time, which are penalized.
.stackedSupport <- function(coefficients, prior, unpenalized) {
  if (prior$type != "laplace") {
    return(names(coefficients))
  }
  penalized <- !names(coefficients) %in% c(unpenalized, prior$exclude)
  return(names(coefficients)[!penalized | coefficients != 0])
}

# The scale s_j (the largest value over the model's risk set, the normalization of FeatureExtraction's
# tidying) of each covariate column of `x`.
.covariateScales <- function(x) {
  return(apply(x, 2, function(values) if (any(values != 0)) max(values[values != 0]) else 1))
}

# The prior of a model's parameters (NULL without one): its type, variance, and per parameter the scale of
# its covariate (`scale`, NA for an unpenalized parameter). A model fit with a prior keeps it when none of
# its parameters is penalized (a Laplace prior that set every covariate to 0 is still conditional on
# that support).
.stackedPenalty <- function(prior, variance, scale) {
  if (prior$type == "none" || is.na(variance)) {
    return(NULL)
  }
  return(list(type = prior$type, variance = variance, scale = scale))
}

# A block from its population (the `personSeqId` and `rowId` of each row), covariates `x` (with the
# intercept column), time functions, parameters (`parameterX`, `parameterTime`: the column of each), fitted `probability`, rows
# `atRisk` in the fit and responses `y` (persons x intervals), whether each person is in the fit
# (`inFit`), the prior (`penalty`, .stackedPenalty()) and the `terms` (see the top of the file). The
# scores, centered over the persons in the fit, and the generalized inverse of the information are
# computed here, with the model's diagnostics; the probabilities, rows and responses are not kept.
.stackedBlock <- function(strategyId, numerator, component, personSeqId, rowId, x, timeFunctions, parameterX, parameterTime,
                          probability, atRisk, y, terms, inFit = rep(TRUE, length(personSeqId)), penalty = NULL) {
  residuals <- (y - probability) * atRisk
  scores <- x[, parameterX, drop = FALSE] * (residuals %*% timeFunctions)[, parameterTime, drop = FALSE]
  # The information sum_il c_il z_il z_il', assembled by kind of parameter (no persons x parameters
  # product per interval). The parameters of the intercept column x[, 1] = 1 (the intercept and the time
  # basis) have z_il = timeFunctions[l, t]; the others (covariates and their interactions) are covariate
  # columns times a time function.
  curvature <- probability * (1 - probability) * atRisk
  information <- matrix(0, length(parameterX), length(parameterX))
  time <- which(parameterX == 1)
  timeColumns <- timeFunctions[, parameterTime[time], drop = FALSE]
  information[time, time] <- crossprod(timeColumns, timeColumns * colSums(curvature))
  covariate <- which(parameterX != 1)
  if (length(covariate) > 0) {
    # sum_i x_ij c_il per covariate column and interval, times the time functions:
    columns <- unique(parameterX[covariate])
    covariateCurvature <- crossprod(x[, columns, drop = FALSE], curvature)[match(parameterX[covariate], columns), , drop = FALSE]
    information[covariate, time] <- (covariateCurvature * t(timeFunctions[, parameterTime[covariate], drop = FALSE])) %*% timeColumns
    information[time, covariate] <- t(information[covariate, time, drop = FALSE])
    # Per pair of time functions (s, t) of the covariate parameters, X_s' diag(sum_l c_il f_ls f_lt) X_t
    # (one X' diag(sum_l c_il) X for the covariates):
    functions <- unique(parameterTime[covariate])
    for (s in functions) {
      for (t in functions) {
        a <- covariate[parameterTime[covariate] == s]
        b <- covariate[parameterTime[covariate] == t]
        weight <- as.vector(curvature %*% (timeFunctions[, s] * timeFunctions[, t]))
        information[a, b] <- crossprod(x[, parameterX[a], drop = FALSE], x[, parameterX[b], drop = FALSE] * weight)
      }
    }
  }
  diagnostics <- list(prior = "none", supportSize = length(parameterX), conditionalOnSupport = FALSE, maxKktGap = NA_real_)
  if (!is.null(penalty)) {
    penalized <- !is.na(penalty$scale)
    diagnostics$prior <- penalty$type
    if (penalty$type == "normal") {
      information <- information + diag(ifelse(penalized, penalty$scale^2 / penalty$variance, 0), length(parameterX))
    } else {
      # Laplace, conditional on the support: |sum_i U_ij| / s_j equals lambda at the optimum (on the
      # normalized scale of the fit) for the penalized parameters in it (none: maxKktGap NA).
      diagnostics$conditionalOnSupport <- TRUE
      if (any(penalized)) {
        lambda <- sqrt(2 / penalty$variance)
        gap <- abs(colSums(scores)[penalized]) / penalty$scale[penalized] - lambda
        diagnostics$maxKktGap <- max(abs(gap))
      }
    }
  }
  psi <- matrix(0, nrow(scores), ncol(scores))
  if (any(inFit)) {
    psi[inFit, ] <- sweep(scores[inFit, , drop = FALSE], 2, colMeans(scores[inFit, , drop = FALSE]))
  }
  inverse <- .generalizedInverse(information)
  diagnostics$rank <- inverse$rank
  diagnostics$conditionNumber <- .conditionNumber(information)
  diagnostics$sampled <- !all(inFit)
  return(list(
    strategyId = as.integer(strategyId),
    numerator = numerator,
    component = if (numerator) paste(component, "numerator") else component,
    personSeqId = personSeqId,
    rowId = rowId,
    x = x,
    timeFunctions = timeFunctions,
    parameterX = parameterX,
    parameterTime = parameterTime,
    terms = terms,
    psi = psi,
    inverse = inverse$inverse,
    diagnostics = diagnostics
  ))
}

# A term (see the top of the file); non-finite derivatives (of persons whose weights it does not give)
# are 0.
.stackedTerm <- function(kappa, cumulativeLog, lastInterval) {
  kappa[!is.finite(kappa)] <- 0
  return(list(kappa = kappa, cumulativeLog = cumulativeLog, lastInterval = lastInterval))
}

# Cumulative sums of the columns of a matrix, along each row.
.rowCumsum <- function(m) {
  for (l in seq_len(ncol(m))[-1]) {
    m[, l] <- m[, l] + m[, l - 1]
  }
  return(m)
}

# log(1 + exp(eta)), without overflow.
.logOnePlusExp <- function(eta) {
  return(ifelse(eta > 30, eta, log1p(exp(eta))))
}

# A binary model at G (fit by CohortMethod::createPs() on the whole risk set) and its numerator, whose
# weights apply on the days `lastInterval` gives 0: with `complement = FALSE` (initiateWithin), W = 1 /
# pi(X) (pi-bar / pi(X) stabilized), log weight -log(pi) with derivative -(1 - pi) x and log numerator
# log(pi-bar) with derivative 1 - pi-bar; with `complement = TRUE` (notInitiateWithin, postGraceModel =
# "binary"), W = 1 / (1 - pi(X)), log weight -log(1 - pi) with derivative pi x and log numerator
# log(1 - pi-bar) with derivative -pi-bar. The numerator pi-bar is the share of starters of the model's
# risk set, or of `numeratorPopulation` (personSeqId, y) when given.
.binaryStackedBlocks <- function(model, strategyId, component, cloneCensorWeightData, context, lastInterval, complement,
                                 numeratorPopulation = NULL) {
  predictions <- model$predictions
  n <- nrow(predictions)
  names <- .stackedSupport(model$coefficients, context$prior, "(Intercept)")
  covariateIds <- setdiff(names, "(Intercept)")
  covariates <- .stackedCovariateMatrix(cloneCensorWeightData, predictions$rowId, covariateIds)
  x <- cbind(1, covariates)
  probability <- matrix(predictions$propensityScore)
  one <- matrix(1, 1, 1)
  scale <- c(NA, ifelse(covariateIds %in% context$prior$exclude, NA, .covariateScales(covariates)))
  term <- if (complement) {
    .stackedTerm(probability, -log(1 - probability), lastInterval)
  } else {
    .stackedTerm(-(1 - probability), -log(probability), lastInterval)
  }
  denominator <- .stackedBlock(strategyId, FALSE, component, predictions$personSeqId, predictions$rowId, x, one, seq_len(ncol(x)), rep(1L, ncol(x)),
                               probability, matrix(1, n, 1), matrix(predictions$y), list(term),
                               penalty = .stackedPenalty(context$prior, model$priorVariance, scale))
  if (is.null(numeratorPopulation)) {
    numeratorPopulation <- predictions
    marginal <- model$numerator$propensityScore
  } else {
    marginal <- mean(numeratorPopulation$y)
  }
  m <- nrow(numeratorPopulation)
  marginal <- matrix(marginal, m, 1)
  term <- if (complement) {
    .stackedTerm(-marginal, log(1 - marginal), lastInterval)
  } else {
    .stackedTerm(1 - marginal, log(marginal), lastInterval)
  }
  numerator <- .stackedBlock(strategyId, TRUE, component, numeratorPopulation$personSeqId, numeratorPopulation$rowId, matrix(1, m, 1), one, 1L, 1L,
                             marginal, matrix(1, m, 1), matrix(numeratorPopulation$y), list(term))
  return(list(denominator, numerator))
}

# The pooled logistic initiation hazard model of a strategy and its numerator model (intercept and time
# basis, the same for every person, fit without a prior). Rows: the intervals of the fit (not those with a
# hazard of 0 or 1) up to each person's last interval, for the persons in the fit. Terms:
# - `prefix` (windowEnd, and stopAt: no prefix from the day evalDay reaches it): W(d) = 1 / S(e | X) with
#   e = min(d - offset, windowEnd) (S-bar / S stabilized): log weight -sum log(1 - h) with derivative
#   sum h z over the intervals that end by e, and log numerator sum log(1 - h-bar) with derivative
#   -sum h-bar z;
# - `following` (firstDay, numeratorFirstDay, cutDay): after cutDay, the probability of starting on days
#   firstDay..G of the persons at risk of deviating at G (their last day at risk is on or after firstDay),
#   1 - Q with Q = S(G) / S(firstDay - 1): log weight -log(1 - Q) with derivative -Q / (1 - Q) h z over
#   the intervals that end on days firstDay..G, and log numerator log(1 - Q-bar) (from numeratorFirstDay)
#   with derivative Q-bar / (1 - Q-bar) h-bar z.
.logisticStackedBlocks <- function(model, strategyId, component, cloneCensorWeightData, context, prefix = NULL, following = NULL) {
  predictions <- model$predictions
  n <- nrow(predictions)
  nIntervals <- as.integer(model$timeBasis$nIntervals)
  intervalLength <- model$timeBasis$intervalLength
  horizon <- model$horizon
  basis <- model$timeBasis$basis
  interactionBasis <- model$interactions$basis
  k <- seq_len(nIntervals) - 1L
  endDay <- pmin((k + 1L) * intervalLength - 1L, horizon)
  fitted <- is.finite(model$intervals$baseline)
  fittedRows <- matrix(fitted, n, nIntervals, byrow = TRUE)
  inFit <- predictions$inFit %in% TRUE
  personLast <- predictions$lastDay %/% intervalLength
  reached <- outer(personLast, k, ">=")
  initiated <- outer(personLast, k, "==") * (predictions$y == 1)
  if (!all(inFit) && identical(model$timeBasis$settings$type, "indicator")) {
    # An indicator interval with rows in the risk set but a hazard of 0 or 1 in the sample gets a plug-in
    # hazard of the whole risk set, which is not an estimating equation of the sample:
    sampleInitiations <- colSums(initiated[inFit, , drop = FALSE])
    sampleAtRisk <- colSums(reached[inFit, , drop = FALSE])
    if (any(fitted & !(sampleInitiations > 0 & sampleInitiations < sampleAtRisk))) {
      stop(sprintf(
        "The %s of strategy %d was fit on a sample (maxCohortSizeForFitting) that leaves indicator intervals without initiations (or with only initiations), whose hazards are plug-ins from the whole risk set; varianceEstimator = \"stackedSandwich\" does not support these: set maxCohortSizeForFitting = 0, use a smooth time basis, or use varianceEstimator = \"bootstrap\"",
        component, strategyId
      ))
    }
  }
  atRisk <- reached * fittedRows * inFit
  y <- atRisk * initiated
  offset <- context$offset

  # The covariate model: intercept, time basis, covariates and covariate-by-time interactions.
  names <- .stackedSupport(model$coefficients, context$prior, c("(Intercept)", colnames(basis)))
  timeColumns <- names[names %in% colnames(basis)]
  isInteraction <- grepl(":", names, fixed = TRUE)
  covariateIds <- names[!names %in% c("(Intercept)", timeColumns) & !isInteraction]
  interactionCovariates <- sub(":.*$", "", names[isInteraction])
  interactionColumns <- sub("^[^:]*:", "", names[isInteraction])
  usedColumns <- unique(interactionColumns)
  allCovariateIds <- unique(c(covariateIds, interactionCovariates))
  covariates <- .stackedCovariateMatrix(cloneCensorWeightData, predictions$rowId, allCovariateIds)
  x <- cbind(1, covariates)
  timeFunctions <- cbind(1, basis[, timeColumns, drop = FALSE], interactionBasis[, match(usedColumns, colnames(interactionBasis)), drop = FALSE])
  parameterX <- c(1L, rep(1L, length(timeColumns)), 1L + match(covariateIds, allCovariateIds), 1L + match(interactionCovariates, allCovariateIds))
  parameterTime <- c(1L, 1L + seq_along(timeColumns), rep(1L, length(covariateIds)), 1L + length(timeColumns) + match(interactionColumns, usedColumns))
  # The interactions of a covariate the prior excludes are penalized (.stackedSupport()), at its scale:
  covariateScale <- .covariateScales(covariates)
  scale <- c(
    rep(NA, 1 + length(timeColumns)),
    ifelse(covariateIds %in% context$prior$exclude, NA, covariateScale[match(covariateIds, allCovariateIds)]),
    covariateScale[match(interactionCovariates, allCovariateIds)]
  )
  eta <- outer(predictions$linearPredictor, model$intervals$baseline, "+")
  interactions <- model$interactions$predictions
  if (!is.null(interactions) && nrow(interactions) > 0) {
    values <- matrix(0, n, ncol(interactionBasis))
    values[cbind(match(interactions$rowId, predictions$rowId), interactions$columnIndex)] <- interactions$value
    eta <- eta + values %*% t(interactionBasis)
  }
  probability <- stats::plogis(eta)
  numeratorEta <- matrix(model$numerator$intervals$baseline, n, nIntervals, byrow = TRUE)
  numeratorProbability <- stats::plogis(numeratorEta)

  terms <- list()
  numeratorTerms <- list()
  if (!is.null(prefix)) {
    # The interval of the weight of day `day` (.evalDaySql() and .computeSurvival()), -1 for none:
    lastInterval <- function(day) {
      evalDay <- pmin(day - offset, prefix$windowEnd)
      interval <- as.integer(ifelse(evalDay >= horizon, nIntervals - 1L, floor((evalDay + 1) / intervalLength) - 1L))
      if (!is.null(prefix$stopAt)) {
        interval[day - offset >= prefix$stopAt] <- -1L
      }
      return(interval)
    }
    terms[[1]] <- .stackedTerm(probability * fittedRows, .rowCumsum(.logOnePlusExp(eta)), lastInterval)
    numeratorTerms[[1]] <- .stackedTerm(-numeratorProbability * fittedRows, -.rowCumsum(.logOnePlusExp(numeratorEta)), lastInterval)
  }
  if (!is.null(following)) {
    lastInterval <- function(day) ifelse(day > following$cutDay, nIntervals - 1L, -1L)
    startTerm <- function(eta, probability, intervals, persons, sign) {
      logQ <- -rowSums(.logOnePlusExp(eta[, intervals, drop = FALSE]))
      logOneMinusQ <- log(-expm1(logQ))
      kappa <- matrix(0, n, nIntervals)
      kappa[, intervals] <- sign * exp(logQ - logOneMinusQ) * probability[, intervals, drop = FALSE]
      kappa <- kappa * fittedRows * persons
      return(.stackedTerm(kappa, matrix(sign * logOneMinusQ * persons, n, nIntervals), lastInterval))
    }
    terms[[length(terms) + 1]] <- startTerm(eta, probability, endDay >= following$firstDay, predictions$lastDay >= following$firstDay, -1)
    numeratorTerms[[length(numeratorTerms) + 1]] <- startTerm(numeratorEta, numeratorProbability, endDay >= following$numeratorFirstDay, rep(TRUE, n), 1)
  }
  denominator <- .stackedBlock(strategyId, FALSE, component, predictions$personSeqId, predictions$rowId, x, timeFunctions, parameterX, parameterTime,
                               probability, atRisk, y, terms, inFit, .stackedPenalty(context$prior, model$priorVariance, scale))
  numeratorColumns <- names(model$numerator$coefficients)
  numeratorColumns <- numeratorColumns[numeratorColumns %in% colnames(basis)]
  numerator <- .stackedBlock(strategyId, TRUE, component, predictions$personSeqId, predictions$rowId, x[, 1, drop = FALSE],
                             cbind(1, basis[, numeratorColumns, drop = FALSE]), rep(1L, 1 + length(numeratorColumns)), seq_len(1 + length(numeratorColumns)),
                             numeratorProbability, atRisk, y, numeratorTerms, inFit)
  return(list(denominator, numerator))
}

# Per model: its strategy, component, prior, support size, whether it is conditional on the support, the
# largest KKT gap, the rank and condition number of its information, whether it was fit on a sample, and
# its strategy's truncated shares on both weight scales (computeCensoringWeights()).
.stackedModelDiagnostics <- function(blocks, caps) {
  value <- function(name, type) vapply(blocks, function(block) block$diagnostics[[name]], type)
  strategyId <- vapply(blocks, function(block) block$strategyId, integer(1))
  pick <- function(column) {
    if (is.null(caps[[column]])) {
      return(rep(NA_real_, length(strategyId)))
    }
    return(as.numeric(caps[[column]][match(strategyId, caps$strategyId)]))
  }
  return(tibble(
    strategyId = strategyId,
    component = vapply(blocks, function(block) block$component, character(1)),
    prior = value("prior", character(1)),
    supportSize = value("supportSize", integer(1)),
    conditionalOnSupport = value("conditionalOnSupport", logical(1)),
    maxKktGap = value("maxKktGap", numeric(1)),
    rank = value("rank", integer(1)),
    conditionNumber = value("conditionNumber", numeric(1)),
    sampled = value("sampled", logical(1)),
    truncatedShare = pick("truncatedShare"),
    otherScaleTruncatedShare = pick("otherScaleTruncatedShare")
  ))
}

# Per block and term, the cells of the `segments` (rowId of the entry, strategyId, lastDay) of the
# block's strategy (`rows`): the block's row of each one's entry (`person`) and its `interval`; a segment
# is `used` when the term gives part of its weight (interval >= 0) and its entry is in the block's
# population (the weights of the others do not depend on the block, which .checkStackedWeights()
# confirms). Blocks of a strategy with the same rows (a model and its numerator) share the match of the
# segments to their rows, and terms with the same lastInterval share the intervals. `stacked$segmentCells`,
# when set (.stackedVcov()), holds them for these segments.
.stackedSegmentCells <- function(stacked, segments) {
  if (!is.null(stacked$segmentCells)) {
    return(stacked$segmentCells)
  }
  matches <- list()
  intervals <- list()
  cells <- vector("list", length(stacked$blocks))
  for (b in seq_along(stacked$blocks)) {
    block <- stacked$blocks[[b]]
    m <- Position(function(entry) entry$strategyId == block$strategyId && identical(entry$rowId, block$rowId), matches)
    if (is.na(m)) {
      rows <- which(segments$strategyId == block$strategyId)
      matches[[length(matches) + 1]] <- list(strategyId = block$strategyId, rowId = block$rowId, rows = rows, person = match(segments$rowId[rows], block$rowId))
      m <- length(matches)
    }
    matched <- matches[[m]]
    cells[[b]] <- vector("list", length(block$terms))
    for (t in seq_along(block$terms)) {
      lastInterval <- block$terms[[t]]$lastInterval
      i <- Position(function(entry) entry$strategyId == block$strategyId && identical(entry$lastInterval, lastInterval), intervals)
      if (is.na(i)) {
        intervals[[length(intervals) + 1]] <- list(strategyId = block$strategyId, lastInterval = lastInterval, interval = .byDay(lastInterval, segments$lastDay[matched$rows]))
        i <- length(intervals)
      }
      interval <- intervals[[i]]$interval
      cells[[b]][[t]] <- list(rows = matched$rows, person = matched$person, interval = interval, used = interval >= 0 & !is.na(matched$person))
    }
  }
  return(cells)
}

# f(days) for a function f of the day (a lastInterval), evaluated once per day of their range.
.byDay <- function(f, days) {
  if (length(days) == 0) {
    return(f(days))
  }
  first <- min(days)
  return(f(seq(first, max(days)))[days - first + 1])
}

# The unstabilized log weight `logWeight` and the log numerator `logNumerator` of each segment
# (rowId, strategyId, lastDay) from the blocks.
.stackedSegmentLogs <- function(stacked, segments) {
  logWeight <- numeric(nrow(segments))
  logNumerator <- numeric(nrow(segments))
  cells <- .stackedSegmentCells(stacked, segments)
  for (b in seq_along(stacked$blocks)) {
    block <- stacked$blocks[[b]]
    for (t in seq_along(block$terms)) {
      cell <- cells[[b]][[t]]
      used <- cell$rows[cell$used]
      values <- block$terms[[t]]$cumulativeLog[cbind(cell$person[cell$used], cell$interval[cell$used] + 1L)]
      if (block$numerator) {
        logNumerator[used] <- logNumerator[used] + values
      } else {
        logWeight[used] <- logWeight[used] + values
      }
    }
  }
  return(list(logWeight = logWeight, logNumerator = logNumerator))
}

# The cap of each segment's strategy (NA: none) from the caps table, column `column`.
.segmentCaps <- function(caps, strategyId, column) {
  if (is.null(caps[[column]])) {
    return(rep(NA_real_, length(strategyId)))
  }
  return(as.numeric(caps[[column]][match(strategyId, caps$strategyId)]))
}

# Per segment, the multipliers of the derivatives of the unstabilized log weight (`denominator`) and of
# the log numerator (`numerator`) for estimates whose weights are stabilized when `stabilized` (the weight
# scale they use; see the top of the file).
.stackedScaleMultipliers <- function(stacked, segments, stabilized) {
  logs <- .stackedSegmentLogs(stacked, segments)
  weights <- stacked$weights
  sigma <- as.numeric(weights$stabilized)
  sameScale <- isTRUE(stabilized) == weights$stabilized
  notTruncated <- function(logWeight, cap) is.na(cap) | logWeight <= log(cap) + 1e-10
  if (sameScale || !weights$truncateMsmWeights) {
    untruncated <- notTruncated(logs$logWeight + sigma * logs$logNumerator, .segmentCaps(weights$caps, segments$strategyId, "cap"))
    numerator <- sigma * untruncated + if (sameScale) 0 else 1 - 2 * sigma
  } else {
    untruncated <- notTruncated(logs$logWeight + (1 - sigma) * logs$logNumerator, .segmentCaps(weights$caps, segments$strategyId, "otherScaleCap"))
    numerator <- (1 - sigma) * untruncated
  }
  return(list(denominator = as.numeric(untruncated), numerator = numerator))
}

# Keys of the segments (rowId, strategyId, lastDay) on which their log weights, multipliers and cells
# depend: the entry (rowId) and a code of the strategy and day class, the days on which every term of the
# strategy's blocks takes the same interval (lastInterval) sharing a class. A list of one key (the entry
# and the code in one number, when that is exact) or two.
.stackedEntryKeys <- function(stacked, rowId, strategyId, lastDay) {
  days <- seq(min(lastDay), max(lastDay))
  strategies <- seq(min(strategyId), max(strategyId))
  codes <- list()
  offset <- 0L
  for (strategy in strategies) {
    intervals <- list(rep(0L, length(days)))
    for (block in stacked$blocks[vapply(stacked$blocks, function(block) block$strategyId == strategy, logical(1))]) {
      for (term in block$terms) {
        intervals[[length(intervals) + 1]] <- term$lastInterval(days)
      }
    }
    key <- do.call(paste, intervals)
    class <- match(key, unique(key))
    codes[[length(codes) + 1]] <- offset + class
    offset <- offset + max(class)
  }
  code <- unlist(codes)[(strategyId - strategies[1]) * length(days) + (lastDay - days[1] + 1)]
  if ((max(rowId) + 1) * offset < 2^53) {
    return(list(rowId * as.numeric(offset) + code))
  }
  return(list(rowId, code))
}

# Classes of consecutive segments with the same value of every key (vectors of the same length; NULL ones
# are skipped; a missing value starts a class): per segment, its `class` (1, 2, ..., in order) and the
# `first` segment of each class.
.stackedClasses <- function(keys) {
  keys <- keys[!vapply(keys, is.null, logical(1))]
  n <- length(keys[[1]])
  change <- logical(max(n - 1, 0))
  for (key in keys) {
    change <- change | key[-1L] != key[-n]
  }
  first <- c(TRUE, change)[seq_len(n)]
  first[is.na(first)] <- TRUE
  return(list(class = cumsum(first), first = which(first)))
}

# Stops when the weights that the blocks give differ from those of the segments (with positive weights)
# of `data`, on either scale (with the caps), which would make the derivatives wrong. Segments of an entry
# with the same keys (.stackedEntryKeys()) and the same weights have the same differences, so each
# run of consecutive such segments is checked once: every segment's difference is checked.
.checkStackedWeights <- function(stacked, data) {
  segments <- data$segments
  lastDay <- segments$stopTime + data$riskStart - 1
  first <- .stackedClasses(c(
    .stackedEntryKeys(stacked, segments$rowId, segments$strategyId, lastDay),
    list(segments$weight, segments$count, segments$stabilizationNumerator, segments$otherScaleWeight)
  ))$first
  first <- first[segments$weight[first] > 0]
  segments <- data.frame(
    rowId = segments$rowId[first], strategyId = segments$strategyId[first], lastDay = lastDay[first], weight = segments$weight[first],
    count = segments$count[first], stabilizationNumerator = segments$stabilizationNumerator[first]
  )
  if (!is.null(data$segments$otherScaleWeight)) {
    segments$otherScaleWeight <- data$segments$otherScaleWeight[first]
  }
  logs <- .stackedSegmentLogs(stacked, segments)
  weights <- stacked$weights
  sigma <- as.numeric(weights$stabilized)
  truncate <- function(logWeight, cap) ifelse(is.na(cap), logWeight, pmin(logWeight, log(cap)))
  logWeight <- truncate(logs$logWeight + sigma * logs$logNumerator, .segmentCaps(weights$caps, segments$strategyId, "cap"))
  differences <- c(logWeight - log(segments$weight / segments$count), logs$logNumerator - log(segments$stabilizationNumerator))
  if (!is.null(segments$otherScaleWeight)) {
    logOther <- if (weights$truncateMsmWeights) {
      truncate(logs$logWeight + (1 - sigma) * logs$logNumerator, .segmentCaps(weights$caps, segments$strategyId, "otherScaleCap"))
    } else {
      logWeight + (1 - 2 * sigma) * logs$logNumerator
    }
    differences <- c(differences, logOther - log(segments$otherScaleWeight))
  }
  if (max(abs(differences)) > 1e-6) {
    stop("The stacked sandwich cannot reproduce the weights from the censoring models: were the weights computed with these censoring models?")
  }
  invisible(NULL)
}

# Derivatives and variance -------------------------------------------------------------------------

# Influences (segments or rows x estimates) are a matrix, or factored (.pooledRowInfluence()): a list
# of the `cell` of each row and `parts`, per model the `residual` of each row and the `loading` (cells x
# estimates), the influence of row r being sum over the parts of residual[r] * loading[cell[r], ]. A
# factored influence is never formed as a rows x estimates array, only by chunks of rows.

# The number of rows and of estimates of an influence.
.influenceCount <- function(influence) {
  return(if (is.matrix(influence)) nrow(influence) else length(influence$cell))
}
.influenceColumns <- function(influence) {
  return(if (is.matrix(influence)) ncol(influence) else ncol(influence$parts[[1]]$loading))
}

# The influence of rows `rows` (a rows x estimates matrix).
.influenceRows <- function(influence, rows) {
  if (is.matrix(influence)) {
    return(influence[rows, , drop = FALSE])
  }
  cell <- influence$cell[rows]
  value <- NULL
  for (part in influence$parts) {
    # (only the columns the part loads on, those of its model)
    columns <- which(colSums(part$loading != 0) > 0)
    term <- part$residual[rows] * part$loading[cell, columns, drop = FALSE]
    if (is.null(value)) {
      value <- matrix(0, length(rows), .influenceColumns(influence))
      value[, columns] <- term
    } else {
      value[, columns] <- value[, columns] + term
    }
  }
  return(value)
}

# The sums of the rows `rows` of an influence matrix by `group` (one per row of `rows`): the groups (`key`,
# in order of appearance) and their sums (`sums`, one row per key).
.influenceGroupSums <- function(influence, rows, group) {
  # (reorder = FALSE: the sums are in the order of unique(group), without parsing their names)
  return(list(key = unique(group), sums = rowsum(influence[rows, , drop = FALSE], group, reorder = FALSE)))
}

# The sums of an influence (rows x estimates, or factored) over the classes `class` (.stackedClasses(),
# one per row) of its rows: a classes x estimates matrix. Each chunk of about `maxEntries` entries of the
# influence is formed once; a class split between chunks gets the sums of its chunks in order.
.influenceClassSums <- function(influence, class, maxEntries = 1e7) {
  n <- length(class)
  if (is.matrix(influence) && class[n] == n) {
    return(influence)
  }
  sums <- matrix(0, class[n], .influenceColumns(influence))
  size <- max(1, floor(maxEntries / ncol(sums)))
  for (start in seq(1, n, by = size)) {
    chunk <- start:min(start + size - 1, n)
    classes <- class[start]:class[chunk[length(chunk)]]
    sums[classes, ] <- sums[classes, ] + rowsum(.influenceRows(influence, chunk), class[chunk], reorder = FALSE)
  }
  return(sums)
}

# Per term of a block, the sums of the segment influences `influence` (a segments x estimates matrix)
# times the scale's `multiplier` of each segment by cell (person and interval) of the term, with the
# block's `cells` (.stackedSegmentCells()): the `sums` (cells x estimates) and the `index` of their
# cells (person + nPersons * interval), or NULL when no segment has a part of its weight from the term.
.stackedCellSums <- function(block, influence, multiplier, cells) {
  nPersons <- length(block$personSeqId)
  return(lapply(seq_along(block$terms), function(k) {
    cell <- cells[[k]]
    if (!any(cell$used)) {
      return(NULL)
    }
    used <- cell$rows[cell$used]
    nCells <- nPersons * ncol(block$terms[[k]]$kappa)
    # The sums of the influence of the used segments with a nonzero multiplier by cell and multiplier
    # (the rows of other strategies and those not used are not visited). The multipliers of the weight
    # scales are 0 or +-1, so scaling the sums by them is exact.
    weight <- multiplier[used]
    nonZero <- weight != 0
    if (!any(nonZero)) {
      return(NULL)
    }
    levels <- unique(weight[nonZero])
    group <- cell$person[cell$used][nonZero] + nPersons * cell$interval[cell$used][nonZero] + nCells * (match(weight[nonZero], levels) - 1)
    groupSums <- .influenceGroupSums(influence, used[nonZero], group)
    key <- groupSums$key
    groupSums <- groupSums$sums * levels[(key - 1) %/% nCells + 1]
    if (length(levels) == 1) {
      # (one group per cell)
      return(list(sums = groupSums, index = key))
    }
    cellOfGroup <- (key - 1) %% nCells + 1
    return(list(sums = rowsum(groupSums, cellOfGroup, reorder = FALSE), index = unique(cellOfGroup)))
  }))
}

# D of one block for the estimates with segment influences `influence` (a segments x estimates matrix)
# times the scale's `multiplier` of each segment, on the `segments` (rowId, strategyId, lastDay) with the
# block's `cells` (.stackedSegmentCells()), from the sums of the influence by cell of each term
# (`cellSums`, .stackedCellSums()): sum_k IJ_k dlog w_k / d alpha. Per term, person and interval l, R_il =
# sum of IJ_k over the person's segments whose weight runs through interval l, so that D = sum_terms
# sum_il kappa_il R_il z_il, assembled by kind of parameter as the information (.stackedBlock()).
.stackedBlockDerivative <- function(block, segments, influence, multiplier = rep(1, .influenceCount(influence)),
                                    cells = .stackedSegmentCells(list(blocks = list(block)), segments)[[1]],
                                    cellSums = .stackedCellSums(block, influence, multiplier, cells)) {
  nEstimates <- .influenceColumns(influence)
  derivative <- matrix(0, length(block$parameterX), nEstimates)
  nPersons <- length(block$personSeqId)
  time <- which(block$parameterX == 1)
  covariate <- which(block$parameterX != 1)
  functions <- unique(block$parameterTime[covariate])
  for (k in seq_along(block$terms)) {
    if (is.null(cellSums[[k]])) {
      next
    }
    term <- block$terms[[k]]
    nIntervals <- ncol(term$kappa)
    for (e in seq_len(nEstimates)) {
      sums <- matrix(0, nPersons, nIntervals)
      sums[cellSums[[k]]$index] <- cellSums[[k]]$sums[, e]
      # Reverse cumulative sums over the intervals:
      for (l in rev(seq_len(nIntervals - 1))) {
        sums[, l] <- sums[, l] + sums[, l + 1]
      }
      weighted <- term$kappa * sums
      derivative[time, e] <- derivative[time, e] + as.vector(colSums(weighted) %*% block$timeFunctions[, block$parameterTime[time], drop = FALSE])
      for (t in functions) {
        a <- covariate[block$parameterTime[covariate] == t]
        derivative[a, e] <- derivative[a, e] + as.vector(crossprod(block$x[, block$parameterX[a], drop = FALSE], weighted %*% block$timeFunctions[, t]))
      }
    }
  }
  return(derivative)
}

# The stacked covariance of estimates with segment (or row) influences `influence` (segments x estimates,
# or factored, see .influenceRows()) on `segments` (rowId, personSeqId, strategyId, stopTime), whose
# weights are stabilized when `stabilized` (the weight scale they use, see .stackedScaleMultipliers()):
# V = IJ'IJ + C'G + G'C + C'BC (see the top of the file).
.stackedVcov <- function(stacked, segments, influence, riskStart, stabilized) {
  if (!is.list(influence) || is.data.frame(influence)) {
    influence <- as.matrix(influence)
  }
  # Consecutive segments of a person's entry with the same keys (.stackedEntryKeys()) have the same
  # multipliers and cells, so the influence is summed over these classes once (each chunk of the
  # influence is formed once) and every group sum is taken from the class sums. (An entry, rowId, is a row
  # of the population and so has one personSeqId: a class is within one person.)
  lastDay <- segments$stopTime + riskStart - 1
  classes <- .stackedClasses(.stackedEntryKeys(stacked, segments$rowId, segments$strategyId, lastDay))
  first <- classes$first
  segments <- data.frame(rowId = segments$rowId[first], personSeqId = segments$personSeqId[first], strategyId = segments$strategyId[first], lastDay = lastDay[first])
  influence <- .influenceClassSums(influence, classes$class)
  # The cells of the classes, matched once for the multipliers and the derivatives:
  stacked$segmentCells <- .stackedSegmentCells(stacked, segments)
  multipliers <- .stackedScaleMultipliers(stacked, segments, stabilized)
  # The cell sums of a block are those of an earlier block with the same cells and multipliers (a model
  # and its numerator, when their multipliers are equal), so they are summed once:
  computed <- list()
  derivatives <- lapply(seq_along(stacked$blocks), function(b) {
    block <- stacked$blocks[[b]]
    multiplier <- if (block$numerator) multipliers$numerator else multipliers$denominator
    cells <- stacked$segmentCells[[b]]
    shape <- c(length(block$personSeqId), vapply(block$terms, function(term) ncol(term$kappa), integer(1)))
    same <- Position(function(entry) identical(entry$shape, shape) && identical(entry$cells, cells) && identical(entry$multiplier, multiplier), computed)
    cellSums <- if (is.na(same)) .stackedCellSums(block, influence, multiplier, cells) else computed[[same]]$cellSums
    computed[[length(computed) + 1]] <<- list(shape = shape, cells = cells, multiplier = multiplier, cellSums = cellSums)
    return(.stackedBlockDerivative(block, segments, influence, multiplier, cells, cellSums))
  })
  coefficients <- stacked$inverse %*% do.call(rbind, derivatives)
  personInfluence <- rowsum(influence, segments$personSeqId, reorder = FALSE)
  person <- match(unique(segments$personSeqId), stacked$persons)
  inPopulation <- !is.na(person)
  cross <- crossprod(stacked$psi[person[inPopulation], , drop = FALSE], personInfluence[inPopulation, , drop = FALSE])
  vcov <- crossprod(personInfluence) + crossprod(coefficients, cross) + crossprod(cross, coefficients) +
    crossprod(coefficients, stacked$B %*% coefficients)
  return(unname(vcov))
}

# Rows ---------------------------------------------------------------------------------------------

# The Kaplan-Meier or Aalen-Johansen rows of one horizon (as in .fitKm()) from the stacked covariance
# `vcov` of (risk, area under the risk) of the target and then of the comparator.
.stackedKmRows <- function(horizon, tau, risks, auc, vcov, competing) {
  se <- function(gradient) .standardError(sum(gradient * (vcov %*% gradient)), vcov, gradient)
  rows <- list(.deltaContrastRows(horizon, risks, rbind(c(1, 0, 0, 0), c(0, 0, 1, 0)), vcov, competing))
  if (!competing) {
    rows[[2]] <- dplyr::bind_rows(
      .amountRow(sprintf("rmst_target_%s", horizon), tau - auc[1], se(c(0, 1, 0, 0))),
      .amountRow(sprintf("rmst_comparator_%s", horizon), tau - auc[2], se(c(0, 0, 0, 1))),
      .differenceRow(sprintf("rmst_%s", horizon), auc[2] - auc[1], se(c(0, -1, 0, 1)))
    )
  }
  rows[[length(rows) + 1]] <- dplyr::bind_rows(
    .amountRow(sprintf("rmtl_target_%s", horizon), auc[1], se(c(0, 1, 0, 0))),
    .amountRow(sprintf("rmtl_comparator_%s", horizon), auc[2], se(c(0, 0, 0, 1))),
    .differenceRow(sprintf("rmtl_%s", horizon), auc[1] - auc[2], se(c(0, 1, 0, -1)))
  )
  return(dplyr::bind_rows(rows))
}

# The variance diagnostics of an outcome model: per estimand (`estimands`), the robust standard error (of
# `robustEstimates`, the same estimates with the robust variance), the stacked one, and their ratio (on
# the log scale for ratios); and per censoring model (`models`) the diagnostics of `stacked`
# (.stackedModelDiagnostics()).
.stackedVarianceDiagnostics <- function(estimates, robustEstimates, stacked) {
  se <- function(rows) ifelse(is.na(rows$seLogEstimate), rows$se, rows$seLogEstimate)
  seRobust <- se(robustEstimates)[match(estimates$estimand, robustEstimates$estimand)]
  seStacked <- se(estimates)
  return(list(
    estimands = tibble(
      estimand = estimates$estimand,
      seRobust = seRobust,
      seStacked = seStacked,
      ratio = seStacked / seRobust
    ),
    models = stacked$diagnostics
  ))
}

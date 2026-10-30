#' Fit a weighted clone-censor-weight outcome model
#'
#' Fits a weighted Cox or pooled logistic model with an appropriate variance
#' estimator to account for within-subject clone correlation.
#'
#' @param personTime Weighted person-time data frame from
#'   [createCensoringWeights()].
#' @param modelType `"cox"` or `"logistic"`.  Ignored when
#'   `fitOutcomeModelArgs` is supplied.
#' @param fitOutcomeModelArgs An optional `FitOutcomeModelArgs` object from
#'   [createFitOutcomeModelArgs()].  When supplied, `modelType` is taken from
#'   the args object.
#' @return A `CloneCensorWeightOutcomeModel` list with elements `summary`
#'   (data frame), `model` (fitted model object), `estimand`, `armLevels`,
#'   `clusterVariable`, and `varianceEstimator`.
#' @export
fitOutcomeModel <- function(personTime,
                            modelType         = "cox",
                            fitOutcomeModelArgs = NULL) {
  checkmate::assertDataFrame(personTime)
  checkmate::assertChoice(modelType, c("cox", "logistic"))
  checkmate::assertSubset(c("arm", "weight", "event", "tstop", "subjectId"), names(personTime))

  if (!is.null(fitOutcomeModelArgs)) {
    checkmate::assertR6(fitOutcomeModelArgs, "FitOutcomeModelArgs")
    modelType <- fitOutcomeModelArgs$modelType
  }

  varEst      <- if (!is.null(fitOutcomeModelArgs)) fitOutcomeModelArgs$varianceEstimator else "robust"
  nBootstraps <- if (!is.null(fitOutcomeModelArgs)) fitOutcomeModelArgs$nBootstraps       else 200L
  seed        <- if (!is.null(fitOutcomeModelArgs)) fitOutcomeModelArgs$seed              else 1L

  armLevels <- sort(unique(as.character(personTime$arm)))
  if (length(armLevels) != 2) stop("personTime must contain exactly two arm levels.", call. = FALSE)

  pd <- personTime
  pd$armValue <- as.integer(as.character(pd$arm) == armLevels[2])

  result <- switch(varEst,
    conventional = .fitConventional(pd, modelType, armLevels),
    robust       = .fitRobust(pd, modelType, armLevels),
    bootstrap    = .fitBootstrap(pd, modelType, armLevels, nBootstraps, seed),
    stop("Unknown varianceEstimator: ", varEst, call. = FALSE)
  )

  result$armLevels        <- armLevels
  result$clusterVariable  <- "subjectId"
  result$varianceEstimator <- varEst
  structure(result, class = "CloneCensorWeightOutcomeModel")
}

# ---------------------------------------------------------------------------
# Conventional path (Cyclops with profile-likelihood CIs)
# ---------------------------------------------------------------------------

.fitConventional <- function(pd, modelType, armLevels) {
  estimand <- if (modelType == "cox") "hazard ratio" else "odds ratio"
  outcomes <- data.frame(
    rowId = seq_len(nrow(pd)),
    y     = as.integer(pd$event),
    time  = as.numeric(pd$tstop),
    stringsAsFactors = FALSE
  )
  covariates <- data.frame(
    rowId = seq_len(nrow(pd)), covariateId = 1L,
    covariateValue = pd$armValue, stringsAsFactors = FALSE
  )
  if (modelType == "logistic") {
    covariates <- rbind(covariates, data.frame(
      rowId = seq_len(nrow(pd)), covariateId = 2L,
      covariateValue = as.numeric(pd$tstop), stringsAsFactors = FALSE
    ))
  }
  cyclopsData <- Cyclops::convertToCyclopsData(
    outcomes   = outcomes,
    covariates = covariates,
    modelType  = if (modelType == "cox") "cox" else "lr",
    addIntercept = modelType == "logistic",
    checkRowIds = FALSE, normalize = NULL, quiet = TRUE
  )
  fit <- Cyclops::fitCyclopsModel(
    cyclopsData,
    prior   = Cyclops::createPrior("none"),
    control = Cyclops::createControl(),
    weights = pd[["weight"]]
  )
  if (!identical(fit[["return_flag"]], "SUCCESS")) {
    stop("Cyclops model did not converge: ", fit[["return_flag"]], call. = FALSE)
  }
  coef <- unname(stats::coef(fit)["1"])
  ci   <- stats::confint(fit, parm = 1, includePenalty = TRUE)
  se   <- (ci[1, 3] - ci[1, 2]) / (2 * stats::qnorm(0.975))
  summary <- data.frame(
    modelType = modelType, term = "arm",
    estimate = exp(coef), lower95 = exp(ci[1, 2]), upper95 = exp(ci[1, 3]),
    pValue = 2 * stats::pnorm(-abs(coef / se)),
    logLikelihood = as.numeric(fit[["log_likelihood"]]),
    stringsAsFactors = FALSE
  )
  list(summary = summary, model = fit, estimand = estimand)
}

# ---------------------------------------------------------------------------
# Robust path (sandwich / Lin-Wei-Ying variance)
# ---------------------------------------------------------------------------

.fitRobust <- function(pd, modelType, armLevels) {
  estimand <- if (modelType == "cox") "hazard ratio" else "odds ratio"

  if (modelType == "cox") {
    fit  <- survival::coxph(
      survival::Surv(tstop, event) ~ armValue + survival::cluster(subjectId),
      data    = pd,
      weights = pd$weight,
      ties    = "efron"
    )
    coef  <- unname(stats::coef(fit)["armValue"])
    # vcov() returns the robust (Lin-Wei-Ying) variance when cluster() is present
    rse   <- sqrt(diag(stats::vcov(fit)))["armValue"]
    lower <- coef - stats::qnorm(0.975) * rse
    upper <- coef + stats::qnorm(0.975) * rse
    ll    <- as.numeric(stats::logLik(fit))
  } else {
    fit   <- stats::glm(
      event ~ armValue + tstop,
      data    = pd,
      family  = stats::binomial(),
      weights = pd$weight
    )
    rse   <- .clusteredSE(fit, pd$subjectId)["armValue"]
    coef  <- unname(stats::coef(fit)["armValue"])
    lower <- coef - stats::qnorm(0.975) * rse
    upper <- coef + stats::qnorm(0.975) * rse
    ll    <- as.numeric(stats::logLik(fit))
  }

  pval <- 2 * stats::pnorm(-abs(coef / rse))
  summary <- data.frame(
    modelType = modelType, term = "arm",
    estimate = exp(coef), lower95 = exp(lower), upper95 = exp(upper),
    pValue = pval, logLikelihood = ll, stringsAsFactors = FALSE
  )
  list(summary = summary, model = fit, estimand = estimand)
}

# Lin-Wei-Ying sandwich SE for a GLM, clustering on `cluster`.
# Returns a named vector of SEs (one per coefficient).
.clusteredSE <- function(fit, cluster) {
  X        <- stats::model.matrix(fit)
  rawResid <- fit$y - stats::fitted(fit)
  scores   <- rawResid * X
  # Sum score contributions within each cluster
  clusterMat <- rowsum(scores, cluster)
  bread      <- tryCatch(
    solve(crossprod(X, X * fit$weights)),
    error = function(e) solve(crossprod(X, X * fit$weights) + diag(1e-8, ncol(X)))
  )
  meat   <- crossprod(clusterMat)
  vcov_cl <- bread %*% meat %*% bread
  sqrt(diag(vcov_cl))
}

# ---------------------------------------------------------------------------
# Bootstrap path (cluster bootstrap over subjects)
# ---------------------------------------------------------------------------

.fitBootstrap <- function(pd, modelType, armLevels, nBootstraps, seed) {
  # Point estimate via fast base-R fit (same as robust path, minus clustering)
  base <- .fitRobust(pd, modelType, armLevels)
  coefPoint <- log(base$summary$estimate)

  subjectIds <- unique(pd$subjectId)
  n          <- length(subjectIds)
  set.seed(seed)

  bootCoefs <- vapply(seq_len(nBootstraps), function(b) {
    sampled  <- sample(subjectIds, n, replace = TRUE)
    bootData <- .buildBootstrapData(pd, sampled)
    tryCatch({
      if (modelType == "cox") {
        fit <- survival::coxph(
          survival::Surv(tstop, event) ~ armValue,
          data    = bootData,
          weights = bootData$weight,
          ties    = "efron"
        )
        unname(stats::coef(fit)["armValue"])
      } else {
        fit <- stats::glm(
          event ~ armValue + tstop,
          data    = bootData,
          family  = stats::binomial(),
          weights = bootData$weight
        )
        unname(stats::coef(fit)["armValue"])
      }
    }, error = function(e) NA_real_)
  }, numeric(1))

  bootCoefs <- bootCoefs[!is.na(bootCoefs)]
  if (length(bootCoefs) < 10L) {
    warning("Fewer than 10 successful bootstrap replicates; CIs may be unreliable.",
            call. = FALSE)
  }
  bootSE <- stats::sd(bootCoefs)
  lower  <- coefPoint - stats::qnorm(0.975) * bootSE
  upper  <- coefPoint + stats::qnorm(0.975) * bootSE
  pval   <- 2 * stats::pnorm(-abs(coefPoint / bootSE))

  estimand <- if (modelType == "cox") "hazard ratio" else "odds ratio"
  summary <- data.frame(
    modelType = modelType, term = "arm",
    estimate = exp(coefPoint), lower95 = exp(lower), upper95 = exp(upper),
    pValue = pval, logLikelihood = base$summary$logLikelihood,
    stringsAsFactors = FALSE
  )
  list(summary = summary, model = base$model, estimand = estimand)
}

.buildBootstrapData <- function(pd, sampledIds) {
  do.call(rbind, lapply(seq_along(sampledIds), function(i) {
    rows <- pd[pd$subjectId == sampledIds[i], , drop = FALSE]
    rows$subjectId <- i  # renumber to make clusters distinct
    rows
  }))
}

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

#' Compute covariate means by clone arm
#'
#' Computes unweighted covariate means separately for each arm, using the
#' FeatureExtraction sparse covariate format (`rowId`, `covariateId`,
#' `covariateValue`) stored in the data object.
#'
#' @param cloneCensorWeightData A `CloneCensorWeightData` object (after
#'   cloning, so the `cohorts` table has an `arm` column).
#' @param weights Optional weighted person-time data frame. When supplied, arm
#'   membership is taken from its `(subjectId, arm)` pairs rather than from
#'   the cohorts table.
#' @return A data frame with columns `covariateId`, `covariate`, `arm`, and
#'   `mean`.  Returns an empty data frame when no covariates are present.
#' @export
computeCovariateBalance <- function(cloneCensorWeightData, weights = NULL) {
  checkmate::assertClass(cloneCensorWeightData, "CloneCensorWeightData")
  covariates   <- .readCloneCensorWeightTable(cloneCensorWeightData, "covariates")
  covariateRef <- .readCloneCensorWeightTable(cloneCensorWeightData, "covariateRef")
  cohorts      <- .readCloneCensorWeightTable(cloneCensorWeightData, "cohorts")

  if (!is.null(weights)) cohorts <- unique(weights[c("subjectId", "arm")])
  if (!"arm" %in% names(cohorts)) cohorts$arm <- "eligible"

  empty <- data.frame(covariateId = integer(), covariate = character(),
                      arm = character(), mean = numeric(), stringsAsFactors = FALSE)

  isSparse <- all(c("rowId", "covariateId", "covariateValue") %in% names(covariates))
  if (!isSparse || nrow(covariates) == 0) return(empty)

  joined <- merge(
    covariates[, c("rowId", "covariateId", "covariateValue")],
    cohorts[, c("subjectId", "arm"), drop = FALSE],
    by.x = "rowId", by.y = "subjectId",
    all.x = FALSE
  )
  if (nrow(joined) == 0) return(empty)

  covariateIds <- sort(unique(joined$covariateId))
  do.call(rbind, lapply(covariateIds, function(cid) {
    sub   <- joined[joined$covariateId == cid, , drop = FALSE]
    means <- tapply(sub$covariateValue, sub$arm, mean, na.rm = TRUE)
    cname <- if (nrow(covariateRef) > 0 && cid %in% covariateRef$covariateId)
      covariateRef$covariateName[covariateRef$covariateId == cid][1]
    else paste0("cov_", cid)
    data.frame(covariateId = cid, covariate = cname,
               arm = names(means), mean = as.numeric(means),
               stringsAsFactors = FALSE)
  }))
}

#' Plot covariate means by arm
#' @param balance Output from [computeCovariateBalance()].
#' @return A `ggplot2` plot.
#' @export
plotCovariateBalance <- function(balance) {
  checkmate::assertDataFrame(balance)
  ggplot2::ggplot(balance, ggplot2::aes(x = rlang::.data$mean,
                                        y = rlang::.data$covariate,
                                        colour = rlang::.data$arm)) +
    ggplot2::geom_point() + ggplot2::theme_minimal()
}

#' Plot weighted Kaplan-Meier curves
#' @param personTime Weighted person-time data frame.
#' @return A `ggplot2` plot.
#' @export
plotKaplanMeier <- function(personTime) {
  checkmate::assertDataFrame(personTime)
  weightVector <- personTime[["weight"]]
  fit <- survival::survfit(survival::Surv(tstop, event) ~ arm, data = personTime,
                           weights = weightVector)
  plotData <- data.frame(time = fit$time, survival = fit$surv,
                         arm = rep(names(fit$strata), fit$strata))
  ggplot2::ggplot(plotData, ggplot2::aes(x = rlang::.data$time,
                                         y = rlang::.data$survival,
                                         colour = rlang::.data$arm)) +
    ggplot2::geom_step() + ggplot2::theme_minimal() +
    ggplot2::labs(x = "Days since trial start", y = "Survival probability")
}

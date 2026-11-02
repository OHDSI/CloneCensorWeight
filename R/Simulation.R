#' Simulate a clone-censor-weight dataset with confounding
#'
#' Generates a synthetic `CloneCensorWeightData` object using a
#' piecewise-exponential data-generating process with user-controlled
#' confounding, treatment prevalence, and true effect size.  Intended for
#' Monte Carlo power studies and validation of the CCW pipeline.
#'
#' **Data-generating process**
#' \enumerate{
#'   \item Generate `nCovariates` standard-normal confounders per subject.
#'   \item Compute a linear predictor `LP = confoundingStrength × rowMean(Z)`.
#'   \item Assign treatment via `Bernoulli(plogis(LP + logit(treatmentPrevalence)))`.
#'   \item Draw initiation timing `Uniform(1, gracePeriod)` days for treated subjects.
#'   \item Compute individual baseline hazard `h0 = baselineHazard × exp(confoundingStrength × LP)`.
#'   \item Draw event time from a piecewise exponential: `h0` before initiation,
#'     `h0 × trueHazardRatio` after.
#'   \item Draw loss-to-follow-up censoring time from `Exponential(censorRate / intervalLength)`.
#'   \item Set observed follow-up to `min(event, censor, maxFollowUpDays)`.
#' }
#'
#' `baselineHazard` and `censorRate` are expressed **per interval** (e.g. per
#' 30-day block).  Internally they are converted to per-day rates.
#'
#' @param n Number of subjects.
#' @param nCovariates Number of continuous baseline confounders.
#' @param treatmentPrevalence Marginal probability of treatment initiation.
#' @param trueHazardRatio True causal hazard ratio for the treatment effect.
#' @param baselineHazard Baseline hazard rate per `intervalLength`-day interval
#'   for an untreated subject with LP = 0.
#' @param censorRate Loss-to-follow-up rate per `intervalLength`-day interval.
#' @param confoundingStrength Coefficient linking the confounder LP to both
#'   treatment propensity and the outcome hazard.
#' @param gracePeriod Grace period in days used later in [createClonesAndCensor()].
#' @param maxFollowUpDays Maximum observation window in days.
#' @param intervalLength Discretisation interval in days (should match the
#'   `intervalLength` passed to [createCreateClonesAndCensorArgs()]).
#' @param seed Integer RNG seed.  Set to `NULL` for no seeding.
#'
#' @return A `CloneCensorWeightData` object ready for the CCW pipeline.
#' @export
simulateCcwData <- function(n                    = 1000,
                            nCovariates          = 5,
                            treatmentPrevalence  = 0.5,
                            trueHazardRatio      = 1.0,
                            baselineHazard       = 0.02,
                            censorRate           = 0.01,
                            confoundingStrength  = 0.5,
                            gracePeriod          = 30,
                            maxFollowUpDays      = 365,
                            intervalLength       = 30,
                            seed                 = NULL) {

  checkmate::assertCount(n, positive = TRUE)
  checkmate::assertCount(nCovariates, positive = TRUE)
  checkmate::assertNumber(treatmentPrevalence, lower = 0, upper = 1)
  checkmate::assertNumber(trueHazardRatio, lower = 0)
  checkmate::assertNumber(baselineHazard, lower = 0)
  checkmate::assertNumber(censorRate, lower = 0)
  checkmate::assertNumber(confoundingStrength)
  checkmate::assertCount(gracePeriod, positive = TRUE)
  checkmate::assertCount(maxFollowUpDays, positive = TRUE)
  checkmate::assertCount(intervalLength, positive = TRUE)
  if (!is.null(seed)) checkmate::assertIntegerish(seed, len = 1)

  if (!is.null(seed)) set.seed(seed)

  # --- Confounders -----------------------------------------------------------
  Z  <- matrix(stats::rnorm(n * nCovariates), nrow = n, ncol = nCovariates)
  LP <- confoundingStrength * rowMeans(Z)

  # --- Treatment assignment --------------------------------------------------
  pTreat  <- stats::plogis(LP + stats::qlogis(treatmentPrevalence))
  treated <- stats::rbinom(n, 1L, pTreat) == 1L

  initDay <- ifelse(
    treated,
    ceiling(stats::runif(n, 0, gracePeriod)),
    NA_real_
  )

  # --- Time-to-event (piecewise exponential) ---------------------------------
  hPerDay <- baselineHazard / intervalLength
  cPerDay <- censorRate     / intervalLength
  h0      <- hPerDay * exp(confoundingStrength * LP)  # individual baseline hazard

  eventDay  <- numeric(n)
  for (i in seq_len(n)) {
    if (!is.na(initDay[i])) {
      # Pre-treatment phase: exponential until initiation
      tPre <- stats::rexp(1, rate = h0[i])
      if (tPre < initDay[i]) {
        eventDay[i] <- tPre        # event before treatment
      } else {
        # Post-treatment phase: start from initiation day
        tPost       <- stats::rexp(1, rate = h0[i] * trueHazardRatio)
        eventDay[i] <- initDay[i] + tPost
      }
    } else {
      eventDay[i] <- stats::rexp(1, rate = h0[i])
    }
  }

  censorDay <- stats::rexp(n, rate = cPerDay)
  obsDay    <- pmin(eventDay, censorDay, maxFollowUpDays)
  event     <- eventDay <= censorDay & eventDay <= maxFollowUpDays

  # --- Cohort table ----------------------------------------------------------
  indexDate <- as.Date("2020-01-01")
  cohorts   <- data.frame(
    subjectId       = seq_len(n),
    cohortStartDate = indexDate,
    cohortEndDate   = indexDate + ceiling(obsDay),
    initiationDate  = as.Date(ifelse(
      !is.na(initDay),
      as.character(indexDate + ceiling(initDay)),
      NA_character_
    )),
    stringsAsFactors = FALSE
  )

  # --- Outcomes table --------------------------------------------------------
  eventIdx <- which(event)
  outcomes <- data.frame(
    subjectId   = eventIdx,
    outcomeDate = indexDate + ceiling(eventDay[eventIdx]),
    outcomeId   = rep(1L, length(eventIdx)),
    stringsAsFactors = FALSE
  )

  # --- Covariates in FeatureExtraction sparse format -------------------------
  # Each continuous confounder gets a unique covariateId = 1000 + column index
  rows <- lapply(seq_len(nCovariates), function(j) {
    data.frame(
      rowId          = seq_len(n),
      covariateId    = 1000L + j,
      covariateValue = Z[, j],
      stringsAsFactors = FALSE
    )
  })
  covariates <- do.call(rbind, rows)

  covariateRef <- data.frame(
    covariateId   = 1000L + seq_len(nCovariates),
    covariateName = paste0("confounder_", seq_len(nCovariates)),
    analysisId    = rep(1L, nCovariates),
    conceptId     = rep(0L, nCovariates),
    stringsAsFactors = FALSE
  )

  analysisRef <- data.frame(
    analysisId       = 1L,
    analysisName     = "Simulated confounders",
    domainId         = "Simulation",
    startDay         = -365L,
    endDay           = 0L,
    isBinary         = "N",
    missingMeansZero = "N",
    stringsAsFactors = FALSE
  )

  createCloneCensorWeightData(cohorts, outcomes, covariates, covariateRef, analysisRef)
}

#' Run a Monte Carlo simulation study
#'
#' Repeatedly calls [simulateCcwData()] and the full CCW pipeline to estimate
#' operating characteristics (bias, RMSE, CI coverage, power/type-I-error rate)
#' across `nSims` replications.
#'
#' @param nSims Number of Monte Carlo replications.
#' @param n Subjects per replication.
#' @param trueHazardRatio True causal hazard ratio.  When 1.0 the returned
#'   `power` column is the empirical type-I error rate.
#' @param baselineHazard Baseline hazard per interval.  See [simulateCcwData()].
#' @param censorRate Loss-to-follow-up rate per interval.  See
#'   [simulateCcwData()].
#' @param nCovariates Number of confounders to simulate.
#' @param treatmentPrevalence Marginal treatment probability.
#' @param confoundingStrength Strength of confounding.  See [simulateCcwData()].
#' @param gracePeriod Grace period in days.
#' @param maxFollowUpDays Maximum observation window.
#' @param intervalLength Discretisation interval in days.
#' @param createCensoringWeightsArgs An optional
#'   [createCreateCensoringWeightsArgs()] object.
#' @param fitOutcomeModelArgs An optional [createFitOutcomeModelArgs()] object.
#'   Defaults to `"cox"` with `"robust"` variance.
#' @param seed Base integer seed; replication `k` uses `seed + k`.
#'
#' @return A data frame with one row per replication and columns:
#'   `simId`, `converged`, `logRr`, `seLogRr`, `estimate`, `lower95`,
#'   `upper95`, `pValue`, `covered`, `bias`, `bias2`.
#' @export
runSimulation <- function(nSims                    = 100,
                          n                        = 500,
                          trueHazardRatio          = 1.0,
                          baselineHazard           = 0.02,
                          censorRate               = 0.01,
                          nCovariates              = 5,
                          treatmentPrevalence      = 0.5,
                          confoundingStrength      = 0.5,
                          gracePeriod              = 30,
                          maxFollowUpDays          = 365,
                          intervalLength           = 30,
                          createCensoringWeightsArgs = NULL,
                          fitOutcomeModelArgs        = NULL,
                          seed                     = 1L) {

  checkmate::assertCount(nSims, positive = TRUE)
  checkmate::assertCount(n, positive = TRUE)
  checkmate::assertNumber(trueHazardRatio, lower = 0)
  checkmate::assertIntegerish(seed, len = 1)

  if (is.null(fitOutcomeModelArgs)) {
    fitOutcomeModelArgs <- createFitOutcomeModelArgs(
      modelType          = "cox",
      varianceEstimator  = "robust"
    )
  }

  cloneArgs <- createCreateClonesAndCensorArgs(
    arms = list(
      createArmDefinition("initiate", "initiateWithin", gracePeriod = gracePeriod),
      createArmDefinition("never",    "deferBeyond",    gracePeriod = gracePeriod)
    ),
    intervalLength = intervalLength
  )

  popArgs <- createCreateStudyPopulationArgs(
    riskWindowEnd             = maxFollowUpDays,
    endAnchor                 = "cohort start",
    removeSubjectsWithPriorOutcome = FALSE
  )

  results <- lapply(seq_len(nSims), function(k) {
    simSeed <- seed + k
    .runOneSim(
      simSeed                    = simSeed,
      n                          = n,
      trueHazardRatio            = trueHazardRatio,
      baselineHazard             = baselineHazard,
      censorRate                 = censorRate,
      nCovariates                = nCovariates,
      treatmentPrevalence        = treatmentPrevalence,
      confoundingStrength        = confoundingStrength,
      gracePeriod                = gracePeriod,
      maxFollowUpDays            = maxFollowUpDays,
      intervalLength             = intervalLength,
      cloneArgs                  = cloneArgs,
      popArgs                    = popArgs,
      createCensoringWeightsArgs = createCensoringWeightsArgs,
      fitOutcomeModelArgs        = fitOutcomeModelArgs
    )
  })

  do.call(rbind, results)
}

#' Summarize Monte Carlo simulation results
#'
#' Computes operating characteristics from the output of [runSimulation()].
#'
#' @param simResults A data frame returned by [runSimulation()].
#' @param trueHazardRatio The true hazard ratio used in the simulation.  When
#'   `NULL`, it is inferred as `exp(mean(simResults$logRr))` — useful when the
#'   simulation is run under the null.
#' @param alpha Significance level for power/type-I-error computation.
#'
#' @return A named list with elements: `nSims`, `convergenceRate`, `trueLogRr`,
#'   `meanLogRr`, `bias`, `rmse`, `coverage`, and `power`.
#' @export
summarizeSimulation <- function(simResults, trueHazardRatio = NULL, alpha = 0.05) {
  checkmate::assertDataFrame(simResults)
  checkmate::assertSubset(c("converged", "logRr", "lower95", "upper95", "pValue"),
                          names(simResults))
  checkmate::assertNumber(alpha, lower = 0, upper = 1)

  conv <- simResults[simResults$converged, , drop = FALSE]

  trueLogRr <- if (!is.null(trueHazardRatio)) {
    log(trueHazardRatio)
  } else {
    log(trueHazardRatio)  # will be NULL — use NULL-safe below
  }

  # Infer from caller argument
  if (is.null(trueHazardRatio)) {
    trueLogRr <- NA_real_
  } else {
    trueLogRr <- log(trueHazardRatio)
  }

  bias     <- if (!is.na(trueLogRr)) mean(conv$logRr - trueLogRr, na.rm = TRUE) else NA_real_
  rmse     <- if (!is.na(trueLogRr)) sqrt(mean((conv$logRr - trueLogRr)^2, na.rm = TRUE)) else NA_real_
  coverage <- if (!is.na(trueLogRr)) {
    trueHr <- exp(trueLogRr)
    mean(conv$lower95 <= trueHr & conv$upper95 >= trueHr, na.rm = TRUE)
  } else NA_real_

  power <- mean(conv$pValue < alpha, na.rm = TRUE)

  list(
    nSims           = nrow(simResults),
    convergenceRate = mean(simResults$converged),
    trueLogRr       = trueLogRr,
    meanLogRr       = mean(conv$logRr, na.rm = TRUE),
    bias            = bias,
    rmse            = rmse,
    coverage        = coverage,
    power           = power
  )
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.runOneSim <- function(simSeed, n, trueHazardRatio, baselineHazard, censorRate,
                       nCovariates, treatmentPrevalence, confoundingStrength,
                       gracePeriod, maxFollowUpDays, intervalLength,
                       cloneArgs, popArgs, createCensoringWeightsArgs,
                       fitOutcomeModelArgs) {

  failed <- data.frame(
    simId     = simSeed,
    converged = FALSE,
    logRr = NA_real_, seLogRr = NA_real_,
    estimate = NA_real_, lower95 = NA_real_, upper95 = NA_real_,
    pValue = NA_real_, covered = NA, bias = NA_real_, bias2 = NA_real_,
    stringsAsFactors = FALSE
  )

  tryCatch({
    data <- simulateCcwData(
      n                   = n,
      nCovariates         = nCovariates,
      treatmentPrevalence = treatmentPrevalence,
      trueHazardRatio     = trueHazardRatio,
      baselineHazard      = baselineHazard,
      censorRate          = censorRate,
      confoundingStrength = confoundingStrength,
      gracePeriod         = gracePeriod,
      maxFollowUpDays     = maxFollowUpDays,
      intervalLength      = intervalLength,
      seed                = simSeed
    )

    createStudyPopulation(data, outcomeId = 1L, popArgs)
    cloned   <- createClonesAndCensor(data, cloneArgs)
    pt       <- discretizePersonTime(cloned, intervalLength = intervalLength)
    weighted <- createCensoringWeights(pt, cloned, createCensoringWeightsArgs)
    model    <- fitOutcomeModel(weighted, fitOutcomeModelArgs = fitOutcomeModelArgs)

    s      <- model$summary
    logRr  <- log(s$estimate)
    seLogRr <- (log(s$upper95) - log(s$lower95)) / (2 * stats::qnorm(0.975))
    trueLogRr <- log(trueHazardRatio)

    data.frame(
      simId     = simSeed,
      converged = TRUE,
      logRr     = logRr,
      seLogRr   = seLogRr,
      estimate  = s$estimate,
      lower95   = s$lower95,
      upper95   = s$upper95,
      pValue    = s$pValue,
      covered   = s$lower95 <= trueHazardRatio & s$upper95 >= trueHazardRatio,
      bias      = logRr - trueLogRr,
      bias2     = (logRr - trueLogRr)^2,
      stringsAsFactors = FALSE
    )
  },
  error = function(e) failed)
}

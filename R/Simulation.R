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

# The DGP is a daily discrete-time process on days 0..maxDays. Each person has three latent
# processes that are independent given X and the start day A:
# - initiation of the treatment, with a constant daily hazard (A is geometric on 0, 1, ...), or, with
#   day-0 initiation parameters, a day-0 hazard and a constant daily hazard after day 0;
# - the outcome, whose daily hazard depends on X and on A (dgp "a" or "b");
# - death, whose daily hazard depends on X only.
# A latent event day is drawn by inversion: with E ~ Exp(1), the event day is the first day d on
# which the cumulative hazard sum_{k <= d} -log(1 - h(k)) reaches E. Then P(T > d) =
# prod_{k <= d} (1 - h(k)), so this is exact, and the same E can be reused under every regime
# (common random numbers). On one day the outcome comes before death ("eventsFirst"): the outcome
# counts when its day is on or before the death day.
#
# An outcome hazard path is described per person by (switchDay, startDay): on days k <= switchDay
# the hazard is the untreated one; on days k > switchDay it is the treated one, with
# s = k - startDay days since initiation (s >= 1).

.simulationSettingsFields <- c(
  "dgp", "nCovariates", "graceEnd", "maxDays", "binaryPrevalence", "initiationHazard", "outcomeHazard",
  "deathHazard", "observationEndHazard", "initiationCoefficients", "outcomeCoefficients",
  "deathCoefficients", "treatmentEffect", "effectRampDays"
)

.simulationTargetId <- 1
.simulationTreatmentId <- 2
.simulationOutcomeId <- 3

# -log(1 - p), the daily increment of the cumulative hazard for a daily probability p. abs() makes it
# +0 when p = 0 (-log1p(-0) is -0, and e / -0 would be -Inf).
.hazardIncrement <- function(probability) {
  return(abs(log1p(-probability)))
}

.dailyCumulativeHazard <- function(eta) {
  return(.hazardIncrement(stats::plogis(eta)))
}

# Covariate j is binary (Bernoulli(binaryPrevalence)) when j is odd and standard normal when j is
# even.
.covariateIsBinary <- function(nCovariates) {
  return(seq_len(nCovariates) %% 2 == 1)
}

.simulateCovariates <- function(n, settings) {
  isBinary <- .covariateIsBinary(settings$nCovariates)
  x <- matrix(0, nrow = n, ncol = settings$nCovariates)
  for (j in seq_len(settings$nCovariates)) {
    if (isBinary[j]) {
      x[, j] <- stats::rbinom(n, 1, settings$binaryPrevalence)
    } else {
      x[, j] <- stats::rnorm(n)
    }
  }
  return(x)
}

# Linear predictors (logit of the daily hazard) for the untreated person. `initiationDay0` is that of
# day 0, which equals `initiation` unless the settings have day-0 initiation parameters.
.linearPredictors <- function(x, settings) {
  linearPredictor <- function(hazard, coefficients) {
    return(stats::qlogis(hazard) + as.vector(x %*% coefficients))
  }
  initiation <- linearPredictor(settings$initiationHazard, settings$initiationCoefficients)
  initiationDay0 <- initiation
  if (.hasDay0Initiation(settings)) {
    initiationDay0 <- linearPredictor(
      if (is.null(settings$initiationHazardDay0)) settings$initiationHazard else settings$initiationHazardDay0,
      if (is.null(settings$initiationCoefficientsDay0)) settings$initiationCoefficients else settings$initiationCoefficientsDay0
    )
  }
  return(list(
    initiation = initiation,
    initiationDay0 = initiationDay0,
    outcome = linearPredictor(settings$outcomeHazard, settings$outcomeCoefficients),
    death = linearPredictor(settings$deathHazard, settings$deathCoefficients)
  ))
}

# Does the initiation hazard of day 0 have its own parameters?
.hasDay0Initiation <- function(settings) {
  return(!is.null(settings$initiationHazardDay0) || !is.null(settings$initiationCoefficientsDay0))
}

# The natural-course start day: the first day on which the cumulative initiation hazard reaches `e`,
# with the day-0 hazard on day 0 and the constant daily hazard after it (Inf if after maxDays).
.initiationEventDay <- function(e, eta, settings, maxDays) {
  rate <- .dailyCumulativeHazard(eta$initiation)
  if (!.hasDay0Initiation(settings)) {
    return(.constantHazardEventDay(e, rate, maxDays))
  }
  day0Rate <- .dailyCumulativeHazard(eta$initiationDay0)
  day <- ceiling((e - day0Rate) / rate)
  day[e <= day0Rate] <- 0
  day[day > maxDays] <- Inf
  return(day)
}

# First day d >= 0 on which a constant daily cumulative hazard `rate` reaches `e`; Inf if after
# maxDays.
.constantHazardEventDay <- function(e, rate, maxDays) {
  day <- ceiling(e / rate) - 1
  day[day > maxDays] <- Inf
  return(day)
}

# The outcome hazard path of a person who starts the treatment on day `startDay` (Inf if never).
.regularPath <- function(startDay, settings) {
  if (settings$dgp == "a") {
    switchDay <- ifelse(startDay <= settings$graceEnd, settings$graceEnd, Inf)
  } else {
    switchDay <- startDay
  }
  return(list(switchDay = switchDay, startDay = startDay))
}

# Days since initiation over which the effect of dgp "b" builds up (1 for dgp "a", whose effect is
# constant).
.rampDays <- function(settings) {
  if (settings$dgp == "b") {
    return(settings$effectRampDays)
  }
  return(1)
}

# The effect on the logit of the daily outcome hazard, s days after initiation (s >= 1).
.treatmentEffect <- function(s, settings) {
  rampDays <- .rampDays(settings)
  return(settings$treatmentEffect * pmin(s, rampDays) / rampDays)
}

# First day on which the cumulative outcome hazard along the path reaches `e` (Inf if after maxDays).
.outcomeEventDay <- function(e, etaOutcome, path, settings) {
  maxDays <- settings$maxDays
  rampDays <- .rampDays(settings)
  untreatedRate <- .dailyCumulativeHazard(etaOutcome)
  eventDay <- .constantHazardEventDay(e, untreatedRate, maxDays)
  eventDay[eventDay > path$switchDay] <- Inf

  # Persons who have not had the event by the switch day continue on the treated hazard:
  active <- which(is.infinite(eventDay) & path$switchDay < maxDays)
  cumulativeHazard <- (path$switchDay[active] + 1) * untreatedRate[active]
  day <- path$switchDay[active] + 1
  # Days while the effect is still building up (dgp "b"), one day at a time:
  for (step in seq_len(rampDays - 1)) {
    s <- day - path$startDay[active]
    inRamp <- which(s < rampDays & day <= maxDays)
    if (length(inRamp) == 0) {
      break
    }
    cumulativeHazard[inRamp] <- cumulativeHazard[inRamp] +
      .dailyCumulativeHazard(etaOutcome[active[inRamp]] + .treatmentEffect(s[inRamp], settings))
    crossed <- inRamp[e[active[inRamp]] <= cumulativeHazard[inRamp]]
    eventDay[active[crossed]] <- day[crossed]
    day[inRamp] <- day[inRamp] + 1
    if (length(crossed) > 0) {
      active <- active[-crossed]
      cumulativeHazard <- cumulativeHazard[-crossed]
      day <- day[-crossed]
    }
  }
  # From `day` on the hazard is constant at the full effect:
  treatedRate <- .dailyCumulativeHazard(etaOutcome[active] + settings$treatmentEffect)
  tailDay <- day - 1 + ceiling((e[active] - cumulativeHazard) / treatedRate)
  tailDay[tailDay > maxDays] <- Inf
  eventDay[active] <- tailDay
  return(eventDay)
}

# Draws, for persons who have not started by day graceEnd, a start day a in 0..graceEnd from
# P(A = a | A <= graceEnd, no outcome or death on or before graceEnd, X), which is proportional to
# P(A = a | X) * P(no outcome on or before graceEnd | A = a, X). Death does not depend on the
# treatment, so its survival cancels. Two passes (total, then inverse CDF), to avoid an
# n x (graceEnd + 1) matrix. With a day-0 initiation hazard h0 (etaInitiationDay0), P(A = 0 | X) = h0
# and P(A = a | X) = (1 - h0) (1 - h)^(a - 1) h for a >= 1.
.drawGraceStartDay <- function(u, etaInitiation, etaOutcome, graceEnd, settings, etaInitiationDay0 = NULL) {
  initiationHazard <- stats::plogis(etaInitiation)
  if (!is.null(etaInitiationDay0)) {
    day0Hazard <- stats::plogis(etaInitiationDay0)
  }
  untreatedRate <- .dailyCumulativeHazard(etaOutcome)
  if (settings$dgp == "a") {
    treatedRate <- .dailyCumulativeHazard(etaOutcome + settings$treatmentEffect)
  }
  startDay <- rep(NA_real_, length(u))
  total <- 0
  for (pass in 1:2) {
    rampSum <- 0
    cumulativeWeight <- 0
    for (daysTreated in 0:graceEnd) {
      a <- graceEnd - daysTreated
      # Cumulative outcome hazard through day graceEnd for a start on day a:
      if (settings$dgp == "a") {
        if (a <= settings$graceEnd && settings$graceEnd < graceEnd) {
          cumulativeHazard <- (settings$graceEnd + 1) * untreatedRate + (graceEnd - settings$graceEnd) * treatedRate
        } else {
          cumulativeHazard <- (graceEnd + 1) * untreatedRate
        }
      } else {
        if (daysTreated > 0) {
          rampSum <- rampSum + .dailyCumulativeHazard(etaOutcome + .treatmentEffect(daysTreated, settings))
        }
        cumulativeHazard <- (a + 1) * untreatedRate + rampSum
      }
      if (is.null(etaInitiationDay0)) {
        pStart <- (1 - initiationHazard)^a * initiationHazard
      } else if (a == 0) {
        pStart <- day0Hazard
      } else {
        pStart <- (1 - day0Hazard) * (1 - initiationHazard)^(a - 1) * initiationHazard
      }
      weight <- pStart * exp(-cumulativeHazard)
      if (pass == 1) {
        total <- total + weight
      } else {
        cumulativeWeight <- cumulativeWeight + weight
        pick <- is.na(startDay) & cumulativeWeight >= u * total
        startDay[pick] <- a
      }
    }
  }
  return(startDay)
}

# The outcome hazard path under a strategy, given the natural-course start day.
.strategyPath <- function(strategy, graceRegime, naturalStartDay, laterStartE, u, eta, settings) {
  if (strategy$type == "neverInitiate") {
    return(.regularPath(rep(Inf, length(naturalStartDay)), settings))
  }
  graceEnd <- strategy$graceEnd
  if (strategy$type == "notInitiateWithin") {
    # No start on or before graceEnd, then the natural initiation hazard from day graceEnd + 1 on:
    startDay <- ifelse(
      naturalStartDay > graceEnd,
      naturalStartDay,
      graceEnd + ceiling(laterStartE / .dailyCumulativeHazard(eta$initiation))
    )
    return(.regularPath(startDay, settings))
  }
  if (graceRegime == "naturalGrace") {
    return(.regularPath(pmin(naturalStartDay, graceEnd), settings))
  }
  # binaryAtG: early starters keep their start day. The others follow the untreated hazard through
  # graceEnd, and then the hazard of a start on a day drawn from the starters by graceEnd.
  path <- .regularPath(naturalStartDay, settings)
  late <- which(naturalStartDay > graceEnd)
  etaInitiationDay0 <- if (.hasDay0Initiation(settings)) eta$initiationDay0[late] else NULL
  drawnStartDay <- .drawGraceStartDay(u[late], eta$initiation[late], eta$outcome[late], graceEnd, settings, etaInitiationDay0)
  latePath <- .regularPath(drawnStartDay, settings)
  path$switchDay[late] <- pmax(latePath$switchDay, graceEnd)
  path$startDay[late] <- drawnStartDay
  return(path)
}

#' Simulate clone-censor-weight data
#'
#' @description
#' Simulates a [CloneCensorWeightData] object from a daily discrete-time data-generating process
#' (DGP) with confounded treatment initiation, a treatment-dependent outcome, competing death, and
#' loss to follow-up. [computeTrueRisks()] gives the true risks of treatment strategies under the
#' same DGP.
#'
#' @details
#' Days are counted from the index date (day 0) up to `maxDays`. Each person has `nCovariates`
#' baseline covariates X: covariate j is binary (with prevalence `binaryPrevalence`) when j is odd and
#' standard normal when j is even. On each day, the daily probabilities (hazards) are, on the logit
#' scale:
#'
#' - initiation of the treatment, for persons who have not started: `qlogis(initiationHazard) + X %*%
#'   initiationCoefficients`, except on day 0 when `initiationHazardDay0` or
#'   `initiationCoefficientsDay0` is given: `qlogis(initiationHazardDay0) + X %*%
#'   initiationCoefficientsDay0` (each `NULL` one takes the value of the other days);
#' - the outcome: `qlogis(outcomeHazard) + X %*% outcomeCoefficients`, plus the treatment effect below;
#' - death: `qlogis(deathHazard) + X %*% deathCoefficients` (not affected by the treatment);
#' - end of observation: `observationEndHazard`, independent of X.
#'
#' The treatment effect on the logit of the daily outcome hazard depends on `dgp`:
#'
#' - `"a"`: `treatmentEffect` on every day after `graceEnd`, for persons who started the treatment on
#'   or before day `graceEnd`; no effect otherwise. The effect acts only through "started by day
#'   `graceEnd`", so the truth of an `initiateWithin(graceEnd)` strategy does not depend on when
#'   within the grace period a person starts.
#' - `"b"`: `treatmentEffect * min(s, effectRampDays) / effectRampDays` on day `A + s`, `s >= 1`, for a
#'   person who started on day `A`. The effect starts the day after initiation and builds up with time
#'   since initiation, so the start day within a grace period changes later risk.
#'
#' Without the day-0 parameters, the daily initiation hazard is constant over time, so the
#' probability of starting by day `G`, `1 - (1 - plogis(qlogis(initiationHazard) + X %*%
#' initiationCoefficients))^(G + 1)`, is not linear in X on the logit scale. The binary censoring model
#' of an `initiateWithin(G)` strategy, a logistic regression on X fitted by [fitCensoringModels()], is
#' therefore misspecified on these data. Its bias is negligible when initiation is weakly confounded
#' or rare, but it can be detected with strong covariate effects on initiation and a high daily
#' initiation hazard. The pooled logistic initiation hazard models of `"notInitiateWithin"` and
#' `"neverInitiate"` strategies are correctly specified. With `initiationCoefficientsDay0` different
#' from `initiationCoefficients`, the covariates change *when* persons start (on day 0 or later), and
#' the initiation hazard model is correctly specified only with
#' `createFitCensoringModelsArgs(covariateTimeInteraction = "day0")` (or `"interval"`).
#'
#' Because the daily initiation hazard is the same on every day, the probability of not starting on
#' days 0 to `d` is `(1 - h)^(d + 1)` for a daily hazard `h`: with the default `initiationHazard =
#' 0.02`, about 0.06% for persons with X = 0 by day 365, and less for those with a higher hazard. A
#' `"neverInitiate"` strategy, or a `"notInitiateWithin"` strategy with a long grace period, then has
#' almost no compliant persons late in follow-up, and its weights are extreme or undefined (a
#' positivity violation). For such strategies, use a lower `initiationHazard` or a shorter `maxDays`.
#'
#' On one day, the outcome comes before death. The outcome, death, and the treatment start are
#' recorded when they fall on or before the end of observation, which is the first of the observation
#' end, death, and `maxDays`. Each person has one index date, one treatment era (from the start to the
#' end of observation), and at most one outcome. The covariates are stored in FeatureExtraction's sparse
#' format: binary covariates (analysis ID 1) only where they are 1, normal covariates (analysis ID 2)
#' for every person, with covariate ID `1000 * j + analysisId` for covariate j.
#'
#' The target cohort ID is 1, the treatment cohort ID is 2 (role `"treatment"`), and the outcome ID
#' is 3. All arguments except `n` and `seed` are stored in `metaData$simulationSettings`, which is the
#' `settings` argument of [computeTrueRisks()].
#'
#' @param n                      Number of persons.
#' @param nCovariates            Number of baseline covariates.
#' @param dgp                    The DGP family: `"a"` or `"b"` (see details).
#' @param graceEnd               The last day of the grace period of DGP `"a"`; the treatment only acts
#'                               after this day. Not used by DGP `"b"`.
#' @param maxDays                The last day (since index) that is simulated.
#' @param binaryPrevalence       Prevalence of the binary covariates.
#' @param initiationHazard       Daily probability of starting the treatment when X = 0.
#' @param outcomeHazard          Daily probability of the outcome when X = 0 and untreated.
#' @param deathHazard            Daily probability of death when X = 0.
#' @param observationEndHazard   Daily probability that observation ends.
#' @param initiationCoefficients Log odds ratios of the daily initiation hazard for the covariates.
#' @param outcomeCoefficients    Log odds ratios of the daily outcome hazard for the covariates.
#' @param deathCoefficients      Log odds ratios of the daily death hazard for the covariates.
#' @param treatmentEffect        The (full) effect of the treatment on the logit of the daily outcome
#'                               hazard.
#' @param effectRampDays         For DGP `"b"`, the number of days since initiation it takes the effect
#'                               to reach `treatmentEffect`.
#' @param initiationHazardDay0   Daily probability of starting the treatment on day 0 when X = 0, or
#'                               `NULL` for `initiationHazard`.
#' @param initiationCoefficientsDay0 Log odds ratios of the initiation hazard of day 0 for the
#'                               covariates, or `NULL` for `initiationCoefficients`.
#' @param seed                   The seed for the random number generator, or `NULL` to leave it as it
#'                               is.
#'
#' @return
#' An object of type [CloneCensorWeightData], with `metaData$simulationSettings`.
#'
#' @examples
#' data <- simulateCloneCensorWeightData(n = 1000, seed = 1)
#' summary(data)
#' close(data)
#'
#' @export
simulateCloneCensorWeightData <- function(n = 5000,
                                          nCovariates = 5,
                                          dgp = "a",
                                          graceEnd = 30,
                                          maxDays = 365,
                                          binaryPrevalence = 0.3,
                                          initiationHazard = 0.02,
                                          outcomeHazard = 5e-04,
                                          deathHazard = 2e-04,
                                          observationEndHazard = 0.001,
                                          initiationCoefficients = rep(0.4, nCovariates),
                                          outcomeCoefficients = rep(0.4, nCovariates),
                                          deathCoefficients = rep(0.2, nCovariates),
                                          treatmentEffect = log(0.5),
                                          effectRampDays = 30,
                                          initiationHazardDay0 = NULL,
                                          initiationCoefficientsDay0 = NULL,
                                          seed = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertCount(n, positive = TRUE, add = errorMessages)
  checkmate::assertCount(nCovariates, add = errorMessages)
  checkmate::assertInt(seed, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  settings <- list(
    dgp = dgp,
    nCovariates = nCovariates,
    graceEnd = graceEnd,
    maxDays = maxDays,
    binaryPrevalence = binaryPrevalence,
    initiationHazard = initiationHazard,
    outcomeHazard = outcomeHazard,
    deathHazard = deathHazard,
    observationEndHazard = observationEndHazard,
    initiationCoefficients = initiationCoefficients,
    outcomeCoefficients = outcomeCoefficients,
    deathCoefficients = deathCoefficients,
    treatmentEffect = treatmentEffect,
    effectRampDays = effectRampDays,
    initiationHazardDay0 = initiationHazardDay0,
    initiationCoefficientsDay0 = initiationCoefficientsDay0
  )
  .assertSimulationSettings(settings)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  x <- .simulateCovariates(n, settings)
  eta <- .linearPredictors(x, settings)
  maxDays <- settings$maxDays
  startDay <- .initiationEventDay(stats::rexp(n), eta, settings, maxDays)
  outcomeDay <- .outcomeEventDay(stats::rexp(n), eta$outcome, .regularPath(startDay, settings), settings)
  deathDay <- .constantHazardEventDay(stats::rexp(n), .dailyCumulativeHazard(eta$death), maxDays)
  observationEndDay <- .constantHazardEventDay(stats::rexp(n), .hazardIncrement(settings$observationEndHazard), maxDays)
  daysToObsEnd <- pmin(observationEndDay, deathDay, maxDays)

  rowIds <- seq_len(n)
  cohorts <- data.frame(
    rowId = rowIds,
    personSeqId = rowIds,
    personId = as.character(rowIds),
    cohortStartDate = as.Date("2020-01-01") + (rowIds - 1) %% 365,
    daysFromObsStart = 365L,
    daysToCohortEnd = 0L,
    daysToObsEnd = as.integer(daysToObsEnd),
    daysToDeath = ifelse(deathDay <= daysToObsEnd, deathDay, NA_integer_)
  )
  idx <- which(outcomeDay <= daysToObsEnd)
  outcomes <- data.frame(rowId = idx, outcomeId = rep(.simulationOutcomeId, length(idx)), daysToEvent = outcomeDay[idx])
  idx <- which(startDay <= daysToObsEnd)
  treatments <- data.frame(
    rowId = idx,
    treatmentId = rep(.simulationTreatmentId, length(idx)),
    daysToStart = startDay[idx],
    daysToEnd = daysToObsEnd[idx]
  )

  isBinary <- .covariateIsBinary(nCovariates)
  analysisIds <- ifelse(isBinary, 1, 2)
  covariateIds <- 1000 * seq_len(nCovariates) + analysisIds
  covariates <- data.frame(
    rowId = rep(rowIds, nCovariates),
    covariateId = rep(covariateIds, each = n),
    covariateValue = as.vector(x)
  )
  covariates <- covariates[covariates$covariateValue != 0, ]
  covariateRef <- data.frame(
    covariateId = covariateIds,
    covariateName = sprintf("Simulated %s covariate %d", ifelse(isBinary, "binary", "normal"), seq_len(nCovariates)),
    analysisId = analysisIds,
    conceptId = seq_len(nCovariates)
  )
  analysisRef <- data.frame(
    analysisId = c(1, 2),
    analysisName = c("SimulatedBinary", "SimulatedNormal"),
    domainId = "Simulated",
    startDay = -365,
    endDay = 0,
    isBinary = c("Y", "N"),
    missingMeansZero = c("Y", "N")
  )
  analysisRef <- analysisRef[analysisRef$analysisId %in% analysisIds, ]

  return(createCloneCensorWeightData(
    cohorts = cohorts,
    outcomes = outcomes,
    treatments = treatments,
    covariates = covariates,
    covariateRef = covariateRef,
    analysisRef = analysisRef,
    metaData = list(
      targetId = .simulationTargetId,
      treatmentIds = c(treatment = .simulationTreatmentId),
      outcomeIds = .simulationOutcomeId,
      simulationSettings = settings
    )
  ))
}

.assertSimulationSettings <- function(settings) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertList(settings, add = errorMessages)
  checkmate::assertNames(names(settings), must.include = .simulationSettingsFields, .var.name = "names(settings)", add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertChoice(settings$dgp, c("a", "b"), .var.name = "dgp", add = errorMessages)
  checkmate::assertCount(settings$nCovariates, .var.name = "nCovariates", add = errorMessages)
  checkmate::assertInt(settings$graceEnd, lower = 0, .var.name = "graceEnd", add = errorMessages)
  checkmate::assertInt(settings$maxDays, lower = 0, .var.name = "maxDays", add = errorMessages)
  checkmate::assertNumber(settings$binaryPrevalence, lower = 0, upper = 1, .var.name = "binaryPrevalence", add = errorMessages)
  checkmate::assertNumber(settings$initiationHazard, lower = 0, upper = 1, .var.name = "initiationHazard", add = errorMessages)
  if (isTRUE(settings$initiationHazard %in% c(0, 1))) {
    errorMessages$push("initiationHazard must be strictly between 0 and 1")
  }
  for (name in c("outcomeHazard", "deathHazard", "observationEndHazard")) {
    checkmate::assertNumber(settings[[name]], lower = 0, upper = 1, .var.name = name, add = errorMessages)
    if (isTRUE(settings[[name]] == 1)) {
      errorMessages$push(sprintf("%s must be smaller than 1", name))
    }
  }
  for (name in c("initiationCoefficients", "outcomeCoefficients", "deathCoefficients")) {
    checkmate::assertNumeric(settings[[name]], finite = TRUE, any.missing = FALSE, len = settings$nCovariates, .var.name = name, add = errorMessages)
  }
  checkmate::assertNumber(settings$treatmentEffect, finite = TRUE, .var.name = "treatmentEffect", add = errorMessages)
  checkmate::assertInt(settings$effectRampDays, lower = 1, .var.name = "effectRampDays", add = errorMessages)
  # Optional (NULL: the same as on the other days), so settings created before they existed are valid:
  checkmate::assertNumber(settings$initiationHazardDay0, lower = 0, upper = 1, null.ok = TRUE, .var.name = "initiationHazardDay0", add = errorMessages)
  if (isTRUE(settings$initiationHazardDay0 %in% c(0, 1))) {
    errorMessages$push("initiationHazardDay0 must be strictly between 0 and 1")
  }
  checkmate::assertNumeric(settings$initiationCoefficientsDay0, finite = TRUE, any.missing = FALSE, len = settings$nCovariates, null.ok = TRUE, .var.name = "initiationCoefficientsDay0", add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  invisible(settings)
}

#' Compute the true risks of treatment strategies
#'
#' @description
#' Computes the true risk of the outcome (the cause-specific cumulative incidence, with death as a
#' competing event) under each treatment strategy, for the data-generating process of
#' [simulateCloneCensorWeightData()], by simulating each strategy directly on `nTruth` new persons.
#'
#' @details
#' The risks are always the cumulative incidence of the outcome with death as a competing event, the
#' estimand of an analysis with `competingEvents = "death"`. An analysis with `competingEvents =
#' "none"` treats death as censoring, so it estimates the risk in a world without death. Its truth is
#' given by this function with `settings$deathHazard <- 0`. Death is not weighted for, so the analysis
#' identifies this risk only when death is non-informative censoring: in this DGP, when death does not
#' depend on X (`deathCoefficients = 0`) or X does not change the outcome (`outcomeCoefficients = 0`).
#'
#' The same random numbers are used for every strategy (and for every value of `graceRegime` given the
#' same `seed`), so contrasts between strategies are estimated more precisely than each risk. There is
#' no loss to follow-up. Each risk is a proportion over `nTruth` independent persons, so its Monte Carlo
#' standard error is `sqrt(risk * (1 - risk) / nTruth)`.
#'
#' The strategies are simulated as follows, with A the start day of the treatment and G the strategy's
#' `graceEnd`:
#'
#' - `neverInitiate`: nobody starts.
#' - `notInitiateWithin(G)`: nobody starts on or before day G; after day G, the natural initiation
#'   hazard applies.
#' - `initiateWithin(G)` with `graceRegime = "binaryAtG"`: the natural course through day G. Each
#'   person with no outcome or death on or before day G who has not started gets a start day a drawn
#'   from P(A = a | A <= G, no outcome or death on or before day G, X); early starters keep their own
#'   start day. After day G, the outcome hazard follows that start day. This is the regime the
#'   binary-at-G censoring weight identifies.
#' - `initiateWithin(G)` with `graceRegime = "naturalGrace"`: the natural course through day G - 1;
#'   persons still untreated on day G start on day G.
#'
#' Under DGP `"a"`, with G equal to the DGP's `graceEnd`, the two grace regimes give the same truth.
#' Under DGP `"b"`, the binary-at-G regime moves the start of persons who would not start by G to
#' earlier days, where the effect has had longer to build up.
#'
#' The start-day probabilities P(A = a | X) use the day-0 initiation hazard when the settings have
#' day-0 initiation parameters. Under DGP `"a"`, the truths of strategies with G equal to the DGP's
#' `graceEnd` (and of `neverInitiate`) do not depend on the initiation process at all, and those of
#' strategies with G up to `graceEnd` do not depend on the day-0 parameters. Under DGP `"b"`, the
#' `initiateWithin` truths depend on them under both grace regimes.
#'
#' With `riskStart = 0` the risks count outcomes from day 0 (the index date), the estimands of analyses
#' with `riskStart = 0`. An analysis with `riskStart = r > 0` (with CohortMethod's default removal of
#' persons with the outcome before `riskStart`, and with the clones of persons who die before
#' `riskStart` dropped) estimates the risk on days r to h among persons with no outcome or death on days
#' 0 to r - 1. With `riskStart = r`, the risk at horizon h is therefore
#' `(risk(h) - risk(r - 1)) / (1 - p)`, where `risk()` is the risk counted from day 0 and `p` the
#' probability of the outcome or death on days 0 to r - 1, the mean over X of
#' `1 - ((1 - hY(X)) (1 - hD(X)))^r`, with `hY` and `hD` the untreated daily outcome and death hazards.
#' This holds when the treatment changes no hazard before day r under any strategy: always for
#' `riskStart = 1`, and for `riskStart <= graceEnd + 1` under DGP `"a"`. Other values stop with an
#' error.
#'
#' @param settings    The simulation settings: `metaData$simulationSettings` of an object created by
#'                    [simulateCloneCensorWeightData()].
#' @param strategies  A list of objects of type `TreatmentStrategy` as created using
#'                    [createTreatmentStrategy()], all referring to the `"treatment"` role, with
#'                    `graceStart = 0` and no `prohibitedTreatments`. Their position in the list is the
#'                    `strategyId` of the result.
#' @param horizons    The days (since index) at which to compute the risks, from `riskStart` to
#'                    `maxDays`. The risk at day h includes outcomes on day h.
#' @param graceRegime How the `initiateWithin` strategies are defined: `"binaryAtG"` or `"naturalGrace"`
#'                    (see details). Other strategy types ignore it.
#' @param nTruth      The number of persons to simulate.
#' @param seed        The seed for the random number generator, or `NULL` to leave it as it is.
#' @param riskStart   The first day of the risk window of the analysis (`riskWindowStart`): the risks
#'                    are those of outcomes from day `riskStart` on, among persons with no outcome or
#'                    death before it (see details).
#'
#' @return
#' A tibble with columns `strategyId` (the position in `strategies`), `day` (the horizon), and `risk`.
#'
#' @examples
#' data <- simulateCloneCensorWeightData(n = 100, seed = 1)
#' settings <- attr(data, "metaData")$simulationSettings
#' close(data)
#' strategies <- list(
#'   createTreatmentStrategy(1, "Initiate within 30 days", "initiateWithin", graceEnd = 30),
#'   createTreatmentStrategy(2, "Do not initiate within 30 days", "notInitiateWithin", graceEnd = 30)
#' )
#' computeTrueRisks(settings, strategies, horizons = c(90, 365), nTruth = 1e5, seed = 1)
#'
#' @export
computeTrueRisks <- function(settings,
                             strategies,
                             horizons,
                             graceRegime = "binaryAtG",
                             nTruth = 1e6,
                             seed = NULL,
                             riskStart = 0) {
  .assertSimulationSettings(settings)
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertList(strategies, min.len = 1, add = errorMessages)
  for (i in seq_along(strategies)) {
    checkmate::assertR6(strategies[[i]], "TreatmentStrategy", .var.name = sprintf("strategies[[%d]]", i), add = errorMessages)
  }
  checkmate::assertIntegerish(horizons, lower = 0, upper = settings$maxDays, any.missing = FALSE, min.len = 1, unique = TRUE, add = errorMessages)
  checkmate::assertChoice(graceRegime, c("binaryAtG", "naturalGrace"), add = errorMessages)
  checkmate::assertCount(nTruth, positive = TRUE, add = errorMessages)
  checkmate::assertInt(seed, null.ok = TRUE, add = errorMessages)
  checkmate::assertInt(riskStart, lower = 0, upper = settings$maxDays, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (any(horizons < riskStart)) {
    stop("Horizons must be on or after riskStart")
  }
  # Day of the first hazard the treatment can change, under any strategy:
  firstTreatedDay <- if (settings$dgp == "a") settings$graceEnd + 1 else 1
  if (riskStart > firstTreatedDay) {
    stop(sprintf(
      "riskStart = %d is not supported: the treatment can change hazards from day %d under DGP '%s', so persons at risk on day riskStart differ between strategies",
      riskStart, firstTreatedDay, settings$dgp
    ))
  }
  for (strategy in strategies) {
    if (strategy$graceStart > 0) {
      stop("Not yet supported in this version: graceStart > 0")
    }
    if (length(strategy$prohibitedTreatments) > 0) {
      stop("Not yet supported in this version: prohibitedTreatments")
    }
    if (strategy$treatment != "treatment") {
      stop("The simulation has only the 'treatment' role, but a strategy refers to '", strategy$treatment, "'")
    }
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }
  # All random numbers are drawn up front, in a fixed order, for common random numbers:
  x <- .simulateCovariates(nTruth, settings)
  initiationE <- stats::rexp(nTruth)
  outcomeE <- stats::rexp(nTruth)
  deathE <- stats::rexp(nTruth)
  laterStartE <- stats::rexp(nTruth)
  u <- stats::runif(nTruth)

  eta <- .linearPredictors(x, settings)
  maxDays <- settings$maxDays
  naturalStartDay <- .initiationEventDay(initiationE, eta, settings, Inf)
  deathDay <- .constantHazardEventDay(deathE, .dailyCumulativeHazard(eta$death), maxDays)

  # Probability of the outcome or death on days 0..riskStart - 1 (untreated hazards):
  pBefore <- 0
  if (riskStart > 0) {
    noEventDay <- (1 - stats::plogis(eta$outcome)) * (1 - stats::plogis(eta$death))
    pBefore <- mean(1 - noEventDay^riskStart)
  }
  results <- list()
  for (i in seq_along(strategies)) {
    path <- .strategyPath(strategies[[i]], graceRegime, naturalStartDay, laterStartE, u, eta, settings)
    outcomeDay <- .outcomeEventDay(outcomeE, eta$outcome, path, settings)
    outcomeDay <- outcomeDay[outcomeDay <= deathDay & is.finite(outcomeDay)]
    risk <- cumsum(tabulate(outcomeDay + 1, nbins = maxDays + 1)) / nTruth
    riskBefore <- if (riskStart > 0) risk[riskStart] else 0
    results[[i]] <- tibble(
      strategyId = as.integer(i),
      day = as.integer(horizons),
      risk = (risk[horizons + 1] - riskBefore) / (1 - pBefore)
    )
  }
  return(dplyr::bind_rows(results))
}

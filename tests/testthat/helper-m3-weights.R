# Fixtures and hand-written references for the M3 (censoring models and weights) tests.

# Persons (contract section 5.2 `persons`) and baseline covariates for `n` persons. Three binary
# covariates (each in its own analysis, so FeatureExtraction's redundancy rule never removes one).
# Initiation follows a daily hazard that depends on the covariates; outcomes and deaths happen on
# random days of the risk window. With `untiedInitiation = TRUE`, initiation days are distinct.
createM3Fixture <- function(n = 150,
                            riskStart = 0L,
                            seed = 1,
                            nominalRiskEnd = 40L,
                            maxTreatmentDay = 30L,
                            untiedInitiation = FALSE,
                            outcomeFraction = 0.3,
                            deathFraction = 0.1) {
  set.seed(seed)
  x <- cbind(rbinom(n, 1, 0.4), rbinom(n, 1, 0.5), rbinom(n, 1, 0.3))
  hazard <- plogis(-2.3 + 0.9 * x[, 1] - 0.7 * x[, 2] + 0.5 * x[, 3])
  treatmentDay <- rgeom(n, hazard)
  if (untiedInitiation) {
    treatmentDay <- sample(0:(5 * n), n)
    treatmentDay <- as.integer(round(treatmentDay * (1 + x[, 1]) / 6))
    treatmentDay[duplicated(treatmentDay)] <- NA
  }
  treatmentDay[!is.na(treatmentDay) & treatmentDay > maxTreatmentDay] <- NA
  obsEnd <- sample(5:(nominalRiskEnd + 30L), n, replace = TRUE)
  riskEnd <- pmin(nominalRiskEnd, obsEnd)
  drawDay <- function(fraction) {
    day <- riskStart + floor(runif(n) * (riskEnd - riskStart + 1))
    day[runif(n) > fraction] <- NA
    as.integer(day)
  }
  outcomeDay <- drawDay(outcomeFraction)
  deathDay <- drawDay(deathFraction)
  deathDay[!is.na(deathDay) & !is.na(outcomeDay) & deathDay == outcomeDay] <- NA
  persons <- data.frame(
    rowId = seq_len(n),
    personSeqId = seq_len(n),
    riskStart = as.integer(riskStart),
    riskEnd = as.integer(riskEnd),
    nominalRiskEnd = as.integer(nominalRiskEnd),
    outcomeDay = outcomeDay,
    deathDay = deathDay,
    treatmentDay = as.integer(treatmentDay),
    comparatorTreatmentDay = NA_integer_
  )
  covariates <- data.frame(
    rowId = rep(seq_len(n), 3),
    covariateId = rep(c(1001, 1002, 1003), each = n),
    covariateValue = as.vector(x)
  )
  covariates <- covariates[covariates$covariateValue != 0, ]
  return(list(persons = persons, covariates = covariates, x = x))
}

createM3Data <- function(fixture) {
  n <- nrow(fixture$persons)
  cohorts <- data.frame(
    rowId = fixture$persons$rowId,
    personSeqId = fixture$persons$personSeqId,
    personId = as.character(1000 + fixture$persons$rowId),
    cohortStartDate = as.Date("2020-01-01"),
    daysFromObsStart = 400L,
    daysToCohortEnd = 0L,
    daysToObsEnd = fixture$persons$riskEnd,
    daysToDeath = fixture$persons$deathDay
  )
  covariateRef <- data.frame(
    covariateId = c(1001, 1002, 1003),
    covariateName = paste("Covariate", 1:3),
    analysisId = c(1, 2, 3),
    conceptId = 0
  )
  analysisRef <- data.frame(
    analysisId = c(1, 2, 3),
    analysisName = paste("Analysis", 1:3),
    domainId = "Condition",
    startDay = -365,
    endDay = 0,
    isBinary = "Y",
    missingMeansZero = "Y"
  )
  return(createCloneCensorWeightData(
    cohorts = cohorts,
    outcomes = NULL,
    treatments = NULL,
    covariates = fixture$covariates,
    covariateRef = covariateRef,
    analysisRef = analysisRef,
    metaData = list(targetId = 1, treatmentIds = c(treatment = 21))
  ))
}

# The clone rules of contract section 5.2, written out per person.
buildM3CloneTable <- function(persons, strategies, sameDayOrdering = "eventsFirst") {
  rows <- list()
  for (strategyId in seq_along(strategies)) {
    strategy <- strategies[[strategyId]]
    for (i in seq_len(nrow(persons))) {
      p <- persons[i, ]
      a <- if (strategy$treatment == "treatment") p$treatmentDay else p$comparatorTreatmentDay
      g <- strategy$graceEnd
      deviationDay <- switch(strategy$type,
        initiateWithin = if (is.na(a) || a > g) g else NA,
        notInitiateWithin = if (!is.na(a) && a <= g) a else NA,
        neverInitiate = a
      )
      candidates <- c(outcome = p$outcomeDay, death = p$deathDay, deviation = deviationDay, end = p$riskEnd)
      candidates <- candidates[!is.na(candidates)]
      endDay <- min(candidates)
      tied <- names(candidates)[candidates == endDay]
      order <- if (sameDayOrdering == "eventsFirst") {
        c("outcome", "death", "deviation", "end")
      } else {
        c("deviation", "outcome", "death", "end")
      }
      reason <- order[order %in% tied][1]
      if (reason == "end") {
        reason <- if (p$riskEnd < p$nominalRiskEnd) "observationEnd" else "riskEnd"
      }
      if (endDay < p$riskStart) {
        next
      }
      rows[[length(rows) + 1]] <- data.frame(
        rowId = p$rowId,
        personSeqId = p$personSeqId,
        strategyId = strategyId,
        riskStart = p$riskStart,
        deviationDay = as.integer(deviationDay),
        endDay = as.integer(endDay),
        survivalTime = as.integer(endDay - p$riskStart + 1),
        endReason = reason,
        outcomeCount = as.integer(reason == "outcome")
      )
    }
  }
  clones <- do.call(rbind, rows)
  clones <- clones[order(clones$strategyId, clones$rowId), ]
  clones <- cbind(cloneRowId = seq_len(nrow(clones)), clones)
  rownames(clones) <- NULL
  clones$strategyId <- as.integer(clones$strategyId)
  clones$cloneRowId <- as.integer(clones$cloneRowId)
  return(clones)
}

buildM3Clones <- function(persons, strategies, sameDayOrdering = "eventsFirst", competingEvents = "death") {
  clones <- buildM3CloneTable(persons, strategies, sameDayOrdering)
  metaData <- list(
    strategies = lapply(strategies, function(strategy) strategy$toList()),
    outcomeId = 3,
    treatmentIds = c(treatment = 21),
    sameDayOrdering = sameDayOrdering,
    competingEvents = competingEvents,
    attrition = dplyr::tibble(
      description = "Original cohorts",
      targetPersons = nrow(persons),
      comparatorPersons = 0L,
      targetExposures = nrow(persons),
      comparatorExposures = 0L
    ),
    targetId = 1
  )
  return(.createCcwClones(persons, clones, metaData))
}

m3Strategies <- function(first = "initiateWithin", second = "notInitiateWithin", graceEnd = 7) {
  make <- function(id, type) {
    createTreatmentStrategy(
      strategyId = id,
      label = paste(type, id),
      type = type,
      graceEnd = if (type == "neverInitiate") NULL else graceEnd
    )
  }
  return(list(make(1, first), make(2, second)))
}

# Unpenalized fits with tight convergence, as plan section 6.4 requires for exact comparisons. The M3
# references give notInitiateWithin clones W = 1 / S(G | X) after G (postGraceModel = "hazard"); the
# binary post-grace weights (contract section 11.2) are tested in test-m3-weights-postGrace.R.
m3Args <- function(hazardModelType = "cox", useCovariates = TRUE, timeBasis = createTimeBasis(), prior = Cyclops::createPrior("none"), postGraceModel = "hazard") {
  return(createFitCensoringModelsArgs(
    hazardModelType = hazardModelType,
    timeBasis = timeBasis,
    postGraceModel = postGraceModel,
    useCovariates = useCovariates,
    errorOnHighCorrelation = FALSE,
    prior = prior,
    control = Cyclops::createControl(tolerance = 1e-10, noiseLevel = "silent", seed = 1, cvType = "auto")
  ))
}

m3CoxControl <- function() {
  return(survival::coxph.control(eps = 1e-12, toler.chol = 1e-13))
}

# Hand-written risk sets ----------------------------------------------------------------------------

# Persons under follow-up at day g (methods section 2.4), with y = started by g.
refBinaryRiskSet <- function(persons, g, sameDayOrdering) {
  keep <- logical(nrow(persons))
  y <- integer(nrow(persons))
  for (i in seq_len(nrow(persons))) {
    p <- persons[i, ]
    eventDays <- c(p$outcomeDay, p$deathDay)
    eventDays <- eventDays[!is.na(eventDays)]
    eventBefore <- if (sameDayOrdering == "eventsFirst") any(eventDays <= g) else any(eventDays < g)
    keep[i] <- !eventBefore && p$riskEnd >= g
    y[i] <- as.integer(!is.na(p$treatmentDay) && p$treatmentDay <= g)
  }
  return(data.frame(rowId = persons$rowId[keep], y = y[keep]))
}

# Day-by-day risk set of initiation through horizon h: one row per person and day at risk, with
# initiated = 1 on the day of initiation (methods section 2.4, contract section 5.3).
refHazardDays <- function(persons, h, sameDayOrdering) {
  rows <- list()
  for (i in seq_len(nrow(persons))) {
    p <- persons[i, ]
    for (day in 0:h) {
      if (day > p$riskEnd) break
      eventToday <- any(c(p$outcomeDay, p$deathDay) == day, na.rm = TRUE)
      eventBefore <- any(c(p$outcomeDay, p$deathDay) < day, na.rm = TRUE)
      if (eventBefore) break
      if (eventToday && sameDayOrdering == "eventsFirst") break
      initiated <- !is.na(p$treatmentDay) && p$treatmentDay == day
      rows[[length(rows) + 1]] <- data.frame(rowId = p$rowId, day = day, initiated = as.integer(initiated))
      if (initiated || eventToday) break
    }
  }
  return(do.call(rbind, rows))
}

# One row per person (Cox): time = last day at risk + 1.
refCoxRows <- function(persons, h, sameDayOrdering) {
  days <- refHazardDays(persons, h, sameDayOrdering)
  last <- days[!duplicated(days$rowId, fromLast = TRUE), ]
  return(data.frame(rowId = last$rowId, time = last$day + 1, initiated = last$initiated))
}

# One row per person and interval of length `intervalLength` (days counted from 0) the person enters.
refIntervalRows <- function(persons, h, sameDayOrdering, intervalLength) {
  days <- refHazardDays(persons, h, sameDayOrdering)
  days$k <- days$day %/% intervalLength
  rows <- stats::aggregate(initiated ~ rowId + k, data = days, FUN = max)
  rows <- rows[order(rows$rowId, rows$k), ]
  rownames(rows) <- NULL
  return(rows)
}

# Reference survival and weights ------------------------------------------------------------------

# S(d | X) for days `days` (a matrix person x day) from a coxph fit and survfit(stype = 2, ctype = 1),
# read at t = d + 1 (S(-1) = 1).
refCoxSurvival <- function(coxFit, newdata, days) {
  sf <- survival::survfit(coxFit, newdata = newdata, stype = 2, ctype = 1)
  out <- matrix(1, nrow = nrow(newdata), ncol = length(days))
  for (j in seq_along(days)) {
    if (days[j] >= 0) {
      out[, j] <- summary(sf, times = days[j] + 1, extend = TRUE)$surv
    }
  }
  return(out)
}

# Matrix of covariates in fixture order.
m3CovariateMatrix <- function(fixture, rowIds) {
  x <- fixture$x[rowIds, , drop = FALSE]
  colnames(x) <- c("x1", "x2", "x3")
  return(x)
}

# The clones' follow-up in the outcome risk set: with deviationFirst, a clone that deviates on its last
# day does so before any outcome that day, so that day is left out (and a clone left with no follow-up
# has no segments).
refOutcomeFollowUp <- function(clones, sameDayOrdering) {
  if (sameDayOrdering == "deviationFirst") {
    dropLastDay <- !is.na(clones$deviationDay) & clones$deviationDay == clones$endDay
    clones$survivalTime[dropLastDay] <- clones$survivalTime[dropLastDay] - 1L
  }
  return(clones[clones$survivalTime > 0, ])
}

# Segments by survSplit: each clone split at its own t-scale cuts.
refSegments <- function(clones, cutsByClone) {
  # survSplit needs the bare name Surv in its formula:
  Surv <- survival::Surv
  out <- list()
  for (i in seq_len(nrow(clones))) {
    clone <- clones[i, ]
    data <- data.frame(
      cloneRowId = clone$cloneRowId,
      startTime = 0,
      stopTime = clone$survivalTime,
      status = as.integer(clone$endReason %in% c("outcome", "death")),
      endReason = clone$endReason
    )
    cuts <- cutsByClone[[as.character(clone$cloneRowId)]]
    split <- survival::survSplit(Surv(startTime, stopTime, status) ~ ., data = data, cut = cuts, start = "startTime", end = "stopTime")
    out[[i]] <- split
  }
  segments <- do.call(rbind, out)
  segments$outcome <- as.integer(segments$status == 1 & segments$endReason == "outcome")
  segments$competing <- as.integer(segments$status == 1 & segments$endReason == "death")
  segments <- segments[order(segments$cloneRowId, segments$startTime), c("cloneRowId", "startTime", "stopTime", "outcome", "competing")]
  rownames(segments) <- NULL
  return(segments)
}

# The t-scale cuts of each clone (contract section 5.4): G - riskStart + 1 for initiateWithin (G - riskStart
# with deviationFirst), and the distinct outcome and competing-event times of all clones inside the
# window where W varies for the hazard strategies. Contract section 9, amendment 5: also the end of day
# tau + 1 (tau the largest graceEnd; riskStart if later), and with `cutIntervalLength` every multiple of
# it (the outcome model's grid).
refCuts <- function(clones, strategies, sameDayOrdering, horizons, cutIntervalLength = NULL) {
  eventTimes <- sort(unique(clones$survivalTime[clones$endReason %in% c("outcome", "death")]))
  graceEnds <- unlist(lapply(strategies, function(strategy) strategy$graceEnd))
  cuts <- list()
  for (i in seq_len(nrow(clones))) {
    clone <- clones[i, ]
    strategy <- strategies[[clone$strategyId]]
    if (strategy$type == "initiateWithin") {
      cutDay <- if (sameDayOrdering == "eventsFirst") strategy$graceEnd else strategy$graceEnd - 1
      c <- cutDay - clone$riskStart + 1
    } else {
      windowEnd <- horizons[[clone$strategyId]]
      c <- eventTimes[eventTimes + clone$riskStart - 1 <= windowEnd]
    }
    if (length(graceEnds) > 0) {
      c <- c(c, max(max(graceEnds) + 1, clone$riskStart) - clone$riskStart + 1)
    }
    if (!is.null(cutIntervalLength)) {
      c <- c(c, seq(cutIntervalLength, clone$survivalTime, by = cutIntervalLength))
    }
    c <- sort(unique(c))
    cuts[[as.character(clone$cloneRowId)]] <- c[c > 0 & c < clone$survivalTime]
  }
  return(cuts)
}

collectSegments <- function(weighted) {
  segments <- as.data.frame(dplyr::collect(weighted$segments))
  segments <- segments[order(segments$cloneRowId, segments$startTime), ]
  rownames(segments) <- NULL
  return(segments)
}

collectClones <- function(clones) {
  x <- as.data.frame(dplyr::collect(clones$clones))
  x <- x[order(x$cloneRowId), ]
  rownames(x) <- NULL
  return(x)
}

# Fixtures and hand-written references for the outcome-model tests (module M4).

# Weighted clones for initiateWithin(grace) (strategy 1) vs notInitiateWithin(grace) (strategy 2),
# built directly from persons following Phase1Contracts sections 5.2 and 5.4 ("eventsFirst"):
# - persons with outcomes, competing deaths (when competingEvents = "death"), early and late
#   initiation, and the last `sharedPersons` rows sharing a personSeqId with the first rows;
# - initiateWithin clones have weight 1 up to day `grace` and 1 / pi after it;
# - notInitiateWithin clones are split at every outcome or death day of all clones within
#   [riskStart, grace] and at day `grace`, with weights growing over the segments;
# - a few extra splits that keep the weight (`extraSplits`).
# `outcomeRate` and `deathRate` are the shares of persons with an outcome or death day drawn uniformly
# from [riskStart, riskStart + maxDays] (kept only within the risk window).
createM4Fixture <- function(n = 120,
                            riskStart = 0,
                            competingEvents = "death",
                            seed = 1,
                            grace = 4,
                            maxDays = 25,
                            sharedPersons = 10,
                            extraSplits = TRUE,
                            outcomeRate = 0.45,
                            deathRate = 0.25) {
  set.seed(seed)
  rowId <- seq_len(n)
  personSeqId <- rowId
  personSeqId[(n - sharedPersons + 1):n] <- seq_len(sharedPersons)
  followUp <- sample(4:maxDays, n, replace = TRUE)
  riskEnd <- riskStart + followUp
  nominalRiskEnd <- riskStart + maxDays
  outcomeDay <- ifelse(runif(n) < outcomeRate, riskStart + sample(0:maxDays, n, replace = TRUE), NA)
  deathDay <- ifelse(runif(n) < deathRate, riskStart + sample(0:maxDays, n, replace = TRUE), NA)
  if (competingEvents == "none") {
    deathDay[] <- NA
  }
  outcomeDay[!is.na(outcomeDay) & outcomeDay > riskEnd] <- NA
  deathDay[!is.na(deathDay) & deathDay > riskEnd] <- NA
  outcomeDay[!is.na(deathDay) & !is.na(outcomeDay) & deathDay < outcomeDay] <- NA
  deathDay[!is.na(deathDay) & !is.na(outcomeDay) & deathDay >= outcomeDay] <- NA
  treatmentDay <- ifelse(runif(n) < 0.6, sample(0:(grace + 4), n, replace = TRUE), NA)
  persons <- data.frame(
    rowId = rowId,
    personSeqId = as.integer(personSeqId),
    riskStart = as.integer(riskStart),
    riskEnd = as.integer(riskEnd),
    nominalRiskEnd = as.integer(nominalRiskEnd),
    outcomeDay = as.integer(outcomeDay),
    deathDay = as.integer(deathDay),
    treatmentDay = as.integer(treatmentDay),
    comparatorTreatmentDay = NA_integer_
  )

  clones <- list()
  for (strategyId in 1:2) {
    if (strategyId == 1) {
      deviationDay <- ifelse(is.na(treatmentDay) | treatmentDay > grace, grace, NA)
    } else {
      deviationDay <- ifelse(!is.na(treatmentDay) & treatmentDay <= grace, treatmentDay, NA)
    }
    candidates <- cbind(outcomeDay, deathDay, deviationDay, riskEnd)
    endDay <- apply(candidates, 1, min, na.rm = TRUE)
    first <- apply(candidates == endDay, 1, function(x) which(x %in% TRUE)[1])
    endReason <- c("outcome", "death", "deviation", "riskEnd")[first]
    endReason[endReason == "riskEnd" & riskEnd < nominalRiskEnd] <- "observationEnd"
    clones[[strategyId]] <- data.frame(
      rowId = rowId,
      personSeqId = as.integer(personSeqId),
      strategyId = as.integer(strategyId),
      riskStart = as.integer(riskStart),
      deviationDay = as.integer(deviationDay),
      endDay = as.integer(endDay),
      endReason = endReason
    )
  }
  clones <- do.call(rbind, clones)
  clones <- clones[clones$endDay >= riskStart, ]
  clones <- clones[order(clones$rowId, clones$strategyId), ]
  clones$cloneRowId <- seq_len(nrow(clones))
  clones$survivalTime <- clones$endDay - clones$riskStart + 1L
  clones$outcomeCount <- as.integer(clones$endReason == "outcome")
  clones <- clones[, c("cloneRowId", "rowId", "personSeqId", "strategyId", "riskStart", "deviationDay", "endDay", "survivalTime", "endReason", "outcomeCount")]
  rownames(clones) <- NULL

  # Segments on the survival time scale t = day - riskStart + 1.
  eventTimes <- sort(unique(c(clones$survivalTime[clones$endReason %in% c("outcome", "death")])))
  graceTime <- grace - riskStart + 1
  piHat <- runif(n, 0.3, 0.9)
  slope <- runif(n, 0, 0.15)
  segments <- list()
  for (i in seq_len(nrow(clones))) {
    clone <- clones[i, ]
    endTime <- clone$survivalTime
    if (clone$strategyId == 1) {
      cuts <- graceTime
    } else {
      cuts <- c(eventTimes[eventTimes <= graceTime], graceTime)
    }
    if (extraSplits && clone$rowId %% 3 == 0) {
      cuts <- c(cuts, round(endTime / 2), endTime - 1)
    }
    cuts <- sort(unique(cuts[cuts > 0 & cuts < endTime]))
    startTime <- c(0, cuts)
    stopTime <- c(cuts, endTime)
    if (clone$strategyId == 1) {
      weight <- ifelse(startTime >= graceTime, 1 / piHat[clone$rowId], 1)
    } else {
      changes <- cumsum(startTime %in% c(0, eventTimes[eventTimes <= graceTime], graceTime))
      weight <- exp(slope[clone$rowId] * (changes - 1)) * (1 + 0.2 * (startTime >= graceTime))
    }
    last <- seq_along(startTime) == length(startTime)
    segments[[i]] <- data.frame(
      cloneRowId = clone$cloneRowId,
      rowId = clone$rowId,
      personSeqId = clone$personSeqId,
      strategyId = clone$strategyId,
      startTime = as.numeric(startTime),
      stopTime = as.numeric(stopTime),
      weight = weight,
      stabilizationNumerator = 1,
      outcome = as.integer(last & clone$endReason == "outcome"),
      competing = as.integer(last & clone$endReason == "death")
    )
  }
  segments <- do.call(rbind, segments)
  segments <- cbind(segmentId = seq_len(nrow(segments)), segments)
  personWeights <- m4PersonWeights(segments)
  strategies <- list(
    createTreatmentStrategy(1, sprintf("Initiate within %d days", grace), "initiateWithin", graceEnd = grace),
    createTreatmentStrategy(2, sprintf("Do not initiate within %d days", grace), "notInitiateWithin", graceEnd = grace)
  )
  metaData <- list(
    strategies = lapply(strategies, function(strategy) strategy$toList()),
    outcomeId = 3,
    treatmentIds = c(treatment = 21),
    sameDayOrdering = "eventsFirst",
    competingEvents = competingEvents,
    attrition = dplyr::tibble(
      description = "Original cohorts",
      targetPersons = as.integer(n),
      comparatorPersons = 0L,
      targetExposures = as.integer(n),
      comparatorExposures = 0L
    ),
    targetId = 1
  )
  return(list(persons = persons, clones = clones, segments = segments, personWeights = personWeights, metaData = metaData))
}

# The personWeights table of segments sorted by cloneRowId and startTime.
m4PersonWeights <- function(segments) {
  last <- !duplicated(segments$cloneRowId, fromLast = TRUE)
  return(data.frame(
    cloneRowId = segments$cloneRowId[last],
    strategyId = segments$strategyId[last],
    maxWeight = as.vector(tapply(segments$weight, segments$cloneRowId, max)),
    weightAtEnd = segments$weight[last]
  ))
}

# Splits segments (sorted by cloneRowId and startTime) at every multiple of `intervalLength` on the t
# scale, keeping their weights, as computeCensoringWeights(outcomeTimeBasis = ) does (Phase1Contracts
# section 9, amendment 5). The outcome and competing-event flags stay on the last piece.
m4SplitOnGrid <- function(segments, intervalLength) {
  # Grid points m L strictly inside (startTime, stopTime]:
  first <- floor(segments$startTime / intervalLength)
  nPieces <- pmax(0, ceiling(segments$stopTime / intervalLength) - 1 - first) + 1
  index <- rep(seq_len(nrow(segments)), nPieces)
  piece <- sequence(nPieces) - 1
  last <- piece == nPieces[index] - 1
  pieces <- segments[index, ]
  pieces$startTime <- ifelse(piece == 0, segments$startTime[index], (first[index] + piece) * intervalLength)
  pieces$stopTime <- ifelse(last, segments$stopTime[index], (first[index] + piece + 1) * intervalLength)
  pieces$outcome <- as.integer(pieces$outcome * last)
  pieces$competing <- as.integer(pieces$competing * last)
  pieces$segmentId <- seq_len(nrow(pieces))
  rownames(pieces) <- NULL
  return(pieces)
}

# Builds the weighted CcwClones object of a fixture (the caller closes it). With `cutIntervalLength`,
# the segments are split on that interval grid and the weights metaData records it, as for clones
# weighted with computeCensoringWeights(outcomeTimeBasis = createTimeBasis(intervalLength = )).
createM4WeightedClones <- function(fixture, cutIntervalLength = NA_integer_) {
  ccwClones <- .createCcwClones(fixture$persons, fixture$clones, fixture$metaData)
  weightsMetaData <- list(stabilized = FALSE, truncation = createWeightTruncation("none")$toList(), cutIntervalLength = as.integer(cutIntervalLength))
  segments <- fixture$segments
  if (!is.na(cutIntervalLength)) {
    segments <- m4SplitOnGrid(segments, cutIntervalLength)
  }
  return(.addSegments(ccwClones, segments, fixture$personWeights, weightsMetaData))
}

# Evaluates `expr` without the warnings of a pooled logistic outcome model with a grace period (a smooth
# time basis, a grid with no interval ending at the end of the grace period (Phase1Contracts section
# 11.6), or clones censored at starts inside the intervals), for tests of other properties of such models. test-integration-warnings.R tests the warnings.
withoutGracePeriodWarnings <- function(expr) {
  return(withCallingHandlers(expr, warning = function(w) {
    if (grepl("^(No interval of the pooled logistic outcome model's grid|The pooled logistic outcome model's smooth time basis|Clones are censored inside the intervals of the pooled logistic outcome model's grid)", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  }))
}

# Fits the outcome model on a fixture (closing the Andromeda object). The pooled logistic model gets
# segments split on its interval grid.
fitM4 <- function(fixture, ...) {
  args <- createFitOutcomeModelArgs(..., control = m4Control)
  cutIntervalLength <- if (args$modelType == "logistic") args$timeBasis$intervalLength else NA_integer_
  weightedClones <- createM4WeightedClones(fixture, cutIntervalLength)
  on.exit(close(weightedClones))
  return(withoutGracePeriodWarnings(suppressMessages(fitOutcomeModel(weightedClones, args))))
}

# Newton steps with a tight standardized-gradient criterion (coordinate descent converges slowly on
# raw polynomial bases).
m4Control <- Cyclops::createControl(algorithm = "newton", convergenceType = "kkt", tolerance = 1e-10, noiseLevel = "silent")

getM4Estimate <- function(model, estimand, column = "estimate") {
  return(model$estimates[[column]][model$estimates$estimand == estimand])
}

# Segments with the survival-package columns: target indicator x and the multi-state factor (censoring
# first).
toSurvivalData <- function(segments) {
  segments$x <- as.integer(segments$strategyId == 1)
  segments$state <- factor(
    ifelse(segments$outcome == 1, "outcome", ifelse(segments$competing == 1, "death", "censor")),
    levels = c("censor", "outcome", "death")
  )
  return(segments)
}

# Sums row-level influence (rows of survfit residuals) by person; returns a vector indexed by the sorted
# unique personSeqIds of `allPersons` (0 for persons absent from the fit).
personInfluence <- function(rowInfluence, personSeqId, allPersons) {
  sums <- rowsum(rowInfluence, personSeqId)
  result <- numeric(length(allPersons))
  result[match(as.integer(rownames(sums)), allPersons)] <- sums[, 1]
  return(result)
}

# Clone x interval rows built independently of the package: one row per clone and interval of
# `intervalLength` days entered, with the weight of the segment covering the clone's last day in the
# interval (Phase1Contracts section 9, amendment 6).
expandM4Rows <- function(segments, intervalLength) {
  rows <- list()
  for (cloneRowId in unique(segments$cloneRowId)) {
    clone <- segments[segments$cloneRowId == cloneRowId, ]
    endTime <- max(clone$stopTime)
    nIntervals <- ceiling(endTime / intervalLength)
    k <- seq_len(nIntervals) - 1
    lastDay <- pmin((k + 1) * intervalLength, endTime)
    weight <- vapply(lastDay, function(t) clone$weight[clone$startTime < t & clone$stopTime >= t], numeric(1))
    last <- k == nIntervals - 1
    rows[[length(rows) + 1]] <- data.frame(
      cloneRowId = cloneRowId,
      personSeqId = clone$personSeqId[1],
      x = as.integer(clone$strategyId[1] == 1),
      k = k,
      weight = weight,
      y = as.integer(last & any(clone$outcome == 1)),
      death = as.integer(last & any(clone$competing == 1))
    )
  }
  return(do.call(rbind, rows))
}

# Young et al. (2020) nonparametric discrete-time hazards and cumulative incidence of one strategy, from
# clone x interval rows (written from the formulas in the paper): the outcome hazard among rows with no
# death in the interval, the death hazard among all rows, and
# CIF(K) = sum_{k < K} h_Y,k (1 - h_D,k) prod_{j < k} (1 - h_Y,j) (1 - h_D,j).
youngNonparametricCif <- function(rows, nIntervals) {
  hazardY <- numeric(nIntervals)
  hazardD <- numeric(nIntervals)
  for (k in seq_len(nIntervals) - 1) {
    atRisk <- rows[rows$k == k, ]
    noDeath <- atRisk[atRisk$death == 0, ]
    hazardD[k + 1] <- sum(atRisk$weight * atRisk$death) / sum(atRisk$weight)
    hazardY[k + 1] <- sum(noDeath$weight * noDeath$y) / sum(noDeath$weight)
  }
  return(youngCif(hazardY, hazardD))
}

youngCif <- function(hazardY, hazardD) {
  survival <- cumprod((1 - hazardY) * (1 - hazardD))
  return(cumsum(hazardY * (1 - hazardD) * c(1, utils::head(survival, -1))))
}

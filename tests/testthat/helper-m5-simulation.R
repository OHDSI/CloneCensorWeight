# Fixtures and hand-written references for the M5 (simulation and diagnostics) tests.

# Simulation settings from simulateCloneCensorWeightData(), for computeTrueRisks().
m5Settings <- function(...) {
  data <- simulateCloneCensorWeightData(n = 10, ..., seed = 1)
  settings <- attr(data, "metaData")$simulationSettings
  close(data)
  return(settings)
}

# A tiny DGP with large daily hazards: one binary and one normal covariate, days 0..8.
m5TinySettings <- function(dgp, graceEnd, treatmentEffect = log(0.3)) {
  return(m5Settings(
    nCovariates = 2,
    dgp = dgp,
    graceEnd = graceEnd,
    maxDays = 8,
    binaryPrevalence = 0.4,
    initiationHazard = 0.3,
    outcomeHazard = 0.06,
    deathHazard = 0.03,
    observationEndHazard = 0.01,
    initiationCoefficients = c(0.8, -0.5),
    outcomeCoefficients = c(0.6, 0.4),
    deathCoefficients = c(-0.3, 0.5),
    treatmentEffect = treatmentEffect,
    effectRampDays = 3
  ))
}

# Truth references from the DGP definition --------------------------------------------------------

# Daily outcome hazards on days 0..maxDays for a person with outcome linear predictor `eta` who
# starts on day `a` (Inf: never), written from the DGP definition in ?simulateCloneCensorWeightData.
# On days up to `untreatedThrough` the person is untreated whatever `a` is (the binary-at-G regime).
m5OutcomeHazards <- function(eta, a, settings, untreatedThrough = a) {
  days <- 0:settings$maxDays
  if (settings$dgp == "a") {
    effect <- ifelse(days > settings$graceEnd & a <= settings$graceEnd, settings$treatmentEffect, 0)
  } else {
    s <- days - a
    effect <- ifelse(s >= 1, settings$treatmentEffect * pmin(s, settings$effectRampDays) / settings$effectRampDays, 0)
  }
  effect[days <= untreatedThrough] <- 0
  return(stats::plogis(eta + effect))
}

# Cause-specific cumulative incidence of the outcome on days 0..maxDays, with the outcome before death
# on the same day: CIF(h) = sum_{d <= h} hY(d) prod_{k < d} (1 - hY(k)) (1 - hD).
m5Cif <- function(outcomeHazards, deathHazard) {
  survivalBefore <- c(1, cumprod((1 - outcomeHazards) * (1 - deathHazard)))[seq_along(outcomeHazards)]
  return(cumsum(outcomeHazards * survivalBefore))
}

# The CIF under a strategy for one person (scalar linear predictors), by enumerating start days with
# their probabilities (contract section 5.6).
m5StrategyCif <- function(etaInitiation, etaOutcome, etaDeath, strategy, graceRegime, settings) {
  hA <- stats::plogis(etaInitiation)
  hD <- stats::plogis(etaDeath)
  maxDays <- settings$maxDays
  cif <- function(a, untreatedThrough = a) {
    return(m5Cif(m5OutcomeHazards(etaOutcome, a, settings, untreatedThrough), hD))
  }
  pStart <- function(a) (1 - hA)^a * hA
  if (strategy$type == "neverInitiate") {
    return(cif(Inf))
  }
  g <- strategy$graceEnd
  if (strategy$type == "notInitiateWithin") {
    # Starts after g follow the natural hazard; starts after maxDays do not matter.
    later <- if (g < maxDays) (g + 1):maxDays else integer(0)
    result <- (1 - hA)^length(later) * cif(Inf)
    for (a in later) {
      result <- result + pStart(a - g - 1) * cif(a)
    }
    return(result)
  }
  if (graceRegime == "naturalGrace") {
    result <- (1 - hA)^g * cif(g)
    for (a in seq_len(g) - 1) {
      result <- result + pStart(a) * cif(a)
    }
    return(result)
  }
  # binaryAtG: early starters keep their day; the others get a ~ P(A = a | A <= g, no outcome by g).
  early <- 0:g
  result <- 0
  for (a in early) {
    result <- result + pStart(a) * cif(a)
  }
  q <- sapply(early, function(a) pStart(a) * prod(1 - m5OutcomeHazards(etaOutcome, a, settings)[1:(g + 1)]))
  q <- q / sum(q)
  for (i in seq_along(early)) {
    result <- result + (1 - hA)^(g + 1) * q[i] * cif(early[i], untreatedThrough = g)
  }
  return(result)
}

# Gauss-Hermite nodes and weights for the standard normal density (Golub-Welsch: eigen decomposition of
# the Jacobi matrix of the probabilists' Hermite polynomials).
m5GaussHermite <- function(m = 40) {
  jacobi <- matrix(0, m, m)
  jacobi[cbind(1:(m - 1), 2:m)] <- sqrt(seq_len(m - 1))
  jacobi[cbind(2:m, 1:(m - 1))] <- sqrt(seq_len(m - 1))
  decomposition <- eigen(jacobi, symmetric = TRUE)
  return(list(nodes = decomposition$values, weights = decomposition$vectors[1, ]^2))
}

# The population risk at each day 0..maxDays under a strategy, for settings with one binary and one
# normal covariate: a sum over the binary value and Gauss-Hermite quadrature over the normal one.
m5ExactRisks <- function(settings, strategy, graceRegime = "binaryAtG") {
  stopifnot(settings$nCovariates == 2)
  eta <- function(hazard, coefficients, x1, x2) stats::qlogis(hazard) + coefficients[1] * x1 + coefficients[2] * x2
  quadrature <- m5GaussHermite()
  risks <- numeric(settings$maxDays + 1)
  for (x1 in 0:1) {
    px1 <- ifelse(x1 == 1, settings$binaryPrevalence, 1 - settings$binaryPrevalence)
    for (k in seq_along(quadrature$nodes)) {
      x2 <- quadrature$nodes[k]
      cif <- m5StrategyCif(
        etaInitiation = eta(settings$initiationHazard, settings$initiationCoefficients, x1, x2),
        etaOutcome = eta(settings$outcomeHazard, settings$outcomeCoefficients, x1, x2),
        etaDeath = eta(settings$deathHazard, settings$deathCoefficients, x1, x2),
        strategy = strategy,
        graceRegime = graceRegime,
        settings = settings
      )
      risks <- risks + px1 * quadrature$weights[k] * cif
    }
  }
  return(risks)
}

# Monte Carlo standard error of a risk estimated from n independent persons.
m5MonteCarloSe <- function(risk, n) {
  return(sqrt(risk * (1 - risk) / n))
}

# Clones and weights built by hand (contract sections 5.2 and 5.4) ------------------------------------

# Dense covariate matrix (rows in rowId order 1..n, columns in covariateId order).
m5CovariateMatrix <- function(data) {
  covariates <- dplyr::collect(data$covariates)
  covariateIds <- sort(dplyr::collect(data$covariateRef)$covariateId)
  n <- attr(data, "metaData")$populationSize
  x <- matrix(0, nrow = n, ncol = length(covariateIds), dimnames = list(NULL, covariateIds))
  x[cbind(covariates$rowId, match(covariates$covariateId, covariateIds))] <- covariates$covariateValue
  return(x)
}

# Persons and clones under initiateWithin(graceEnd) (strategy 1) vs notInitiateWithin(graceEnd)
# (strategy 2), "eventsFirst", from a CohortMethod study population, following contract section 5.2.
m5ReferenceClones <- function(data, graceEnd, riskStart, riskWindowEnd = 365) {
  population <- CohortMethod::createStudyPopulation(
    cohortMethodData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
      riskWindowStart = riskStart,
      riskWindowEnd = riskWindowEnd,
      endAnchor = "cohort start"
    )
  )
  treatments <- dplyr::collect(data$treatments)
  treatments <- treatments[treatments$daysToStart >= 0, ]
  firstStart <- tapply(treatments$daysToStart, treatments$rowId, min)
  outcomeDay <- as.integer(population$daysToEvent)
  deathDay <- population$daysToDeath
  deathDay[!is.na(deathDay) & (deathDay < population$riskStart | deathDay > population$riskEnd)] <- NA
  deathDay[!is.na(deathDay) & !is.na(outcomeDay) & deathDay == outcomeDay] <- NA
  persons <- data.frame(
    rowId = population$rowId,
    personSeqId = population$personSeqId,
    riskStart = as.integer(population$riskStart),
    riskEnd = as.integer(population$riskEnd),
    nominalRiskEnd = as.integer(riskWindowEnd),
    outcomeDay = outcomeDay,
    deathDay = as.integer(deathDay),
    treatmentDay = as.integer(firstStart[as.character(population$rowId)]),
    comparatorTreatmentDay = NA_integer_
  )
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within G", "initiateWithin", graceEnd = graceEnd),
    createTreatmentStrategy(2, "Do not initiate within G", "notInitiateWithin", graceEnd = graceEnd)
  )
  a <- persons$treatmentDay
  deviationDays <- list(
    ifelse(is.na(a) | a > graceEnd, graceEnd, NA_integer_),
    ifelse(!is.na(a) & a <= graceEnd, a, NA_integer_)
  )
  clones <- list()
  for (strategyId in 1:2) {
    # Candidate end days in the "eventsFirst" order:
    candidates <- cbind(persons$outcomeDay, persons$deathDay, deviationDays[[strategyId]], persons$riskEnd)
    endDay <- apply(candidates, 1, min, na.rm = TRUE)
    first <- apply(candidates == endDay, 1, function(x) which(!is.na(x) & x)[1])
    endReason <- c("outcome", "death", "deviation", "riskEnd")[first]
    endReason[endReason == "riskEnd" & persons$riskEnd < persons$nominalRiskEnd] <- "observationEnd"
    clones[[strategyId]] <- data.frame(
      rowId = persons$rowId,
      personSeqId = persons$personSeqId,
      strategyId = strategyId,
      riskStart = persons$riskStart,
      deviationDay = as.integer(deviationDays[[strategyId]]),
      endDay = as.integer(endDay),
      endReason = endReason
    )
  }
  clones <- do.call(rbind, clones)
  clones <- clones[clones$endDay >= clones$riskStart, ]
  clones <- clones[order(clones$rowId, clones$strategyId), ]
  clones$cloneRowId <- seq_len(nrow(clones))
  clones$survivalTime <- clones$endDay - clones$riskStart + 1L
  clones$outcomeCount <- as.integer(clones$endReason == "outcome")
  rownames(clones) <- NULL
  attrition <- CohortMethod::getAttritionTable(population)
  attrition <- rbind(attrition, dplyr::tibble(
    description = "Clones with follow-up",
    targetPersons = sum(clones$strategyId == 1),
    comparatorPersons = sum(clones$strategyId == 2),
    targetExposures = sum(clones$strategyId == 1),
    comparatorExposures = sum(clones$strategyId == 2)
  ))
  metaData <- list(
    strategies = lapply(strategies, function(strategy) strategy$toList()),
    outcomeId = 3,
    treatmentIds = c(treatment = 2),
    sameDayOrdering = "eventsFirst",
    competingEvents = "death",
    attrition = attrition,
    targetId = 1
  )
  return(list(persons = persons, clones = clones, metaData = metaData))
}

# Segments split at day graceEnd, with weight `before` on days up to graceEnd and `after` on later
# days (both per clone). Outcome and competing event on the last segment.
m5ReferenceSegments <- function(clones, graceEnd, before, after) {
  splitTime <- graceEnd - clones$riskStart + 1
  split <- splitTime >= 1 & clones$survivalTime > splitTime
  segment <- function(idx, startTime, stopTime, weight, last) {
    return(data.frame(
      cloneRowId = clones$cloneRowId[idx],
      rowId = clones$rowId[idx],
      personSeqId = clones$personSeqId[idx],
      strategyId = clones$strategyId[idx],
      startTime = as.numeric(startTime),
      stopTime = as.numeric(stopTime),
      weight = weight,
      outcome = as.integer(last & clones$endReason[idx] == "outcome"),
      competing = as.integer(last & clones$endReason[idx] == "death")
    ))
  }
  all <- seq_len(nrow(clones))
  segments <- rbind(
    segment(all, 0, ifelse(split, splitTime, clones$survivalTime), ifelse(splitTime < 1, after, before), !split),
    segment(which(split), splitTime[split], clones$survivalTime[split], after[split], TRUE)
  )
  segments <- segments[order(segments$cloneRowId, segments$startTime), ]
  segments <- cbind(segmentId = seq_len(nrow(segments)), segments)
  rownames(segments) <- NULL
  personWeights <- data.frame(
    cloneRowId = clones$cloneRowId,
    strategyId = clones$strategyId,
    maxWeight = as.numeric(tapply(segments$weight, segments$cloneRowId, max)[as.character(clones$cloneRowId)]),
    weightAtEnd = as.numeric(tapply(segments$weight, segments$cloneRowId, function(w) w[length(w)])[as.character(clones$cloneRowId)])
  )
  return(list(segments = segments, personWeights = personWeights))
}

# Weighted clones for simulated data under initiateWithin(graceEnd) vs notInitiateWithin(graceEnd),
# with the true DGP weights after the grace period: 1 / P(A <= G | X) for the initiators and
# 1 / P(A > G | X) for the deferrers. Before the grace period ends, deferrers get the weight
# 1 / sqrt(P(A > G | X)) (any value would do), so that the two segments differ.
m5WeightedClones <- function(data, graceEnd, riskStart) {
  reference <- m5ReferenceClones(data, graceEnd, riskStart)
  clones <- reference$clones
  settings <- attr(data, "metaData")$simulationSettings
  x <- m5CovariateMatrix(data)
  initiationHazard <- stats::plogis(stats::qlogis(settings$initiationHazard) + as.vector(x %*% settings$initiationCoefficients))
  notStarted <- ((1 - initiationHazard)^(graceEnd + 1))[clones$rowId]
  isTarget <- clones$strategyId == 1
  before <- ifelse(isTarget, 1, 1 / sqrt(notStarted))
  after <- ifelse(isTarget, 1 / (1 - notStarted), 1 / notStarted)
  weighted <- m5ReferenceSegments(clones, graceEnd, before, after)
  ccwClones <- .createCcwClones(reference$persons, clones, reference$metaData)
  ccwClones <- .addSegments(ccwClones, weighted$segments, weighted$personWeights, list(stabilized = FALSE, truncation = "none"))
  return(list(ccwClones = ccwClones, clones = clones, after = after, x = x))
}

# Standardized differences of means computed by hand, with CohortMethod's conventions: the unweighted
# SD divides by n; the weighted SD uses normalized weights w and divides by 1 - sum(w^2); the
# denominator is sqrt((sdTarget^2 + sdComparator^2) / 2).
m5HandBalance <- function(x, treatment, weight = rep(1, length(treatment)), weighted = TRUE) {
  moments <- function(values, w) {
    if (weighted) {
      w <- w / sum(w)
      mean <- sum(w * values)
      sd <- sqrt((sum(w * values^2) - mean^2) / (1 - sum(w^2)))
    } else {
      mean <- mean(values)
      sd <- sqrt(mean(values^2) - mean^2)
    }
    return(c(mean = mean, sd = sd))
  }
  result <- lapply(colnames(x), function(covariateId) {
    target <- moments(x[treatment == 1, covariateId], weight[treatment == 1])
    comparator <- moments(x[treatment == 0, covariateId], weight[treatment == 0])
    data.frame(
      covariateId = as.numeric(covariateId),
      meanTarget = target[["mean"]],
      meanComparator = comparator[["mean"]],
      sdTarget = target[["sd"]],
      sdComparator = comparator[["sd"]],
      stdDiff = (target[["mean"]] - comparator[["mean"]]) / sqrt((target[["sd"]]^2 + comparator[["sd"]]^2) / 2)
    )
  })
  return(do.call(rbind, result))
}

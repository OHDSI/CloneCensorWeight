# Small fixtures for the Foundation tests.

# Data frames for a CloneCensorWeightData with `n` persons: ~10% die, ~30% have outcome 3 (some
# before index), and about half start treatment 21 in the first 30 days (more often with covariate
# 1002). Five binary baseline covariates.
createFoundationTables <- function(n = 200, seed = 1) {
  set.seed(seed)
  rowIds <- seq_len(n)
  cohorts <- data.frame(
    rowId = rowIds,
    personSeqId = rowIds,
    personId = as.character(1000 + rowIds),
    cohortStartDate = as.Date("2020-01-01") + sample.int(365, n, replace = TRUE),
    daysFromObsStart = 365L + sample.int(1000, n, replace = TRUE),
    daysToCohortEnd = 0L,
    daysToObsEnd = sample(30:730, n, replace = TRUE),
    daysToDeath = NA_integer_
  )
  dead <- sample(rowIds, n %/% 10)
  cohorts$daysToDeath[dead] <- as.integer(pmin(cohorts$daysToObsEnd[dead], sample(1:400, length(dead), replace = TRUE)))

  covariates <- data.frame(
    rowId = rep(rowIds, 5),
    covariateId = rep(1001:1005, each = n),
    covariateValue = 1
  )
  covariates <- covariates[runif(nrow(covariates)) < 0.3, ]
  covariates$covariateId <- as.numeric(covariates$covariateId)
  covariateRef <- data.frame(
    covariateId = as.numeric(1001:1005),
    covariateName = paste("Covariate", 1:5),
    analysisId = 1,
    conceptId = 0
  )
  analysisRef <- data.frame(
    analysisId = 1,
    analysisName = "Test",
    domainId = "Condition",
    startDay = -365,
    endDay = 0,
    isBinary = "Y",
    missingMeansZero = "Y"
  )

  hasCovariate <- rowIds %in% covariates$rowId[covariates$covariateId == 1002]
  treated <- rowIds[runif(n) < ifelse(hasCovariate, 0.7, 0.35)]
  treatments <- data.frame(
    rowId = treated,
    treatmentId = 21,
    daysToStart = sample(0:30, length(treated), replace = TRUE)
  )
  treatments$daysToEnd <- treatments$daysToStart + 30L

  withOutcome <- sample(rowIds, round(0.3 * n))
  outcomes <- data.frame(
    rowId = withOutcome,
    outcomeId = 3,
    daysToEvent = sample(-50:400, length(withOutcome), replace = TRUE)
  )

  return(list(
    cohorts = cohorts,
    outcomes = outcomes,
    treatments = treatments,
    covariates = covariates,
    covariateRef = covariateRef,
    analysisRef = analysisRef,
    metaData = list(targetId = 1, treatmentIds = c(treatment = 21))
  ))
}

createFoundationData <- function(n = 200, seed = 1) {
  return(do.call(createCloneCensorWeightData, createFoundationTables(n, seed)))
}

# Valid persons, clones, and metaData for three persons under initiateWithin(7) (strategy 1) vs
# notInitiateWithin(7) (strategy 2).
createFoundationClones <- function() {
  persons <- data.frame(
    rowId = 1:3,
    personSeqId = 1:3,
    riskStart = 0L,
    riskEnd = c(100L, 50L, 200L),
    nominalRiskEnd = c(365L, 365L, 365L),
    outcomeDay = c(NA, 20L, NA),
    deathDay = NA_integer_,
    treatmentDay = c(3L, NA, 10L),
    comparatorTreatmentDay = NA_integer_
  )
  endDay <- c(100L, 3L, 7L, 20L, 7L, 200L)
  clones <- data.frame(
    cloneRowId = 1:6,
    rowId = rep(1:3, each = 2),
    personSeqId = rep(1:3, each = 2),
    strategyId = rep(1:2, 3),
    riskStart = 0L,
    deviationDay = c(NA, 3L, 7L, NA, 7L, NA),
    endDay = endDay,
    survivalTime = endDay + 1L,
    endReason = c("observationEnd", "deviation", "deviation", "outcome", "deviation", "observationEnd"),
    outcomeCount = c(0L, 0L, 0L, 1L, 0L, 0L)
  )
  strategies <- list(
    createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
    createTreatmentStrategy(2, "Do not initiate within 7 days", "notInitiateWithin", graceEnd = 7)
  )
  metaData <- list(
    strategies = lapply(strategies, function(strategy) strategy$toList()),
    outcomeId = 3,
    treatmentIds = c(treatment = 21),
    sameDayOrdering = "eventsFirst",
    competingEvents = "death",
    attrition = dplyr::tibble(
      description = "Original cohorts",
      targetPersons = 3L,
      comparatorPersons = 0L,
      targetExposures = 3L,
      comparatorExposures = 0L
    ),
    targetId = 1
  )
  return(list(persons = persons, clones = clones, metaData = metaData))
}

# One segment per clone of createFoundationClones(), with weight 1.
createFoundationSegments <- function(clones) {
  segments <- data.frame(
    segmentId = clones$cloneRowId,
    cloneRowId = clones$cloneRowId,
    rowId = clones$rowId,
    personSeqId = clones$personSeqId,
    strategyId = clones$strategyId,
    startTime = 0,
    stopTime = as.numeric(clones$survivalTime),
    weight = 1,
    outcome = clones$outcomeCount,
    competing = 0L
  )
  personWeights <- data.frame(
    cloneRowId = clones$cloneRowId,
    strategyId = clones$strategyId,
    maxWeight = 1,
    weightAtEnd = 1
  )
  return(list(segments = segments, personWeights = personWeights))
}

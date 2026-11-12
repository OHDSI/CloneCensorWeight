# Fixtures and a slow pure-R reference implementation of the clone rules (module M2, clones).

# Hand-built 12-person fixture -----------------------------------------------------------------------
#
# Treatment 21 has role "treatment", 22 "comparatorTreatment". Outcome 3 is the study outcome. With
# riskWindowEnd = 30 (anchored at cohort start) and G = 7, the persons are:
#
#  1: starts on day 0; observed past day 30 (outcome 4 on day 2 is another outcome)
#  2: starts on day G = 7 (and had an earlier era before index); outcome on day 25
#  3: starts after G, on day 12 (and again on day 40)
#  4: never starts; outcome on day G = 7, the day of death and of observation end
#  5: outcome on day 3 during the grace period, before starting on day 5
#  6: never starts; death on day 4 during the grace period (observation ends at death)
#  7: never starts; observation ends on day 5 during the grace period; comparator starts on day 1
#  8: outcome on the initiation day (day 4)
#  9: never starts; observation ends on day G = 7; comparator starts on day 2
# 10: prior outcome (day -10), removed by createStudyPopulation()
# 11: outcome on day 0 (a prior outcome when riskStart = 1); never starts
# 12: starts on day 3, the day of death and of observation end
createM2FixtureTables <- function() {
  cohorts <- data.frame(
    rowId = 1:12,
    personSeqId = 101:112,
    personId = as.character(1001:1012),
    cohortStartDate = as.Date("2021-01-01") + 0:11,
    daysFromObsStart = 500L,
    daysToCohortEnd = 0L,
    daysToObsEnd = c(100L, 100L, 100L, 7L, 100L, 4L, 5L, 100L, 7L, 100L, 100L, 3L),
    daysToDeath = c(NA, NA, NA, 7L, NA, 4L, NA, NA, NA, NA, NA, 3L)
  )
  outcomes <- data.frame(
    rowId = c(1L, 2L, 4L, 5L, 5L, 8L, 10L, 11L),
    outcomeId = c(4, 3, 3, 3, 3, 3, 3, 3),
    daysToEvent = c(2L, 25L, 7L, 3L, 25L, 4L, -10L, 0L)
  )
  treatments <- data.frame(
    rowId = c(1L, 2L, 2L, 3L, 3L, 5L, 7L, 8L, 9L, 10L, 12L),
    treatmentId = c(21, 21, 21, 21, 21, 21, 22, 21, 22, 21, 21),
    daysToStart = c(0L, -30L, 7L, 12L, 40L, 5L, 1L, 4L, 2L, 2L, 3L)
  )
  treatments$daysToEnd <- treatments$daysToStart + 10L
  return(list(
    cohorts = cohorts,
    outcomes = outcomes,
    treatments = treatments,
    metaData = list(targetId = 1, treatmentIds = c(treatment = 21, comparatorTreatment = 22))
  ))
}

createM2FixtureData <- function() {
  return(do.call(createCloneCensorWeightData, createM2FixtureTables()))
}

createM2FixturePopulation <- function(data, riskStart) {
  return(suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
      riskWindowStart = riskStart,
      riskWindowEnd = 30,
      endAnchor = "cohort start"
    )
  )))
}

m2Strategy <- function(strategyId, type, graceEnd = NULL, treatment = "treatment") {
  return(createTreatmentStrategy(
    strategyId = strategyId,
    label = paste(type, graceEnd, treatment),
    type = type,
    treatment = treatment,
    graceEnd = graceEnd
  ))
}

# Runs createClones() and returns the persons and clones tables as data frames ordered by row, and the
# metaData.
runM2Clones <- function(population, data, strategies, sameDayOrdering = "eventsFirst", competingEvents = "death") {
  clones <- suppressMessages(createClones(
    population = population,
    cloneCensorWeightData = data,
    createClonesArgs = createCreateClonesArgs(
      strategies = strategies,
      sameDayOrdering = sameDayOrdering,
      competingEvents = competingEvents
    )
  ))
  on.exit(close(clones))
  persons <- as.data.frame(dplyr::arrange(dplyr::collect(clones$persons), .data$rowId))
  cloneTable <- as.data.frame(dplyr::arrange(dplyr::collect(clones$clones), .data$cloneRowId))
  return(list(persons = persons, clones = cloneTable, metaData = attr(clones, "metaData")))
}

# Completes a hand-written expected clones table (rowId, strategyId, deviationDay, endDay, endReason)
# with the columns that follow from the definitions: cloneRowId (ordered by rowId and strategyId),
# personSeqId (100 + rowId in the fixture), riskStart, survivalTime, and outcomeCount.
completeM2ExpectedClones <- function(expected, riskStart) {
  expected <- as.data.frame(expected)
  expected <- expected[order(expected$rowId, expected$strategyId), ]
  result <- data.frame(
    cloneRowId = seq_len(nrow(expected)),
    rowId = as.integer(expected$rowId),
    personSeqId = as.integer(100 + expected$rowId),
    strategyId = as.integer(expected$strategyId),
    riskStart = as.integer(riskStart),
    deviationDay = as.integer(expected$deviationDay),
    endDay = as.integer(expected$endDay),
    survivalTime = as.integer(expected$endDay - riskStart + 1),
    endReason = expected$endReason,
    outcomeCount = as.integer(expected$endReason == "outcome")
  )
  rownames(result) <- NULL
  return(result)
}

# Replaces rows of a hand-written expected table (matched on rowId, and strategyId if present).
replaceM2Rows <- function(expected, replacements) {
  if (nrow(replacements) == 0) {
    return(expected)
  }
  keyColumns <- intersect(c("rowId", "strategyId"), colnames(expected))
  key <- do.call(paste, as.list(as.data.frame(expected)[keyColumns]))
  replacementKey <- do.call(paste, as.list(as.data.frame(replacements)[keyColumns]))
  stopifnot(all(replacementKey %in% key))
  expected <- expected[!key %in% replacementKey, ]
  return(dplyr::bind_rows(expected, replacements))
}

# Random data -------------------------------------------------------------------------------------

# Data frames for a CloneCensorWeightData with `n` persons on a short day grid, so that outcomes,
# deaths, starts, grace-period ends, and observation ends often fall on the same day.
simulateM2Tables <- function(n = 1000, seed = 1) {
  set.seed(seed)
  rowIds <- seq_len(n)
  daysToObsEnd <- sample(0:40, n, replace = TRUE)
  cohorts <- data.frame(
    rowId = rowIds,
    personSeqId = rowIds + 5000L,
    personId = as.character(rowIds),
    cohortStartDate = as.Date("2020-01-01"),
    daysFromObsStart = 400L,
    daysToCohortEnd = sample(0:15, n, replace = TRUE),
    daysToObsEnd = daysToObsEnd,
    daysToDeath = NA_integer_
  )
  dead <- runif(n) < 0.3
  atObsEnd <- runif(n) < 0.5
  cohorts$daysToDeath[dead] <- ifelse(
    atObsEnd[dead],
    daysToObsEnd[dead],
    vapply(daysToObsEnd[dead], function(x) x + 1L - sample.int(x + 1L, 1), integer(1))
  )

  nOutcomes <- sample(0:3, n, replace = TRUE, prob = c(0.4, 0.3, 0.2, 0.1))
  outcomes <- data.frame(
    rowId = rep(rowIds, nOutcomes),
    outcomeId = sample(c(3, 3, 3, 4), sum(nOutcomes), replace = TRUE),
    daysToEvent = sample(-5:40, sum(nOutcomes), replace = TRUE)
  )
  # Some outcomes on the day of death.
  onDeath <- which(dead & runif(n) < 0.3)
  outcomes <- rbind(outcomes, data.frame(rowId = onDeath, outcomeId = 3, daysToEvent = cohorts$daysToDeath[onDeath]))

  nTreatments <- sample(0:2, n, replace = TRUE)
  nComparators <- sample(0:2, n, replace = TRUE, prob = c(0.6, 0.3, 0.1))
  treatments <- data.frame(
    rowId = c(rep(rowIds, nTreatments), rep(rowIds, nComparators)),
    treatmentId = c(rep(21, sum(nTreatments)), rep(22, sum(nComparators))),
    daysToStart = sample(-10:30, sum(nTreatments) + sum(nComparators), replace = TRUE)
  )
  treatments$daysToEnd <- treatments$daysToStart + sample(0:20, nrow(treatments), replace = TRUE)
  return(list(
    cohorts = cohorts,
    outcomes = outcomes,
    treatments = treatments,
    metaData = list(targetId = 1, treatmentIds = c(treatment = 21, comparatorTreatment = 22))
  ))
}

# Slow pure-R reference -------------------------------------------------------------------------------

# Builds the persons and clones tables from the raw data frames (`tables`, as from
# simulateM2Tables()), the study population (only its rows and risk windows are used), and the
# study-population and clone arguments. Each clone is followed by walking its days one at a time
# from the index day: on each day, the events of that day are checked in the same-day order, and the
# first that occurs ends follow-up. A clone that deviates before riskStart has no follow-up.
referenceM2Clones <- function(tables, population, createStudyPopulationArgs, createClonesArgs, outcomeId = 3) {
  treatmentIds <- tables$metaData$treatmentIds
  strategies <- createClonesArgs$strategies
  if (createClonesArgs$sameDayOrdering == "eventsFirst") {
    order <- c("outcome", "death", "deviation")
  } else {
    order <- c("deviation", "outcome", "death")
  }
  population <- population[order(population$rowId), ]
  persons <- vector("list", nrow(population))
  clones <- list()
  for (i in seq_len(nrow(population))) {
    rowId <- population$rowId[i]
    cohort <- tables$cohorts[tables$cohorts$rowId == rowId, ]
    riskStart <- as.integer(population$riskStart[i])
    riskEnd <- as.integer(population$riskEnd[i])

    nominalRiskEnd <- createStudyPopulationArgs$riskWindowEnd
    if (createStudyPopulationArgs$endAnchor == "cohort end") {
      nominalRiskEnd <- nominalRiskEnd + cohort$daysToCohortEnd
    }
    nominalRiskEnd <- min(nominalRiskEnd, riskStart + createStudyPopulationArgs$maxDaysAtRisk)

    outcomeDays <- tables$outcomes$daysToEvent[tables$outcomes$rowId == rowId &
      tables$outcomes$outcomeId == outcomeId &
      tables$outcomes$daysToEvent >= riskStart &
      tables$outcomes$daysToEvent <= riskEnd]
    outcomeDay <- if (length(outcomeDays) > 0) min(outcomeDays) else NA_integer_

    deathDay <- NA_integer_
    death <- cohort$daysToDeath
    if (createClonesArgs$competingEvents == "death" && !is.na(death) && death >= riskStart && death <= riskEnd) {
      if (is.na(outcomeDay) || death != outcomeDay) {
        deathDay <- death
      }
    }

    firstStart <- function(role) {
      if (!role %in% names(treatmentIds)) {
        return(NA_integer_)
      }
      starts <- tables$treatments$daysToStart[tables$treatments$rowId == rowId &
        tables$treatments$treatmentId == treatmentIds[[role]] &
        tables$treatments$daysToStart >= 0]
      if (length(starts) > 0) min(starts) else NA_integer_
    }
    starts <- list(treatment = firstStart("treatment"), comparatorTreatment = firstStart("comparatorTreatment"))

    persons[[i]] <- data.frame(
      rowId = as.integer(rowId),
      personSeqId = as.integer(population$personSeqId[i]),
      riskStart = riskStart,
      riskEnd = riskEnd,
      nominalRiskEnd = as.integer(nominalRiskEnd),
      outcomeDay = as.integer(outcomeDay),
      deathDay = as.integer(deathDay),
      treatmentDay = as.integer(starts$treatment),
      comparatorTreatmentDay = as.integer(starts$comparatorTreatment)
    )

    for (strategyId in 1:2) {
      strategy <- strategies[[strategyId]]
      a <- starts[[strategy$treatment]]
      violatedOn <- function(day) {
        switch(strategy$type,
          initiateWithin = day == strategy$graceEnd && (is.na(a) || a > strategy$graceEnd),
          notInitiateWithin = !is.na(a) && a == day && day <= strategy$graceEnd,
          neverInitiate = !is.na(a) && a == day
        )
      }
      deviationDay <- NA_integer_
      endDay <- NA_integer_
      endReason <- NA_character_
      days <- if (riskEnd >= 0) 0:riskEnd else integer(0)
      for (day in days) {
        violated <- violatedOn(day)
        if (violated && is.na(deviationDay)) {
          deviationDay <- day
        }
        if (!is.na(endDay)) {
          next
        }
        if (day < riskStart) {
          if (violated) {
            endDay <- day
            endReason <- "deviation"
          }
          next
        }
        happens <- c(
          outcome = !is.na(outcomeDay) && outcomeDay == day,
          death = !is.na(deathDay) && deathDay == day,
          deviation = violated
        )
        for (event in order) {
          if (happens[[event]]) {
            endDay <- day
            endReason <- event
            break
          }
        }
        if (is.na(endDay) && day == riskEnd) {
          endDay <- day
          endReason <- if (riskEnd < nominalRiskEnd) "observationEnd" else "riskEnd"
        }
      }
      if (is.na(endDay)) {
        # riskEnd < riskStart: no day of follow-up.
        endDay <- riskEnd
      }
      if (endDay >= riskStart) {
        clones[[length(clones) + 1]] <- data.frame(
          rowId = as.integer(rowId),
          personSeqId = as.integer(population$personSeqId[i]),
          strategyId = as.integer(strategyId),
          riskStart = riskStart,
          deviationDay = as.integer(deviationDay),
          endDay = as.integer(endDay),
          survivalTime = as.integer(endDay - riskStart + 1),
          endReason = endReason,
          outcomeCount = as.integer(endReason == "outcome")
        )
      }
    }
  }
  persons <- do.call(rbind, persons)
  clones <- do.call(rbind, clones)
  clones <- cbind(cloneRowId = seq_len(nrow(clones)), clones)
  rownames(persons) <- NULL
  rownames(clones) <- NULL
  return(list(persons = persons, clones = clones))
}

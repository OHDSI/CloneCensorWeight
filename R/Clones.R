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

#' Create clones
#'
#' @details
#' Clones every person in the study population once per treatment strategy, and follows each clone
#' from the start of the risk window (`riskStart`) until the first of: the outcome, death (when
#' `competingEvents = "death"`), deviation from the strategy, the end of observation, or the end of the
#' risk window. All days are relative to the index date (day 0).
#'
#' Let `A` be the first start (on or after day 0) of the strategy's treatment, `G` its `graceEnd` and
#' `S` its `graceStart`. The strategy is violated (`deviationDay`) on:
#'
#' - `"initiateWithin"`: day `A` if `A < S` (a start before the grace period, only when `S > 0`), and
#'   day `G` if the treatment does not start on or before day `G`.
#' - `"notInitiateWithin"`: day `A` if `A <= G`.
#' - `"neverInitiate"`: day `A`.
#'
#' Observation ends at death: `riskEnd` is the end of the risk window of the study population, capped
#' at the day of death. A violation after `riskEnd` is not recorded (`deviationDay` is `NA`).
#' Follow-up ends at the first of `deviationDay`, the first outcome, death (when
#' `competingEvents = "death"`), and `riskEnd`. When several of these fall on the same day,
#' `sameDayOrdering` decides the end reason (see [createCreateClonesArgs()]). With
#' `"eventsFirst"`, the clone is still under follow-up on the day it deviates, so `endDay` is that
#' day. With `"deviationFirst"`, the deviation comes before everything else on that day, so a clone
#' that ends because of a deviation is last under follow-up on the day before (`endDay =
#' deviationDay - 1`) and is not in the risk set of that day's outcomes.
#'
#' A clone that ends at `riskEnd` ends because of the end of observation (`"observationEnd"`) if
#' `riskEnd` is the end of observation (the earlier of `daysToObsEnd` and `daysToDeath`) and is before
#' the nominal end of the risk window (`nominalRiskEnd`). Otherwise it ends because of the end of the
#' risk window (`"riskEnd"`), which includes a risk window shortened at the start of the next risk
#' window of the same person (`censorAtNewRiskWindow`). With `competingEvents = "none"`, death is not
#' an event: it ends follow-up as the end of observation.
#'
#' Clones that end before `riskStart` (for example a start on day 0 in a `"neverInitiate"` strategy
#' when `riskStart = 1`) have no follow-up and are dropped. In the attrition, the rows of the study
#' population count every person in both the target and the comparator columns (every person is
#' cloned into both strategies), and two rows are added, `"Clones created"` and
#' `"Clones with follow-up"`, with the target strategy counted in the target columns and the comparator
#' strategy in the comparator columns, so that [drawAttritionDiagram()] shows the clones without
#' follow-up as removed from each strategy.
#'
#' This version supports strategies without `prohibitedTreatments`.
#'
#' @param population            A study population as created using [createStudyPopulation()].
#' @param cloneCensorWeightData An object of type [CloneCensorWeightData], the one used to create the
#'                              `population`.
#' @param createClonesArgs      An object of type `CreateClonesArgs` as created using
#'                              [createCreateClonesArgs()].
#'
#' @return
#' An object of type [CcwClones] with the tables
#'
#' - `persons`: one row per row of the population, with `rowId`, `personSeqId`, `riskStart`,
#'   `riskEnd` (capped at death), `nominalRiskEnd`, `outcomeDay` (first outcome in the risk window),
#'   `deathDay` (death in the risk window, `NA` when `competingEvents = "none"` or when it falls on the
#'   outcome day), `treatmentDay`, and `comparatorTreatmentDay` (first start on or after day 0 of each
#'   role).
#' - `clones`: one row per person and strategy with follow-up, with `cloneRowId`, `rowId`,
#'   `personSeqId`, `strategyId` (1 = target, 2 = comparator), `riskStart`, `deviationDay`, `endDay`,
#'   `survivalTime` (`endDay - riskStart + 1`), `endReason` (`"outcome"`, `"death"`, `"deviation"`,
#'   `"observationEnd"`, or `"riskEnd"`), and `outcomeCount`.
#'
#' The `metaData` attribute holds the strategies, outcome ID, treatment IDs, same-day ordering,
#' competing events, attrition, and target ID.
#'
#' @examples
#' data <- createCloneCensorWeightData(
#'   cohorts = data.frame(
#'     rowId = 1:3,
#'     personSeqId = 1:3,
#'     personId = c("101", "102", "103"),
#'     cohortStartDate = as.Date("2020-01-01"),
#'     daysFromObsStart = 400,
#'     daysToCohortEnd = 0,
#'     daysToObsEnd = c(365, 200, 100),
#'     daysToDeath = c(NA, 150, NA)
#'   ),
#'   outcomes = data.frame(rowId = 1, outcomeId = 3, daysToEvent = 50),
#'   treatments = data.frame(rowId = 2, treatmentId = 21, daysToStart = 5, daysToEnd = 35),
#'   metaData = list(targetId = 1, treatmentIds = c(treatment = 21))
#' )
#' studyPop <- createStudyPopulation(
#'   cloneCensorWeightData = data,
#'   outcomeId = 3,
#'   createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
#'     riskWindowEnd = 180,
#'     endAnchor = "cohort start"
#'   )
#' )
#' clones <- createClones(
#'   population = studyPop,
#'   cloneCensorWeightData = data,
#'   createClonesArgs = createCreateClonesArgs(
#'     strategies = list(
#'       createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
#'       createTreatmentStrategy(2, "Do not initiate within 7 days", "notInitiateWithin", graceEnd = 7)
#'     )
#'   )
#' )
#' clones$clones
#' close(clones)
#' close(data)
#'
#' @export
createClones <- function(population, cloneCensorWeightData, createClonesArgs) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertStudyPopulation(population, add = errorMessages)
  .assertCloneCensorWeightData(cloneCensorWeightData, add = errorMessages)
  checkmate::assertR6(createClonesArgs, "CreateClonesArgs", add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  populationMetaData <- attr(population, "metaData")
  if (is.null(populationMetaData$nominalRiskEnd) || is.null(populationMetaData$outcomeId)) {
    stop("The population must be created using createStudyPopulation() from CloneCensorWeight")
  }
  dataMetaData <- attr(cloneCensorWeightData, "metaData")
  treatmentIds <- dataMetaData$treatmentIds
  strategies <- createClonesArgs$strategies
  for (strategy in strategies) {
    if (length(strategy$prohibitedTreatments) > 0) {
      stop(sprintf("Not yet supported in this version: prohibitedTreatments (strategy '%s')", strategy$label))
    }
    if (!strategy$treatment %in% names(treatmentIds)) {
      stop(sprintf(
        "Strategy '%s' refers to treatment role '%s', which has no cohort ID in metaData$treatmentIds",
        strategy$label, strategy$treatment
      ))
    }
  }

  message(sprintf("Creating clones of %d persons for 2 strategies", nrow(population)))
  work <- Andromeda::andromeda(
    population = as.data.frame(population)[, c(
      "rowId", "personSeqId", "riskStart", "riskEnd", "daysToCohortEnd", "daysToObsEnd", "daysToDeath", "daysToEvent"
    )],
    strategies = data.frame(
      strategyId = 1:2,
      strategyType = vapply(strategies, function(strategy) strategy$type, character(1)),
      role = vapply(strategies, function(strategy) strategy$treatment, character(1)),
      graceStart = vapply(strategies, function(strategy) as.integer(strategy$graceStart), integer(1)),
      graceEnd = vapply(strategies, function(strategy) {
        if (is.null(strategy$graceEnd)) NA_integer_ else as.integer(strategy$graceEnd)
      }, integer(1))
    ),
    roles = data.frame(role = names(treatmentIds), treatmentId = as.numeric(treatmentIds))
  )
  on.exit(close(work))
  treatmentIdValues <- as.numeric(treatmentIds)
  work$firstStarts <- cloneCensorWeightData$treatments |>
    filter(.data$treatmentId %in% treatmentIdValues & .data$daysToStart >= 0) |>
    dplyr::group_by(.data$rowId, .data$treatmentId) |>
    dplyr::summarise(startDay = min(.data$daysToStart, na.rm = TRUE), .groups = "drop")

  .andromedaExecute(work, .createPersonsSql(populationMetaData$nominalRiskEnd, createClonesArgs$competingEvents))
  .andromedaExecute(work, .createClonesSql(createClonesArgs$sameDayOrdering))
  .andromedaExecute(work, "
    CREATE TABLE persons AS
    SELECT rowId,
      personSeqId,
      riskStart,
      riskEnd,
      nominalRiskEnd,
      outcomeDay,
      deathDay,
      treatmentDay,
      comparatorTreatmentDay
    FROM workPersons
    ORDER BY rowId;")

  metaData <- list(
    strategies = lapply(strategies, function(strategy) strategy$toList()),
    outcomeId = populationMetaData$outcomeId,
    treatmentIds = treatmentIds,
    sameDayOrdering = createClonesArgs$sameDayOrdering,
    competingEvents = createClonesArgs$competingEvents,
    attrition = dplyr::bind_rows(.cloneAttritionRows(populationMetaData$attrition), .getCloneCounts(work)),
    targetId = dataMetaData$targetId
  )
  return(.createCcwClones(work$persons, work$clones, metaData))
}

# SQL creating the `workPersons` table from the `population`, `roles`, and `firstStarts` tables: the
# columns of the `persons` table plus `observationEnd`, the earlier of the end of observation and death.
# Observation ends at death, so `riskEnd` is capped at death and outcomes after death are ignored.
.createPersonsSql <- function(nominalRiskEnd, competingEvents) {
  nominalRiskEndSql <- sprintf("%d", as.integer(nominalRiskEnd$riskWindowEnd))
  if (grepl("end$", nominalRiskEnd$endAnchor, ignore.case = TRUE)) {
    nominalRiskEndSql <- paste(nominalRiskEndSql, "+ pop.daysToCohortEnd")
  }
  if (!is.null(nominalRiskEnd$maxDaysAtRisk)) {
    nominalRiskEndSql <- sprintf(
      "LEAST(%s, pop.riskStart + %d)", nominalRiskEndSql, as.integer(nominalRiskEnd$maxDaysAtRisk)
    )
  }
  sql <- "
    CREATE TABLE workPersons AS
    WITH typed AS (
      SELECT CAST(rowId AS INTEGER) AS rowId,
        CAST(personSeqId AS INTEGER) AS personSeqId,
        CAST(riskStart AS INTEGER) AS riskStart,
        CAST(riskEnd AS INTEGER) AS riskEnd,
        CAST(daysToCohortEnd AS INTEGER) AS daysToCohortEnd,
        CAST(daysToObsEnd AS INTEGER) AS daysToObsEnd,
        CAST(daysToDeath AS INTEGER) AS daysToDeath,
        CAST(daysToEvent AS INTEGER) AS daysToEvent
      FROM population
    ),
    capped AS (
      SELECT rowId,
        personSeqId,
        riskStart,
        LEAST(riskEnd, daysToDeath) AS riskEnd,
        LEAST(daysToObsEnd, daysToDeath) AS observationEnd,
        daysToCohortEnd,
        daysToDeath,
        daysToEvent
      FROM typed
    ),
    pop AS (
      SELECT *,
        CASE WHEN daysToEvent <= riskEnd THEN daysToEvent END AS outcomeDay
      FROM capped
    ),
    roleStarts AS (
      SELECT firstStarts.rowId,
        roles.role,
        CAST(MIN(firstStarts.startDay) AS INTEGER) AS startDay
      FROM firstStarts
      INNER JOIN roles
        ON firstStarts.treatmentId = roles.treatmentId
      GROUP BY firstStarts.rowId,
        roles.role
    )
    SELECT pop.rowId,
      pop.personSeqId,
      pop.riskStart,
      pop.riskEnd,
      CAST(%s AS INTEGER) AS nominalRiskEnd,
      pop.outcomeDay,
      CASE
        WHEN %s
          AND pop.daysToDeath >= pop.riskStart
          AND pop.daysToDeath <= pop.riskEnd
          AND (pop.outcomeDay IS NULL OR pop.daysToDeath <> pop.outcomeDay)
        THEN pop.daysToDeath
      END AS deathDay,
      treatmentStarts.startDay AS treatmentDay,
      comparatorStarts.startDay AS comparatorTreatmentDay,
      pop.observationEnd
    FROM pop
    LEFT JOIN roleStarts treatmentStarts
      ON pop.rowId = treatmentStarts.rowId
        AND treatmentStarts.role = 'treatment'
    LEFT JOIN roleStarts comparatorStarts
      ON pop.rowId = comparatorStarts.rowId
        AND comparatorStarts.role = 'comparatorTreatment'
    ORDER BY pop.rowId;"
  return(sprintf(sql, nominalRiskEndSql, if (competingEvents == "death") "TRUE" else "FALSE"))
}

# SQL creating the `clones` table from the `workPersons` and `strategies` tables. `exitDay` is the day
# follow-up stops. With "deviationFirst", a clone that stops because of a deviation leaves before that
# day's events, so its last day of follow-up (`endDay`) is the day before.
.createClonesSql <- function(sameDayOrdering) {
  eventOrder <- list(
    outcome = "WHEN outcomeDay = exitDay THEN 'outcome'",
    death = "WHEN deathDay = exitDay THEN 'death'",
    deviation = "WHEN deviationDay = exitDay THEN 'deviation'"
  )
  if (sameDayOrdering == "eventsFirst") {
    eventOrder <- eventOrder[c("outcome", "death", "deviation")]
    deviationOffset <- 0L
  } else {
    eventOrder <- eventOrder[c("deviation", "outcome", "death")]
    deviationOffset <- 1L
  }
  sql <- "
    CREATE TABLE clones AS
    WITH cloned AS (
      SELECT workPersons.rowId,
        workPersons.personSeqId,
        workPersons.riskStart,
        workPersons.riskEnd,
        workPersons.nominalRiskEnd,
        workPersons.observationEnd,
        workPersons.outcomeDay,
        workPersons.deathDay,
        strategies.strategyId,
        strategies.strategyType,
        strategies.graceStart,
        strategies.graceEnd,
        CASE
          WHEN strategies.role = 'treatment' THEN workPersons.treatmentDay
          ELSE workPersons.comparatorTreatmentDay
        END AS startDay
      FROM workPersons
      CROSS JOIN strategies
    ),
    violations AS (
      SELECT *,
        CASE
          WHEN strategyType = 'initiateWithin' AND startDay < graceStart THEN startDay
          WHEN strategyType = 'initiateWithin' AND (startDay IS NULL OR startDay > graceEnd) THEN graceEnd
          WHEN strategyType = 'notInitiateWithin' AND startDay <= graceEnd THEN startDay
          WHEN strategyType = 'neverInitiate' THEN startDay
        END AS violationDay
      FROM cloned
    ),
    deviations AS (
      SELECT *,
        CASE WHEN violationDay <= riskEnd THEN violationDay END AS deviationDay
      FROM violations
    ),
    followUp AS (
      SELECT *,
        LEAST(deviationDay, outcomeDay, deathDay, riskEnd) AS exitDay
      FROM deviations
    ),
    endReasons AS (
      SELECT *,
        CASE
          %s
          WHEN riskEnd = observationEnd AND riskEnd < nominalRiskEnd THEN 'observationEnd'
          ELSE 'riskEnd'
        END AS endReason
      FROM followUp
    ),
    endDays AS (
      SELECT *,
        CASE WHEN endReason = 'deviation' THEN exitDay - %d ELSE exitDay END AS endDay
      FROM endReasons
    )
    SELECT CAST(ROW_NUMBER() OVER (ORDER BY rowId, strategyId) AS INTEGER) AS cloneRowId,
      rowId,
      personSeqId,
      CAST(strategyId AS INTEGER) AS strategyId,
      riskStart,
      CAST(deviationDay AS INTEGER) AS deviationDay,
      CAST(endDay AS INTEGER) AS endDay,
      CAST(endDay - riskStart + 1 AS INTEGER) AS survivalTime,
      endReason,
      CAST(CASE WHEN endReason = 'outcome' THEN 1 ELSE 0 END AS INTEGER) AS outcomeCount
    FROM endDays
    WHERE endDay >= riskStart
    ORDER BY rowId,
      strategyId;"
  return(sprintf(sql, paste(unlist(eventOrder), collapse = "\n          "), deviationOffset))
}

# The population attrition (CohortMethod counts everyone in the target columns) with the target counts
# also in the comparator columns, because every person is cloned into both strategies. The steps of
# CohortMethod::drawAttritionDiagram() (previous minus current row) then never remove a negative number.
.cloneAttritionRows <- function(attrition) {
  attrition$comparatorPersons <- attrition$targetPersons
  attrition$comparatorExposures <- attrition$targetExposures
  return(attrition)
}

# Attrition rows for the clones: clones created, and clones with follow-up (endDay >= riskStart). The
# target strategy is counted in the target columns, the comparator strategy in the comparator columns.
.getCloneCounts <- function(work) {
  created <- .andromedaSql(work, "
    SELECT COUNT(DISTINCT personSeqId) AS persons,
      COUNT(*) AS clones
    FROM persons;")
  withFollowUp <- .andromedaSql(work, "
    SELECT strategyId,
      COUNT(DISTINCT personSeqId) AS persons,
      COUNT(*) AS clones
    FROM clones
    GROUP BY strategyId;")
  idx <- match(1:2, withFollowUp$strategyId)
  persons <- as.integer(ifelse(is.na(idx), 0, withFollowUp$persons[idx]))
  clones <- as.integer(ifelse(is.na(idx), 0, withFollowUp$clones[idx]))
  return(tibble(
    description = c("Clones created", "Clones with follow-up"),
    targetPersons = as.integer(c(created$persons, persons[1])),
    comparatorPersons = as.integer(c(created$persons, persons[2])),
    targetExposures = as.integer(c(created$clones, clones[1])),
    comparatorExposures = as.integer(c(created$clones, clones[2]))
  ))
}

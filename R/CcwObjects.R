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

#' Clones of a clone-censor-weight analysis
#'
#' @description
#' `CcwClones` is an S4 class that inherits from [Andromeda][Andromeda::Andromeda]. It holds the clones
#' created by `createClones()`:
#'
#' - `persons`: one row per row of the study population (uncloned), with the days (since index) of the
#'   risk window, the first outcome, death, and the first start of each treatment.
#' - `clones`: one row per person and strategy with follow-up, with the day the strategy is violated,
#'   the last day of follow-up, and the reason follow-up ended.
#'
#' The weighted clones that `computeCensoringWeights()` returns (a copy in a new Andromeda) also hold
#' `segments` (weighted pieces of each clone's follow-up) and `personWeights` (weight summaries per
#' clone). The `metaData` attribute holds the
#' strategies, outcome ID, treatment IDs, same-day ordering, competing events, attrition, and target ID.
#'
#' @name CcwClones-class
#' @aliases CcwClones
#' @export
setClass("CcwClones", contains = "Andromeda")

.personsColumns <- c(
  rowId = "integer",
  personSeqId = "integer",
  riskStart = "integer",
  riskEnd = "integer",
  nominalRiskEnd = "integer",
  outcomeDay = "integer",
  deathDay = "integer",
  treatmentDay = "integer",
  comparatorTreatmentDay = "integer"
)
.personsNullable <- c("outcomeDay", "deathDay", "treatmentDay", "comparatorTreatmentDay")

.clonesColumns <- c(
  cloneRowId = "integer",
  rowId = "integer",
  personSeqId = "integer",
  strategyId = "integer",
  riskStart = "integer",
  deviationDay = "integer",
  endDay = "integer",
  survivalTime = "integer",
  endReason = "character",
  outcomeCount = "integer"
)
.clonesNullable <- "deviationDay"

.endReasons <- c("outcome", "death", "deviation", "observationEnd", "riskEnd")

.segmentsColumns <- c(
  segmentId = "integer",
  cloneRowId = "integer",
  rowId = "integer",
  personSeqId = "integer",
  strategyId = "integer",
  startTime = "numeric",
  stopTime = "numeric",
  weight = "numeric",
  outcome = "integer",
  competing = "integer"
)

.personWeightsColumns <- c(
  cloneRowId = "integer",
  strategyId = "integer",
  maxWeight = "numeric",
  weightAtEnd = "numeric"
)

.ccwClonesMetaDataFields <- c(
  "strategies", "outcomeId", "treatmentIds", "sameDayOrdering", "competingEvents", "attrition", "targetId"
)

# Runs a query returning one row of violation counts, and stops naming the violated rules.
.checkRules <- function(andromeda, sql, tableName) {
  counts <- .andromedaSql(andromeda, sql)
  violated <- names(counts)[!is.na(unlist(counts)) & unlist(counts) > 0]
  if (length(violated) > 0) {
    stop(sprintf("Table '%s' violates: %s", tableName, paste(violated, collapse = ", ")))
  }
}

# Builds a CcwClones object from `persons` and `clones` (data frames or Andromeda tables), after
# checking their columns, types, and the rules of the clones contract.
.createCcwClones <- function(persons, clones, metaData) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assert(checkmate::checkDataFrame(persons), checkmate::checkClass(persons, "tbl_dbi"), .var.name = "persons", add = errorMessages)
  checkmate::assert(checkmate::checkDataFrame(clones), checkmate::checkClass(clones, "tbl_dbi"), .var.name = "clones", add = errorMessages)
  checkmate::assertList(metaData, add = errorMessages)
  checkmate::assertNames(names(metaData), must.include = .ccwClonesMetaDataFields, .var.name = "names(metaData)", add = errorMessages)
  checkmate::assertList(metaData$strategies, len = 2, .var.name = "metaData$strategies", add = errorMessages)
  checkmate::assertChoice(metaData$sameDayOrdering, c("eventsFirst", "deviationFirst"), .var.name = "metaData$sameDayOrdering", add = errorMessages)
  checkmate::assertChoice(metaData$competingEvents, c("death", "none"), .var.name = "metaData$competingEvents", add = errorMessages)
  checkmate::assertDataFrame(metaData$attrition, .var.name = "metaData$attrition", add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  .checkColumns(persons, .personsColumns, "persons", nullable = .personsNullable)
  .checkColumns(clones, .clonesColumns, "clones", nullable = .clonesNullable)

  ccwClones <- Andromeda::andromeda()
  ccwClones$persons <- persons
  ccwClones$clones <- clones
  tryCatch(
    {
      .checkRules(ccwClones, "
        SELECT COUNT(*) - COUNT(DISTINCT rowId) AS duplicateRowId
        FROM persons;", "persons")
      .checkRules(ccwClones, sprintf("
        SELECT COUNT(*) - COUNT(DISTINCT cloneRowId) AS duplicateCloneRowId,
          CASE WHEN COUNT(*) > 0 AND (MIN(cloneRowId) <> 1 OR MAX(cloneRowId) <> COUNT(*)) THEN 1 ELSE 0 END AS cloneRowIdNotOneToN,
          SUM(CASE WHEN strategyId NOT IN (1, 2) THEN 1 ELSE 0 END) AS strategyIdNotOneOrTwo,
          SUM(CASE WHEN survivalTime <> endDay - riskStart + 1 OR survivalTime < 1 THEN 1 ELSE 0 END) AS badSurvivalTime,
          SUM(CASE WHEN endReason NOT IN (%s) THEN 1 ELSE 0 END) AS badEndReason,
          SUM(CASE WHEN outcomeCount <> CASE WHEN endReason = 'outcome' THEN 1 ELSE 0 END THEN 1 ELSE 0 END) AS badOutcomeCount
        FROM clones;", paste0("'", .endReasons, "'", collapse = ", ")), "clones")
    },
    error = function(e) {
      close(ccwClones)
      stop(conditionMessage(e), call. = FALSE)
    }
  )
  attr(ccwClones, "metaData") <- metaData
  class(ccwClones) <- "CcwClones"
  attr(class(ccwClones), "package") <- "CloneCensorWeight"
  return(ccwClones)
}

# A new CcwClones object, in its own Andromeda, with copies of the `persons` and `clones` tables and of
# the metaData (without `metaData$weights`) of `clones`. computeCensoringWeights() weights such a copy,
# so weighting the same clones twice leaves the earlier weighted object unchanged.
.copyCcwClones <- function(clones) {
  copy <- Andromeda::andromeda()
  copy$persons <- clones$persons
  copy$clones <- clones$clones
  metaData <- attr(clones, "metaData")
  metaData$weights <- NULL
  attr(copy, "metaData") <- metaData
  class(copy) <- "CcwClones"
  attr(class(copy), "package") <- "CloneCensorWeight"
  return(copy)
}

# Adds the `segments` and `personWeights` tables (data frames or Andromeda tables) to a CcwClones
# object, after checking their columns, types and rules, and stores `weightsMetaData` as
# `metaData$weights`. The tables are added in place (the caller's object shares the connection), but
# `metaData$weights` is set only on the returned object, so callers must use the returned object.
# .assertCcwClones(x, weighted = TRUE) rejects the caller's original object.
.addSegments <- function(clones, segments, personWeights, weightsMetaData) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(clones, add = errorMessages)
  checkmate::assert(checkmate::checkDataFrame(segments), checkmate::checkClass(segments, "tbl_dbi"), .var.name = "segments", add = errorMessages)
  checkmate::assert(checkmate::checkDataFrame(personWeights), checkmate::checkClass(personWeights, "tbl_dbi"), .var.name = "personWeights", add = errorMessages)
  checkmate::assertList(weightsMetaData, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  .checkColumns(segments, .segmentsColumns, "segments")
  .checkColumns(personWeights, .personWeightsColumns, "personWeights")

  clones$segments <- segments
  clones$personWeights <- personWeights
  tryCatch(
    .checkRules(clones, "
      SELECT COUNT(*) - COUNT(DISTINCT segmentId) AS duplicateSegmentId,
        SUM(CASE WHEN stopTime <= startTime THEN 1 ELSE 0 END) AS stopTimeNotAfterStartTime,
        SUM(CASE WHEN weight < 0 THEN 1 ELSE 0 END) AS negativeWeight,
        SUM(CASE WHEN outcome NOT IN (0, 1) OR competing NOT IN (0, 1) OR outcome + competing > 1 THEN 1 ELSE 0 END) AS badEventIndicator
      FROM segments;", "segments"),
    error = function(e) {
      clones$segments <- NULL
      clones$personWeights <- NULL
      stop(conditionMessage(e), call. = FALSE)
    }
  )
  metaData <- attr(clones, "metaData")
  metaData$weights <- weightsMetaData
  attr(clones, "metaData") <- metaData
  return(clones)
}

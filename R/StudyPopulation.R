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

#' Create a study population
#'
#' @details
#' Creates the study population of a clone-censor-weight analysis by calling
#' [CohortMethod::createStudyPopulation()] on the eligible cohort: it removes persons with a prior
#' outcome, defines the risk window (`riskStart`, `riskEnd`, in days since index, with `riskEnd`
#' capped at the end of observation), removes persons with too few days at risk, and adds the first
#' outcome in the risk window (`daysToEvent`).
#'
#' The risk window must start relative to the cohort start (`startAnchor = "cohort start"`), because
#' the treatment strategies are defined in days since index, and it cannot start before the index date.
#'
#' The rules that define the risk-window end before it is capped at the end of observation
#' (`riskWindowEnd`, `endAnchor` and `maxDaysAtRisk`) are stored in `metaData$nominalRiskEnd`, so
#' [createClones()] can tell a clone that reached the end of observation (`"observationEnd"`) from one
#' that reached the end of the risk window (`"riskEnd"`). The outcome ID is stored in
#' `metaData$outcomeId`.
#'
#' @param cloneCensorWeightData     An object of type [CloneCensorWeightData] as generated using
#'                                  `getDbCloneCensorWeightData()` or [createCloneCensorWeightData()].
#' @param population                If specified, this population will be used as the starting point
#'                                  instead of the cohorts in the `cloneCensorWeightData` object.
#' @param outcomeId                 The ID of the outcome.
#' @param createStudyPopulationArgs An object of type `CreateStudyPopulationArgs` as created using
#'                                  [CohortMethod::createCreateStudyPopulationArgs()].
#'
#' @return
#' A data frame with one row per person-index, as returned by [CohortMethod::createStudyPopulation()].
#' Columns used downstream: `rowId`, `personSeqId`, `cohortStartDate`, `daysFromObsStart`,
#' `daysToCohortEnd`, `daysToObsEnd`, `daysToDeath`, `riskStart`, `riskEnd`, and `daysToEvent`. The
#' `metaData` attribute holds the attrition, `outcomeId`, and `nominalRiskEnd`.
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
#'   outcomes = data.frame(rowId = c(1, 3), outcomeId = 3, daysToEvent = c(50, -20)),
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
#' studyPop
#' close(data)
#'
#' @export
createStudyPopulation <- function(cloneCensorWeightData,
                                  population = NULL,
                                  outcomeId,
                                  createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs()) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCloneCensorWeightData(cloneCensorWeightData, add = errorMessages)
  checkmate::assertDataFrame(population, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumber(outcomeId, add = errorMessages)
  checkmate::assertR6(createStudyPopulationArgs, "CreateStudyPopulationArgs", add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (outcomeId %% 1 != 0) {
    stop("outcomeId must be a (64-bit) integer")
  }
  .checkRiskWindowStart(createStudyPopulationArgs)

  population <- CohortMethod::createStudyPopulation(
    cohortMethodData = cloneCensorWeightData,
    population = population,
    outcomeId = outcomeId,
    createStudyPopulationArgs = createStudyPopulationArgs
  )
  return(.addStudyPopulationMetaData(population, outcomeId, createStudyPopulationArgs))
}

.checkRiskWindowStart <- function(createStudyPopulationArgs) {
  if (!grepl("start$", createStudyPopulationArgs$startAnchor, ignore.case = TRUE)) {
    stop("The risk window must start relative to the cohort start (startAnchor = 'cohort start')")
  }
  if (createStudyPopulationArgs$riskWindowStart < 0) {
    stop("The risk window cannot start before the index date (riskWindowStart must be >= 0)")
  }
}

.addStudyPopulationMetaData <- function(population, outcomeId, createStudyPopulationArgs) {
  metaData <- attr(population, "metaData")
  metaData$outcomeId <- outcomeId
  metaData$nominalRiskEnd <- list(
    riskWindowEnd = createStudyPopulationArgs$riskWindowEnd,
    endAnchor = createStudyPopulationArgs$endAnchor,
    maxDaysAtRisk = createStudyPopulationArgs$maxDaysAtRisk
  )
  attr(population, "metaData") <- metaData
  return(population)
}

# The study population without an outcome, for censoring models shared across outcomes
# (runCcwAnalyses() with refitCensoringModelsForEveryOutcome = FALSE): the persons and risk windows of
# createStudyPopulation(), without removing persons with a prior outcome, and with no outcomes
# (daysToEvent NA, metaData$outcomeId NA), so clones of it end only at deviation, death, or the end of the
# risk window. Every study population created with the same arguments for an outcome is a subset of it.
.createStudyPopulationWithoutOutcome <- function(cloneCensorWeightData, createStudyPopulationArgs) {
  .checkRiskWindowStart(createStudyPopulationArgs)
  population <- CohortMethod::createStudyPopulation(
    cohortMethodData = cloneCensorWeightData,
    outcomeId = NULL,
    createStudyPopulationArgs = createStudyPopulationArgs
  )
  population$daysToEvent <- rep(NA_integer_, nrow(population))
  return(.addStudyPopulationMetaData(population, NA_real_, createStudyPopulationArgs))
}

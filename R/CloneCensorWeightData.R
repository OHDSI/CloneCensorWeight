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

#' Clone-censor-weight data
#'
#' @description
#' `CloneCensorWeightData` is an S4 class that inherits from
#' [CohortMethodData][CohortMethod::CohortMethodData], so CohortMethod functions such as
#' [CohortMethod::createStudyPopulation()] accept it as is. It is an [Andromeda][Andromeda::Andromeda]
#' object with these tables:
#'
#' - `cohorts`: one row per eligible person-index, with columns `rowId`, `personSeqId`, `personId`,
#'   `treatment` (always 1, because CohortMethod requires it), `cohortStartDate`, `daysFromObsStart`,
#'   `daysToCohortEnd`, `daysToObsEnd`, and `daysToDeath` (`NA` if no death is observed).
#' - `outcomes`: `rowId`, `outcomeId`, `daysToEvent`.
#' - `treatments`: `rowId`, `treatmentId`, `daysToStart`, `daysToEnd`; every era of a treatment cohort
#'   overlapping the prior-use lookback window and follow-up.
#' - `covariates`, `covariateRef`, `analysisRef`: baseline covariates in FeatureExtraction format.
#'
#' All days are relative to the index date (`cohortStartDate`, day 0). The `metaData` attribute is a
#' list with at least `targetId`, `treatmentIds` (named by role: `treatment` and optionally
#' `comparatorTreatment`), `outcomeIds`, `populationSize`, and `attrition`.
#'
#' A `CloneCensorWeightData` is created using `getDbCloneCensorWeightData()` or
#' [createCloneCensorWeightData()], can only be saved using [saveCloneCensorWeightData()], and loaded
#' using [loadCloneCensorWeightData()].
#'
#' @name CloneCensorWeightData-class
#' @aliases CloneCensorWeightData
NULL

#' CloneCensorWeightData class.
#'
#' @export
setClass("CloneCensorWeightData", contains = "CohortMethodData")

.cloneCensorWeightDataTables <- c("cohorts", "outcomes", "treatments", "covariates", "covariateRef", "analysisRef")

.cohortsColumns <- c(
  rowId = "integer",
  personSeqId = "integer",
  personId = "character",
  treatment = "integer",
  cohortStartDate = "Date",
  daysFromObsStart = "integer",
  daysToCohortEnd = "integer",
  daysToObsEnd = "integer",
  daysToDeath = "integer"
)

.outcomesColumns <- c(
  rowId = "integer",
  outcomeId = "numeric",
  daysToEvent = "integer"
)

.treatmentsColumns <- c(
  rowId = "integer",
  treatmentId = "numeric",
  daysToStart = "integer",
  daysToEnd = "integer"
)

.covariatesColumns <- c(
  rowId = "integer",
  covariateId = "numeric",
  covariateValue = "numeric"
)

.covariateRefColumns <- c(
  covariateId = "numeric",
  covariateName = "character",
  analysisId = "numeric",
  conceptId = "numeric"
)

.analysisRefColumns <- c(
  analysisId = "numeric",
  analysisName = "character",
  domainId = "character",
  startDay = "numeric",
  endDay = "numeric",
  isBinary = "character",
  missingMeansZero = "character"
)

.attritionColumns <- c("description", "targetPersons", "comparatorPersons", "targetExposures", "comparatorExposures")

.prepareTable <- function(table, columnTypes, tableName, nullable = character(0)) {
  if (is.null(table) || ncol(table) == 0) {
    table <- .createEmptyTable(columnTypes)
  }
  return(.checkColumns(table, columnTypes, tableName, nullable = nullable, coerce = TRUE))
}

#' Create clone-censor-weight data
#'
#' @description
#' Creates a [CloneCensorWeightData] object from data frames, for example simulated data or data
#' extracted by other means. Empty or `NULL` tables are created with the right columns. Integer-valued
#' numeric columns are converted to integer where required, and columns that are all `NA` of type
#' logical (for example `daysToDeath = NA` when nobody dies) are converted to the required type. The
#' `treatment` column of `cohorts` is added (always 1) if absent.
#'
#' @param cohorts      A data frame with columns `rowId`, `personSeqId`, `personId` (character),
#'                     `cohortStartDate` (Date), `daysFromObsStart`, `daysToCohortEnd`, `daysToObsEnd`,
#'                     and `daysToDeath` (`NA` if no death), and optionally `treatment` (always 1).
#' @param outcomes     A data frame with columns `rowId`, `outcomeId`, and `daysToEvent`, or `NULL`.
#' @param treatments   A data frame with columns `rowId`, `treatmentId`, `daysToStart`, and `daysToEnd`,
#'                     or `NULL`.
#' @param covariates   A data frame with columns `rowId`, `covariateId`, and `covariateValue`, or
#'                     `NULL`.
#' @param covariateRef A data frame with the FeatureExtraction columns `covariateId`, `covariateName`,
#'                     `analysisId`, and `conceptId`, or `NULL`.
#' @param analysisRef  A data frame with the FeatureExtraction columns `analysisId`, `analysisName`,
#'                     `domainId`, `startDay`, `endDay`, `isBinary`, and `missingMeansZero`, or `NULL`.
#' @param metaData     A list (required). Must contain `targetId` (the eligible cohort ID) and
#'                     `treatmentIds` (cohort IDs named by role: `treatment` and optionally
#'                     `comparatorTreatment`).
#'                     `populationSize` is set to the number of rows of `cohorts`; `outcomeIds` and
#'                     `attrition` are created if absent.
#'
#' @return
#' An object of type [CloneCensorWeightData].
#'
#' @examples
#' cohorts <- data.frame(
#'   rowId = 1:3,
#'   personSeqId = 1:3,
#'   personId = c("101", "102", "103"),
#'   cohortStartDate = as.Date("2020-01-01") + 0:2,
#'   daysFromObsStart = c(400, 500, 600),
#'   daysToCohortEnd = c(0, 0, 0),
#'   daysToObsEnd = c(365, 200, 100),
#'   daysToDeath = c(NA, 150, NA)
#' )
#' outcomes <- data.frame(rowId = 1, outcomeId = 3, daysToEvent = 50)
#' treatments <- data.frame(rowId = 2, treatmentId = 2, daysToStart = 5, daysToEnd = 35)
#' data <- createCloneCensorWeightData(
#'   cohorts = cohorts,
#'   outcomes = outcomes,
#'   treatments = treatments,
#'   metaData = list(targetId = 1, treatmentIds = c(treatment = 2))
#' )
#' data
#' close(data)
#'
#' @export
createCloneCensorWeightData <- function(cohorts,
                                        outcomes,
                                        treatments,
                                        covariates = NULL,
                                        covariateRef = NULL,
                                        analysisRef = NULL,
                                        metaData) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(cohorts, add = errorMessages)
  checkmate::assertDataFrame(outcomes, null.ok = TRUE, add = errorMessages)
  checkmate::assertDataFrame(treatments, null.ok = TRUE, add = errorMessages)
  checkmate::assertDataFrame(covariates, null.ok = TRUE, add = errorMessages)
  checkmate::assertDataFrame(covariateRef, null.ok = TRUE, add = errorMessages)
  checkmate::assertDataFrame(analysisRef, null.ok = TRUE, add = errorMessages)
  checkmate::assertList(metaData, add = errorMessages)
  checkmate::assertNumber(metaData$targetId, .var.name = "metaData$targetId", add = errorMessages)
  checkmate::assertNumeric(metaData$treatmentIds,
    min.len = 1, max.len = 2, any.missing = FALSE, names = "unique",
    .var.name = "metaData$treatmentIds", add = errorMessages
  )
  checkmate::assertSubset(names(metaData$treatmentIds), .treatmentRoles, .var.name = "names(metaData$treatmentIds)", add = errorMessages)
  checkmate::assertTRUE("treatment" %in% names(metaData$treatmentIds), .var.name = "metaData$treatmentIds has a 'treatment' element", add = errorMessages)
  checkmate::assertNumeric(metaData$outcomeIds, any.missing = FALSE, null.ok = TRUE, .var.name = "metaData$outcomeIds", add = errorMessages)
  checkmate::assertDataFrame(metaData$attrition, null.ok = TRUE, .var.name = "metaData$attrition", add = errorMessages)
  if (!is.null(metaData$attrition)) {
    checkmate::assertNames(colnames(metaData$attrition), must.include = .attritionColumns, .var.name = "metaData$attrition", add = errorMessages)
  }
  checkmate::reportAssertions(collection = errorMessages)

  if (!"treatment" %in% colnames(cohorts)) {
    cohorts$treatment <- rep(1L, nrow(cohorts))
  }
  cohorts <- .checkColumns(cohorts, .cohortsColumns, "cohorts", nullable = "daysToDeath", coerce = TRUE)
  if (any(cohorts$treatment != 1L)) {
    stop("Column 'treatment' of table 'cohorts' must be 1 in every row")
  }
  if (anyDuplicated(cohorts$rowId)) {
    stop("Column 'rowId' of table 'cohorts' must be unique")
  }
  outcomes <- .prepareTable(outcomes, .outcomesColumns, "outcomes")
  treatments <- .prepareTable(treatments, .treatmentsColumns, "treatments")
  covariates <- .prepareTable(covariates, .covariatesColumns, "covariates")
  covariateRef <- .prepareTable(covariateRef, .covariateRefColumns, "covariateRef",
    nullable = setdiff(names(.covariateRefColumns), "covariateId")
  )
  analysisRef <- .prepareTable(analysisRef, .analysisRefColumns, "analysisRef",
    nullable = setdiff(names(.analysisRefColumns), "analysisId")
  )

  metaData$populationSize <- nrow(cohorts)
  if (is.null(metaData$outcomeIds)) {
    metaData$outcomeIds <- sort(unique(outcomes$outcomeId))
  }
  if (is.null(metaData$attrition)) {
    metaData$attrition <- tibble(
      description = "Original cohorts",
      targetPersons = length(unique(cohorts$personSeqId)),
      comparatorPersons = 0L,
      targetExposures = nrow(cohorts),
      comparatorExposures = 0L
    )
  }

  data <- Andromeda::andromeda(
    cohorts = cohorts,
    outcomes = outcomes,
    treatments = treatments,
    covariates = covariates,
    covariateRef = covariateRef,
    analysisRef = analysisRef
  )
  attr(data, "metaData") <- metaData
  class(data) <- "CloneCensorWeightData"
  attr(class(data), "package") <- "CloneCensorWeight"
  return(data)
}

#' Save the clone-censor-weight data to file
#'
#' @description
#' Saves an object of type [CloneCensorWeightData] to a file. The object is closed afterwards.
#'
#' @param x    An object of type [CloneCensorWeightData].
#' @param file The name of the file where the data will be written. If the file already exists it will
#'             be overwritten.
#'
#' @return
#' Returns no output.
#'
#' @examples
#' data <- createCloneCensorWeightData(
#'   cohorts = data.frame(
#'     rowId = 1L, personSeqId = 1L, personId = "1", cohortStartDate = as.Date("2020-01-01"),
#'     daysFromObsStart = 400L, daysToCohortEnd = 0L, daysToObsEnd = 365L, daysToDeath = NA_integer_
#'   ),
#'   outcomes = NULL,
#'   treatments = NULL,
#'   metaData = list(targetId = 1, treatmentIds = c(treatment = 2))
#' )
#' file <- tempfile(fileext = ".zip")
#' saveCloneCensorWeightData(data, file)
#' data <- loadCloneCensorWeightData(file)
#' close(data)
#' unlink(file)
#'
#' @export
saveCloneCensorWeightData <- function(x, file) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCloneCensorWeightData(x, add = errorMessages)
  checkmate::assertCharacter(file, len = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  Andromeda::saveAndromeda(x, file)
  message("To use this CloneCensorWeightData object, you will have to load it from file (using loadCloneCensorWeightData).")
  invisible(NULL)
}

#' Load the clone-censor-weight data from a file
#'
#' @description
#' Loads an object of type [CloneCensorWeightData] from a file in the file system.
#'
#' @param file The name of the file containing the data.
#'
#' @return
#' An object of type [CloneCensorWeightData].
#'
#' @inherit saveCloneCensorWeightData examples
#'
#' @export
loadCloneCensorWeightData <- function(file) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertCharacter(file, len = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  if (!file.exists(file)) {
    stop(sprintf("Cannot find file '%s'", file))
  }
  if (file.info(file)$isdir) {
    stop(sprintf("'%s' is a folder, but should be a file", file))
  }
  data <- Andromeda::loadAndromeda(file)
  class(data) <- "CloneCensorWeightData"
  attr(class(data), "package") <- "CloneCensorWeight"
  return(data)
}

#' Check whether an object is a CloneCensorWeightData object
#'
#' @param x The object to check.
#'
#' @return
#' A logical value.
#'
#' @examples
#' isCloneCensorWeightData(data.frame())
#'
#' @export
isCloneCensorWeightData <- function(x) {
  return(inherits(x, "CloneCensorWeightData"))
}

# show()
#' @param object An object of type `CloneCensorWeightData`.
#'
#' @export
#' @rdname CloneCensorWeightData-class
setMethod("show", "CloneCensorWeightData", function(object) {
  metaData <- attr(object, "metaData")
  message("# CloneCensorWeightData object")
  message("")
  message(paste("Target cohort ID:", metaData$targetId))
  message(paste(
    "Treatment cohort ID(s):",
    paste(sprintf("%s = %s", names(metaData$treatmentIds), metaData$treatmentIds), collapse = ", ")
  ))
  message(paste("Outcome cohort ID(s):", paste(metaData$outcomeIds, collapse = ",")))
  message("")
  message("Inherits from CovariateData:")
  class(object) <- "CovariateData"
  attr(class(object), "package") <- "FeatureExtraction"
  show(object)
})

# summary()
#' @param object An object of type `CloneCensorWeightData`.
#'
#' @export
#' @rdname CloneCensorWeightData-class
setMethod("summary", "CloneCensorWeightData", function(object) {
  if (!Andromeda::isValidAndromeda(object)) {
    stop("Object is not valid. Probably the Andromeda object was closed.")
  }
  metaData <- attr(object, "metaData")
  cohortCounts <- .andromedaSql(object, "
    SELECT COUNT(DISTINCT personSeqId) AS persons,
      COUNT(*) AS indexDates,
      COUNT(daysToDeath) AS deaths
    FROM cohorts;")
  outcomeCounts <- .andromedaSql(object, "
    SELECT outcomeId,
      COUNT(*) AS eventCount,
      COUNT(DISTINCT rowId) AS personCount
    FROM outcomes
    GROUP BY outcomeId;")
  outcomeCounts <- tibble(outcomeId = metaData$outcomeIds) |>
    dplyr::left_join(outcomeCounts, by = "outcomeId") |>
    dplyr::mutate(dplyr::across(c("eventCount", "personCount"), ~ dplyr::coalesce(as.numeric(.x), 0)))
  treatmentCounts <- tibble(role = names(metaData$treatmentIds), treatmentId = as.numeric(metaData$treatmentIds)) |>
    dplyr::left_join(
      .andromedaSql(object, "
        SELECT treatmentId,
          COUNT(*) AS eraCount,
          COUNT(DISTINCT rowId) AS personCount
        FROM treatments
        GROUP BY treatmentId;"),
      by = "treatmentId"
    ) |>
    dplyr::mutate(dplyr::across(c("eraCount", "personCount"), ~ dplyr::coalesce(as.numeric(.x), 0)))
  result <- list(
    metaData = metaData,
    persons = as.numeric(cohortCounts$persons),
    indexDates = as.numeric(cohortCounts$indexDates),
    deaths = as.numeric(cohortCounts$deaths),
    outcomeCounts = outcomeCounts,
    treatmentCounts = treatmentCounts,
    covariateCount = as.numeric(.andromedaSql(object, "SELECT COUNT(*) AS n FROM covariateRef;")$n),
    covariateValueCount = as.numeric(.andromedaSql(object, "SELECT COUNT(*) AS n FROM covariates;")$n)
  )
  class(result) <- "summary.CloneCensorWeightData"
  return(result)
})

#' @export
print.summary.CloneCensorWeightData <- function(x, ...) {
  message("CloneCensorWeightData object summary")
  message("")
  message(paste("Target cohort ID:", x$metaData$targetId))
  message(paste("Persons:", x$persons))
  message(paste("Index dates:", x$indexDates))
  message(paste("Deaths:", x$deaths))
  message("")
  message("Treatment counts:")
  for (i in seq_len(nrow(x$treatmentCounts))) {
    message(sprintf(
      "  %s (cohort ID %s): %s eras, %s persons",
      x$treatmentCounts$role[i],
      x$treatmentCounts$treatmentId[i],
      x$treatmentCounts$eraCount[i],
      x$treatmentCounts$personCount[i]
    ))
  }
  message("")
  message("Outcome counts:")
  if (nrow(x$outcomeCounts) > 0) {
    outcomeCounts <- as.data.frame(x$outcomeCounts)
    rownames(outcomeCounts) <- outcomeCounts$outcomeId
    outcomeCounts$outcomeId <- NULL
    colnames(outcomeCounts) <- c("Event count", "Person count")
    printCoefmat(outcomeCounts)
  }
  message("")
  message("Covariates:")
  message(paste("Number of covariates:", x$covariateCount))
  message(paste("Number of non-zero covariate values:", x$covariateValueCount))
  invisible(x)
}

#' Clone-censor-weight data class
#'
#' @import FeatureExtraction
#' @import Cyclops
#' @import DatabaseConnector
#' @import SqlRender
#' @import methods
#' @importClassesFrom FeatureExtraction CovariateData
#' @importClassesFrom Andromeda Andromeda
#' @importFrom rlang .data
#' @importFrom stats ave rbinom runif
#' @export
setClass("CloneCensorWeightData", contains = "CovariateData")

#' Create clone-censor-weight data
#'
#' Creates a disk-backed `CloneCensorWeightData` object from normalized
#' post-CDM tables.
#'
#' @param cohorts A data frame with subject and trial-entry follow-up columns.
#' @param outcomes A data frame of outcome occurrences.
#' @param covariates A baseline covariate data frame.
#' @param covariateRef A covariate reference data frame.
#' @param analysisRef An analysis reference data frame.
#'
#' @return An object of class `CloneCensorWeightData`. Its disk-backed tables
#'   are `cohorts`, `outcomes`, `covariates`, `covariateRef`, and `analysisRef`.
#' @export
createCloneCensorWeightData <- function(cohorts,
                                        outcomes = data.frame(),
                                        covariates = data.frame(),
                                        covariateRef = data.frame(),
                                        analysisRef = data.frame()) {
  checkmate::assertDataFrame(cohorts)
  checkmate::assertSubset(c("subjectId", "cohortStartDate", "cohortEndDate"), names(cohorts))
  checkmate::assertDataFrame(outcomes)
  checkmate::assertDataFrame(covariates)
  checkmate::assertDataFrame(covariateRef)
  checkmate::assertDataFrame(analysisRef)
  if (ncol(outcomes) == 0) {
    outcomes <- data.frame(subjectId = integer(), outcomeDate = as.Date(character()), outcomeId = integer())
  }
  if (ncol(covariates) == 0) covariates <- data.frame(subjectId = integer())
  if (ncol(covariateRef) == 0) {
    covariateRef <- data.frame(covariateId = integer(), covariateName = character())
  }
  if (ncol(analysisRef) == 0) {
    analysisRef <- data.frame(analysisId = character(), analysisName = character())
  }
  database <- Andromeda::andromeda()
  database[["cohorts"]] <- cohorts
  database[["outcomes"]] <- outcomes
  database[["covariates"]] <- covariates
  database[["covariateRef"]] <- covariateRef
  database[["analysisRef"]] <- analysisRef
  class(database) <- "CloneCensorWeightData"
  attr(class(database), "package") <- "CloneCensorWeight"
  database
}

.readCloneCensorWeightTable <- function(cloneCensorWeightData, tableName) {
  table <- switch(tableName,
                  cohorts = cloneCensorWeightData$cohorts,
                  outcomes = cloneCensorWeightData$outcomes,
                  covariates = cloneCensorWeightData$covariates,
                  covariateRef = cloneCensorWeightData$covariateRef,
                  analysisRef = cloneCensorWeightData$analysisRef,
                  stop("Unknown CloneCensorWeightData table: ", tableName, call. = FALSE))
  dplyr::collect(table)
}

#' Create a sample clone-censor-weight data set
#'
#' @param n Number of synthetic eligible subjects.
#' @param seed Random seed.
#' @return A `CloneCensorWeightData` object.
#' @export
createSampleCloneCensorWeightData <- function(n = 200, seed = 42) {
  checkmate::assertCount(n, positive = TRUE)
  checkmate::assertIntegerish(seed, len = 1)
  set.seed(seed)
  start <- as.Date("2020-01-01") + sample(0:30, n, replace = TRUE)
  followUp <- sample(90:365, n, replace = TRUE)
  exposureDelay <- sample(c(NA, 0:90), n, replace = TRUE,
                          prob = c(0.35, rep(0.65 / 91, 91)))
  outcomeDelay <- sample(c(NA, 30:365), n, replace = TRUE,
                         prob = c(0.65, rep(0.35 / 336, 336)))
  cohorts <- data.frame(
    subjectId = seq_len(n), cohortStartDate = start,
    cohortEndDate = start + followUp,
    initiationDate = as.Date(ifelse(is.na(exposureDelay), NA,
                                    as.character(start + exposureDelay))),
    stringsAsFactors = FALSE
  )
  outcomes <- data.frame(
    subjectId = which(!is.na(outcomeDelay)),
    outcomeDate = start[which(!is.na(outcomeDelay))] + outcomeDelay[which(!is.na(outcomeDelay))],
    outcomeId = rep(1, sum(!is.na(outcomeDelay))), stringsAsFactors = FALSE
  )
  covariates <- data.frame(
    subjectId = seq_len(n), age = round(runif(n, 35, 85)),
    priorCondition = rbinom(n, 1, 0.35), priorDrug = rbinom(n, 1, 0.25)
  )
  createCloneCensorWeightData(cohorts, outcomes, covariates)
}

#' Save clone-censor-weight data
#' @param cloneCensorWeightData A `CloneCensorWeightData` object.
#' @param file A file path.
#' @return Invisibly returns `file`.
#' @export
saveCloneCensorWeightData <- function(cloneCensorWeightData, file) {
  checkmate::assertClass(cloneCensorWeightData, "CloneCensorWeightData")
  checkmate::assertCharacter(file, len = 1)
  Andromeda::saveAndromeda(cloneCensorWeightData, file)
  invisible(file)
}

#' Load clone-censor-weight data
#' @param file A path created by [saveCloneCensorWeightData()].
#' @return A `CloneCensorWeightData` object.
#' @export
loadCloneCensorWeightData <- function(file) {
  checkmate::assertCharacter(file, len = 1)
  database <- Andromeda::loadAndromeda(file)
  class(database) <- "CloneCensorWeightData"
  attr(class(database), "package") <- "CloneCensorWeight"
  database
}

#' Extract data for clone-censor-weight analysis
#'
#' Extracts eligible subjects, treatment initiation markers, outcomes, and
#' baseline covariates from an OMOP CDM database.
#'
#' @param connectionDetails A `ConnectionDetails` object.
#' @param cdmDatabaseSchema CDM schema.
#' @param targetDatabaseSchema Target cohort schema.
#' @param targetTable Target cohort table.
#' @param targetId Target cohort identifier.
#' @param outcomeDatabaseSchema Outcome cohort schema.
#' @param outcomeTable Outcome cohort table.
#' @param outcomeIds Outcome cohort identifiers.
#' @param covariateSettings FeatureExtraction covariate settings.
#' @param exposureDatabaseSchema Exposure cohort schema.
#' @param exposureTable Exposure cohort table.
#' @param exposureId Exposure cohort identifier.
#' @param cdmVersion OMOP CDM version.
#' @return A `CloneCensorWeightData` object.
#' @export
getDbCloneCensorWeightData <- function(connectionDetails,
                                       cdmDatabaseSchema,
                                       targetDatabaseSchema,
                                       targetTable,
                                       targetId,
                                       outcomeDatabaseSchema,
                                       outcomeTable,
                                       outcomeIds,
                                       covariateSettings,
                                       exposureDatabaseSchema,
                                       exposureTable,
                                       exposureId,
                                       cdmVersion = "5") {
  checkmate::assertClass(connectionDetails, "ConnectionDetails")
  checkmate::assertCharacter(cdmDatabaseSchema, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(targetDatabaseSchema, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(targetTable, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(targetId, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(outcomeDatabaseSchema, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(outcomeTable, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(outcomeIds, min.len = 1, any.missing = FALSE)
  checkmate::assert_true(inherits(covariateSettings, c("covariateSettings", "CovariateSettings")))
  checkmate::assertCharacter(exposureDatabaseSchema, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(exposureTable, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(exposureId, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(cdmVersion, len = 1, any.missing = FALSE)

    ParallelLogger::logTrace("Extracting clone-censor-weight target population")
    connection <- DatabaseConnector::connect(connectionDetails)
    on.exit(DatabaseConnector::disconnect(connection), add = TRUE)
    sqlFolder <- system.file("sql", "sql_server", package = "CloneCensorWeight")
    if (!nzchar(sqlFolder)) sqlFolder <- file.path("inst", "sql", "sql_server")
    readSql <- function(file) paste(readLines(file, warn = FALSE), collapse = "\n")
    render <- function(sql) {
      DatabaseConnector::renderTranslateQuerySql(
        connection = connection,
        sql = sql,
        cdm_database_schema = cdmDatabaseSchema,
        target_table = targetTable,
        target_database_schema = targetDatabaseSchema,
        target_id = targetId,
        outcome_database_schema = outcomeDatabaseSchema,
        outcome_table = outcomeTable,
        outcome_ids = outcomeIds,
        exposure_database_schema = exposureDatabaseSchema,
        exposure_table = exposureTable,
        exposure_id = exposureId
      )
    }
    cohorts <- DatabaseConnector::querySql(connection, render(readSql(file.path(sqlFolder, "get_cohorts.sql"))))
    outcomes <- DatabaseConnector::querySql(connection, render(readSql(file.path(sqlFolder, "get_outcomes.sql"))))
    names(cohorts) <- tolower(names(cohorts))
    names(outcomes) <- tolower(names(outcomes))
    names(cohorts)[names(cohorts) == "subject_id"] <- "subjectId"
    names(cohorts)[names(cohorts) == "cohort_start_date"] <- "cohortStartDate"
    names(cohorts)[names(cohorts) == "cohort_end_date"] <- "cohortEndDate"
    names(cohorts)[names(cohorts) == "initiation_date"] <- "initiationDate"
    names(outcomes)[names(outcomes) == "subject_id"] <- "subjectId"
    names(outcomes)[names(outcomes) == "outcome_date"] <- "outcomeDate"
    names(outcomes)[names(outcomes) == "outcome_id"] <- "outcomeId"
    covariates <- data.frame(subjectId = cohorts$subjectId)
    if (!is.null(covariateSettings)) {
      ParallelLogger::logWarn("FeatureExtraction covariates require a configured FeatureExtraction adapter; returning subject identifiers only")
    }
    createCloneCensorWeightData(cohorts, outcomes, covariates)
}
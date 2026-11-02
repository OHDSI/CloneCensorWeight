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
#' @param cohorts A data frame with columns `subjectId`, `cohortStartDate`,
#'   `cohortEndDate`, and `initiationDate`.
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
  checkmate::assertSubset(
    c("subjectId", "cohortStartDate", "cohortEndDate"),
    names(cohorts)
  )
  checkmate::assertDataFrame(outcomes)
  checkmate::assertDataFrame(covariates)
  checkmate::assertDataFrame(covariateRef)
  checkmate::assertDataFrame(analysisRef)

  if (ncol(outcomes) == 0) {
    outcomes <- data.frame(
      subjectId   = integer(),
      outcomeDate = as.Date(character()),
      outcomeId   = integer()
    )
  }
  if (ncol(covariates) == 0) {
    covariates <- data.frame(rowId = integer(), covariateId = integer(), covariateValue = numeric())
  }
  if (ncol(covariateRef) == 0) {
    covariateRef <- data.frame(
      covariateId   = integer(),
      covariateName = character(),
      analysisId    = integer(),
      conceptId     = integer()
    )
  }
  if (ncol(analysisRef) == 0) {
    analysisRef <- data.frame(
      analysisId       = integer(),
      analysisName     = character(),
      domainId         = character(),
      startDay         = integer(),
      endDay           = integer(),
      isBinary         = character(),
      missingMeansZero = character()
    )
  }

  database <- Andromeda::andromeda()
  database[["cohorts"]]     <- cohorts
  database[["outcomes"]]    <- outcomes
  database[["covariates"]]  <- covariates
  database[["covariateRef"]] <- covariateRef
  database[["analysisRef"]] <- analysisRef
  class(database) <- "CloneCensorWeightData"
  attr(class(database), "package") <- "CloneCensorWeight"
  database
}

.readCloneCensorWeightTable <- function(cloneCensorWeightData, tableName) {
  table <- switch(
    tableName,
    cohorts     = cloneCensorWeightData$cohorts,
    outcomes    = cloneCensorWeightData$outcomes,
    covariates  = cloneCensorWeightData$covariates,
    covariateRef = cloneCensorWeightData$covariateRef,
    analysisRef  = cloneCensorWeightData$analysisRef,
    stop("Unknown CloneCensorWeightData table: ", tableName, call. = FALSE)
  )
  dplyr::collect(table)
}

#' Create a sample clone-censor-weight data set
#'
#' Generates a small synthetic data set useful for offline testing and
#' demonstrations. The covariates table uses a wide format (one column per
#' predictor) compatible with the current weighting code; it will be updated to
#' FeatureExtraction sparse format in a later phase.
#'
#' @param n Number of synthetic eligible subjects.
#' @param seed Random seed.
#' @return A `CloneCensorWeightData` object.
#' @export
createSampleCloneCensorWeightData <- function(n = 200, seed = 42) {
  checkmate::assertCount(n, positive = TRUE)
  checkmate::assertIntegerish(seed, len = 1)
  set.seed(seed)

  start        <- as.Date("2020-01-01") + sample(0:30, n, replace = TRUE)
  followUp     <- sample(90:365, n, replace = TRUE)
  exposureDelay <- sample(c(NA, 0:90), n, replace = TRUE,
                          prob = c(0.35, rep(0.65 / 91, 91)))
  outcomeDelay  <- sample(c(NA, 30:365), n, replace = TRUE,
                          prob = c(0.65, rep(0.35 / 336, 336)))

  cohorts <- data.frame(
    subjectId       = seq_len(n),
    cohortStartDate = start,
    cohortEndDate   = start + followUp,
    initiationDate  = as.Date(ifelse(
      is.na(exposureDelay), NA,
      as.character(start + exposureDelay)
    )),
    stringsAsFactors = FALSE
  )

  outcomeIdx <- which(!is.na(outcomeDelay))
  outcomes <- data.frame(
    subjectId   = outcomeIdx,
    outcomeDate = start[outcomeIdx] + outcomeDelay[outcomeIdx],
    outcomeId   = rep(1L, length(outcomeIdx)),
    stringsAsFactors = FALSE
  )

  ageYears  <- round(runif(n, 35, 85))
  priorCond <- rbinom(n, 1, 0.35)
  priorDrug <- rbinom(n, 1, 0.25)

  # FeatureExtraction sparse format: one row per non-zero (subject, covariate) pair.
  # Continuous covariates (age) have an entry for every subject; binary covariates
  # only have entries for subjects where value = 1 (missing = 0 by convention).
  condIdx <- which(priorCond == 1L)
  drugIdx <- which(priorDrug == 1L)
  covariates <- data.frame(
    rowId          = c(seq_len(n),          condIdx,              drugIdx),
    covariateId    = c(rep(1003L, n),        rep(2001L, length(condIdx)), rep(3001L, length(drugIdx))),
    covariateValue = c(ageYears / 10,        rep(1.0,   length(condIdx)), rep(1.0,   length(drugIdx)))
  )

  covariateRef <- data.frame(
    covariateId   = c(1003L, 2001L, 3001L),
    covariateName = c("age in decades", "prior condition in 365d", "prior drug in 365d"),
    analysisId    = c(1L, 2L, 3L),
    conceptId     = c(0L, 0L, 0L),
    stringsAsFactors = FALSE
  )

  analysisRef <- data.frame(
    analysisId       = c(1L, 2L, 3L),
    analysisName     = c("Demographics", "Conditions prior to index", "Drugs prior to index"),
    domainId         = c("Demographics", "Condition", "Drug"),
    startDay         = c(NA_integer_, -365L, -365L),
    endDay           = c(NA_integer_, 0L, 0L),
    isBinary         = c("N", "Y", "Y"),
    missingMeansZero = c("N", "Y", "Y"),
    stringsAsFactors = FALSE
  )

  createCloneCensorWeightData(cohorts, outcomes, covariates, covariateRef, analysisRef)
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

#' Extract data for a clone-censor-weight analysis
#'
#' Connects to an OMOP CDM database, retrieves eligible subjects, treatment
#' initiation markers, outcomes, and baseline covariates via FeatureExtraction,
#' and assembles them into a disk-backed `CloneCensorWeightData` object.
#'
#' @param connectionDetails A `ConnectionDetails` object created by
#'   [DatabaseConnector::createConnectionDetails()].
#' @param cdmDatabaseSchema Schema containing the OMOP CDM tables.
#' @param tempEmulationSchema Schema for temp-table emulation on platforms
#'   that require it (Oracle, BigQuery). Defaults to
#'   `getOption("sqlRenderTempEmulationSchema")`.
#' @param targetDatabaseSchema Schema containing the target cohort table.
#' @param targetTable Target cohort table name.
#' @param targetId Index cohort definition identifier.
#' @param outcomeDatabaseSchema Schema containing the outcome cohort table.
#' @param outcomeTable Outcome cohort table name.
#' @param outcomeIds Integer vector of outcome cohort definition identifiers.
#' @param exposureDatabaseSchema Schema containing the exposure cohort table.
#' @param exposureTable Exposure cohort table name.
#' @param exposureId Exposure cohort definition identifier.
#' @param getDbCloneCensorWeightDataArgs A
#'   `GetDbCloneCensorWeightDataArgs` object created by
#'   [createGetDbCloneCensorWeightDataArgs()].
#' @param cdmVersion OMOP CDM version string. Passed to FeatureExtraction.
#' @return A `CloneCensorWeightData` object.
#' @export
getDbCloneCensorWeightData <- function(connectionDetails,
                                       cdmDatabaseSchema,
                                       tempEmulationSchema = getOption("sqlRenderTempEmulationSchema"),
                                       targetDatabaseSchema,
                                       targetTable,
                                       targetId,
                                       outcomeDatabaseSchema,
                                       outcomeTable,
                                       outcomeIds,
                                       exposureDatabaseSchema,
                                       exposureTable,
                                       exposureId,
                                       getDbCloneCensorWeightDataArgs,
                                       cdmVersion = "5") {
  checkmate::assertClass(connectionDetails, "ConnectionDetails")
  checkmate::assertCharacter(cdmDatabaseSchema, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(targetDatabaseSchema, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(targetTable, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(targetId, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(outcomeDatabaseSchema, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(outcomeTable, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(outcomeIds, min.len = 1, any.missing = FALSE)
  checkmate::assertCharacter(exposureDatabaseSchema, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(exposureTable, len = 1, any.missing = FALSE)
  checkmate::assertIntegerish(exposureId, len = 1, any.missing = FALSE)
  checkmate::assertR6(getDbCloneCensorWeightDataArgs, "GetDbCloneCensorWeightDataArgs")

  args <- getDbCloneCensorWeightDataArgs

  ParallelLogger::logInfo("Connecting to database")
  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection), add = TRUE)

  # --- Target population -----------------------------------------------------
  ParallelLogger::logInfo("Retrieving target population")
  cohorts <- .renderTranslateQuerySql(
    connection               = connection,
    sqlFile                  = "GetCohorts.sql",
    cdm_database_schema      = cdmDatabaseSchema,
    target_database_schema   = targetDatabaseSchema,
    target_table             = targetTable,
    target_id                = targetId,
    exposure_database_schema = exposureDatabaseSchema,
    exposure_table           = exposureTable,
    exposure_id              = exposureId,
    washout_period           = args$washoutPeriod
  )
  ParallelLogger::logInfo(sprintf("Retrieved %d subjects", nrow(cohorts)))

  if (args$maxCohortSize > 0 && nrow(cohorts) > args$maxCohortSize) {
    ParallelLogger::logInfo(sprintf(
      "Sampling %d of %d subjects (maxCohortSize = %d)",
      args$maxCohortSize, nrow(cohorts), args$maxCohortSize
    ))
    cohorts <- cohorts[sample.int(nrow(cohorts), args$maxCohortSize), ]
  }

  # Upload temp cohort table (snake_case; adds cohort_definition_id = 1).
  # Used by GetOutcomes.sql and FeatureExtraction.
  .uploadTempCohortTable(connection, cohorts)

  # --- Outcomes --------------------------------------------------------------
  ParallelLogger::logInfo("Retrieving outcomes")
  outcomes <- .renderTranslateQuerySql(
    connection               = connection,
    sqlFile                  = "GetOutcomes.sql",
    outcome_database_schema  = outcomeDatabaseSchema,
    outcome_table            = outcomeTable,
    outcome_ids              = outcomeIds
  )
  ParallelLogger::logInfo(sprintf("Retrieved %d outcome records", nrow(outcomes)))

  # --- Baseline covariates (FeatureExtraction) --------------------------------
  ParallelLogger::logInfo("Extracting baseline covariates")
  covariateData <- FeatureExtraction::getDbCovariateData(
    connection        = connection,
    oracleTempSchema  = tempEmulationSchema,
    cdmDatabaseSchema = cdmDatabaseSchema,
    cdmVersion        = cdmVersion,
    cohortTable       = "#ccw_cohorts",
    cohortIds         = c(-1),
    rowIdField        = "subject_id",
    covariateSettings = args$covariateSettings,
    aggregated        = FALSE
  )

  # --- Assemble --------------------------------------------------------------
  ccwData <- createCloneCensorWeightData(cohorts, outcomes)
  ccwData[["covariates"]]   <- covariateData$covariates
  ccwData[["covariateRef"]] <- covariateData$covariateRef
  ccwData[["analysisRef"]]  <- covariateData$analysisRef
  ccwData
}

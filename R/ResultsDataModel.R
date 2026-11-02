#' Get the results data model specification
#'
#' Returns the column-level schema for all CloneCensorWeight results tables as a
#' data frame.  Each row describes one column: its table, name, SQL data type,
#' whether it is required, and whether it is part of the primary key.
#'
#' The spec is also available as `inst/csv/resultsDataModel.csv` for use with
#' `ResultModelManager`.
#'
#' @return A data frame with columns `tableName`, `columnName`, `dataType`,
#'   `isRequired`, and `primaryKey`.
#' @export
getResultsDataModelSpec <- function() {
  path <- system.file("csv", "resultsDataModel.csv", package = "CloneCensorWeight")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

#' Create the results data model tables
#'
#' Executes the DDL that creates all `ccw_*` results tables in the target
#' schema.  Existing tables with the same names are **dropped and recreated**.
#'
#' @param connection A `DatabaseConnector` connection object.
#' @param databaseSchema Schema in which to create tables, e.g.
#'   `"my_results"`.  Pass `""` to use the default schema (not recommended for
#'   production).
#' @param tablePrefix An optional string prepended to every table name, e.g.
#'   `"study1_"` yields `study1_ccw_analysis`, etc.  Default `""`.
#'
#' @return Invisibly `NULL`.  Called for its side-effect.
#' @export
createResultsDataModel <- function(connection, databaseSchema, tablePrefix = "") {
  checkmate::assertString(databaseSchema)
  checkmate::assertString(tablePrefix)

  schemaPrefix <- if (nchar(databaseSchema) > 0) paste0(databaseSchema, ".") else ""

  sql <- SqlRender::loadRenderTranslateSql(
    sqlFilename    = "CreateResultsDataModel.sql",
    packageName    = "CloneCensorWeight",
    database_schema = schemaPrefix,
    table_prefix   = tablePrefix,
    dbms           = DatabaseConnector::dbms(connection)
  )
  DatabaseConnector::executeSql(connection, sql, progressBar = FALSE, reportOverallTime = FALSE)
  invisible(NULL)
}

#' Drop the results data model tables
#'
#' Removes all `ccw_*` results tables (with optional prefix) from the target
#' schema.  Tables that do not exist are silently ignored.
#'
#' @inheritParams createResultsDataModel
#' @return Invisibly `NULL`.  Called for its side-effect.
#' @export
dropResultsDataModel <- function(connection, databaseSchema, tablePrefix = "") {
  checkmate::assertString(databaseSchema)
  checkmate::assertString(tablePrefix)

  spec         <- getResultsDataModelSpec()
  tableNames   <- unique(spec$tableName)
  schemaPrefix <- if (nchar(databaseSchema) > 0) paste0(databaseSchema, ".") else ""
  dbms         <- DatabaseConnector::dbms(connection)

  for (tbl in tableNames) {
    fullName <- paste0(schemaPrefix, tablePrefix, tbl)
    sql      <- SqlRender::translate(
      paste0("DROP TABLE IF EXISTS ", fullName, ";"),
      targetDialect = dbms
    )
    DatabaseConnector::executeSql(connection, sql,
                                  progressBar = FALSE, reportOverallTime = FALSE)
  }
  invisible(NULL)
}

#' Upload analysis results to the results data model
#'
#' Inserts a named list of result data frames — typically produced by
#' [runCcwAnalyses()] — into the corresponding `ccw_*` tables.  The `analysisId`
#' is stamped onto every row before insertion so results from multiple studies
#' can coexist in the same schema.
#'
#' Supported list names (all optional; omitted elements are skipped):
#' \describe{
#'   \item{`analysis`}{Single-row data frame with fields matching
#'     `ccw_analysis` (minus `analysis_id`).}
#'   \item{`cohortDefinition`}{Data frame matching `ccw_cohort_definition`.}
#'   \item{`cohortCount`}{Data frame matching `ccw_cohort_count`.}
#'   \item{`outcomeModel`}{`outcomeModel$summary` from [fitOutcomeModel()], or
#'     the calibrated summary from [calibrateEstimates()].}
#'   \item{`standardizedDifference`}{Data frame from
#'     [computeStandardizedDifferences()].}
#'   \item{`kaplanMeier`}{Data frame from [computeKaplanMeier()].}
#'   \item{`powerMdrr`}{Data frame from [computeMdrr()].}
#'   \item{`negativeControl`}{Data frame with `logRr` and `seLogRr`, and
#'     optionally a `controlId` column.}
#' }
#'
#' @param connection A `DatabaseConnector` connection object.
#' @param databaseSchema Schema containing the results tables.
#' @param analysisId A character string uniquely identifying this analysis run.
#' @param results A named list of data frames as described above.
#' @param tablePrefix Optional table-name prefix matching the one used in
#'   [createResultsDataModel()].
#'
#' @return Invisibly `NULL`.  Called for its side-effect.
#' @export
uploadResults <- function(connection, databaseSchema, analysisId, results,
                          tablePrefix = "") {
  checkmate::assertString(databaseSchema)
  checkmate::assertString(analysisId, min.chars = 1)
  checkmate::assertList(results, names = "unique")
  checkmate::assertString(tablePrefix)

  schemaPrefix <- if (nchar(databaseSchema) > 0) paste0(databaseSchema, ".") else ""

  .insert <- function(tableName, df) {
    df$analysis_id <- analysisId
    DatabaseConnector::insertTable(
      connection     = connection,
      tableName      = paste0(schemaPrefix, tablePrefix, tableName),
      data           = df,
      dropTableIfExists = FALSE,
      createTable    = FALSE,
      progressBar    = FALSE
    )
  }

  if (!is.null(results$analysis)) {
    df <- results$analysis
    checkmate::assertDataFrame(df)
    .insert("ccw_analysis", df)
  }

  if (!is.null(results$cohortDefinition)) {
    .insert("ccw_cohort_definition", results$cohortDefinition)
  }

  if (!is.null(results$cohortCount)) {
    .insert("ccw_cohort_count", results$cohortCount)
  }

  if (!is.null(results$outcomeModel)) {
    df <- results$outcomeModel
    checkmate::assertDataFrame(df)
    # Rename camelCase summary columns to snake_case for the DB table
    df <- .renameOutcomeModelCols(df)
    .insert("ccw_outcome_model", df)
  }

  if (!is.null(results$standardizedDifference)) {
    df <- results$standardizedDifference
    checkmate::assertDataFrame(df)
    df <- .renameSmdCols(df)
    .insert("ccw_standardized_difference", df)
  }

  if (!is.null(results$kaplanMeier)) {
    df <- results$kaplanMeier
    checkmate::assertDataFrame(df)
    df <- .renameKmCols(df)
    .insert("ccw_kaplan_meier", df)
  }

  if (!is.null(results$powerMdrr)) {
    df <- results$powerMdrr
    checkmate::assertDataFrame(df)
    df <- .renameMdrrCols(df)
    .insert("ccw_power_mdrr", df)
  }

  if (!is.null(results$negativeControl)) {
    df <- results$negativeControl
    checkmate::assertDataFrame(df)
    if (!"control_id" %in% names(df) && !"controlId" %in% names(df)) {
      df$control_id <- as.character(seq_len(nrow(df)))
    }
    df <- .renameNcCols(df)
    .insert("ccw_negative_control", df)
  }

  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Internal column-rename helpers (camelCase R → snake_case DB)
# ---------------------------------------------------------------------------

.renameOutcomeModelCols <- function(df) {
  map <- c(
    modelType          = "model_type",
    estimand           = "estimand",
    varianceEstimator  = "variance_estimator",
    arm1               = "arm1",
    arm2               = "arm2",
    estimate           = "estimate",
    lower95            = "lower_95",
    upper95            = "upper_95",
    pValue             = "p_value",
    logLikelihood      = "log_likelihood",
    calibratedP        = "calibrated_p",
    calibratedLower95  = "calibrated_lower_95",
    calibratedUpper95  = "calibrated_upper_95"
  )
  .applyRenameMap(df, map)
}

.renameSmdCols <- function(df) {
  map <- c(
    covariateId = "covariate_id",
    covariate   = "covariate",
    weighted    = "weighted",
    meanArm1    = "mean_arm1",
    sdArm1      = "sd_arm1",
    meanArm2    = "mean_arm2",
    sdArm2      = "sd_arm2",
    smd         = "smd",
    absSmd      = "abs_smd"
  )
  .applyRenameMap(df, map)
}

.renameKmCols <- function(df) {
  map <- c(
    arm      = "arm",
    time     = "time",
    survival = "survival",
    lower    = "lower",
    upper    = "upper",
    nRisk    = "n_risk",
    nEvent   = "n_event"
  )
  .applyRenameMap(df, map)
}

.renameMdrrCols <- function(df) {
  map <- c(
    arm1            = "arm1",
    arm2            = "arm2",
    alpha           = "alpha",
    power           = "power",
    nEventsArm1     = "n_events_arm1",
    nEventsArm2     = "n_events_arm2",
    nEventsEffArm1  = "n_events_eff_arm1",
    nEventsEffArm2  = "n_events_eff_arm2",
    mdrr            = "mdrr"
  )
  .applyRenameMap(df, map)
}

.renameNcCols <- function(df) {
  map <- c(
    controlId = "control_id",
    logRr     = "log_rr",
    seLogRr   = "se_log_rr"
  )
  .applyRenameMap(df, map)
}

# Apply a rename map (names = old, values = new) to a data frame.
# Columns not in the map are kept as-is.
.applyRenameMap <- function(df, map) {
  presentOld <- intersect(names(map), names(df))
  idx        <- match(presentOld, names(df))
  names(df)[idx] <- map[presentOld]
  df
}

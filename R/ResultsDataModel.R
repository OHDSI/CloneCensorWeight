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

# The results data model follows CohortMethod's (R/ResultsDataModel.R, Copyright Observational Health
# Data Sciences and Informatics, Apache 2.0), with the ccw_ table prefix.

#' Get the specifications of the results data model
#'
#' @details
#' The tables of the results data model all start with `ccw_`. `database_id` is part of the primary key
#' of every table that holds site-specific results; the definition tables `ccw_analysis`,
#' `ccw_strategy`, `ccw_target_treatment` and `ccw_target_treatment_outcome` have none and are the same
#' for every site. Tables per strategy have `strategy_id`, and tables that compare two strategies
#' have `target_strategy_id` and `comparator_strategy_id`. Counts in the columns with
#' `minCellCount = "Yes"` below the minimum cell count of [exportToCsv()] are written as
#' `-minCellCount`.
#'
#' @return
#' A tibble with the columns `tableName`, `columnName`, `dataType`, `isRequired`, `primaryKey`,
#' `minCellCount`, `deprecated` and `description`, one row per column of each table.
#'
#' @examples
#' specifications <- getResultsDataModelSpecifications()
#' unique(specifications$tableName)
#'
#' @export
getResultsDataModelSpecifications <- function() {
  pathToCsv <- system.file("csv", "resultsDataModelSpecification.csv", package = "CloneCensorWeight", mustWork = TRUE)
  specifications <- readr::read_csv(file = pathToCsv, col_types = readr::cols(.default = readr::col_character()))
  colnames(specifications) <- SqlRender::snakeCaseToCamelCase(colnames(specifications))
  return(specifications)
}

#' Create the results data model tables on a database server
#'
#' @details
#' Creates the tables of [getResultsDataModelSpecifications()] with the DDL of the first version of the
#' model, then applies the migrations to the current version with [migrateDataModel()]. Only
#' PostgreSQL and SQLite servers are supported.
#'
#' @param connectionDetails An object of type `ConnectionDetails` as created using
#'                          [DatabaseConnector::createConnectionDetails()].
#' @param databaseSchema    The schema on the server where the tables will be created. On SQLite, this
#'                          must be `"main"`.
#' @param tablePrefix       (Optional) A string to insert before the table names (which already start
#'                          with `ccw_`).
#'
#' @return
#' Does not return anything. Is called for its side effect of creating the tables.
#'
#' @examples
#' \dontrun{
#' connectionDetails <- DatabaseConnector::createConnectionDetails(
#'   dbms = "sqlite",
#'   server = tempfile(fileext = ".sqlite")
#' )
#' createResultsDataModel(connectionDetails, databaseSchema = "main")
#' }
#'
#' @export
createResultsDataModel <- function(connectionDetails = NULL,
                                   databaseSchema,
                                   tablePrefix = "") {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertClass(connectionDetails, "ConnectionDetails", add = errorMessages)
  checkmate::assertCharacter(databaseSchema, len = 1, add = errorMessages)
  checkmate::assertCharacter(tablePrefix, len = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (connectionDetails$dbms == "sqlite" && databaseSchema != "main") {
    stop("Invalid schema for sqlite, use databaseSchema = 'main'")
  }
  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection))
  # The first version of the results model:
  sql <- .loadRenderTranslateSql(
    sqlFilename = "CreateResultsDataModel.sql",
    dbms = DatabaseConnector::dbms(connection),
    tempEmulationSchema = NULL,
    database_schema = databaseSchema,
    table_prefix = tablePrefix
  )
  DatabaseConnector::executeSql(connection, sql, progressBar = FALSE, reportOverallTime = FALSE)
  # Migrate to the current version:
  migrateDataModel(
    connectionDetails = connectionDetails,
    databaseSchema = databaseSchema,
    tablePrefix = tablePrefix
  )
  invisible(NULL)
}

#' Upload results to the database server
#'
#' @description
#' Requires that the results data model tables have been created using [createResultsDataModel()].
#'
#' @details
#' Uses [ResultModelManager::uploadResults()]. Rows of the site-specific tables replace the rows of the
#' same site with the same primary key. Rows of the definition tables whose primary key is already in
#' the database are not uploaded, unless `forceOverWriteOfSpecifications = TRUE`.
#'
#' @param connectionDetails              An object of type `ConnectionDetails` as created using
#'                                       [DatabaseConnector::createConnectionDetails()].
#' @param schema                         The schema on the server where the tables have been created.
#' @param zipFileName                    The name of the zip file created by [exportToCsv()].
#' @param forceOverWriteOfSpecifications If `TRUE`, the rows of the definition tables (analyses,
#'                                       strategies, and target-treatment-outcome combinations) that
#'                                       already exist in the database are overwritten. Only use this if
#'                                       these definitions have changed since the last upload.
#' @param purgeSiteDataBeforeUploading   If `TRUE`, all the data of the site (`database_id`) is deleted
#'                                       before uploading. This assumes the zip file contains the full
#'                                       data of that site, and needs the `databaseIdentifierFile` of
#'                                       [ResultModelManager::uploadResults()] (passed through `...`).
#' @param tempFolder                     A folder on the local file system where the zip file is extracted
#'                                       to. It is cleaned up when the function is finished.
#' @param tablePrefix                    (Optional) A string to insert before the table names.
#' @param ...                            Other arguments of [ResultModelManager::uploadResults()].
#'
#' @return
#' Does not return anything. Is called for its side effect of uploading the results.
#'
#' @examples
#' \dontrun{
#' uploadResults(
#'   connectionDetails = connectionDetails,
#'   schema = "main",
#'   zipFileName = "export/Results_MyDatabase.zip",
#'   purgeSiteDataBeforeUploading = FALSE
#' )
#' }
#'
#' @export
uploadResults <- function(connectionDetails,
                          schema,
                          zipFileName,
                          forceOverWriteOfSpecifications = FALSE,
                          purgeSiteDataBeforeUploading = TRUE,
                          tempFolder = tempdir(),
                          tablePrefix = "",
                          ...) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertClass(connectionDetails, "ConnectionDetails", add = errorMessages)
  checkmate::assertCharacter(schema, len = 1, add = errorMessages)
  checkmate::assertFileExists(zipFileName, add = errorMessages)
  checkmate::assertLogical(forceOverWriteOfSpecifications, len = 1, add = errorMessages)
  checkmate::assertLogical(purgeSiteDataBeforeUploading, len = 1, add = errorMessages)
  checkmate::assertCharacter(tablePrefix, len = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  rlang::check_installed("zip", reason = "to unzip the results")
  unzipFolder <- tempfile("unzipTempFolder", tmpdir = tempFolder)
  dir.create(path = unzipFolder, recursive = TRUE)
  on.exit(unlink(unzipFolder, recursive = TRUE), add = TRUE)

  ParallelLogger::logInfo("Unzipping ", zipFileName)
  zip::unzip(zipFileName, exdir = unzipFolder)

  ResultModelManager::uploadResults(
    connectionDetails = connectionDetails,
    schema = schema,
    resultsFolder = unzipFolder,
    tablePrefix = tablePrefix,
    forceOverWriteOfSpecifications = forceOverWriteOfSpecifications,
    purgeSiteDataBeforeUploading = purgeSiteDataBeforeUploading,
    runCheckAndFixCommands = FALSE,
    specifications = getResultsDataModelSpecifications(),
    warnOnMissingTable = FALSE,
    ...
  )
  invisible(NULL)
}

#' Find the results zip file in an export folder
#'
#' @description
#' Finds the `Results_<databaseId>.zip` file that [exportToCsv()] made in `exportFolder`, as the
#' Strategus module does before [uploadResults()].
#'
#' @param exportFolder The folder that [exportToCsv()] exported to.
#'
#' @return
#' The full path of the one `Results_*.zip` file in `exportFolder`, or `NULL` when there is none. Stops
#' with an error naming them when there is more than one (for example the zip of an earlier run, with
#' another `databaseId`), because which one is meant is then ambiguous.
#'
#' @examples
#' \dontrun{
#' zipFileName <- findResultsZipFile("export")
#' }
#'
#' @export
findResultsZipFile <- function(exportFolder) {
  checkmate::assertDirectoryExists(exportFolder)
  exportFolder <- normalizePath(exportFolder)
  zipFiles <- list.files(path = exportFolder, pattern = "^Results_.*\\.zip$", full.names = TRUE)
  if (length(zipFiles) > 1) {
    stop(
      "The results folder ", exportFolder, " has more than one results zip file (",
      paste(basename(zipFiles), collapse = ", "), "), so which one this run made is ambiguous. ",
      "Remove those of earlier runs and upload again."
    )
  }
  if (length(zipFiles) == 0) {
    return(NULL)
  }
  return(zipFiles)
}

#' Migrate the results data model
#'
#' @description
#' Migrates the results data model tables from their current version to the version of this package,
#' applying the migrations in `inst/sql/sql_server/migrations` that have not been applied yet. The
#' applied migrations are recorded in the `ccw_migration` table, which [createResultsDataModel()]
#' creates. This version of the package has no migrations: its tables are those of the first version.
#'
#' It is strongly advised to back up the data (the SQLite file, a backup of the PostgreSQL database, or
#' the zip files of [exportToCsv()]) before migrating.
#'
#' @inheritParams getDataMigrator
#'
#' @return
#' Does not return anything. Is called for its side effect of migrating the tables.
#'
#' @export
migrateDataModel <- function(connectionDetails, databaseSchema, tablePrefix = "") {
  ParallelLogger::logInfo("Migrating data set")
  migrator <- getDataMigrator(
    connectionDetails = connectionDetails,
    databaseSchema = databaseSchema,
    tablePrefix = tablePrefix
  )
  on.exit(migrator$closeConnection())
  # ResultModelManager 0.6.2 fails when a package has no migrations yet, as this version has none:
  if (nrow(migrator$getStatus()) > 0) {
    migrator$executeMigrations()
  }
  invisible(NULL)
}

#' Get the data migrator of the results data model
#'
#' @description
#' Returns a `ResultModelManager::DataMigrationManager` for the results data model tables, which can
#' report the migrations already applied (`getStatus()`) and apply the others (`executeMigrations()`).
#'
#' @param connectionDetails An object of type `ConnectionDetails` as created using
#'                          [DatabaseConnector::createConnectionDetails()].
#' @param databaseSchema    The schema on the server where the tables are.
#' @param tablePrefix       (Optional) A string inserted before the table names.
#'
#' @return
#' An instance of `ResultModelManager::DataMigrationManager`.
#'
#' @export
getDataMigrator <- function(connectionDetails, databaseSchema, tablePrefix = "") {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertClass(connectionDetails, "ConnectionDetails", add = errorMessages)
  checkmate::assertCharacter(databaseSchema, len = 1, add = errorMessages)
  checkmate::assertCharacter(tablePrefix, len = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  return(ResultModelManager::DataMigrationManager$new(
    connectionDetails = connectionDetails,
    databaseSchema = databaseSchema,
    tablePrefix = tablePrefix,
    packageTablePrefix = "ccw_",
    migrationPath = "migrations",
    packageName = "CloneCensorWeight"
  ))
}

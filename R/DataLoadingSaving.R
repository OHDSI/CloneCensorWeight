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

#' Get the clone-censor-weight data from the server
#'
#' @description
#' Executes a set of SQL statements against a database in OMOP CDM format to extract the data needed
#' for a clone-censor-weight analysis: the eligible cohort, the eras of the treatment cohorts, the
#' outcomes, death, and baseline covariates.
#'
#' @details
#' The eligible ("target") cohort defines eligibility and time zero: its `cohort_start_date` is the
#' index date (day 0). The eligible cohort is restricted in this order:
#'
#' 1. If `firstOnly`, only the first cohort entry of each person is kept.
#' 2. The index date must fall in an observation period with at least `washoutPeriod` days of
#'    observation before it.
#' 3. The index date must fall between `studyStartDate` and `studyEndDate`. As in CohortMethod, the
#'    end of the cohort and of the observation period are truncated at `studyEndDate`.
#' 4. If `removeSubjectsWithPriorTreatment`, index dates with prior use of either treatment are
#'    removed. Prior use is an era of the treatment or comparator treatment cohort that overlaps days
#'    `-priorUseLookback` to `-1` relative to the index date, so a start on day 0 is allowed.
#' 5. If the cohort is larger than `maxCohortSize` (and `maxCohortSize > 0`), a random sample of that
#'    size is drawn.
#'
#' The number of persons and index dates after each step is recorded in the attrition table, in
#' CohortMethod's format (the counts are in the target columns; the comparator columns are 0).
#'
#' The treatment cohorts are not restricted to new users: every era of the treatment (and comparator
#' treatment) cohort that overlaps days `-priorUseLookback` to the end of observation is stored in the
#' `treatments` table, so that strategies can be evaluated on the first start on or after day 0. The
#' treatment cohorts must therefore hold every era of each person, not only the first. A cohort that
#' keeps one era per person is accepted without a warning, but it misses the starts after index of
#' persons with an earlier era (who then count as not starting) and their prior use. Cohort
#' definitions often keep only the first event by default: with Capr, for example, set
#' `entry(primaryCriteriaLimit = "All", qualifiedLimit = "All")` and
#' `attrition(expressionLimit = "All")`.
#' Outcomes are stored as in CohortMethod: every outcome cohort start in the index observation period,
#' including those before the index date. `daysToDeath` is the number of days from the index date to the
#' person's earliest death date in the CDM `death` table (`NA` if none); it is not restricted to the
#' observation period. Covariates are constructed by FeatureExtraction at the index date.
#'
#' All cohorts must be in tables with the format of the CDM `cohort` table: `cohort_definition_id`,
#' `subject_id`, `cohort_start_date`, and `cohort_end_date`.
#'
#' Either `connectionDetails` or `connection` must be provided. An open `connection` is used as is and
#' is not closed, which allows, for example, a DuckDB database attached read-only to the connection.
#'
#' As in CohortMethod, the eligible cohort defaults to the table of the treatment cohorts
#' (`targetDatabaseSchema = exposureDatabaseSchema`, `targetTable = exposureTable`), but the outcome
#' cohorts default to the `cohort` table of the CDM schema (`outcomeDatabaseSchema =
#' cdmDatabaseSchema`, `outcomeTable = "cohort"`). When all cohorts are in one table in another
#' schema, also set `outcomeDatabaseSchema` and `outcomeTable` (see the second example).
#'
#' @param connectionDetails              An R object of type `connectionDetails` created using the
#'                                       [DatabaseConnector::createConnectionDetails()] function. Either
#'                                       `connectionDetails` or `connection` must be provided.
#' @param connection                     An open connection created using
#'                                       [DatabaseConnector::connect()]. It is not closed when the
#'                                       function returns. Either `connectionDetails` or `connection`
#'                                       must be provided.
#' @param cdmDatabaseSchema              The name of the database schema that contains the OMOP CDM
#'                                       instance. Requires read permissions to this database. On SQL
#'                                       Server, this should specify both the database and the schema,
#'                                       so for example 'cdm_instance.dbo'.
#' @param tempEmulationSchema            Some database platforms like Oracle and Impala do not truly
#'                                       support temp tables. To emulate temp tables, provide a schema
#'                                       with write privileges where temp tables can be created.
#' @param targetId                       The cohort definition ID of the eligible ("target") cohort in
#'                                       `targetTable`.
#' @param treatmentId                    The cohort definition ID of the treatment cohort (role
#'                                       `"treatment"`) in `exposureTable`.
#' @param comparatorTreatmentId          The cohort definition ID of the comparator treatment cohort
#'                                       (role `"comparatorTreatment"`) in `exposureTable`, or `NULL`
#'                                       if the strategies refer to one treatment only.
#' @param outcomeIds                     The cohort definition IDs of the outcome cohorts in
#'                                       `outcomeTable`.
#' @param exposureDatabaseSchema         The name of the database schema that contains the treatment
#'                                       cohorts.
#' @param exposureTable                  The name of the table that contains the treatment cohorts.
#' @param outcomeDatabaseSchema          The name of the database schema that contains the outcome
#'                                       cohorts.
#' @param outcomeTable                   The name of the table that contains the outcome cohorts.
#' @param targetDatabaseSchema           The name of the database schema that contains the eligible
#'                                       cohort.
#' @param targetTable                    The name of the table that contains the eligible cohort.
#' @param cdmVersion                     The version of the Common Data Model used. Currently only
#'                                       version 5 is supported.
#' @param getDbCloneCensorWeightDataArgs An object of type `GetDbCloneCensorWeightDataArgs` as created
#'                                       using the [createGetDbCloneCensorWeightDataArgs()] function.
#'
#' @return
#' An object of type [CloneCensorWeightData].
#'
#' @examples
#' \dontrun{
#' connectionDetails <- Eunomia::getEunomiaConnectionDetails()
#' Eunomia::createCohorts(connectionDetails)
#' # Eligible: new users of NSAIDs; treatments: celecoxib and diclofenac; outcome: GI bleed
#' data <- getDbCloneCensorWeightData(
#'   connectionDetails = connectionDetails,
#'   cdmDatabaseSchema = "main",
#'   targetId = 4,
#'   treatmentId = 1,
#'   comparatorTreatmentId = 2,
#'   outcomeIds = 3,
#'   getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs()
#' )
#' summary(data)
#'
#' # A DuckDB CDM attached read-only to an in-memory database, with all cohorts in one table of the
#' # in-memory database (generated there beforehand), so the eligible, treatment and outcome cohort
#' # tables are all set:
#' connection <- DatabaseConnector::connect(dbms = "duckdb", server = ":memory:")
#' DatabaseConnector::executeSql(connection, "ATTACH 'cdm.duckdb' AS cdm (READ_ONLY);")
#' data <- getDbCloneCensorWeightData(
#'   connection = connection,
#'   cdmDatabaseSchema = "cdm.main",
#'   targetId = 1,
#'   treatmentId = 2,
#'   outcomeIds = 3,
#'   exposureDatabaseSchema = "main",
#'   exposureTable = "my_cohorts",
#'   outcomeDatabaseSchema = "main",
#'   outcomeTable = "my_cohorts",
#'   targetDatabaseSchema = "main",
#'   targetTable = "my_cohorts",
#'   getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs()
#' )
#' DatabaseConnector::disconnect(connection)
#' }
#'
#' @export
getDbCloneCensorWeightData <- function(connectionDetails = NULL,
                                       connection = NULL,
                                       cdmDatabaseSchema,
                                       tempEmulationSchema = getOption("sqlRenderTempEmulationSchema"),
                                       targetId,
                                       treatmentId,
                                       comparatorTreatmentId = NULL,
                                       outcomeIds,
                                       exposureDatabaseSchema = cdmDatabaseSchema,
                                       exposureTable = "cohort",
                                       outcomeDatabaseSchema = cdmDatabaseSchema,
                                       outcomeTable = "cohort",
                                       targetDatabaseSchema = exposureDatabaseSchema,
                                       targetTable = exposureTable,
                                       cdmVersion = "5",
                                       getDbCloneCensorWeightDataArgs) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertClass(connectionDetails, "ConnectionDetails", null.ok = TRUE, add = errorMessages)
  checkmate::assertClass(connection, "DatabaseConnectorConnection", null.ok = TRUE, add = errorMessages)
  checkmate::assertCharacter(cdmDatabaseSchema, len = 1, add = errorMessages)
  checkmate::assertCharacter(tempEmulationSchema, len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumber(targetId, add = errorMessages)
  checkmate::assertNumber(treatmentId, add = errorMessages)
  checkmate::assertNumber(comparatorTreatmentId, null.ok = TRUE, add = errorMessages)
  checkmate::assertNumeric(outcomeIds, any.missing = FALSE, min.len = 1, add = errorMessages)
  ids <- c(targetId, treatmentId, comparatorTreatmentId, outcomeIds)
  if (is.numeric(ids)) {
    checkmate::assertTRUE(all(ids %% 1 == 0), .var.name = "all cohort IDs are integers", add = errorMessages)
  }
  checkmate::assertCharacter(exposureDatabaseSchema, len = 1, add = errorMessages)
  checkmate::assertCharacter(exposureTable, len = 1, add = errorMessages)
  checkmate::assertCharacter(outcomeDatabaseSchema, len = 1, add = errorMessages)
  checkmate::assertCharacter(outcomeTable, len = 1, add = errorMessages)
  checkmate::assertCharacter(targetDatabaseSchema, len = 1, add = errorMessages)
  checkmate::assertCharacter(targetTable, len = 1, add = errorMessages)
  checkmate::assertCharacter(cdmVersion, len = 1, add = errorMessages)
  checkmate::assertR6(getDbCloneCensorWeightDataArgs, "GetDbCloneCensorWeightDataArgs", add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (is.null(connectionDetails) == is.null(connection)) {
    stop("Provide either connectionDetails or connection, but not both")
  }
  if (!is.null(comparatorTreatmentId) && comparatorTreatmentId == treatmentId) {
    stop("treatmentId and comparatorTreatmentId must be different")
  }

  if (is.null(connection)) {
    connection <- DatabaseConnector::connect(connectionDetails)
    on.exit(DatabaseConnector::disconnect(connection))
  }
  dbms <- DatabaseConnector::dbms(connection)
  args <- getDbCloneCensorWeightDataArgs
  treatmentIds <- c(treatment = as.numeric(treatmentId), comparatorTreatment = as.numeric(comparatorTreatmentId))
  removePriorUsers <- args$removeSubjectsWithPriorTreatment && args$priorUseLookback > 0

  message("Constructing eligible cohort")
  cohortSqlArgs <- list(
    sqlFilename = "CreateOrCountCohorts.sql",
    dbms = dbms,
    tempEmulationSchema = tempEmulationSchema,
    cdm_database_schema = cdmDatabaseSchema,
    target_database_schema = targetDatabaseSchema,
    target_table = targetTable,
    exposure_database_schema = exposureDatabaseSchema,
    exposure_table = exposureTable,
    target_id = targetId,
    treatment_ids = treatmentIds,
    study_start_date = args$studyStartDate,
    study_end_date = args$studyEndDate,
    first_only = args$firstOnly,
    washout_period = args$washoutPeriod,
    remove_subjects_with_prior_treatment = removePriorUsers,
    prior_use_lookback = args$priorUseLookback
  )
  sql <- do.call(.loadRenderTranslateSql, c(cohortSqlArgs, action = "CREATE"))
  DatabaseConnector::executeSql(connection, sql, progressBar = FALSE, reportOverallTime = FALSE)
  sql <- do.call(.loadRenderTranslateSql, c(cohortSqlArgs, action = "COUNT"))
  counts <- DatabaseConnector::querySql(connection, sql, snakeCaseToCamelCase = TRUE)
  counts <- counts[order(counts$seqId), ]
  attrition <- tibble(
    description = counts$description,
    targetPersons = as.numeric(counts$persons),
    comparatorPersons = 0,
    targetExposures = as.numeric(counts$exposures),
    comparatorExposures = 0
  )

  cohortSize <- attrition$targetExposures[nrow(attrition)]
  sampled <- args$maxCohortSize != 0 && cohortSize > args$maxCohortSize
  if (sampled) {
    message("Downsampling eligible cohort from ", cohortSize, " to ", args$maxCohortSize)
    sql <- .loadRenderTranslateSql(
      sqlFilename = "SampleCohorts.sql",
      dbms = dbms,
      tempEmulationSchema = tempEmulationSchema,
      max_cohort_size = args$maxCohortSize
    )
    DatabaseConnector::executeSql(connection, sql, progressBar = FALSE, reportOverallTime = FALSE)
    cohortTable <- "#ccw_sample"
  } else {
    cohortTable <- "#ccw_cohort"
  }
  if (cohortSize == 0) {
    warning("Eligible cohort is empty")
  }

  message("Constructing covariates")
  covariateSettings <- .setCovariateCohortTable(args$covariateSettings, exposureDatabaseSchema, exposureTable)
  data <- FeatureExtraction::getDbCovariateData(
    connection = connection,
    tempEmulationSchema = tempEmulationSchema,
    cdmDatabaseSchema = cdmDatabaseSchema,
    cdmVersion = cdmVersion,
    cohortTable = cohortTable,
    cohortTableIsTemp = TRUE,
    rowIdField = "row_id",
    covariateSettings = covariateSettings
  )

  message("Fetching cohorts, treatments, and outcomes from server")
  start <- Sys.time()
  sql <- .loadRenderTranslateSql(
    sqlFilename = "GetCohorts.sql",
    dbms = dbms,
    tempEmulationSchema = tempEmulationSchema,
    cdm_database_schema = cdmDatabaseSchema,
    cohort_table = cohortTable
  )
  DatabaseConnector::querySqlToAndromeda(connection, sql, andromeda = data, andromedaTableName = "cohorts", snakeCaseToCamelCase = TRUE)
  sql <- .loadRenderTranslateSql(
    sqlFilename = "GetTreatments.sql",
    dbms = dbms,
    tempEmulationSchema = tempEmulationSchema,
    exposure_database_schema = exposureDatabaseSchema,
    exposure_table = exposureTable,
    treatment_ids = treatmentIds,
    prior_use_lookback = args$priorUseLookback,
    cohort_table = cohortTable
  )
  DatabaseConnector::querySqlToAndromeda(connection, sql, andromeda = data, andromedaTableName = "treatments", snakeCaseToCamelCase = TRUE)
  sql <- .loadRenderTranslateSql(
    sqlFilename = "GetOutcomes.sql",
    dbms = dbms,
    tempEmulationSchema = tempEmulationSchema,
    outcome_database_schema = outcomeDatabaseSchema,
    outcome_table = outcomeTable,
    outcome_ids = outcomeIds,
    cohort_table = cohortTable
  )
  DatabaseConnector::querySqlToAndromeda(connection, sql, andromeda = data, andromedaTableName = "outcomes", snakeCaseToCamelCase = TRUE)
  delta <- Sys.time() - start
  message("Fetching data took ", signif(delta, 3), " ", attr(delta, "units"))

  sql <- .loadRenderTranslateSql(
    sqlFilename = "RemoveCohortTempTables.sql",
    dbms = dbms,
    tempEmulationSchema = tempEmulationSchema,
    sampled = sampled
  )
  DatabaseConnector::executeSql(connection, sql, progressBar = FALSE, reportOverallTime = FALSE)

  .castAndromedaTable(data, "cohorts", .cohortsColumns, nullable = "daysToDeath")
  .castAndromedaTable(data, "treatments", .treatmentsColumns)
  .castAndromedaTable(data, "outcomes", .outcomesColumns)
  .castAndromedaTable(data, "covariates", .covariatesColumns)
  .castAndromedaTable(data, "covariateRef", .covariateRefColumns, nullable = setdiff(names(.covariateRefColumns), "covariateId"))
  .castAndromedaTable(data, "analysisRef", .analysisRefColumns, nullable = setdiff(names(.analysisRefColumns), "analysisId"))

  populationCounts <- .andromedaSql(data, "SELECT COUNT(DISTINCT personSeqId) AS persons, COUNT(*) AS exposures FROM cohorts;")
  if (sampled) {
    attrition <- dplyr::bind_rows(attrition, tibble(
      description = "Random sample",
      targetPersons = as.numeric(populationCounts$persons),
      comparatorPersons = 0,
      targetExposures = as.numeric(populationCounts$exposures),
      comparatorExposures = 0
    ))
  }
  metaData <- list(
    targetId = as.numeric(targetId),
    treatmentIds = treatmentIds,
    outcomeIds = as.numeric(outcomeIds),
    studyStartDate = args$studyStartDate,
    studyEndDate = args$studyEndDate,
    priorUseLookback = args$priorUseLookback,
    populationSize = as.numeric(populationCounts$exposures),
    attrition = attrition
  )
  covariateMetaData <- attr(data, "metaData")
  attr(data, "metaData") <- c(covariateMetaData[setdiff(names(covariateMetaData), names(metaData))], metaData)
  class(data) <- "CloneCensorWeightData"
  attr(class(data), "package") <- "CloneCensorWeight"
  return(data)
}

# Reads a parameterized SQL file from inst/sql/sql_server, renders it, and translates it to `dbms`.
# Uses system.file() from this package's namespace, so it also works under devtools::load_all().
.loadRenderTranslateSql <- function(sqlFilename, dbms, tempEmulationSchema, ...) {
  pathToSql <- system.file("sql", "sql_server", sqlFilename, package = "CloneCensorWeight", mustWork = TRUE)
  sql <- SqlRender::readSql(pathToSql)
  sql <- SqlRender::render(sql, ...)
  return(SqlRender::translate(sql, targetDialect = dbms, tempEmulationSchema = tempEmulationSchema))
}

# Cohort-based covariates (FeatureExtraction::createCohortBasedCovariateSettings) without their own
# cohort table use the exposure cohort table, as in CohortMethod.
.setCovariateCohortTable <- function(covariateSettings, exposureDatabaseSchema, exposureTable) {
  if (inherits(covariateSettings, "covariateSettings")) {
    covariateSettings <- list(covariateSettings)
  }
  for (i in seq_along(covariateSettings)) {
    if ("covariateCohorts" %in% names(covariateSettings[[i]]) && is.null(covariateSettings[[i]]$covariateCohortTable)) {
      covariateSettings[[i]]$covariateCohortDatabaseSchema <- exposureDatabaseSchema
      covariateSettings[[i]]$covariateCohortTable <- exposureTable
    }
  }
  return(covariateSettings)
}

# Casts the columns of an Andromeda table to the types in `columnTypes` (as in the foundation column
# specifications), keeping any other columns, and checks the result. An empty table is replaced by an
# empty table with the right columns, because the column types of an empty query result are unknown.
.castAndromedaTable <- function(andromeda, tableName, columnTypes, nullable = character(0)) {
  rowCount <- .andromedaSql(andromeda, sprintf("SELECT COUNT(*) AS n FROM %s;", tableName))$n
  if (rowCount == 0) {
    andromeda[[tableName]] <- .createEmptyTable(columnTypes)
  } else {
    sqlTypes <- c(integer = "INTEGER", numeric = "DOUBLE", character = "VARCHAR", Date = "DATE")
    otherColumns <- setdiff(colnames(andromeda[[tableName]]), names(columnTypes))
    columns <- c(
      sprintf("CAST(%s AS %s) AS %s", names(columnTypes), sqlTypes[columnTypes], names(columnTypes)),
      otherColumns
    )
    .andromedaExecute(andromeda, sprintf(
      "CREATE TABLE %s_typed AS SELECT %s FROM %s;",
      tableName,
      paste(columns, collapse = ", "),
      tableName
    ))
    .andromedaExecute(andromeda, sprintf("DROP TABLE %s;", tableName))
    .andromedaExecute(andromeda, sprintf("ALTER TABLE %s_typed RENAME TO %s;", tableName, tableName))
  }
  .checkColumns(andromeda[[tableName]], columnTypes, tableName, nullable = nullable)
  invisible(NULL)
}

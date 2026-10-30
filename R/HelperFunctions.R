# Internal utility functions — none of these are exported.

`%||%` <- function(left, right) if (!is.null(left)) left else right

# Returns the path to inst/sql/sql_server/, working for both installed and
# development layouts.
.sqlFolder <- function() {
  folder <- system.file("sql", "sql_server", package = "CloneCensorWeight")
  if (!nzchar(folder)) folder <- file.path("inst", "sql", "sql_server")
  folder
}

# Renders, translates, and queries a package SQL file, returning a data frame.
# `...` passes named render parameters to DatabaseConnector.
.renderTranslateQuerySql <- function(connection, sqlFile, ...) {
  sql <- readLines(file.path(.sqlFolder(), sqlFile), warn = FALSE)
  sql <- paste(sql, collapse = "\n")
  DatabaseConnector::renderTranslateQuerySql(
    connection           = connection,
    sql                  = sql,
    snakeCaseToCamelCase = TRUE,
    ...
  )
}

# Renders, translates, and executes a package SQL file (no rows returned).
# `...` passes named render parameters to DatabaseConnector.
.renderTranslateExecuteSql <- function(connection, sqlFile, ...) {
  sql <- readLines(file.path(.sqlFolder(), sqlFile), warn = FALSE)
  sql <- paste(sql, collapse = "\n")
  DatabaseConnector::renderTranslateExecuteSql(
    connection         = connection,
    sql                = sql,
    progressBar        = FALSE,
    reportOverallTime  = FALSE,
    ...
  )
}

# Uploads cohorts to the #ccw_cohorts temp table (snake_case columns).
# Expected input columns (camelCase): subjectId, cohortStartDate, cohortEndDate,
# initiationDate.  cohortDefinitionId is added as 1L if absent.
# The resulting table is FeatureExtraction-compatible (has cohort_definition_id).
.uploadTempCohortTable <- function(connection, cohorts) {
  if (!"cohortDefinitionId" %in% names(cohorts)) {
    cohorts$cohortDefinitionId <- 1L
  }
  DatabaseConnector::insertTable(
    connection           = connection,
    tableName            = "#ccw_cohorts",
    data                 = cohorts,
    dropTableIfExists    = TRUE,
    createTable          = TRUE,
    tempTable            = TRUE,
    camelCaseToSnakeCase = TRUE,
    progressBar          = FALSE
  )
  invisible("#ccw_cohorts")
}

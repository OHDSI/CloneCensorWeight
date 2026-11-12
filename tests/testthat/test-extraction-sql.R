# Tests of the extraction SQL (render and translate for every supported platform) and of the argument
# checks of getDbCloneCensorWeightData(). No database needed.

extractionDialects <- c("postgresql", "sql server", "oracle", "bigquery", "spark", "snowflake", "duckdb", "sqlite", "redshift")

# Every SQL file with parameter sets that cover all of its conditional blocks.
extractionSqlParameters <- function() {
  cohortParameters <- expand.grid(
    action = c("CREATE", "COUNT"),
    first_only = c(TRUE, FALSE),
    study_start_date = c("", "20100101"),
    study_end_date = c("", "20191231"),
    remove_subjects_with_prior_treatment = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  cohortParameters <- lapply(seq_len(nrow(cohortParameters)), function(i) {
    c(as.list(cohortParameters[i, ]), list(
      cdm_database_schema = "cdm",
      target_database_schema = "results",
      target_table = "eligible_cohort",
      exposure_database_schema = "results",
      exposure_table = "cohort",
      target_id = 1,
      treatment_ids = c(2, 3),
      washout_period = 365,
      prior_use_lookback = 365
    ))
  })
  cohortTables <- c("#ccw_cohort", "#ccw_sample")
  return(list(
    CreateOrCountCohorts.sql = cohortParameters,
    SampleCohorts.sql = list(list(max_cohort_size = 1000)),
    GetCohorts.sql = lapply(cohortTables, function(table) list(cdm_database_schema = "cdm", cohort_table = table)),
    GetTreatments.sql = lapply(cohortTables, function(table) {
      list(exposure_database_schema = "results", exposure_table = "cohort", treatment_ids = c(2, 3), prior_use_lookback = 365, cohort_table = table)
    }),
    GetOutcomes.sql = lapply(cohortTables, function(table) {
      list(outcome_database_schema = "results", outcome_table = "cohort", outcome_ids = c(4, 5), cohort_table = table)
    }),
    RemoveCohortTempTables.sql = list(list(sampled = TRUE), list(sampled = FALSE))
  ))
}

test_that("Every SQL file is covered by the render tests", {
  sqlFolder <- system.file("sql", "sql_server", package = "CloneCensorWeight")
  expect_setequal(list.files(sqlFolder, pattern = "\\.sql$"), names(extractionSqlParameters()))
})

test_that("Every SQL file renders and translates for every supported platform", {
  skip_if_not_installed("SqlRender")
  parameters <- extractionSqlParameters()
  for (sqlFilename in names(parameters)) {
    for (parameterSet in parameters[[sqlFilename]]) {
      # Rendering leaves no parameters, and warns about none
      rendered <- expect_no_warning(do.call(
        CloneCensorWeight:::.loadRenderTranslateSql,
        c(list(sqlFilename = sqlFilename, dbms = "sql server", tempEmulationSchema = NULL), parameterSet)
      ))
      code <- gsub("(?s)/\\*.*?\\*/", "", rendered, perl = TRUE)
      expect_false(grepl("@", code), info = sqlFilename)
      for (dbms in extractionDialects) {
        translated <- expect_no_error(do.call(
          CloneCensorWeight:::.loadRenderTranslateSql,
          c(list(sqlFilename = sqlFilename, dbms = dbms, tempEmulationSchema = "scratch"), parameterSet)
        ))
        expect_true(nchar(translated) > 0, info = paste(sqlFilename, dbms))
        if (!dbms %in% c("sql server", "redshift")) {
          # Temp tables are emulated or translated everywhere else (Redshift has # temp tables)
          expect_false(grepl("#ccw_", gsub("(?s)/\\*.*?\\*/", "", translated, perl = TRUE)), info = paste(sqlFilename, dbms))
        }
      }
    }
  }
})

test_that("The cohort SQL includes each step only when requested", {
  skip_if_not_installed("SqlRender")
  render <- function(...) {
    parameters <- extractionSqlParameters()$CreateOrCountCohorts.sql[[1]]
    parameters[names(list(...))] <- list(...)
    do.call(CloneCensorWeight:::.loadRenderTranslateSql, c(
      list(sqlFilename = "CreateOrCountCohorts.sql", dbms = "sql server", tempEmulationSchema = NULL),
      parameters
    ))
  }
  sql <- render(action = "COUNT", first_only = TRUE, study_start_date = "20100101", remove_subjects_with_prior_treatment = TRUE)
  expect_true(grepl("First cohort entry only", sql))
  expect_true(grepl("365 days of prior observation", sql))
  expect_true(grepl("Restrict to study period", sql))
  expect_true(grepl("No prior treatment use in the 365 days before index", sql))
  sql <- render(action = "COUNT", first_only = FALSE, study_start_date = "", study_end_date = "", remove_subjects_with_prior_treatment = FALSE)
  expect_false(grepl("First cohort entry only", sql))
  expect_true(grepl("365 days of prior observation", sql))
  expect_false(grepl("Restrict to study period", sql))
  expect_false(grepl("No prior treatment use", sql))
  sql <- render(action = "CREATE", remove_subjects_with_prior_treatment = TRUE)
  expect_true(grepl("INTO #ccw_cohort\\s+FROM no_prior_treatment", sql))
})

test_that("getDbCloneCensorWeightData checks its arguments", {
  skip_if_not_installed("DatabaseConnector")
  skip_if_not_installed("SqlRender")
  args <- createGetDbCloneCensorWeightDataArgs()
  connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "sqlite", server = tempfile())
  call <- function(...) {
    arguments <- list(
      connectionDetails = connectionDetails,
      cdmDatabaseSchema = "main",
      targetId = 1,
      treatmentId = 2,
      comparatorTreatmentId = 3,
      outcomeIds = 4,
      getDbCloneCensorWeightDataArgs = args
    )
    arguments[names(list(...))] <- list(...)
    do.call(getDbCloneCensorWeightData, arguments)
  }
  expect_error(call(connectionDetails = NULL), "either connectionDetails or connection")
  expect_error(call(connection = "x"), "connection")
  expect_error(call(connectionDetails = "x"), "connectionDetails")
  expect_error(call(comparatorTreatmentId = 2), "must be different")
  expect_error(call(targetId = 1.5), "integers")
  expect_error(call(outcomeIds = numeric(0)), "outcomeIds")
  expect_error(call(targetId = c(1, 2)), "targetId")
  expect_error(call(cdmDatabaseSchema = 1), "cdmDatabaseSchema")
  expect_error(call(getDbCloneCensorWeightDataArgs = list()), "getDbCloneCensorWeightDataArgs")
})

test_that("Cohort-based covariate settings without a cohort table use the exposure table", {
  settings <- list(
    FeatureExtraction::createCovariateSettings(useDemographicsGender = TRUE),
    FeatureExtraction::createCohortBasedCovariateSettings(
      analysisId = 999,
      covariateCohorts = data.frame(cohortId = 1, cohortName = "A"),
      valueType = "binary"
    ),
    FeatureExtraction::createCohortBasedCovariateSettings(
      analysisId = 998,
      covariateCohorts = data.frame(cohortId = 1, cohortName = "A"),
      covariateCohortDatabaseSchema = "other",
      covariateCohortTable = "other_cohort",
      valueType = "binary"
    )
  )
  result <- CloneCensorWeight:::.setCovariateCohortTable(settings, "results", "cohort")
  expect_null(result[[1]]$covariateCohortTable)
  expect_equal(result[[2]]$covariateCohortDatabaseSchema, "results")
  expect_equal(result[[2]]$covariateCohortTable, "cohort")
  expect_equal(result[[3]]$covariateCohortTable, "other_cohort")
  single <- CloneCensorWeight:::.setCovariateCohortTable(settings[[2]], "results", "cohort")
  expect_length(single, 1)
  expect_equal(single[[1]]$covariateCohortTable, "cohort")
  expect_identical(attributes(single[[1]]), attributes(result[[2]]))
})

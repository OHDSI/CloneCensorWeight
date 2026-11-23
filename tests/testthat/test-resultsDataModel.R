library(testthat)

# The results data model (plan section 3.5, Phase 2 item 3): the specification, the frozen DDL, the
# migrator, exportToCsv() with minimum cell counts, and the round trip export -> createResultsDataModel
# -> uploadResults into SQLite, on the shipped fixture and on a fresh Eunomia run.

specifications <- getResultsDataModelSpecifications()
definitionTables <- c("ccw_analysis", "ccw_strategy", "ccw_target_treatment", "ccw_target_treatment_outcome")

quietly <- function(expr) {
  value <- NULL
  utils::capture.output(value <- suppressMessages(expr))
  return(value)
}

newSqlite <- function() {
  # The file name is fixed first: DatabaseConnector evaluates the arguments at every connect.
  sqliteFile <- tempfile(fileext = ".sqlite")
  return(DatabaseConnector::createConnectionDetails(dbms = "sqlite", server = sqliteFile))
}

tableColumns <- function(connection, tableName) {
  return(tolower(colnames(DatabaseConnector::querySql(connection, sprintf("SELECT * FROM main.%s WHERE 1 = 0;", tableName)))))
}

countRows <- function(connection, tableName) {
  return(DatabaseConnector::querySql(connection, sprintf("SELECT COUNT(*) FROM main.%s;", tableName))[1, 1])
}

csvColumns <- function(fileName) {
  return(colnames(readr::read_csv(fileName, n_max = 0, col_types = readr::cols())))
}

test_that("The specification follows the rules of Strategus and its consumers", {
  expect_setequal(
    colnames(specifications),
    c("tableName", "columnName", "dataType", "isRequired", "primaryKey", "minCellCount", "deprecated", "description")
  )
  expect_true(all(startsWith(specifications$tableName, "ccw_")))
  expect_true(all(grepl("^[a-z][a-z0-9_]*$", specifications$tableName)))
  expect_true(all(grepl("^[a-z][a-z0-9_]*$", specifications$columnName)))
  expect_false(anyDuplicated(paste(specifications$tableName, specifications$columnName)) > 0)
  expect_true(all(sub("\\(.*\\)", "", specifications$dataType) %in% c("int", "bigint", "float", "varchar")))
  for (column in c("isRequired", "primaryKey", "minCellCount", "deprecated")) {
    expect_true(all(specifications[[column]] %in% c("Yes", "No")))
  }
  expect_true(all(specifications$isRequired[specifications$primaryKey == "Yes"] == "Yes"))
  expect_true(all(specifications$dataType[specifications$minCellCount == "Yes"] %in% c("int", "float")))
  for (tableName in unique(specifications$tableName)) {
    table <- specifications[specifications$tableName == tableName, ]
    keys <- table$columnName[table$primaryKey == "Yes"]
    expect_true(length(keys) > 0, info = tableName)
    # database_id is in the primary key of every site-specific table, and absent from the definitions:
    if (tableName %in% definitionTables) {
      expect_false("database_id" %in% table$columnName, info = tableName)
    } else {
      expect_true("database_id" %in% keys, info = tableName)
    }
    # Strategy-aware keys:
    if (any(grepl("strategy_id$", table$columnName)) && tableName != "ccw_strategy") {
      strategyColumns <- table$columnName[grepl("strategy_id$", table$columnName)]
      expect_true(all(strategyColumns %in% keys), info = tableName)
      expect_true(
        setequal(strategyColumns, "strategy_id") || setequal(strategyColumns, c("target_strategy_id", "comparator_strategy_id")),
        info = tableName
      )
    }
  }
  expect_true(all(c("target_strategy_id", "comparator_strategy_id", "estimand") %in%
    specifications$columnName[specifications$tableName == "ccw_result" & specifications$primaryKey == "Yes"]))
  # The columns EvidenceSynthesis and the viewers read, with CohortMethod's names:
  expect_true(all(c("log_rr", "se_log_rr", "ci_95_lb", "ci_95_ub", "p", "one_sided_p", "calibrated_log_rr", "calibrated_se_log_rr") %in%
    specifications$columnName[specifications$tableName == "ccw_result"]))
  expect_true(all(c("mdrr", "unblind", "unblind_for_evidence_synthesis") %in%
    specifications$columnName[specifications$tableName == "ccw_diagnostics_summary"]))
})

test_that("The target-treatment ID is a hash of the three cohort IDs", {
  ids <- .createTargetTreatmentId(c(4, 4, 4), c(1, 1, 2), c(2, NA, 1))
  expect_type(ids, "integer")
  expect_equal(length(unique(ids)), 3)
  expected <- as.integer(as.numeric(paste0("0x", digest::digest("4 1 2", algo = "murmur32", serialize = FALSE))) - 2^31)
  expect_equal(ids[1], expected)
  expect_equal(.createTargetTreatmentId(4, 1, NA), ids[2])
})

test_that(".enforceMinCellValue censors counts below the minimum, as CohortMethod", {
  data <- dplyr::tibble(x = c(0, 1, 4, 5, 100, NA))
  censored <- .enforceMinCellValue(data, "x", 5, silent = TRUE)
  expect_equal(censored$x, c(0, -5, -5, 5, 100, NA))
  censored <- .enforceMinCellValue(data, "x", c(1, 2, 5, 6, 1, 1), silent = TRUE)
  expect_equal(censored$x, c(0, -2, -5, -6, 100, NA))
  expect_message(.enforceMinCellValue(data, "x", 5), "censoring 2 values")
})

test_that("Censoring-model densities are keyed by horizon", {
  # An initiateWithin strategy with graceStart > 0 has two hazard models, at two horizons:
  positivity <- dplyr::tibble(strategyId = 1L, model = c("hazard", "hazard"), day = c(6L, 30L))
  attr(positivity, "predictions") <- dplyr::tibble(
    strategyId = 1L, model = "hazard", day = rep(c(6L, 30L), each = 4), probability = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9)
  )
  densities <- .censoringModelDensities(positivity)
  keys <- specifications$columnName[specifications$tableName == "ccw_censoring_model_dist" & specifications$primaryKey == "Yes"]
  keys <- SqlRender::snakeCaseToCamelCase(setdiff(keys, c("database_id", "analysis_id", "target_treatment_id", "outcome_id")))
  expect_true("horizon" %in% keys)
  expect_true(anyDuplicated(densities[, setdiff(keys, "horizon")]) > 0)
  expect_false(anyDuplicated(densities[, keys]) > 0)
  # Each horizon's density uses its own predictions:
  expect_false(isTRUE(all.equal(densities$density[densities$horizon == 6], densities$density[densities$horizon == 30])))
})

test_that("Likelihood profiles are exported as CohortMethod's, without contrasts null by design", {
  outputFolder <- tempfile()
  dir.create(outputFolder)
  on.exit(unlink(outputFolder, recursive = TRUE))
  profile <- data.frame(point = c(-1, 0, 1), value = c(-3, -1, -2), derivative = c(2, 0, -1))
  estimates <- function(se) tibble(estimand = "hr", seLogEstimate = se)
  saveRDS(list(modelType = "cox", estimates = estimates(0.2)), file.path(outputFolder, "none.rds"))
  saveRDS(list(modelType = "cox", estimates = estimates(0), logLikelihoodProfile = profile), file.path(outputFolder, "null.rds"))
  # Cyclops' robust variance leaves rounding noise in the standard error of a contrast null by design:
  saveRDS(list(modelType = "cox", estimates = estimates(2.4e-7), logLikelihoodProfile = profile), file.path(outputFolder, "noise.rds"))
  saveRDS(list(modelType = "cox", estimates = estimates(0.2), logLikelihoodProfile = profile), file.path(outputFolder, "profile.rds"))
  reference <- tibble(
    analysisId = 1:5, targetTreatmentId = 5, outcomeId = 3, targetStrategyId = 1L, comparatorStrategyId = 2L,
    outcomeModelFile = c("none.rds", "null.rds", "noise.rds", "", "profile.rds")
  )
  # (a) no profile, (b) a contrast null by design, with a standard error of 0 or of rounding noise, (c) no
  # outcome model:
  fileName <- quietly(.exportLikelihoodProfiles(outputFolder, reference[1:4, ], outputFolder, "db"))
  expect_equal(csvColumns(file.path(outputFolder, fileName)), specifications$columnName[specifications$tableName == "ccw_likelihood_profile"])
  expect_equal(nrow(readr::read_csv(file.path(outputFolder, fileName), col_types = readr::cols())), 0)
  # A profile, with the log likelihood normalized to a maximum of 0 and the gradient as is:
  quietly(.exportLikelihoodProfiles(outputFolder, reference, outputFolder, "db"))
  exported <- readr::read_csv(file.path(outputFolder, fileName), col_types = readr::cols())
  expect_equal(exported$analysis_id, rep(5, 3))
  expect_equal(exported$estimand, rep("hr", 3))
  expect_equal(exported$log_rr, profile$point)
  expect_equal(exported$log_likelihood, c(-2, 0, -1))
  expect_equal(exported$gradient, profile$derivative)
})

test_that("createResultsDataModel creates the specified tables, and migrations run from an empty schema", {
  skip_if_not_installed("RSQLite")
  connectionDetails <- newSqlite()
  quietly(createResultsDataModel(connectionDetails, databaseSchema = "main"))
  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection))
  for (tableName in unique(specifications$tableName)) {
    expect_equal(
      tableColumns(connection, tableName),
      specifications$columnName[specifications$tableName == tableName],
      info = tableName
    )
  }
  # The migrations table exists, and this version has no migrations to apply:
  expect_equal(tableColumns(connection, "ccw_migration"), c("migration_file", "migration_order"))
  expect_equal(countRows(connection, "ccw_migration"), 0)
  migrator <- quietly(getDataMigrator(connectionDetails, databaseSchema = "main"))
  expect_s3_class(migrator, "DataMigrationManager")
  expect_equal(nrow(migrator$getStatus()), 0)
  expect_true(migrator$check())
  migrator$closeConnection()
  expect_silent(quietly(migrateDataModel(connectionDetails, databaseSchema = "main")))
  expect_equal(countRows(connection, "ccw_migration"), 0)

  # Migrating an empty schema does not fail:
  emptyConnectionDetails <- newSqlite()
  expect_silent(quietly(migrateDataModel(emptyConnectionDetails, databaseSchema = "main")))

  # A table prefix, and the SQLite schema check:
  prefixed <- newSqlite()
  quietly(createResultsDataModel(prefixed, databaseSchema = "main", tablePrefix = "pre_"))
  prefixedConnection <- DatabaseConnector::connect(prefixed)
  on.exit(DatabaseConnector::disconnect(prefixedConnection), add = TRUE)
  expect_true(DatabaseConnector::existsTable(prefixedConnection, "main", "pre_ccw_result"))
  expect_error(createResultsDataModel(newSqlite(), databaseSchema = "results"), "Invalid schema for sqlite")
})

test_that("The shipped Results_Eunomia.zip has the specified columns and uploads", {
  skip_if_not_installed("RSQLite")
  skip_if_not_installed("zip")
  zipFile <- system.file("Results_Eunomia.zip", package = "CloneCensorWeight")
  expect_true(file.exists(zipFile))
  unzipFolder <- tempfile()
  on.exit(unlink(unzipFolder, recursive = TRUE))
  zip::unzip(zipFile, exdir = unzipFolder)
  expect_setequal(list.files(unzipFolder), paste0(unique(specifications$tableName), ".csv"))
  for (tableName in unique(specifications$tableName)) {
    expect_equal(
      csvColumns(file.path(unzipFolder, paste0(tableName, ".csv"))),
      specifications$columnName[specifications$tableName == tableName],
      info = tableName
    )
  }
  # The ratio estimands carry the EASE of their null distribution; the summary keeps the largest per analysis:
  results <- readr::read_csv(file.path(unzipFolder, "ccw_result.csv"), col_types = readr::cols())
  diagnostics <- readr::read_csv(file.path(unzipFolder, "ccw_diagnostics_summary.csv"), col_types = readr::cols())
  expect_true(all(!is.na(results$ease[results$estimand %in% c("hr", "rr_365")])))
  expect_true(all(is.na(results$ease[!results$estimand %in% c("hr", "rr_365")])))
  for (analysisId in unique(results$analysis_id)) {
    expect_equal(diagnostics$ease[diagnostics$analysis_id == analysisId],
      rep(max(results$ease[results$analysis_id == analysisId], na.rm = TRUE), sum(diagnostics$analysis_id == analysisId)),
      tolerance = 1e-6, info = analysisId
    )
  }

  connectionDetails <- newSqlite()
  quietly(createResultsDataModel(connectionDetails, databaseSchema = "main"))
  quietly(uploadResults(connectionDetails, schema = "main", zipFileName = zipFile, purgeSiteDataBeforeUploading = FALSE))
  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection), add = TRUE)
  for (tableName in unique(specifications$tableName)) {
    expect_true(countRows(connection, tableName) > 0, info = tableName)
  }
  expect_equal(
    DatabaseConnector::querySql(connection, "SELECT DISTINCT database_id FROM main.ccw_result;")[, 1],
    "Eunomia"
  )
})

# Eunomia ----------------------------------------------------------------------------------------

skip_on_cran()
skip_if_not_installed("Eunomia")
skip_if_not_installed("RSQLite")
skip_if_not_installed("zip")

databaseFile <- tempfile(fileext = ".sqlite")
eunomiaConnectionDetails <- tryCatch(
  suppressMessages(Eunomia::getEunomiaConnectionDetails(databaseFile = databaseFile)),
  error = function(e) NULL
)
if (is.null(eunomiaConnectionDetails)) {
  skip("Eunomia GiBleed data could not be loaded")
}
utils::capture.output(suppressMessages(Eunomia::createCohorts(eunomiaConnectionDetails)))
local({
  connection <- DatabaseConnector::connect(eunomiaConnectionDetails)
  on.exit(DatabaseConnector::disconnect(connection))
  # A negative-control outcome cohort: every occurrence of a common condition.
  DatabaseConnector::executeSql(connection, "
    INSERT INTO main.cohort (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
    SELECT 101, person_id, condition_start_date, condition_start_date
    FROM main.condition_occurrence
    WHERE condition_concept_id = 4112343;", progressBar = FALSE, reportOverallTime = FALSE)
})

# Strategy IDs other than 1 and 2, to check they are exported (not the positions of the strategies):
exportAnalysis <- createCcwAnalysis(
  analysisId = 1,
  description = "Initiate vs not within 30 days",
  getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(
    covariateSettings = FeatureExtraction::createCovariateSettings(
      useDemographicsGender = TRUE,
      useDemographicsAgeGroup = TRUE,
      endDays = -1
    )
  ),
  createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = TRUE,
    riskWindowEnd = 365,
    endAnchor = "cohort start"
  ),
  createClonesArgs = createCreateClonesArgs(strategies = list(
    createTreatmentStrategy(7, "Celecoxib within 30 days", "initiateWithin", graceEnd = 30),
    createTreatmentStrategy(8, "No celecoxib within 30 days", "notInitiateWithin", graceEnd = 30)
  )),
  fitCensoringModelsArgs = createFitCensoringModelsArgs(
    prior = Cyclops::createPrior("laplace", exclude = c(0), variance = 1)
  ),
  computeWeightsArgs = createComputeWeightsArgs(),
  computeCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(),
  fitOutcomeModelArgs = createFitOutcomeModelArgs(horizons = 365)
)
exportSpecifications <- createCcwAnalysesSpecifications(
  ccwAnalysisList = list(exportAnalysis),
  targetTreatmentOutcomesList = list(createTargetTreatmentOutcomes(4, 1, 2, outcomes = list(
    CohortMethod::createOutcome(3),
    CohortMethod::createOutcome(101, outcomeOfInterest = FALSE, trueEffectSize = 1)
  )))
)
outputFolder <- tempfile()
quietly(runCcwAnalyses(
  connectionDetails = eunomiaConnectionDetails,
  cdmDatabaseSchema = "main",
  outputFolder = outputFolder,
  ccwAnalysesSpecifications = exportSpecifications
))
exportFolder <- file.path(outputFolder, "export")
zipFile <- quietly(exportToCsv(outputFolder, exportFolder = exportFolder, databaseId = "Eunomia", minCellCount = 5))

readExported <- function(folder, tableName) {
  return(readr::read_csv(file.path(folder, paste0(tableName, ".csv")), col_types = readr::cols()))
}

test_that("exportToCsv writes every table with exactly the specified columns", {
  expect_equal(zipFile, file.path(normalizePath(exportFolder), "Results_Eunomia.zip"))
  expect_true(file.exists(zipFile))
  expect_setequal(zip::zip_list(zipFile)$filename, paste0(unique(specifications$tableName), ".csv"))
  for (tableName in unique(specifications$tableName)) {
    expect_equal(
      csvColumns(file.path(exportFolder, paste0(tableName, ".csv"))),
      specifications$columnName[specifications$tableName == tableName],
      info = tableName
    )
    data <- readExported(exportFolder, tableName)
    expect_true(nrow(data) > 0, info = tableName)
    # Primary keys are unique and complete:
    keys <- specifications$columnName[specifications$tableName == tableName & specifications$primaryKey == "Yes"]
    expect_false(anyDuplicated(data[, keys]) > 0, info = tableName)
    expect_false(anyNA(data[, keys]), info = tableName)
  }

  # Strategy IDs, the target-treatment ID, and the results against getResultsSummary():
  targetTreatmentId <- .createTargetTreatmentId(4, 1, 2)
  expect_equal(readExported(exportFolder, "ccw_target_treatment")$target_treatment_id, targetTreatmentId)
  strategies <- readExported(exportFolder, "ccw_strategy")
  expect_equal(strategies$strategy_id, c(7, 8))
  expect_equal(strategies$type, c("initiateWithin", "notInitiateWithin"))
  expect_equal(strategies$grace_end, c(30, 30))
  results <- readExported(exportFolder, "ccw_result")
  expect_equal(unique(results$target_strategy_id), 7)
  expect_equal(unique(results$comparator_strategy_id), 8)
  expect_equal(unique(results$target_treatment_id), targetTreatmentId)
  summary <- getResultsSummary(outputFolder)
  expect_equal(results$estimate, summary$estimate, tolerance = 1e-10)
  expect_equal(results$se_log_rr, summary$seLogRr, tolerance = 1e-10)
  expect_equal(results$se, summary$se, tolerance = 1e-10)
  expect_equal(as.numeric(results$ease), summary$ease, tolerance = 1e-10)
  # A missing standard error is NA, never 0:
  expect_true(all(is.na(results$se_log_rr[!results$estimand %in% c("hr", "rr_365")])))
  diagnostics <- readExported(exportFolder, "ccw_diagnostics_summary")
  expect_equal(diagnostics$mdrr, getDiagnosticsSummary(outputFolder)$mdrr, tolerance = 1e-10)
  expect_true(all(diagnostics$unblind %in% c(0, 1)))

  # The tables of the outcome of interest only, per strategy:
  for (tableName in c("ccw_attrition", "ccw_follow_up_dist", "ccw_censoring_model", "ccw_censoring_model_dist", "ccw_weight_dist", "ccw_cumulative_incidence")) {
    data <- readExported(exportFolder, tableName)
    expect_equal(unique(data$outcome_id), 3, info = tableName)
    expect_setequal(unique(data$strategy_id), c(7, 8))
  }
  attrition <- readExported(exportFolder, "ccw_attrition")
  expect_equal(tail(attrition$description[attrition$strategy_id == 7], 2), c("Clones created", "Clones with follow-up"))
  incidence <- readExported(exportFolder, "ccw_cumulative_incidence")
  outcomeModel <- readRDS(file.path(outputFolder, getFileReference(outputFolder)$outcomeModelFile[1]))
  expect_equal(incidence$risk, outcomeModel$curves$risk, tolerance = 1e-10)
  # Rows every 30 days from day 0 and at the horizon (365), through the horizon:
  gridDays <- c(seq(0, 360, by = 30), 365)
  for (strategyId in c(7, 8)) {
    days <- incidence$time_day[incidence$strategy_id == strategyId]
    expect_true(all(gridDays[gridDays > min(days)] %in% days))
    expect_true(max(days) <= 365)
  }
  censoringModel <- readExported(exportFolder, "ccw_censoring_model")
  expect_setequal(unique(censoringModel$model), c("binary", "hazard"))
  expect_true(all(censoringModel$coefficient != 0))
  balance <- readExported(exportFolder, "ccw_covariate_balance")
  expect_true(all(balance$covariate_id %in% readExported(exportFolder, "ccw_covariate")$covariate_id))
  expect_equal(unique(balance$target_strategy_id), 7)

  # The likelihood profiles of the hazard ratio, of every outcome (controls included), normalized to a
  # maximum of 0 per key, for keys of ccw_result:
  profiles <- readExported(exportFolder, "ccw_likelihood_profile")
  expect_equal(unique(profiles$estimand), "hr")
  expect_setequal(unique(profiles$outcome_id), c(3, 101))
  profileKeys <- c("database_id", "analysis_id", "target_treatment_id", "outcome_id", "target_strategy_id", "comparator_strategy_id", "estimand")
  maxima <- dplyr::summarise(dplyr::group_by(profiles, dplyr::across(dplyr::all_of(profileKeys))), maximum = max(.data$log_likelihood), .groups = "drop")
  expect_equal(maxima$maximum, rep(0, nrow(maxima)))
  expect_equal(nrow(dplyr::anti_join(maxima, results, by = profileKeys)), 0)
})

test_that("exportToCsv enforces the minimum cell count", {
  censoredFolder <- tempfile()
  on.exit(unlink(censoredFolder, recursive = TRUE))
  minCellCount <- 200
  quietly(exportToCsv(outputFolder, exportFolder = censoredFolder, databaseId = "Eunomia", minCellCount = minCellCount))
  uncensored <- readExported(exportFolder, "ccw_result")
  censored <- readExported(censoredFolder, "ccw_result")
  # Eunomia has between 1 and 199 comparator outcomes of GI bleed:
  expect_true(any(uncensored$comparator_outcomes > 0 & uncensored$comparator_outcomes < minCellCount))
  countColumns <- specifications[specifications$minCellCount == "Yes" & !grepl("_mean_|^mean_", specifications$columnName), ]
  censoredAny <- FALSE
  for (tableName in unique(countColumns$tableName)) {
    data <- readExported(censoredFolder, tableName)
    for (column in countColumns$columnName[countColumns$tableName == tableName]) {
      values <- data[[column]]
      ok <- is.na(values) | values == 0 | values >= minCellCount | values == -minCellCount
      expect_true(all(ok), info = paste(tableName, column))
      censoredAny <- censoredAny || any(values == -minCellCount, na.rm = TRUE)
    }
  }
  expect_true(censoredAny)
  expect_equal(
    censored$comparator_outcomes,
    ifelse(uncensored$comparator_outcomes > 0 & uncensored$comparator_outcomes < minCellCount, -minCellCount, uncensored$comparator_outcomes)
  )
  # Means of the covariate balance are censored at minCellCount over the number of clones, so a
  # censored mean is negative:
  balance <- readExported(censoredFolder, "ccw_covariate_balance")
  uncensoredBalance <- readExported(exportFolder, "ccw_covariate_balance")
  expect_true(sum(balance$target_mean_before < 0) > sum(uncensoredBalance$target_mean_before < 0))
  # Other columns are unchanged:
  expect_equal(censored$estimate, uncensored$estimate)
})

test_that("Export, createResultsDataModel and uploadResults round trip into SQLite", {
  connectionDetails <- newSqlite()
  quietly(createResultsDataModel(connectionDetails, databaseSchema = "main"))
  quietly(uploadResults(connectionDetails, schema = "main", zipFileName = zipFile, purgeSiteDataBeforeUploading = FALSE))
  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection))
  for (tableName in unique(specifications$tableName)) {
    expect_equal(countRows(connection, tableName), nrow(readExported(exportFolder, tableName)), info = tableName)
  }
  uploaded <- DatabaseConnector::querySql(connection, "
    SELECT estimand, estimate, se_log_rr, target_strategy_id, comparator_strategy_id, database_id
    FROM main.ccw_result
    WHERE outcome_id = 3;", snakeCaseToCamelCase = TRUE)
  summary <- getResultsSummary(outputFolder)
  summary <- summary[summary$outcomeId == 3, ]
  uploaded <- uploaded[match(summary$estimand, uploaded$estimand), ]
  expect_equal(uploaded$estimate, summary$estimate, tolerance = 1e-10)
  expect_equal(uploaded$seLogRr, summary$seLogRr, tolerance = 1e-10)
  expect_equal(unique(uploaded$targetStrategyId), 7)
  expect_equal(unique(uploaded$databaseId), "Eunomia")
  definition <- DatabaseConnector::querySql(connection, "SELECT definition FROM main.ccw_analysis;")[, 1]
  expect_equal(jsonlite::fromJSON(definition), jsonlite::fromJSON(exportAnalysis$toJson()))

  # Uploading the same site again replaces its rows:
  quietly(uploadResults(connectionDetails, schema = "main", zipFileName = zipFile, purgeSiteDataBeforeUploading = FALSE))
  for (tableName in unique(specifications$tableName)) {
    expect_equal(countRows(connection, tableName), nrow(readExported(exportFolder, tableName)), info = tableName)
  }
  # A second site adds its rows to the site-specific tables only:
  secondFolder <- tempfile()
  on.exit(unlink(secondFolder, recursive = TRUE), add = TRUE)
  secondZip <- quietly(exportToCsv(outputFolder, exportFolder = secondFolder, databaseId = "Second", minCellCount = 5))
  quietly(uploadResults(connectionDetails, schema = "main", zipFileName = secondZip, purgeSiteDataBeforeUploading = FALSE))
  for (tableName in unique(specifications$tableName)) {
    expected <- nrow(readExported(exportFolder, tableName)) * if (tableName %in% definitionTables) 1 else 2
    expect_equal(countRows(connection, tableName), expected, info = tableName)
  }
})

test_that("exportToCsv checks its arguments", {
  expect_error(exportToCsv(tempfile(), databaseId = "x"), "Directory")
  emptyFolder <- tempfile()
  dir.create(emptyFolder)
  on.exit(unlink(emptyFolder, recursive = TRUE))
  expect_error(exportToCsv(emptyFolder, databaseId = "x"), "ccwAnalysesSpecifications.rds")
  expect_error(exportToCsv(outputFolder, databaseId = "x", minCellCount = -1), "minCellCount")
})

unlink(outputFolder, recursive = TRUE)

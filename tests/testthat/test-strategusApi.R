library(testthat)

# The functions a Strategus module calls (plan section 3.8, module pattern A: the module is an adapted copy
# of Strategus's CohortMethodModule). The tests replay what the module's execute, createResultsDataModel
# and uploadResults do, with a job context built by hand, so they need no Strategus (plan section 3.8 (f)).
# The replay is a copy of those calls, not the module: the module file itself
# (extras/strategus/CloneCensorWeightModule.R) is checked by extras/strategus/RunModulePrototype.R, which
# runs it with Strategus's own job context, directly and under Strategus::execute.

# The CCW functions the module calls, and the arguments it passes to them:
moduleCalls <- list(
  createDefaultCcwMultiThreadingSettings = "maxCores",
  convertUntypedListToCcwAnalysesSpecifications = "untypedList",
  runCcwAnalyses = c(
    "connectionDetails", "cdmDatabaseSchema", "tempEmulationSchema", "exposureDatabaseSchema", "exposureTable",
    "outcomeDatabaseSchema", "outcomeTable", "outputFolder", "multiThreadingSettings", "ccwAnalysesSpecifications"
  ),
  exportToCsv = c("outputFolder", "exportFolder", "databaseId", "minCellCount", "maxCores"),
  createResultsDataModel = c("connectionDetails", "databaseSchema", "tablePrefix"),
  uploadResults = c("connectionDetails", "schema", "zipFileName", "purgeSiteDataBeforeUploading")
)

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

# The module specifications as CloneCensorWeightModule$createModuleSpecifications() makes them, and as
# Strategus hands them to the module after ParallelLogger::saveSettingsToJson() and loadSettingsFromJson().
createModuleSettings <- function(specifications) {
  return(list(ccwAnalysesSpecifications = specifications$toList()))
}

loadedFromJson <- function(settings) {
  fileName <- tempfile(fileext = ".json")
  on.exit(unlink(fileName))
  ParallelLogger::saveSettingsToJson(settings, fileName)
  return(ParallelLogger::loadSettingsFromJson(fileName))
}

# What CloneCensorWeightModule$execute() does once the job context is built.
replayExecute <- function(jobContext, connectionDetails) {
  executionSettings <- jobContext$moduleExecutionSettings
  multiThreadingSettings <- createDefaultCcwMultiThreadingSettings(executionSettings$maxCores)
  if (is.null(jobContext$settings$ccwAnalysesSpecifications)) {
    stop("The CloneCensorWeightModule specification is missing the required `ccwAnalysesSpecifications` setting.")
  }
  args <- jobContext$settings
  args$connectionDetails <- connectionDetails
  args$cdmDatabaseSchema <- executionSettings$cdmDatabaseSchema
  args$tempEmulationSchema <- executionSettings$tempEmulationSchema
  args$exposureDatabaseSchema <- executionSettings$workDatabaseSchema
  args$exposureTable <- executionSettings$cohortTableNames$cohortTable
  args$outcomeDatabaseSchema <- executionSettings$workDatabaseSchema
  args$outcomeTable <- executionSettings$cohortTableNames$cohortTable
  args$outputFolder <- executionSettings$workSubFolder
  args$multiThreadingSettings <- multiThreadingSettings
  args$ccwAnalysesSpecifications <- convertUntypedListToCcwAnalysesSpecifications(jobContext$settings$ccwAnalysesSpecifications)
  do.call(runCcwAnalyses, args)

  exportFolder <- executionSettings$resultsSubFolder
  exportToCsv(
    outputFolder = executionSettings$workSubFolder,
    exportFolder = exportFolder,
    databaseId = executionSettings$cdmDatabaseMetaData$databaseId,
    minCellCount = executionSettings$minCellCount,
    maxCores = executionSettings$maxCores
  )
  resultsDataModel <- getResultsDataModelSpecifications()
  resultsDataModel <- resultsDataModel[file.exists(file.path(exportFolder, paste0(resultsDataModel$tableName, ".csv"))), ]
  if (any(!startsWith(resultsDataModel$tableName, "ccw_"))) {
    stop("Table names do not have required prefix")
  }
  return(resultsDataModel)
}

# What CloneCensorWeightModule$uploadResults() does: the first zip file in the results folder.
replayUploadResults <- function(resultsSubFolder, resultsConnectionDetails) {
  zipFiles <- list.files(path = normalizePath(resultsSubFolder), pattern = "\\.zip$", full.names = TRUE)
  uploadResults(
    connectionDetails = resultsConnectionDetails,
    schema = "main",
    zipFileName = zipFiles[1],
    purgeSiteDataBeforeUploading = FALSE
  )
  return(zipFiles[1])
}

test_that("The functions the module calls are exported and take the module's arguments", {
  for (name in names(moduleCalls)) {
    expect_true(name %in% getNamespaceExports("CloneCensorWeight"), info = name)
    expect_true(all(moduleCalls[[name]] %in% names(formals(get(name)))), info = name)
  }
  # The thread counts the module's options strategus.CloneCensorWeightModule.* override:
  expect_true(all(c("getDbCcwDataThreads", "fitOutcomeModelThreads") %in% names(createDefaultCcwMultiThreadingSettings(4))))
})

test_that("The module settings are formals of runCcwAnalyses and convert in both stored forms", {
  specifications <- createCcwAnalysesSpecifications(
    ccwAnalysisList = list(createCcwAnalysis(
      analysisId = 1,
      getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(
        covariateSettings = list(
          FeatureExtraction::createCovariateSettings(useDemographicsGender = TRUE, endDays = -1),
          # Without a cohort table, which Strategus fills in (plan section 3.8 (b)):
          FeatureExtraction::createCohortBasedCovariateSettings(
            analysisId = 999,
            covariateCohorts = data.frame(cohortId = 1, cohortName = "Treatment"),
            valueType = "binary"
          )
        )
      ),
      createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 365),
      createClonesArgs = createCreateClonesArgs(strategies = list(
        createTreatmentStrategy(1, "Initiate within 30 days", "initiateWithin", graceEnd = 30),
        createTreatmentStrategy(2, "Do not initiate within 30 days", "notInitiateWithin", graceEnd = 30)
      )),
      fitCensoringModelsArgs = createFitCensoringModelsArgs(),
      computeWeightsArgs = createComputeWeightsArgs(),
      fitOutcomeModelArgs = createFitOutcomeModelArgs()
    )),
    targetTreatmentOutcomesList = list(createTargetTreatmentOutcomes(4, 1, 2, outcomes = list(CohortMethod::createOutcome(3))))
  )
  inSession <- createModuleSettings(specifications)
  # do.call(runCcwAnalyses, settings) passes every setting as an argument:
  expect_true(all(names(inSession) %in% names(formals(runCcwAnalyses))))
  # Built in the same session, the settings objects carry their classes as attr_class members, so
  # Strategus's .replaceCovariateSettings, which looks for class covariateSettings with attribute fun
  # "getDbCohortBasedCovariatesData", leaves the cohort table NULL, and CCW defaults it to the exposure
  # table (test-extraction-sql.R):
  sessionCovariateSettings <- inSession$ccwAnalysesSpecifications$ccwAnalysisList[[1]]$getDbCloneCensorWeightDataArgs$covariateSettings
  expect_equal(sessionCovariateSettings[[1]]$attr_class, "covariateSettings")
  expect_false(inherits(sessionCovariateSettings[[2]], "covariateSettings"))
  expect_null(sessionCovariateSettings[[2]]$covariateCohortTable)
  fromSession <- convertUntypedListToCcwAnalysesSpecifications(inSession$ccwAnalysesSpecifications)
  expect_equal(fromSession$toList(), specifications$toList())
  # Loaded from JSON, ParallelLogger restores the classes and the fun attribute, so the hook fills in
  # the cohort table:
  loaded <- loadedFromJson(inSession)
  jsonCovariateSettings <- loaded$ccwAnalysesSpecifications$ccwAnalysisList[[1]]$getDbCloneCensorWeightDataArgs$covariateSettings
  expect_s3_class(jsonCovariateSettings[[1]], "covariateSettings")
  expect_s3_class(jsonCovariateSettings[[2]], "covariateSettings")
  expect_equal(attr(jsonCovariateSettings[[2]], "fun"), "getDbCohortBasedCovariatesData")
  expect_null(jsonCovariateSettings[[2]]$covariateCohortTable)
  fromJson <- convertUntypedListToCcwAnalysesSpecifications(loaded$ccwAnalysesSpecifications)
  expect_equal(fromJson$toList(), specifications$toList())

  # The likely mistakes stop with a message that names them:
  expect_error(convertUntypedListToCcwAnalysesSpecifications(specifications), "pass its `\\$toList\\(\\)`")
  expect_error(convertUntypedListToCcwAnalysesSpecifications(NULL), "has no ccwAnalysisList or targetTreatmentOutcomesList")
  expect_error(
    convertUntypedListToCcwAnalysesSpecifications(list(ccwAnalysisList = inSession$ccwAnalysesSpecifications$ccwAnalysisList)),
    "has no targetTreatmentOutcomesList"
  )
})

# Eunomia ------------------------------------------------------------------------------------------------

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

moduleSpecifications <- createCcwAnalysesSpecifications(
  ccwAnalysisList = list(createCcwAnalysis(
    analysisId = 1,
    description = "Celecoxib within 30 days vs not",
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
      createTreatmentStrategy(1, "Celecoxib within 30 days", "initiateWithin", graceEnd = 30),
      createTreatmentStrategy(2, "No celecoxib within 30 days", "notInitiateWithin", graceEnd = 30)
    )),
    fitCensoringModelsArgs = createFitCensoringModelsArgs(
      prior = Cyclops::createPrior("laplace", exclude = c(0), variance = 1)
    ),
    computeWeightsArgs = createComputeWeightsArgs(),
    computeCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(),
    fitOutcomeModelArgs = createFitOutcomeModelArgs(horizons = 365)
  )),
  targetTreatmentOutcomesList = list(createTargetTreatmentOutcomes(4, 1, 2, outcomes = list(CohortMethod::createOutcome(3))))
)

# The job context of a CdmExecutionSettings, with the fields the module reads, as StrategusModule's
# .createJobContext() builds it. databaseId is an integer hash under Strategus::execute().
studyFolder <- tempfile("strategusApi")
handBuiltContext <- list(
  settings = createModuleSettings(moduleSpecifications),
  moduleExecutionSettings = list(
    workDatabaseSchema = "main",
    cdmDatabaseSchema = "main",
    cohortTableNames = list(cohortTable = "cohort"),
    tempEmulationSchema = NULL,
    workFolder = file.path(studyFolder, "work_folder"),
    resultsFolder = file.path(studyFolder, "results_folder"),
    minCellCount = 5,
    maxCores = 1,
    workSubFolder = file.path(studyFolder, "work_folder", "CloneCensorWeightModule"),
    resultsSubFolder = file.path(studyFolder, "results_folder", "CloneCensorWeightModule"),
    cdmDatabaseMetaData = list(databaseId = 123456789L)
  )
)
dir.create(handBuiltContext$moduleExecutionSettings$workSubFolder, recursive = TRUE)
dir.create(handBuiltContext$moduleExecutionSettings$resultsSubFolder, recursive = TRUE)
resultsDataModel <- quietly(replayExecute(handBuiltContext, eunomiaConnectionDetails))
# The files of the stages (runCcwAnalyses rewrites the specifications and the reference table on every run):
workFiles <- list.files(handBuiltContext$moduleExecutionSettings$workSubFolder, recursive = TRUE, full.names = TRUE)
workFiles <- workFiles[!basename(workFiles) %in% c("ccwAnalysesSpecifications.rds", "outcomeModelReference.rds")]
workFilesModified <- file.mtime(workFiles)

test_that("execute: runCcwAnalyses and exportToCsv write every table of the module's specification", {
  exportFolder <- handBuiltContext$moduleExecutionSettings$resultsSubFolder
  csvFiles <- list.files(exportFolder, pattern = "\\.csv$")
  expect_setequal(paste0(unique(resultsDataModel$tableName), ".csv"), csvFiles)
  expect_setequal(unique(resultsDataModel$tableName), unique(getResultsDataModelSpecifications()$tableName))
  expect_true(file.exists(file.path(exportFolder, "Results_123456789.zip")))
  result <- readr::read_csv(file.path(exportFolder, "ccw_result.csv"), col_types = readr::cols(database_id = "c"))
  expect_equal(unique(result$database_id), "123456789")
  expect_true("hr" %in% result$estimand)
  expect_true(all(is.finite(result$estimate[result$estimand == "hr"])))
})

test_that("execute again resumes from the work folder without recomputing", {
  quietly(replayExecute(handBuiltContext, eunomiaConnectionDetails))
  expect_identical(file.mtime(workFiles), workFilesModified)
})

test_that("createResultsDataModel and uploadResults: the first zip file uploads into SQLite", {
  resultsConnectionDetails <- newSqlite()
  quietly(createResultsDataModel(connectionDetails = resultsConnectionDetails, databaseSchema = "main", tablePrefix = ""))
  zipFile <- quietly(replayUploadResults(handBuiltContext$moduleExecutionSettings$resultsSubFolder, resultsConnectionDetails))
  expect_equal(basename(zipFile), "Results_123456789.zip")
  connection <- DatabaseConnector::connect(resultsConnectionDetails)
  on.exit(DatabaseConnector::disconnect(connection))
  exportFolder <- handBuiltContext$moduleExecutionSettings$resultsSubFolder
  for (tableName in unique(resultsDataModel$tableName)) {
    csvRows <- nrow(readr::read_csv(file.path(exportFolder, paste0(tableName, ".csv")), col_types = readr::cols()))
    uploadedRows <- DatabaseConnector::querySql(connection, sprintf("SELECT COUNT(*) FROM main.%s;", tableName))[1, 1]
    expect_equal(uploadedRows, csvRows, info = tableName)
  }
})

unlink(studyFolder, recursive = TRUE)

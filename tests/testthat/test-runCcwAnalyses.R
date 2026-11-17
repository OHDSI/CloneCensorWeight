library(testthat)

# runCcwAnalyses() and its summary (plan section 3.4, Phase 2 item 2). The unit tests at the top need no
# database. The Eunomia tests run the NSAID design of test-integration-eunomia.R (eligible cohort 4,
# treatment cohort 1 celecoxib, comparator treatment cohort 2 diclofenac, outcome cohort 3 GI bleed)
# with five negative-control outcome cohorts built from condition occurrences (Eunomia has none), and
# check: the reference table and the files shared between analyses, that a censoring model shared by
# analyses is fit once, the estimates against the single-study functions, resume after deleted files,
# the non-interactive handling of changed specifications, refitting the censoring models per outcome,
# a run with an open connection (cross-validated censoring models, a pooled logistic outcome model and
# competingEvents = "none"), and that a multi-threaded run gives the results of a single-threaded one.

# Settings ---------------------------------------------------------------------------------------

test_that("createCcwMultiThreadingSettings and the defaults", {
  settings <- createCcwMultiThreadingSettings(fitOutcomeModelThreads = 3)
  expect_s3_class(settings, "CcwMultiThreadingSettings")
  expect_equal(settings$fitOutcomeModelThreads, 3)
  expect_equal(settings$getDbCcwDataThreads, 1)
  expect_error(createCcwMultiThreadingSettings(calibrationThreads = 0), "calibrationThreads")
  defaults <- createDefaultCcwMultiThreadingSettings(20)
  expect_equal(unclass(defaults), list(
    getDbCcwDataThreads = 3,
    createStudyPopThreads = 3,
    fitCensoringModelThreads = 2,
    censoringModelCvThreads = 10,
    computeBalanceThreads = 5,
    fitOutcomeModelThreads = 5,
    outcomeCvThreads = 4,
    calibrationThreads = 4
  ))
  expect_true(all(unlist(createDefaultCcwMultiThreadingSettings(1)) == 1))
  expect_error(createDefaultCcwMultiThreadingSettings(0))
})

orchestrationAnalysis <- function(analysisId,
                                  strategies,
                                  fitOutcomeModelArgs = createFitOutcomeModelArgs(horizons = 365),
                                  computeWeightsArgs = createComputeWeightsArgs(),
                                  sameDayOrdering = "eventsFirst",
                                  competingEvents = "death",
                                  # A fixed prior variance keeps the fits fast (no cross-validation):
                                  fitCensoringModelsArgs = createFitCensoringModelsArgs(
                                    prior = Cyclops::createPrior("laplace", exclude = c(0), variance = 1)
                                  ),
                                  covariateSettings = FeatureExtraction::createCovariateSettings(
                                    useDemographicsGender = TRUE,
                                    useDemographicsAgeGroup = TRUE,
                                    endDays = -1
                                  )) {
  return(createCcwAnalysis(
    analysisId = analysisId,
    getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(covariateSettings = covariateSettings),
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
      removeSubjectsWithPriorOutcome = TRUE,
      riskWindowEnd = 365,
      endAnchor = "cohort start"
    ),
    createClonesArgs = createCreateClonesArgs(strategies = strategies, sameDayOrdering = sameDayOrdering, competingEvents = competingEvents),
    fitCensoringModelsArgs = fitCensoringModelsArgs,
    computeWeightsArgs = computeWeightsArgs,
    computeCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(),
    fitOutcomeModelArgs = fitOutcomeModelArgs
  ))
}

initiate30 <- createTreatmentStrategy(1, "Celecoxib within 30 days", "initiateWithin", graceEnd = 30)
notInitiate30 <- createTreatmentStrategy(2, "No celecoxib within 30 days", "notInitiateWithin", graceEnd = 30)
neverInitiate <- createTreatmentStrategy(2, "Never celecoxib", "neverInitiate")

test_that("The reference table shares files between analyses with the same settings", {
  outputFolder <- tempfile()
  on.exit(unlink(outputFolder, recursive = TRUE))
  dir.create(outputFolder)
  analyses <- list(
    orchestrationAnalysis(1, list(initiate30, notInitiate30)),
    # Differs from 1 only in the outcome model:
    orchestrationAnalysis(2, list(initiate30, notInitiate30), fitOutcomeModelArgs = createFitOutcomeModelArgs(modelType = "km", horizons = 365)),
    # Differs from 1 in the comparator strategy:
    orchestrationAnalysis(3, list(initiate30, neverInitiate)),
    # Differs from 1 in the weights:
    orchestrationAnalysis(4, list(initiate30, notInitiate30), computeWeightsArgs = createComputeWeightsArgs(stabilized = FALSE)),
    # Differs from 1 in sameDayOrdering, which changes the risk sets:
    orchestrationAnalysis(5, list(initiate30, notInitiate30), sameDayOrdering = "deviationFirst"),
    # Differs from 1 in only the strategy labels and IDs, which the censoring models do not depend on:
    orchestrationAnalysis(6, list(
      createTreatmentStrategy(7, "Start", "initiateWithin", graceEnd = 30),
      createTreatmentStrategy(8, "Defer", "notInitiateWithin", graceEnd = 30)
    ))
  )
  ttos <- list(
    createTargetTreatmentOutcomes(4, 1, 2, outcomes = list(
      CohortMethod::createOutcome(3),
      CohortMethod::createOutcome(101, outcomeOfInterest = FALSE, trueEffectSize = 1),
      # Its own risk window, so its own study population (but the same censoring population):
      CohortMethod::createOutcome(102, outcomeOfInterest = TRUE, priorOutcomeLookback = 30)
    )),
    createTargetTreatmentOutcomes(4, 2, NULL, outcomes = list(CohortMethod::createOutcome(3)))
  )
  specifications <- createCcwAnalysesSpecifications(
    ccwAnalysisList = analyses,
    targetTreatmentOutcomesList = ttos,
    analysesToExclude = data.frame(analysisId = 6, treatmentId = 2)
  )
  referenceTable <- NULL
  expect_message(
    referenceTable <- .createCcwReferenceTable(specifications, outputFolder),
    "Removed 1 of the 24"
  )
  expect_equal(nrow(referenceTable), 23)
  expect_true(all(dir.exists(file.path(outputFolder, sprintf("Analysis_%d", 1:6)))))
  row <- function(analysisId, outcomeId, treatmentId = 1) {
    referenceTable[referenceTable$analysisId == analysisId & referenceTable$outcomeId == outcomeId & referenceTable$treatmentId == treatmentId, ]
  }
  expect_equal(row(1, 3)$ccwDataFile, "CcwData_l1_t4_tr1_c2.zip")
  expect_equal(row(1, 3, treatmentId = 2)$ccwDataFile, "CcwData_l1_t4_tr2.zip")
  expect_equal(row(1, 3)$studyPopFile, "StudyPop_l1_s1_t4_tr1_c2_o3.rds")
  expect_equal(row(1, 101)$studyPopFile, "")
  expect_equal(row(1, 102)$studyPopFile, "StudyPop_l1_s2_t4_tr1_c2_o102.rds")
  expect_equal(row(1, 3)$outcomeModelFile, file.path("Analysis_1", "om_t4_tr1_c2_o3.rds"))
  # Shared censoring models across outcomes, and across analyses 1, 2, 4 and 6:
  expect_equal(unique(referenceTable$censoringModelFile[referenceTable$analysisId %in% c(1, 2, 4, 6) & referenceTable$treatmentId == 1]), "CensModel_l1_s1_m1_t4_tr1_c2.rds")
  expect_equal(row(3, 3)$censoringModelFile, "CensModel_l1_s1_m2_t4_tr1_c2.rds")
  expect_equal(row(5, 3)$censoringModelFile, "CensModel_l1_s1_m3_t4_tr1_c2.rds")
  # Balance is shared by all outcomes, and by analyses 1, 2 and 6, but not 4 (other weights):
  expect_equal(row(1, 101)$balanceFile, row(1, 3)$balanceFile)
  expect_equal(row(2, 3)$balanceFile, row(1, 3)$balanceFile)
  expect_equal(row(6, 3)$balanceFile, row(1, 3)$balanceFile)
  expect_false(row(4, 3)$balanceFile == row(1, 3)$balanceFile)
  # One censoring task per censoring population and sameDayOrdering:
  expect_equal(length(unique(referenceTable$censoringGroup[referenceTable$treatmentId == 1])), 2)

  # Refitting per outcome: the censoring files are per outcome, balance only for outcomes of interest.
  specifications$refitCensoringModelsForEveryOutcome <- TRUE
  specifications$analysesToExclude <- NULL
  referenceTable <- .createCcwReferenceTable(specifications, outputFolder)
  expect_equal(row(1, 3)$censoringModelFile, "CensModel_l1_s1_m1_t4_tr1_c2_o3.rds")
  expect_equal(row(1, 102)$censoringModelFile, "CensModel_l1_s2_m1_t4_tr1_c2_o102.rds")
  expect_equal(row(1, 101)$balanceFile, "")
  expect_equal(row(1, 3)$balanceFile, "Balance_l1_s1_m1_w1_b1_t4_tr1_c2_o3.rds")
})

# Summary ----------------------------------------------------------------------------------------

test_that("summarizeResults applies the thresholds, calibrates with the eligible negative controls and unblinds", {
  outputFolder <- tempfile()
  on.exit(unlink(outputFolder, recursive = TRUE))
  dir.create(file.path(outputFolder, "Analysis_1"), recursive = TRUE)
  outcomeIds <- c(1, 11:16)
  ncLogRr <- c(0.05, -0.1, 0.12, 0.02, -0.04, 2)
  ncSe <- c(0.1, 0.12, 0.15, 0.1, 0.2, 0.1)
  referenceTable <- dplyr::tibble(
    analysisId = 1L,
    targetId = 4,
    treatmentId = 1,
    comparatorTreatmentId = 2,
    outcomeId = outcomeIds,
    outcomeOfInterest = outcomeIds == 1,
    trueEffectSize = ifelse(outcomeIds == 1, NA, 1),
    balanceFile = "Balance.rds",
    outcomeModelFile = file.path("Analysis_1", sprintf("om_o%d.rds", outcomeIds))
  )
  saveRDS(data.frame(afterMatchingStdDiff = c(0.05, -0.08, NA)), file.path(outputFolder, "Balance.rds"))
  logRr <- c(0.4, ncLogRr)
  seLogRr <- c(0.15, ncSe)
  for (i in seq_along(outcomeIds)) {
    outcomeModel <- list(
      estimates = dplyr::tibble(
        estimand = c("hr", "rd_365"),
        estimate = c(exp(logRr[i]), 0.01),
        logEstimate = c(logRr[i], NA),
        seLogEstimate = c(seLogRr[i], NA),
        se = c(NA, 0.005),
        lb95 = c(exp(logRr[i] - 1.96 * seLogRr[i]), 0.0002),
        ub95 = c(exp(logRr[i] + 1.96 * seLogRr[i]), 0.0198),
        p = c(2 * stats::pnorm(-abs(logRr[i] / seLogRr[i])), 0.045)
      ),
      counts = dplyr::tibble(
        strategyId = 1:2, clones = c(100L, 90L), persons = c(100L, 90L), outcomes = c(5L, 4L),
        competingEvents = c(1L, 0L), personDays = c(1000, 900), weightedPersonDays = c(1010, 880),
        weightedOutcomes = c(5.5, 3.9)
      ),
      diagnostics = list(values = dplyr::tibble(
        # The outcome of interest fails the MDRR, and the last negative control the positivity:
        mdrr = if (outcomeIds[i] == 1) 20 else 2,
        positivity = if (outcomeIds[i] == 16) 0.001 else 0.2,
        positivityShareExtreme = if (outcomeIds[i] == 16) 0.3 else 0,
        essFraction = if (outcomeIds[i] == 1) 0.5 else NA,
        maxWeight = if (outcomeIds[i] == 1) 4 else NA,
        maxWeightBeforeTruncation = if (outcomeIds[i] == 1) 9 else NA
      ))
    )
    saveRDS(outcomeModel, file.path(outputFolder, referenceTable$outcomeModelFile[i]))
  }
  set.seed(123)
  summarizeResults(
    referenceTable = referenceTable,
    outputFolder = outputFolder,
    mainFileName = file.path(outputFolder, "resultsSummary.rds"),
    diagnosticsSummaryFileName = file.path(outputFolder, "diagnosticsSummary.rds"),
    calibrationThreads = 1,
    ccwDiagnosticThresholds = createCcwDiagnosticThresholds()
  ) |> expect_message("Calibrating")
  results <- getResultsSummary(outputFolder)
  diagnostics <- getDiagnosticsSummary(outputFolder)
  expect_equal(nrow(results), 14)
  expect_equal(nrow(diagnostics), 7)

  diagnostic <- function(outcomeId, column) diagnostics[[column]][diagnostics$outcomeId == outcomeId]
  expect_equal(diagnostic(1, "maxSdm"), 0.08)
  expect_equal(diagnostic(1, "balanceDiagnostic"), "PASS")
  expect_equal(diagnostic(1, "mdrrDiagnostic"), "FAIL")
  expect_equal(diagnostic(11, "mdrrDiagnostic"), "PASS")
  expect_equal(diagnostic(1, "essFractionDiagnostic"), "PASS")
  expect_equal(diagnostic(11, "essFractionDiagnostic"), "NOT EVALUATED")
  expect_equal(diagnostic(1, "maxWeightDiagnostic"), "PASS")
  expect_equal(diagnostic(1, "maxWeightBeforeTruncation"), 9)
  expect_equal(diagnostic(16, "positivityDiagnostic"), "FAIL")
  expect_equal(diagnostic(16, "positivityShareExtreme"), 0.3)

  # Calibration uses the five negative controls that pass, not outcome 16:
  set.seed(123)
  null <- EmpiricalCalibration::fitMcmcNull(logRr = ncLogRr[1:5], seLogRr = ncSe[1:5])
  ease <- EmpiricalCalibration::computeExpectedAbsoluteSystematicError(null)$ease
  expect_equal(unique(diagnostics$ease), ease)
  expect_equal(unique(diagnostics$easeDiagnostic), if (ease < 0.25) "PASS" else "FAIL")
  hr <- results[results$outcomeId == 1 & results$estimand == "hr", ]
  expect_equal(hr$calibratedP, EmpiricalCalibration::calibrateP(null, 0.4, 0.15)$p)
  calibratedCi <- EmpiricalCalibration::calibrateConfidenceInterval(
    logRr = 0.4, seLogRr = 0.15, model = EmpiricalCalibration::convertNullToErrorModel(null)
  )
  expect_equal(hr$calibratedRr, exp(calibratedCi$logRr))
  expect_equal(hr$calibratedCi95Lb, exp(calibratedCi$logLb95Rr))
  expect_equal(hr$oneSidedP, EmpiricalCalibration::computeTraditionalP(0.4, 0.15, twoSided = FALSE, upper = TRUE))
  expect_equal(hr$targetSubjects, 100)
  expect_equal(hr$comparatorDays, 900)
  expect_equal(hr$targetWeightedOutcomes, 5.5)
  expect_equal(hr$comparatorWeightedDays, 880)
  rd <- results[results$outcomeId == 1 & results$estimand == "rd_365", ]
  expect_true(is.na(rd$calibratedP) && is.na(rd$logRr) && is.na(rd$oneSidedP))
  expect_equal(rd$se, 0.005)

  # Unblinding, as in CohortMethod:
  expect_equal(diagnostics$unblindForEvidenceSynthesis, diagnostics$positivityDiagnostic != "FAIL" & diagnostics$easeDiagnostic != "FAIL")
  expect_equal(diagnostics$unblind, diagnostics$unblindForEvidenceSynthesis & diagnostics$mdrrDiagnostic != "FAIL")
  expect_false(diagnostic(1, "unblind"))
  expect_false(diagnostic(16, "unblindForEvidenceSynthesis"))

  # Fewer than 5 eligible negative controls: no calibration, EASE not evaluated.
  file.remove(file.path(outputFolder, referenceTable$outcomeModelFile[2]))
  summarizeResults(
    referenceTable = referenceTable[-2, ],
    outputFolder = outputFolder,
    mainFileName = file.path(outputFolder, "resultsSummary.rds"),
    diagnosticsSummaryFileName = file.path(outputFolder, "diagnosticsSummary.rds"),
    calibrationThreads = 1,
    ccwDiagnosticThresholds = createCcwDiagnosticThresholds()
  ) |> suppressMessages()
  expect_true(all(is.na(getResultsSummary(outputFolder)$calibratedP)))
  expect_true(all(getDiagnosticsSummary(outputFolder)$easeDiagnostic == "NOT EVALUATED"))
})

test_that("Positive controls calibrate the confidence intervals of the hazard ratio only", {
  # Synthetic positive controls have their true effect size on the hazard scale (plan section 3.8 (e)).
  trueEffectSize <- c(rep(1, 6), rep(2, 5), NA)
  logRr <- c(0.05, -0.1, 0.12, 0.02, -0.04, 0.08, log(2) + c(0.1, 0.2, 0.05, 0.15, 0.12), 0.3)
  seLogRr <- c(0.1, 0.12, 0.15, 0.1, 0.2, 0.1, 0.1, 0.12, 0.15, 0.1, 0.2, 0.15)
  group <- function(estimand) {
    dplyr::tibble(estimand = estimand, trueEffectSize = trueEffectSize, logRr = logRr, seLogRr = seLogRr, unblindForCalibration = TRUE)
  }
  set.seed(1)
  hr <- calibrateGroup(group("hr"))
  set.seed(1)
  rr <- calibrateGroup(group("rr_365"))
  controls <- !is.na(trueEffectSize)
  negativeControls <- trueEffectSize %in% 1
  set.seed(1)
  null <- EmpiricalCalibration::fitMcmcNull(logRr = logRr[negativeControls], seLogRr = seLogRr[negativeControls])
  model <- EmpiricalCalibration::fitSystematicErrorModel(
    logRr = logRr[controls],
    seLogRr = seLogRr[controls],
    trueLogRr = log(trueEffectSize[controls]),
    estimateCovarianceMatrix = FALSE
  )
  expected <- EmpiricalCalibration::calibrateConfidenceInterval(logRr = logRr, seLogRr = seLogRr, model = model)
  expect_equal(hr$calibratedCi95Lb, exp(expected$logLb95Rr))
  expect_equal(hr$calibratedCi95Ub, exp(expected$logUb95Rr))
  expected <- EmpiricalCalibration::calibrateConfidenceInterval(
    logRr = logRr, seLogRr = seLogRr, model = EmpiricalCalibration::convertNullToErrorModel(null)
  )
  expect_equal(rr$calibratedCi95Lb, exp(expected$logLb95Rr))
  expect_equal(rr$calibratedCi95Ub, exp(expected$logUb95Rr))
  expect_false(isTRUE(all.equal(hr$calibratedCi95Lb, rr$calibratedCi95Lb)))
  # The p-values use the negative controls for every estimand:
  expect_equal(rr$calibratedP, hr$calibratedP)
  expect_equal(hr$calibratedP, EmpiricalCalibration::calibrateP(null, logRr, seLogRr)$p)
})

test_that("Specifications with conflicting target-treatment-outcomes, or with nothing to run, stop", {
  analysis <- orchestrationAnalysis(1, list(initiate30, notInitiate30))
  outcome3 <- list(CohortMethod::createOutcome(3))
  outcome101 <- list(CohortMethod::createOutcome(101, outcomeOfInterest = FALSE, trueEffectSize = 1))
  # The same cohorts in two entries share the data, so they need the same covariate concept IDs:
  expect_error(createCcwAnalysesSpecifications(list(analysis), list(
    createTargetTreatmentOutcomes(4, 1, 2, outcomes = outcome3, excludedCovariateConceptIds = 1118084),
    createTargetTreatmentOutcomes(4, 1, 2, outcomes = outcome101)
  )), "same includedCovariateConceptIds")
  expect_error(createCcwAnalysesSpecifications(list(analysis), list(
    createTargetTreatmentOutcomes(4, 1, 2, outcomes = outcome3),
    createTargetTreatmentOutcomes(4, 1, 2, outcomes = c(outcome3, outcome101))
  )), "more than once")
  # Outcomes of the same cohorts split over two entries, or other cohorts, are allowed:
  expect_no_error(createCcwAnalysesSpecifications(list(analysis), list(
    createTargetTreatmentOutcomes(4, 1, 2, outcomes = outcome3, excludedCovariateConceptIds = 1118084),
    createTargetTreatmentOutcomes(4, 1, 2, outcomes = outcome101, excludedCovariateConceptIds = 1118084),
    createTargetTreatmentOutcomes(4, 1, NULL, outcomes = outcome3),
    createTargetTreatmentOutcomes(4, 2, 1, outcomes = outcome3)
  )))

  specifications <- createCcwAnalysesSpecifications(
    list(analysis),
    list(createTargetTreatmentOutcomes(4, 1, 2, outcomes = outcome3)),
    analysesToExclude = data.frame(analysisId = 1)
  )
  outputFolder <- tempfile()
  on.exit(unlink(outputFolder, recursive = TRUE))
  expect_error(
    suppressMessages(runCcwAnalyses(
      connectionDetails = DatabaseConnector::createConnectionDetails(dbms = "sqlite", server = tempfile()),
      cdmDatabaseSchema = "main",
      outputFolder = outputFolder,
      ccwAnalysesSpecifications = specifications
    )),
    "nothing to run"
  )
})

# Eunomia ----------------------------------------------------------------------------------------

skip_on_cran()
skip_if_not_installed("Eunomia")
skip_if_not_installed("RSQLite")
skip_if_not_installed("EmpiricalCalibration")

databaseFile <- tempfile(fileext = ".sqlite")
connectionDetails <- tryCatch(
  suppressMessages(Eunomia::getEunomiaConnectionDetails(databaseFile = databaseFile)),
  error = function(e) NULL
)
if (is.null(connectionDetails)) {
  skip("Eunomia GiBleed data could not be loaded")
}
utils::capture.output(suppressMessages(Eunomia::createCohorts(connectionDetails)))
# Negative-control outcome cohorts 101 to 105: every occurrence of a common condition.
negativeControlConceptIds <- c(4112343, 260139, 28060, 4283893, 4294548)
local({
  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection))
  for (i in seq_along(negativeControlConceptIds)) {
    DatabaseConnector::executeSql(connection, sprintf("
      INSERT INTO main.cohort (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
      SELECT %d, person_id, condition_start_date, condition_start_date
      FROM main.condition_occurrence
      WHERE condition_concept_id = %d;", 100 + i, negativeControlConceptIds[i]), progressBar = FALSE, reportOverallTime = FALSE)
  }
})
negativeControls <- lapply(101:105, function(outcomeId) {
  CohortMethod::createOutcome(outcomeId, outcomeOfInterest = FALSE, trueEffectSize = 1)
})

runQuietly <- function(...) {
  messages <- character()
  value <- NULL
  utils::capture.output(value <- withCallingHandlers(
    runCcwAnalyses(connectionDetails = connectionDetails, cdmDatabaseSchema = "main", ...),
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  ))
  return(list(referenceTable = value, messages = messages))
}

# The messages of the stage whose heading contains `stage`.
stageMessages <- function(messages, stage) {
  starts <- grep("^\\*\\*\\* ", messages)
  start <- starts[grepl(stage, messages[starts], fixed = TRUE)]
  if (length(start) == 0) {
    return(character())
  }
  end <- c(starts[starts > start], length(messages) + 1)[1]
  return(messages[seq_len(end - start - 1) + start])
}

estimateColumns <- c(
  "analysisId", "targetId", "treatmentId", "comparatorTreatmentId", "outcomeId", "estimand", "targetSubjects",
  "comparatorSubjects", "targetClones", "comparatorClones", "targetDays", "comparatorDays", "targetOutcomes",
  "comparatorOutcomes", "estimate", "ci95Lb", "ci95Ub", "p", "oneSidedP", "logRr", "seLogRr", "se"
)

mainAnalyses <- list(
  orchestrationAnalysis(1, list(initiate30, notInitiate30)),
  orchestrationAnalysis(2, list(initiate30, notInitiate30), fitOutcomeModelArgs = createFitOutcomeModelArgs(modelType = "km", horizons = 365)),
  orchestrationAnalysis(3, list(initiate30, neverInitiate)),
  orchestrationAnalysis(4, list(initiate30, notInitiate30), fitOutcomeModelArgs = createFitOutcomeModelArgs(
    modelType = "km",
    horizons = 365,
    varianceEstimator = "bootstrap",
    bootstrapArgs = createBootstrapArgs(replicates = 2)
  ))
)
mainSpecifications <- createCcwAnalysesSpecifications(
  ccwAnalysisList = mainAnalyses,
  targetTreatmentOutcomesList = list(createTargetTreatmentOutcomes(4, 1, 2, outcomes = c(list(CohortMethod::createOutcome(3)), negativeControls))),
  # Only analysis 1 runs the negative controls:
  analysesToExclude = data.frame(analysisId = rep(2:4, each = 5), outcomeId = rep(101:105, 3))
)
mainFolder <- tempfile()
mainRun <- runQuietly(outputFolder = mainFolder, ccwAnalysesSpecifications = mainSpecifications)

test_that("runCcwAnalyses runs several analyses on Eunomia and fits a shared censoring model once", {
  referenceTable <- mainRun$referenceTable
  expect_equal(referenceTable, getFileReference(mainFolder))
  expect_equal(nrow(referenceTable), 6 + 3)
  expect_setequal(
    names(referenceTable),
    c("analysisId", "targetId", "treatmentId", "comparatorTreatmentId", "outcomeId", "outcomeOfInterest", "trueEffectSize",
      "ccwDataFile", "studyPopFile", "censoringModelFile", "balanceFile", "outcomeModelFile")
  )
  files <- unlist(referenceTable[, c("ccwDataFile", "studyPopFile", "censoringModelFile", "balanceFile", "outcomeModelFile")])
  expect_true(all(file.exists(file.path(mainFolder, files[files != ""]))))
  expect_equal(unique(referenceTable$ccwDataFile), "CcwData_l1_t4_tr1_c2.zip")
  expect_equal(unique(referenceTable$studyPopFile[referenceTable$outcomeId == 3]), "StudyPop_l1_s1_t4_tr1_c2_o3.rds")
  expect_equal(length(unique(referenceTable$censoringModelFile)), 2)
  expect_true(file.exists(file.path(mainFolder, "ccwAnalysesSpecifications.rds")))
  expect_false(any(grepl("\\.tmp", list.files(mainFolder, recursive = TRUE))))

  # Each censoring model is fit once: the binary model of "initiate within 30 days" is shared by the
  # three analyses that compare it with "not within 30 days" (one file) and reused for the analysis
  # that compares it with "never".
  censoringMessages <- stageMessages(mainRun$messages, "Fitting censoring models")
  expect_equal(sum(grepl("^Fitting binary censoring model", censoringMessages)), 1)
  expect_equal(sum(grepl("^Reusing the fitted binary censoring model", censoringMessages)), 1)
  expect_equal(sum(grepl("^Fitting logistic initiation hazard model through day 30 ", censoringMessages)), 1)
  expect_equal(sum(grepl("^Fitting logistic initiation hazard model through day 365 ", censoringMessages)), 1)
  expect_equal(sum(grepl("^Computing covariate balance took", stageMessages(mainRun$messages, "Computing covariate balance"))), 2)

  results <- getResultsSummary(mainFolder)
  diagnostics <- getDiagnosticsSummary(mainFolder)
  expect_equal(nrow(diagnostics), nrow(referenceTable))
  expect_equal(
    unique(results[, c("analysisId", "outcomeId")]),
    referenceTable[, c("analysisId", "outcomeId")],
    ignore_attr = TRUE
  )
  expect_true(all(c("hr", "risk_target_365", "risk_comparator_365", "rd_365", "rr_365") %in% results$estimand[results$analysisId == 1 & results$outcomeId == 3]))
  expect_true(all(is.finite(results$estimate)))
  # Five negative controls, so the ratio estimands of analysis 1 are calibrated:
  calibrated <- results[results$analysisId == 1 & results$estimand %in% c("hr", "rr_365"), ]
  expect_true(all(!is.na(calibrated$calibratedP)))
  expect_true(all(is.na(results$calibratedP[!results$estimand %in% c("hr", "rr_365")])))
  expect_true(all(!is.na(diagnostics$ease[diagnostics$analysisId == 1])))
  expect_true(all(is.na(diagnostics$ease[diagnostics$analysisId != 1])))
  # Diagnostics of the weights for every outcome, because they decide which negative controls
  # calibrate; balance shared by all outcomes:
  expect_true(all(!is.na(diagnostics$essFraction)))
  expect_true(all(!is.na(diagnostics$maxWeight)))
  expect_true(all(diagnostics$maxWeightBeforeTruncation >= diagnostics$maxWeight))
  expect_equal(length(unique(diagnostics$maxSdm[diagnostics$analysisId == 1])), 1)
  expect_true(all(diagnostics$balanceDiagnostic == "PASS"))
  expect_equal(diagnostics$unblind, diagnostics$unblindForEvidenceSynthesis & diagnostics$mdrrDiagnostic != "FAIL")

  # The outcome model file holds the model, its diagnostics and the curves:
  outcomeModel <- readRDS(file.path(mainFolder, referenceTable$outcomeModelFile[referenceTable$analysisId == 1 & referenceTable$outcomeId == 3]))
  expect_s3_class(outcomeModel, "CcwOutcomeModel")
  expect_s3_class(outcomeModel$diagnostics$weightSummary, "CcwWeightSummary")
  expect_s3_class(outcomeModel$diagnostics$positivity, "CcwPositivity")
  expect_true(nrow(outcomeModel$curves) > 0)
  negativeControlModel <- readRDS(file.path(mainFolder, referenceTable$outcomeModelFile[referenceTable$analysisId == 1 & referenceTable$outcomeId == 101]))
  expect_null(negativeControlModel$diagnostics$weightSummary)
  expect_null(negativeControlModel$curves)
  bootstrapModel <- readRDS(file.path(mainFolder, referenceTable$outcomeModelFile[referenceTable$analysisId == 4]))
  expect_equal(bootstrapModel$varianceEstimator, "bootstrap")
  expect_equal(bootstrapModel$metaData$bootstrap$replicates, 2)
})

test_that("runCcwAnalyses gives the estimates of the single-study functions with censoring models shared across outcomes", {
  referenceTable <- mainRun$referenceTable
  data <- loadCloneCensorWeightData(file.path(mainFolder, referenceTable$ccwDataFile[1]))
  on.exit(close(data))
  analysis <- mainAnalyses[[1]]
  # The shared censoring models are fit on the population without an outcome:
  populationWithoutOutcome <- suppressMessages(.createStudyPopulationWithoutOutcome(data, analysis$createStudyPopulationArgs))
  clonesWithoutOutcome <- suppressMessages(createClones(populationWithoutOutcome, data, analysis$createClonesArgs))
  on.exit(close(clonesWithoutOutcome), add = TRUE)
  censoringModels <- suppressMessages(fitCensoringModels(clonesWithoutOutcome, data, analysis$fitCensoringModelsArgs))
  expect_true(is.na(attr(censoringModels, "metaData")$outcomeId))
  savedModels <- readRDS(file.path(mainFolder, referenceTable$censoringModelFile[referenceTable$analysisId == 1][1]))
  expect_equal(savedModels[[1]]$coefficients, censoringModels[[1]]$coefficients, tolerance = 1e-10)
  expect_equal(savedModels[[2]]$coefficients, censoringModels[[2]]$coefficients, tolerance = 1e-10)

  results <- getResultsSummary(mainFolder)
  for (outcomeId in c(3, 101)) {
    population <- suppressMessages(createStudyPopulation(data, outcomeId = outcomeId, createStudyPopulationArgs = analysis$createStudyPopulationArgs))
    clones <- suppressMessages(createClones(population, data, analysis$createClonesArgs))
    weightedClones <- suppressMessages(computeCensoringWeights(clones, censoringModels, analysis$computeWeightsArgs))
    outcomeModel <- suppressMessages(fitOutcomeModel(weightedClones, analysis$fitOutcomeModelArgs))
    expected <- outcomeModel$estimates
    actual <- results[results$analysisId == 1 & results$outcomeId == outcomeId, ]
    expect_equal(actual$estimand, expected$estimand)
    expect_equal(actual$estimate, expected$estimate, tolerance = 1e-10)
    expect_equal(actual$ci95Ub, expected$ub95, tolerance = 1e-10)
    expect_equal(unique(actual$targetOutcomes), outcomeModel$counts$outcomes[1])
    mdrr <- computeMdrr(clones)$mdrr
    expect_equal(getDiagnosticsSummary(mainFolder)$mdrr[getDiagnosticsSummary(mainFolder)$analysisId == 1 & getDiagnosticsSummary(mainFolder)$outcomeId == outcomeId], mdrr)
    # The outcome of interest also keeps the attrition and the follow-up of the clones:
    savedModel <- readRDS(file.path(mainFolder, referenceTable$outcomeModelFile[referenceTable$analysisId == 1 & referenceTable$outcomeId == outcomeId]))
    if (outcomeId == 3) {
      expect_equal(savedModel$attrition, attr(clones, "metaData")$attrition)
      survivalTime <- dplyr::collect(clones$clones)
      expect_equal(
        savedModel$followUpDistribution$medianDays,
        as.numeric(tapply(survivalTime$survivalTime, survivalTime$strategyId, stats::median))
      )
      expect_equal(savedModel$followUpDistribution$maxDays, as.numeric(tapply(survivalTime$survivalTime, survivalTime$strategyId, max)))
    } else {
      expect_null(savedModel$attrition)
      expect_null(savedModel$followUpDistribution)
    }
    close(weightedClones)
    close(clones)
  }
})

test_that("runCcwAnalyses resumes after deleted outputs and recomputes only those", {
  referenceTable <- mainRun$referenceTable
  oldResults <- getResultsSummary(mainFolder)
  row <- referenceTable[referenceTable$analysisId == 1 & referenceTable$outcomeId == 3, ]
  neverModelFile <- referenceTable$censoringModelFile[referenceTable$analysisId == 3]
  oldNeverModels <- readRDS(file.path(mainFolder, neverModelFile))
  oldOutcomeModel <- readRDS(file.path(mainFolder, row$outcomeModelFile))
  unlink(file.path(mainFolder, c(row$outcomeModelFile, neverModelFile, row$balanceFile)))

  rerun <- runQuietly(outputFolder = mainFolder, ccwAnalysesSpecifications = mainSpecifications)
  expect_false(any(grepl("Creating CloneCensorWeightData|Creating study populations", rerun$messages)))
  # The censoring models of analysis 3 are refit, reusing the shared binary model of the other file:
  censoringMessages <- stageMessages(rerun$messages, "Fitting censoring models")
  expect_equal(sum(grepl("^Fitting binary censoring model", censoringMessages)), 0)
  expect_equal(sum(grepl("^Reusing the fitted binary censoring model", censoringMessages)), 1)
  expect_equal(sum(grepl("^Fitting logistic initiation hazard model", censoringMessages)), 1)
  expect_equal(sum(grepl("^Computing covariate balance took", stageMessages(rerun$messages, "Computing covariate balance"))), 1)
  expect_equal(sum(grepl("^Fitting outcome model took", stageMessages(rerun$messages, "Fitting outcome models"))), 1)
  expect_true(any(grepl("Summarizing results", rerun$messages)))

  newNeverModels <- readRDS(file.path(mainFolder, neverModelFile))
  expect_equal(newNeverModels[[2]]$coefficients, oldNeverModels[[2]]$coefficients, tolerance = 1e-10)
  newOutcomeModel <- readRDS(file.path(mainFolder, row$outcomeModelFile))
  expect_equal(newOutcomeModel$estimates, oldOutcomeModel$estimates, tolerance = 1e-10)
  expect_equal(getResultsSummary(mainFolder)[, estimateColumns], oldResults[, estimateColumns], tolerance = 1e-10)

  # Nothing to do: no stage runs, and the summary is kept.
  summaryTime <- file.mtime(file.path(mainFolder, "resultsSummary.rds"))
  rerun <- runQuietly(outputFolder = mainFolder, ccwAnalysesSpecifications = mainSpecifications)
  expect_false(any(grepl("^\\*\\*\\*", rerun$messages)))
  expect_equal(file.mtime(file.path(mainFolder, "resultsSummary.rds")), summaryTime)
})

# A smaller design for the remaining tests: two target-treatment combinations (celecoxib and diclofenac
# in the treatment role), so the data and censoring stages have two tasks each.
smallSpecifications <- function(refit = FALSE, thresholds = createCcwDiagnosticThresholds()) {
  outcomes <- list(CohortMethod::createOutcome(3), negativeControls[[1]])
  return(createCcwAnalysesSpecifications(
    ccwAnalysisList = list(orchestrationAnalysis(1, list(initiate30, notInitiate30))),
    targetTreatmentOutcomesList = list(
      createTargetTreatmentOutcomes(4, 1, 2, outcomes = outcomes),
      createTargetTreatmentOutcomes(4, 2, 1, outcomes = outcomes)
    ),
    refitCensoringModelsForEveryOutcome = refit,
    ccwDiagnosticThresholds = thresholds
  ))
}
smallFolder <- tempfile()
smallRun <- runQuietly(outputFolder = smallFolder, ccwAnalysesSpecifications = smallSpecifications())

test_that("runCcwAnalyses deletes the outputs of changed specifications outside an interactive session", {
  changedFolder <- tempfile()
  on.exit(unlink(changedFolder, recursive = TRUE))
  dir.create(changedFolder)
  file.copy(list.files(smallFolder, full.names = TRUE), changedFolder, recursive = TRUE)
  writeLines("stale", file.path(changedFolder, "stale.txt"))
  # Unchanged specifications keep the files:
  rerun <- runQuietly(outputFolder = changedFolder, ccwAnalysesSpecifications = smallSpecifications())
  expect_true(file.exists(file.path(changedFolder, "stale.txt")))
  expect_false(any(grepl("have changed", rerun$messages)))

  # Changed specifications (with one combination, to keep the rerun short):
  changed <- smallSpecifications(thresholds = createCcwDiagnosticThresholds(mdrrThreshold = 1))
  changed$analysesToExclude <- data.frame(treatmentId = 2)
  rerun <- runQuietly(outputFolder = changedFolder, ccwAnalysesSpecifications = changed)
  expect_true(any(grepl("the analysis settings have changed", rerun$messages)))
  expect_true(any(grepl("Deleting the old files", rerun$messages)))
  expect_false(file.exists(file.path(changedFolder, "stale.txt")))
  expect_true(any(grepl("Creating CloneCensorWeightData", rerun$messages)))
  expect_equal(readRDS(file.path(changedFolder, "ccwAnalysesSpecifications.rds"))$toList(), changed$toList())
  diagnostics <- getDiagnosticsSummary(changedFolder)
  expect_true(all(diagnostics$mdrrDiagnostic == "FAIL"))
  smallResults <- getResultsSummary(smallFolder)
  expect_equal(
    getResultsSummary(changedFolder)[, estimateColumns],
    smallResults[smallResults$treatmentId == 1, estimateColumns],
    tolerance = 1e-10,
    ignore_attr = TRUE
  )
})

test_that("runCcwAnalyses refits the censoring models for every outcome when asked", {
  refitFolder <- tempfile()
  on.exit(unlink(refitFolder, recursive = TRUE))
  specifications <- smallSpecifications(refit = TRUE)
  specifications$analysesToExclude <- data.frame(treatmentId = 1)
  run <- runQuietly(outputFolder = refitFolder, ccwAnalysesSpecifications = specifications)
  referenceTable <- run$referenceTable
  expect_equal(referenceTable$censoringModelFile, sprintf("CensModel_l1_s1_m1_t4_tr2_c1_o%d.rds", c(3, 101)))
  expect_equal(referenceTable$balanceFile[referenceTable$outcomeId == 101], "")
  # One censoring task per outcome:
  expect_equal(sum(grepl("^Fitting binary censoring model", stageMessages(run$messages, "Fitting censoring models"))), 2)
  results <- getResultsSummary(refitFolder)
  analysis <- smallSpecifications()$ccwAnalysisList[[1]]
  data <- loadCloneCensorWeightData(file.path(refitFolder, referenceTable$ccwDataFile[1]))
  on.exit(close(data), add = TRUE)
  for (outcomeId in c(3, 101)) {
    # The standard single-study pipeline:
    population <- suppressMessages(createStudyPopulation(data, outcomeId = outcomeId, createStudyPopulationArgs = analysis$createStudyPopulationArgs))
    clones <- suppressMessages(createClones(population, data, analysis$createClonesArgs))
    censoringModels <- suppressMessages(fitCensoringModels(clones, data, analysis$fitCensoringModelsArgs))
    weightedClones <- suppressMessages(computeCensoringWeights(clones, censoringModels, analysis$computeWeightsArgs))
    expected <- suppressMessages(fitOutcomeModel(weightedClones, analysis$fitOutcomeModelArgs))$estimates
    actual <- results[results$treatmentId == 2 & results$outcomeId == outcomeId, ]
    expect_equal(actual$estimate, expected$estimate, tolerance = 1e-10)
    close(weightedClones)
    close(clones)
  }
  # Persons with an earlier outcome 101 are in the risk sets of the shared models but not in those of
  # outcome 101's own models, so the estimates differ:
  sharedResults <- getResultsSummary(smallFolder)
  expect_false(isTRUE(all.equal(
    sharedResults$estimate[sharedResults$treatmentId == 2 & sharedResults$outcomeId == 101],
    results$estimate[results$outcomeId == 101],
    tolerance = 1e-12
  )))
})

test_that("runCcwAnalyses runs with an open connection, cross-validated censoring models, a pooled logistic outcome model and competingEvents = 'none'", {
  analysis <- orchestrationAnalysis(
    1,
    list(initiate30, notInitiate30),
    fitOutcomeModelArgs = createFitOutcomeModelArgs(
      modelType = "logistic",
      horizons = 340,
      # 31-day intervals have a boundary at the end of the grace period:
      timeBasis = createTimeBasis(type = "indicator", intervalLength = 31)
    ),
    competingEvents = "none",
    fitCensoringModelsArgs = createFitCensoringModelsArgs()
  )
  specifications <- createCcwAnalysesSpecifications(
    ccwAnalysisList = list(analysis),
    targetTreatmentOutcomesList = list(createTargetTreatmentOutcomes(4, 1, 2, outcomes = list(CohortMethod::createOutcome(3))))
  )
  outputFolder <- tempfile()
  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit({
    DatabaseConnector::disconnect(connection)
    unlink(outputFolder, recursive = TRUE)
  })
  # Cross-validation shrinks the coefficients of the uninformative demographics to zero:
  muffleZeroCoefficients <- function(expr) {
    withCallingHandlers(expr, warning = function(w) {
      if (grepl("All coefficients (except maybe the intercept) are zero", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    })
  }
  referenceTable <- NULL
  utils::capture.output(suppressMessages(muffleZeroCoefficients(referenceTable <- runCcwAnalyses(
    connection = connection,
    cdmDatabaseSchema = "main",
    outputFolder = outputFolder,
    # Cyclops threads do not change the cross-validated fits:
    multiThreadingSettings = createCcwMultiThreadingSettings(censoringModelCvThreads = 2),
    ccwAnalysesSpecifications = specifications
  ))))
  # The connection stays open:
  expect_equal(nrow(DatabaseConnector::querySql(connection, "SELECT 1 AS one;")), 1)

  # The single-study pipeline, with the censoring models fit on the population without an outcome:
  data <- loadCloneCensorWeightData(file.path(outputFolder, referenceTable$ccwDataFile))
  on.exit(close(data), add = TRUE)
  population <- suppressMessages(.createStudyPopulationWithoutOutcome(data, analysis$createStudyPopulationArgs))
  clones <- suppressMessages(createClones(population, data, analysis$createClonesArgs))
  censoringModels <- suppressMessages(muffleZeroCoefficients(fitCensoringModels(clones, data, analysis$fitCensoringModelsArgs)))
  close(clones)
  savedModels <- readRDS(file.path(outputFolder, referenceTable$censoringModelFile))
  for (i in 1:2) {
    expect_equal(savedModels[[i]]$priorVariance, censoringModels[[i]]$priorVariance, tolerance = 1e-5)
    expect_equal(savedModels[[i]]$coefficients, censoringModels[[i]]$coefficients, tolerance = 1e-5)
  }
  population <- suppressMessages(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = analysis$createStudyPopulationArgs))
  clones <- suppressMessages(createClones(population, data, analysis$createClonesArgs))
  on.exit(close(clones), add = TRUE)
  weightedClones <- suppressMessages(computeCensoringWeights(
    clones, censoringModels, analysis$computeWeightsArgs,
    outcomeTimeBasis = analysis$fitOutcomeModelArgs$timeBasis
  ))
  on.exit(close(weightedClones), add = TRUE)
  expected <- suppressMessages(fitOutcomeModel(weightedClones, analysis$fitOutcomeModelArgs))$estimates
  results <- getResultsSummary(outputFolder)
  expect_equal(results$estimand, expected$estimand)
  expect_true("or" %in% results$estimand)
  expect_equal(results$estimate, expected$estimate, tolerance = 1e-5)
  # Deaths end the clones as observationEnd, so there are no competing events:
  expect_true(all(results$targetCompetingEvents == 0 & results$comparatorCompetingEvents == 0))
  expect_equal(getDiagnosticsSummary(outputFolder)$mdrr, computeMdrr(clones, modelType = "logistic")$mdrr)
})

test_that("A multi-threaded run gives the results of a single-threaded one", {
  # Worker processes load the installed package, which must hold this version.
  installed <- find.package("CloneCensorWeight", lib.loc = .libPaths(), quiet = TRUE)
  skip_if(
    length(installed) == 0 || !any(grepl("runCcwAnalyses", readLines(file.path(installed[1], "NAMESPACE")))),
    "The installed CloneCensorWeight has no runCcwAnalyses()"
  )
  threadedFolder <- tempfile()
  on.exit(unlink(threadedFolder, recursive = TRUE))
  runQuietly(
    outputFolder = threadedFolder,
    ccwAnalysesSpecifications = smallSpecifications(),
    multiThreadingSettings = createCcwMultiThreadingSettings(
      getDbCcwDataThreads = 2,
      createStudyPopThreads = 2,
      fitCensoringModelThreads = 2,
      censoringModelCvThreads = 2,
      computeBalanceThreads = 2,
      fitOutcomeModelThreads = 2,
      outcomeCvThreads = 2,
      calibrationThreads = 2
    )
  )
  expect_equal(getFileReference(threadedFolder), getFileReference(smallFolder))
  expect_equal(getResultsSummary(threadedFolder), getResultsSummary(smallFolder), tolerance = 1e-10)
  expect_equal(getDiagnosticsSummary(threadedFolder), getDiagnosticsSummary(smallFolder), tolerance = 1e-10)
})

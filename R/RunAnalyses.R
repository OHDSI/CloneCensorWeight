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

# runCcwAnalyses() follows the design of CohortMethod's runCmAnalyses() (R/RunAnalyses.R, Copyright
# Observational Health Data Sciences and Informatics, Apache 2.0): a reference table, stages that share
# files between analyses with the same settings, resume by file, and ParallelLogger clusters.

# The CloneCensorWeightData object last loaded by this process, reused by the tasks of a stage.
cache <- new.env(parent = emptyenv())

#' Create clone-censor-weight multi-threading settings
#'
#' @description
#' The numbers of parallel threads [runCcwAnalyses()] uses in each of its stages.
#'
#' @param getDbCcwDataThreads      The number of parallel threads to use for building the
#'                                 `CloneCensorWeightData` objects.
#' @param createStudyPopThreads    The number of parallel threads to use for creating the study
#'                                 populations.
#' @param fitCensoringModelThreads The number of parallel threads to use for fitting the censoring
#'                                 models. Each thread fits the models of one study population.
#' @param censoringModelCvThreads  The number of parallel threads to use for the cross-validation of
#'                                 the censoring models. The total number of threads at one time can
#'                                 be `fitCensoringModelThreads * censoringModelCvThreads`.
#' @param computeBalanceThreads    The number of parallel threads to use for computing covariate
#'                                 balance.
#' @param fitOutcomeModelThreads   The number of parallel threads to use for fitting the outcome
#'                                 models.
#' @param outcomeCvThreads         The number of threads Cyclops uses in each outcome model fit. The
#'                                 total number of threads at one time can be
#'                                 `fitOutcomeModelThreads * outcomeCvThreads`.
#' @param calibrationThreads       The number of parallel threads to use for empirical calibration.
#'
#' @return
#' An object of type `CcwMultiThreadingSettings`.
#'
#' @seealso [createDefaultCcwMultiThreadingSettings()]
#'
#' @examples
#' settings <- createCcwMultiThreadingSettings(fitOutcomeModelThreads = 2)
#'
#' @export
createCcwMultiThreadingSettings <- function(getDbCcwDataThreads = 1,
                                            createStudyPopThreads = 1,
                                            fitCensoringModelThreads = 1,
                                            censoringModelCvThreads = 1,
                                            computeBalanceThreads = 1,
                                            fitOutcomeModelThreads = 1,
                                            outcomeCvThreads = 1,
                                            calibrationThreads = 1) {
  settings <- list()
  for (name in names(formals(createCcwMultiThreadingSettings))) {
    settings[[name]] <- get(name)
  }
  errorMessages <- checkmate::makeAssertCollection()
  for (name in names(settings)) {
    checkmate::assertInt(settings[[name]], lower = 1, .var.name = name, add = errorMessages)
  }
  checkmate::reportAssertions(collection = errorMessages)
  class(settings) <- "CcwMultiThreadingSettings"
  return(settings)
}

#' Create default clone-censor-weight multi-threading settings
#'
#' @description
#' Create clone-censor-weight multi-threading settings based on the maximum number of cores to be used,
#' as CohortMethod's `createDefaultMultiThreadingSettings()` does.
#'
#' @param maxCores Maximum number of CPU cores to use.
#'
#' @return
#' An object of type `CcwMultiThreadingSettings`.
#'
#' @seealso [createCcwMultiThreadingSettings()]
#'
#' @examples
#' settings <- createDefaultCcwMultiThreadingSettings(10)
#'
#' @export
createDefaultCcwMultiThreadingSettings <- function(maxCores) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertInt(maxCores, lower = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  return(createCcwMultiThreadingSettings(
    getDbCcwDataThreads = min(3, maxCores),
    createStudyPopThreads = min(3, maxCores),
    fitCensoringModelThreads = max(1, floor(maxCores / 10)),
    censoringModelCvThreads = min(10, maxCores),
    computeBalanceThreads = min(5, maxCores),
    fitOutcomeModelThreads = max(1, floor(maxCores / 4)),
    outcomeCvThreads = min(4, maxCores),
    calibrationThreads = min(4, maxCores)
  ))
}

#' Run a list of clone-censor-weight analyses
#'
#' @details
#' Runs every analysis in `ccwAnalysesSpecifications$ccwAnalysisList` against every
#' target-treatment-outcome combination in `ccwAnalysesSpecifications$targetTreatmentOutcomesList`,
#' except those in `analysesToExclude` (see [createCcwAnalysesSpecifications()]). The work is done in
#' stages, each of which writes files to `outputFolder`. A file is computed once and shared by every
#' analysis whose settings for that stage (and the stages before it) are the same:
#'
#' 1. **Data** (`CcwData_l<l>_t<t>_tr<tr>[_c<c>].zip`): [getDbCloneCensorWeightData()] per
#'    `getDbCloneCensorWeightDataArgs` and target, treatment and comparator treatment cohort, with all
#'    outcomes of the combinations that use it. The `excludedCovariateConceptIds` and
#'    `includedCovariateConceptIds` of the target-treatment-outcomes are added to the covariate
#'    settings.
#' 2. **Study populations** (`StudyPop_l<l>_s<s>_..._o<o>.rds`): [createStudyPopulation()] per data,
#'    study population settings and outcome, for outcomes of interest. The `priorOutcomeLookback`,
#'    `riskWindowStart`, `startAnchor`, `riskWindowEnd` and `endAnchor` of an outcome (see
#'    [CohortMethod::createOutcome()]) override those of the analysis, as in CohortMethod. Populations
#'    of other outcomes (negative controls) are created when they are needed and not saved.
#' 3. **Censoring models** (`CensModel_l<l>_s<s>_m<m>_..[_o<o>].rds`): [fitCensoringModels()] per data,
#'    censoring population, censoring settings (`m`: the strategies' types, treatments and grace
#'    periods, `sameDayOrdering`, `competingEvents` and `fitCensoringModelsArgs`). With
#'    `refitCensoringModelsForEveryOutcome = FALSE` (the default), the censoring population is the
#'    study population without an outcome: persons with a prior outcome are not removed, and the risk
#'    sets do not end at the outcome. The models are then fit once and shared by every outcome.
#'    With `refitCensoringModelsForEveryOutcome = TRUE`, the models are fit on the clones of each
#'    outcome's study population (`_o<o>`). The models of one censoring population are fit by one
#'    thread, and a model shared by strategies of different analyses is fit once: for example the
#'    binary model of an `initiateWithin(G)` strategy, whether the comparator is
#'    `notInitiateWithin(G)` or `neverInitiate`.
#' 4. **Weights** are computed by [computeCensoringWeights()] when needed and not saved.
#' 5. **Covariate balance** (`Balance_l<l>_s<s>_m<m>_w<w>_b<b>_..[_o<o>].rds`):
#'    [computeCovariateBalance()] at the end of the grace periods, on the weighted clones of the
#'    censoring population, per censoring models, `computeWeightsArgs` (`w`) and
#'    `computeCovariateBalanceArgs` (`b`). Only for analyses with `computeCovariateBalanceArgs` and a
#'    strategy with a grace period. With shared censoring models it is shared by all outcomes;
#'    otherwise it is computed for outcomes of interest only.
#' 6. **Outcome models** (`Analysis_<a>/om_t<t>_tr<tr>[_c<c>]_o<o>.rds`): [fitOutcomeModel()], or
#'    [bootstrapCcwAnalysis()] when `fitOutcomeModelArgs$varianceEstimator = "bootstrap"`, on the
#'    weighted clones of the outcome's study population. The saved `CcwOutcomeModel` also holds
#'    `diagnostics`: the [computeMdrr()] result, [computePositivity()] of the censoring models, and,
#'    for outcomes of interest, [computeWeightSummary()] of the weighted clones and `curves`, from
#'    [computeCumulativeIncidence()] at every event time, `attrition` (that of the clones) and
#'    `followUpDistribution` (quantiles of the days of follow-up of the clones of each strategy).
#'    `diagnostics$values` holds the values that the summary compares with the diagnostic thresholds.
#' 7. **Summary** (`resultsSummary.rds`, `diagnosticsSummary.rds`): see [getResultsSummary()] and
#'    [getDiagnosticsSummary()].
#'
#' File names hold the IDs of the settings (`l`, `s`, `m`, `w`, `b`, numbered in order of first use)
#' and the cohort IDs (`t` target, `tr` treatment, `c` comparator treatment, `o` outcome). Every stage
#' computes only the files that do not exist yet, so an interrupted run resumes where it stopped, and
#' deleting a file recomputes it. Files are written under a temporary name and then renamed, so a file
#' that exists is complete. The summary is recomputed when an outcome model was (re)computed or when
#' it is missing.
#'
#' The specifications are saved in `outputFolder` (`ccwAnalysesSpecifications.rds`). If the folder
#' already holds outputs of different specifications, then in an interactive session the function asks
#' whether to delete them (answering no keeps them, and cancelling returns without running); outside an
#' interactive session, as when called by a Strategus module, it deletes them with a message and
#' recomputes everything. CohortMethod's `runCmAnalyses()` asks with [utils::askYesNo()] in both cases,
#' whose default answer outside an interactive session is also to delete.
#'
#' Thread counts come only from `multiThreadingSettings`. `censoringModelCvThreads` and
#' `outcomeCvThreads` set the `threads` of the Cyclops control of the censoring and outcome models,
#' which does not change the fits. Stages run in clusters made by [ParallelLogger::makeCluster()]; a
#' stage with 1 thread runs in the main R session. With an open `connection`, the data are extracted in
#' the main R session whatever `getDbCcwDataThreads` is.
#'
#' After completion, a tibble with the files of each combination is available from
#' [getFileReference()], the results from [getResultsSummary()], and the diagnostics from
#' [getDiagnosticsSummary()].
#'
#' @param connectionDetails         An R object of type `connectionDetails` created using the
#'                                  [DatabaseConnector::createConnectionDetails()] function. Either
#'                                  `connectionDetails` or `connection` must be provided.
#' @param cdmDatabaseSchema         The name of the database schema that contains the OMOP CDM instance.
#'                                  Requires read permissions to this database. On SQL Server, this
#'                                  should specify both the database and the schema, so for example
#'                                  'cdm_instance.dbo'.
#' @param tempEmulationSchema       Some database platforms like Oracle and Impala do not truly support
#'                                  temp tables. To emulate temp tables, provide a schema with write
#'                                  privileges where temp tables can be created.
#' @param exposureDatabaseSchema    The name of the database schema that contains the eligible (target)
#'                                  and treatment cohorts.
#' @param exposureTable             The name of the table that contains the eligible (target) and
#'                                  treatment cohorts, in the format of the CDM `cohort` table.
#' @param outcomeDatabaseSchema     The name of the database schema that contains the outcome cohorts.
#' @param outcomeTable              The name of the table that contains the outcome cohorts, in the
#'                                  format of the CDM `cohort` table.
#' @param outputFolder              Name of the folder where all the outputs will be written to.
#' @param multiThreadingSettings    An object of type `CcwMultiThreadingSettings` as created using
#'                                  [createCcwMultiThreadingSettings()] or
#'                                  [createDefaultCcwMultiThreadingSettings()].
#' @param ccwAnalysesSpecifications An object of type `CcwAnalysesSpecifications` as created using
#'                                  [createCcwAnalysesSpecifications()].
#' @param connection                An open connection created using [DatabaseConnector::connect()],
#'                                  instead of `connectionDetails`. It is not closed when the function
#'                                  returns.
#'
#' @return
#' A tibble with, for each analysis and target-treatment-outcome combination, the files of each stage
#' relative to `outputFolder` (`""` when the stage does not apply), as [getFileReference()] returns.
#'
#' @examples
#' \dontrun{
#' connectionDetails <- Eunomia::getEunomiaConnectionDetails()
#' Eunomia::createCohorts(connectionDetails)
#' analysis <- createCcwAnalysis(
#'   analysisId = 1,
#'   getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(
#'     covariateSettings = FeatureExtraction::createCovariateSettings(
#'       useDemographicsGender = TRUE,
#'       useDemographicsAgeGroup = TRUE,
#'       endDays = -1
#'     )
#'   ),
#'   createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
#'     riskWindowEnd = 365,
#'     endAnchor = "cohort start"
#'   ),
#'   createClonesArgs = createCreateClonesArgs(
#'     strategies = list(
#'       createTreatmentStrategy(1, "Celecoxib within 30 days", "initiateWithin", graceEnd = 30),
#'       createTreatmentStrategy(2, "No celecoxib within 30 days", "notInitiateWithin", graceEnd = 30)
#'     )
#'   ),
#'   fitCensoringModelsArgs = createFitCensoringModelsArgs(),
#'   computeWeightsArgs = createComputeWeightsArgs(),
#'   computeCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(),
#'   fitOutcomeModelArgs = createFitOutcomeModelArgs()
#' )
#' specifications <- createCcwAnalysesSpecifications(
#'   ccwAnalysisList = list(analysis),
#'   targetTreatmentOutcomesList = list(createTargetTreatmentOutcomes(
#'     targetId = 4,
#'     treatmentId = 1,
#'     comparatorTreatmentId = 2,
#'     outcomes = list(CohortMethod::createOutcome(outcomeId = 3))
#'   ))
#' )
#' runCcwAnalyses(
#'   connectionDetails = connectionDetails,
#'   cdmDatabaseSchema = "main",
#'   outputFolder = "output",
#'   multiThreadingSettings = createDefaultCcwMultiThreadingSettings(4),
#'   ccwAnalysesSpecifications = specifications
#' )
#' getResultsSummary("output")
#' }
#'
#' @export
runCcwAnalyses <- function(connectionDetails = NULL,
                           cdmDatabaseSchema,
                           tempEmulationSchema = getOption("sqlRenderTempEmulationSchema"),
                           exposureDatabaseSchema = cdmDatabaseSchema,
                           exposureTable = "cohort",
                           outcomeDatabaseSchema = cdmDatabaseSchema,
                           outcomeTable = "cohort",
                           outputFolder = "./CloneCensorWeightOutput",
                           multiThreadingSettings = createCcwMultiThreadingSettings(),
                           ccwAnalysesSpecifications,
                           connection = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertClass(connectionDetails, "ConnectionDetails", null.ok = TRUE, add = errorMessages)
  checkmate::assertClass(connection, "DatabaseConnectorConnection", null.ok = TRUE, add = errorMessages)
  checkmate::assertCharacter(cdmDatabaseSchema, len = 1, add = errorMessages)
  checkmate::assertCharacter(tempEmulationSchema, len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::assertCharacter(exposureDatabaseSchema, len = 1, add = errorMessages)
  checkmate::assertCharacter(exposureTable, len = 1, add = errorMessages)
  checkmate::assertCharacter(outcomeDatabaseSchema, len = 1, add = errorMessages)
  checkmate::assertCharacter(outcomeTable, len = 1, add = errorMessages)
  checkmate::assertCharacter(outputFolder, len = 1, add = errorMessages)
  checkmate::assertClass(multiThreadingSettings, "CcwMultiThreadingSettings", add = errorMessages)
  checkmate::assertR6(ccwAnalysesSpecifications, "CcwAnalysesSpecifications", add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (is.null(connectionDetails) == is.null(connection)) {
    stop("Provide either connectionDetails or connection, but not both")
  }
  specifications <- ccwAnalysesSpecifications
  on.exit(.clearCache())

  outputFolder <- normalizePath(outputFolder, mustWork = FALSE)
  specificationsFile <- file.path(outputFolder, "ccwAnalysesSpecifications.rds")
  if (file.exists(specificationsFile)) {
    oldSpecifications <- readRDS(specificationsFile)
    if (!isTRUE(all.equal(oldSpecifications$toList(), specifications$toList()))) {
      .clearCache()
      message(sprintf("Output files already exist in '%s', but the analysis settings have changed.", outputFolder))
      if (interactive()) {
        response <- utils::askYesNo("Do you want to delete the old files before proceeding?")
        if (is.na(response)) {
          return(invisible(NULL))
        } else if (response) {
          unlink(outputFolder, recursive = TRUE)
        }
      } else {
        message("Deleting the old files, because this is not an interactive session, and recomputing all outputs.")
        unlink(outputFolder, recursive = TRUE)
      }
    }
  }
  if (!dir.exists(outputFolder)) {
    dir.create(outputFolder, recursive = TRUE)
  }
  saveRDS(specifications, specificationsFile)
  referenceTable <- .createCcwReferenceTable(specifications, outputFolder)
  saveRDS(.publicReferenceTable(referenceTable), file.path(outputFolder, "outcomeModelReference.rds"))
  toCompute <- function(files) {
    files != "" & !file.exists(file.path(outputFolder, files))
  }
  getAnalysis <- function(analysisId) {
    ParallelLogger::matchInList(specifications$ccwAnalysisList, list(analysisId = analysisId))[[1]]
  }
  path <- function(file) {
    if (file == "") "" else file.path(outputFolder, file)
  }

  # Data ------------------------------------------------------------------------------------------
  subset <- referenceTable[!duplicated(referenceTable$ccwDataFile) & toCompute(referenceTable$ccwDataFile), ]
  if (nrow(subset) > 0) {
    message("*** Creating CloneCensorWeightData objects ***")
    tasks <- lapply(seq_len(nrow(subset)), function(i) {
      refRow <- subset[i, ]
      getDbArgs <- getAnalysis(refRow$analysisId)$getDbCloneCensorWeightDataArgs$clone()
      getDbArgs$covariateSettings <- .addCovariateConceptIds(getDbArgs$covariateSettings, refRow)
      outcomeIds <- unique(referenceTable$outcomeId[referenceTable$ccwDataFile == refRow$ccwDataFile])
      list(
        args = list(
          connectionDetails = connectionDetails,
          cdmDatabaseSchema = cdmDatabaseSchema,
          tempEmulationSchema = tempEmulationSchema,
          targetId = refRow$targetId,
          treatmentId = refRow$treatmentId,
          comparatorTreatmentId = if (is.na(refRow$comparatorTreatmentId)) NULL else refRow$comparatorTreatmentId,
          outcomeIds = outcomeIds,
          exposureDatabaseSchema = exposureDatabaseSchema,
          exposureTable = exposureTable,
          outcomeDatabaseSchema = outcomeDatabaseSchema,
          outcomeTable = outcomeTable,
          getDbCloneCensorWeightDataArgs = getDbArgs
        ),
        ccwDataFile = path(refRow$ccwDataFile)
      )
    })
    if (is.null(connection)) {
      .runTasks(tasks, .doCreateCcwData, multiThreadingSettings$getDbCcwDataThreads)
    } else {
      for (task in tasks) {
        task$args$connectionDetails <- NULL
        task$args$connection <- connection
        .doCreateCcwData(task)
      }
    }
  }

  # Study populations -------------------------------------------------------------------------------
  subset <- referenceTable[!duplicated(referenceTable$studyPopFile) & toCompute(referenceTable$studyPopFile), ]
  if (nrow(subset) > 0) {
    message("*** Creating study populations ***")
    tasks <- lapply(seq_len(nrow(subset)), function(i) {
      refRow <- subset[i, ]
      list(
        ccwDataFile = path(refRow$ccwDataFile),
        outcomeId = refRow$outcomeId,
        createStudyPopulationArgs = .getStudyPopulationArgs(specifications, getAnalysis(refRow$analysisId), refRow),
        studyPopFile = path(refRow$studyPopFile)
      )
    })
    .runTasks(tasks, .doCreateStudyPopulation, multiThreadingSettings$createStudyPopThreads)
  }

  # Censoring models -----------------------------------------------------------------------------
  # One task per censoring population, which fits the models of all its analyses, so a model shared by
  # strategies of several analyses is fit once.
  subset <- referenceTable[!duplicated(referenceTable$censoringModelFile) & toCompute(referenceTable$censoringModelFile), ]
  if (nrow(subset) > 0) {
    message("*** Fitting censoring models ***")
    groups <- split(seq_len(nrow(subset)), subset$censoringGroup)
    tasks <- lapply(groups, function(rows) {
      refRow <- subset[rows[1], ]
      groupFiles <- unique(referenceTable$censoringModelFile[referenceTable$censoringGroup == refRow$censoringGroup])
      list(
        ccwDataFile = path(refRow$ccwDataFile),
        outcomeId = if (specifications$refitCensoringModelsForEveryOutcome) refRow$outcomeId else NA_real_,
        createStudyPopulationArgs = .getStudyPopulationArgs(specifications, getAnalysis(refRow$analysisId), refRow),
        studyPopFile = if (specifications$refitCensoringModelsForEveryOutcome) path(refRow$studyPopFile) else "",
        existingModelFiles = vapply(groupFiles[!toCompute(groupFiles)], path, character(1)),
        models = lapply(rows, function(row) {
          analysis <- getAnalysis(subset$analysisId[row])
          fitCensoringModelsArgs <- analysis$fitCensoringModelsArgs$clone()
          fitCensoringModelsArgs$control$threads <- multiThreadingSettings$censoringModelCvThreads
          list(
            createClonesArgs = analysis$createClonesArgs,
            fitCensoringModelsArgs = fitCensoringModelsArgs,
            censoringModelFile = path(subset$censoringModelFile[row])
          )
        })
      )
    })
    .runTasks(unname(tasks), .doFitCensoringModels, multiThreadingSettings$fitCensoringModelThreads)
  }

  # Covariate balance ----------------------------------------------------------------------------
  subset <- referenceTable[!duplicated(referenceTable$balanceFile) & toCompute(referenceTable$balanceFile), ]
  if (nrow(subset) > 0) {
    message("*** Computing covariate balance ***")
    tasks <- lapply(seq_len(nrow(subset)), function(i) {
      refRow <- subset[i, ]
      analysis <- getAnalysis(refRow$analysisId)
      list(
        ccwDataFile = path(refRow$ccwDataFile),
        outcomeId = if (specifications$refitCensoringModelsForEveryOutcome) refRow$outcomeId else NA_real_,
        createStudyPopulationArgs = .getStudyPopulationArgs(specifications, analysis, refRow),
        studyPopFile = if (specifications$refitCensoringModelsForEveryOutcome) path(refRow$studyPopFile) else "",
        createClonesArgs = analysis$createClonesArgs,
        censoringModelFile = path(refRow$censoringModelFile),
        computeWeightsArgs = analysis$computeWeightsArgs,
        computeCovariateBalanceArgs = analysis$computeCovariateBalanceArgs,
        balanceFile = path(refRow$balanceFile)
      )
    })
    .runTasks(tasks, .doComputeBalance, multiThreadingSettings$computeBalanceThreads)
  }

  # Outcome models ---------------------------------------------------------------------------------
  subset <- referenceTable[toCompute(referenceTable$outcomeModelFile), ]
  subset <- subset[order(!subset$outcomeOfInterest), ]
  if (nrow(subset) > 0) {
    message("*** Fitting outcome models ***")
    tasks <- lapply(seq_len(nrow(subset)), function(i) {
      refRow <- subset[i, ]
      analysis <- getAnalysis(refRow$analysisId)
      fitOutcomeModelArgs <- analysis$fitOutcomeModelArgs$clone()
      fitOutcomeModelArgs$control$threads <- multiThreadingSettings$outcomeCvThreads
      fitCensoringModelsArgs <- analysis$fitCensoringModelsArgs$clone()
      fitCensoringModelsArgs$control$threads <- multiThreadingSettings$censoringModelCvThreads
      list(
        ccwDataFile = path(refRow$ccwDataFile),
        outcomeId = refRow$outcomeId,
        outcomeOfInterest = refRow$outcomeOfInterest,
        createStudyPopulationArgs = .getStudyPopulationArgs(specifications, analysis, refRow),
        studyPopFile = path(refRow$studyPopFile),
        createClonesArgs = analysis$createClonesArgs,
        fitCensoringModelsArgs = fitCensoringModelsArgs,
        censoringModelFile = path(refRow$censoringModelFile),
        computeWeightsArgs = analysis$computeWeightsArgs,
        fitOutcomeModelArgs = fitOutcomeModelArgs,
        positivityThreshold = specifications$ccwDiagnosticThresholds$positivityThreshold,
        outcomeModelFile = path(refRow$outcomeModelFile)
      )
    })
    .runTasks(tasks, .doFitOutcomeModel, multiThreadingSettings$fitOutcomeModelThreads)
  }

  # Summary --------------------------------------------------------------------------------------
  mainFileName <- file.path(outputFolder, "resultsSummary.rds")
  diagnosticsSummaryFileName <- file.path(outputFolder, "diagnosticsSummary.rds")
  if (nrow(subset) > 0 || !file.exists(mainFileName) || !file.exists(diagnosticsSummaryFileName)) {
    message("*** Summarizing results ***")
    summarizeResults(
      referenceTable = referenceTable,
      outputFolder = outputFolder,
      mainFileName = mainFileName,
      diagnosticsSummaryFileName = diagnosticsSummaryFileName,
      calibrationThreads = multiThreadingSettings$calibrationThreads,
      ccwDiagnosticThresholds = specifications$ccwDiagnosticThresholds
    )
  }
  invisible(.publicReferenceTable(referenceTable))
}

# Reference table --------------------------------------------------------------------------------

.f <- function(x) {
  return(format(x, scientific = FALSE, trim = TRUE))
}

# "_t<targetId>_tr<treatmentId>[_c<comparatorTreatmentId>]"
.createTtSubstring <- function(targetId, treatmentId, comparatorTreatmentId) {
  return(paste0(
    sprintf("_t%s_tr%s", .f(targetId), .f(treatmentId)),
    ifelse(is.na(comparatorTreatmentId), "", sprintf("_c%s", .f(comparatorTreatmentId)))
  ))
}

# The strategy settings the censoring models depend on (not the labels or strategy IDs).
.censoringSettingsJson <- function(analysis) {
  createClonesArgs <- analysis$createClonesArgs
  fitCensoringModelsArgs <- analysis$fitCensoringModelsArgs$toList()
  fitCensoringModelsArgs$control$threads <- NULL
  return(.settingsJson(list(
    strategies = lapply(createClonesArgs$strategies, function(strategy) {
      list(
        type = strategy$type,
        treatment = strategy$treatment,
        graceStart = strategy$graceStart,
        graceEnd = strategy$graceEnd,
        prohibitedTreatments = strategy$prohibitedTreatments
      )
    }),
    sameDayOrdering = createClonesArgs$sameDayOrdering,
    competingEvents = createClonesArgs$competingEvents,
    fitCensoringModelsArgs = fitCensoringModelsArgs
  )))
}

# The study population settings of an analysis for an outcome: those of the analysis, with the
# outcome's own priorOutcomeLookback and risk window when it has them (as in CohortMethod).
.getStudyPopulationArgs <- function(specifications, analysis, refRow) {
  createStudyPopulationArgs <- analysis$createStudyPopulationArgs$clone()
  outcome <- specifications$targetTreatmentOutcomesList[[refRow$ttoIndex]]$outcomes[[refRow$outcomeIndex]]
  for (field in c("priorOutcomeLookback", "riskWindowStart", "startAnchor", "riskWindowEnd", "endAnchor")) {
    if (!is.null(outcome[[field]])) {
      createStudyPopulationArgs[[field]] <- outcome[[field]]
    }
  }
  return(createStudyPopulationArgs)
}

.addCovariateConceptIds <- function(covariateSettings, refRow) {
  isList <- !is(covariateSettings, "covariateSettings")
  if (!isList) {
    covariateSettings <- list(covariateSettings)
  }
  conceptIds <- function(x) {
    x <- as.character(x)
    if (x == "") numeric(0) else as.numeric(strsplit(x, ",")[[1]])
  }
  for (i in seq_along(covariateSettings)) {
    covariateSettings[[i]]$excludedCovariateConceptIds <- unique(c(
      conceptIds(refRow$excludedCovariateConceptIds),
      covariateSettings[[i]]$excludedCovariateConceptIds
    ))
    covariateSettings[[i]]$includedCovariateConceptIds <- unique(c(
      conceptIds(refRow$includedCovariateConceptIds),
      covariateSettings[[i]]$includedCovariateConceptIds
    ))
  }
  if (!isList) {
    covariateSettings <- covariateSettings[[1]]
  }
  return(covariateSettings)
}

.createCcwReferenceTable <- function(specifications, outputFolder) {
  analysisList <- specifications$ccwAnalysisList
  refit <- specifications$refitCensoringModelsForEveryOutcome
  analyses <- dplyr::bind_rows(lapply(analysisList, function(analysis) {
    tibble(
      analysisId = as.integer(analysis$analysisId),
      analysisFolder = sprintf("Analysis_%d", analysis$analysisId),
      sameDayOrdering = analysis$createClonesArgs$sameDayOrdering,
      competingEvents = analysis$createClonesArgs$competingEvents,
      loadArgs = .settingsJson(analysis$getDbCloneCensorWeightDataArgs$toList()),
      censoringArgs = .censoringSettingsJson(analysis),
      weightsArgs = .settingsJson(analysis$computeWeightsArgs$toList()),
      balanceArgs = if (is.null(analysis$computeCovariateBalanceArgs) ||
        all(vapply(analysis$createClonesArgs$strategies, function(strategy) is.null(strategy$graceEnd), logical(1)))) {
        ""
      } else {
        .settingsJson(analysis$computeCovariateBalanceArgs$toList())
      }
    )
  }))
  ttos <- dplyr::bind_rows(lapply(seq_along(specifications$targetTreatmentOutcomesList), function(ttoIndex) {
    tto <- specifications$targetTreatmentOutcomesList[[ttoIndex]]
    dplyr::bind_rows(lapply(seq_along(tto$outcomes), function(outcomeIndex) {
      outcome <- tto$outcomes[[outcomeIndex]]
      tibble(
        targetId = as.numeric(tto$targetId),
        treatmentId = as.numeric(tto$treatmentId),
        comparatorTreatmentId = if (is.null(tto$comparatorTreatmentId)) NA_real_ else as.numeric(tto$comparatorTreatmentId),
        outcomeId = as.numeric(outcome$outcomeId),
        outcomeOfInterest = outcome$outcomeOfInterest,
        trueEffectSize = if (is.null(outcome$trueEffectSize)) NA_real_ else as.numeric(outcome$trueEffectSize),
        includedCovariateConceptIds = paste(tto$includedCovariateConceptIds, collapse = ","),
        excludedCovariateConceptIds = paste(tto$excludedCovariateConceptIds, collapse = ","),
        ttoIndex = ttoIndex,
        outcomeIndex = outcomeIndex
      )
    }))
  }))
  referenceTable <- dplyr::cross_join(analyses, ttos)

  # The study population settings of each row, and those of its censoring population: without an
  # outcome, the settings about prior outcomes do not matter.
  studyPopArgs <- character(nrow(referenceTable))
  censoringPopArgs <- character(nrow(referenceTable))
  for (i in seq_len(nrow(referenceTable))) {
    analysis <- analysisList[[match(referenceTable$analysisId[i], analyses$analysisId)]]
    argsList <- .getStudyPopulationArgs(specifications, analysis, referenceTable[i, ])$toList()
    studyPopArgs[i] <- .settingsJson(argsList)
    argsList$removeSubjectsWithPriorOutcome <- NULL
    argsList$priorOutcomeLookback <- NULL
    censoringPopArgs[i] <- .settingsJson(argsList)
  }
  id <- function(x) match(x, unique(x))
  referenceTable$loadId <- id(referenceTable$loadArgs)
  referenceTable$studyPopId <- id(studyPopArgs)
  referenceTable$censoringPopId <- if (refit) referenceTable$studyPopId else id(censoringPopArgs)
  referenceTable$censoringId <- id(referenceTable$censoringArgs)
  referenceTable$weightsId <- id(referenceTable$weightsArgs)
  referenceTable$balanceId <- id(referenceTable$balanceArgs)

  tt <- .createTtSubstring(referenceTable$targetId, referenceTable$treatmentId, referenceTable$comparatorTreatmentId)
  o <- sprintf("_o%s", .f(referenceTable$outcomeId))
  censoringO <- if (refit) o else ""
  referenceTable$ccwDataFile <- sprintf("CcwData_l%d%s.zip", referenceTable$loadId, tt)
  referenceTable$studyPopFile <- ifelse(
    referenceTable$outcomeOfInterest,
    sprintf("StudyPop_l%d_s%d%s%s.rds", referenceTable$loadId, referenceTable$studyPopId, tt, o),
    ""
  )
  referenceTable$censoringModelFile <- sprintf(
    "CensModel_l%d_s%d_m%d%s%s.rds",
    referenceTable$loadId, referenceTable$censoringPopId, referenceTable$censoringId, tt, censoringO
  )
  referenceTable$balanceFile <- ifelse(
    referenceTable$balanceArgs == "" | (refit & !referenceTable$outcomeOfInterest),
    "",
    sprintf(
      "Balance_l%d_s%d_m%d_w%d_b%d%s%s.rds",
      referenceTable$loadId, referenceTable$censoringPopId, referenceTable$censoringId, referenceTable$weightsId,
      referenceTable$balanceId, tt, censoringO
    )
  )
  referenceTable$outcomeModelFile <- file.path(referenceTable$analysisFolder, sprintf("om%s%s.rds", tt, o))
  # The censoring models of one censoring population are fit by one task, because they share fits.
  # sameDayOrdering and competingEvents change the risk sets, so they are part of it:
  referenceTable$censoringGroup <- paste(
    referenceTable$ccwDataFile, referenceTable$censoringPopId, censoringO, referenceTable$sameDayOrdering,
    referenceTable$competingEvents
  )

  referenceTable <- referenceTable[order(
    referenceTable$analysisId,
    referenceTable$targetId,
    referenceTable$treatmentId,
    referenceTable$comparatorTreatmentId,
    referenceTable$outcomeId
  ), ]
  analysesToExclude <- specifications$analysesToExclude
  if (!is.null(analysesToExclude)) {
    matchingColumns <- intersect(
      colnames(analysesToExclude),
      c("targetId", "treatmentId", "comparatorTreatmentId", "outcomeId", "analysisId")
    )
    analysesToExclude <- as.data.frame(analysesToExclude)[, matchingColumns, drop = FALSE]
    analysesToExclude[] <- lapply(analysesToExclude, as.numeric)
    keys <- as.data.frame(referenceTable[, matchingColumns])
    keys[] <- lapply(keys, as.numeric)
    countBefore <- nrow(referenceTable)
    referenceTable <- referenceTable[is.na(dplyr::left_join(
      keys,
      dplyr::mutate(dplyr::distinct(analysesToExclude), excluded = TRUE),
      by = matchingColumns
    )$excluded), ]
    message(sprintf(
      "Removed %d of the %d target-treatment-outcome-analysis combinations as specified by the user.",
      countBefore - nrow(referenceTable),
      countBefore
    ))
  }
  folders <- file.path(outputFolder, unique(referenceTable$analysisFolder))
  for (folder in folders[!dir.exists(folders)]) {
    dir.create(folder)
  }
  return(referenceTable)
}

.publicReferenceTable <- function(referenceTable) {
  return(referenceTable[, c(
    "analysisId",
    "targetId",
    "treatmentId",
    "comparatorTreatmentId",
    "outcomeId",
    "outcomeOfInterest",
    "trueEffectSize",
    "ccwDataFile",
    "studyPopFile",
    "censoringModelFile",
    "balanceFile",
    "outcomeModelFile"
  )])
}

# Tasks ------------------------------------------------------------------------------------------

.runTasks <- function(tasks, fun, threads) {
  cluster <- ParallelLogger::makeCluster(min(length(tasks), threads))
  on.exit(ParallelLogger::stopCluster(cluster))
  ParallelLogger::clusterRequire(cluster, "CloneCensorWeight")
  dummy <- ParallelLogger::clusterApply(cluster, tasks, fun)
  invisible(NULL)
}

# Writes under a temporary name and renames, so that a file that exists is complete.
.saveFile <- function(object, file, save = saveRDS) {
  tempFile <- paste0(file, ".tmp", Sys.getpid(), ".", sub(".*\\.", "", basename(file)))
  save(object, tempFile)
  if (!file.rename(tempFile, file)) {
    unlink(tempFile)
    stop(sprintf("Could not write '%s'", file))
  }
  invisible(NULL)
}

.clearCache <- function() {
  if (!is.null(cache$ccwData) && Andromeda::isValidAndromeda(cache$ccwData)) {
    close(cache$ccwData)
  }
  rm(list = ls(envir = cache), envir = cache)
}

.getCcwData <- function(ccwDataFile) {
  if (!identical(cache$ccwDataFile, ccwDataFile) || !Andromeda::isValidAndromeda(cache$ccwData)) {
    .clearCache()
    cache$ccwData <- loadCloneCensorWeightData(ccwDataFile)
    cache$ccwDataFile <- ccwDataFile
  }
  return(cache$ccwData)
}

# The study population of an outcome (outcomeId NA: without an outcome), from its file when it has one.
.getPopulation <- function(data, outcomeId, createStudyPopulationArgs, studyPopFile) {
  if (studyPopFile != "" && file.exists(studyPopFile)) {
    return(readRDS(studyPopFile))
  }
  if (is.na(outcomeId)) {
    return(.createStudyPopulationWithoutOutcome(data, createStudyPopulationArgs))
  }
  return(createStudyPopulation(data, outcomeId = outcomeId, createStudyPopulationArgs = createStudyPopulationArgs))
}

.doCreateCcwData <- function(params) {
  ParallelLogger::logDebug(sprintf(
    "Calling getDbCloneCensorWeightData() for targetId %s, treatmentId %s",
    params$args$targetId, params$args$treatmentId
  ))
  data <- do.call(getDbCloneCensorWeightData, params$args)
  # Saving closes the data.
  on.exit(if (Andromeda::isValidAndromeda(data)) close(data))
  .saveFile(data, params$ccwDataFile, save = function(x, file) suppressMessages(saveCloneCensorWeightData(x, file)))
  return(NULL)
}

.doCreateStudyPopulation <- function(params) {
  data <- .getCcwData(params$ccwDataFile)
  population <- createStudyPopulation(data, outcomeId = params$outcomeId, createStudyPopulationArgs = params$createStudyPopulationArgs)
  .saveFile(population, params$studyPopFile)
  return(NULL)
}

# Fits the censoring models of the analyses of one censoring population, reusing the fits of the models
# already saved for it and those of earlier analyses of the task.
.doFitCensoringModels <- function(params) {
  data <- .getCcwData(params$ccwDataFile)
  population <- .getPopulation(data, params$outcomeId, params$createStudyPopulationArgs, params$studyPopFile)
  fits <- list()
  for (file in params$existingModelFiles) {
    existingFits <- .censoringFitsFromModels(readRDS(file))
    fits[names(existingFits)] <- existingFits
  }
  for (model in params$models) {
    clones <- createClones(population, data, model$createClonesArgs)
    result <- tryCatch(
      .fitCensoringModelsWithFits(clones, data, model$fitCensoringModelsArgs, fits = fits),
      finally = close(clones)
    )
    fits <- result$fits
    .saveFile(result$models, model$censoringModelFile)
  }
  return(NULL)
}

.doComputeBalance <- function(params) {
  data <- .getCcwData(params$ccwDataFile)
  population <- .getPopulation(data, params$outcomeId, params$createStudyPopulationArgs, params$studyPopFile)
  clones <- createClones(population, data, params$createClonesArgs)
  on.exit(close(clones))
  weightedClones <- computeCensoringWeights(clones, readRDS(params$censoringModelFile), params$computeWeightsArgs)
  on.exit(close(weightedClones), add = TRUE)
  balance <- computeCovariateBalance(weightedClones, data, params$computeCovariateBalanceArgs)
  .saveFile(balance, params$balanceFile)
  return(NULL)
}

.doFitOutcomeModel <- function(params) {
  data <- .getCcwData(params$ccwDataFile)
  createStudyPopulationArgs <- params$createStudyPopulationArgs
  population <- .getPopulation(data, params$outcomeId, createStudyPopulationArgs, params$studyPopFile)
  clones <- createClones(population, data, params$createClonesArgs)
  on.exit(close(clones))
  censoringModels <- readRDS(params$censoringModelFile)
  fitOutcomeModelArgs <- params$fitOutcomeModelArgs
  outcomeTimeBasis <- if (fitOutcomeModelArgs$modelType == "logistic") fitOutcomeModelArgs$timeBasis else NULL
  weightedClones <- computeCensoringWeights(clones, censoringModels, params$computeWeightsArgs, outcomeTimeBasis = outcomeTimeBasis)
  on.exit(close(weightedClones), add = TRUE)
  if (fitOutcomeModelArgs$varianceEstimator == "bootstrap") {
    outcomeModel <- bootstrapCcwAnalysis(
      cloneCensorWeightData = data,
      outcomeId = params$outcomeId,
      createStudyPopulationArgs = createStudyPopulationArgs,
      createClonesArgs = params$createClonesArgs,
      fitCensoringModelsArgs = params$fitCensoringModelsArgs,
      computeWeightsArgs = params$computeWeightsArgs,
      fitOutcomeModelArgs = fitOutcomeModelArgs,
      censoringModels = censoringModels
    )
  } else {
    outcomeModel <- fitOutcomeModel(weightedClones, fitOutcomeModelArgs)
  }
  outcomeModel$diagnostics <- .computeOutcomeDiagnostics(
    clones = clones,
    weightedClones = weightedClones,
    censoringModels = censoringModels,
    fitOutcomeModelArgs = fitOutcomeModelArgs,
    positivityThreshold = params$positivityThreshold,
    summarizeWeights = params$outcomeOfInterest
  )
  if (params$outcomeOfInterest) {
    outcomeModel$curves <- computeCumulativeIncidence(weightedClones)
    outcomeModel$attrition <- attr(clones, "metaData")$attrition
    outcomeModel$followUpDistribution <- .computeFollowUpDistribution(clones)
  }
  .saveFile(outcomeModel, params$outcomeModelFile)
  return(NULL)
}

# The quantiles of the days of follow-up (survivalTime) of the clones of each strategy.
.computeFollowUpDistribution <- function(clones) {
  followUp <- clones$clones |>
    select("strategyId", "survivalTime") |>
    collect()
  return(dplyr::bind_rows(lapply(sort(unique(followUp$strategyId)), function(strategyId) {
    quantiles <- stats::quantile(followUp$survivalTime[followUp$strategyId == strategyId], c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1), names = FALSE)
    tibble(
      strategyId = as.integer(strategyId),
      minDays = quantiles[1],
      p10Days = quantiles[2],
      p25Days = quantiles[3],
      medianDays = quantiles[4],
      p75Days = quantiles[5],
      p90Days = quantiles[6],
      maxDays = quantiles[7]
    )
  })))
}

# The diagnostics of one outcome model, and the values compared with the CcwDiagnosticThresholds:
# - mdrr: from computeMdrr() on the clones;
# - positivity: the smallest predicted probability of following a strategy (the probability of
#   starting by the end of the horizon for "initiateWithin", of not starting for the others), over the
#   censoring models of both strategies, and positivityShareExtreme, the largest shareExtreme of
#   computePositivity();
# - with summarizeWeights, from computeWeightSummary(): essFraction, the smallest effectiveSampleShare,
#   and maxWeight, the largest weight, after truncation, and maxWeightBeforeTruncation, over the strategies,
#   bands of days, and the scales of the weights the estimates use (that of computeCensoringWeights(),
#   used by Kaplan-Meier and Aalen-Johansen, and msmWeights for the Cox and pooled logistic models).
.computeOutcomeDiagnostics <- function(clones,
                                       weightedClones,
                                       censoringModels,
                                       fitOutcomeModelArgs,
                                       positivityThreshold,
                                       summarizeWeights) {
  mdrr <- computeMdrr(clones, modelType = if (fitOutcomeModelArgs$modelType == "logistic") "logistic" else "cox")
  positivity <- computePositivity(clones, censoringModels, threshold = positivityThreshold)
  types <- vapply(attr(clones, "metaData")$strategies, function(strategy) strategy$type, character(1))
  following <- ifelse(types[positivity$strategyId] == "initiateWithin", positivity$min, 1 - positivity$max)
  shareExtreme <- if (nrow(positivity) == 0) NA_real_ else max(positivity$shareExtreme)
  values <- tibble(
    mdrr = as.numeric(mdrr$mdrr),
    positivity = if (length(following) == 0) NA_real_ else min(following),
    positivityShareExtreme = shareExtreme,
    essFraction = NA_real_,
    maxWeight = NA_real_,
    maxWeightBeforeTruncation = NA_real_
  )
  weightSummary <- NULL
  if (summarizeWeights) {
    weightSummary <- computeWeightSummary(weightedClones)
    weightsMetaData <- attr(weightedClones, "metaData")$weights
    scales <- if (weightsMetaData$stabilized) "stabilized" else "unstabilized"
    if (fitOutcomeModelArgs$modelType != "km") {
      scales <- unique(c(scales, fitOutcomeModelArgs$msmWeights))
    }
    used <- weightSummary[weightSummary$scale %in% scales, ]
    after <- used[used$truncation == "after", ]
    before <- used[used$truncation == "before", ]
    extreme <- function(x, f) {
      x <- x[!is.na(x)]
      if (length(x) == 0) NA_real_ else f(x)
    }
    values$essFraction <- extreme(after$effectiveSampleShare, min)
    values$maxWeight <- extreme(after$maxWeight, max)
    values$maxWeightBeforeTruncation <- extreme(before$maxWeight, max)
  }
  return(list(values = values, mdrr = mdrr, positivity = positivity, weightSummary = weightSummary))
}

# Getters ----------------------------------------------------------------------------------------

#' Get file reference
#'
#' @param outputFolder Name of the folder where all the outputs have been written to.
#'
#' @return
#' A tibble with, for each analysis and target-treatment-outcome combination, the files generated by
#' [runCcwAnalyses()] relative to `outputFolder` (`""` when a stage does not apply): `ccwDataFile`,
#' `studyPopFile`, `censoringModelFile`, `balanceFile` and `outcomeModelFile`, with the columns
#' `analysisId`, `targetId`, `treatmentId`, `comparatorTreatmentId` (`NA` without one), `outcomeId`,
#' `outcomeOfInterest` and `trueEffectSize`.
#'
#' @export
getFileReference <- function(outputFolder) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertCharacter(outputFolder, len = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  outputFolder <- normalizePath(outputFolder)
  return(readRDS(file.path(outputFolder, "outcomeModelReference.rds")))
}

#' Get a summary report of the analyses results
#'
#' @details
#' One row per analysis, target-treatment-outcome combination, and estimand of the outcome model
#' (`"hr"`, `"or"`, `"risk_target_<h>"`, `"risk_comparator_<h>"`, `"rd_<h>"`, `"rr_<h>"`,
#' `"rmst_<h>"`, `"rmtl_<h>"`; see [fitOutcomeModel()]), with the counts per strategy (the target
#' strategy is strategy 1, the comparator strategy 2): `targetSubjects` (persons), `targetClones`,
#' `targetDays` (days of follow-up), `targetOutcomes`, `targetCompetingEvents`, `targetWeightedDays`
#' and `targetWeightedOutcomes` (weighted by the weights of [computeCensoringWeights()]), and the same
#' for the comparator; and `estimate`, `ci95Lb`, `ci95Ub`, `p`, and `se` (additive estimands), or `logRr` and
#' `seLogRr` (ratio estimands), as in [fitOutcomeModel()]. Ratio estimands (`"hr"`, `"or"`, `"rr_<h>"`)
#' also have `oneSidedP` (for a ratio above 1) and, when there are at least 5 negative controls (outcomes
#' with `trueEffectSize = 1`) with a standard error whose diagnostics allow calibration (see
#' [getDiagnosticsSummary()]), the empirically calibrated `calibratedRr`, `calibratedCi95Lb`,
#' `calibratedCi95Ub`, `calibratedP`, `calibratedOneSidedP`, `calibratedLogRr` and
#' `calibratedSeLogRr`, as in CohortMethod: each estimand of each analysis and target-treatment
#' combination is calibrated with its own negative controls, using
#' [EmpiricalCalibration::fitMcmcNull()], and, with at least 5 positive controls, a systematic error
#' model for the confidence intervals.
#'
#' @param outputFolder Name of the folder where all the outputs have been written to.
#'
#' @return
#' A tibble containing the results of each target-treatment-outcome-analysis combination and estimand.
#'
#' @export
getResultsSummary <- function(outputFolder) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertCharacter(outputFolder, len = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  outputFolder <- normalizePath(outputFolder)
  return(readRDS(file.path(outputFolder, "resultsSummary.rds")))
}

#' Get a summary report of the analyses diagnostics
#'
#' @details
#' One row per analysis and target-treatment-outcome combination, with each diagnostic and its
#' `"PASS"`, `"FAIL"` or `"NOT EVALUATED"` against the thresholds of [createCcwDiagnosticThresholds()]:
#'
#' - `maxSdm`, `balanceDiagnostic`: the largest absolute standardized difference of mean after
#'   weighting, from [computeCovariateBalance()] (`< sdmThreshold`). Not evaluated without
#'   `computeCovariateBalanceArgs`, and, with `refitCensoringModelsForEveryOutcome = TRUE`, for outcomes
#'   not of interest.
#' - `mdrr`, `mdrrDiagnostic`: the minimum detectable relative risk from [computeMdrr()]
#'   (`< mdrrThreshold`).
#' - `ease`, `easeDiagnostic`: the expected absolute systematic error from the negative controls, the
#'   largest over the ratio estimands of the analysis and target-treatment combination
#'   (`< easeThreshold`). Not evaluated with fewer than 5 negative controls.
#' - `essFraction`, `essFractionDiagnostic`: the smallest effective sample size of a strategy as a
#'   fraction of its clones under follow-up, averaged over the days of a band, from
#'   [computeWeightSummary()] (`>= essFractionThreshold`).
#' - `maxWeight`, `maxWeightDiagnostic`: the largest weight the estimates use, after truncation
#'   (`< maxWeightThreshold`). `maxWeightBeforeTruncation` is the largest one before truncation.
#'   `essFraction` and `maxWeight` are taken over the strategies, the bands of days of
#'   [computeWeightSummary()], and the scales of the weights the estimates use. They are not evaluated
#'   for outcomes not of interest.
#' - `positivity`, `positivityDiagnostic`: the smallest predicted probability of following a strategy
#'   from the censoring models, from [computePositivity()] (`>= positivityThreshold`).
#'   `positivityShareExtreme` is the largest share of persons with a probability below the threshold.
#'
#' As in CohortMethod, the negative controls used for calibration are those whose balance, ESS
#' fraction, maximum weight and positivity diagnostics do not fail; `unblindForEvidenceSynthesis` is
#' that, and the EASE diagnostic not failing, and `unblind` also requires the MDRR diagnostic not to
#' fail.
#'
#' @param outputFolder Name of the folder where all the outputs have been written to.
#'
#' @return
#' A tibble containing the diagnostics of each target-treatment-outcome-analysis combination.
#'
#' @export
getDiagnosticsSummary <- function(outputFolder) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertCharacter(outputFolder, len = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  outputFolder <- normalizePath(outputFolder)
  return(readRDS(file.path(outputFolder, "diagnosticsSummary.rds")))
}

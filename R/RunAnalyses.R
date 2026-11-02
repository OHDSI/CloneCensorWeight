#' Run a full suite of CCW target trial emulation analyses
#'
#' Iterates over every combination of `TargetComparatorOutcome` × `CcwAnalysis` ×
#' outcome in `ccwAnalysesSpecifications`, executes the full clone-censor-weight
#' pipeline for each, and appends results to CSV files in `outputFolder`.
#'
#' **Pipeline per combination**
#' \enumerate{
#'   \item Extract data from the CDM via [getDbCloneCensorWeightData()].
#'   \item Apply study population restrictions via [createStudyPopulation()].
#'   \item Clone and censor via [createClonesAndCensor()].
#'   \item Discretize person-time via [discretizePersonTime()].
#'   \item Compute IPCW weights via [createCensoringWeights()].
#'   \item Fit the outcome model via [fitOutcomeModel()].
#'   \item Compute covariate balance via [computeStandardizedDifferences()].
#'   \item Compute weighted Kaplan-Meier via [computeKaplanMeier()].
#'   \item Compute MDRR via [computeMdrr()].
#'   \item Calibrate using negative control outcomes (when supplied).
#' }
#'
#' **Resume support** — If `outputFolder` already contains a
#' `ccw_outcome_model.csv` file, any composite analysis ID already present in
#' that file is skipped.
#'
#' **Composite analysis ID** — Each combination is uniquely identified by
#' `"{ccwAnalysis$analysisId}_t{targetId}_o{outcomeId}"`.
#'
#' @param connectionDetails A `ConnectionDetails` object from
#'   [DatabaseConnector::createConnectionDetails()].
#' @param cdmDatabaseSchema Schema containing the OMOP CDM tables.
#' @param tempEmulationSchema Schema for temp-table emulation (Oracle,
#'   BigQuery).  Defaults to `getOption("sqlRenderTempEmulationSchema")`.
#' @param targetDatabaseSchema Schema containing the target cohort table.
#' @param targetTable Target cohort table name.
#' @param outcomeDatabaseSchema Schema containing the outcome cohort table.
#' @param outcomeTable Outcome cohort table name.
#' @param exposureDatabaseSchema Schema containing the exposure cohort table.
#' @param exposureTable Exposure cohort table name.
#' @param ccwAnalysesSpecifications A [createCcwAnalysesSpecifications()]
#'   object defining the analyses to run.
#' @param outputFolder Path to the directory that will receive CSV results.
#'   Created if absent.
#' @param databaseId Short identifier for the data source stamped onto every
#'   result row.
#' @param databaseName Human-readable data source name.
#' @param negativeControlOutcomeIds Integer vector of outcome IDs (a subset of
#'   the `outcomeIds` in each `TargetComparatorOutcome`) to treat as negative controls
#'   for empirical calibration.  Leave empty to skip calibration.
#' @param cdmVersion OMOP CDM version string passed to FeatureExtraction.
#'
#' @return Invisibly, the path to `outputFolder`.
#' @export
runCcwAnalyses <- function(connectionDetails,
                           cdmDatabaseSchema,
                           tempEmulationSchema       = getOption("sqlRenderTempEmulationSchema"),
                           targetDatabaseSchema,
                           targetTable,
                           outcomeDatabaseSchema,
                           outcomeTable,
                           exposureDatabaseSchema,
                           exposureTable,
                           ccwAnalysesSpecifications,
                           outputFolder,
                           databaseId                = "unknown",
                           databaseName              = databaseId,
                           negativeControlOutcomeIds = integer(0),
                           cdmVersion                = "5") {
  checkmate::assertClass(connectionDetails, "ConnectionDetails")
  checkmate::assertString(cdmDatabaseSchema)
  checkmate::assertString(targetDatabaseSchema)
  checkmate::assertString(targetTable)
  checkmate::assertString(outcomeDatabaseSchema)
  checkmate::assertString(outcomeTable)
  checkmate::assertString(exposureDatabaseSchema)
  checkmate::assertString(exposureTable)
  checkmate::assertR6(ccwAnalysesSpecifications, "CcwAnalysesSpecifications")
  checkmate::assertString(outputFolder)
  checkmate::assertString(databaseId)
  checkmate::assertIntegerish(negativeControlOutcomeIds, any.missing = FALSE)

  if (!dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

  specs <- ccwAnalysesSpecifications
  nTao  <- length(specs$targetArmOutcomeList)
  nAna  <- length(specs$ccwAnalysisList)
  ParallelLogger::logInfo(sprintf(
    "Running %d target-arm-outcome(s) x %d analysis spec(s)",
    nTao, nAna
  ))

  for (taoIdx in seq_along(specs$targetArmOutcomeList)) {
    tao <- specs$targetArmOutcomeList[[taoIdx]]

    for (anaIdx in seq_along(specs$ccwAnalysisList)) {
      analysis <- specs$ccwAnalysisList[[anaIdx]]

      ParallelLogger::logInfo(sprintf(
        "TAO %d/%d, analysis '%s': extracting data",
        taoIdx, nTao, analysis$analysisId
      ))

      data <- tryCatch(
        getDbCloneCensorWeightData(
          connectionDetails       = connectionDetails,
          cdmDatabaseSchema       = cdmDatabaseSchema,
          tempEmulationSchema     = tempEmulationSchema,
          targetDatabaseSchema    = targetDatabaseSchema,
          targetTable             = targetTable,
          targetId                = tao$indexCohortId,
          outcomeDatabaseSchema   = outcomeDatabaseSchema,
          outcomeTable            = outcomeTable,
          outcomeIds              = tao$outcomeIds,
          exposureDatabaseSchema  = exposureDatabaseSchema,
          exposureTable           = exposureTable,
          exposureId              = tao$targetCohortId,
          getDbCloneCensorWeightDataArgs = analysis$getDbCloneCensorWeightDataArgs,
          cdmVersion              = cdmVersion
        ),
        error = function(e) {
          ParallelLogger::logError(sprintf(
            "Data extraction failed for TAO %d, analysis '%s': %s",
            taoIdx, analysis$analysisId, conditionMessage(e)
          ))
          NULL
        }
      )

      if (is.null(data)) next

      .runForData(
        data                      = data,
        tao                       = tao,
        analysis                  = analysis,
        outputFolder              = outputFolder,
        databaseId                = databaseId,
        databaseName              = databaseName,
        negativeControlOutcomeIds = negativeControlOutcomeIds
      )

      Andromeda::close(data)
    }
  }

  ParallelLogger::logInfo("runCcwAnalyses complete. Results written to: ", outputFolder)
  invisible(outputFolder)
}

#' Run CCW analyses on a pre-loaded data object
#'
#' A convenience wrapper around the internal pipeline that accepts an already
#' loaded `CloneCensorWeightData` object instead of database connection
#' details.  Useful for vignettes, simulation studies with [simulateCcwData()],
#' and integration testing without a live CDM connection.
#'
#' @param data A `CloneCensorWeightData` object, e.g. from
#'   [createSampleCloneCensorWeightData()] or [simulateCcwData()].
#' @inheritParams runCcwAnalyses
#'
#' @return Invisibly, the path to `outputFolder`.
#' @export
runCcwAnalysesOnData <- function(data,
                                 ccwAnalysesSpecifications,
                                 outputFolder,
                                 databaseId                = "unknown",
                                 databaseName              = databaseId,
                                 negativeControlOutcomeIds = integer(0)) {
  checkmate::assertClass(data, "CloneCensorWeightData")
  checkmate::assertR6(ccwAnalysesSpecifications, "CcwAnalysesSpecifications")
  checkmate::assertString(outputFolder)
  checkmate::assertIntegerish(negativeControlOutcomeIds, any.missing = FALSE)

  if (!dir.exists(outputFolder)) dir.create(outputFolder, recursive = TRUE)

  specs <- ccwAnalysesSpecifications
  for (tao in specs$targetArmOutcomeList) {
    for (analysis in specs$ccwAnalysisList) {
      .runForData(
        data                      = data,
        tao                       = tao,
        analysis                  = analysis,
        outputFolder              = outputFolder,
        databaseId                = databaseId,
        databaseName              = databaseName,
        negativeControlOutcomeIds = negativeControlOutcomeIds
      )
    }
  }
  invisible(outputFolder)
}

# ---------------------------------------------------------------------------
# Internal: run the full analysis pipeline on a pre-loaded data object.
# ---------------------------------------------------------------------------

.runForData <- function(data,
                        tao,
                        analysis,
                        outputFolder,
                        databaseId                = "unknown",
                        databaseName              = databaseId,
                        negativeControlOutcomeIds = integer(0)) {

  originalCohorts <- .readCloneCensorWeightTable(data, "cohorts")
  allOutcomes     <- .readCloneCensorWeightTable(data, "outcomes")

  # Accumulate results in memory so we can calibrate at the end of the block.
  outcomeResults <- list()

  for (outcomeId in tao$outcomeIds) {
    compositeId <- .makeCompositeId(analysis$analysisId, tao$indexCohortId, outcomeId)

    if (.alreadyRun(compositeId, outputFolder)) {
      ParallelLogger::logInfo(sprintf("Skipping '%s' (already in output)", compositeId))
      next
    }

    ParallelLogger::logInfo(sprintf("Running '%s'", compositeId))

    # Restore per-outcome state
    data[["cohorts"]]  <- originalCohorts
    data[["outcomes"]] <- allOutcomes[allOutcomes$outcomeId == outcomeId, , drop = FALSE]

    result <- tryCatch(
      .runOnePipeline(data, outcomeId, analysis),
      error = function(e) {
        ParallelLogger::logError(sprintf("Pipeline failed for '%s': %s",
                                         compositeId, conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(result)) {
      result$compositeId <- compositeId
      outcomeResults[[as.character(outcomeId)]] <- result
    }
  }

  # Restore the outcomes table to its original full state
  data[["outcomes"]] <- allOutcomes

  if (length(outcomeResults) == 0L) return(invisible(NULL))

  # Empirical calibration -------------------------------------------------------
  ncIds <- intersect(tao$outcomeIds, negativeControlOutcomeIds)
  ncDf  <- .buildNcDf(outcomeResults, ncIds)

  for (outcomeIdStr in names(outcomeResults)) {
    r           <- outcomeResults[[outcomeIdStr]]
    summaryDf   <- r$model$summary
    isNc        <- as.integer(outcomeIdStr) %in% negativeControlOutcomeIds

    if (!isNc && !is.null(ncDf) && nrow(ncDf) >= 5L) {
      summaryDf <- tryCatch(
        calibrateEstimates(r$model, negativeControlResults = ncDf),
        error = function(e) {
          ParallelLogger::logWarn(sprintf(
            "Calibration failed for '%s': %s", r$compositeId, conditionMessage(e)
          ))
          r$model$summary
        }
      )
    }

    analysisRow <- data.frame(
      analysis_id          = r$compositeId,
      description          = analysis$description,
      created_date         = Sys.Date(),
      package_version      = as.character(utils::packageVersion("CloneCensorWeight")),
      database_id          = databaseId,
      database_name        = databaseName,
      stringsAsFactors     = FALSE
    )

    exportToCsv(
      results = list(
        analysis               = analysisRow,
        outcomeModel           = summaryDf,
        standardizedDifference = r$smd,
        kaplanMeier            = r$km,
        powerMdrr              = r$mdrr
      ),
      analysisId   = r$compositeId,
      exportFolder = outputFolder,
      append       = TRUE
    )
  }

  invisible(NULL)
}

# Run the full pipeline for a single outcome on pre-prepared data.
.runOnePipeline <- function(data, outcomeId, analysis) {
  createStudyPopulation(data, outcomeId, analysis$createStudyPopulationArgs)

  intervalLength <- analysis$createClonesAndCensorArgs$intervalLength %||% 30L
  cloned   <- createClonesAndCensor(data, analysis$createClonesAndCensorArgs)
  pt       <- discretizePersonTime(cloned, intervalLength = intervalLength)
  weighted <- createCensoringWeights(pt, cloned, analysis$createCensoringWeightsArgs)

  model  <- fitOutcomeModel(weighted, fitOutcomeModelArgs = analysis$fitOutcomeModelArgs)

  smdBefore          <- computeStandardizedDifferences(data)
  smdAfter           <- computeStandardizedDifferences(data, weights = weighted)
  smdBefore$weighted <- 0L
  smdAfter$weighted  <- 1L
  smd <- rbind(smdBefore, smdAfter)

  km   <- computeKaplanMeier(weighted)
  mdrr <- computeMdrr(weighted)

  list(model = model, smd = smd, km = km, mdrr = mdrr)
}

# Build a negative-control data frame (logRr, seLogRr) from the results list.
.buildNcDf <- function(outcomeResults, ncIds) {
  rows <- lapply(ncIds, function(oid) {
    r <- outcomeResults[[as.character(oid)]]
    if (is.null(r)) return(NULL)
    s <- r$model$summary
    se <- (log(s$upper95) - log(s$lower95)) / (2 * stats::qnorm(0.975))
    data.frame(logRr = log(s$estimate), seLogRr = se, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(NULL)
  do.call(rbind, rows)
}

# Build the composite analysis ID string.
.makeCompositeId <- function(analysisId, indexCohortId, outcomeId) {
  paste0(analysisId, "_t", indexCohortId, "_o", outcomeId)
}

# TRUE if the composite ID already appears in ccw_outcome_model.csv.
.alreadyRun <- function(compositeId, outputFolder) {
  path <- file.path(outputFolder, "ccw_outcome_model.csv")
  if (!file.exists(path)) return(FALSE)
  # Read only the analysis_id column to stay fast on large files
  ids <- tryCatch(
    readr::read_csv(path, col_select = "analysis_id",
                    show_col_types = FALSE)$analysis_id,
    error = function(e) character(0)
  )
  compositeId %in% ids
}

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

# The export follows CohortMethod's exportToCsv() (R/Export.R, Copyright Observational Health Data
# Sciences and Informatics, Apache 2.0): the target-comparator hash and the tidying of the covariate
# balance are adapted from there, and the minimum cell count rule is its enforceMinCellValue().

#' Export clone-censor-weight results to CSV files
#'
#' @details
#' This requires that [runCcwAnalyses()] has been executed first. It exports the results in
#' `outputFolder` to one CSV file per table of the results data model (see
#' [getResultsDataModelSpecifications()]), with the columns of the specification in its order, and
#' adds them to `Results_<databaseId>.zip` in `exportFolder`, for sharing with other sites and for
#' [uploadResults()].
#'
#' - The definition tables: `ccw_analysis` and `ccw_strategy` (the JSON of each analysis and strategy),
#'   `ccw_target_treatment` (one row per target-treatment-comparator treatment combination, whose
#'   `target_treatment_id` is a hash of the three cohort IDs, as CohortMethod's
#'   `target_comparator_id`), and `ccw_target_treatment_outcome`.
#' - For every analysis and target-treatment-outcome combination: `ccw_result` (one row per
#'   estimand, see [getResultsSummary()]) and `ccw_diagnostics_summary` (see [getDiagnosticsSummary()]).
#'   Also `ccw_likelihood_profile`, as CohortMethod's `cm_likelihood_profile`: the likelihood profile of
#'   the log `hr`, `shr` or `or` (`logLikelihoodProfile` of [fitOutcomeModel()]), with `log_likelihood`
#'   normalized to a maximum of 0 and `gradient` as is, also for fits that did not converge, but not for
#'   contrasts that are null by design (standard error 0 up to rounding: below 1e-6). No minimum cell
#'   count applies to it.
#' - For outcomes of interest: `ccw_attrition`, `ccw_follow_up_dist`, `ccw_covariate_balance` (at the
#'   end of the grace periods, see [computeCovariateBalance()]), `ccw_censoring_model` (the non-zero
#'   coefficients) and `ccw_censoring_model_dist` (the density of the predicted probability of
#'   starting the treatment, from [computePositivity()]), `ccw_weight_dist` (from
#'   [computeWeightSummary()]) and `ccw_cumulative_incidence` (from [computeCumulativeIncidence()]).
#' - `ccw_covariate` and `ccw_covariate_analysis`: the covariates of each analysis.
#'
#' Strategies are identified by their `strategyId` (see [createTreatmentStrategy()]): tables per
#' strategy have `strategy_id`, and tables that compare the two strategies of an analysis have
#' `target_strategy_id` and `comparator_strategy_id`.
#'
#' To preserve privacy, counts in the columns flagged `minCellCount` in the specifications that are
#' greater than 0 and smaller than `minCellCount` are written as `-minCellCount`, as in CohortMethod.
#' Means of the covariate balance are censored likewise, at `minCellCount` divided by the number of
#' clones they are computed on.
#'
#' @param outputFolder The folder where [runCcwAnalyses()] generated all results.
#' @param exportFolder The folder where the CSV files and the zip file will be written.
#' @param databaseId   A unique ID for the database. It is written in the `database_id` column of every
#'                     site-specific table.
#' @param minCellCount To preserve privacy: the minimum count that can be included in the results. A
#'                     count below this threshold (and above 0) is written as `-minCellCount`.
#' @param maxCores     How many parallel cores may be used. Kept for the arguments of CohortMethod's
#'                     `exportToCsv()`; the export currently runs in one thread, because the curves
#'                     and distributions are computed by [runCcwAnalyses()].
#'
#' @return
#' The path of the zip file, invisibly. Is called for its side effect of writing the CSV files and the
#' zip file to `exportFolder`.
#'
#' @examples
#' \dontrun{
#' exportToCsv(
#'   outputFolder = "output",
#'   databaseId = "Eunomia",
#'   minCellCount = 5
#' )
#' }
#'
#' @export
exportToCsv <- function(outputFolder,
                        exportFolder = file.path(outputFolder, "export"),
                        databaseId,
                        minCellCount = 5,
                        maxCores = 1) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertCharacter(outputFolder, len = 1, add = errorMessages)
  checkmate::assertDirectoryExists(outputFolder, add = errorMessages)
  checkmate::assertFileExists(file.path(outputFolder, "ccwAnalysesSpecifications.rds"), add = errorMessages)
  checkmate::assertFileExists(file.path(outputFolder, "resultsSummary.rds"), add = errorMessages)
  checkmate::assertCharacter(exportFolder, len = 1, add = errorMessages)
  checkmate::assertAtomic(databaseId, len = 1, any.missing = FALSE, add = errorMessages)
  checkmate::assertInt(minCellCount, lower = 0, add = errorMessages)
  checkmate::assertInt(maxCores, lower = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (!dir.exists(exportFolder)) {
    dir.create(exportFolder, recursive = TRUE)
  }
  start <- Sys.time()
  message("Exporting results to CSV")
  databaseId <- as.character(databaseId)
  specifications <- readRDS(file.path(outputFolder, "ccwAnalysesSpecifications.rds"))
  targetTreatment <- .getTargetTreatment(specifications)
  strategyIds <- .getStrategyIds(specifications)
  reference <- getFileReference(outputFolder) |>
    dplyr::inner_join(targetTreatment, by = c("targetId", "treatmentId", "comparatorTreatmentId")) |>
    dplyr::inner_join(strategyIds, by = "analysisId")

  files <- c(
    .exportDefinitions(specifications, targetTreatment, exportFolder),
    .exportCovariates(outputFolder, reference, exportFolder, databaseId),
    .exportResults(outputFolder, targetTreatment, strategyIds, exportFolder, databaseId, minCellCount),
    .exportDiagnosticsSummary(outputFolder, targetTreatment, strategyIds, exportFolder, databaseId),
    .exportCovariateBalance(outputFolder, reference, specifications, exportFolder, databaseId, minCellCount),
    .exportCensoringModels(outputFolder, reference, exportFolder, databaseId),
    .exportFromOutcomeModels(outputFolder, reference, exportFolder, databaseId, minCellCount),
    .exportLikelihoodProfiles(outputFolder, reference, exportFolder, databaseId)
  )

  message("Adding results to zip file")
  zipName <- file.path(normalizePath(exportFolder), sprintf("Results_%s.zip", databaseId))
  if (file.exists(zipName)) {
    unlink(zipName)
  }
  DatabaseConnector::createZipFile(zipFile = zipName, files = file.path(exportFolder, files), rootFolder = exportFolder)
  delta <- Sys.time() - start
  message("Exporting to CSV took ", signif(delta, 3), " ", attr(delta, "units"))
  message("Results are ready for sharing at: ", zipName)
  invisible(zipName)
}

# Helpers --------------------------------------------------------------------------------------------

# A murmur32 hash of the three cohort IDs, as CohortMethod's target_comparator_id.
.createTargetTreatmentId <- function(targetId, treatmentId, comparatorTreatmentId) {
  hashString <- paste(targetId, treatmentId, comparatorTreatmentId)
  return(vapply(hashString, function(string) {
    as.integer(as.numeric(paste0("0x", digest::digest(string, algo = "murmur32", serialize = FALSE))) - 2^31)
  }, integer(1), USE.NAMES = FALSE))
}

.getTargetTreatment <- function(specifications) {
  targetTreatment <- dplyr::bind_rows(lapply(specifications$targetTreatmentOutcomesList, function(tto) {
    tibble(
      targetId = as.numeric(tto$targetId),
      treatmentId = as.numeric(tto$treatmentId),
      comparatorTreatmentId = if (is.null(tto$comparatorTreatmentId)) NA_real_ else as.numeric(tto$comparatorTreatmentId)
    )
  })) |>
    dplyr::distinct()
  targetTreatment$targetTreatmentId <- .createTargetTreatmentId(
    targetTreatment$targetId,
    targetTreatment$treatmentId,
    targetTreatment$comparatorTreatmentId
  )
  return(targetTreatment)
}

# The strategyId of the target (first) and comparator (second) strategy of each analysis.
.getStrategyIds <- function(specifications) {
  return(dplyr::bind_rows(lapply(specifications$ccwAnalysisList, function(analysis) {
    strategies <- analysis$createClonesArgs$strategies
    tibble(
      analysisId = as.integer(analysis$analysisId),
      targetStrategyId = as.integer(strategies[[1]]$strategyId),
      comparatorStrategyId = as.integer(strategies[[2]]$strategyId)
    )
  })))
}

# The strategyId of the strategies at positions 1 (target) and 2 (comparator) of a reference row.
.strategyIdOf <- function(refRow, position) {
  return(ifelse(position == 1, refRow$targetStrategyId, refRow$comparatorStrategyId))
}

# Writes a table with exactly the columns of its specification, in their order, and returns the file
# name. Integer-valued columns are written without scientific notation and logicals as 0/1.
.writeToCsv <- function(data, tableName, exportFolder) {
  specifications <- getResultsDataModelSpecifications()
  columns <- specifications$columnName[specifications$tableName == tableName]
  if (ncol(data) == 0) {
    data <- dplyr::as_tibble(stats::setNames(rep(list(character()), length(columns)), columns))
  } else {
    colnames(data) <- SqlRender::camelCaseToSnakeCase(colnames(data))
  }
  if (!setequal(colnames(data), columns) || anyDuplicated(colnames(data))) {
    stop(sprintf(
      "The columns of %s (%s) differ from its specification (%s)",
      tableName, paste(colnames(data), collapse = ", "), paste(columns, collapse = ", ")
    ))
  }
  data <- data[, columns]
  data[] <- lapply(data, function(column) if (is.logical(column)) as.integer(column) else column)
  fileName <- paste0(tableName, ".csv")
  readr::write_csv(x = data, file = file.path(exportFolder, fileName), na = "")
  return(fileName)
}

# Counts greater than 0 and below minValues are written as -minValues, by
# CohortMethod::enforceMinCellValue(). A table without the column (an empty table) is returned as is.
.enforceMinCellValue <- function(data, fieldName, minValues, silent = FALSE) {
  if (is.null(data[[fieldName]])) {
    return(data)
  }
  return(CohortMethod::enforceMinCellValue(data, fieldName, minValues, silent = silent))
}

.enforceMinCellValues <- function(data, tableName, minCellCount, silent = TRUE) {
  specifications <- getResultsDataModelSpecifications()
  fields <- specifications$columnName[specifications$tableName == tableName & specifications$minCellCount == "Yes"]
  for (field in SqlRender::snakeCaseToCamelCase(fields)) {
    data <- .enforceMinCellValue(data, field, minCellCount, silent = silent)
  }
  return(data)
}

# Definition tables ----------------------------------------------------------------------------------

.exportDefinitions <- function(specifications, targetTreatment, exportFolder) {
  message("- ccw_analysis, ccw_strategy, ccw_target_treatment and ccw_target_treatment_outcome tables")
  analyses <- dplyr::bind_rows(lapply(specifications$ccwAnalysisList, function(analysis) {
    tibble(
      analysisId = as.integer(analysis$analysisId),
      description = analysis$description,
      definition = as.character(analysis$toJson())
    )
  }))
  strategies <- dplyr::bind_rows(lapply(specifications$ccwAnalysisList, function(analysis) {
    dplyr::bind_rows(lapply(analysis$createClonesArgs$strategies, function(strategy) {
      tibble(
        analysisId = as.integer(analysis$analysisId),
        strategyId = as.integer(strategy$strategyId),
        label = strategy$label,
        type = strategy$type,
        treatment = strategy$treatment,
        graceStart = as.integer(strategy$graceStart),
        graceEnd = if (is.null(strategy$graceEnd)) NA_integer_ else as.integer(strategy$graceEnd),
        prohibitedTreatments = if (length(strategy$prohibitedTreatments) == 0) NA_character_ else paste(strategy$prohibitedTreatments, collapse = ","),
        definition = as.character(strategy$toJson())
      )
    }))
  }))
  targetTreatmentOutcome <- dplyr::bind_rows(lapply(specifications$targetTreatmentOutcomesList, function(tto) {
    tibble(
      targetId = as.numeric(tto$targetId),
      treatmentId = as.numeric(tto$treatmentId),
      comparatorTreatmentId = if (is.null(tto$comparatorTreatmentId)) NA_real_ else as.numeric(tto$comparatorTreatmentId),
      outcomeId = vapply(tto$outcomes, function(outcome) as.numeric(outcome$outcomeId), numeric(1)),
      outcomeOfInterest = vapply(tto$outcomes, function(outcome) as.integer(outcome$outcomeOfInterest), integer(1)),
      trueEffectSize = vapply(tto$outcomes, function(outcome) {
        if (is.null(outcome$trueEffectSize)) NA_real_ else as.numeric(outcome$trueEffectSize)
      }, numeric(1))
    )
  })) |>
    dplyr::inner_join(targetTreatment, by = c("targetId", "treatmentId", "comparatorTreatmentId")) |>
    dplyr::select("targetTreatmentId", "outcomeId", "outcomeOfInterest", "trueEffectSize") |>
    dplyr::distinct(.data$targetTreatmentId, .data$outcomeId, .keep_all = TRUE)
  return(c(
    .writeToCsv(analyses, "ccw_analysis", exportFolder),
    .writeToCsv(strategies, "ccw_strategy", exportFolder),
    .writeToCsv(targetTreatment[, c("targetTreatmentId", "targetId", "treatmentId", "comparatorTreatmentId")], "ccw_target_treatment", exportFolder),
    .writeToCsv(targetTreatmentOutcome, "ccw_target_treatment_outcome", exportFolder)
  ))
}

# Covariates -----------------------------------------------------------------------------------------

.exportCovariates <- function(outputFolder, reference, exportFolder, databaseId) {
  message("- ccw_covariate_analysis and ccw_covariate tables")
  covariateAnalysis <- list()
  covariates <- list()
  for (ccwDataFile in unique(reference$ccwDataFile)) {
    data <- loadCloneCensorWeightData(file.path(outputFolder, ccwDataFile))
    analysisIds <- tibble(analysisId = unique(reference$analysisId[reference$ccwDataFile == ccwDataFile]))
    covariateAnalysis[[length(covariateAnalysis) + 1]] <- data$analysisRef |>
      dplyr::select(covariateAnalysisId = "analysisId", covariateAnalysisName = "analysisName") |>
      dplyr::collect() |>
      dplyr::cross_join(analysisIds)
    covariates[[length(covariates) + 1]] <- data$covariateRef |>
      dplyr::select("covariateId", "covariateName", covariateAnalysisId = "analysisId") |>
      dplyr::collect() |>
      dplyr::cross_join(analysisIds)
    close(data)
  }
  covariateAnalysis <- dplyr::bind_rows(covariateAnalysis)
  if (nrow(covariateAnalysis) > 0) {
    covariateAnalysis <- covariateAnalysis |>
      dplyr::distinct(.data$covariateAnalysisId, .data$analysisId, .keep_all = TRUE) |>
      dplyr::mutate(
        covariateAnalysisId = as.integer(.data$covariateAnalysisId),
        analysisId = as.integer(.data$analysisId),
        databaseId = !!databaseId
      )
  }
  covariates <- dplyr::bind_rows(covariates)
  if (nrow(covariates) > 0) {
    covariates <- covariates |>
      dplyr::distinct(.data$covariateId, .data$analysisId, .keep_all = TRUE) |>
      dplyr::mutate(
        covariateAnalysisId = as.integer(.data$covariateAnalysisId),
        analysisId = as.integer(.data$analysisId),
        databaseId = !!databaseId
      )
  }
  return(c(
    .writeToCsv(covariateAnalysis, "ccw_covariate_analysis", exportFolder),
    .writeToCsv(covariates, "ccw_covariate", exportFolder)
  ))
}

# Results and diagnostics ----------------------------------------------------------------------------

.exportResults <- function(outputFolder, targetTreatment, strategyIds, exportFolder, databaseId, minCellCount) {
  message("- ccw_result table")
  results <- getResultsSummary(outputFolder) |>
    dplyr::inner_join(targetTreatment, by = c("targetId", "treatmentId", "comparatorTreatmentId")) |>
    dplyr::inner_join(strategyIds, by = "analysisId") |>
    dplyr::mutate(databaseId = !!databaseId) |>
    dplyr::select(
      "databaseId",
      "analysisId",
      "targetTreatmentId",
      "outcomeId",
      "targetStrategyId",
      "comparatorStrategyId",
      "estimand",
      "estimate",
      "ci95Lb",
      "ci95Ub",
      "p",
      "oneSidedP",
      "logRr",
      "seLogRr",
      "se",
      "calibratedRr",
      "calibratedCi95Lb",
      "calibratedCi95Ub",
      "calibratedP",
      "calibratedOneSidedP",
      "calibratedLogRr",
      "calibratedSeLogRr",
      "ease",
      "targetSubjects",
      "comparatorSubjects",
      "targetClones",
      "comparatorClones",
      "targetDays",
      "comparatorDays",
      "targetOutcomes",
      "comparatorOutcomes",
      "targetCompetingEvents",
      "comparatorCompetingEvents",
      "targetWeightedDays",
      "comparatorWeightedDays",
      "targetWeightedOutcomes",
      "comparatorWeightedOutcomes"
    ) |>
    .enforceMinCellValues("ccw_result", minCellCount, silent = FALSE)
  return(.writeToCsv(results, "ccw_result", exportFolder))
}

.exportDiagnosticsSummary <- function(outputFolder, targetTreatment, strategyIds, exportFolder, databaseId) {
  message("- ccw_diagnostics_summary table")
  diagnostics <- getDiagnosticsSummary(outputFolder) |>
    dplyr::inner_join(targetTreatment, by = c("targetId", "treatmentId", "comparatorTreatmentId")) |>
    dplyr::inner_join(strategyIds, by = "analysisId") |>
    dplyr::mutate(databaseId = !!databaseId) |>
    dplyr::select(-"targetId", -"treatmentId", -"comparatorTreatmentId")
  return(.writeToCsv(diagnostics, "ccw_diagnostics_summary", exportFolder))
}

# Covariate balance ----------------------------------------------------------------------------------

.exportCovariateBalance <- function(outputFolder, reference, specifications, exportFolder, databaseId, minCellCount) {
  message("- ccw_covariate_balance table")
  sdmThreshold <- specifications$ccwDiagnosticThresholds$sdmThreshold
  reference <- reference[reference$outcomeOfInterest & reference$balanceFile != "", ]
  balances <- lapply(unique(reference$balanceFile), function(balanceFile) {
    rows <- reference[reference$balanceFile == balanceFile, ]
    balance <- .tidyBalance(readRDS(file.path(outputFolder, balanceFile)), minCellCount, sdmThreshold)
    rows |>
      dplyr::transmute(
        databaseId = !!databaseId,
        analysisId = as.integer(.data$analysisId),
        .data$targetTreatmentId,
        .data$outcomeId,
        .data$targetStrategyId,
        .data$comparatorStrategyId
      ) |>
      dplyr::cross_join(balance)
  })
  return(.writeToCsv(dplyr::bind_rows(balances), "ccw_covariate_balance", exportFolder))
}

# The balance of computeCovariateBalance(), in the columns of CohortMethod's cm_covariate_balance, with
# the balanced flags at the sdmThreshold of the diagnostics, without covariates absent from both
# strategies, and with the means censored at minCellCount divided by the number of clones (inferred from
# the sums and means, as in CohortMethod).
.tidyBalance <- function(balance, minCellCount, sdmThreshold) {
  inferredSize <- function(sum, mean) {
    size <- mean(sum / mean, na.rm = TRUE)
    if (is.nan(size)) NA_real_ else size
  }
  targetBeforeSize <- inferredSize(balance$beforeMatchingSumTarget, balance$beforeMatchingMeanTarget)
  comparatorBeforeSize <- inferredSize(balance$beforeMatchingSumComparator, balance$beforeMatchingMeanComparator)
  targetAfterSize <- inferredSize(balance$afterMatchingSumTarget, balance$afterMatchingMeanTarget)
  comparatorAfterSize <- inferredSize(balance$afterMatchingSumComparator, balance$afterMatchingMeanComparator)
  zeroIfNa <- function(x) ifelse(is.na(x), 0, x)
  balance <- tibble(
    covariateId = balance$covariateId,
    targetMeanBefore = zeroIfNa(balance$beforeMatchingMeanTarget),
    comparatorMeanBefore = zeroIfNa(balance$beforeMatchingMeanComparator),
    meanBefore = zeroIfNa(balance$beforeMatchingMean),
    stdDiffBefore = balance$beforeMatchingStdDiff,
    stdDiffVarBefore = balance$beforeMatchingSdmVariance,
    balancedBefore = as.integer(abs(balance$beforeMatchingStdDiff) < sdmThreshold),
    targetMeanAfter = zeroIfNa(balance$afterMatchingMeanTarget),
    comparatorMeanAfter = zeroIfNa(balance$afterMatchingMeanComparator),
    meanAfter = zeroIfNa(balance$afterMatchingMean),
    stdDiffAfter = balance$afterMatchingStdDiff,
    stdDiffVarAfter = balance$afterMatchingSdmVariance,
    balancedAfter = as.integer(abs(balance$afterMatchingStdDiff) < sdmThreshold),
    targetStdDiff = balance$targetStdDiff,
    comparatorStdDiff = balance$comparatorStdDiff,
    targetComparatorStdDiff = balance$targetComparatorStdDiff
  )
  balance <- balance[!(round(balance$targetMeanBefore, 3) == 0 & round(balance$comparatorMeanBefore, 3) == 0 &
    round(balance$targetMeanAfter, 3) == 0 & round(balance$comparatorMeanAfter, 3) == 0), ]
  balance <- balance |>
    .enforceMinCellValue("targetMeanBefore", minCellCount / targetBeforeSize, silent = TRUE) |>
    .enforceMinCellValue("comparatorMeanBefore", minCellCount / comparatorBeforeSize, silent = TRUE) |>
    .enforceMinCellValue("meanBefore", minCellCount / (targetBeforeSize + comparatorBeforeSize), silent = TRUE) |>
    .enforceMinCellValue("targetMeanAfter", minCellCount / targetAfterSize, silent = TRUE) |>
    .enforceMinCellValue("comparatorMeanAfter", minCellCount / comparatorAfterSize, silent = TRUE) |>
    .enforceMinCellValue("meanAfter", minCellCount / (targetAfterSize + comparatorAfterSize), silent = TRUE)
  rounded <- c(
    "targetMeanBefore", "comparatorMeanBefore", "meanBefore", "stdDiffBefore", "targetMeanAfter", "comparatorMeanAfter",
    "meanAfter", "stdDiffAfter", "targetStdDiff", "comparatorStdDiff", "targetComparatorStdDiff"
  )
  balance[rounded] <- lapply(balance[rounded], round, digits = 3)
  return(balance)
}

# Censoring models -----------------------------------------------------------------------------------

.exportCensoringModels <- function(outputFolder, reference, exportFolder, databaseId) {
  message("- ccw_censoring_model table")
  reference <- reference[reference$outcomeOfInterest, ]
  models <- lapply(unique(reference$censoringModelFile), function(censoringModelFile) {
    coefficients <- dplyr::bind_rows(lapply(readRDS(file.path(outputFolder, censoringModelFile)), function(model) {
      dplyr::bind_rows(lapply(c(list(model), if (is.null(model$postGrace)) list() else list(model$postGrace)), function(fit) {
        coefficients <- fit$coefficients[!is.na(fit$coefficients) & fit$coefficients != 0]
        tibble(
          position = rep(as.integer(model$strategyId), length(coefficients)),
          model = rep(fit$kind, length(coefficients)),
          term = names(coefficients),
          covariateId = suppressWarnings(as.numeric(sub(":.*$", "", names(coefficients)))),
          coefficient = as.numeric(coefficients)
        )
      }))
    }))
    rows <- reference[reference$censoringModelFile == censoringModelFile, ]
    dplyr::bind_rows(lapply(seq_len(nrow(rows)), function(i) {
      coefficients |>
        dplyr::mutate(
          databaseId = !!databaseId,
          analysisId = as.integer(rows$analysisId[i]),
          targetTreatmentId = rows$targetTreatmentId[i],
          outcomeId = rows$outcomeId[i],
          strategyId = .strategyIdOf(rows[i, ], .data$position)
        ) |>
        dplyr::select(-"position")
    }))
  })
  return(.writeToCsv(dplyr::bind_rows(models), "ccw_censoring_model", exportFolder))
}

# Tables from the outcome model files ----------------------------------------------------------------

# ccw_attrition, ccw_follow_up_dist, ccw_censoring_model_dist, ccw_weight_dist and
# ccw_cumulative_incidence, from the outcome models of the outcomes of interest.
.exportFromOutcomeModels <- function(outputFolder, reference, exportFolder, databaseId, minCellCount) {
  message("- ccw_attrition, ccw_follow_up_dist, ccw_censoring_model_dist, ccw_weight_dist and ccw_cumulative_incidence tables")
  reference <- reference[reference$outcomeOfInterest, ]
  tables <- list(attrition = list(), followUp = list(), censoringDist = list(), weightDist = list(), incidence = list())
  for (i in seq_len(nrow(reference))) {
    refRow <- reference[i, ]
    outcomeModel <- readRDS(file.path(outputFolder, refRow$outcomeModelFile))
    keys <- function(data, position) {
      data <- dplyr::mutate(
        data,
        databaseId = !!databaseId,
        analysisId = as.integer(refRow$analysisId),
        targetTreatmentId = refRow$targetTreatmentId,
        outcomeId = refRow$outcomeId,
        strategyId = .strategyIdOf(refRow, position)
      )
      return(data)
    }
    if (!is.null(outcomeModel$attrition)) {
      attrition <- outcomeModel$attrition
      for (position in 1:2) {
        tables$attrition[[length(tables$attrition) + 1]] <- tibble(
          sequenceNumber = seq_len(nrow(attrition)),
          description = attrition$description,
          subjects = as.integer(if (position == 1) attrition$targetPersons else attrition$comparatorPersons)
        ) |>
          keys(position)
      }
    }
    if (!is.null(outcomeModel$followUpDistribution)) {
      followUp <- outcomeModel$followUpDistribution
      tables$followUp[[length(tables$followUp) + 1]] <- keys(dplyr::select(followUp, -"strategyId"), followUp$strategyId)
    }
    positivity <- outcomeModel$diagnostics$positivity
    densities <- if (is.null(positivity)) NULL else .censoringModelDensities(positivity)
    if (!is.null(densities) && nrow(densities) > 0) {
      tables$censoringDist[[length(tables$censoringDist) + 1]] <- keys(dplyr::select(densities, -"strategyId"), densities$strategyId)
    }
    weightSummary <- outcomeModel$diagnostics$weightSummary
    if (!is.null(weightSummary)) {
      weightDist <- tibble(
        scale = weightSummary$scale,
        truncation = weightSummary$truncation,
        bandStart = as.integer(weightSummary$bandStart),
        dayBand = weightSummary$dayBand,
        cloneDays = weightSummary$cloneDays,
        days = as.integer(weightSummary$days),
        meanWeight = weightSummary$meanWeight,
        minWeight = weightSummary$minWeight,
        p25Weight = weightSummary$p25Weight,
        medianWeight = weightSummary$medianWeight,
        p75Weight = weightSummary$p75Weight,
        p95Weight = weightSummary$p95Weight,
        p99Weight = weightSummary$p99Weight,
        maxWeight = weightSummary$maxWeight,
        effectiveSampleSize = weightSummary$effectiveSampleSize,
        effectiveSampleShare = weightSummary$effectiveSampleShare,
        shareTruncated = weightSummary$shareTruncated,
        cap = weightSummary$cap
      )
      tables$weightDist[[length(tables$weightDist) + 1]] <- keys(weightDist, weightSummary$strategyId)
    }
    curves <- outcomeModel$curves
    if (!is.null(curves)) {
      incidence <- tibble(
        timeDay = as.integer(curves$day),
        risk = curves$risk,
        riskLb = curves$lb95,
        riskUb = curves$ub95,
        atRiskWeighted = curves$atRiskWeighted,
        atRisk = as.integer(round(curves$atRiskCount))
      )
      tables$incidence[[length(tables$incidence) + 1]] <- keys(incidence, curves$strategyId)
    }
  }
  bind <- function(table, tableName) {
    .enforceMinCellValues(dplyr::bind_rows(table), tableName, minCellCount)
  }
  return(c(
    .writeToCsv(bind(tables$attrition, "ccw_attrition"), "ccw_attrition", exportFolder),
    .writeToCsv(bind(tables$followUp, "ccw_follow_up_dist"), "ccw_follow_up_dist", exportFolder),
    .writeToCsv(bind(tables$censoringDist, "ccw_censoring_model_dist"), "ccw_censoring_model_dist", exportFolder),
    .writeToCsv(bind(tables$weightDist, "ccw_weight_dist"), "ccw_weight_dist", exportFolder),
    .writeToCsv(bind(tables$incidence, "ccw_cumulative_incidence"), "ccw_cumulative_incidence", exportFolder)
  ))
}

# ccw_likelihood_profile, as CohortMethod's cm_likelihood_profile (CohortMethod R/Export.R,
# exportLikelihoodProfiles(), Apache 2.0): the likelihood profile of the log hr, shr or or of every outcome
# model (controls included, and fits that did not converge), with the log likelihood normalized to a
# maximum of 0 and the gradient as is, without a minimum cell count. Unlike CohortMethod, a contrast that
# is null by design (standard error 0 up to rounding) is left out, and so are non-finite points.
.exportLikelihoodProfiles <- function(outputFolder, reference, exportFolder, databaseId) {
  message("- ccw_likelihood_profile table")
  reference <- reference[reference$outcomeModelFile != "", ]
  profiles <- list()
  for (i in seq_len(nrow(reference))) {
    refRow <- reference[i, ]
    outcomeModel <- readRDS(file.path(outputFolder, refRow$outcomeModelFile))
    profile <- outcomeModel$logLikelihoodProfile
    if (is.null(profile)) {
      next
    }
    estimand <- .profiledEstimands[[outcomeModel$modelType]]
    # A contrast that is null by design has a peaked pseudo-likelihood but no information. Its standard
    # error is 0, or rounding noise (of the order of 1e-7 with Cyclops' robust variance); a log ratio with a
    # standard error below 1e-6 would need of the order of 1e12 events:
    if (isTRUE(outcomeModel$estimates$seLogEstimate[outcomeModel$estimates$estimand == estimand] < 1e-6)) {
      next
    }
    profile <- profile[is.finite(profile$point) & is.finite(profile$value), ]
    if (nrow(profile) == 0) {
      next
    }
    profiles[[length(profiles) + 1]] <- tibble(
      databaseId = databaseId,
      analysisId = as.integer(refRow$analysisId),
      targetTreatmentId = refRow$targetTreatmentId,
      outcomeId = refRow$outcomeId,
      targetStrategyId = refRow$targetStrategyId,
      comparatorStrategyId = refRow$comparatorStrategyId,
      estimand = estimand,
      logRr = profile$point,
      logLikelihood = profile$value - max(profile$value),
      gradient = profile$derivative
    )
  }
  return(.writeToCsv(dplyr::bind_rows(profiles), "ccw_likelihood_profile", exportFolder))
}

# The density of the predicted probability of starting the treatment of each censoring model, on 100
# points from 0 to 1, as CohortMethod's preference score distribution. Models whose predictions do not
# vary are left out.
.censoringModelDensities <- function(positivity) {
  predictions <- attr(positivity, "predictions")
  densities <- lapply(seq_len(nrow(positivity)), function(i) {
    probability <- predictions$probability[predictions$strategyId == positivity$strategyId[i] & predictions$model == positivity$model[i] &
      predictions$day == positivity$day[i]]
    probability <- probability[!is.na(probability)]
    if (length(probability) < 2 || min(probability) == max(probability)) {
      return(NULL)
    }
    density <- stats::density(probability, from = 0, to = 1, n = 100)
    tibble(
      strategyId = positivity$strategyId[i],
      model = positivity$model[i],
      horizon = as.integer(positivity$day[i]),
      probability = density$x,
      density = density$y
    )
  })
  return(dplyr::bind_rows(densities))
}

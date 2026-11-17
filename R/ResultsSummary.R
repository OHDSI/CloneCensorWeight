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

# The summary of runCcwAnalyses(), modelled on CohortMethod's summarizeResults(), addBalance() and
# calibrateEstimates() (R/RunAnalyses.R, Copyright Observational Health Data Sciences and
# Informatics, Apache 2.0). See getResultsSummary() and getDiagnosticsSummary() for the columns.

# Estimands that are ratios, which are calibrated (plan section 9 item 8 (ii)).
.ratioEstimandPattern <- "^(hr|or|rr_[0-9]+)$"

.diagnostic <- function(value, pass) {
  return(dplyr::case_when(
    is.na(value) ~ "NOT EVALUATED",
    pass ~ "PASS",
    TRUE ~ "FAIL"
  ))
}

summarizeResults <- function(referenceTable,
                             outputFolder,
                             mainFileName,
                             diagnosticsSummaryFileName,
                             calibrationThreads,
                             ccwDiagnosticThresholds) {
  keyColumns <- c("analysisId", "targetId", "treatmentId", "comparatorTreatmentId", "outcomeId")
  thresholds <- ccwDiagnosticThresholds
  subset <- referenceTable[referenceTable$outcomeModelFile != "", ]

  balanceFiles <- unique(subset$balanceFile[subset$balanceFile != ""])
  maxSdms <- vapply(balanceFiles, function(balanceFile) {
    balance <- readRDS(file.path(outputFolder, balanceFile))
    sdm <- abs(balance$afterMatchingStdDiff)
    sdm <- sdm[!is.na(sdm)]
    if (length(sdm) == 0) NA_real_ else max(sdm)
  }, numeric(1))

  results <- vector("list", nrow(subset))
  diagnostics <- vector("list", nrow(subset))
  for (i in seq_len(nrow(subset))) {
    refRow <- subset[i, ]
    outcomeModel <- readRDS(file.path(outputFolder, refRow$outcomeModelFile))
    counts <- outcomeModel$counts
    count <- function(strategyId, column) {
      value <- counts[[column]][counts$strategyId == strategyId]
      if (length(value) == 0) 0 else value
    }
    estimates <- outcomeModel$estimates
    isRatio <- grepl(.ratioEstimandPattern, estimates$estimand)
    oneSidedP <- rep(NA_real_, nrow(estimates))
    ok <- isRatio & !is.na(estimates$logEstimate) & !is.na(estimates$seLogEstimate) & estimates$seLogEstimate > 0
    if (any(ok)) {
      oneSidedP[ok] <- EmpiricalCalibration::computeTraditionalP(
        logRr = estimates$logEstimate[ok],
        seLogRr = estimates$seLogEstimate[ok],
        twoSided = FALSE,
        upper = TRUE
      )
    }
    results[[i]] <- tibble(
      analysisId = refRow$analysisId,
      targetId = refRow$targetId,
      treatmentId = refRow$treatmentId,
      comparatorTreatmentId = refRow$comparatorTreatmentId,
      outcomeId = refRow$outcomeId,
      trueEffectSize = refRow$trueEffectSize,
      estimand = estimates$estimand,
      targetSubjects = count(1, "persons"),
      comparatorSubjects = count(2, "persons"),
      targetClones = count(1, "clones"),
      comparatorClones = count(2, "clones"),
      targetDays = count(1, "personDays"),
      comparatorDays = count(2, "personDays"),
      targetOutcomes = count(1, "outcomes"),
      comparatorOutcomes = count(2, "outcomes"),
      targetCompetingEvents = count(1, "competingEvents"),
      comparatorCompetingEvents = count(2, "competingEvents"),
      targetWeightedDays = count(1, "weightedPersonDays"),
      comparatorWeightedDays = count(2, "weightedPersonDays"),
      targetWeightedOutcomes = count(1, "weightedOutcomes"),
      comparatorWeightedOutcomes = count(2, "weightedOutcomes"),
      estimate = estimates$estimate,
      ci95Lb = estimates$lb95,
      ci95Ub = estimates$ub95,
      p = estimates$p,
      oneSidedP = oneSidedP,
      logRr = ifelse(isRatio, estimates$logEstimate, NA_real_),
      seLogRr = ifelse(isRatio, estimates$seLogEstimate, NA_real_),
      se = ifelse(isRatio, NA_real_, estimates$se)
    )
    values <- outcomeModel$diagnostics$values
    maxSdm <- if (refRow$balanceFile == "") NA_real_ else maxSdms[[refRow$balanceFile]]
    diagnostics[[i]] <- tibble(
      analysisId = refRow$analysisId,
      targetId = refRow$targetId,
      treatmentId = refRow$treatmentId,
      comparatorTreatmentId = refRow$comparatorTreatmentId,
      outcomeId = refRow$outcomeId,
      maxSdm = maxSdm,
      balanceDiagnostic = .diagnostic(maxSdm, maxSdm < thresholds$sdmThreshold),
      mdrr = values$mdrr,
      mdrrDiagnostic = .diagnostic(values$mdrr, values$mdrr < thresholds$mdrrThreshold),
      essFraction = values$essFraction,
      essFractionDiagnostic = .diagnostic(values$essFraction, values$essFraction >= thresholds$essFractionThreshold),
      maxWeight = values$maxWeight,
      maxWeightBeforeTruncation = values$maxWeightBeforeTruncation,
      maxWeightDiagnostic = .diagnostic(values$maxWeight, values$maxWeight < thresholds$maxWeightThreshold),
      positivity = values$positivity,
      positivityShareExtreme = values$positivityShareExtreme,
      positivityDiagnostic = .diagnostic(values$positivity, values$positivity >= thresholds$positivityThreshold)
    )
  }
  results <- dplyr::bind_rows(results)
  diagnostics <- dplyr::bind_rows(diagnostics)
  diagnostics$unblindForCalibration <- diagnostics$balanceDiagnostic != "FAIL" &
    diagnostics$essFractionDiagnostic != "FAIL" &
    diagnostics$maxWeightDiagnostic != "FAIL" &
    diagnostics$positivityDiagnostic != "FAIL"

  results <- calibrateEstimates(
    results = dplyr::left_join(results, diagnostics[, c(keyColumns, "unblindForCalibration")], by = keyColumns),
    calibrationThreads = calibrationThreads
  )
  withEase <- results[!is.na(results$ease), ]
  if (nrow(withEase) == 0) {
    diagnostics$ease <- rep(NA_real_, nrow(diagnostics))
  } else {
    ease <- withEase |>
      dplyr::group_by(dplyr::across(dplyr::all_of(keyColumns))) |>
      dplyr::summarise(ease = max(abs(.data$ease)), .groups = "drop")
    diagnostics <- dplyr::left_join(diagnostics, ease, by = keyColumns)
  }
  diagnostics$easeDiagnostic <- .diagnostic(diagnostics$ease, abs(diagnostics$ease) < thresholds$easeThreshold)
  diagnostics$unblindForEvidenceSynthesis <- diagnostics$unblindForCalibration & diagnostics$easeDiagnostic != "FAIL"
  diagnostics$unblind <- diagnostics$unblindForEvidenceSynthesis & diagnostics$mdrrDiagnostic != "FAIL"

  resultsSummary <- results[, c(
    keyColumns,
    "estimand",
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
    "comparatorWeightedOutcomes",
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
    "calibratedSeLogRr"
  )]
  resultsSummary <- resultsSummary[order(
    resultsSummary$analysisId, resultsSummary$targetId, resultsSummary$treatmentId,
    resultsSummary$comparatorTreatmentId, resultsSummary$outcomeId
  ), ]
  diagnosticsSummary <- diagnostics[, c(
    keyColumns,
    "maxSdm",
    "balanceDiagnostic",
    "mdrr",
    "mdrrDiagnostic",
    "ease",
    "easeDiagnostic",
    "essFraction",
    "essFractionDiagnostic",
    "maxWeight",
    "maxWeightBeforeTruncation",
    "maxWeightDiagnostic",
    "positivity",
    "positivityShareExtreme",
    "positivityDiagnostic",
    "unblindForEvidenceSynthesis",
    "unblind"
  )]
  .saveFile(diagnosticsSummary, diagnosticsSummaryFileName)
  .saveFile(resultsSummary, mainFileName)
  invisible(NULL)
}

# Calibrates the ratio estimands per analysis, target-treatment combination and estimand, using the
# negative controls whose diagnostics allow calibration (unblindForCalibration).
calibrateEstimates <- function(results, calibrationThreads) {
  results$rowIndex <- seq_len(nrow(results))
  isRatio <- grepl(.ratioEstimandPattern, results$estimand)
  groups <- split(
    results[isRatio, ],
    paste(results$analysisId, results$targetId, results$treatmentId, results$comparatorTreatmentId, results$estimand)[isRatio]
  )
  calibrated <- list()
  if (length(groups) > 0) {
    message("Calibrating estimates")
    cluster <- ParallelLogger::makeCluster(min(length(groups), calibrationThreads))
    on.exit(ParallelLogger::stopCluster(cluster))
    ParallelLogger::clusterRequire(cluster, "CloneCensorWeight")
    calibrated <- ParallelLogger::clusterApply(cluster, unname(groups), calibrateGroup)
  }
  calibrated <- dplyr::bind_rows(c(calibrated, list(.addEmptyCalibration(results[!isRatio, ]))))
  calibrated <- calibrated[order(calibrated$rowIndex), ]
  calibrated$rowIndex <- NULL
  return(calibrated)
}

.addEmptyCalibration <- function(group) {
  group$calibratedRr <- rep(NA_real_, nrow(group))
  group$calibratedCi95Lb <- rep(NA_real_, nrow(group))
  group$calibratedCi95Ub <- rep(NA_real_, nrow(group))
  group$calibratedP <- rep(NA_real_, nrow(group))
  group$calibratedOneSidedP <- rep(NA_real_, nrow(group))
  group$calibratedLogRr <- rep(NA_real_, nrow(group))
  group$calibratedSeLogRr <- rep(NA_real_, nrow(group))
  group$ease <- rep(NA_real_, nrow(group))
  return(group)
}

calibrateGroup <- function(group) {
  usable <- !is.na(group$trueEffectSize) & !is.na(group$seLogRr) & group$seLogRr > 0 & group$unblindForCalibration
  ncs <- group[usable & group$trueEffectSize == 1, ]
  pcs <- group[usable & group$trueEffectSize != 1, ]
  if (nrow(ncs) < 5) {
    return(.addEmptyCalibration(group))
  }
  null <- EmpiricalCalibration::fitMcmcNull(logRr = ncs$logRr, seLogRr = ncs$seLogRr)
  ease <- EmpiricalCalibration::computeExpectedAbsoluteSystematicError(null)
  calibratedP <- EmpiricalCalibration::calibrateP(null = null, logRr = group$logRr, seLogRr = group$seLogRr)
  calibratedOneSidedP <- EmpiricalCalibration::calibrateP(
    null = null,
    logRr = group$logRr,
    seLogRr = group$seLogRr,
    twoSided = FALSE,
    upper = TRUE
  )
  # Synthetic positive controls have their true effect size on the hazard scale, so they calibrate
  # only the hazard ratio (plan section 3.8 (e)); negative controls are null on every scale.
  if (nrow(pcs) >= 5 && group$estimand[1] == "hr") {
    model <- EmpiricalCalibration::fitSystematicErrorModel(
      logRr = c(ncs$logRr, pcs$logRr),
      seLogRr = c(ncs$seLogRr, pcs$seLogRr),
      trueLogRr = log(c(ncs$trueEffectSize, pcs$trueEffectSize)),
      estimateCovarianceMatrix = FALSE
    )
  } else {
    model <- EmpiricalCalibration::convertNullToErrorModel(null)
  }
  calibratedCi <- EmpiricalCalibration::calibrateConfidenceInterval(model = model, logRr = group$logRr, seLogRr = group$seLogRr)
  group$calibratedRr <- exp(calibratedCi$logRr)
  group$calibratedCi95Lb <- exp(calibratedCi$logLb95Rr)
  group$calibratedCi95Ub <- exp(calibratedCi$logUb95Rr)
  group$calibratedP <- calibratedP$p
  group$calibratedOneSidedP <- calibratedOneSidedP$p
  group$calibratedLogRr <- calibratedCi$logRr
  group$calibratedSeLogRr <- calibratedCi$seLogRr
  group$ease <- ease$ease
  return(group)
}

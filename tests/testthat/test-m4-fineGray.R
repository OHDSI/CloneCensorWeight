# Tier 1 reference tests of the weighted Fine-Gray MSM (modelType = "fineGray"; plan sections 2.6 and
# 6.4). Fixtures run with riskStart = 0 and 1.

library(testthat)

fineGrayCoxControl <- function() survival::coxph.control(eps = 1e-12, toler.chol = 1e-14, iter.max = 100)

# The weighted Fine-Gray rows built independently of the package: each clone that dies at T stays in
# the risk set of every later outcome time s (of either strategy) with weight w(T) G(s-) / G(T-), on the
# grid of the outcome times (one row per outcome time), where G is a reverse Kaplan-Meier (survfit) of
# the clone's strategy, with outcomes and deaths moved just before the censorings at their time:
# - "natural": of the ends at observationEnd and riskEnd, weighted by the segment weights;
# - "finegray": of every end that is not an outcome or a death, unweighted.
handBuiltFineGrayRows <- function(fixture, censoring) {
  data <- toSurvivalData(fixture$segments)
  data$endReason <- fixture$clones$endReason[match(data$cloneRowId, fixture$clones$cloneRowId)]
  data$last <- !duplicated(data$cloneRowId, fromLast = TRUE)
  data$event <- data$outcome == 1 | data$competing == 1
  outcomeTimes <- sort(unique(data$stopTime[data$outcome == 1]))
  rows <- data
  rows$status <- rows$outcome
  tails <- list()
  for (strategyId in 1:2) {
    strategy <- data[data$strategyId == strategyId, ]
    if (censoring == "natural") {
      censored <- strategy$last & !strategy$event & strategy$endReason %in% c("observationEnd", "riskEnd")
      g <- strategy$weight
    } else {
      censored <- strategy$last & !strategy$event
      g <- rep(1, nrow(strategy))
    }
    fit <- survival::survfit(
      survival::Surv(startTime, stopTime - 0.2 * event, as.integer(censored)) ~ 1,
      data = strategy, weights = g
    )
    before <- function(t) {
      jumps <- fit$time[fit$n.event > 0]
      surv <- fit$surv[fit$n.event > 0]
      vapply(t, function(u) if (any(jumps < u)) surv[max(which(jumps < u))] else 1, numeric(1))
    }
    dead <- strategy[strategy$competing == 1, ]
    for (i in seq_len(nrow(dead))) {
      later <- outcomeTimes[outcomeTimes > dead$stopTime[i]]
      if (length(later) == 0) {
        next
      }
      tail <- dead[rep(i, length(later)), ]
      tail$startTime <- c(dead$stopTime[i], later[-length(later)])
      tail$stopTime <- later
      tail$weight <- dead$weight[i] * before(later) / before(dead$stopTime[i])
      tail$status <- 0L
      tails[[length(tails) + 1]] <- tail
    }
  }
  return(do.call(rbind, c(list(rows), tails)))
}

# The weight of each clone in the risk set of each outcome time (t, cloneRowId, weight) of rows
# (startTime, stopTime, weight, cloneRowId).
riskSetWeights <- function(rows, times) {
  cells <- do.call(rbind, lapply(times, function(t) {
    atRisk <- rows$startTime < t & rows$stopTime >= t
    data.frame(t = rep(t, sum(atRisk)), cloneRowId = rows$cloneRowId[atRisk], weight = rows$weight[atRisk])
  }))
  return(stats::aggregate(weight ~ t + cloneRowId, cells, sum))
}

test_that("The \"finegray\" censoring distribution gives the weights of survival::finegray", {
  skip_if_not_installed("survival")
  fixture <- createM4Fixture(riskStart = 0, competingEvents = "death", seed = 3)
  data <- toSurvivalData(fixture$segments)
  weightedClones <- createM4WeightedClones(fixture)
  outcomeData <- .getOutcomeData(weightedClones)
  close(weightedClones)
  rows <- .fineGraySegments(outcomeData, "finegray")
  times <- sort(unique(data$stopTime[data$outcome == 1]))
  for (strategyId in 1:2) {
    # finegray() is called per strategy: it extends a clone that dies over the outcome times of its own
    # stratum only (as for a Cox model stratified by strategy), and with strata() and weights together
    # survival 3.8-6 misaligns the weights (user.weights[split$row], with split$row indexing the
    # stratum's rows). So it is compared on its own (time, clone) cells:
    strategy <- data[data$strategyId == strategyId, ]
    fgData <- survival::finegray(
      survival::Surv(startTime, stopTime, state) ~ x + cloneRowId,
      data = strategy, weights = weight, id = cloneRowId, etype = "outcome"
    )
    reference <- riskSetWeights(data.frame(startTime = fgData$fgstart, stopTime = fgData$fgstop, weight = fgData$fgwt, cloneRowId = fgData$cloneRowId), times)
    ours <- riskSetWeights(rows[rows$strategyId == strategyId, ], times)
    matched <- merge(reference, ours, by = c("t", "cloneRowId"), all.x = TRUE, suffixes = c(".reference", ".ours"))
    expect_gt(nrow(matched), 100)
    expect_equal(matched$weight.ours, matched$weight.reference, tolerance = 1e-10)
  }
})

test_that("Fine-Gray MSM matches coxph on hand-built Fine-Gray rows, for both censoring distributions", {
  skip_if_not_installed("survival")
  for (riskStart in 0:1) {
    fixture <- createM4Fixture(riskStart = riskStart, competingEvents = "death", seed = 3 + riskStart)
    estimates <- list()
    for (censoring in c("natural", "finegray")) {
      rows <- handBuiltFineGrayRows(fixture, censoring)
      reference <- survival::coxph(
        survival::Surv(startTime, stopTime, status) ~ x,
        data = rows, weights = weight, cluster = personSeqId, ties = "breslow", control = fineGrayCoxControl()
      )
      for (engine in c("cyclops", "ccw")) {
        robust <- fitM4(fixture, modelType = "fineGray", fineGrayCensoring = censoring, horizons = 10, robustVarianceEngine = engine)
        expect_equal(getM4Estimate(robust, "shr", "logEstimate"), unname(coef(reference)), tolerance = 1e-6)
        expect_equal(getM4Estimate(robust, "shr", "seLogEstimate"), sqrt(vcov(reference)[1, 1]), tolerance = 1e-5)
      }
      modelBased <- fitM4(fixture, modelType = "fineGray", fineGrayCensoring = censoring, horizons = 10, varianceEstimator = "modelBased")
      expect_equal(getM4Estimate(modelBased, "shr", "seLogEstimate"), sqrt(reference$naive.var[1, 1]), tolerance = 1e-6)
      expect_equal(names(robust$coefficients), "target")
      estimates[[censoring]] <- getM4Estimate(robust, "shr", "logEstimate")
    }
    # The two censoring distributions differ (finegray's also counts the deviations):
    expect_gt(abs(estimates$natural - estimates$finegray), 1e-6)
  }
})

test_that("Fine-Gray MSM reports the Aalen-Johansen rows of \"km\" at the horizons", {
  fixture <- createM4Fixture(riskStart = 0, competingEvents = "death", seed = 5)
  fineGray <- fitM4(fixture, modelType = "fineGray", horizons = c(5, 15))
  km <- fitM4(fixture, modelType = "km", horizons = c(5, 15))
  expect_equal(fineGray$estimates$estimand[1], "shr")
  expect_equal(fineGray$estimates[-1, ], km$estimates, ignore_attr = TRUE)
  expect_false(is.na(getM4Estimate(fineGray, "shr")))
})

test_that("The subdistribution hazard ratio is exported as a ratio", {
  fixture <- createM4Fixture(riskStart = 0, competingEvents = "death", seed = 5)
  outcomeModel <- fitM4(fixture, modelType = "fineGray", horizons = 15)
  outcomeModel$diagnostics <- list(values = dplyr::tibble(
    mdrr = 2, positivity = 0.2, positivityShareExtreme = 0, essFraction = 0.5, maxWeight = 4, maxWeightBeforeTruncation = 4
  ))
  outputFolder <- tempfile()
  on.exit(unlink(outputFolder, recursive = TRUE))
  dir.create(outputFolder)
  saveRDS(outcomeModel, file.path(outputFolder, "om.rds"))
  referenceTable <- dplyr::tibble(
    analysisId = 1L, targetId = 4, treatmentId = 1, comparatorTreatmentId = 2, outcomeId = 3,
    trueEffectSize = NA_real_, balanceFile = "", outcomeModelFile = "om.rds"
  )
  suppressMessages(summarizeResults(
    referenceTable = referenceTable,
    outputFolder = outputFolder,
    mainFileName = file.path(outputFolder, "resultsSummary.rds"),
    diagnosticsSummaryFileName = file.path(outputFolder, "diagnosticsSummary.rds"),
    calibrationThreads = 1,
    ccwDiagnosticThresholds = createCcwDiagnosticThresholds()
  ))
  targetTreatment <- dplyr::tibble(targetId = 4, treatmentId = 1, comparatorTreatmentId = 2, targetTreatmentId = 1L)
  strategyIds <- dplyr::tibble(analysisId = 1L, targetStrategyId = 1L, comparatorStrategyId = 2L)
  suppressMessages(.exportResults(outputFolder, targetTreatment, strategyIds, outputFolder, "db", minCellCount = 0))
  results <- readr::read_csv(file.path(outputFolder, "ccw_result.csv"), col_types = readr::cols())
  shr <- results[results$estimand == "shr", ]
  expect_equal(shr$log_rr, getM4Estimate(outcomeModel, "shr", "logEstimate"))
  expect_equal(shr$se_log_rr, getM4Estimate(outcomeModel, "shr", "seLogEstimate"))
  expect_false(is.na(shr$se_log_rr) || is.na(shr$one_sided_p))
  expect_true(is.na(shr$se))

  # The bootstrap summarizes the replicates of the ratio on the log scale:
  t0 <- dplyr::tibble(estimand = "shr", estimate = 1.25, logEstimate = log(1.25))
  replicates <- cbind(shr = c(1.1, 1.3, 1.2, 1.5))
  estimates <- .bootstrapEstimates(t0, replicates, "percentile")
  expect_equal(estimates$seLogEstimate, sd(log(replicates[, "shr"])))
  expect_true(is.na(estimates$se))
})

test_that("Fine-Gray MSM needs competing death and a valid censoring distribution", {
  fixture <- createM4Fixture(riskStart = 0, competingEvents = "none")
  expect_error(fitM4(fixture, modelType = "fineGray"), "competingEvents = \"death\"")
  expect_error(createFitOutcomeModelArgs(modelType = "fineGray", fineGrayCensoring = "all"), "fineGrayCensoring")
  args <- createFitOutcomeModelArgs(modelType = "fineGray", fineGrayCensoring = "finegray")
  restored <- FitOutcomeModelArgs$new(json = args$toJson())
  expect_equal(restored$fineGrayCensoring, "finegray")
  expect_equal(restored$modelType, "fineGray")
})

test_that("The Fine-Gray likelihood profile is the same on native and expanded rows, and is coxph's on hand-built rows", {
  skip_if_not_installed("survival")
  fixture <- createM4Fixture(n = 150, riskStart = 0, competingEvents = "death", seed = 3)
  grid <- c(-0.3, 0, 0.3)
  native <- fitM4(fixture, modelType = "fineGray", horizons = 10, profileGrid = grid, profileBounds = NULL)
  # The ccw robust engine fits the expanded rows:
  expanded <- fitM4(fixture, modelType = "fineGray", horizons = 10, profileGrid = grid, profileBounds = NULL, robustVarianceEngine = "ccw")
  rows <- handBuiltFineGrayRows(fixture, "natural")
  reference <- vapply(grid, function(point) {
    fit <- survival::coxph(
      survival::Surv(startTime, stopTime, status) ~ x,
      data = rows, weights = weight, ties = "breslow", init = point, control = survival::coxph.control(iter.max = 0)
    )
    return(tail(fit$loglik, 1))
  }, numeric(1))
  expect_equal(native$logLikelihoodProfile$value, reference, tolerance = 1e-8)
  expect_equal(expanded$logLikelihoodProfile$value, reference, tolerance = 1e-8)
  expect_equal(native$logLikelihoodProfile$derivative, expanded$logLikelihoodProfile$derivative, tolerance = 1e-8)
})

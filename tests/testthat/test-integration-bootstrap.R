library(testthat)

# Tests of the bootstrap driver, bootstrapCcwAnalysis() (Phase1Contracts section 9, amendment 9; plan
# sections 4.4 and 6.4): reproducibility, replay of given counts against a copy-based run of the
# analysis, confidence limits against boot::boot.ci, all model types, fixed hyperparameters, and
# failed replicates.

bootstrapData <- function(n = 120, seed = 1, initiationCoefficient = 0.4) {
  return(simulateCloneCensorWeightData(
    n = n, maxDays = 15, graceEnd = 3, initiationHazard = 0.1, outcomeHazard = 0.02, deathHazard = 0.005,
    initiationCoefficients = rep(initiationCoefficient, 5), seed = seed
  ))
}

# The analysis settings: initiateWithin(3) vs notInitiateWithin(3), without competing events.
bootstrapSettings <- function(modelType = "km",
                              useCovariates = FALSE,
                              prior = Cyclops::createPrior("none"),
                              replicates = 3,
                              seed = 1,
                              ci = "percentile",
                              fixHyperparameters = TRUE) {
  return(list(
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 15, endAnchor = "cohort start"),
    createClonesArgs = createCreateClonesArgs(
      strategies = list(
        createTreatmentStrategy(1, "Initiate within 3 days", "initiateWithin", graceEnd = 3),
        createTreatmentStrategy(2, "Do not initiate within 3 days", "notInitiateWithin", graceEnd = 3)
      ),
      competingEvents = "none"
    ),
    fitCensoringModelsArgs = createFitCensoringModelsArgs(
      useCovariates = useCovariates,
      prior = prior,
      control = Cyclops::createControl(noiseLevel = "silent", cvType = "auto", seed = 1, resetCoefficients = TRUE, tolerance = 2e-07, cvRepetitions = 1, startingVariance = 0.01)
    ),
    computeWeightsArgs = createComputeWeightsArgs(),
    fitOutcomeModelArgs = createFitOutcomeModelArgs(
      modelType = modelType,
      horizons = c(8, 15),
      bootstrapArgs = createBootstrapArgs(replicates = replicates, seed = seed, ci = ci, fixHyperparameters = fixHyperparameters)
    )
  ))
}

runBootstrap <- function(data, settings, ...) {
  return(suppressMessages(do.call(bootstrapCcwAnalysis, c(list(cloneCensorWeightData = data, outcomeId = 3), settings, list(...)))))
}

# The analysis of a data set step by step, as a user would run it; returns the outcome model.
runAnalysisByHand <- function(data, settings, fixedPriorVariances = NULL) {
  population <- createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = settings$createStudyPopulationArgs)
  clones <- createClones(population, data, settings$createClonesArgs)
  on.exit(close(clones))
  models <- fitCensoringModels(clones, data, settings$fitCensoringModelsArgs, .fixedPriorVariances = fixedPriorVariances)
  outcomeTimeBasis <- if (settings$fitOutcomeModelArgs$modelType == "logistic") settings$fitOutcomeModelArgs$timeBasis else NULL
  weighted <- computeCensoringWeights(clones, models, settings$computeWeightsArgs, outcomeTimeBasis = outcomeTimeBasis)
  args <- settings$fitOutcomeModelArgs
  args$varianceEstimator <- "robust"
  return(fitOutcomeModel(weighted, args))
}

# A new CloneCensorWeightData with person i (in increasing personSeqId) copied counts[i] times, built
# in R independently of the driver: copies are numbered by person and copy, rows by copy and rowId.
copyPersons <- function(data, counts) {
  tables <- list()
  for (table in c("cohorts", "outcomes", "treatments", "covariates", "covariateRef", "analysisRef")) {
    tables[[table]] <- as.data.frame(dplyr::collect(data[[table]]))
  }
  personSeqIds <- sort(unique(tables$cohorts$personSeqId))
  cohorts <- list()
  rowIdMap <- list()
  nextRowId <- 0
  nextPersonSeqId <- 0
  for (i in seq_along(personSeqIds)) {
    rows <- tables$cohorts[tables$cohorts$personSeqId == personSeqIds[i], ]
    rows <- rows[order(rows$rowId), ]
    for (copy in seq_len(counts[i])) {
      nextPersonSeqId <- nextPersonSeqId + 1
      newRowIds <- nextRowId + seq_len(nrow(rows))
      nextRowId <- nextRowId + nrow(rows)
      rowIdMap[[length(rowIdMap) + 1]] <- data.frame(rowId = rows$rowId, newRowId = newRowIds)
      copied <- rows
      copied$rowId <- newRowIds
      copied$personSeqId <- nextPersonSeqId
      cohorts[[length(cohorts) + 1]] <- copied
    }
  }
  rowIdMap <- do.call(rbind, rowIdMap)
  remap <- function(table) {
    table <- merge(table, rowIdMap, by = "rowId")
    table$rowId <- table$newRowId
    table$newRowId <- NULL
    return(table)
  }
  return(createCloneCensorWeightData(
    cohorts = do.call(rbind, cohorts),
    outcomes = remap(tables$outcomes),
    treatments = remap(tables$treatments),
    covariates = remap(tables$covariates),
    covariateRef = tables$covariateRef,
    analysisRef = tables$analysisRef,
    metaData = list(targetId = 1, treatmentIds = c(treatment = 2))
  ))
}

# A boot object holding the full-data estimate t0 and the replicates t of one statistic.
bootObject <- function(t0, t) {
  return(structure(list(
    t0 = t0,
    t = matrix(t, ncol = 1),
    R = length(t),
    data = seq_along(t),
    seed = 1,
    statistic = function(data, indices) 0,
    sim = "ordinary",
    call = quote(boot::boot(data = 1, statistic = function(data, indices) 0, R = 1)),
    stype = "i",
    strata = rep(1, length(t)),
    weights = rep(1 / length(t), length(t))
  ), class = "boot"))
}

# Checks the standard errors and limits of `estimates` (from .bootstrapEstimates()) against sd() and
# boot::boot.ci() on the replicates, on the log scale for ratios.
expectBootCiLimits <- function(estimates, replicates, ci) {
  type <- if (ci == "percentile") "perc" else "norm"
  for (i in seq_len(nrow(estimates))) {
    estimand <- estimates$estimand[i]
    isRatio <- grepl("^(hr|or|rr_)", estimand)
    transform <- if (isRatio) log else identity
    t0 <- transform(estimates$estimate[i])
    t <- transform(replicates[, estimand])
    finite <- t[is.finite(t)]
    if (length(finite) < 2 || !is.finite(t0)) {
      next
    }
    limits <- suppressWarnings(boot::boot.ci(bootObject(t0, t), type = type))
    limits <- if (type == "perc") limits$percent[4:5] else limits$normal[2:3]
    if (isRatio) {
      expect_equal(c(estimates$lb95[i], estimates$ub95[i]), exp(limits), tolerance = 1e-12, info = estimand)
      expect_equal(estimates$seLogEstimate[i], sd(finite), tolerance = 1e-12, info = estimand)
    } else {
      expect_equal(c(estimates$lb95[i], estimates$ub95[i]), limits, tolerance = 1e-12, info = estimand)
      expect_equal(estimates$se[i], sd(finite), tolerance = 1e-12, info = estimand)
    }
  }
}

test_that("The bootstrap is reproducible with a seed and replays the counts it returns", {
  data <- bootstrapData()
  on.exit(close(data))
  settings <- bootstrapSettings(replicates = 2, seed = 123)
  result <- runBootstrap(data, settings, returnCounts = TRUE)

  counts <- result$bootstrapCounts
  expect_equal(dim(counts), c(120, 2))
  expect_equal(rownames(counts), as.character(1:120))
  expect_equal(unname(colSums(counts)), rep(120L, 2))
  # The counts are multinomial draws made with the seed before any model is fitted:
  set.seed(123)
  expect_equal(unname(counts), stats::rmultinom(2, size = 120, prob = rep(1, 120)))

  # Replaying the counts gives the same replicates, so the same seed gives the same result:
  replay <- runBootstrap(data, settings, counts = counts)
  expect_equal(replay$replicates, result$replicates)
  expect_equal(replay$estimates, result$estimates)
  expect_null(replay$bootstrapCounts)
  expect_equal(dim(result$replicates), c(2, nrow(result$estimates)))
  expect_equal(colnames(result$replicates), result$estimates$estimand)
  expect_false(anyNA(result$replicates[, grepl("^(risk|rd)_", colnames(result$replicates))]))
  expect_false(isTRUE(all.equal(result$replicates[1, ], result$replicates[2, ])))
  expect_equal(result$metaData$bootstrap$replicates, 2)
  expect_equal(result$metaData$bootstrap$failedReplicates, 0)
})

test_that("Replicates with boot.array counts equal a copy-based analysis of the resampled persons", {
  skip_if_not_installed("boot")
  data <- bootstrapData(n = 100, seed = 2)
  on.exit(close(data))
  settings <- bootstrapSettings(modelType = "cox", useCovariates = TRUE)
  set.seed(5)
  resamples <- boot::boot(data.frame(person = 1:100), function(persons, frequencies) sum(frequencies), R = 2, stype = "f")
  counts <- t(boot::boot.array(resamples))
  result <- runBootstrap(data, settings, counts = counts, returnCounts = TRUE)
  expect_equal(unname(result$bootstrapCounts), unname(counts))
  for (replicate in 1:2) {
    copied <- copyPersons(data, counts[, replicate])
    byHand <- suppressMessages(runAnalysisByHand(copied, settings))
    close(copied)
    expect_equal(unname(result$replicates[replicate, ]), byHand$estimates$estimate, tolerance = 1e-8)
  }
})

test_that("Bootstrap limits equal boot::boot.ci (perc and norm)", {
  skip_if_not_installed("boot")
  set.seed(1)
  t0 <- dplyr::tibble(
    estimand = c("hr", "risk_target_30", "rd_30", "rr_30"),
    estimate = c(0.8, 0.1, 0.02, 1.25),
    logEstimate = c(log(0.8), NA, NA, log(1.25))
  )
  for (nReplicates in c(19, 39, 50, 200)) {
    replicates <- cbind(
      hr = exp(stats::rnorm(nReplicates, log(0.8), 0.2)),
      risk_target_30 = stats::rbeta(nReplicates, 10, 90),
      rd_30 = stats::rnorm(nReplicates, 0.02, 0.01),
      rr_30 = exp(stats::rnorm(nReplicates, log(1.25), 0.3))
    )
    # Failed or infinite replicates are left out, as boot.ci does:
    replicates[3, "hr"] <- NA
    replicates[5, "rr_30"] <- Inf
    for (ci in c("percentile", "normal")) {
      estimates <- .bootstrapEstimates(t0, replicates, ci)
      expect_equal(estimates$estimate, t0$estimate)
      expectBootCiLimits(estimates, replicates, ci)
    }
  }
  # The p-value of a difference or ratio uses the bootstrap standard error:
  estimates <- .bootstrapEstimates(t0, replicates, "percentile")
  expect_equal(estimates$p[3], 2 * pnorm(-abs(0.02 / estimates$se[3])))
  expect_equal(estimates$p[1], 2 * pnorm(-abs(log(0.8) / estimates$seLogEstimate[1])))
  expect_true(is.na(estimates$p[2]))
})

test_that("The bootstrap runs for km, cox and logistic outcome models", {
  skip_if_not_installed("boot")
  data <- bootstrapData(seed = 3)
  on.exit(close(data))
  results <- list()
  for (modelType in c("km", "cox", "logistic")) {
    settings <- bootstrapSettings(modelType = modelType, replicates = 2)
    result <- runBootstrap(data, settings)
    results[[modelType]] <- result
    full <- suppressMessages(runAnalysisByHand(data, settings))
    expect_s3_class(result, "CcwOutcomeModel")
    expect_equal(result$modelType, modelType)
    expect_equal(result$varianceEstimator, "bootstrap")
    expect_equal(result$estimates$estimand, full$estimates$estimand)
    expect_equal(colnames(result$estimates), colnames(full$estimates))
    expect_equal(result$estimates$estimate, full$estimates$estimate, tolerance = 1e-10)
    expect_equal(result$coefficients, full$coefficients, tolerance = 1e-10)
    expect_equal(result$counts, full$counts)
    expect_equal(result$metaData$riskStart, 0)
    expect_equal(result$metaData$outcomeId, 3)
    expect_equal(result$metaData$args$modelType, modelType)
    expect_equal(nrow(result$replicates), 2)
    expect_false(anyNA(result$replicates[, grepl("^risk_", colnames(result$replicates))]))
    expectBootCiLimits(result$estimates, result$replicates, "percentile")
    normal <- .bootstrapEstimates(result$estimates[, c("estimand", "estimate", "logEstimate")], result$replicates, "normal")
    expectBootCiLimits(normal, result$replicates, "normal")
    if (modelType == "km") {
      expect_null(result$coefficients)
      expect_null(result$vcov)
    } else {
      expect_equal(names(result$coefficients), names(full$coefficients))
      expect_equal(dimnames(result$vcov), list(names(full$coefficients), names(full$coefficients)))
      # The natural spline basis of the default pooled logistic model is rebuilt in every replicate, so
      # only the target coefficient has a bootstrap covariance:
      invariant <- if (modelType == "logistic") "target" else names(full$coefficients)
      expect_false(anyNA(result$vcov[invariant, invariant]))
      expect_true(all(is.na(result$vcov[!rownames(result$vcov) %in% invariant, ])))
    }
    # print, coef and confint work on the result:
    expect_equal(coef(result), stats::setNames(result$estimates$estimate, result$estimates$estimand))
    expect_equal(unname(confint(result)[, 1]), result$estimates$lb95)
  }
  # The same seed draws the same persons, so the cox bootstrap bootstraps the hazard ratio and exactly
  # the Kaplan-Meier rows of the km bootstrap:
  cox <- results$cox
  km <- results$km
  expect_equal(colnames(cox$replicates), c("hr", colnames(km$replicates)))
  expect_false(anyNA(cox$replicates))
  expect_identical(cox$replicates[, -1], km$replicates)
  expect_identical(cox$estimates[-1, ], km$estimates)
  expect_false(is.na(cox$estimates$seLogEstimate[1]))
})

test_that("A cox replicate without outcomes keeps its Kaplan-Meier rows and leaves out its hazard ratio", {
  data <- bootstrapData(seed = 5)
  on.exit(close(data))
  settings <- bootstrapSettings(modelType = "cox", replicates = 2)
  # The second replicate draws only the persons without the outcome, so its Cox model has no events:
  cohorts <- dplyr::collect(data$cohorts)
  withOutcome <- unique(cohorts$personSeqId[cohorts$rowId %in% dplyr::collect(data$outcomes)$rowId])
  personSeqIds <- sort(unique(cohorts$personSeqId))
  counts <- cbind(rep(1L, 120), as.integer(!personSeqIds %in% withOutcome))
  expect_gt(length(withOutcome), 0)
  warnings <- capture_warnings(result <- runBootstrap(data, settings, counts = counts))
  expect_length(warnings, 1)
  # Its risk ratios (0 / 0) are missing too:
  expect_match(warnings, "missing or infinite in bootstrap replicates that did not fail.*\\(finite replicates of 2\\): hr \\(1\\), rr_8 \\(1\\), rr_15 \\(1\\)$")
  expect_equal(result$metaData$bootstrap$failedReplicates, 0)
  expect_true(is.na(result$replicates[2, "hr"]))
  risks <- grepl("^risk_", colnames(result$replicates))
  expect_equal(unname(result$replicates[2, risks]), rep(0, sum(risks)))
  expect_false(anyNA(result$replicates[2, grepl("^(risk|rd|rmst|rmtl)", colnames(result$replicates))]))
})

test_that("With fixHyperparameters, replicates use the full-data prior variances without cross-validation", {
  data <- bootstrapData(seed = 4, initiationCoefficient = 0.8)
  on.exit(close(data))
  settings <- bootstrapSettings(
    modelType = "km",
    useCovariates = TRUE,
    prior = Cyclops::createPrior("laplace", exclude = c(0), useCrossValidation = TRUE),
    replicates = 2
  )
  # The full-data censoring models, with variances chosen by cross-validation:
  population <- suppressMessages(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = settings$createStudyPopulationArgs))
  clones <- suppressMessages(createClones(population, data, settings$createClonesArgs))
  fullModels <- suppressMessages(fitCensoringModels(clones, data, settings$fitCensoringModelsArgs))
  close(clones)
  # With postGraceModel = "binary" (contract section 11.2), strategy 2 also has the binary model of its
  # weights after the grace period, shared with strategy 1:
  variances <- list("1" = fullModels[[1]]$priorVariance, "2" = fullModels[[2]]$priorVariance, "2.postGrace" = fullModels[[2]]$postGrace$priorVariance)
  expect_equal(variances[["2.postGrace"]], variances[["1"]])
  expect_false(anyNA(unlist(variances)))
  expect_true(all(fullModels[[1]]$coefficients[-1] != 0))

  realFitCensoringModels <- fitCensoringModels
  calls <- list()
  local_mocked_bindings(fitCensoringModels = function(clones, cloneCensorWeightData, fitCensoringModelsArgs, .fixedPriorVariances = NULL) {
    models <- realFitCensoringModels(clones, cloneCensorWeightData, fitCensoringModelsArgs, .fixedPriorVariances)
    calls[[length(calls) + 1]] <<- list(fixed = .fixedPriorVariances, fitted = lapply(models, function(model) model$priorVariance))
    return(models)
  })
  result <- runBootstrap(data, settings, censoringModels = fullModels)
  expect_equal(result$metaData$bootstrap$fixedPriorVariances, variances)
  # Only the replicates fit censoring models (the full-data models were given):
  expect_length(calls, 2)
  for (call in calls) {
    expect_equal(call$fixed, variances)
    expect_equal(unname(unlist(call$fitted)), unname(unlist(variances[c("1", "2")])))
  }

  calls <- list()
  settings$fitOutcomeModelArgs$bootstrapArgs <- createBootstrapArgs(replicates = 1, fixHyperparameters = FALSE)
  result <- runBootstrap(data, settings, censoringModels = fullModels)
  expect_null(result$metaData$bootstrap$fixedPriorVariances)
  expect_length(calls, 1)
  expect_null(calls[[1]]$fixed)
})

test_that("Failed replicates are NA, left out with a warning, and inputs are checked", {
  data <- bootstrapData(seed = 5)
  on.exit(close(data))
  settings <- bootstrapSettings(replicates = 2)
  # The second replicate draws nobody, so its analysis fails:
  counts <- cbind(rep(1L, 120), rep(0L, 120))
  expect_warning(result <- runBootstrap(data, settings, counts = counts), "failed in 1 of 2 bootstrap replicates")
  expect_true(all(is.na(result$replicates[2, ])))
  expect_equal(result$metaData$bootstrap$failedReplicates, 1)
  # Drawing everyone once reproduces the full data:
  expect_equal(unname(result$replicates[1, ]), result$estimates$estimate, tolerance = 1e-10)
  expect_true(all(is.na(result$estimates$se)))

  expect_error(runBootstrap(data, settings, counts = matrix(1L, 119, 2)), "one row per person")
  rownamed <- matrix(1L, 120, 1, dimnames = list(as.character(120:1), NULL))
  expect_error(runBootstrap(data, settings, counts = rownamed), "row names")
  expect_error(runBootstrap(data, settings, counts = matrix(-1L, 120, 1)), "counts")
  expect_error(runBootstrap(data, settings[-1]), "createStudyPopulationArgs")
  weighted <- createM4WeightedClones(createM4Fixture(n = 20))
  on.exit(close(weighted), add = TRUE)
  expect_error(fitOutcomeModel(weighted, createFitOutcomeModelArgs(varianceEstimator = "bootstrap")), "bootstrapCcwAnalysis")
})

test_that("Coefficients whose meaning changes between replicates get no bootstrap covariance", {
  names <- c("(Intercept)", "target", "time1", "target:time1", "death:(Intercept)", "death:target", "death:time1", "death:target:time1")
  args <- function(...) createFitOutcomeModelArgs(modelType = "logistic", ...)
  expect_equal(.replicateInvariantCoefficients(names, args()), c("target", "death:target"))
  expect_equal(.replicateInvariantCoefficients(names, args(armByTimeInteraction = TRUE)), c("target", "death:target"))
  for (timeBasis in list(createTimeBasis("indicator"), createTimeBasis("linear"), createTimeBasis("polynomial", raw = TRUE))) {
    expect_equal(.replicateInvariantCoefficients(names, args(timeBasis = timeBasis, armByTimeInteraction = TRUE)), names)
  }
  orthogonal <- createTimeBasis("polynomial", raw = FALSE)
  expect_equal(.replicateInvariantCoefficients(names, args(timeBasis = orthogonal)), c("target", "death:target"))
  expect_equal(.replicateInvariantCoefficients(names, args(timeBasis = orthogonal, armByTimeInteraction = TRUE)), character(0))
  expect_equal(.replicateInvariantCoefficients("target", createFitOutcomeModelArgs(modelType = "cox")), "target")
})

test_that("A replicate whose outcome model cannot be fitted counts as failed", {
  data <- bootstrapData(seed = 5)
  on.exit(close(data))
  settings <- bootstrapSettings(modelType = "logistic", replicates = 2)
  # The second replicate draws only the persons without the outcome, so its pooled logistic outcome
  # model has no events and gives no estimate, without an error:
  cohorts <- dplyr::collect(data$cohorts)
  withOutcome <- unique(cohorts$personSeqId[cohorts$rowId %in% dplyr::collect(data$outcomes)$rowId])
  personSeqIds <- sort(unique(cohorts$personSeqId))
  counts <- cbind(rep(1L, 120), as.integer(!personSeqIds %in% withOutcome))
  expect_gt(length(withOutcome), 0)
  warnings <- capture_warnings(result <- runBootstrap(data, settings, counts = counts))
  expect_length(warnings, 1)
  expect_match(warnings, "failed in 1 of 2 bootstrap replicates, which are left out: replicate 2: no finite estimate \\(No outcome events")
  expect_equal(result$metaData$bootstrap$failedReplicates, 1)
  expect_true(all(is.na(result$replicates[2, ])))
  expect_equal(result$metaData$bootstrap$finiteReplicates, stats::setNames(rep(1, nrow(result$estimates)), result$estimates$estimand))
})

test_that("Estimates missing in replicates that did not fail are counted, with a warning", {
  data <- bootstrapData(seed = 6)
  on.exit(close(data))
  settings <- bootstrapSettings(replicates = 3)
  realPoint <- .fitOutcomeModelPoint
  calls <- 0
  # The second replicate (third fit) loses its risks, as when its death model cannot be fitted:
  local_mocked_bindings(.fitOutcomeModelPoint = function(weightedClones, fitOutcomeModelArgs, ...) {
    calls <<- calls + 1
    result <- realPoint(weightedClones, fitOutcomeModelArgs, ...)
    if (calls == 3) {
      result$estimates$estimate[grepl("^(risk|rd|rr)_", result$estimates$estimand)] <- NA
      warning("Fitting the death model failed: test")
    }
    return(result)
  })
  warnings <- capture_warnings(result <- runBootstrap(data, settings))
  # The replicate's own warning is not shown; the summary is:
  expect_length(warnings, 1)
  expect_match(warnings, "missing or infinite in bootstrap replicates that did not fail.*\\(finite replicates of 3\\): risk_target_8 \\(2\\)")
  expect_equal(result$metaData$bootstrap$failedReplicates, 0)
  missing <- grepl("^(risk|rd|rr)_", result$estimates$estimand)
  expect_equal(unname(result$metaData$bootstrap$finiteReplicates), ifelse(missing, 2, 3))
  expect_true(all(is.na(result$replicates[2, missing])))
  expect_false(anyNA(result$replicates[2, !missing]))
})

test_that("Replicates fit the censoring models without the high-correlation check", {
  data <- bootstrapData(seed = 7)
  on.exit(close(data))
  settings <- bootstrapSettings(useCovariates = TRUE, replicates = 2)
  expect_true(settings$fitCensoringModelsArgs$errorOnHighCorrelation)
  realFitCensoringModels <- fitCensoringModels
  checks <- logical(0)
  local_mocked_bindings(fitCensoringModels = function(clones, cloneCensorWeightData, fitCensoringModelsArgs, .fixedPriorVariances = NULL) {
    checks <<- c(checks, fitCensoringModelsArgs$errorOnHighCorrelation)
    return(realFitCensoringModels(clones, cloneCensorWeightData, fitCensoringModelsArgs, .fixedPriorVariances))
  })
  result <- runBootstrap(data, settings)
  # The full data first, with the check; then the replicates, without it:
  expect_equal(checks, c(TRUE, FALSE, FALSE))
  expect_true(settings$fitCensoringModelsArgs$errorOnHighCorrelation)
})

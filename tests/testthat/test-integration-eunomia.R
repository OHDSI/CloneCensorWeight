library(testthat)

# End-to-end test on Eunomia (GiBleed, SQLite; a Phase 1 exit item, plan section 8): the NSAID design of
# test-extraction.R (eligible cohort 4, new users of NSAIDs; treatment cohort 1, celecoxib; comparator
# treatment cohort 2, diclofenac; outcome cohort 3, GI bleed) through getDbCloneCensorWeightData(
# connectionDetails = ) -> createStudyPopulation -> createClones -> fitCensoringModels ->
# computeCensoringWeights -> computeCovariateBalance -> fitOutcomeModel ("cox" and "km") ->
# bootstrapCcwAnalysis (5 replicates), for "start celecoxib within 30 days" vs "do not start celecoxib
# within 30 days", with risks at 180 and 365 days.
#
# In Eunomia every NSAID new user starts celecoxib or diclofenac on the index day, so the grace period is
# degenerate (celecoxib-first vs diclofenac-first, with 30 days of shared follow-up). It still exercises
# every step: the binary model at G, the initiation-hazard model on days 0 to G, and the split at G.
# Eunomia has no deaths, so the competing-event code runs with none. The covariates are gender and age
# group, as in extractionCovariateSettings() (helper-extraction.R) but with windows ending on day -1, before
# the day-0 initiation decision (Phase1Contracts section 11.6). They are fitted without a prior: they are
# only weakly associated with starting celecoxib, and cross-validation shrinks them all to zero
# (CohortMethod::createPs then warns).
#
# Asserted: structure, finiteness, the identity of the default report's risks with "km", counts that follow
# from the design, and that a rerun gives the same results: a second extraction gives the same tables, and
# the bootstrap driver, which first runs the whole analysis on the full data, gives the same estimates from
# it. Rerun estimates are compared to 1e-12 (relative), allowing for the order of floating-point sums in
# multithreaded DuckDB queries. Runtime is about 50 s, mostly the fixed cost of each bootstrap replicate.

skip_on_cran()
skip_if_not_installed("Eunomia")
skip_if_not_installed("DatabaseConnector")
skip_if_not_installed("SqlRender")
skip_if_not_installed("RSQLite")

test_that("The Eunomia NSAID design runs end to end and reproduces on rerun", {
  databaseFile <- tempfile(fileext = ".sqlite")
  connectionDetails <- tryCatch(
    suppressMessages(Eunomia::getEunomiaConnectionDetails(databaseFile = databaseFile)),
    error = function(e) NULL
  )
  if (is.null(connectionDetails)) {
    skip("Eunomia GiBleed data could not be loaded")
  }
  on.exit(unlink(databaseFile))
  suppressMessages(Eunomia::createCohorts(connectionDetails))
  horizons <- c(180, 365)
  graceEnd <- 30

  # Extraction ------------------------------------------------------------------------------------
  extract <- function() {
    data <- NULL
    # FeatureExtraction prints a progress bar.
    utils::capture.output(data <- suppressMessages(getDbCloneCensorWeightData(
      connectionDetails = connectionDetails,
      cdmDatabaseSchema = "main",
      targetId = 4,
      treatmentId = 1,
      comparatorTreatmentId = 2,
      outcomeIds = 3,
      getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(
        covariateSettings = FeatureExtraction::createCovariateSettings(
          useDemographicsGender = TRUE,
          useDemographicsAgeGroup = TRUE,
          endDays = -1
        )
      )
    )))
    return(data)
  }
  data <- extract()
  on.exit(close(data), add = TRUE)
  expect_true(isCloneCensorWeightData(data))
  expect_equal(attr(data, "metaData")$treatmentIds, c(treatment = 1, comparatorTreatment = 2))
  nCohort <- dplyr::pull(dplyr::count(data$cohorts))
  expect_gt(nCohort, 1000)

  # Study population and clones -------------------------------------------------------------------
  createStudyPopulationArgs <- CohortMethod::createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = TRUE,
    riskWindowStart = 0,
    riskWindowEnd = max(horizons),
    endAnchor = "cohort start"
  )
  population <- suppressMessages(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = createStudyPopulationArgs))
  n <- nrow(population)
  expect_gt(n, 1000)
  expect_lte(n, nCohort)
  expect_true(all(population$riskStart == 0 & population$riskEnd <= max(horizons)))
  # One entry per person (firstOnly):
  expect_equal(length(unique(population$personSeqId)), n)

  createClonesArgs <- createCreateClonesArgs(strategies = list(
    createTreatmentStrategy(1, "Start celecoxib within 30 days", "initiateWithin", graceEnd = graceEnd),
    createTreatmentStrategy(2, "Do not start celecoxib within 30 days", "notInitiateWithin", graceEnd = graceEnd)
  ))
  clones <- suppressMessages(createClones(population, data, createClonesArgs))
  on.exit(close(clones), add = TRUE)
  persons <- as.data.frame(dplyr::collect(clones$persons))
  cloneTable <- as.data.frame(dplyr::collect(clones$clones))
  expect_equal(nrow(persons), n)
  # Every person starts exactly one of the two NSAIDs, on the index day:
  celecoxib <- persons$rowId[persons$treatmentDay %in% 0]
  diclofenac <- persons$rowId[persons$comparatorTreatmentDay %in% 0]
  expect_setequal(c(celecoxib, diclofenac), persons$rowId)
  expect_length(intersect(celecoxib, diclofenac), 0)
  expect_gt(length(celecoxib), 0)
  expect_gt(length(diclofenac), 0)
  # Every clone has follow-up (riskStart = 0, "eventsFirst"), and there are no deaths:
  expect_equal(as.vector(table(cloneTable$strategyId)), c(n, n))
  expect_true(all(cloneTable$survivalTime >= 1 & cloneTable$endDay <= max(horizons)))
  expect_true(all(cloneTable$endReason %in% c("outcome", "deviation", "observationEnd", "riskEnd")))
  expect_equal(cloneTable$outcomeCount, as.integer(cloneTable$endReason == "outcome"))
  attrition <- attr(clones, "metaData")$attrition
  expect_equal(utils::tail(attrition$description, 2), c("Clones created", "Clones with follow-up"))
  expect_equal(utils::tail(attrition$targetPersons, 2), c(n, n))
  # Strategy 2 censors celecoxib starters on day 0 (unless the outcome comes first that day), strategy 1
  # censors diclofenac starters at day 30, and each person follows one strategy to the end:
  deviators2 <- cloneTable[cloneTable$strategyId == 2 & cloneTable$rowId %in% celecoxib, ]
  expect_true(all(deviators2$endDay == 0 & deviators2$endReason %in% c("deviation", "outcome")))
  deviators1 <- cloneTable[cloneTable$strategyId == 1 & cloneTable$rowId %in% diclofenac, ]
  expect_true(all(deviators1$endDay <= graceEnd))
  expect_true(all(deviators1$endDay[deviators1$endReason == "deviation"] == graceEnd))
  expect_gt(sum(deviators1$endReason == "deviation"), 0)
  followers <- cloneTable[(cloneTable$strategyId == 1 & cloneTable$rowId %in% celecoxib) | (cloneTable$strategyId == 2 & cloneTable$rowId %in% diclofenac), ]
  expect_equal(nrow(followers), n)
  expect_false(any(followers$endReason == "deviation"))
  # So every outcome in the risk window is counted in at least one strategy:
  expect_setequal(unique(cloneTable$rowId[cloneTable$outcomeCount == 1]), population$rowId[population$outcomeCount > 0])
  expect_gt(sum(cloneTable$outcomeCount[cloneTable$strategyId == 1]), 0)
  expect_gt(sum(cloneTable$outcomeCount[cloneTable$strategyId == 2]), 0)

  # Censoring models, weights and balance ---------------------------------------------------------
  fitCensoringModelsArgs <- createFitCensoringModelsArgs(prior = Cyclops::createPrior("none"))
  models <- suppressMessages(fitCensoringModels(clones, data, fitCensoringModelsArgs))
  expect_s3_class(models, "CcwCensoringModels")
  expect_equal(vapply(models, function(model) model$modelType, character(1)), c("binary", "logistic"))
  expect_equal(vapply(models, function(model) model$horizon, numeric(1)), c(graceEnd, graceEnd))
  expect_equal(vapply(models, function(model) model$diagnostics$status, character(1)), c("OK", "OK"))
  expect_true(all(is.finite(unlist(lapply(models, function(model) model$coefficients)))))

  computeWeightsArgs <- createComputeWeightsArgs()
  weighted <- suppressMessages(computeCensoringWeights(clones, models, computeWeightsArgs))
  on.exit(close(weighted), add = TRUE)
  segments <- as.data.frame(dplyr::collect(weighted$segments))
  expect_true(all(is.finite(segments$weight) & segments$weight > 0))
  expect_setequal(unique(segments$cloneRowId), cloneTable$cloneRowId)
  expect_equal(sum(segments$outcome), sum(cloneTable$outcomeCount))
  expect_equal(sum(segments$competing), 0)
  weightsMetaData <- attr(weighted, "metaData")$weights
  expect_true(weightsMetaData$stabilized)
  expect_true(all(segments$weight <= max(weightsMetaData$caps$cap)))

  balance <- suppressMessages(computeCovariateBalance(weighted, data))
  expect_gt(nrow(balance), 0)
  expect_true(all(c("covariateId", "beforeMatchingStdDiff", "afterMatchingStdDiff") %in% colnames(balance)))
  expect_true(all(is.finite(balance$beforeMatchingStdDiff) & is.finite(balance$afterMatchingStdDiff)))

  # Outcome models --------------------------------------------------------------------------------
  coxArgs <- createFitOutcomeModelArgs(horizons = horizons)
  cox <- suppressMessages(fitOutcomeModel(weighted, coxArgs))
  km <- suppressMessages(fitOutcomeModel(weighted, createFitOutcomeModelArgs(modelType = "km", horizons = horizons)))
  expect_s3_class(cox, "CcwOutcomeModel")
  expect_equal(cox$estimates$estimand[1], "hr")
  for (model in list(cox, km)) {
    estimates <- model$estimates
    expect_true(all(is.finite(estimates$estimate) & is.finite(estimates$lb95) & is.finite(estimates$ub95)))
    expect_true(all(estimates$lb95 < estimates$estimate & estimates$estimate < estimates$ub95))
    se <- ifelse(is.na(estimates$se), estimates$seLogEstimate, estimates$se)
    expect_true(all(is.finite(se) & se > 0))
    for (h in horizons) {
      risks <- estimates$estimate[estimates$estimand %in% sprintf("risk_%s_%d", c("target", "comparator"), h)]
      expect_length(risks, 2)
      expect_true(all(risks > 0 & risks < 1))
    }
    counts <- model$counts
    expect_equal(counts$strategyId, 1:2)
    expect_equal(counts$clones, c(n, n))
    expect_equal(counts$persons, c(n, n))
    expect_equal(counts$outcomes, as.vector(tapply(cloneTable$outcomeCount, cloneTable$strategyId, sum)))
    expect_equal(counts$competingEvents, c(0L, 0L))
    expect_equal(counts$personDays, as.vector(tapply(as.numeric(cloneTable$survivalTime), cloneTable$strategyId, sum)))
  }
  # The default report's risks, RD, RR and RMTL are those of "km" (Phase1Contracts section 10.5):
  expect_identical(cox$estimates[-1, ], km$estimates)
  expect_identical(cox$counts, km$counts)

  # Rerun and bootstrap ---------------------------------------------------------------------------
  # A second extraction gives the same data:
  rerunData <- extract()
  on.exit(close(rerunData), add = TRUE)
  for (table in c("cohorts", "outcomes", "treatments", "covariates", "covariateRef")) {
    expect_equal(sortExtractionRows(dplyr::collect(rerunData[[table]])), sortExtractionRows(dplyr::collect(data[[table]])), info = table)
  }
  # The bootstrap driver runs the whole analysis on the full (re-extracted) data first: a rerun of every
  # step from the study population to the outcome model.
  bootstrapArgs <- createFitOutcomeModelArgs(horizons = horizons, bootstrapArgs = createBootstrapArgs(replicates = 5))
  bootstrap <- suppressMessages(bootstrapCcwAnalysis(
    cloneCensorWeightData = rerunData,
    outcomeId = 3,
    createStudyPopulationArgs = createStudyPopulationArgs,
    createClonesArgs = createClonesArgs,
    fitCensoringModelsArgs = fitCensoringModelsArgs,
    computeWeightsArgs = computeWeightsArgs,
    fitOutcomeModelArgs = bootstrapArgs,
    returnCounts = TRUE
  ))
  expect_s3_class(bootstrap, "CcwOutcomeModel")
  expect_equal(bootstrap$varianceEstimator, "bootstrap")
  expect_equal(bootstrap$metaData$bootstrap$failedReplicates, 0)
  expect_equal(dim(bootstrap$bootstrapCounts), c(nCohort, 5))
  expect_equal(unname(colSums(bootstrap$bootstrapCounts)), rep(nCohort, 5))
  expect_equal(dim(bootstrap$replicates), c(5, nrow(cox$estimates)))
  expect_equal(colnames(bootstrap$replicates), cox$estimates$estimand)
  expect_true(all(is.finite(bootstrap$replicates)))
  expect_true(all(is.finite(bootstrap$estimates$lb95) & is.finite(bootstrap$estimates$ub95)))
  bootstrapSe <- ifelse(is.na(bootstrap$estimates$se), bootstrap$estimates$seLogEstimate, bootstrap$estimates$se)
  expect_true(all(is.finite(bootstrapSe) & bootstrapSe > 0))
  # The rerun gives the same estimates:
  expect_equal(bootstrap$estimates$estimand, cox$estimates$estimand)
  expect_equal(bootstrap$estimates$estimate, cox$estimates$estimate, tolerance = 1e-12)
  expect_equal(bootstrap$coefficients, cox$coefficients, tolerance = 1e-12)
  expect_identical(bootstrap$counts, cox$counts)
})

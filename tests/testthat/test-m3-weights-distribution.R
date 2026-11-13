library(testthat)

# Tests of computeWeightDistribution() (Phase1Contracts section 11.7): the weight it reports for a clone
# on day d is the exact W(d), also on the days inside a segment that no event splits (a segment carries
# the weight of its last day), for every strategy type, post-grace model, stabilization and truncation.
# The references write W(d) out from the fitted censoring models (refDayWeights(), helper-m3-weights.R).

# Simulated data, clones and censoring models (with covariates and no prior, so the weights vary by
# person) of initiateWithin(5) vs `comparatorType`, with few events, so that some segments span several
# days. The caller closes `data` and `clones`.
distributionAnalysis <- function(comparatorType = "notInitiateWithin", postGraceModel = "binary", sameDayOrdering = "eventsFirst", hazardModelType = "logistic") {
  data <- simulateCloneCensorWeightData(
    n = 300, maxDays = 20, graceEnd = 5, initiationHazard = 0.05, outcomeHazard = 0.0005, deathHazard = 0.0002, seed = 1
  )
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 20, endAnchor = "cohort start")
  ))
  strategies <- list(
    createTreatmentStrategy(1, "Initiate", "initiateWithin", graceEnd = 5),
    createTreatmentStrategy(2, "Comparator", comparatorType, graceEnd = if (comparatorType == "neverInitiate") NULL else 5)
  )
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = strategies, sameDayOrdering = sameDayOrdering)))
  models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(
    hazardModelType = hazardModelType, prior = Cyclops::createPrior("none"), postGraceModel = postGraceModel
  )))
  return(list(data = data, clones = clones, models = models))
}

# The distribution computeWeightDistribution() should report on `days`, from W(d) written out by hand,
# stabilized or not, and truncated at the caps of the weighted clones.
handDistribution <- function(analysis, weighted, days, stabilized) {
  sameDayOrdering <- attr(analysis$clones, "metaData")$sameDayOrdering
  cloneDays <- refDayWeights(analysis$models, collectClones(analysis$clones), days, sameDayOrdering)
  w <- if (stabilized) cloneDays$weight * cloneDays$numerator else cloneDays$weight
  caps <- attr(weighted, "metaData")$weights$caps
  cap <- caps$cap[match(cloneDays$strategyId, caps$strategyId)]
  w <- ifelse(is.na(cap), w, pmin(w, cap))
  return(list(distribution = refWeightDistribution(cloneDays, w), cloneDays = cloneDays, w = w))
}

asPlainDataFrame <- function(distribution) {
  distribution <- as.data.frame(distribution)
  attr(distribution, "strategyLabels") <- NULL
  return(distribution)
}

test_that("computeWeightDistribution reports the exact W(d) for every strategy type, post-grace model, stabilization and truncation", {
  configurations <- list(
    list(comparatorType = "notInitiateWithin", postGraceModel = "binary"),
    list(comparatorType = "notInitiateWithin", postGraceModel = "hazard"),
    list(comparatorType = "neverInitiate"),
    list(comparatorType = "notInitiateWithin", postGraceModel = "binary", sameDayOrdering = "deviationFirst", hazardModelType = "cox")
  )
  truncations <- list(none = createWeightTruncation("none"), percentile = createWeightTruncation("percentile", 0.9))
  # Days 0 to tau + 1 = 6, and later days (where the neverInitiate weight still changes):
  days <- c(0:6, 10, 20)
  for (configuration in configurations) {
    analysis <- do.call(distributionAnalysis, configuration)
    eventsFirst <- is.null(configuration$sameDayOrdering)
    for (stabilized in c(FALSE, TRUE)) {
      for (truncation in names(truncations)) {
        info <- paste(c(unlist(configuration), if (stabilized) "stabilized" else "unstabilized", truncation), collapse = " ")
        weighted <- suppressMessages(computeCensoringWeights(
          analysis$clones, analysis$models, createComputeWeightsArgs(stabilized = stabilized, truncation = truncations[[truncation]])
        ))
        distribution <- computeWeightDistribution(weighted, days = days)
        hand <- handDistribution(analysis, weighted, days, stabilized)
        expect_equal(asPlainDataFrame(distribution), hand$distribution, tolerance = 1e-10, info = info)
        expect_equal(sort(unique(distribution$day)), days, info = info)
        if (truncation == "percentile") {
          untruncated <- if (stabilized) hand$cloneDays$weight * hand$cloneDays$numerator else hand$cloneDays$weight
          expect_true(any(hand$w < untruncated - 1e-12), info = info)
        }
        if (eventsFirst) {
          # W(0) = 1 in both strategies:
          dayZero <- distribution[distribution$day == 0, ]
          expect_equal(dayZero$strategyId, 1:2, info = info)
          expect_equal(dayZero$minWeight, c(1, 1), info = info)
          expect_equal(dayZero$maxWeight, c(1, 1), info = info)
          expect_equal(dayZero$effectiveSampleSize, dayZero$clones, info = info)
        }
        if (truncation == "none") {
          # The weight of the segment that contains a day is the weight of the segment's last day, which
          # is not W(d) on some of these days:
          segments <- collectSegments(weighted)
          t <- hand$cloneDays$day + 1
          segmentWeight <- vapply(seq_along(t), function(i) {
            segments$weight[segments$cloneRowId == hand$cloneDays$cloneRowId[i] & segments$startTime < t[i] & t[i] <= segments$stopTime]
          }, numeric(1))
          expect_true(any(abs(segmentWeight - hand$w) > 1e-6), info = info)
        }
        close(weighted)
      }
    }
    close(analysis$clones)
    close(analysis$data)
  }
})

test_that("computeWeightDistribution reports the exact W(d) on its default days and leaves no tables behind", {
  analysis <- distributionAnalysis("neverInitiate")
  on.exit({
    close(analysis$clones)
    close(analysis$data)
  })
  weighted <- suppressMessages(computeCensoringWeights(analysis$clones, analysis$models, createComputeWeightsArgs()))
  on.exit(close(weighted), add = TRUE)
  tables <- names(weighted)
  distribution <- computeWeightDistribution(weighted)
  segments <- collectSegments(weighted)
  days <- sort(unique(as.integer(segments$stopTime - 1)))
  expect_equal(sort(unique(distribution$day)), days)
  hand <- handDistribution(analysis, weighted, days, stabilized = TRUE)
  expect_equal(asPlainDataFrame(distribution), hand$distribution, tolerance = 1e-10)
  expect_equal(names(weighted), tables)
  expect_equal(attr(distribution, "strategyLabels"), c("Initiate", "Comparator"))

  # Weighted clones without the censoring models of their weights:
  withoutModels <- weighted
  attr(withoutModels, "metaData")$weights$censoringModels <- NULL
  expect_error(computeWeightDistribution(withoutModels), "censoring models")
})

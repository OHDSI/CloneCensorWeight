library(testthat)

# Contract section 11.4: computeCensoringModelCalibration() compares observed and predicted
# P(start by G | X) by decile of the prediction, for the binary model and for 1 - S(G | X) of the
# initiation hazard model, per strategy with a grace period.

test_that("computeCensoringModelCalibration equals a hand computation", {
  fixture <- createM3Fixture(150, riskStart = 0L, seed = 51)
  data <- createM3Data(fixture)
  on.exit(close(data))
  configurations <- list(
    list(sameDayOrdering = "eventsFirst", postGraceModel = "binary", hazardModelType = "logistic"),
    list(sameDayOrdering = "deviationFirst", postGraceModel = "binary", hazardModelType = "cox"),
    list(sameDayOrdering = "deviationFirst", postGraceModel = "hazard", hazardModelType = "logistic")
  )
  for (configuration in configurations) {
    sameDayOrdering <- configuration$sameDayOrdering
    postGraceModel <- configuration$postGraceModel
    hazardModelType <- configuration$hazardModelType
    clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = 7), sameDayOrdering)
    reference <- refBinaryRiskSet(fixture$persons, 7, sameDayOrdering)
    asPlain <- function(x) data.frame(group = x$group, persons = x$persons, predicted = x$predicted, observed = x$observed)
    # Groups of about equal size by rank of the prediction; ties share the group of the first of them:
    handCalibration <- function(predicted, numberOfGroups) {
      group <- ceiling(rank(predicted, ties.method = "min") * numberOfGroups / length(predicted))
      group <- match(group, sort(unique(group)))
      data.frame(
        group = seq_len(max(group)),
        persons = as.integer(table(group)),
        predicted = as.vector(tapply(predicted, group, mean)),
        observed = as.vector(tapply(reference$y, group, mean))
      )
    }
    models <- suppressMessages(fitCensoringModels(clones, data, m3Args(hazardModelType, postGraceModel = postGraceModel)))
    calibration <- computeCensoringModelCalibration(clones, models, numberOfGroups = 5)
    expect_equal(colnames(calibration), c("strategyId", "model", "covariateTimeInteraction", "grouping", "covariateId", "group", "persons", "predicted", "observed"))
    expect_equal(attr(calibration, "strategyLabels"), c("initiateWithin 1", "notInitiateWithin 2"))
    expect_equal(attr(calibration, "graceEnds"), c(7L, 7L))
    expect_equal(unique(calibration$model[calibration$strategyId == 2]), if (postGraceModel == "binary") c("hazard", "binary") else "hazard")

    binary <- models[[1]]$predictions$propensityScore[match(reference$rowId, models[[1]]$predictions$rowId)]
    expected <- handCalibration(binary, 5)
    actual <- calibration[calibration$strategyId == 1 & calibration$model == "binary" & calibration$grouping == "prediction", ]
    expect_equal(asPlain(actual), expected)
    # The intercept of the maximum likelihood fit reproduces the observed share overall:
    expect_equal(sum(actual$predicted * actual$persons), sum(reference$y), tolerance = 1e-6)
    if (postGraceModel == "binary") {
      actual <- calibration[calibration$strategyId == 2 & calibration$model == "binary" & calibration$grouping == "prediction", ]
      expect_equal(asPlain(actual), expected)
    }

    # 1 - S(7 | X) by hand: logistic, prod of 1 - h over the daily intervals; Cox, exp(-Lambda0(8) e^lp).
    hazard <- models[[2]]
    linearPredictor <- hazard$predictions$linearPredictor[match(reference$rowId, hazard$predictions$rowId)]
    if (hazardModelType == "logistic") {
      survival <- sapply(linearPredictor, function(lp) prod(1 - plogis(hazard$intervals$baseline + lp)))
    } else {
      lambda <- hazard$baselineHazard$cumHazard[max(which(hazard$baselineHazard$time <= 8))]
      survival <- exp(-lambda * exp(linearPredictor))
    }
    expected <- handCalibration(1 - survival, 5)
    actual <- calibration[calibration$strategyId == 2 & calibration$model == "hazard" & calibration$grouping == "prediction", ]
    expect_equal(asPlain(actual), expected, tolerance = 1e-10)
    close(clones)
  }

  # Strategies without a grace period are left out:
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "neverInitiate", graceEnd = 7), "eventsFirst")
  on.exit(close(clones), add = TRUE)
  models <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic")))
  calibration <- computeCensoringModelCalibration(clones, models)
  expect_equal(unique(calibration$strategyId), 1L)
  expect_equal(attr(calibration, "graceEnds"), c(7L, NA_integer_))
  groups <- calibration[calibration$grouping == "prediction", ]
  expect_equal(sum(groups$persons), nrow(refBinaryRiskSet(fixture$persons, 7, "eventsFirst")))
  expect_equal(calibration$persons[calibration$grouping == "overall"], sum(groups$persons))
  expect_lte(nrow(groups), 10)

  # Input checks:
  expect_error(computeCensoringModelCalibration(clones, models[1]), "CcwCensoringModels")
  expect_error(computeCensoringModelCalibration(clones, models, numberOfGroups = 0), "numberOfGroups")
  expect_error(computeCensoringModelCalibration(data, models), "CcwClones")
})

test_that("The calibration diagnostic shows a miscalibrated initiation hazard model", {
  # One binary and one normal covariate; the binary one makes persons start on day 0 (log odds ratio 3)
  # and less often later (-1.5). Without covariate-by-time interactions the hazard model is
  # misspecified, and its 1 - S(G | X) is miscalibrated; the binary model at G is close to calibrated.
  data <- simulateCloneCensorWeightData(
    n = 20000, nCovariates = 2, dgp = "a", graceEnd = 7, maxDays = 20, binaryPrevalence = 0.4,
    initiationHazard = 0.06, initiationCoefficients = c(-1.5, 0.3), initiationHazardDay0 = 0.1,
    initiationCoefficientsDay0 = c(3, 0.3), outcomeHazard = 0.005, deathHazard = 0.002,
    outcomeCoefficients = c(0.5, 0.2), deathCoefficients = c(0.3, 0.1), seed = 2
  )
  on.exit(close(data))
  population <- suppressMessages(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
    riskWindowEnd = 20, endAnchor = "cohort start"
  )))
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = list(
    createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
    createTreatmentStrategy(2, "Do not initiate within 7 days", "notInitiateWithin", graceEnd = 7)
  ))))
  on.exit(close(clones), add = TRUE)
  zScores <- function(calibration) {
    se <- sqrt(calibration$observed * (1 - calibration$observed) / calibration$persons)
    return((calibration$observed - calibration$predicted) / se)
  }
  maxAbsZ <- list()
  for (covariateTimeInteraction in c("none", "day0")) {
    models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(
      prior = Cyclops::createPrior("none"),
      covariateTimeInteraction = covariateTimeInteraction
    )))
    calibration <- computeCensoringModelCalibration(clones, models)
    expect_equal(sort(unique(calibration$group)), 1:10)
    hazard <- calibration[calibration$model == "hazard" & calibration$grouping == "prediction", ]
    maxAbsZ[[covariateTimeInteraction]] <- max(abs(zScores(hazard)))
    binary <- calibration[calibration$model == "binary" & calibration$grouping == "prediction", ]
    expect_lt(max(abs(zScores(binary))), 4)
  }
  expect_gt(maxAbsZ$none, 6)
  expect_lt(maxAbsZ$day0, 4)
})

test_that("plotCensoringModelCalibration returns a ggplot and saves it", {
  skip_if_not_installed("ggplot2")
  calibration <- data.frame(
    strategyId = c(1, 1, 2, 2, 2, 2),
    model = c("binary", "binary", "hazard", "hazard", "binary", "binary"),
    group = c(1, 2, 1, 2, 1, 2),
    persons = 100,
    predicted = c(0.2, 0.6, 0.3, 0.7, 0.2, 0.6),
    observed = c(0.21, 0.59, 0.2, 0.62, 0.21, 0.59)
  )
  attr(calibration, "strategyLabels") <- c("Initiate", "Defer")
  plot <- plotCensoringModelCalibration(calibration)
  expect_s3_class(plot, "ggplot")
  expect_equal(levels(plot$data$strategy), c("Initiate", "Defer"))
  plot <- plotCensoringModelCalibration(calibration, labels = c("A", "B"))
  expect_equal(levels(plot$data$strategy), c("A", "B"))
  fileName <- tempfile(fileext = ".png")
  on.exit(unlink(fileName))
  plotCensoringModelCalibration(calibration, fileName = fileName)
  expect_true(file.exists(fileName))
  expect_error(plotCensoringModelCalibration(calibration[, -2]), "model")
  expect_error(plotCensoringModelCalibration(calibration, labels = "A"), "labels")
})

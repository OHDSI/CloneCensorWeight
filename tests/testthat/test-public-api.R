test_that("data can be saved and loaded", {
  data <- createSampleCloneCensorWeightData(n = 5, seed = 8)
  file <- tempfile()
  expect_invisible(saveCloneCensorWeightData(data, file))
  loaded <- loadCloneCensorWeightData(file)
  expect_s4_class(loaded, "CloneCensorWeightData")
  expect_equal(nrow(dplyr::collect(loaded$cohorts)), 5)
})

test_that("zero-outcome samples retain valid schemas", {
  data <- createSampleCloneCensorWeightData(n = 1, seed = 99)
  expect_true(all(c("subjectId", "outcomeDate", "outcomeId") %in%
                    names(dplyr::collect(data$outcomes))))
})

test_that("logistic models and diagnostics run", {
  data <- createSampleCloneCensorWeightData(n = 40, seed = 11)
  arms <- list(list(arm = "initiate", strategy = "initiate"),
               list(arm = "never", strategy = "never"))
  cloned <- createClonesAndCensor(data, arms, 30, 0)
  weighted <- createCensoringWeights(discretizePersonTime(cloned, 30), cloned)
  logisticModel <- fitOutcomeModel(weighted, "logistic")
  expect_s3_class(logisticModel, "CloneCensorWeightOutcomeModel")
  expect_s3_class(logisticModel$model, "cyclopsFit")
  expect_true(is.data.frame(computeWeightDistribution(weighted)))
  expect_true(is.data.frame(computeCovariateBalance(cloned, weighted)))
  expect_s3_class(plotWeightDistribution(weighted), "ggplot")
  expect_s3_class(plotKaplanMeier(weighted), "ggplot")
  expect_s3_class(plotCovariateBalance(computeCovariateBalance(cloned, weighted)), "ggplot")
  expect_true(is.data.frame(computeCovariateBalance(cloned)))
  expect_true(nrow(discretizePersonTime(cloned, cutPoints = c(0, 30, 60, 120))) > 0)
  analysis <- createCloneCensorWeightAnalysis("api", arms, 30)
  expect_equal(analysis$toList()[["analysisId"]], "api")
  outputFolder <- tempfile("ccw-api-")
  results <- runCloneCensorWeightAnalyses(data, list(analysis), outputFolder)
  expect_true("api" %in% names(results))
  expect_true(file.exists(file.path(outputFolder, "api", "results.csv")))
})

test_that("invalid arguments fail early", {
  data <- createSampleCloneCensorWeightData(n = 5)
  expect_error(createClonesAndCensor(data, list(list(arm = "one")), 30, 0))
  expect_error(discretizePersonTime(data, intervalLength = 0))
  expect_error(createCloneCensorWeightAnalysis("bad", list(), 30))
})
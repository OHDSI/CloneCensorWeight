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

# Shared arms for this test block
.testClonesArgs <- function() {
  createCreateClonesAndCensorArgs(
    arms = list(
      createArmDefinition("initiate", "initiateWithin", gracePeriod = 30),
      createArmDefinition("never",    "deferBeyond",    gracePeriod = 30)
    ),
    intervalLength = 30
  )
}

test_that("logistic models and diagnostics run", {
  data   <- createSampleCloneCensorWeightData(n = 40, seed = 11)
  cloned <- createClonesAndCensor(data, .testClonesArgs())
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
})

test_that("invalid arguments fail early", {
  data <- createSampleCloneCensorWeightData(n = 5)
  # createClonesAndCensor now requires CreateClonesAndCensorArgs R6 object
  expect_error(createClonesAndCensor(data, list(arm = "one")))
  expect_error(discretizePersonTime(data, intervalLength = 0))
  expect_error(createArmDefinition(armId = "a", strategyType = "invalid"))
})
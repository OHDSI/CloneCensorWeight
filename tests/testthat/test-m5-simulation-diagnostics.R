library(testthat)

test_that("computeMdrr equals CohortMethod::computeMdrr and the Schoenfeld formula on the clones", {
  data <- simulateCloneCensorWeightData(n = 200, graceEnd = 30, seed = 1)
  on.exit(close(data))
  for (riskStart in 0:1) {
    reference <- m5ReferenceClones(data, graceEnd = 30, riskStart = riskStart)
    ccwClones <- .createCcwClones(reference$persons, reference$clones, reference$metaData)
    clones <- reference$clones
    population <- data.frame(
      personSeqId = clones$personSeqId,
      treatment = as.integer(clones$strategyId == 1),
      outcomeCount = clones$outcomeCount,
      timeAtRisk = clones$survivalTime
    )
    for (twoSided in c(TRUE, FALSE)) {
      for (modelType in c("cox", "logistic")) {
        mdrr <- computeMdrr(ccwClones, alpha = 0.1, power = 0.9, twoSided = twoSided, modelType = modelType)
        expected <- CohortMethod::computeMdrr(population, alpha = 0.1, power = 0.9, twoSided = twoSided, modelType = modelType)
        expect_equal(mdrr, expected)
      }
    }

    mdrr <- computeMdrr(ccwClones)
    events <- sum(clones$outcomeCount)
    pTarget <- mean(clones$strategyId == 1)
    expect_gt(events, 0)
    expect_equal(mdrr$mdrr, exp(sqrt((stats::qnorm(0.8) + stats::qnorm(0.975))^2 / (events * pTarget * (1 - pTarget)))), tolerance = 1e-12)
    expect_equal(mdrr$se, 1 / sqrt(events * pTarget * (1 - pTarget)), tolerance = 1e-12)
    expect_equal(mdrr$targetPersons, length(unique(clones$personSeqId[clones$strategyId == 1])))
    expect_equal(mdrr$comparatorExposures, sum(clones$strategyId == 2))
    expect_equal(mdrr$targetDays, sum(clones$survivalTime[clones$strategyId == 1]))
    expect_equal(mdrr$totalOutcomes, events)
    close(ccwClones)
  }
})

test_that("computeMdrr accepts weighted clones and ignores the weights", {
  data <- simulateCloneCensorWeightData(n = 200, graceEnd = 30, seed = 2)
  on.exit(close(data))
  fixture <- m5WeightedClones(data, graceEnd = 30, riskStart = 0)
  on.exit(close(fixture$ccwClones), add = TRUE)
  reference <- m5ReferenceClones(data, graceEnd = 30, riskStart = 0)
  unweighted <- .createCcwClones(reference$persons, reference$clones, reference$metaData)
  on.exit(close(unweighted), add = TRUE)
  expect_equal(computeMdrr(fixture$ccwClones), computeMdrr(unweighted))
  expect_error(computeMdrr(data), "CcwClones")
})

test_that("the attrition of clones works with CohortMethod's attrition functions", {
  data <- simulateCloneCensorWeightData(n = 200, graceEnd = 30, seed = 1)
  on.exit(close(data))
  reference <- m5ReferenceClones(data, graceEnd = 30, riskStart = 1)
  ccwClones <- .createCcwClones(reference$persons, reference$clones, reference$metaData)
  on.exit(close(ccwClones), add = TRUE)
  expect_equal(CohortMethod::getAttritionTable(ccwClones), reference$metaData$attrition)
  expect_equal(CohortMethod::getAttritionTable(data), attr(data, "metaData")$attrition)
  skip_if_not_installed("ggplot2")
  diagram <- CohortMethod::drawAttritionDiagram(ccwClones, targetLabel = "Initiate", comparatorLabel = "Defer")
  expect_s3_class(diagram, "ggplot")
})

test_that("getAttritionTable and drawAttritionDiagram are re-exported from CohortMethod", {
  # The re-exports take effect once NAMESPACE is regenerated with devtools::document().
  skip_if_not("getAttritionTable" %in% getNamespaceExports("CloneCensorWeight"), "NAMESPACE not yet regenerated")
  expect_identical(getAttritionTable, CohortMethod::getAttritionTable)
  expect_identical(drawAttritionDiagram, CohortMethod::drawAttritionDiagram)
})

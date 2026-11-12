library(testthat)

# Compares CohortMethod's balance columns with the hand computation, matched by covariate ID.
expectBalanceEqual <- function(balance, hand, prefix) {
  balance <- balance[match(hand$covariateId, balance$covariateId), ]
  for (column in c("MeanTarget", "MeanComparator", "SdTarget", "SdComparator", "StdDiff")) {
    expect_equal(balance[[paste0(prefix, column)]], hand[[paste0(tolower(substr(column, 1, 1)), substring(column, 2))]],
      tolerance = 1e-8, label = paste0(prefix, column)
    )
  }
}

test_that("computeCovariateBalance equals a hand computation on compliant clones", {
  data <- simulateCloneCensorWeightData(n = 200, graceEnd = 30, seed = 1)
  on.exit(close(data))
  for (riskStart in 0:1) {
    fixture <- m5WeightedClones(data, graceEnd = 30, riskStart = riskStart)
    balance <- computeCovariateBalance(fixture$ccwClones, data, CohortMethod::createComputeCovariateBalanceArgs(maxCohortSize = 0))
    close(fixture$ccwClones)

    # Clones still following their strategy on day 31, with the weight of days after the grace period:
    clones <- fixture$clones
    compliant <- clones$endDay >= 31
    expect_true(any(clones$endDay < 31 & clones$endReason == "deviation"))
    x <- fixture$x[clones$rowId[compliant], , drop = FALSE]
    treatment <- as.integer(clones$strategyId[compliant] == 1)
    expect_equal(nrow(balance), ncol(x))

    # Before: the crude comparison of compliant clones.
    before <- m5HandBalance(x, treatment, weighted = FALSE)
    expectBalanceEqual(balance, before, "beforeMatching")
    binary <- balance$covariateId %in% c(1001, 3001, 5001)
    expect_equal(
      balance$beforeMatchingSumTarget[binary],
      colSums(x[treatment == 1, as.character(balance$covariateId[binary])]),
      ignore_attr = TRUE
    )

    # After: weighted by W(31).
    after <- m5HandBalance(x, treatment, fixture$after[compliant])
    expectBalanceEqual(balance, after, "afterMatching")
    expect_false(isTRUE(all.equal(before$stdDiff, after$stdDiff)))
  }
})

test_that("computeCovariateBalance reads the weight of the segment that contains day tau + 1", {
  data <- simulateCloneCensorWeightData(n = 200, graceEnd = 30, seed = 1)
  on.exit(close(data))
  x <- m5CovariateMatrix(data)
  for (riskStart in 0:1) {
    # riskWindowEnd = 31 makes every clone still under follow-up on day 31 end on that day.
    for (riskWindowEnd in c(365, 31)) {
      reference <- m5ReferenceClones(data, graceEnd = 30, riskStart = riskStart, riskWindowEnd = riskWindowEnd)
      clones <- reference$clones
      # tau = 30. Day 31 is time t31 = 31 - riskStart + 1; cut every clone at t31 - 1 and t31, so that
      # (t31 - 1, t31] is day 31 alone, and give each of the three pieces a different weight:
      t31 <- 31 - riskStart + 1
      weights <- cbind(1 + clones$cloneRowId %% 2, 1 + clones$cloneRowId %% 5, 1 + clones$cloneRowId %% 3)
      weighted <- m5CutSegments(clones, c(t31 - 1, t31), weights)
      ccwClones <- .createCcwClones(reference$persons, clones, reference$metaData)
      ccwClones <- .addSegments(ccwClones, weighted$segments, weighted$personWeights, list(stabilized = FALSE, truncation = "none"))
      balance <- computeCovariateBalance(ccwClones, data)
      close(ccwClones)

      compliant <- clones$endDay >= 31
      if (riskWindowEnd == 31) {
        expect_true(all(clones$endDay[compliant] == 31))
      }
      xCompliant <- x[clones$rowId[compliant], , drop = FALSE]
      treatment <- as.integer(clones$strategyId[compliant] == 1)
      expectBalanceEqual(balance, m5HandBalance(xCompliant, treatment, weighted = FALSE), "beforeMatching")
      expectBalanceEqual(balance, m5HandBalance(xCompliant, treatment, weights[compliant, 2]), "afterMatching")
      # The weights of the neighbouring pieces give a different answer:
      for (k in c(1, 3)) {
        other <- m5HandBalance(xCompliant, treatment, weights[compliant, k])
        expect_false(isTRUE(all.equal(other$stdDiff, m5HandBalance(xCompliant, treatment, weights[compliant, 2])$stdDiff)))
      }
    }
  }
})

test_that("computeCovariateBalance uses the largest graceEnd of the two strategies", {
  data <- simulateCloneCensorWeightData(n = 200, graceEnd = 30, seed = 2)
  on.exit(close(data))
  fixture <- m5WeightedClones(data, graceEnd = 30, riskStart = 0)
  on.exit(close(fixture$ccwClones), add = TRUE)
  # Make the comparator neverInitiate and the target initiateWithin(20): tau = 20, so W(21) is read
  # from each clone's segment containing day 21, which is the first segment when the clone is split at
  # day 30.
  metaData <- attr(fixture$ccwClones, "metaData")
  metaData$strategies[[1]]$graceEnd <- 20
  metaData$strategies[[2]] <- createTreatmentStrategy(2, "Never initiate", "neverInitiate")$toList()
  attr(fixture$ccwClones, "metaData") <- metaData
  balance <- computeCovariateBalance(fixture$ccwClones, data)

  clones <- fixture$clones
  segments <- dplyr::collect(fixture$ccwClones$segments)
  firstSegment <- segments[!duplicated(segments$cloneRowId), ]
  compliant <- clones$endDay >= 21
  x <- fixture$x[clones$rowId[compliant], , drop = FALSE]
  treatment <- as.integer(clones$strategyId[compliant] == 1)
  weight <- firstSegment$weight[match(clones$cloneRowId[compliant], firstSegment$cloneRowId)]
  expectBalanceEqual(balance, m5HandBalance(x, treatment, weighted = FALSE), "beforeMatching")
  expectBalanceEqual(balance, m5HandBalance(x, treatment, weight), "afterMatching")
})

test_that("weighting with the true DGP weights balances the compliant clones", {
  data <- simulateCloneCensorWeightData(n = 5000, graceEnd = 30, seed = 1)
  on.exit(close(data))
  fixture <- m5WeightedClones(data, graceEnd = 30, riskStart = 1)
  on.exit(close(fixture$ccwClones), add = TRUE)
  balance <- computeCovariateBalance(fixture$ccwClones, data)
  # Initiators by day 30 have higher covariate values than deferrers:
  expect_true(all(balance$beforeMatchingStdDiff > 0.15))
  # After weighting, every covariate is within CohortMethod's default threshold, and within sampling
  # error of zero:
  expect_true(all(abs(balance$afterMatchingStdDiff) < 0.1))
  expect_true(all(abs(balance$afterMatchingStdDiff) < 3 * sqrt(balance$afterMatchingSdmVariance)))
})

test_that("computeCovariateBalance checks its input", {
  data <- simulateCloneCensorWeightData(n = 100, seed = 1)
  on.exit(close(data))
  fixture <- m5WeightedClones(data, graceEnd = 30, riskStart = 0)
  on.exit(close(fixture$ccwClones), add = TRUE)
  ccwClones <- fixture$ccwClones

  expect_error(computeCovariateBalance(data, data), "CcwClones")
  expect_error(computeCovariateBalance(ccwClones, ccwClones), "CloneCensorWeightData")
  expect_error(computeCovariateBalance(ccwClones, data, CohortMethod::createCreateStudyPopulationArgs()), "ComputeCovariateBalanceArgs")
  unweighted <- ccwClones
  attr(unweighted, "metaData")$weights <- NULL
  expect_error(computeCovariateBalance(unweighted, data), "metaData\\$weights")

  noGrace <- ccwClones
  never <- createTreatmentStrategy(1, "Never initiate", "neverInitiate")$toList()
  attr(noGrace, "metaData")$strategies <- list(never, never)
  expect_error(computeCovariateBalance(noGrace, data), "graceEnd")

  lateGrace <- ccwClones
  attr(lateGrace, "metaData")$strategies[[1]]$graceEnd <- 400
  expect_error(computeCovariateBalance(lateGrace, data), "No clones of strategy 1 are under follow-up on day 401")
})

library(testthat)

test_that(".createCcwClones accepts valid tables", {
  fixture <- createFoundationClones()
  ccwClones <- .createCcwClones(fixture$persons, fixture$clones, fixture$metaData)
  on.exit(close(ccwClones))

  expect_true(is(ccwClones, "CcwClones"))
  expect_true(is(ccwClones, "Andromeda"))
  expect_true(checkmate::testTRUE(.checkCcwClones(ccwClones)))
  expect_equal(as.data.frame(dplyr::collect(ccwClones$persons)), fixture$persons)
  expect_equal(as.data.frame(dplyr::collect(ccwClones$clones)), fixture$clones)
  expect_equal(attr(ccwClones, "metaData"), fixture$metaData)
})

test_that(".createCcwClones accepts Andromeda tables", {
  fixture <- createFoundationClones()
  source <- Andromeda::andromeda(persons = fixture$persons, clones = fixture$clones)
  on.exit(close(source))
  ccwClones <- .createCcwClones(source$persons, source$clones, fixture$metaData)
  on.exit(close(ccwClones), add = TRUE)
  expect_equal(as.data.frame(dplyr::collect(ccwClones$clones)), fixture$clones)

  clones <- fixture$clones
  clones$endDay[1] <- NA
  source$badClones <- clones
  expect_error(.createCcwClones(source$persons, source$badClones, fixture$metaData), "'endDay' of table 'clones' must not contain NA")
})

test_that(".createCcwClones rejects invalid tables", {
  fixture <- createFoundationClones()
  expectRejected <- function(persons = fixture$persons, clones = fixture$clones, metaData = fixture$metaData, regexp) {
    expect_error(.createCcwClones(persons, clones, metaData), regexp)
  }

  persons <- fixture$persons
  persons$nominalRiskEnd <- NULL
  expectRejected(persons = persons, regexp = "missing column\\(s\\): nominalRiskEnd")

  persons <- fixture$persons
  persons$riskEnd <- as.numeric(persons$riskEnd)
  expectRejected(persons = persons, regexp = "'riskEnd' of table 'persons' must be of type integer")

  persons <- fixture$persons
  persons$rowId[2] <- 1L
  expectRejected(persons = persons, regexp = "duplicateRowId")

  clones <- fixture$clones
  clones$endReason <- NULL
  expectRejected(clones = clones, regexp = "missing column\\(s\\): endReason")

  clones <- fixture$clones
  clones$strategyId <- as.numeric(clones$strategyId)
  expectRejected(clones = clones, regexp = "'strategyId' of table 'clones' must be of type integer")

  clones <- fixture$clones
  clones$endReason <- factor(clones$endReason)
  expectRejected(clones = clones, regexp = "'endReason' of table 'clones' must be of type character")

  clones <- fixture$clones
  clones$survivalTime[1] <- NA
  expectRejected(clones = clones, regexp = "must not contain NA")

  clones <- fixture$clones
  clones$cloneRowId[2] <- 1L
  expectRejected(clones = clones, regexp = "duplicateCloneRowId")

  clones <- fixture$clones
  clones$cloneRowId <- clones$cloneRowId + 1L
  expectRejected(clones = clones, regexp = "cloneRowIdNotOneToN")

  clones <- fixture$clones
  clones$strategyId[1] <- 3L
  expectRejected(clones = clones, regexp = "strategyIdNotOneOrTwo")

  clones <- fixture$clones
  clones$survivalTime[1] <- clones$survivalTime[1] + 1L
  expectRejected(clones = clones, regexp = "badSurvivalTime")

  clones <- fixture$clones
  clones$endReason[1] <- "censored"
  expectRejected(clones = clones, regexp = "badEndReason")

  clones <- fixture$clones
  clones$outcomeCount[1] <- 1L
  expectRejected(clones = clones, regexp = "badOutcomeCount")

  metaData <- fixture$metaData
  metaData$sameDayOrdering <- NULL
  expectRejected(metaData = metaData, regexp = "sameDayOrdering")

  metaData <- fixture$metaData
  metaData$strategies <- metaData$strategies[1]
  expectRejected(metaData = metaData, regexp = "strategies")

  expectRejected(persons = list(), regexp = "persons")
})

test_that(".addSegments accepts valid tables and rejects invalid ones", {
  fixture <- createFoundationClones()
  ccwClones <- .createCcwClones(fixture$persons, fixture$clones, fixture$metaData)
  on.exit(close(ccwClones))
  weighted <- createFoundationSegments(fixture$clones)

  expect_error(.addSegments(ccwClones, weighted$segments[, -1], weighted$personWeights, list()), "missing column\\(s\\): segmentId")
  segments <- weighted$segments
  segments$outcome <- as.numeric(segments$outcome)
  expect_error(.addSegments(ccwClones, segments, weighted$personWeights, list()), "'outcome' of table 'segments' must be of type integer")
  segments <- weighted$segments
  segments$weight[1] <- NA
  expect_error(.addSegments(ccwClones, segments, weighted$personWeights, list()), "must not contain NA")
  segments <- weighted$segments
  segments$stopTime[1] <- 0
  expect_error(.addSegments(ccwClones, segments, weighted$personWeights, list()), "stopTimeNotAfterStartTime")
  segments <- weighted$segments
  segments$weight[1] <- -1
  expect_error(.addSegments(ccwClones, segments, weighted$personWeights, list()), "negativeWeight")
  segments <- weighted$segments
  segments$competing[4] <- 1L
  expect_error(.addSegments(ccwClones, segments, weighted$personWeights, list()), "badEventIndicator")
  segments <- weighted$segments
  segments$segmentId[2] <- 1L
  expect_error(.addSegments(ccwClones, segments, weighted$personWeights, list()), "duplicateSegmentId")
  expect_false(any(c("segments", "personWeights") %in% names(ccwClones)))
  personWeights <- weighted$personWeights
  personWeights$maxWeight <- NULL
  expect_error(.addSegments(ccwClones, weighted$segments, personWeights, list()), "missing column\\(s\\): maxWeight")
  expect_error(.addSegments(fixture$clones, weighted$segments, weighted$personWeights, list()), "CcwClones")
  expect_error(.addSegments(ccwClones, weighted$segments, weighted$personWeights, "stabilized"), "weightsMetaData")

  weightsMetaData <- list(stabilized = TRUE, truncation = createWeightTruncation()$toList(), cap = c(`1` = 5, `2` = 4))
  result <- .addSegments(ccwClones, weighted$segments, weighted$personWeights, weightsMetaData)
  expect_true(is(result, "CcwClones"))
  expect_true(checkmate::testTRUE(.checkCcwClones(result, weighted = TRUE)))
  # The caller's original object shares the tables but not metaData$weights, so it is rejected:
  expect_true(all(c("segments", "personWeights") %in% names(ccwClones)))
  expect_null(attr(ccwClones, "metaData")$weights)
  expect_error(.assertCcwClones(ccwClones, weighted = TRUE), "metaData\\$weights")
  expect_equal(as.data.frame(dplyr::collect(result$segments)), weighted$segments)
  expect_equal(as.data.frame(dplyr::collect(result$personWeights)), weighted$personWeights)
  expect_equal(attr(result, "metaData")$weights, weightsMetaData)
  expect_equal(attr(result, "metaData")$strategies, fixture$metaData$strategies)
})

test_that(".assertCcwClones", {
  fixture <- createFoundationClones()
  ccwClones <- .createCcwClones(fixture$persons, fixture$clones, fixture$metaData)

  expect_silent(.assertCcwClones(ccwClones))
  expect_error(.assertCcwClones(ccwClones, weighted = TRUE), "segments, personWeights")
  expect_error(.assertCcwClones(fixture$clones), "CcwClones")
  collection <- checkmate::makeAssertCollection()
  .assertCcwClones(fixture$clones, add = collection)
  expect_false(collection$isEmpty())
  close(ccwClones)
  expect_error(.assertCcwClones(ccwClones), "closed")
})

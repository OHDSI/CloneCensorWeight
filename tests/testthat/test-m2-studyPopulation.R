library(testthat)

test_that("createStudyPopulation equals CohortMethod's and records the nominal risk end", {
  data <- createM2FixtureData()
  on.exit(close(data))
  for (riskStart in 0:1) {
    args <- CohortMethod::createCreateStudyPopulationArgs(
      riskWindowStart = riskStart,
      riskWindowEnd = 30,
      endAnchor = "cohort start"
    )
    population <- suppressMessages(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = args))
    cmPopulation <- suppressMessages(CohortMethod::createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = args))

    metaData <- attr(population, "metaData")
    expect_equal(metaData$outcomeId, 3)
    expect_equal(metaData$nominalRiskEnd, list(riskWindowEnd = 30, endAnchor = "cohort start", maxDaysAtRisk = 99999))
    metaData$outcomeId <- NULL
    metaData$nominalRiskEnd <- NULL
    attr(population, "metaData") <- metaData
    expect_equal(population, cmPopulation)
    expect_true(checkmate::testTRUE(.checkStudyPopulation(population)))

    # Prior outcomes are removed: person 10 (day -10) and, when riskStart = 1, person 11 (day 0).
    if (riskStart == 0) {
      expect_equal(population$rowId, c(1:9, 11:12))
      expect_equal(population$daysToEvent, c(NA, 25, NA, 7, 3, NA, NA, 4, NA, 0, NA))
      expect_equal(population$riskEnd, c(30, 30, 30, 7, 30, 4, 5, 30, 7, 30, 3))
      expect_equal(metaData$attrition$targetPersons, c(12, 11, 11))
    } else {
      expect_equal(population$rowId, c(1:9, 12))
      expect_equal(population$daysToEvent, c(NA, 25, NA, 7, 3, NA, NA, 4, NA, NA))
      expect_equal(population$riskEnd, c(30, 30, 30, 7, 30, 4, 5, 30, 7, 3))
      expect_equal(metaData$attrition$targetPersons, c(12, 10, 10))
    }
    expect_true(all(population$riskStart == riskStart))
  }
})

test_that("createStudyPopulation starts from a given population", {
  data <- createM2FixtureData()
  on.exit(close(data))
  full <- suppressMessages(createStudyPopulation(data, outcomeId = 3))
  subset <- full[full$rowId %in% c(1, 2, 3, 10), ]
  population <- suppressMessages(createStudyPopulation(
    data,
    population = subset,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 10, endAnchor = "cohort start")
  ))
  expect_equal(population$rowId, 1:3)
  expect_equal(population$riskEnd, c(10, 10, 10))
  expect_equal(attr(population, "metaData")$nominalRiskEnd$riskWindowEnd, 10)
  expect_equal(attr(population, "metaData")$nominalRiskEnd$endAnchor, "cohort start")
})

test_that("createStudyPopulation rejects invalid input", {
  data <- createM2FixtureData()
  on.exit(close(data))
  expect_error(
    createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(startAnchor = "cohort end")),
    "startAnchor = 'cohort start'"
  )
  expect_error(
    createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowStart = -1)),
    "riskWindowStart must be >= 0"
  )
  expect_error(createStudyPopulation(data, outcomeId = 3.5), "integer")
  expect_error(createStudyPopulation(data, outcomeId = c(3, 4)), "outcomeId")
  expect_error(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = list()), "createStudyPopulationArgs")
  expect_error(createStudyPopulation(list(), outcomeId = 3), "CloneCensorWeightData")
})

library(testthat)

# Hand-written expected tables for the 12-person fixture of helper-m2-clones.R (riskWindowEnd = 30,
# anchored at cohort start). Person 10 has a prior outcome and is never in the study population;
# person 11 (outcome on day 0) is a prior outcome when riskStart = 1.

expectedM2Persons <- list(
  riskStart0 = dplyr::tribble(
    ~rowId, ~riskEnd, ~nominalRiskEnd, ~outcomeDay, ~deathDay, ~treatmentDay, ~comparatorTreatmentDay,
    1,      30,       30,              NA,          NA,        0,             NA,
    2,      30,       30,              25,          NA,        7,             NA,
    3,      30,       30,              NA,          NA,        12,            NA,
    4,      7,        30,              7,           NA,        NA,            NA,
    5,      30,       30,              3,           NA,        5,             NA,
    6,      4,        30,              NA,          4,         NA,            NA,
    7,      5,        30,              NA,          NA,        NA,            1,
    8,      30,       30,              4,           NA,        4,             NA,
    9,      7,        30,              NA,          NA,        NA,            2,
    11,     30,       30,              0,           NA,        NA,            NA,
    12,     3,        30,              NA,          3,         3,             NA
  ),
  riskStart1 = dplyr::tribble(
    ~rowId, ~riskEnd, ~nominalRiskEnd, ~outcomeDay, ~deathDay, ~treatmentDay, ~comparatorTreatmentDay,
    1,      30,       30,              NA,          NA,        0,             NA,
    2,      30,       30,              25,          NA,        7,             NA,
    3,      30,       30,              NA,          NA,        12,            NA,
    4,      7,        30,              7,           NA,        NA,            NA,
    5,      30,       30,              3,           NA,        5,             NA,
    6,      4,        30,              NA,          4,         NA,            NA,
    7,      5,        30,              NA,          NA,        NA,            1,
    8,      30,       30,              4,           NA,        4,             NA,
    9,      7,        30,              NA,          NA,        NA,            2,
    12,     3,        30,              NA,          3,         3,             NA
  )
)

# Expected clones per strategy, with sameDayOrdering = "eventsFirst" and competingEvents = "death".
# `deviationFirst` lists the rows that differ under sameDayOrdering = "deviationFirst".
expectedM2Clones <- list(
  initiateWithin7 = list(
    riskStart0 = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      1,      NA,            30,      "riskEnd",
      2,      NA,            25,      "outcome",
      3,      7,             7,       "deviation",
      4,      7,             7,       "outcome",
      5,      NA,            3,       "outcome",
      6,      NA,            4,       "death",
      7,      NA,            5,       "observationEnd",
      8,      NA,            4,       "outcome",
      9,      7,             7,       "deviation",
      11,     7,             0,       "outcome",
      12,     NA,            3,       "death"
    ),
    riskStart1 = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      1,      NA,            30,      "riskEnd",
      2,      NA,            25,      "outcome",
      3,      7,             7,       "deviation",
      4,      7,             7,       "outcome",
      5,      NA,            3,       "outcome",
      6,      NA,            4,       "death",
      7,      NA,            5,       "observationEnd",
      8,      NA,            4,       "outcome",
      9,      7,             7,       "deviation",
      12,     NA,            3,       "death"
    ),
    deviationFirst = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      4,      7,             7,       "deviation"
    )
  ),
  notInitiateWithin7 = list(
    riskStart0 = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      1,      0,             0,       "deviation",
      2,      7,             7,       "deviation",
      3,      NA,            30,      "riskEnd",
      4,      NA,            7,       "outcome",
      5,      5,             3,       "outcome",
      6,      NA,            4,       "death",
      7,      NA,            5,       "observationEnd",
      8,      4,             4,       "outcome",
      9,      NA,            7,       "observationEnd",
      11,     NA,            0,       "outcome",
      12,     3,             3,       "death"
    ),
    riskStart1 = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      2,      7,             7,       "deviation",
      3,      NA,            30,      "riskEnd",
      4,      NA,            7,       "outcome",
      5,      5,             3,       "outcome",
      6,      NA,            4,       "death",
      7,      NA,            5,       "observationEnd",
      8,      4,             4,       "outcome",
      9,      NA,            7,       "observationEnd",
      12,     3,             3,       "death"
    ),
    deviationFirst = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      8,      4,             4,       "deviation",
      12,     3,             3,       "deviation"
    )
  ),
  neverInitiate = list(
    riskStart0 = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      1,      0,             0,       "deviation",
      2,      7,             7,       "deviation",
      3,      12,            12,      "deviation",
      4,      NA,            7,       "outcome",
      5,      5,             3,       "outcome",
      6,      NA,            4,       "death",
      7,      NA,            5,       "observationEnd",
      8,      4,             4,       "outcome",
      9,      NA,            7,       "observationEnd",
      11,     NA,            0,       "outcome",
      12,     3,             3,       "death"
    ),
    riskStart1 = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      2,      7,             7,       "deviation",
      3,      12,            12,      "deviation",
      4,      NA,            7,       "outcome",
      5,      5,             3,       "outcome",
      6,      NA,            4,       "death",
      7,      NA,            5,       "observationEnd",
      8,      4,             4,       "outcome",
      9,      NA,            7,       "observationEnd",
      12,     3,             3,       "death"
    ),
    deviationFirst = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      8,      4,             4,       "deviation",
      12,     3,             3,       "deviation"
    )
  ),
  initiateWithin0 = list(
    riskStart0 = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      1,      NA,            30,      "riskEnd",
      2,      0,             0,       "deviation",
      3,      0,             0,       "deviation",
      4,      0,             0,       "deviation",
      5,      0,             0,       "deviation",
      6,      0,             0,       "deviation",
      7,      0,             0,       "deviation",
      8,      0,             0,       "deviation",
      9,      0,             0,       "deviation",
      11,     0,             0,       "outcome",
      12,     0,             0,       "deviation"
    ),
    riskStart1 = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      1,      NA,            30,      "riskEnd"
    ),
    deviationFirst = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      11,     0,             0,       "deviation"
    )
  ),
  notInitiateWithin0 = list(
    riskStart0 = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      1,      0,             0,       "deviation",
      2,      NA,            25,      "outcome",
      3,      NA,            30,      "riskEnd",
      4,      NA,            7,       "outcome",
      5,      NA,            3,       "outcome",
      6,      NA,            4,       "death",
      7,      NA,            5,       "observationEnd",
      8,      NA,            4,       "outcome",
      9,      NA,            7,       "observationEnd",
      11,     NA,            0,       "outcome",
      12,     NA,            3,       "death"
    ),
    riskStart1 = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      2,      NA,            25,      "outcome",
      3,      NA,            30,      "riskEnd",
      4,      NA,            7,       "outcome",
      5,      NA,            3,       "outcome",
      6,      NA,            4,       "death",
      7,      NA,            5,       "observationEnd",
      8,      NA,            4,       "outcome",
      9,      NA,            7,       "observationEnd",
      12,     NA,            3,       "death"
    ),
    deviationFirst = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason
    )
  )
)

m2Strategies <- list(
  initiateWithin7 = function(strategyId) m2Strategy(strategyId, "initiateWithin", 7),
  notInitiateWithin7 = function(strategyId) m2Strategy(strategyId, "notInitiateWithin", 7),
  neverInitiate = function(strategyId) m2Strategy(strategyId, "neverInitiate"),
  initiateWithin0 = function(strategyId) m2Strategy(strategyId, "initiateWithin", 0),
  notInitiateWithin0 = function(strategyId) m2Strategy(strategyId, "notInitiateWithin", 0)
)

expectedM2Strategy <- function(name, strategyId, riskStart, sameDayOrdering) {
  expected <- expectedM2Clones[[name]][[paste0("riskStart", riskStart)]]
  if (sameDayOrdering == "deviationFirst") {
    replacements <- expectedM2Clones[[name]]$deviationFirst
    replacements <- replacements[replacements$rowId %in% expected$rowId, ]
    expected <- replaceM2Rows(expected, replacements)
  }
  expected$strategyId <- rep(strategyId, nrow(expected))
  return(expected)
}

toM2Persons <- function(expected, riskStart) {
  expected <- as.data.frame(expected)
  result <- data.frame(
    rowId = as.integer(expected$rowId),
    personSeqId = as.integer(100 + expected$rowId),
    riskStart = as.integer(riskStart)
  )
  for (column in setdiff(colnames(expected), "rowId")) {
    result[[column]] <- as.integer(expected[[column]])
  }
  return(result)
}

for (pair in list(
  c("initiateWithin7", "notInitiateWithin7"),
  c("initiateWithin7", "neverInitiate"),
  c("initiateWithin0", "notInitiateWithin0")
)) {
  for (riskStart in 0:1) {
    for (sameDayOrdering in c("eventsFirst", "deviationFirst")) {
      test_that(sprintf("Fixture clones for %s vs %s, riskStart = %d, %s", pair[1], pair[2], riskStart, sameDayOrdering), {
        data <- createM2FixtureData()
        on.exit(close(data))
        population <- createM2FixturePopulation(data, riskStart)
        result <- runM2Clones(
          population = population,
          data = data,
          strategies = list(m2Strategies[[pair[1]]](10), m2Strategies[[pair[2]]](20)),
          sameDayOrdering = sameDayOrdering
        )
        expected <- dplyr::bind_rows(
          expectedM2Strategy(pair[1], 1, riskStart, sameDayOrdering),
          expectedM2Strategy(pair[2], 2, riskStart, sameDayOrdering)
        )
        expect_equal(result$clones, completeM2ExpectedClones(expected, riskStart))
        expect_equal(result$persons, toM2Persons(expectedM2Persons[[paste0("riskStart", riskStart)]], riskStart))
      })
    }
  }
}

test_that("Fixture clones for a strategy on the comparator treatment role", {
  data <- createM2FixtureData()
  on.exit(close(data))
  population <- createM2FixturePopulation(data, 0)
  result <- runM2Clones(
    population = population,
    data = data,
    strategies = list(
      m2Strategy(1, "initiateWithin", 7),
      m2Strategy(2, "neverInitiate", treatment = "comparatorTreatment")
    )
  )
  neverComparator <- dplyr::tribble(
    ~rowId, ~deviationDay, ~endDay, ~endReason,
    1,      NA,            30,      "riskEnd",
    2,      NA,            25,      "outcome",
    3,      NA,            30,      "riskEnd",
    4,      NA,            7,       "outcome",
    5,      NA,            3,       "outcome",
    6,      NA,            4,       "death",
    7,      1,             1,       "deviation",
    8,      NA,            4,       "outcome",
    9,      2,             2,       "deviation",
    11,     NA,            0,       "outcome",
    12,     NA,            3,       "death"
  )
  neverComparator$strategyId <- 2
  expected <- dplyr::bind_rows(expectedM2Strategy("initiateWithin7", 1, 0, "eventsFirst"), neverComparator)
  expect_equal(result$clones, completeM2ExpectedClones(expected, 0))
})

test_that("Fixture clones with competingEvents = 'none'", {
  data <- createM2FixtureData()
  on.exit(close(data))
  for (riskStart in 0:1) {
    population <- createM2FixturePopulation(data, riskStart)
    result <- runM2Clones(
      population = population,
      data = data,
      strategies = list(m2Strategy(1, "initiateWithin", 7), m2Strategy(2, "notInitiateWithin", 7)),
      competingEvents = "none"
    )
    # Deaths (persons 6 and 12) are no longer events: follow-up runs to the end of observation, which
    # is the death day. Person 12 starts on day 3, the day observation ends, so the deviation comes
    # first in the notInitiateWithin arm.
    noneReplacements <- list(
      dplyr::tribble(
        ~rowId, ~deviationDay, ~endDay, ~endReason,
        6,      NA,            4,       "observationEnd",
        12,     NA,            3,       "observationEnd"
      ),
      dplyr::tribble(
        ~rowId, ~deviationDay, ~endDay, ~endReason,
        6,      NA,            4,       "observationEnd",
        12,     3,             3,       "deviation"
      )
    )
    expected <- dplyr::bind_rows(
      replaceM2Rows(expectedM2Clones$initiateWithin7[[paste0("riskStart", riskStart)]], noneReplacements[[1]]) |>
        dplyr::mutate(strategyId = 1),
      replaceM2Rows(expectedM2Clones$notInitiateWithin7[[paste0("riskStart", riskStart)]], noneReplacements[[2]]) |>
        dplyr::mutate(strategyId = 2)
    )
    expect_equal(result$clones, completeM2ExpectedClones(expected, riskStart))
    expectedPersons <- expectedM2Persons[[paste0("riskStart", riskStart)]]
    expectedPersons$deathDay <- NA
    expect_equal(result$persons, toM2Persons(expectedPersons, riskStart))
  }
})

test_that("Nominal risk end follows endAnchor and maxDaysAtRisk", {
  tables <- createM2FixtureTables()
  tables$cohorts$daysToCohortEnd <- c(2L, 10L, 50L, 9L, 0L, 9L, 9L, 0L, 9L, 0L, 0L, 9L)
  data <- do.call(createCloneCensorWeightData, tables)
  on.exit(close(data))
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
      riskWindowEnd = 0,
      endAnchor = "cohort end",
      maxDaysAtRisk = 20
    )
  ))
  result <- runM2Clones(
    population = population,
    data = data,
    strategies = list(m2Strategy(1, "initiateWithin", 7), m2Strategy(2, "notInitiateWithin", 7))
  )
  expectedPersons <- data.frame(
    rowId = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 11L, 12L),
    riskEnd = c(2L, 10L, 20L, 7L, 0L, 4L, 5L, 0L, 7L, 0L, 3L),
    nominalRiskEnd = c(2L, 10L, 20L, 9L, 0L, 9L, 9L, 0L, 9L, 0L, 9L)
  )
  expect_equal(result$persons[, c("rowId", "riskEnd", "nominalRiskEnd")], expectedPersons)
  notInitiate <- result$clones[result$clones$strategyId == 2, c("rowId", "endDay", "endReason")]
  rownames(notInitiate) <- NULL
  expectedNotInitiate <- data.frame(
    rowId = c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 11L, 12L),
    endDay = c(0L, 7L, 20L, 7L, 0L, 4L, 5L, 0L, 7L, 0L, 3L),
    endReason = c(
      "deviation", "deviation", "riskEnd", "outcome", "riskEnd", "death", "observationEnd", "riskEnd",
      "observationEnd", "outcome", "death"
    )
  )
  expect_equal(notInitiate, expectedNotInitiate)
})

test_that("Clone attrition and metaData", {
  data <- createM2FixtureData()
  on.exit(close(data))
  population <- createM2FixturePopulation(data, 1)
  strategies <- list(m2Strategy(5, "initiateWithin", 0), m2Strategy(3, "notInitiateWithin", 0))
  result <- runM2Clones(population = population, data = data, strategies = strategies)
  metaData <- result$metaData
  expect_equal(metaData$strategies, lapply(strategies, function(strategy) strategy$toList()))
  expect_equal(metaData$outcomeId, 3)
  expect_equal(metaData$treatmentIds, c(treatment = 21, comparatorTreatment = 22))
  expect_equal(metaData$sameDayOrdering, "eventsFirst")
  expect_equal(metaData$competingEvents, "death")
  expect_equal(metaData$targetId, 1)

  attrition <- metaData$attrition
  populationAttrition <- attr(population, "metaData")$attrition
  expect_equal(attrition[seq_len(nrow(populationAttrition)), ], populationAttrition)
  cloneAttrition <- as.data.frame(attrition[nrow(populationAttrition) + 1:2, ])
  expect_equal(cloneAttrition, data.frame(
    description = c("Clones created", "Clones with follow-up"),
    targetPersons = c(10L, 1L),
    comparatorPersons = c(10L, 9L),
    targetExposures = c(10L, 1L),
    comparatorExposures = c(10L, 9L)
  ))
})

test_that("createClones stops for unsupported strategies and invalid input", {
  data <- createM2FixtureData()
  on.exit(close(data))
  population <- createM2FixturePopulation(data, 0)
  clonesArgs <- function(...) {
    createCreateClonesArgs(strategies = list(createTreatmentStrategy(1, "Initiate", "initiateWithin", ...), m2Strategy(2, "neverInitiate")))
  }
  expect_error(
    createClones(population, data, clonesArgs(graceStart = 2, graceEnd = 7)),
    "Not yet supported in this version: graceStart > 0"
  )
  expect_error(
    createClones(population, data, clonesArgs(graceEnd = 7, prohibitedTreatments = "comparatorTreatment")),
    "Not yet supported in this version: prohibitedTreatments"
  )

  tables <- createM2FixtureTables()
  tables$metaData$treatmentIds <- c(treatment = 21)
  dataWithoutComparator <- do.call(createCloneCensorWeightData, tables)
  on.exit(close(dataWithoutComparator), add = TRUE)
  comparatorArgs <- createCreateClonesArgs(strategies = list(
    m2Strategy(1, "initiateWithin", 7),
    m2Strategy(2, "neverInitiate", treatment = "comparatorTreatment")
  ))
  expect_error(
    createClones(population, dataWithoutComparator, comparatorArgs),
    "role 'comparatorTreatment', which has no cohort ID"
  )

  cmPopulation <- suppressMessages(CohortMethod::createStudyPopulation(data, outcomeId = 3))
  expect_error(
    createClones(cmPopulation, data, clonesArgs(graceEnd = 7)),
    "createStudyPopulation\\(\\) from CloneCensorWeight"
  )
  expect_error(createClones(population[, -1], data, clonesArgs(graceEnd = 7)), "rowId")
  expect_error(createClones(population, data, list()), "createClonesArgs")
})

test_that("createClones on an empty population", {
  data <- createM2FixtureData()
  on.exit(close(data))
  population <- createM2FixturePopulation(data, 0)
  empty <- population[0, ]
  attr(empty, "metaData") <- attr(population, "metaData")
  result <- runM2Clones(
    population = empty,
    data = data,
    strategies = list(m2Strategy(1, "initiateWithin", 7), m2Strategy(2, "notInitiateWithin", 7))
  )
  expect_equal(nrow(result$persons), 0)
  expect_equal(nrow(result$clones), 0)
  expect_equal(colnames(result$clones), names(.clonesColumns))
})

test_that("createClones equals the slow reference on 1,000 random persons", {
  tables <- simulateM2Tables(n = 1000, seed = 123)
  data <- do.call(createCloneCensorWeightData, tables)
  on.exit(close(data))
  settings <- list(
    list(riskWindowStart = 0, riskWindowEnd = 30, endAnchor = "cohort start", maxDaysAtRisk = 99999, minDaysAtRisk = 1, removeSubjectsWithPriorOutcome = TRUE,
         strategies = list(m2Strategy(1, "initiateWithin", 7), m2Strategy(2, "notInitiateWithin", 7)), sameDayOrdering = "eventsFirst", competingEvents = "death"),
    list(riskWindowStart = 1, riskWindowEnd = 30, endAnchor = "cohort start", maxDaysAtRisk = 99999, minDaysAtRisk = 0, removeSubjectsWithPriorOutcome = FALSE,
         strategies = list(m2Strategy(1, "initiateWithin", 7), m2Strategy(2, "notInitiateWithin", 7)), sameDayOrdering = "deviationFirst", competingEvents = "death"),
    list(riskWindowStart = 0, riskWindowEnd = 10, endAnchor = "cohort end", maxDaysAtRisk = 99999, minDaysAtRisk = 1, removeSubjectsWithPriorOutcome = TRUE,
         strategies = list(m2Strategy(1, "initiateWithin", 3), m2Strategy(2, "neverInitiate")), sameDayOrdering = "eventsFirst", competingEvents = "none"),
    list(riskWindowStart = 1, riskWindowEnd = 25, endAnchor = "cohort start", maxDaysAtRisk = 15, minDaysAtRisk = 1, removeSubjectsWithPriorOutcome = TRUE,
         strategies = list(m2Strategy(1, "initiateWithin", 0), m2Strategy(2, "notInitiateWithin", 0)), sameDayOrdering = "deviationFirst", competingEvents = "death"),
    list(riskWindowStart = 0, riskWindowEnd = 40, endAnchor = "cohort start", maxDaysAtRisk = 99999, minDaysAtRisk = 0, removeSubjectsWithPriorOutcome = FALSE,
         strategies = list(m2Strategy(1, "neverInitiate", treatment = "comparatorTreatment"), m2Strategy(2, "initiateWithin", 5, treatment = "comparatorTreatment")),
         sameDayOrdering = "eventsFirst", competingEvents = "death"),
    list(riskWindowStart = 1, riskWindowEnd = 5, endAnchor = "cohort end", maxDaysAtRisk = 12, minDaysAtRisk = 0, removeSubjectsWithPriorOutcome = TRUE,
         strategies = list(m2Strategy(1, "notInitiateWithin", 10, treatment = "comparatorTreatment"), m2Strategy(2, "neverInitiate")),
         sameDayOrdering = "deviationFirst", competingEvents = "none")
  )
  endReasons <- character(0)
  droppedClones <- 0
  tiedClones <- 0
  for (setting in settings) {
    studyPopulationArgs <- CohortMethod::createCreateStudyPopulationArgs(
      removeSubjectsWithPriorOutcome = setting$removeSubjectsWithPriorOutcome,
      minDaysAtRisk = setting$minDaysAtRisk,
      maxDaysAtRisk = setting$maxDaysAtRisk,
      riskWindowStart = setting$riskWindowStart,
      riskWindowEnd = setting$riskWindowEnd,
      endAnchor = setting$endAnchor
    )
    population <- suppressMessages(createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = studyPopulationArgs))
    clonesArgs <- createCreateClonesArgs(
      strategies = setting$strategies,
      sameDayOrdering = setting$sameDayOrdering,
      competingEvents = setting$competingEvents
    )
    result <- runM2Clones(
      population = population,
      data = data,
      strategies = setting$strategies,
      sameDayOrdering = setting$sameDayOrdering,
      competingEvents = setting$competingEvents
    )
    reference <- referenceM2Clones(tables, population, studyPopulationArgs, clonesArgs)
    expect_equal(result$persons, reference$persons)
    expect_equal(result$clones, reference$clones)

    endReasons <- c(endReasons, result$clones$endReason)
    droppedClones <- droppedClones + 2 * nrow(result$persons) - nrow(result$clones)
    person <- match(result$clones$rowId, result$persons$rowId)
    tied <- result$clones$deviationDay == result$clones$endDay &
      (result$persons$outcomeDay[person] == result$clones$endDay | result$persons$deathDay[person] == result$clones$endDay)
    tiedClones <- tiedClones + sum(tied, na.rm = TRUE)
  }
  # The random data exercise every end reason, dropped clones, and deviations on an event day.
  expect_setequal(unique(endReasons), .endReasons)
  expect_gt(droppedClones, 0)
  expect_gt(tiedClones, 0)
})

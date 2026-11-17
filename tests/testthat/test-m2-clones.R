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

# Expected clones per strategy, ordering, and riskStart, with competingEvents = "death". With
# "deviationFirst", a clone that ends because of a deviation on day d is last under follow-up on day
# d - 1, and is dropped when d - 1 < riskStart.
expectedM2Clones <- list(
  initiateWithin7 = list(
    eventsFirst = list(
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
      )
    ),
    deviationFirst = list(
      riskStart0 = dplyr::tribble(
        ~rowId, ~deviationDay, ~endDay, ~endReason,
        1,      NA,            30,      "riskEnd",
        2,      NA,            25,      "outcome",
        3,      7,             6,       "deviation",
        4,      7,             6,       "deviation",
        5,      NA,            3,       "outcome",
        6,      NA,            4,       "death",
        7,      NA,            5,       "observationEnd",
        8,      NA,            4,       "outcome",
        9,      7,             6,       "deviation",
        11,     7,             0,       "outcome",
        12,     NA,            3,       "death"
      ),
      riskStart1 = dplyr::tribble(
        ~rowId, ~deviationDay, ~endDay, ~endReason,
        1,      NA,            30,      "riskEnd",
        2,      NA,            25,      "outcome",
        3,      7,             6,       "deviation",
        4,      7,             6,       "deviation",
        5,      NA,            3,       "outcome",
        6,      NA,            4,       "death",
        7,      NA,            5,       "observationEnd",
        8,      NA,            4,       "outcome",
        9,      7,             6,       "deviation",
        12,     NA,            3,       "death"
      )
    )
  ),
  notInitiateWithin7 = list(
    eventsFirst = list(
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
      )
    ),
    deviationFirst = list(
      riskStart0 = dplyr::tribble(
        ~rowId, ~deviationDay, ~endDay, ~endReason,
        2,      7,             6,       "deviation",
        3,      NA,            30,      "riskEnd",
        4,      NA,            7,       "outcome",
        5,      5,             3,       "outcome",
        6,      NA,            4,       "death",
        7,      NA,            5,       "observationEnd",
        8,      4,             3,       "deviation",
        9,      NA,            7,       "observationEnd",
        11,     NA,            0,       "outcome",
        12,     3,             2,       "deviation"
      ),
      riskStart1 = dplyr::tribble(
        ~rowId, ~deviationDay, ~endDay, ~endReason,
        2,      7,             6,       "deviation",
        3,      NA,            30,      "riskEnd",
        4,      NA,            7,       "outcome",
        5,      5,             3,       "outcome",
        6,      NA,            4,       "death",
        7,      NA,            5,       "observationEnd",
        8,      4,             3,       "deviation",
        9,      NA,            7,       "observationEnd",
        12,     3,             2,       "deviation"
      )
    )
  ),
  neverInitiate = list(
    eventsFirst = list(
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
      )
    ),
    deviationFirst = list(
      riskStart0 = dplyr::tribble(
        ~rowId, ~deviationDay, ~endDay, ~endReason,
        2,      7,             6,       "deviation",
        3,      12,            11,      "deviation",
        4,      NA,            7,       "outcome",
        5,      5,             3,       "outcome",
        6,      NA,            4,       "death",
        7,      NA,            5,       "observationEnd",
        8,      4,             3,       "deviation",
        9,      NA,            7,       "observationEnd",
        11,     NA,            0,       "outcome",
        12,     3,             2,       "deviation"
      ),
      riskStart1 = dplyr::tribble(
        ~rowId, ~deviationDay, ~endDay, ~endReason,
        2,      7,             6,       "deviation",
        3,      12,            11,      "deviation",
        4,      NA,            7,       "outcome",
        5,      5,             3,       "outcome",
        6,      NA,            4,       "death",
        7,      NA,            5,       "observationEnd",
        8,      4,             3,       "deviation",
        9,      NA,            7,       "observationEnd",
        12,     3,             2,       "deviation"
      )
    )
  ),
  initiateWithin0 = list(
    eventsFirst = list(
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
      )
    ),
    deviationFirst = list(
      riskStart0 = dplyr::tribble(
        ~rowId, ~deviationDay, ~endDay, ~endReason,
        1,      NA,            30,      "riskEnd"
      ),
      riskStart1 = dplyr::tribble(
        ~rowId, ~deviationDay, ~endDay, ~endReason,
        1,      NA,            30,      "riskEnd"
      )
    )
  ),
  notInitiateWithin0 = list(
    eventsFirst = list(
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
      )
    ),
    deviationFirst = list(
      riskStart0 = dplyr::tribble(
        ~rowId, ~deviationDay, ~endDay, ~endReason,
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
      )
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
  expected <- expectedM2Clones[[name]][[sameDayOrdering]][[paste0("riskStart", riskStart)]]
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

test_that("Fixture clones of initiateWithin with graceStart > 0: a start before graceStart is a deviation", {
  data <- createM2FixtureData()
  on.exit(close(data))
  population <- createM2FixturePopulation(data, 0)
  strategies <- list(m2Strategy(1, "initiateWithin", 7, graceStart = 4), m2Strategy(2, "notInitiateWithin", 7))
  # Person 1 starts on day 0 and person 12 on day 3 (the day of death), before graceStart = 4; persons 2
  # (day 7), 5 (day 5, after its outcome on day 3) and 8 (day 4) start within days 4..7; the others
  # deviate on day 7 as with graceStart = 0 (expectedM2Clones$initiateWithin7).
  expected <- list(
    eventsFirst = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      1,      0,             0,       "deviation",
      2,      NA,            25,      "outcome",
      3,      7,             7,       "deviation",
      4,      7,             7,       "outcome",
      5,      NA,            3,       "outcome",
      6,      NA,            4,       "death",
      7,      NA,            5,       "observationEnd",
      8,      NA,            4,       "outcome",
      9,      7,             7,       "deviation",
      11,     7,             0,       "outcome",
      12,     3,             3,       "death"
    ),
    # Person 1 deviates on day 0 before any follow-up, so its clone is dropped; person 12 deviates on
    # day 3 before its death that day:
    deviationFirst = dplyr::tribble(
      ~rowId, ~deviationDay, ~endDay, ~endReason,
      2,      NA,            25,      "outcome",
      3,      7,             6,       "deviation",
      4,      7,             6,       "deviation",
      5,      NA,            3,       "outcome",
      6,      NA,            4,       "death",
      7,      NA,            5,       "observationEnd",
      8,      NA,            4,       "outcome",
      9,      7,             6,       "deviation",
      11,     7,             0,       "outcome",
      12,     3,             2,       "deviation"
    )
  )
  for (sameDayOrdering in names(expected)) {
    result <- runM2Clones(population, data, strategies, sameDayOrdering)
    initiate <- result$clones[result$clones$strategyId == 1, c("rowId", "deviationDay", "endDay", "endReason")]
    rownames(initiate) <- NULL
    expect_equal(initiate, as.data.frame(dplyr::mutate(expected[[sameDayOrdering]], dplyr::across(c("rowId", "deviationDay", "endDay"), as.integer))))
    expect_equal(result$metaData$strategies[[1]]$graceStart, 4)
  }
})

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
      replaceM2Rows(
        expectedM2Clones$initiateWithin7$eventsFirst[[paste0("riskStart", riskStart)]],
        noneReplacements[[1]]
      ) |>
        dplyr::mutate(strategyId = 1),
      replaceM2Rows(
        expectedM2Clones$notInitiateWithin7$eventsFirst[[paste0("riskStart", riskStart)]],
        noneReplacements[[2]]
      ) |>
        dplyr::mutate(strategyId = 2)
    )
    expect_equal(result$clones, completeM2ExpectedClones(expected, riskStart))
    expectedPersons <- expectedM2Persons[[paste0("riskStart", riskStart)]]
    expectedPersons$deathDay <- NA
    expect_equal(result$persons, toM2Persons(expectedPersons, riskStart))
  }
})

test_that("Follow-up ends at death when observation continues, with and without competing death", {
  # Observation (daysToObsEnd) continues after death for persons 1 to 5:
  # 1: death on day 5, starts on day 2
  # 2: death on day 12, after G = 7; never starts
  # 3: death on day 0 (before riskStart when riskStart = 1); never starts
  # 4: death on day 6; an outcome on day 9, after death, is ignored; never starts
  # 5: death on day 3, the day of the start
  # 6: no death; starts on day 10; observation and the risk window end on day 30
  data <- createCloneCensorWeightData(
    cohorts = data.frame(
      rowId = 1:6,
      personSeqId = 101:106,
      personId = as.character(1001:1006),
      cohortStartDate = as.Date("2021-01-01"),
      daysFromObsStart = 500L,
      daysToCohortEnd = 0L,
      daysToObsEnd = c(30L, 40L, 30L, 30L, 20L, 30L),
      daysToDeath = c(5L, 12L, 0L, 6L, 3L, NA)
    ),
    outcomes = data.frame(rowId = 4L, outcomeId = 3, daysToEvent = 9L),
    treatments = data.frame(
      rowId = c(1L, 5L, 6L),
      treatmentId = 21,
      daysToStart = c(2L, 3L, 10L),
      daysToEnd = c(12L, 13L, 20L)
    ),
    metaData = list(targetId = 1, treatmentIds = c(treatment = 21))
  )
  on.exit(close(data))
  expectedClones <- list(
    # Death is not an event: it ends follow-up as the end of observation.
    none = dplyr::tribble(
      ~rowId, ~strategyId, ~deviationDay, ~endDay, ~endReason,
      1,      1,           NA,            5,       "observationEnd",
      1,      2,           2,             2,       "deviation",
      2,      1,           7,             7,       "deviation",
      2,      2,           NA,            12,      "observationEnd",
      3,      1,           NA,            0,       "observationEnd",
      3,      2,           NA,            0,       "observationEnd",
      4,      1,           NA,            6,       "observationEnd",
      4,      2,           NA,            6,       "observationEnd",
      5,      1,           NA,            3,       "observationEnd",
      5,      2,           3,             3,       "deviation",
      6,      1,           7,             7,       "deviation",
      6,      2,           NA,            30,      "riskEnd"
    ),
    death = dplyr::tribble(
      ~rowId, ~strategyId, ~deviationDay, ~endDay, ~endReason,
      1,      1,           NA,            5,       "death",
      1,      2,           2,             2,       "deviation",
      2,      1,           7,             7,       "deviation",
      2,      2,           NA,            12,      "death",
      3,      1,           NA,            0,       "death",
      3,      2,           NA,            0,       "death",
      4,      1,           NA,            6,       "death",
      4,      2,           NA,            6,       "death",
      5,      1,           NA,            3,       "death",
      5,      2,           3,             3,       "death",
      6,      1,           7,             7,       "deviation",
      6,      2,           NA,            30,      "riskEnd"
    )
  )
  for (riskStart in 0:1) {
    population <- suppressMessages(createStudyPopulation(
      cloneCensorWeightData = data,
      outcomeId = 3,
      createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
        riskWindowStart = riskStart,
        riskWindowEnd = 30,
        endAnchor = "cohort start"
      )
    ))
    for (competingEvents in c("none", "death")) {
      result <- runM2Clones(
        population = population,
        data = data,
        strategies = list(m2Strategy(1, "initiateWithin", 7), m2Strategy(2, "notInitiateWithin", 7)),
        competingEvents = competingEvents
      )
      deathDay <- if (competingEvents == "death") c(5, 12, 0, 6, 3, NA) else rep(NA, 6)
      if (riskStart == 1) {
        deathDay[3] <- NA
      }
      expectedPersons <- data.frame(
        rowId = 1:6,
        riskEnd = c(5, 12, 0, 6, 3, 30),
        nominalRiskEnd = 30,
        outcomeDay = NA,
        deathDay = deathDay,
        treatmentDay = c(2, NA, NA, NA, 3, 10),
        comparatorTreatmentDay = NA
      )
      expect_equal(result$persons, toM2Persons(expectedPersons, riskStart), info = competingEvents)
      # Person 3 has no follow-up when riskStart = 1.
      expected <- expectedClones[[competingEvents]]
      expected <- expected[expected$rowId != 3 | riskStart == 0, ]
      expect_equal(result$clones, completeM2ExpectedClones(expected, riskStart), info = competingEvents)
    }
  }
})

test_that("A risk window shortened at the next risk window ends with 'riskEnd', not 'observationEnd'", {
  # Rows 1 and 2 are two index dates of one person, 10 days apart, so censorAtNewRiskWindow ends the
  # risk window of row 1 on day 9 while observation continues. Observation ends on day 12 for row 3.
  data <- createCloneCensorWeightData(
    cohorts = data.frame(
      rowId = 1:3,
      personSeqId = c(101L, 101L, 102L),
      personId = c("1001", "1001", "1002"),
      cohortStartDate = as.Date("2021-01-01") + c(0, 10, 0),
      daysFromObsStart = 500L,
      daysToCohortEnd = 0L,
      daysToObsEnd = c(100L, 90L, 12L),
      daysToDeath = NA_integer_
    ),
    outcomes = data.frame(rowId = integer(0), outcomeId = numeric(0), daysToEvent = integer(0)),
    treatments = data.frame(
      rowId = integer(0),
      treatmentId = numeric(0),
      daysToStart = integer(0),
      daysToEnd = integer(0)
    ),
    metaData = list(targetId = 1, treatmentIds = c(treatment = 21))
  )
  on.exit(close(data))
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
      riskWindowEnd = 30,
      endAnchor = "cohort start",
      censorAtNewRiskWindow = TRUE
    )
  ))
  result <- runM2Clones(
    population = population,
    data = data,
    strategies = list(m2Strategy(1, "notInitiateWithin", 7), m2Strategy(2, "neverInitiate"))
  )
  expect_equal(result$persons$riskEnd, c(9L, 30L, 12L))
  expect_equal(result$persons$nominalRiskEnd, c(30L, 30L, 30L))
  expect_equal(result$clones[, c("rowId", "strategyId", "endDay", "endReason")], data.frame(
    rowId = rep(1:3, each = 2),
    strategyId = rep(1:2, 3),
    endDay = rep(c(9L, 30L, 12L), each = 2),
    endReason = rep(c("riskEnd", "riskEnd", "observationEnd"), each = 2)
  ))
})

test_that("Clones and the contract weights give the true risk with no confounding and no effect", {
  # The initiation day A and the outcome day Y are exactly independent: each pair (A, Y) appears
  # aCount * yCount times (108 persons), and Y does not depend on A. The risk of every strategy is then
  # P(riskStart <= Y <= h | Y >= riskStart). The censoring weights W(d) of contract section 5.4,
  # computed from the true P(A > d), are the same for all clones of a strategy on a day, so the
  # weighted Kaplan-Meier over the clones table equals that risk exactly under both orderings only if a
  # clone is in the risk set of a day exactly when its weight says so. With "deviationFirst", a clone
  # that deviates on day d must not be in the risk set of day d.
  aDays <- c(0L, 1L, 2L, 3L, NA)
  aCount <- c(2, 1, 2, 1, 3)
  yDays <- c(0L, 1L, 2L, 3L, 4L, 5L, NA)
  yCount <- c(1, 2, 1, 1, 2, 1, 4)
  pairs <- expand.grid(a = seq_along(aDays), y = seq_along(yDays))
  pairs <- pairs[rep(seq_len(nrow(pairs)), aCount[pairs$a] * yCount[pairs$y]), ]
  a <- aDays[pairs$a]
  y <- yDays[pairs$y]
  n <- length(a)
  expect_equal(n, 108)
  data <- createCloneCensorWeightData(
    cohorts = data.frame(
      rowId = seq_len(n),
      personSeqId = seq_len(n),
      personId = as.character(seq_len(n)),
      cohortStartDate = as.Date("2020-01-01"),
      daysFromObsStart = 400L,
      daysToCohortEnd = 0L,
      daysToObsEnd = 10L,
      daysToDeath = NA_integer_
    ),
    outcomes = data.frame(rowId = which(!is.na(y)), outcomeId = 3, daysToEvent = y[!is.na(y)]),
    treatments = data.frame(
      rowId = which(!is.na(a)),
      treatmentId = 21,
      daysToStart = a[!is.na(a)],
      daysToEnd = a[!is.na(a)] + 5L
    ),
    metaData = list(targetId = 1, treatmentIds = c(treatment = 21))
  )
  on.exit(close(data))
  graceEnd <- 2
  # S(d) = P(A > d), no initiation on days 0..d (S(-1) = 1).
  survivalA <- function(d) sum(aCount[is.na(aDays) | aDays > d]) / sum(aCount)
  # Contract section 5.4: S(d - 1) with "eventsFirst", S(d) with "deviationFirst"; for initiateWithin,
  # 1/pi from day G + 1 ("eventsFirst") or day G ("deviationFirst"), with pi = P(A <= G).
  weight <- function(type, d, sameDayOrdering) {
    lag <- if (sameDayOrdering == "eventsFirst") 1 else 0
    switch(type,
      initiateWithin = if (d < graceEnd + lag) 1 else 1 / (1 - survivalA(graceEnd)),
      notInitiateWithin = 1 / survivalA(min(d - lag, graceEnd)),
      neverInitiate = 1 / survivalA(d - lag)
    )
  }
  for (riskStart in 0:1) {
    # Persons with the outcome before riskStart are removed as prior outcomes.
    population <- suppressMessages(createStudyPopulation(
      cloneCensorWeightData = data,
      outcomeId = 3,
      createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
        riskWindowStart = riskStart,
        riskWindowEnd = 10,
        endAnchor = "cohort start"
      )
    ))
    atRiskY <- is.na(yDays) | yDays >= riskStart
    trueRisk <- function(h) sum(yCount[!is.na(yDays) & yDays >= riskStart & yDays <= h]) / sum(yCount[atRiskY])
    for (pair in list(c("initiateWithin", "notInitiateWithin"), c("initiateWithin", "neverInitiate"))) {
      strategies <- lapply(1:2, function(i) m2Strategy(i, pair[i], if (pair[i] == "neverInitiate") NULL else graceEnd))
      for (sameDayOrdering in c("eventsFirst", "deviationFirst")) {
        result <- runM2Clones(population, data, strategies, sameDayOrdering = sameDayOrdering)
        for (strategyId in 1:2) {
          clones <- result$clones[result$clones$strategyId == strategyId, ]
          survival <- 1
          for (d in riskStart:6) {
            w <- weight(pair[strategyId], d, sameDayOrdering)
            atRisk <- w * sum(clones$endDay >= d)
            events <- w * sum(clones$endDay == d & clones$endReason == "outcome")
            survival <- survival * (1 - events / atRisk)
            expect_equal(1 - survival, trueRisk(d), tolerance = 1e-12, info = sprintf(
              "%s, %s, riskStart = %d, day %d", pair[strategyId], sameDayOrdering, riskStart, d
            ))
          }
        }
      }
    }
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
  # Every person is cloned into both strategies, so the population rows count them in both columns:
  expectedPopulationAttrition <- populationAttrition
  expectedPopulationAttrition$comparatorPersons <- populationAttrition$targetPersons
  expectedPopulationAttrition$comparatorExposures <- populationAttrition$targetExposures
  expect_equal(attrition[seq_len(nrow(populationAttrition)), ], expectedPopulationAttrition)
  cloneAttrition <- as.data.frame(attrition[nrow(populationAttrition) + 1:2, ])
  expect_equal(cloneAttrition, data.frame(
    description = c("Clones created", "Clones with follow-up"),
    targetPersons = c(10L, 1L),
    comparatorPersons = c(10L, 9L),
    targetExposures = c(10L, 1L),
    comparatorExposures = c(10L, 9L)
  ))

  # drawAttritionDiagram() removes previous minus current persons at each step: never a negative number,
  # nothing when the clones are created, and the clones without follow-up of each strategy after that.
  removed <- -diff(as.matrix(attrition[, c("targetPersons", "comparatorPersons")]))
  expect_true(all(removed >= 0))
  expect_equal(unname(removed[nrow(removed) - 1, ]), c(0L, 0L))
  expect_equal(unname(removed[nrow(removed), ]), c(9L, 1L))
  skip_if_not_installed("ggplot2")
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = strategies)))
  on.exit(close(clones), add = TRUE)
  diagram <- drawAttritionDiagram(clones)
  labels <- unlist(lapply(ggplot2::ggplot_build(diagram)$data, function(layer) as.character(layer$label)))
  expect_true(any(grepl("Comparator: n = 9", labels, fixed = TRUE)))
  expect_false(any(grepl("n = -", labels, fixed = TRUE)))
})

test_that("createClones stops for unsupported strategies and invalid input", {
  data <- createM2FixtureData()
  on.exit(close(data))
  population <- createM2FixturePopulation(data, 0)
  clonesArgs <- function(...) {
    createCreateClonesArgs(strategies = list(createTreatmentStrategy(1, "Initiate", "initiateWithin", ...), m2Strategy(2, "neverInitiate")))
  }
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

test_that("createClones equals the slow reference on 1,000 random persons (500 outside the full tier)", {
  # (500 persons still exercise every case counted at the end.)
  tables <- simulateM2Tables(n = testConfigurations(1000, 500), seed = 123)
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
         sameDayOrdering = "deviationFirst", competingEvents = "none"),
    list(riskWindowStart = 0, riskWindowEnd = 30, endAnchor = "cohort start", maxDaysAtRisk = 99999, minDaysAtRisk = 1, removeSubjectsWithPriorOutcome = TRUE,
         censorAtNewRiskWindow = TRUE, strategies = list(m2Strategy(1, "initiateWithin", 7), m2Strategy(2, "neverInitiate")),
         sameDayOrdering = "deviationFirst", competingEvents = "none"),
    # graceStart > 0 (a start before day graceStart is a deviation), under both orderings and with
    # riskStart after some early starts:
    list(riskWindowStart = 0, riskWindowEnd = 30, endAnchor = "cohort start", maxDaysAtRisk = 99999, minDaysAtRisk = 1, removeSubjectsWithPriorOutcome = TRUE,
         strategies = list(m2Strategy(1, "initiateWithin", 7, graceStart = 3), m2Strategy(2, "notInitiateWithin", 7)), sameDayOrdering = "eventsFirst", competingEvents = "death"),
    list(riskWindowStart = 2, riskWindowEnd = 30, endAnchor = "cohort start", maxDaysAtRisk = 99999, minDaysAtRisk = 0, removeSubjectsWithPriorOutcome = FALSE,
         strategies = list(m2Strategy(1, "initiateWithin", 9, graceStart = 4, treatment = "comparatorTreatment"), m2Strategy(2, "initiateWithin", 4)),
         sameDayOrdering = "deviationFirst", competingEvents = "death")
  )
  endReasons <- character(0)
  droppedClones <- 0
  tiedClones <- c(eventsFirst = 0, deviationFirst = 0)
  shortenedRiskWindows <- 0
  deathsEndingObservation <- 0
  for (setting in settings) {
    studyPopulationArgs <- CohortMethod::createCreateStudyPopulationArgs(
      removeSubjectsWithPriorOutcome = setting$removeSubjectsWithPriorOutcome,
      minDaysAtRisk = setting$minDaysAtRisk,
      maxDaysAtRisk = setting$maxDaysAtRisk,
      riskWindowStart = setting$riskWindowStart,
      riskWindowEnd = setting$riskWindowEnd,
      endAnchor = setting$endAnchor,
      censorAtNewRiskWindow = isTRUE(setting$censorAtNewRiskWindow)
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
    clones <- result$clones
    person <- match(clones$rowId, result$persons$rowId)
    eventOnDeviationDay <- clones$deviationDay == result$persons$outcomeDay[person] |
      clones$deviationDay == result$persons$deathDay[person]
    if (setting$sameDayOrdering == "eventsFirst") {
      tied <- eventOnDeviationDay & clones$endDay == clones$deviationDay & clones$endReason != "deviation"
    } else {
      tied <- eventOnDeviationDay & clones$endDay == clones$deviationDay - 1 & clones$endReason == "deviation"
    }
    tiedClones[setting$sameDayOrdering] <- tiedClones[setting$sameDayOrdering] + sum(tied, na.rm = TRUE)
    shortenedRiskWindows <- shortenedRiskWindows +
      sum(clones$endReason == "riskEnd" & clones$endDay < result$persons$nominalRiskEnd[person])
    cohort <- tables$cohorts[match(clones$rowId, tables$cohorts$rowId), ]
    deathsEndingObservation <- deathsEndingObservation + sum(
      setting$competingEvents == "none" & clones$endReason == "observationEnd" &
        clones$endDay == cohort$daysToDeath & cohort$daysToDeath < cohort$daysToObsEnd,
      na.rm = TRUE
    )
  }
  # The random data exercise every end reason, dropped clones, deviations on an event day under both
  # orderings, risk windows shortened at the next risk window, and deaths before the end of
  # observation with competingEvents = "none".
  expect_setequal(unique(endReasons), .endReasons)
  expect_gt(droppedClones, 0)
  expect_gt(tiedClones[["eventsFirst"]], 0)
  expect_gt(tiedClones[["deviationFirst"]], 0)
  expect_gt(shortenedRiskWindows, 0)
  expect_gt(deathsEndingObservation, 0)
})

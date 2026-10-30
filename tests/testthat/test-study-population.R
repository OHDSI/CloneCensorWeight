make_data <- function() {
  # 10 subjects, cohort start 2020-01-01, cohort end varies
  cohorts <- data.frame(
    subjectId       = 1:10,
    cohortStartDate = as.Date("2020-01-01"),
    cohortEndDate   = as.Date("2020-01-01") + c(rep(365L, 8), 0L, 3L),
    initiationDate  = as.Date(NA),
    stringsAsFactors = FALSE
  )
  # Subjects 3 and 7 have a prior outcome (before cohort start)
  # Subjects 1 and 5 have an outcome during follow-up
  outcomes <- data.frame(
    subjectId   = c(3L, 7L, 1L, 5L),
    outcomeDate = as.Date(c("2019-06-01", "2019-12-31", "2020-06-01", "2020-08-01")),
    outcomeId   = rep(1L, 4),
    stringsAsFactors = FALSE
  )
  createCloneCensorWeightData(cohorts, outcomes)
}

test_that("createStudyPopulation removes subjects with prior outcome", {
  data <- make_data()
  args <- createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = TRUE,
    priorOutcomeLookback           = 99999,
    minDaysAtRisk                  = 1,
    riskWindowEnd                  = 365,
    endAnchor                      = "cohort start"
  )
  pop <- createStudyPopulation(data, outcomeId = 1L, createStudyPopulationArgs = args)
  cohorts <- dplyr::collect(pop$cohorts)

  # Subjects 3 and 7 had prior outcomes; 9 and 10 have 0/3 days (< 1 min)
  expect_false(3L %in% cohorts$subjectId)
  expect_false(7L %in% cohorts$subjectId)
  # Subjects 1 and 5 had outcomes AFTER cohort start — should be retained
  expect_true(1L %in% cohorts$subjectId)
  expect_true(5L %in% cohorts$subjectId)
})

test_that("prior outcome lookback window is respected", {
  # Subject 7's prior outcome is on 2019-12-31 (1 day before cohort start).
  # Each call gets a fresh Andromeda — createStudyPopulation modifies in place.
  args_0d <- createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = TRUE,
    priorOutcomeLookback           = 0,
    minDaysAtRisk                  = 0,
    riskWindowEnd                  = 365,
    endAnchor                      = "cohort start"
  )
  args_1d <- createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = TRUE,
    priorOutcomeLookback           = 1,
    minDaysAtRisk                  = 0,
    riskWindowEnd                  = 365,
    endAnchor                      = "cohort start"
  )
  pop_0d <- createStudyPopulation(make_data(), outcomeId = 1L,
                                  createStudyPopulationArgs = args_0d)
  pop_1d <- createStudyPopulation(make_data(), outcomeId = 1L,
                                  createStudyPopulationArgs = args_1d)

  expect_true(7L  %in% dplyr::collect(pop_0d$cohorts)$subjectId)
  expect_false(7L %in% dplyr::collect(pop_1d$cohorts)$subjectId)
})

test_that("risk window trims cohortEndDate", {
  data <- make_data()
  args <- createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = FALSE,
    minDaysAtRisk                  = 0,
    riskWindowEnd                  = 180,
    endAnchor                      = "cohort start"
  )
  pop    <- createStudyPopulation(data, outcomeId = 1L, createStudyPopulationArgs = args)
  cohorts <- dplyr::collect(pop$cohorts)

  # All retained subjects should have follow-up capped at day 180
  endDays <- as.integer(
    as.Date(cohorts$cohortEndDate) - as.Date(cohorts$cohortStartDate)
  )
  expect_true(all(endDays <= 180))
})

test_that("minDaysAtRisk removes short-follow-up subjects", {
  data <- make_data()
  args <- createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = FALSE,
    minDaysAtRisk                  = 30,
    riskWindowEnd                  = 365,
    endAnchor                      = "cohort start"
  )
  pop     <- createStudyPopulation(data, outcomeId = 1L, createStudyPopulationArgs = args)
  cohorts <- dplyr::collect(pop$cohorts)

  # Subjects 9 (0 days) and 10 (3 days) should be excluded
  expect_false(9L %in% cohorts$subjectId)
  expect_false(10L %in% cohorts$subjectId)
})

test_that("getAttritionTable returns correct counts", {
  data <- make_data()
  args <- createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = TRUE,
    priorOutcomeLookback           = 99999,
    minDaysAtRisk                  = 1,
    riskWindowEnd                  = 365,
    endAnchor                      = "cohort start"
  )
  pop      <- createStudyPopulation(data, outcomeId = 1L, createStudyPopulationArgs = args)
  attrition <- getAttritionTable(pop)

  expect_s3_class(attrition, "data.frame")
  expect_true(all(c("description", "nSubjects", "nExcluded") %in% names(attrition)))
  # First step: all 10 subjects, 0 excluded
  expect_equal(attrition$nSubjects[1], 10)
  expect_equal(attrition$nExcluded[1], 0)
  # Prior outcome step removes subjects 3 and 7
  prior_row <- attrition[grepl("prior", attrition$description, ignore.case = TRUE), ]
  expect_equal(prior_row$nExcluded, 2)
  # Excluded counts should be non-negative
  expect_true(all(attrition$nExcluded >= 0))
})

test_that("getAttritionTable errors if createStudyPopulation not run", {
  data <- make_data()
  expect_error(getAttritionTable(data), "createStudyPopulation")
})

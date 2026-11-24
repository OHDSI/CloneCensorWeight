library(testthat)

test_that("createCloneCensorWeightData builds the object with the contracted types", {
  data <- createFoundationData()
  on.exit(close(data))

  expect_true(isCloneCensorWeightData(data))
  expect_true(is(data, "CohortMethodData"))
  expect_true(is(data, "CovariateData"))
  expect_true(all(c("cohorts", "outcomes", "treatments", "covariates", "covariateRef", "analysisRef") %in% names(data)))

  cohorts <- dplyr::collect(data$cohorts)
  expect_equal(nrow(cohorts), 200)
  expect_type(cohorts$rowId, "integer")
  expect_type(cohorts$personSeqId, "integer")
  expect_type(cohorts$personId, "character")
  expect_true(all(cohorts$treatment == 1L))
  expect_type(cohorts$treatment, "integer")
  expect_s3_class(cohorts$cohortStartDate, "Date")
  for (column in c("daysFromObsStart", "daysToCohortEnd", "daysToObsEnd", "daysToDeath")) {
    expect_type(cohorts[[column]], "integer")
  }
  expect_equal(sum(!is.na(cohorts$daysToDeath)), 20)

  outcomes <- dplyr::collect(data$outcomes)
  expect_type(outcomes$rowId, "integer")
  expect_type(outcomes$outcomeId, "double")
  expect_type(outcomes$daysToEvent, "integer")
  treatments <- dplyr::collect(data$treatments)
  expect_type(treatments$treatmentId, "double")
  expect_type(treatments$daysToStart, "integer")
  covariates <- dplyr::collect(data$covariates)
  expect_type(covariates$rowId, "integer")
  expect_type(covariates$covariateId, "double")
  expect_type(covariates$covariateValue, "double")

  metaData <- attr(data, "metaData")
  expect_equal(metaData$targetId, 1)
  expect_equal(metaData$treatmentIds, c(treatment = 21))
  expect_equal(metaData$outcomeIds, 3)
  expect_equal(metaData$populationSize, 200)
  expect_equal(
    colnames(metaData$attrition),
    c("description", "targetPersons", "comparatorPersons", "targetExposures", "comparatorExposures")
  )
  expect_equal(metaData$attrition$targetPersons, 200)
})

test_that("createCloneCensorWeightData fills empty tables with the right columns", {
  tables <- createFoundationTables(n = 5)
  data <- createCloneCensorWeightData(
    cohorts = tables$cohorts,
    outcomes = NULL,
    treatments = data.frame(),
    metaData = tables$metaData
  )
  on.exit(close(data))

  outcomes <- dplyr::collect(data$outcomes)
  expect_equal(nrow(outcomes), 0)
  expect_equal(colnames(outcomes), c("rowId", "outcomeId", "daysToEvent"))
  expect_type(outcomes$daysToEvent, "integer")
  expect_equal(colnames(dplyr::collect(data$treatments)), c("rowId", "treatmentId", "daysToStart", "daysToEnd"))
  expect_equal(colnames(dplyr::collect(data$covariates)), c("rowId", "covariateId", "covariateValue"))
  expect_equal(colnames(dplyr::collect(data$covariateRef)), c("covariateId", "covariateName", "analysisId", "conceptId"))
  expect_equal(
    colnames(dplyr::collect(data$analysisRef)),
    c("analysisId", "analysisName", "domainId", "startDay", "endDay", "isBinary", "missingMeansZero")
  )
  expect_equal(attr(data, "metaData")$outcomeIds, numeric(0))
  expect_true(all(dplyr::pull(data$cohorts, "treatment") == 1L))
})

test_that("createCloneCensorWeightData checks columns, types and metaData", {
  tables <- createFoundationTables(n = 10)
  build <- function(...) {
    args <- tables
    args[names(list(...))] <- list(...)
    data <- do.call(createCloneCensorWeightData, args)
    close(data)
  }

  cohorts <- tables$cohorts
  cohorts$daysToDeath <- NULL
  expect_error(build(cohorts = cohorts), "missing column\\(s\\): daysToDeath")

  cohorts <- tables$cohorts
  cohorts$personId <- as.numeric(cohorts$personId)
  expect_error(build(cohorts = cohorts), "'personId' of table 'cohorts' must be of type character")

  cohorts <- tables$cohorts
  cohorts$cohortStartDate <- as.character(cohorts$cohortStartDate)
  expect_error(build(cohorts = cohorts), "'cohortStartDate' of table 'cohorts' must be of type Date")

  cohorts <- tables$cohorts
  cohorts$daysToObsEnd <- cohorts$daysToObsEnd + 0.5
  expect_error(build(cohorts = cohorts), "'daysToObsEnd' of table 'cohorts' must be of type integer")

  cohorts <- tables$cohorts
  cohorts$daysToObsEnd[1] <- NA
  expect_error(build(cohorts = cohorts), "must not contain NA")

  cohorts <- tables$cohorts
  cohorts$treatment <- 0L
  expect_error(build(cohorts = cohorts), "must be 1 in every row")

  cohorts <- tables$cohorts
  cohorts$rowId[2] <- 1L
  expect_error(build(cohorts = cohorts), "must be unique")

  outcomes <- tables$outcomes
  outcomes$outcomeId <- as.character(outcomes$outcomeId)
  expect_error(build(outcomes = outcomes), "'outcomeId' of table 'outcomes' must be of type numeric")

  expect_error(build(treatments = data.frame(rowId = 1L, treatmentId = 21)), "missing column\\(s\\): daysToStart, daysToEnd")
  expect_error(build(covariates = "covariates"), "covariates")

  expect_error(createCloneCensorWeightData(tables$cohorts, tables$outcomes, tables$treatments), "metaData")
  expect_error(build(metaData = list(treatmentIds = c(treatment = 21))), "targetId")
  expect_error(build(metaData = list(targetId = 1)), "treatmentIds")
  expect_error(build(metaData = list(targetId = 1, treatmentIds = c(drug = 21))), "treatmentIds")
  expect_error(build(metaData = list(targetId = 1, treatmentIds = c(comparatorTreatment = 22))), "'treatment' element")

  # Integer-valued doubles are accepted and stored as integers:
  cohorts <- tables$cohorts
  cohorts$rowId <- as.numeric(cohorts$rowId)
  cohorts$daysToObsEnd <- as.numeric(cohorts$daysToObsEnd)
  data <- createCloneCensorWeightData(cohorts, tables$outcomes, tables$treatments, metaData = tables$metaData)
  on.exit(close(data))
  expect_type(dplyr::pull(data$cohorts, "daysToObsEnd"), "integer")
  expect_type(dplyr::pull(data$cohorts, "rowId"), "integer")
})

test_that("createCloneCensorWeightData accepts all-NA logical columns where NA is allowed", {
  tables <- createFoundationTables(n = 10)
  tables$cohorts$daysToDeath <- NA
  tables$covariateRef$conceptId <- NA
  data <- do.call(createCloneCensorWeightData, tables)
  on.exit(close(data))
  daysToDeath <- dplyr::pull(data$cohorts, "daysToDeath")
  expect_type(daysToDeath, "integer")
  expect_true(all(is.na(daysToDeath)))
  expect_type(dplyr::pull(data$covariateRef, "conceptId"), "double")
  expect_equal(summary(data)$deaths, 0)

  # Where NA is not allowed, the error is about the NA values, not the type:
  tables$cohorts$daysToObsEnd <- NA
  expect_error(do.call(createCloneCensorWeightData, tables), "'daysToObsEnd' of table 'cohorts' must not contain NA")
})

test_that("Saving and loading CloneCensorWeightData", {
  data <- createFoundationData()
  cohorts <- dplyr::collect(data$cohorts)
  covariates <- dplyr::collect(data$covariates)
  metaData <- attr(data, "metaData")
  file <- tempfile(fileext = ".zip")
  on.exit(unlink(file))

  messages <- capture_messages(saveCloneCensorWeightData(data, file))
  expect_true(any(grepl("loadCloneCensorWeightData", messages)))
  expect_false(Andromeda::isValidAndromeda(data))
  loaded <- loadCloneCensorWeightData(file)
  on.exit(close(loaded), add = TRUE)

  expect_true(isCloneCensorWeightData(loaded))
  expect_equal(dplyr::collect(loaded$cohorts), cohorts)
  expect_equal(dplyr::collect(loaded$covariates), covariates)
  expect_equal(attr(loaded, "metaData"), metaData)

  expect_error(loadCloneCensorWeightData(tempfile()), "Cannot find file")
  expect_error(loadCloneCensorWeightData(tempdir()), "is a folder")
  expect_error(saveCloneCensorWeightData(data.frame(), file), "CloneCensorWeightData")
})

test_that("show and summary", {
  data <- createFoundationData()
  on.exit(close(data))

  messages <- capture_messages(capture.output(show(data)))
  expect_true(any(grepl("CloneCensorWeightData object", messages)))
  expect_true(any(grepl("treatment = 21", messages)))
  s <- summary(data)
  expect_s3_class(s, "summary.CloneCensorWeightData")
  expect_equal(s$persons, 200)
  expect_equal(s$indexDates, 200)
  expect_equal(s$deaths, 20)
  tables <- createFoundationTables()
  expect_equal(s$outcomeCounts$eventCount, nrow(tables$outcomes))
  expect_equal(s$treatmentCounts$personCount, length(unique(tables$treatments$rowId)))
  expect_equal(s$covariateCount, 5)
  expect_equal(s$covariateValueCount, nrow(tables$covariates))
  messages <- capture_messages(capture.output(print(s)))
  expect_true(any(grepl("CloneCensorWeightData object summary", messages)))
})

test_that("CohortMethod::createStudyPopulation runs on CloneCensorWeightData and keeps daysToDeath", {
  tables <- createFoundationTables()
  data <- do.call(createCloneCensorWeightData, tables)
  on.exit(close(data))

  studyPop <- suppressMessages(CohortMethod::createStudyPopulation(
    cohortMethodData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
      removeSubjectsWithPriorOutcome = TRUE,
      priorOutcomeLookback = 99999,
      minDaysAtRisk = 1,
      riskWindowStart = 0,
      startAnchor = "cohort start",
      riskWindowEnd = 365,
      endAnchor = "cohort start"
    )
  ))
  expect_true(all(c(
    "rowId", "personSeqId", "cohortStartDate", "daysFromObsStart", "daysToCohortEnd", "daysToObsEnd",
    "daysToDeath", "riskStart", "riskEnd", "daysToEvent", "outcomeCount", "survivalTime"
  ) %in% colnames(studyPop)))

  priorOutcomeRowIds <- tables$outcomes$rowId[tables$outcomes$daysToEvent < 0]
  expect_equal(sort(studyPop$rowId), setdiff(tables$cohorts$rowId, priorOutcomeRowIds))
  idx <- match(studyPop$rowId, tables$cohorts$rowId)
  expect_equal(studyPop$daysToDeath, tables$cohorts$daysToDeath[idx])
  expect_equal(studyPop$riskEnd, pmin(365L, tables$cohorts$daysToObsEnd[idx]))
  expect_equal(
    CohortMethod::getAttritionTable(studyPop)$targetPersons,
    c(200, nrow(studyPop), nrow(studyPop))
  )
})

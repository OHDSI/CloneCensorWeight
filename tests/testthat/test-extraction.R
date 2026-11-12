# Tests of getDbCloneCensorWeightData() on Eunomia, against a reference implementation in R
# (helper-extraction.R) applied to the raw CDM tables fetched by direct SQL.

skip_if_not_installed("Eunomia")
skip_if_not_installed("DatabaseConnector")
skip_if_not_installed("SqlRender")
skip_if_not_installed("RSQLite")

connectionDetails <- createExtractionTestDatabase()
edgePersons <- attr(connectionDetails, "edgePersons")
connection <- DatabaseConnector::connect(connectionDetails)
rawTables <- fetchExtractionTables(connection)

extractQuietly <- function(args,
                           connectionDetails = NULL,
                           connection = NULL,
                           ids = extractionIds) {
  data <- NULL
  # FeatureExtraction prints a progress bar.
  utils::capture.output(data <- suppressMessages(getDbCloneCensorWeightData(
    connectionDetails = connectionDetails,
    connection = connection,
    cdmDatabaseSchema = "main",
    targetId = ids$targetId,
    treatmentId = ids$treatmentId,
    comparatorTreatmentId = ids$comparatorTreatmentId,
    outcomeIds = ids$outcomeIds,
    getDbCloneCensorWeightDataArgs = args
  )))
  return(data)
}

collectCohorts <- function(data) {
  cohorts <- as.data.frame(collect(data$cohorts))
  cohorts <- cohorts[order(cohorts$rowId), ]
  rownames(cohorts) <- NULL
  return(cohorts)
}

expectMatchesReference <- function(data, expected) {
  expect_equal(collectCohorts(data), expected$cohorts)
  expect_equal(sortExtractionRows(collect(data$treatments)), sortExtractionRows(expected$treatments))
  expect_equal(sortExtractionRows(collect(data$outcomes)), sortExtractionRows(expected$outcomes))
  expect_equal(as.data.frame(attr(data, "metaData")$attrition), expected$attrition)
}

# Whether each hand-made entry of extractionBoundary is in the cohorts of `data`, named by case.
boundaryKept <- function(data) {
  entries <- extractionBoundary$entries
  return(stats::setNames(as.character(entries$personId) %in% collectCohorts(data)$personId, entries$case))
}

# The expected kept flags of extractionBoundary in column `column`, named by case.
expectedBoundaryKept <- function(column) {
  return(stats::setNames(extractionBoundary$entries[[column]], extractionBoundary$entries$case))
}

# The cohorts row(s) of a hand-made case.
boundaryCohort <- function(data, case) {
  cohorts <- collectCohorts(data)
  return(cohorts[cohorts$personId == as.character(boundaryPersonId(case)), ])
}

# The rows of a collected table for one rowId, as a data frame.
rowsOf <- function(table, rowId) {
  rows <- as.data.frame(collect(table))
  rows <- rows[rows$rowId == rowId, ]
  rownames(rows) <- NULL
  return(rows)
}

defaultArgs <- createGetDbCloneCensorWeightDataArgs(covariateSettings = extractionCovariateSettings())
defaultData <- extractQuietly(defaultArgs, connectionDetails = connectionDetails)
defaultExpected <- computeExpectedExtraction(rawTables, defaultArgs)

test_that("Extraction matches the reference with the default arguments", {
  expect_true(isCloneCensorWeightData(defaultData))
  expect_s4_class(defaultData, "CohortMethodData")
  expectMatchesReference(defaultData, defaultExpected)

  # The fixture exercises every step
  attrition <- defaultExpected$attrition
  expect_equal(attrition$description, c(
    "Original cohorts", "First cohort entry only", "365 days of prior observation",
    "No prior treatment use in the 365 days before index"
  ))
  expect_true(all(diff(attrition$targetExposures) < 0))
  expect_gt(nrow(defaultExpected$cohorts), 50)
  expect_gt(nrow(defaultExpected$treatments), 50)
  # Prior users are removed, so the remaining eras overlapping the lookback window start on day 0 or later
  expect_true(all(defaultExpected$treatments$daysToStart >= 0))
  expect_true(any(defaultExpected$outcomes$daysToEvent < 0))
  expect_true(any(defaultExpected$outcomes$daysToEvent > 0))
})

test_that("Extraction stores the metaData of the contract", {
  metaData <- attr(defaultData, "metaData")
  expect_equal(metaData$targetId, 101)
  expect_identical(metaData$treatmentIds, c(treatment = 102, comparatorTreatment = 103))
  expect_identical(metaData$outcomeIds, c(3, 1))
  expect_equal(metaData$populationSize, nrow(defaultExpected$cohorts))
  expect_equal(metaData$priorUseLookback, 365)
  expect_named(metaData$attrition, c("description", "targetPersons", "comparatorPersons", "targetExposures", "comparatorExposures"))
  expect_setequal(names(defaultData), c("cohorts", "outcomes", "treatments", "covariates", "covariateRef", "analysisRef"))
  # Column types as in the contract (checks NAs too)
  expect_silent(CloneCensorWeight:::.checkColumns(defaultData$cohorts, CloneCensorWeight:::.cohortsColumns, "cohorts", nullable = "daysToDeath"))
  expect_silent(CloneCensorWeight:::.checkColumns(defaultData$treatments, CloneCensorWeight:::.treatmentsColumns, "treatments"))
  expect_silent(CloneCensorWeight:::.checkColumns(defaultData$outcomes, CloneCensorWeight:::.outcomesColumns, "outcomes"))
  expect_silent(CloneCensorWeight:::.checkColumns(defaultData$covariates, CloneCensorWeight:::.covariatesColumns, "covariates"))
})

test_that("Death days are the earliest death record relative to index", {
  cohorts <- collectCohorts(defaultData)
  personId <- as.numeric(cohorts$personId)
  # Death on the first day of observation (person_id divisible by 30), else 100 days before the end of
  # observation (divisible by 40), else at the end of observation (divisible by 20), else none.
  expected <- rep(NA_integer_, nrow(cohorts))
  dies <- personId %% 20 == 0
  expected[dies] <- cohorts$daysToObsEnd[dies]
  dies <- personId %% 40 == 0
  expected[dies] <- cohorts$daysToObsEnd[dies] - 100L
  dies <- personId %% 30 == 0
  expected[dies] <- -cohorts$daysFromObsStart[dies]
  expect_equal(cohorts$daysToDeath, expected)
  expect_true(any(personId %% 40 == 0))
  expect_true(any(personId %% 30 == 0))
  expect_true(any(is.na(cohorts$daysToDeath)))
})

test_that("Covariates match the person table", {
  cohorts <- collectCohorts(defaultData)
  covariates <- as.data.frame(collect(defaultData$covariates))
  expect_true(all(covariates$rowId %in% cohorts$rowId))
  persons <- rawTables$persons[match(as.numeric(cohorts$personId), rawTables$persons$personId), ]
  expectedGender <- data.frame(
    rowId = cohorts$rowId,
    covariateId = persons$genderConceptId * 1000 + 1,
    covariateValue = 1
  )
  expectedAgeGroup <- data.frame(
    rowId = cohorts$rowId,
    covariateId = floor((as.numeric(format(cohorts$cohortStartDate, "%Y")) - persons$yearOfBirth) / 5) * 1000 + 3,
    covariateValue = 1
  )
  expect_equal(
    sortExtractionRows(covariates),
    sortExtractionRows(rbind(expectedGender, expectedAgeGroup))
  )
})

test_that("An open connection is used, left open, and its temp tables are dropped", {
  data <- extractQuietly(defaultArgs, connection = connection)
  on.exit(close(data))
  expectMatchesReference(data, defaultExpected)
  expect_equal(
    sortExtractionRows(collect(data$covariates)),
    sortExtractionRows(collect(defaultData$covariates))
  )
  expect_equal(DatabaseConnector::querySql(connection, "SELECT 1 AS x;")[[1]], 1)
  # sqlite_stat1 and sqlite_stat4 come from ANALYZE
  tempTablesSql <- "SELECT name FROM sqlite_temp_master WHERE name NOT LIKE 'sqlite%';"
  expect_equal(nrow(DatabaseConnector::querySql(connection, tempTablesSql)), 0)

  # A second extraction on the same connection
  args <- createGetDbCloneCensorWeightDataArgs(maxCohortSize = 20, covariateSettings = extractionCovariateSettings())
  sampled <- extractQuietly(args, connection = connection)
  on.exit(close(sampled), add = TRUE)
  expect_equal(nrow(collect(sampled$cohorts)), 20)
  expect_equal(nrow(DatabaseConnector::querySql(connection, tempTablesSql)), 0)
})

test_that("Prior use is an era overlapping days -priorUseLookback to -1", {
  # extractionEdgeEras: eras on days -400..-365 (prior use), -500..-366 (not), 0..10 (not), -1..5
  # (prior use), and -10..0 (prior use)
  cohorts <- collectCohorts(defaultData)
  present <- as.character(edgePersons) %in% cohorts$personId
  expect_equal(present, !extractionEdgeEras$priorUse)
  treatments <- collect(defaultData$treatments)
  rowIdDay0 <- cohorts$rowId[cohorts$personId == as.character(edgePersons[3])]
  expect_equal(nrow(dplyr::filter(treatments, .data$rowId == rowIdDay0, .data$daysToStart == 0, .data$daysToEnd == 10)), 1)

  args <- createGetDbCloneCensorWeightDataArgs(removeSubjectsWithPriorTreatment = FALSE, covariateSettings = extractionCovariateSettings())
  data <- extractQuietly(args, connectionDetails = connectionDetails)
  on.exit(close(data))
  expected <- computeExpectedExtraction(rawTables, args)
  expectMatchesReference(data, expected)
  expect_true(any(expected$treatments$daysToStart < 0))
  cohorts <- collectCohorts(data)
  expect_true(all(as.character(edgePersons) %in% cohorts$personId))
  treatments <- dplyr::inner_join(collect(data$treatments), cohorts, by = "rowId")
  edgeTreatments <- dplyr::filter(treatments, .data$personId %in% as.character(edgePersons), .data$treatmentId == 102)
  # The era ending on day -366 is outside the lookback window; the one ending on day -365 is inside
  expect_equal(sum(edgeTreatments$daysToStart == -400 & edgeTreatments$daysToEnd == -365), 1)
  expect_equal(sum(edgeTreatments$daysToStart == -500), 0)

  # With priorUseLookback = 0 nobody is removed for prior use
  args <- createGetDbCloneCensorWeightDataArgs(priorUseLookback = 0, covariateSettings = extractionCovariateSettings())
  dataNoLookback <- extractQuietly(args, connectionDetails = connectionDetails)
  on.exit(close(dataNoLookback), add = TRUE)
  expected <- computeExpectedExtraction(rawTables, args)
  expectMatchesReference(dataNoLookback, expected)
  expect_equal(nrow(attr(dataNoLookback, "metaData")$attrition), 3)
  expect_true(all(as.character(edgePersons) %in% collectCohorts(dataNoLookback)$personId))
  expect_true(all(collect(dataNoLookback$treatments)$daysToEnd >= 0))
})

test_that("Extraction matches the reference with all entries, no washout, and a study period", {
  args <- createGetDbCloneCensorWeightDataArgs(
    studyStartDate = "19800101",
    studyEndDate = "20101231",
    washoutPeriod = 0,
    firstOnly = FALSE,
    priorUseLookback = 180,
    covariateSettings = extractionCovariateSettings()
  )
  data <- extractQuietly(args, connectionDetails = connectionDetails)
  on.exit(close(data))
  expected <- computeExpectedExtraction(rawTables, args)
  expectMatchesReference(data, expected)
  expect_equal(expected$attrition$description, c(
    "Original cohorts", "0 days of prior observation", "Restrict to study period",
    "No prior treatment use in the 180 days before index"
  ))
  cohorts <- collectCohorts(data)
  expect_true(anyDuplicated(cohorts$personId) > 0)
  expect_true(all(cohorts$cohortStartDate >= as.Date("1980-01-01") & cohorts$cohortStartDate <= as.Date("2010-12-31")))
  # Observation end and cohort end are truncated at the study end date
  expect_true(all(cohorts$cohortStartDate + cohorts$daysToObsEnd <= as.Date("2010-12-31")))
  expect_true(any(cohorts$cohortStartDate + cohorts$daysToObsEnd == as.Date("2010-12-31")))
  expect_true(all(cohorts$cohortStartDate + cohorts$daysToCohortEnd <= as.Date("2010-12-31")))

  args <- createGetDbCloneCensorWeightDataArgs(studyStartDate = "20000101", covariateSettings = extractionCovariateSettings())
  dataStart <- extractQuietly(args, connectionDetails = connectionDetails)
  on.exit(close(dataStart), add = TRUE)
  expectMatchesReference(dataStart, computeExpectedExtraction(rawTables, args))
})

test_that("Hand-made entries on the washout, observation start, and study period boundaries are kept or dropped", {
  # Default arguments: 365 days of washout, no study period
  expect_equal(boundaryKept(defaultData), expectedBoundaryKept("keptDefault"))
  expect_equal(boundaryCohort(defaultData, "washoutExact")$daysFromObsStart, 365)
  # Without a study end date, the cohort end and the observation end are not truncated
  cohort <- boundaryCohort(defaultData, "cohortEndAfterStudyEnd")
  expect_equal(cohort$daysToCohortEnd, 273)
  expect_equal(cohort$daysToObsEnd, 2039)

  args <- createGetDbCloneCensorWeightDataArgs(
    studyStartDate = "19800101",
    studyEndDate = "20101231",
    washoutPeriod = 0,
    firstOnly = FALSE,
    priorUseLookback = 180,
    covariateSettings = extractionCovariateSettings()
  )
  data <- extractQuietly(args, connectionDetails = connectionDetails)
  on.exit(close(data))
  expect_equal(boundaryKept(data), expectedBoundaryKept("keptStudy"))
  expect_equal(boundaryCohort(data, "obsStart")$daysFromObsStart, 0)
  expect_equal(boundaryCohort(data, "studyStart")$cohortStartDate, as.Date("1980-01-01"))
  # The cohort end (2011-03-01) and the observation end are truncated at the study end date, 213 days
  # after the index date
  cohort <- boundaryCohort(data, "cohortEndAfterStudyEnd")
  expect_equal(cohort$daysToCohortEnd, 213)
  expect_equal(cohort$daysToObsEnd, 213)
  cohort <- boundaryCohort(data, "studyEnd")
  expect_equal(cohort$daysToCohortEnd, 0)
  expect_equal(cohort$daysToObsEnd, 0)
})

test_that("Outcomes on the first and last day of the index observation period are stored", {
  outcomeAtObsEnd <- boundaryCohort(defaultData, "outcomeAtObsEnd")
  outcomeAtObsStart <- boundaryCohort(defaultData, "outcomeAtObsStart")
  expect_equal(outcomeAtObsEnd$daysToObsEnd, 2191)
  expect_equal(outcomeAtObsStart$daysFromObsStart, 730)
  # Outcomes on the last day of observation and the day after: only the first is stored
  expect_equal(
    rowsOf(defaultData$outcomes, outcomeAtObsEnd$rowId),
    data.frame(rowId = outcomeAtObsEnd$rowId, outcomeId = 3, daysToEvent = 2191L)
  )
  # Outcomes on the first day of observation and the day before: only the first is stored
  expect_equal(
    rowsOf(defaultData$outcomes, outcomeAtObsStart$rowId),
    data.frame(rowId = outcomeAtObsStart$rowId, outcomeId = 3, daysToEvent = -730L)
  )

  # In CohortMethod's study population, the outcome on the last day of observation is counted on riskEnd,
  # and the outcome on the first day of observation is a prior outcome
  for (riskStart in c(0, 1)) {
    population <- suppressMessages(CohortMethod::createStudyPopulation(
      defaultData,
      outcomeId = 3,
      createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
        removeSubjectsWithPriorOutcome = TRUE,
        riskWindowStart = riskStart,
        riskWindowEnd = 99999
      )
    ))
    row <- population[population$rowId == outcomeAtObsEnd$rowId, ]
    expect_equal(nrow(row), 1)
    expect_equal(row$riskEnd, 2191)
    expect_equal(row$daysToEvent, 2191)
    expect_equal(row$outcomeCount, 1)
    expect_equal(row$survivalTime, 2191 - riskStart + 1)
    expect_false(outcomeAtObsStart$rowId %in% population$rowId)
  }
})

test_that("The index observation period is the one that contains the index date", {
  # secondObsPeriod: the index date is 100 days into the second of two observation periods; the first
  # holds a treatment era on days -251 to -176 and an outcome on day -220
  args <- createGetDbCloneCensorWeightDataArgs(
    washoutPeriod = 100,
    removeSubjectsWithPriorTreatment = FALSE,
    covariateSettings = extractionCovariateSettings()
  )
  data <- extractQuietly(args, connectionDetails = connectionDetails)
  on.exit(close(data))
  expectMatchesReference(data, computeExpectedExtraction(rawTables, args))
  expect_equal(boundaryKept(data), expectedBoundaryKept("keptWashout100"))
  cohort <- boundaryCohort(data, "secondObsPeriod")
  expect_equal(nrow(cohort), 1)
  expect_equal(cohort$daysFromObsStart, 100)
  expect_equal(cohort$daysToObsEnd, 2762)
  # The era in the first period overlaps the lookback window, so it is stored
  expect_equal(
    rowsOf(data$treatments, cohort$rowId),
    data.frame(rowId = cohort$rowId, treatmentId = 102, daysToStart = -251L, daysToEnd = -176L)
  )
  # The outcome in the first period is not stored; the one in the second period before index is
  expect_equal(
    rowsOf(data$outcomes, cohort$rowId),
    data.frame(rowId = cohort$rowId, outcomeId = 3, daysToEvent = -69L)
  )

  # The era in the first period is prior use
  args <- createGetDbCloneCensorWeightDataArgs(washoutPeriod = 100, covariateSettings = extractionCovariateSettings())
  dataPriorUse <- extractQuietly(args, connectionDetails = connectionDetails)
  on.exit(close(dataPriorUse), add = TRUE)
  expectMatchesReference(dataPriorUse, computeExpectedExtraction(rawTables, args))
  expected <- expectedBoundaryKept("keptWashout100")
  expected["secondObsPeriod"] <- FALSE
  expect_equal(boundaryKept(dataPriorUse), expected)
})

test_that("A single treatment role extracts only that treatment", {
  ids <- extractionIds
  ids$comparatorTreatmentId <- NULL
  args <- createGetDbCloneCensorWeightDataArgs(covariateSettings = extractionCovariateSettings())
  data <- extractQuietly(args, connectionDetails = connectionDetails, ids = ids)
  on.exit(close(data))
  expected <- computeExpectedExtraction(rawTables, args, treatmentIds = 102)
  expectMatchesReference(data, expected)
  expect_identical(attr(data, "metaData")$treatmentIds, c(treatment = 102))
  expect_equal(unique(collect(data$treatments)$treatmentId), 102)
  # Fewer persons are removed for prior use when only one treatment counts
  expect_gt(nrow(expected$cohorts), nrow(defaultExpected$cohorts))
})

test_that("maxCohortSize draws a random sample of the eligible cohort", {
  args <- createGetDbCloneCensorWeightDataArgs(maxCohortSize = 50, covariateSettings = extractionCovariateSettings())
  data <- extractQuietly(args, connectionDetails = connectionDetails)
  on.exit(close(data))
  cohorts <- collectCohorts(data)
  expect_equal(nrow(cohorts), 50)
  # Row IDs are those of the unsampled cohort; personSeqId is renumbered over the sample
  expectedCohorts <- defaultExpected$cohorts[match(cohorts$rowId, defaultExpected$cohorts$rowId), ]
  rownames(expectedCohorts) <- NULL
  columns <- setdiff(colnames(cohorts), "personSeqId")
  expect_equal(cohorts[, columns], expectedCohorts[, columns])
  expect_equal(cohorts$personSeqId, match(as.numeric(cohorts$personId), sort(unique(as.numeric(cohorts$personId)))))
  expect_equal(
    sortExtractionRows(collect(data$treatments)),
    sortExtractionRows(defaultExpected$treatments[defaultExpected$treatments$rowId %in% cohorts$rowId, ])
  )
  expect_equal(
    sortExtractionRows(collect(data$outcomes)),
    sortExtractionRows(defaultExpected$outcomes[defaultExpected$outcomes$rowId %in% cohorts$rowId, ])
  )
  expect_true(all(collect(data$covariates)$rowId %in% cohorts$rowId))
  attrition <- attr(data, "metaData")$attrition
  expect_equal(as.data.frame(attrition[seq_len(nrow(attrition) - 1), ]), defaultExpected$attrition)
  expect_equal(attrition$description[nrow(attrition)], "Random sample")
  expect_equal(attrition$targetExposures[nrow(attrition)], 50)
  expect_equal(attrition$targetPersons[nrow(attrition)], length(unique(cohorts$personId)))
  expect_equal(attr(data, "metaData")$populationSize, 50)
})

test_that("The Eunomia NSAID design (cohorts 4, 1, 2, 3) matches the reference", {
  ids <- list(targetId = 4, treatmentId = 1, comparatorTreatmentId = 2, outcomeIds = 3)
  args <- createGetDbCloneCensorWeightDataArgs(maxCohortSize = 150, covariateSettings = extractionCovariateSettings())
  data <- extractQuietly(args, connectionDetails = connectionDetails, ids = ids)
  on.exit(close(data))
  expected <- computeExpectedExtraction(fetchExtractionTables(connection, ids), args, treatmentIds = c(1, 2), outcomeIds = 3)
  cohorts <- collectCohorts(data)
  expect_equal(nrow(cohorts), 150)
  expectedCohorts <- expected$cohorts[match(cohorts$rowId, expected$cohorts$rowId), ]
  rownames(expectedCohorts) <- NULL
  columns <- setdiff(colnames(cohorts), "personSeqId")
  expect_equal(cohorts[, columns], expectedCohorts[, columns])
  expect_equal(
    sortExtractionRows(collect(data$treatments)),
    sortExtractionRows(expected$treatments[expected$treatments$rowId %in% cohorts$rowId, ])
  )
  expect_equal(
    sortExtractionRows(collect(data$outcomes)),
    sortExtractionRows(expected$outcomes[expected$outcomes$rowId %in% cohorts$rowId, ])
  )
})

test_that("The result feeds CohortMethod::createStudyPopulation", {
  for (riskStart in c(0, 1)) {
    population <- suppressMessages(CohortMethod::createStudyPopulation(
      defaultData,
      outcomeId = 3,
      createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
        removeSubjectsWithPriorOutcome = TRUE,
        riskWindowStart = riskStart,
        riskWindowEnd = 365
      )
    ))
    cohorts <- collectCohorts(defaultData)
    expect_true(all(c("daysToDeath", "riskStart", "riskEnd", "daysToEvent", "survivalTime") %in% colnames(population)))
    expect_true(all(population$rowId %in% cohorts$rowId))
    expect_true(all(population$riskStart == riskStart))
    idx <- match(population$rowId, cohorts$rowId)
    expect_equal(population$daysToDeath, cohorts$daysToDeath[idx])
    expect_equal(population$personSeqId, cohorts$personSeqId[idx])
    # First outcome on or after riskStart and within the risk window
    outcomes <- defaultExpected$outcomes
    outcomes <- outcomes[outcomes$outcomeId == 3 & outcomes$daysToEvent >= riskStart, ]
    firstOutcome <- stats::aggregate(daysToEvent ~ rowId, data = outcomes, FUN = min)
    expectedDays <- firstOutcome$daysToEvent[match(population$rowId, firstOutcome$rowId)]
    expectedDays[!is.na(expectedDays) & expectedDays > population$riskEnd] <- NA
    expect_equal(population$daysToEvent, expectedDays)
    attrition <- CohortMethod::getAttritionTable(population)
    expect_equal(as.data.frame(attrition[1:4, ]), defaultExpected$attrition)
  }
})

test_that("An empty eligible cohort gives empty tables with the right columns", {
  args <- createGetDbCloneCensorWeightDataArgs(studyStartDate = "20990101", covariateSettings = extractionCovariateSettings())
  warnings <- character(0)
  data <- withCallingHandlers(
    extractQuietly(args, connectionDetails = connectionDetails),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  on.exit(close(data))
  expect_true("Eligible cohort is empty" %in% warnings)
  cohorts <- collect(data$cohorts)
  expect_equal(nrow(cohorts), 0)
  expect_equal(nrow(collect(data$treatments)), 0)
  expect_equal(nrow(collect(data$outcomes)), 0)
  expect_equal(nrow(collect(data$covariates)), 0)
  expect_equal(
    vapply(cohorts, function(x) class(x)[1], character(1)),
    CloneCensorWeight:::.cohortsColumns
  )
  attrition <- attr(data, "metaData")$attrition
  expect_equal(attrition$targetExposures[nrow(attrition)], 0)
  expect_equal(attr(data, "metaData")$populationSize, 0)
})

test_that("Saved data loads as CloneCensorWeightData", {
  file <- tempfile(fileext = ".zip")
  on.exit(unlink(file))
  data <- extractQuietly(defaultArgs, connectionDetails = connectionDetails)
  suppressMessages(saveCloneCensorWeightData(data, file))
  loaded <- loadCloneCensorWeightData(file)
  on.exit(close(loaded), add = TRUE)
  expect_true(isCloneCensorWeightData(loaded))
  expect_equal(collectCohorts(loaded), defaultExpected$cohorts)
  expect_equal(attr(loaded, "metaData")$treatmentIds, c(treatment = 102, comparatorTreatment = 103))
})

close(defaultData)
DatabaseConnector::disconnect(connection)
unlink(connectionDetails$server())

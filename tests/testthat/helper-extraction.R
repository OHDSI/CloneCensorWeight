# Fixtures and the reference implementation for the extraction tests (test-extraction*.R).

# Cohort IDs of the Eunomia design built by createExtractionTestDatabase().
extractionIds <- list(
  targetId = 101, # all ankle sprains (concept 81151) of persons with person_id divisible by 10
  treatmentId = 102, # ibuprofen (1177480) drug eras, plus hand-made eras for prior-use edge cases
  comparatorTreatmentId = 103, # acetaminophen (1125315) drug eras
  outcomeIds = c(3, 1) # Eunomia's GI bleed and celecoxib cohorts (Eunomia::createCohorts)
)

# Offsets (days relative to the first index date) of hand-made treatment eras, one per person, and
# whether they are prior use when the lookback is 365 days: an era is prior use when it overlaps days
# -365 to -1.
extractionEdgeEras <- data.frame(
  startOffset = c(-400, -500, 0, -1, -10),
  endOffset = c(-365, -366, 10, 5, 0),
  priorUse = c(TRUE, FALSE, FALSE, TRUE, TRUE)
)

# Creates a copy of the Eunomia GiBleed database with Eunomia's cohorts 1-4, the cohorts in
# `extractionIds`, and death records (none in Eunomia): death on the last day of observation for
# persons with person_id divisible by 20, a second, earlier death record (100 days before the end of
# observation) for those divisible by 40, and a death on the first day of observation (before most
# index dates) for those divisible by 30. Hand-made treatment eras (extractionEdgeEras) are added for
# five persons whose first index date has 365 days of prior observation and no other prior use.
# Returns the connection details, with the persons of the edge cases as attribute "edgePersons".
createExtractionTestDatabase <- function() {
  databaseFile <- tempfile(fileext = ".sqlite")
  connectionDetails <- tryCatch(
    suppressMessages(Eunomia::getEunomiaConnectionDetails(databaseFile = databaseFile)),
    error = function(e) NULL
  )
  if (is.null(connectionDetails)) {
    testthat::skip("Eunomia GiBleed data could not be loaded")
  }
  suppressMessages(Eunomia::createCohorts(connectionDetails))
  connection <- DatabaseConnector::connect(connectionDetails)
  on.exit(DatabaseConnector::disconnect(connection))
  sql <- "
    INSERT INTO main.cohort (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
    SELECT @target_id, person_id, condition_start_date, MAX(COALESCE(condition_end_date, condition_start_date))
    FROM main.condition_occurrence
    WHERE condition_concept_id = 81151 AND person_id % 10 = 0
    GROUP BY person_id, condition_start_date;

    INSERT INTO main.cohort (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
    SELECT @treatment_id, person_id, drug_era_start_date, drug_era_end_date
    FROM main.drug_era
    WHERE drug_concept_id = 1177480;

    INSERT INTO main.cohort (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
    SELECT @comparator_treatment_id, person_id, drug_era_start_date, drug_era_end_date
    FROM main.drug_era
    WHERE drug_concept_id = 1125315;

    INSERT INTO main.death (person_id, death_date, death_type_concept_id)
    SELECT person_id, observation_period_end_date, 32817
    FROM main.observation_period
    WHERE person_id % 20 = 0;

    INSERT INTO main.death (person_id, death_date, death_type_concept_id)
    SELECT person_id, DATEADD(DAY, -100, observation_period_end_date), 32817
    FROM main.observation_period
    WHERE person_id % 40 = 0;

    INSERT INTO main.death (person_id, death_date, death_type_concept_id)
    SELECT person_id, observation_period_start_date, 32817
    FROM main.observation_period
    WHERE person_id % 30 = 0;"
  DatabaseConnector::renderTranslateExecuteSql(
    connection = connection,
    sql = sql,
    target_id = extractionIds$targetId,
    treatment_id = extractionIds$treatmentId,
    comparator_treatment_id = extractionIds$comparatorTreatmentId,
    progressBar = FALSE,
    reportOverallTime = FALSE
  )

  # Persons for the hand-made eras: first index date with 365 days of prior observation and no
  # treatment era overlapping days -365 to -1, and not dying on the first day of observation.
  tables <- fetchExtractionTables(connection)
  firstIndex <- tables$target[order(tables$target$subjectId, tables$target$cohortStartDate), ]
  firstIndex <- firstIndex[!duplicated(firstIndex$subjectId), ]
  firstIndex <- merge(firstIndex, tables$observationPeriods, by.x = "subjectId", by.y = "personId")
  firstIndex <- firstIndex[
    firstIndex$cohortStartDate >= firstIndex$observationPeriodStartDate + 365 &
      firstIndex$cohortStartDate <= firstIndex$observationPeriodEndDate &
      firstIndex$subjectId %% 30 != 0,
  ]
  priorUse <- merge(firstIndex, tables$exposures, by = "subjectId", suffixes = c("", "Exposure"))
  priorUse <- priorUse[
    priorUse$cohortStartDateExposure < priorUse$cohortStartDate &
      priorUse$cohortEndDateExposure >= priorUse$cohortStartDate - 365,
  ]
  candidates <- sort(setdiff(firstIndex$subjectId, priorUse$subjectId))
  edgePersons <- candidates[seq_len(nrow(extractionEdgeEras))]
  for (i in seq_along(edgePersons)) {
    DatabaseConnector::renderTranslateExecuteSql(
      connection = connection,
      sql = "
        INSERT INTO main.cohort (cohort_definition_id, subject_id, cohort_start_date, cohort_end_date)
        SELECT @treatment_id, subject_id, DATEADD(DAY, @start_offset, first_index), DATEADD(DAY, @end_offset, first_index)
        FROM (
          SELECT subject_id, MIN(cohort_start_date) AS first_index
          FROM main.cohort
          WHERE cohort_definition_id = @target_id AND subject_id = @subject_id
          GROUP BY subject_id
        ) tmp;",
      treatment_id = extractionIds$treatmentId,
      target_id = extractionIds$targetId,
      subject_id = edgePersons[i],
      start_offset = extractionEdgeEras$startOffset[i],
      end_offset = extractionEdgeEras$endOffset[i],
      progressBar = FALSE,
      reportOverallTime = FALSE
    )
  }
  attr(connectionDetails, "edgePersons") <- edgePersons
  return(connectionDetails)
}

# Fetches the raw tables the reference implementation needs, by direct SQL on the CDM.
fetchExtractionTables <- function(connection, ids = extractionIds) {
  query <- function(sql, ...) {
    result <- DatabaseConnector::renderTranslateQuerySql(connection, sql, snakeCaseToCamelCase = TRUE, ...)
    for (column in grep("Date$", colnames(result), value = TRUE)) {
      result[[column]] <- as.Date(result[[column]])
    }
    return(result)
  }
  return(list(
    target = query(
      "SELECT subject_id, cohort_start_date, cohort_end_date FROM main.cohort WHERE cohort_definition_id = @id;",
      id = ids$targetId
    ),
    exposures = query(
      "SELECT cohort_definition_id, subject_id, cohort_start_date, cohort_end_date FROM main.cohort WHERE cohort_definition_id IN (@ids);",
      ids = c(ids$treatmentId, ids$comparatorTreatmentId)
    ),
    outcomes = query(
      "SELECT cohort_definition_id, subject_id, cohort_start_date FROM main.cohort WHERE cohort_definition_id IN (@ids);",
      ids = ids$outcomeIds
    ),
    observationPeriods = query(
      "SELECT person_id, observation_period_start_date, observation_period_end_date FROM main.observation_period;"
    ),
    deaths = query("SELECT person_id, death_date FROM main.death;"),
    persons = query("SELECT person_id, gender_concept_id, year_of_birth FROM main.person;")
  ))
}

# Reference implementation of the extraction (in R, from the raw tables): returns the expected
# cohorts, treatments, outcomes and attrition (before sampling).
computeExpectedExtraction <- function(tables,
                                      args,
                                      treatmentIds = c(extractionIds$treatmentId, extractionIds$comparatorTreatmentId),
                                      outcomeIds = extractionIds$outcomeIds) {
  counts <- function(description, x) {
    data.frame(
      description = description,
      targetPersons = as.numeric(length(unique(x$subjectId))),
      comparatorPersons = 0,
      targetExposures = as.numeric(nrow(x)),
      comparatorExposures = 0
    )
  }
  lookback <- args$priorUseLookback
  cohorts <- tables$target
  attrition <- counts("Original cohorts", cohorts)
  if (args$firstOnly) {
    cohorts <- cohorts[order(cohorts$subjectId, cohorts$cohortStartDate), ]
    cohorts <- cohorts[!duplicated(cohorts$subjectId), ]
    attrition <- rbind(attrition, counts("First cohort entry only", cohorts))
  }
  cohorts <- merge(cohorts, tables$observationPeriods, by.x = "subjectId", by.y = "personId")
  cohorts <- cohorts[
    cohorts$cohortStartDate >= cohorts$observationPeriodStartDate &
      cohorts$cohortStartDate <= cohorts$observationPeriodEndDate &
      as.numeric(cohorts$cohortStartDate - cohorts$observationPeriodStartDate) >= args$washoutPeriod,
  ]
  attrition <- rbind(attrition, counts(sprintf("%d days of prior observation", args$washoutPeriod), cohorts))
  toDate <- function(x) as.Date(x, format = "%Y%m%d")
  if (args$studyStartDate != "" || args$studyEndDate != "") {
    if (args$studyStartDate != "") {
      cohorts <- cohorts[cohorts$cohortStartDate >= toDate(args$studyStartDate), ]
    }
    if (args$studyEndDate != "") {
      cohorts <- cohorts[cohorts$cohortStartDate <= toDate(args$studyEndDate), ]
    }
    attrition <- rbind(attrition, counts("Restrict to study period", cohorts))
  }
  exposures <- tables$exposures[tables$exposures$cohortDefinitionId %in% treatmentIds, ]
  if (args$removeSubjectsWithPriorTreatment && lookback > 0) {
    priorUse <- merge(cohorts, exposures, by = "subjectId", suffixes = c("", "Exposure"))
    priorUse <- priorUse[
      priorUse$cohortStartDateExposure <= priorUse$cohortStartDate - 1 &
        priorUse$cohortEndDateExposure >= priorUse$cohortStartDate - lookback,
    ]
    key <- paste(cohorts$subjectId, cohorts$cohortStartDate)
    cohorts <- cohorts[!key %in% paste(priorUse$subjectId, priorUse$cohortStartDate), ]
    attrition <- rbind(attrition, counts(sprintf("No prior treatment use in the %d days before index", lookback), cohorts))
  }

  cohorts <- cohorts[order(cohorts$subjectId, cohorts$cohortStartDate), ]
  cohortEnd <- cohorts$cohortEndDate
  observationEnd <- cohorts$observationPeriodEndDate
  if (args$studyEndDate != "") {
    cohortEnd <- pmin(cohortEnd, toDate(args$studyEndDate))
    observationEnd <- pmin(observationEnd, toDate(args$studyEndDate))
  }
  deaths <- stats::aggregate(deathDate ~ personId, data = tables$deaths, FUN = min)
  deathDate <- deaths$deathDate[match(cohorts$subjectId, deaths$personId)]
  expectedCohorts <- data.frame(
    rowId = seq_len(nrow(cohorts)),
    personSeqId = match(cohorts$subjectId, sort(unique(cohorts$subjectId))),
    personId = as.character(cohorts$subjectId),
    treatment = 1L,
    cohortStartDate = cohorts$cohortStartDate,
    daysFromObsStart = as.integer(cohorts$cohortStartDate - cohorts$observationPeriodStartDate),
    daysToCohortEnd = as.integer(cohortEnd - cohorts$cohortStartDate),
    daysToObsEnd = as.integer(observationEnd - cohorts$cohortStartDate),
    daysToDeath = as.integer(deathDate - cohorts$cohortStartDate)
  )
  index <- data.frame(
    rowId = expectedCohorts$rowId,
    subjectId = cohorts$subjectId,
    indexDate = cohorts$cohortStartDate,
    daysFromObsStart = expectedCohorts$daysFromObsStart,
    daysToObsEnd = expectedCohorts$daysToObsEnd
  )

  treatments <- merge(index, exposures, by = "subjectId")
  treatments$daysToStart <- as.integer(treatments$cohortStartDate - treatments$indexDate)
  treatments$daysToEnd <- as.integer(treatments$cohortEndDate - treatments$indexDate)
  treatments <- treatments[treatments$daysToStart <= treatments$daysToObsEnd & treatments$daysToEnd >= -lookback, ]
  expectedTreatments <- unique(data.frame(
    rowId = treatments$rowId,
    treatmentId = as.numeric(treatments$cohortDefinitionId),
    daysToStart = treatments$daysToStart,
    daysToEnd = treatments$daysToEnd
  ))

  outcomes <- merge(index, tables$outcomes[tables$outcomes$cohortDefinitionId %in% outcomeIds, ], by = "subjectId")
  outcomes$daysToEvent <- as.integer(outcomes$cohortStartDate - outcomes$indexDate)
  outcomes <- outcomes[-outcomes$daysToEvent <= outcomes$daysFromObsStart & outcomes$daysToEvent <= outcomes$daysToObsEnd, ]
  expectedOutcomes <- unique(data.frame(
    rowId = outcomes$rowId,
    outcomeId = as.numeric(outcomes$cohortDefinitionId),
    daysToEvent = outcomes$daysToEvent
  ))

  return(list(
    cohorts = expectedCohorts,
    treatments = expectedTreatments,
    outcomes = expectedOutcomes,
    attrition = attrition
  ))
}

# Sorts the rows of a data frame by all its columns and drops row names, for comparisons.
sortExtractionRows <- function(x) {
  x <- as.data.frame(x)
  x <- x[do.call(order, unname(as.list(x))), , drop = FALSE]
  rownames(x) <- NULL
  return(x)
}

# Small covariate settings for the extraction tests: gender and 5-year age group.
extractionCovariateSettings <- function() {
  return(FeatureExtraction::createCovariateSettings(
    useDemographicsGender = TRUE,
    useDemographicsAgeGroup = TRUE
  ))
}

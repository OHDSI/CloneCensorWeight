library(testthat)

# Treatment indicator for the view: started treatment 21 by day 7.
createViewCohorts <- function(tables) {
  startedBy7 <- unique(tables$treatments$rowId[tables$treatments$daysToStart <= 7])
  return(data.frame(
    rowId = tables$cohorts$rowId,
    personSeqId = tables$cohorts$personSeqId,
    treatment = as.integer(tables$cohorts$rowId %in% startedBy7)
  ))
}

sortCovariates <- function(covariates) {
  covariates <- as.data.frame(covariates)[, c("rowId", "covariateId", "covariateValue")]
  covariates <- covariates[order(covariates$rowId, covariates$covariateId), ]
  rownames(covariates) <- NULL
  return(covariates)
}

test_that(".createCmDataView without rowIdMap keeps the covariates of the given rows", {
  tables <- createFoundationTables()
  data <- do.call(createCloneCensorWeightData, tables)
  on.exit(close(data))
  cohorts <- createViewCohorts(tables)[seq(1, 200, by = 2), ]

  view <- .createCmDataView(data, cohorts)
  on.exit(close(view), add = TRUE)

  expect_true(is(view, "CohortMethodData"))
  expect_false(isCloneCensorWeightData(view))
  viewCohorts <- dplyr::collect(view$cohorts)
  expect_equal(viewCohorts$treatment, cohorts$treatment)
  expect_equal(viewCohorts$personId, tables$cohorts$personId[cohorts$rowId])
  expect_equal(viewCohorts$cohortStartDate, tables$cohorts$cohortStartDate[cohorts$rowId])
  expect_equal(viewCohorts$daysToObsEnd, tables$cohorts$daysToObsEnd[cohorts$rowId])
  expected <- tables$covariates[tables$covariates$rowId %in% cohorts$rowId, ]
  expect_equal(sortCovariates(dplyr::collect(view$covariates)), sortCovariates(expected))
  expect_equal(nrow(dplyr::collect(view$covariateRef)), 5)
  expect_equal(nrow(dplyr::collect(view$analysisRef)), 1)
  expect_equal(attr(view, "metaData")$populationSize, 100)
  expect_equal(sort(names(view)), c("analysisRef", "cohorts", "covariateRef", "covariates"))

  # The source is unchanged:
  expect_equal(sort(names(data)), sort(.cloneCensorWeightDataTables))
  expect_equal(attr(data, "metaData")$populationSize, 200)
})

test_that(".createCmDataView with rowIdMap remaps the covariates", {
  tables <- createFoundationTables()
  data <- do.call(createCloneCensorWeightData, tables)
  on.exit(close(data))
  # Two clone rows per person, with new row IDs in reverse order:
  n <- nrow(tables$cohorts)
  rowIdMap <- data.frame(newRowId = (2 * n):1, rowId = rep(tables$cohorts$rowId, 2))
  cohorts <- data.frame(
    rowId = rowIdMap$newRowId,
    personSeqId = rowIdMap$rowId,
    treatment = rep(c(1L, 0L), each = n),
    cohortStartDate = as.Date("2000-01-01")
  )

  view <- .createCmDataView(data, cohorts, rowIdMap)
  on.exit(close(view), add = TRUE)

  viewCohorts <- dplyr::collect(view$cohorts)
  expect_equal(viewCohorts$rowId, cohorts$rowId)
  expect_equal(viewCohorts$personId, tables$cohorts$personId[rowIdMap$rowId])
  # Given columns are kept, missing ones filled from the source row:
  expect_true(all(viewCohorts$cohortStartDate == as.Date("2000-01-01")))
  expect_equal(viewCohorts$daysFromObsStart, tables$cohorts$daysFromObsStart[rowIdMap$rowId])

  expected <- merge(tables$covariates, rowIdMap, by = "rowId")
  expected$rowId <- expected$newRowId
  covariates <- dplyr::collect(view$covariates)
  expect_equal(nrow(covariates), 2 * nrow(tables$covariates))
  expect_equal(sortCovariates(covariates), sortCovariates(expected))
  expect_type(covariates$rowId, "integer")
  expect_equal(attr(view, "metaData")$populationSize, 2 * n)

  # A subset of the map's rows:
  view2 <- .createCmDataView(data, cohorts[1:10, ], rowIdMap)
  on.exit(close(view2), add = TRUE)
  expect_equal(sort(unique(dplyr::pull(view2$covariates, "rowId"))), sort(intersect(cohorts$rowId[1:10], expected$rowId)))

  expect_error(.createCmDataView(data, cohorts, rowIdMap[-1, ]), "newRowId in rowIdMap")
  expect_error(.createCmDataView(data, cohorts[, c("rowId", "personSeqId")]), "treatment")
  expect_error(.createCmDataView(data, cohorts, data.frame(newRowId = 1, oldRowId = 1)), "rowId")
  expect_error(.createCmDataView(view, cohorts), "CloneCensorWeightData")
})

test_that("CohortMethod::createPs runs on a view", {
  tables <- createFoundationTables()
  data <- do.call(createCloneCensorWeightData, tables)
  on.exit(close(data))
  cohorts <- createViewCohorts(tables)
  view <- .createCmDataView(data, cohorts)
  on.exit(close(view), add = TRUE)

  ps <- suppressMessages(CohortMethod::createPs(
    cohortMethodData = view,
    createPsArgs = CohortMethod::createCreatePsArgs(prior = Cyclops::createPrior("none"))
  ))
  expect_equal(nrow(ps), nrow(cohorts))
  expect_true(all(ps$propensityScore > 0 & ps$propensityScore < 1))
  model <- attr(ps, "metaData")$psModelCoef
  # Covariate 1002 raises the probability of starting treatment:
  expect_gt(model[["1002"]], 0)

  # And with remapped rows (one clone per person):
  rowIdMap <- data.frame(newRowId = cohorts$rowId + 1000L, rowId = cohorts$rowId)
  cohorts2 <- cohorts
  cohorts2$rowId <- rowIdMap$newRowId
  view2 <- .createCmDataView(data, cohorts2, rowIdMap)
  on.exit(close(view2), add = TRUE)
  ps2 <- suppressMessages(CohortMethod::createPs(
    cohortMethodData = view2,
    createPsArgs = CohortMethod::createCreatePsArgs(prior = Cyclops::createPrior("none"))
  ))
  expect_equal(ps2$propensityScore[order(ps2$rowId)], ps$propensityScore[order(ps$rowId)], tolerance = 1e-8)
})

test_that(".andromedaSql and .andromedaExecute", {
  andromeda <- Andromeda::andromeda(x = data.frame(a = 1:3, b = c("p", "q", "r")))
  on.exit(close(andromeda))
  result <- .andromedaSql(andromeda, "SELECT SUM(a) AS total FROM x;")
  expect_equal(result$total, 6)
  .andromedaExecute(andromeda, "CREATE TABLE y AS SELECT a * 2 AS a2 FROM x WHERE b <> 'q';")
  expect_true("y" %in% names(andromeda))
  expect_equal(dplyr::pull(andromeda$y, "a2"), c(2, 6))
})

test_that(".assertCloneCensorWeightData", {
  data <- createFoundationData(n = 5)
  expect_silent(.assertCloneCensorWeightData(data))
  expect_error(.assertCloneCensorWeightData(data.frame()), "CloneCensorWeightData")
  collection <- checkmate::makeAssertCollection()
  .assertCloneCensorWeightData(list(), add = collection)
  expect_false(collection$isEmpty())
  data$treatments <- NULL
  expect_error(.assertCloneCensorWeightData(data), "treatments")
  close(data)
  expect_error(.assertCloneCensorWeightData(data), "closed")
})

test_that(".checkColumns coerces only when asked", {
  table <- data.frame(a = c(1, 2), b = 1:2)
  columnTypes <- c(a = "integer", b = "numeric")
  expect_error(.checkColumns(table, columnTypes, "t"), "'a' of table 't' must be of type integer")
  coerced <- .checkColumns(table, columnTypes, "t", coerce = TRUE)
  expect_type(coerced$a, "integer")
  expect_type(coerced$b, "double")
  expect_equal(colnames(.createEmptyTable(c(x = "Date", y = "character"))), c("x", "y"))

  # All-NA logical columns take the specified type:
  table <- data.frame(i = NA, n = NA, c = NA, d = NA)
  columnTypes <- c(i = "integer", n = "numeric", c = "character", d = "Date")
  expect_error(.checkColumns(table, columnTypes, "t", nullable = names(columnTypes)), "must be of type integer, but is logical")
  coerced <- .checkColumns(table, columnTypes, "t", nullable = names(columnTypes), coerce = TRUE)
  expect_type(coerced$i, "integer")
  expect_type(coerced$n, "double")
  expect_type(coerced$c, "character")
  expect_s3_class(coerced$d, "Date")
  expect_true(all(is.na(unlist(coerced))))
  # A logical column with non-NA values is still rejected:
  expect_error(.checkColumns(data.frame(i = c(TRUE, NA)), c(i = "integer"), "t", nullable = "i", coerce = TRUE), "must be of type integer")
})

test_that(".assertStudyPopulation", {
  data <- createFoundationData()
  on.exit(close(data))
  studyPop <- suppressMessages(CohortMethod::createStudyPopulation(
    cohortMethodData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 365, endAnchor = "cohort start")
  ))
  expect_silent(.assertStudyPopulation(studyPop))
  expect_error(.assertStudyPopulation(studyPop[, setdiff(colnames(studyPop), "daysToDeath")]), "daysToDeath")
  expect_error(.assertStudyPopulation(list(rowId = 1)), "data.frame")
  collection <- checkmate::makeAssertCollection()
  .assertStudyPopulation(studyPop[, c("rowId", "personSeqId")], add = collection)
  expect_false(collection$isEmpty())
})

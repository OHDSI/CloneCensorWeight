# Helpers ------------------------------------------------------------------

make_run_data <- function(n = 60, seed = 7) {
  createSampleCloneCensorWeightData(n = n, seed = seed)
}

make_analysis <- function(id = "main") {
  createCcwAnalysis(
    analysisId                 = id,
    description                = paste("Analysis", id),
    createStudyPopulationArgs  = createCreateStudyPopulationArgs(
      riskWindowEnd = 365, endAnchor = "cohort start"
    ),
    createClonesAndCensorArgs  = createCreateClonesAndCensorArgs(
      arms = list(
        createArmDefinition("initiate", "initiateWithin", gracePeriod = 30),
        createArmDefinition("never",    "deferBeyond",    gracePeriod = 30)
      ),
      intervalLength = 30
    ),
    createCensoringWeightsArgs = createCreateCensoringWeightsArgs(),
    fitOutcomeModelArgs        = createFitOutcomeModelArgs("cox", "robust")
  )
}

make_tao <- function() {
  createTargetComparatorOutcome(
    indexCohortId      = 1L,
    targetCohortId  = 2L,
    comparatorCohortId = 2L,
    outcomeIds         = c(10L, 11L, 12L)
  )
}

# Tests -------------------------------------------------------------------

test_that(".makeCompositeId produces correct string", {
  id <- CloneCensorWeight:::.makeCompositeId("main", 1L, 10L)
  expect_equal(id, "main_t1_o10")
})

test_that(".alreadyRun returns FALSE when folder is empty", {
  tmp <- withr::local_tempdir()
  expect_false(CloneCensorWeight:::.alreadyRun("main_t1_o10", tmp))
})

test_that(".alreadyRun returns TRUE after result is exported", {
  tmp <- withr::local_tempdir()
  result <- list(
    outcomeModel = data.frame(modelType = "cox", term = "arm",
                              estimate = 1.1, lower95 = 0.8, upper95 = 1.5,
                              pValue = 0.3, logLikelihood = -90,
                              stringsAsFactors = FALSE)
  )
  exportToCsv(result, "main_t1_o10", tmp, append = FALSE)
  expect_true(CloneCensorWeight:::.alreadyRun("main_t1_o10", tmp))
  expect_false(CloneCensorWeight:::.alreadyRun("main_t1_o99", tmp))
})

test_that(".buildNcDf returns NULL when no NC ids match results", {
  res <- list(
    "10" = list(model = list(summary = data.frame(
      estimate = 1.0, lower95 = 0.7, upper95 = 1.4,
      stringsAsFactors = FALSE)))
  )
  out <- CloneCensorWeight:::.buildNcDf(res, integer(0))
  expect_null(out)
})

test_that(".buildNcDf extracts logRr and seLogRr for matching NC ids", {
  makeMod <- function(est, lo, hi) {
    structure(list(summary = data.frame(
      modelType = "cox", term = "arm",
      estimate = est, lower95 = lo, upper95 = hi,
      pValue = 0.5, logLikelihood = -80, stringsAsFactors = FALSE
    )), class = "CloneCensorWeightOutcomeModel")
  }
  res <- list(
    "11" = list(model = makeMod(1.0, 0.75, 1.33)),
    "12" = list(model = makeMod(0.9, 0.65, 1.25))
  )
  out <- CloneCensorWeight:::.buildNcDf(res, c(11L, 12L))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)
  expect_true(all(c("logRr", "seLogRr") %in% names(out)))
  expect_true(all(is.finite(out$logRr)))
  expect_true(all(out$seLogRr > 0))
})

test_that(".runForData writes CSV files for each outcome", {
  tmp      <- withr::local_tempdir()
  data     <- make_run_data()
  analysis <- make_analysis("a1")
  tao      <- make_tao()

  CloneCensorWeight:::.runForData(
    data         = data,
    tao          = tao,
    analysis     = analysis,
    outputFolder = tmp,
    databaseId   = "synpuf"
  )

  files <- list.files(tmp)
  expect_true("ccw_outcome_model.csv"    %in% files)
  expect_true("ccw_kaplan_meier.csv"     %in% files)
  expect_true("ccw_power_mdrr.csv"       %in% files)
})

test_that(".runForData stamps analysis_id correctly", {
  tmp      <- withr::local_tempdir()
  data     <- make_run_data()
  analysis <- make_analysis("beta")
  tao      <- make_tao()

  CloneCensorWeight:::.runForData(
    data         = data,
    tao          = tao,
    analysis     = analysis,
    outputFolder = tmp,
    databaseId   = "synpuf"
  )

  om <- readr::read_csv(file.path(tmp, "ccw_outcome_model.csv"), show_col_types = FALSE)
  expect_true(all(grepl("^beta_t1_o", om$analysis_id)))
  expect_equal(nrow(om), length(tao$outcomeIds))
})

test_that(".runForData skips already-completed analyses", {
  tmp      <- withr::local_tempdir()
  data     <- make_run_data()
  analysis <- make_analysis("skip_test")
  tao      <- make_tao()

  # First run
  CloneCensorWeight:::.runForData(data, tao, analysis, tmp, "synpuf")
  om1 <- readr::read_csv(file.path(tmp, "ccw_outcome_model.csv"), show_col_types = FALSE)

  # Second run with same ids — should not duplicate
  CloneCensorWeight:::.runForData(data, tao, analysis, tmp, "synpuf")
  om2 <- readr::read_csv(file.path(tmp, "ccw_outcome_model.csv"), show_col_types = FALSE)

  expect_equal(nrow(om1), nrow(om2))
})

test_that(".runForData accumulates results from multiple analyses", {
  tmp       <- withr::local_tempdir()
  data      <- make_run_data()
  analysis1 <- make_analysis("sens1")
  analysis2 <- make_analysis("sens2")
  tao       <- make_tao()

  CloneCensorWeight:::.runForData(data, tao, analysis1, tmp, "db1")
  CloneCensorWeight:::.runForData(data, tao, analysis2, tmp, "db1")

  om <- readr::read_csv(file.path(tmp, "ccw_outcome_model.csv"), show_col_types = FALSE)
  expect_equal(nrow(om), 2L * length(tao$outcomeIds))
  expect_true("sens1_t1_o10" %in% om$analysis_id)
  expect_true("sens2_t1_o10" %in% om$analysis_id)
})

test_that("getResultsDataModelSpec returns a data frame with required columns", {
  spec <- getResultsDataModelSpec()
  expect_s3_class(spec, "data.frame")
  expect_true(all(c("tableName", "columnName", "dataType", "isRequired", "primaryKey")
                  %in% names(spec)))
  expect_gt(nrow(spec), 0)
})

test_that("spec covers all expected tables", {
  spec   <- getResultsDataModelSpec()
  tables <- unique(spec$tableName)
  expected <- c(
    "ccw_analysis", "ccw_cohort_definition", "ccw_cohort_count",
    "ccw_outcome_model", "ccw_standardized_difference",
    "ccw_kaplan_meier", "ccw_power_mdrr", "ccw_negative_control"
  )
  expect_true(all(expected %in% tables))
})

test_that("every table has at least one primary key column", {
  spec <- getResultsDataModelSpec()
  pkCounts <- tapply(spec$primaryKey, spec$tableName, function(x) sum(x == TRUE))
  expect_true(all(pkCounts >= 1))
})

test_that(".applyRenameMap renames camelCase to snake_case correctly", {
  df  <- data.frame(modelType = "cox", lower95 = 0.8, extra = 1)
  map <- c(modelType = "model_type", lower95 = "lower_95")
  out <- CloneCensorWeight:::.applyRenameMap(df, map)
  expect_true("model_type" %in% names(out))
  expect_true("lower_95"   %in% names(out))
  expect_true("extra"      %in% names(out))
  expect_false("modelType" %in% names(out))
})

test_that(".renameOutcomeModelCols maps all summary columns", {
  df <- data.frame(
    modelType = "cox", term = "arm", estimate = 1.2,
    lower95 = 0.9, upper95 = 1.6, pValue = 0.05,
    logLikelihood = -100, stringsAsFactors = FALSE
  )
  out <- CloneCensorWeight:::.renameOutcomeModelCols(df)
  expect_true("model_type"    %in% names(out))
  expect_true("lower_95"      %in% names(out))
  expect_true("upper_95"      %in% names(out))
  expect_true("p_value"       %in% names(out))
  expect_true("log_likelihood" %in% names(out))
})

test_that(".renameSmdCols maps all SMD columns", {
  df <- data.frame(
    covariateId = 1L, covariate = "age", weighted = 0L,
    meanArm1 = 0.5, sdArm1 = 0.1, meanArm2 = 0.6, sdArm2 = 0.1,
    smd = 0.2, absSmd = 0.2, stringsAsFactors = FALSE
  )
  out <- CloneCensorWeight:::.renameSmdCols(df)
  expect_true("covariate_id" %in% names(out))
  expect_true("mean_arm1"    %in% names(out))
  expect_true("abs_smd"      %in% names(out))
})

test_that(".renameKmCols maps nRisk and nEvent", {
  df <- data.frame(arm = "a", time = 10, survival = 0.9,
                   lower = 0.8, upper = 1.0, nRisk = 50, nEvent = 2)
  out <- CloneCensorWeight:::.renameKmCols(df)
  expect_true("n_risk"  %in% names(out))
  expect_true("n_event" %in% names(out))
})

test_that(".renameMdrrCols maps all MDRR columns", {
  df <- data.frame(arm1 = "a", arm2 = "b", alpha = 0.05, power = 0.8,
                   nEventsArm1 = 10L, nEventsArm2 = 12L,
                   nEventsEffArm1 = 9.5, nEventsEffArm2 = 11.0, mdrr = 1.8)
  out <- CloneCensorWeight:::.renameMdrrCols(df)
  expect_true("n_events_arm1"      %in% names(out))
  expect_true("n_events_eff_arm1"  %in% names(out))
  expect_true("mdrr"               %in% names(out))
})

test_that("CSV is internally consistent: PK columns are all isRequired", {
  spec <- getResultsDataModelSpec()
  pkRows <- spec[spec$primaryKey == TRUE, ]
  expect_true(all(pkRows$isRequired == TRUE))
})

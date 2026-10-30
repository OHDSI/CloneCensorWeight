# Helpers ------------------------------------------------------------------

make_minimal_results <- function(analysis_id = "test_01") {
  list(
    analysis = data.frame(
      database_id  = "synpuf",
      database_name = "SynPUF",
      created_date  = Sys.Date(),
      stringsAsFactors = FALSE
    ),
    outcomeModel = data.frame(
      modelType = "cox", term = "arm",
      estimate = 1.2, lower95 = 0.9, upper95 = 1.6,
      pValue = 0.12, logLikelihood = -100,
      stringsAsFactors = FALSE
    ),
    negativeControl = data.frame(
      logRr = c(-0.1, 0.05), seLogRr = c(0.2, 0.15),
      stringsAsFactors = FALSE
    )
  )
}

# Tests -------------------------------------------------------------------

test_that("exportToCsv creates the export folder if absent", {
  tmp <- file.path(tempdir(), paste0("ccw_export_", Sys.getpid()))
  on.exit(unlink(tmp, recursive = TRUE))
  exportToCsv(make_minimal_results(), "a01", tmp)
  expect_true(dir.exists(tmp))
})

test_that("exportToCsv writes expected CSV files", {
  tmp <- withr::local_tempdir()
  exportToCsv(make_minimal_results(), "a01", tmp)
  files <- list.files(tmp)
  expect_true("ccw_analysis.csv"         %in% files)
  expect_true("ccw_outcome_model.csv"    %in% files)
  expect_true("ccw_negative_control.csv" %in% files)
})

test_that("exportToCsv stamps analysis_id onto every file", {
  tmp <- withr::local_tempdir()
  exportToCsv(make_minimal_results(), "RUN42", tmp)
  om  <- readr::read_csv(file.path(tmp, "ccw_outcome_model.csv"), show_col_types = FALSE)
  nc  <- readr::read_csv(file.path(tmp, "ccw_negative_control.csv"), show_col_types = FALSE)
  expect_true(all(om$analysis_id == "RUN42"))
  expect_true(all(nc$analysis_id == "RUN42"))
})

test_that("exportToCsv renames camelCase outcome model columns to snake_case", {
  tmp <- withr::local_tempdir()
  exportToCsv(make_minimal_results(), "a01", tmp)
  om <- readr::read_csv(file.path(tmp, "ccw_outcome_model.csv"), show_col_types = FALSE)
  expect_true("model_type"    %in% names(om))
  expect_true("lower_95"      %in% names(om))
  expect_true("p_value"       %in% names(om))
  expect_false("modelType"    %in% names(om))
  expect_false("lower95"      %in% names(om))
})

test_that("exportToCsv auto-assigns control_id when absent", {
  tmp <- withr::local_tempdir()
  exportToCsv(make_minimal_results(), "a01", tmp)
  nc <- readr::read_csv(file.path(tmp, "ccw_negative_control.csv"), show_col_types = FALSE)
  expect_true("control_id" %in% names(nc))
  expect_equal(nrow(nc), 2L)
})

test_that("exportToCsv returns paths of written files invisibly", {
  tmp <- withr::local_tempdir()
  paths <- exportToCsv(make_minimal_results(), "a01", tmp)
  expect_type(paths, "character")
  expect_true(length(paths) >= 3)
  expect_true(all(file.exists(paths)))
})

test_that("append = TRUE accumulates rows across calls", {
  tmp <- withr::local_tempdir()
  exportToCsv(make_minimal_results(), "run1", tmp, append = FALSE)
  exportToCsv(make_minimal_results(), "run2", tmp, append = TRUE)
  nc <- readr::read_csv(file.path(tmp, "ccw_negative_control.csv"), show_col_types = FALSE)
  expect_equal(nrow(nc), 4L)
  expect_true("run1" %in% nc$analysis_id)
  expect_true("run2" %in% nc$analysis_id)
})

test_that("append = FALSE overwrites existing file", {
  tmp <- withr::local_tempdir()
  exportToCsv(make_minimal_results(), "run1", tmp, append = FALSE)
  exportToCsv(make_minimal_results(), "run2", tmp, append = FALSE)
  nc <- readr::read_csv(file.path(tmp, "ccw_negative_control.csv"), show_col_types = FALSE)
  expect_equal(nrow(nc), 2L)
  expect_false("run1" %in% nc$analysis_id)
})

test_that("importFromCsv reads back what exportToCsv wrote", {
  tmp <- withr::local_tempdir()
  exportToCsv(make_minimal_results(), "a01", tmp)
  imported <- importFromCsv(tmp)
  expect_type(imported, "list")
  expect_true("ccw_outcome_model"    %in% names(imported))
  expect_true("ccw_negative_control" %in% names(imported))
  expect_s3_class(imported$ccw_outcome_model, "data.frame")
})

test_that("importFromCsv returns empty list when folder has no ccw_ CSVs", {
  tmp <- withr::local_tempdir()
  result <- importFromCsv(tmp)
  expect_equal(length(result), 0L)
})

test_that("column order in exported CSV matches schema spec for outcome_model", {
  tmp <- withr::local_tempdir()
  exportToCsv(make_minimal_results(), "a01", tmp)
  om   <- readr::read_csv(file.path(tmp, "ccw_outcome_model.csv"), show_col_types = FALSE)
  spec <- getResultsDataModelSpec()
  specCols <- spec$columnName[spec$tableName == "ccw_outcome_model"]
  presentSpec <- intersect(specCols, names(om))
  expect_equal(names(om)[seq_along(presentSpec)], presentSpec)
})

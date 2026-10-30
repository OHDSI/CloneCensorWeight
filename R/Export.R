#' Export analysis results to CSV files
#'
#' Writes each element of a results list to a CSV file named after its
#' corresponding `ccw_*` results table.  An `analysis_id` column is stamped
#' onto every row.  Files are written with column order matching the results
#' data model specification so that CSVs from multiple databases can be safely
#' appended and later bulk-loaded into a database.
#'
#' The `results` list accepts the same named elements as [uploadResults()]:
#' `analysis`, `cohortDefinition`, `cohortCount`, `outcomeModel`,
#' `standardizedDifference`, `kaplanMeier`, `powerMdrr`, `negativeControl`.
#' Unrecognised or `NULL` elements are silently skipped.
#'
#' @param results A named list of result data frames (see [uploadResults()]).
#' @param analysisId A character string that uniquely identifies this analysis
#'   run and is stamped onto every exported row.
#' @param exportFolder Path to the directory that will receive the CSV files.
#'   Created automatically if it does not exist.
#' @param append Logical.  When `TRUE`, rows are appended to any existing CSV
#'   file of the same name (useful when accumulating results from multiple
#'   databases in the same folder).  When `FALSE` (default), existing files
#'   are overwritten.
#'
#' @return Invisibly, a character vector of the file paths that were written.
#' @export
exportToCsv <- function(results, analysisId, exportFolder, append = FALSE) {
  checkmate::assertList(results, names = "unique")
  checkmate::assertString(analysisId, min.chars = 1)
  checkmate::assertString(exportFolder)
  checkmate::assertFlag(append)

  if (!dir.exists(exportFolder)) {
    dir.create(exportFolder, recursive = TRUE)
  }

  spec     <- getResultsDataModelSpec()
  written  <- character(0)

  .write <- function(tableName, df) {
    df$analysis_id <- analysisId
    df <- .orderToSpec(df, tableName, spec)
    path <- file.path(exportFolder, paste0(tableName, ".csv"))
    if (append && file.exists(path)) {
      readr::write_csv(df, path, append = TRUE)
    } else {
      readr::write_csv(df, path, append = FALSE)
    }
    written <<- c(written, path)
  }

  if (!is.null(results$analysis)) {
    checkmate::assertDataFrame(results$analysis)
    .write("ccw_analysis", results$analysis)
  }

  if (!is.null(results$cohortDefinition)) {
    checkmate::assertDataFrame(results$cohortDefinition)
    .write("ccw_cohort_definition", results$cohortDefinition)
  }

  if (!is.null(results$cohortCount)) {
    checkmate::assertDataFrame(results$cohortCount)
    .write("ccw_cohort_count", results$cohortCount)
  }

  if (!is.null(results$outcomeModel)) {
    checkmate::assertDataFrame(results$outcomeModel)
    df <- .renameOutcomeModelCols(results$outcomeModel)
    .write("ccw_outcome_model", df)
  }

  if (!is.null(results$standardizedDifference)) {
    checkmate::assertDataFrame(results$standardizedDifference)
    df <- .renameSmdCols(results$standardizedDifference)
    .write("ccw_standardized_difference", df)
  }

  if (!is.null(results$kaplanMeier)) {
    checkmate::assertDataFrame(results$kaplanMeier)
    df <- .renameKmCols(results$kaplanMeier)
    .write("ccw_kaplan_meier", df)
  }

  if (!is.null(results$powerMdrr)) {
    checkmate::assertDataFrame(results$powerMdrr)
    df <- .renameMdrrCols(results$powerMdrr)
    .write("ccw_power_mdrr", df)
  }

  if (!is.null(results$negativeControl)) {
    checkmate::assertDataFrame(results$negativeControl)
    df <- results$negativeControl
    if (!"control_id" %in% names(df) && !"controlId" %in% names(df)) {
      df$control_id <- as.character(seq_len(nrow(df)))
    }
    df <- .renameNcCols(df)
    .write("ccw_negative_control", df)
  }

  invisible(written)
}

#' Import previously exported CSV results
#'
#' Reads all `ccw_*.csv` files found in `exportFolder` and returns them as a
#' named list.  The list names are the bare table names without the `.csv`
#' suffix (e.g., `"ccw_outcome_model"`).
#'
#' @param exportFolder Path to the directory produced by [exportToCsv()].
#'
#' @return A named list of data frames, one per CSV file found.  Returns an
#'   empty list when no matching files exist.
#' @export
importFromCsv <- function(exportFolder) {
  checkmate::assertDirectoryExists(exportFolder)
  files <- list.files(exportFolder, pattern = "^ccw_.*\\.csv$", full.names = TRUE)
  if (length(files) == 0) return(list())

  out <- lapply(files, function(f) {
    readr::read_csv(f, show_col_types = FALSE)
  })
  names(out) <- sub("\\.csv$", "", basename(files))
  out
}

# Reorder a data frame's columns to match the schema spec for tableName.
# Columns present in the spec come first in spec order; extra columns follow.
.orderToSpec <- function(df, tableName, spec) {
  specCols <- spec$columnName[spec$tableName == tableName]
  present  <- intersect(specCols, names(df))
  extra    <- setdiff(names(df), specCols)
  df[, c(present, extra), drop = FALSE]
}

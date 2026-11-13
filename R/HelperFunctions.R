# Copyright 2026 Observational Health Data Sciences and Informatics
#
# This file is part of CloneCensorWeight
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Copied from CohortMethod (R/HelperFunctions.R), which does not export it.
.assertCovariateId <- function(covariateId, len = NULL, min.len = NULL, null.ok = FALSE, add = NULL) {
  checkmate::assertNumeric(covariateId, null.ok = null.ok, len = len, min.len = 1, add = add)
  if (!is.null(covariateId)) {
    message <- sprintf(
      "Variable '%s' is a (64-bit) integer",
      paste0(deparse(eval.parent(substitute(substitute(covariateId))), width.cutoff = 500L), collapse = "\n")
    )
    checkmate::assertTRUE(all(covariateId == round(covariateId)), .var.name = message, add = add)
  }
}

# Table validation ---------------------------------------------------------------------------------

.emptyColumn <- list(
  integer = integer(0),
  numeric = numeric(0),
  character = character(0),
  Date = as.Date(character(0))
)

# Creates a zero-row tibble with the given columns. `columnTypes` is a named character vector:
# names are column names, values are "integer", "numeric", "character", or "Date".
.createEmptyTable <- function(columnTypes) {
  return(dplyr::as_tibble(lapply(columnTypes, function(type) .emptyColumn[[type]])))
}

# Checks that `table` (a data frame or an Andromeda table) has the columns in `columnTypes` with the
# right types, and no NA values outside the `nullable` columns. With `coerce = TRUE` (data frames
# only), integer-valued numeric columns are converted to integer, integer columns to numeric, and
# all-NA logical columns (such as `daysToDeath = NA`) to the type the specification asks for.
# Returns the (possibly coerced) table.
.checkColumns <- function(table, columnTypes, tableName, nullable = character(0), coerce = FALSE) {
  isLazy <- inherits(table, "tbl_dbi")
  prototype <- if (isLazy) collect(utils::head(table, 0)) else table
  missingColumns <- setdiff(names(columnTypes), colnames(prototype))
  if (length(missingColumns) > 0) {
    stop(sprintf("Table '%s' is missing column(s): %s", tableName, paste(missingColumns, collapse = ", ")))
  }
  for (column in names(columnTypes)) {
    type <- columnTypes[[column]]
    values <- prototype[[column]]
    if (coerce && !isLazy) {
      if (is.logical(values) && all(is.na(values))) {
        values <- switch(type,
          integer = as.integer(values),
          numeric = as.numeric(values),
          character = as.character(values),
          Date = as.Date(values)
        )
      } else if (type == "integer" && !is.integer(values) && is.numeric(values) && checkmate::testIntegerish(values)) {
        values <- as.integer(values)
      } else if (type == "numeric" && is.integer(values)) {
        values <- as.numeric(values)
      }
      table[[column]] <- values
    }
    typeOk <- switch(type,
      integer = is.integer(values),
      numeric = is.numeric(values),
      character = is.character(values),
      Date = inherits(values, "Date")
    )
    if (!typeOk) {
      stop(sprintf("Column '%s' of table '%s' must be of type %s, but is %s", column, tableName, type, class(values)[1]))
    }
    if (!column %in% nullable) {
      if (isLazy) {
        hasNa <- table |>
          filter(is.na(.data[[column]])) |>
          count() |>
          pull() > 0
      } else {
        hasNa <- anyNA(values)
      }
      if (hasNa) {
        stop(sprintf("Column '%s' of table '%s' must not contain NA values", column, tableName))
      }
    }
  }
  return(table)
}

# Andromeda SQL ------------------------------------------------------------------------------------

# An Andromeda object is a DBI connection to its DuckDB database, so raw DuckDB SQL can be run on it
# directly. Table and column names are case-insensitive in DuckDB.
.andromedaSql <- function(andromeda, sql) {
  return(dplyr::as_tibble(DBI::dbGetQuery(andromeda, sql)))
}

.andromedaExecute <- function(andromeda, sql) {
  DBI::dbExecute(andromeda, sql)
  invisible(NULL)
}

# CohortMethodData view ----------------------------------------------------------------------------

.cmCohortColumns <- c(
  "rowId", "personSeqId", "personId", "treatment", "cohortStartDate", "daysFromObsStart",
  "daysToCohortEnd", "daysToObsEnd"
)

# Returns a new CohortMethodData whose cohorts are `cohorts` and whose covariates are the source
# covariates of those rows. With `rowIdMap` (columns newRowId and rowId), `cohorts$rowId` holds new row
# IDs, and each new row gets the covariates of the source row it maps to. CohortMethod columns missing
# from `cohorts` are filled from the source row. The caller must close the result.
.createCmDataView <- function(cloneCensorWeightData, cohorts, rowIdMap = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCloneCensorWeightData(cloneCensorWeightData, add = errorMessages)
  checkmate::assertDataFrame(cohorts, add = errorMessages)
  checkmate::assertNames(colnames(cohorts), must.include = c("rowId", "personSeqId", "treatment"), add = errorMessages)
  checkmate::assertIntegerish(cohorts$rowId, any.missing = FALSE, unique = TRUE, add = errorMessages)
  checkmate::assertDataFrame(rowIdMap, null.ok = TRUE, add = errorMessages)
  if (!is.null(rowIdMap)) {
    checkmate::assertNames(colnames(rowIdMap), must.include = c("newRowId", "rowId"), add = errorMessages)
    checkmate::assertIntegerish(rowIdMap$newRowId, any.missing = FALSE, unique = TRUE, add = errorMessages)
  }
  checkmate::reportAssertions(collection = errorMessages)

  if (is.null(rowIdMap)) {
    sourceRowIds <- cohorts$rowId
  } else {
    sourceRowIds <- rowIdMap$rowId[match(cohorts$rowId, rowIdMap$newRowId)]
    if (anyNA(sourceRowIds)) {
      stop("Every rowId in cohorts must appear as newRowId in rowIdMap")
    }
  }
  missingColumns <- setdiff(.cmCohortColumns, colnames(cohorts))
  if (length(missingColumns) > 0) {
    sourceCohorts <- cloneCensorWeightData$cohorts |>
      select(dplyr::all_of(c("rowId", missingColumns))) |>
      collect()
    idx <- match(sourceRowIds, sourceCohorts$rowId)
    if (anyNA(idx)) {
      stop("Some rows of cohorts have no matching row in the source cohorts")
    }
    cohorts[missingColumns] <- sourceCohorts[idx, missingColumns]
  }

  view <- Andromeda::andromeda(
    cohorts = cohorts,
    rowIdMap = tibble(newRowId = cohorts$rowId, rowId = sourceRowIds)
  )
  view$sourceCovariates <- cloneCensorWeightData$covariates
  .andromedaExecute(view, "
    CREATE TABLE covariates AS
    SELECT rowIdMap.newRowId AS rowId,
      sourceCovariates.covariateId,
      sourceCovariates.covariateValue
    FROM sourceCovariates
    INNER JOIN rowIdMap
      ON sourceCovariates.rowId = rowIdMap.rowId;")
  view$sourceCovariates <- NULL
  view$rowIdMap <- NULL
  view$covariateRef <- cloneCensorWeightData$covariateRef
  view$analysisRef <- cloneCensorWeightData$analysisRef

  metaData <- attr(cloneCensorWeightData, "metaData")
  metaData$populationSize <- nrow(cohorts)
  attr(view, "metaData") <- metaData
  class(view) <- "CohortMethodData"
  attr(class(view), "package") <- "CohortMethod"
  return(view)
}

# Assertions ---------------------------------------------------------------------------------------

.checkCloneCensorWeightData <- function(x) {
  if (!inherits(x, "CloneCensorWeightData")) {
    return("Must be an object of class 'CloneCensorWeightData'")
  }
  if (!Andromeda::isValidAndromeda(x)) {
    return("Must be an open Andromeda object, but it was closed")
  }
  missingTables <- setdiff(.cloneCensorWeightDataTables, names(x))
  if (length(missingTables) > 0) {
    return(sprintf("Must contain table(s) %s", paste(missingTables, collapse = ", ")))
  }
  return(TRUE)
}

.assertCloneCensorWeightData <- function(x, add = NULL, .var.name = checkmate::vname(x)) {
  return(checkmate::makeAssertion(x, .checkCloneCensorWeightData(x), .var.name, add))
}

.checkCcwClones <- function(x, weighted = FALSE) {
  if (!inherits(x, "CcwClones")) {
    return("Must be an object of class 'CcwClones'")
  }
  if (!Andromeda::isValidAndromeda(x)) {
    return("Must be an open Andromeda object, but it was closed")
  }
  tables <- c("persons", "clones", if (weighted) c("segments", "personWeights"))
  missingTables <- setdiff(tables, names(x))
  if (length(missingTables) > 0) {
    return(sprintf(
      "Must contain table(s) %s%s", paste(missingTables, collapse = ", "),
      if (weighted) ": use the weighted clones returned by computeCensoringWeights()" else ""
    ))
  }
  if (weighted && is.null(attr(x, "metaData")$weights)) {
    return("Must have metaData$weights: use the weighted clones returned by computeCensoringWeights()")
  }
  return(TRUE)
}

# With weighted = TRUE, the clones must also have the segments and personWeights tables, and
# metaData$weights.
.assertCcwClones <- function(x, add = NULL, weighted = FALSE, .var.name = checkmate::vname(x)) {
  return(checkmate::makeAssertion(x, .checkCcwClones(x, weighted), .var.name, add))
}

# Columns of the study population (contract section 5.1) that the downstream modules use.
.studyPopulationColumns <- c(
  "rowId", "personSeqId", "cohortStartDate", "daysFromObsStart", "daysToCohortEnd", "daysToObsEnd",
  "daysToDeath", "riskStart", "riskEnd", "daysToEvent"
)

.checkStudyPopulation <- function(x) {
  check <- checkmate::checkDataFrame(x)
  if (!isTRUE(check)) {
    return(check)
  }
  return(checkmate::checkNames(colnames(x), must.include = .studyPopulationColumns))
}

.assertStudyPopulation <- function(x, add = NULL, .var.name = checkmate::vname(x)) {
  return(checkmate::makeAssertion(x, .checkStudyPopulation(x), .var.name, add))
}

# Time basis ---------------------------------------------------------------------------------------

# The time-basis columns (Phase1Contracts section 4, without the intercept) of a pooled model, shared by
# the censoring models and the pooled logistic outcome model. The basis is defined by the interval
# indices `interval` (0-based) of the rows the model is fitted to, each row counting `multiplicity`
# times (a person's frequency weight, so frequency weights give the basis of the data with persons
# copied; rows counting 0 times are left out):
# - "indicator": one column per distinct interval except the first (the reference);
# - "linear": the interval index;
# - "polynomial": stats::poly(k, degree, raw = raw), whose centring and scaling (orthogonal) depend on
#   the distribution of the rows' intervals;
# - "naturalSpline": splines::ns(k, df), with knots at quantiles of the rows' intervals.
# The degree and df are capped at the number of distinct intervals minus 1, and a single distinct
# interval gives no columns. Columns that are zero at every distinct interval are dropped. Returns
# `build(k)`, the columns at the interval indices k, and `labels`: the interval of each indicator
# column, or the column number.
.timeBasisColumns <- function(timeBasis, interval, multiplicity = rep(1, length(interval))) {
  interval <- interval[multiplicity > 0]
  multiplicity <- multiplicity[multiplicity > 0]
  intervals <- sort(unique(interval))
  nDistinct <- length(intervals)
  if (nDistinct <= 1) {
    return(list(build = function(k) matrix(numeric(0), nrow = length(k), ncol = 0), labels = integer(0)))
  }
  if (timeBasis$type == "indicator") {
    levels <- intervals[-1]
    build <- function(k) outer(k, levels, "==") * 1
    labels <- levels
  } else if (timeBasis$type == "linear") {
    build <- function(k) matrix(as.numeric(k), ncol = 1)
    labels <- 1L
  } else {
    sample <- rep(intervals, .sumBy(multiplicity, match(interval, intervals), nDistinct))
    if (timeBasis$type == "polynomial") {
      basis <- stats::poly(sample, degree = min(timeBasis$degree, nDistinct - 1), raw = timeBasis$raw)
    } else {
      basis <- splines::ns(sample, df = min(timeBasis$df, nDistinct - 1))
    }
    build <- function(k) matrix(unclass(stats::predict(basis, k)), nrow = length(k))
    labels <- seq_len(ncol(basis))
  }
  keep <- colSums(build(intervals) != 0) > 0
  return(list(
    build = function(k) build(k)[, keep, drop = FALSE],
    labels = labels[keep]
  ))
}

#' Apply study population inclusion/exclusion criteria
#'
#' Filters a `CloneCensorWeightData` object to the analysis-eligible study
#' population.  Filters are applied in this order:
#'
#' 1. **Risk window**: `cohortStartDate` and `cohortEndDate` are adjusted to
#'    the observation window defined by `riskWindowStart`/`startAnchor` and
#'    `riskWindowEnd`/`endAnchor`.
#' 2. **Prior outcome**: subjects with an occurrence of `outcomeId` in the
#'    window `[cohortStartDate - priorOutcomeLookback, cohortStartDate)` are
#'    removed when `removeSubjectsWithPriorOutcome = TRUE`.
#' 3. **Minimum days at risk**: subjects whose adjusted follow-up is shorter
#'    than `minDaysAtRisk` are removed.
#'
#' The resulting attrition table is stored inside the object and can be
#' retrieved with [getAttritionTable()].
#'
#' @param cloneCensorWeightData A `CloneCensorWeightData` object, usually
#'   returned by [getDbCloneCensorWeightData()].
#' @param outcomeId Integer.  Outcome cohort definition identifier used to
#'   test for prior occurrence.
#' @param createStudyPopulationArgs A `CreateStudyPopulationArgs` object
#'   created by [createCreateStudyPopulationArgs()].
#'
#' @return The input object with its `cohorts` table filtered and an
#'   `attrition` table added.  The modified object is returned invisibly.
#' @export
createStudyPopulation <- function(cloneCensorWeightData,
                                  outcomeId,
                                  createStudyPopulationArgs) {
  checkmate::assertClass(cloneCensorWeightData, "CloneCensorWeightData")
  checkmate::assertIntegerish(outcomeId, len = 1, any.missing = FALSE)
  checkmate::assertR6(createStudyPopulationArgs, "CreateStudyPopulationArgs")

  args     <- createStudyPopulationArgs
  cohorts  <- .readCloneCensorWeightTable(cloneCensorWeightData, "cohorts")
  outcomes <- .readCloneCensorWeightTable(cloneCensorWeightData, "outcomes")

  attrition <- data.frame(
    description = character(),
    nSubjects   = integer(),
    nExcluded   = integer(),
    stringsAsFactors = FALSE
  )

  .addStep <- function(desc) {
    n    <- nrow(cohorts)
    prev <- if (nrow(attrition) == 0L) n else attrition$nSubjects[nrow(attrition)]
    attrition <<- rbind(attrition, data.frame(
      description = desc,
      nSubjects   = n,
      nExcluded   = prev - n,
      stringsAsFactors = FALSE
    ))
  }

  .addStep("Initial target population after washout")

  # ── 1. Risk window ──────────────────────────────────────────────────────────
  cohortStart <- as.Date(cohorts$cohortStartDate)
  cohortEnd   <- as.Date(cohorts$cohortEndDate)

  obsStart <- if (tolower(args$startAnchor) == "cohort start") {
    cohortStart + args$riskWindowStart
  } else {
    cohortEnd + args$riskWindowStart
  }
  obsEnd <- if (tolower(args$endAnchor) == "cohort start") {
    cohortStart + args$riskWindowEnd
  } else {
    cohortEnd + args$riskWindowEnd
  }

  cohorts$cohortStartDate <- pmax(cohortStart, obsStart)
  cohorts$cohortEndDate   <- pmin(cohortEnd,   obsEnd)

  # ── 2. Prior outcome ────────────────────────────────────────────────────────
  if (args$removeSubjectsWithPriorOutcome) {
    # Join outcomes to cohort start dates for vectorised comparison
    relevant <- outcomes[outcomes$outcomeId == outcomeId, , drop = FALSE]
    if (nrow(relevant) > 0) {
      byStart <- cohorts[, c("subjectId", "cohortStartDate")]
      joined  <- merge(relevant, byStart, by = "subjectId")
      oDate   <- as.Date(joined$outcomeDate)
      cStart  <- as.Date(joined$cohortStartDate)
      prior   <- joined[oDate >= cStart - args$priorOutcomeLookback & oDate < cStart, ]
      cohorts <- cohorts[!cohorts$subjectId %in% prior$subjectId, ]
    }
    .addStep("After removing subjects with prior outcome")
  }

  # ── 3. Minimum days at risk ─────────────────────────────────────────────────
  daysAtRisk <- as.numeric(
    as.Date(cohorts$cohortEndDate) - as.Date(cohorts$cohortStartDate)
  )
  cohorts <- cohorts[daysAtRisk >= args$minDaysAtRisk, ]
  .addStep(sprintf("After requiring >= %d days at risk", args$minDaysAtRisk))

  # ── Write back ───────────────────────────────────────────────────────────────
  cloneCensorWeightData[["cohorts"]]   <- cohorts
  cloneCensorWeightData[["attrition"]] <- attrition
  invisible(cloneCensorWeightData)
}

#' Retrieve the attrition table
#'
#' Returns the per-step subject counts recorded by [createStudyPopulation()].
#'
#' @param cloneCensorWeightData A `CloneCensorWeightData` object that has been
#'   passed through [createStudyPopulation()].
#'
#' @return A data frame with columns `description`, `nSubjects`, and
#'   `nExcluded`.
#' @export
getAttritionTable <- function(cloneCensorWeightData) {
  checkmate::assertClass(cloneCensorWeightData, "CloneCensorWeightData")
  if (is.null(cloneCensorWeightData$attrition)) {
    stop(
      "No attrition table found. Run createStudyPopulation() first.",
      call. = FALSE
    )
  }
  dplyr::collect(cloneCensorWeightData$attrition)
}

#' Create treatment-strategy clones and apply artificial censoring
#'
#' Duplicates each eligible subject across the supplied treatment strategies,
#' identifies protocol deviation relative to the grace period, and truncates
#' follow-up at the first deviation. The resulting cohorts are retained in a
#' disk-backed object and contain one row per subject-clone follow-up.
#'
#' @param cloneCensorWeightData An object of class
#'   `CloneCensorWeightData`, usually returned by
#'   [getDbCloneCensorWeightData()].
#' @param arms A list defining the treatment strategies. Each arm must provide
#'   a unique `arm` identifier and a strategy definition understood by the
#'   package.
#' @param gracePeriod A non-negative numeric scalar giving the number of days
#'   after `cohortStartDate` during which initiation is permitted.
#' @param washoutPeriod A non-negative numeric scalar giving the required
#'   baseline lookback period in days.
#'
#' @return The input object with its `cohorts` table replaced by the cloned,
#'   artificially censored follow-up table. The table contains `subjectId`,
#'   `cloneId`, `arm`, `cohortStartDate`, `cohortEndDate`, `survivalTime`,
#'   `event`, and `artificiallyCensored`.
#' @export
createClonesAndCensor <- function(cloneCensorWeightData,
                                  arms,
                                  gracePeriod,
                                  washoutPeriod) {
  checkmate::assertClass(cloneCensorWeightData, "CloneCensorWeightData")
  checkmate::assertList(arms, min.len = 2, any.missing = FALSE)
  checkmate::assertNumeric(gracePeriod, len = 1, lower = 0, any.missing = FALSE)
  checkmate::assertNumeric(washoutPeriod, len = 1, lower = 0, any.missing = FALSE)
  cohorts <- .readCloneCensorWeightTable(cloneCensorWeightData, "cohorts")
  outcomes <- .readCloneCensorWeightTable(cloneCensorWeightData, "outcomes")
  checkmate::assertSubset(c("subjectId", "cohortStartDate", "cohortEndDate", "initiationDate"),
                          names(cohorts))
  # Support both ArmDefinition R6 objects (armId / strategyType) and legacy
  # plain lists (arm / strategy).  Per-arm gracePeriod overrides the scalar.
  armIds <- vapply(arms, function(arm) {
    id <- arm$armId %||% arm$arm
    if (is.null(id)) stop("Each arm must have an 'armId' or 'arm' field.", call. = FALSE)
    as.character(id)
  }, character(1))
  checkmate::assert_true(!anyNA(armIds) && !anyDuplicated(armIds))

  cloneRows <- lapply(seq_along(arms), function(index) {
    arm <- arms[[index]]
    strategy <- tolower(as.character(
      arm$strategyType %||% arm$strategy %||% arm$type %||% arm$arm
    ))
    armGracePeriod <- arm$gracePeriod %||% gracePeriod
    initiationDay <- as.numeric(as.Date(cohorts$initiationDate) - as.Date(cohorts$cohortStartDate))
    isInitiateArm <- strategy %in% c("initiatewithin", "initiate", "within", "treated", "arm1")
    followsStrategy <- if (isInitiateArm) {
      !is.na(initiationDay) & initiationDay >= 0 & initiationDay <= armGracePeriod
    } else {
      is.na(initiationDay) | initiationDay > armGracePeriod
    }
    censorDay <- if (isInitiateArm) {
      # Censor at gracePeriod whenever the subject did NOT initiate within the window:
      # either never initiated (NA) or initiated after the grace period.
      ifelse(is.na(initiationDay) | initiationDay > armGracePeriod, armGracePeriod, NA_real_)
    } else {
      ifelse(!is.na(initiationDay) & initiationDay <= armGracePeriod, pmax(initiationDay, 0), NA_real_)
    }
    originalEnd  <- as.Date(cohorts$cohortEndDate)
    censorDate   <- as.Date(cohorts$cohortStartDate) + censorDay
    censorDate[is.na(censorDay)] <- originalEnd[is.na(censorDay)]
    endDate <- pmin(originalEnd, censorDate)
    event <- vapply(seq_len(nrow(cohorts)), function(row) {
      subjectOutcomes <- outcomes[outcomes$subjectId == cohorts$subjectId[row], , drop = FALSE]
      any(as.Date(subjectOutcomes$outcomeDate) >= as.Date(cohorts$cohortStartDate[row]) &
            as.Date(subjectOutcomes$outcomeDate) <= endDate[row])
    }, logical(1))
    data.frame(subjectId = cohorts$subjectId,
               cloneId = paste(cohorts$subjectId, armIds[index], sep = "_"),
               arm = armIds[index],
               cohortStartDate = as.Date(cohorts$cohortStartDate),
               cohortEndDate = endDate,
               survivalTime = pmax(as.numeric(endDate - cohorts$cohortStartDate), 0),
               event = as.integer(event & followsStrategy),
               artificiallyCensored = as.integer(!followsStrategy),
               initiationDate = as.Date(cohorts$initiationDate),
               stringsAsFactors = FALSE)
  })
  cloneCensorWeightData[["cohorts"]] <- do.call(rbind, cloneRows)
  cloneCensorWeightData
}

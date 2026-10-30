#' Create treatment-strategy clones and apply artificial censoring
#'
#' Duplicates each eligible subject across both treatment-strategy arms,
#' applies artificial censoring when observed behaviour deviates from the
#' arm's protocol, and stores the resulting cloned cohort table back into
#' the data object.
#'
#' The two censoring rules are determined by `strategyType` in each
#' `ArmDefinition`:
#'
#' * **`"initiateWithin"`** — a clone is artificially censored at day
#'   `gracePeriod` if treatment has **not** started by then.
#' * **`"deferBeyond"`** — a clone is artificially censored at the
#'   initiation date if treatment starts **within** `gracePeriod` days.
#'
#' The `event` flag is set when any row in the `outcomes` table falls within
#' the clone's (possibly truncated) follow-up window.  For multi-outcome
#' analyses the caller is responsible for filtering the outcomes table to the
#' outcome of interest before calling this function.
#'
#' @param cloneCensorWeightData A `CloneCensorWeightData` object, usually
#'   the output of [createStudyPopulation()].
#' @param createClonesAndCensorArgs A `CreateClonesAndCensorArgs` object
#'   created by [createCreateClonesAndCensorArgs()].
#'
#' @return The input object with its `cohorts` table replaced by the cloned,
#'   artificially censored follow-up table.  The table contains `subjectId`,
#'   `cloneId`, `arm`, `cohortStartDate`, `cohortEndDate`, `survivalTime`,
#'   `event`, `artificiallyCensored`, and `initiationDate`.
#' @export
createClonesAndCensor <- function(cloneCensorWeightData,
                                  createClonesAndCensorArgs) {
  checkmate::assertClass(cloneCensorWeightData, "CloneCensorWeightData")
  checkmate::assertR6(createClonesAndCensorArgs, "CreateClonesAndCensorArgs")

  cohorts  <- .readCloneCensorWeightTable(cloneCensorWeightData, "cohorts")
  outcomes <- .readCloneCensorWeightTable(cloneCensorWeightData, "outcomes")

  checkmate::assertSubset(
    c("subjectId", "cohortStartDate", "cohortEndDate", "initiationDate"),
    names(cohorts)
  )

  arms    <- createClonesAndCensorArgs$arms
  armIds  <- vapply(arms, function(a) a$armId, character(1))
  if (anyNA(armIds) || anyDuplicated(armIds)) {
    stop("Each arm must have a unique, non-missing armId.", call. = FALSE)
  }

  cloneRows <- lapply(seq_along(arms), function(index) {
    arm            <- arms[[index]]
    armId          <- arm$armId
    isInitiateArm  <- tolower(arm$strategyType) == "initiatewithin"
    armGracePeriod <- arm$gracePeriod

    initiationDay <- as.numeric(
      as.Date(cohorts$initiationDate) - as.Date(cohorts$cohortStartDate)
    )

    followsStrategy <- if (isInitiateArm) {
      !is.na(initiationDay) & initiationDay >= 0 & initiationDay <= armGracePeriod
    } else {
      is.na(initiationDay) | initiationDay > armGracePeriod
    }

    censorDay <- if (isInitiateArm) {
      # Censor at gracePeriod when no initiation occurred within the window.
      ifelse(is.na(initiationDay) | initiationDay > armGracePeriod,
             armGracePeriod, NA_real_)
    } else {
      ifelse(!is.na(initiationDay) & initiationDay <= armGracePeriod,
             pmax(initiationDay, 0), NA_real_)
    }

    originalEnd <- as.Date(cohorts$cohortEndDate)
    censorDate  <- as.Date(cohorts$cohortStartDate) + censorDay
    censorDate[is.na(censorDay)] <- originalEnd[is.na(censorDay)]
    endDate <- pmin(originalEnd, censorDate)

    event <- vapply(seq_len(nrow(cohorts)), function(i) {
      sub <- outcomes[outcomes$subjectId == cohorts$subjectId[i], , drop = FALSE]
      any(
        as.Date(sub$outcomeDate) >= as.Date(cohorts$cohortStartDate[i]) &
        as.Date(sub$outcomeDate) <= endDate[i]
      )
    }, logical(1))

    data.frame(
      subjectId            = cohorts$subjectId,
      cloneId              = paste(cohorts$subjectId, armId, sep = "_"),
      arm                  = armId,
      cohortStartDate      = as.Date(cohorts$cohortStartDate),
      cohortEndDate        = endDate,
      survivalTime         = pmax(as.numeric(endDate - as.Date(cohorts$cohortStartDate)), 0),
      event                = as.integer(event),
      artificiallyCensored = as.integer(!followsStrategy),
      initiationDate       = as.Date(cohorts$initiationDate),
      stringsAsFactors     = FALSE
    )
  })

  cloneCensorWeightData[["cohorts"]] <- do.call(rbind, cloneRows)
  cloneCensorWeightData
}

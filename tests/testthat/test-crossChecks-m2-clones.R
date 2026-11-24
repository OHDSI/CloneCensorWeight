library(testthat)

# Tier 2 cross-check (plan section 6.4, `ccwr` row): two-arm cloning, initiateWithin(G) vs
# notInitiateWithin(G), against ccwr's clone_arms() -> apply_logics(create_policy_A(),
# create_censoring_logics_A()) -> create_final_data().
#
# Mapping, with riskStart = 0 and sameDayOrdering = "eventsFirst": ccwr time = day + 1, so
# `followup` = min(outcome day, riskEnd) + 1, `grace_period` = G + 1, and `treatment` = 1 with
# `time_to_treatment` = A + 1 only when A + 1 <= followup. "Surgery" is initiateWithin(G) and "Control"
# is notInitiateWithin(G). ccwr puts the deviation first in its Control arm and the event first at G in
# its Surgery arm, and it does not tell observation end from deviation at G, so persons with the
# outcome on their initiation day A <= G, or with observation end on day G without having started, are
# dropped. ccwr has no death: deaths are observation ends (competingEvents = "none").

skipM2CrossChecks <- function() {
  skip_on_cran()
  skip_if(Sys.getenv("CCW_RUN_CROSSCHECKS") != "true")
  skip_if_not_installed("ccwr", minimum_version = "0.0.2")
}

runCcwrClones <- function(ccwrData, gracePeriod) {
  arms <- c("Control", "Surgery")
  clones <- ccwr::clone_arms(ccwrData, arms)
  policies <- ccwr::create_policy_A(
    arms, "treatment", "time_to_treatment", gracePeriod, "outcome", "followup",
    clone_outcome = "cloneOutcome", clone_followup = "fup"
  )
  clones <- ccwr::apply_logics(clones, policies)
  censoringLogics <- ccwr::create_censoring_logics_A(
    arms, "treatment", "time_to_treatment", gracePeriod, "followup",
    clone_censoring = "censoring", clone_uncensored_followup = "fupUncensored"
  )
  clones <- ccwr::apply_logics(clones, censoringLogics)
  final <- ccwr::create_final_data(
    clones,
    clone_followup = "fup", clone_outcome = "cloneOutcome", clone_censoring = "censoring", col_ids = "id"
  )
  return(list(clones = clones, final = final))
}

# Compares our clones (persons identified by rowId = ccwr `id`) with ccwr's, per clone and on the
# summed person-time, events and censorings of ccwr's split rows.
expectClonesEqualCcwr <- function(clones, ccwr) {
  for (arm in c("Surgery", "Control")) {
    strategyId <- if (arm == "Surgery") 1L else 2L
    ours <- clones[clones$strategyId == strategyId, ]
    theirs <- as.data.frame(ccwr$clones[[arm]])
    theirs <- theirs[match(ours$rowId, theirs$id), ]
    expect_equal(nrow(ours), nrow(ccwr$clones[[arm]]))
    expect_equal(ours$survivalTime, as.integer(theirs$fup), info = arm)
    expect_equal(ours$outcomeCount, as.integer(theirs$cloneOutcome), info = arm)
    expect_equal(as.integer(ours$endReason == "deviation"), as.integer(theirs$censoring), info = arm)

    final <- ccwr$final[[arm]]
    expect_equal(sum(final$Tstop - final$Tstart), sum(ours$survivalTime), info = arm)
    expect_equal(sum(final$cloneOutcome), sum(ours$outcomeCount), info = arm)
    expect_equal(sum(final$censoring), sum(ours$endReason == "deviation"), info = arm)
  }
}

test_that("Clones equal ccwr on patients13 (Maringe 2020, Figure 2)", {
  skipM2CrossChecks()
  patients <- as.data.frame(ccwr::patients13)
  # Times are whole days, and every surgery is within follow-up.
  expect_true(all(patients$followup %% 1 == 0))
  expect_true(all(patients$time_to_surgery %% 1 == 0, na.rm = TRUE))
  expect_true(all(patients$time_to_surgery <= patients$followup, na.rm = TRUE))
  gracePeriod <- 180

  # Our data on the day scale: day = ccwr time - 1. Death is the outcome, and the end of follow-up is
  # the end of observation. The risk window (30 days) is widened to cover all follow-up.
  n <- nrow(patients)
  surgery <- which(patients$surgery == 1)
  died <- which(patients$death == 1)
  data <- createCloneCensorWeightData(
    cohorts = data.frame(
      rowId = seq_len(n),
      personSeqId = seq_len(n),
      personId = patients$id,
      cohortStartDate = as.Date("2015-01-01"),
      daysFromObsStart = 365L,
      daysToCohortEnd = 0L,
      daysToObsEnd = as.integer(patients$followup - 1),
      daysToDeath = NA
    ),
    outcomes = data.frame(rowId = died, outcomeId = 3, daysToEvent = as.integer(patients$followup[died] - 1)),
    treatments = data.frame(
      rowId = surgery,
      treatmentId = 21,
      daysToStart = as.integer(patients$time_to_surgery[surgery] - 1),
      daysToEnd = as.integer(patients$time_to_surgery[surgery] - 1)
    ),
    metaData = list(targetId = 1, treatmentIds = c(treatment = 21))
  )
  on.exit(close(data))
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 9999, endAnchor = "cohort start")
  ))
  ours <- runM2Clones(
    population = population,
    data = data,
    strategies = list(m2Strategy(1, "initiateWithin", gracePeriod - 1), m2Strategy(2, "notInitiateWithin", gracePeriod - 1)),
    sameDayOrdering = "eventsFirst",
    competingEvents = "none"
  )

  ccwrData <- data.frame(
    id = seq_len(n),
    treatment = patients$surgery,
    time_to_treatment = patients$time_to_surgery,
    outcome = patients$death,
    followup = patients$followup
  )
  expectClonesEqualCcwr(ours$clones, runCcwrClones(ccwrData, gracePeriod))
})

test_that("Clones equal ccwr on the 12-person fixture", {
  skipM2CrossChecks()
  graceEnd <- 7
  tables <- createM2FixtureTables()
  data <- do.call(createCloneCensorWeightData, tables)
  on.exit(close(data))
  ours <- runM2Clones(
    population = createM2FixturePopulation(data, riskStart = 0),
    data = data,
    strategies = list(m2Strategy(1, "initiateWithin", graceEnd), m2Strategy(2, "notInitiateWithin", graceEnd)),
    sameDayOrdering = "eventsFirst",
    competingEvents = "none"
  )

  # ccwr input built from the raw fixture tables for the persons in our study population.
  rowIds <- ours$persons$rowId
  cohorts <- tables$cohorts[match(rowIds, tables$cohorts$rowId), ]
  riskEnd <- pmin(cohorts$daysToObsEnd, 30L)
  firstDay <- function(days) {
    days <- days[!is.na(days)]
    if (length(days) == 0) NA_integer_ else min(days)
  }
  outcomeDay <- vapply(seq_along(rowIds), function(i) {
    with(tables$outcomes, firstDay(daysToEvent[rowId == rowIds[i] & outcomeId == 3 & daysToEvent >= 0 & daysToEvent <= riskEnd[i]]))
  }, integer(1))
  startDay <- vapply(rowIds, function(r) {
    with(tables$treatments, firstDay(daysToStart[rowId == r & treatmentId == 21 & daysToStart >= 0]))
  }, integer(1))
  followup <- pmin(outcomeDay, riskEnd, na.rm = TRUE) + 1
  treated <- !is.na(startDay) & startDay + 1 <= followup
  ccwrData <- data.frame(
    id = rowIds,
    treatment = as.numeric(treated),
    time_to_treatment = ifelse(treated, startDay + 1, NA),
    outcome = as.numeric(!is.na(outcomeDay)),
    followup = followup
  )

  outcomeOnStartDay <- !is.na(outcomeDay) & !is.na(startDay) & outcomeDay == startDay & startDay <= graceEnd
  observationEndAtGraceEnd <- riskEnd == graceEnd & (is.na(startDay) | startDay > graceEnd) &
    (is.na(outcomeDay) | outcomeDay != graceEnd)
  keep <- !outcomeOnStartDay & !observationEndAtGraceEnd
  # Persons 8 (outcome on the initiation day) and 9 (observation ends on day G) are dropped.
  expect_equal(rowIds[!keep], c(8L, 9L))

  expectClonesEqualCcwr(ours$clones[ours$clones$rowId %in% rowIds[keep], ], runCcwrClones(ccwrData[keep, ], graceEnd + 1))
})

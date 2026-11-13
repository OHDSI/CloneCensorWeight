library(testthat)

# Tier 2 cross-check (plan section 6.4, the `ccwr` outcome-stage row; a Phase 1 exit item, plan section
# 8): ccwr::lungcancer (200 simulated patients from the Maringe et al. 2020 supplement) through ccwr's own
# clone_arms() -> apply_logics(create_policy_A(), create_censoring_logics_A()) -> create_final_data() ->
# estimate_censoring(method = "Cox") -> weight_cases(), with the 7 Maringe covariates as censoring
# predictors. Given ccwr's clones and weights, our outcome models run on its rows mapped to segments:
# startTime = Tstart, stopTime = Tstop, weight = weight_Cox, personSeqId = id, one cloneRowId per id x arm,
# "Surgery" = strategy 1 (initiateWithin, the target) and "Control" = strategy 2 (notInitiateWithin).
# riskStart = 0, so ccwr time t is day t - 1, and the 1-year horizon h = 364 is read at t_h = 365 on
# ccwr's scale. Checks:
# - our weighted KM risks equal ccwr::emul_estimate(method = "KM") to 1e-10;
# - our Cox MSM log HR and robust SE equal ccwr::emul_estimate(method = "Cox"), the weighted coxph with
#   cluster(id) and robust = TRUE, refitted with Breslow ties, to 1e-6. ccwr hard-codes Efron ties, and the
#   clones of the two arms share death times, so the two tie methods differ (log HR -0.49875 with Efron vs
#   -0.49844 with Breslow).
# Nightly only (CCW_RUN_CROSSCHECKS=true).

skipLungcancerCrossChecks <- function() {
  skip_on_cran()
  skip_if(Sys.getenv("CCW_RUN_CROSSCHECKS") != "true", "Set CCW_RUN_CROSSCHECKS=true to run the cross-checks")
  skip_if_not_installed("ccwr", minimum_version = "0.0.2")
  skip_if_not_installed("survival")
}

lungcancerArms <- c("Control", "Surgery")
lungcancerPredictors <- c("sex", "charlson", "perf", "stage", "emergency", "age", "deprivation")

# ccwr's clones and Cox censoring weights of `lungcancer`, one data frame per arm.
runCcwrLungcancer <- function(lungcancer, gracePeriod) {
  arms <- lungcancerArms
  clones <- ccwr::clone_arms(lungcancer, arms)
  clones <- ccwr::apply_logics(clones, ccwr::create_policy_A(
    arms, "surgery", "timetosurgery", gracePeriod, "death", "fup_obs",
    clone_outcome = "outcome", clone_followup = "fup"
  ))
  clones <- ccwr::apply_logics(clones, ccwr::create_censoring_logics_A(
    arms, "surgery", "timetosurgery", gracePeriod, "fup_obs",
    clone_censoring = "censoring", clone_uncensored_followup = "fupUncensored"
  ))
  final <- ccwr::create_final_data(
    clones,
    clone_followup = "fup", clone_outcome = "outcome", clone_censoring = "censoring",
    col_ids = c("id", lungcancerPredictors)
  )
  return(ccwr::weight_cases(ccwr::estimate_censoring(final, predictors = lungcancerPredictors, method = "Cox")))
}

# Our weighted clones (the caller closes them) built from ccwr's weighted rows, with createM4WeightedClones()
# (helper-m4-outcomes.R). ccwr's times need not be whole days (administrative censoring is at 365.24), so
# the clones table, which only feeds the counts, ends each clone on the day that holds its last time.
ccwrRowsToWeightedClones <- function(weighted, lungcancer, graceEnd) {
  rows <- dplyr::bind_rows(lapply(weighted, as.data.frame), .id = "arm")
  rows$strategyId <- ifelse(rows$arm == "Surgery", 1L, 2L)
  ids <- sort(unique(rows$id))
  rows$personSeqId <- match(rows$id, ids)
  rows <- rows[order(rows$personSeqId, rows$strategyId, rows$Tstart), ]
  cloneKey <- paste(rows$personSeqId, rows$strategyId)
  rows$cloneRowId <- match(cloneKey, unique(cloneKey))
  segments <- data.frame(
    segmentId = seq_len(nrow(rows)),
    cloneRowId = rows$cloneRowId,
    rowId = rows$personSeqId,
    personSeqId = rows$personSeqId,
    strategyId = rows$strategyId,
    startTime = rows$Tstart,
    stopTime = rows$Tstop,
    weight = unname(rows$weight_Cox),
    outcome = as.integer(rows$outcome),
    competing = 0L
  )
  last <- !duplicated(segments$cloneRowId, fromLast = TRUE)
  endTime <- segments$stopTime[last]
  clones <- data.frame(
    cloneRowId = segments$cloneRowId[last],
    rowId = segments$rowId[last],
    personSeqId = segments$personSeqId[last],
    strategyId = segments$strategyId[last],
    riskStart = 0L,
    deviationDay = NA_integer_,
    endDay = as.integer(ceiling(endTime) - 1),
    survivalTime = as.integer(ceiling(endTime)),
    endReason = ifelse(segments$outcome[last] == 1, "outcome", ifelse(rows$censoring[last] == 1, "deviation", "riskEnd")),
    outcomeCount = segments$outcome[last]
  )
  patients <- lungcancer[match(ids, lungcancer$id), ]
  persons <- data.frame(
    rowId = seq_along(ids),
    personSeqId = seq_along(ids),
    riskStart = 0L,
    riskEnd = as.integer(ceiling(patients$fup_obs) - 1),
    nominalRiskEnd = 365L,
    outcomeDay = as.integer(ifelse(patients$death == 1, ceiling(patients$fup_obs) - 1, NA)),
    deathDay = NA_integer_,
    treatmentDay = as.integer(ifelse(patients$surgery == 1, patients$timetosurgery - 1, NA)),
    comparatorTreatmentDay = NA_integer_
  )
  strategies <- list(
    createTreatmentStrategy(1, "Surgery", "initiateWithin", graceEnd = graceEnd),
    createTreatmentStrategy(2, "Control", "notInitiateWithin", graceEnd = graceEnd)
  )
  fixture <- list(
    persons = persons,
    clones = clones,
    segments = segments,
    personWeights = m4PersonWeights(segments),
    metaData = list(
      strategies = lapply(strategies, function(strategy) strategy$toList()),
      outcomeId = 1,
      treatmentIds = c(treatment = 21),
      sameDayOrdering = "eventsFirst",
      competingEvents = "none",
      attrition = dplyr::tibble(
        description = "ccwr::lungcancer",
        targetPersons = length(ids),
        comparatorPersons = 0L,
        targetExposures = length(ids),
        comparatorExposures = 0L
      ),
      targetId = 1
    )
  )
  return(createM4WeightedClones(fixture))
}

# ccwr's rows of both arms, bound as ccwr::emul_estimate() binds them.
bindCcwrArms <- function(weighted) {
  rows <- dplyr::bind_rows(weighted, .id = "arms")
  rows$arms <- factor(rows$arms, levels = lungcancerArms)
  return(rows)
}

test_that("Outcome stage: KM risks and Cox MSM on ccwr's lungcancer clones and weights equal ccwr's", {
  skipLungcancerCrossChecks()
  lungcancer <- as.data.frame(ccwr::lungcancer)
  # ccwr's vignette grace period (182.62 days, so surgery on day t <= 182 on ccwr's scale, our A <= 181).
  weighted <- runCcwrLungcancer(lungcancer, gracePeriod = 182.62)
  for (arm in lungcancerArms) {
    rows <- weighted[[arm]]
    expect_equal(length(unique(rows$id)), nrow(lungcancer), info = arm)
    expect_true(all(rows$Tstop > rows$Tstart), info = arm)
    expect_true(all(is.finite(rows$weight_Cox) & rows$weight_Cox >= 1), info = arm)
  }
  weightedClones <- ccwrRowsToWeightedClones(weighted, lungcancer, graceEnd = 181)
  on.exit(close(weightedClones))
  km <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(modelType = "km", horizons = 364)))
  cox <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(horizons = 364)))
  expect_equal(km$counts$clones, c(200L, 200L))
  expect_equal(km$counts$outcomes, c(sum(weighted$Surgery$outcome), sum(weighted$Control$outcome)))

  # KM risks at t = 365 on ccwr's scale (h = 364):
  ccwrKm <- summary(ccwr::emul_estimate(weighted, method = "KM", weights = "weight_Cox"), times = 365, extend = TRUE)
  ccwrRisk <- stats::setNames(1 - ccwrKm$surv, sub("arms=", "", ccwrKm$strata))
  for (model in list(km, cox)) {
    expect_equal(getM4Estimate(model, "risk_target_364"), unname(ccwrRisk["Surgery"]), tolerance = 1e-10)
    expect_equal(getM4Estimate(model, "risk_comparator_364"), unname(ccwrRisk["Control"]), tolerance = 1e-10)
    expect_equal(getM4Estimate(model, "rd_364"), unname(ccwrRisk["Surgery"] - ccwrRisk["Control"]), tolerance = 1e-10)
  }

  # The Cox MSM. ccwr's model is coxph(Surv(Tstart, Tstop, outcome) ~ arms + cluster(id), weights,
  # robust = TRUE, ties = "efron"); the reference is that model with Breslow ties.
  ccwrCox <- ccwr::emul_estimate(weighted, method = "Cox", weights = "weight_Cox")
  rows <- bindCcwrArms(weighted)
  control <- survival::coxph.control(eps = 1e-12, toler.chol = 1e-14, iter.max = 100)
  efron <- survival::coxph(
    survival::Surv(Tstart, Tstop, outcome) ~ arms + cluster(id),
    data = rows, weights = weight_Cox, robust = TRUE, ties = "efron", control = control
  )
  expect_equal(unname(coef(efron)), unname(coef(ccwrCox)), tolerance = 1e-8)
  expect_equal(sqrt(stats::vcov(efron)[1, 1]), sqrt(stats::vcov(ccwrCox)[1, 1]), tolerance = 1e-8)
  breslow <- survival::coxph(
    survival::Surv(Tstart, Tstop, outcome) ~ arms + cluster(id),
    data = rows, weights = weight_Cox, robust = TRUE, ties = "breslow", control = control
  )
  expect_equal(getM4Estimate(cox, "hr", "logEstimate"), unname(coef(breslow)), tolerance = 1e-6)
  expect_equal(getM4Estimate(cox, "hr", "seLogEstimate"), sqrt(stats::vcov(breslow)[1, 1]), tolerance = 1e-6)
})

# The full pipelines (plan section 6.4, Tier 3 `ccwr` full-pipeline row): reported, not gated. Times are
# rounded to whole days for both packages (only the administrative censoring time changes, 365.24 -> 365),
# ccwr's grace period is 182 days (our G = 181), and both use the 7 Maringe covariates as censoring
# predictors without a prior. ccwr: Cox censoring models in both arms (Efron), unstabilized weights,
# emul_estimate(method = "KM") and "Cox" with predictors = NULL. Ours: initiateWithin(181) (binary-at-G
# model) vs notInitiateWithin(181) (Cox initiation-hazard model, Breslow), unstabilized weights without
# truncation, "km" and the default "cox" report. Exact agreement is impossible: ccwr's initiate arm uses a
# Cox model with every censoring event at G instead of a binary model at G, its ties are Efron, and its
# same-day rules differ (plan section 6.4).
#
# Offline run, 2026-11-13, with 500 bootstrap replicates of the same person resamples on both sides (ours
# through bootstrapCcwAnalysis(counts = ), ccwr's by rerunning the steps above), no failed replicate:
#
#                 ours      ccwr   difference  ours SD  ccwr SD  difference / SD  SD ratio
#   riskSurgery   0.1828    0.1840  -1.25e-03   0.0330   0.0320   -0.038           1.03
#   riskControl   0.3034    0.3047  -1.32e-03   0.0494   0.0514   -0.027           0.96
#   rd           -0.1206   -0.1207   6.6e-05    0.0550   0.0561    0.001           0.98
#   log HR       -0.5144   -0.4988  -1.57e-02   0.2472   0.2432   -0.063           1.02
#
# The 1-year RD and log HR agree to well within one bootstrap SD, and the bootstrap SDs within +-10% (the
# plan's reporting targets). ccwr's own emul_estimate_bootstrap(n_bootstrap = 500, seed = 1) gives SD(log
# HR) = 0.2499 (ratio 0.99). The HRs (0.598 ours, 0.607 ccwr) and the 1-year RD (-0.12) are within the
# Maringe gist's sanity targets (HR about 0.607 +- 0.01, survival difference about 0.12).
test_that("Full pipeline vs ccwr on lungcancer (reported, not gated)", {
  skipLungcancerCrossChecks()
  lungcancer <- as.data.frame(ccwr::lungcancer)
  lungcancer$fup_obs <- round(lungcancer$fup_obs)
  gracePeriod <- 182
  graceEnd <- gracePeriod - 1

  weighted <- runCcwrLungcancer(lungcancer, gracePeriod)
  ccwrKm <- summary(ccwr::emul_estimate(weighted, method = "KM", weights = "weight_Cox"), times = 365, extend = TRUE)
  ccwrRisk <- stats::setNames(1 - ccwrKm$surv, sub("arms=", "", ccwrKm$strata))
  ccwrLogHr <- unname(coef(ccwr::emul_estimate(weighted, method = "Cox", weights = "weight_Cox")))

  # Our data on the day scale (day = ccwr time - 1). Death is the outcome; there is no competing death.
  n <- nrow(lungcancer)
  died <- which(lungcancer$death == 1)
  operated <- which(lungcancer$surgery == 1)
  # Linear terms as in ccwr (sex and stage, coded 1/2, as indicators):
  covariateMatrix <- cbind(
    female = as.numeric(lungcancer$sex == 2),
    charlson = lungcancer$charlson,
    perf = lungcancer$perf,
    stage2 = as.numeric(lungcancer$stage == 2),
    emergency = lungcancer$emergency,
    age = lungcancer$age,
    deprivation = lungcancer$deprivation
  )
  covariateIds <- seq_len(ncol(covariateMatrix)) * 1000 + 1
  covariates <- data.frame(
    rowId = rep(seq_len(n), ncol(covariateMatrix)),
    covariateId = rep(covariateIds, each = n),
    covariateValue = as.vector(covariateMatrix)
  )
  data <- createCloneCensorWeightData(
    cohorts = data.frame(
      rowId = seq_len(n),
      personSeqId = seq_len(n),
      personId = lungcancer$id,
      cohortStartDate = as.Date("2015-01-01"),
      daysFromObsStart = 365L,
      daysToCohortEnd = 0L,
      daysToObsEnd = as.integer(lungcancer$fup_obs - 1),
      daysToDeath = NA_integer_
    ),
    outcomes = data.frame(rowId = died, outcomeId = 1, daysToEvent = as.integer(lungcancer$fup_obs[died] - 1)),
    treatments = data.frame(
      rowId = operated,
      treatmentId = 21,
      daysToStart = as.integer(lungcancer$timetosurgery[operated] - 1),
      daysToEnd = as.integer(lungcancer$timetosurgery[operated] - 1)
    ),
    covariates = covariates[covariates$covariateValue != 0, ],
    covariateRef = data.frame(covariateId = covariateIds, covariateName = colnames(covariateMatrix), analysisId = seq_along(covariateIds), conceptId = 0),
    analysisRef = data.frame(
      analysisId = seq_along(covariateIds),
      analysisName = colnames(covariateMatrix),
      domainId = "Demographics",
      startDay = 0,
      endDay = 0,
      isBinary = ifelse(colnames(covariateMatrix) %in% c("female", "stage2", "emergency"), "Y", "N"),
      missingMeansZero = "Y"
    ),
    metaData = list(targetId = 1, treatmentIds = c(treatment = 21))
  )
  on.exit(close(data))
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 1,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowStart = 0, riskWindowEnd = 364, endAnchor = "cohort start")
  ))
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(
    strategies = list(
      createTreatmentStrategy(1, "Surgery within 182 days", "initiateWithin", graceEnd = graceEnd),
      createTreatmentStrategy(2, "No surgery within 182 days", "notInitiateWithin", graceEnd = graceEnd)
    ),
    competingEvents = "none"
  )))
  on.exit(close(clones), add = TRUE)
  models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(
    hazardModelType = "cox",
    prior = Cyclops::createPrior("none")
  )))
  weightedClones <- suppressMessages(computeCensoringWeights(
    clones, models, createComputeWeightsArgs(stabilized = FALSE, truncation = createWeightTruncation("none"))
  ))
  on.exit(close(weightedClones), add = TRUE)
  cox <- suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(horizons = 364)))
  ours <- c(
    riskSurgery = getM4Estimate(cox, "risk_target_364"),
    riskControl = getM4Estimate(cox, "risk_comparator_364"),
    rd = getM4Estimate(cox, "rd_364"),
    logHr = getM4Estimate(cox, "hr", "logEstimate")
  )
  theirs <- c(riskSurgery = unname(ccwrRisk["Surgery"]), riskControl = unname(ccwrRisk["Control"]), rd = unname(ccwrRisk["Surgery"] - ccwrRisk["Control"]), logHr = ccwrLogHr)
  robustSe <- c(
    riskSurgery = getM4Estimate(cox, "risk_target_364", "se"),
    riskControl = getM4Estimate(cox, "risk_comparator_364", "se"),
    rd = getM4Estimate(cox, "rd_364", "se"),
    logHr = getM4Estimate(cox, "hr", "seLogEstimate")
  )
  expect_true(all(is.finite(ours)) && all(is.finite(theirs)) && all(robustSe > 0))
  report <- data.frame(ours = ours, ccwr = theirs, difference = ours - theirs, ourRobustSe = robustSe, differenceInSe = (ours - theirs) / robustSe)
  message(
    "Full pipeline vs ccwr on lungcancer, 1-year horizon (reported, not gated):\n",
    paste(utils::capture.output(print(signif(report, 4))), collapse = "\n")
  )
})

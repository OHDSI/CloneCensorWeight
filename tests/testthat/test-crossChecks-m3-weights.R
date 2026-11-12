library(testthat)

# Tier 2 cross-checks for the censoring weights (plan section 6.4): weight products from ipw::ipwtm
# (pooled logistic) and the Cox initiation-hazard weights of ccwr. Nightly only.

skipCrossChecks <- function() {
  skip_on_cran()
  skip_if(Sys.getenv("CCW_RUN_CROSSCHECKS") != "true", "Set CCW_RUN_CROSSCHECKS=true to run the cross-checks")
}

test_that("Pooled logistic weights equal ipw::ipwtm weight products", {
  skipCrossChecks()
  skip_if_not_installed("ipw", minimum_version = "1.3.2")
  graceEnd <- 8
  fixture <- createM3Fixture(150, riskStart = 0L, seed = 31)
  data <- createM3Data(fixture)
  on.exit(close(data))
  clones <- buildM3Clones(fixture$persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = graceEnd), "eventsFirst")
  on.exit(close(clones), add = TRUE)
  # A linear time basis, so the ipw glm has the same design (and no separation in empty days):
  models <- suppressMessages(fitCensoringModels(clones, data, m3Args("logistic", timeBasis = createTimeBasis("linear"))))

  # One arm's K_g rows, up to and including the initiation row:
  rows <- refIntervalRows(fixture$persons, graceEnd, "eventsFirst", 1)
  rows <- cbind(rows, m3CovariateMatrix(fixture, rows$rowId))
  rows <- rows[order(rows$rowId, rows$k), ]
  for (stabilized in c(FALSE, TRUE)) {
    weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = stabilized, truncation = createWeightTruncation("none"))))
    segments <- collectSegments(weighted)
    segments <- segments[segments$strategyId == 2, ]
    ipwFit <- ipw::ipwtm(
      exposure = initiated,
      family = "binomial",
      link = "logit",
      numerator = if (stabilized) ~k else NULL,
      denominator = ~ k + x1 + x2 + x3,
      id = rowId,
      timevar = k,
      type = "first",
      data = rows
    )
    # Our W on day d equals the ipw weight of the row of day d - 1 (lag 1 within id), and W(0) = 1.
    lastDay <- pmin(segments$stopTime - 1 - 1, graceEnd)
    expected <- ipwFit$ipw.weights[match(paste(segments$rowId, lastDay), paste(rows$rowId, rows$k))]
    expected[lastDay < 0] <- 1
    expect_false(anyNA(expected))
    expect_equal(segments$weight, expected, tolerance = 1e-6)
  }
})

test_that("Cox initiation-hazard weights equal those of ccwr on untied data", {
  skipCrossChecks()
  skip_if_not_installed("ccwr", minimum_version = "0.0.2")
  skip_if_not_installed("survival")
  graceEnd <- 20
  fixture <- createM3Fixture(80, riskStart = 0L, seed = 32, untiedInitiation = TRUE, deathFraction = 0, maxTreatmentDay = 40L)
  # ccwr uses Efron ties and orders events differently from us, so keep only untied cases: distinct
  # initiation days, and no outcome on a day someone initiates within the grace period.
  persons <- fixture$persons
  initiationDays <- persons$treatmentDay[!is.na(persons$treatmentDay) & persons$treatmentDay <= graceEnd]
  expect_false(anyDuplicated(initiationDays) > 0)
  keep <- is.na(persons$outcomeDay) | !(persons$outcomeDay %in% c(initiationDays, persons$treatmentDay))
  fixture$persons <- persons[keep, ]
  fixture$persons$rowId <- seq_len(nrow(fixture$persons))
  fixture$persons$personSeqId <- fixture$persons$rowId
  fixture$x <- fixture$x[keep, ]
  fixture$covariates <- data.frame(
    rowId = rep(fixture$persons$rowId, 3),
    covariateId = rep(c(1001, 1002, 1003), each = nrow(fixture$x)),
    covariateValue = as.vector(fixture$x)
  )
  fixture$covariates <- fixture$covariates[fixture$covariates$covariateValue != 0, ]
  persons <- fixture$persons
  data <- createM3Data(fixture)
  on.exit(close(data))
  clones <- buildM3Clones(persons, m3Strategies("initiateWithin", "notInitiateWithin", graceEnd = graceEnd), "eventsFirst")
  on.exit(close(clones), add = TRUE)
  args <- m3Args("cox")
  # The binary model of the other strategy may separate on this fixture; it is not compared:
  args$stopOnError <- FALSE
  models <- suppressWarnings(suppressMessages(fitCensoringModels(clones, data, args)))
  weighted <- suppressMessages(computeCensoringWeights(clones, models, createComputeWeightsArgs(stabilized = FALSE, truncation = createWeightTruncation("none"))))
  segments <- collectSegments(weighted)
  segments <- segments[segments$strategyId == 2, ]

  # ccwr input on its time scale (day d is d + 1; riskStart = 0): followup = min(outcome, risk end) + 1,
  # time_to_treatment = A + 1 when A + 1 <= followup, grace period G + 1.
  followup <- pmin(persons$outcomeDay, persons$riskEnd, na.rm = TRUE) + 1
  treated <- !is.na(persons$treatmentDay) & persons$treatmentDay + 1 <= followup
  input <- data.frame(
    id = persons$rowId,
    treatment = as.integer(treated),
    time_to_treatment = ifelse(treated, persons$treatmentDay + 1, NA),
    outcome = as.integer(!is.na(persons$outcomeDay) & persons$outcomeDay + 1 == followup),
    followup = followup,
    x1 = fixture$x[, 1],
    x2 = fixture$x[, 2],
    x3 = fixture$x[, 3]
  )
  arms <- c("Control", "Surgery")
  ccwrClones <- ccwr::clone_arms(input, arms)
  ccwrClones <- ccwr::apply_logics(ccwrClones, ccwr::create_policy_A(arms, "treatment", "time_to_treatment", graceEnd + 1, "outcome", "followup"))
  ccwrClones <- ccwr::apply_logics(ccwrClones, ccwr::create_censoring_logics_A(arms, "treatment", "time_to_treatment", graceEnd + 1, "followup"))
  finalData <- ccwr::create_final_data(ccwrClones, ".fup", ".outcome", ".censoring", c("id", "x1", "x2", "x3"))
  censoring <- ccwr::estimate_censoring(finalData, predictors = c("x1", "x2", "x3"), method = "Cox", id = "id", censoring = ".censoring")
  control <- ccwr::weight_cases(censoring)$Control

  # ccwr's Cox model equals ours (Efron = Breslow without tied initiations):
  expect_equal(unname(models[[2]]$coefficients), unname(coef(survival::coxph(
    survival::Surv(Tstart, Tstop, .censoring) ~ x1 + x2 + x3,
    data = finalData$Control, ties = "efron", control = m3CoxControl()
  ))), tolerance = 1e-6)
  # ccwr splits at every follow-up time, so each of our segments ends where one of its rows ends:
  expected <- unname(control$weight_Cox)[match(paste(segments$rowId, segments$stopTime), paste(control$id, control$Tstop))]
  expect_false(anyNA(expected))
  expect_equal(segments$weight, expected, tolerance = 1e-6)
})

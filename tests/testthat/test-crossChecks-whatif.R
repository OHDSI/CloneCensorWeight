library(testthat)

# The What If G = 0 reduction (plan section 6.4, Tier 3 row "G = 0 reduction"; a Phase 2 exit item):
# Program 17.3 of Hernan & Robins, What If (survival curves from an IP-weighted pooled logistic hazards
# model on NHEFS), against initiateWithin(0) vs notInitiateWithin(0). With G = 0, K = {0} in both arms
# and both censoring models reduce to one logistic model at day 0: the book's propensity model. Months
# are coded as days with riskStart = 1, so interval k (0-based) is the book's `time`. Gate: 1e-6 on the
# outcome-model coefficients (signs flipped: the book models event == 0) and on the risks at months
# 1..120, against a fresh run of the book code.
# NHEFS is not shipped and never committed: the test reads it from the file named by the environment
# variable CCW_NHEFS_PATH and skips without it. Source, date and checksum: extras/data/nhefs/SOURCE.md
# (the owner's checkout); the reviewed transcription of the book code it follows is
# extras/comparison/whatif/program17_3.R.

skipWhatIf <- function() {
  skip_on_cran()
  path <- Sys.getenv("CCW_NHEFS_PATH")
  skip_if(path == "" || !file.exists(path), "Set CCW_NHEFS_PATH to the NHEFS csv to run the What If check")
  return(path)
}

# Program 17.3 (book code; splitstackshape::expandRows() replaced by rep(), plots left out).
whatIfProgram17_3 <- function(nhefs) {
  nhefs$survtime <- ifelse(nhefs$death == 0, 120, (nhefs$yrdth - 83) * 12 + nhefs$modth)
  p.denom <- glm(qsmk ~ sex + race + age + I(age * age) + as.factor(education) +
    smokeintensity + I(smokeintensity * smokeintensity) + smokeyrs + I(smokeyrs * smokeyrs) +
    as.factor(exercise) + as.factor(active) + wt71 + I(wt71 * wt71), data = nhefs, family = binomial())
  pd <- predict(p.denom, nhefs, type = "response")
  pn <- predict(glm(qsmk ~ 1, data = nhefs, family = binomial()), nhefs, type = "response")
  nhefs$sw.a <- ifelse(nhefs$qsmk == 1, pn / pd, (1 - pn) / (1 - pd))
  nhefs.ipw <- nhefs[rep(seq_len(nrow(nhefs)), nhefs$survtime), ]
  nhefs.ipw$time <- sequence(rle(nhefs.ipw$seqn)$lengths) - 1
  nhefs.ipw$event <- ifelse(nhefs.ipw$time == nhefs.ipw$survtime - 1 & nhefs.ipw$death == 1, 1, 0)
  nhefs.ipw$timesq <- nhefs.ipw$time^2
  ipw.model <- suppressWarnings(glm(event == 0 ~ qsmk + I(qsmk * time) + I(qsmk * timesq) + time + timesq,
    family = binomial(), weights = sw.a, data = nhefs.ipw
  ))
  grid <- data.frame(time = 0:119, timesq = (0:119)^2)
  survival <- function(qsmk) cumprod(predict(ipw.model, cbind(grid, qsmk = qsmk), type = "response"))
  return(list(model = ipw.model, weights = nhefs$sw.a, risk1 = 1 - survival(1), risk0 = 1 - survival(0)))
}

whatIfCcwData <- function(nhefs) {
  survtime <- ifelse(nhefs$death == 0, 120, (nhefs$yrdth - 83) * 12 + nhefs$modth)
  x <- data.frame(
    sex = nhefs$sex, race = nhefs$race, age = nhefs$age, ageSq = nhefs$age^2,
    education2 = nhefs$education == 2, education3 = nhefs$education == 3,
    education4 = nhefs$education == 4, education5 = nhefs$education == 5,
    smokeintensity = nhefs$smokeintensity, smokeintensitySq = nhefs$smokeintensity^2,
    smokeyrs = nhefs$smokeyrs, smokeyrsSq = nhefs$smokeyrs^2,
    exercise1 = nhefs$exercise == 1, exercise2 = nhefs$exercise == 2,
    active1 = nhefs$active == 1, active2 = nhefs$active == 2,
    wt71 = nhefs$wt71, wt71Sq = nhefs$wt71^2
  )
  n <- nrow(nhefs)
  covariates <- data.frame(
    rowId = rep(seq_len(n), ncol(x)),
    covariateId = rep(1000 * seq_along(x) + 2, each = n),
    covariateValue = as.numeric(unlist(x, use.names = FALSE))
  )
  died <- which(nhefs$death == 1)
  quit <- which(nhefs$qsmk == 1)
  return(createCloneCensorWeightData(
    cohorts = data.frame(
      rowId = seq_len(n), personSeqId = seq_len(n), personId = as.character(nhefs$seqn),
      cohortStartDate = as.Date("1971-01-01"), daysFromObsStart = 365L, daysToCohortEnd = 0L,
      daysToObsEnd = 120L, daysToDeath = NA_integer_
    ),
    outcomes = data.frame(rowId = died, outcomeId = 3, daysToEvent = as.integer(survtime[died])),
    treatments = data.frame(rowId = quit, treatmentId = 21, daysToStart = 0L, daysToEnd = 120L),
    covariates = covariates[covariates$covariateValue != 0, ],
    covariateRef = data.frame(covariateId = 1000 * seq_along(x) + 2, covariateName = names(x), analysisId = 2, conceptId = 0),
    analysisRef = data.frame(
      analysisId = 2, analysisName = "NHEFS", domainId = "NHEFS", startDay = -365, endDay = -1,
      isBinary = "N", missingMeansZero = "Y"
    ),
    metaData = list(targetId = 1, treatmentIds = c(treatment = 21), outcomeIds = 3)
  ))
}

test_that("initiateWithin(0) vs notInitiateWithin(0) reproduces What If Program 17.3 on NHEFS", {
  path <- skipWhatIf()
  nhefs <- utils::read.csv(path)
  expect_equal(nrow(nhefs), 1629L)
  book <- whatIfProgram17_3(nhefs)

  data <- whatIfCcwData(nhefs)
  on.exit(close(data))
  population <- createStudyPopulation(data, outcomeId = 3, createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
    removeSubjectsWithPriorOutcome = FALSE, minDaysAtRisk = 1, riskWindowStart = 1,
    startAnchor = "cohort start", riskWindowEnd = 120, endAnchor = "cohort start"
  ))
  clones <- createClones(population, data, createCreateClonesArgs(strategies = list(
    createTreatmentStrategy(1, "Quit", "initiateWithin", graceEnd = 0),
    createTreatmentStrategy(2, "Do not quit", "notInitiateWithin", graceEnd = 0)
  ), competingEvents = "none"))
  on.exit(close(clones), add = TRUE)
  models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(
    covariateTimeInteraction = "none", prior = Cyclops::createPrior("none"),
    control = Cyclops::createControl(algorithm = "newton", convergenceType = "gradient", tolerance = 1e-10, noiseLevel = "silent")
  )))
  timeBasis <- createTimeBasis(type = "polynomial", intervalLength = 1, degree = 2, raw = TRUE)
  weighted <- suppressMessages(computeCensoringWeights(
    clones, models, createComputeWeightsArgs(truncation = createWeightTruncation(type = "none")),
    outcomeTimeBasis = timeBasis
  ))
  on.exit(close(weighted), add = TRUE)
  fit <- suppressMessages(fitOutcomeModel(weighted, createFitOutcomeModelArgs(
    modelType = "logistic", horizons = 1:120, timeBasis = timeBasis, armByTimeInteraction = TRUE, msmWeights = "stabilized"
  )))

  # One clone per person: quitters in the initiate arm, the others in the comparator arm, with the
  # book's person-months and its stabilized weights.
  expect_equal(fit$counts$clones, c(sum(nhefs$qsmk == 1), sum(nhefs$qsmk == 0)))
  expect_equal(sum(fit$counts$personDays), sum(ifelse(nhefs$death == 0, 120, (nhefs$yrdth - 83) * 12 + nhefs$modth)))
  personWeights <- dplyr::collect(weighted$personWeights)
  cloneRows <- dplyr::collect(weighted$clones)
  rowIds <- cloneRows$rowId[match(personWeights$cloneRowId, cloneRows$cloneRowId)]
  expect_equal(personWeights$maxWeight, book$weights[rowIds], tolerance = 1e-6)

  bookCoefficients <- coef(book$model)[c("(Intercept)", "qsmk", "time", "timesq", "I(qsmk * time)", "I(qsmk * timesq)")]
  ours <- fit$coefficients[c("(Intercept)", "target", "time1", "time2", "target:time1", "target:time2")]
  expect_lt(max(abs(unname(ours) + unname(bookCoefficients))), 1e-6)
  risk <- function(arm) fit$estimates$estimate[match(sprintf("risk_%s_%d", arm, 1:120), fit$estimates$estimand)]
  expect_lt(max(abs(risk("target") - book$risk1)), 1e-6)
  expect_lt(max(abs(risk("comparator") - book$risk0)), 1e-6)
})

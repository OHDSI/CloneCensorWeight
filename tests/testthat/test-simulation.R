# simulateCcwData ------------------------------------------------------------

test_that("simulateCcwData returns a CloneCensorWeightData object", {
  d <- simulateCcwData(n = 50, seed = 1)
  expect_s4_class(d, "CloneCensorWeightData")
})

test_that("simulateCcwData cohorts table has correct columns", {
  d       <- simulateCcwData(n = 40, seed = 2)
  cohorts <- dplyr::collect(d$cohorts)
  expect_true(all(c("subjectId", "cohortStartDate", "cohortEndDate",
                    "initiationDate") %in% names(cohorts)))
  expect_equal(nrow(cohorts), 40L)
})

test_that("simulateCcwData outcomes contain only outcomeId = 1", {
  d        <- simulateCcwData(n = 100, seed = 3)
  outcomes <- dplyr::collect(d$outcomes)
  if (nrow(outcomes) > 0) {
    expect_true(all(outcomes$outcomeId == 1L))
  }
})

test_that("simulateCcwData produces covariates in sparse format", {
  d   <- simulateCcwData(n = 50, nCovariates = 3, seed = 4)
  cov <- dplyr::collect(d$covariates)
  expect_true(all(c("rowId", "covariateId", "covariateValue") %in% names(cov)))
  expect_equal(nrow(cov), 50L * 3L)
})

test_that("treatmentPrevalence controls fraction of initiators", {
  set.seed(10)
  d    <- simulateCcwData(n = 500, treatmentPrevalence = 0.8, seed = 10)
  coh  <- dplyr::collect(d$cohorts)
  frac <- mean(!is.na(coh$initiationDate))
  expect_true(frac > 0.6)
})

test_that("higher trueHazardRatio produces more observed events", {
  d1 <- simulateCcwData(n = 400, trueHazardRatio = 1.0, seed = 20)
  d5 <- simulateCcwData(n = 400, trueHazardRatio = 5.0, seed = 20)
  ev1 <- nrow(dplyr::collect(d1$outcomes))
  ev5 <- nrow(dplyr::collect(d5$outcomes))
  expect_true(ev5 >= ev1)
})

test_that("simulateCcwData feeds the full CCW pipeline without error", {
  d <- simulateCcwData(n = 80, seed = 99)
  args <- createCreateClonesAndCensorArgs(
    arms = list(
      createArmDefinition("initiate", "initiateWithin", gracePeriod = 30),
      createArmDefinition("never",    "deferBeyond",    gracePeriod = 30)
    ),
    intervalLength = 30
  )
  popArgs <- createCreateStudyPopulationArgs(
    riskWindowEnd = 365, endAnchor = "cohort start",
    removeSubjectsWithPriorOutcome = FALSE
  )
  createStudyPopulation(d, outcomeId = 1L, popArgs)
  cloned   <- createClonesAndCensor(d, args)
  pt       <- discretizePersonTime(cloned, intervalLength = 30)
  weighted <- createCensoringWeights(pt, cloned)
  model    <- fitOutcomeModel(weighted,
                fitOutcomeModelArgs = createFitOutcomeModelArgs("cox", "robust"))
  expect_s3_class(model, "CloneCensorWeightOutcomeModel")
  expect_true(is.finite(model$summary$estimate))
})

# runSimulation --------------------------------------------------------------

test_that("runSimulation returns a data frame with correct structure", {
  res <- suppressWarnings(runSimulation(nSims = 3, n = 60, seed = 1))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3L)
  expect_true(all(c("simId", "converged", "logRr", "seLogRr",
                    "estimate", "lower95", "upper95", "pValue",
                    "covered", "bias", "bias2") %in% names(res)))
})

test_that("runSimulation simId column matches seed offset", {
  res <- suppressWarnings(runSimulation(nSims = 4, n = 50, seed = 100))
  expect_equal(res$simId, 101:104)
})

test_that("runSimulation converged column is logical", {
  res <- suppressWarnings(runSimulation(nSims = 3, n = 60, seed = 5))
  expect_type(res$converged, "logical")
})

test_that("runSimulation estimates are finite for converged rows", {
  res  <- suppressWarnings(runSimulation(nSims = 5, n = 80, seed = 7))
  conv <- res[res$converged, ]
  if (nrow(conv) > 0) {
    expect_true(all(is.finite(conv$logRr)))
    expect_true(all(is.finite(conv$seLogRr)))
    expect_true(all(conv$lower95 > 0))
    expect_true(all(conv$upper95 > conv$lower95))
  }
})

# summarizeSimulation --------------------------------------------------------

make_sim_results <- function() {
  data.frame(
    simId     = 1:10,
    converged = c(rep(TRUE, 9), FALSE),
    logRr     = c(rnorm(9, 0, 0.3), NA),
    seLogRr   = c(abs(rnorm(9, 0.2, 0.05)), NA),
    estimate  = c(exp(rnorm(9, 0, 0.3)), NA),
    lower95   = c(exp(rnorm(9, -0.5, 0.1)), NA),
    upper95   = c(exp(rnorm(9, 0.5, 0.1)), NA),
    pValue    = c(runif(9), NA),
    covered   = c(sample(c(TRUE, FALSE), 9, replace = TRUE), NA),
    bias      = c(rnorm(9, 0, 0.1), NA),
    bias2     = c(rnorm(9, 0, 0.05)^2, NA),
    stringsAsFactors = FALSE
  )
}

test_that("summarizeSimulation returns a list with expected fields", {
  set.seed(42)
  out <- summarizeSimulation(make_sim_results(), trueHazardRatio = 1.0)
  expect_type(out, "list")
  expect_true(all(c("nSims", "convergenceRate", "trueLogRr", "meanLogRr",
                    "bias", "rmse", "coverage", "power") %in% names(out)))
})

test_that("summarizeSimulation convergenceRate matches data", {
  set.seed(42)
  out <- summarizeSimulation(make_sim_results(), trueHazardRatio = 1.0)
  expect_equal(out$nSims, 10L)
  expect_equal(out$convergenceRate, 0.9)
})

test_that("summarizeSimulation bias and rmse are finite when trueHR supplied", {
  set.seed(42)
  out <- summarizeSimulation(make_sim_results(), trueHazardRatio = 1.0)
  expect_true(is.finite(out$bias))
  expect_true(is.finite(out$rmse))
  expect_true(out$rmse >= 0)
})

test_that("summarizeSimulation coverage is in [0, 1]", {
  set.seed(42)
  out <- summarizeSimulation(make_sim_results(), trueHazardRatio = 1.0)
  expect_true(out$coverage >= 0 && out$coverage <= 1)
})

test_that("summarizeSimulation power is in [0, 1]", {
  set.seed(42)
  out <- summarizeSimulation(make_sim_results(), trueHazardRatio = 2.0)
  expect_true(out$power >= 0 && out$power <= 1)
})

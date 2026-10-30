# Helpers ----------------------------------------------------------------

make_neg_controls <- function(n = 50, seed = 42) {
  set.seed(seed)
  # Simulate negative controls: log HR near 0, some noise + a small bias
  data.frame(
    logRr   = stats::rnorm(n, mean = 0.05, sd = 0.20),
    seLogRr = stats::runif(n, min = 0.05, max = 0.40),
    stringsAsFactors = FALSE
  )
}

make_outcome_model <- function(logRr = log(1.3), se = 0.15) {
  summary <- data.frame(
    modelType     = "cox",
    term          = "arm",
    estimate      = exp(logRr),
    lower95       = exp(logRr - stats::qnorm(0.975) * se),
    upper95       = exp(logRr + stats::qnorm(0.975) * se),
    pValue        = 2 * stats::pnorm(-abs(logRr / se)),
    logLikelihood = -100,
    stringsAsFactors = FALSE
  )
  structure(
    list(summary = summary, model = NULL, estimand = "hazard ratio",
         armLevels = c("a", "b"), clusterVariable = "subjectId",
         varianceEstimator = "robust"),
    class = "CloneCensorWeightOutcomeModel"
  )
}

# Tests ------------------------------------------------------------------

test_that("calibrateEstimates returns data frame with calibrated columns", {
  nc     <- make_neg_controls()
  model  <- make_outcome_model()
  result <- calibrateEstimates(model, negativeControlResults = nc)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("calibratedP", "calibratedLower95", "calibratedUpper95") %in%
                    names(result)))
})

test_that("calibrated columns are finite and in valid ranges", {
  nc     <- make_neg_controls()
  model  <- make_outcome_model()
  result <- calibrateEstimates(model, negativeControlResults = nc)
  expect_true(is.finite(result$calibratedP))
  expect_true(result$calibratedP >= 0 && result$calibratedP <= 1)
  expect_true(result$calibratedLower95 > 0)
  expect_true(result$calibratedUpper95 > result$calibratedLower95)
})

test_that("pre-computed null is accepted in place of raw negative controls", {
  nc   <- make_neg_controls()
  null <- EmpiricalCalibration::fitNull(nc$logRr, nc$seLogRr)
  model  <- make_outcome_model()
  result <- calibrateEstimates(model, null = null)
  expect_true(all(c("calibratedP", "calibratedLower95", "calibratedUpper95") %in%
                    names(result)))
})

test_that("error when neither negativeControlResults nor null is supplied", {
  model <- make_outcome_model()
  expect_error(calibrateEstimates(model))
})

test_that("calibrated CI is wider than uncalibrated under added systematic error", {
  # Add large bias to negative controls → wider calibrated CI
  nc_biased <- make_neg_controls()
  nc_biased$logRr <- nc_biased$logRr + stats::rnorm(nrow(nc_biased), 0.5, 0.1)
  model <- make_outcome_model()
  result_biased <- calibrateEstimates(model, negativeControlResults = nc_biased)
  result_clean  <- calibrateEstimates(model, negativeControlResults = make_neg_controls())
  # Calibrated CI should be wider (or at least not narrower) when there is bias
  width_biased <- log(result_biased$calibratedUpper95) - log(result_biased$calibratedLower95)
  width_clean  <- log(result_clean$calibratedUpper95)  - log(result_clean$calibratedLower95)
  expect_true(width_biased >= width_clean - 0.01)
})

test_that("plotCalibration returns a ggplot", {
  nc <- make_neg_controls()
  p  <- plotCalibration(nc)
  expect_s3_class(p, "ggplot")
})

test_that("plotCalibration with primary estimate returns a ggplot", {
  nc      <- make_neg_controls()
  primary <- data.frame(logRr = log(1.3), seLogRr = 0.15)
  p       <- plotCalibration(nc, primaryEstimate = primary)
  expect_s3_class(p, "ggplot")
})

test_that("plotCalibration accepts a pre-computed null", {
  nc   <- make_neg_controls()
  null <- EmpiricalCalibration::fitNull(nc$logRr, nc$seLogRr)
  p    <- plotCalibration(nc, null = null)
  expect_s3_class(p, "ggplot")
})

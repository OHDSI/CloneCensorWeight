library(testthat)

# pooledLogisticEngine: the pooled logistic hazard model fit by Cyclops's pooled logistic regression
# ("plr") is the one fit on the expanded rows ("lr").

# Nobody initiates on day 3, so the indicator basis leaves that day out (explicit plr rows).
fixture <- createM3Fixture(50, seed = 3)
fixture$persons$treatmentDay[fixture$persons$treatmentDay %in% 3L] <- 4L
data <- createM3Data(fixture)
clones <- buildM3Clones(fixture$persons, list(
  createTreatmentStrategy(1, "Defer", "notInitiateWithin", graceEnd = 7),
  createTreatmentStrategy(2, "Defer too", "notInitiateWithin", graceEnd = 7)
))

test_that("plr and lr fit the same hazard model", {
  skip_if_not(.cyclopsHasPooledLogistic(), "The loaded Cyclops has no plr")
  fitBoth <- function(args) {
    lapply(c(plr = "plr", lr = "lr"), function(engine) {
      args$pooledLogisticEngine <- engine
      suppressMessages(fitCensoringModels(clones, data, args))[[1]]
    })
  }
  # Indicator basis, a day left out (explicit rows), day-0 interactions (an extra fixed basis column):
  args <- m3Args("logistic", prior = Cyclops::createPrior("laplace", variance = 0.5, exclude = 0), covariateTimeInteraction = "day0")
  models <- fitBoth(args)
  expect_equal(models$plr$diagnostics$engine, "plr")
  expect_true(any(grepl(":day0$", names(models$plr$coefficients))))
  expect_equal(models$plr$coefficients, models$lr$coefficients[names(models$plr$coefficients)], tolerance = 1e-6)
  # Spline basis without covariates (implicit rows):
  models <- fitBoth(m3Args("logistic", useCovariates = FALSE, timeBasis = createTimeBasis("naturalSpline", intervalLength = 2, df = 2)))
  expect_equal(models$plr$coefficients, models$lr$coefficients, tolerance = 1e-6)
})

test_that("plr stops before any fit when the loaded Cyclops has no plr", {
  local_mocked_bindings(.cyclopsHasPooledLogistic = function() FALSE)
  args <- m3Args("logistic", useCovariates = FALSE)
  args$pooledLogisticEngine <- "plr"
  expect_error(fitCensoringModels(clones, data, args), "pooledLogisticEngine = \"lr\"")
})

close(clones)
close(data)

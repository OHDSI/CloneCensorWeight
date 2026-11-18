library(testthat)

# pooledLogisticEngine: the pooled logistic hazard model fit by Cyclops's pooled logistic regression
# ("plr") is the one fit on the expanded rows ("lr").

# Two strategies that share one hazard model. Nobody initiates on day 3, so the indicator basis leaves
# that day out and the persons get one plr row per day (explicit rows).
fixture <- createM3Fixture(100, seed = 3)
fixture$persons$treatmentDay[fixture$persons$treatmentDay %in% 3L] <- 4L
data <- createM3Data(fixture)
clones <- buildM3Clones(fixture$persons, list(
  createTreatmentStrategy(1, "Defer", "notInitiateWithin", graceEnd = 7),
  createTreatmentStrategy(2, "Defer too", "notInitiateWithin", graceEnd = 7)
))

# The hazard models of both engines; the order of the interaction coefficients of "lr" depends on the
# order of the covariates in Andromeda.
expectSameFit <- function(args) {
  args$pooledLogisticEngine <- "plr"
  plr <- suppressMessages(fitCensoringModels(clones, data, args))
  args$pooledLogisticEngine <- "lr"
  lr <- suppressMessages(fitCensoringModels(clones, data, args))
  expect_equal(c(plr[[1]]$diagnostics$engine, lr[[1]]$diagnostics$engine), c("plr", "lr"))
  expect_equal(plr[[1]]$coefficients, lr[[1]]$coefficients[names(plr[[1]]$coefficients)], tolerance = 1e-6)
  expect_equal(plr[[1]]$numerator$coefficients, lr[[1]]$numerator$coefficients, tolerance = 1e-6)
  expect_equal(plr[[1]]$priorVariance, lr[[1]]$priorVariance, tolerance = 1e-6)
  return(list(plr = plr, lr = lr))
}

test_that("plr and lr give the same model and weights: indicator basis, day-0 interactions, a day left out, a sample", {
  skip_if_not(.cyclopsHasPooledLogistic(), "The loaded Cyclops has no plr")
  args <- m3Args("logistic", prior = Cyclops::createPrior("laplace", variance = 0.5, exclude = 0), covariateTimeInteraction = "day0")
  args$maxCohortSizeForFitting <- 80
  models <- expectSameFit(args)
  expect_true(models$plr[[1]]$diagnostics$sampled)
  expect_true(any(grepl(":day0$", names(models$plr[[1]]$coefficients))))
  weights <- lapply(models, function(models) {
    weighted <- suppressMessages(computeCensoringWeights(clones, models))
    on.exit(close(weighted))
    collectSegments(weighted)$weight
  })
  expect_equal(weights$plr, weights$lr, tolerance = 1e-6)
})

test_that("plr and lr give the same model: spline basis, interactions by interval", {
  skip_if_not(.cyclopsHasPooledLogistic(), "The loaded Cyclops has no plr")
  args <- m3Args("logistic", timeBasis = createTimeBasis("naturalSpline", intervalLength = 2, df = 2), prior = Cyclops::createPrior("laplace", variance = 0.5, exclude = 0), covariateTimeInteraction = "interval")
  models <- expectSameFit(args)
  expect_true(any(grepl(":timeBasis", names(models$plr[[1]]$coefficients))))
})

test_that("plr stops with a clear error when the loaded Cyclops has no plr", {
  local_mocked_bindings(.cyclopsHasPooledLogistic = function() FALSE)
  args <- m3Args("logistic", useCovariates = FALSE)
  args$pooledLogisticEngine <- "plr"
  expect_error(suppressMessages(fitCensoringModels(clones, data, args)), "pooledLogisticEngine = \"lr\"")
})

close(clones)
close(data)

library(testthat)

# maxCohortSizeForFitting for the initiation hazard models: a uniform sample of whole persons.

# Two strategies (the clones need them) that share one hazard model:
subsampleStrategies <- function() {
  return(list(
    createTreatmentStrategy(1, "Defer", "notInitiateWithin", graceEnd = 7),
    createTreatmentStrategy(2, "Defer too", "notInitiateWithin", graceEnd = 7)
  ))
}

# One data set for both tests:
fixture <- createM3Fixture(100, seed = 3)
data <- createM3Data(fixture)
clones <- buildM3Clones(fixture$persons, subsampleStrategies())

test_that("A hazard model whose risk set is not larger than maxCohortSizeForFitting is not sampled", {
  args <- m3Args("cox")
  args$maxCohortSizeForFitting <- 100
  model <- suppressMessages(fitCensoringModels(clones, data, args))[[1]]
  args$maxCohortSizeForFitting <- 0
  expect_identical(suppressMessages(fitCensoringModels(clones, data, args))[[1]], model)
  expect_false(model$diagnostics$sampled)
  expect_equal(model$diagnostics$sampleSize, model$diagnostics$riskSetSize)
})

test_that("A larger hazard risk set is fit on a seeded sample of whole persons and predicts for everyone", {
  args <- m3Args("logistic", prior = Cyclops::createPrior("normal", variance = 1, exclude = c(0)))
  args$maxCohortSizeForFitting <- 30
  messages <- capture_messages(model <- fitCensoringModels(clones, data, args)[[1]])
  expect_true(any(grepl("uniform sample of 30 of its", messages)))
  expect_true(model$diagnostics$sampled)
  expect_equal(model$diagnostics$sampleSize, 30)
  expect_gt(model$diagnostics$riskSetSize, 30)
  expect_equal(nrow(model$predictions), model$diagnostics$riskSetSize)
  expect_false(any(model$intervals$baseline == Inf))
  # Two intervals without initiations in the sample but some in the risk set get the risk set's
  # hazard; only the last interval (none in the risk set either) stays at -Inf:
  expect_equal(sum(is.infinite(model$intervals$baseline)), 1L)
  expect_equal(sum(is.infinite(model$numerator$intervals$baseline)), 1L)

  # The same fit as on the persons drawn with the seed of control alone (whole persons, numerator and
  # denominator alike):
  rowIds <- sort(model$predictions$rowId)
  sampled <- rowIds[sort(withr::with_seed(1, sample.int(length(rowIds), 30)))]
  subset <- fixture
  subset$persons <- fixture$persons[fixture$persons$rowId %in% sampled, ]
  subset$covariates <- fixture$covariates[fixture$covariates$rowId %in% sampled, ]
  subsetData <- createM3Data(subset)
  subsetClones <- buildM3Clones(subset$persons, subsampleStrategies())
  args$maxCohortSizeForFitting <- 0
  reference <- suppressMessages(fitCensoringModels(subsetClones, subsetData, args))[[1]]
  expect_equal(model$coefficients, reference$coefficients, tolerance = 1e-8)
  expect_equal(model$numerator$coefficients, reference$numerator$coefficients, tolerance = 1e-8)
  expect_equal(model$diagnostics$sampleEvents, reference$diagnostics$events)

  # The draw leaves the caller's random number stream as it was:
  set.seed(5)
  suppressMessages(.sampleRiskSet(data.frame(y = 1:5), 2, 1, "model"))
  expect_identical(runif(1), withr::with_seed(5, runif(1)))
  # A sample without initiations, from a risk set with some, is not used:
  expect_message(whole <- .sampleRiskSet(data.frame(y = c(1, rep(0, 9))), 2, 1, "model"), "no initiations")
  expect_equal(nrow(whole), 10)
  close(subsetClones)
  close(subsetData)
})

close(clones)
close(data)

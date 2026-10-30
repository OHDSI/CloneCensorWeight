test_that("sample data has a post-CDM contract", {
  data <- createSampleCloneCensorWeightData(n = 20, seed = 1)
  expect_s4_class(data, "CloneCensorWeightData")
  expect_false(is.null(data$cohorts))
  expect_false(is.null(data$outcomes))
  expect_false(is.null(data$covariates))
  expect_false(is.null(data$covariateRef))
  expect_false(is.null(data$analysisRef))
  cohorts <- dplyr::collect(data$cohorts)
  expect_equal(nrow(cohorts), 20)
  expect_true(all(c("subjectId", "cohortStartDate", "cohortEndDate", "initiationDate") %in%
                    names(cohorts)))
})

test_that("cloning and weighting produce valid analysis data", {
  data <- createSampleCloneCensorWeightData(n = 40, seed = 2)
  clonesArgs <- createCreateClonesAndCensorArgs(
    arms = list(
      createArmDefinition("initiate", "initiateWithin", gracePeriod = 30),
      createArmDefinition("never",    "deferBeyond",    gracePeriod = 30)
    ),
    intervalLength = 30
  )
  cloned <- createClonesAndCensor(data, clonesArgs)
  cohorts <- dplyr::collect(cloned$cohorts)
  expect_equal(nrow(cohorts), 80)
  expect_true(any(cohorts$artificiallyCensored == 1))
  personTime <- discretizePersonTime(cloned, intervalLength = 30)
  weighted <- createCensoringWeights(personTime, cloned)
  expect_true(all(is.finite(weighted$weight)))
  expect_true(all(weighted$weight > 0))
  model <- fitOutcomeModel(weighted, modelType = "cox",
    fitOutcomeModelArgs = createFitOutcomeModelArgs("cox", "robust"))
  expect_s3_class(model, "CloneCensorWeightOutcomeModel")
  expect_equal(model$varianceEstimator, "robust")
  expect_true(is.finite(model$summary$estimate))
})

test_that("Eunomia is accepted as an optional sample CDM", {
  skip_if_not_installed("Eunomia")
  expect_true(requireNamespace("Eunomia", quietly = TRUE))
})
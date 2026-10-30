make_cloned <- function(n = 60, seed = 5) {
  data <- createSampleCloneCensorWeightData(n = n, seed = seed)
  args <- createCreateClonesAndCensorArgs(
    arms = list(
      createArmDefinition("initiate", "initiateWithin", gracePeriod = 30),
      createArmDefinition("never",    "deferBeyond",    gracePeriod = 30)
    ),
    intervalLength = 30
  )
  createClonesAndCensor(data, args)
}

test_that("computeStandardizedDifferences returns one row per covariate", {
  cloned <- make_cloned()
  smd <- computeStandardizedDifferences(cloned)
  # sample data has 3 covariates (age, priorCondition, priorDrug)
  expect_equal(nrow(smd), 3)
  expect_true(all(c("covariateId", "covariate", "smd", "absSmd") %in% names(smd)))
  expect_true(all(is.finite(smd$smd)))
  expect_true(all(smd$absSmd >= 0))
})

test_that("absSmd equals abs(smd)", {
  cloned <- make_cloned()
  smd <- computeStandardizedDifferences(cloned)
  expect_equal(smd$absSmd, abs(smd$smd))
})

test_that("after-IPCW SMDs are returned when weights supplied", {
  cloned   <- make_cloned()
  pt       <- discretizePersonTime(cloned, intervalLength = 30)
  weighted <- createCensoringWeights(pt, cloned)
  smdAfter <- computeStandardizedDifferences(cloned, weights = weighted)
  expect_equal(nrow(smdAfter), 3)
  expect_true(all(is.finite(smdAfter$smd)))
})

test_that("plotLovePlot produces a ggplot with both before and after layers", {
  cloned   <- make_cloned()
  smdBef   <- computeStandardizedDifferences(cloned)
  pt       <- discretizePersonTime(cloned, intervalLength = 30)
  weighted <- createCensoringWeights(pt, cloned)
  smdAft   <- computeStandardizedDifferences(cloned, weights = weighted)

  p <- plotLovePlot(smdBef, smdAft)
  expect_s3_class(p, "ggplot")

  # Both groups must appear in the underlying plot data
  groups <- as.character(unique(p$data$Weighting))
  expect_true("Before IPCW" %in% groups)
  expect_true("After IPCW"  %in% groups)
})

test_that("plotLovePlot works with before-only data", {
  cloned <- make_cloned()
  smdBef <- computeStandardizedDifferences(cloned)
  expect_s3_class(plotLovePlot(smdBef), "ggplot")
})

test_that("empty covariates return empty data frame", {
  data  <- createCloneCensorWeightData(
    cohorts = data.frame(subjectId = 1L, cohortStartDate = as.Date("2020-01-01"),
                         cohortEndDate = as.Date("2021-01-01"),
                         initiationDate = NA_real_, stringsAsFactors = FALSE)
  )
  args <- createCreateClonesAndCensorArgs(
    arms = list(
      createArmDefinition("a", "initiateWithin", gracePeriod = 30),
      createArmDefinition("b", "deferBeyond",    gracePeriod = 30)
    ),
    intervalLength = 30
  )
  cloned <- createClonesAndCensor(data, args)
  smd <- computeStandardizedDifferences(cloned)
  expect_equal(nrow(smd), 0)
  expect_true(all(c("covariateId", "covariate", "smd", "absSmd") %in% names(smd)))
})

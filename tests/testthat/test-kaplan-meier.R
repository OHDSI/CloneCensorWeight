make_weighted <- function(n = 60, seed = 9) {
  data <- createSampleCloneCensorWeightData(n = n, seed = seed)
  args <- createCreateClonesAndCensorArgs(
    arms = list(
      createArmDefinition("initiate", "initiateWithin", gracePeriod = 30),
      createArmDefinition("never",    "deferBeyond",    gracePeriod = 30)
    ),
    intervalLength = 30
  )
  cloned   <- createClonesAndCensor(data, args)
  personTime <- discretizePersonTime(cloned, intervalLength = 30)
  createCensoringWeights(personTime, cloned)
}

test_that("computeKaplanMeier returns correct structure", {
  weighted <- make_weighted()
  km <- computeKaplanMeier(weighted)
  expect_s3_class(km, "data.frame")
  expect_true(all(c("arm", "time", "survival", "lower", "upper", "nRisk", "nEvent") %in%
                    names(km)))
  expect_equal(sort(unique(km$arm)), c("initiate", "never"))
})

test_that("survival starts at 1 and stays in [0, 1]", {
  weighted <- make_weighted()
  km <- computeKaplanMeier(weighted)
  # origin rows
  origin <- km[km$time == 0, ]
  expect_equal(nrow(origin), 2)  # one per arm
  expect_true(all(origin$survival == 1))
  # all values in range
  expect_true(all(km$survival >= 0 & km$survival <= 1))
  expect_true(all(km$lower >= 0 & km$lower <= 1, na.rm = TRUE))
  expect_true(all(km$upper >= 0 & km$upper <= 1, na.rm = TRUE))
})

test_that("survival is non-increasing within each arm", {
  weighted <- make_weighted()
  km <- computeKaplanMeier(weighted)
  for (a in unique(km$arm)) {
    sub <- km[km$arm == a, ]
    sub <- sub[order(sub$time), ]
    expect_true(all(diff(sub$survival) <= 1e-9),
                info = paste("non-increasing check failed for arm:", a))
  }
})

test_that("CI lower <= survival <= upper", {
  weighted <- make_weighted()
  km <- computeKaplanMeier(weighted)
  # exclude origin row (lower/upper = 1 there)
  km2 <- km[km$time > 0, ]
  expect_true(all(km2$lower <= km2$survival + 1e-9))
  expect_true(all(km2$survival <= km2$upper + 1e-9))
})

test_that("plotKaplanMeier returns a ggplot", {
  weighted <- make_weighted()
  p <- plotKaplanMeier(weighted)
  expect_s3_class(p, "ggplot")
  # Plot data has both arms
  arms <- unique(p$data$arm)
  expect_true(length(arms) == 2)
})

test_that("no-event data returns flat survival of 1", {
  # All subjects have no event
  cohorts <- data.frame(
    subjectId = 1:10, cohortStartDate = as.Date("2020-01-01"),
    cohortEndDate = as.Date("2020-04-10"),
    initiationDate = as.Date(rep(c("2020-01-15", NA), 5)),
    stringsAsFactors = FALSE
  )
  data   <- createCloneCensorWeightData(cohorts)
  args   <- createCreateClonesAndCensorArgs(
    arms = list(
      createArmDefinition("a", "initiateWithin", gracePeriod = 30),
      createArmDefinition("b", "deferBeyond",    gracePeriod = 30)
    ),
    intervalLength = 30
  )
  cloned <- createClonesAndCensor(data, args)
  pt     <- discretizePersonTime(cloned, intervalLength = 30)
  wt     <- createCensoringWeights(pt, cloned)
  km     <- computeKaplanMeier(wt)
  expect_true(all(km$survival == 1))
})

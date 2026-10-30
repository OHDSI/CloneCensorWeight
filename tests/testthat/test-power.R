make_weighted_power <- function(n = 80, seed = 13) {
  data <- createSampleCloneCensorWeightData(n = n, seed = seed)
  args <- createCreateClonesAndCensorArgs(
    arms = list(
      createArmDefinition("initiate", "initiateWithin", gracePeriod = 30),
      createArmDefinition("never",    "deferBeyond",    gracePeriod = 30)
    ),
    intervalLength = 30
  )
  cloned   <- createClonesAndCensor(data, args)
  pt       <- discretizePersonTime(cloned, intervalLength = 30)
  createCensoringWeights(pt, cloned)
}

test_that("computeMdrr returns correct structure", {
  weighted <- make_weighted_power()
  result   <- computeMdrr(weighted)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true(all(c("arm1", "arm2", "alpha", "power", "nEventsArm1", "nEventsArm2",
                    "nEventsEffArm1", "nEventsEffArm2", "mdrr") %in% names(result)))
})

test_that("MDRR is >= 1 when events are present", {
  weighted <- make_weighted_power()
  result   <- computeMdrr(weighted)
  if (is.finite(result$mdrr)) {
    expect_true(result$mdrr >= 1)
  }
})

test_that("MDRR is Inf when an arm has no events", {
  # Force zero events in one arm by using very short follow-up
  cohorts <- data.frame(
    subjectId = 1:20, cohortStartDate = as.Date("2020-01-01"),
    cohortEndDate = as.Date("2020-02-01"),
    initiationDate = as.Date(rep(c("2020-01-10", NA), 10)),
    stringsAsFactors = FALSE
  )
  # No outcomes → zero events in both arms
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
  result <- computeMdrr(wt)
  expect_true(is.infinite(result$mdrr))
})

test_that("more events → smaller MDRR (more power)", {
  # Compare MDRR for n=40 vs n=200; more events should give lower MDRR
  r40  <- computeMdrr(make_weighted_power(n = 40,  seed = 1))
  r200 <- computeMdrr(make_weighted_power(n = 200, seed = 1))
  if (is.finite(r40$mdrr) && is.finite(r200$mdrr)) {
    expect_true(r200$mdrr <= r40$mdrr)
  }
})

test_that("alpha and power parameters are respected", {
  weighted <- make_weighted_power()
  r80  <- computeMdrr(weighted, power = 0.80)
  r90  <- computeMdrr(weighted, power = 0.90)
  # Higher power requirement → larger MDRR needed
  if (is.finite(r80$mdrr) && is.finite(r90$mdrr)) {
    expect_true(r90$mdrr >= r80$mdrr)
  }
  expect_equal(r80$alpha, 0.05)
  expect_equal(r80$power, 0.80)
})

test_that("plotMdrr returns a ggplot", {
  weighted <- make_weighted_power()
  result   <- computeMdrr(weighted)
  p <- plotMdrr(result)
  expect_s3_class(p, "ggplot")
})

test_that("plotMdrr works with multiple stacked results", {
  r1 <- computeMdrr(make_weighted_power(n = 60,  seed = 1))
  r2 <- computeMdrr(make_weighted_power(n = 120, seed = 2))
  r1$analysisId <- "A"
  r2$analysisId <- "B"
  p <- plotMdrr(rbind(r1, r2))
  expect_s3_class(p, "ggplot")
  expect_equal(length(unique(p$data$analysisId)), 2)
})

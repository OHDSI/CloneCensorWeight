library(testthat)

# Tests of the weights of the outcome MSMs (Phase1Contracts section 11.1): createFitOutcomeModelArgs(
# msmWeights = ), the stabilization numerator computeCensoringWeights() stores on every segment, and the
# Cox, pooled logistic, Kaplan-Meier and bootstrap paths that use them.

# Simulated data, clones and censoring models (with covariates and no prior, so the weights vary by
# person) of initiateWithin(5) vs `comparatorType`. The caller closes `data` and `clones`.
msmAnalysis <- function(comparatorType = "neverInitiate", sameDayOrdering = "eventsFirst", hazardModelType = "logistic", postGraceModel = "binary", seed = 1) {
  data <- simulateCloneCensorWeightData(
    n = 300, maxDays = 20, graceEnd = 5, initiationHazard = 0.05, outcomeHazard = 0.01, deathHazard = 0.005, seed = seed
  )
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 20, endAnchor = "cohort start")
  ))
  strategies <- list(
    createTreatmentStrategy(1, "Initiate", "initiateWithin", graceEnd = 5),
    createTreatmentStrategy(2, "Comparator", comparatorType, graceEnd = if (comparatorType == "neverInitiate") NULL else 5)
  )
  clones <- suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = strategies, sameDayOrdering = sameDayOrdering)))
  models <- suppressMessages(fitCensoringModels(clones, data, createFitCensoringModelsArgs(
    hazardModelType = hazardModelType, prior = Cyclops::createPrior("none"), postGraceModel = postGraceModel,
    # The hand-computed references have no covariate-by-time interactions:
    covariateTimeInteraction = "none"
  )))
  return(list(data = data, clones = clones, models = models))
}

closeMsmAnalysis <- function(analysis) {
  close(analysis$clones)
  close(analysis$data)
}

# truncateMsmWeights = FALSE: these tests check the section 11.1 behaviour, which it reproduces (section
# 12.2); TRUE is tested in test-integration-designChoices.R.
msmComputeWeightsArgs <- function(stabilized, truncation = createWeightTruncation("none"), truncateMsmWeights = FALSE) {
  return(createComputeWeightsArgs(stabilized = stabilized, truncation = truncation, truncateMsmWeights = truncateMsmWeights))
}

weightMsm <- function(analysis, computeWeightsArgs, outcomeTimeBasis = NULL) {
  return(suppressMessages(computeCensoringWeights(analysis$clones, analysis$models, computeWeightsArgs, outcomeTimeBasis = outcomeTimeBasis)))
}

collectMsmSegments <- function(weightedClones) {
  segments <- as.data.frame(dplyr::collect(weightedClones$segments))
  segments <- segments[order(segments$cloneRowId, segments$startTime), ]
  rownames(segments) <- NULL
  segments$x <- as.integer(segments$strategyId == 1)
  return(segments)
}

fitMsm <- function(weightedClones, ...) {
  return(withoutGracePeriodWarnings(suppressMessages(fitOutcomeModel(weightedClones, createFitOutcomeModelArgs(...)))))
}

msmCoxControl <- function() survival::coxph.control(eps = 1e-12, toler.chol = 1e-14, iter.max = 100)

# The unstabilized W(d) of each segment (d its last day) and its stabilization numerator, written out
# from the fitted censoring models ("eventsFirst", daily pooled logistic hazard model; contract section
# 5.4): 1 up to day G and 1 / pi after it, with numerator pi-bar after G (binary); 1 / prod over days
# k < d, k <= H of (1 - h_k), with the same product of the numerator model as numerator (hazard), except
# after G for a hazard strategy with the binary model at G (postGraceModel = "binary", contract section
# 11.2): 1 / (1 - pi), with numerator 1 - pi-bar.
handMsmWeights <- function(models, segments, riskStart) {
  day <- segments$stopTime + riskStart - 1
  weight <- rep(NA_real_, nrow(segments))
  numerator <- rep(NA_real_, nrow(segments))
  for (strategyId in 1:2) {
    model <- models[[strategyId]]
    idx <- which(segments$strategyId == strategyId)
    if (model$kind == "binary") {
      pi <- model$predictions$propensityScore[match(segments$rowId[idx], model$predictions$rowId)]
      after <- day[idx] > model$horizon
      weight[idx] <- ifelse(after, 1 / pi, 1)
      numerator[idx] <- ifelse(after, model$numerator$propensityScore, 1)
    } else {
      linearPredictor <- model$predictions$linearPredictor[match(segments$rowId[idx], model$predictions$rowId)]
      baseline <- model$intervals$baseline[match(seq_len(model$horizon + 1) - 1, model$intervals$k)]
      numeratorBaseline <- model$numerator$intervals$baseline[match(seq_len(model$horizon + 1) - 1, model$numerator$intervals$k)]
      for (j in seq_along(idx)) {
        k <- seq_len(min(day[idx[j]], model$horizon + 1))
        weight[idx[j]] <- 1 / prod(1 - stats::plogis(baseline[k] + linearPredictor[j]))
        numerator[idx[j]] <- prod(1 - stats::plogis(numeratorBaseline[k]))
      }
      if (!is.null(model$postGrace)) {
        after <- idx[day[idx] > model$horizon]
        pi <- model$postGrace$predictions$propensityScore[match(segments$rowId[after], model$postGrace$predictions$rowId)]
        weight[after] <- 1 / (1 - pi)
        numerator[after] <- 1 - model$postGrace$numerator$propensityScore
      }
    }
  }
  return(list(weight = weight, numerator = numerator))
}

test_that("The Cox MSM uses the msmWeights weights: unstabilized and stabilized both match coxph", {
  skip_if_not_installed("survival")
  analysis <- msmAnalysis("neverInitiate")
  on.exit(closeMsmAnalysis(analysis))
  # The weights of computeCensoringWeights() are stabilized (its default):
  weighted <- weightMsm(analysis, msmComputeWeightsArgs(stabilized = TRUE))
  on.exit(close(weighted), add = TRUE)
  segments <- collectMsmSegments(weighted)
  hand <- handMsmWeights(analysis$models, segments, riskStart = 0)
  expect_equal(segments$stabilizationNumerator, hand$numerator, tolerance = 1e-10)
  expect_equal(segments$weight, hand$weight * hand$numerator, tolerance = 1e-10)
  segments$unstabilizedWeight <- hand$weight

  # (1) The default, msmWeights = "unstabilized": coxph with the unstabilized weights formed by hand.
  unstabilizedReference <- survival::coxph(
    survival::Surv(startTime, stopTime, outcome) ~ x + cluster(personSeqId),
    data = segments, weights = unstabilizedWeight, ties = "breslow", control = msmCoxControl()
  )
  # (2) msmWeights = "stabilized" is the earlier behaviour: coxph with the stored (stabilized) weights.
  stabilizedReference <- survival::coxph(
    survival::Surv(startTime, stopTime, outcome) ~ x + cluster(personSeqId),
    data = segments, weights = weight, ties = "breslow", control = msmCoxControl()
  )
  for (msmWeights in c("unstabilized", "stabilized")) {
    reference <- if (msmWeights == "unstabilized") unstabilizedReference else stabilizedReference
    args <- if (msmWeights == "unstabilized") list() else list(msmWeights = msmWeights)
    robust <- do.call(fitMsm, c(list(weighted, modelType = "cox", horizons = c(10, 20)), args))
    modelBased <- do.call(fitMsm, c(list(weighted, modelType = "cox", horizons = c(10, 20), varianceEstimator = "modelBased"), args))
    expect_equal(robust$metaData$args$msmWeights, msmWeights)
    expect_equal(robust$estimates$logEstimate[1], unname(coef(reference)), tolerance = 1e-8, info = msmWeights)
    expect_equal(robust$estimates$seLogEstimate[1], sqrt(reference$var[1, 1]), tolerance = 1e-8, info = msmWeights)
    expect_equal(unname(robust$vcov), unname(reference$var), tolerance = 1e-8, info = msmWeights)
    expect_equal(modelBased$estimates$logEstimate[1], unname(coef(reference)), tolerance = 1e-8, info = msmWeights)
    expect_equal(modelBased$estimates$seLogEstimate[1], sqrt(reference$naive.var[1, 1]), tolerance = 1e-8, info = msmWeights)
  }
  # The stabilization numerators are arm- and day-specific, so the two hazard ratios differ:
  expect_gt(abs(coef(unstabilizedReference) - coef(stabilizedReference)), 0.05)
})

test_that("The Kaplan-Meier and Aalen-Johansen rows do not depend on msmWeights or on the stabilization", {
  analysis <- msmAnalysis("notInitiateWithin")
  on.exit(closeMsmAnalysis(analysis))
  horizons <- c(3, 10, 20)
  for (truncation in list(createWeightTruncation("none"), createWeightTruncation("percentile", 0.9))) {
    weighted <- weightMsm(analysis, msmComputeWeightsArgs(stabilized = TRUE, truncation = truncation))
    km <- fitMsm(weighted, modelType = "km", horizons = horizons)
    for (variance in c("robust", "modelBased")) {
      unstabilized <- fitMsm(weighted, modelType = "cox", horizons = horizons, varianceEstimator = variance)
      stabilized <- fitMsm(weighted, modelType = "cox", horizons = horizons, varianceEstimator = variance, msmWeights = "stabilized")
      expect_false(isTRUE(all.equal(unstabilized$estimates[1, ], stabilized$estimates[1, ])))
      # (3) The rows at the horizons are identical under both, and (robust) identical to "km":
      expect_identical(unstabilized$estimates[-1, ], stabilized$estimates[-1, ])
      if (variance == "robust") {
        expect_identical(unstabilized$estimates[-1, ], km$estimates)
      }
    }
    close(weighted)
  }
  # Without truncation the stabilization cancels in the Kaplan-Meier and Aalen-Johansen estimates:
  kmEstimates <- lapply(c(TRUE, FALSE), function(stabilized) {
    weighted <- weightMsm(analysis, msmComputeWeightsArgs(stabilized = stabilized))
    on.exit(close(weighted))
    fitMsm(weighted, modelType = "km", horizons = horizons)$estimates
  })
  expect_equal(kmEstimates[[2]], kmEstimates[[1]], tolerance = 1e-10)
})

test_that("The pooled logistic MSM uses the msmWeights weights, as glm does", {
  skip_if_not_installed("sandwich")
  analysis <- msmAnalysis("neverInitiate")
  on.exit(closeMsmAnalysis(analysis))
  timeBasis <- createTimeBasis("linear", intervalLength = 1)
  weighted <- weightMsm(analysis, msmComputeWeightsArgs(stabilized = TRUE), outcomeTimeBasis = timeBasis)
  on.exit(close(weighted), add = TRUE)
  # Daily segments are the clone x day rows, each with the weight of its day:
  rows <- collectMsmSegments(weighted)
  expect_true(all(rows$stopTime - rows$startTime == 1))
  rows$k <- rows$stopTime - 1
  noDeath <- rows$competing == 0
  glmControl <- glm.control(epsilon = 1e-14, maxit = 100)
  odds <- c()
  for (msmWeights in c("unstabilized", "stabilized")) {
    rows$w <- if (msmWeights == "unstabilized") rows$weight / rows$stabilizationNumerator else rows$weight
    model <- fitMsm(weighted, modelType = "logistic", horizons = 20, timeBasis = timeBasis, msmWeights = msmWeights)
    # (4) Young et al. (2020): the outcome model on the rows without death, the death model on all rows.
    outcomeFit <- glm(outcome ~ x + k, family = quasibinomial(), data = rows[noDeath, ], weights = w, control = glmControl)
    deathFit <- glm(competing ~ x + k, family = quasibinomial(), data = rows, weights = w, control = glmControl)
    expect_equal(unname(model$coefficients), unname(c(coef(outcomeFit), coef(deathFit))), tolerance = 1e-8, info = msmWeights)
    robust <- sandwich::vcovCL(outcomeFit, cluster = rows$personSeqId[noDeath], type = "HC0", cadjust = FALSE)
    expect_equal(model$vcov["target", "target"], robust["x", "x"], tolerance = 1e-6, info = msmWeights)
    modelBased <- fitMsm(weighted, modelType = "logistic", horizons = 20, timeBasis = timeBasis, msmWeights = msmWeights, varianceEstimator = "modelBased")
    information <- crossprod(model.matrix(outcomeFit), model.matrix(outcomeFit) * (outcomeFit$prior.weights * fitted(outcomeFit) * (1 - fitted(outcomeFit))))
    expect_equal(modelBased$vcov["target", "target"], solve(information)[2, 2], tolerance = 1e-6, info = msmWeights)
    odds[msmWeights] <- model$coefficients[["target"]]
  }
  expect_gt(abs(odds[["unstabilized"]] - odds[["stabilized"]]), 0.05)
})

test_that("The stored numerator reconstructs the other weight for every strategy type, with and without truncation", {
  skip_if_not_installed("survival")
  configurations <- list(
    list(comparatorType = "notInitiateWithin", sameDayOrdering = "eventsFirst", hazardModelType = "logistic"),
    list(comparatorType = "neverInitiate", sameDayOrdering = "deviationFirst", hazardModelType = "cox")
  )
  structure <- c("segmentId", "cloneRowId", "strategyId", "startTime", "stopTime", "outcome", "competing")
  for (configuration in configurations) {
    label <- paste(unlist(configuration), collapse = " ")
    analysis <- do.call(msmAnalysis, configuration)
    weighted <- list(
      stabilized = weightMsm(analysis, msmComputeWeightsArgs(stabilized = TRUE)),
      unstabilized = weightMsm(analysis, msmComputeWeightsArgs(stabilized = FALSE))
    )
    s <- collectMsmSegments(weighted$stabilized)
    u <- collectMsmSegments(weighted$unstabilized)
    expect_equal(s[, structure], u[, structure], info = label)
    # The numerator is stored whatever computeWeightsArgs$stabilized is, and varies in both strategies:
    expect_equal(s$stabilizationNumerator, u$stabilizationNumerator, tolerance = 1e-12, info = label)
    expect_true(all(s$stabilizationNumerator > 0 & s$stabilizationNumerator <= 1), info = label)
    expect_true(all(tapply(s$stabilizationNumerator < 1 - 1e-6, s$strategyId, any)), info = label)
    # (5) Each weight is the other times or divided by the numerator:
    expect_equal(s$weight / s$stabilizationNumerator, u$weight, tolerance = 1e-12, info = label)
    expect_equal(u$weight * u$stabilizationNumerator, s$weight, tolerance = 1e-12, info = label)
    # The binary model's numerator is pi-bar after G and 1 up to G:
    initiate <- s$strategyId == 1
    afterG <- s$stopTime - 1 > if (configuration$sameDayOrdering == "eventsFirst") 5 else 4
    expect_equal(unique(s$stabilizationNumerator[initiate & afterG]), analysis$models[[1]]$numerator$propensityScore, info = label)
    expect_equal(unique(s$stabilizationNumerator[initiate & !afterG]), 1, info = label)
    # So the outcome models fitted with either msmWeights are the same from either weighted clones:
    for (msmWeights in c("unstabilized", "stabilized")) {
      fromStabilized <- fitMsm(weighted$stabilized, modelType = "cox", horizons = 20, msmWeights = msmWeights)
      fromUnstabilized <- fitMsm(weighted$unstabilized, modelType = "cox", horizons = 20, msmWeights = msmWeights)
      expect_equal(fromStabilized$estimates, fromUnstabilized$estimates, tolerance = 1e-8, info = paste(label, msmWeights))
    }

    # With truncation, the same segments are truncated under both weights: the other weight is formed
    # from the truncated weight.
    for (stabilized in c(TRUE, FALSE)) {
      truncatedClones <- weightMsm(analysis, msmComputeWeightsArgs(stabilized = stabilized, truncation = createWeightTruncation("percentile", 0.9)))
      truncated <- collectMsmSegments(truncatedClones)
      untruncated <- if (stabilized) s else u
      other <- if (stabilized) u else s
      expect_equal(truncated$stabilizationNumerator, untruncated$stabilizationNumerator, tolerance = 1e-12, info = label)
      isTruncated <- truncated$weight < untruncated$weight * (1 - 1e-12)
      expect_true(any(isTruncated), info = label)
      otherWeight <- if (stabilized) truncated$weight / truncated$stabilizationNumerator else truncated$weight * truncated$stabilizationNumerator
      expect_equal(otherWeight[!isTruncated], other$weight[!isTruncated], tolerance = 1e-12, info = label)
      expect_true(all(otherWeight[isTruncated] < other$weight[isTruncated]), info = label)
      # The Cox MSM with the other weight uses exactly these weights:
      otherMsmWeights <- if (stabilized) "unstabilized" else "stabilized"
      model <- fitMsm(truncatedClones, modelType = "cox", horizons = 20, msmWeights = otherMsmWeights)
      truncated$otherWeight <- otherWeight
      reference <- survival::coxph(
        survival::Surv(startTime, stopTime, outcome) ~ x + cluster(personSeqId),
        data = truncated, weights = otherWeight, ties = "breslow", control = msmCoxControl()
      )
      expect_equal(model$estimates$logEstimate[1], unname(coef(reference)), tolerance = 1e-8, info = label)
      expect_equal(model$estimates$seLogEstimate[1], sqrt(reference$var[1, 1]), tolerance = 1e-8, info = label)
      close(truncatedClones)
    }
    close(weighted$stabilized)
    close(weighted$unstabilized)
    closeMsmAnalysis(analysis)
  }
})

test_that("Stabilized and unstabilized weights reconstruct each other for every strategy type, post-grace model and truncation basis", {
  # Contract sections 11.1, 11.2 and 11.5: after G, the numerator of a notInitiateWithin(G) clone is
  # 1 - pi-bar of the binary model at G (postGraceModel = "binary") or S-bar(G) ("hazard"); either weight
  # is the other times or divided by the numerator, and under truncation (with either basis) the other
  # weight is formed from the truncated weight, so its cap is the cap divided or multiplied by the
  # numerator.
  configurations <- list(
    list(comparatorType = "notInitiateWithin", postGraceModel = "binary"),
    list(comparatorType = "notInitiateWithin", postGraceModel = "hazard"),
    list(comparatorType = "neverInitiate", postGraceModel = "binary"),
    list(comparatorType = "notInitiateWithin", postGraceModel = "binary", sameDayOrdering = "deviationFirst", hazardModelType = "cox")
  )
  truncations <- list(
    none = createWeightTruncation("none"),
    varying = createWeightTruncation("percentile", 0.9, basis = "varying"),
    allClones = createWeightTruncation("percentile", 0.9, basis = "allClones")
  )
  kinds <- c("stabilized", "unstabilized")
  structure <- c("segmentId", "cloneRowId", "rowId", "strategyId", "startTime", "stopTime", "outcome", "competing")
  for (configuration in configurations) {
    label <- paste(unlist(configuration), collapse = " ")
    analysis <- do.call(msmAnalysis, configuration)
    # The segments, caps, the weights each msmWeights gives the MSMs, and (untruncated) the Cox MSMs:
    results <- list()
    for (truncation in names(truncations)) {
      for (kind in kinds) {
        weighted <- weightMsm(analysis, msmComputeWeightsArgs(stabilized = kind == "stabilized", truncation = truncations[[truncation]]))
        outcomeData <- .getOutcomeData(weighted)
        result <- list(
          segments = collectMsmSegments(weighted),
          caps = attr(weighted, "metaData")$weights$caps,
          msm = lapply(stats::setNames(kinds, kinds), function(msmWeights) .withMsmWeights(outcomeData, msmWeights)$segments$weight)
        )
        if (truncation == "none") {
          result$cox <- lapply(stats::setNames(kinds, kinds), function(msmWeights) fitMsm(weighted, modelType = "cox", horizons = 20, msmWeights = msmWeights)$estimates)
        }
        close(weighted)
        results[[truncation]][[kind]] <- result
      }
    }
    closeMsmAnalysis(analysis)

    s <- results$none$stabilized$segments
    u <- results$none$unstabilized$segments
    expect_equal(s[, structure], u[, structure], info = label)
    numerator <- u$stabilizationNumerator
    expect_equal(s$stabilizationNumerator, numerator, tolerance = 1e-12, info = label)
    expect_true(all(numerator > 0 & numerator <= 1), info = label)
    expect_true(all(tapply(numerator < 1 - 1e-6, u$strategyId, any)), info = label)
    # Without truncation, each weight is exactly the other times or divided by the numerator, the MSMs
    # get the same weights from either weighted clones, and so the same fits:
    expect_equal(s$weight, u$weight * numerator, tolerance = 1e-12, info = label)
    expect_equal(u$weight, s$weight / numerator, tolerance = 1e-12, info = label)
    for (kind in kinds) {
      for (msmWeights in kinds) {
        expect_equal(results$none[[kind]]$msm[[msmWeights]], results$none[[msmWeights]]$segments$weight, tolerance = 1e-12, info = paste(label, kind, msmWeights))
      }
    }
    for (msmWeights in kinds) {
      expect_equal(results$none$stabilized$cox[[msmWeights]], results$none$unstabilized$cox[[msmWeights]], tolerance = 1e-8, info = paste(label, msmWeights))
    }

    # The numerators: pi-bar after G and 1 up to G (initiateWithin); after G, 1 - pi-bar of the binary
    # model at G or S-bar(G) (notInitiateWithin); with "eventsFirst" and the daily logistic hazard model,
    # every weight and numerator written out by hand.
    deviationFirst <- identical(configuration$sameDayOrdering, "deviationFirst")
    afterG <- if (deviationFirst) u$stopTime - 1 >= 5 else u$stopTime - 1 > 5
    initiate <- u$strategyId == 1
    expect_equal(unique(numerator[initiate & afterG]), analysis$models[[1]]$numerator$propensityScore, info = label)
    expect_equal(unique(numerator[initiate & !afterG]), 1, info = label)
    postGrace <- analysis$models[[2]]$postGrace
    if (configuration$comparatorType == "notInitiateWithin") {
      comparatorAfterG <- !initiate & afterG
      expect_true(any(comparatorAfterG), info = label)
      if (configuration$postGraceModel == "binary") {
        expect_equal(unique(numerator[comparatorAfterG]), 1 - postGrace$numerator$propensityScore, info = label)
        pi <- postGrace$predictions$propensityScore[match(u$rowId[comparatorAfterG], postGrace$predictions$rowId)]
        expect_equal(u$weight[comparatorAfterG], 1 / (1 - pi), tolerance = 1e-12, info = label)
      } else {
        expect_null(postGrace)
        expect_length(unique(round(numerator[comparatorAfterG], 12)), 1)
      }
    } else {
      expect_null(postGrace)
    }
    if (!deviationFirst) {
      hand <- handMsmWeights(analysis$models, u, riskStart = 0)
      expect_equal(numerator, hand$numerator, tolerance = 1e-10, info = label)
      expect_equal(u$weight, hand$weight, tolerance = 1e-10, info = label)
    }

    # With truncation, under either basis: the weight is the untruncated weight capped at the strategy's
    # cap, and the other weight the other untruncated weight capped at the cap divided (unstabilized) or
    # multiplied (stabilized) by the numerator.
    for (truncation in c("varying", "allClones")) {
      for (kind in kinds) {
        other <- setdiff(kinds, kind)
        result <- results[[truncation]][[kind]]
        info <- paste(label, truncation, kind)
        expect_equal(result$segments[, structure], u[, structure], info = info)
        expect_equal(result$segments$stabilizationNumerator, numerator, tolerance = 1e-12, info = info)
        cap <- result$caps$cap[match(result$segments$strategyId, result$caps$strategyId)]
        raw <- results$none[[kind]]$segments$weight
        expect_equal(result$segments$weight, pmin(raw, cap), tolerance = 1e-12, info = info)
        expect_true(any(result$segments$weight < raw * (1 - 1e-12)), info = info)
        otherCap <- if (kind == "stabilized") cap / numerator else cap * numerator
        expect_equal(result$msm[[other]], pmin(results$none[[other]]$segments$weight, otherCap), tolerance = 1e-12, info = info)
        expect_equal(result$msm[[kind]], result$segments$weight, info = info)
      }
    }
    # The two bases give different caps here, so both are exercised:
    for (kind in kinds) {
      expect_false(isTRUE(all.equal(results$varying[[kind]]$caps$cap, results$allClones[[kind]]$caps$cap)), info = paste(label, kind))
    }
  }
})

test_that("Bootstrap replicates fit the Cox MSM with the msmWeights weights", {
  skip_if_not_installed("survival")
  data <- simulateCloneCensorWeightData(
    n = 120, maxDays = 15, graceEnd = 3, initiationHazard = 0.1, outcomeHazard = 0.02, deathHazard = 0.005, seed = 2
  )
  on.exit(close(data))
  createStudyPopulationArgs <- CohortMethod::createCreateStudyPopulationArgs(riskWindowEnd = 15, endAnchor = "cohort start")
  createClonesArgs <- createCreateClonesArgs(strategies = list(
    createTreatmentStrategy(1, "Initiate within 3 days", "initiateWithin", graceEnd = 3),
    createTreatmentStrategy(2, "Never initiate", "neverInitiate")
  ))
  # Censoring models without covariates keep this test fast: the stabilized weights are then 1 and the
  # unstabilized ones 1 / stabilizationNumerator.
  fitCensoringModelsArgs <- createFitCensoringModelsArgs(useCovariates = FALSE, prior = Cyclops::createPrior("none"))
  # Stabilized and truncated weights (the defaults of computeCensoringWeights()):
  computeWeightsArgs <- createComputeWeightsArgs()
  # The weighted segments of the analysis of `cloneCensorWeightData` (with the censoring models
  # `censoringModels`, or fitted), and the censoring models:
  analyze <- function(cloneCensorWeightData, censoringModels = NULL) {
    population <- suppressMessages(createStudyPopulation(cloneCensorWeightData, outcomeId = 3, createStudyPopulationArgs = createStudyPopulationArgs))
    clones <- suppressMessages(createClones(population, cloneCensorWeightData, createClonesArgs))
    on.exit(close(clones))
    if (is.null(censoringModels)) {
      censoringModels <- suppressMessages(fitCensoringModels(clones, cloneCensorWeightData, fitCensoringModelsArgs))
    }
    weighted <- suppressMessages(computeCensoringWeights(clones, censoringModels, computeWeightsArgs))
    on.exit(close(weighted), add = TRUE)
    segments <- collectMsmSegments(weighted)
    segments$unstabilizedWeight <- segments$weight / segments$stabilizationNumerator
    return(list(segments = segments, censoringModels = censoringModels))
  }
  # The log hazard ratio of coxph on the segments with the msmWeights weights:
  referenceLogHr <- function(segments, msmWeights) {
    weights <- if (msmWeights == "unstabilized") segments$unstabilizedWeight else segments$weight
    reference <- survival::coxph(
      survival::Surv(startTime, stopTime, outcome) ~ x,
      data = segments, weights = weights, ties = "breslow", control = msmCoxControl()
    )
    return(unname(coef(reference)))
  }

  full <- analyze(data)
  set.seed(4)
  counts <- stats::rmultinom(1, size = 120, prob = rep(1, 120))
  personSeqIds <- sort(unique(dplyr::pull(data$cohorts, "personSeqId")))
  resampled <- .resampleCloneCensorWeightData(data, personSeqIds, counts[, 1])
  replicate <- analyze(resampled)
  close(resampled)
  results <- list()
  for (msmWeights in c("unstabilized", "stabilized")) {
    result <- suppressMessages(bootstrapCcwAnalysis(
      cloneCensorWeightData = data,
      outcomeId = 3,
      createStudyPopulationArgs = createStudyPopulationArgs,
      createClonesArgs = createClonesArgs,
      fitCensoringModelsArgs = fitCensoringModelsArgs,
      computeWeightsArgs = computeWeightsArgs,
      fitOutcomeModelArgs = createFitOutcomeModelArgs(modelType = "cox", horizons = c(8, 15), msmWeights = msmWeights),
      censoringModels = full$censoringModels,
      counts = counts
    ))
    expect_equal(result$metaData$args$msmWeights, msmWeights)
    # (6) The full-data and the replicate's hazard ratios are those of coxph with the chosen weights:
    expect_equal(result$estimates$logEstimate[result$estimates$estimand == "hr"], referenceLogHr(full$segments, msmWeights), tolerance = 1e-8, info = msmWeights)
    expect_equal(log(result$replicates[[1, "hr"]]), referenceLogHr(replicate$segments, msmWeights), tolerance = 1e-8, info = msmWeights)
    results[[msmWeights]] <- result
  }
  expect_gt(abs(log(results$unstabilized$replicates[[1, "hr"]] / results$stabilized$replicates[[1, "hr"]])), 0.01)
  # The Kaplan-Meier replicates are the same under both:
  km <- colnames(results$unstabilized$replicates) != "hr"
  expect_equal(results$unstabilized$replicates[, km], results$stabilized$replicates[, km], tolerance = 1e-12)
})

test_that("Weighted clones without stabilization numerators stop with a message to rerun computeCensoringWeights", {
  # Review of contract 11.1: weighted clones made before the numerator was stored gave a raw dplyr error.
  analysis <- msmAnalysis("notInitiateWithin")
  on.exit(closeMsmAnalysis(analysis))
  weighted <- weightMsm(analysis, msmComputeWeightsArgs(stabilized = TRUE))
  on.exit(close(weighted), add = TRUE)
  segments <- dplyr::collect(weighted$segments)
  segments$stabilizationNumerator <- NULL
  weighted$segments <- segments
  for (modelType in c("cox", "km")) {
    expect_error(
      fitOutcomeModel(weighted, createFitOutcomeModelArgs(modelType = modelType, horizons = 10)),
      "earlier version of computeCensoringWeights\\(\\). Rerun computeCensoringWeights\\(\\)"
    )
  }
})

test_that("The cap bounds the weights of computeCensoringWeights, and the MSM weights only when of the same kind", {
  # Contract 11.1 as documented at review: the other weight is formed from the truncated weight, so the
  # cap computed on the stabilized scale does not bound the unstabilized MSM weights of neverInitiate.
  analysis <- msmAnalysis("neverInitiate")
  on.exit(closeMsmAnalysis(analysis))
  truncation <- createWeightTruncation("percentile", 0.9)
  for (stabilized in c(TRUE, FALSE)) {
    weighted <- weightMsm(analysis, msmComputeWeightsArgs(stabilized = stabilized, truncation = truncation))
    segments <- collectMsmSegments(weighted)
    caps <- attr(weighted, "metaData")$weights$caps
    close(weighted)
    cap <- caps$cap[segments$strategyId]
    expect_true(all(segments$weight <= cap + 1e-12))
    unstabilized <- if (stabilized) segments$weight / segments$stabilizationNumerator else segments$weight
    if (stabilized) {
      expect_gt(max(unstabilized[segments$strategyId == 2]), caps$cap[2])
    } else {
      expect_true(all(unstabilized <= cap + 1e-12))
    }
  }
  # With truncateMsmWeights = TRUE (contract section 12.2) the unstabilized MSM weights have a cap of their own:
  weighted <- weightMsm(analysis, msmComputeWeightsArgs(stabilized = TRUE, truncation = truncation, truncateMsmWeights = TRUE))
  outcomeData <- .withMsmWeights(.getOutcomeData(weighted), "unstabilized")
  caps <- attr(weighted, "metaData")$weights$caps
  close(weighted)
  expect_true(all(outcomeData$segments$weight <= caps$otherScaleCap[outcomeData$segments$strategyId] + 1e-12))
})

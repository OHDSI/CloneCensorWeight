library(testthat)

# Round trip through JSON (toJson() then json =) and through an untyped list.
expectRoundTrip <- function(settings, generator) {
  json <- settings$toJson()
  fromJson <- generator$new(json = json)
  expect_equal(fromJson, settings)
  expect_equal(fromJson$toJson(), json)
  untypedList <- jsonlite::fromJSON(json, simplifyDataFrame = FALSE)
  fromList <- generator$new(untypedList = untypedList)
  expect_equal(fromList, settings)
  expect_equal(fromList$toJson(), json)
}

createTestStrategies <- function() {
  return(list(
    createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
    createTreatmentStrategy(2, "Do not initiate within 7 days", "notInitiateWithin", graceEnd = 7)
  ))
}

createTestAnalysis <- function(analysisId = 1, computeCovariateBalanceArgs = NULL) {
  return(createCcwAnalysis(
    analysisId = analysisId,
    description = "Test analysis",
    getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(),
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
      riskWindowEnd = 365,
      endAnchor = "cohort start"
    ),
    createClonesArgs = createCreateClonesArgs(strategies = createTestStrategies()),
    fitCensoringModelsArgs = createFitCensoringModelsArgs(),
    computeWeightsArgs = createComputeWeightsArgs(),
    computeCovariateBalanceArgs = computeCovariateBalanceArgs,
    fitOutcomeModelArgs = createFitOutcomeModelArgs()
  ))
}

createTestTargetTreatmentOutcomes <- function(targetId = 1) {
  return(createTargetTreatmentOutcomes(
    targetId = targetId,
    treatmentId = 21,
    comparatorTreatmentId = 22,
    outcomes = list(
      CohortMethod::createOutcome(outcomeId = 3),
      CohortMethod::createOutcome(outcomeId = 4, outcomeOfInterest = FALSE, trueEffectSize = 1),
      CohortMethod::createOutcome(outcomeId = 5, riskWindowEnd = 30, endAnchor = "cohort start")
    ),
    excludedCovariateConceptIds = c(100, 200)
  ))
}

test_that("TreatmentStrategy round trip and validation", {
  expectRoundTrip(createTestStrategies()[[1]], TreatmentStrategy)
  expectRoundTrip(createTestStrategies()[[2]], TreatmentStrategy)
  expectRoundTrip(createTreatmentStrategy(3, "Never", "neverInitiate"), TreatmentStrategy)
  expectRoundTrip(
    createTreatmentStrategy(4, "Head to head", "initiateWithin",
      treatment = "comparatorTreatment", graceStart = 2, graceEnd = 14,
      prohibitedTreatments = "treatment"
    ),
    TreatmentStrategy
  )

  expect_error(createTreatmentStrategy(1, "x", "initiateWithin"), "graceEnd is required")
  expect_error(createTreatmentStrategy(1, "x", "neverInitiate", graceEnd = 7), "must be NULL")
  expect_error(createTreatmentStrategy(1, "x", "initiateWithin", graceStart = 10, graceEnd = 7), "greater than or equal")
  expect_error(createTreatmentStrategy(1, "x", "initiateWithin", graceEnd = -1), "graceEnd")
  expect_error(createTreatmentStrategy(1, "x", "startSometime", graceEnd = 7), "type")
  expect_error(createTreatmentStrategy(1, "x", "initiateWithin", treatment = "drugA", graceEnd = 7), "treatment")
  expect_error(createTreatmentStrategy(1.5, "x", "initiateWithin", graceEnd = 7), "strategyId")
  expect_error(
    createTreatmentStrategy(1, "x", "initiateWithin", graceEnd = 7, prohibitedTreatments = "treatment"),
    "own treatment"
  )
  expect_error(
    createTreatmentStrategy(1, "x", "initiateWithin", graceEnd = 7, prohibitedTreatments = "drugB"),
    "prohibitedTreatments"
  )
})

test_that("GetDbCloneCensorWeightDataArgs round trip and validation", {
  expectRoundTrip(createGetDbCloneCensorWeightDataArgs(), GetDbCloneCensorWeightDataArgs)
  expectRoundTrip(
    createGetDbCloneCensorWeightDataArgs(
      studyStartDate = "20100101",
      studyEndDate = "20201231",
      covariateSettings = FeatureExtraction::createDefaultCovariateSettings(
        excludedCovariateConceptIds = c(1, 2),
        addDescendantsToExclude = TRUE
      )
    ),
    GetDbCloneCensorWeightDataArgs
  )
  complexCovariateSettings <- list(
    FeatureExtraction::createCovariateSettings(useDemographicsAgeGroup = TRUE, useConditionGroupEraLongTerm = TRUE),
    FeatureExtraction::createCohortBasedCovariateSettings(
      analysisId = 999,
      covariateCohorts = data.frame(cohortId = 3, cohortName = "Feature cohort")
    )
  )
  expectRoundTrip(
    createGetDbCloneCensorWeightDataArgs(covariateSettings = complexCovariateSettings),
    GetDbCloneCensorWeightDataArgs
  )

  expect_error(createGetDbCloneCensorWeightDataArgs(studyStartDate = "20101301"), "not valid")
  expect_error(createGetDbCloneCensorWeightDataArgs(washoutPeriod = -1), "washoutPeriod")
  expect_error(createGetDbCloneCensorWeightDataArgs(firstOnly = "yes"), "firstOnly")
  expect_error(createGetDbCloneCensorWeightDataArgs(covariateSettings = "default"), "covariateSettings")
})

test_that("CreateClonesArgs round trip and validation", {
  expectRoundTrip(createCreateClonesArgs(strategies = createTestStrategies()), CreateClonesArgs)
  expectRoundTrip(
    createCreateClonesArgs(
      strategies = list(
        createTreatmentStrategy(1, "Initiate within 30 days", "initiateWithin", graceEnd = 30),
        createTreatmentStrategy(2, "Never initiate", "neverInitiate")
      ),
      sameDayOrdering = "deviationFirst",
      competingEvents = "none"
    ),
    CreateClonesArgs
  )

  strategies <- createTestStrategies()
  expect_error(createCreateClonesArgs(strategies = strategies[1]), "length 2")
  expect_error(createCreateClonesArgs(strategies = c(strategies, strategies[1])), "length 2")
  expect_error(createCreateClonesArgs(strategies = list(strategies[[1]], strategies[[1]])), "different strategyId")
  expect_error(createCreateClonesArgs(strategies = list(strategies[[1]], list(strategyId = 2))), "strategies.*R6")
  expect_error(createCreateClonesArgs(strategies = strategies, sameDayOrdering = "random"), "sameDayOrdering")
  expect_error(createCreateClonesArgs(strategies = strategies, competingEvents = "censor"), "competingEvents")
})

test_that("TimeBasis round trip and validation", {
  expectRoundTrip(createTimeBasis(), TimeBasis)
  expectRoundTrip(createTimeBasis(type = "naturalSpline", intervalLength = 7, df = 3), TimeBasis)

  expect_error(createTimeBasis(type = "cubic"), "type")
  expect_error(createTimeBasis(intervalLength = 0), "intervalLength")
  expect_error(createTimeBasis(df = 0), "df")
  expectRoundTrip(createTimeBasis(type = "polynomial", intervalLength = 7, degree = 3, raw = FALSE), TimeBasis)
  expect_error(createTimeBasis(type = "polynomial", degree = 0), "degree")
  expect_error(createTimeBasis(raw = NA), "raw")
})

test_that("FitCensoringModelsArgs round trip and validation", {
  expectRoundTrip(createFitCensoringModelsArgs(), FitCensoringModelsArgs)
  # Contract section 9, amendment 4: the default hazard model is pooled logistic with daily indicators.
  defaults <- createFitCensoringModelsArgs()
  expect_equal(defaults$hazardModelType, "logistic")
  expect_equal(defaults$timeBasis$type, "indicator")
  expect_equal(defaults$timeBasis$intervalLength, 1)
  expectRoundTrip(
    createFitCensoringModelsArgs(
      hazardModelType = "logistic",
      timeBasis = createTimeBasis(type = "linear", intervalLength = 7),
      useCovariates = FALSE,
      excludeCovariateIds = c(1002, 1003),
      includeCovariateIds = c(1001),
      prior = Cyclops::createPrior("normal", variance = 2),
      control = Cyclops::createControl(seed = 123, threads = 2)
    ),
    FitCensoringModelsArgs
  )

  expect_error(createFitCensoringModelsArgs(hazardModelType = "weibull"), "hazardModelType")
  expect_error(createFitCensoringModelsArgs(censoringMechanisms = c("artificial", "observationEnd")), "censoringMechanisms")
  expect_error(createFitCensoringModelsArgs(timeBasis = list(type = "indicator")), "timeBasis")
  expect_error(createFitCensoringModelsArgs(excludeCovariateIds = 1.5), "integer")
  expect_error(createFitCensoringModelsArgs(prior = list(priorType = "laplace")), "cyclopsPrior")
  expect_error(createFitCensoringModelsArgs(control = list()), "cyclopsControl")
})

test_that("FitCensoringModelsArgs postGraceModel and covariateTimeInteraction (contract sections 11.2-11.3 and 12.1)", {
  defaults <- createFitCensoringModelsArgs()
  expect_equal(defaults$postGraceModel, "hazard")
  expect_equal(defaults$covariateTimeInteraction, "day0")
  # The Cox model has no interactions, so its default is "none":
  expect_equal(createFitCensoringModelsArgs(hazardModelType = "cox")$covariateTimeInteraction, "none")
  # With longer intervals "day0" is the first interval (contract section 12.1):
  expect_equal(createFitCensoringModelsArgs(timeBasis = createTimeBasis(intervalLength = 7))$covariateTimeInteraction, "day0")
  for (covariateTimeInteraction in c("none", "day0", "interval")) {
    for (postGraceModel in c("binary", "hazard")) {
      expectRoundTrip(
        createFitCensoringModelsArgs(covariateTimeInteraction = covariateTimeInteraction, postGraceModel = postGraceModel),
        FitCensoringModelsArgs
      )
    }
  }
  expectRoundTrip(createFitCensoringModelsArgs(hazardModelType = "cox", postGraceModel = "hazard"), FitCensoringModelsArgs)
  expect_error(createFitCensoringModelsArgs(postGraceModel = "logistic"), "postGraceModel")
  expect_error(createFitCensoringModelsArgs(covariateTimeInteraction = "week1"), "covariateTimeInteraction")
  for (covariateTimeInteraction in c("day0", "interval")) {
    expect_error(
      createFitCensoringModelsArgs(hazardModelType = "cox", covariateTimeInteraction = covariateTimeInteraction),
      "covariateTimeInteraction must be \"none\""
    )
  }
  # Settings saved before these fields existed load with the defaults ("none" for the Cox model):
  for (hazardModelType in c("logistic", "cox")) {
    untypedList <- jsonlite::fromJSON(createFitCensoringModelsArgs(hazardModelType = hazardModelType, postGraceModel = "binary")$toJson(), simplifyDataFrame = FALSE)
    untypedList$postGraceModel <- NULL
    untypedList$covariateTimeInteraction <- NULL
    expect_message(old <- FitCensoringModelsArgs$new(untypedList = untypedList), "saved without covariateTimeInteraction, postGraceModel")
    expect_equal(old$postGraceModel, "hazard")
    expect_equal(old$covariateTimeInteraction, if (hazardModelType == "cox") "none" else "day0")
    expect_message(old <- FitCensoringModelsArgs$new(json = jsonlite::toJSON(untypedList, auto_unbox = TRUE, null = "null", digits = NA)), sprintf("covariateTimeInteraction = %s", if (hazardModelType == "cox") "none" else "day0"))
    expect_equal(old$covariateTimeInteraction, if (hazardModelType == "cox") "none" else "day0")
  }
})

test_that("WeightTruncation basis (contract section 11.5)", {
  expect_equal(createWeightTruncation()$basis, "varying")
  expectRoundTrip(createWeightTruncation(basis = "allClones"), WeightTruncation)
  expectRoundTrip(createWeightTruncation(type = "percentile", upper = 0.95, basis = "varying"), WeightTruncation)
  expectRoundTrip(createComputeWeightsArgs(truncation = createWeightTruncation(basis = "allClones")), ComputeWeightsArgs)
  expect_error(createWeightTruncation(basis = "all"), "basis")
  untypedList <- jsonlite::fromJSON(createWeightTruncation(basis = "allClones")$toJson(), simplifyDataFrame = FALSE)
  untypedList$basis <- NULL
  expect_message(expect_equal(WeightTruncation$new(untypedList = untypedList)$basis, "varying"), "basis = varying")
  # otherScaleUpper (review of section 12): the cap of a fixed rule on the other scale
  expectRoundTrip(createWeightTruncation(type = "fixed", upper = 5, otherScaleUpper = 50), WeightTruncation)
  expect_null(createWeightTruncation(type = "fixed", upper = 5)$otherScaleUpper)
  expect_error(createWeightTruncation(type = "fixed", upper = 5, otherScaleUpper = 0.5), "otherScaleUpper")
})

test_that("WeightTruncation and ComputeWeightsArgs round trip and validation", {
  expectRoundTrip(createWeightTruncation(), WeightTruncation)
  expectRoundTrip(createWeightTruncation(type = "fixed", upper = 20), WeightTruncation)
  expectRoundTrip(createWeightTruncation(type = "none", upper = NULL), WeightTruncation)
  expectRoundTrip(createComputeWeightsArgs(), ComputeWeightsArgs)
  expectRoundTrip(
    createComputeWeightsArgs(stabilized = FALSE, truncation = createWeightTruncation(type = "fixed", upper = 10)),
    ComputeWeightsArgs
  )

  expect_error(createWeightTruncation(type = "percentile", upper = 1), "in \\(0, 1\\)")
  expect_error(createWeightTruncation(type = "percentile", upper = 0), "in \\(0, 1\\)")
  expect_error(createWeightTruncation(type = "fixed", upper = 0.5), "greater than 1")
  expect_error(createWeightTruncation(type = "trim"), "type")
  expect_error(createComputeWeightsArgs(stabilized = NA), "stabilized")
  expect_error(createComputeWeightsArgs(truncation = 0.99), "truncation")
  # truncateMsmWeights (contract section 12.2): TRUE by default, FALSE available, and settings saved
  # before it existed load with TRUE:
  expect_true(createComputeWeightsArgs()$truncateMsmWeights)
  expectRoundTrip(createComputeWeightsArgs(truncateMsmWeights = FALSE), ComputeWeightsArgs)
  expect_error(createComputeWeightsArgs(truncateMsmWeights = NA), "truncateMsmWeights")
  untypedList <- jsonlite::fromJSON(createComputeWeightsArgs(truncateMsmWeights = FALSE)$toJson(), simplifyDataFrame = FALSE)
  untypedList$truncateMsmWeights <- NULL
  expect_message(expect_true(ComputeWeightsArgs$new(untypedList = untypedList)$truncateMsmWeights), "ComputeWeightsArgs saved without truncateMsmWeights: filled with the current default truncateMsmWeights = TRUE")
  # Complete settings load without a message:
  expect_no_message(ComputeWeightsArgs$new(json = createComputeWeightsArgs()$toJson()))
  expect_no_message(FitCensoringModelsArgs$new(json = createFitCensoringModelsArgs()$toJson()))
})

test_that("BootstrapArgs and FitOutcomeModelArgs round trip and validation", {
  expectRoundTrip(createBootstrapArgs(), BootstrapArgs)
  expectRoundTrip(createBootstrapArgs(replicates = 50, seed = 42, ci = "normal", fixHyperparameters = FALSE), BootstrapArgs)
  expectRoundTrip(createFitOutcomeModelArgs(), FitOutcomeModelArgs)
  # Contract section 9, amendments 6 and 7: daily natural-spline basis and Newton steps until the relative
  # change of the objective is below 1e-10.
  defaults <- createFitOutcomeModelArgs()
  expect_equal(defaults$timeBasis$toList()[c("type", "intervalLength", "df")], list(type = "naturalSpline", intervalLength = 1, df = 4))
  expect_equal(defaults$control[c("algorithm", "convergenceType", "tolerance")], list(algorithm = "newton", convergenceType = "gradient", tolerance = 1e-10))
  expect_equal(FitOutcomeModelArgs$new(json = defaults$toJson())$control$algorithm, "newton")
  # Contract section 11.1: the MSMs use unstabilized weights by default.
  expect_equal(defaults$msmWeights, "unstabilized")
  expectRoundTrip(
    createFitOutcomeModelArgs(
      modelType = "logistic",
      horizons = c(30, 180, 365),
      timeBasis = createTimeBasis(type = "naturalSpline", intervalLength = 7, df = 5),
      armByTimeInteraction = TRUE,
      msmWeights = "stabilized",
      varianceEstimator = "bootstrap",
      bootstrapArgs = createBootstrapArgs(replicates = 100),
      prior = Cyclops::createPrior("laplace", variance = 0.5),
      control = Cyclops::createControl(tolerance = 1e-08)
    ),
    FitOutcomeModelArgs
  )

  expect_error(createBootstrapArgs(replicates = 0), "replicates")
  expect_error(createBootstrapArgs(ci = "bca"), "ci")
  expect_error(createFitOutcomeModelArgs(modelType = "poisson"), "modelType")
  expect_error(createFitOutcomeModelArgs(horizons = c(30.5)), "horizons")
  expect_error(createFitOutcomeModelArgs(horizons = c(-1)), "horizons")
  expect_error(createFitOutcomeModelArgs(varianceEstimator = "jackknife"), "varianceEstimator")
  # Settings saved before msmWeights existed load with the default:
  untypedList <- jsonlite::fromJSON(createFitOutcomeModelArgs(msmWeights = "stabilized")$toJson(), simplifyDataFrame = FALSE)
  untypedList$msmWeights <- NULL
  expect_message(expect_equal(FitOutcomeModelArgs$new(untypedList = untypedList)$msmWeights, "unstabilized"), "msmWeights = unstabilized")
  expect_message(expect_equal(FitOutcomeModelArgs$new(json = jsonlite::toJSON(untypedList, auto_unbox = TRUE, null = "null", digits = NA))$msmWeights, "unstabilized"), "msmWeights")
  expect_error(createFitOutcomeModelArgs(msmWeights = "truncated"), "msmWeights")
  expect_error(createFitOutcomeModelArgs(msmWeights = TRUE), "msmWeights")
  expect_error(createFitOutcomeModelArgs(useCovariates = TRUE), "Not yet supported")
  expect_error(createFitOutcomeModelArgs(stratified = TRUE), "Not yet supported")
})

test_that("CcwAnalysis round trip with nested CohortMethod, Cyclops and FeatureExtraction settings", {
  expectRoundTrip(createTestAnalysis(), CcwAnalysis)
  expectRoundTrip(
    createTestAnalysis(computeCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs()),
    CcwAnalysis
  )
  expectRoundTrip(
    createTestAnalysis(computeCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(
      covariateFilter = FeatureExtraction::getDefaultTable1Specifications()
    )),
    CcwAnalysis
  )
  analysis <- CcwAnalysis$new(json = createTestAnalysis()$toJson())
  expect_true(is(analysis$createStudyPopulationArgs, "CreateStudyPopulationArgs"))
  expect_equal(analysis$createStudyPopulationArgs$riskWindowEnd, 365)
  expect_true(is(analysis$createClonesArgs$strategies[[2]], "TreatmentStrategy"))
  expect_s3_class(analysis$fitCensoringModelsArgs$prior, "cyclopsPrior")
  expect_s3_class(analysis$fitOutcomeModelArgs$control, "cyclopsControl")
  expect_true(is.list(analysis$getDbCloneCensorWeightDataArgs$covariateSettings))
  expect_s3_class(analysis$getDbCloneCensorWeightDataArgs$covariateSettings, "covariateSettings")

  expect_error(createTestAnalysis(analysisId = "a"), "analysisId")
  expect_error(
    createCcwAnalysis(
      analysisId = 1,
      getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(),
      createStudyPopulationArgs = list(),
      createClonesArgs = createCreateClonesArgs(strategies = createTestStrategies()),
      fitCensoringModelsArgs = createFitCensoringModelsArgs(),
      computeWeightsArgs = createComputeWeightsArgs(),
      fitOutcomeModelArgs = createFitOutcomeModelArgs()
    ),
    "createStudyPopulationArgs.*R6"
  )
  expect_error(createTestAnalysis(computeCovariateBalanceArgs = createComputeWeightsArgs()), "ComputeCovariateBalanceArgs")
  # CohortMethod's class of the same name is rejected:
  expect_error(
    createCcwAnalysis(
      analysisId = 1,
      getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(),
      createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(),
      createClonesArgs = createCreateClonesArgs(strategies = createTestStrategies()),
      fitCensoringModelsArgs = createFitCensoringModelsArgs(),
      computeWeightsArgs = createComputeWeightsArgs(),
      fitOutcomeModelArgs = CohortMethod::createFitOutcomeModelArgs(modelType = "cox")
    ),
    "fitOutcomeModelArgs.*horizons"
  )
})

test_that("Settings keep full numeric precision through JSON", {
  truncation <- createWeightTruncation(type = "percentile", upper = 0.99999)
  expectRoundTrip(truncation, WeightTruncation)
  expect_identical(WeightTruncation$new(json = truncation$toJson())$upper, 0.99999)

  args <- createFitCensoringModelsArgs(
    prior = Cyclops::createPrior("laplace", variance = 0.123456789, exclude = 0),
    control = Cyclops::createControl(startingVariance = 0.0012345678, tolerance = 1e-9)
  )
  expectRoundTrip(args, FitCensoringModelsArgs)
  reloaded <- FitCensoringModelsArgs$new(json = args$toJson())
  expect_identical(reloaded$prior$variance, 0.123456789)
  expect_identical(reloaded$control$startingVariance, 0.0012345678)
  expect_identical(reloaded$control$tolerance, 1e-9)

  thresholds <- createCcwDiagnosticThresholds(easeThreshold = 0.123456)
  expect_identical(CcwDiagnosticThresholds$new(json = thresholds$toJson())$easeThreshold, 0.123456)

  # Through the list save and load functions:
  tto <- createTargetTreatmentOutcomes(1, 21, outcomes = list(CohortMethod::createOutcome(3, trueEffectSize = 1.23456789)))
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file))
  saveTargetTreatmentOutcomesList(list(tto), file)
  expect_identical(loadTargetTreatmentOutcomesList(file)[[1]]$outcomes[[1]]$trueEffectSize, 1.23456789)
})

test_that("TargetTreatmentOutcomes round trip and validation", {
  expectRoundTrip(createTestTargetTreatmentOutcomes(), TargetTreatmentOutcomes)
  expectRoundTrip(
    createTargetTreatmentOutcomes(targetId = 1, treatmentId = 21, outcomes = list(CohortMethod::createOutcome(outcomeId = 3))),
    TargetTreatmentOutcomes
  )
  tto <- TargetTreatmentOutcomes$new(json = createTestTargetTreatmentOutcomes()$toJson())
  expect_true(is(tto$outcomes[[1]], "Outcome"))
  expect_true(is.na(tto$outcomes[[1]]$trueEffectSize))
  expect_equal(tto$outcomes[[2]]$trueEffectSize, 1)
  expect_false(tto$outcomes[[2]]$outcomeOfInterest)

  expect_error(
    createTargetTreatmentOutcomes(targetId = 1, treatmentId = 21, outcomes = list()),
    "outcomes"
  )
  expect_error(
    createTargetTreatmentOutcomes(targetId = 1, treatmentId = 21, outcomes = list(list(outcomeId = 3))),
    "outcomes.*R6"
  )
  expect_error(
    createTargetTreatmentOutcomes(
      targetId = 1, treatmentId = 21,
      outcomes = list(CohortMethod::createOutcome(outcomeId = 3), CohortMethod::createOutcome(outcomeId = 3))
    ),
    "duplicate outcome IDs"
  )
  expect_error(
    createTargetTreatmentOutcomes(targetId = 1.5, treatmentId = 21, outcomes = list(CohortMethod::createOutcome(outcomeId = 3))),
    "integer"
  )
})

test_that("CcwDiagnosticThresholds round trip and validation", {
  expectRoundTrip(createCcwDiagnosticThresholds(), CcwDiagnosticThresholds)
  expectRoundTrip(createCcwDiagnosticThresholds(mdrrThreshold = 5, maxWeightThreshold = 20), CcwDiagnosticThresholds)

  expect_error(createCcwDiagnosticThresholds(sdmThreshold = -0.1), "sdmThreshold")
  expect_error(createCcwDiagnosticThresholds(essFractionThreshold = 2), "essFractionThreshold")
})

test_that("CcwAnalysesSpecifications round trip, validation and conversion", {
  specifications <- createCcwAnalysesSpecifications(
    ccwAnalysisList = list(
      createTestAnalysis(1),
      createTestAnalysis(2, CohortMethod::createComputeCovariateBalanceArgs())
    ),
    targetTreatmentOutcomesList = list(createTestTargetTreatmentOutcomes(1), createTestTargetTreatmentOutcomes(2)),
    analysesToExclude = data.frame(targetId = 2, analysisId = 1),
    refitCensoringModelsForEveryOutcome = TRUE,
    ccwDiagnosticThresholds = createCcwDiagnosticThresholds(mdrrThreshold = 4)
  )
  expectRoundTrip(specifications, CcwAnalysesSpecifications)
  expectRoundTrip(
    createCcwAnalysesSpecifications(list(createTestAnalysis()), list(createTestTargetTreatmentOutcomes())),
    CcwAnalysesSpecifications
  )

  untypedList <- jsonlite::fromJSON(specifications$toJson(), simplifyDataFrame = FALSE)
  converted <- convertUntypedListToCcwAnalysesSpecifications(untypedList)
  expect_equal(converted, specifications)
  expect_true(is(converted, "CcwAnalysesSpecifications"))

  expect_error(
    createCcwAnalysesSpecifications(list(createTestAnalysis(1), createTestAnalysis(1)), list(createTestTargetTreatmentOutcomes())),
    "Duplicate analysis IDs"
  )
  expect_error(
    createCcwAnalysesSpecifications(list(createTestAnalysis()), list(createTestTargetTreatmentOutcomes()), analysesToExclude = data.frame(x = 1)),
    "at least one of these columns"
  )
  expect_error(createCcwAnalysesSpecifications(list(), list(createTestTargetTreatmentOutcomes())), "ccwAnalysisList")
  expect_error(createCcwAnalysesSpecifications(list(createTestAnalysis()), list(createTestAnalysis())), "TargetTreatmentOutcomes")
})

test_that("Saving and loading analysis and target-treatment-outcomes lists", {
  analysisList <- list(createTestAnalysis(1), createTestAnalysis(2, CohortMethod::createComputeCovariateBalanceArgs()))
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file))
  saveCcwAnalysisList(analysisList, file)
  expect_equal(loadCcwAnalysisList(file), analysisList)

  ttoList <- list(createTestTargetTreatmentOutcomes(1), createTestTargetTreatmentOutcomes(2))
  saveTargetTreatmentOutcomesList(ttoList, file)
  expect_equal(loadTargetTreatmentOutcomesList(file), ttoList)

  expect_error(saveCcwAnalysisList(ttoList, file), "CcwAnalysis")
  expect_error(saveTargetTreatmentOutcomesList(analysisList, file), "TargetTreatmentOutcomes")
})

test_that("Constructors without arguments fail", {
  expect_error(TimeBasis$new(), "Must provide either")
})

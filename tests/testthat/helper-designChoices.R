# Helpers of the design-choice tests (contract section 12): test-integration-designChoices*.R.

designStrategies <- function(second = "notInitiateWithin", graceEnd = 7) {
  list(
    createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = graceEnd),
    createTreatmentStrategy(2, "Comparator", second, graceEnd = if (second == "neverInitiate") NULL else graceEnd)
  )
}

designClones <- function(data, second = "notInitiateWithin", sameDayOrdering = "eventsFirst") {
  population <- suppressMessages(createStudyPopulation(
    cloneCensorWeightData = data,
    outcomeId = 3,
    createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(riskWindowStart = 0, riskWindowEnd = 60, endAnchor = "cohort start")
  ))
  return(suppressMessages(createClones(population, data, createCreateClonesArgs(strategies = designStrategies(second), sameDayOrdering = sameDayOrdering))))
}

# A fixed prior: unpenalized fits of the covariate-by-day-0 interactions can separate in small samples.
designArgs <- function(..., prior = Cyclops::createPrior("laplace", exclude = c(0), variance = 1)) {
  createFitCensoringModelsArgs(prior = prior, ...)
}

# The grace-period DGP of the statistical review in which starting the treatment cuts the outcome hazard
# by 80% from the next day ("bStrong"), and its null version.
selectionData <- function(treatmentEffect, n = 5000, seed = 1) {
  return(simulateCloneCensorWeightData(
    n = n, seed = seed, nCovariates = 1, dgp = "b", graceEnd = 7, maxDays = 60, binaryPrevalence = 0.4,
    initiationHazard = 0.08, outcomeHazard = 0.03, deathHazard = 0.002, observationEndHazard = 0.003,
    initiationCoefficients = 1.5, outcomeCoefficients = 1.5, deathCoefficients = 0.5,
    treatmentEffect = treatmentEffect, effectRampDays = 1
  ))
}

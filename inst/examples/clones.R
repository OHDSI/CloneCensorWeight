# Simulated data, its study population, and the clones of two strategies:
data <- simulateCloneCensorWeightData(n = 200, maxDays = 60, graceEnd = 7, seed = 1)
population <- createStudyPopulation(
  cloneCensorWeightData = data,
  outcomeId = 3,
  createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
    riskWindowEnd = 60,
    endAnchor = "cohort start"
  )
)
clones <- createClones(
  population = population,
  cloneCensorWeightData = data,
  createClonesArgs = createCreateClonesArgs(
    strategies = list(
      createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
      createTreatmentStrategy(2, "Do not initiate within 7 days", "notInitiateWithin", graceEnd = 7)
    )
  )
)

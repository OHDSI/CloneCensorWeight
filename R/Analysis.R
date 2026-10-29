#' Create a clone-censor-weight analysis specification
#' @param analysisId Analysis identifier.
#' @param arms Treatment strategy definitions.
#' @param gracePeriod Grace period in days.
#' @param washoutPeriod Washout period in days.
#' @param intervalLength Person-time interval length in days.
#' @param modelType Outcome model type.
#' @return A `CloneCensorWeightAnalysis` object.
#' @import R6
#' @export
CloneCensorWeightAnalysis <- R6::R6Class(
  "CloneCensorWeightAnalysis",
  public = list(
    analysisId = NULL, arms = NULL, gracePeriod = NULL, washoutPeriod = NULL,
    intervalLength = NULL, modelType = NULL,
    initialize = function(analysisId, arms, gracePeriod, washoutPeriod = 0,
                          intervalLength = 30, modelType = "cox") {
      checkmate::assertCharacter(analysisId, len = 1, any.missing = FALSE)
      checkmate::assertList(arms, min.len = 2)
      checkmate::assertNumeric(gracePeriod, len = 1, lower = 0)
      checkmate::assertNumeric(washoutPeriod, len = 1, lower = 0)
      checkmate::assertNumeric(intervalLength, len = 1, lower = 1)
      checkmate::assertChoice(modelType, c("cox", "logistic"))
      self$analysisId <- analysisId
      self$arms <- arms
      self$gracePeriod <- gracePeriod
      self$washoutPeriod <- washoutPeriod
      self$intervalLength <- intervalLength
      self$modelType <- modelType
    },
    toList = function() as.list(self)
  )
)

#' Create a clone-censor-weight analysis specification
#' @param analysisId Analysis identifier.
#' @param arms Treatment strategy definitions.
#' @param gracePeriod Grace period in days.
#' @param washoutPeriod Washout period in days.
#' @param intervalLength Person-time interval length in days.
#' @param modelType Outcome model type.
#' @return A `CloneCensorWeightAnalysis` object.
#' @export
createCloneCensorWeightAnalysis <- function(analysisId, arms, gracePeriod,
                                            washoutPeriod = 0, intervalLength = 30,
                                            modelType = "cox") {
  CloneCensorWeightAnalysis$new(analysisId, arms, gracePeriod, washoutPeriod,
                                intervalLength, modelType)
}

#' Run clone-censor-weight analyses
#' @param cloneCensorWeightData A `CloneCensorWeightData` object.
#' @param analyses A list of `CloneCensorWeightAnalysis` objects.
#' @param outputFolder Optional output folder.
#' @return A named list of analysis results.
#' @export
runCloneCensorWeightAnalyses <- function(cloneCensorWeightData, analyses, outputFolder = NULL) {
  checkmate::assertClass(cloneCensorWeightData, "CloneCensorWeightData")
  checkmate::assertList(analyses, min.len = 1)
  if (!is.null(outputFolder)) dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
  results <- lapply(analyses, function(analysis) {
    checkmate::assertR6(analysis, "CloneCensorWeightAnalysis")
    cloned <- createClonesAndCensor(cloneCensorWeightData, analysis$arms,
                                    analysis$gracePeriod, analysis$washoutPeriod)
    personTime <- discretizePersonTime(cloned, analysis$intervalLength)
    weighted <- createCensoringWeights(personTime, cloned)
    outcomeModel <- fitOutcomeModel(weighted, analysis$modelType)
    result <- list(analysis = analysis, data = cloned, personTime = weighted,
                   outcomeModel = outcomeModel,
                   weightDistribution = computeWeightDistribution(weighted))
    if (!is.null(outputFolder)) {
      folder <- file.path(outputFolder, analysis$analysisId)
      dir.create(folder, recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(outcomeModel$summary, file.path(folder, "results.csv"), row.names = FALSE)
      utils::write.csv(result$weightDistribution,
                       file.path(folder, "weightDistribution.csv"), row.names = FALSE)
      saveCloneCensorWeightData(cloned, file.path(folder, "cloneCensorWeightData"))
    }
    result
  })
  names(results) <- vapply(analyses, function(analysis) analysis$analysisId, character(1))
  results
}
# -------------------------------------------------------------------------
# Internal helpers for covariate settings serialization
# -------------------------------------------------------------------------

# Preserves the S3/S4 class of a covariateSettings object as a list member
# so that jsonlite can serialize it. Mirrors CohortMethod's approach.
.convertAttrToMember <- function(object) {
  result <- lapply(object, function(x) x)
  result[["attr_class"]] <- class(object)
  result
}

# Restores the class attribute after deserializing from JSON.
.convertMemberToAttr <- function(object) {
  if (!is.null(object[["attr_class"]])) {
    cls <- object[["attr_class"]]
    object[["attr_class"]] <- NULL
    class(object) <- cls
  }
  object
}

# Reconstruct a plain list into a typed settings object. Used by load functions.
.reconstructSettings <- function(settingsList, className) {
  obj <- get(className)$new()
  nms <- intersect(names(settingsList), ls(obj))
  for (nm in nms) {
    if (!is.function(obj[[nm]])) obj[[nm]] <- settingsList[[nm]]
  }
  obj
}


# -------------------------------------------------------------------------
# AbstractSerializableSettings
# -------------------------------------------------------------------------

#' Abstract base class for serializable settings
#'
#' Provides `toList()`, `toJson()`, and `validate()` to all settings classes.
#' Users should not instantiate this class directly; use the concrete subclasses.
#'
#' @export
AbstractSerializableSettings <- R6::R6Class(
  "AbstractSerializableSettings",
  public = list(
    #' @description Convert settings to a plain R list.
    #' @return A named list.
    toList = function() {
      nms  <- ls(self)
      vals <- lapply(nms, function(n) self[[n]])
      names(vals) <- nms
      isFunc <- vapply(vals, is.function, logical(1))
      vals[!isFunc]
    },
    #' @description Convert settings to a JSON string.
    #' @return A JSON character string.
    toJson = function() {
      jsonlite::toJSON(self$toList(), auto_unbox = TRUE, null = "null",
                       pretty = TRUE)
    },
    #' @description Validate settings. Throws informative errors on invalid input.
    #' @return Invisibly returns `self`.
    validate = function() invisible(self)
  )
)


# -------------------------------------------------------------------------
# ArmDefinition
# -------------------------------------------------------------------------

#' @title ArmDefinition class
#' @description R6 class representing a treatment strategy arm. Use
#'   [createArmDefinition()] to construct.
#' @export
ArmDefinition <- R6::R6Class(
  "ArmDefinition",
  inherit = AbstractSerializableSettings,
  public = list(
    #' @field armId Unique arm identifier.
    armId        = NULL,
    #' @field strategyType One of `"initiateWithin"` or `"deferBeyond"`.
    strategyType = NULL,
    #' @field gracePeriod Days from cohort start defining the protocol window.
    gracePeriod  = NULL,
    #' @field label Optional display label.
    label        = NULL,

    #' @description Create an `ArmDefinition`.
    #' @param armId Unique arm identifier (character).
    #' @param strategyType `"initiateWithin"` or `"deferBeyond"`.
    #' @param gracePeriod Days from cohort start defining the protocol window.
    #' @param label Optional display label (character or `NULL`).
    initialize = function(armId        = NULL,
                          strategyType = NULL,
                          gracePeriod  = 0,
                          label        = NULL) {
      self$armId        <- armId
      self$strategyType <- strategyType
      self$gracePeriod  <- gracePeriod
      self$label        <- label
      self$validate()
      invisible(self)
    },

    #' @description Validate the arm definition.
    #' @return Invisibly returns `self`.
    validate = function() {
      errs <- checkmate::makeAssertCollection()
      checkmate::assertCharacter(self$armId, len = 1, min.chars = 1,
                                 any.missing = FALSE, add = errs)
      checkmate::assertChoice(self$strategyType,
                              choices = c("initiateWithin", "deferBeyond"),
                              add = errs)
      checkmate::assertNumeric(self$gracePeriod, len = 1, lower = 0,
                               any.missing = FALSE, add = errs)
      checkmate::assertCharacter(self$label, len = 1, null.ok = TRUE, add = errs)
      checkmate::reportAssertions(errs)
      invisible(self)
    }
  )
)

#' Create an arm definition
#'
#' Defines one of the two treatment strategies in a clone-censor-weight
#' analysis.
#'
#' @param armId      A unique character identifier for this arm.
#' @param strategyType One of `"initiateWithin"` or `"deferBeyond"`.
#'   `"initiateWithin"`: the clone is artificially censored if treatment has
#'   **not** occurred by day `gracePeriod`.
#'   `"deferBeyond"`: the clone is artificially censored if treatment occurs
#'   **within** `gracePeriod` days; treatment after that window is permitted.
#' @param gracePeriod Number of days from cohort start defining the protocol
#'   window. Must be `>= 0`.
#' @param label Optional display label used in plots and tables.
#'
#' @return An `ArmDefinition` object.
#'
#' @examples
#' armVaccinate <- createArmDefinition(
#'   armId        = "rzv_within90",
#'   strategyType = "initiateWithin",
#'   gracePeriod  = 90,
#'   label        = "RZV within 90 days"
#' )
#' armDeferred <- createArmDefinition(
#'   armId        = "no_rzv_within90",
#'   strategyType = "deferBeyond",
#'   gracePeriod  = 90,
#'   label        = "No RZV within 90 days"
#' )
#'
#' @export
createArmDefinition <- function(armId, strategyType, gracePeriod = 0,
                                label = NULL) {
  ArmDefinition$new(armId = armId, strategyType = strategyType,
                    gracePeriod = gracePeriod, label = label)
}


# -------------------------------------------------------------------------
# GetDbCloneCensorWeightDataArgs
# -------------------------------------------------------------------------

#' @title GetDbCloneCensorWeightDataArgs class
#' @description R6 class for data extraction arguments. Use
#'   [createGetDbCloneCensorWeightDataArgs()] to construct.
#' @export
GetDbCloneCensorWeightDataArgs <- R6::R6Class(
  "GetDbCloneCensorWeightDataArgs",
  inherit = AbstractSerializableSettings,
  public = list(
    #' @field washoutPeriod Required days of prior observation.
    washoutPeriod   = NULL,
    #' @field covariateSettings A FeatureExtraction covariate settings object.
    covariateSettings = NULL,
    #' @field maxCohortSize Maximum eligible subjects; 0 means no limit.
    maxCohortSize   = NULL,

    #' @description Create `GetDbCloneCensorWeightDataArgs`.
    #' @param washoutPeriod Required days of prior observation.
    #' @param covariateSettings A FeatureExtraction covariate settings object.
    #' @param maxCohortSize Maximum eligible subjects (0 = no limit).
    initialize = function(washoutPeriod   = 365,
                          covariateSettings = NULL,
                          maxCohortSize   = 0) {
      self$washoutPeriod    <- washoutPeriod
      self$covariateSettings <- covariateSettings
      self$maxCohortSize    <- maxCohortSize
      invisible(self)
    },

    #' @description Convert to list, with special handling for covariateSettings.
    #' @return A named list.
    toList = function() {
      result <- list(
        washoutPeriod = self$washoutPeriod,
        maxCohortSize = self$maxCohortSize
      )
      if (!is.null(self$covariateSettings)) {
        result[["covariateSettings"]] <-
          .convertAttrToMember(self$covariateSettings)
      }
      result
    },

    #' @description Validate.
    #' @return Invisibly returns `self`.
    validate = function() {
      errs <- checkmate::makeAssertCollection()
      checkmate::assertNumeric(self$washoutPeriod, len = 1, lower = 0,
                               any.missing = FALSE, add = errs)
      checkmate::assertNumeric(self$maxCohortSize, len = 1, lower = 0,
                               any.missing = FALSE, add = errs)
      checkmate::reportAssertions(errs)
      invisible(self)
    }
  )
)

#' Create data extraction arguments
#'
#' @param washoutPeriod Required days of prior observation before cohort entry.
#'   Subjects without this lookback are excluded.
#' @param covariateSettings A FeatureExtraction covariate settings object
#'   (e.g. from [FeatureExtraction::createDefaultCovariateSettings()]). Used
#'   to build the baseline covariate matrix for IPCW models.
#' @param maxCohortSize Maximum number of eligible subjects to include. If the
#'   eligible population exceeds this, a random sample is drawn. `0` means no
#'   limit.
#'
#' @return A `GetDbCloneCensorWeightDataArgs` object.
#'
#' @examples
#' \dontrun{
#' args <- createGetDbCloneCensorWeightDataArgs(
#'   washoutPeriod     = 365,
#'   covariateSettings = FeatureExtraction::createDefaultCovariateSettings()
#' )
#' }
#'
#' @export
createGetDbCloneCensorWeightDataArgs <- function(washoutPeriod   = 365,
                                                 covariateSettings = NULL,
                                                 maxCohortSize   = 0) {
  obj <- GetDbCloneCensorWeightDataArgs$new(
    washoutPeriod   = washoutPeriod,
    covariateSettings = covariateSettings,
    maxCohortSize   = maxCohortSize
  )
  obj$validate()
  obj
}


# -------------------------------------------------------------------------
# CreateStudyPopulationArgs
# -------------------------------------------------------------------------

#' @title CreateStudyPopulationArgs class
#' @description R6 class for study population arguments. Use
#'   [createCreateStudyPopulationArgs()] to construct.
#' @export
CreateStudyPopulationArgs <- R6::R6Class(
  "CreateStudyPopulationArgs",
  inherit = AbstractSerializableSettings,
  public = list(
    #' @field removeSubjectsWithPriorOutcome Remove subjects with prior outcome.
    removeSubjectsWithPriorOutcome = NULL,
    #' @field priorOutcomeLookback Days to look back for prior outcome.
    priorOutcomeLookback           = NULL,
    #' @field minDaysAtRisk Minimum days of follow-up required.
    minDaysAtRisk                  = NULL,
    #' @field riskWindowStart Days from anchor to start of risk window.
    riskWindowStart                = NULL,
    #' @field startAnchor Anchor for risk window start.
    startAnchor                    = NULL,
    #' @field riskWindowEnd Days from anchor to end of risk window.
    riskWindowEnd                  = NULL,
    #' @field endAnchor Anchor for risk window end.
    endAnchor                      = NULL,

    #' @description Create `CreateStudyPopulationArgs`.
    #' @param removeSubjectsWithPriorOutcome Logical.
    #' @param priorOutcomeLookback Days.
    #' @param minDaysAtRisk Minimum follow-up days.
    #' @param riskWindowStart Days from anchor.
    #' @param startAnchor `"cohort start"` or `"cohort end"`.
    #' @param riskWindowEnd Days from anchor.
    #' @param endAnchor `"cohort start"` or `"cohort end"`.
    initialize = function(removeSubjectsWithPriorOutcome = TRUE,
                          priorOutcomeLookback           = 99999,
                          minDaysAtRisk                  = 1,
                          riskWindowStart                = 0,
                          startAnchor                    = "cohort start",
                          riskWindowEnd                  = 0,
                          endAnchor                      = "cohort end") {
      self$removeSubjectsWithPriorOutcome <- removeSubjectsWithPriorOutcome
      self$priorOutcomeLookback           <- priorOutcomeLookback
      self$minDaysAtRisk                  <- minDaysAtRisk
      self$riskWindowStart                <- riskWindowStart
      self$startAnchor                    <- startAnchor
      self$riskWindowEnd                  <- riskWindowEnd
      self$endAnchor                      <- endAnchor
      invisible(self)
    },

    #' @description Validate.
    #' @return Invisibly returns `self`.
    validate = function() {
      errs <- checkmate::makeAssertCollection()
      checkmate::assertLogical(self$removeSubjectsWithPriorOutcome, len = 1,
                               add = errs)
      checkmate::assertNumeric(self$priorOutcomeLookback, len = 1, lower = 0,
                               add = errs)
      checkmate::assertNumeric(self$minDaysAtRisk, len = 1, lower = 0,
                               add = errs)
      checkmate::assertNumeric(self$riskWindowStart, len = 1, add = errs)
      checkmate::assertChoice(self$startAnchor,
                              c("cohort start", "cohort end"), add = errs)
      checkmate::assertNumeric(self$riskWindowEnd, len = 1, add = errs)
      checkmate::assertChoice(self$endAnchor,
                              c("cohort start", "cohort end"), add = errs)
      checkmate::reportAssertions(errs)
      invisible(self)
    }
  )
)

#' Create study population arguments
#'
#' @param removeSubjectsWithPriorOutcome If `TRUE`, subjects with the outcome
#'   in the `priorOutcomeLookback` window before cohort entry are excluded.
#' @param priorOutcomeLookback Days before cohort entry to search for prior
#'   outcome. Only used when `removeSubjectsWithPriorOutcome = TRUE`.
#' @param minDaysAtRisk Subjects with fewer days of follow-up are excluded.
#' @param riskWindowStart Days from `startAnchor` to the start of the risk
#'   window. Use `0` to start at cohort entry.
#' @param startAnchor `"cohort start"` or `"cohort end"`.
#' @param riskWindowEnd Days from `endAnchor` to the end of the risk window.
#'   Use `0` to end at the natural observation end.
#' @param endAnchor `"cohort start"` or `"cohort end"`.
#'
#' @return A `CreateStudyPopulationArgs` object.
#'
#' @examples
#' args <- createCreateStudyPopulationArgs(
#'   removeSubjectsWithPriorOutcome = TRUE,
#'   priorOutcomeLookback           = 365,
#'   minDaysAtRisk                  = 1,
#'   riskWindowEnd                  = 730,
#'   endAnchor                      = "cohort start"
#' )
#'
#' @export
createCreateStudyPopulationArgs <- function(
    removeSubjectsWithPriorOutcome = TRUE,
    priorOutcomeLookback           = 99999,
    minDaysAtRisk                  = 1,
    riskWindowStart                = 0,
    startAnchor                    = "cohort start",
    riskWindowEnd                  = 0,
    endAnchor                      = "cohort end") {
  obj <- CreateStudyPopulationArgs$new(
    removeSubjectsWithPriorOutcome = removeSubjectsWithPriorOutcome,
    priorOutcomeLookback           = priorOutcomeLookback,
    minDaysAtRisk                  = minDaysAtRisk,
    riskWindowStart                = riskWindowStart,
    startAnchor                    = startAnchor,
    riskWindowEnd                  = riskWindowEnd,
    endAnchor                      = endAnchor
  )
  obj$validate()
  obj
}


# -------------------------------------------------------------------------
# CreateClonesAndCensorArgs
# -------------------------------------------------------------------------

#' @title CreateClonesAndCensorArgs class
#' @description R6 class for cloning arguments. Use
#'   [createCreateClonesAndCensorArgs()] to construct.
#' @export
CreateClonesAndCensorArgs <- R6::R6Class(
  "CreateClonesAndCensorArgs",
  inherit = AbstractSerializableSettings,
  public = list(
    #' @field arms List of exactly two `ArmDefinition` objects.
    arms           = NULL,
    #' @field intervalLength Person-time interval length in days.
    intervalLength = NULL,

    #' @description Create `CreateClonesAndCensorArgs`.
    #' @param arms List of exactly two `ArmDefinition` objects.
    #' @param intervalLength Person-time interval length in days.
    initialize = function(arms = NULL, intervalLength = 30) {
      self$arms           <- arms
      self$intervalLength <- intervalLength
      invisible(self)
    },

    #' @description Convert to list, serializing each arm.
    #' @return A named list.
    toList = function() {
      list(
        arms           = lapply(self$arms, function(a) a$toList()),
        intervalLength = self$intervalLength
      )
    },

    #' @description Validate.
    #' @return Invisibly returns `self`.
    validate = function() {
      errs <- checkmate::makeAssertCollection()
      checkmate::assertList(self$arms, len = 2, any.missing = FALSE, add = errs)
      checkmate::assertNumeric(self$intervalLength, len = 1, lower = 1,
                               add = errs)
      checkmate::reportAssertions(errs)
      if (!all(vapply(self$arms, inherits, logical(1), "ArmDefinition"))) {
        stop("Each element of `arms` must be an ArmDefinition object.",
             call. = FALSE)
      }
      armIds <- vapply(self$arms, function(a) a$armId, character(1))
      if (anyDuplicated(armIds)) {
        stop("Both arms must have distinct `armId` values.", call. = FALSE)
      }
      invisible(self)
    }
  )
)

#' Create cloning and censoring arguments
#'
#' @param arms A list of exactly two [createArmDefinition()] objects defining
#'   the two treatment strategies. The arms must have distinct `armId` values.
#' @param intervalLength Length in days of each person-time interval used in
#'   the pooled logistic IPCW models. Default is 30 days.
#'
#' @return A `CreateClonesAndCensorArgs` object.
#'
#' @examples
#' arm1 <- createArmDefinition("treated",  "initiateWithin", gracePeriod = 90)
#' arm2 <- createArmDefinition("deferred", "deferBeyond",    gracePeriod = 90)
#' args <- createCreateClonesAndCensorArgs(
#'   arms           = list(arm1, arm2),
#'   intervalLength = 30
#' )
#'
#' @export
createCreateClonesAndCensorArgs <- function(arms, intervalLength = 30) {
  obj <- CreateClonesAndCensorArgs$new(arms = arms, intervalLength = intervalLength)
  obj$validate()
  obj
}


# -------------------------------------------------------------------------
# CreateCensoringWeightsArgs
# -------------------------------------------------------------------------

#' @title CreateCensoringWeightsArgs class
#' @description R6 class for IPCW weight arguments. Use
#'   [createCreateCensoringWeightsArgs()] to construct.
#' @export
CreateCensoringWeightsArgs <- R6::R6Class(
  "CreateCensoringWeightsArgs",
  inherit = AbstractSerializableSettings,
  public = list(
    #' @field artificialCensoringFormula Formula for the artificial-censoring
    #'   IPCW model, or `NULL` to build automatically from covariates.
    artificialCensoringFormula = NULL,
    #' @field classicalCensoringFormula Formula for the classical-censoring
    #'   IPCW model, or `NULL` to build automatically.
    classicalCensoringFormula  = NULL,
    #' @field trim Two-element probability vector for weight trimming.
    trim      = NULL,
    #' @field maxWeight Hard cap applied after trimming.
    maxWeight = NULL,

    #' @description Create `CreateCensoringWeightsArgs`.
    #' @param artificialCensoringFormula Formula or `NULL`.
    #' @param classicalCensoringFormula Formula or `NULL`.
    #' @param trim Two-element probability vector.
    #' @param maxWeight Hard cap on cumulative weights.
    initialize = function(artificialCensoringFormula = NULL,
                          classicalCensoringFormula  = NULL,
                          trim      = c(0.01, 0.99),
                          maxWeight = 10) {
      self$artificialCensoringFormula <- artificialCensoringFormula
      self$classicalCensoringFormula  <- classicalCensoringFormula
      self$trim      <- trim
      self$maxWeight <- maxWeight
      invisible(self)
    },

    #' @description Validate.
    #' @return Invisibly returns `self`.
    validate = function() {
      errs <- checkmate::makeAssertCollection()
      checkmate::assertNumeric(self$trim, len = 2, lower = 0, upper = 1,
                               sorted = TRUE, add = errs)
      checkmate::assertNumeric(self$maxWeight, len = 1, lower = 1,
                               add = errs)
      checkmate::reportAssertions(errs)
      invisible(self)
    }
  )
)

#' Create censoring weight arguments
#'
#' @param artificialCensoringFormula An optional [stats::formula()] for the
#'   pooled logistic model predicting artificial censoring. If `NULL`, the
#'   formula is built automatically from all baseline covariates.
#' @param classicalCensoringFormula An optional [stats::formula()] for the
#'   pooled logistic model predicting classical (loss-to-follow-up) censoring.
#'   If `NULL`, built automatically.
#' @param trim A two-element numeric vector `c(lower, upper)` of quantile
#'   probabilities at which cumulative weights are trimmed. Default
#'   `c(0.01, 0.99)`.
#' @param maxWeight Hard upper cap on cumulative weights applied after
#'   trimming. Default `10`.
#'
#' @return A `CreateCensoringWeightsArgs` object.
#'
#' @examples
#' args <- createCreateCensoringWeightsArgs(
#'   trim      = c(0.01, 0.99),
#'   maxWeight = 10
#' )
#'
#' @export
createCreateCensoringWeightsArgs <- function(artificialCensoringFormula = NULL,
                                             classicalCensoringFormula  = NULL,
                                             trim      = c(0.01, 0.99),
                                             maxWeight = 10) {
  obj <- CreateCensoringWeightsArgs$new(
    artificialCensoringFormula = artificialCensoringFormula,
    classicalCensoringFormula  = classicalCensoringFormula,
    trim      = trim,
    maxWeight = maxWeight
  )
  obj$validate()
  obj
}


# -------------------------------------------------------------------------
# FitOutcomeModelArgs
# -------------------------------------------------------------------------

#' @title FitOutcomeModelArgs class
#' @description R6 class for outcome model arguments. Use
#'   [createFitOutcomeModelArgs()] to construct.
#' @export
FitOutcomeModelArgs <- R6::R6Class(
  "FitOutcomeModelArgs",
  inherit = AbstractSerializableSettings,
  public = list(
    #' @field modelType `"cox"` or `"logistic"`.
    modelType         = NULL,
    #' @field varianceEstimator `"robust"`, `"bootstrap"`, or `"conventional"`.
    varianceEstimator = NULL,
    #' @field nBootstraps Number of bootstrap replicates.
    nBootstraps       = NULL,
    #' @field seed Random seed for bootstrap.
    seed              = NULL,

    #' @description Create `FitOutcomeModelArgs`.
    #' @param modelType Model type.
    #' @param varianceEstimator Variance estimator.
    #' @param nBootstraps Bootstrap replicates.
    #' @param seed Random seed.
    initialize = function(modelType         = "cox",
                          varianceEstimator = "robust",
                          nBootstraps       = 200,
                          seed              = 1) {
      self$modelType         <- modelType
      self$varianceEstimator <- varianceEstimator
      self$nBootstraps       <- nBootstraps
      self$seed              <- seed
      invisible(self)
    },

    #' @description Validate.
    #' @return Invisibly returns `self`.
    validate = function() {
      errs <- checkmate::makeAssertCollection()
      checkmate::assertChoice(self$modelType, c("cox", "logistic"), add = errs)
      checkmate::assertChoice(self$varianceEstimator,
                              c("robust", "bootstrap", "conventional"),
                              add = errs)
      checkmate::assertCount(self$nBootstraps, positive = TRUE, add = errs)
      checkmate::assertIntegerish(self$seed, len = 1, add = errs)
      checkmate::reportAssertions(errs)
      invisible(self)
    }
  )
)

#' Create outcome model arguments
#'
#' @param modelType `"cox"` for a weighted Cox proportional hazards model, or
#'   `"logistic"` for a weighted pooled logistic model.
#' @param varianceEstimator Method for computing standard errors that account
#'   for within-subject clone correlation:
#'   - `"robust"`: Lin-Wei-Ying sandwich estimator, clustering on subject.
#'     Fast and appropriate for most studies.
#'   - `"bootstrap"`: Subject-level bootstrap. Slower but makes fewer
#'     distributional assumptions.
#'   - `"conventional"`: Ignores clustering. Not recommended for publication.
#' @param nBootstraps Number of bootstrap replicates. Only used when
#'   `varianceEstimator = "bootstrap"`. Default `200`.
#' @param seed Random seed for the bootstrap. Default `1`.
#'
#' @return A `FitOutcomeModelArgs` object.
#'
#' @examples
#' args <- createFitOutcomeModelArgs(
#'   modelType         = "cox",
#'   varianceEstimator = "robust"
#' )
#'
#' @export
createFitOutcomeModelArgs <- function(modelType         = "cox",
                                      varianceEstimator = "robust",
                                      nBootstraps       = 200,
                                      seed              = 1) {
  obj <- FitOutcomeModelArgs$new(
    modelType         = modelType,
    varianceEstimator = varianceEstimator,
    nBootstraps       = nBootstraps,
    seed              = seed
  )
  obj$validate()
  obj
}


# -------------------------------------------------------------------------
# CcwAnalysis
# -------------------------------------------------------------------------

#' @title CcwAnalysis class
#' @description R6 class representing a complete CCW analysis specification.
#'   Use [createCcwAnalysis()] to construct.
#' @export
CcwAnalysis <- R6::R6Class(
  "CcwAnalysis",
  inherit = AbstractSerializableSettings,
  public = list(
    #' @field analysisId Unique analysis identifier.
    analysisId                     = NULL,
    #' @field description Human-readable description.
    description                    = NULL,
    #' @field getDbCloneCensorWeightDataArgs Data extraction args.
    getDbCloneCensorWeightDataArgs = NULL,
    #' @field createStudyPopulationArgs Study population args.
    createStudyPopulationArgs      = NULL,
    #' @field createClonesAndCensorArgs Cloning and censoring args.
    createClonesAndCensorArgs      = NULL,
    #' @field createCensoringWeightsArgs IPCW weight args.
    createCensoringWeightsArgs     = NULL,
    #' @field fitOutcomeModelArgs Outcome model args.
    fitOutcomeModelArgs            = NULL,

    #' @description Create a `CcwAnalysis`.
    #' @param analysisId Unique character identifier.
    #' @param description Human-readable description.
    #' @param getDbCloneCensorWeightDataArgs A
    #'   `GetDbCloneCensorWeightDataArgs` object.
    #' @param createStudyPopulationArgs A `CreateStudyPopulationArgs` object.
    #' @param createClonesAndCensorArgs A `CreateClonesAndCensorArgs` object.
    #' @param createCensoringWeightsArgs A `CreateCensoringWeightsArgs` object.
    #' @param fitOutcomeModelArgs A `FitOutcomeModelArgs` object.
    initialize = function(analysisId                     = NULL,
                          description                    = "",
                          getDbCloneCensorWeightDataArgs = NULL,
                          createStudyPopulationArgs      = NULL,
                          createClonesAndCensorArgs      = NULL,
                          createCensoringWeightsArgs     = NULL,
                          fitOutcomeModelArgs            = NULL) {
      self$analysisId                     <- analysisId
      self$description                    <- description
      self$getDbCloneCensorWeightDataArgs <- getDbCloneCensorWeightDataArgs
      self$createStudyPopulationArgs      <- createStudyPopulationArgs
      self$createClonesAndCensorArgs      <- createClonesAndCensorArgs
      self$createCensoringWeightsArgs     <- createCensoringWeightsArgs
      self$fitOutcomeModelArgs            <- fitOutcomeModelArgs
      invisible(self)
    },

    #' @description Convert to list, recursively serializing nested args.
    #' @return A named list.
    toList = function() {
      asList <- function(x) if (is.null(x)) NULL else x$toList()
      list(
        analysisId                     = self$analysisId,
        description                    = self$description,
        getDbCloneCensorWeightDataArgs = asList(self$getDbCloneCensorWeightDataArgs),
        createStudyPopulationArgs      = asList(self$createStudyPopulationArgs),
        createClonesAndCensorArgs      = asList(self$createClonesAndCensorArgs),
        createCensoringWeightsArgs     = asList(self$createCensoringWeightsArgs),
        fitOutcomeModelArgs            = asList(self$fitOutcomeModelArgs)
      )
    },

    #' @description Validate.
    #' @return Invisibly returns `self`.
    validate = function() {
      errs <- checkmate::makeAssertCollection()
      checkmate::assertCharacter(self$analysisId, len = 1, min.chars = 1,
                                 any.missing = FALSE, add = errs)
      checkmate::reportAssertions(errs)
      for (field in c("getDbCloneCensorWeightDataArgs",
                      "createStudyPopulationArgs",
                      "createClonesAndCensorArgs",
                      "createCensoringWeightsArgs",
                      "fitOutcomeModelArgs")) {
        if (!is.null(self[[field]])) self[[field]]$validate()
      }
      invisible(self)
    }
  )
)

#' Create a CCW analysis specification
#'
#' Combines all step-level argument objects into a single, serializable
#' analysis specification that can be passed to [runCcwAnalyses()].
#'
#' @param analysisId A unique character identifier for this analysis.
#' @param description A human-readable description.
#' @param getDbCloneCensorWeightDataArgs A
#'   [createGetDbCloneCensorWeightDataArgs()] object.
#' @param createStudyPopulationArgs A [createCreateStudyPopulationArgs()]
#'   object.
#' @param createClonesAndCensorArgs A [createCreateClonesAndCensorArgs()]
#'   object.
#' @param createCensoringWeightsArgs A [createCreateCensoringWeightsArgs()]
#'   object.
#' @param fitOutcomeModelArgs A [createFitOutcomeModelArgs()] object.
#'
#' @return A `CcwAnalysis` object.
#'
#' @examples
#' arm1 <- createArmDefinition("treated",  "initiateWithin", gracePeriod = 90)
#' arm2 <- createArmDefinition("deferred", "deferBeyond",    gracePeriod = 90)
#'
#' analysis <- createCcwAnalysis(
#'   analysisId                 = "main",
#'   description                = "Primary analysis, 90-day grace period",
#'   createClonesAndCensorArgs  = createCreateClonesAndCensorArgs(list(arm1, arm2)),
#'   createStudyPopulationArgs  = createCreateStudyPopulationArgs(riskWindowEnd = 730,
#'                                                                endAnchor = "cohort start"),
#'   createCensoringWeightsArgs = createCreateCensoringWeightsArgs(),
#'   fitOutcomeModelArgs        = createFitOutcomeModelArgs()
#' )
#'
#' @export
createCcwAnalysis <- function(analysisId,
                              description                    = "",
                              getDbCloneCensorWeightDataArgs = NULL,
                              createStudyPopulationArgs      = NULL,
                              createClonesAndCensorArgs      = NULL,
                              createCensoringWeightsArgs     = NULL,
                              fitOutcomeModelArgs            = NULL) {
  obj <- CcwAnalysis$new(
    analysisId                     = analysisId,
    description                    = description,
    getDbCloneCensorWeightDataArgs = getDbCloneCensorWeightDataArgs,
    createStudyPopulationArgs      = createStudyPopulationArgs,
    createClonesAndCensorArgs      = createClonesAndCensorArgs,
    createCensoringWeightsArgs     = createCensoringWeightsArgs,
    fitOutcomeModelArgs            = fitOutcomeModelArgs
  )
  obj$validate()
  obj
}

#' Save a list of CCW analyses to a JSON file
#'
#' @param ccwAnalysisList A list of `CcwAnalysis` objects.
#' @param file Path to the output JSON file.
#'
#' @return Invisibly returns `file`.
#'
#' @examples
#' \dontrun{
#' saveCcwAnalysisList(analysisList, "analyses.json")
#' }
#'
#' @export
saveCcwAnalysisList <- function(ccwAnalysisList, file) {
  checkmate::assertList(ccwAnalysisList, min.len = 1)
  checkmate::assertCharacter(file, len = 1)
  json <- jsonlite::toJSON(lapply(ccwAnalysisList, function(a) a$toList()),
                           auto_unbox = TRUE, null = "null", pretty = TRUE)
  writeLines(json, file)
  invisible(file)
}

#' Load a list of CCW analyses from a JSON file
#'
#' @param file Path to a JSON file created by [saveCcwAnalysisList()].
#'
#' @return A list of `CcwAnalysis` objects.
#'
#' @examples
#' \dontrun{
#' analysisList <- loadCcwAnalysisList("analyses.json")
#' }
#'
#' @export
loadCcwAnalysisList <- function(file) {
  checkmate::assertFileExists(file)
  rawList <- jsonlite::fromJSON(file, simplifyVector = FALSE)
  lapply(rawList, .reconstructCcwAnalysis)
}

.reconstructCcwAnalysis <- function(x) {
  arms <- if (!is.null(x$createClonesAndCensorArgs$arms)) {
    lapply(x$createClonesAndCensorArgs$arms, function(a) {
      createArmDefinition(armId        = a$armId,
                          strategyType = a$strategyType,
                          gracePeriod  = a$gracePeriod %||% 0,
                          label        = a$label)
    })
  } else NULL

  createCcwAnalysis(
    analysisId = x$analysisId,
    description = x$description %||% "",
    getDbCloneCensorWeightDataArgs = if (!is.null(x$getDbCloneCensorWeightDataArgs)) {
      cs <- x$getDbCloneCensorWeightDataArgs$covariateSettings
      createGetDbCloneCensorWeightDataArgs(
        washoutPeriod     = x$getDbCloneCensorWeightDataArgs$washoutPeriod %||% 365,
        covariateSettings = if (!is.null(cs)) .convertMemberToAttr(cs) else NULL,
        maxCohortSize     = x$getDbCloneCensorWeightDataArgs$maxCohortSize %||% 0
      )
    } else NULL,
    createStudyPopulationArgs = if (!is.null(x$createStudyPopulationArgs)) {
      do.call(createCreateStudyPopulationArgs, x$createStudyPopulationArgs)
    } else NULL,
    createClonesAndCensorArgs = if (!is.null(x$createClonesAndCensorArgs)) {
      createCreateClonesAndCensorArgs(
        arms           = arms,
        intervalLength = x$createClonesAndCensorArgs$intervalLength %||% 30
      )
    } else NULL,
    createCensoringWeightsArgs = if (!is.null(x$createCensoringWeightsArgs)) {
      createCreateCensoringWeightsArgs(
        trim      = unlist(x$createCensoringWeightsArgs$trim) %||% c(0.01, 0.99),
        maxWeight = x$createCensoringWeightsArgs$maxWeight %||% 10
      )
    } else NULL,
    fitOutcomeModelArgs = if (!is.null(x$fitOutcomeModelArgs)) {
      do.call(createFitOutcomeModelArgs, x$fitOutcomeModelArgs)
    } else NULL
  )
}


# -------------------------------------------------------------------------
# TargetArmOutcome
# -------------------------------------------------------------------------

#' @title TargetArmOutcome class
#' @description R6 class representing a target cohort, exposure cohort pair,
#'   and associated outcomes. Use [createTargetArmOutcome()] to construct.
#' @export
TargetArmOutcome <- R6::R6Class(
  "TargetArmOutcome",
  inherit = AbstractSerializableSettings,
  public = list(
    #' @field targetId Eligible population cohort definition ID.
    targetId                    = NULL,
    #' @field arm1CohortId Exposure cohort ID for arm 1.
    arm1CohortId                = NULL,
    #' @field arm2CohortId Exposure cohort ID for arm 2.
    arm2CohortId                = NULL,
    #' @field outcomeIds Outcome cohort definition IDs.
    outcomeIds                  = NULL,
    #' @field excludedCovariateConceptIds Concept IDs to exclude from covariates.
    excludedCovariateConceptIds = NULL,
    #' @field includedCovariateConceptIds Concept IDs to include in covariates.
    includedCovariateConceptIds = NULL,

    #' @description Create a `TargetArmOutcome`.
    #' @param targetId Eligible population cohort definition ID.
    #' @param arm1CohortId Arm 1 exposure cohort ID.
    #' @param arm2CohortId Arm 2 exposure cohort ID.
    #' @param outcomeIds Outcome cohort IDs.
    #' @param excludedCovariateConceptIds Excluded concept IDs.
    #' @param includedCovariateConceptIds Included concept IDs.
    initialize = function(targetId                    = NULL,
                          arm1CohortId                = NULL,
                          arm2CohortId                = NULL,
                          outcomeIds                  = NULL,
                          excludedCovariateConceptIds = integer(0),
                          includedCovariateConceptIds = integer(0)) {
      self$targetId                    <- targetId
      self$arm1CohortId                <- arm1CohortId
      self$arm2CohortId                <- arm2CohortId
      self$outcomeIds                  <- outcomeIds
      self$excludedCovariateConceptIds <- excludedCovariateConceptIds %||% integer(0)
      self$includedCovariateConceptIds <- includedCovariateConceptIds %||% integer(0)
      invisible(self)
    },

    #' @description Validate.
    #' @return Invisibly returns `self`.
    validate = function() {
      errs <- checkmate::makeAssertCollection()
      checkmate::assertNumeric(self$targetId, len = 1, any.missing = FALSE,
                               add = errs)
      checkmate::assertNumeric(self$arm1CohortId, len = 1, any.missing = FALSE,
                               add = errs)
      checkmate::assertNumeric(self$arm2CohortId, len = 1, any.missing = FALSE,
                               add = errs)
      checkmate::assertNumeric(self$outcomeIds, min.len = 1, any.missing = FALSE,
                               add = errs)
      checkmate::assertNumeric(self$excludedCovariateConceptIds,
                               any.missing = FALSE, add = errs)
      checkmate::assertNumeric(self$includedCovariateConceptIds,
                               any.missing = FALSE, add = errs)
      checkmate::reportAssertions(errs)
      invisible(self)
    }
  )
)

#' Create a target arm outcome specification
#'
#' Defines the cohorts involved in a CCW analysis: the eligible target
#' population, the two exposure cohorts (one per arm), and the outcome cohorts.
#'
#' @param targetId Cohort definition ID for the eligible population (e.g. new
#'   JAK inhibitor users).
#' @param arm1CohortId Cohort definition ID for the arm 1 exposure (e.g. RZV
#'   recipients).
#' @param arm2CohortId Cohort definition ID for the arm 2 exposure. For a
#'   `"deferBeyond"` arm this is the same exposure cohort as arm 1 — the arm
#'   strategy determines how initiation events are interpreted.
#' @param outcomeIds One or more outcome cohort definition IDs. Include both
#'   the primary outcome and negative control outcomes.
#' @param excludedCovariateConceptIds Concept IDs to exclude from the covariate
#'   set (e.g. the treatment drug ingredients).
#' @param includedCovariateConceptIds Concept IDs to restrict the covariate set
#'   to. Usually left empty.
#'
#' @return A `TargetArmOutcome` object.
#'
#' @examples
#' \dontrun{
#' tao <- createTargetArmOutcome(
#'   targetId     = 101,
#'   arm1CohortId = 102,
#'   arm2CohortId = 102,
#'   outcomeIds   = c(201, 301:320)
#' )
#' }
#'
#' @export
createTargetArmOutcome <- function(targetId,
                                   arm1CohortId,
                                   arm2CohortId,
                                   outcomeIds,
                                   excludedCovariateConceptIds = integer(0),
                                   includedCovariateConceptIds = integer(0)) {
  obj <- TargetArmOutcome$new(
    targetId                    = targetId,
    arm1CohortId                = arm1CohortId,
    arm2CohortId                = arm2CohortId,
    outcomeIds                  = outcomeIds,
    excludedCovariateConceptIds = excludedCovariateConceptIds,
    includedCovariateConceptIds = includedCovariateConceptIds
  )
  obj$validate()
  obj
}

#' Save a list of target arm outcomes to a JSON file
#'
#' @param targetArmOutcomeList A list of `TargetArmOutcome` objects.
#' @param file Path to the output JSON file.
#'
#' @return Invisibly returns `file`.
#'
#' @examples
#' \dontrun{
#' saveTargetArmOutcomeList(taoList, "targetArmOutcomes.json")
#' }
#'
#' @export
saveTargetArmOutcomeList <- function(targetArmOutcomeList, file) {
  checkmate::assertList(targetArmOutcomeList, min.len = 1)
  checkmate::assertCharacter(file, len = 1)
  json <- jsonlite::toJSON(lapply(targetArmOutcomeList, function(t) t$toList()),
                           auto_unbox = TRUE, null = "null", pretty = TRUE)
  writeLines(json, file)
  invisible(file)
}

#' Load a list of target arm outcomes from a JSON file
#'
#' @param file Path to a JSON file created by [saveTargetArmOutcomeList()].
#'
#' @return A list of `TargetArmOutcome` objects.
#'
#' @examples
#' \dontrun{
#' taoList <- loadTargetArmOutcomeList("targetArmOutcomes.json")
#' }
#'
#' @export
loadTargetArmOutcomeList <- function(file) {
  checkmate::assertFileExists(file)
  rawList <- jsonlite::fromJSON(file, simplifyVector = FALSE)
  lapply(rawList, function(x) {
    createTargetArmOutcome(
      targetId                    = x$targetId,
      arm1CohortId                = x$arm1CohortId,
      arm2CohortId                = x$arm2CohortId,
      outcomeIds                  = unlist(x$outcomeIds),
      excludedCovariateConceptIds = unlist(x$excludedCovariateConceptIds) %||% integer(0),
      includedCovariateConceptIds = unlist(x$includedCovariateConceptIds) %||% integer(0)
    )
  })
}


# -------------------------------------------------------------------------
# CcwAnalysesSpecifications
# -------------------------------------------------------------------------

#' @title CcwAnalysesSpecifications class
#' @description Top-level R6 class combining analyses and target-arm-outcome
#'   combinations into a single serializable specification. Compatible with
#'   the Strategus module interface without requiring the Strategus package.
#'   Use [createCcwAnalysesSpecifications()] to construct.
#' @export
CcwAnalysesSpecifications <- R6::R6Class(
  "CcwAnalysesSpecifications",
  inherit = AbstractSerializableSettings,
  public = list(
    #' @field ccwAnalysisList List of `CcwAnalysis` objects.
    ccwAnalysisList      = NULL,
    #' @field targetArmOutcomeList List of `TargetArmOutcome` objects.
    targetArmOutcomeList = NULL,

    #' @description Create a `CcwAnalysesSpecifications`.
    #' @param ccwAnalysisList List of `CcwAnalysis` objects.
    #' @param targetArmOutcomeList List of `TargetArmOutcome` objects.
    initialize = function(ccwAnalysisList = NULL, targetArmOutcomeList = NULL) {
      self$ccwAnalysisList      <- ccwAnalysisList
      self$targetArmOutcomeList <- targetArmOutcomeList
      invisible(self)
    },

    #' @description Convert to list.
    #' @return A named list.
    toList = function() {
      list(
        ccwAnalysisList      = lapply(self$ccwAnalysisList, function(a) a$toList()),
        targetArmOutcomeList = lapply(self$targetArmOutcomeList, function(t) t$toList())
      )
    },

    #' @description Validate.
    #' @return Invisibly returns `self`.
    validate = function() {
      errs <- checkmate::makeAssertCollection()
      checkmate::assertList(self$ccwAnalysisList, min.len = 1, add = errs)
      checkmate::assertList(self$targetArmOutcomeList, min.len = 1, add = errs)
      checkmate::reportAssertions(errs)
      analysisIds <- vapply(self$ccwAnalysisList,
                            function(a) a$analysisId, character(1))
      if (anyDuplicated(analysisIds)) {
        stop("All CcwAnalysis objects must have unique analysisId values.",
             call. = FALSE)
      }
      invisible(self)
    }
  )
)

#' Create a CCW analyses specifications object
#'
#' Combines a list of [createCcwAnalysis()] objects and a list of
#' [createTargetArmOutcome()] objects into a top-level specification that
#' can be passed to [runCcwAnalyses()] or serialized for a Strategus module.
#'
#' @param ccwAnalysisList A list of `CcwAnalysis` objects.
#' @param targetArmOutcomeList A list of `TargetArmOutcome` objects.
#'
#' @return A `CcwAnalysesSpecifications` object.
#'
#' @examples
#' \dontrun{
#' specs <- createCcwAnalysesSpecifications(
#'   ccwAnalysisList      = list(analysis1, analysis2),
#'   targetArmOutcomeList = list(tao1)
#' )
#' specs$toJson()
#' }
#'
#' @export
createCcwAnalysesSpecifications <- function(ccwAnalysisList, targetArmOutcomeList) {
  obj <- CcwAnalysesSpecifications$new(
    ccwAnalysisList      = ccwAnalysisList,
    targetArmOutcomeList = targetArmOutcomeList
  )
  obj$validate()
  obj
}

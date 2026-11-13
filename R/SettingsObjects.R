# Copyright 2026 Observational Health Data Sciences and Informatics
#
# This file is part of CloneCensorWeight
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# The code from here to the end of .convertEmptyListToLogical() is copied from CohortMethod
# (R/SettingsObjects.R, Copyright Observational Health Data Sciences and Informatics, Apache 2.0),
# because CohortMethod does not export it. Replace it with the CohortMethod exports once they exist.

.recurseToList <- function(x) {
  if ("toList" %in% names(x)) {
    return(x$toList())
  } else if (is.list(x)) {
    for (i in seq_along(x)) {
      # If member is NULL then assigning NULL will remove it from the list altogether.
      if (!is.null(x[[i]])) {
        x[[i]] <- .recurseToList(x[[i]])
      }
    }
    if (is.object(x)) {
      # jsonlite::toJSON throws error on S3 objects
      class(x) <- NULL
    }
    return(x)
  } else {
    return(x)
  }
}

AbstractSerializableSettings <- R6Class(
  "AbstractSerializableSettings",
  public = list(
    initialize = function(typedList, untypedList, json) {
      if (!missing(typedList)) {
        # Initialize with list where objects already have correct type
        args <- typedList
        requireTyping <- FALSE
      } else if (!missing(untypedList)) {
        # Initialize with list where objects are still lists, amd need to be converted to right type.
        args <- untypedList
        requireTyping <- TRUE
      } else if (!missing(json)) {
        # Initialize with a JSON string
        args <- jsonlite::fromJSON(json, simplifyDataFrame = FALSE)
        requireTyping <- TRUE
      } else {
        stop("Must provide either typedList, untypedList, or json argument to constructor")
      }
      self$fromList(args, requireTyping)
      self$validate()
    },
    toList = function() {
      # Using is(x, "environment") instead of is.environment(x) because former is FALSE for R6 objects.
      asList <- Filter(function(x) !is.function(x) & !is(x, "environment"), as.list(self))
      asList <- .recurseToList(asList)
      return(asList)
    },
    toJson  = function() {
      # Unlike CohortMethod, digits = NA (15 significant digits) instead of jsonlite's default of 4
      # decimal places, so that reloaded settings keep their values.
      jsonlite::toJSON(
        self$toList(),
        auto_unbox = TRUE,
        pretty = TRUE,
        null = "null",
        digits = NA
      )
    },
    fromList = function(list, requireTyping) {
      for (name in names(list)) {
        if (name %in% names(self)) {
          self[[name]] <- list[[name]]
        }
      }
    }
  )
)

.dateCheck <- function(date) {
  if (date != "") {
    tryCatch({
      dateFormatted <- paste0(
        substr(x = date, start = 1, stop = 4), "-",
        substr(x = date, start = 5, stop = 6), "-",
        substr(x = date, start = 7, stop = 8)
      )

      date <- as.Date(dateFormatted)
    }, error = function(e) {
      stop(sprintf("Date: %s (%s) is not valid", date, dateFormatted))
    })
    checkmate::assertDate(date, lower = "1000-01-01", upper = "2999-12-31")
  }
}

.convertAttrToMember <- function(object) {
  if (is.function(object)) {
    return(list(serialized_code = as.character(serialize(object, NULL))))
  } else if (is.list(object)) {
    if (length(object) > 0) {
      for (i in 1:length(object)) {
        if (!is.null(object[[i]])) {
          object[[i]] <- .convertAttrToMember(object[[i]])
        }
      }
    }
    a <- names(attributes(object))
    a <- a[!a %in% c("names", "class")]
    class <- class(object)
    if (length(class) > 1 || class != "list") {
      class(object) <- "list"
      object$attr_class <- class
    }
    if (length(a) > 0) {
      object[paste("attr", a, sep = "_")] <- attributes(object)[a]
    }
  }
  return(object)
}

.convertMemberToAttr <- function(object) {
  if (is.list(object)) {
    if (length(object) > 0) {
      if (length(object) == 1 && !is.null(names(object)) && names(object) == "serialized_code") {
        return(unserialize(as.raw(sapply(object$serialized_code, strtoi, base = 16L))))
      }
      for (i in 1:length(object)) {
        if (!is.null(object[[i]])) {
          object[[i]] <- .convertMemberToAttr(object[[i]])
        }
      }
      attrNames <- names(object)[grep("^attr_", names(object))]
      cleanNames <- gsub("^attr_", "", attrNames)
      if (any(cleanNames == "class")) {
        class(object) <- object$attr_class
        object$attr_class <- NULL
        attrNames <- attrNames[attrNames != "attr_class"]
        cleanNames <- cleanNames[cleanNames != "class"]
      }
      if (any(cleanNames == "row.names") &&  length(object$attr_row.names) ==  0) {
        object$attr_row.names <- as.character(c())
      }
      attributes(object)[cleanNames] <- object[attrNames]
      object[attrNames] <- NULL
    }
  }
  return(object)
}

.convertEmptyListToLogical <- function(object) {
  # To mimic FeatureExtraction covariate settings behavior: empty vectors are logical
  if (is.list(object)) {
    if (length(object) > 0) {
      for (i in 1:length(object)) {
        if (!is.null(object[[i]])) {
          object[[i]] <- .convertEmptyListToLogical(object[[i]])
        }
      }

      idx <-sapply(object, function(x) is.list(x) && length(x) == 0)
      object[idx] <- sapply(object[idx], as.logical)
    }
  }
  return(object)
}

# End of code copied from CohortMethod --------------------------------------------------------------

.treatmentRoles <- c("treatment", "comparatorTreatment")

# Re-types a CohortMethod settings object from an untyped list by calling its public factory.
# Fields the factory does not know are dropped; missing fields take the factory's defaults.
.typeCohortMethodSettings <- function(untypedList, factory) {
  args <- untypedList[intersect(names(untypedList), names(formals(factory)))]
  return(do.call(factory, args))
}

# Treatment strategy -------------------------------------------------------------------------------

#' Create a treatment strategy
#'
#' @details
#' A treatment strategy is one arm of the emulated target trial. Strategies refer to treatments by
#' role (`"treatment"` or `"comparatorTreatment"`), never by cohort ID. The roles are mapped to cohort
#' IDs by [createTargetTreatmentOutcomes()].
#'
#' Strategy types:
#'
#' - `"initiateWithin"`: start the treatment on a day in `graceStart..graceEnd` (days since index).
#' - `"notInitiateWithin"`: do not start the treatment on or before day `graceEnd`.
#' - `"neverInitiate"`: never start the treatment during follow-up.
#'
#' This version only supports `graceStart = 0` and no `prohibitedTreatments` when creating clones.
#'
#' @param strategyId           An integer identifying the strategy. Used for labels and serialization;
#'                             within an analysis the strategies are numbered by their position.
#' @param label                A short label for the strategy, used in output.
#' @param type                 The strategy type: `"initiateWithin"`, `"notInitiateWithin"`, or
#'                             `"neverInitiate"`.
#' @param treatment            The role of the treatment the strategy refers to: `"treatment"` or
#'                             `"comparatorTreatment"`.
#' @param graceStart           First day (since index) of the grace period.
#' @param graceEnd             Last day (since index) of the grace period. Required for
#'                             `"initiateWithin"` and `"notInitiateWithin"`; must be `NULL` for
#'                             `"neverInitiate"`.
#' @param prohibitedTreatments Roles of treatments whose start is a deviation from the strategy, or
#'                             `NULL`.
#'
#' @return
#' An object of type `TreatmentStrategy`.
#'
#' @examples
#' initiate <- createTreatmentStrategy(
#'   strategyId = 1,
#'   label = "Initiate within 7 days",
#'   type = "initiateWithin",
#'   graceEnd = 7
#' )
#' defer <- createTreatmentStrategy(
#'   strategyId = 2,
#'   label = "Do not initiate within 7 days",
#'   type = "notInitiateWithin",
#'   graceEnd = 7
#' )
#'
#' @export
createTreatmentStrategy <- function(strategyId,
                                    label,
                                    type,
                                    treatment = "treatment",
                                    graceStart = 0,
                                    graceEnd = NULL,
                                    prohibitedTreatments = NULL) {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(TreatmentStrategy$new(typedList = args))
}

TreatmentStrategy <- R6Class(
  "TreatmentStrategy",
  inherit = AbstractSerializableSettings,
  public = list(
    strategyId = NULL,
    label = NULL,
    type = NULL,
    treatment = NULL,
    graceStart = NULL,
    graceEnd = NULL,
    prohibitedTreatments = NULL,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertInt(self$strategyId, add = errorMessages)
      checkmate::assertString(self$label, add = errorMessages)
      checkmate::assertChoice(self$type, c("initiateWithin", "notInitiateWithin", "neverInitiate"), add = errorMessages)
      checkmate::assertChoice(self$treatment, .treatmentRoles, add = errorMessages)
      checkmate::assertInt(self$graceStart, lower = 0, add = errorMessages)
      checkmate::assertInt(self$graceEnd, lower = 0, null.ok = TRUE, add = errorMessages)
      checkmate::assertCharacter(self$prohibitedTreatments, any.missing = FALSE, null.ok = TRUE, add = errorMessages)
      checkmate::assertSubset(self$prohibitedTreatments, .treatmentRoles, add = errorMessages)
      checkmate::reportAssertions(collection = errorMessages)
      if (self$type == "neverInitiate") {
        if (!is.null(self$graceEnd)) {
          stop("graceEnd must be NULL for strategy type 'neverInitiate'")
        }
      } else {
        if (is.null(self$graceEnd)) {
          stop(sprintf("graceEnd is required for strategy type '%s'", self$type))
        }
        if (self$graceEnd < self$graceStart) {
          stop("graceEnd must be greater than or equal to graceStart")
        }
      }
      if (self$treatment %in% self$prohibitedTreatments) {
        stop("A strategy cannot prohibit its own treatment")
      }
    }
  )
)

# Data extraction ----------------------------------------------------------------------------------

#' Create a parameter object for extracting clone-censor-weight data
#'
#' @details
#' Create an object defining the parameter values for extracting a [CloneCensorWeightData] object from
#' the database.
#'
#' @param studyStartDate                    A calendar date specifying the minimum date that a cohort
#'                                          index date can appear. Date format is 'yyyymmdd', or `""`
#'                                          for no restriction.
#' @param studyEndDate                      A calendar date specifying the maximum date that a cohort
#'                                          index date can appear. Date format is 'yyyymmdd', or `""`
#'                                          for no restriction.
#' @param washoutPeriod                     The minimum required continuous observation time (days)
#'                                          prior to index date.
#' @param firstOnly                         Keep only the first index date per person?
#' @param priorUseLookback                  How many days before the index date to look for prior use
#'                                          of the treatments. For "no prior use ever", use a number
#'                                          of days longer than any observation, such as `99999`.
#' @param removeSubjectsWithPriorTreatment  Remove persons who used one of the treatments in the
#'                                          `priorUseLookback` window before index?
#' @param maxCohortSize                     If the eligible cohort is larger than this number it will
#'                                          be sampled to this size. `0` means no maximum.
#' @param covariateSettings                 An object of type `covariateSettings` as created using
#'                                          [FeatureExtraction::createCovariateSettings()], or a list of
#'                                          such objects. FeatureExtraction's default windows end on
#'                                          the index date (`endDays = 0`), so the covariates include
#'                                          records of day 0, the first day on which a strategy can
#'                                          start the treatment: for example drugs started and
#'                                          diagnoses recorded on the index date. Exclude the treatment
#'                                          from the covariates: its concepts with their descendants
#'                                          (`excludedCovariateConceptIds` with
#'                                          `addDescendantsToExclude = TRUE`) and the drug classes that
#'                                          contain it, such as its ATC groups, which are ancestors of
#'                                          the ingredient rather than descendants. Otherwise the
#'                                          censoring models can stop with a high-correlation error.
#'
#' @return
#' An object of type `GetDbCloneCensorWeightDataArgs`.
#'
#' @examples
#' args <- createGetDbCloneCensorWeightDataArgs(washoutPeriod = 180)
#'
#' @export
createGetDbCloneCensorWeightDataArgs <- function(studyStartDate = "",
                                                 studyEndDate = "",
                                                 washoutPeriod = 365,
                                                 firstOnly = TRUE,
                                                 priorUseLookback = 365,
                                                 removeSubjectsWithPriorTreatment = TRUE,
                                                 maxCohortSize = 0,
                                                 covariateSettings = FeatureExtraction::createDefaultCovariateSettings()) {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(GetDbCloneCensorWeightDataArgs$new(typedList = args))
}

GetDbCloneCensorWeightDataArgs <- R6Class(
  "GetDbCloneCensorWeightDataArgs",
  inherit = AbstractSerializableSettings,
  public = list(
    studyStartDate = NULL,
    studyEndDate = NULL,
    washoutPeriod = NULL,
    firstOnly = NULL,
    priorUseLookback = NULL,
    removeSubjectsWithPriorTreatment = NULL,
    maxCohortSize = NULL,
    covariateSettings = NULL,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertCharacter(self$studyStartDate, len = 1, add = errorMessages)
      checkmate::assertCharacter(self$studyEndDate, len = 1, add = errorMessages)
      checkmate::assertInt(self$washoutPeriod, lower = 0, add = errorMessages)
      checkmate::assertFlag(self$firstOnly, add = errorMessages)
      checkmate::assertInt(self$priorUseLookback, lower = 0, add = errorMessages)
      checkmate::assertFlag(self$removeSubjectsWithPriorTreatment, add = errorMessages)
      checkmate::assertInt(self$maxCohortSize, lower = 0, add = errorMessages)
      checkmate::assertList(self$covariateSettings, add = errorMessages)
      checkmate::reportAssertions(collection = errorMessages)
      .dateCheck(self$studyStartDate)
      .dateCheck(self$studyEndDate)
    },
    toList = function() {
      asList <- super$toList()
      asList$covariateSettings <- .convertAttrToMember(self$covariateSettings)
      return(asList)
    },
    fromList = function(list, requireTyping) {
      super$fromList(list)
      if (requireTyping) {
        self$covariateSettings <- .convertEmptyListToLogical(.convertMemberToAttr(list$covariateSettings))
      }
    }
  )
)

# Cloning ------------------------------------------------------------------------------------------

#' Create a parameter object for creating clones
#'
#' @details
#' Create an object defining the parameter values for cloning. Each eligible person is cloned once
#' per strategy and followed until the first of outcome, death, deviation from the strategy,
#' observation end, or risk-window end.
#'
#' Same-day ordering (`sameDayOrdering`):
#'
#' - `"eventsFirst"`: on one day the order is outcome, death, deviation, observation end, risk-window
#'   end. An outcome on the day of a deviation counts.
#' - `"deviationFirst"`: deviation, outcome, death, observation end, risk-window end.
#'
#' @param strategies      A list of exactly two objects of type `TreatmentStrategy` as created using
#'                        [createTreatmentStrategy()]. The first is the target strategy, the second
#'                        the comparator (reference) strategy.
#' @param sameDayOrdering How to order events that happen on the same day: `"eventsFirst"` or
#'                        `"deviationFirst"`. See details.
#' @param competingEvents Treat death as a competing event (`"death"`) or not (`"none"`).
#'
#' @return
#' An object of type `CreateClonesArgs`.
#'
#' @examples
#' args <- createCreateClonesArgs(
#'   strategies = list(
#'     createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
#'     createTreatmentStrategy(2, "Do not initiate within 7 days", "notInitiateWithin", graceEnd = 7)
#'   )
#' )
#'
#' @export
createCreateClonesArgs <- function(strategies,
                                   sameDayOrdering = "eventsFirst",
                                   competingEvents = "death") {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(CreateClonesArgs$new(typedList = args))
}

CreateClonesArgs <- R6Class(
  "CreateClonesArgs",
  inherit = AbstractSerializableSettings,
  public = list(
    strategies = NULL,
    sameDayOrdering = NULL,
    competingEvents = NULL,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertList(self$strategies, len = 2, add = errorMessages)
      for (i in seq_along(self$strategies)) {
        checkmate::assertR6(self$strategies[[i]], "TreatmentStrategy", add = errorMessages)
      }
      checkmate::assertChoice(self$sameDayOrdering, c("eventsFirst", "deviationFirst"), add = errorMessages)
      checkmate::assertChoice(self$competingEvents, c("death", "none"), add = errorMessages)
      checkmate::reportAssertions(collection = errorMessages)
      if (self$strategies[[1]]$strategyId == self$strategies[[2]]$strategyId) {
        stop("The two strategies must have different strategyId values")
      }
    },
    fromList = function(list, requireTyping) {
      super$fromList(list)
      if (requireTyping) {
        for (i in seq_along(self$strategies)) {
          self$strategies[[i]] <- TreatmentStrategy$new(untypedList = self$strategies[[i]])
        }
      }
    }
  )
)

# Censoring models ---------------------------------------------------------------------------------

#' Create a time basis
#'
#' @details
#' Defines how time enters a pooled (discrete-time) model. Every model that uses a basis also has an
#' intercept, so `"indicator"` has one column per interval except the first (the reference), and a
#' single interval gives no basis columns.
#'
#' @param type           The basis: `"indicator"` (one indicator per interval), `"naturalSpline"`,
#'                       `"linear"`, or `"polynomial"`.
#' @param intervalLength The length of each interval in days.
#' @param df             Degrees of freedom of the natural spline.
#' @param degree         Degree of the polynomial (`stats::poly(k, degree, raw = raw)` of the 0-based
#'                       interval index `k`).
#' @param raw            Use raw (not orthogonal) polynomials?
#'
#' @return
#' An object of type `TimeBasis`.
#'
#' @examples
#' timeBasis <- createTimeBasis(type = "naturalSpline", intervalLength = 7, df = 3)
#'
#' @export
createTimeBasis <- function(type = "indicator",
                            intervalLength = 1,
                            df = 4,
                            degree = 2,
                            raw = TRUE) {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(TimeBasis$new(typedList = args))
}

TimeBasis <- R6Class(
  "TimeBasis",
  inherit = AbstractSerializableSettings,
  public = list(
    type = NULL,
    intervalLength = NULL,
    df = NULL,
    degree = 2,
    raw = TRUE,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertChoice(self$type, c("indicator", "naturalSpline", "linear", "polynomial"), add = errorMessages)
      checkmate::assertInt(self$intervalLength, lower = 1, add = errorMessages)
      checkmate::assertInt(self$df, lower = 1, add = errorMessages)
      checkmate::assertInt(self$degree, lower = 1, add = errorMessages)
      checkmate::assertLogical(self$raw, len = 1, any.missing = FALSE, add = errorMessages)
      checkmate::reportAssertions(collection = errorMessages)
    }
  )
)

#' Create a parameter object for fitting the censoring models
#'
#' @details
#' Create an object defining the parameter values for fitting the models of artificial censoring
#' (deviation from a strategy). For `"initiateWithin"` strategies a binary model at the end of the grace
#' period is always used. For `"notInitiateWithin"` and `"neverInitiate"` strategies the hazard of
#' initiation is modelled with `hazardModelType`.
#'
#' The default hazard model, pooled logistic regression with a daily indicator time basis, is the
#' exact discrete-time hazard of initiation. The Cox model treats the initiation days as continuous
#' times with Breslow ties, so its survival of non-initiation, `exp(-Lambda0(d + 1) exp(X beta))`, is
#' biased when persons initiate on the same day: `exp(-h)` is larger than `1 - h` for a daily
#' probability `h`, and the error builds up over the days of the horizon. The bias of the Cox weights
#' therefore grows with the horizon and with each person's daily probability of initiating; it is
#' largest for `"neverInitiate"` strategies. For long horizons, prefer the pooled logistic model, with
#' a longer `intervalLength` if daily intervals are too costly.
#'
#' `maxCohortSizeForFitting` applies only to the binary models at the end of the grace period
#' (through [CohortMethod::createPs()]). In this version the hazard models are always fitted to their
#' whole risk set.
#'
#' **Weights of `"notInitiateWithin"` after the grace period.** With `postGraceModel = "binary"` (the
#' default), the weight of a `"notInitiateWithin(G)"` clone after day `G` is `1 / (1 - pi(X))`, where
#' `pi(X)` is the probability of starting by day `G` from the binary model at `G`, the same model the
#' `"initiateWithin(G)"` strategy uses (the fit is shared when the analysis has such a strategy for the
#' same treatment, and made otherwise). The initiation hazard model still gives the weights on days up
#' to `G`, and all weights of `"neverInitiate"` strategies. The binary model is fitted on the persons
#' under follow-up at `G`, so `1 - pi(X)` is the probability of not starting by `G` among persons with
#' no outcome or death by `G`. It equals the probability `S(G | X)` of not initiating on days `0..G` of
#' the hazard model when starting the treatment during the grace period does not change the risk of
#' the outcome or death during the grace period. Unlike `S(G | X)`, it does not assume that the
#' covariate effects on initiation are the same on every day of the grace period (with the default
#' `covariateTimeInteraction = "none"`): when a covariate changes *when* persons start (for example
#' on day 0 rather than later), `1 - S(G | X)` is miscalibrated and the weights `1 / S(G | X)` do not
#' balance the covariates after `G`; see [computeCensoringModelCalibration()]. With
#' `postGraceModel = "hazard"`, the weight after `G` is `1 / S(G | X)`.
#'
#' **Covariate-by-time interactions.** With `covariateTimeInteraction = "day0"`, the pooled logistic
#' initiation hazard model also has, for each covariate, its product with an indicator of the first
#' interval of the time basis (day 0 with the default daily intervals), so the covariates can change
#' initiation on day 0 differently from later days. With `"interval"`, it has each covariate's product
#' with each time-basis column (with the indicator basis: with the indicator of each interval except
#' the reference), so the covariate effects can differ in every interval. The interaction columns are
#' penalized by the prior like the covariates; the time-basis columns are never penalized. The Cox
#' model has no interactions.
#'
#' @param hazardModelType          The hazard model: `"logistic"` (pooled logistic, the default) or
#'                                 `"cox"` (see details for its ties caveat).
#' @param timeBasis                An object of type `TimeBasis` as created using
#'                                 [createTimeBasis()], used by the pooled logistic model. The default
#'                                 is a daily indicator basis.
#' @param covariateTimeInteraction Covariate-by-time interactions in the pooled logistic initiation
#'                                 hazard model: `"none"`, `"day0"` (each covariate times an indicator
#'                                 of the first interval, day 0 with daily intervals), or `"interval"`
#'                                 (each covariate times each time-basis column). See details. Must be
#'                                 `"none"` with `hazardModelType = "cox"`.
#' @param postGraceModel           The model that gives the weights of `"notInitiateWithin"` strategies
#'                                 after the grace period: `"binary"` (`1 / (1 - pi(X))` from the
#'                                 binary model at the end of the grace period) or `"hazard"`
#'                                 (`1 / S(G | X)` from the initiation hazard model). See details.
#' @param censoringMechanisms      The censoring mechanisms to weight for. This version supports only
#'                                 `"artificial"` (deviation from a strategy). Censoring at the end of
#'                                 observation is not weighted, so it is assumed to be independent of
#'                                 the outcome. Risks are biased in data where the end of observation
#'                                 depends on the outcome, for example when persons who die are
#'                                 observed until death but survivors only until their last record.
#' @param useCovariates            Use the baseline covariates in the censoring models?
#' @param excludeCovariateIds      Exclude these covariates from the censoring models.
#' @param includeCovariateIds      Include only these covariates in the censoring models.
#' @param maxCohortSizeForFitting  If the risk set of the binary model of an `"initiateWithin"`
#'                                 strategy is larger than this number, it is sampled to this size
#'                                 before fitting. `0` means no sampling. The hazard models do not
#'                                 sample in this version.
#' @param errorOnHighCorrelation   Check each covariate for high correlation with the response, and
#'                                 throw an error if found?
#' @param stopOnError              If an error occurs while fitting, should the function stop?
#' @param prior                    The prior used to fit the models. See [Cyclops::createPrior()].
#' @param control                  The control object used to fit the models. See
#'                                 [Cyclops::createControl()].
#'
#' @return
#' An object of type `FitCensoringModelsArgs`.
#'
#' @examples
#' args <- createFitCensoringModelsArgs(hazardModelType = "cox")
#'
#' @export
createFitCensoringModelsArgs <- function(hazardModelType = "logistic",
                                         timeBasis = createTimeBasis(),
                                         covariateTimeInteraction = "none",
                                         postGraceModel = "binary",
                                         censoringMechanisms = "artificial",
                                         useCovariates = TRUE,
                                         excludeCovariateIds = c(),
                                         includeCovariateIds = c(),
                                         maxCohortSizeForFitting = 250000,
                                         errorOnHighCorrelation = TRUE,
                                         stopOnError = TRUE,
                                         prior = Cyclops::createPrior("laplace",
                                                                      exclude = c(0),
                                                                      useCrossValidation = TRUE),
                                         control = Cyclops::createControl(noiseLevel = "silent",
                                                                          cvType = "auto",
                                                                          seed = 1,
                                                                          resetCoefficients = TRUE,
                                                                          tolerance = 2e-07,
                                                                          cvRepetitions = 10,
                                                                          startingVariance = 0.01)) {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(FitCensoringModelsArgs$new(typedList = args))
}

FitCensoringModelsArgs <- R6Class(
  "FitCensoringModelsArgs",
  inherit = AbstractSerializableSettings,
  public = list(
    hazardModelType = NULL,
    timeBasis = NULL,
    # Defaults, so settings saved before these fields existed load:
    covariateTimeInteraction = "none",
    postGraceModel = "binary",
    censoringMechanisms = NULL,
    useCovariates = NULL,
    excludeCovariateIds = NULL,
    includeCovariateIds = NULL,
    maxCohortSizeForFitting = NULL,
    errorOnHighCorrelation = NULL,
    stopOnError = NULL,
    prior = NULL,
    control = NULL,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertChoice(self$hazardModelType, c("cox", "logistic"), add = errorMessages)
      checkmate::assertR6(self$timeBasis, "TimeBasis", add = errorMessages)
      checkmate::assertChoice(self$covariateTimeInteraction, c("none", "day0", "interval"), add = errorMessages)
      checkmate::assertChoice(self$postGraceModel, c("binary", "hazard"), add = errorMessages)
      if (identical(self$hazardModelType, "cox") && !identical(self$covariateTimeInteraction, "none")) {
        errorMessages$push("covariateTimeInteraction must be \"none\" with hazardModelType = \"cox\": covariate-by-time interactions are available only in the pooled logistic initiation hazard model (hazardModelType = \"logistic\")")
      }
      checkmate::assertCharacter(self$censoringMechanisms, min.len = 1, add = errorMessages)
      checkmate::assertSubset(self$censoringMechanisms, "artificial", add = errorMessages)
      checkmate::assertFlag(self$useCovariates, add = errorMessages)
      .assertCovariateId(self$excludeCovariateIds, null.ok = TRUE, add = errorMessages)
      .assertCovariateId(self$includeCovariateIds, null.ok = TRUE, add = errorMessages)
      checkmate::assertInt(self$maxCohortSizeForFitting, lower = 0, add = errorMessages)
      checkmate::assertFlag(self$errorOnHighCorrelation, add = errorMessages)
      checkmate::assertFlag(self$stopOnError, add = errorMessages)
      checkmate::assertClass(self$prior, "cyclopsPrior", add = errorMessages)
      checkmate::assertClass(self$control, "cyclopsControl", add = errorMessages)
      checkmate::reportAssertions(collection = errorMessages)
    },
    fromList = function(list, requireTyping) {
      super$fromList(list)
      if (requireTyping) {
        self$timeBasis <- TimeBasis$new(untypedList = self$timeBasis)
        class(self$control) <- "cyclopsControl"
        class(self$prior) <- "cyclopsPrior"
      }
    }
  )
)

# Weights ------------------------------------------------------------------------------------------

#' Create a weight truncation rule
#'
#' @details
#' Weights are truncated per strategy: every segment weight above the strategy's cap is set to the cap.
#' With `type = "fixed"`, the cap is `upper`. With `type = "percentile"`, the cap is the `upper`
#' percentile (`quantile` type 7) of the clones' maximum weights, taken over:
#'
#' - `basis = "varying"` (the default): the clones under follow-up on day `tau + 1` (or on their
#'   `riskStart`, if later), where `tau` is the largest `graceEnd` of the strategies. These are the
#'   clones whose weight varies. Clones censored during a grace period (for example the non-starters
#'   of an `"initiateWithin"` strategy, whose weight is 1 throughout) do not lower the cap. Without a
#'   grace period, or when no clone of a strategy is under follow-up on that day, all clones of the
#'   strategy are used.
#' - `basis = "allClones"`: all clones of the strategy.
#'
#' When few persons start the treatment by the end of a grace period, most clones of an
#' `"initiateWithin"` strategy keep a weight of 1, and a percentile over all clones lies inside the
#' distribution of the weights that vary: it truncates far more than `1 - upper` of them.
#'
#' @param type  The truncation rule: `"none"`, `"percentile"`, or `"fixed"`.
#' @param upper For `"percentile"`, the percentile in (0, 1). For `"fixed"`, the cap (greater than 1).
#'              Ignored for `"none"`.
#' @param basis For `"percentile"`, the clones over which the percentile is taken: `"varying"` or
#'              `"allClones"` (see details). Ignored for other types.
#'
#' @return
#' An object of type `WeightTruncation`.
#'
#' @examples
#' truncation <- createWeightTruncation(type = "fixed", upper = 20)
#' truncation <- createWeightTruncation(type = "percentile", upper = 0.99, basis = "allClones")
#'
#' @export
createWeightTruncation <- function(type = "percentile",
                                   upper = 0.99,
                                   basis = "varying") {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(WeightTruncation$new(typedList = args))
}

WeightTruncation <- R6Class(
  "WeightTruncation",
  inherit = AbstractSerializableSettings,
  public = list(
    type = NULL,
    upper = NULL,
    # Default, so settings saved before this field existed load:
    basis = "varying",
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertChoice(self$type, c("none", "percentile", "fixed"), add = errorMessages)
      checkmate::assertNumber(self$upper, null.ok = TRUE, add = errorMessages)
      checkmate::assertChoice(self$basis, c("varying", "allClones"), add = errorMessages)
      checkmate::reportAssertions(collection = errorMessages)
      if (self$type == "percentile" && (is.null(self$upper) || self$upper <= 0 || self$upper >= 1)) {
        stop("For percentile truncation, upper must be in (0, 1)")
      }
      if (self$type == "fixed" && (is.null(self$upper) || self$upper <= 1)) {
        stop("For fixed truncation, upper must be greater than 1")
      }
    }
  )
)

#' Create a parameter object for computing the censoring weights
#'
#' @param stabilized Stabilize the weights with a numerator model of time only?
#' @param truncation An object of type `WeightTruncation` as created using
#'                   [createWeightTruncation()].
#'
#' @return
#' An object of type `ComputeWeightsArgs`.
#'
#' @examples
#' args <- createComputeWeightsArgs(stabilized = FALSE)
#'
#' @export
createComputeWeightsArgs <- function(stabilized = TRUE,
                                     truncation = createWeightTruncation()) {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(ComputeWeightsArgs$new(typedList = args))
}

ComputeWeightsArgs <- R6Class(
  "ComputeWeightsArgs",
  inherit = AbstractSerializableSettings,
  public = list(
    stabilized = NULL,
    truncation = NULL,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertFlag(self$stabilized, add = errorMessages)
      checkmate::assertR6(self$truncation, "WeightTruncation", add = errorMessages)
      checkmate::reportAssertions(collection = errorMessages)
    },
    fromList = function(list, requireTyping) {
      super$fromList(list)
      if (requireTyping) {
        self$truncation <- WeightTruncation$new(untypedList = self$truncation)
      }
    }
  )
)

# Outcome model ------------------------------------------------------------------------------------

#' Create a parameter object for the bootstrap
#'
#' @details
#' The bootstrap ([bootstrapCcwAnalysis()]) resamples persons (so all clones of a person stay
#' together) and re-estimates the censoring models and weights in every replicate.
#'
#' @param replicates         The number of bootstrap replicates.
#' @param seed               The random seed.
#' @param ci                 How to compute confidence intervals: `"percentile"` or `"normal"`
#'                           (bias-corrected), as the `"perc"` and `"norm"` types of [boot::boot.ci()].
#' @param fixHyperparameters Reuse the prior variances selected by cross-validation on the full data
#'                           in every replicate?
#'
#' @return
#' An object of type `BootstrapArgs`.
#'
#' @examples
#' args <- createBootstrapArgs(replicates = 500)
#'
#' @export
createBootstrapArgs <- function(replicates = 200,
                                seed = 1,
                                ci = "percentile",
                                fixHyperparameters = TRUE) {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(BootstrapArgs$new(typedList = args))
}

BootstrapArgs <- R6Class(
  "BootstrapArgs",
  inherit = AbstractSerializableSettings,
  public = list(
    replicates = NULL,
    seed = NULL,
    ci = NULL,
    fixHyperparameters = NULL,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertInt(self$replicates, lower = 1, add = errorMessages)
      checkmate::assertInt(self$seed, add = errorMessages)
      checkmate::assertChoice(self$ci, c("percentile", "normal"), add = errorMessages)
      checkmate::assertFlag(self$fixHyperparameters, add = errorMessages)
      checkmate::reportAssertions(collection = errorMessages)
    }
  )
)

#' Create a parameter object for fitting the outcome model
#'
#' @details
#' Create an object defining the parameter values for the weighted (marginal structural) outcome
#' model fitted to the clones.
#'
#' This version requires `useCovariates = FALSE` and `stratified = FALSE`.
#'
#' The pooled logistic model needs weighted clones cut on its interval grid
#' (`computeCensoringWeights(outcomeTimeBasis = timeBasis)`). Its default, daily intervals with a
#' natural spline of time, uses the exact daily weights. With longer intervals, a clone whose
#' follow-up ends inside an interval without the outcome (for example at the end of a grace period, or
#' when it deviates) counts as at risk, without the outcome, for the whole interval, which biases that
#' interval's hazard down. This bias is large when many clones end inside an interval, as at the end of
#' a grace period that does not end an interval of the grid.
#'
#' The default `control` (Newton steps until the relative change of the log likelihood is below
#' `1e-10`) fits the models to the precision of their reference implementations (`survival::coxph()`
#' and `stats::glm()`), also on ill-scaled time bases, where coordinate descent stops early. The
#' criterion is relative, so it is met on large data (millions of weighted rows) and on models with
#' quasi-separated hazards, where a criterion on the gradient can fail to be met at the optimum. If a
#' fit with `control` fails, [fitOutcomeModel()] refits it with Cyclops' default algorithm and
#' convergence criterion, with a warning.
#'
#' The default, `"cox"`, reports the hazard ratio of the weighted Cox model together with the weighted
#' Kaplan-Meier (or, with competing death, Aalen-Johansen) risks, risk differences and risk ratios of
#' `"km"` on the same weighted clones. Risks from the Cox model would assume proportional hazards,
#' which a grace period violates by design: during the grace period the clones of both strategies
#' follow the same course, so their true hazards are equal there. Such risks are biased, also at
#' horizons within the grace period, so they are not reported. The hazard ratio is kept because
#' meta-analysis and empirical calibration use it.
#'
#' The pooled logistic (`"logistic"`) risks are model-based. Without `armByTimeInteraction` they assume
#' an odds ratio that is constant over follow-up and are biased under a grace period for the same
#' reason, so use `armByTimeInteraction = TRUE` with a grace period. With daily intervals, an indicator
#' (saturated) time basis and the interaction give the weighted Kaplan-Meier or Aalen-Johansen risks.
#' A smooth time basis, such as the default natural spline, can still give biased risks with the
#' interaction: each strategy's hazard can change abruptly after the grace period, when the clones of
#' the two strategies start to follow different treatments, and a smooth function of time cannot follow
#' that change.
#'
#' The Cox and pooled logistic models are fitted with the weights chosen by `msmWeights`, whatever
#' [createComputeWeightsArgs()] `stabilized` is ([computeCensoringWeights()] stores what is needed to
#' form both). The default is the unstabilized weights. Under a grace period the stabilization
#' numerators depend on the strategy and the day: after the grace period of an `"initiateWithin"`
#' strategy the weight is multiplied by the share of persons who started by its end, and during the
#' grace period of a strategy of not initiating by the share who have not initiated yet, which falls
#' from day to day. Within a strategy, all clones at risk on a day share the numerator, so it cancels
#' in the weighted Kaplan-Meier and Aalen-Johansen risks. It does not cancel in the Cox partial
#' likelihood, which compares the clones of both strategies in each risk set: it changes how much each
#' strategy, and each period of follow-up, counts. Because a grace period makes the hazards
#' non-proportional (during the grace period the clones of both strategies follow the same course),
#' the hazard ratio fitted with stabilized weights is then a different average over follow-up, not the
#' hazard ratio of the regimes. In a simulation with a 7-day grace period and a regime log hazard ratio
#' of -0.513, the Cox model estimated a log hazard ratio of about -0.355 with stabilized weights and
#' -0.516 with unstabilized weights (means over replicates). The same holds for a pooled logistic model
#' whose hazards are shared across strategies or days. `"stabilized"` reproduces the earlier
#' behaviour, which fitted these models with the stabilized weights of [createComputeWeightsArgs()]
#' `stabilized = TRUE`; its standard errors are usually smaller, but it estimates that other average.
#'
#' @param modelType            The outcome model: `"cox"` (the hazard ratio of a weighted Cox model,
#'                             with the risks of `"km"`), `"logistic"` (pooled logistic), or `"km"`
#'                             (weighted Kaplan-Meier, or Aalen-Johansen with competing events).
#' @param horizons             Days since index at which risks, risk differences and risk ratios are
#'                             reported.
#' @param timeBasis            An object of type `TimeBasis` as created using [createTimeBasis()], used
#'                             by the pooled logistic model.
#' @param armByTimeInteraction Add strategy-by-time interactions to the pooled logistic model? With the
#'                             interactions, the odds ratio reported as `"or"` is the odds ratio where
#'                             the time basis is 0: in the first interval for the indicator, natural
#'                             spline, linear and raw polynomial bases (where clones of both strategies
#'                             often have the same weights and hazards, so it is close to 1), and at no
#'                             particular interval for an orthogonal polynomial basis. It is not a
#'                             summary of the effect over follow-up.
#' @param msmWeights           The weights of the Cox and pooled logistic models: `"unstabilized"` or
#'                             `"stabilized"` (see Details). The weighted Kaplan-Meier and
#'                             Aalen-Johansen estimates (of `"km"`, and the rows at the horizons of
#'                             `"cox"`) do not depend on it.
#' @param varianceEstimator    The variance estimator: `"robust"` (sandwich clustered on person),
#'                             `"bootstrap"` (computed by [bootstrapCcwAnalysis()]), or `"modelBased"`.
#' @param bootstrapArgs        An object of type `BootstrapArgs` as created using
#'                             [createBootstrapArgs()], used by [bootstrapCcwAnalysis()].
#' @param useCovariates        Add the baseline covariates to the outcome model?
#' @param stratified           Fit a stratified model?
#' @param prior                The prior used to fit the model. See [Cyclops::createPrior()]. The target
#'                             indicator, the time basis and their interactions are never regularized
#'                             (they are added to the prior's `exclude`), as CohortMethod excludes the
#'                             treatment variable. This version fits no other covariates, so the prior
#'                             has no effect on the estimates.
#' @param control              The control object used to fit the model. See
#'                             [Cyclops::createControl()].
#'
#' @return
#' An object of type `FitOutcomeModelArgs`.
#'
#' @examples
#' args <- createFitOutcomeModelArgs(modelType = "km", horizons = c(30, 365))
#'
#' @export
createFitOutcomeModelArgs <- function(modelType = "cox",
                                      horizons = c(365),
                                      timeBasis = createTimeBasis(type = "naturalSpline",
                                                                  intervalLength = 1,
                                                                  df = 4),
                                      armByTimeInteraction = FALSE,
                                      msmWeights = "unstabilized",
                                      varianceEstimator = "robust",
                                      bootstrapArgs = createBootstrapArgs(),
                                      useCovariates = FALSE,
                                      stratified = FALSE,
                                      prior = Cyclops::createPrior("none"),
                                      control = Cyclops::createControl(algorithm = "newton",
                                                                       convergenceType = "gradient",
                                                                       tolerance = 1e-10,
                                                                       noiseLevel = "silent")) {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(FitOutcomeModelArgs$new(typedList = args))
}

FitOutcomeModelArgs <- R6Class(
  "FitOutcomeModelArgs",
  inherit = AbstractSerializableSettings,
  public = list(
    modelType = NULL,
    horizons = NULL,
    timeBasis = NULL,
    armByTimeInteraction = NULL,
    msmWeights = NULL,
    varianceEstimator = NULL,
    bootstrapArgs = NULL,
    useCovariates = NULL,
    stratified = NULL,
    prior = NULL,
    control = NULL,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertChoice(self$modelType, c("cox", "logistic", "km"), add = errorMessages)
      checkmate::assertIntegerish(self$horizons, lower = 0, min.len = 1, any.missing = FALSE, add = errorMessages)
      checkmate::assertR6(self$timeBasis, "TimeBasis", add = errorMessages)
      checkmate::assertFlag(self$armByTimeInteraction, add = errorMessages)
      checkmate::assertChoice(self$msmWeights, c("unstabilized", "stabilized"), add = errorMessages)
      checkmate::assertChoice(self$varianceEstimator, c("robust", "bootstrap", "modelBased"), add = errorMessages)
      checkmate::assertR6(self$bootstrapArgs, "BootstrapArgs", add = errorMessages)
      checkmate::assertFlag(self$useCovariates, add = errorMessages)
      checkmate::assertFlag(self$stratified, add = errorMessages)
      checkmate::assertClass(self$prior, "cyclopsPrior", add = errorMessages)
      checkmate::assertClass(self$control, "cyclopsControl", add = errorMessages)
      checkmate::reportAssertions(collection = errorMessages)
      if (self$useCovariates) {
        stop("Not yet supported in this version: useCovariates = TRUE")
      }
      if (self$stratified) {
        stop("Not yet supported in this version: stratified = TRUE")
      }
    },
    fromList = function(list, requireTyping) {
      super$fromList(list)
      if (requireTyping) {
        self$timeBasis <- TimeBasis$new(untypedList = self$timeBasis)
        self$bootstrapArgs <- BootstrapArgs$new(untypedList = self$bootstrapArgs)
        class(self$control) <- "cyclopsControl"
        class(self$prior) <- "cyclopsPrior"
      }
    }
  )
)

# Analysis -----------------------------------------------------------------------------------------

#' Create a clone-censor-weight analysis specification
#'
#' @details
#' Create a set of analysis choices. In this version the specification can be saved and loaded (see
#' [saveCcwAnalysisList()]); the driver that runs many analyses, `runCcwAnalyses()`, comes in a later
#' version.
#'
#' @param analysisId                     An integer that will be used later to refer to this specific
#'                                       set of analysis choices.
#' @param description                    A short description of the analysis.
#' @param getDbCloneCensorWeightDataArgs An object of type `GetDbCloneCensorWeightDataArgs` as created
#'                                       using [createGetDbCloneCensorWeightDataArgs()].
#' @param createStudyPopulationArgs      An object of type `CreateStudyPopulationArgs` as created using
#'                                       [CohortMethod::createCreateStudyPopulationArgs()].
#' @param createClonesArgs               An object of type `CreateClonesArgs` as created using
#'                                       [createCreateClonesArgs()].
#' @param fitCensoringModelsArgs         An object of type `FitCensoringModelsArgs` as created using
#'                                       [createFitCensoringModelsArgs()].
#' @param computeWeightsArgs             An object of type `ComputeWeightsArgs` as created using
#'                                       [createComputeWeightsArgs()].
#' @param computeCovariateBalanceArgs    An object of type `ComputeCovariateBalanceArgs` as created using
#'                                       [CohortMethod::createComputeCovariateBalanceArgs()], or `NULL`
#'                                       to skip covariate balance.
#' @param fitOutcomeModelArgs            An object of type `FitOutcomeModelArgs` as created using
#'                                       [createFitOutcomeModelArgs()].
#'
#' @return
#' An object of type `CcwAnalysis`.
#'
#' @examples
#' analysis <- createCcwAnalysis(
#'   analysisId = 1,
#'   description = "Initiate vs do not initiate within 7 days",
#'   getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(),
#'   createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(
#'     riskWindowEnd = 365,
#'     endAnchor = "cohort start"
#'   ),
#'   createClonesArgs = createCreateClonesArgs(
#'     strategies = list(
#'       createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
#'       createTreatmentStrategy(2, "Do not initiate within 7 days", "notInitiateWithin", graceEnd = 7)
#'     )
#'   ),
#'   fitCensoringModelsArgs = createFitCensoringModelsArgs(),
#'   computeWeightsArgs = createComputeWeightsArgs(),
#'   fitOutcomeModelArgs = createFitOutcomeModelArgs()
#' )
#'
#' @export
createCcwAnalysis <- function(analysisId,
                              description = "",
                              getDbCloneCensorWeightDataArgs,
                              createStudyPopulationArgs,
                              createClonesArgs,
                              fitCensoringModelsArgs,
                              computeWeightsArgs,
                              computeCovariateBalanceArgs = NULL,
                              fitOutcomeModelArgs) {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(CcwAnalysis$new(typedList = args))
}

CcwAnalysis <- R6Class(
  "CcwAnalysis",
  inherit = AbstractSerializableSettings,
  public = list(
    analysisId = NULL,
    description = NULL,
    getDbCloneCensorWeightDataArgs = NULL,
    createStudyPopulationArgs = NULL,
    createClonesArgs = NULL,
    fitCensoringModelsArgs = NULL,
    computeWeightsArgs = NULL,
    computeCovariateBalanceArgs = NULL,
    fitOutcomeModelArgs = NULL,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertInt(self$analysisId, add = errorMessages)
      checkmate::assertString(self$description, add = errorMessages)
      checkmate::assertR6(self$getDbCloneCensorWeightDataArgs, "GetDbCloneCensorWeightDataArgs", add = errorMessages)
      checkmate::assertR6(self$createStudyPopulationArgs, "CreateStudyPopulationArgs", add = errorMessages)
      checkmate::assertR6(self$createClonesArgs, "CreateClonesArgs", add = errorMessages)
      checkmate::assertR6(self$fitCensoringModelsArgs, "FitCensoringModelsArgs", add = errorMessages)
      checkmate::assertR6(self$computeWeightsArgs, "ComputeWeightsArgs", add = errorMessages)
      checkmate::assertR6(self$computeCovariateBalanceArgs, "ComputeCovariateBalanceArgs", null.ok = TRUE, add = errorMessages)
      # CohortMethod has an R6 class with the same name; require members only this package's class has.
      checkmate::assertR6(self$fitOutcomeModelArgs, "FitOutcomeModelArgs",
        public = c("horizons", "timeBasis", "varianceEstimator", "bootstrapArgs"), add = errorMessages
      )
      checkmate::reportAssertions(collection = errorMessages)
    },
    fromList = function(list, requireTyping) {
      super$fromList(list)
      if (requireTyping) {
        self$getDbCloneCensorWeightDataArgs <- GetDbCloneCensorWeightDataArgs$new(untypedList = self$getDbCloneCensorWeightDataArgs)
        self$createStudyPopulationArgs <- .typeCohortMethodSettings(
          self$createStudyPopulationArgs,
          CohortMethod::createCreateStudyPopulationArgs
        )
        self$createClonesArgs <- CreateClonesArgs$new(untypedList = self$createClonesArgs)
        self$fitCensoringModelsArgs <- FitCensoringModelsArgs$new(untypedList = self$fitCensoringModelsArgs)
        self$computeWeightsArgs <- ComputeWeightsArgs$new(untypedList = self$computeWeightsArgs)
        if (!is.null(self$computeCovariateBalanceArgs)) {
          self$computeCovariateBalanceArgs <- .typeCohortMethodSettings(
            self$computeCovariateBalanceArgs,
            CohortMethod::createComputeCovariateBalanceArgs
          )
        }
        self$fitOutcomeModelArgs <- FitOutcomeModelArgs$new(untypedList = self$fitOutcomeModelArgs)
      }
    }
  )
)

#' Create target-treatment-outcomes combinations
#'
#' @details
#' Maps the treatment roles used by the strategies (`"treatment"` and `"comparatorTreatment"`) to
#' cohort IDs, for one eligible (target) cohort and a set of outcomes.
#'
#' @param targetId                    The cohort ID of the eligible (target) cohort. Its start date is
#'                                    time zero.
#' @param treatmentId                 The cohort ID of the treatment with role `"treatment"`.
#' @param comparatorTreatmentId       The cohort ID of the treatment with role `"comparatorTreatment"`,
#'                                    or `NULL`.
#' @param outcomes                    A list of objects of type `Outcome` as created using
#'                                    [CohortMethod::createOutcome()].
#' @param excludedCovariateConceptIds Concept IDs that cannot be used to construct covariates, specific
#'                                    to this combination.
#' @param includedCovariateConceptIds Concept IDs that must be used to construct covariates, specific to
#'                                    this combination.
#'
#' @return
#' An object of type `TargetTreatmentOutcomes`.
#'
#' @examples
#' tto <- createTargetTreatmentOutcomes(
#'   targetId = 1,
#'   treatmentId = 2,
#'   outcomes = list(CohortMethod::createOutcome(outcomeId = 3))
#' )
#'
#' @export
createTargetTreatmentOutcomes <- function(targetId,
                                          treatmentId,
                                          comparatorTreatmentId = NULL,
                                          outcomes,
                                          excludedCovariateConceptIds = c(),
                                          includedCovariateConceptIds = c()) {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(TargetTreatmentOutcomes$new(typedList = args))
}

TargetTreatmentOutcomes <- R6Class(
  "TargetTreatmentOutcomes",
  inherit = AbstractSerializableSettings,
  public = list(
    targetId = NULL,
    treatmentId = NULL,
    comparatorTreatmentId = NULL,
    outcomes = NULL,
    excludedCovariateConceptIds = NULL,
    includedCovariateConceptIds = NULL,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      .assertCovariateId(self$targetId, len = 1, add = errorMessages)
      .assertCovariateId(self$treatmentId, len = 1, add = errorMessages)
      .assertCovariateId(self$comparatorTreatmentId, len = 1, null.ok = TRUE, add = errorMessages)
      checkmate::assertList(self$outcomes, min.len = 1, add = errorMessages)
      for (i in seq_along(self$outcomes)) {
        checkmate::assertR6(self$outcomes[[i]], "Outcome", add = errorMessages)
      }
      checkmate::assertIntegerish(self$excludedCovariateConceptIds, null.ok = TRUE, add = errorMessages)
      checkmate::assertIntegerish(self$includedCovariateConceptIds, null.ok = TRUE, add = errorMessages)
      checkmate::reportAssertions(collection = errorMessages)
      outcomeIds <- sapply(self$outcomes, function(outcome) outcome$outcomeId)
      duplicatedIds <- outcomeIds[duplicated(outcomeIds)]
      if (length(duplicatedIds) > 0) {
        stop(sprintf("Found duplicate outcome IDs: %s", paste(duplicatedIds, collapse = ", ")))
      }
    },
    fromList = function(list, requireTyping) {
      super$fromList(list)
      if (requireTyping) {
        for (i in seq_along(self$outcomes)) {
          outcome <- self$outcomes[[i]]
          # As CohortMethod's Outcome$fromList(): NA is serialized as null.
          if (is.null(outcome$trueEffectSize)) {
            outcome$trueEffectSize <- NA
          }
          self$outcomes[[i]] <- .typeCohortMethodSettings(outcome, CohortMethod::createOutcome)
        }
      }
    }
  )
)

#' Create clone-censor-weight diagnostics thresholds
#'
#' @description
#' Thresholds used to decide whether an analysis passes its diagnostics.
#'
#' @param mdrrThreshold        What is the maximum allowed minimum detectable relative risk (MDRR)?
#' @param easeThreshold        What is the maximum allowed expected absolute systematic error (EASE)?
#' @param sdmThreshold         What is the maximum allowed standardized difference of mean (SDM) of the
#'                             weighted compliant clones at the end of the grace period?
#' @param essFractionThreshold What is the minimum allowed effective sample size of a strategy, as a
#'                             fraction of its number of clones?
#' @param maxWeightThreshold   What is the maximum allowed censoring weight?
#' @param positivityThreshold  What is the minimum allowed probability of following a strategy
#'                             (for example the fitted probability of initiating by the end of the grace
#'                             period)?
#'
#' @return
#' An object of type `CcwDiagnosticThresholds`.
#'
#' @examples
#' thresholds <- createCcwDiagnosticThresholds(maxWeightThreshold = 50)
#'
#' @export
createCcwDiagnosticThresholds <- function(mdrrThreshold = 10,
                                          easeThreshold = 0.25,
                                          sdmThreshold = 0.1,
                                          essFractionThreshold = 0.1,
                                          maxWeightThreshold = 100,
                                          positivityThreshold = 0.01) {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(CcwDiagnosticThresholds$new(typedList = args))
}

CcwDiagnosticThresholds <- R6Class(
  "CcwDiagnosticThresholds",
  inherit = AbstractSerializableSettings,
  public = list(
    mdrrThreshold = NULL,
    easeThreshold = NULL,
    sdmThreshold = NULL,
    essFractionThreshold = NULL,
    maxWeightThreshold = NULL,
    positivityThreshold = NULL,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertNumber(self$mdrrThreshold, lower = 0, add = errorMessages)
      checkmate::assertNumber(self$easeThreshold, lower = 0, add = errorMessages)
      checkmate::assertNumber(self$sdmThreshold, lower = 0, add = errorMessages)
      checkmate::assertNumber(self$essFractionThreshold, lower = 0, upper = 1, add = errorMessages)
      checkmate::assertNumber(self$maxWeightThreshold, lower = 0, add = errorMessages)
      checkmate::assertNumber(self$positivityThreshold, lower = 0, upper = 1, add = errorMessages)
      checkmate::reportAssertions(collection = errorMessages)
    }
  )
)

#' Create full clone-censor-weight analysis specifications
#'
#' @details
#' ## Analyses to Exclude
#'
#' The driver that runs many analyses, `runCcwAnalyses()`, comes in a later version; in this version
#' the specifications can be created, saved and loaded. It will run all combinations of
#' target-treatment-outcomes and analyses. The
#' `analysesToExclude` argument removes some of them. It is a data frame with at least one of the
#' columns `targetId`, `treatmentId`, `comparatorTreatmentId`, `outcomeId`, and `analysisId`. Rows of
#' the reference table that match a row of `analysesToExclude` are skipped.
#'
#' @param ccwAnalysisList                     A list of objects of type `CcwAnalysis` as created using
#'                                            [createCcwAnalysis()].
#' @param targetTreatmentOutcomesList         A list of objects of type `TargetTreatmentOutcomes` as
#'                                            created using [createTargetTreatmentOutcomes()].
#' @param analysesToExclude                   Analyses to exclude. See the Analyses to Exclude section
#'                                            for details.
#' @param refitCensoringModelsForEveryOutcome Should the censoring models be fitted for every outcome
#'                                            (after removing persons with that outcome from the risk
#'                                            sets)? If `FALSE`, the models are shared across outcomes.
#' @param ccwDiagnosticThresholds             An object of type `CcwDiagnosticThresholds` as created
#'                                            using [createCcwDiagnosticThresholds()].
#'
#' @return
#' An object of type `CcwAnalysesSpecifications`.
#'
#' @examples
#' analysis <- createCcwAnalysis(
#'   analysisId = 1,
#'   getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(),
#'   createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(),
#'   createClonesArgs = createCreateClonesArgs(
#'     strategies = list(
#'       createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
#'       createTreatmentStrategy(2, "Never initiate", "neverInitiate")
#'     )
#'   ),
#'   fitCensoringModelsArgs = createFitCensoringModelsArgs(),
#'   computeWeightsArgs = createComputeWeightsArgs(),
#'   fitOutcomeModelArgs = createFitOutcomeModelArgs()
#' )
#' tto <- createTargetTreatmentOutcomes(
#'   targetId = 1,
#'   treatmentId = 2,
#'   outcomes = list(CohortMethod::createOutcome(outcomeId = 3))
#' )
#' specifications <- createCcwAnalysesSpecifications(
#'   ccwAnalysisList = list(analysis),
#'   targetTreatmentOutcomesList = list(tto)
#' )
#'
#' @export
createCcwAnalysesSpecifications <- function(ccwAnalysisList,
                                            targetTreatmentOutcomesList,
                                            analysesToExclude = NULL,
                                            refitCensoringModelsForEveryOutcome = FALSE,
                                            ccwDiagnosticThresholds = createCcwDiagnosticThresholds()) {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(CcwAnalysesSpecifications$new(typedList = args))
}

CcwAnalysesSpecifications <- R6Class(
  "CcwAnalysesSpecifications",
  inherit = AbstractSerializableSettings,
  public = list(
    ccwAnalysisList = NULL,
    targetTreatmentOutcomesList = NULL,
    analysesToExclude = NULL,
    refitCensoringModelsForEveryOutcome = NULL,
    ccwDiagnosticThresholds = NULL,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertList(self$ccwAnalysisList, min.len = 1, add = errorMessages)
      for (i in seq_along(self$ccwAnalysisList)) {
        checkmate::assertR6(self$ccwAnalysisList[[i]], "CcwAnalysis", add = errorMessages)
      }
      checkmate::assertList(self$targetTreatmentOutcomesList, min.len = 1, add = errorMessages)
      for (i in seq_along(self$targetTreatmentOutcomesList)) {
        checkmate::assertR6(self$targetTreatmentOutcomesList[[i]], "TargetTreatmentOutcomes", add = errorMessages)
      }
      checkmate::assertDataFrame(self$analysesToExclude, null.ok = TRUE, add = errorMessages)
      checkmate::assertFlag(self$refitCensoringModelsForEveryOutcome, add = errorMessages)
      checkmate::assertR6(self$ccwDiagnosticThresholds, "CcwDiagnosticThresholds", add = errorMessages)
      checkmate::reportAssertions(collection = errorMessages)
      excludeColumns <- c("targetId", "treatmentId", "comparatorTreatmentId", "outcomeId", "analysisId")
      if (!is.null(self$analysesToExclude) && !any(excludeColumns %in% colnames(self$analysesToExclude))) {
        stop(sprintf("analysesToExclude should have at least one of these columns: %s", paste(excludeColumns, collapse = ", ")))
      }
      analysisIds <- sapply(self$ccwAnalysisList, function(analysis) analysis$analysisId)
      if (any(duplicated(analysisIds))) {
        stop("Duplicate analysis IDs are not allowed")
      }
    },
    fromList = function(list, requireTyping) {
      super$fromList(list)
      if (requireTyping) {
        for (i in seq_along(self$ccwAnalysisList)) {
          self$ccwAnalysisList[[i]] <- CcwAnalysis$new(untypedList = self$ccwAnalysisList[[i]])
        }
        for (i in seq_along(self$targetTreatmentOutcomesList)) {
          self$targetTreatmentOutcomesList[[i]] <- TargetTreatmentOutcomes$new(untypedList = self$targetTreatmentOutcomesList[[i]])
        }
        if (!is.null(self$analysesToExclude)) {
          self$analysesToExclude <- as.data.frame(self$analysesToExclude)
        }
        self$ccwDiagnosticThresholds <- CcwDiagnosticThresholds$new(untypedList = self$ccwDiagnosticThresholds)
      }
    }
  )
)

# Loading, saving, conversion ----------------------------------------------------------------------

.saveSettingsList <- function(settingsList, className, file) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertList(settingsList, min.len = 1, add = errorMessages)
  for (i in seq_along(settingsList)) {
    checkmate::assertR6(settingsList[[i]], className, add = errorMessages)
  }
  checkmate::assertCharacter(file, len = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  for (i in seq_along(settingsList)) {
    settingsList[[i]] <- settingsList[[i]]$toList()
  }
  json <- jsonlite::toJSON(
    settingsList,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    digits = NA
  )
  file <- normalizePath(file, mustWork = FALSE)
  write(json, file)
}

.loadSettingsList <- function(file, generator) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertCharacter(file, len = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  file <- normalizePath(file)
  json <- readChar(file, file.info(file)$size)
  settingsList <- jsonlite::fromJSON(json, simplifyDataFrame = FALSE)
  for (i in seq_along(settingsList)) {
    settingsList[[i]] <- generator$new(untypedList = settingsList[[i]])
  }
  return(settingsList)
}

#' Save a list of CcwAnalysis to file
#'
#' @description
#' Write a list of objects of type `CcwAnalysis` to file. The file is in JSON format.
#'
#' @param ccwAnalysisList A list of objects of type `CcwAnalysis` as created using the
#'                        [createCcwAnalysis()] function.
#' @param file            The name of the file where the results will be written.
#'
#' @return
#' Returns no output.
#'
#' @examples
#' analysis <- createCcwAnalysis(
#'   analysisId = 1,
#'   getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(),
#'   createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(),
#'   createClonesArgs = createCreateClonesArgs(
#'     strategies = list(
#'       createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
#'       createTreatmentStrategy(2, "Never initiate", "neverInitiate")
#'     )
#'   ),
#'   fitCensoringModelsArgs = createFitCensoringModelsArgs(),
#'   computeWeightsArgs = createComputeWeightsArgs(),
#'   fitOutcomeModelArgs = createFitOutcomeModelArgs()
#' )
#' file <- tempfile(fileext = ".json")
#' saveCcwAnalysisList(list(analysis), file)
#' analysisList <- loadCcwAnalysisList(file)
#' unlink(file)
#'
#' @export
saveCcwAnalysisList <- function(ccwAnalysisList, file) {
  .saveSettingsList(ccwAnalysisList, "CcwAnalysis", file)
  invisible(NULL)
}

#' Load a list of CcwAnalysis from file
#'
#' @description
#' Load a list of objects of type `CcwAnalysis` from file. The file is in JSON format.
#'
#' @param file The name of the file.
#'
#' @return
#' A list of objects of type `CcwAnalysis`.
#'
#' @inherit saveCcwAnalysisList examples
#'
#' @export
loadCcwAnalysisList <- function(file) {
  return(.loadSettingsList(file, CcwAnalysis))
}

#' Save a list of TargetTreatmentOutcomes to file
#'
#' @description
#' Write a list of objects of type `TargetTreatmentOutcomes` to file. The file is in JSON format.
#'
#' @param targetTreatmentOutcomesList A list of objects of type `TargetTreatmentOutcomes` as created
#'                                    using the [createTargetTreatmentOutcomes()] function.
#' @param file                        The name of the file where the results will be written.
#'
#' @return
#' Returns no output.
#'
#' @examples
#' tto <- createTargetTreatmentOutcomes(
#'   targetId = 1,
#'   treatmentId = 2,
#'   outcomes = list(CohortMethod::createOutcome(outcomeId = 3))
#' )
#' file <- tempfile(fileext = ".json")
#' saveTargetTreatmentOutcomesList(list(tto), file)
#' ttoList <- loadTargetTreatmentOutcomesList(file)
#' unlink(file)
#'
#' @export
saveTargetTreatmentOutcomesList <- function(targetTreatmentOutcomesList, file) {
  .saveSettingsList(targetTreatmentOutcomesList, "TargetTreatmentOutcomes", file)
  invisible(NULL)
}

#' Load a list of TargetTreatmentOutcomes from file
#'
#' @description
#' Load a list of objects of type `TargetTreatmentOutcomes` from file. The file is in JSON format.
#'
#' @param file The name of the file.
#'
#' @return
#' A list of objects of type `TargetTreatmentOutcomes`.
#'
#' @inherit saveTargetTreatmentOutcomesList examples
#'
#' @export
loadTargetTreatmentOutcomesList <- function(file) {
  return(.loadSettingsList(file, TargetTreatmentOutcomes))
}

#' Convert untyped list to CcwAnalysesSpecifications
#'
#' @param untypedList A list of untyped objects. For example, these could be objects from a call
#'                    to `jsonlite::fromJSON()`. Importantly, `simplifyDataFrame` must be set to
#'                    `FALSE` when doing so.
#'
#' @return
#' An object of type `CcwAnalysesSpecifications`.
#'
#' @examples
#' analysis <- createCcwAnalysis(
#'   analysisId = 1,
#'   getDbCloneCensorWeightDataArgs = createGetDbCloneCensorWeightDataArgs(),
#'   createStudyPopulationArgs = CohortMethod::createCreateStudyPopulationArgs(),
#'   createClonesArgs = createCreateClonesArgs(
#'     strategies = list(
#'       createTreatmentStrategy(1, "Initiate within 7 days", "initiateWithin", graceEnd = 7),
#'       createTreatmentStrategy(2, "Never initiate", "neverInitiate")
#'     )
#'   ),
#'   fitCensoringModelsArgs = createFitCensoringModelsArgs(),
#'   computeWeightsArgs = createComputeWeightsArgs(),
#'   fitOutcomeModelArgs = createFitOutcomeModelArgs()
#' )
#' specifications <- createCcwAnalysesSpecifications(
#'   ccwAnalysisList = list(analysis),
#'   targetTreatmentOutcomesList = list(createTargetTreatmentOutcomes(
#'     targetId = 1,
#'     treatmentId = 2,
#'     outcomes = list(CohortMethod::createOutcome(outcomeId = 3))
#'   ))
#' )
#' untypedList <- jsonlite::fromJSON(specifications$toJson(), simplifyDataFrame = FALSE)
#' specifications2 <- convertUntypedListToCcwAnalysesSpecifications(untypedList)
#'
#' @export
convertUntypedListToCcwAnalysesSpecifications <- function(untypedList) {
  return(CcwAnalysesSpecifications$new(untypedList = untypedList))
}

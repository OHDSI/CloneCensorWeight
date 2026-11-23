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

# The base class of the settings objects is CohortMethod's AbstractSerializableSettings, which
# provides toList() and fromList(). This package's version also says which fields a list or JSON saved
# by an earlier version did not have, and writes JSON with 15 significant digits.

AbstractSerializableSettings <- R6Class(
  "AbstractSerializableSettings",
  inherit = CohortMethod::AbstractSerializableSettings,
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
      if (requireTyping) {
        .messageFilledFields(self, names(args))
      }
      self$validate()
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
    }
  )
)

# Settings loaded from a list or JSON without some fields (saved by an earlier version) take the current
# class-level defaults of those fields, which can differ from the behaviour of the version that saved
# them, so say which fields were filled.
.messageFilledFields <- function(settings, savedNames) {
  fields <- Filter(function(name) {
    value <- settings[[name]]
    !is.function(value) && !is(value, "environment") && !is.null(value)
  }, setdiff(ls(settings), c(".__enclos_env__", savedNames)))
  if (length(fields) > 0) {
    values <- vapply(fields, function(name) paste(format(settings[[name]]), collapse = ", "), character(1))
    message(sprintf(
      "%s saved without %s: filled with the current default%s %s, which can differ from the behaviour of the version that saved them (see NEWS)",
      class(settings)[1],
      paste(fields, collapse = ", "),
      if (length(fields) > 1) "s" else "",
      paste(sprintf("%s = %s", fields, values), collapse = ", ")
    ))
  }
}

# .dateCheck() is copied from CohortMethod (R/SettingsObjects.R, Copyright Observational Health Data
# Sciences and Informatics, Apache 2.0), because CohortMethod does not export it.

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
#' - `"initiateWithin"`: start the treatment on a day in `graceStart..graceEnd` (days since index). With
#'   `graceStart > 0`, a start before day `graceStart` is a deviation.
#' - `"notInitiateWithin"`: do not start the treatment on or before day `graceEnd`.
#' - `"neverInitiate"`: never start the treatment during follow-up.
#'
#' `graceStart` applies to `"initiateWithin"` only and must be 0 for the other types. This version does
#' not support `prohibitedTreatments` when creating clones.
#'
#' @param strategyId           An integer identifying the strategy. Used for labels and serialization;
#'                             within an analysis the strategies are numbered by their position.
#' @param label                A short label for the strategy, used in output.
#' @param type                 The strategy type: `"initiateWithin"`, `"notInitiateWithin"`, or
#'                             `"neverInitiate"`.
#' @param treatment            The role of the treatment the strategy refers to: `"treatment"` or
#'                             `"comparatorTreatment"`.
#' @param graceStart           First day (since index) of the grace period of an `"initiateWithin"`
#'                             strategy (0 for the other types).
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
      if (self$type != "initiateWithin" && self$graceStart != 0) {
        stop(sprintf("graceStart must be 0 for strategy type '%s': it applies to 'initiateWithin' strategies only", self$type))
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
#'                                          diagnoses recorded on the index date. [fitCensoringModels()]
#'                                          then warns, because its models predict initiation on day 0;
#'                                          use `endDays = -1` to end the windows the day before. Exclude
#'                                          the treatment from the covariates: its concepts with their
#'                                          descendants (`excludedCovariateConceptIds` with
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
      asList$covariateSettings <- CohortMethod::.convertAttrToMember(self$covariateSettings)
      return(asList)
    },
    fromList = function(list, requireTyping) {
      super$fromList(list)
      if (requireTyping) {
        self$covariateSettings <- CohortMethod::.convertEmptyListToLogical(CohortMethod::.convertMemberToAttr(list$covariateSettings))
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
#' period is used by default (`initiateWithinModel = "binary"`), or the initiation hazard model
#' (`initiateWithinModel = "hazard"`). For `"notInitiateWithin"` and `"neverInitiate"` strategies the
#' hazard of initiation is modelled with `hazardModelType`.
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
#' `pooledLogisticEngine` chooses how the pooled logistic hazard model is fit. `"lr"`
#' expands the risk set to one row per person and interval, with the person's covariates, the time
#' basis and the interactions on every row. `"plr"` (the default) gives Cyclops's pooled logistic regression
#' (`modelType = "plr"`) the covariates once per person, each person's intervals (as a first and last
#' interval, or as a list when the indicator basis leaves out an interval), the time basis per interval
#' and the interactions as covariate-by-basis products, so its memory stays near that of one row per
#' person. A `"day0"` interaction column that is not a basis column enters with its main effect fixed
#' at 0. With a fixed prior variance (or no prior, when the maximum likelihood estimate exists) the two
#' engines agree to about 1e-6 in the
#' coefficients and weights. With cross-validation, the plr of Cyclops develop (from `2b211c6a`)
#' flags separation in the time-basis and interaction columns as `"lr"` does,
#' and starts its data-driven `startingVariance` (the default `-1`) from the person-interval rows, so
#' on a test case it selected the variance `"lr"` selects; with a `"day0"` interaction its start
#' differs slightly from `"lr"`'s (about 1e-4 relative), because the fixed main effect counts as a
#' column. Earlier plr builds could pick a much smaller variance (for example with `fold = 5`), and
#' before Cyclops' `fitter = "super"` covered `modelType = "plr"`, the plr hazard model's search
#' under the default control could select a different variance than `"lr"`'s. `"plr"` needs much
#' less memory and is faster (10^6 persons, G = 30, fixed prior: 235 s and 14 GB against 353 s and
#' 23 GB).
#'
#' `maxCohortSizeForFitting` applies to every censoring model: the binary models at the end of the
#' grace period sample their starters and non-starters (through [CohortMethod::createPs()]), and a
#' hazard model whose risk set holds more persons is fit on a uniform sample of that many persons and
#' then predicts for every person (see [fitCensoringModels()]).
#'
#' **Weights of `"notInitiateWithin"` after the grace period.** With `postGraceModel = "hazard"` (the
#' default), the weight of a `"notInitiateWithin(G)"` clone after day `G` is `1 / S(G | X)`, where
#' `S(G | X)` is the probability of not initiating on days `0..G` from the initiation hazard model, the
#' same model that gives the weights on days up to `G`. With `postGraceModel = "binary"`, it is
#' `1 / (1 - pi(X))`, where `pi(X)` is the probability of starting by day `G` from the binary model at
#' `G`, the model the `"initiateWithin(G)"` strategy uses (the fit is shared when the analysis has such
#' a strategy for the same treatment, and made otherwise, with cross-validation when the prior uses
#' it). The initiation hazard model gives the weights on days up to `G`, and all weights of
#' `"neverInitiate"` strategies, under both options.
#'
#' When each option is valid:
#'
#' - **`"binary"`** needs two assumptions. (1) *No effect of grace-period initiation on grace-period
#'   risk.* The binary model is fitted on the persons still under follow-up at `G`, so `1 - pi(X)` is
#'   the probability of not starting by `G` among persons with no outcome or death by `G`. It equals
#'   `S(G | X)` only when starting the treatment during the grace period does not change the risk of
#'   the outcome or death during the grace period. When starting lowers that risk (a treatment that
#'   acts within days), the persons who started are over-represented among those still under
#'   follow-up at `G`: `pi(X)` is too large, the weights `1 / (1 - pi(X))` are too large in a
#'   covariate-dependent way, and the risks of the `"notInitiateWithin"` strategy are biased (in a
#'   simulation in which starting cut the outcome hazard by 80% from the next day, the risk at day 60
#'   was biased by +0.020, z = 12, coverage 0.72, while `"hazard"` was unbiased).
#'   [computeGracePeriodSelection()] checks this assumption; [computeCensoringModelCalibration()]
#'   cannot, because the binary model reproduces the observed share among the persons under follow-up
#'   at `G` by construction. (2) *A correctly specified binary model:* a logistic regression of
#'   starting by `G` on the main effects of the covariates. With continuous covariates that change
#'   *when* persons start, the probability of starting by `G` is not logit-linear in them.
#' - **`"hazard"`** follows each person only until the outcome or death, so it needs neither
#'   assumption (1) nor a model of the selected survivors. It needs a correctly specified hazard model.
#'   Without interactions (`covariateTimeInteraction = "none"`) the covariate effects on initiation are
#'   the same on every day of the grace period; when covariates change *when* persons start (for
#'   example who starts on day 0, the day of the index visit), `1 - S(G | X)` is then miscalibrated.
#'   With `covariateTimeInteraction = "day0"` (the default) the covariate effects on day 0 differ from
#'   those on later days, which handles covariates that change start timing on day 0; `"interval"`
#'   lets them differ in every interval of the time basis.
#'
#' Compare the options on your data with [computeCcwDiagnostics()] (calibration, the grace-period
#' selection check, weights, balance and positivity per option) and [computeCcwSensitivity()] (the
#' estimates per option).
#'
#' **Weights of `"initiateWithin"` after the grace period.** With `initiateWithinModel = "binary"` (the
#' default), the weight of an `"initiateWithin(G)"` clone after day `G` is `1 / pi(X)`, where `pi(X)` is
#' the probability of starting by day `G` from the binary model at `G`. With
#' `initiateWithinModel = "hazard"`, it is `1 / (1 - S(G | X))`, where `S(G | X)` is the probability of
#' not initiating on days `0..G` from the initiation hazard model through `G` (with
#' `hazardModelType`, `timeBasis` and `covariateTimeInteraction`), and the stabilized weight multiplies
#' by `1 - S-bar(G)` of its numerator model. That hazard model is the one of a `"notInitiateWithin(G)"`
#' strategy of the same treatment (the fit is shared when the analysis has one). When each option is
#' valid:
#'
#' - **`"binary"`** is fitted on the persons under follow-up at `G`, whose probability of starting by
#'   `G` is the one the weight needs, so it is valid whether or not starting during the grace period
#'   changes the risk during the grace period. It needs a correctly specified binary model: a logistic
#'   regression of starting by `G` on the main effects of the covariates. With covariates that change
#'   *when* persons start (for example who starts on day 0), the probability of starting by `G` is not
#'   logit-linear in them, and the weights are biased (in simulations with five covariates, and with
#'   covariates that change who starts on day 0, the risk of the `"initiateWithin"` strategy at day 60
#'   was biased by -0.26 and +0.32 empirical standard errors).
#' - **`"hazard"`** models the timing of initiation, with covariate effects that can differ on day 0
#'   (`covariateTimeInteraction = "day0"`) or in every interval (`"interval"`). But `1 - S(G | X)` is
#'   the probability of starting by `G` of all persons at day 0, followed until the outcome or death,
#'   not of the persons still under follow-up at `G`. The two agree only when starting the treatment
#'   during the grace period does not change the risk of the outcome or death during the grace period:
#'   the same assumption (1) as `postGraceModel = "binary"` above, the other way round. When starting
#'   lowers that risk, the persons who started are over-represented among those under follow-up at `G`,
#'   so `1 - S(G | X)` is too small for them and the weights after `G` too large. It also needs a
#'   correctly specified hazard model.
#'
#' [computeGracePeriodSelection()] checks the assumption of `"hazard"`; a flag there means that
#' `"binary"` is the valid choice. [computeCensoringModelCalibration()] compares `pi(X)` and
#' `1 - S(G | X)` with the observed share who started by `G` (with the same caveat: under grace-period
#' selection the observed share is that of the selected survivors).
#'
#' **The `"initiateWithin"` regime (`graceRegime`).** "Start within the grace period" does not say on
#' which day of it a person starts, and two regimes are in use (plan section 2.2). They agree up to
#' `G` and differ after it when the start day within the grace period changes later risk.
#'
#' - **`"binaryAtG"`** (the default): the clones still followed after `G` (those that started by `G`)
#'   are weighted by `1 / pi(X)`, the inverse probability of starting by `G` of the persons under
#'   follow-up at `G`. The weight is the same for a start on day 0 and on day `G`, so each person still
#'   followed at `G` is represented by the starters of the same covariates, pooled over their start
#'   days: after `G`, the regime follows the risk of a start on a day drawn from the starters' start
#'   days. This is standard clone-censor-weight practice and gives stable weights.
#' - **`"naturalGrace"`** (Cain et al. 2010; Wanis et al. 2024): persons who would start during the grace
#'   period start on their own day, and persons still untreated on day `G` start on day `G`. Its
#'   weights are the last-interval weights: 1 for the clones that started before day `G`, and
#'   `1 / P(start on day G | not started before G, X)` for those that start on day `G`. With
#'   `initiateWithinModel = "binary"` that probability comes from a binary model at `G` fit by
#'   [CohortMethod::createPs()] on the persons who reach day `G` untreated and under follow-up, with the
#'   response "started on day `G`"; with `"hazard"` it is the initiation hazard of day `G`,
#'   `1 - S(G | X) / S(G - 1 | X)`, whose risk set on day `G` is the same persons (so the grace-period
#'   selection that biases `initiateWithinModel = "hazard"` under `"binaryAtG"` does not arise). With
#'   the pooled logistic hazard model the last interval must be day `G` alone (`graceEnd` a multiple of
#'   `intervalLength`). The weights of few clones (those that start on day `G`) vary, so they are less
#'   stable than those of `"binaryAtG"`. Only those clones stand for the persons who reach day `G`
#'   untreated, so the regime needs starts on day `G`: [fitCensoringModels()] warns when nobody starts
#'   on day `G`, and [computePositivity()] shows covariate strata where few do. When few day-`G`
#'   starters carry the weight after `G`, the standard errors of risks at horizons shortly after `G`
#'   fall with the estimate and their intervals can under-cover (simulation of plan section 9 item 7,
#'   `G = 7`, 6000 persons with about 140 day-7 starters and 14 of their outcomes by day 14: coverage
#'   0.89 to 0.90 at day 14, 0.95 to 0.96 at day 60; with 24000 persons, 0.94 to 0.95 at day 14).
#'   The stabilized weights after `G` multiply by the same quantity for every clone (the marginal
#'   probability of starting by `G` of the persons under follow-up at `G`), so the weighted
#'   Kaplan-Meier risks do not depend on stabilization.
#'
#' [computeTrueRisks()] gives the truth of either regime (its argument `graceRegime`). The choice is a
#' choice of estimand, not of model: [computeCcwSensitivity()] reports both.
#'
#' **`"initiateWithin"` strategies with `graceStart = S > 0`.** A start before day `S` is a deviation,
#' so the weights of days up to `S` are those of a `"notInitiateWithin(S - 1)"` strategy, `1 / S(d - 1 |
#' X)` from the initiation hazard model through `S - 1` (with `initiateWithinModel = "binary"`; the
#' fit is shared with such a strategy in the analysis) or through `G` (with `"hazard"`). The binary model
#' at `G` (and the hazard model's probability, `1 - S(G | X) / S(S - 1 | X)`) is then that of starting
#' on days `S..G` of the persons who have not started before day `S`, and multiplies `S(S - 1 | X)`
#' after `G`. With the pooled logistic hazard model and `initiateWithinModel = "hazard"`, `S` must be a
#' multiple of `intervalLength`.
#'
#' **Covariate-by-time interactions.** With `covariateTimeInteraction = "day0"` (the default), the
#' pooled logistic initiation hazard model also has, for each covariate, its product with an indicator
#' of the first interval of the time basis. With daily intervals (`timeBasis` with
#' `intervalLength = 1`, the default) that is day 0, so the covariates can change initiation on day 0
#' differently from later days. With longer intervals the model has one row per person and interval,
#' so it cannot tell day 0 from the other days of the first interval: `"day0"` is then the product
#' with the indicator of the first interval (days `0..intervalLength - 1`). Use daily intervals when
#' initiation on day 0 itself depends on the covariates differently. With `"interval"`, the model has
#' each covariate's product with each time-basis column (with the indicator basis: with the indicator
#' of each interval except the reference), so the covariate effects can differ in every interval.
#' With `"none"` it has no interactions. The interaction columns are penalized by the prior like the
#' covariates; the time-basis columns are never penalized. The Cox model has no interactions, so with
#' `hazardModelType = "cox"` the default is `"none"`, and `"day0"` and `"interval"` are not available.
#'
#' Settings saved before `covariateTimeInteraction`, `postGraceModel`, `initiateWithinModel` and
#' `graceRegime` existed load with `"day0"` (`"none"` for the Cox model), `"hazard"`, `"binary"` and
#' `"binaryAtG"`.
#'
#' @param hazardModelType          The hazard model: `"logistic"` (pooled logistic, the default) or
#'                                 `"cox"` (see details for its ties caveat).
#' @param timeBasis                An object of type `TimeBasis` as created using
#'                                 [createTimeBasis()], used by the pooled logistic model. The default
#'                                 is a daily indicator basis.
#' @param covariateTimeInteraction Covariate-by-time interactions in the pooled logistic initiation
#'                                 hazard model: `"day0"` (each covariate times an indicator of day 0,
#'                                 or of the first interval with longer intervals), `"interval"` (each
#'                                 covariate times each time-basis column), or `"none"`. See details.
#'                                 Must be `"none"` with `hazardModelType = "cox"`, and is `"none"` by
#'                                 default then.
#' @param pooledLogisticEngine     How the pooled logistic initiation hazard model is fit: `"lr"`
#'                                 (logistic regression on the rows expanded to one per person and
#'                                 interval, with their covariates) or `"plr"` (Cyclops's pooled
#'                                 logistic regression, which keeps the covariates one row per person
#'                                 and builds the person-interval rows, time basis and interactions on
#'                                 the fly, so it needs much less memory). Both fit the same model; see
#'                                 details for when they agree. `"plr"`, the default, needs a Cyclops with
#'                                 pooled logistic regression and stops with an error otherwise. The
#'                                 model's `diagnostics$engine` says which was used.
#' @param postGraceModel           The model that gives the weights of `"notInitiateWithin"` strategies
#'                                 after the grace period: `"hazard"` (`1 / S(G | X)` from the
#'                                 initiation hazard model) or `"binary"` (`1 / (1 - pi(X))` from the
#'                                 binary model at the end of the grace period). See details for when
#'                                 each is valid.
#' @param initiateWithinModel      The model that gives the weights of `"initiateWithin"` strategies
#'                                 after the grace period: `"binary"` (`1 / pi(X)` from the binary
#'                                 model at the end of the grace period) or `"hazard"`
#'                                 (`1 / (1 - S(G | X))` from the initiation hazard model). See details
#'                                 for when each is valid.
#' @param graceRegime              The regime of `"initiateWithin"` strategies, which their weights after
#'                                 the grace period identify: `"binaryAtG"` (`1 / pi(X)` for every clone
#'                                 that started by `G`) or `"naturalGrace"` (the last-interval weight:
#'                                 1 for clones that started before `G`, the inverse probability of
#'                                 starting on day `G` for those that start on day `G`). See details.
#' @param censoringMechanisms      The censoring mechanisms to weight for. This version supports only
#'                                 `"artificial"` (deviation from a strategy). Censoring at the end of
#'                                 observation is not weighted, so it is assumed to be independent of
#'                                 the outcome. Risks are biased in data where the end of observation
#'                                 depends on the outcome, for example when persons who die are
#'                                 observed until death but survivors only until their last record.
#' @param useCovariates            Use the baseline covariates in the censoring models?
#' @param excludeCovariateIds      Exclude these covariates from the censoring models.
#' @param includeCovariateIds      Include only these covariates in the censoring models.
#' @param maxCohortSizeForFitting  The largest number of persons (for the binary model: of starters
#'                                 and of non-starters each) a censoring model is fit on. The
#'                                 binary model at the end of the grace period samples its starters
#'                                 and its non-starters to at most this many each (through
#'                                 [CohortMethod::createPs()]). A hazard model (pooled logistic or Cox)
#'                                 whose risk set holds more persons is fit on a uniform sample of this
#'                                 many persons, with all their rows, drawn with the seed of `control`;
#'                                 a uniform sample needs no correction of the intercept, and the
#'                                 linear predictor is computed for every person. `0` means no
#'                                 sampling.
#' @param errorOnHighCorrelation   Check each covariate for high correlation with the response, and
#'                                 throw an error if found?
#' @param stopOnError              If an error occurs while fitting, should the function stop?
#' @param prior                    The prior used to fit the models. See [Cyclops::createPrior()].
#' @param control                  The control object used to fit the models. See
#'                                 [Cyclops::createControl()]. For the Cox initiation hazard model,
#'                                 its `newtonMemoryBudget` is raised by 12 bytes per covariate value of
#'                                 the risk set (about 2.7 GB at 10^6 persons with ~240 values each), the
#'                                 row-wise copy of the design that Cyclops' prox-Newton steps make, so
#'                                 that large Laplace fits are not left on coordinate descent. The budget
#'                                 holds per fit: cross-validation folds on several threads may each use it.
#'                                 The default `resetCoefficients = FALSE` starts each cross-validation fold
#'                                 fit from that fold's solution at the nearest prior variance already tried
#'                                 (Cyclops' `cvWarmStart`): the same fits as cold starts to rounding, and
#'                                 independent of `threads`. The default `fitter = "super"` (an interim
#'                                 choice, to be revisited by the stopping-rule study of the plan) changes
#'                                 the variance search of the logistic models (both engines): it steps down from
#'                                 `startingVariance` first, never below the variance at which every
#'                                 penalized coefficient is 0, and refines its result. It selects a smaller
#'                                 variance than `fitter = "default"` and is faster. The fits still stop by
#'                                 Cyclops' `"gap"` criterion under its default `algorithm = "auto"`. Cox models
#'                                 are not affected. No `tolerance` is set, since under `"super"` a logistic
#'                                 fit ignores (and warns about) a tolerance given without a `convergenceType`.
#'
#' @return
#' An object of type `FitCensoringModelsArgs`.
#'
#' @examples
#' args <- createFitCensoringModelsArgs(hazardModelType = "cox")
#' args <- createFitCensoringModelsArgs(covariateTimeInteraction = "none", postGraceModel = "binary")
#' args <- createFitCensoringModelsArgs(initiateWithinModel = "hazard")
#' args <- createFitCensoringModelsArgs(graceRegime = "naturalGrace")
#'
#' @export
createFitCensoringModelsArgs <- function(hazardModelType = "logistic",
                                         timeBasis = createTimeBasis(),
                                         covariateTimeInteraction = if (identical(hazardModelType, "cox")) "none" else "day0",
                                         pooledLogisticEngine = "plr",
                                         postGraceModel = "hazard",
                                         initiateWithinModel = "binary",
                                         graceRegime = "binaryAtG",
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
                                                                          resetCoefficients = FALSE,
                                                                          fitter = "super",
                                                                          cvRepetitions = 3,
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
    # Defaults, so settings saved before these fields existed load (with the engine they were fit with):
    covariateTimeInteraction = "day0",
    pooledLogisticEngine = "lr",
    postGraceModel = "hazard",
    initiateWithinModel = "binary",
    graceRegime = "binaryAtG",
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
      checkmate::assertChoice(self$pooledLogisticEngine, c("plr", "lr"), add = errorMessages)
      checkmate::assertChoice(self$postGraceModel, c("binary", "hazard"), add = errorMessages)
      checkmate::assertChoice(self$initiateWithinModel, c("binary", "hazard"), add = errorMessages)
      checkmate::assertChoice(self$graceRegime, c("binaryAtG", "naturalGrace"), add = errorMessages)
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
      # The Cox model has no interactions, so its settings saved without the field load with "none":
      if (is.null(list$covariateTimeInteraction) && identical(self$hazardModelType, "cox")) {
        self$covariateTimeInteraction <- "none"
      }
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
#' The rule is applied on the scale of the weights of [computeCensoringWeights()], stabilized or not as
#' set by `stabilized` in [createComputeWeightsArgs()]. With `truncateMsmWeights = TRUE` in
#' [createComputeWeightsArgs()] (the default), the same rule (type, `upper` and basis) is also applied,
#' separately, to the weights on the other scale, on their own scale: the percentile is then taken
#' over the clones' maximum weights on that scale. So the weights of the Cox and pooled logistic
#' outcome models (`msmWeights` in [createFitOutcomeModelArgs()]) are bounded by a cap of their own
#' scale whichever scale they are on. With `truncateMsmWeights = FALSE`, the weights on the other scale
#' are formed from the truncated weights and are not bounded by any cap.
#'
#' A percentile cap adapts to each scale, but a fixed cap does not: unstabilized weights are larger than
#' stabilized ones (by the inverse of the stabilization numerator, the share of persons still following
#' the strategy), so a cap meant for stabilized weights truncates far more of the unstabilized weights,
#' which are the default weights of the Cox and pooled logistic outcome models. For `"fixed"`,
#' `otherScaleUpper` sets the cap on the other scale; by default (`NULL`) it is `upper`, the same number
#' on both scales, and [computeCensoringWeights()] warns when it truncates a much larger share of the
#' clones on the other scale than on the scale of `stabilized`.
#'
#' @param type  The truncation rule: `"none"`, `"percentile"`, or `"fixed"`.
#' @param upper For `"percentile"`, the percentile in (0, 1). For `"fixed"`, the cap (greater than 1).
#'              Ignored for `"none"`.
#' @param basis For `"percentile"`, the clones over which the percentile is taken: `"varying"` or
#'              `"allClones"` (see details). Ignored for other types.
#' @param otherScaleUpper For `"fixed"`, the cap (greater than 1) on the other scale (unstabilized when
#'              `stabilized = TRUE` in [createComputeWeightsArgs()]), used with `truncateMsmWeights =
#'              TRUE`. If `NULL`, `upper` (see details). Ignored for other types.
#'
#' @return
#' An object of type `WeightTruncation`.
#'
#' @examples
#' truncation <- createWeightTruncation(type = "fixed", upper = 20)
#' truncation <- createWeightTruncation(type = "fixed", upper = 5, otherScaleUpper = 50)
#' truncation <- createWeightTruncation(type = "percentile", upper = 0.99, basis = "allClones")
#'
#' @export
createWeightTruncation <- function(type = "percentile",
                                   upper = 0.99,
                                   basis = "varying",
                                   otherScaleUpper = NULL) {
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
    otherScaleUpper = NULL,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertChoice(self$type, c("none", "percentile", "fixed"), add = errorMessages)
      checkmate::assertNumber(self$upper, null.ok = TRUE, add = errorMessages)
      checkmate::assertNumber(self$otherScaleUpper, lower = 1, finite = TRUE, null.ok = TRUE, add = errorMessages)
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
#' @details
#' `stabilized` sets the weights `weight` that [computeCensoringWeights()] stores on the segments. These
#' weights give the weighted Kaplan-Meier and Aalen-Johansen risks (which, without truncation, are the
#' same either way) and the weighted counts of [fitOutcomeModel()]. The Cox and pooled logistic outcome
#' models do not use `stabilized`: they use the weights chosen by `msmWeights` in
#' [createFitOutcomeModelArgs()] (unstabilized by default). [computeCensoringWeights()] therefore also
#' stores the weight on the other scale (`otherScaleWeight`: unstabilized when `stabilized = TRUE`,
#' stabilized otherwise).
#'
#' **Truncation scales.** The truncation rule (`truncation`) is applied to `weight`, with a cap computed
#' on its scale. `truncateMsmWeights` sets how the weights on the other scale are truncated:
#'
#' - `TRUE` (the default): the same rule (type, `upper` and basis) is applied separately to the weights
#'   on the other scale, with a cap computed on that scale from the untruncated weights on that scale (a
#'   fixed cap is the same number on both scales unless `otherScaleUpper` is set in
#'   [createWeightTruncation()], see there).
#'   So the weights of the Cox and pooled logistic models are capped on their own scale whatever
#'   `stabilized` and `msmWeights` are. The two scales can truncate different segments.
#' - `FALSE`: the weights on the other scale are formed from the truncated `weight` (divided or
#'   multiplied by the segment's stabilization numerator), so the same segments are truncated on both
#'   scales, but the cap does not bound the weights on the other scale: their cap on a day is the cap
#'   divided (or multiplied) by that day's numerator. For a `"neverInitiate"` strategy the numerator is
#'   the share of persons who have not initiated yet, which can be small at long horizons: in a
#'   simulation with a 365-day horizon and the defaults otherwise, the 99th percentile cap was 5.7 and
#'   the largest unstabilized weight of the `"neverInitiate"` strategy was 369. This was the behaviour
#'   before `truncateMsmWeights` existed; settings saved without it load with `TRUE`.
#'
#' The argument sits here rather than in [createFitOutcomeModelArgs()] because the cap of a percentile
#' rule is computed from the untruncated weights of all clones (and, with `basis = "varying"`, from
#' which clones are under follow-up after the grace periods), which only [computeCensoringWeights()]
#' has: the outcome models receive weights that are already truncated. [computeCensoringWeights()]
#' reports the caps and the shares truncated on both scales, and [fitOutcomeModel()] reports the cap,
#' the largest weight and the share truncated of the weights its Cox or pooled logistic model used.
#' [computeWeightSummary()] summarizes both scales before and after truncation.
#'
#' @param stabilized         Stabilize the weights `weight` with a numerator model of time only? It
#'                           does not change which weights the Cox and pooled logistic outcome models
#'                           use (see details).
#' @param truncation         An object of type `WeightTruncation` as created using
#'                           [createWeightTruncation()].
#' @param truncateMsmWeights Also apply the truncation rule to the weights on the other scale (the
#'                           unstabilized weights when `stabilized = TRUE`), on their own scale? If
#'                           `FALSE`, they are formed from the truncated weights (see details).
#'
#' @return
#' An object of type `ComputeWeightsArgs`.
#'
#' @examples
#' args <- createComputeWeightsArgs(stabilized = FALSE)
#' args <- createComputeWeightsArgs(truncateMsmWeights = FALSE)
#'
#' @export
createComputeWeightsArgs <- function(stabilized = TRUE,
                                     truncation = createWeightTruncation(),
                                     truncateMsmWeights = TRUE) {
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
    # Default, so settings saved before this field existed load:
    truncateMsmWeights = TRUE,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertFlag(self$stabilized, add = errorMessages)
      checkmate::assertR6(self$truncation, "WeightTruncation", add = errorMessages)
      checkmate::assertFlag(self$truncateMsmWeights, add = errorMessages)
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
#' With `method = "reweight"` (the default), a replicate gives each person a frequency weight, the
#' number of times it was drawn, and refits the censoring models of the full data with these weights
#' (reusing their Cyclops data where it can), then the outcome model with each clone's weight times its
#' person's frequency weight. With `method = "copy"`, a replicate copies each drawn person as many times
#' as drawn into new data and repeats the whole analysis on it. With the same counts and
#' `fixHyperparameters = TRUE`, the two give the same replicates (see [bootstrapCcwAnalysis()] for
#' when they differ).
#'
#' @param replicates         The number of bootstrap replicates.
#' @param seed               The random seed.
#' @param ci                 How to compute confidence intervals: `"percentile"` or `"normal"`
#'                           (bias-corrected), as the `"perc"` and `"norm"` types of [boot::boot.ci()].
#' @param fixHyperparameters Reuse the prior variances selected by cross-validation on the full data
#'                           in every replicate?
#' @param method             How a replicate is analyzed: `"reweight"` (frequency weights on the full
#'                           data) or `"copy"` (copies of the drawn persons).
#' @param warmStart          With `method = "reweight"`: start the refits of the censoring models from
#'                           the coefficients of the full-data models? With `FALSE`, they start as new
#'                           fits do.
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
                                fixHyperparameters = TRUE,
                                method = "reweight",
                                warmStart = TRUE) {
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
    # Defaults, so settings saved before these fields existed load:
    method = "reweight",
    warmStart = TRUE,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertInt(self$replicates, lower = 1, add = errorMessages)
      checkmate::assertInt(self$seed, add = errorMessages)
      checkmate::assertChoice(self$ci, c("percentile", "normal"), add = errorMessages)
      checkmate::assertFlag(self$fixHyperparameters, add = errorMessages)
      checkmate::assertChoice(self$method, c("reweight", "copy"), add = errorMessages)
      checkmate::assertFlag(self$warmStart, add = errorMessages)
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
#' `"fineGray"` (an option, not the default) needs clones created with `competingEvents = "death"`. It
#' reports the subdistribution hazard ratio (`"shr"`) of a weighted Fine-Gray model: a weighted Cox
#' model of the outcome in which a clone that dies stays in the risk set after its death, with its weight
#' at death times G(t-) / G(T-), where T is the time of death and G the censoring survival of its
#' strategy (`fineGrayCensoring`). Like `"cox"`, it reports the Aalen-Johansen risks, risk differences
#' and risk ratios of `"km"` at the horizons, not the model-based cumulative incidence, which would
#' assume proportional subdistribution hazards. Its robust variance treats the weights and G as fixed.
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
#' When the hazard ratio of the regimes is constant over follow-up (for example without a grace
#' period, and with an effect that does not change over time), the numerators scale the events and the
#' risk set of a strategy alike on each day, and both weights estimate the same hazard ratio; the
#' stabilized one then usually has the smaller standard error.
#'
#' Truncation: with `truncateMsmWeights = TRUE` in [createComputeWeightsArgs()] (the default), the
#' weights of `msmWeights` are truncated on their own scale with the rule of
#' [createWeightTruncation()], whatever scale `stabilized` sets for the weights of
#' [computeCensoringWeights()]. With `truncateMsmWeights = FALSE`, they are formed from the truncated
#' weights of [computeCensoringWeights()]; when the two scales differ (as with `stabilized = TRUE` and
#' `msmWeights = "unstabilized"`), the cap then does not bound them, and they can be much larger,
#' for example for a `"neverInitiate"` strategy at long horizons; see [createComputeWeightsArgs()].
#' [fitOutcomeModel()] reports the cap, largest weight and share truncated of the weights it used.
#'
#' Robust variance: with `robustVarianceEngine = "cyclops"` (the default), the Cox model is fitted with
#' each clone's person as its Cyclops subject, and the pooled logistic models (the outcome model and,
#' with competing death, the death model, stacked into one model so that their joint covariance is
#' clustered on person too) are refitted, from their estimates, to one row per person,
#' strategy-by-interval cell and response, per model (cells with a fixed hazard left out); the variance
#' is Cyclops' `vcov(fit, type = "robust")`, at the estimates of that refit. `"ccw"` computes the same
#' sandwich in R from the score residuals summed by person, at the reported estimates, without the
#' refit, so it takes less time and memory for the pooled logistic models. The two agree to rounding
#' when the outcome model is converged, as with the default `control`; with a loose convergence
#' criterion the refit moves the estimates, and the two differ accordingly.
#'
#' @param modelType            The outcome model: `"cox"` (the hazard ratio of a weighted Cox model,
#'                             with the risks of `"km"`), `"logistic"` (pooled logistic), `"km"`
#'                             (weighted Kaplan-Meier, or Aalen-Johansen with competing events), or
#'                             `"fineGray"` (the subdistribution hazard ratio of a weighted Fine-Gray
#'                             model, with the risks of `"km"`; needs competing death).
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
#'                             summary of the effect over follow-up, and it is `NA` when either
#'                             strategy has no outcomes, or only outcomes, in the first interval.
#' @param msmWeights           The weights of the Cox and pooled logistic models: `"unstabilized"` or
#'                             `"stabilized"` (see Details). The weighted Kaplan-Meier and
#'                             Aalen-Johansen estimates (of `"km"`, and the rows at the horizons of
#'                             `"cox"`) do not depend on it.
#' @param varianceEstimator    The variance estimator: `"robust"` (sandwich clustered on person, the
#'                             default), `"bootstrap"` (computed by [bootstrapCcwAnalysis()]),
#'                             `"modelBased"`, or `"stackedSandwich"` (the robust sandwich stacked with
#'                             the estimating equations of the censoring models, so that it accounts
#'                             for the estimation of the weights; [fitOutcomeModel()] then needs
#'                             `cloneCensorWeightData`, see its Details). `"stackedSandwich"` is not
#'                             available for `modelType = "fineGray"`.
#' @param robustVarianceEngine What computes the robust variance of the Cox and pooled logistic models
#'                             (with `varianceEstimator = "robust"`): `"cyclops"` (Cyclops'
#'                             `vcov(fit, type = "robust")`, fitted with the persons as subjects) or
#'                             `"ccw"` (this package's own sandwich in R, the earlier behaviour). Both
#'                             are the HC0 sandwich clustered on person and agree to rounding when the
#'                             outcome model is converged (see Details). The standard errors of the
#'                             Kaplan-Meier and Aalen-Johansen estimates (of `"km"`, and the rows at
#'                             the horizons of `"cox"`) are this package's own either way.
#' @param fineGrayCensoring    The censoring survival G of `"fineGray"`, per strategy: `"natural"` (the
#'                             default) estimates it from the clones whose follow-up ended at the end of
#'                             observation or of the risk window, with the clones weighted by their
#'                             weights; clones censored at a deviation from their strategy leave its risk
#'                             set without counting as censored, because the weights already account for
#'                             that (artificial) censoring. `"finegray"` estimates it as
#'                             [survival::finegray()] does: unweighted, from every end that is not an
#'                             outcome or a death, deviations included, so artificial censoring counts in
#'                             both the weights and G; it reproduces [survival::finegray()] and is meant
#'                             for comparisons. Censorings at a time follow the outcomes and deaths at
#'                             that time.
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
#' args <- createFitOutcomeModelArgs(horizons = c(30, 365), msmWeights = "stabilized")
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
                                      robustVarianceEngine = "cyclops",
                                      fineGrayCensoring = "natural",
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
    # Default, so settings saved before this field existed load:
    msmWeights = "unstabilized",
    varianceEstimator = NULL,
    # Default, so settings saved before this field existed load:
    robustVarianceEngine = "cyclops",
    # Default, so settings saved before this field existed load:
    fineGrayCensoring = "natural",
    bootstrapArgs = NULL,
    useCovariates = NULL,
    stratified = NULL,
    prior = NULL,
    control = NULL,
    validate = function() {
      errorMessages <- checkmate::makeAssertCollection()
      checkmate::assertChoice(self$modelType, c("cox", "logistic", "km", "fineGray"), add = errorMessages)
      checkmate::assertIntegerish(self$horizons, lower = 0, min.len = 1, any.missing = FALSE, add = errorMessages)
      checkmate::assertR6(self$timeBasis, "TimeBasis", add = errorMessages)
      checkmate::assertFlag(self$armByTimeInteraction, add = errorMessages)
      checkmate::assertChoice(self$msmWeights, c("unstabilized", "stabilized"), add = errorMessages)
      checkmate::assertChoice(self$varianceEstimator, c("robust", "bootstrap", "modelBased", "stackedSandwich"), add = errorMessages)
      checkmate::assertChoice(self$robustVarianceEngine, c("cyclops", "ccw"), add = errorMessages)
      checkmate::assertChoice(self$fineGrayCensoring, c("natural", "finegray"), add = errorMessages)
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
      if (self$varianceEstimator == "stackedSandwich" && self$modelType == "fineGray") {
        stop("varianceEstimator = \"stackedSandwich\" is not available for modelType = \"fineGray\", whose censoring survival G also depends on the weights")
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
#' Create a set of analysis choices. The specification can be saved and loaded (see
#' [saveCcwAnalysisList()]), and [runCcwAnalyses()] runs a list of them.
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
      if (self$fitOutcomeModelArgs$varianceEstimator == "stackedSandwich" && self$fitCensoringModelsArgs$hazardModelType == "cox") {
        stop("varianceEstimator = \"stackedSandwich\" does not yet support the Cox censoring model: use hazardModelType = \"logistic\"")
      }
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
#' [runCcwAnalyses()] runs all combinations of target-treatment-outcomes and analyses. The
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
#' @param refitCensoringModelsForEveryOutcome Should the censoring models be fitted for every outcome,
#'                                            on the clones of that outcome's study population (whose
#'                                            risk sets end at the outcome)? If `FALSE`, the models are
#'                                            fitted once, on the study population without an outcome
#'                                            (no persons with a prior outcome removed, risk sets that
#'                                            do not end at an outcome), and shared across outcomes. See
#'                                            [runCcwAnalyses()].
#' @param ccwDiagnosticThresholds             An object of type `CcwDiagnosticThresholds` as created
#'                                            using [createCcwDiagnosticThresholds()].
#' @param calibrationSeed                     The random seed of the empirical calibration: the null
#'                                            distribution of each calibration group is fitted by
#'                                            [EmpiricalCalibration::fitMcmcNull()], an MCMC, with this
#'                                            seed and R's default random number kinds (whatever the
#'                                            caller's [RNGkind()]), so that calibrated p-values and
#'                                            confidence intervals are the same from run to run (and for
#'                                            any number of calibration threads). The random number
#'                                            stream of the caller is restored afterwards. `NULL` runs
#'                                            the MCMC unseeded, on the current random number stream, as
#'                                            CohortMethod does. Specifications saved without it (by an
#'                                            earlier version) load with `NULL`. A change of it alone
#'                                            makes [runCcwAnalyses()] recompute only the summaries.
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
                                            ccwDiagnosticThresholds = createCcwDiagnosticThresholds(),
                                            calibrationSeed = 1) {
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
    calibrationSeed = NULL,
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
      checkmate::assertInt(self$calibrationSeed, null.ok = TRUE, add = errorMessages)
      checkmate::reportAssertions(collection = errorMessages)
      excludeColumns <- c("targetId", "treatmentId", "comparatorTreatmentId", "outcomeId", "analysisId")
      if (!is.null(self$analysesToExclude) && !any(excludeColumns %in% colnames(self$analysesToExclude))) {
        stop(sprintf("analysesToExclude should have at least one of these columns: %s", paste(excludeColumns, collapse = ", ")))
      }
      analysisIds <- sapply(self$ccwAnalysisList, function(analysis) analysis$analysisId)
      if (any(duplicated(analysisIds))) {
        stop("Duplicate analysis IDs are not allowed")
      }
      # Target-treatment-outcomes with the same cohorts share the CloneCensorWeightData (with the union
      # of their outcomes) and the names of the outcome model files:
      ttos <- self$targetTreatmentOutcomesList
      cohortKeys <- vapply(ttos, function(tto) {
        paste(tto$targetId, tto$treatmentId, if (is.null(tto$comparatorTreatmentId)) NA else tto$comparatorTreatmentId)
      }, character(1))
      conceptIds <- vapply(ttos, function(tto) {
        paste(paste(sort(tto$includedCovariateConceptIds), collapse = ","), paste(sort(tto$excludedCovariateConceptIds), collapse = ","))
      }, character(1))
      if (any(tapply(conceptIds, cohortKeys, function(x) length(unique(x))) > 1)) {
        stop("Target-treatment-outcomes with the same targetId, treatmentId and comparatorTreatmentId must have the same includedCovariateConceptIds and excludedCovariateConceptIds, because they share the data")
      }
      outcomeKeys <- unlist(lapply(seq_along(ttos), function(i) {
        paste(cohortKeys[i], vapply(ttos[[i]]$outcomes, function(outcome) as.numeric(outcome$outcomeId), numeric(1)))
      }))
      if (any(duplicated(outcomeKeys))) {
        stop("An outcome appears more than once for the same targetId, treatmentId and comparatorTreatmentId")
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
#' @details
#' A Strategus module stores the specifications as their `$toList()` and calls this function on what
#' it is given back: that list itself when the analysis specifications were built in the same
#' session, or the list that `ParallelLogger::loadSettingsFromJson()` returns for it, in which the
#' settings objects have their classes again. Both convert to the same specifications. The function
#' stops when given the `CcwAnalysesSpecifications` object itself instead of its `$toList()`, or a
#' list without `ccwAnalysisList` or `targetTreatmentOutcomesList`.
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
  # The two likely mistakes, which would otherwise fail with unrelated messages:
  if (inherits(untypedList, "CcwAnalysesSpecifications")) {
    stop("`untypedList` is a CcwAnalysesSpecifications object: pass its `$toList()` instead")
  }
  missingNames <- setdiff(c("ccwAnalysisList", "targetTreatmentOutcomesList"), names(untypedList))
  if (length(missingNames) > 0) {
    stop(sprintf(
      "`untypedList` has no %s: it should be the `$toList()` of an object created by createCcwAnalysesSpecifications()",
      paste(missingNames, collapse = " or ")
    ))
  }
  return(CcwAnalysesSpecifications$new(untypedList = untypedList))
}

#' Discretize clone follow-up into person-time intervals
#'
#' Expands each clone's follow-up into a sequence of discrete intervals.  The
#' resulting data frame drives all pooled logistic IPCW models and the weighted
#' outcome analysis.
#'
#' @param cloneCensorWeightData A `CloneCensorWeightData` object after cloning.
#' @param intervalLength Length of each interval in days.
#' @param cutPoints Optional sorted vector of interval cut points (overrides
#'   `intervalLength`).
#' @return A data frame with one row per clone-interval containing columns
#'   `subjectId`, `cloneId`, `arm`, `interval`, `tstart`, `tstop`,
#'   `personTime`, `acCensorInInterval`, `classicalCensorInInterval`, `event`.
#' @export
discretizePersonTime <- function(cloneCensorWeightData, intervalLength = 30, cutPoints = NULL) {
  checkmate::assertClass(cloneCensorWeightData, "CloneCensorWeightData")
  checkmate::assertNumeric(intervalLength, len = 1, lower = 1, any.missing = FALSE)
  if (!is.null(cutPoints)) checkmate::assertNumeric(cutPoints, sorted = TRUE, unique = TRUE)
  cohorts <- .readCloneCensorWeightTable(cloneCensorWeightData, "cohorts")
  checkmate::assertSubset(c("cloneId", "survivalTime", "event", "artificiallyCensored", "arm"),
                          names(cohorts))
  result <- lapply(seq_len(nrow(cohorts)), function(row) {
    end    <- max(0, cohorts$survivalTime[row])
    points <- if (is.null(cutPoints)) seq(0, end, by = intervalLength) else cutPoints[cutPoints <= end]
    points <- sort(unique(c(0, points, end)))
    if (length(points) < 2) points <- c(0, end)
    isFinal <- cohorts$survivalTime[row] > points[-length(points)] &
               cohorts$survivalTime[row] <= points[-1]
    data.frame(
      subjectId  = cohorts$subjectId[row],
      cloneId    = cohorts$cloneId[row],
      arm        = cohorts$arm[row],
      interval   = seq_len(length(points) - 1),
      tstart     = points[-length(points)],
      tstop      = points[-1],
      personTime = diff(points),
      acCensorInInterval = as.integer(
        cohorts$artificiallyCensored[row] == 1 & isFinal
      ),
      classicalCensorInInterval = as.integer(
        cohorts$event[row] == 0 &
        cohorts$artificiallyCensored[row] == 0 &
        isFinal
      ),
      event = as.integer(cohorts$event[row] == 1 & isFinal),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}

# Pivot sparse covariates (rowId/covariateId/covariateValue) to wide format.
# Subjects absent for a given covariate are assigned 0 (missing = zero).
.pivotCovariates <- function(covariates, allSubjectIds) {
  if (nrow(covariates) == 0) return(data.frame(subjectId = allSubjectIds))
  wide <- data.frame(subjectId = allSubjectIds)
  for (cid in sort(unique(covariates$covariateId))) {
    sub     <- covariates[covariates$covariateId == cid, , drop = FALSE]
    idx     <- match(allSubjectIds, sub$rowId)
    wide[[paste0("cov_", cid)]] <- ifelse(is.na(idx), 0, sub$covariateValue[idx])
  }
  wide
}

#' Calculate stabilized inverse probability of censoring weights
#'
#' Fits separate pooled logistic regression models for artificial censoring
#' (per-protocol non-adherence) and classical censoring (loss to follow-up),
#' then multiplies the resulting interval weights to form a cumulative combined
#' weight per clone-interval.
#'
#' Stabilized weights use a marginal numerator
#' (`P(uncensored | arm, interval)`) and a covariate-adjusted denominator
#' (`P(uncensored | arm, interval, X)`).  When no covariates are present, both
#' numerator and denominator are equivalent and all weights equal 1.
#'
#' @param personTime A person-time data frame from [discretizePersonTime()].
#' @param cloneCensorWeightData The source data object containing baseline
#'   covariates.
#' @param createCensoringWeightsArgs An optional `CreateCensoringWeightsArgs`
#'   object from [createCreateCensoringWeightsArgs()]. If `NULL`, default
#'   trimming (`c(0.01, 0.99)`) and no hard cap are applied.
#' @return A person-time data frame augmented with `acIntervalWeight`,
#'   `ccIntervalWeight`, `acWeight`, `ccWeight`, and the final trimmed/capped
#'   `weight` column.
#' @export
createCensoringWeights <- function(personTime,
                                   cloneCensorWeightData,
                                   createCensoringWeightsArgs = NULL) {
  checkmate::assertDataFrame(personTime)
  checkmate::assertClass(cloneCensorWeightData, "CloneCensorWeightData")
  if (!is.null(createCensoringWeightsArgs)) {
    checkmate::assertR6(createCensoringWeightsArgs, "CreateCensoringWeightsArgs")
  }

  args      <- createCensoringWeightsArgs
  trim      <- if (!is.null(args)) args$trim      else c(0.01, 0.99)
  maxWeight <- if (!is.null(args)) args$maxWeight else Inf
  acFormula <- if (!is.null(args)) args$artificialCensoringFormula else NULL
  ccFormula <- if (!is.null(args)) args$classicalCensoringFormula  else NULL

  covariates <- .readCloneCensorWeightTable(cloneCensorWeightData, "covariates")
  allIds     <- sort(unique(personTime$subjectId))
  wideCov    <- .pivotCovariates(covariates, allIds)
  modelData  <- merge(personTime, wideCov, by = "subjectId", all.x = TRUE)
  covNames   <- setdiff(names(wideCov), "subjectId")
  rhs        <- paste(c("factor(arm)", "interval", covNames), collapse = " + ")

  # --- Artificial censoring IPCW ---
  # Numerator: marginal (arm + interval only); denominator: covariate-adjusted.
  modelData$remainUncensored_ac <- 1L - modelData$acCensorInInterval
  numAC <- fitPooledLogistic(remainUncensored_ac ~ factor(arm) + interval, modelData)
  denAC <- fitPooledLogistic(
    if (!is.null(acFormula)) acFormula
    else stats::as.formula(paste("remainUncensored_ac ~", rhs)),
    modelData
  )
  pNA <- pmin(pmax(as.numeric(numAC$probability), 1e-6), 1 - 1e-6)
  pDA <- pmin(pmax(as.numeric(denAC$probability), 1e-6), 1 - 1e-6)
  modelData$acIntervalWeight <- pNA / pDA

  # --- Classical censoring IPCW ---
  modelData$remainUncensored_cc <- 1L - modelData$classicalCensorInInterval
  numCC <- fitPooledLogistic(remainUncensored_cc ~ factor(arm) + interval, modelData)
  denCC <- fitPooledLogistic(
    if (!is.null(ccFormula)) ccFormula
    else stats::as.formula(paste("remainUncensored_cc ~", rhs)),
    modelData
  )
  pNC <- pmin(pmax(as.numeric(numCC$probability), 1e-6), 1 - 1e-6)
  pDC <- pmin(pmax(as.numeric(denCC$probability), 1e-6), 1 - 1e-6)
  modelData$ccIntervalWeight <- pNC / pDC

  # Cumulative products within clone; combined weight is their product
  modelData$acWeight <- ave(modelData$acIntervalWeight, modelData$cloneId, FUN = cumprod)
  modelData$ccWeight <- ave(modelData$ccIntervalWeight, modelData$cloneId, FUN = cumprod)
  modelData$weight   <- modelData$acWeight * modelData$ccWeight

  modelData <- trimWeights(modelData, trim)
  if (is.finite(maxWeight)) modelData$weight <- pmin(modelData$weight, maxWeight)

  attr(modelData, "acNumeratorModel")   <- numAC$model
  attr(modelData, "acDenominatorModel") <- denAC$model
  attr(modelData, "ccNumeratorModel")   <- numCC$model
  attr(modelData, "ccDenominatorModel") <- denCC$model
  modelData
}

fitPooledLogistic <- function(formula, data) {
  outcomeName <- deparse(formula[[2]])
  engine <- "glm"
  if (requireNamespace("Cyclops", quietly = TRUE) &&
      exists("convertToCyclopsData", asNamespace("Cyclops"))) {
    cyclopsFit <- tryCatch({
      design <- stats::model.matrix(formula, data = data)
      design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
      covDF  <- data.frame(
        rowId          = rep(seq_len(nrow(data)), times = ncol(design)),
        covariateId    = rep(seq_len(ncol(design)), each = nrow(data)),
        covariateValue = as.vector(design),
        stringsAsFactors = FALSE
      )
      outDF <- data.frame(rowId = seq_len(nrow(data)), outcome = data[[outcomeName]])
      cyclopsData <- Cyclops::convertToCyclopsData(
        outcomes = outDF, covariates = covDF, addIntercept = TRUE,
        modelType = "logistic", checkRowIds = FALSE, normalize = NULL, quiet = TRUE
      )
      fit    <- Cyclops::fitCyclopsModel(
        cyclopsData,
        prior   = Cyclops::createPrior("laplace", variance = 0.01),
        control = Cyclops::createControl()
      )
      coefs  <- stats::coef(fit)
      cnames <- names(coefs)
      if (length(coefs) != ncol(design) + 1L) stop("Unexpected Cyclops coefficient shape")
      iIdx      <- which(cnames %in% c("0", "(Intercept)"))[1]
      if (is.na(iIdx)) iIdx <- 1L
      intercept <- coefs[iIdx]
      slopes    <- coefs[-iIdx]
      list(model = fit,
           probability = as.numeric(stats::plogis(as.numeric(intercept) + design %*% slopes)),
           engine = "Cyclops")
    }, error = function(e) NULL)
    if (!is.null(cyclopsFit)) return(cyclopsFit)
  }
  fit <- stats::glm(formula, data = data, family = stats::binomial())
  list(model = fit, probability = stats::predict(fit, type = "response"), engine = engine)
}

#' Trim cumulative weights to percentile limits
#' @param personTime Weighted person-time data frame.
#' @param trim Two-element probability vector.
#' @return The data frame with a trimmed `weight` column.
#' @export
trimWeights <- function(personTime, trim = c(0.01, 0.99)) {
  checkmate::assertDataFrame(personTime)
  checkmate::assertSubset("weight", names(personTime))
  limits <- stats::quantile(personTime$weight, probs = trim, na.rm = TRUE, names = FALSE)
  personTime$weight <- pmin(pmax(personTime$weight, limits[1]), limits[2])
  attr(personTime, "weightLimits") <- limits
  personTime
}

#' Summarize the weight distribution
#' @param personTime Weighted person-time data frame.
#' @return A data frame of distribution statistics by arm.
#' @export
computeWeightDistribution <- function(personTime) {
  checkmate::assertDataFrame(personTime)
  do.call(rbind, lapply(split(personTime$weight, personTime$arm), function(w) {
    data.frame(n = length(w), mean = mean(w), sd = stats::sd(w),
               min = min(w), q01 = stats::quantile(w, .01),
               median = stats::median(w), q99 = stats::quantile(w, .99),
               max = max(w))
  }))
}

#' Plot the weight distribution
#' @param personTime Weighted person-time data frame.
#' @return A `ggplot2` plot.
#' @export
plotWeightDistribution <- function(personTime) {
  ggplot2::ggplot(personTime, ggplot2::aes(x = rlang::.data$weight, colour = rlang::.data$arm)) +
    ggplot2::geom_density() + ggplot2::theme_minimal()
}

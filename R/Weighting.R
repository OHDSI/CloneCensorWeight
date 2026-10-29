#' Discretize clone follow-up into person-time intervals
#' @param cloneCensorWeightData A `CloneCensorWeightData` object after cloning.
#' @param intervalLength Length of each interval in days.
#' @param cutPoints Optional sorted vector of interval cut points.
#' @return A data frame with one row per clone-interval.
#' @export
discretizePersonTime <- function(cloneCensorWeightData, intervalLength = 30, cutPoints = NULL) {
  checkmate::assertClass(cloneCensorWeightData, "CloneCensorWeightData")
  checkmate::assertNumeric(intervalLength, len = 1, lower = 1, any.missing = FALSE)
  if (!is.null(cutPoints)) checkmate::assertNumeric(cutPoints, sorted = TRUE, unique = TRUE)
  cohorts <- .readCloneCensorWeightTable(cloneCensorWeightData, "cohorts")
  checkmate::assertSubset(c("cloneId", "survivalTime", "event", "artificiallyCensored", "arm"), names(cohorts))
  result <- lapply(seq_len(nrow(cohorts)), function(row) {
    end <- max(0, cohorts$survivalTime[row])
    points <- if (is.null(cutPoints)) seq(0, end, by = intervalLength) else cutPoints[cutPoints <= end]
    points <- sort(unique(c(0, points, end)))
    if (length(points) < 2) points <- c(0, end)
    data.frame(subjectId = cohorts$subjectId[row], cloneId = cohorts$cloneId[row],
               arm = cohorts$arm[row], interval = seq_len(length(points) - 1),
               tstart = points[-length(points)], tstop = points[-1],
               personTime = diff(points),
               censorInInterval = as.integer(cohorts$artificiallyCensored[row] == 1 &
                 cohorts$survivalTime[row] > points[-length(points)] &
                 cohorts$survivalTime[row] <= points[-1]),
               event = as.integer(cohorts$event[row] == 1 &
                 cohorts$survivalTime[row] > points[-length(points)] &
                 cohorts$survivalTime[row] <= points[-1]),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, result)
}

#' Calculate stabilized inverse probability of censoring weights
#' @param personTime A person-time data frame from [discretizePersonTime()].
#' @param cloneCensorWeightData The source data object containing covariates.
#' @param trim Optional two-element probability vector for trimming.
#' @return A person-time data frame with interval and cumulative weights.
#' @export
createCensoringWeights <- function(personTime, cloneCensorWeightData, trim = c(0.01, 0.99)) {
  checkmate::assertDataFrame(personTime)
  checkmate::assertClass(cloneCensorWeightData, "CloneCensorWeightData")
  checkmate::assertNumeric(trim, len = 2, lower = 0, upper = 1, sorted = TRUE)
  covariates <- .readCloneCensorWeightTable(cloneCensorWeightData, "covariates")
  modelData <- merge(personTime, covariates, by = "subjectId", all.x = TRUE)
  modelData$remainUncensored <- 1L - modelData$censorInInterval
  numerator <- fitPooledLogistic(remainUncensored ~ factor(arm) + interval, modelData)
  covariateNames <- setdiff(names(covariates), "subjectId")
  rhs <- paste(c("factor(arm)", "interval", covariateNames), collapse = " + ")
  denominator <- fitPooledLogistic(stats::as.formula(paste("remainUncensored ~", rhs)), modelData)
  numeratorProbability <- pmin(pmax(numerator$probability, 1e-06), 1 - 1e-06)
  denominatorProbability <- pmin(pmax(denominator$probability, 1e-06), 1 - 1e-06)
  modelData$intervalWeight <- numeratorProbability / denominatorProbability
  modelData$weight <- ave(modelData$intervalWeight, modelData$cloneId, FUN = cumprod)
  modelData <- trimWeights(modelData, trim)
  attr(modelData, "numeratorModel") <- numerator$model
  attr(modelData, "denominatorModel") <- denominator$model
  attr(modelData, "weightModelEngine") <- c(numerator$engine, denominator$engine)
  modelData
}

fitPooledLogistic <- function(formula, data) {
  probability <- NULL
  fit <- NULL
  engine <- "glm"
  if (requireNamespace("Cyclops", quietly = TRUE) &&
      exists("convertToCyclopsData", asNamespace("Cyclops"))) {
    cyclopsFit <- tryCatch({
      design <- stats::model.matrix(formula, data = data)
      design <- design[, colnames(design) != "(Intercept)", drop = FALSE]
      covariates <- data.frame(
        rowId = rep(seq_len(nrow(data)), times = ncol(design)),
        covariateId = rep(seq_len(ncol(design)), each = nrow(data)),
        covariateValue = as.vector(design),
        stringsAsFactors = FALSE
      )
      outcomes <- data.frame(rowId = seq_len(nrow(data)), outcome = data$remainUncensored)
      cyclopsData <- Cyclops::convertToCyclopsData(
        outcomes = outcomes, covariates = covariates, addIntercept = TRUE,
        modelType = "logistic", checkRowIds = FALSE, normalize = NULL, quiet = TRUE
      )
      prior <- Cyclops::createPrior("laplace", variance = 0.01)
      control <- Cyclops::createControl()
      fit <- Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control)
      coefficients <- stats::coef(fit)
      coefficientNames <- names(coefficients)
      if (length(coefficients) != ncol(design) + 1L) stop("Unexpected Cyclops coefficient shape")
      intercept <- coefficients[coefficientNames %in% c("0", "(Intercept)")][1]
      if (is.na(intercept)) intercept <- coefficients[1]
      slopes <- coefficients[setdiff(seq_along(coefficients), match(intercept, coefficients))]
      probability <- stats::plogis(as.numeric(intercept) + design %*% slopes)
      list(model = fit, probability = as.numeric(probability), engine = "Cyclops")
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
  do.call(rbind, lapply(split(personTime$weight, personTime$arm), function(weights) {
    data.frame(n = length(weights), mean = mean(weights), sd = stats::sd(weights),
               min = min(weights), q01 = stats::quantile(weights, .01),
               median = stats::median(weights), q99 = stats::quantile(weights, .99),
               max = max(weights))
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
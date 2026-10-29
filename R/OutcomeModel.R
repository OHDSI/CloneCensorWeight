#' Fit a weighted clone-censor-weight outcome model
#' @param personTime Weighted person-time data frame.
#' @param modelType Either `"cox"` or `"logistic"`.
#' @return A `CloneCensorWeightOutcomeModel` with estimates and model metadata.
#' @export
fitOutcomeModel <- function(personTime, modelType = "cox") {
  checkmate::assertDataFrame(personTime)
  checkmate::assertChoice(modelType, c("cox", "logistic"))
  checkmate::assertSubset(c("arm", "weight", "event", "tstop", "subjectId"), names(personTime))
  armLevels <- sort(unique(as.character(personTime$arm)))
  checkmate::assert_true(length(armLevels) == 2)
  armValue <- as.integer(as.character(personTime$arm) == armLevels[2])
  outcomes <- data.frame(
    rowId = seq_len(nrow(personTime)),
    y = as.integer(personTime$event),
    time = as.numeric(personTime$tstop),
    stringsAsFactors = FALSE
  )
  covariates <- data.frame(
    rowId = seq_len(nrow(personTime)),
    covariateId = 1,
    covariateValue = armValue,
    stringsAsFactors = FALSE
  )
  if (modelType == "logistic") {
    timeCovariate <- data.frame(
      rowId = seq_len(nrow(personTime)), covariateId = 2,
      covariateValue = as.numeric(personTime$tstop), stringsAsFactors = FALSE
    )
    covariates <- rbind(covariates, timeCovariate)
  }
  cyclopsData <- Cyclops::convertToCyclopsData(
    outcomes = outcomes, covariates = covariates,
    modelType = if (modelType == "cox") "cox" else "lr",
    addIntercept = modelType == "logistic", checkRowIds = FALSE,
    normalize = NULL, quiet = TRUE
  )
  fit <- Cyclops::fitCyclopsModel(
    cyclopsData, prior = Cyclops::createPrior("none"),
    control = Cyclops::createControl(), weights = personTime[["weight"]]
  )
  if (!identical(fit[["return_flag"]], "SUCCESS")) {
    stop("Cyclops outcome model did not converge: ", fit[["return_flag"]], call. = FALSE)
  }
  coefficient <- unname(stats::coef(fit)["1"])
  confidence <- stats::confint(fit, parm = 1, includePenalty = TRUE)
  standardError <- (confidence[1, 3] - confidence[1, 2]) / (2 * stats::qnorm(0.975))
  estimand <- if (modelType == "cox") "hazard ratio" else "odds ratio"
  summary <- data.frame(
    modelType = modelType, term = "arm", estimate = exp(coefficient),
    lower95 = exp(confidence[1, 2]), upper95 = exp(confidence[1, 3]),
    pValue = 2 * stats::pnorm(-abs(coefficient / standardError)),
    logLikelihood = as.numeric(fit[["log_likelihood"]]), stringsAsFactors = FALSE
  )
  structure(list(summary = summary, model = fit, estimand = estimand,
                 armLevels = armLevels, clusterVariable = "subjectId"),
            class = "CloneCensorWeightOutcomeModel")
}

#' Compute covariate means by clone arm
#' @param cloneCensorWeightData A `CloneCensorWeightData` object.
#' @param weights Optional weighted person-time data frame.
#' @return A data frame of covariate means by arm.
#' @export
computeCovariateBalance <- function(cloneCensorWeightData, weights = NULL) {
  checkmate::assertClass(cloneCensorWeightData, "CloneCensorWeightData")
  covariates <- .readCloneCensorWeightTable(cloneCensorWeightData, "covariates")
  cohorts <- .readCloneCensorWeightTable(cloneCensorWeightData, "cohorts")
  if (!is.null(weights)) cohorts <- unique(weights[c("subjectId", "arm")])
  if (!"arm" %in% names(cohorts)) cohorts$arm <- "eligible"
  numericNames <- setdiff(names(covariates)[vapply(covariates, is.numeric, logical(1))], "subjectId")
  do.call(rbind, lapply(numericNames, function(name) {
    means <- tapply(covariates[[name]][match(cohorts$subjectId, covariates$subjectId)],
                    cohorts$arm, mean, na.rm = TRUE)
    data.frame(covariate = name, arm = names(means), mean = as.numeric(means), stringsAsFactors = FALSE)
  }))
}

#' Plot covariate means by arm
#' @param balance Output from [computeCovariateBalance()].
#' @return A `ggplot2` plot.
#' @export
plotCovariateBalance <- function(balance) {
  checkmate::assertDataFrame(balance)
  ggplot2::ggplot(balance, ggplot2::aes(x = rlang::.data$mean,
                                        y = rlang::.data$covariate,
                                        colour = rlang::.data$arm)) +
    ggplot2::geom_point() + ggplot2::theme_minimal()
}

#' Plot weighted Kaplan-Meier curves
#' @param personTime Weighted person-time data frame.
#' @return A `ggplot2` plot.
#' @export
plotKaplanMeier <- function(personTime) {
  checkmate::assertDataFrame(personTime)
  weightVector <- personTime[["weight"]]
  fit <- survival::survfit(survival::Surv(tstop, event) ~ arm, data = personTime,
                           weights = weightVector)
  plotData <- data.frame(time = fit$time, survival = fit$surv,
                         arm = rep(names(fit$strata), fit$strata))
  ggplot2::ggplot(plotData, ggplot2::aes(x = rlang::.data$time,
                                         y = rlang::.data$survival,
                                         colour = rlang::.data$arm)) +
    ggplot2::geom_step() + ggplot2::theme_minimal() +
    ggplot2::labs(x = "Days since trial start", y = "Survival probability")
}
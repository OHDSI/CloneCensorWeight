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

#' Compute the calibration of the censoring models at the end of the grace period
#'
#' @description
#' Compares, for each strategy with a grace period (`"initiateWithin"` and `"notInitiateWithin"`), the
#' predicted and the observed probability of starting the treatment by the last day `G` of the grace
#' period, by group (decile) of the prediction.
#'
#' @details
#' The persons compared are those of the binary model's risk set at `G`: no outcome or death on or
#' before day `G` (before day `G` with `sameDayOrdering = "deviationFirst"`) and `riskEnd >= G`. The
#' observed value is the share of them who started the treatment on or before day `G`. The
#' predictions are:
#'
#' - `model = "binary"`: `pi(X)`, the probability of starting by `G` from the binary model at `G`
#'   (the model of an `"initiateWithin"` strategy, and of the weights after `G` of a
#'   `"notInitiateWithin"` strategy fit with `postGraceModel = "binary"`). Fitted by maximum
#'   likelihood with an intercept, it reproduces the observed share overall.
#' - `model = "hazard"`: `1 - S(G | X)` from the initiation hazard model of a `"notInitiateWithin"`
#'   strategy, the probability of initiating on days `0..G`, whose complement gives the weights with
#'   `postGraceModel = "hazard"`.
#'
#' The persons are ranked by their prediction and split into `numberOfGroups` groups of about the same
#' size (deciles by default); persons with tied predictions are kept in one group, so there can be
#' fewer groups (for example two with a single binary covariate). In a well-calibrated model the
#' observed share is close to the mean prediction in every group.
#'
#' A miscalibrated `1 - S(G | X)` means that `1 / S(G | X)` does not reweight the compliant clones of
#' the `"notInitiateWithin"` strategy after `G` to the persons under follow-up at `G`, so the
#' covariates are not balanced after the grace period (see [computeCovariateBalance()]). It can come
#' from a hazard model that assumes the same covariate effects on every day of the grace period when
#' covariates change *when* persons start (see `covariateTimeInteraction` in
#' [createFitCensoringModelsArgs()]), or from starting the treatment during the grace period changing
#' the risk of the outcome or death during the grace period (then the observed share among persons
#' still under follow-up at `G` differs from the probability of initiating by `G`).
#'
#' @param clones          An object of type [CcwClones] as created by `createClones()` (or the
#'                        weighted clones returned by [computeCensoringWeights()]).
#' @param censoringModels An object of class `CcwCensoringModels` as created by [fitCensoringModels()]
#'                        on the same clones.
#' @param numberOfGroups  The number of groups (quantiles) of the prediction.
#'
#' @return
#' A tibble with columns `strategyId`, `model` (`"binary"` or `"hazard"`), `group` (1 for the lowest
#' predictions), `persons`, `predicted` (the mean prediction in the group) and `observed` (the share
#' who started by `G`). The attribute `strategyLabels` holds the strategy labels, and `graceEnds` the
#' `G` of each strategy (`NA` without a grace period).
#'
#' @example inst/examples/clones.R
#' @examples
#' censoringModels <- fitCensoringModels(
#'   clones = clones,
#'   cloneCensorWeightData = data,
#'   fitCensoringModelsArgs = createFitCensoringModelsArgs(prior = Cyclops::createPrior("none"))
#' )
#' calibration <- computeCensoringModelCalibration(clones, censoringModels)
#' calibration
#' close(clones)
#' close(data)
#'
#' @export
computeCensoringModelCalibration <- function(clones, censoringModels, numberOfGroups = 10) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(clones, add = errorMessages)
  checkmate::assertClass(censoringModels, "CcwCensoringModels", add = errorMessages)
  checkmate::assertInt(numberOfGroups, lower = 1, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  metaData <- attr(clones, "metaData")
  .checkModelsMatchClones(censoringModels, metaData)

  strategies <- metaData$strategies
  persons <- clones$persons |>
    collect()
  results <- list()
  for (strategyId in seq_along(strategies)) {
    strategy <- strategies[[strategyId]]
    if (is.null(strategy$graceEnd)) {
      next
    }
    graceEnd <- as.integer(strategy$graceEnd)
    treatmentDay <- if (strategy$treatment == "treatment") persons$treatmentDay else persons$comparatorTreatmentDay
    riskSet <- .getBinaryRiskSet(persons, treatmentDay, graceEnd, metaData$sameDayOrdering)
    model <- censoringModels[[strategyId]]
    predictions <- list()
    if (model$kind == "binary") {
      predictions$binary <- model$predictions$propensityScore[match(riskSet$rowId, model$predictions$rowId)]
    } else {
      predictions$hazard <- 1 - .computeHazardSurvival(model, riskSet$rowId, graceEnd)
      if (!is.null(model$postGrace)) {
        predictions$binary <- model$postGrace$predictions$propensityScore[match(riskSet$rowId, model$postGrace$predictions$rowId)]
      }
    }
    for (modelName in names(predictions)) {
      predicted <- predictions[[modelName]]
      if (anyNA(predicted)) {
        stop(sprintf("The censoring model of strategy %d has no prediction for some persons under follow-up at day %d. Were the censoring models fit on these clones?", strategyId, graceEnd))
      }
      groups <- .groupByRank(predicted, numberOfGroups)
      results[[length(results) + 1]] <- tibble(
        strategyId = as.integer(strategyId),
        model = modelName,
        group = seq_len(max(groups)),
        persons = tabulate(groups),
        predicted = as.vector(tapply(predicted, groups, mean)),
        observed = as.vector(tapply(riskSet$y, groups, mean))
      )
    }
  }
  calibration <- dplyr::bind_rows(c(
    list(tibble(strategyId = integer(0), model = character(0), group = integer(0), persons = integer(0), predicted = numeric(0), observed = numeric(0))),
    results
  ))
  attr(calibration, "strategyLabels") <- vapply(strategies, function(strategy) as.character(strategy$label), character(1))
  attr(calibration, "graceEnds") <- vapply(strategies, function(strategy) if (is.null(strategy$graceEnd)) NA_integer_ else as.integer(strategy$graceEnd), integer(1))
  return(calibration)
}

# Groups 1, 2, ... of `x` of about equal size, in increasing order of `x`; tied values share the group
# of the first of them, so there can be fewer than numberOfGroups groups.
.groupByRank <- function(x, numberOfGroups) {
  group <- ceiling(rank(x, ties.method = "min") * numberOfGroups / length(x))
  return(match(group, sort(unique(group))))
}

#' Plot the calibration of the censoring models at the end of the grace period
#'
#' @description
#' Plots, per strategy with a grace period, the observed share of persons who started the treatment
#' by the end of the grace period against the mean predicted probability, per group of the
#' prediction, for each censoring model. Points on the dashed diagonal are well calibrated. Requires
#' the `ggplot2` package.
#'
#' @param calibration A tibble as created by [computeCensoringModelCalibration()].
#' @param labels      The labels of strategies 1 (target) and 2 (comparator). If `NULL`, the
#'                    `strategyLabels` attribute of `calibration` is used.
#' @param fileName    Name of the file where the plot should be saved, for example 'plot.png'. See
#'                    [ggplot2::ggsave()] for supported file formats. If `NULL`, the plot is not saved.
#'
#' @return
#' A ggplot object. Use the [ggplot2::ggsave()] function to save to file in a different format.
#'
#' @examples
#' calibration <- data.frame(
#'   strategyId = c(1, 1, 2, 2, 2, 2),
#'   model = c("binary", "binary", "hazard", "hazard", "binary", "binary"),
#'   group = c(1, 2, 1, 2, 1, 2),
#'   persons = 100,
#'   predicted = c(0.2, 0.6, 0.3, 0.7, 0.2, 0.6),
#'   observed = c(0.21, 0.59, 0.2, 0.62, 0.21, 0.59)
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plotCensoringModelCalibration(calibration)
#' }
#'
#' @export
plotCensoringModelCalibration <- function(calibration,
                                          labels = NULL,
                                          fileName = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(calibration, add = errorMessages)
  checkmate::assertNames(colnames(calibration), must.include = c("strategyId", "model", "group", "predicted", "observed"), add = errorMessages)
  checkmate::assertCharacter(labels, len = 2, any.missing = FALSE, null.ok = TRUE, add = errorMessages)
  checkmate::assertCharacter(fileName, len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required to plot the calibration. Please install it.")
  }
  if (is.null(labels)) {
    labels <- attr(calibration, "strategyLabels")
    if (is.null(labels)) {
      labels <- c("Target", "Comparator")
    }
  }

  data <- as.data.frame(calibration)
  data$strategy <- factor(labels[data$strategyId], levels = labels)
  modelLabels <- c(binary = "Binary model at the end of the grace period", hazard = "1 - S(G | X) from the initiation hazard model")
  data$modelLabel <- factor(modelLabels[data$model], levels = modelLabels)
  colors <- c("#2a78d6", "#eb6834")
  names(colors) <- modelLabels
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$predicted, y = .data$observed, color = .data$modelLabel)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#888888") +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~strategy) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Predicted probability of starting by the end of the grace period",
      y = "Observed share",
      color = "Model"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", legend.direction = "vertical", panel.grid.minor = ggplot2::element_blank())
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 7, height = 4.5, dpi = 400)
  }
  return(plot)
}

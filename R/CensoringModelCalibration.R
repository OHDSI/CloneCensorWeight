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
#' observed value is the share of them who started the treatment on or before day `G` (a strategy
#' with nobody under follow-up at `G` has no rows). The predictions are:
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
#' The observed share is taken among the persons still under follow-up at `G`. It estimates the
#' probability of starting by `G` only when starting the treatment during the grace period does not
#' change the risk of the outcome or death during the grace period. A gap between `1 - S(G | X)` and
#' the observed share therefore has two possible causes, which this diagnostic cannot tell apart:
#'
#' - The hazard model is misspecified, for example it assumes the same covariate effects on every day
#'   of the grace period when covariates change *when* persons start (see `covariateTimeInteraction`
#'   in [createFitCensoringModelsArgs()]). Then `1 / S(G | X)` does not reweight the compliant clones
#'   of the `"notInitiateWithin"` strategy after `G` to the persons under follow-up at `G`. A gap that
#'   closes when interactions are added points to this cause.
#' - Starting the treatment changes the risk of the outcome or death during the grace period. Then the
#'   persons under follow-up at `G` are a selected group, the observed share differs from the
#'   probability of starting by `G`, and the hazard model, which follows each person only until the
#'   outcome or death, is the correct one: it is the binary model, calibrated to the observed share by
#'   construction, whose weights after `G` are biased (see `postGraceModel` in
#'   [createFitCensoringModelsArgs()]).
#'
#' So good calibration of the binary model does not show that `postGraceModel = "binary"` is
#' appropriate, and a gap for the hazard model does not show that it is not. To tell the causes apart,
#' the result also holds the grace-period selection check of [computeGracePeriodSelection()] (the
#' outcome and death rates during the grace period after and without initiation), and the hazard
#' models with and without covariate-by-time interactions can be compared by calling this function on
#' each fit (as [computeCcwDiagnostics()] does).
#'
#' Besides the groups of the prediction, the result has the calibration in the large (all persons in
#' one group) and, with `strataCovariateIds`, the persons with and without each of these covariates.
#'
#' @param clones             An object of type [CcwClones] as created by `createClones()` (or the
#'                           weighted clones returned by [computeCensoringWeights()]).
#' @param censoringModels    An object of class `CcwCensoringModels` as created by
#'                           [fitCensoringModels()] on the same clones.
#' @param numberOfGroups     The number of groups (quantiles) of the prediction.
#' @param cloneCensorWeightData The object of type [CloneCensorWeightData] the clones were created
#'                           from. Needed only with `strataCovariateIds`.
#' @param strataCovariateIds IDs of covariates by which the calibration is also reported: per
#'                           covariate, group 1 holds the persons with a non-zero value and group 0 the
#'                           others.
#'
#' @return
#' A tibble of class `CcwCensoringModelCalibration` with columns `strategyId`, `model` (`"binary"` or
#' `"hazard"`), `covariateTimeInteraction` (the setting of the hazard model, `NA` for the binary
#' model), `grouping` (`"prediction"`: groups of the prediction, `group` 1 holding the lowest
#' predictions; `"overall"`: all persons, the calibration in the large; `"covariate"`: by
#' `covariateId`), `covariateId`, `group`, `persons`, `predicted` (the mean prediction in the group)
#' and `observed` (the share who started by `G`). The attribute `strategyLabels` holds the strategy
#' labels, `graceEnds` the `G` of each strategy (`NA` without a grace period), and
#' `gracePeriodSelection` the result of [computeGracePeriodSelection()].
#'
#' @example inst/examples/clones.R
#' @examples
#' censoringModels <- fitCensoringModels(
#'   clones = clones,
#'   cloneCensorWeightData = data,
#'   fitCensoringModelsArgs = createFitCensoringModelsArgs(
#'     prior = Cyclops::createPrior("laplace", exclude = c(0), variance = 1)
#'   )
#' )
#' calibration <- computeCensoringModelCalibration(
#'   clones = clones,
#'   censoringModels = censoringModels,
#'   cloneCensorWeightData = data,
#'   strataCovariateIds = 1001
#' )
#' calibration
#' close(clones)
#' close(data)
#'
#' @export
computeCensoringModelCalibration <- function(clones,
                                             censoringModels,
                                             numberOfGroups = 10,
                                             cloneCensorWeightData = NULL,
                                             strataCovariateIds = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(clones, add = errorMessages)
  checkmate::assertClass(censoringModels, "CcwCensoringModels", add = errorMessages)
  checkmate::assertInt(numberOfGroups, lower = 1, add = errorMessages)
  .assertCovariateId(strataCovariateIds, null.ok = TRUE, add = errorMessages)
  if (length(strataCovariateIds) > 0) {
    .assertCloneCensorWeightData(cloneCensorWeightData, add = errorMessages)
  }
  checkmate::reportAssertions(collection = errorMessages)
  metaData <- attr(clones, "metaData")
  .checkModelsMatchClones(censoringModels, metaData)
  modelArgs <- attr(censoringModels, "metaData")$fitCensoringModelsArgs
  interaction <- if (identical(modelArgs$hazardModelType, "cox")) "none" else modelArgs$covariateTimeInteraction
  if (is.null(interaction)) {
    interaction <- "none"
  }
  strata <- NULL
  if (length(strataCovariateIds) > 0) {
    strataCovariateIds <- as.numeric(strataCovariateIds)
    strata <- cloneCensorWeightData$covariates |>
      filter(.data$covariateId %in% strataCovariateIds, .data$covariateValue != 0) |>
      select("rowId", "covariateId") |>
      collect()
  }

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
    if (nrow(riskSet) == 0) {
      next
    }
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
      rows <- list(
        .calibrationRows(predicted, riskSet$y, .groupByRank(predicted, numberOfGroups), "prediction"),
        .calibrationRows(predicted, riskSet$y, rep(1L, length(predicted)), "overall")
      )
      for (covariateId in strataCovariateIds) {
        hasCovariate <- as.integer(riskSet$rowId %in% strata$rowId[strata$covariateId == covariateId])
        rows[[length(rows) + 1]] <- .calibrationRows(predicted, riskSet$y, hasCovariate, "covariate", covariateId)
      }
      results[[length(results) + 1]] <- dplyr::bind_rows(rows) |>
        dplyr::mutate(
          strategyId = as.integer(strategyId),
          model = modelName,
          covariateTimeInteraction = if (modelName == "hazard") interaction else NA_character_,
          .before = 1
        )
    }
  }
  calibration <- dplyr::bind_rows(c(
    list(tibble(
      strategyId = integer(0), model = character(0), covariateTimeInteraction = character(0), grouping = character(0),
      covariateId = numeric(0), group = integer(0), persons = integer(0), predicted = numeric(0), observed = numeric(0)
    )),
    results
  ))
  attr(calibration, "strategyLabels") <- vapply(strategies, function(strategy) as.character(strategy$label), character(1))
  attr(calibration, "graceEnds") <- vapply(strategies, function(strategy) if (is.null(strategy$graceEnd)) NA_integer_ else as.integer(strategy$graceEnd), integer(1))
  attr(calibration, "gracePeriodSelection") <- computeGracePeriodSelection(clones, censoringModels, numberOfGroups = numberOfGroups)
  class(calibration) <- c("CcwCensoringModelCalibration", class(calibration))
  return(calibration)
}

# The calibration of `predicted` against `y` in each group of `groups` (sorted).
.calibrationRows <- function(predicted, y, groups, grouping, covariateId = NA_real_) {
  levels <- sort(unique(groups))
  index <- match(groups, levels)
  return(tibble(
    grouping = grouping,
    covariateId = as.numeric(covariateId),
    group = as.integer(levels),
    persons = tabulate(index, nbins = length(levels)),
    predicted = as.vector(tapply(predicted, index, mean)),
    observed = as.vector(tapply(y, index, mean))
  ))
}

#' @export
print.CcwCensoringModelCalibration <- function(x, ...) {
  writeLines("Calibration of the censoring models at the end of the grace period")
  data <- as.data.frame(x)
  overall <- data[data$grouping == "overall", ]
  if (nrow(overall) > 0) {
    groups <- data[data$grouping == "prediction", ]
    key <- function(d) paste(d$strategyId, d$model, d$covariateTimeInteraction)
    gaps <- tapply(abs(groups$observed - groups$predicted), key(groups), max)
    summary <- data.frame(
      strategyId = overall$strategyId,
      model = overall$model,
      interaction = overall$covariateTimeInteraction,
      persons = overall$persons,
      observed = round(overall$observed, 4),
      predicted = round(overall$predicted, 4),
      largestGroupGap = round(as.vector(gaps[key(overall)]), 4)
    )
    writeLines("Calibration in the large, and the largest |observed - predicted| over the groups of the prediction:")
    print(summary, row.names = FALSE)
  }
  selection <- attr(x, "gracePeriodSelection")
  if (!is.null(selection)) {
    .printSelectionFlags(selection)
  }
  writeLines("")
  NextMethod()
  invisible(x)
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
#' prediction, for each censoring model (each hazard model fit, for example with and without
#' covariate-by-time interactions, gets its own line). Points on the dashed diagonal are well
#' calibrated. Requires the `ggplot2` package.
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
  if (!is.null(data$grouping)) {
    data <- data[data$grouping == "prediction", ]
  }
  data$strategy <- factor(labels[data$strategyId], levels = labels)
  # The hazard models of several fits (with and without interactions) get one label each:
  interaction <- if (is.null(data$covariateTimeInteraction)) rep(NA_character_, nrow(data)) else data$covariateTimeInteraction
  data$modelLabel <- ifelse(
    data$model == "binary",
    "Binary model at the end of the grace period",
    ifelse(is.na(interaction), "1 - S(G | X) from the initiation hazard model", sprintf("1 - S(G | X) from the initiation hazard model (interaction: %s)", interaction))
  )
  modelLabels <- unique(data$modelLabel[order(data$model != "binary", interaction)])
  data$modelLabel <- factor(data$modelLabel, levels = modelLabels)
  colors <- rep(c("#2a78d6", "#eb6834", "#1f9e6e", "#8a5cc2", "#c9a227"), length.out = length(modelLabels))
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

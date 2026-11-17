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

# Diagnostics that help choose the design options for the data at hand (Phase1Contracts section 12.5).

# Grace-period selection -----------------------------------------------------------------------------

#' Check whether starting the treatment changes the outcome and death rates during the grace period
#'
#' @description
#' Compares, for each strategy with a grace period (days `0..G`), the rates of the outcome and of death
#' during the grace period in the person-time after initiation with those in the person-time before or
#' without initiation, overall and by group (decile) of the predicted probability of starting by `G`,
#' and flags a difference.
#'
#' @details
#' `postGraceModel = "binary"` (see [createFitCensoringModelsArgs()]) assumes that starting the
#' treatment during the grace period does not change the risk of the outcome or death during the grace
#' period: otherwise the persons still under follow-up at `G`, on whom the binary model is fitted, are
#' a selected group and its weights after `G` are biased. This function checks that assumption.
#'
#' Each person is followed from day 0 to the first of day `G`, `riskEnd`, the outcome and death. The
#' days from the day of initiation on (from the day after it with `sameDayOrdering = "eventsFirst"`,
#' under which an event on the day of initiation comes first) are *initiated* person-time; the earlier
#' days, and all days of persons who do not start by then, are *not initiated* person-time. With
#' `"eventsFirst"` no person-time of day 0 is initiated, so day 0 is left out (with its events), and
#' follow-up starts on day 1. Initiators are the persons who start while followed.
#'
#' Persons who start differ from those who do not in their covariates, so the crude rates can differ
#' when starting has no effect. The rates are therefore also compared within groups of the predicted
#' probability of starting by `G`, `1 - S(G | X)` from an initiation hazard model of the same
#' treatment in `censoringModels` (of a `"notInitiateWithin(G)"` or `"neverInitiate"` strategy). The
#' summary is the Mantel-Haenszel rate ratio (initiated vs not initiated) over the strata, with the
#' Greenland-Robins variance. With `stratifyByDay = TRUE` (the default) the strata are the groups
#' crossed with the day of the grace period, so initiated and not-initiated person-time are compared on
#' the same day. Initiated person-time falls later in the grace period than not-initiated person-time,
#' so without the day strata a baseline risk that changes over the grace period (for example a risk
#' that is highest just after the index date) shows as a difference between them when starting has no
#' effect. `stratifyByDay = FALSE` compares them within the groups only. Without such a model the groups
#' are left out (one group). The flag is raised when the rate ratio of the events that end follow-up
#' (outcome or death) differs from 1 at level `alpha`. Then `"binary"` is likely biased, and
#' `postGraceModel = "hazard"` is the valid choice. A flag can also come from confounding that the
#' groups do not remove (or, with `stratifyByDay = FALSE`, from a time trend in the baseline risk), and
#' the absence of a flag does not prove the assumption when there are few events: read the rate ratio
#' and its interval, not only the flag. `postGraceModel` concerns the `"notInitiateWithin"`
#' strategies; the check of an `"initiateWithin"` strategy of the same treatment and `G` uses the same
#' person-time, and concerns `initiateWithinModel`: a flag there means that
#' `initiateWithinModel = "hazard"` is likely biased and `"binary"` is the valid choice (see
#' [createFitCensoringModelsArgs()]). This holds for `graceRegime = "binaryAtG"` (the default). With
#' `graceRegime = "naturalGrace"` both models give the probability of starting on day `G` of the
#' persons who reach day `G` untreated, the risk set of the hazard model on day `G`, so the check does
#' not choose between them.
#'
#' @param clones          An object of type [CcwClones] as created by `createClones()` (or the
#'                        weighted clones returned by [computeCensoringWeights()]).
#' @param censoringModels `NULL`, or an object of class `CcwCensoringModels` as created by
#'                        [fitCensoringModels()] on the same clones, which gives the groups.
#' @param numberOfGroups  The number of groups (quantiles) of the predicted probability of starting.
#' @param alpha           The level of the test behind the flag.
#' @param stratifyByDay   Also stratify the Mantel-Haenszel rate ratio by the day of the grace period?
#'
#' @return
#' A tibble of class `CcwGracePeriodSelection` with columns `strategyId`, `grouping` (`"overall"` or
#' `"prediction"`), `group`, `persons`, `initiators`, `initiatedDays`, `notInitiatedDays`,
#' `initiatedOutcomes`, `notInitiatedOutcomes`, `initiatedDeaths`, `notInitiatedDeaths`, the rates per
#' person-day `initiatedOutcomeRate`, `notInitiatedOutcomeRate`, `initiatedDeathRate`,
#' `notInitiatedDeathRate`, and `rateRatio`, the ratio of the rates of the outcome or death. The
#' attribute `summary` holds, per strategy and event (`"outcomeOrDeath"`, `"outcome"`, `"death"`),
#' whether the rate ratio is `stratified` over the groups, whether it is stratified `byDay`,
#' `rateRatio`, `lb95`, `ub95`, `p` and `flag`; `strategyLabels` and `strategyTypes` hold the strategy
#' labels and types.
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
#' selection <- computeGracePeriodSelection(clones, censoringModels, numberOfGroups = 5)
#' selection
#' close(clones)
#' close(data)
#'
#' @export
computeGracePeriodSelection <- function(clones, censoringModels = NULL, numberOfGroups = 10, alpha = 0.05, stratifyByDay = TRUE) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(clones, add = errorMessages)
  checkmate::assertClass(censoringModels, "CcwCensoringModels", null.ok = TRUE, add = errorMessages)
  checkmate::assertInt(numberOfGroups, lower = 1, add = errorMessages)
  checkmate::assertNumber(alpha, lower = 0, upper = 1, add = errorMessages)
  checkmate::assertFlag(stratifyByDay, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  metaData <- attr(clones, "metaData")
  if (!is.null(censoringModels)) {
    .checkModelsMatchClones(censoringModels, metaData)
  }

  strategies <- metaData$strategies
  persons <- clones$persons |>
    collect()
  eventsFirst <- metaData$sameDayOrdering == "eventsFirst"
  rows <- list()
  summaries <- list()
  for (strategyId in seq_along(strategies)) {
    strategy <- strategies[[strategyId]]
    if (is.null(strategy$graceEnd)) {
      next
    }
    graceEnd <- as.integer(strategy$graceEnd)
    treatmentDay <- if (strategy$treatment == "treatment") persons$treatmentDay else persons$comparatorTreatmentDay
    # With "eventsFirst" no person-time of day 0 is initiated, so day 0 is left out: it would compare the
    # day-0 events of all persons with none.
    firstDay <- if (eventsFirst) 1L else 0L
    personTime <- .getGracePeriodPersonTime(persons, treatmentDay, graceEnd, firstDay)
    if (nrow(personTime) == 0) {
      next
    }
    predicted <- .predictStartByGraceEnd(censoringModels, strategy, graceEnd, personTime$rowId)
    hasPrediction <- !is.na(predicted)
    groups <- rep(NA_integer_, nrow(personTime))
    if (any(hasPrediction)) {
      groups[hasPrediction] <- .groupByRank(predicted[hasPrediction], numberOfGroups)
    }
    strategyRows <- dplyr::bind_rows(
      .selectionRows(personTime, rep(1L, nrow(personTime)), "overall"),
      .selectionRows(personTime[hasPrediction, ], groups[hasPrediction], "prediction")
    )
    rows[[length(rows) + 1]] <- dplyr::mutate(strategyRows, strategyId = as.integer(strategyId), .before = 1)
    stratified <- any(hasPrediction)
    strata <- if (stratified) {
      .selectionStrata(personTime[hasPrediction, ], groups[hasPrediction], firstDay, graceEnd, stratifyByDay)
    } else {
      .selectionStrata(personTime, rep(1L, nrow(personTime)), firstDay, graceEnd, stratifyByDay)
    }
    events <- list(
      outcomeOrDeath = list(a = strata$initiatedOutcomes + strata$initiatedDeaths, b = strata$notInitiatedOutcomes + strata$notInitiatedDeaths),
      outcome = list(a = strata$initiatedOutcomes, b = strata$notInitiatedOutcomes),
      death = list(a = strata$initiatedDeaths, b = strata$notInitiatedDeaths)
    )
    for (event in names(events)) {
      ratio <- .mantelHaenszelRateRatio(events[[event]]$a, events[[event]]$b, strata$initiatedDays, strata$notInitiatedDays)
      summaries[[length(summaries) + 1]] <- tibble(
        strategyId = as.integer(strategyId),
        event = event,
        stratified = stratified,
        byDay = stratifyByDay,
        rateRatio = ratio$rateRatio,
        lb95 = ratio$lb95,
        ub95 = ratio$ub95,
        p = ratio$p,
        flag = !is.na(ratio$p) && ratio$p < alpha
      )
    }
  }
  selection <- dplyr::bind_rows(c(list(.emptySelectionRows()), rows))
  attr(selection, "summary") <- dplyr::bind_rows(c(
    list(tibble(strategyId = integer(0), event = character(0), stratified = logical(0), byDay = logical(0), rateRatio = numeric(0), lb95 = numeric(0), ub95 = numeric(0), p = numeric(0), flag = logical(0))),
    summaries
  ))
  attr(selection, "strategyLabels") <- vapply(strategies, function(strategy) as.character(strategy$label), character(1))
  attr(selection, "strategyTypes") <- vapply(strategies, function(strategy) as.character(strategy$type), character(1))
  class(selection) <- c("CcwGracePeriodSelection", class(selection))
  return(selection)
}

# Per person at risk on the first day (day 1 with "eventsFirst", else day 0; followed through day G,
# riskEnd, the outcome or death): the days not
# initiated and initiated, whether the person initiated while followed, and the event that ended follow-up
# in the grace period, with whether it fell in initiated person-time.
.getGracePeriodPersonTime <- function(persons, treatmentDay, graceEnd, firstDay) {
  eventsFirst <- firstDay == 1L
  eventDay <- pmin(persons$outcomeDay, persons$deathDay, na.rm = TRUE)
  lastDay <- pmin(graceEnd, persons$riskEnd, eventDay, na.rm = TRUE)
  # The first day of initiated person-time: an event on the day of initiation comes first with
  # "eventsFirst".
  initiatedFrom <- ifelse(is.na(treatmentDay), graceEnd + 1L, treatmentDay + if (eventsFirst) 1L else 0L)
  hasEvent <- !is.na(eventDay) & eventDay <= lastDay
  isOutcome <- hasEvent & !is.na(persons$outcomeDay) & persons$outcomeDay == eventDay
  initiatedEvent <- hasEvent & eventDay >= initiatedFrom
  personTime <- tibble(
    rowId = persons$rowId,
    initiator = !is.na(treatmentDay) & treatmentDay <= lastDay,
    notInitiatedDays = pmax(0, pmin(lastDay, initiatedFrom - 1L) - firstDay + 1),
    initiatedDays = pmax(0, lastDay - initiatedFrom + 1),
    initiatedOutcome = isOutcome & initiatedEvent,
    notInitiatedOutcome = isOutcome & !initiatedEvent,
    initiatedDeath = hasEvent & !isOutcome & initiatedEvent,
    notInitiatedDeath = hasEvent & !isOutcome & !initiatedEvent,
    # The first initiated day and the last day followed, for the day strata:
    initiatedFrom = as.integer(initiatedFrom),
    lastDay = as.integer(lastDay)
  )
  return(personTime[lastDay >= firstDay, ])
}

# 1 - S(G | X) from an initiation hazard model of the strategy's treatment whose horizon reaches G (NA
# for persons outside its risk set, or for all without such a model).
.predictStartByGraceEnd <- function(censoringModels, strategy, graceEnd, rowIds) {
  for (model in censoringModels) {
    if (model$kind == "hazard" && model$treatment == strategy$treatment && model$horizon >= graceEnd) {
      return(1 - .computeHazardSurvival(model, rowIds, graceEnd))
    }
  }
  return(rep(NA_real_, length(rowIds)))
}

.selectionRows <- function(personTime, groups, grouping) {
  if (nrow(personTime) == 0) {
    return(.emptySelectionRows()[, -1])
  }
  levels <- sort(unique(groups))
  index <- factor(match(groups, levels), levels = seq_along(levels))
  total <- function(x) as.vector(tapply(as.numeric(x), index, sum))
  rows <- tibble(
    grouping = grouping,
    group = as.integer(levels),
    persons = as.integer(table(index)),
    initiators = as.integer(total(personTime$initiator)),
    initiatedDays = total(personTime$initiatedDays),
    notInitiatedDays = total(personTime$notInitiatedDays),
    initiatedOutcomes = as.integer(total(personTime$initiatedOutcome)),
    notInitiatedOutcomes = as.integer(total(personTime$notInitiatedOutcome)),
    initiatedDeaths = as.integer(total(personTime$initiatedDeath)),
    notInitiatedDeaths = as.integer(total(personTime$notInitiatedDeath))
  )
  rate <- function(events, days) ifelse(days > 0, events / days, NA_real_)
  rows$initiatedOutcomeRate <- rate(rows$initiatedOutcomes, rows$initiatedDays)
  rows$notInitiatedOutcomeRate <- rate(rows$notInitiatedOutcomes, rows$notInitiatedDays)
  rows$initiatedDeathRate <- rate(rows$initiatedDeaths, rows$initiatedDays)
  rows$notInitiatedDeathRate <- rate(rows$notInitiatedDeaths, rows$notInitiatedDays)
  rows$rateRatio <- rate(rows$initiatedOutcomes + rows$initiatedDeaths, rows$initiatedDays) /
    rate(rows$notInitiatedOutcomes + rows$notInitiatedDeaths, rows$notInitiatedDays)
  return(rows)
}

# Per stratum (group, crossed with the day of the grace period when `byDay`): the initiated and
# not-initiated person-days and their outcomes and deaths. With "eventsFirst" no person-time of day 0
# is initiated, so day 0 (firstDay = 1) is left out.
.selectionStrata <- function(personTime, groups, firstDay, graceEnd, byDay) {
  groupLevels <- sort(unique(groups))
  group <- match(groups, groupLevels)
  nGroups <- length(groupLevels)
  nDays <- graceEnd - firstDay + 1L
  # Day index 1..nDays, plus one row for the day after the last:
  key <- function(g, day) (g - 1L) * (nDays + 1L) + (day - firstDay + 1L)
  size <- nGroups * (nDays + 1L)
  countDays <- function(start, end) {
    valid <- end >= start
    change <- tabulate(key(group[valid], start[valid]), size) - tabulate(key(group[valid], end[valid] + 1L), size)
    counts <- apply(matrix(change, nDays + 1L, nGroups), 2, cumsum)
    return(matrix(counts, nDays + 1L, nGroups)[seq_len(nDays), , drop = FALSE])
  }
  countEvents <- function(event) {
    counts <- tabulate(key(group[event], personTime$lastDay[event]), size)
    return(matrix(counts, nDays + 1L, nGroups)[seq_len(nDays), , drop = FALSE])
  }
  notInitiatedEnd <- pmin(personTime$lastDay, personTime$initiatedFrom - 1L)
  strata <- list(
    initiatedDays = countDays(pmax(personTime$initiatedFrom, firstDay), personTime$lastDay),
    notInitiatedDays = countDays(rep(firstDay, nrow(personTime)), notInitiatedEnd),
    initiatedOutcomes = countEvents(personTime$initiatedOutcome),
    notInitiatedOutcomes = countEvents(personTime$notInitiatedOutcome),
    initiatedDeaths = countEvents(personTime$initiatedDeath),
    notInitiatedDeaths = countEvents(personTime$notInitiatedDeath)
  )
  # Sum over days within each group:
  collapse <- if (byDay) as.vector else colSums
  return(lapply(strata, function(counts) as.numeric(collapse(counts))))
}

.emptySelectionRows <- function() {
  return(tibble(
    strategyId = integer(0), grouping = character(0), group = integer(0), persons = integer(0), initiators = integer(0),
    initiatedDays = numeric(0), notInitiatedDays = numeric(0), initiatedOutcomes = integer(0), notInitiatedOutcomes = integer(0),
    initiatedDeaths = integer(0), notInitiatedDeaths = integer(0), initiatedOutcomeRate = numeric(0),
    notInitiatedOutcomeRate = numeric(0), initiatedDeathRate = numeric(0), notInitiatedDeathRate = numeric(0), rateRatio = numeric(0)
  ))
}

# The Mantel-Haenszel rate ratio over strata, with events a (exposed) and b (unexposed) in person-time t1
# and t0, and the Greenland-Robins variance of its log (Rothman, Greenland and Lash, Modern
# Epidemiology, 3rd edition, chapter 15).
.mantelHaenszelRateRatio <- function(a, b, t1, t0) {
  t <- t1 + t0
  use <- t1 > 0 & t0 > 0
  a <- a[use]
  b <- b[use]
  t1 <- t1[use]
  t0 <- t0[use]
  t <- t[use]
  numerator <- sum(a * t0 / t)
  denominator <- sum(b * t1 / t)
  if (numerator <= 0 || denominator <= 0) {
    return(list(rateRatio = if (denominator > 0) 0 else NA_real_, lb95 = NA_real_, ub95 = NA_real_, p = NA_real_))
  }
  rateRatio <- numerator / denominator
  se <- sqrt(sum(t1 * t0 * (a + b) / t^2) / (numerator * denominator))
  z <- log(rateRatio) / se
  return(list(
    rateRatio = rateRatio,
    lb95 = exp(log(rateRatio) - stats::qnorm(0.975) * se),
    ub95 = exp(log(rateRatio) + stats::qnorm(0.975) * se),
    p = 2 * stats::pnorm(-abs(z))
  ))
}

.printSelectionFlags <- function(selection) {
  summary <- attr(selection, "summary")
  summary <- summary[summary$event == "outcomeOrDeath", ]
  if (nrow(summary) == 0) {
    return(invisible(NULL))
  }
  types <- attr(selection, "strategyTypes")
  writeLines("Grace-period selection check (rate ratio of outcome or death, initiated vs not initiated person-time):")
  for (i in seq_len(nrow(summary))) {
    row <- summary[i, ]
    strata <- c(if (row$stratified) "group", if (isTRUE(row$byDay)) "day")
    flag <- ""
    if (isTRUE(row$flag)) {
      flag <- " FLAG: initiated and not-initiated person-time differ in risk"
      if (is.null(types) || types[row$strategyId] == "notInitiateWithin") {
        flag <- paste0(flag, "; postGraceModel = \"binary\" is likely biased, \"hazard\" is valid")
      } else if (types[row$strategyId] == "initiateWithin") {
        flag <- paste0(flag, "; initiateWithinModel = \"hazard\" is likely biased, \"binary\" is valid (with graceRegime = \"binaryAtG\"; with \"naturalGrace\" both are fit on the persons who reach day G untreated)")
      }
    }
    writeLines(sprintf(
      "  Strategy %d: %s %.2f (%.2f-%.2f), p = %s%s",
      row$strategyId,
      if (length(strata) == 0) "crude" else paste("stratified by", paste(strata, collapse = " and ")),
      row$rateRatio, row$lb95, row$ub95,
      format.pval(row$p, digits = 2),
      flag
    ))
  }
  invisible(NULL)
}

#' @export
print.CcwGracePeriodSelection <- function(x, ...) {
  writeLines("Outcome and death rates during the grace period, initiated vs not initiated person-time")
  .printSelectionFlags(x)
  writeLines("")
  NextMethod()
  invisible(x)
}

#' Plot the grace-period selection check
#'
#' @description
#' Plots, per strategy with a grace period, the rate of the outcome or death during the grace period
#' in initiated and in not-initiated person-time, per group of the predicted probability of starting
#' (overall when there are no groups). Requires the `ggplot2` package.
#'
#' @param gracePeriodSelection A tibble as created by [computeGracePeriodSelection()].
#' @param labels               The labels of strategies 1 (target) and 2 (comparator). If `NULL`, the
#'                             `strategyLabels` attribute of `gracePeriodSelection` is used.
#' @param fileName             Name of the file where the plot should be saved, for example
#'                             'plot.png'. See [ggplot2::ggsave()] for supported file formats. If
#'                             `NULL`, the plot is not saved.
#'
#' @return
#' A ggplot object. Use the [ggplot2::ggsave()] function to save to file in a different format.
#'
#' @examples
#' selection <- data.frame(
#'   strategyId = 2,
#'   grouping = "prediction",
#'   group = 1:3,
#'   initiatedDays = c(100, 300, 600),
#'   notInitiatedDays = c(2000, 1500, 900),
#'   initiatedOutcomes = c(0, 1, 2),
#'   notInitiatedOutcomes = c(8, 9, 10),
#'   initiatedDeaths = 0,
#'   notInitiatedDeaths = c(2, 1, 1)
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plotGracePeriodSelection(selection)
#' }
#'
#' @export
plotGracePeriodSelection <- function(gracePeriodSelection, labels = NULL, fileName = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(gracePeriodSelection, add = errorMessages)
  checkmate::assertNames(colnames(gracePeriodSelection), must.include = c(
    "strategyId", "grouping", "group", "initiatedDays", "notInitiatedDays", "initiatedOutcomes",
    "notInitiatedOutcomes", "initiatedDeaths", "notInitiatedDeaths"
  ), add = errorMessages)
  checkmate::assertCharacter(labels, len = 2, any.missing = FALSE, null.ok = TRUE, add = errorMessages)
  checkmate::assertCharacter(fileName, len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required to plot the grace-period selection check. Please install it.")
  }
  labels <- .strategyLabels(labels, gracePeriodSelection)
  data <- as.data.frame(gracePeriodSelection)
  strategyIds <- unique(data$strategyId)
  keep <- unlist(lapply(strategyIds, function(id) {
    rows <- which(data$strategyId == id)
    grouped <- rows[data$grouping[rows] == "prediction"]
    if (length(grouped) > 0) grouped else rows
  }))
  data <- data[keep, ]
  rate <- function(events, days) ifelse(days > 0, events / days, NA_real_)
  long <- data.frame(
    strategy = factor(rep(labels[data$strategyId], 2), levels = labels),
    group = rep(data$group, 2),
    personTime = factor(rep(c("Initiated", "Not initiated"), each = nrow(data)), levels = c("Initiated", "Not initiated")),
    rate = 1000 * c(
      rate(data$initiatedOutcomes + data$initiatedDeaths, data$initiatedDays),
      rate(data$notInitiatedOutcomes + data$notInitiatedDeaths, data$notInitiatedDays)
    )
  )
  colors <- c(Initiated = "#2a78d6", `Not initiated` = "#eb6834")
  plot <- ggplot2::ggplot(long[!is.na(long$rate), ], ggplot2::aes(x = .data$group, y = .data$rate, color = .data$personTime)) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~strategy) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      x = "Group of the predicted probability of starting by the end of the grace period",
      y = "Outcome or death per 1,000 person-days",
      color = "Person-time"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 7, height = 4.5, dpi = 400)
  }
  return(plot)
}

# The strategy labels: `labels`, else the strategyLabels attribute, else Target and Comparator.
.strategyLabels <- function(labels, x) {
  if (is.null(labels)) {
    labels <- attr(x, "strategyLabels")
    if (is.null(labels)) {
      labels <- c("Target", "Comparator")
    }
  }
  return(labels)
}

# Positivity ---------------------------------------------------------------------------------------

#' Summarize the predicted probabilities of the censoring models (positivity)
#'
#' @description
#' Summarizes, for each strategy and censoring model, the distribution of the predicted probability of
#' starting the treatment by the end of the model's horizon, and the share of persons whose predicted
#' probability of following the strategy is extreme.
#'
#' @details
#' For a binary model at `G` (of an `"initiateWithin"` strategy, and of the weights after `G` of a
#' `"notInitiateWithin"` strategy with `postGraceModel = "binary"`), the prediction is `pi(X)`, the
#' probability of starting by `G`, over the persons under follow-up at `G`. For an initiation hazard
#' model it is `1 - S(H | X)`, the probability of starting on days `0..H` for its horizon `H` (`G` of a
#' `"notInitiateWithin"` strategy, or of an `"initiateWithin"` strategy with
#' `initiateWithinModel = "hazard"`, the largest `riskEnd` of a `"neverInitiate"` strategy), over the
#' persons at risk on day 0. The probability of following the strategy is `pi(X)` (or `1 - S(G | X)`)
#' for `"initiateWithin"` and `S(H | X)` (or `1 - pi(X)`) for the strategies of not initiating; its
#' inverse is the weight, so persons with a probability below `threshold` get weights above
#' `1 / threshold`.
#'
#' For an `"initiateWithin"` strategy with `graceStart = S > 0` the binary model at `G` is fit on the
#' persons who have not started before day `S`, and the hazard model's prediction at `G` is
#' `1 - S(G | X) / S(S - 1 | X)`, the probability of starting on days `S..G` of those persons; the
#' result also has the probability of starting before day `S`, `1 - S(S - 1 | X)` (`day = S - 1`),
#' from the hazard model of the days before `S` (or the hazard model through `G`), whose complement is
#' the probability of following the strategy. With `graceRegime = "naturalGrace"` (see
#' [createFitCensoringModelsArgs()]) the prediction at `G` is the probability of starting on day `G` of
#' the persons who reach it untreated (binary model fit on them, or `1 - S(G | X) / S(G - 1 | X)`).
#'
#' @param clones          An object of type [CcwClones] as created by `createClones()`.
#' @param censoringModels An object of class `CcwCensoringModels` as created by [fitCensoringModels()]
#'                        on the same clones.
#' @param threshold       Probabilities of following a strategy below this value count as extreme
#'                        (see `positivityThreshold` in [createCcwDiagnosticThresholds()]).
#'
#' @return
#' A tibble of class `CcwPositivity` with columns `strategyId`, `model` (`"binary"` or `"hazard"`),
#' `graceRegime` (that of the models of an `"initiateWithin"` strategy, `NA` for the other strategies),
#' `day` (the horizon of the prediction), `persons`, the `mean`, `min`, `p01`, `p05`, `p25`, `median`,
#' `p75`, `p95`, `p99` and `max` of the predicted probability of starting by `day`, `shareBelow` and
#' `shareAbove` (the shares of predictions below `threshold` and above `1 - threshold`), and
#' `shareExtreme`, the share whose probability of following the strategy is below `threshold`. The
#' attribute `predictions` holds the predictions (`strategyId`, `model`, `graceRegime`, `day`, `rowId`,
#' `probability`), and
#' `strategyLabels` the strategy labels.
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
#' positivity <- computePositivity(clones, censoringModels)
#' positivity
#' close(clones)
#' close(data)
#'
#' @export
computePositivity <- function(clones, censoringModels, threshold = 0.01) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(clones, add = errorMessages)
  checkmate::assertClass(censoringModels, "CcwCensoringModels", add = errorMessages)
  checkmate::assertNumber(threshold, lower = 0, upper = 0.5, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  metaData <- attr(clones, "metaData")
  .checkModelsMatchClones(censoringModels, metaData)

  strategies <- metaData$strategies
  rows <- list()
  predictions <- list()
  for (strategyId in seq_along(strategies)) {
    model <- censoringModels[[strategyId]]
    strategy <- strategies[[strategyId]]
    isInitiateWithin <- strategy$type == "initiateWithin"
    graceRegime <- NA_character_
    firstDay <- 0L
    if (isInitiateWithin) {
      graceRegime <- if (is.null(model$graceRegime)) "binaryAtG" else model$graceRegime
      firstDay <- .binaryModelFirstDay(strategy, graceRegime)
    }
    # Each entry: the model, its day, the predictions and whether the probability of following the
    # strategy is the prediction (starting) or its complement (not starting):
    entries <- list()
    if (model$kind == "binary") {
      entries[[1]] <- list(model = "binary", day = model$horizon, rowId = model$predictions$rowId, probability = model$predictions$propensityScore, starting = isInitiateWithin)
    } else {
      probability <- if (isInitiateWithin) {
        .hazardStartProbability(model, model$predictions$rowId, firstDay, model$horizon)
      } else {
        1 - .computeHazardSurvival(model, model$predictions$rowId, model$horizon)
      }
      entries[[1]] <- list(model = "hazard", day = model$horizon, rowId = model$predictions$rowId, probability = probability, starting = isInitiateWithin)
      if (!is.null(model$postGrace)) {
        entries[[2]] <- list(model = "binary", day = model$postGrace$horizon, rowId = model$postGrace$predictions$rowId, probability = model$postGrace$predictions$propensityScore, starting = FALSE)
      }
    }
    # initiateWithin with graceStart > 0: also the probability of starting before graceStart, from the
    # hazard model of those days (or the hazard model through G):
    graceStart <- if (is.null(model$graceStart)) 0L else as.integer(model$graceStart)
    if (isInitiateWithin && graceStart > 0) {
      earlyModel <- if (model$kind == "hazard") model else model$early
      entries[[length(entries) + 1]] <- list(
        model = "hazard", day = graceStart - 1L, rowId = earlyModel$predictions$rowId,
        probability = 1 - .computeHazardSurvival(earlyModel, earlyModel$predictions$rowId, graceStart - 1L), starting = FALSE
      )
    }
    for (entry in entries) {
      p <- entry$probability
      following <- if (entry$starting) p else 1 - p
      quantiles <- stats::quantile(p, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), names = FALSE, type = 7)
      rows[[length(rows) + 1]] <- tibble(
        strategyId = as.integer(strategyId),
        model = entry$model,
        graceRegime = graceRegime,
        day = as.integer(entry$day),
        persons = length(p),
        mean = mean(p),
        min = min(p),
        p01 = quantiles[1],
        p05 = quantiles[2],
        p25 = quantiles[3],
        median = quantiles[4],
        p75 = quantiles[5],
        p95 = quantiles[6],
        p99 = quantiles[7],
        max = max(p),
        shareBelow = mean(p < threshold),
        shareAbove = mean(p > 1 - threshold),
        shareExtreme = mean(following < threshold)
      )
      predictions[[length(predictions) + 1]] <- tibble(strategyId = as.integer(strategyId), model = entry$model, graceRegime = graceRegime, day = as.integer(entry$day), rowId = entry$rowId, probability = p)
    }
  }
  positivity <- dplyr::bind_rows(rows)
  attr(positivity, "predictions") <- dplyr::bind_rows(predictions)
  attr(positivity, "threshold") <- threshold
  attr(positivity, "strategyLabels") <- vapply(strategies, function(strategy) as.character(strategy$label), character(1))
  class(positivity) <- c("CcwPositivity", class(positivity))
  return(positivity)
}

#' @export
print.CcwPositivity <- function(x, ...) {
  writeLines(sprintf("Predicted probability of starting by the end of each model's horizon (extreme: probability of following the strategy below %s)", format(attr(x, "threshold"))))
  NextMethod()
  invisible(x)
}

#' Plot the predicted probabilities of the censoring models
#'
#' @description
#' Plots, per strategy and censoring model, the distribution (histogram) of the predicted probability
#' of starting the treatment by the end of the model's horizon. Requires the `ggplot2` package.
#'
#' @param positivity A tibble as created by [computePositivity()] (with its `predictions` attribute).
#' @param labels     The labels of strategies 1 (target) and 2 (comparator). If `NULL`, the
#'                   `strategyLabels` attribute of `positivity` is used.
#' @param fileName   Name of the file where the plot should be saved, for example 'plot.png'. See
#'                   [ggplot2::ggsave()] for supported file formats. If `NULL`, the plot is not saved.
#'
#' @return
#' A ggplot object. Use the [ggplot2::ggsave()] function to save to file in a different format.
#'
#' @examples
#' positivity <- data.frame(strategyId = 1:2)
#' attr(positivity, "predictions") <- data.frame(
#'   strategyId = rep(1:2, each = 50),
#'   model = rep(c("binary", "hazard"), each = 50),
#'   probability = c(stats::runif(50, 0.1, 0.6), stats::runif(50, 0.2, 0.7))
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plotPositivity(positivity)
#' }
#'
#' @export
plotPositivity <- function(positivity, labels = NULL, fileName = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(positivity, add = errorMessages)
  checkmate::assertDataFrame(attr(positivity, "predictions"), .var.name = "attr(positivity, \"predictions\")", add = errorMessages)
  checkmate::assertCharacter(labels, len = 2, any.missing = FALSE, null.ok = TRUE, add = errorMessages)
  checkmate::assertCharacter(fileName, len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required to plot the predicted probabilities. Please install it.")
  }
  labels <- .strategyLabels(labels, positivity)
  data <- as.data.frame(attr(positivity, "predictions"))
  data$strategy <- factor(labels[data$strategyId], levels = labels)
  modelLabels <- c(binary = "pi(X), binary model", hazard = "1 - S(H | X), hazard model")
  data$modelLabel <- modelLabels[data$model]
  # Several predictions of one model kind for a strategy (an initiateWithin strategy with graceStart > 0
  # also has the probability of starting before graceStart) are told apart by their day:
  if (!is.null(data$day) && anyDuplicated(unique(data[, c("strategyId", "model", "day")])[, c("strategyId", "model")]) > 0) {
    data$modelLabel <- sprintf("%s, day %d", data$modelLabel, as.integer(data$day))
  }
  levels <- unique(data$modelLabel[if (is.null(data$day)) order(data$model) else order(data$model, data$day)])
  data$modelLabel <- factor(data$modelLabel, levels = levels)
  colors <- rep(c("#2a78d6", "#eb6834", "#1f9e6e", "#8a5cc2"), length.out = length(levels))
  names(colors) <- levels
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$probability, fill = .data$modelLabel)) +
    ggplot2::geom_histogram(bins = 40, boundary = 0, alpha = 0.6, position = "identity") +
    ggplot2::facet_wrap(~strategy, scales = "free_y") +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Predicted probability of starting by the end of the horizon", y = "Persons", fill = "Model") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 7, height = 4.5, dpi = 400)
  }
  return(plot)
}

# Weights -----------------------------------------------------------------------------------------

#' Summarize the censoring weights on both scales, before and after truncation
#'
#' @description
#' Summarizes, per strategy, scale (stabilized and unstabilized), truncation (before and after) and
#' band of days, the weights `W(d)` of the clones under follow-up on the days of the band: quantiles,
#' maximum, mean, effective sample size, the share truncated, and the cap.
#'
#' @details
#' The weight of a clone on a day is evaluated as in [computeWeightDistribution()], exactly, from the
#' censoring models stored in the weighted clones, on both scales: the scale of the weights of
#' [computeCensoringWeights()] (`stabilized` in [createComputeWeightsArgs()]) and the other one. After
#' truncation, the weights on the scale of [computeCensoringWeights()] are capped at its cap, and those
#' on the other scale at their own cap with `truncateMsmWeights = TRUE`, or are formed from the
#' truncated weights with `truncateMsmWeights = FALSE` (`cap` is then `NA`, and `maxWeight` is the
#' effective cap). The Cox and pooled logistic outcome models use the weights after truncation on the
#' scale of `msmWeights` in [createFitOutcomeModelArgs()]; the Kaplan-Meier and Aalen-Johansen risks use
#' those on the scale of [computeCensoringWeights()].
#'
#' The statistics of a band pool the clone-days of the evaluated days in it (each clone counts once per
#' day). The mean of the stabilized weights should be close to 1 on every day. The effective sample
#' size, `sum(w)^2 / sum(w^2)`, is computed per day and averaged over the band's days;
#' `effectiveSampleShare` is the average of its ratio to the number of clones. `shareTruncated` is the
#' share of the clone-days whose weight truncation lowered.
#'
#' @param weightedClones An object of type [CcwClones] as returned by [computeCensoringWeights()].
#' @param days           The days (since index) at which the weights are evaluated. If `NULL`, the
#'                       last days of all segments (the days at which the outcome models use the
#'                       weights), as in [computeWeightDistribution()].
#' @param dayBandStarts  The first days of the bands. If `NULL`, the grace periods (from day 0), the day
#'                       after them (`tau + 1`, where `tau` is the largest `graceEnd`, the day of
#'                       [computeCovariateBalance()]) and the days after that; without a grace period,
#'                       one band.
#'
#' @return
#' A tibble of class `CcwWeightSummary` with columns `strategyId`, `scale` (`"stabilized"` or
#' `"unstabilized"`), `truncation` (`"before"` or `"after"`), `dayBand` (a label), `bandStart`,
#' `cloneDays`, `days`, `meanWeight`, `minWeight`, `p25Weight`, `medianWeight`, `p75Weight`,
#' `p95Weight`, `p99Weight`, `maxWeight`, `effectiveSampleSize`, `effectiveSampleShare`,
#' `shareTruncated` (`NA` before truncation) and `cap` (the cap on that scale; `NA` before truncation and
#' for a scale not truncated on its own scale). Attributes: `strategyLabels`, `stabilized` (the scale
#' of the weights of [computeCensoringWeights()]) and `truncateMsmWeights`.
#'
#' @example inst/examples/clones.R
#' @example inst/examples/weightedClones.R
#' @examples
#' weightSummary <- computeWeightSummary(weightedClones)
#' weightSummary
#' close(weightedClones)
#' close(clones)
#' close(data)
#'
#' @export
computeWeightSummary <- function(weightedClones, days = NULL, dayBandStarts = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(weightedClones, weighted = TRUE, add = errorMessages)
  checkmate::assertIntegerish(days, any.missing = FALSE, null.ok = TRUE, add = errorMessages)
  checkmate::assertIntegerish(dayBandStarts, any.missing = FALSE, min.len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  metaData <- attr(weightedClones, "metaData")
  graceEnds <- unlist(lapply(metaData$strategies, function(strategy) strategy$graceEnd))
  if (is.null(dayBandStarts)) {
    dayBandStarts <- if (length(graceEnds) == 0) 0L else c(0L, max(graceEnds) + 1L, max(graceEnds) + 2L)
  }
  dayBandStarts <- sort(unique(as.integer(dayBandStarts)))
  on.exit({
    .dropTables(weightedClones, "ccwD")
    .dropTables(weightedClones, "ccwW")
  })
  .createCloneDayTables(weightedClones, days)
  weightedClones$ccwDBands <- tibble(
    bandStart = dayBandStarts,
    bandEnd = c(dayBandStarts[-1] - 1L, NA_integer_)
  )
  stabilized <- metaData$weights$stabilized
  truncateMsmWeights <- isTRUE(metaData$weights$truncateMsmWeights)
  weight <- .weightSql(stabilized)
  otherWeight <- .weightSql(!stabilized)
  truncated <- .truncatedWeightSql(weight)
  otherTruncated <- if (truncateMsmWeights) {
    .truncatedWeightSql(otherWeight, "ccwWOtherCaps")
  } else {
    sprintf("(%s) %s ccwWSurvival.numeratorSurvival", truncated, if (stabilized) "/" else "*")
  }
  .andromedaExecute(weightedClones, sprintf("
    CREATE TABLE ccwDWeights AS
    SELECT ccwDCloneDays.strategyId,
      ccwDCloneDays.day,
      ccwDBands.bandStart,
      %s AS weightBefore,
      %s AS weightAfter,
      %s AS otherBefore,
      %s AS otherAfter
    FROM ccwDCloneDays
    INNER JOIN ccwWSurvival
      ON ccwDCloneDays.strategyId = ccwWSurvival.strategyId
        AND ccwDCloneDays.rowId = ccwWSurvival.rowId
        AND ccwDCloneDays.evalDay = ccwWSurvival.evalDay
    INNER JOIN ccwDBands
      ON ccwDCloneDays.day >= ccwDBands.bandStart
        AND (ccwDBands.bandEnd IS NULL OR ccwDCloneDays.day <= ccwDBands.bandEnd)
    LEFT JOIN ccwWCaps
      ON ccwDCloneDays.strategyId = ccwWCaps.strategyId
    LEFT JOIN ccwWOtherCaps
      ON ccwDCloneDays.strategyId = ccwWOtherCaps.strategyId;", weight, truncated, otherWeight, otherTruncated))

  scales <- if (stabilized) c("stabilized", "unstabilized") else c("unstabilized", "stabilized")
  caps <- metaData$weights$caps
  combinations <- list(
    list(scale = scales[1], truncation = "before", column = "weightBefore", before = NULL, cap = NULL),
    list(scale = scales[1], truncation = "after", column = "weightAfter", before = "weightBefore", cap = caps$cap),
    list(scale = scales[2], truncation = "before", column = "otherBefore", before = NULL, cap = NULL),
    list(scale = scales[2], truncation = "after", column = "otherAfter", before = "otherBefore", cap = if (truncateMsmWeights) caps$otherScaleCap else NULL)
  )
  summary <- dplyr::bind_rows(lapply(combinations, function(combination) {
    rows <- .summarizeWeightColumn(weightedClones, combination$column, combination$before)
    rows$scale <- combination$scale
    rows$truncation <- combination$truncation
    rows$cap <- if (is.null(combination$cap)) NA_real_ else as.numeric(combination$cap[match(rows$strategyId, caps$strategyId)])
    rows
  }))
  bandEnds <- c(dayBandStarts[-1] - 1L, max(c(summary$lastDay, dayBandStarts)))
  bandLabels <- ifelse(dayBandStarts == bandEnds, as.character(dayBandStarts), sprintf("%d-%d", dayBandStarts, bandEnds))
  summary$dayBand <- bandLabels[match(summary$bandStart, dayBandStarts)]
  summary <- summary[order(summary$strategyId, match(summary$scale, scales), summary$truncation != "before", summary$bandStart), c(
    "strategyId", "scale", "truncation", "dayBand", "bandStart", "cloneDays", "days", "meanWeight", "minWeight", "p25Weight",
    "medianWeight", "p75Weight", "p95Weight", "p99Weight", "maxWeight", "effectiveSampleSize", "effectiveSampleShare",
    "shareTruncated", "cap"
  )]
  attr(summary, "strategyLabels") <- vapply(metaData$strategies, function(strategy) as.character(strategy$label), character(1))
  attr(summary, "stabilized") <- stabilized
  attr(summary, "truncateMsmWeights") <- truncateMsmWeights
  class(summary) <- c("CcwWeightSummary", class(summary))
  return(summary)
}

# Per strategy and band, the statistics of `column` of ccwDWeights, with the share of clone-days whose
# weight is below `before` (the column before truncation; NULL: NA).
.summarizeWeightColumn <- function(clones, column, before) {
  shareTruncated <- if (is.null(before)) "CAST(NULL AS DOUBLE)" else sprintf("AVG(CASE WHEN %s > %s * (1 + 1e-9) THEN 1.0 ELSE 0.0 END)", before, column)
  statistics <- .andromedaSql(clones, sprintf("
    SELECT strategyId,
      bandStart,
      COUNT(*) AS cloneDays,
      COUNT(DISTINCT day) AS days,
      MAX(day) AS lastDay,
      AVG(%1$s) AS meanWeight,
      MIN(%1$s) AS minWeight,
      quantile_cont(%1$s, 0.25) AS p25Weight,
      quantile_cont(%1$s, 0.5) AS medianWeight,
      quantile_cont(%1$s, 0.75) AS p75Weight,
      quantile_cont(%1$s, 0.95) AS p95Weight,
      quantile_cont(%1$s, 0.99) AS p99Weight,
      MAX(%1$s) AS maxWeight,
      %2$s AS shareTruncated
    FROM ccwDWeights
    GROUP BY strategyId,
      bandStart;", column, shareTruncated))
  effective <- .andromedaSql(clones, sprintf("
    SELECT strategyId,
      bandStart,
      AVG(ess) AS effectiveSampleSize,
      AVG(ess / n) AS effectiveSampleShare
    FROM (
      SELECT strategyId,
        bandStart,
        day,
        POWER(SUM(%1$s), 2) / SUM(POWER(%1$s, 2)) AS ess,
        COUNT(*) AS n
      FROM ccwDWeights
      GROUP BY strategyId,
        bandStart,
        day
    ) dailyEss
    GROUP BY strategyId,
      bandStart;", column))
  statistics <- dplyr::left_join(statistics, effective, by = c("strategyId", "bandStart"))
  statistics$strategyId <- as.integer(statistics$strategyId)
  statistics$bandStart <- as.integer(statistics$bandStart)
  statistics$cloneDays <- as.integer(statistics$cloneDays)
  statistics$days <- as.integer(statistics$days)
  statistics$shareTruncated <- as.numeric(statistics$shareTruncated)
  return(dplyr::as_tibble(statistics))
}

#' @export
print.CcwWeightSummary <- function(x, ...) {
  writeLines(sprintf(
    "Weights per strategy, scale and band of days (weights of computeCensoringWeights(): %s; other scale truncated on its own scale: %s)",
    if (isTRUE(attr(x, "stabilized"))) "stabilized" else "unstabilized",
    if (isTRUE(attr(x, "truncateMsmWeights"))) "yes" else "no"
  ))
  NextMethod()
  invisible(x)
}

#' Plot the summary of the censoring weights
#'
#' @description
#' Plots, per strategy and scale, the median, 99th percentile and maximum weight in each band of days,
#' before and after truncation (log scale). Requires the `ggplot2` package.
#'
#' @param weightSummary A tibble as created by [computeWeightSummary()].
#' @param labels        The labels of strategies 1 (target) and 2 (comparator). If `NULL`, the
#'                      `strategyLabels` attribute of `weightSummary` is used.
#' @param fileName      Name of the file where the plot should be saved, for example 'plot.png'. See
#'                      [ggplot2::ggsave()] for supported file formats. If `NULL`, the plot is not
#'                      saved.
#'
#' @return
#' A ggplot object. Use the [ggplot2::ggsave()] function to save to file in a different format.
#'
#' @examples
#' weightSummary <- data.frame(
#'   strategyId = rep(1:2, each = 4),
#'   scale = rep(c("stabilized", "unstabilized"), each = 2, times = 2),
#'   truncation = rep(c("before", "after"), 4),
#'   dayBand = "8-60",
#'   bandStart = 8,
#'   medianWeight = c(1, 1, 2, 2, 1, 1, 1.5, 1.5),
#'   p99Weight = c(4, 3, 9, 7, 3, 3, 5, 4),
#'   maxWeight = c(9, 3, 20, 7, 6, 3, 9, 4)
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plotWeightSummary(weightSummary)
#' }
#'
#' @export
plotWeightSummary <- function(weightSummary, labels = NULL, fileName = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(weightSummary, add = errorMessages)
  checkmate::assertNames(colnames(weightSummary), must.include = c("strategyId", "scale", "truncation", "dayBand", "bandStart", "medianWeight", "p99Weight", "maxWeight"), add = errorMessages)
  checkmate::assertCharacter(labels, len = 2, any.missing = FALSE, null.ok = TRUE, add = errorMessages)
  checkmate::assertCharacter(fileName, len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required to plot the weight summary. Please install it.")
  }
  labels <- .strategyLabels(labels, weightSummary)
  data <- as.data.frame(weightSummary)
  bands <- unique(data$dayBand[order(data$bandStart)])
  statistics <- c(medianWeight = "Median", p99Weight = "99th percentile", maxWeight = "Maximum")
  long <- data.frame(
    strategy = factor(rep(labels[data$strategyId], length(statistics)), levels = labels),
    scale = rep(data$scale, length(statistics)),
    truncation = factor(rep(data$truncation, length(statistics)), levels = c("before", "after")),
    dayBand = factor(rep(data$dayBand, length(statistics)), levels = bands),
    statistic = factor(rep(statistics, each = nrow(data)), levels = statistics),
    weight = unlist(data[names(statistics)], use.names = FALSE)
  )
  colors <- c(before = "#eb6834", after = "#2a78d6")
  plot <- ggplot2::ggplot(long, ggplot2::aes(x = .data$dayBand, y = .data$weight, color = .data$truncation, shape = .data$statistic)) +
    ggplot2::geom_point(size = 2.2, position = ggplot2::position_dodge(width = 0.4)) +
    ggplot2::facet_grid(scale ~ strategy) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Days since index", y = "Weight (log scale)", color = "Truncation", shape = "Statistic") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 7, height = 5, dpi = 400)
  }
  return(plot)
}

# Sensitivity --------------------------------------------------------------------------------------

#' Refit the analysis under alternative design choices and tabulate the estimates
#'
#' @description
#' Fits the censoring models, the weights and the outcome model under the specified settings and under
#' each variant of them, and returns the estimates side by side. It refits, so it costs up to one
#' analysis per variant; censoring models and weights are shared between variants that do not change
#' them.
#'
#' @details
#' By default the variants change one choice at a time from the specified settings, to each of its
#' other values: `postGraceModel`, `initiateWithinModel`, `graceRegime` and `covariateTimeInteraction`
#' (pooled logistic hazard model only) of [createFitCensoringModelsArgs()] (only when a strategy uses
#' them), `truncateMsmWeights` (with
#' truncation, not for `modelType = "km"`, which does not use the weights on the MSM scale) and the
#' truncation `basis` of [createComputeWeightsArgs()] (with percentile truncation), and `msmWeights` of
#' [createFitOutcomeModelArgs()] (not for `modelType = "km"`). An `initiateWithinModel` or `graceRegime`
#' variant that the pooled logistic hazard model's `intervalLength` does not allow (see
#' [createFitCensoringModelsArgs()]) is left out. Pass `variants` for other combinations.
#'
#' `graceRegime` is not a modelling choice: `"binaryAtG"` and `"naturalGrace"` define different
#' `"initiateWithin"` regimes (see [createFitCensoringModelsArgs()]), so their estimates differ when the
#' start day within the grace period changes later risk. A large difference shows that the choice of
#' regime matters for the data at hand, not that one of them is biased.
#'
#' @param clones                 An object of type [CcwClones] as created by `createClones()`.
#' @param cloneCensorWeightData  The object of type [CloneCensorWeightData] the clones were created
#'                               from.
#' @param fitCensoringModelsArgs An object of type `FitCensoringModelsArgs` as created using
#'                               [createFitCensoringModelsArgs()]: the specified settings.
#' @param computeWeightsArgs     An object of type `ComputeWeightsArgs` as created using
#'                               [createComputeWeightsArgs()].
#' @param fitOutcomeModelArgs    An object of type `FitOutcomeModelArgs` as created using
#'                               [createFitOutcomeModelArgs()].
#' @param variants               `NULL` (the default variants, see details), or a named list of
#'                               variants. Each variant is a list with optional elements
#'                               `fitCensoringModelsArgs`, `computeWeightsArgs` and
#'                               `fitOutcomeModelArgs`, each a named list of the settings to change, for
#'                               example `list("binary, no interactions" = list(fitCensoringModelsArgs =
#'                               list(postGraceModel = "binary", covariateTimeInteraction = "none")))`.
#' @param .censoringModelFits    For internal use ([computeCcwDiagnostics()]): `NULL`, or censoring
#'                               models fitted earlier, which are reused instead of refitted.
#'
#' @return
#' A tibble of class `CcwSensitivity` with columns `variant` (`"specified"` for the specified
#' settings), `changes` (a description of the changed settings) and the columns of the estimates of
#' [fitOutcomeModel()] (`estimand`, `estimate`, `logEstimate`, `seLogEstimate`, `se`, `lb95`, `ub95`,
#' `p`). The attribute `variants` holds the settings of each variant. Each censoring model is fit once
#' and shared by the variants that use it.
#'
#' @example inst/examples/clones.R
#' @examples
#' # With the Cox hazard model, which has no covariate-by-time interactions, there are fewer variants
#' # than with the default pooled logistic model, which keeps the example fast:
#' sensitivity <- computeCcwSensitivity(
#'   clones = clones,
#'   cloneCensorWeightData = data,
#'   fitCensoringModelsArgs = createFitCensoringModelsArgs(
#'     hazardModelType = "cox",
#'     prior = Cyclops::createPrior("laplace", exclude = c(0), variance = 1)
#'   ),
#'   fitOutcomeModelArgs = createFitOutcomeModelArgs(modelType = "km", horizons = 30)
#' )
#' sensitivity
#' close(clones)
#' close(data)
#'
#' @export
computeCcwSensitivity <- function(clones,
                                  cloneCensorWeightData,
                                  fitCensoringModelsArgs = createFitCensoringModelsArgs(),
                                  computeWeightsArgs = createComputeWeightsArgs(),
                                  fitOutcomeModelArgs = createFitOutcomeModelArgs(),
                                  variants = NULL,
                                  .censoringModelFits = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(clones, add = errorMessages)
  .assertCloneCensorWeightData(cloneCensorWeightData, add = errorMessages)
  checkmate::assertR6(fitCensoringModelsArgs, "FitCensoringModelsArgs", add = errorMessages)
  checkmate::assertR6(computeWeightsArgs, "ComputeWeightsArgs", add = errorMessages)
  checkmate::assertR6(fitOutcomeModelArgs, "FitOutcomeModelArgs", add = errorMessages)
  checkmate::assertList(variants, names = "unique", null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  strategies <- attr(clones, "metaData")$strategies
  if (is.null(variants)) {
    variants <- .defaultSensitivityVariants(strategies, fitCensoringModelsArgs, computeWeightsArgs, fitOutcomeModelArgs)
  }
  variants <- c(list(specified = list()), variants)
  base <- list(fitCensoringModelsArgs = fitCensoringModelsArgs, computeWeightsArgs = computeWeightsArgs, fitOutcomeModelArgs = fitOutcomeModelArgs)

  censoringModels <- list()
  # Each censoring model is fit once, also when variants (or computeCcwDiagnostics()) share it:
  fits <- if (is.null(.censoringModelFits)) list() else .censoringModelFits
  warnDay0CovariateWindows <- is.null(.censoringModelFits)
  warnNaturalGrace <- is.null(.censoringModelFits)
  weightedClones <- list()
  on.exit(for (weighted in weightedClones) close(weighted))
  estimates <- list()
  settings <- list()
  for (name in names(variants)) {
    args <- .applyVariant(base, variants[[name]], name)
    censoringKey <- as.character(args$fitCensoringModelsArgs$toJson())
    if (is.null(censoringModels[[censoringKey]])) {
      fitted <- .fitCensoringModelsWithFits(clones, cloneCensorWeightData, args$fitCensoringModelsArgs, fits = fits, warnDay0CovariateWindows = warnDay0CovariateWindows, warnNaturalGrace = warnNaturalGrace)
      warnDay0CovariateWindows <- FALSE
      if (args$fitCensoringModelsArgs$graceRegime == "naturalGrace") {
        warnNaturalGrace <- FALSE
      }
      fits <- fitted$fits
      censoringModels[[censoringKey]] <- fitted$models
    }
    outcomeTimeBasis <- if (args$fitOutcomeModelArgs$modelType == "logistic") args$fitOutcomeModelArgs$timeBasis else NULL
    weightsKey <- paste(censoringKey, as.character(args$computeWeightsArgs$toJson()), if (is.null(outcomeTimeBasis)) "" else as.character(outcomeTimeBasis$toJson()))
    if (is.null(weightedClones[[weightsKey]])) {
      weightedClones[[weightsKey]] <- computeCensoringWeights(clones, censoringModels[[censoringKey]], args$computeWeightsArgs, outcomeTimeBasis)
    }
    model <- fitOutcomeModel(weightedClones[[weightsKey]], args$fitOutcomeModelArgs)
    changes <- .describeVariant(variants[[name]])
    estimates[[name]] <- dplyr::mutate(model$estimates, variant = name, changes = changes, .before = 1)
    settings[[name]] <- tibble(
      variant = name,
      changes = changes,
      postGraceModel = args$fitCensoringModelsArgs$postGraceModel,
      initiateWithinModel = args$fitCensoringModelsArgs$initiateWithinModel,
      graceRegime = args$fitCensoringModelsArgs$graceRegime,
      covariateTimeInteraction = args$fitCensoringModelsArgs$covariateTimeInteraction,
      stabilized = args$computeWeightsArgs$stabilized,
      truncationType = args$computeWeightsArgs$truncation$type,
      truncationBasis = args$computeWeightsArgs$truncation$basis,
      truncateMsmWeights = args$computeWeightsArgs$truncateMsmWeights,
      modelType = args$fitOutcomeModelArgs$modelType,
      msmWeights = args$fitOutcomeModelArgs$msmWeights
    )
  }
  sensitivity <- dplyr::bind_rows(estimates)
  attr(sensitivity, "variants") <- dplyr::bind_rows(settings)
  class(sensitivity) <- c("CcwSensitivity", class(sensitivity))
  return(sensitivity)
}

# One variant per other value of each choice that the analysis uses.
.defaultSensitivityVariants <- function(strategies, fitCensoringModelsArgs, computeWeightsArgs, fitOutcomeModelArgs) {
  types <- vapply(strategies, function(strategy) strategy$type, character(1))
  variants <- list()
  add <- function(argsName, field, value) {
    change <- list()
    change[[argsName]] <- stats::setNames(list(value), field)
    variants[[sprintf("%s = %s", field, format(value))]] <<- change
  }
  if (any(types == "notInitiateWithin")) {
    for (value in setdiff(c("hazard", "binary"), fitCensoringModelsArgs$postGraceModel)) {
      add("fitCensoringModelsArgs", "postGraceModel", value)
    }
  }
  if (any(types == "initiateWithin")) {
    # Choices the hazard model's grid does not allow (.checkInitiateWithinHazardGrid()) are left out:
    for (value in setdiff(c("binary", "hazard"), fitCensoringModelsArgs$initiateWithinModel)) {
      if (.hazardGridAllows(strategies, fitCensoringModelsArgs, initiateWithinModel = value)) {
        add("fitCensoringModelsArgs", "initiateWithinModel", value)
      }
    }
    for (value in setdiff(c("binaryAtG", "naturalGrace"), fitCensoringModelsArgs$graceRegime)) {
      if (.hazardGridAllows(strategies, fitCensoringModelsArgs, graceRegime = value)) {
        add("fitCensoringModelsArgs", "graceRegime", value)
      }
    }
  }
  usesHazardModel <- .usesHazardModel(strategies, fitCensoringModelsArgs$initiateWithinModel)
  if (usesHazardModel && fitCensoringModelsArgs$hazardModelType == "logistic") {
    for (value in setdiff(c("day0", "none", "interval"), fitCensoringModelsArgs$covariateTimeInteraction)) {
      add("fitCensoringModelsArgs", "covariateTimeInteraction", value)
    }
  }
  truncation <- computeWeightsArgs$truncation
  # The Kaplan-Meier analysis does not use the weights on the MSM scale:
  if (truncation$type != "none" && fitOutcomeModelArgs$modelType != "km") {
    add("computeWeightsArgs", "truncateMsmWeights", !computeWeightsArgs$truncateMsmWeights)
  }
  if (truncation$type == "percentile") {
    otherBasis <- setdiff(c("varying", "allClones"), truncation$basis)
    variants[[sprintf("basis = %s", otherBasis)]] <- list(computeWeightsArgs = list(truncation = createWeightTruncation(type = "percentile", upper = truncation$upper, basis = otherBasis)))
  }
  if (fitOutcomeModelArgs$modelType != "km") {
    add("fitOutcomeModelArgs", "msmWeights", setdiff(c("unstabilized", "stabilized"), fitOutcomeModelArgs$msmWeights))
  }
  return(variants)
}

# Does the hazard model's interval grid allow these choices of initiateWithinModel and graceRegime
# (.checkInitiateWithinHazardGrid())?
.hazardGridAllows <- function(strategies,
                              fitCensoringModelsArgs,
                              initiateWithinModel = fitCensoringModelsArgs$initiateWithinModel,
                              graceRegime = fitCensoringModelsArgs$graceRegime) {
  args <- fitCensoringModelsArgs$clone(deep = TRUE)
  args$initiateWithinModel <- initiateWithinModel
  args$graceRegime <- graceRegime
  return(!inherits(tryCatch(.checkInitiateWithinHazardGrid(strategies, args), error = function(e) e), "error"))
}

# Does an analysis use an initiation hazard model (whose covariateTimeInteraction then matters)? Every
# strategy but initiateWithin does, and initiateWithin with initiateWithinModel = "hazard" or with
# graceStart > 0 (the days before graceStart).
.usesHazardModel <- function(strategies, initiateWithinModel) {
  return(any(vapply(strategies, function(strategy) {
    strategy$type != "initiateWithin" || strategy$graceStart > 0
  }, logical(1))) || "hazard" %in% initiateWithinModel)
}

# Copies of the settings objects in `base` with the changes of `variant`.
.applyVariant <- function(base, variant, name) {
  unknown <- setdiff(names(variant), names(base))
  if (length(unknown) > 0) {
    stop(sprintf("Variant '%s' changes unknown settings: %s", name, paste(unknown, collapse = ", ")))
  }
  args <- lapply(base, function(settings) settings$clone(deep = TRUE))
  for (argsName in names(variant)) {
    for (field in names(variant[[argsName]])) {
      if (!field %in% names(args[[argsName]]) || is.function(args[[argsName]][[field]])) {
        stop(sprintf("Variant '%s': %s has no setting '%s'", name, argsName, field))
      }
      args[[argsName]][[field]] <- variant[[argsName]][[field]]
    }
    args[[argsName]]$validate()
  }
  return(args)
}

.describeVariant <- function(variant) {
  changes <- unlist(lapply(names(variant), function(argsName) {
    vapply(names(variant[[argsName]]), function(field) {
      value <- variant[[argsName]][[field]]
      if (inherits(value, "WeightTruncation")) {
        value <- sprintf("%s(%s, basis = %s)", value$type, format(value$upper), value$basis)
      }
      sprintf("%s = %s", field, paste(format(value), collapse = ", "))
    }, character(1))
  }))
  return(if (length(changes) == 0) "" else paste(changes, collapse = "; "))
}

#' @export
print.CcwSensitivity <- function(x, ...) {
  writeLines("Estimates under the specified settings and the variants:")
  data <- as.data.frame(x)
  cell <- ifelse(is.na(data$lb95), sprintf("%.3f", data$estimate), sprintf("%.3f (%.3f, %.3f)", data$estimate, data$lb95, data$ub95))
  variants <- unique(data$variant)
  estimands <- unique(data$estimand)
  table <- matrix("", nrow = length(estimands), ncol = length(variants), dimnames = list(estimands, variants))
  table[cbind(match(data$estimand, estimands), match(data$variant, variants))] <- cell
  print(noquote(table))
  invisible(x)
}

#' Plot the estimates under alternative design choices
#'
#' @description
#' Plots the estimates and 95% confidence intervals of [computeCcwSensitivity()] per estimand, one row
#' per variant. Requires the `ggplot2` package.
#'
#' @param sensitivity A tibble as created by [computeCcwSensitivity()].
#' @param estimands   The estimands to plot. If `NULL`, all.
#' @param fileName    Name of the file where the plot should be saved, for example 'plot.png'. See
#'                    [ggplot2::ggsave()] for supported file formats. If `NULL`, the plot is not saved.
#'
#' @return
#' A ggplot object. Use the [ggplot2::ggsave()] function to save to file in a different format.
#'
#' @examples
#' sensitivity <- data.frame(
#'   variant = c("specified", "postGraceModel = binary"),
#'   estimand = "rd_60",
#'   estimate = c(-0.020, -0.031),
#'   lb95 = c(-0.035, -0.046),
#'   ub95 = c(-0.005, -0.016)
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plotCcwSensitivity(sensitivity)
#' }
#'
#' @export
plotCcwSensitivity <- function(sensitivity, estimands = NULL, fileName = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(sensitivity, add = errorMessages)
  checkmate::assertNames(colnames(sensitivity), must.include = c("variant", "estimand", "estimate", "lb95", "ub95"), add = errorMessages)
  checkmate::assertCharacter(estimands, min.len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::assertCharacter(fileName, len = 1, null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required to plot the sensitivity analysis. Please install it.")
  }
  data <- as.data.frame(sensitivity)
  if (!is.null(estimands)) {
    data <- data[data$estimand %in% estimands, ]
  }
  data$variant <- factor(data$variant, levels = rev(unique(data$variant)))
  data$estimand <- factor(data$estimand, levels = unique(data$estimand))
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = .data$estimate, y = .data$variant)) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lb95, xmax = .data$ub95), width = 0.2, orientation = "y", color = "#2a78d6", na.rm = TRUE) +
    ggplot2::geom_point(size = 2, color = "#2a78d6") +
    ggplot2::facet_wrap(~estimand, scales = "free_x") +
    ggplot2::labs(x = "Estimate (95% confidence interval)", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  if (!is.null(fileName)) {
    ggplot2::ggsave(fileName, plot, width = 7, height = 4.5, dpi = 400)
  }
  return(plot)
}

# All diagnostics -----------------------------------------------------------------------------------

#' Compute the diagnostics for choosing the design options of a clone-censor-weight analysis
#'
#' @description
#' Computes, for the specified settings and for the alternative choices of `postGraceModel`,
#' `initiateWithinModel`, `graceRegime` and `covariateTimeInteraction`, the diagnostics that help choose between them for the data at hand:
#' calibration of the censoring models, the grace-period selection check, the weights on both scales,
#' covariate balance after the grace period, positivity and, optionally, the estimates under each
#' choice.
#'
#' @details
#' The choices compared are `postGraceModel` `"hazard"` and `"binary"` (when a strategy is
#' `"notInitiateWithin"`), crossed with `covariateTimeInteraction` `"none"`, `"day0"` and the
#' specified one (with the pooled logistic hazard model, when a strategy uses it; the Cox model has
#' none), crossed with `initiateWithinModel` `"binary"` and `"hazard"` and with `graceRegime`
#' `"binaryAtG"` and `"naturalGrace"` (when a strategy is `"initiateWithin"`; combinations that the
#' hazard model's interval grid does not allow are left out). `graceRegime` chooses the regime (the
#' estimand) rather than a model, so its diagnostics describe the weights of each regime and do not
#' rank them. Choices that give the same censoring models (for example the interaction
#' settings when no strategy uses the hazard model) are compared once. Each censoring model is fit once
#' and shared by the choices that use it (for example the hazard model through `G` serves both
#' `postGraceModel` choices, and both strategies with `initiateWithinModel = "hazard"`), and the models
#' in `censoringModels` are reused, so only the models they lack are fit. For each choice the result
#' holds:
#'
#' - `calibration`: [computeCensoringModelCalibration()], for the binary model at `G` and for the
#'   hazard model of each interaction setting, by group of the prediction, in the large, and by the
#'   covariates in `strataCovariateIds`;
#' - `gracePeriodSelection`: [computeGracePeriodSelection()] for the specified settings, which tells
#'   whether `postGraceModel = "binary"` and `initiateWithinModel = "hazard"` are valid;
#' - `weightSummary`: [computeWeightSummary()] for each choice (column `variant`);
#' - `balance`: per choice, the largest absolute standardized difference of mean and the number of
#'   covariates above `diagnosticThresholds$sdmThreshold`, before and after weighting, of
#'   [computeCovariateBalance()] at day `tau + 1`;
#' - `positivity`: [computePositivity()] for the specified interaction setting, with the binary model
#'   at `G` and the hazard model of each strategy with a grace period;
#' - `sensitivity`: [computeCcwSensitivity()] with its default variants, only when `sensitivity = TRUE`,
#'   because it refits the outcome model under every variant.
#'
#' Read the balance with the grace-period selection check. The clones compared at `tau + 1` are the
#' survivors of two regimes. When starting the treatment changes the risk during the grace period, the
#' persons still followed differ between the strategies, and the valid `"hazard"` weights, which
#' reweight to all persons at day 0, need not balance them, while the `"binary"` weights, fitted on the
#' survivors, do (in the simulation of [createFitCensoringModelsArgs()], the largest absolute
#' standardized difference was 0.27 with `"hazard"` and 0.002 with `"binary"`, whose risks were the
#' biased ones). So when the selection check flags, balance cannot choose `postGraceModel`.
#'
#' The printed result flags what exceeds `diagnosticThresholds`.
#'
#' @param clones                      An object of type [CcwClones] as created by `createClones()`.
#' @param cloneCensorWeightData       The object of type [CloneCensorWeightData] the clones were
#'                                    created from.
#' @param fitCensoringModelsArgs      An object of type `FitCensoringModelsArgs` as created using
#'                                    [createFitCensoringModelsArgs()]: the specified settings.
#' @param computeWeightsArgs          An object of type `ComputeWeightsArgs` as created using
#'                                    [createComputeWeightsArgs()].
#' @param computeCovariateBalanceArgs An object of type `ComputeCovariateBalanceArgs` as created using
#'                                    [CohortMethod::createComputeCovariateBalanceArgs()].
#' @param diagnosticThresholds        An object of type `CcwDiagnosticThresholds` as created using
#'                                    [createCcwDiagnosticThresholds()].
#' @param numberOfGroups              The number of groups (quantiles) of the predictions.
#' @param strataCovariateIds          IDs of covariates by which the calibration is also reported.
#' @param sensitivity                 Also refit the outcome model under alternative choices (see
#'                                    [computeCcwSensitivity()])?
#' @param fitOutcomeModelArgs         An object of type `FitOutcomeModelArgs` as created using
#'                                    [createFitOutcomeModelArgs()], used with `sensitivity = TRUE`.
#' @param censoringModels             `NULL`, or an object of class `CcwCensoringModels` as created by
#'                                    [fitCensoringModels()] on the same clones (without the internal
#'                                    `.fixedPriorVariances`): its models are reused, for the specified
#'                                    settings and for every other choice that uses them, instead of
#'                                    refitted. When `fitCensoringModelsArgs` is not given, the
#'                                    specified settings are those the models were fit with. When they
#'                                    were fit with settings other than `fitCensoringModelsArgs` (besides
#'                                    the design choices compared, including `graceRegime`), the models
#'                                    that depend on those
#'                                    settings are fit again. For example, the binary model at G is
#'                                    reused when only the hazard model's settings (`timeBasis` or
#'                                    `hazardModelType`) differ, and the hazard models when only
#'                                    `maxCohortSizeForFitting` differs.
#'
#' @return
#' An object of class `CcwDiagnostics`: a list with `variants` (a tibble with `variant`,
#' `postGraceModel`, `covariateTimeInteraction`, `initiateWithinModel`, `graceRegime` and
#' `specified`; the `variant` label names `graceRegime` only when a strategy is `"initiateWithin"`),
#' `calibration`,
#' `gracePeriodSelection`, `weightSummary`, `balance`, `positivity`, `sensitivity` (`NULL` unless
#' requested) and `thresholds`.
#'
#' @example inst/examples/clones.R
#' @examples
#' # With the Cox hazard model, which has no covariate-by-time interactions, 8 variants are compared
#' # rather than the 16 of the default pooled logistic model, which keeps the example fast:
#' censoringModels <- fitCensoringModels(
#'   clones = clones,
#'   cloneCensorWeightData = data,
#'   fitCensoringModelsArgs = createFitCensoringModelsArgs(
#'     hazardModelType = "cox",
#'     prior = Cyclops::createPrior("laplace", exclude = c(0), variance = 1)
#'   )
#' )
#' # Reuses the fitted models, and fits only those of the other choices:
#' diagnostics <- computeCcwDiagnostics(
#'   clones = clones,
#'   cloneCensorWeightData = data,
#'   numberOfGroups = 5,
#'   censoringModels = censoringModels
#' )
#' diagnostics
#' close(clones)
#' close(data)
#'
#' @export
computeCcwDiagnostics <- function(clones,
                                  cloneCensorWeightData,
                                  fitCensoringModelsArgs = createFitCensoringModelsArgs(),
                                  computeWeightsArgs = createComputeWeightsArgs(),
                                  computeCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs(),
                                  diagnosticThresholds = createCcwDiagnosticThresholds(),
                                  numberOfGroups = 10,
                                  strataCovariateIds = NULL,
                                  sensitivity = FALSE,
                                  fitOutcomeModelArgs = createFitOutcomeModelArgs(),
                                  censoringModels = NULL) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(clones, add = errorMessages)
  .assertCloneCensorWeightData(cloneCensorWeightData, add = errorMessages)
  checkmate::assertR6(fitCensoringModelsArgs, "FitCensoringModelsArgs", add = errorMessages)
  checkmate::assertR6(computeWeightsArgs, "ComputeWeightsArgs", add = errorMessages)
  checkmate::assertR6(diagnosticThresholds, "CcwDiagnosticThresholds", add = errorMessages)
  checkmate::assertInt(numberOfGroups, lower = 1, add = errorMessages)
  checkmate::assertFlag(sensitivity, add = errorMessages)
  checkmate::assertR6(fitOutcomeModelArgs, "FitOutcomeModelArgs", add = errorMessages)
  checkmate::assertClass(censoringModels, "CcwCensoringModels", null.ok = TRUE, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)
  # The fitted censoring models, reused by every choice that uses them:
  fits <- list()
  if (!is.null(censoringModels)) {
    .checkModelsMatchClones(censoringModels, attr(clones, "metaData"))
    modelsArgs <- attr(censoringModels, "metaData")$fitCensoringModelsArgs
    if (missing(fitCensoringModelsArgs)) {
      fitCensoringModelsArgs <- FitCensoringModelsArgs$new(untypedList = modelsArgs)
    } else if (.withoutDesignChoices(modelsArgs) != .withoutDesignChoices(fitCensoringModelsArgs$toList())) {
      message("The censoring models given were fit with other settings than fitCensoringModelsArgs (besides postGraceModel, initiateWithinModel, graceRegime and covariateTimeInteraction), so the models that depend on them are fit again")
    }
    fits <- .censoringFitsFromModels(censoringModels)
  }
  strategies <- attr(clones, "metaData")$strategies
  types <- vapply(strategies, function(strategy) strategy$type, character(1))

  postGraceModels <- if (any(types == "notInitiateWithin")) c("hazard", "binary") else fitCensoringModelsArgs$postGraceModel
  initiateWithinModels <- if (any(types == "initiateWithin")) c("binary", "hazard") else fitCensoringModelsArgs$initiateWithinModel
  graceRegimes <- if (any(types == "initiateWithin")) c("binaryAtG", "naturalGrace") else fitCensoringModelsArgs$graceRegime
  interactions <- fitCensoringModelsArgs$covariateTimeInteraction
  usesHazardModel <- .usesHazardModel(strategies, initiateWithinModels)
  if (fitCensoringModelsArgs$hazardModelType == "logistic" && usesHazardModel) {
    interactions <- unique(c(interactions, "none", "day0"))
  }
  variants <- expand.grid(
    postGraceModel = postGraceModels,
    covariateTimeInteraction = interactions,
    initiateWithinModel = initiateWithinModels,
    graceRegime = graceRegimes,
    stringsAsFactors = FALSE
  )
  # Combinations the hazard model's grid does not allow (.checkInitiateWithinHazardGrid()) are left out:
  allowed <- vapply(seq_len(nrow(variants)), function(i) {
    .hazardGridAllows(strategies, fitCensoringModelsArgs, variants$initiateWithinModel[i], variants$graceRegime[i])
  }, logical(1))
  variants <- variants[allowed, ]
  variants$variant <- sprintf(
    "postGraceModel = %s, covariateTimeInteraction = %s, initiateWithinModel = %s",
    variants$postGraceModel, variants$covariateTimeInteraction, variants$initiateWithinModel
  )
  # The regime is named only when the analysis has an initiateWithin strategy, so the labels of other
  # analyses stay as they were:
  if (any(types == "initiateWithin")) {
    variants$variant <- sprintf("%s, graceRegime = %s", variants$variant, variants$graceRegime)
  }
  variants$specified <- variants$postGraceModel == fitCensoringModelsArgs$postGraceModel &
    variants$covariateTimeInteraction == fitCensoringModelsArgs$covariateTimeInteraction &
    variants$initiateWithinModel == fitCensoringModelsArgs$initiateWithinModel &
    variants$graceRegime == fitCensoringModelsArgs$graceRegime
  variants <- dplyr::as_tibble(variants[order(!variants$specified), c("variant", "postGraceModel", "covariateTimeInteraction", "initiateWithinModel", "graceRegime", "specified")])

  # The censoring models of each choice, each model fit once (or taken from censoringModels):
  variantModels <- list()
  for (i in seq_len(nrow(variants))) {
    args <- fitCensoringModelsArgs$clone(deep = TRUE)
    args$postGraceModel <- variants$postGraceModel[i]
    args$covariateTimeInteraction <- variants$covariateTimeInteraction[i]
    args$initiateWithinModel <- variants$initiateWithinModel[i]
    args$graceRegime <- variants$graceRegime[i]
    # Warn about day-0 covariate windows once, and about a natural-grace regime without starts on day G
    # once (with the first variant of that regime):
    fitted <- .fitCensoringModelsWithFits(
      clones, cloneCensorWeightData, args,
      fits = fits,
      warnDay0CovariateWindows = i == 1,
      warnNaturalGrace = isTRUE(i == match("naturalGrace", variants$graceRegime))
    )
    fits <- fitted$fits
    variantModels[[i]] <- fitted$models
  }
  # Choices with the same censoring models are compared once (the specified one, or the first):
  duplicate <- vapply(seq_along(variantModels), function(i) {
    any(vapply(seq_len(i - 1), function(j) .sameCensoringModels(variantModels[[i]], variantModels[[j]]), logical(1)))
  }, logical(1))
  variants <- variants[!duplicate, ]
  variantModels <- variantModels[!duplicate]
  specifiedModels <- variantModels[[1]]

  calibration <- dplyr::bind_rows(lapply(variantModels, function(models) {
    computeCensoringModelCalibration(clones, models, numberOfGroups, cloneCensorWeightData, strataCovariateIds)
  }))
  calibration <- calibration[!duplicated(calibration[, c("strategyId", "model", "covariateTimeInteraction", "graceRegime", "grouping", "covariateId", "group")]), ]

  balance <- list()
  weightSummaries <- list()
  for (i in seq_len(nrow(variants))) {
    weightedClones <- computeCensoringWeights(clones, variantModels[[i]], computeWeightsArgs)
    tryCatch(
      {
        covariateBalance <- computeCovariateBalance(weightedClones, cloneCensorWeightData, computeCovariateBalanceArgs)
        balance[[i]] <- .summarizeBalance(covariateBalance, diagnosticThresholds$sdmThreshold) |>
          dplyr::mutate(variant = variants$variant[i], .before = 1)
        weightSummaries[[i]] <- computeWeightSummary(weightedClones) |>
          dplyr::as_tibble() |>
          dplyr::mutate(variant = variants$variant[i], .before = 1)
      },
      finally = close(weightedClones)
    )
  }

  # Positivity of the specified interaction setting, with the binary model at G of notInitiateWithin's
  # weights after G when there is one, and both models of initiateWithin:
  postGraceModel <- if ("binary" %in% postGraceModels) "binary" else fitCensoringModelsArgs$postGraceModel
  positivityModels <- variantModels[variants$covariateTimeInteraction == fitCensoringModelsArgs$covariateTimeInteraction & variants$postGraceModel == postGraceModel]
  if (length(positivityModels) == 0) {
    positivityModels <- list(specifiedModels)
  }
  positivity <- .combinePositivity(lapply(positivityModels, function(models) {
    computePositivity(clones, models, diagnosticThresholds$positivityThreshold)
  }))

  diagnostics <- list(
    variants = variants,
    calibration = calibration,
    gracePeriodSelection = computeGracePeriodSelection(clones, specifiedModels, numberOfGroups),
    weightSummary = dplyr::bind_rows(weightSummaries),
    balance = dplyr::bind_rows(balance),
    positivity = positivity,
    sensitivity = if (sensitivity) computeCcwSensitivity(clones, cloneCensorWeightData, fitCensoringModelsArgs, computeWeightsArgs, fitOutcomeModelArgs, .censoringModelFits = fits) else NULL,
    thresholds = diagnosticThresholds$toList()
  )
  class(diagnostics) <- "CcwDiagnostics"
  return(diagnostics)
}

# The settings of the censoring models (FitCensoringModelsArgs$toList()) without the design choices that
# computeCcwDiagnostics() compares, as JSON.
.withoutDesignChoices <- function(argsList) {
  argsList$postGraceModel <- NULL
  argsList$initiateWithinModel <- NULL
  argsList$graceRegime <- NULL
  argsList$covariateTimeInteraction <- NULL
  return(.settingsJson(argsList))
}

# Do two CcwCensoringModels hold the same models (whatever settings they record)?
.sameCensoringModels <- function(a, b) {
  a <- unclass(a)
  b <- unclass(b)
  attr(a, "metaData") <- NULL
  attr(b, "metaData") <- NULL
  return(identical(a, b))
}

# One CcwPositivity from several, each model of each strategy once (in strategy order).
.combinePositivity <- function(positivities) {
  positivity <- dplyr::bind_rows(lapply(positivities, dplyr::as_tibble))
  keep <- !duplicated(positivity[, c("strategyId", "model", "graceRegime", "day")])
  positivity <- positivity[keep, ]
  positivity <- positivity[order(positivity$strategyId), ]
  predictions <- dplyr::bind_rows(lapply(positivities, function(x) attr(x, "predictions")))
  predictions <- predictions[!duplicated(predictions[, c("strategyId", "model", "graceRegime", "day", "rowId")]), ]
  attr(positivity, "predictions") <- predictions[order(predictions$strategyId), ]
  attr(positivity, "threshold") <- attr(positivities[[1]], "threshold")
  attr(positivity, "strategyLabels") <- attr(positivities[[1]], "strategyLabels")
  class(positivity) <- c("CcwPositivity", class(positivity))
  return(positivity)
}

.summarizeBalance <- function(balance, threshold) {
  before <- abs(balance$beforeMatchingStdDiff)
  after <- abs(balance$afterMatchingStdDiff)
  return(tibble(
    covariates = nrow(balance),
    maxAbsStdDiffBefore = if (all(is.na(before))) NA_real_ else max(before, na.rm = TRUE),
    maxAbsStdDiffAfter = if (all(is.na(after))) NA_real_ else max(after, na.rm = TRUE),
    aboveThresholdBefore = sum(before > threshold, na.rm = TRUE),
    aboveThresholdAfter = sum(after > threshold, na.rm = TRUE)
  ))
}

#' @export
print.CcwDiagnostics <- function(x, ...) {
  thresholds <- x$thresholds
  writeLines("Clone-censor-weight diagnostics")
  writeLines("")
  writeLines(sprintf("Covariate balance at tau + 1 (SDM threshold %s):", format(thresholds$sdmThreshold)))
  balance <- as.data.frame(x$balance)
  balance$specified <- x$variants$specified[match(balance$variant, x$variants$variant)]
  print(balance, row.names = FALSE, digits = 3)
  writeLines("")
  .printSelectionFlags(x$gracePeriodSelection)
  writeLines("")
  calibration <- as.data.frame(x$calibration)
  overall <- calibration[calibration$grouping == "overall", c("strategyId", "model", "covariateTimeInteraction", "graceRegime", "persons", "observed", "predicted")]
  groups <- calibration[calibration$grouping == "prediction", ]
  key <- function(d) paste(d$strategyId, d$model, d$covariateTimeInteraction, d$graceRegime)
  overall$largestGroupGap <- as.vector(tapply(abs(groups$observed - groups$predicted), key(groups), max)[key(overall)])
  writeLines("Calibration at the end of the grace period (in the large, and the largest gap over groups of the prediction):")
  print(overall, row.names = FALSE, digits = 3)
  writeLines("")
  weights <- as.data.frame(x$weightSummary)
  weights <- weights[weights$truncation == "after", ]
  worst <- do.call(rbind, lapply(split(weights, list(weights$variant, weights$strategyId, weights$scale), drop = TRUE), function(d) {
    data.frame(variant = d$variant[1], strategyId = d$strategyId[1], scale = d$scale[1], maxWeight = max(d$maxWeight), minEffectiveSampleShare = min(d$effectiveSampleShare))
  }))
  worst$flag <- ifelse(worst$maxWeight > thresholds$maxWeightThreshold | worst$minEffectiveSampleShare < thresholds$essFractionThreshold, "FLAG", "")
  writeLines(sprintf("Weights after truncation (flagged above %s or with an effective sample share below %s):", format(thresholds$maxWeightThreshold), format(thresholds$essFractionThreshold)))
  print(worst[order(worst$variant, worst$strategyId, worst$scale), ], row.names = FALSE, digits = 3)
  writeLines("")
  positivity <- as.data.frame(x$positivity)[, c("strategyId", "model", "graceRegime", "day", "persons", "min", "median", "max", "shareExtreme")]
  writeLines(sprintf("Positivity (share of persons whose probability of following the strategy is below %s):", format(thresholds$positivityThreshold)))
  print(positivity, row.names = FALSE, digits = 3)
  if (!is.null(x$sensitivity)) {
    writeLines("")
    print(x$sensitivity)
  }
  invisible(x)
}

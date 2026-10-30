#' Compute weighted Kaplan-Meier survival estimates
#'
#' Fits a weighted KM estimator using the counting-process formulation
#' (`Surv(tstart, tstop, event)`) so that interval-format person-time data is
#' handled correctly.  Each clone-interval contributes risk only from `tstart`
#' to `tstop`, avoiding the immortal-time inflation that would occur if every
#' interval row were treated as entering at time 0.
#'
#' @param personTime Weighted person-time data frame from
#'   [createCensoringWeights()].  Must contain `tstart`, `tstop`, `event`,
#'   `arm`, and `weight`.
#' @return A data frame with columns `arm`, `time`, `survival`, `lower`,
#'   `upper` (95 % CI), `nRisk`, and `nEvent`.  A row with `time = 0` and
#'   `survival = 1` is prepended for each arm so that the step function
#'   originates at the y-axis.
#' @export
computeKaplanMeier <- function(personTime) {
  checkmate::assertDataFrame(personTime)
  checkmate::assertSubset(c("tstart", "tstop", "event", "arm", "weight"),
                          names(personTime))

  # Extract weights before the survfit call to avoid name shadowing: the data
  # frame has a column named "personTime", so personTime$weight evaluated inside
  # survfit's model.frame would resolve "personTime" to that column, not this arg.
  wts <- personTime[["weight"]]
  fit <- survival::survfit(
    survival::Surv(tstart, tstop, event) ~ arm,
    data    = personTime,
    weights = wts
  )

  s <- summary(fit)

  km <- if (length(s$time) > 0) {
    data.frame(
      arm      = sub("^[^=]+=", "", as.character(s$strata)),
      time     = s$time,
      survival = s$surv,
      lower    = s$lower,
      upper    = s$upper,
      nRisk    = s$n.risk,
      nEvent   = s$n.event,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(arm = character(), time = numeric(), survival = numeric(),
               lower = numeric(), upper = numeric(),
               nRisk = numeric(), nEvent = numeric(),
               stringsAsFactors = FALSE)
  }

  # Prepend time=0 (survival=1) for each arm so curves start at the origin
  armLevels <- sort(unique(c(km$arm, sub("^[^=]+=", "", names(fit$strata)))))
  origin <- data.frame(
    arm = armLevels, time = 0, survival = 1, lower = 1, upper = 1,
    nRisk = NA_real_, nEvent = 0, stringsAsFactors = FALSE
  )

  out <- rbind(origin, km)
  out[order(out$arm, out$time), ]
}

#' Plot weighted Kaplan-Meier survival curves
#'
#' Calls [computeKaplanMeier()] and renders a step-function plot with a
#' shaded 95 % confidence band.
#'
#' @param personTime Weighted person-time data frame from
#'   [createCensoringWeights()].
#' @return A `ggplot2` object.
#' @export
plotKaplanMeier <- function(personTime) {
  checkmate::assertDataFrame(personTime)
  kmData <- computeKaplanMeier(personTime)

  ggplot2::ggplot(kmData,
    ggplot2::aes(x = time, y = survival, colour = arm, fill = arm)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         alpha = 0.15, colour = NA) +
    ggplot2::geom_step() +
    ggplot2::scale_y_continuous(limits = c(NA, 1)) +
    ggplot2::labs(
      x      = "Days since trial start",
      y      = "Survival probability",
      colour = NULL,
      fill   = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

#' Calibrate estimates using negative control outcomes
#'
#' Fits a systematic error null distribution from negative control outcomes and
#' uses it to produce calibrated p-values and 95 % confidence intervals for a
#' primary outcome estimate.
#'
#' The null distribution characterizes the residual bias and imprecision
#' present in the study design.  By default, the frequentist
#' [EmpiricalCalibration::fitNull()] estimator is used; for publication-quality
#' analyses, pre-compute the null with
#' [EmpiricalCalibration::fitMcmcNull()] and pass it via the `null` argument.
#'
#' @param outcomeModel A `CloneCensorWeightOutcomeModel` from
#'   [fitOutcomeModel()].
#' @param negativeControlResults A data frame with columns `logRr` and
#'   `seLogRr` — one row per negative control outcome.  Exactly one of
#'   `negativeControlResults` and `null` must be supplied.
#' @param null A pre-computed null distribution object from
#'   [EmpiricalCalibration::fitNull()] or
#'   [EmpiricalCalibration::fitMcmcNull()]. Takes precedence over
#'   `negativeControlResults` when both are supplied.
#'
#' @return A data frame identical to `outcomeModel$summary` with three
#'   additional columns: `calibratedP`, `calibratedLower95`, and
#'   `calibratedUpper95`.
#' @export
calibrateEstimates <- function(outcomeModel,
                               negativeControlResults = NULL,
                               null                   = NULL) {
  checkmate::assertClass(outcomeModel, "CloneCensorWeightOutcomeModel")
  if (is.null(null) && is.null(negativeControlResults)) {
    stop("Exactly one of 'negativeControlResults' or 'null' must be supplied.",
         call. = FALSE)
  }
  if (!is.null(negativeControlResults)) {
    checkmate::assertDataFrame(negativeControlResults)
    checkmate::assertSubset(c("logRr", "seLogRr"), names(negativeControlResults))
  }

  if (is.null(null)) {
    null <- EmpiricalCalibration::fitNull(
      logRr   = negativeControlResults$logRr,
      seLogRr = negativeControlResults$seLogRr
    )
  }

  summary <- outcomeModel$summary
  logRr   <- log(summary$estimate)
  se      <- (log(summary$upper95) - log(summary$lower95)) / (2 * stats::qnorm(0.975))

  calP <- EmpiricalCalibration::calibrateP(
    null    = null,
    logRr   = logRr,
    seLogRr = se
  )
  # calibrateP returns a scalar for fitNull and a data.frame for fitMcmcNull
  pValue <- if (is.data.frame(calP)) calP$p else calP

  em    <- EmpiricalCalibration::convertNullToErrorModel(null)
  calCI <- EmpiricalCalibration::calibrateConfidenceInterval(
    logRr   = logRr,
    seLogRr = se,
    model   = em
  )

  summary$calibratedP        <- pValue
  summary$calibratedLower95  <- exp(calCI$logLb95Rr)
  summary$calibratedUpper95  <- exp(calCI$logUb95Rr)
  summary
}

#' Plot the empirical calibration effect
#'
#' Displays negative control estimates as a scatter plot of log ratio vs.
#' standard error, with a shaded funnel showing the 95 % expected region under
#' the fitted null distribution.  A primary estimate can be overlaid as a
#' distinct symbol.
#'
#' @param negativeControlResults A data frame with columns `logRr` and
#'   `seLogRr`.
#' @param primaryEstimate An optional single-row data frame with columns
#'   `logRr` and `seLogRr` representing the primary outcome estimate.
#' @param null An optional pre-computed null distribution.  When `NULL`, the
#'   null is fitted internally with [EmpiricalCalibration::fitNull()].
#'
#' @return A `ggplot2` object.
#' @export
plotCalibration <- function(negativeControlResults,
                            primaryEstimate = NULL,
                            null            = NULL) {
  checkmate::assertDataFrame(negativeControlResults)
  checkmate::assertSubset(c("logRr", "seLogRr"), names(negativeControlResults))

  if (is.null(null)) {
    null <- EmpiricalCalibration::fitNull(
      logRr   = negativeControlResults$logRr,
      seLogRr = negativeControlResults$seLogRr
    )
  }

  nullMean <- null["mean"]
  nullSd   <- null["sd"]

  # Funnel polygon: 95 % band of null distribution as a function of SE.
  # SE runs from 0 (tip) to maxSe (widest).
  maxSe    <- max(negativeControlResults$seLogRr, na.rm = TRUE) * 1.1
  seGrid   <- seq(0, maxSe, length.out = 200)
  lower    <- nullMean - stats::qnorm(0.975) * sqrt(nullSd^2 + seGrid^2)
  upper    <- nullMean + stats::qnorm(0.975) * sqrt(nullSd^2 + seGrid^2)

  funnelDf <- data.frame(
    logRr = c(lower, rev(upper)),
    se    = c(seGrid, rev(seGrid))
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = funnelDf,
      ggplot2::aes(x = logRr, y = se),
      fill = "grey85", colour = NA
    ) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_point(
      data = negativeControlResults,
      ggplot2::aes(x = logRr, y = seLogRr),
      colour = "steelblue", alpha = 0.7, size = 1.8
    ) +
    ggplot2::scale_y_reverse(limits = c(maxSe, 0)) +
    ggplot2::labs(
      x = "Log ratio",
      y = "Standard error"
    ) +
    ggplot2::theme_minimal()

  if (!is.null(primaryEstimate)) {
    checkmate::assertDataFrame(primaryEstimate)
    checkmate::assertSubset(c("logRr", "seLogRr"), names(primaryEstimate))
    p <- p + ggplot2::geom_point(
      data   = primaryEstimate,
      ggplot2::aes(x = logRr, y = seLogRr),
      colour = "firebrick", size = 3.5, shape = 17
    )
  }

  p
}

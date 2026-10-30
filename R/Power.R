#' Compute the minimum detectable relative risk
#'
#' Estimates the smallest hazard ratio detectable at the requested power level
#' given the observed number of events in each arm.  IPCW-weighted event counts
#' are adjusted downward using Kish's effective-sample-size formula to account
#' for variance inflation from heterogeneous weights.
#'
#' The formula follows Schoenfeld's (1981) approximation for the log-rank test:
#' \deqn{
#'   \mathrm{MDRR} = \exp\!\left(
#'     \frac{z_{1-\alpha/2} + z_{\mathrm{power}}}{\sqrt{d_1^\star d_2^\star / (d_1^\star + d_2^\star)}}
#'   \right)
#' }
#' where \eqn{d_k^\star} is the Kish effective event count for arm \eqn{k}.
#'
#' @param personTime Weighted person-time data frame from
#'   [createCensoringWeights()].  Must contain `arm`, `event`, and `weight`.
#' @param alpha Two-sided significance level.  Default `0.05`.
#' @param power Desired statistical power.  Default `0.80`.
#'
#' @return A single-row data frame with columns `arm1`, `arm2`, `alpha`,
#'   `power`, `nEventsArm1`, `nEventsArm2`, `nEventsEffArm1`, `nEventsEffArm2`,
#'   and `mdrr`.  `nEvents*` is the raw observed count; `nEventsEff*` is the
#'   Kish-adjusted effective count.  `mdrr` is `Inf` when either arm has zero
#'   events.
#' @export
computeMdrr <- function(personTime, alpha = 0.05, power = 0.80) {
  checkmate::assertDataFrame(personTime)
  checkmate::assertSubset(c("arm", "event", "weight"), names(personTime))
  checkmate::assertNumber(alpha, lower = 0, upper = 1)
  checkmate::assertNumber(power, lower = 0, upper = 1)

  armLevels <- sort(unique(as.character(personTime$arm)))
  if (length(armLevels) != 2) {
    stop("personTime must contain exactly two arm levels.", call. = FALSE)
  }

  .armStats <- function(armLabel) {
    sub       <- personTime[as.character(personTime$arm) == armLabel & personTime$event == 1, ]
    nEvents   <- nrow(sub)
    nEff      <- if (nEvents == 0) 0 else sum(sub$weight)^2 / sum(sub$weight^2)
    list(nEvents = nEvents, nEff = nEff)
  }

  s1 <- .armStats(armLevels[1])
  s2 <- .armStats(armLevels[2])

  d1 <- s1$nEff
  d2 <- s2$nEff
  mdrr <- if (d1 == 0 || d2 == 0) {
    Inf
  } else {
    z <- stats::qnorm(1 - alpha / 2) + stats::qnorm(power)
    exp(z / sqrt(d1 * d2 / (d1 + d2)))
  }

  data.frame(
    arm1            = armLevels[1],
    arm2            = armLevels[2],
    alpha           = alpha,
    power           = power,
    nEventsArm1     = s1$nEvents,
    nEventsArm2     = s2$nEvents,
    nEventsEffArm1  = d1,
    nEventsEffArm2  = d2,
    mdrr            = mdrr,
    stringsAsFactors = FALSE
  )
}

#' Plot minimum detectable relative risk
#'
#' Displays the MDRR from one or more [computeMdrr()] results as a dot plot on
#' a linear scale, with a dashed reference line at `MDRR = 2`.  Useful for
#' comparing power across multiple phenotype definitions or databases when the
#' rows of `mdrrTable` come from different analyses.
#'
#' @param mdrrTable A data frame — either the single-row output of
#'   [computeMdrr()] or several such rows stacked with `rbind()`.  Must contain
#'   a `mdrr` column.  An optional `analysisId` column labels each row on the
#'   y-axis; if absent, rows are numbered 1, 2, …
#' @return A `ggplot2` object.
#' @export
plotMdrr <- function(mdrrTable) {
  checkmate::assertDataFrame(mdrrTable)
  checkmate::assertSubset("mdrr", names(mdrrTable))

  if (!"analysisId" %in% names(mdrrTable)) {
    mdrrTable$analysisId <- factor(seq_len(nrow(mdrrTable)))
  } else {
    mdrrTable$analysisId <- factor(mdrrTable$analysisId)
  }

  ggplot2::ggplot(mdrrTable,
    ggplot2::aes(x = mdrr, y = analysisId)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_vline(xintercept = 2, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_x_continuous(limits = c(1, NA)) +
    ggplot2::labs(
      x       = "Minimum detectable relative risk",
      y       = NULL,
      caption = "Dashed line: MDRR = 2"
    ) +
    ggplot2::theme_minimal()
}

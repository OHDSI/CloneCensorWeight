#' Compute standardized mean differences for covariate balance
#'
#' For each baseline covariate, computes the standardized mean difference
#' (SMD) between the two arms.  Call twice — once with `weights = NULL` to
#' get pre-IPCW SMDs among artificially-uncensored subjects, and once with
#' the weighted person-time data frame to get post-IPCW SMDs — then compare
#' them with [plotLovePlot()].
#'
#' The SMD formula is
#' \eqn{(\bar{x}_1 - \bar{x}_2) / \sqrt{(s_1^2 + s_2^2) / 2}},
#' using weighted means and variances when IPCW weights are supplied.
#' Missing values in binary covariates are treated as zero (the
#' FeatureExtraction convention for `missingMeansZero = "Y"` covariates).
#'
#' @param cloneCensorWeightData A `CloneCensorWeightData` object after cloning
#'   (so the `cohorts` table has an `arm` column).
#' @param weights Optional weighted person-time data frame from
#'   [createCensoringWeights()]. When `NULL`, all subjects who were **not**
#'   artificially censored receive equal weight 1.  When supplied, the final
#'   cumulative IPCW weight for each clone is used.
#'
#' @return A data frame with one row per covariate and columns `covariateId`,
#'   `covariate`, `meanArm1`, `sdArm1`, `meanArm2`, `sdArm2`, `smd`, and
#'   `absSmd`.  Column names for arm 1 and arm 2 reflect the alphabetically
#'   first and second arm labels respectively.
#' @export
computeStandardizedDifferences <- function(cloneCensorWeightData, weights = NULL) {
  checkmate::assertClass(cloneCensorWeightData, "CloneCensorWeightData")
  if (!is.null(weights)) checkmate::assertDataFrame(weights)

  covariates   <- .readCloneCensorWeightTable(cloneCensorWeightData, "covariates")
  covariateRef <- .readCloneCensorWeightTable(cloneCensorWeightData, "covariateRef")
  cohorts      <- .readCloneCensorWeightTable(cloneCensorWeightData, "cohorts")

  empty <- data.frame(
    covariateId = integer(), covariate = character(),
    meanArm1 = numeric(), sdArm1 = numeric(),
    meanArm2 = numeric(), sdArm2 = numeric(),
    smd = numeric(), absSmd = numeric(),
    stringsAsFactors = FALSE
  )

  isSparse <- all(c("rowId", "covariateId", "covariateValue") %in% names(covariates))
  if (!isSparse || nrow(covariates) == 0) return(empty)

  if (!is.null(weights)) {
    armWeights <- .finalCloneWeights(weights)
  } else {
    checkmate::assertSubset("artificiallyCensored", names(cohorts))
    uncensored <- cohorts[cohorts$artificiallyCensored == 0, c("subjectId", "arm"),
                          drop = FALSE]
    armWeights <- uncensored
    armWeights$ipcwWeight <- 1.0
  }

  armLevels <- sort(unique(as.character(armWeights$arm)))
  if (length(armLevels) != 2) {
    stop("Expected exactly two arm levels; found: ",
         paste(armLevels, collapse = ", "), call. = FALSE)
  }

  covariateIds <- sort(unique(covariates$covariateId))
  do.call(rbind, lapply(covariateIds, function(cid) {
    sub   <- covariates[covariates$covariateId == cid, c("rowId", "covariateValue"),
                        drop = FALSE]
    stats1 <- .weightedStats(sub, armWeights, armLevels[1])
    stats2 <- .weightedStats(sub, armWeights, armLevels[2])

    pooledSD <- sqrt((stats1$sd^2 + stats2$sd^2) / 2)
    smd <- if (!is.na(pooledSD) && pooledSD > 0)
      (stats1$mean - stats2$mean) / pooledSD else 0

    cname <- if (nrow(covariateRef) > 0 && cid %in% covariateRef$covariateId)
      covariateRef$covariateName[covariateRef$covariateId == cid][1]
    else paste0("cov_", cid)

    data.frame(
      covariateId = cid, covariate = cname,
      meanArm1 = stats1$mean, sdArm1 = stats1$sd,
      meanArm2 = stats2$mean, sdArm2 = stats2$sd,
      smd = smd, absSmd = abs(smd),
      stringsAsFactors = FALSE
    )
  }))
}

# Final cumulative IPCW weight per clone (last interval per cloneId).
.finalCloneWeights <- function(personTime) {
  checkmate::assertSubset(c("subjectId", "cloneId", "arm", "interval", "weight"),
                          names(personTime))
  # For each clone, take the row with the highest interval number
  maxInterval <- tapply(personTime$interval, personTime$cloneId, max)
  keep <- mapply(function(cid, maxI) {
    rows <- which(personTime$cloneId == cid & personTime$interval == maxI)
    rows[1]
  }, names(maxInterval), maxInterval)
  out <- personTime[unlist(keep), c("subjectId", "arm", "weight"), drop = FALSE]
  names(out)[names(out) == "weight"] <- "ipcwWeight"
  rownames(out) <- NULL
  out
}

# Weighted mean and population-weighted SD of a covariate for one arm.
# Subjects absent from `covSubset` receive value 0 (binary missing = zero).
.weightedStats <- function(covSubset, armWeights, armLabel) {
  arm    <- armWeights[armWeights$arm == armLabel, c("subjectId", "ipcwWeight"),
                       drop = FALSE]
  merged <- merge(arm, covSubset, by.x = "subjectId", by.y = "rowId", all.x = TRUE)
  merged$covariateValue[is.na(merged$covariateValue)] <- 0

  w  <- merged$ipcwWeight
  x  <- merged$covariateValue
  sw <- sum(w)
  if (sw == 0) return(list(mean = NA_real_, sd = NA_real_))

  wMean <- sum(w * x) / sw
  wVar  <- sum(w * (x - wMean)^2) / sw
  list(mean = wMean, sd = sqrt(wVar))
}

#' Plot a love plot of covariate balance
#'
#' Displays absolute standardized mean differences (SMD) per covariate,
#' optionally overlaying pre- and post-IPCW values.  A dashed vertical line
#' marks the conventional balance threshold (`threshold`, default 0.1).
#'
#' @param smdBefore A data frame from [computeStandardizedDifferences()] with
#'   `weights = NULL`, representing pre-IPCW balance.
#' @param smdAfter An optional data frame from [computeStandardizedDifferences()]
#'   with IPCW weights, representing post-IPCW balance.
#' @param threshold Absolute SMD threshold to mark with a dashed reference line.
#'   Default `0.1`.
#' @return A `ggplot2` object.
#' @export
plotLovePlot <- function(smdBefore, smdAfter = NULL, threshold = 0.1) {
  checkmate::assertDataFrame(smdBefore)
  checkmate::assertSubset(c("covariate", "absSmd"), names(smdBefore))
  checkmate::assertNumber(threshold, lower = 0)

  plotData          <- smdBefore[, c("covariate", "absSmd"), drop = FALSE]
  plotData$Weighting <- "Before IPCW"

  if (!is.null(smdAfter)) {
    checkmate::assertDataFrame(smdAfter)
    checkmate::assertSubset(c("covariate", "absSmd"), names(smdAfter))
    after          <- smdAfter[, c("covariate", "absSmd"), drop = FALSE]
    after$Weighting <- "After IPCW"
    plotData       <- rbind(plotData, after)
    plotData$Weighting <- factor(plotData$Weighting,
                                 levels = c("Before IPCW", "After IPCW"))
  }

  ggplot2::ggplot(plotData,
    ggplot2::aes(x = absSmd, y = covariate, colour = Weighting, shape = Weighting)) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", colour = "grey50") +
    ggplot2::scale_x_continuous(limits = c(0, NA)) +
    ggplot2::labs(
      x      = "Absolute standardized mean difference",
      y      = NULL,
      colour = NULL,
      shape  = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}

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

#' Compute covariate balance between the compliant clones of two strategies
#'
#' @description
#' Compares the baseline covariates of the clones that still follow each strategy after the grace
#' period, before and after weighting, using [CohortMethod::computeCovariateBalance()].
#'
#' @details
#' Let τ be the largest `graceEnd` of the two strategies, and d = τ + 1 (or the clone's `riskStart`, if
#' later) the first day after the grace period. The clones compared are those of either strategy still
#' under follow-up and uncensored on day d (`endDay >= d`); they are the target (strategy 1) and
#' comparator (strategy 2) groups. Each clone has the covariates of its person. The "before" columns
#' of the result compare these clones crudely; the "after" columns weight each clone by W(d), the
#' weight of its segment that contains day d (time `d - riskStart + 1`), which is the weight the outcome
#' models use for risk-set membership on day d. For `initiateWithin` and `notInitiateWithin` strategies
#' W is constant after the grace period, so this is exactly W(τ + 1). For a `neverInitiate` strategy W
#' keeps changing, and the segment's weight is that of the first event time on or after day d.
#'
#' @param weightedClones              An object of type [CcwClones] with segments and weights, as created
#'                                    by `computeCensoringWeights()`.
#' @param cloneCensorWeightData       The object of type [CloneCensorWeightData] the clones were created
#'                                    from.
#' @param computeCovariateBalanceArgs An object of type `ComputeCovariateBalanceArgs` as created using
#'                                    [CohortMethod::createComputeCovariateBalanceArgs()].
#'
#' @return
#' The tibble returned by [CohortMethod::computeCovariateBalance()]: one row per covariate, where
#' "target" is strategy 1 and "comparator" is strategy 2, "before" is the unweighted comparison of the
#' compliant clones, and "after" the weighted comparison.
#'
#' @examples
#' \dontrun{
#' balance <- computeCovariateBalance(
#'   weightedClones = weightedClones,
#'   cloneCensorWeightData = cloneCensorWeightData
#' )
#' }
#'
#' @export
computeCovariateBalance <- function(weightedClones,
                                    cloneCensorWeightData,
                                    computeCovariateBalanceArgs = CohortMethod::createComputeCovariateBalanceArgs()) {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(weightedClones, weighted = TRUE, add = errorMessages)
  .assertCloneCensorWeightData(cloneCensorWeightData, add = errorMessages)
  checkmate::assertR6(computeCovariateBalanceArgs, "ComputeCovariateBalanceArgs", add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  strategies <- attr(weightedClones, "metaData")$strategies
  graceEnds <- unlist(lapply(strategies, function(strategy) strategy$graceEnd))
  if (length(graceEnds) == 0) {
    stop("Covariate balance needs a strategy with a grace period (graceEnd)")
  }
  balanceDay <- max(graceEnds) + 1

  sql <- sprintf("
    SELECT clones.cloneRowId,
      clones.rowId,
      clones.personSeqId,
      CAST(CASE WHEN clones.strategyId = 1 THEN 1 ELSE 0 END AS INTEGER) AS treatment,
      segments.weight
    FROM clones
    INNER JOIN segments
      ON segments.cloneRowId = clones.cloneRowId
    WHERE clones.endDay >= GREATEST(%1$d, clones.riskStart)
      AND segments.startTime < GREATEST(%1$d, clones.riskStart) - clones.riskStart + 1
      AND segments.stopTime >= GREATEST(%1$d, clones.riskStart) - clones.riskStart + 1
    ORDER BY clones.cloneRowId;", as.integer(balanceDay))
  compliant <- .andromedaSql(weightedClones, sql)
  for (treatment in c(1, 0)) {
    if (!any(compliant$treatment == treatment)) {
      stop(sprintf("No clones of strategy %d are under follow-up on day %d", 2 - treatment, balanceDay))
    }
  }

  cohortMethodData <- .createCmDataView(
    cloneCensorWeightData = cloneCensorWeightData,
    cohorts = data.frame(
      rowId = compliant$cloneRowId,
      personSeqId = compliant$personSeqId,
      treatment = compliant$treatment
    ),
    rowIdMap = data.frame(newRowId = compliant$cloneRowId, rowId = compliant$rowId)
  )
  on.exit(close(cohortMethodData))
  population <- data.frame(
    rowId = compliant$cloneRowId,
    treatment = compliant$treatment,
    iptw = compliant$weight
  )
  balance <- CohortMethod::computeCovariateBalance(
    population = population,
    cohortMethodData = cohortMethodData,
    computeCovariateBalanceArgs = computeCovariateBalanceArgs
  )
  return(balance)
}

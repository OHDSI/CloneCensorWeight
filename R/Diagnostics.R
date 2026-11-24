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

#' @importFrom CohortMethod getAttritionTable
#' @export
CohortMethod::getAttritionTable

#' @importFrom CohortMethod drawAttritionDiagram
#' @export
CohortMethod::drawAttritionDiagram

#' Compute the minimum detectable relative risk of a clone-censor-weight analysis
#'
#' @description
#' Computes the minimum detectable relative risk (MDRR) from the clones, using
#' [CohortMethod::computeMdrr()].
#'
#' @details
#' Each clone is one row of the population passed to [CohortMethod::computeMdrr()]: `treatment` is 1
#' for clones of the target strategy (strategy 1) and 0 for the comparator (strategy 2),
#' `outcomeCount` is the clone's outcome count, and `timeAtRisk` its `survivalTime`. The weights are
#' not used, so the MDRR describes the unweighted clone population, and a person contributes a clone to
#' both arms. An outcome during a grace period, when the person's clones still follow both strategies,
#' is therefore counted once in each arm, so `totalOutcomes` can exceed the number of persons with the
#' outcome, and the MDRR is optimistic.
#'
#' @param clones     An object of type [CcwClones], as created by `createClones()` or
#'                   `computeCensoringWeights()`.
#' @param alpha      Type I error.
#' @param power      1 - beta, where beta is the type II error.
#' @param twoSided   Consider a two-sided test?
#' @param modelType  The type of outcome model that will be used: `"cox"` or `"logistic"`.
#'
#' @return
#' The data frame returned by [CohortMethod::computeMdrr()]: counts of persons, clones (exposures),
#' days at risk and outcomes per arm, the MDRR, and the expected standard error.
#'
#' @example inst/examples/clones.R
#' @examples
#' computeMdrr(clones)
#' close(clones)
#' close(data)
#'
#' @export
computeMdrr <- function(clones, alpha = 0.05, power = 0.8, twoSided = TRUE, modelType = "cox") {
  errorMessages <- checkmate::makeAssertCollection()
  .assertCcwClones(clones, add = errorMessages)
  checkmate::reportAssertions(collection = errorMessages)

  population <- clones$clones |>
    select("personSeqId", "strategyId", "outcomeCount", "survivalTime") |>
    collect()
  population <- data.frame(
    personSeqId = population$personSeqId,
    treatment = as.integer(population$strategyId == 1),
    outcomeCount = population$outcomeCount,
    timeAtRisk = population$survivalTime
  )
  return(CohortMethod::computeMdrr(
    population = population,
    alpha = alpha,
    power = power,
    twoSided = twoSided,
    modelType = modelType
  ))
}

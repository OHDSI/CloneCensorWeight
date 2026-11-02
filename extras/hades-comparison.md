# CloneCensorWeight and CohortMethod

## Shared HADES conventions

Both packages are organized as standard R packages with roxygen documentation,
`R/` function modules, SQL under `inst/sql`, `Andromeda`-backed analysis data,
`checkmate` validation, `ParallelLogger` messaging, `DatabaseConnector` SQL
execution, and `testthat` coverage. `CloneCensorWeightData` follows the
CohortMethod pattern by extending FeatureExtraction's S4 `CovariateData`
class, which itself extends `Andromeda`. Analysis specifications use R6
objects so settings can be validated, serialized, and reused across analyses.

## Important differences

`CohortMethod` compares two observed treatment cohorts. Its primary exposure is
fixed at cohort entry, and propensity-score design choices include trimming,
matching, stratification, and IPTW. `CloneCensorWeight` emulates a trial with
multiple treatment strategies by duplicating each eligible subject, censoring a
clone when its observed initiation violates its assigned strategy, and using
stabilized IPCW to recover the strategy-specific target population.

Consequently, a CCW row is a subject-strategy clone, not an independent person.
Outcome inference must account for repeated clones from the same subject. The
grace period defines the protocol and the artificial censoring time; it is not
an ordinary CohortMethod risk-window parameter.

## Where computation belongs

* SQL should select CDM cohorts, treatment initiation dates, outcomes, and
  baseline features. `SqlRender` should translate parameterized templates and
  `DatabaseConnector` should execute them. Large joins and date restrictions
  should remain on the database server.
* `FeatureExtraction` should build high-dimensional covariates in the database
  and expose them as `CovariateData`. `Andromeda` should hold the resulting
  sparse tables; use `collect()` only for small summaries, model metadata, or
  plotting samples.
* `Cyclops` should fit regularized pooled logistic models for numerator and
  denominator censoring hazards and high-dimensional outcome models. Ordinary
  `glm` is appropriate only as a small-data fallback or for validation fixtures.
  The current weighting implementation attempts Cyclops first and records the
   selected engine. Both the pooled logistic and Cox outcome paths use Cyclops;
   the package requires the current CohortMethod-compatible Cyclops release
   (3.6 or newer).
* Base R and `survival` are appropriate for date arithmetic, interval assembly,
  small diagnostic summaries, and plotting inputs. They should not be used to
  materialize millions of CDM covariates or person-time rows in production.

## Optimization priorities

1. Materialize the eligible clone cohort and initiation markers once as a
   database temporary table, then extract only the required rows and covariates.
2. Expand person-time in SQL or DuckDB-backed Andromeda operations rather than
   constructing one large in-memory data frame in R.
3. Fit numerator and denominator models once per analysis specification and
   reuse them across outcome definitions. Parallelize independent analyses with
   `ParallelLogger`/cluster settings, following CohortMethod's multi-analysis
   runner.
4. Trim or truncate weights only after recording untrimmed diagnostics, and
   retain effective sample size and clone-level clustering diagnostics.
5. Keep output files deterministic and reference-based: results summaries are
   small CSV/RDS files, while clone, covariate, and person-time stores remain
   Andromeda caches.
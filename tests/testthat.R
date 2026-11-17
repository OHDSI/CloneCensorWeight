library(testthat)
library(CloneCensorWeight)

# The test files run in parallel (Config/testthat/parallel in DESCRIPTION). testthat starts 2
# processes unless TESTTHAT_CPUS or options(Ncpus) says otherwise; off CRAN (NOT_CRAN=true, as set by
# devtools::check()) use up to 8, leaving 2 cores free.
if (Sys.getenv("NOT_CRAN") == "true" && Sys.getenv("TESTTHAT_CPUS") == "") {
  cores <- parallel::detectCores()
  Sys.setenv(TESTTHAT_CPUS = max(2, min(8, if (is.na(cores)) 2 else cores - 2)))
}

test_check("CloneCensorWeight")

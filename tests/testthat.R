library(testthat)
library(CloneCensorWeight)

# The test files run in parallel (Config/testthat/parallel in DESCRIPTION). testthat starts 2
# processes unless TESTTHAT_CPUS or options(Ncpus) says otherwise; off CRAN (NOT_CRAN=true, as set by
# devtools::check()) use 4.
if (Sys.getenv("NOT_CRAN") == "true" && Sys.getenv("TESTTHAT_CPUS") == "") {
  Sys.setenv(TESTTHAT_CPUS = "4")
}

test_check("CloneCensorWeight")

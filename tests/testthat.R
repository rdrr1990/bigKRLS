library(testthat)
library(bigKRLS)

# have to disable R_TESTS environment variable or parallel::makePSOCKcluster hangs.
Sys.setenv("R_TESTS" = "")


test_check("bigKRLS")


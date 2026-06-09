## Author: Mikhailova Anna
## First created: 2025-08-25
## Description: formal testing of sg-gof-tp function and its helpers
## Keywords: SimuRg, sg-gof-tp, goodness-of-fit


test_that("sg-tp-gof output is correct", {
  x <- sg_gof_tp(obj1)
  expect_true(inherits(x[[1]], "ggplot"))
  expect_snapshot(ggplot2::layer_data(x[[1]]))
})
test_that("sg-gof-tp does not work", {
  expect_error(sg_gof_tp(data.frame()))
  expect_error(sg_gof_tp(obj1, tsld = T))
})
test_that("sg-tp-gof file load", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  x <- sg_gof_tp(fpath_i, log_y = T, cov_cols = c("AGE"), col_i = "AGE",
                 sort_by = c("AGE"), desc = T)
  expect_true(inherits(x[[1]], "ggplot"))
  expect_snapshot(ggplot2::layer_data(x[[1]]))
})

## Author: Ugolkov Yaroslav
## First created: 2025-09-05
## Description: formal testing of sg-gof-obpr function and its helpers
## Keywords: SimuRg, sg-gof-obpr, goodness-of-fit

test_that("sg-gof-obpr output is correct", {
  x <- sg_gof_obpr(obj1,col_i = 'AGE',cov_cols = c('AGE'), facet_i = 'AGE')
  expect_true(inherits(x, "ggplot"))
  expect_snapshot(ggplot2::layer_data(x))
})

test_that("sg-gof-obpr does not work", {
  expect_error(sg_gof_obpr(data.frame()))
  # expect_error(sg_gof_obpr(obj1, f_scales = "invalid_scale"))
})

test_that("sg-gof-obpr file load", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  x <- sg_gof_obpr(fpath_i, log_axes = TRUE, no_leg = T)
  expect_true(inherits(x, "ggplot"))
  expect_snapshot(ggplot2::layer_data(x))
})

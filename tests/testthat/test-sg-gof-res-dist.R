## Author:Ugolkov Yaroslav
## First created: 2025-09-05
## Description: formal testing of sg-gof-res-dist function and its helpers
## Keywords: SimuRg, sg-gof-res-dist, goodness-of-fit

test_that("sg-gof-res-dist output is correct", {
  x <- sg_gof_res_dist(obj1, ndist = F)
  expect_true(inherits(x, "ggplot"))
  expect_snapshot(ggplot2::layer_data(x))
})

test_that("sg-gof-res-dist does not work", {
  expect_error(sg_gof_res_dist(data.frame()))
  # expect_error(sg_gof_res_dist(obj1, plot_type = "invalid_type"))
})

test_that("sg-gof-res-dist file load", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  x <- sg_gof_res_dist(fpath_i, res_type = 'WRES', plot_type = 'QQ')
  expect_true(inherits(x, "ggplot"))
  # expect_snapshot(ggplot2::layer_data(x))
})

## Author: Mikhailova Anna
## First created: 2025-09-05
## Description: formal testing of sg-modcomp function and its helpers
## Keywords: SimuRg, sg-modcomp, model building, comparison

test_that("sg-modcomp works", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK_cov.RData", package = "SimuRg")
  sum_tab <- sg_modcomp(fpath_i, 3)
  expect_true(inherits(sum_tab, "data.frame"))
  expect_snapshot(sum_tab)
  sum_tab <- sg_modcomp(obj1, 4)
  expect_true(inherits(sum_tab, "data.frame"))
  expect_snapshot(sum_tab)
})
test_that("sg-modcomp does not work", {
  expect_error(sg_modcomp(data.frame()))
})

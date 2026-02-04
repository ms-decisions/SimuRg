## Author: Mikhailova Anna
## First created: 2025-09-05
## Description: formal testing of sg-parsum function and its helpers
## Keywords: SimuRg, sg-parsum, goodness-of-fit

test_that("sg-parsum file load works", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  sum_tab <- sg_parsum(fpath_i)
  expect_true(inherits(sum_tab, "data.frame"))
  expect_snapshot(sum_tab)
  sum_tab <- sg_parsum(obj1)
  expect_true(inherits(sum_tab, "data.frame"))
  expect_snapshot(sum_tab)
})
test_that("sg-parsum file load works dataset with covariates", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK_cov.RData", package = "SimuRg")
  sum_tab <- sg_parsum(fpath_i)
  expect_true(inherits(sum_tab, "data.frame"))
  expect_snapshot(sum_tab)
  sum_tab <- sg_parsum(obj1)
  expect_true(inherits(sum_tab, "data.frame"))
  expect_snapshot(sum_tab)
})
test_that("sg-parsum does not work", {
  expect_error(sg_parsum(data.frame()))
})

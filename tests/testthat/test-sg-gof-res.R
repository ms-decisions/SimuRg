## Author: Ugolkov Yaroslav
## First created: 2025-09-05
## Description: formal testing of sg-gof-res function and its helpers
## Keywords: SimuRg, sg-gof-res, goodness-of-fit

test_that("sg-gof-res output is correct", {
  x <- sg_gof_res(obj1, res_type = "IWRES", cov_cols = 'AGE', col_i = 'AGE', log_x = T, facet_i = 'AGE')
  expect_true(inherits(x, "ggplot"))
  expect_snapshot(ggplot2::layer_data(x))
})

test_that("sg-gof-res does not work", {
  expect_error(sg_gof_res(data.frame(), res_type = "IWRES"))
  # expect_error(sg_gof_res(obj1, res_type = "IWRES", f_scales = "invalid_scale"))
})

test_that("sg-gof-res file load", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  x <- sg_gof_res(fpath_i, res_type = "IWRES", log_x = TRUE, addline = T, dens =T,
                  no_leg = T, cov_cols = 'AGE', col_i = 'AGE')
  expect_true(inherits(x, "ggplot"))
  expect_snapshot(ggplot2::layer_data(x))
})

# Дополнительные тесты для 100% покрытия
test_that("sg-gof-res works with different residual types", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  # Test with WRES
  x <- sg_gof_res(fpath_i, res_type = "WRES", vs_time = TRUE)
  expect_true(inherits(x, "ggplot"))

  # Test with CWRES
  x <- sg_gof_res(fpath_i, res_type = "CWRES", vs_time = FALSE)
  expect_true(inherits(x, "ggplot"))

  # Test with RES
  x <- sg_gof_res(fpath_i, res_type = "RES", vs_time = TRUE)
  expect_true(inherits(x, "ggplot"))
})

test_that("sg-gof-res works with different plot configurations", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  # Test without covariates
  x <- sg_gof_res(fpath_i, res_type = "IWRES", vs_time = TRUE, addline = FALSE, smooth = FALSE)
  expect_true(inherits(x, "ggplot"))

  # Test with PRED instead of IPRED
  x <- sg_gof_res(fpath_i, res_type = "IWRES", indiv = FALSE, vs_time = FALSE)
  expect_true(inherits(x, "ggplot"))

  # Test density plot without color
  x <- sg_gof_res(fpath_i, res_type = "IWRES", dens = TRUE, addline = FALSE)
  expect_true(inherits(x, "ggplot"))

  # Test density plot with color
  x <- sg_gof_res(fpath_i, res_type = "IWRES", dens = TRUE, cov_cols = 'AGE', col_i = 'AGE')
  expect_true(inherits(x, "ggplot"))
})

test_that("sg-gof-res works with different scale factors", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  x <- sg_gof_res(fpath_i, res_type = "IWRES", sc_factor = 2)
  expect_true(inherits(x, "ggplot"))
})

test_that("sg-gof-res works with list input", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  smrg_obj <- get(load(fpath_i))

  # Test with list input instead of file path
  x <- sg_gof_res(smrg_obj, res_type = "IWRES")
  expect_true(inherits(x, "ggplot"))
})

test_that("sg-gof-res handles axis limits", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  x <- sg_gof_res(fpath_i, res_type = "IWRES", min_y = -5, max_y = 5, min_x = 0, max_x = 100)
  expect_true(inherits(x, "ggplot"))
})

test_that("sg-gof-res works with different facet scales", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  x <- sg_gof_res(fpath_i, res_type = "IWRES", cov_cols = 'AGE', facet_i = 'AGE', f_scales = "free")
  expect_true(inherits(x, "ggplot"))

  x <- sg_gof_res(fpath_i, res_type = "IWRES", cov_cols = 'AGE', facet_i = 'AGE', f_scales = "free_x")
  expect_true(inherits(x, "ggplot"))

  x <- sg_gof_res(fpath_i, res_type = "IWRES", cov_cols = 'AGE', facet_i = 'AGE', f_scales = "free_y")
  expect_true(inherits(x, "ggplot"))
})

test_that("sg-gof-res works with different quantile bins", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  x <- sg_gof_res(fpath_i, res_type = "IWRES", cov_cols = 'AGE', n_quantiles = 4)
  expect_true(inherits(x, "ggplot"))
})

test_that("sg-gof-res works with discrete covariates", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  # Test with SEX which should be discrete
  x <- sg_gof_res(fpath_i, res_type = "IWRES", cov_cols = 'SEX', col_i = 'SEX')
  expect_true(inherits(x, "ggplot"))
})

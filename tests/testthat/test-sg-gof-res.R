## Author: Ugolkov Yaroslav
## First created: 2025-09-05
## Description: formal testing of sg-gof-res function and its helpers
## Keywords: SimuRg, sg-gof-res, goodness-of-fit

test_that("sg-gof-res output is correct", {
  skip_on_ci()
  x <- sg_gof_res(obj1, cov_cols = 'AGE', col_i = 'AGE', log_x = T, facet_i = 'AGE', vs_time = F, indiv = F)
  expect_true(inherits(x, "ggplot"))
  expect_snapshot(ggplot2::layer_data(x))
})

test_that("sg-gof-res does not work", {
  expect_error(sg_gof_res(data.frame()))
  # expect_error(sg_gof_res(obj1, f_scales = "invalid_scale"))
})

test_that("sg-gof-res file load", {
  skip_on_ci()
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  x <- sg_gof_res(fpath_i, log_x = TRUE, addline = T, legend_fl = T,
                  cov_cols = 'AGE', col_i = 'AGE', n_quantiles = 10)
  expect_true(inherits(x, "ggplot"))
  expect_snapshot(ggplot2::layer_data(x))
})

# Дополнительные тесты для 100% покрытия
test_that("sg-gof-res works with different residual types", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  # Test with WRES
  x <- sg_gof_res(fpath_i, indiv = F, vs_time = TRUE)
  expect_true(inherits(x, "ggplot"))

  # Test with RES
  x <- sg_gof_res(fpath_i, vs_time = TRUE, weighted = F, f_scales = "free")
  expect_true(inherits(x, "ggplot"))
})

test_that("sg-gof-res works with different plot configurations", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  # Test without covariates
  x <- sg_gof_res(fpath_i, vs_time = TRUE, addline = FALSE, smooth = FALSE)
  expect_true(inherits(x, "ggplot"))

  # Test with PRED instead of IPRED
  x <- sg_gof_res(fpath_i, indiv = FALSE, vs_time = FALSE)
  expect_true(inherits(x, "ggplot"))

})

test_that("sg-gof-res works with list input", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  smrg_obj <- get(load(fpath_i))

  # Test with list input instead of file path
  x <- sg_gof_res(smrg_obj)
  expect_true(inherits(x, "ggplot"))
})

test_that("sg-gof-res handles axis limits and labels", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  x <- sg_gof_res(fpath_i, min_y = -0.5, max_y = 0.5, min_x = 0, max_x = 100,
                  weighted = F, lab_y = "Individual residuals, uM", lab_x = "Time, h")
  expect_true(inherits(x, "ggplot"))
})

test_that("sg-gof-res works with different facet scales", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  x <- sg_gof_res(fpath_i, vs_time = F, cov_cols = 'AGE', facet_i = 'AGE', f_scales = "free")
  expect_true(inherits(x, "ggplot"))

  x <- sg_gof_res(fpath_i, vs_time = F, cov_cols = 'AGE', facet_i = 'AGE', f_scales = "free_x")
  expect_true(inherits(x, "ggplot"))

  x <- sg_gof_res(fpath_i, cov_cols = 'AGE', facet_i = 'AGE', f_scales = "free_y")
  expect_true(inherits(x, "ggplot"))
})

test_that("sg-gof-res works with different quantile bins", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  x <- sg_gof_res(fpath_i, cov_cols = 'AGE', col_i = 'AGE', n_quantiles = 4, legend_fl = T)
  expect_true(inherits(x, "ggplot"))
})

test_that("sg-gof-res works with discrete covariates", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  # Test with SEX which should be discrete
  x <- sg_gof_res(fpath_i, cov_cols = 'SEX', col_i = 'SEX', legend_fl = T)
  expect_true(inherits(x, "ggplot"))
})

test_that("sg-gof-res scales x, exponentiates, labels quantiles, and prefers CWRES", {
  set.seed(1)
  mock_obj <- list(
    SDTAB = data.frame(
      ID = rep(1:4, each = 5),
      TIME = rep(c(1, 2, 4, 8, 12), 4),
      DV = rep(c(10, 20, 40, 80, 100), 4),
      PRED = rep(c(12, 22, 38, 70, 90), 4),
      IPRED = rep(log(c(10, 20, 40, 80, 100)), 4),
      IWRES = rnorm(20),
      IRES = rnorm(20),
      WRES = rnorm(20),
      CWRES = rnorm(20),
      RES = rnorm(20),
      MDV = 0
    ),
    COTAB = data.frame(ID = 1:4, AGE = c(20, 40, 60, 80), WT = c(50, 70, 90, 110)),
    CATAB = data.frame(ID = 1:4, SEX = c("F", "M", "F", "M"))
  )

  point_x <- function(p) {
    idx <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1L)))
    ggplot2::layer_data(p, idx[[1]])$x
  }

  p_time <- sg_gof_res(mock_obj, vs_time = TRUE, sc_factor = 2, smooth = FALSE, addline = FALSE)
  expect_equal(range(point_x(p_time)), c(0.5, 6))

  p_pred <- sg_gof_res(
    mock_obj, vs_time = FALSE, indiv = TRUE, exp_dv = TRUE, sc_factor = 10,
    smooth = FALSE, addline = FALSE
  )
  expect_equal(range(point_x(p_pred)), c(1, 10), tolerance = 1e-8)

  p_cwres <- sg_gof_res(mock_obj, indiv = FALSE, weighted = TRUE, smooth = FALSE, addline = FALSE)
  expect_identical(p_cwres$labels$y, "CWRES")

  p_leg <- sg_gof_res(
    mock_obj, cov_cols = "AGE", col_i = "AGE", legend_fl = TRUE,
    levels_discrete = 2, smooth = FALSE, addline = FALSE
  )
  col_scale <- Filter(function(s) "colour" %in% s$aesthetics, p_leg$scales$scales)[[1]]
  expect_true(any(grepl("^AGE: Q1, ", col_scale$labels)))

  p_num <- sg_gof_res(mock_obj, log_x = TRUE, log_labs = "number", smooth = FALSE, addline = FALSE)
  p_math <- sg_gof_res(mock_obj, log_x = TRUE, log_labs = "math", smooth = FALSE, addline = FALSE)
  num_labs <- as.character(p_num$scales$get_scales("x")$labels(c(1, 10, 100)))
  math_labs <- as.character(p_math$scales$get_scales("x")$labels(c(1, 10, 100)))
  expect_false(any(grepl("10\\^", num_labs)))
  expect_true(any(grepl("10", math_labs)))
})

test_that("sg-gof-res works with multiple covariates", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

  x <- sg_gof_res(fpath_i, cov_cols = c('AGE','WEIGHT'), col_i = c('AGE'), legend_fl = T, facet_i = 'WEIGHT')
  expect_true(inherits(x, "ggplot"))
})

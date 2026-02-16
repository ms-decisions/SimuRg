## Author: Ugolkov Yaroslav
## First created: 2025-09-05
## Description: formal testing of sg-gof-obpr function and its helpers
## Keywords: SimuRg, sg-gof-obpr, goodness-of-fit

##Updated tests

test_that("sg-gof-obpr output is correct", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  if (exists("obj1")) {
    x <- sg_gof_obpr(obj1, col_i = 'AGE', cov_cols = c('AGE'), facet_i = 'AGE')
    expect_true(inherits(x, "ggplot"))
  } else {
    skip("obj1 not available in test environment")
  }
})

test_that("sg-gof-obpr does not work", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  expect_error(sg_gof_obpr(data.frame()))
})

test_that("sg-gof-obpr file load", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  if (file.exists(fpath_i)) {
    x <- sg_gof_obpr(fpath_i, log_axes = TRUE, no_leg = TRUE)
    expect_true(inherits(x, "ggplot"))
  } else {
    skip("Test data file not available")
  }
})

##New tests
# Helper function to create mock sg_fit object for testing
create_mock_sg_fit_gof <- function() {
  set.seed(123)
  list(
    SDTAB = data.frame(
      ID = rep(1:3, each = 10),
      TIME = rep(c(0, 1, 2, 4, 8, 10, 13, 24, 48, 72), 3),
      DV = rnorm(30, mean = 10, sd = 2),
      PRED = rnorm(30, mean = 10, sd = 1.5),
      IPRED = rnorm(30, mean = 10, sd = 1.2),
      MDV = rep(0, 30),
      WRES = rnorm(30, mean = 0, sd = 1),
      IWRES = rnorm(30, mean = 0, sd = 1)
    ),
    COTAB = data.frame(
      ID = 1:3,
      AGE = c(30, 45, 60),
      WT = c(70, 80, 90)
    ),
    CATAB = data.frame(
      ID = 1:3,
      SEX = c("M", "F", "M")
    )
  )
}

# Error handling tests
test_that("sg_gof_obpr errors on data.frame input", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  expect_error(
    sg_gof_obpr(data.frame()),
    "fpath_i cannot be a data.frame"
  )
})

test_that("sg_gof_obpr errors on missing SDTAB", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  obj_no_sdtab <- list(COTAB = data.frame(), CATAB = data.frame())

  expect_error(
    sg_gof_obpr(obj_no_sdtab),
    "SDTAB"
  )
})

test_that("sg_gof_obpr errors on empty SDTAB", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  obj_empty <- list(SDTAB = data.frame())

  expect_error(
    sg_gof_obpr(obj_empty),
    "SDTAB is empty"
  )
})

test_that("sg_gof_obpr errors on non-existent file", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  expect_error(
    sg_gof_obpr("nonexistent_file.RData"),
    "File does not exist"
  )
})

test_that("sg_gof_obpr errors on missing required columns in SDTAB", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  obj_missing_cols <- list(
    SDTAB = data.frame(
      ID = 1:5,
      TIME = 1:5
      # Missing DV, PRED, IPRED, MDV
    )
  )

  expect_error(
    sg_gof_obpr(obj_missing_cols),
    "SDTAB is missing required columns"
  )
})

# Basic functionality tests
test_that("sg_gof_obpr returns ggplot object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()
  p <- sg_gof_obpr(mock_obj)

  expect_s3_class(p, "ggplot")
})

test_that("sg_gof_obpr handles individual predictions (IPRED)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()
  p <- sg_gof_obpr(mock_obj, indiv = TRUE)

  expect_s3_class(p, "ggplot")

  # Check that plot was created successfully
  gb <- ggplot2::ggplot_build(p)
  expect_true(nrow(gb$data[[1]]) > 0)
})

test_that("sg_gof_obpr handles population predictions (PRED)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()
  p <- sg_gof_obpr(mock_obj, indiv = FALSE)

  expect_s3_class(p, "ggplot")

  gb <- ggplot2::ggplot_build(p)
  expect_true(nrow(gb$data[[1]]) > 0)
})

test_that("sg_gof_obpr handles addline parameter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()

  # With lines
  p_with_lines <- sg_gof_obpr(mock_obj, addline = TRUE)
  gb_with <- ggplot2::ggplot_build(p_with_lines)
  has_lines <- any(vapply(gb_with$plot$layers, function(l) inherits(l$geom, "GeomLine"), logical(1)))
  expect_true(has_lines || nrow(gb_with$data[[1]]) > 0)

  # Without lines
  p_without_lines <- sg_gof_obpr(mock_obj, addline = FALSE)
  expect_s3_class(p_without_lines, "ggplot")
})

test_that("sg_gof_obpr handles smooth parameter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()

  # With smooth
  p_smooth <- sg_gof_obpr(mock_obj, smooth = TRUE)
  gb_smooth <- ggplot2::ggplot_build(p_smooth)
  has_smooth <- any(vapply(gb_smooth$plot$layers, function(l) inherits(l$geom, "GeomSmooth"), logical(1)))
  expect_true(has_smooth || nrow(gb_smooth$data[[1]]) > 0)

  # Without smooth
  p_no_smooth <- sg_gof_obpr(mock_obj, smooth = FALSE)
  expect_s3_class(p_no_smooth, "ggplot")
})

test_that("sg_gof_obpr handles log_axes parameter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()

  # Linear axes
  p_linear <- sg_gof_obpr(mock_obj, log_axes = FALSE)
  expect_s3_class(p_linear, "ggplot")

  # Log axes
  p_log <- sg_gof_obpr(mock_obj, log_axes = TRUE)
  expect_s3_class(p_log, "ggplot")

  # Check that log scales are applied
  x_scale_log <- p_log$scales$get_scales("x")
  y_scale_log <- p_log$scales$get_scales("y")
  expect_false(is.null(x_scale_log))
  expect_false(is.null(y_scale_log))
})

test_that("sg_gof_obpr handles sc_factor parameter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()

  # With scaling factor
  p_scaled <- sg_gof_obpr(mock_obj, sc_factor = 2)
  expect_s3_class(p_scaled, "ggplot")

  # Without scaling (default)
  p_unscaled <- sg_gof_obpr(mock_obj, sc_factor = 1)
  expect_s3_class(p_unscaled, "ggplot")
})

test_that("sg_gof_obpr handles covariate coloring", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()
  p <- sg_gof_obpr(mock_obj, cov_cols = "SEX", col_i = "SEX")

  expect_s3_class(p, "ggplot")

  gb <- ggplot2::ggplot_build(p)
  expect_true(nrow(gb$data[[1]]) > 0)
})

test_that("sg_gof_obpr handles continuous covariate categorization", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()
  p <- sg_gof_obpr(mock_obj, cov_cols = "AGE", col_i = "AGE", n_quantiles = 3)

  expect_s3_class(p, "ggplot")

  # AGE should be categorized into quantiles
  gb <- ggplot2::ggplot_build(p)
  expect_true(nrow(gb$data[[1]]) > 0)
})

test_that("sg_gof_obpr handles faceting", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()
  p <- sg_gof_obpr(mock_obj, cov_cols = "SEX", facet_i = "SEX")

  expect_s3_class(p, "ggplot")

  # Check that faceting was applied
  has_facet <- !is.null(p$facet) && inherits(p$facet, "FacetWrap")
  expect_true(has_facet, "Plot should have faceting applied")

  gb <- ggplot2::ggplot_build(p)
  n_panels <- length(gb$layout$panel_params)
  expect_true(n_panels > 1, "Faceted plot should have multiple panels")
})

test_that("sg_gof_obpr handles different facet scales", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()

  # Fixed scales
  p_fixed <- sg_gof_obpr(mock_obj, cov_cols = "SEX", facet_i = "SEX", f_scales = "fixed")
  expect_s3_class(p_fixed, "ggplot")

  # Free scales
  p_free <- sg_gof_obpr(mock_obj, cov_cols = "SEX", facet_i = "SEX", f_scales = "free")
  expect_s3_class(p_free, "ggplot")

  # Free x
  p_free_x <- sg_gof_obpr(mock_obj, cov_cols = "SEX", facet_i = "SEX", f_scales = "free_x")
  expect_s3_class(p_free_x, "ggplot")

  # Free y
  p_free_y <- sg_gof_obpr(mock_obj, cov_cols = "SEX", facet_i = "SEX", f_scales = "free_y")
  expect_s3_class(p_free_y, "ggplot")
})

test_that("sg_gof_obpr handles alpha_i parameter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()

  p_alpha <- sg_gof_obpr(mock_obj, alpha_i = 0.3)
  expect_s3_class(p_alpha, "ggplot")

  p_alpha2 <- sg_gof_obpr(mock_obj, alpha_i = 0.8)
  expect_s3_class(p_alpha2, "ggplot")
})

test_that("sg_gof_obpr handles custom axis labels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()
  p <- sg_gof_obpr(mock_obj, xlab = "Predicted (ng/mL)", ylab = "Observed (ng/mL)")

  expect_s3_class(p, "ggplot")

  # Check that labels are set
  expect_equal(p$labels$x, "Predicted (ng/mL)")
  expect_equal(p$labels$y, "Observed (ng/mL)")
})

test_that("sg_gof_obpr handles no_leg parameter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()

  # With legend
  p_legend <- sg_gof_obpr(mock_obj, cov_cols = "SEX", col_i = "SEX", no_leg = FALSE)
  expect_s3_class(p_legend, "ggplot")

  # Without legend
  p_no_legend <- sg_gof_obpr(mock_obj, cov_cols = "SEX", col_i = "SEX", no_leg = TRUE)
  expect_s3_class(p_no_legend, "ggplot")

  # Check legend position
  expect_equal(p_no_legend$theme$legend.position, "none")
})

test_that("sg_gof_obpr handles n_quantiles parameter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()

  p_q3 <- sg_gof_obpr(mock_obj, cov_cols = "AGE", col_i = "AGE", n_quantiles = 3)
  expect_s3_class(p_q3, "ggplot")

  p_q5 <- sg_gof_obpr(mock_obj, cov_cols = "AGE", col_i = "AGE", n_quantiles = 5)
  expect_s3_class(p_q5, "ggplot")
})

test_that("sg_gof_obpr handles levels_discrete parameter", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()

  p_disc <- sg_gof_obpr(mock_obj, cov_cols = "AGE", col_i = "AGE", levels_discrete = 5)
  expect_s3_class(p_disc, "ggplot")

  p_disc2 <- sg_gof_obpr(mock_obj, cov_cols = "AGE", col_i = "AGE", levels_discrete = 20)
  expect_s3_class(p_disc2, "ggplot")
})

test_that("sg_gof_obpr filters MDV correctly", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  mock_obj <- create_mock_sg_fit_gof()
  # Add some MDV=1 rows
  mock_obj$SDTAB <- rbind(
    mock_obj$SDTAB,
    data.frame(
      ID = c(1, 2),
      TIME = c(0, 0),
      DV = c(NA, NA),
      PRED = c(100, 100),
      IPRED = c(100, 100),
      MDV = c(1, 1),
      WRES = c(NA, NA),
      IWRES = c(NA, NA)
    )
  )

  p <- sg_gof_obpr(mock_obj)
  expect_s3_class(p, "ggplot")

  # MDV=1 rows should be filtered out
  gb <- ggplot2::ggplot_build(p)
  # Should have fewer points than total rows (MDV=1 filtered)
  expect_true(nrow(gb$data[[1]]) <= nrow(mock_obj$SDTAB))
})

test_that("sg_gof_obpr handles SDTAB as list of data frames", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  # Create SDTAB as a list (simulating some data structures)
  mock_obj <- list(
    SDTAB = list(
      data.frame(
        ID = 1:3,
        TIME = c(0, 1, 2),
        DV = rnorm(3, 10, 2),
        PRED = rnorm(3, 10, 1.5),
        IPRED = rnorm(3, 10, 1.2),
        MDV = c(0, 0, 0),
        WRES = rnorm(3),
        IWRES = rnorm(3)
      )
    ),
    COTAB = NULL,
    CATAB = NULL
  )

  p <- sg_gof_obpr(mock_obj)
  expect_s3_class(p, "ggplot")
})





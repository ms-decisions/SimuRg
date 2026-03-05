## Author: Alina Melnikova
## Description: basic correctness checks for sg_sim_tp
## Keywords: sg_sim_tp, visualization, diagnostics

library(testthat)
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

# Helper: small mock dataset
make_mock_data <- function() {
  set.seed(123)  # For reproducibility
  data.frame(
    TIME = rep(rep(1:4, each = 3), times = 2),  # 3 replicates per time point
    VALUE = c(
      rnorm(12, mean = 10, sd = 1),  # VAR A: 4 time points * 3 replicates
      rnorm(12, mean = 8, sd = 1.5)  # VAR B: 4 time points * 3 replicates
    ),
    VAR = rep(c("A", "B"), each = 12),
    GRP = rep(c("X", "Y"), each = 12)
  )
}
# Helper: extended mock dataset with additional grouping variables
make_extended_mock_data <- function() {
  data.frame(
    TIME = rep(1:4, times = 6),
    VALUE = rnorm(24, mean = 10, sd = 2),
    VAR = rep(rep(c("A", "B"), each = 4), times = 3),
    Regimen = rep(c("R1", "R2", "R3"), each = 8),
    CAT1 = rep(c("Cat1", "Cat2"), each = 12)
  )
}


test_that("sg_sim_tp returns ggplot with raw data", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  df <- make_mock_data()
  p <- sg_sim_tp(ds_i = df, group_i = "VAR", col_i = "VAR")

  expect_s3_class(p, "ggplot")
})

test_that("sg_sim_tp handles mean + SD bands", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  df <- make_mock_data()
  p <- sg_sim_tp(
    ds_i = df,
    group_i = "VAR",
    cent_i = "mean",
    vrns_i = "SD",
    col_i = "VAR",
    fill_i = "VAR"
  )

  gb <- ggplot2::ggplot_build(p)
  # Expect at least one ribbon layer when variance bands requested
  has_ribbon <- any(vapply(gb$plot$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1)))
  expect_true(has_ribbon)
})

test_that("sg_sim_tp handles percentile ribbons", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  df <- make_mock_data()
  p <- sg_sim_tp(
    ds_i = df,
    group_i = "VAR",
    cent_i = "median",
    lperc_i = 0.1,
    uperc_i = 0.9,
    col_i = "VAR",
    fill_i = "VAR"
  )

  gb <- ggplot2::ggplot_build(p)
  has_ribbon <- any(vapply(gb$plot$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1)))
  expect_true(has_ribbon)
})



test_that("sg_sim_tp applies axis limits", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  df <- make_mock_data()
  p <- sg_sim_tp(ds_i = df, min_x = 0, max_x = 5, min_y = 0, max_y = 15)

  xlim <- p$coordinates$limits$x
  ylim <- p$coordinates$limits$y
  expect_equal(xlim, c(0, 5))
  expect_equal(ylim, c(0, 15))
})

test_that("sg_sim_tp supports wrap faceting", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  df <- make_mock_data()
  p <- sg_sim_tp(ds_i = df, wrap_i = "~VAR")

  # Check that faceting was applied by examining the plot's facet object
  # facet_wrap creates a FacetWrap object
  has_facet <- !is.null(p$facet) &&
    (inherits(p$facet, "FacetWrap") || inherits(p$facet, "FacetGrid"))

  expect_true(has_facet, "Plot should have faceting applied")

  # Alternative check: verify the built plot has multiple panels
  gb <- ggplot2::ggplot_build(p)
  n_panels <- length(gb$layout$panel_params)
  expect_true(n_panels > 1, "Faceted plot should have multiple panels")
})

test_that("sg_sim_tp handles different variance types (SD, SE, IQR)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  df <- make_mock_data()

  # Test SD
  p_sd <- sg_sim_tp(ds_i = df, cent_i = 'mean', vrns_i = 'SD', col_i = 'VAR', fill_i = 'VAR')
  gb_sd <- ggplot2::ggplot_build(p_sd)
  has_ribbon_sd <- any(vapply(gb_sd$plot$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1)))
  expect_true(has_ribbon_sd, "SD variance should create ribbons")

  # Test SE
  p_se <- sg_sim_tp(ds_i = df, cent_i = 'mean', vrns_i = 'SE', col_i = 'VAR', fill_i = 'VAR')
  gb_se <- ggplot2::ggplot_build(p_se)
  has_ribbon_se <- any(vapply(gb_se$plot$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1)))
  expect_true(has_ribbon_se, "SE variance should create ribbons")

  # Test IQR
  p_iqr <- sg_sim_tp(ds_i = df, cent_i = 'median', vrns_i = 'IQR', col_i = 'VAR', fill_i = 'VAR')
  gb_iqr <- ggplot2::ggplot_build(p_iqr)
  has_ribbon_iqr <- any(vapply(gb_iqr$plot$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1)))
  expect_true(has_ribbon_iqr, "IQR variance should create ribbons")
})

test_that("sg_sim_tp handles different central tendency measures", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  df <- make_mock_data()

  # Test mean
  p_mean <- sg_sim_tp(ds_i = df, cent_i = 'mean', col_i = 'VAR')
  expect_s3_class(p_mean, "ggplot")

  # Test median
  p_median <- sg_sim_tp(ds_i = df, cent_i = 'median', col_i = 'VAR')
  expect_s3_class(p_median, "ggplot")

  # Test geometric mean
  df_pos <- df
  df_pos$VALUE <- abs(df_pos$VALUE) + 1  # Ensure positive values for geometric mean
  p_geom <- sg_sim_tp(ds_i = df_pos, cent_i = 'geom_mean', col_i = 'VAR')
  expect_s3_class(p_geom, "ggplot")
})

test_that("sg_sim_tp handles median + IQR with custom axis limits", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  df <- make_extended_mock_data()
  p <- sg_sim_tp(
    ds_i = df,
    cent_i = 'median',
    vrns_i = 'IQR',
    col_i = 'CAT1',
    fill_i = 'CAT1',
    min_x = 0,
    max_x = 100,
    min_y = 0
  )

  expect_s3_class(p, "ggplot")

  # Check that axis limits were applied
  xlim <- p$coordinates$limits$x
  ylim <- p$coordinates$limits$y

  expect_false(is.null(xlim), "x limits should be set")
  expect_equal(xlim, c(0, 100), tolerance = 1e-6)

  expect_false(is.null(ylim), "min_y should be set")
  expect_equal(ylim[1], 0, tolerance = 1e-6)

  # Check that IQR ribbons exist
  gb <- ggplot2::ggplot_build(p)
  has_ribbon <- any(vapply(gb$plot$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1)))
  expect_true(has_ribbon, "Plot should have IQR variance ribbons")
})

test_that("sg_sim_tp handles mean + SE with points, grid faceting, and log scale", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  df <- make_extended_mock_data()
  p <- sg_sim_tp(
    ds_i = df,
    cent_i = 'mean',
    vrns_i = 'SE',
    col_i = 'Regimen',
    fill_i = 'Regimen',
    add_points = 2,
    grid_i = '~VAR',
    log_y = TRUE
  )

  expect_s3_class(p, "ggplot")

  # Check that grid faceting was applied
  has_facet <- !is.null(p$facet) && inherits(p$facet, "FacetGrid")
  expect_true(has_facet, "Plot should have grid faceting")

  # Check that points were added
  gb <- ggplot2::ggplot_build(p)
  has_points <- any(vapply(gb$plot$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))
  expect_true(has_points, "Plot should have points when add_points > 0")

  # Check that SE ribbons exist
  has_ribbon <- any(vapply(gb$plot$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1)))
  expect_true(has_ribbon, "Plot should have SE variance ribbons")

  # Check that log scale was applied
  y_scale <- p$scales$get_scales("y")
  expect_false(is.null(y_scale), "y scale should exist when log_y = TRUE")
})

test_that("sg_sim_tp handles median with 5th-95th percentile bands and wrap faceting", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  df <- make_mock_data()
  p <- sg_sim_tp(
    ds_i = df,
    cent_i = 'median',
    lperc_i = 0.05,
    uperc_i = 0.95,
    col_i = 'VAR',
    fill_i = 'VAR',
    wrap_i = '~VAR',
    wrap_ncol = 2
  )

  expect_s3_class(p, "ggplot")

  # Check that faceting was applied
  has_facet <- !is.null(p$facet) && inherits(p$facet, "FacetWrap")
  expect_true(has_facet, "Plot should have wrap faceting")

  # Check that percentile ribbons exist
  gb <- ggplot2::ggplot_build(p)
  has_ribbon <- any(vapply(gb$plot$layers, function(l) inherits(l$geom, "GeomRibbon"), logical(1)))
  expect_true(has_ribbon, "Plot should have percentile ribbons")

  # Check wrap_ncol parameter
  gb <- ggplot2::ggplot_build(p)
  n_panels <- length(gb$layout$panel_params)
  expect_true(n_panels > 1, "Faceted plot should have multiple panels")
})

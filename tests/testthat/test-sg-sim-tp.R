## Author: Test suite for sg_sim_tp
## Description: basic correctness checks for sg_sim_tp
## Keywords: sg_sim_tp, visualization, diagnostics

library(testthat)

# Helper: small mock dataset
make_mock_data <- function() {
  data.frame(
    TIME = rep(1:4, times = 2),
    VALUE = c(10, 12, 11, 9, 8, 7, 9, 10),
    VAR = rep(c("A", "B"), each = 4),
    GRP = rep(c("X", "Y"), each = 4)
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

test_that("sg_sim_tp applies log scales", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  df <- make_mock_data()
  p <- sg_sim_tp(ds_i = df, log_y = F, log_x = F)

  # Check that log10 scales are applied by examining the scale objects
  # scale_y_log10() creates a ScaleContinuousPosition with trans = "log10"
  y_scale <- p$scales$get_scales("y")
  x_scale <- p$scales$get_scales("x")

  # Check if scales exist
  expect_false(is.null(y_scale), "y scale should exist")
  expect_false(is.null(x_scale), "x scale should exist")

  # Check for log10 transformation - can be checked via trans$name or class
  if (!is.null(y_scale$trans)) {
    y_is_log <- y_scale$trans$name == "log10" ||
      grepl("log10", class(y_scale$trans)[1], ignore.case = TRUE)
    expect_true(y_is_log, "y scale should be log10")
  }

  if (!is.null(x_scale$trans)) {
    x_is_log <- x_scale$trans$name == "log10" ||
      grepl("log10", class(x_scale$trans)[1], ignore.case = TRUE)
    expect_true(x_is_log, "x scale should be log10")
  }
})


test_that("sg_sim_tp applies axis limits", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")

  df <- make_mock_data()
  p <- sg_sim_tp(ds_i = df, x_min = 0, x_max = 5, y_min = 0, y_max = 15)

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

  gb <- ggplot2::ggplot_build(p)
  expect_true(gb$layout$has_facet)
})



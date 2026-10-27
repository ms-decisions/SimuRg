## Author:Ugolkov Yaroslav
## First created: 2025-09-05
## Description: formal testing of sg-gof-res-dist function and its helpers
## Keywords: SimuRg, sg-gof-res-dist, goodness-of-fit

test_that("sg-gof-res-dist output is correct", {
  skip_on_ci()
  x <- sg_gof_res_dist(obj1, ndist = F)
  expect_true(inherits(x, "ggplot"))
  expect_snapshot(ggplot2::layer_data(x))
})

test_that("sg-gof-res-dist does not work", {
  expect_error(sg_gof_res_dist(data.frame()))
  # expect_error(sg_gof_res_dist(obj1, plot_type = "invalid_type"))
})

test_that("sg-gof-res-dist file load", {
  skip_on_ci()
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  x <- sg_gof_res_dist(fpath_i, res_type = 'WRES', plot_type = 'QQ')
  expect_true(inherits(x, "ggplot"))
  # expect_snapshot(ggplot2::layer_data(x))
})

test_that("sg-gof-res-dist facets one residual by a covariate", {
  set.seed(1)
  mock <- list(
    SDTAB = data.frame(
      ID = rep(1:4, each = 40),
      DVID = 1,
      IWRES = c(rnorm(40, 0, 0.5), rnorm(40, 2, 0.5), rnorm(80, 0, 0.5)),
      CWRES = rnorm(160)
    ),
    CATAB = data.frame(
      ID = 1:4,
      POP = c("Healthy", "Asthma", "CRNP", "Atopic dermatitis")
    )
  )

  p <- sg_gof_res_dist(
    mock,
    res_type = "IWRES",
    cov_cols = "POP",
    facet_i = "POP",
    ndist = TRUE
  )
  expect_true(inherits(p, "ggplot"))
  expect_equal(p$scales$get_scales("x")$name, "IWRES")
  expect_equal(p$scales$get_scales("y")$name, "Density")
  expect_equal(nrow(ggplot2::ggplot_build(p)$layout$layout), 4)
  expect_false(p$facet$params$free$x)
  expect_false(p$facet$params$free$y)

  line_data <- ggplot2::layer_data(p, i = 2)
  panel_median <- tapply(line_data$x, line_data$PANEL, median)
  expect_true(abs(panel_median[[1]] - panel_median[[2]]) > 0.5)

  p_cw <- sg_gof_res_dist(
    mock,
    res_type = "CWRES",
    cov_cols = "POP",
    facet_i = "POP",
    ndist = FALSE
  )
  expect_equal(p_cw$scales$get_scales("x")$name, "CWRES")
  expect_equal(length(p_cw$layers), 1)
})

test_that("sg-gof-res-dist bins a continuous covariate and rejects an unknown one", {
  mock <- list(
    SDTAB = data.frame(
      ID = rep(1:6, each = 10),
      DVID = 1,
      IWRES = rnorm(60)
    ),
    COTAB = data.frame(
      ID = 1:6,
      AGE = c(20, 30, 40, 50, 60, 70)
    )
  )

  p <- sg_gof_res_dist(
    mock,
    res_type = "IWRES",
    cov_cols = "AGE",
    facet_i = "AGE",
    n_quantiles = 3,
    levels_discrete = 4
  )
  age_panels <- ggplot2::ggplot_build(p)$layout$layout
  expect_equal(nrow(age_panels), 3)
  expect_equal(
    as.character(age_panels$AGE),
    c("Q1 [20, 36.67]", "Q2 (36.67, 53.33]", "Q3 (53.33, 70]")
  )

  wide <- list(
    SDTAB = data.frame(
      ID = 1:12,
      DVID = 1,
      IWRES = rnorm(12)
    ),
    CATAB = data.frame(
      ID = 1:12,
      POP = paste0("G", 1:12)
    )
  )
  p_wide <- sg_gof_res_dist(wide, res_type = "IWRES", facet_i = "POP", ndist = FALSE)
  expect_equal(nrow(ggplot2::ggplot_build(p_wide)$layout$layout), 12)

  expect_error(
    sg_gof_res_dist(mock, res_type = "IWRES", facet_i = "POP"),
    "not found"
  )
  expect_error(
    sg_gof_res_dist(mock, res_type = "CWRES"),
    "CWRES"
  )
})

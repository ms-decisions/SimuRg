## Author: Ugolkov Yaroslav
## First created: 2025-09-05
## Description: formal testing of sg-gof-par-dist function and its helpers
## Keywords: SimuRg, sg-gof-par-dist, goodness-of-fit

test_that("Test_1_1: sg-gof-par-dist output is correct", {
  skip_on_ci()
  x <- sg_gof_par_dist(obj1, tdist =F)
  expect_true(inherits(x, "ggplot"))
  expect_snapshot(ggplot2::layer_data(x))
})

test_that("Test_2_2: sg-gof-par-dist does not work", {
  expect_error(sg_gof_par_dist(data.frame()))
  # expect_error(sg_gof_par_dist(obj1, plot_type = "invalid_type"))
})

test_that("Test_3_3: sg-gof-par-dist file load", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  x <- sg_gof_par_dist(fpath_i, tdist = TRUE, plot_type = 'QQ')
  expect_true(inherits(x, "ggplot"))
  # expect_snapshot(ggplot2::layer_data(x))
})

# ---- New coverage tests for sg_gof_par_dist (mock GFO layout based on run1c_smrg_results.json) ----
make_mock_gfo_for_par_dist <- function(include_shrinkage = TRUE, include_distribution = TRUE) {
  sumtab <- data.frame(
    PAR = c("ka_pop", "V_pop", "Cl_pop", "omega_ka", "omega_V", "omega_Cl", "b"),
    VALUE = c(1.0, 10.0, 1.2, 0.08, 0.15, 0.30, 0.20),
    TYPE = c(
      "Typical values", "Typical values", "Typical values",
      "Random effects", "Random effects", "Random effects", "Residual error model"
    ),
    EST = "ESTIMATED",
    SE = c(0.02, 0.20, 0.06, 0.02, 0.02, 0.03, 0.01),
    RSE = c(2, 2, 5, 25, 13, 10, 5),
    LCI = c(0.96, 9.6, 1.1, 0.04, 0.11, 0.24, 0.18),
    UCI = c(1.04, 10.4, 1.3, 0.12, 0.19, 0.36, 0.22),
    stringsAsFactors = FALSE
  )

  if (include_distribution) {
    sumtab$DISTRIBUTION <- c("lognormal", "normal", "logitnormal", NA, NA, NA, NA)
  }
  if (include_shrinkage) {
    sumtab$ETAshrinkage_var <- c(NA, NA, NA, 80, 20, 5, NA)
    sumtab$ETAshrinkage_sd <- c(NA, NA, NA, 60, 10, 3, NA)
  }

  patab <- data.frame(
    ID = 1:8,
    eta_ka = c(-0.10, -0.05, 0.00, 0.03, 0.08, -0.02, 0.05, -0.01),
    eta_V = c(0.20, 0.10, -0.05, 0.00, 0.15, -0.10, 0.06, -0.03),
    eta_Cl = c(0.30, 0.20, 0.10, -0.20, 0.05, -0.15, 0.12, -0.08),
    ka = c(0.90, 0.95, 1.00, 1.03, 1.08, 0.98, 1.05, 0.99),
    V = c(12.0, 11.2, 9.8, 10.1, 12.8, 9.6, 10.7, 10.3),
    Cl = c(1.7, 1.5, 1.3, 1.0, 1.4, 0.9, 1.6, 1.1),
    stringsAsFactors = FALSE
  )

  gfo <- list(
    PROJNAME = "mock_run1c",
    SUMTAB = sumtab,
    PATAB = patab,
    COVMAT = diag(c(0.01, 0.02, 0.03)),
    CORRMAT = diag(3),
    OFV = data.frame(LL = 100, AIC = 110, BIC = 115, BICc = 116),
    OMEGAMAT = matrix(c(0.08, 0, 0, 0, 0.15, 0, 0, 0, 0.30), nrow = 3, byrow = TRUE),
    SIGMAMAT = matrix(0.2, nrow = 1)
  )

  rownames(gfo$COVMAT) <- c("ka_pop", "V_pop", "Cl_pop")
  colnames(gfo$COVMAT) <- c("ka_pop", "V_pop", "Cl_pop")
  rownames(gfo$CORRMAT) <- c("ka_pop", "V_pop", "Cl_pop")
  colnames(gfo$CORRMAT) <- c("ka_pop", "V_pop", "Cl_pop")
  rownames(gfo$OMEGAMAT) <- c("omega_ka", "omega_V", "omega_Cl")
  colnames(gfo$OMEGAMAT) <- c("omega_ka", "omega_V", "omega_Cl")

  gfo
}

test_that("Test_4_4: sg-gof-par-dist mock GFO DIST without theoretical density returns ggplot", {
  mock_gfo <- make_mock_gfo_for_par_dist()
  x <- sg_gof_par_dist(mock_gfo, tdist = FALSE, plot_type = "DIST")
  expect_true(inherits(x, "ggplot"))
})

test_that("Test_5_5: sg-gof-par-dist mock GFO DIST with theoretical density adds density layer", {
  mock_gfo <- make_mock_gfo_for_par_dist()
  x <- sg_gof_par_dist(mock_gfo, tdist = TRUE, plot_type = "DIST")
  expect_true(inherits(x, "ggplot"))
  expect_true(any(vapply(x$layers, function(z) inherits(z$geom, "GeomDensity"), logical(1))))
})

test_that("Test_6_6: sg-gof-par-dist QQ plot with selected base parameter maps to eta_*", {
  mock_gfo <- make_mock_gfo_for_par_dist()
  x <- sg_gof_par_dist(mock_gfo, plot_type = "QQ", par_seq = "ka")
  expect_true(inherits(x, "ggplot"))
  expect_equal(unique(x$data$PAR), "eta_ka")
})

test_that("Test_7_7: sg-gof-par-dist QQ plot warns for missing eta parameter names", {
  mock_gfo <- make_mock_gfo_for_par_dist()
  expect_warning(
    sg_gof_par_dist(mock_gfo, plot_type = "QQ", par_seq = c("ka", "not_existing")),
    "Parameters not found in PATAB"
  )
})

test_that("Test_8_8: sg-gof-par-dist correlations plot returns recordedplot object", {
  mock_gfo <- make_mock_gfo_for_par_dist()
  x <- sg_gof_par_dist(mock_gfo, plot_type = "correlations", par_seq = c("ka", "V", "Cl"))
  expect_true(inherits(x, "recordedplot"))
})

test_that("Test_9_9: sg-gof-par-dist RE mode produces ggplot for eta parameters", {
  mock_gfo <- make_mock_gfo_for_par_dist()
  x <- sg_gof_par_dist(mock_gfo, par_type = "RE", par_seq = c("eta_ka", "eta_V"), tdist = TRUE, plot_type = "DIST")
  expect_true(inherits(x, "ggplot"))
})

test_that("Test_10_10: sg-gof-par-dist warns on missing DIST parameters but still returns plot for valid ones", {
  mock_gfo <- make_mock_gfo_for_par_dist()
  expect_warning(
    sg_gof_par_dist(mock_gfo, par_seq = c("ka", "bad_parameter"), tdist = FALSE, plot_type = "DIST"),
    "Parameters not found in PATAB"
  )
  x <- suppressWarnings(
    sg_gof_par_dist(mock_gfo, par_seq = c("ka", "bad_parameter"), tdist = FALSE, plot_type = "DIST")
  )
  expect_true(inherits(x, "ggplot"))
})

test_that("Test_11_11: sg-gof-par-dist returns NULL when all requested parameters are absent", {
  mock_gfo <- make_mock_gfo_for_par_dist()
  x <- suppressWarnings(sg_gof_par_dist(mock_gfo, par_seq = c("bad1", "bad2"), tdist = FALSE, plot_type = "DIST"))
  expect_null(x)
})

test_that("Test_12_12: sg-gof-par-dist handles missing ETA shrinkage column with message", {
  mock_gfo <- make_mock_gfo_for_par_dist(include_shrinkage = FALSE)
  expect_message(
    x <- sg_gof_par_dist(mock_gfo, tdist = FALSE, plot_type = "DIST"),
    "Shrinkage was not calculated in GFO object"
  )
  expect_true(inherits(x, "ggplot"))
})

test_that("Test_13_13: sg-gof-par-dist reads JSON containing top-level GFO/GCO layout", {
  mock_gfo <- make_mock_gfo_for_par_dist()
  f_json <- tempfile(fileext = ".json")
  on.exit(unlink(f_json), add = TRUE)
  jsonlite::write_json(list(GFO = mock_gfo, GCO = list(PROJNAME = "mock")), f_json, auto_unbox = TRUE, pretty = TRUE)

  x <- sg_gof_par_dist(f_json, tdist = FALSE, plot_type = "DIST")
  expect_true(inherits(x, "ggplot"))
})

test_that("Test_14_14: sg-gof-par-dist maps OMEGAMAT V1 V2 V3 names from SUMTAB omega order", {
  mock_gfo <- make_mock_gfo_for_par_dist()
  rownames(mock_gfo$OMEGAMAT) <- c("V1", "V2", "V3")
  colnames(mock_gfo$OMEGAMAT) <- c("V1", "V2", "V3")

  x <- sg_gof_par_dist(mock_gfo, tdist = TRUE, plot_type = "DIST")
  expect_true(inherits(x, "ggplot"))
  expect_true(any(vapply(x$layers, function(z) inherits(z$geom, "GeomDensity"), logical(1))))
  expect_setequal(unique(x$data$PAR), c("ka", "V", "Cl"))
})

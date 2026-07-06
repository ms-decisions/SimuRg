## Description: formal testing of sg-predist-vis function
## Keywords: SimuRg, sg-predist-vis, prediction distribution, diagnostics

# Helper: mock simulation results in the long format produced by sg_predist_sim()
# (sg_sim() returns columns sim.id, ID, TIME, VAR, VALUE when omega is used).
create_mock_predist_sim <- function(n_ids = 2, n_times = 5, n_sims = 20, var_name = "Cc") {
  times <- c(0, 1, 2, 4, 8)[seq_len(n_times)]
  data.frame(
    sim.id = rep(seq_len(n_sims), each = n_ids * n_times),
    ID = rep(rep(seq_len(n_ids), each = n_times), n_sims),
    TIME = rep(rep(times, n_ids), n_sims),
    VAR = var_name,
    VALUE = rnorm(n_ids * n_times * n_sims, mean = 10, sd = 2)
  )
}

# Helper: minimal mock sg_fit object providing SDTAB for the observed-data overlay.
# SDTAB carries a column named after the output variable so the dt_obs_fl path works.
create_mock_obj <- function(var_name = "Cc", n_ids = 2, n_times = 5) {
  times <- c(0, 1, 2, 4, 8)[seq_len(n_times)]
  sdtab <- data.frame(
    ID = rep(seq_len(n_ids), each = n_times),
    TIME = rep(times, n_ids),
    MDV = 0
  )
  sdtab[[var_name]] <- rnorm(n_ids * n_times, mean = 10, sd = 2)
  list(SDTAB = sdtab)
}

# ---- Basic functionality ----
test_that("sg_predist_vis returns a named list of ggplot objects", {

  ds_sim <- create_mock_predist_sim()
  obj <- create_mock_obj()

  result <- sg_predist_vis(fpath_i = obj, ds_sim = ds_sim)

  expect_type(result, "list")
  expect_length(result, 1)
  expect_named(result, "Cc")
  expect_true(all(vapply(result, function(x) inherits(x, "ggplot"), logical(1))))
})

test_that("sg_predist_vis produces a plot that builds without error", {

  ds_sim <- create_mock_predist_sim()
  obj <- create_mock_obj()

  p <- sg_predist_vis(fpath_i = obj, ds_sim = ds_sim)[[1]]
  built <- ggplot2::ggplot_build(p)

  expect_s3_class(built, "ggplot_built")
  # Ribbon layer must carry finite ymin/ymax (guards against the summarise()
  # column-naming regression where L_Q/H_Q were not found).
  ribbon <- built$data[[1]]
  expect_true(all(is.finite(ribbon$ymin)))
  expect_true(all(is.finite(ribbon$ymax)))
})

test_that("sg_predist_vis handles multiple output variables", {

  ds_sim <- rbind(
    create_mock_predist_sim(var_name = "Cc"),
    create_mock_predist_sim(var_name = "Cp")
  )
  obj <- create_mock_obj()

  result <- sg_predist_vis(fpath_i = obj, ds_sim = ds_sim)

  expect_length(result, 2)
  expect_named(result, c("Cc", "Cp"))
})

# ---- Prediction interval handling ----
test_that("sg_predist_vis accepts every supported prediction interval", {

  ds_sim <- create_mock_predist_sim()
  obj <- create_mock_obj()

  for (pi in c("95%", "90%", "80%", "50%")) {
    result <- sg_predist_vis(fpath_i = obj, ds_sim = ds_sim, pred_interval = pi)
    expect_true(inherits(result[[1]], "ggplot"))
  }
})

test_that("sg_predist_vis widens the band for a wider prediction interval", {

  ds_sim <- create_mock_predist_sim(n_sims = 200)
  obj <- create_mock_obj()

  b50 <- ggplot2::ggplot_build(
    sg_predist_vis(fpath_i = obj, ds_sim = ds_sim, pred_interval = "50%")[[1]]
  )
  b95 <- ggplot2::ggplot_build(
    sg_predist_vis(fpath_i = obj, ds_sim = ds_sim, pred_interval = "95%")[[1]]
  )

  width_50 <- mean(b50$data[[1]]$ymax - b50$data[[1]]$ymin)
  width_95 <- mean(b95$data[[1]]$ymax - b95$data[[1]]$ymin)

  expect_gt(width_95, width_50)
})

test_that("sg_predist_vis errors on an invalid prediction interval", {

  ds_sim <- create_mock_predist_sim()
  obj <- create_mock_obj()

  expect_error(
    sg_predist_vis(fpath_i = obj, ds_sim = ds_sim, pred_interval = "99%"),
    regexp = "pred_interval"
  )
})

# ---- Input validation ----
test_that("sg_predist_vis errors when ds_sim lacks required columns", {

  ds_sim <- data.frame(ID = 1:5, TIME = 1:5)  # missing VAR and VALUE
  obj <- create_mock_obj()

  expect_error(
    sg_predist_vis(fpath_i = obj, ds_sim = ds_sim),
    regexp = "VAR|VALUE"
  )
})

test_that("sg_predist_vis errors on invalid fpath_i", {

  ds_sim <- create_mock_predist_sim()

  expect_error(sg_predist_vis(fpath_i = 123, ds_sim = ds_sim))
  expect_error(sg_predist_vis(fpath_i = "nonexistent_file.rda", ds_sim = ds_sim))
})

# ---- Display options ----
test_that("sg_predist_vis honours the log-y option", {

  ds_sim <- create_mock_predist_sim()
  obj <- create_mock_obj()

  p_lin <- sg_predist_vis(fpath_i = obj, ds_sim = ds_sim, logy = FALSE)[[1]]
  p_log <- sg_predist_vis(fpath_i = obj, ds_sim = ds_sim, logy = TRUE)[[1]]

  expect_true(inherits(p_lin, "ggplot"))
  expect_true(inherits(p_log, "ggplot"))
  # The log scale registers a log-10 transform on the y aesthetic.
  expect_match(p_log$scales$get_scales("y")$trans$name, "log")
})

test_that("sg_predist_vis honours the legend flag", {

  ds_sim <- create_mock_predist_sim()
  obj <- create_mock_obj()

  p_legend <- sg_predist_vis(fpath_i = obj, ds_sim = ds_sim, legend_fl = TRUE)[[1]]
  p_none   <- sg_predist_vis(fpath_i = obj, ds_sim = ds_sim, legend_fl = FALSE)[[1]]

  expect_equal(p_none$theme$legend.position, "none")
  expect_true(inherits(p_legend, "ggplot"))
})

test_that("sg_predist_vis overlays observed data when dt_obs_fl = TRUE", {

  ds_sim <- create_mock_predist_sim(var_name = "Cc")
  obj <- create_mock_obj(var_name = "Cc")

  p_obs    <- sg_predist_vis(fpath_i = obj, ds_sim = ds_sim, dt_obs_fl = TRUE)[[1]]
  p_no_obs <- sg_predist_vis(fpath_i = obj, ds_sim = ds_sim, dt_obs_fl = FALSE)[[1]]

  # Observed overlay adds an extra (point) layer.
  expect_gt(length(p_obs$layers), length(p_no_obs$layers))
})

test_that("sg_predist_vis honours custom axis labels", {

  ds_sim <- create_mock_predist_sim()
  obj <- create_mock_obj()

  p <- sg_predist_vis(fpath_i = obj, ds_sim = ds_sim,
                      name_x = "Time (h)", name_y = "Concentration")[[1]]

  built <- ggplot2::ggplot_build(p)
  expect_equal(built$plot$scales$get_scales("x")$name, "Time (h)")
  expect_equal(built$plot$scales$get_scales("y")$name, "Concentration")
})

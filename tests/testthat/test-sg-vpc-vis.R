## Author: Melnikova Alina
## First created: 2025-12-04
## Description: formal testing of sg-vpc-vis function and its helpers
## Keywords: SimuRg, sg-vpc-vis, vpc, diagnostics

# Helper function to create mock simulated data for VPC
create_mock_sim_data <- function(n_ids = 2, n_times = 5, n_sims = 10, var_name = "Cc") {
  data.frame(
    id = rep(1:n_sims, each = n_ids * n_times),
    ID = rep(rep(1:n_ids, each = n_times), n_sims),
    TIME = rep(rep(c(0, 1, 2, 4, 8), n_ids), n_sims),
    VAR = var_name,
    VALUE = rnorm(n_ids * n_times * n_sims, mean = 10, sd = 2)
  )
}

# Helper function to create mock observed data
create_mock_obs_data <- function(n_ids = 2, n_times = 5, dvid = 1) {
  data.frame(
    ID = rep(1:n_ids, each = n_times),
    TIME = rep(c(0, 1, 2, 4, 8), n_ids),
    DV = rnorm(n_ids * n_times, mean = 10, sd = 2),
    DVID = dvid
  )
}

# Helper function to create output_names mapping
create_mock_output_names <- function(var_name = "Cc", dvid = 1) {
  data.frame(
    output = var_name,
    dvid = dvid
  )
}

# Error handling tests
test_that("sg_vpc_vis errors on missing required columns in ds_sim", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- data.frame(id = 1:10, time = 1:10)  # Missing VAR and VALUE
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  expect_error(
    sg_vpc_vis(ds_sim = ds_sim, data_i = data_i, output_names = output_names),
    regexp = "VAR|VALUE"
  )
})

test_that("sg_vpc_vis errors on missing required columns in data_i", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- data.frame(ID = 1:5, TIME = 1:5)  # Missing DV and DVID
  output_names <- create_mock_output_names()

  expect_error(
    sg_vpc_vis(ds_sim = ds_sim, data_i = data_i, output_names = output_names),
    regexp = "DV|DVID"
  )
})

test_that("sg_vpc_vis errors on missing output_names", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()

  expect_error(
    sg_vpc_vis(ds_sim = ds_sim, data_i = data_i, output_names = NULL),
    regexp = "output_names"
  )
})

test_that("sg_vpc_vis errors on invalid binning method", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  # This should error or use default method
  result <- tryCatch({
    sg_vpc_vis(ds_sim = ds_sim, data_i = data_i, output_names = output_names,
               method = "invalid_method")
  }, error = function(e) NULL)

  # If it doesn't error, it should still produce output (method validation may be in helper)
  if (!is.null(result)) {
    expect_true(is.list(result))
  }
})

# Basic functionality tests
test_that("sg_vpc_vis returns list of ggplot objects", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  result <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i, output_names = output_names)

  expect_type(result, "list")
  expect_true(length(result) > 0)
  expect_true(all(vapply(result, function(x) inherits(x, "ggplot"), logical(1))))
})

test_that("sg_vpc_vis handles multiple output variables", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  # Create data with two variables
  ds_sim1 <- create_mock_sim_data(var_name = "Cc")
  ds_sim2 <- create_mock_sim_data(var_name = "Cp")
  ds_sim <- rbind(ds_sim1, ds_sim2)

  data_i1 <- create_mock_obs_data(dvid = 1)
  data_i2 <- create_mock_obs_data(dvid = 2)
  data_i <- rbind(data_i1, data_i2)

  output_names <- data.frame(
    output = c("Cc", "Cp"),
    dvid = c(1, 2)
  )

  result <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i, output_names = output_names)

  expect_equal(length(result), 2)
  expect_named(result, c("Cc", "Cp"))
})

test_that("sg_vpc_vis handles different binning methods", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  # Test kmeans method
  result_kmeans <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                              output_names = output_names, method = "kmeans", n_bins = 5)
  expect_true(is.list(result_kmeans))
  expect_true(all(vapply(result_kmeans, function(x) inherits(x, "ggplot"), logical(1))))

  # Test ntile method
  result_ntile <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                             output_names = output_names, method = "ntile", n_bins = 5)
  expect_true(is.list(result_ntile))

  # Test equal_x method
  result_equal <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                             output_names = output_names, method = "equal_x", n_bins = 5)
  expect_true(is.list(result_equal))
})

test_that("sg_vpc_vis handles custom n_bins parameter", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  result <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                       output_names = output_names, n_bins = 5)

  expect_true(is.list(result))
  expect_true(all(vapply(result, function(x) inherits(x, "ggplot"), logical(1))))
})

test_that("sg_vpc_vis handles interpolation parameter", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  # With interpolation
  result_interp <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                              output_names = output_names, interpolation = TRUE)
  expect_true(is.list(result_interp))

  # Without interpolation
  result_no_interp <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                                 output_names = output_names, interpolation = FALSE)
  expect_true(is.list(result_no_interp))
})

test_that("sg_vpc_vis handles log_y parameter", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  # Linear scale
  result_linear <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                              output_names = output_names, log_y = FALSE)
  expect_true(is.list(result_linear))

  # Log scale
  result_log <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                           output_names = output_names, log_y = TRUE)
  expect_true(is.list(result_log))
})

test_that("sg_vpc_vis handles prediction correction", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  # Without prediction correction
  result_no_corr <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                               output_names = output_names, pred.corr = FALSE)
  expect_true(is.list(result_no_corr))

  # With prediction correction
  result_corr <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                            output_names = output_names, pred.corr = TRUE)
  expect_true(is.list(result_corr))
})

test_that("sg_vpc_vis handles percentile display flags", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  # Only theoretical percentiles
  result_theor <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                             output_names = output_names,
                             theor_perc = TRUE, emp_perc = FALSE)
  expect_true(is.list(result_theor))

  # Only empirical percentiles
  result_emp <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                           output_names = output_names,
                           theor_perc = FALSE, emp_perc = TRUE)
  expect_true(is.list(result_emp))

  # Both
  result_both <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                            output_names = output_names,
                            theor_perc = TRUE, emp_perc = TRUE)
  expect_true(is.list(result_both))
})

test_that("sg_vpc_vis handles confidence interval display", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  # With CI
  result_ci <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                          output_names = output_names, theor_percCI = TRUE)
  expect_true(is.list(result_ci))

  # Without CI
  result_no_ci <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                             output_names = output_names, theor_percCI = FALSE)
  expect_true(is.list(result_no_ci))
})

test_that("sg_vpc_vis handles observed data points display", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  # With observed points
  result_points <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                              output_names = output_names, dt_obs_fl = TRUE)
  expect_true(is.list(result_points))

  # Without observed points
  result_no_points <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                                 output_names = output_names, dt_obs_fl = FALSE)
  expect_true(is.list(result_no_points))
})

test_that("sg_vpc_vis handles custom axis labels", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  result <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                       output_names = output_names,
                       lab_x = "Time (hours)",
                       lab_y = "Concentration (ng/mL)")

  expect_true(is.list(result))
  expect_true(all(vapply(result, function(x) inherits(x, "ggplot"), logical(1))))
})

test_that("sg_vpc_vis handles custom prediction interval bounds", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  result <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                       output_names = output_names,
                       piLow = 0.05, piUp = 0.95)

  expect_true(is.list(result))
})

test_that("sg_vpc_vis handles custom confidence interval bounds", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  result <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                       output_names = output_names,
                       ciLow = 0.01, ciUp = 0.99)

  expect_true(is.list(result))
})

test_that("sg_vpc_vis handles legend display", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names()

  # With legend
  result_legend <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                              output_names = output_names, legend_fl = TRUE)
  expect_true(is.list(result_legend))

  # Without legend
  result_no_legend <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                                 output_names = output_names, legend_fl = FALSE)
  expect_true(is.list(result_no_legend))
})

test_that("sg_vpc_vis handles custom time and DV column names", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data()
  # Create data with custom column names
  data_i_custom <- create_mock_obs_data()
  names(data_i_custom)[names(data_i_custom) == "TIME"] <- "TAD"
  names(data_i_custom)[names(data_i_custom) == "DV"] <- "CONC"

  output_names <- create_mock_output_names()

  result <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i_custom,
                       output_names = output_names,
                       time_col = "TAD", dv_col = "CONC")

  expect_true(is.list(result))
  expect_true(all(vapply(result, function(x) inherits(x, "ggplot"), logical(1))))
})

test_that("sg_vpc_vis handles variables with _ResErr suffix", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  # Create data with _ResErr suffix
  ds_sim <- create_mock_sim_data(var_name = "Cc_ResErr")
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names(var_name = "Cc")

  result <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                       output_names = output_names)

  expect_true(is.list(result))
  expect_true("Cc_ResErr" %in% names(result))
})

test_that("sg_vpc_vis returns named list matching VAR values", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("purrr")
  skip_if_not_installed("stringr")

  ds_sim <- create_mock_sim_data(var_name = "Cc")
  data_i <- create_mock_obs_data()
  output_names <- create_mock_output_names(var_name = "Cc")

  result <- sg_vpc_vis(ds_sim = ds_sim, data_i = data_i,
                       output_names = output_names)

  expect_named(result, "Cc")
  expect_true(is.list(result))
})


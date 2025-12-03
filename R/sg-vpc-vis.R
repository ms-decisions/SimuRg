## Author: Mikhailova Anna
## First created: 2025-11-17
## Description: functions for vpc visualisation
## Keywords: SimuRg, vpc, diagnostics

#' Visual Predictive Check (VPC) Function
#'
#' Generates VPC plots to assess model performance by comparing observed data
#' with prediction intervals derived from simulated data.
#'
#' @param ds_sim Data frame with simulated data containing columns:
#'   - id: patient identifier
#'   - sim.id: simulation replicate identifier
#'   - time: time points
#'   - VAR: output variable name
#'   - VALUE: simulated values
#' @param data_i Data frame with observed data containing columns:
#'   - ID: patient identifier
#'   - TIME: time points
#'   - DV: observed dependent variable values
#'   - DVID: dependent variable identifier
#' @param output_names Data frame mapping VAR names to DVID
#' @param x Character. Name of time column in data_i (default: "TIME")
#' @param y Character. Name of DV column in data_i (default: "DV")
#' @param logy Logical. Use log scale for y-axis (default: TRUE)
#' @param piLow Numeric. Lower prediction interval bound (default: 0.10)
#' @param piUp Numeric. Upper prediction interval bound (default: 0.90)
#' @param ciLow Numeric. Lower confidence interval bound (default: 0.025)
#' @param ciUp Numeric. Upper confidence interval bound (default: 0.975)
#' @param pred.corr Logical. Apply prediction correction (default: FALSE)
#' @param name_y Character. Y-axis label
#' @param name_x Character. X-axis label
#' @param theor_perc Logical. Show theoretical percentiles (default: TRUE)
#' @param theor_percCI Logical. Show CI around theoretical percentiles (default: TRUE)
#' @param emp_perc Logical. Show empirical percentiles (default: TRUE)
#' @param dt_obs_fl Logical. Show observed data points (default: FALSE)
#' @param legend_fl Logical. Show legend (default: FALSE)
#' @param bins Integer. Number of bins for time stratification (default: 10)
#' @param method Character. Binning method: "kmeans", "equal", "quantile" (default: "kmeans")
#' @param interpolation Logical. Use line interpolation vs rectangles (default: FALSE)
#' @param strat_by_dose Character. Variable name for dose stratification (default: NULL)
#'
#' @return List of ggplot objects, one for each output variable
#'
#' @examples
#' \dontrun{
#' vpc_plots <- generate_vpc(
#'   ds_sim = sim_data,
#'   data_i = obs_data,
#'   output_names = output_map,
#'   name_x = "Time (hours)",
#'   name_y = "Concentration (ng/mL)",
#'   bins = 8,
#'   pred.corr = TRUE
#' )
#' }
#' @import dplyr tidyr ggplot2 purrr stringr
#' @export

sg_vpc_vis <- function(ds_sim,
                         data_i,
                         output_names,
                         x = "TIME",
                         y = "DV",
                         logy = TRUE,
                         piLow = 0.10,
                         piUp = 0.90,
                         ciLow = 0.025,
                         ciUp = 0.975,
                         pred.corr = FALSE,
                         name_y,
                         name_x,
                         theor_perc = TRUE,
                         theor_percCI = TRUE,
                         emp_perc = TRUE,
                         dt_obs_fl = FALSE,
                         legend_fl = FALSE,
                         bins = 10,
                         method = "kmeans",
                         interpolation = FALSE,
                         strat_by_dose = NULL) {


  # Standardize column names in ds_sim
  ds_sim <- ds_sim %>%
    rename(ID = id, TIME = time)

  # Get unique output variables
  unique_vars <- unique(ds_sim$VAR)

  # Generate VPC for each output variable
  vpc_plots <- purrr::map(unique_vars, function(var) {

    # Filter data for current variable
    ds_sim_var <- ds_sim %>%
      filter(VAR == var)

    # Match observed data to simulation variable
    var_clean <- stringr::str_remove(var, "_ResErr")
    target_dvid <- output_names$dvid[output_names$output == var_clean]

    data_i_var <- data_i %>%
      filter(DVID == target_dvid) %>%
      rename(TIME = all_of(x), DV = all_of(y))

    # Store original observed data for overlay
    data_obs <- data_i_var %>%
      mutate(DT_label = "Data")

    # ===== TIME BINNING =====
    data_i_var <- bin_data(data_i_var, bins = bins, method = method)

    # Extract binning assignments
    bin_assignments <- data_i_var %>%
      select(ID, TIME, BIN, TIME_BIN, TIME_BIN_min, TIME_BIN_max)

    # Apply bins to simulated data
    ds_sim_var <- ds_sim_var %>%
      left_join(bin_assignments, by = c("ID", "TIME"))

    # ===== PREDICTION CORRECTION =====
    if (pred.corr) {

      # Calculate population prediction (PPRED) for each individual observation
      ppred <- ds_sim_var %>%
        group_by(ID, TIME) %>%
        summarise(PPRED = mean(VALUE, na.rm = TRUE), .groups = "drop")

      # Calculate mean PPRED per bin
      ppred_bin <- ds_sim_var %>%
        group_by(BIN) %>%
        summarise(PPRED_bin = mean(VALUE, na.rm = TRUE), .groups = "drop")

      # Join PPREDs to observed data
      data_i_var <- data_i_var %>%
        left_join(ppred, by = c("ID", "TIME")) %>%
        left_join(ppred_bin, by = "BIN")

      # Join PPREDs to simulated data
      ds_sim_var <- ds_sim_var %>%
        left_join(ppred, by = c("ID", "TIME")) %>%
        left_join(ppred_bin, by = "BIN")

      # Apply prediction correction: DV_corr = DV * (mean_PPRED_bin / PPRED_ij)
      data_i_var <- data_i_var %>%
        mutate(DV = DV * PPRED_bin / PPRED)

      ds_sim_var <- ds_sim_var %>%
        mutate(VALUE = VALUE * PPRED_bin / PPRED)
    }

    # ===== CALCULATE PERCENTILES =====

    # Define grouping variables
    vec_group <- if (!is.null(strat_by_dose)) {
      c(strat_by_dose, "TIME_BIN", "TIME_BIN_min", "TIME_BIN_max")
    } else {
      c("TIME_BIN", "TIME_BIN_min", "TIME_BIN_max")
    }

    # Summary functions
    summ_pi <- list(
      median = ~median(., na.rm = TRUE),
      pi_l = ~quantile(., piLow, na.rm = TRUE),
      pi_u = ~quantile(., piUp, na.rm = TRUE)
    )

    summ_ci <- list(
      median = ~median(., na.rm = TRUE),
      ci_l = ~quantile(., ciLow, na.rm = TRUE),
      ci_u = ~quantile(., ciUp, na.rm = TRUE)
    )

    # Empirical percentiles from observed data
    emp_int <- data_i_var %>%
      group_by(across(all_of(vec_group))) %>%
      summarise(across(DV, summ_pi), .groups = "drop")

    # Prediction intervals from simulated data (per simulation)
    pred_int <- ds_sim_var %>%
      mutate(sim.id = as.numeric(sim.id)) %>%
      group_by(across(all_of(c(vec_group, "sim.id")))) %>%
      summarise(across(VALUE, summ_pi), .groups = "drop")

    # Confidence intervals around prediction intervals
    conf_int <- pred_int %>%
      group_by(across(all_of(vec_group))) %>%
      summarise(across(c(median, pi_l, pi_u), summ_ci), .groups = "drop")

    # ===== PREPARE DATA FOR PLOTTING =====

    # Theoretical percentiles (median of prediction intervals)
    pred_int_median <- conf_int %>%
      select(all_of(vec_group), median_median, pi_l_median, pi_u_median) %>%
      pivot_longer(
        cols = c("median_median", "pi_l_median", "pi_u_median"),
        values_to = "theor_median",
        names_to = "PI"
      ) %>%
      mutate(PI = case_when(
        str_detect(PI, "pi_l") ~ "lower",
        str_detect(PI, "median_") ~ "median",
        str_detect(PI, "pi_u") ~ "upper"
      ))

    # Lower confidence interval
    pred_int_l_ci <- conf_int %>%
      select(all_of(vec_group), median_ci_l, pi_l_ci_l, pi_u_ci_l) %>%
      pivot_longer(
        cols = c("median_ci_l", "pi_l_ci_l", "pi_u_ci_l"),
        values_to = "theor_LCI",
        names_to = "PI"
      ) %>%
      mutate(PI = case_when(
        str_detect(PI, "pi_l") ~ "lower",
        str_detect(PI, "median_") ~ "median",
        str_detect(PI, "pi_u") ~ "upper"
      ))

    # Upper confidence interval
    pred_int_u_ci <- conf_int %>%
      select(all_of(vec_group), median_ci_u, pi_l_ci_u, pi_u_ci_u) %>%
      pivot_longer(
        cols = c("median_ci_u", "pi_l_ci_u", "pi_u_ci_u"),
        values_to = "theor_UCI",
        names_to = "PI"
      ) %>%
      mutate(PI = case_when(
        str_detect(PI, "pi_l") ~ "lower",
        str_detect(PI, "median_") ~ "median",
        str_detect(PI, "pi_u") ~ "upper"
      ))

    # Empirical percentiles
    emp_w <- emp_int %>%
      pivot_longer(
        cols = c("median", "pi_l", "pi_u"),
        values_to = "OBS",
        names_to = "PI"
      ) %>%
      mutate(PI = case_when(
        str_detect(PI, "pi_l") ~ "lower",
        str_detect(PI, "median") ~ "median",
        str_detect(PI, "pi_u") ~ "upper"
      ))

    # Combine all elements
    pred_int_ci <- pred_int_median %>%
      left_join(pred_int_l_ci, by = c(vec_group, "PI")) %>%
      left_join(pred_int_u_ci, by = c(vec_group, "PI")) %>%
      left_join(emp_w, by = c(vec_group, "PI"))

    # Prepare datasets for legend distinction
    emp_w_mod <- emp_w %>%
      mutate(TYPE = "empirical percentiles")

    pred_int_median_mod <- pred_int_median %>%
      rename(OBS = theor_median) %>%
      mutate(TYPE = "theoretical percentiles")

    pred_int_median_emp_mod <- bind_rows(pred_int_median_mod, emp_w_mod) %>%
      mutate(TYPE = as.factor(TYPE))

    # Add CI labels
    pred_int_l_ci_mod <- pred_int_l_ci %>%
      mutate(CI_label = ifelse(
        PI == "median",
        "CI for the median",
        paste0("CI for the ", piLow, " and ", piUp, " percentiles")
      ))

    pred_int_u_ci_mod <- pred_int_u_ci %>%
      mutate(CI_label = ifelse(
        PI == "median",
        "CI for the median",
        paste0("CI for the ", piLow, " and ", piUp, " percentiles")
      ))

    pred_int_ci_mod <- left_join(pred_int_l_ci_mod, pred_int_u_ci_mod,
                                 by = c(vec_group, "PI", "CI_label"))

    # ===== CREATE PLOT =====

    vpc_plot <- create_vpc_plot(
      pred_int_median_emp_mod = pred_int_median_emp_mod,
      pred_int_ci_mod = pred_int_ci_mod,
      data_obs = data_obs,
      pred_int_ci = pred_int_ci,
      data_i_var = data_i_var,
      theor_perc = theor_perc,
      emp_perc = emp_perc,
      theor_percCI = theor_percCI,
      interpolation = interpolation,
      dt_obs_fl = dt_obs_fl,
      strat_by_dose = strat_by_dose,
      logy = logy,
      pred.corr = pred.corr,
      name_y = name_y,
      name_x = name_x,
      legend_fl = legend_fl
    )

    return(vpc_plot)
  })

  names(vpc_plots) <- unique_vars
  return(vpc_plots)
}

## Author: Mikhailova Anna
## First created: 2025-12-03
## Description: functions for vpc visualisation
## Keywords: SimuRg, vpc, diagnostics

add_percentile_lines <- function(plot_i, line_data, ci_data, show_ci, interpolation) {

  if (interpolation) {
    # Use smooth lines
    plot_i <- plot_i +
      geom_line(data = line_data,
                aes(x = TIME_BIN, y = OBS,
                    group = interaction(PI, TYPE), lty = TYPE),
                col = "royalblue", linewidth = 0.8)

    if (show_ci) {
      plot_i <- plot_i +
        geom_ribbon(data = ci_data %>% filter(PI != "median"),
                    aes(x = TIME_BIN, ymin = theor_LCI,
                        ymax = theor_UCI, group = PI,
                        fill = CI_label),
                    col = NA, linewidth = 1, alpha = 0.3) +
        geom_ribbon(data = ci_data %>% filter(PI == "median"),
                    aes(x = TIME_BIN, ymin = theor_LCI,
                        ymax = theor_UCI, group = PI,
                        fill = CI_label),
                    col = NA, linewidth = 1, alpha = 0.3)
    }

  } else {
    # Use step rectangles
    plot_i <- plot_i +
      geom_line(data = line_data,
                aes(x = TIME_BIN, y = OBS,
                    group = interaction(PI, TYPE), lty = TYPE),
                col = "royalblue", linewidth = 0.8)

    if (show_ci) {
      plot_i <- plot_i +
        geom_rect(data = ci_data %>% filter(PI != "median"),
                  aes(xmin = TIME_BIN_min, xmax = TIME_BIN_max,
                      ymin = theor_LCI, ymax = theor_UCI,
                      group = PI, fill = CI_label),
                  col = NA, linewidth = 1, alpha = 0.3) +
        geom_rect(data = ci_data %>% filter(PI == "median"),
                  aes(xmin = TIME_BIN_min, xmax = TIME_BIN_max,
                      ymin = theor_LCI, ymax = theor_UCI,
                      group = PI, fill = CI_label),
                  col = NA, linewidth = 1, alpha = 0.3)
    }
  }

  return(plot_i)
}

create_vpc_plot <- function(pred_int_median_emp_mod,
                            pred_int_ci_mod,
                            data_obs,
                            pred_int_ci,
                            data_i_var,
                            theor_perc,
                            emp_perc,
                            theor_percCI,
                            interpolation,
                            dt_obs_fl,
                            strat_by_dose,
                            log_y,
                            pred.corr,
                            lab_y,
                            lab_x,
                            legend_fl) {

  # Base plot aesthetics
  p_char <- list(
    theme_bw(),
    theme(
      panel.background = element_rect(fill = "transparent", colour = "black"),
      strip.text = element_text(size = 18, colour = "black"),
      plot.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 18)
    ),
    scale_fill_manual(values = c("royalblue", "firebrick3")),
    scale_linetype_manual(values = c(
      `empirical percentiles` = "solid",
      `theoretical percentiles` = "dashed"
    )),
    scale_x_continuous(name = lab_x)
  )

  # Legend settings
  p_char_leg_T <- list(
    theme(legend.position = c(0.8, 0.85)),
    theme(legend.box.background = element_rect(fill = "white", color = "black", linetype = "solid")),
    theme(legend.box.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")),
    theme(legend.text = element_text(size = 18)),
    theme(legend.title = element_blank()),
    theme(legend.key.size = unit(1, "cm")),
    theme(legend.margin = margin(0, 0, 0, 0, "cm")),
    #theme(legend.margin = unit(0, "cm")),
    labs(lty = "", fill = "", color = "", shape = "")
  )

  p_char_leg_F <- list(theme(legend.position = "none"))

  # Initialize plot
  vpc_plot <- ggplot()

  # Add percentile lines and confidence intervals
  if (theor_perc && emp_perc) {
    vpc_plot <- add_percentile_lines(vpc_plot, pred_int_median_emp_mod,
                                     pred_int_ci_mod, theor_percCI, interpolation)
  } else if (theor_perc) {
    data_filtered <- pred_int_median_emp_mod %>%
      filter(TYPE == "theoretical percentiles")
    vpc_plot <- add_percentile_lines(vpc_plot, data_filtered,
                                     pred_int_ci_mod, theor_percCI, interpolation)
  } else if (emp_perc) {
    data_filtered <- pred_int_median_emp_mod %>%
      filter(TYPE == "empirical percentiles")
    vpc_plot <- add_percentile_lines(vpc_plot, data_filtered,
                                     pred_int_ci_mod, theor_percCI, interpolation)
  }

  # Add stratification by dose
  if (!is.null(strat_by_dose)) {
    vpc_plot <- vpc_plot + facet_grid(~COV_add)
  }

  # Add observed data points
  if (dt_obs_fl) {
    vpc_plot <- vpc_plot +
      geom_point(data = data_obs, aes(x = TIME, y = DV, shape = DT_label),
                 col = "royalblue4", size = 1.5)
  }

  # Y-axis scaling
  y_label <- paste0(if_else(pred.corr, "Prediction corrected ", ""), lab_y)

  x_range <- c(
    min(c(pred_int_ci$TIME_BIN, data_i_var$TIME)),
    max(c(pred_int_ci$TIME_BIN, data_i_var$TIME))
  )

  if (log_y) {
    vpc_plot <- vpc_plot +
      scale_y_log10(
        name = y_label,
        breaks = 10^seq(-2, 4, 1),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      coord_cartesian(xlim = x_range)
  } else {
    vpc_plot <- vpc_plot +
      scale_y_continuous(name = y_label) +
      coord_cartesian(xlim = x_range)
  }

  # Apply styling
  vpc_plot <- vpc_plot + p_char

  # Apply legend settings
  if (legend_fl) {
    vpc_plot <- vpc_plot + p_char_leg_T
  } else {
    vpc_plot <- vpc_plot + p_char_leg_F
  }

  return(vpc_plot)
}

fun_Bin_smrg <- function(ds, n_bins, method = c("kmeans", "ntile", "equal_x", "custom")){

  #browser()

  if(method == "kmeans"){

    if(n_bins > length(unique(ds$TIME))){
      n_bins <- length(unique(ds$TIME))
    }

    set.seed(123)

    #browser()

    ds_bin <- ds %>%
      mutate(BIN = as.factor(kmeans(TIME, n_bins)$cluster))

    sorted_bins <- ds_bin %>%
      group_by(BIN) %>%
      summarise_at(vars(TIME), list(min_BIN_time = ~min(., na.rm=T))) %>%
      ungroup() %>%
      select(BIN, min_BIN_time) %>%
      unique() %>%
      arrange(min_BIN_time) %>%
      mutate(BIN_sort = 1:nrow(.)) %>%
      select(-min_BIN_time)

    ds_bin <- ds_bin %>%
      left_join(sorted_bins, by='BIN') %>%
      select(-BIN) %>%
      rename(BIN=BIN_sort)

  } else if(method == "ntile"){
    ds_bin <- ds %>%
      mutate(BIN = as.factor(ntile(TIME, n_bins)))
  } else if(method == "equal_x"){


    #browser()

    #browser()
    n_bins <- seq(min(ds$TIME, na.rm = T),
                max(ds$TIME, na.rm = T),
                length.out = n_bins + 1)
    n_bins <- n_bins[-1]

    time_bin <- ds %>%
      select(TIME) %>%
      filter(!is.na(TIME)) %>%
      mutate(BIN = 0) %>%
      unique()

    time_bin <- time_bin %>%
      mutate(nrow = 1:nrow(.)) %>%
      group_by(nrow) %>%
      mutate(BIN = min((1:length(n_bins))[TIME <= n_bins], na.rm = T)) %>%
      ungroup() %>%
      select(-nrow)

    # for (i in 1:nrow(time_bin)){
    #   time_bin[i,]$BIN <- min((1:length(n_bins))[time_bin[i,]$TIME <= n_bins], na.rm = T)
    # }

    time_bin <- time_bin %>%
      mutate(BIN = as.factor(BIN))

    ds_bin <- ds %>%
      left_join(time_bin, by = "TIME")

  } else if(method == "custom"){
    n_bins <- c(n_bins, max(ds$TIME, na.rm = T))

    time_bin <- ds %>%
      select(TIME) %>%
      filter(!is.na(TIME)) %>%
      mutate(BIN = 0)

    for (i in 1:nrow(time_bin)){
      time_bin[i,]$BIN <- min((1:length(n_bins))[time_bin[i,]$TIME <= n_bins], na.rm = T)
    }

    time_bin <- time_bin %>%
      mutate(BIN = as.factor(BIN))

    ds_bin <- ds %>%
      left_join(time_bin, by = "TIME")
  }

  ds_bin <- ds_bin %>%
    group_by(BIN) %>%
    mutate_at(vars(TIME), list(TIME_BIN = ~mean(., na.rm = T))) %>%
    ungroup()

  bin_stats <- ds_bin %>%
    group_by(BIN) %>%
    summarise(TIME_BIN_min = min(TIME),
              TIME_BIN_max = max(TIME)) %>%
    mutate(TIME_BIN_max_prev = lag(TIME_BIN_max),
           TIME_BIN_min_next = lead(TIME_BIN_min)) %>%
    ungroup() %>%
    mutate(TIME_BIN_min = (TIME_BIN_min + TIME_BIN_max_prev) / 2,
           TIME_BIN_max = (TIME_BIN_max + TIME_BIN_min_next) / 2) %>%
    select(-TIME_BIN_max_prev, -TIME_BIN_min_next)

  bin_stats$TIME_BIN_min[is.na(bin_stats$TIME_BIN_min)] <- -Inf
  bin_stats$TIME_BIN_max[is.na(bin_stats$TIME_BIN_max)] <- Inf

  ds_bin <- ds_bin %>%
    left_join(bin_stats, by = "BIN")

  return(ds_bin)
}

#' Visual Predictive Check (VPC) Function
#'
#' Generates VPC plots to assess model performance by comparing observed data
#' with prediction intervals derived from simulated data.
#'
#' @inheritParams sg_dummy
#' @param ds_sim A data frame with simulated data containing columns:
#'   - ID: simulation replicate identifier
#'   - TIME: time points
#'   - VAR: output variable name
#'   - VALUE: simulated values
#' @param data_i A data frame with observed data containing columns:
#'   - ID: patient identifier
#'   - TIME: time points
#'   - DV: observed dependent variable values
#'   - DVID: dependent variable identifier
#' @param output_names A data frame mapping `VAR` values (from `ds_sim`)
#'  to `DVID` values (from `data_i`). Must contain columns:
#'   - output: character, `VAR` values (e.g., "Cc", "Cp")
#'   - dvid: numeric, corresponding DVIDs (e.g., 1, 2)
#' @param lab_x String. X-axis label. Default is "TIME, h".
#' @param lab_y String. Y-axis label. Default is "DV".
#' @param method Character. Binning method: "kmeans", "equal", "quantile" (default: "kmeans").
#' @param n_bins Integer. Number of bins to use in the histogram. Default is 10.
#' @param interpolation Logical. Use line interpolation vs rectangles (default: FALSE).
#' @param strat_by_dose Character. Variable name for dose stratification (default: NULL).
#'
#' @return List of ggplot objects, one for each output variable
#' @details
#' For now, the model in the example is NOT the generalized model object (GMO),
#' as the parameters in generalized fit object are backtransformed, and,
#' therefore, not transformed in the model. This issue will be fixed in the next
#' versions of SimuRg package
#'
#' @examples
#' \donttest{
#' library(rxode2)
#' fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
#' mod <- rxode2::rxode2({
#'   ka_pop = 0.1;
#'   Vd_pop = 10;
#'   Cl_pop = 0.5;
#'
#'   omega_ka = 0;
#'   omega_Vd = 0;
#'   omega_Cl = 0;
#'
#'   Cc_b = 0;
#'   ka_tv = exp(ka_pop);
#'   Vd_tv = exp(Vd_pop);
#'   Cl_tv = exp(Cl_pop);
#'
#'   ka = ka_tv * exp(omega_ka);
#'   Vd = Vd_tv * exp(omega_Vd);
#'   Cl = Cl_tv * exp(omega_Cl);
#'
#'   Cc = Ac / Vd;
#'
#'   Ad(0) = 0;
#'   Ac(0) = 0;
#'
#'   d/dt(Ad) = -ka * Ad;
#'   d/dt(Ac) = ka * Ad - Cl * Cc;
#'
#'   Cc_ResErr = Cc * (1 + Cc_b);
#' })
#'
#' sim_data <- sg_vpc_sim(fpath_i, model=mod, outputs = "Cc_ResErr")
#' outp_nms <- data.frame(dvid = 1, output = "Cc")
#' vpc_plots <- sg_vpc_vis(
#'   ds_sim = sim_data,
#'   data_i = warfarin,
#'   output_names = outp_nms,
#'   lab_x = "Time (hours)",
#'   lab_y = "Concentration (ng/mL)",
#'   n_bins = 8,
#'   pred.corr = TRUE
#' )
#' }
#' @import dplyr tidyr ggplot2 purrr stringr
#' @export
sg_vpc_vis <- function(ds_sim,
                       data_i,
                       output_names,
                       time_col = "TIME",
                       dv_col = "DV",
                       log_y = F,
                       piLow = 0.10,
                       piUp = 0.90,
                       ciLow = 0.025,
                       ciUp = 0.975,
                       pred.corr = FALSE,
                       lab_y = "DV",
                       lab_x = "Time, h",
                       theor_perc = TRUE,
                       theor_percCI = TRUE,
                       emp_perc = TRUE,
                       dt_obs_fl = FALSE,
                       legend_fl = FALSE,
                       n_bins = 10,
                       method = "kmeans",
                       interpolation = T,
                       strat_by_dose = NULL) {

  # Check input dataset for required columns
  required_ds_sim_cols <- c("sim.id", "ID", "TIME", "VAR", "VALUE")
  missing_ds_sim_cols <- setdiff(required_ds_sim_cols, colnames(ds_sim))
  if (length(missing_ds_sim_cols) > 0) {
    stop(
      "ds_sim is missing required columns: ",
      paste(missing_ds_sim_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Validate output_names: must be provided and contain required columns
  if (is.null(output_names)) {
    stop("output_names must be provided and map VAR to DVID", call. = FALSE)
  }
  required_output_cols <- c("output", "dvid")
  missing_output_cols <- setdiff(required_output_cols, colnames(output_names))
  if (length(missing_output_cols) > 0) {
    stop(
      "output_names is missing required columns: ",
      paste(missing_output_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Check if the method is valid
  valid_methods <- c("kmeans", "ntile", "equal_x", "custom")
  if (!method %in% valid_methods) {
    cat("Invalid binning method. Valid methods are: ", paste(valid_methods, collapse = ", "), "\n",
        "Using default method: kmeans\n")
    method <- "kmeans"
  }

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
      rename(TIME = all_of(time_col), DV = all_of(dv_col))

    # Store original observed data for overlay
    data_obs <- data_i_var %>%
      mutate(DT_label = "Data")

    # ===== TIME BINNING =====
    data_i_var <- fun_Bin_smrg(data_i_var, n_bins = n_bins, method = method)

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
    colnames(emp_int) <- str_remove(colnames(emp_int), "DV_")

    # Prediction intervals from simulated data (per simulation)
    pred_int <- ds_sim_var %>%
      mutate(ID = as.numeric(ID)) %>%
      group_by(across(all_of(c(vec_group, "ID")))) %>%
      summarise(across(VALUE, summ_pi), .groups = "drop")
    colnames(pred_int) <- str_remove(colnames(pred_int), "VALUE_")
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
      log_y = log_y,
      pred.corr = pred.corr,
      lab_y = lab_y,
      lab_x = lab_x,
      legend_fl = legend_fl
    )

    return(vpc_plot)
  })

  names(vpc_plots) <- unique_vars
  return(vpc_plots)
}

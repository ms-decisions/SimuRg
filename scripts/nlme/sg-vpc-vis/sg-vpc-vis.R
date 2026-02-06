library(tidyverse)
library(rxode2)
library(devtools)
source("scripts/nlme/5.1.1-sg-vpc-sim/sg-vpc-sim.R")
load("data/warfarin.rda")
#' Add Percentile Lines to VPC Plot
#'
#' Internal helper function
#'
#' @keywords internal

add_percentile_lines <- function(plot, line_data, ci_data, show_ci, interpolation) {

  if (interpolation) {
    # Use smooth lines
    plot <- plot +
      geom_line(data = line_data,
                aes(x = TIME_BIN, y = OBS, group = interaction(PI, TYPE), lty = TYPE),
                col = "royalblue", size = 0.8)

    if (show_ci) {
      plot <- plot +
        geom_ribbon(data = ci_data %>% filter(PI != "median"),
                    aes(x = TIME_BIN, ymin = theor_LCI, ymax = theor_UCI,
                        group = PI, fill = CI_label),
                    col = NA, size = 1, alpha = 0.3) +
        geom_ribbon(data = ci_data %>% filter(PI == "median"),
                    aes(x = TIME_BIN, ymin = theor_LCI, ymax = theor_UCI,
                        group = PI, fill = CI_label),
                    col = NA, size = 1, alpha = 0.3)
    }

  } else {
    # Use step rectangles
    plot <- plot +
      geom_line(data = line_data,
                aes(x = TIME_BIN, y = OBS, group = interaction(PI, TYPE), lty = TYPE),
                col = "royalblue", size = 0.8)

    if (show_ci) {
      plot <- plot +
        geom_rect(data = ci_data %>% filter(PI != "median"),
                  aes(xmin = TIME_BIN_min, xmax = TIME_BIN_max,
                      ymin = theor_LCI, ymax = theor_UCI,
                      group = PI, fill = CI_label),
                  col = NA, size = 1, alpha = 0.3) +
        geom_rect(data = ci_data %>% filter(PI == "median"),
                  aes(xmin = TIME_BIN_min, xmax = TIME_BIN_max,
                      ymin = theor_LCI, ymax = theor_UCI,
                      group = PI, fill = CI_label),
                  col = NA, size = 1, alpha = 0.3)
    }
  }

  return(plot)
}
#' Create VPC Plot
#'
#' Internal function to generate the VPC ggplot object
#'
#' @param ... Multiple data frames and plot parameters
#' @return ggplot object
#' @keywords internal

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
                            logy,
                            pred.corr,
                            name_y,
                            name_x,
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
    scale_x_continuous(name = name_x)
  )

  # Legend settings
  p_char_leg_T <- list(
    theme(legend.position = c(0.8, 0.85)),
    theme(legend.box.background = element_rect(fill = "white", color = "black", linetype = "solid")),
    theme(legend.box.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")),
    theme(legend.text = element_text(size = 18)),
    theme(legend.title = element_blank()),
    theme(legend.key.size = unit(1, "cm")),
    theme(legend.margin = unit(0, "cm")),
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
  y_label <- paste0(if_else(pred.corr, "Prediction corrected ", ""), name_y)

  x_range <- c(
    min(c(pred_int_ci$TIME_BIN, data_i_var$TIME)),
    max(c(pred_int_ci$TIME_BIN, data_i_var$TIME))
  )

  if (logy) {
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
fun_Bin_smrg <- function(ds, bins, method = c("kmeans", "ntile", "equal_x", "custom")){

  #browser()

  if(method == "kmeans"){

    if(bins > length(unique(ds$TIME))){
      bins <- length(unique(ds$TIME))
    }

    set.seed(123)

    #browser()

    ds_bin <- ds %>%
      mutate(BIN = as.factor(kmeans(TIME, bins)$cluster))

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
      mutate(BIN = as.factor(ntile(TIME, bins)))
  } else if(method == "equal_x"){


    #browser()

    print(c(1))
    #browser()
    bins <- seq(min(ds$TIME, na.rm = T),
                max(ds$TIME, na.rm = T),
                length.out = bins + 1)
    print(c(2))
    bins <- bins[-1]
    print(c(3))

    time_bin <- ds %>%
      select(TIME) %>%
      filter(!is.na(TIME)) %>%
      mutate(BIN = 0) %>%
      unique()

    print(c(4))

    time_bin <- time_bin %>%
      mutate(nrow = 1:nrow(.)) %>%
      group_by(nrow) %>%
      mutate(BIN = min((1:length(bins))[TIME <= bins], na.rm = T)) %>%
      ungroup() %>%
      select(-nrow)

    print(c(5))
    # for (i in 1:nrow(time_bin)){
    #   time_bin[i,]$BIN <- min((1:length(bins))[time_bin[i,]$TIME <= bins], na.rm = T)
    # }

    time_bin <- time_bin %>%
      mutate(BIN = as.factor(BIN))

    print(c(6))

    ds_bin <- ds %>%
      left_join(time_bin, by = "TIME")

    print(c(7))
  } else if(method == "custom"){
    bins <- c(bins, max(ds$TIME, na.rm = T))

    time_bin <- ds %>%
      select(TIME) %>%
      filter(!is.na(TIME)) %>%
      mutate(BIN = 0)

    for (i in 1:nrow(time_bin)){
      time_bin[i,]$BIN <- min((1:length(bins))[time_bin[i,]$TIME <= bins], na.rm = T)
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

sg_vpc_vis <- function(ds_sim,
                       data_i,
                       output_names,
                       x = "TIME",
                       y = "DV",
                       logy = F,
                       piLow = 0.10,
                       piUp = 0.90,
                       ciLow = 0.025,
                       ciUp = 0.975,
                       pred.corr = FALSE,
                       name_y = "DV",
                       name_x = "Time, h",
                       theor_perc = TRUE,
                       theor_percCI = TRUE,
                       emp_perc = TRUE,
                       dt_obs_fl = FALSE,
                       legend_fl = FALSE,
                       bins = 10,
                       method = "kmeans",
                       interpolation = T,
                       strat_by_dose = NULL) {

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
    data_i_var <- fun_Bin_smrg(data_i_var, bins = bins, method = method)

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
outp_nms <- data.frame(dvid = 1, output = "Cc")
plt1 <- sg_vpc_vis(res, warfarin, outp_nms)

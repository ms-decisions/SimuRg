#' Visualize Local Sensitivity Analysis Results
#'
#' @description
#' This function creates plots to visualize the results of a local sensitivity analysis:
#' (1) a family of curves plot showing how model outputs vary over time for different parameter perturbations,
#' and (2) a tornado plot summarizing the relative influence of each parameter on selected summary metrics
#' (e.g., output value at a chosen time or Cmax).
#'
#' @param sens_data A data frame containing sensitivity analysis results from sg_localsens_sim()
#' @param tornado_time Numeric value specifying the time point (in the same units as `TIME`)
#' at which to evaluate the model output for the tornado plot. Defaults to the maximum time in `sens_data` if not specified.
#' @param metrics Character vector specifying which metrics to include in the tornado plot.
#'   Supported options are `"value"` (value at `tornado_time`) and `"cmax"` (maximum output). Defaults to `"value"`.
#' @param ref_data Optional data frame providing reference (baseline) output values,
#'   formatted similarly to `sens_data`. If `NULL`, the function uses the mid-range parameter value as the baseline.
#' @param log_scale Logical. If `TRUE`, the y-axis of the family-of-curves plot is displayed on a log scale.
#' @param color_low Character. Color used for the lower parameter bound (default: "#1f77b4").
#' @param color_high Character. Color used for the upper parameter bound (default: "#ff7f0e").
#' @param facet_scales Character string passed to `facet_grid()`, defining how y-axis scales behave across facets (default = `"free"`).
#'
#' @return A named list containing:
#'   \describe{
#'     \item{`family_of_curves`}{A `ggplot2` object showing parameter sensitivity curves over time.}
#'     \item{`tornado`}{A `ggplot2` object showing the tornado plot of relative sensitivity.}
#'   }
#'
#' @examples
#' \donttest{
#' library(rxode2)
#' library(tibble)
#' mod_ex <- RxODE({
#'   # Doses in mg
#'   # Time in hours
#'   ### Parameter values
#'   # Typical
#'   POPCL = 5;
#'   POPVC = 180;
#'   POPQ = 7;
#'   POPVP = 52;
#'   POPKTR = 6;
#'   FBIOPAR = 1;

#'   # Covariate coefficients
#'   POPIGFRCOV = 0.8;
#'   POPPATCOVCLGR1 = 0.85;
#'   POPPATCOVCLGR2 = 0.85;
#'   POPPATCOVCLGR3 = 0.85;
#'   POPPATCOVCLGR4 = 0.85;

#'   # Random effects
#'   PPVCL = 0;
#'   PPVVC = 0;
#'   PPVKTR = 0;

#'   # Residual error
#'   RUV = 0;

#'   ### Covariates
#'   PATCOVCL = 1
#'   if(POPN == 1){PATCOVCL = POPPATCOVCLGR1}
#'   if(POPN == 2){PATCOVCL = POPPATCOVCLGR2}
#'   if(POPN == 3){PATCOVCL = POPPATCOVCLGR3}
#'   if(POPN == 4){PATCOVCL = POPPATCOVCLGR4}
#'   IGFRCOV = (IGFR/112)^POPIGFRCOV

#'   ## Parameters
#'   CL = POPCL * IGFRCOV * PATCOVCL * exp(PPVCL);
#'   VC = POPVC * exp(PPVVC);
#'   Q = POPQ;
#'   VP = POPVP;
#'   KTR = POPKTR * exp(PPVKTR);

#'   ### Explicit functions
#'   Cc = Ac/VC;                 # nmol/L
#'   Cp = Ap/VP;                 # nmol/L

#'   ### Initial conditions
#'   At1(0) = 0;         # mg
#'   At2(0) = 0;         # mg
#'   At3(0) = 0;         # mg
#'   At4(0) = 0;         # mg
#'   At5(0) = 0;         # mg
#'   At6(0) = 0;         # mg
#'   Ad(0) = 0;          # mg
#'   Ac(0) = 0;          # mg
#'   Ap(0) = 0;          # mg

#'   ### ODEs
#'   d/dt(At1) = - KTR*At1;
#'   d/dt(At2) = KTR*At1 - KTR*At2;
#'   d/dt(At3) = KTR*At2 - KTR*At3;
#'   d/dt(At4) = KTR*At3 - KTR*At4;
#'   d/dt(At5) = KTR*At4 - KTR*At5;
#'   d/dt(At6) = KTR*At5 - KTR*At6;
#'   d/dt(Ad) = KTR*At6 - KTR*Ad;
#'   d/dt(Ac) = KTR*Ad - CL*Cc - Q*(Cc - Cp);
#'   d/dt(Ap) = Q*(Cc - Cp);

#'   FBIO = FBIOPAR
#'   f(At1) = FBIO*1000000/505;
#'   CHECKRUV = RUV;
#'   Cc_ResErr = Cc + RUV*Cc;
#' })
#' et_base <- tribble(
#'   ~id, ~time, ~evid, ~cmt, ~amt, ~addl, ~ii, ~IGFR, ~POPN,
#'   1,   0,     1,     1,    10,   2,     24,  112,   1
#' )
#' sens_loc <- sg_localsens_sim(
#'   model = mod_ex,
#'   params = c("POPCL", "POPVC"),
#'   stimes = seq(0, 168, 0.1),
#'   output = "Cc",
#'   perc = 0.9,
#'   cov = c("IGFR", "POPN"),
#'   et = et_base
#' )
#' plots <- sg_localsens_vis(
#'   sens_data = sens_loc,
#'   tornado_time = 168,
#'   metrics = c("value", "cmax"),
#'   log_scale = TRUE,
#'   color_low = "#4575b4",
#'   color_high = "#d73027",
#'   facet_scales = "free"
#' )
#'
#' # Display plots
#' print(plots$family_of_curves)
#' print(plots$tornado)
#' }
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom grid unit
#' @export

sg_localsens_vis <- function(sens_data,
                             tornado_time = NULL,
                             metrics = "value",
                             ref_data = NULL,
                             log_scale = FALSE,
                             color_low = "#1f77b4",
                             color_high = "#ff7f0e",
                             facet_scales = "free") {

  output_var <- unique(sens_data$VAR)

  if(is.null(tornado_time)) {
    tornado_time <- max(sens_data$TIME)
  }

  foc_plot <- ggplot(data = sens_data) +
    geom_line(aes(x = TIME, y = VALUE, group = interaction(PARVAL, ID),
                  color = PARVAL_NORM), alpha = 0.7, linewidth = 0.6) +
    facet_grid(PARNAME ~ ., scales = facet_scales) +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13),
      legend.text = element_text(size = 10),
      strip.text = element_text(size = 11),
      plot.title = element_text(size = 14, hjust = 0.5),
      legend.key.size = unit(0.3, "cm")
    ) +
    labs(
      x = "Time",
      y = output_var,
      color = "Parameter\nvalue"
    ) +
    scale_color_gradient(
      low = color_low,
      high = color_high,
      breaks = c(0, 1),
      labels = c("Lower bound", "Upper bound"),
      guide = guide_colorbar(
        title.theme = element_text(size = 10),
        label.theme = element_text(size = 9),
        barwidth = unit(0.5, "cm"),
        barheight = unit(2, "cm")
      )
    )

  if(log_scale) {
    foc_plot <- foc_plot + scale_y_log10()
  }


  # Tornado plot
  metrics_data <- sens_data %>%
    group_by(PARNAME, PARVAL, PARVAL_NORM) %>%
    summarise(
      value_at_time = VALUE[which.min(abs(TIME - tornado_time))],
      cmax = max(VALUE, na.rm = TRUE),
      .groups = 'drop'
    )

  # Get reference values (either from ref_data or middle parameter value)
  if(!is.null(ref_data)) {
    ref_metrics <- ref_data %>%
      summarise(
        value_ref = VALUE[which.min(abs(TIME - tornado_time))],
        cmax_ref = max(VALUE, na.rm = TRUE)
      )
  } else {
    ref_metrics <- metrics_data %>%
      group_by(PARNAME) %>%
      filter(abs(PARVAL_NORM - 0.5) == min(abs(PARVAL_NORM - 0.5))) %>%
      summarise(
        value_ref = first(value_at_time),
        cmax_ref = first(cmax),
        .groups = 'drop'
      )
  }

  tornado_data <- metrics_data %>%
    group_by(PARNAME) %>%
    filter(PARVAL_NORM %in% c(0, 1)) %>%
    mutate(
      PAR_group = if_else(PARVAL_NORM == 0, "L", "U")
    ) %>%
    select(PARNAME, PARVAL, PAR_group, value_at_time, cmax) %>%
    pivot_longer(
      cols = c(value_at_time, cmax),
      names_to = "METRIC",
      values_to = "VALUE"
    ) %>%
    # Join with reference values
    left_join(ref_metrics, by = "PARNAME") %>%
    mutate(
      value_cfb = if_else(METRIC == "value_at_time",
                          (VALUE - value_ref) * 100 / value_ref,
                          (VALUE - cmax_ref) * 100 / cmax_ref),
      METRIC = case_when(
        METRIC == "value_at_time" ~ paste0("Value at t=", tornado_time),
        METRIC == "cmax" ~ "Cmax"
      )
    )

  metric_labels <- c(
    "value" = paste0("Value at t=", tornado_time),
    "cmax" = "Cmax"
  )

  selected_metrics <- metric_labels[metrics]
  tornado_data <- tornado_data %>%
    filter(METRIC %in% selected_metrics)

  tornado_plot <- ggplot(tornado_data,
                         aes(x = reorder(PARNAME, abs(value_cfb)),
                             y = value_cfb,
                             fill = PAR_group)) +
    geom_bar(stat = "identity", position = "identity", alpha = 0.8) +
    facet_grid(METRIC ~ ., scales = facet_scales) +
    coord_flip() +
    labs(
      y = "Change from baseline, %",
      x = "Parameter"
    ) +
    scale_fill_manual(
      values = c("U" = color_high, "L" = color_low),
      labels = c("U" = "Upper bound", "L" = "Lower bound"),
      name = "Parameter\nvalue"
    ) +
    geom_hline(yintercept = 0, linetype = "longdash", alpha = 0.7) +
    theme_minimal() +
    theme(
      axis.title.y = element_blank(),
      strip.text = element_text(size = 11),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "bottom"
    )

   return(list(
     family_of_curves = foc_plot,
     tornado = tornado_plot
   ))
}

#' Generate prediction distribution visualizations from sg_predist_sim() simulation results
#' @inheritParams sg_dummy
#' @param ds_sim [GSO]. Data frame with simulation results in long format
#'   (columns `ID`, `TIME`, `VAR`, `VALUE`), typically from [sg_predist_sim()].
#' @param dt_obs_fl Logical. If TRUE, observed data points are overlaid on the simulation plots. Default is FALSE.
#' @param logy Logical. If TRUE, y-axis is displayed on a logarithmic scale. Default is FALSE.
#' @param legend_fl Logical. If FALSE, the legend is hidden. Default is TRUE.
#' @param pred_interval Character. Prediction interval to display. Options are "95%", "90%", "80%", "50%". Default is "80%".
#' @param lab_y Character. Y-axis label. If `NULL` (default), each plot is labelled
#' @param lab_x Character. X-axis label. Default is `"TIME"`
#'   with its own output variable (the `VAR` value, with any `_ResErr` suffix removed).
#'
#' @returns A list of ggplot objects, one per output variable in the simulation dataset.
#' @seealso [sg_predist_sim()], [sg_vpc_vis()]
#' @examples
#' \donttest{
#' library(rxode2)
#' fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
#' mod_fin <- rxode2({
#'  # Differential equations
#'  ka = ka_pop*exp(omega_ka)
#'  V = V_pop
#'  Cl = Cl_pop*exp(omega_Cl)
#'  d/dt(Ad) = -ka * Ad
#'  d/dt(Ac) = ka * Ad - Cl/V * Ac
#'  # Concentration calculations
#'  Cc = Ac / V
#' })
#' predist_sim <- sg_predist_sim(fpath_i, model =mod_fin, outputs = "Cc")
#' predist_vis <- sg_predist_vis(fpath_i,predist_sim)
#' predist_vis
#' }
#' @import ggplot2
#' @import dplyr
#' @importFrom purrr map
#' @importFrom scales pretty_breaks trans_format math_format
#' @export
sg_predist_vis <- function(fpath_i,
                           ds_sim,
                           time_col = "TIME",
                           lab_x = "TIME", lab_y = NULL,
                           dt_obs_fl = FALSE, logy = FALSE, legend_fl = TRUE,
                           pred_interval = "80%") {

  MSDcol <- c("#1a1866", "#f2b93b", "#b73b58", "#a2d620", "#5839bb", "#9c4ec7", "#3a6eba", "#efdd3c", "#69686d")

  # ds_sim from sg_predist_sim() is already in long format (ID, TIME, VAR, VALUE);
  # no reshaping is needed (see sg_sim() output specification).
  required_ds_sim_cols <- c("ID", "TIME", "VAR", "VALUE")
  missing_ds_sim_cols <- setdiff(required_ds_sim_cols, colnames(ds_sim))
  if (length(missing_ds_sim_cols) > 0) {
    stop(
      "ds_sim is missing required columns: ",
      paste(missing_ds_sim_cols, collapse = ", "),
      call. = FALSE
    )
  }

  obj <- read_smrg_obj(fpath_i)

  data_fin.noex <-  obj$SDTAB %>% filter(MDV != 1) %>% select(-MDV)
  data_fin.noex$TIME  <- data_fin.noex[[time_col]]

  lower_quantile <- switch(
    pred_interval,
    "95%" = 0.025,
    "90%" = 0.05,
    "80%" = 0.10,
    "50%" = 0.25,
    stop("pred_interval must be one of 95%, 90%, 80%, 50%")
  )
  upper_quantile <- 1 - lower_quantile

  local_funSum <- list(
    mean   = ~mean(.),
    median = ~median(.),
    min    = ~min(.),
    max    = ~max(.),
    sd     = ~sd(.),
    se     = ~sd(.)/sqrt(n()),
    L_Q    = ~quantile(., lower_quantile),
    H_Q    = ~quantile(., upper_quantile)
  )

  # ---- Shared plot styling (from sg_vpc_vis) ----
  p_char <- list(
    theme_bw(),
    theme(panel.background = element_rect(fill = "transparent", colour = "black")),
    theme(strip.text  = element_text(size = 18, colour = "black")),
    theme(plot.title  = element_text(size = 18, face = "bold")),
    theme(axis.text   = element_text(size = 14)),
    theme(axis.title  = element_text(size = 18)),
    scale_color_manual(values = c(MSDcol[1])),
    scale_fill_manual(values = c(MSDcol[1]))
  )

  p_char_leg_T <- list(
    theme(legend.position = c(0.8, 0.85)),
    theme(legend.box.background = element_rect(fill = "white", color = "black", linetype = "solid")),
    theme(legend.box.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")),
    theme(legend.text = element_text(size = 18)),
    theme(legend.title = element_blank()),
    theme(legend.key.size = unit(1, "cm")),
    theme(legend.margin = margin(0, 0, 0, 0, "cm")),
    labs(fill = "", color = "")
  )

  p_char_leg_F <- list(theme(legend.position = "none"))

  p_list <- unique(ds_sim$VAR) %>% map(function(v){

    ds_sim_v <- ds_sim %>% filter(VAR == v)

    # Default the y-axis label to the output variable itself (one plot per VAR);
    # an explicit `lab_y` overrides it.
    lab_y <- if (is.null(lab_y)) sub("_ResErr$", "", v) else lab_y

    ds_sim_sum <- ds_sim_v %>%
      group_by(TIME, VAR) %>%
      summarise(across(VALUE, local_funSum, .names = "{.fn}"), .groups = "drop") %>%
      mutate(
        label_median = "Median",
        label_band   = paste0(pred_interval, " Prediction Interval")
      )

    p <- ggplot() +
      geom_ribbon(
        data = ds_sim_sum,
        aes(x = TIME, ymin = L_Q, ymax = H_Q, fill = label_band),
        col = NA, alpha = 0.2
      ) +
      geom_line(
        data = ds_sim_sum,
        aes(x = TIME, y = median, color = label_median),
        linewidth = 0.8, alpha = 0.8
      ) +
      scale_x_continuous(name = lab_x, breaks = pretty_breaks()) +
      p_char

    if (dt_obs_fl) {
      data_obs_v <- data_fin.noex %>%
        select(ID, TIME, all_of(v)) %>%
        rename(DV = all_of(v)) %>%
        mutate(label = "Data")
      p <- p +
        geom_point(
          data = data_obs_v,
          aes(x = TIME, y = DV, shape = label),
          color = "royalblue4", size = 1.5
        )
    }

    if (logy) {
      p <- p +
        scale_y_log10(
          name = lab_y,
          breaks = 10^seq(-10, 10, 1),
          labels = scales::trans_format("log10", scales::math_format(10^.x))
        )
    } else {
      p <- p +
        scale_y_continuous(name = lab_y, breaks = pretty_breaks())
    }

    if (legend_fl) {
      p <- p + p_char_leg_T
    } else {
      p <- p + p_char_leg_F
    }

    p
  })

  names(p_list) <- unique(ds_sim$VAR)
  return(p_list)
}

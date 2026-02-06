## Author: Ugolkov Yaroslav
## First created: 2025-09-05
## Description: sg-gof-res-dist and its helper functions
## Keywords: SimuRg, sg-gof-res-dist, goodness-of-fit
#'
#' Plot distribution or QQ-plot of residuals from NONMEM/Simurg run
#'@description
#' This function visualizes residuals from estimation results
#' stored in a Simurg output object. The residuals can be plotted either as
#' histograms with overlaid normal density curves, or as QQ-plots to assess
#' normality.
#'
#' @inheritParams sg_dummy
#' @param res_type Character vector. One or several types of residuals to plot.
#'   Values must correspond to column name(s) in `smrg_obj$SDTAB`, e.g.
#'   `"RES"`, `"IWRES"`, `"IRES"`. If multiple residual types are provided,
#'   plots will be generated for each of them.
#' @param ndist Logical. If `TRUE`, overlays the corresponding normal density
#'   curve on the histogram (default: TRUE).

#'
#' @return A `ggplot` object.
#' @examples
#' \dontrun{
#' # Plot the distribution of individual weighted residuals (IWRES)
#' # as histograms with overlaid normal curves for each DVID
#' sg_gof_res_dist(fpath_i = "PK.RData", res_type = "IWRES")
#'
#' # Generate QQ-plots to visually assess normality of standard residuals (RES)
#' # stratified by DVID
#' sg_gof_res_dist(fpath_i = "PK.RData", res_type = "RES", plot_type = "QQ")
#'
#' # Multiple residual types can be specified (e.g., IWRES and IRES);
#' # the function will produce plots for each residual type across all DVID groups
#' sg_gof_res_dist(fpath_i = "PK.RData", res_type = c("IWRES", "IRES"))
#' }
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @export


sg_gof_res_dist <- function(fpath_i, res_type = 'RES', n_bins = 30, ndist = T, plot_type = 'DIST'){
  smrg_obj <- read_smrg_obj(fpath_i)
  if (is.null(smrg_obj$SDTAB)) {
    stop("sg_fit object must contain SDTAB component")
  }

  sdtab <- smrg_obj$SDTAB

  if (is.data.frame(sdtab) && nrow(sdtab) == 0) {
    stop("SDTAB is empty (no rows)")
  }
  if (is.list(sdtab) && length(sdtab) == 0) {
    stop("SDTAB is empty (no elements)")
  }
  if (is.data.frame(sdtab)) {
    sdtab_i <- sdtab
  } else if (is.list(sdtab)) {
    sdtab_i <- as.data.frame(do.call(rbind, sdtab))
  } else {
    stop("SDTAB must be a data frame or a list of data frames")
  }

  # Check for required columns and convert to numeric
  required_cols <- c("DVID", res_type)
  required_cols <- unique(required_cols)
  missing_cols <- setdiff(required_cols, colnames(sdtab_i))
  if (length(missing_cols) > 0) {
    stop("SDTAB is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Convert columns to numeric (only if they exist)
  sdtab_i <- sdtab_i %>%
    mutate(
      DVID = if ("DVID" %in% colnames(.)) as.numeric(DVID) else NA_real_
    )

  # Convert residual columns to numeric
  for (res_col in res_type) {
    if (res_col %in% colnames(sdtab_i)) {
      sdtab_i <- sdtab_i %>% mutate(across(all_of(res_col), as.numeric))
    }
  }

  # Add optional columns if they exist
  if ("WRES" %in% colnames(sdtab_i) && !("WRES" %in% res_type)) {
    sdtab_i <- sdtab_i %>% mutate(WRES = as.numeric(WRES))
  }
  if ("IWRES" %in% colnames(sdtab_i) && !("IWRES" %in% res_type)) {
    sdtab_i <- sdtab_i %>% mutate(IWRES = as.numeric(IWRES))
  }
  if ("CWRES" %in% colnames(sdtab_i) && !("CWRES" %in% res_type)) {
    sdtab_i <- sdtab_i %>% mutate(CWRES = as.numeric(CWRES))
  }
  if ("RES" %in% colnames(sdtab_i) && !("RES" %in% res_type)) {
    sdtab_i <- sdtab_i %>% mutate(RES = as.numeric(RES))
  }
  if ("IRES" %in% colnames(sdtab_i) && !("IRES" %in% res_type)) {
    sdtab_i <- sdtab_i %>% mutate(IRES = as.numeric(IRES))
  }

  res_for_plot <- sdtab_i %>% select(DVID,all_of(res_type))
  res_for_plot2 <- res_for_plot %>%
    gather(key = "residual_type", value = "value", -DVID)

  p_i <- ggplot(data = res_for_plot2, aes(x = value, y = ..density..)) +
    geom_histogram(bins = n_bins, col = "grey25", fill = MSDcol[2]) +
    facet_wrap(DVID~residual_type, scales = "free") +
    scale_y_continuous(name = "Density", breaks = scales::pretty_breaks(7), expand = c(0, 0), lim = c(0, NA)) +
    scale_x_continuous(name = 'Residuals', breaks = scales::pretty_breaks(7), expand = c(0, 0))+
    theme(panel.grid.minor = element_blank()) +
    theme_bw()
  if (ndist) {

    # mean_val <- mean(res_for_plot[[res_type]], na.rm = TRUE)
    # sd_val <- sd(res_for_plot[[res_type]], na.rm = TRUE)
    #
    #
    # x_min <- mean_val - 4 * sd_val
    # x_max <- mean_val + 4 * sd_val
    # x_seq <- seq(x_min, x_max, length.out = 1000)
    # p_i <- p_i + annotate("line", x = x_seq, y = dnorm(x_seq, mean = mean_val, sd = sd_val), linewidth = 0.8, lty = "dashed")


    norm_params <- res_for_plot2 %>%
      group_by(DVID, residual_type) %>%
      summarise(
        mean_val = mean(value, na.rm = TRUE),
        sd_val = sd(value, na.rm = TRUE),
        min_val = min(value, na.rm = TRUE),
        max_val = max(value, na.rm = TRUE),
        .groups = "drop"
      )

    norm_curves <- norm_params %>%
      group_by(DVID, residual_type) %>%
      summarise(
        x = list(seq(min_val - 0.5 * sd_val, max_val + 0.5 * sd_val, length.out = 1000)),
        y = list(dnorm(x[[1]], mean = mean_val, sd = sd_val)),
        .groups = "drop"
      ) %>%
      unnest(c(x, y))

    p_i <- p_i +  geom_line(data = norm_curves, aes(x = x, y = y),
                            color = "black", linewidth = 0.8, linetype = "dashed")
  }

  if (plot_type == 'QQ') {

    # Создаем базовый QQ-plot
    p_i <- ggplot(res_for_plot2, aes(sample = value)) +
      stat_qq(size = 1.75, color = MSDcol[1], alpha = 0.8) +
      stat_qq_line(col = "firebrick") +
      labs(x = "Theoretical quantiles", y = "Sample quantiles") +
      geom_abline(size = 0.5, col = "black", linetype = "dashed") +
      scale_x_continuous(breaks = scales::pretty_breaks(7)) +
      scale_y_continuous(breaks = scales::pretty_breaks(7)) +
      theme(panel.grid.minor = element_blank()) +
      theme_bw()+
      facet_wrap(DVID ~ residual_type, scales = "free")
  }
  return(p_i)



}

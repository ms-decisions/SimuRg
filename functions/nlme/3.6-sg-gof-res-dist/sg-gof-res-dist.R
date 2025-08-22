#' Plot distribution or QQ-plot of residuals from NONMEM/Simurg run
#'
#' This function visualizes residuals from estimation results
#' stored in a Simurg output object. The residuals can be plotted either as
#' histograms with overlaid normal density curves, or as QQ-plots to assess
#' normality.
#'
#' @param fpath_i Path to .RData file containing modeling results (must contain  `$SDTAB`).
#' @param res_type Character vector. One or several types of residuals to plot.  
#'   Values must correspond to column name(s) in `smrg_obj$SDTAB`, e.g.  
#'   `"RES"`, `"IWRES"`, `"IRES"`. If multiple residual types are provided,  
#'   plots will be generated for each of them.
#' @param n_bins Integer. Number of bins in histogram (default: 30).
#' @param ndist Logical. If `TRUE`, overlays the corresponding normal density
#'   curve on the histogram (default: TRUE).
#' @param plot_type Character. Type of plot to produce:  
#'   * `"DIST"` (default) — histogram of residuals with optional normal density,  
#'   * `"QQ"` — QQ-plot of residuals.
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
#' @export


sg_gof_res_dist <- function(fpath_i, res_type = 'RES', n_bins = 30, ndist = T, plot_type = 'DIST'){
  smrg_obj <- get(load(fpath_i))
  sdtab_i <- smrg_obj$SDTAB
  
  MSDcol <- c("#1a1866", "#f2b93b", "#b73b58", "#a2d620", "#14D98E", "#9c4ec7", "#3a6eba", "#efdd3c", "#69686d",'#844538', '#D91477','#F3A9FF')
  
  
  
  res_for_plot <- sdtab_i %>% select(DVID,all_of(res_type))
  res_for_plot2 <- res_for_plot %>% 
    gather(key = "residual_type", value = "value", -DVID)
  
  p_i <- ggplot(data = res_for_plot2, aes(x = value, y = ..density..)) +
    geom_histogram(bins = n_bins, col = "grey25", fill = MSDcol[2]) +
    facet_wrap(DVID~residual_type, scales = "free") +
    scale_y_continuous(name = "Density", breaks = scales::pretty_breaks(7), expand = c(0, 0), lim = c(0, NA)) +
    scale_x_continuous(name = all_of(res_type), breaks = scales::pretty_breaks(7), expand = c(0, 0))+
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
  
  # Создаем данные для нормальных кривых для каждой группы
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
    # Преобразуем в длинный формат
    ds_i_e2 <- ds_i_e %>% 
      gather(key = "residual_type", value = "value", -DVID)
    
    # Создаем базовый QQ-plot
    p_i <- ggplot(res_for_plot2, aes(sample = value)) +
      stat_qq(size = 1.75, color = MSDcol[1], alpha = 0.8) +
      stat_qq_line(col = "firebrick") +
      labs(x = "Theoretical quantiles", y = "Sample quantiles") +
      geom_abline(size = 0.5, col = "black", linetype = "dashed") +
      scale_x_continuous(breaks = scales::pretty_breaks(7)) +
      scale_y_continuous(breaks = scales::pretty_breaks(7)) +
      theme_bw()+ 
      facet_wrap(DVID ~ residual_type, scales = "free")
  }
 
  return(p_i)
  
  
 
}



#' Generate prediction distribution visualizations from sg_predist_sim () simulation results
#' @inheritParams sg_dummy
#' @param ds_sim A data frame containing simulation results, typically obtained from sg_predist_sim().
#' @param dt_obs_fl Logical. If TRUE, observed data points are overlaid on the simulation plots. Default is FALSE.
#' @param logy Logical. If TRUE, y-axis is displayed on a logarithmic scale. Default is FALSE.
#' @param legend_fl Logical. If FALSE, the legend is hidden. Default is TRUE.
#' @param pred_interval Character. Prediction interval to display. Options are "95%", "90%", "80%", "50%". Default is "80%".
#'
#' @returns A list of ggplot objects, one per output variable in the simulation dataset.
#' @import ggplot2
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @importFrom purrr map
#' @importFrom scales pretty_breaks trans_breaks trans_format math_format
#' @export
sg_predist_vis <- function(fpath_i,
                           ds_sim,
                           time_col = "TIME",
                           name_x = "TIME", name_y = 'DV',
                           dt_obs_fl = F, logy = F, legend_fl = T,
                           pred_interval = '80%'){
  
  MSDcol <- c("#1a1866", "#f2b93b", "#b73b58", "#a2d620", "#5839bb", "#9c4ec7", "#3a6eba", "#efdd3c", "#69686d")
  
  if (inherits(fpath_i, "character")) {
    if (file.exists(fpath_i)) {
      obj <- get(load(fpath_i))
    } else {
      stop("File specified by fpath_i does not exist")
    }
  } else if (inherits(fpath_i, "list")) {
    obj <- fpath_i
  } else {
    stop("fpath_i object should be either an sg_fit object, or a path to saved sg_fit object")
  }
  
  data_fin.noex <-  obj$SDTAB %>% filter(MDV != 1) %>% select(-MDV)
  data_fin.noex$TIME  <- data_fin.noex[[time_col]]
  
  ds_sim_l <- ds_sim %>%
    pivot_longer(
      cols = -c(ID, TIME),
      names_to = "VAR",
      values_to = "VALUE"
    )
  
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
  
  p_list <- unique(ds_sim_l$VAR) %>% map(function(v){
    
    ds_sim_v <- ds_sim_l %>% filter(VAR == v)
    
    data_obs_v <- data_fin.noex %>%
      select(ID, TIME, all_of(v)) %>%
      rename(DV = all_of(v))
    
    ds_sim_sum <- ds_sim_v %>%
      group_by(TIME, VAR) %>%
      summarise(across(VALUE, local_funSum), .groups = "drop") %>%
      mutate(
        label_median = "Median",
        label_band   = paste(pred_interval, "Prediction Interval")
      )
    
    p <- ggplot() +
      geom_ribbon(
        data = ds_sim_sum,
        aes(x = TIME, ymin = L_Q, ymax = H_Q, fill = label_band),
        alpha = 0.2
      ) +
      geom_line(
        data = ds_sim_sum,
        aes(x = TIME, y = median, color = label_median),
        size = 0.9
      ) +
      scale_x_continuous(name = name_x, breaks = pretty_breaks()) +
      scale_y_continuous(name = name_y, breaks = pretty_breaks()) +
      scale_color_manual(values = MSDcol[1]) +
      scale_fill_manual(values = MSDcol[1]) +
      theme_bw() +
      theme(
        axis.text  = element_text(size = 12),
        axis.title = element_text(size = 14)
      )
    
    if (dt_obs_fl) {
      p <- p +
        geom_point(
          data = data_obs_v,
          aes(x = TIME, y = DV),
          color = "royalblue4",
          size = 1.5
        )
    }
    
    if (logy) {
      p <- p +
        scale_y_log10(
          name = name_y,
          breaks = scales::trans_breaks("log10", function(x) 10^x),
          labels = scales::trans_format("log10", scales::math_format(10^.x))
        )
    }
    
    if (!legend_fl) {
      p <- p + theme(legend.position = "none")
    }
    
    p
  })
  
  return(p_list)
}

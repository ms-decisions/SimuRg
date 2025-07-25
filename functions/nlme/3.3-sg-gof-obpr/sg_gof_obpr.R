#'Observed vs. Predicted Plot
#'
#' @description
#' Creates a ggplot2 scatter plot of observed (DV) vs predicted (PRED/IPRED) values 
#' with options for faceting, coloring by covariates, and trend lines.
#'
#' @param fpath_i Path to .RData file containing modeling results (must contain `obj1` with `SDTAB`, `COTAB`, `CATAB`)
#' @param cov_cols Optional character vector of column names for covariates to include
#' @param indiv Logical. If `TRUE` (default), uses individual predictions (IPRED); otherwise uses population predictions (PRED)
#' @param addline Logical. Add connecting lines for individual subjects. Default `TRUE`
#' @param alpha_i Transparency level (0-1) for points/lines. Default 0.5
#' @param smooth Logical. Add LOESS smooth line. Default `TRUE`
#' @param log_axes Logical. Use log10 scaling for axes? Default `FALSE`
#' @param sc_factor Scaling factor for DV/PRED/IPRED values. Default 1 (no scaling)
#' @param abreaks Axis breaks function. Default `scales::pretty_breaks(7)`
#' @param xlab X-axis label. Default "Model-predicted values"
#' @param ylab Y-axis label. Default "Observed values"
#' @param col_i Optional column name for color
#' @param col_lab Label for color legend
#' @param facet_i Optional column name for facet
#' @param f_scales Facet scales. Default "fixed"
#' @param no_leg Logical. Remove legend? Default `FALSE`
#' @param n_quantiles Integer. Number of quantile groups for continuous variables in `col_i`. Default 3.
#' @param levels_discrete Integer. Maximum unique values to consider a variable discrete. Default 10.
#'
#' @return A ggplot2 object
#'
#' @details
#' The function automatically:
#' - Filters out MDV != 1 records
#' - Handles continuous covariates by converting to quantile categories
#' - Adds identity line and customizable trend lines
#'
#' @examples
#' \dontrun{
#' # Basic plot
#' p <- fun_ObPr.nm("model_results.RData")
#' 
#' # With covariates
#' p <- fun_ObPr.nm(
#'   "model_results.RData",
#'   cov_cols = "SEX",
#'   col_i = "SEX",
#'   facet_i = "SEX"
#' )
#' }
#'
#' @importFrom dplyr %>% filter mutate_at rename_at select distinct left_join
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_smooth labs 
#' @importFrom scales pretty_breaks number_format
#' @export
NULL
sg_gof_obpr <- function(
  fpath_i, cov_cols = NULL, indiv = T, addline = T, alpha_i = 0.5,
  smooth = T, log_axes = F, sc_factor = 1, abreaks = scales::pretty_breaks(7), 
  xlab = "Model-predicted values", ylab = "Observed values", col_i = NULL, col_lab = NULL, 
  facet_i = NULL, f_scales = "fixed",
  no_leg = F, n_quantiles = 3, levels_discrete = 10
){
  is_discrete <- function(x, max_levels = levels_discrete) {
    n_unique <- length(unique(na.omit(x)))
    n_unique <= max_levels
  }
  continuous_to_categories <- function(x, n_quant = 3) {
    cut(
      x,
      breaks = quantile(x, probs = seq(0, 1, length.out = n_quant + 1), na.rm = TRUE),
      include.lowest = TRUE,
      labels = paste0("Q", 1:n_quant)
    )
  }
  X <- ifelse(indiv, "IPRED", "PRED")
  load(fpath_i)
  ds_i <- obj1$SDTAB %>% 
    #read_table(fpath_i, skip = 1, col_names = T, col_types = cols()) %>%
    filter(MDV != 1) %>% 
    mutate_at(vars(IPRED, PRED, DV), function(s){s/sc_factor}) %>% 
    rename_at(vars(one_of(X)), function(n){n = "X"})
  
  if(!is.null(cov_cols)){
    suppressMessages({
      ds_covs <- left_join(obj1$COTAB,obj1$CATAB)})
    ds_covs_i <- ds_covs %>% select(ID, one_of(cov_cols)) %>% unique()
    stopifnot(n_distinct(ds_covs_i$ID) == nrow(ds_covs_i))
    suppressMessages({
      ds_i <- ds_i %>% left_join(ds_covs_i, by = "ID")})
    if (!is_discrete(ds_i[[cov_cols]])) {
      ds_i[[cov_cols]] <- continuous_to_categories(ds_i[[cov_cols]], n_quant = n_quantiles)
      #message(paste(cov_cols, "was converted into", n_quantiles, "quantile-based categories."))
    } else {
      ds_i[[cov_cols]] <- factor(ds_i[[cov_cols]])
    }
  }
  
  lim_obpr <- c(min(ds_i$X, ds_i$DV), max(ds_i$X, ds_i$DV))
  
  p_char <- list(
    labs(y = ylab, x = xlab),
    geom_abline(size = 0.5, col = "black", linetype = "dashed"),
    scale_color_manual(values = rep(MSDcol[c(2, 3, 4, 5, 7, 9, 6, 8, 1)], 30)),
    theme(legend.justification = c("left", "center"),
          legend.box.just = "left",
          legend.background = element_rect(fill = "white", size = 0.15, linetype = "solid", colour = "black"),
          legend.key.size = unit(0.38, "cm"),
          legend.title = element_text(size = 8), 
          legend.text = element_text(size = 8),
          plot.title = element_text(size = 12))
  )
  
  if(log_axes){
    p_char <- c(p_char, 
                scale_x_log10(
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))),
                scale_y_log10(
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))))
  } else {
    p_char <- c(p_char,
                scale_x_continuous(lim = lim_obpr, breaks = abreaks,
                                   labels = scales::number_format()),
                scale_y_continuous(lim = lim_obpr, breaks = abreaks,
                                   labels = scales::number_format()))
  }
  
  p_ObPr <<- ggplot(data = ds_i, aes(x = X, y = DV)) + p_char
  
  if(!is.null(col_i)){
    p_ObPr <- p_ObPr +
      geom_point(aes(col = !!sym(col_i)), size = 1.5, alpha = alpha_i) +
      guides(color = guide_legend(override.aes = list(alpha = 1))) +
      labs(col = col_lab)
    
    if(addline){ p_ObPr <- p_ObPr + geom_line(aes(col = !!sym(col_i), group = ID), lwd = 0.4, alpha = alpha_i) }
    if(smooth){ p_ObPr <- p_ObPr + geom_smooth(aes(col = !!sym(col_i), group = !!sym(col_i)), formula = "y ~ x", method = "loess", size = 1.2, se = F) }
    
  } else {
    p_ObPr <- p_ObPr +
      geom_point(size = 1.5, alpha = alpha_i) +
      guides(color = guide_legend(override.aes = list(alpha = 1))) +
      labs(col = col_lab)
    
    if(addline){ p_ObPr <- p_ObPr + geom_line(aes(group = ID), col = MSDcol[1], lwd = 0.4, alpha = alpha_i) }
    if(smooth){ p_ObPr <- p_ObPr + geom_smooth(formula = "y ~ x", method = "loess", size = 1.2, se = F, col = MSDcol[3]) }
  }
  
  if(!is.null(facet_i)){
    p_ObPr <- p_ObPr + facet_wrap(as.formula(str_c("~", facet_i)), scales = f_scales)
  }
  
  if(no_leg){
    p_ObPr <- p_ObPr + theme(legend.position = "none")
  }
  
  return(p_ObPr)
}

roxygen2::roxygenise(package.dir = '~/Simurg_function/sg_gof_obpr.R')

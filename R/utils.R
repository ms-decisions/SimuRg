#' The function to store common parameters description
#'
#' @param abreaks Axis breaks function. Default is `scales::pretty_breaks(7)`
#' @param addOFV logical. If `TRUE`, information about OFV and AIC of the model will be added. Default is `TRUE`
#' @param addline logical. If `TRUE`, lines connecting observations of individual subjects will be added. Default is `TRUE`
#' @param alpha_i numeric. Transparency level (from 0 to 1) for points/lines. Default is 0.5
#' @param cap string. Plot caption.
#' @param cov_cols vector of characters. Name of the columns with covariates
#' @param col_i string. Column name for color
#' @param col_lab string. Label for color legend
#' @param dens Logical. If `TRUE`, plot histogram/density of residuals instead of scatter
#' @param ds_covs data.frame. The dataframe with covariates
#' @param eta_seq vector of strings. Character vector of parameter names to be plotted (e.g., `c("ka", "Cl")`). If `NULL`, all parameters be included. Default is `NULL`
#' @param f_scales one of `"fixed"`, `"free"`, `"free_x"`, `"free_y"`. User can specify whether the scales (x and y axes) should be fixed across all panels (`"fixed"`), free for each panel (`"free"`), or free only in one dimension (`"free_x"` or `"free_y"`). Default is `"fixed"`
#' @param facet_i string. Column name for facet
#' @param fpath_i string or sg-fit object. If the string is given, the path to `.Rdata` file with sg-fit object is expected
#' @param indiv logical. If `TRUE` uses individual predictions (`"IPRED"`); otherwise uses population predictions (`"PRED"`). Default is `TRUE`
#' @param lab_x string. X-axis label
#' @param lab_y string. Y-axis label
#' @param levels_discrete integer. Maximum unique values to consider a variable discrete. Default is 10
#' @param log_x logical. If `TRUE`, a logarithmic scale is applied to x-axis. Default is `FALSE`
#' @param log_y logical. If `TRUE`, a logarithmic scale is applied to y-axis. Default is `FALSE`
#' @param log_axes logical. If `TRUE`, a logarithmic scale is applied to all axes. Default is `FALSE`
#' @param max_x numeric. X ax maximum limit. Default is `NULL`
#' @param min_x numeric. X ax minimum limit. Default is `NULL`
#' @param max_y numeric. Y ax maximum limit. Default is `NULL`
#' @param min_y numeric. Y ax minimum limit. Default is `NULL`
#' @param n_bins integer. Number of bins to use in the histogram. Default is 30.
#' @param n_quantiles integer. Number of quantile groups for continuous variables in `col_i`. Default is 3
#' @param no_leg logical. If `TRUE`, no legend will be displayed. Default is `FALSE`
#' @param plot_type Character. Type of plot to produce:
#'   * `"DIST"` (default) - histogram of individual parameters,
#'   * `"QQ"` - QQ-plot of individual parameters
#' @param run_id integer. Tested model ID. Default is 1.
#' @param sc_factor numeric. Scaling factor for DV/PRED/IPRED values. Default is 1 (no scaling)
#' @param smooth logical. Add LOESS smooth line. Default is `TRUE`
#' @param tdist logical. If `TRUE`, overlay theoretical parameter distributions based on population mean and OMEGA matrix. Default is `TRUE`
#' @param tsld logical. If `TRUE`, uses time since last dose instead of time from first dose. Default is `FALSE`
sg_dummy <- function() {}

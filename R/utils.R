#' The function to store common parameters description
#'
#' @param abreaks Axis breaks function. Default is `scales::pretty_breaks(7)`
#' @param addcov logical. If `TRUE`, columns with covariate values will be added to the resulting dataset. Default is `TRUE`
#' @param addOFV logical. If `TRUE`, information about OFV and AIC of the model will be added. Default is `TRUE`
#' @param addline logical. If `TRUE`, lines connecting observations of individual subjects will be added. Default is `TRUE`
#' @param aggr string. Aggregation type. Set to `NULL` for no aggregation, `ID` for aggregation by time, population and ID, `NPOP` for aggregation by population and time, `TIME` for aggregation by time. Default is `NULL`
#' @param alpha_i numeric. Transparency level (from 0 to 1) for points/lines. Default is 0.5
#' @param atol numeric. A numeric absolute tolerance used by the ODE solver to determine if a good solution has been achieved. This is also used in the solved linear model to check if prior doses do no add anything to the solution. Default is 1e-8
#' @param cap string. Plot caption.
#' @param cov_cols vector of characters. Name of the columns with covariates
#' @param covint string. Specifies the interpolation method for time-varying covariates. When solving ODEs it often samples times outside the sampling time specified in event table. When this happens, the time varying covariates are interpolated. Currently this can be:
#'  * `"linear"` interpolation, which interpolates the covariate by solving the line between the observed covariates and extrapolating the new covariate value
#'  * `"locf"` last observation carried forward
#'  * `"NOCB"` next observation carried backward. This is the same method that NONMEM uses
#'  * `"midpoint"` last observation carried forward to midpoint; next observation carried backward to midpoint
#'
#' Default is `"locf"`
#' @param col_i string. Column name for color
#' @param col_lab string. Label for color legend
#' @param data string. Path to the dataset used to fit a model
#' @param dens Logical. If `TRUE`, plot histogram/density of residuals instead of scatter
#' @param ds_covs data.frame. The dataframe with covariates
#' @param et data.frame. Event table
#' @param eta_seq vector of strings. Character vector of parameter names to be plotted (e.g., `c("ka", "Cl")`). If `NULL`, all parameters be included. Default is `NULL`
#' @param f_scales one of `"fixed"`, `"free"`, `"free_x"`, `"free_y"`. User can specify whether the scales (x and y axes) should be fixed across all panels (`"fixed"`), free for each panel (`"free"`), or free only in one dimension (`"free_x"` or `"free_y"`). Default is `"fixed"`
#' @param facet_i string. Column name for facet
#' @param fpath_i string or sg-fit object. If the string is given, the path to `.Rdata` or `.json` file with sg-fit object is expected
#' @param headers list. List with dataframe headers.
#' @param ncores integer. Number of cores used for calculations. Default is 1
#' @param indiv logical. If `TRUE` uses individual predictions (`"IPRED"`); otherwise uses population predictions (`"PRED"`). Default is `TRUE`
#' @param inits named vector. Initial conditions of model variables. Default is `NULL`
#' @param keep vector of strings. Columns of event table to keep in the output dataframe. Default is `NULL`
#' @param lab_x string. X-axis label
#' @param lab_y string. Y-axis label
#' @param levels_discrete integer. Maximum unique values to consider a variable discrete. Default is 10
#' @param log_x logical. If `TRUE`, a logarithmic scale is applied to x-axis. Default is `FALSE`
#' @param log_y logical. If `TRUE`, a logarithmic scale is applied to y-axis. Default is `FALSE`
#' @param log_axes logical. If `TRUE`, a logarithmic scale is applied to all axes. Default is `FALSE`
#' @param max_x numeric. X ax maximum limit. Default is `NULL`
#' @param max_y numeric. Y ax maximum limit. Default is `NULL`
#' @param maxsteps integer. Maximum number of steps allowed during one call to the solver. Default is 70000
#' @param min_x numeric. X ax minimum limit. Default is `NULL`
#' @param min_y numeric. Y ax minimum limit. Default is `NULL`
#' @param model RxODE model. The model to simulate from.
#' @param n_bins integer. Number of bins to use in the histogram. Default is 30.
#' @param n_quantiles integer. Number of quantile groups for continuous variables in `col_i`. Default is 3
#' @param no_leg logical. If `TRUE`, no legend will be displayed. Default is `FALSE`
#' @param npop integer. Number of population replicates. Default is 1
#' @param nsub integer. Number of subjects sampled per population (omega/sigma matrices per ID). Default is 1
#' @param occ interoccasion variability object. Object to set properties of interoccasion variability
#' @param omega named mztrix or vector. Matrix
#' @param outputs vector of strings. Names of the model variabeles to output. If `NULL`, all varaibles returned. Default is `NULL`
#' @param plot_type Character. Type of plot to produce:
#'   * `"DIST"` (default) - histogram of individual parameters,
#'   * `"QQ"` - QQ-plot of individual parameters
#' @param re random effects object. Contains options for random effects in model fit
#' @param rtol numeric. A numberic relative tolerance used by the ODE solver to determine if a good solution has been achieved.This is also used in the solved linear model to check if prior doses do not add anything to the solution. Default is 1e-6
#' @param run_id integer. Tested model ID. Default is 1.
#' @param ruv residual error object. Options for residual error used in model fit
#' @param sc_factor numeric. Scaling factor for DV/PRED/IPRED values. Default is 1 (no scaling)
#' @param smooth logical. Add LOESS smooth line. Default is `TRUE`
#' @param stimes vector of numeric. Sampling time points. Default is `NULL`
#' @param tdist logical. If `TRUE`, overlay theoretical parameter distributions based on population mean and OMEGA matrix. Default is `TRUE`
#' @param theta named vector or data.frame. Values of population parameters to simulate with. Default is `NULL`
#' @param thetamat matrix. Named theta matrix. Default is `NULL`
#' @param tsld logical. If `TRUE`, uses time since last dose instead of time from first dose. Default is `FALSE`
sg_dummy <- function() {}

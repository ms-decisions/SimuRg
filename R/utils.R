utils::globalVariables(c(":=", ".", "..density..", ".x", "95% CI", "ANOVA", "ATSLD",
                         "BIN", "BIN_sort", "COHORTC", "COL", "CORR_PAR", "COV",
                         "COVNAME", "COVVAL", "CV", "CV % (95% CI)", "CWRES",
                         "DT_label", "DV", "DVID", "EST", "Estimate", "ETA",
                         "ETAshrinkage_sd", "ETAshrinkage_var", "EVID", "ID",
                         "INIT", "IPRED", "IRES",
                         "IWRES", "LABEL", "LARGPVAL", "LCI",  "MDV", "METRIC",
                         "NAME", "NAME_INIT", "NTILE", "OCCMAN", "PAR", "PAR1",
                         "PAR2", "PAR_NAME", "Parameter", "PARVAL", "PAR_group",
                         "PARVAL_NORM", "PARNAME", "PEARS", "PI", "PNAME", "POPN",
                         "PPRED", "PPRED_bin", "PRED", "PVAL", "Project_name",
                         "R_diff", "RES", "RSE", "RUVpars", "Run_ID", "SE",
                         "SIGNPEARS", "Shrinkage", "Shrinkage_var", "Shrinkage (var), %",
                         "Source", "TDIST", "TIME", "TIME_BIN", "TIME_BIN_max",
                         "TIME_BIN_max_prev", "TIME_BIN_min", "TIME_BIN_min_next",
                         "TRANS", "TSLD", "TV", "TYPE", "UCI", "VALUE", "VAR",
                         "Var1", "Var2", "Var_epsilon", "Var_mc", "Var_total",
                         "WRES", "X", "Y", "aov", "as.formula", "category", "cmax",
                         "cmax_ref", "combn", "cor", "cor.test", "data", "density",
                         "dnorm", "dvid", "errorModel", "inside", "geom_mean_val",
                         "kmeans", "ks.test", "key", "lperc", "max_val", "mean_val",
                         "median", "median_ci_l", "median_ci_u", "median_median",
                         "median_val", "min_BIN_time", "min_val", "na.omit", "p",
                         "parameter", "pi_l", "pi_l_ci_l", "pi_l_ci_u", "pi_l_median",
                         "pi_u", "pi_u_ci_l", "pi_u_ci_u", "pi_u_median", "popPred",
                         "q25", "q75", "quantile", "r", "read.csv", "reorder",
                         "residual_type", "rnorm", "row_id", "sd", "sd_val",
                         "se_val", "setNames", "theor_median", "time", "toJSON",
                         "type", "use", "uperc","value", "value_at_time",
                         "value_cfb", "value_ref", "var", "x", "y", "y_central",
                         "y_lower", "y_lower_perc", "y_upper", "y_upper_perc"))
#' The function to store common parameters description
#'
#' @param abreaks axis breaks function. Default is `scales::pretty_breaks(7)`
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
#' @param ... other arguments that will be passed to rxSolve function.
#' @param ciLow numeric. Lower confidence interval bound. Default is 0.025
#' @param ciUp numeric. Upper confidence interval bound. Default is 0.975
#' @param col_i string. Column name for color
#' @param col_lab string. Label for color legend
#' @param covs list of covariate structures. List object for covariates specification
#' @param data string. Path to the dataset used to fit a model
#' @param dens logical. If `TRUE`, plot histogram/density of residuals instead of scatter
#' @param ds_covs data.frame. The dataframe with covariates
#' @param ds_i data.frame. The data frame with source data.
#' @param dt_obs_fl logical. Show observed data points. Default is `FALSE`
#' @param dv_col character. Name of DV column in data_i. Default is`DV`
#' @param emp_perc logical. Show empirical percentiles. Default is `TRUE`
#' @param et data.frame. Event table
#' @param eta_seq vector of strings. Character vector of parameter names to be plotted (e.g., `c("ka", "Cl")`). If `NULL`, all parameters be included. Default is `NULL`
#' @param excl_col character vector. Contains column names to exclude from synthesis. Default: \code{NULL}
#' @param f_scales one of `"fixed"`, `"free"`, `"free_x"`, `"free_y"`. User can specify whether the scales (x and y axes) should be fixed across all panels (`"fixed"`), free for each panel (`"free"`), or free only in one dimension (`"free_x"` or `"free_y"`). Default is `"fixed"`
#' @param facet_i string. Column name for facet
#' @param fill_i string. Column name for fill aesthetic. Default is `NULL`
#' @param filt string. Provide a filter to apply. Default is `"T"`
#' @param fit logical. If`TRUE`, fit will be performed. Default is `FALSE`
#' @param fpath_i string or sg-fit object. If the string is given, the path to `.Rdata` or `.json` file with sg-fit object is expected
#' @param free_stat string. Facet scaling option. One of `"free"`, `"free_x"`, `"free_y"`, or `"fixed"`. Default is `'free'`
#' @param group_i string. Primary grouping variable for lines. Default is `'VAR'`
#' @param headers list. List with dataframe headers.
#' @param ncores integer. Number of cores used for calculations. Default is 1
#' @param id_col character string. Specify the name of the identifier column to exclude from synthesis. Default: \code{NULL}
#' @param indiv logical. If `TRUE` uses individual predictions (`"IPRED"`); otherwise uses population predictions (`"PRED"`). Default is `TRUE`
#' @param inits named vector. Initial conditions of model variables. Default is `NULL`
#' @param keep vector of strings. Columns of event table to keep in the output dataframe. Default is `NULL`
#' @param lab_x string. X-axis label
#' @param lab_y string. Y-axis label
#' @param legend_fl logical. Show legend. Default is `FALSE`
#' @param levels_discrete integer. Maximum unique values to consider a variable discrete. Default is 10
#' @param log_x logical. If `TRUE`, a logarithmic scale is applied to x-axis. Default is `FALSE`
#' @param log_y logical. If `TRUE`, a logarithmic scale is applied to y-axis. Default is `FALSE`
#' @param log_axes logical. If `TRUE`, a logarithmic scale is applied to all axes. Default is `FALSE`
#' @param lty_i string. Column name for linetype aesthetic. Default is `NULL`
#' @param max_x numeric. X ax maximum limit. Default is `NA`
#' @param max_y numeric. Y ax maximum limit. Default is `NA`
#' @param method one of c("liblsoda", "lsoda", "dop853", "indLin"). Method for solving ODE.
#' @param min_x numeric. X ax minimum limit. Default is `NA`
#' @param min_y numeric. Y ax minimum limit. Default is `NA`
#' @param model RxODE model. The model to simulate from.
#' @param n_bins integer. Number of bins to use in the histogram. Default is 30.
#' @param n_quantiles integer. Number of quantile groups for continuous variables in `col_i`. Default is 3
#' @param no_leg logical. If `TRUE`, no legend will be displayed. Default is `FALSE`
#' @param npop integer. Number of population replicates. Default is 1
#' @param nsub integer. Number of subjects sampled per population (omega/sigma matrices per ID). Default is 1
#' @param occ interoccasion variability object. Object to set properties of interoccasion variability
#' @param omega named mztrix or vector. Matrix
#' @param opt_name string. Specify the optimizer to use. Currently, only `"Monolix"` is avaliable. Default is `"Monolix"`
#' @param outputs vector of strings. Names of the model variabeles to output. If `NULL`, all varaibles returned. Default is `NULL`
#' @param path_to_save_output string. Path to save fit output files. If `NULL`, current working directroy will be used. Default is `NULL`
#' @param path_to_fitter string. The path to the program fitter. If `NULL`, "C:/ProgramData/Lixoft/MonolixSuite2023R1/bin/monolix.bat" will be used. Defult is `NULL`
#' @param piLow numeric. Lower prediction interval bound. Default is 0.10
#' @param piUp numeric. Upper prediction interval bound Default is 0.90
#' @param plot_type Character. Type of plot to produce:
#'   * `"DIST"` (default) - histogram of individual parameters,
#'   * `"QQ"` - QQ-plot of individual parameters
#' @param pred.corr logical. Apply prediction correction. Default is `FALSE`
#' @param project_name string. The name of the project.
#' @param re random effects object. Contains options for random effects in model fit
#' @param rtol numeric. A numberic relative tolerance used by the ODE solver to determine if a good solution has been achieved. This is also used in the solved linear model to check if prior doses do not add anything to the solution. Default is 1e-6
#' @param run_id integer. Tested model ID. Default is 1.
#' @param ruv residual error object. Options for residual error used in model fit
#' @param sc_factor numeric. Scaling factor for DV/PRED/IPRED values. Default is 1 (no scaling)
#' @param scale numeric named vector. Scaling for ode parameters of the system. The names must correspond to the parameter identifiers in the ODE specification. Each of the ODE variables will be divided by the scaling factor. For example scale=c(center=2) will divide the center ODE variable by 2. Defult is `NULL`
#' @param seed integer. Random seed for synthetic data generation reproducibility.Default is \code{123}.
#' @param shp_i string. Column name for shape aesthetic. Default is `NULL`
#' @param sigma matrix. Named sigma covariance or Cholesky decomposition of a covariance matrix. Defult is `NULL`
#' @param smooth logical. Add LOESS smooth line. Default is `TRUE`
#' @param stimes vector of numeric. Sampling time points. Default is `NULL`
#' @param task_opt string. Options to be passed to Monolix fit. Default is `NULL`
#' @param tdist logical. If `TRUE`, overlay theoretical parameter distributions based on population mean and OMEGA matrix. Default is `TRUE`
#' @param time_col string. The column to use as a time column. Currently, can be only `TIME`. Default is `TIME`
#' @param theor_perc logical. Show theoretical percentiles. Default is `TRUE`
#' @param theor_percCI logical. Show CI around theoretical percentiles. Default is `TRUE`
#' @param theta named vector or data.frame. Values of population parameters to simulate with. Default is `NULL`
#' @param thetamat matrix. Named theta matrix. Default is `NULL`
#' @param tsld logical. If `TRUE`, uses time since last dose instead of time from first dose. Default is `FALSE`
#' @param val_col string. Name of value column. Default is `VALUE`
#' @param wrap_i string. Faceting formula for `facet_wrap` (e.g., `'~VAR'`). Default is `NULL`
#' @param wrap_ncol integer. Number of columns for `facet_wrap`. Default is `NULL`
#' @param wrap_nrow integer. Number of rows for `facet_wrap`. Default is `NULL`
sg_dummy <- function(
    ...,
  abreaks = scales::pretty_breaks(7),
  addcov = TRUE,
  addOFV = TRUE,
  addline = TRUE,
  aggr = NULL,
  alpha_i = 0.5,
  atol = 1e-8,
  cap,
  cov_cols,
  covint = "locf",
  ciLow = 0.025,
  ciUp = 0.975,
  col_i,
  col_lab,
  covs,
  data,
  dens,
  ds_covs,
  ds_i,
  dt_obs_fl = FALSE,
  dv_col = "DV",
  emp_perc = TRUE,
  et,
  eta_seq = NULL,
  excl_col = NULL,
  f_scales = "fixed",
  facet_i,
  fill_i = NULL,
  filt = "T",
  fit = FALSE,
  fpath_i,
  free_stat = "free",
  group_i = "VAR",
  headers,
  ncores = 1,
  id_col = NULL,
  indiv = TRUE,
  inits = NULL,
  keep = NULL,
  lab_x,
  lab_y,
  legend_fl = FALSE,
  levels_discrete = 10,
  log_x = FALSE,
  log_y = FALSE,
  log_axes = FALSE,
  lty_i = NULL,
  max_x = NULL,
  max_y = NULL,
  method = "lsoda",
  min_x = NULL,
  min_y = NULL,
  model,
  n_bins = 30,
  n_quantiles = 3,
  no_leg = FALSE,
  npop = 1,
  nsub = 1,
  occ,
  opt_name = "Monolix",
  omega,
  outputs = NULL,
  path_to_fitter = NULL,
  path_to_save_output = NULL,
  piLow = 0.10,
  piUp = 0.90,
  plot_type = "DIST",
  pred.corr = FALSE,
  project_name,
  re,
  rtol = 1e-6,
  run_id = 1,
  ruv,
  sc_factor = 1,
  scale = NULL,
  seed = 123,
  sigma = NULL,
  shp_i = NULL,
  smooth = TRUE,
  stimes = NULL,
  task_opt = NULL,
  tdist = TRUE,
  time_col = "TIME",
  theor_perc = TRUE,
  theor_percCI = TRUE,
  theta = NULL,
  thetamat = NULL,
  tsld = FALSE,
  val_col = "VALUE",
  wrap_i = NULL,
  wrap_ncol = NULL,
  wrap_nrow = NULL
) {}

# Ensure SimuRg object table components are single data frames
#
# Converts SDTAB, EVTAB, COTAB, CATAB, SUMTAB from list-of-rows to one data frame
# when needed. Idempotent if already data frames.

smrg_ensure_tables_df <- function(obj) {
  table_names <- c("SDTAB", "EVTAB", "COTAB", "CATAB", "SUMTAB")
  for (nm in table_names) {
    x <- obj[[nm]]
    if (is.null(x)) next
    if (is.data.frame(x)) next
    if (!is.list(x) || length(x) == 0) next
    # List of rows (each element is a list or atomic vector)
    if (is.list(x[[1]]) && !is.data.frame(x[[1]])) {
      obj[[nm]] <- as.data.frame(do.call(rbind, lapply(x, as.data.frame)))
    } else {
      obj[[nm]] <- as.data.frame(x)
    }
  }
  obj
}

read_smrg_obj <- function(fpath_i) {
  if (inherits(fpath_i, "character")) {
    if (!file.exists(fpath_i)) {
      stop("File does not exist: ", fpath_i)
    }

    ext <- tools::file_ext(fpath_i)

    if (tolower(ext) == "rdata") {
      result <- get(load(fpath_i))
    } else if (tolower(ext) == "json") {
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Package 'jsonlite' is required for reading JSON files.")
      }
      # simplifyVector = TRUE, simplifyDataFrame = TRUE so SDTAB and other
      # array-of-objects become single data frames
      result <- jsonlite::fromJSON(fpath_i, simplifyVector = TRUE, simplifyDataFrame = TRUE)
      # Ensure known table components are data frames (in case of list of rows)
      result <- smrg_ensure_tables_df(result)
    } else {
      stop("Unsupported file type: ", ext, ". Supported: .RData, .json")
    }

    return(result)
  } else if (inherits(fpath_i, "list")) {
    # fpath_i is already an object (sg_fit object); ensure tables are data frames
    return(smrg_ensure_tables_df(fpath_i))
  } else if (inherits(fpath_i, "data.frame")) {
    stop("fpath_i cannot be a data.frame. Provide either a file path (character) or an sg_fit object (list with SDTAB)")
  } else {
    stop("fpath_i should be either a file path (character) or an sg_fit object (list). Got: ", class(fpath_i)[1])
  }
}

utils::globalVariables(c(":=", ".", "..density..", ".x", "95% CI", "90%CI", "ANOVA", "ATSLD", "BCOVVAL",
                         "BIN", "BIN_sort", "BTR", "CATDES", "COHORTC", "COL", "CORR_PAR", "COV",
                         "COVNAME", "COVVAL", "CV", "CV % (95% CI)", "CWRES",
                         "DT_label", "DV", "DVID", "EST", "Estimate", "ETA",
                         "ETAshrinkage_sd", "ETAshrinkage_var", "EVID", "ID",
                         "First order", "INIT", "IPRED", "IRES",
                         "IWRES", "KEY", "LAB", "LABEL", "LARGPVAL", "LCI", "LP", "MDV", "METRIC",
                         "NAME", "NAME_INIT", "NICEN", "NTILE", "NVALUE", "OCCMAN", "OUT","P05", "P95", "P975", "PAR", "PAR1",
                         "PAR2", "PAR_NAME", "PAR_f", "Parameter", "PARVAL", "PAR_group", "PCH",
                         "PARVAL_NORM", "PARNAME", "PEARS", "PI", "PNAME", "POPN",
                         "PPRED", "PPRED_bin", "PRED", "PVAL", "Project_name",
                         "R_diff", "REF_spec", "REFVAL", "REF", "RES", "Regimen", "RSE", "RUVpars", "Run_ID", "SE",
                         "SIGNPEARS", "Shrinkage", "Shrinkage_var", "Shrinkage (var), %",
                         "Source", "STAT", "TDIST", "TIME", "TIME_BIN", "TIME_BIN_max",
                         "TIME_BIN_max_prev", "TIME_BIN_min", "TIME_BIN_min_next", "TR",
                         "TRANS", "TSLD", "Total order", "TV", "TVALUE", "TYPE", "Type", "UCI", "UP", "VALUE", "VAR",
                         "Var1", "Var2", "Var_epsilon", "Var_mc", "Var_total",
                         "WRES", "X", "Y", "DVNAME", "aov", "as.formula", "category", "cmax",
                         "cmax_ref", "combn", "cor", "cor.test", "data", "density",
                         "dnorm", "dvid", "errorModel", "geom_mean_val", "inside",
                         "kmeans", "ks.test", "key", "lperc", "max_val", "mean_val",
                         "median", "median_ci_l", "median_ci_u", "median_median",
                         "median_val", "min_BIN_time", "min_val", "na.omit", "p",
                         "parameter", "pi_l", "pi_l_ci_l", "pi_l_ci_u", "pi_l_median",
                         "pi_u", "pi_u_ci_l", "pi_u_ci_u", "pi_u_median", "popPred",
                         "q25", "q75", "quantile", "quantiles", "r", "read.csv", "reorder",
                         "residual_type", "rnorm", "row_id", "sd", "sd_val",
                         "se_val", "setNames", "sim.id", "theor_median", "time", "toJSON",
                         "type", "use", "uperc","value", "value_at_time",
                         "value_cfb", "value_ref", "var", "x", "y", "y_central",
                         "y_lower", "y_lower_perc", "y_upper", "y_upper_perc"))
#' The function to store common parameters description
#'
#' @param abreaks A function that generates axis breaks. Default is `scales::pretty_breaks(7)`.
#' @param addcov Logical. If `TRUE`, columns with covariate values will be added
#'  to the resulting dataset. Default is `TRUE`.
#' @param addOFV Logical. If `TRUE`, information about OFV and AIC of the model
#'  will be added. Default is `TRUE`.
#' @param addline Logical. If `TRUE`, lines connecting observations of
#'  individual subjects will be added. Default is `TRUE`.
#' @param aggr A string that specifies the aggregation type. Set to `NULL` for
#'  no aggregation, `ID` for aggregation by time, population and ID, `NPOP` for
#'  aggregation by population and time, `TIME` for aggregation by time.
#'  Default is `NULL`.
#' @param alpha_i Numeric. Transparency level (from 0 to 1) for points/lines.
#'  Default is 0.5.
#' @param atol A numeric absolute tolerance used by the ODE solver to determine
#'  if a good solution has been achieved. This is also used in the solved linear
#'  model to check if prior doses do no add anything to the solution.
#'  Default is 1e-8.
#' @param cap A string that specifies the plot caption.
#' @param cat_cov_l A named list. Each element defines one categorical covariate
#'   and must itself be a list with components:
#'   \describe{
#'     \item{\code{NAME}}{Character. Column name of the covariate
#'       (e.g. \code{"SEX"}).}
#'     \item{\code{NICENAME}}{Character or \code{NULL}. Display label.}
#'     \item{\code{REF}}{Character or \code{NULL}. Reference category value.
#'       If \code{NULL}, the first factor level (alphabetically) is used.}
#'     \item{\code{par_vec}}{Character vector. Model parameter(s) affected by
#'       this covariate (e.g. \code{c("ka")}).}
#'   }
#' @param cont_cov_l A named list. Each element defines one continuous covariate
#'   and must itself be a list with components:
#'   \describe{
#'     \item{\code{NAME}}{Character. Column name of the (transformed)
#'       covariate in the dataset (e.g. \code{"LG_AGE"}).}
#'     \item{\code{UTNAME}}{Character or \code{NULL}. Column name of the
#'       untransformed (back-transformed) covariate
#'       (e.g. \code{"AGE"}).  If \code{NULL} or \code{NA}, defaults to
#'       \code{NAME}.}
#'     \item{\code{REF}}{Character or numeric. Reference value for the
#'       covariate.  Use \code{"median"} to derive from data, or a numeric
#'       value.}
#'     \item{\code{NICENAME}}{Character or \code{NULL}. Display label for
#'       plots and tables (e.g. \code{"Age, years"}).}
#'     \item{\code{par_vec}}{Character vector. Model parameter(s) affected by
#'       this covariate (e.g. \code{c("CL")}).}
#'   }
#' @param cov_cols A character vector specifying the names of the columns with covariates.
#' @param covint String. Specifies the interpolation method for time-varying
#'  covariates. When solving ODEs it often samples times outside the sampling
#'  time specified in event table. When this happens, the time varying
#'  covariates are interpolated. Currently this can be:
#'  * `"linear"` interpolation, which interpolates the covariate by solving the
#'  line between the observed covariates and extrapolating the new covariate value
#'  * `"locf"` last observation carried forward
#'  * `"NOCB"` next observation carried backward.
#'  This is the same method that NONMEM uses
#'  * `"midpoint"` last observation carried forward to midpoint;
#'  next observation carried backward to midpoint
#'
#' Default is `"locf"`
#' @param ... Other arguments that will be passed to rxSolve function.
#' @param ciLow Numeric. Lower confidence interval bound. Default is 0.025
#' @param ciUp Numeric. Upper confidence interval bound. Default is 0.975
#' @param ci_band_alpha Numeric in \eqn{[0,1]}: transparency of the shaded band.
#'   Default \code{0.2}.
#' @param ci_band_col Color for the band fill and dotted limit lines.
#'   Default \code{"firebrick"}.
#' @param ci_limits Numeric vector of length 2: lower and upper bounds of the
#'   shaded acceptance band and of the dotted horizontal guides.  Default
#'   \code{c(0.8, 1.25)} is a common bioequivalence-style window on the ratio
#'   scale.
#' @param ci_quantiles Character vector of length 2: names of the lower and
#'   upper uncertainty columns in the sensitivity table, in that order.
#'   Defaults \code{c("P025", "P975")} to match the default percentiles in
#'   \code{sg_covsens_sim}.  Use other names (e.g. \code{c("P05", "P95")}) if
#'   you changed \code{quantiles} in the simulation and the columns exist.
#' @param col_i String. Column name for color
#' @param col_lab String. Label for color legend
#' @param covs A list of covariate structures. Each element should be a list containing covariate relationship definitions with the following structure:
#'   * `PAR` - string, name of the parameter to which covariate effect is applied
#'   * `COVNAME` - string, name of the covariate column in the dataset
#'   * `FUNC` - string, functional form of the covariate effect
#'   * `TRANS` - string, transformation applied to covariate
#'   * `REF` - numeric, reference value for categorical covariates
#'   * `INIT` - numeric, initial value for the covariate effect parameter
#'   * `EST` - logical, whether to estimate the covariate effect
#'   Can be `NULL`, if no covariates are used. Default is `NULL`
#' @param covsens_res Named list as returned by \code{sg_covsens_sim()}.  Must
#' contain the element selected by \code{type} (\code{PARSENS} and/or
#' \code{EXPSENS} data.frames with columns \code{LAB}, \code{VAR},
#' \code{mean}, \code{Type}, and the interval columns named by
#' \code{ci_quantiles}).
#' @param data String. Path to the dataset used to fit a model. Should be a CSV file containing the pharmacokinetic/pharmacodynamic data with appropriate column structure matching the headers specification
#' @param dens Logical. If `TRUE`, plot histogram/density of residuals instead of scatter
#' @param ds_covs Data.frame. The dataframe with covariates
#' @param ds_i Data.frame. The data frame with source data.
#' @param ds_parest Data.frame. Parameter estimates table with columns
#'   \code{parameter} and \code{value}.
#'   Required when \code{fpath_i} is \code{NULL}; must be provided together
#'   with \code{ds_covs}.  Default is \code{NULL}.
#' @param dt_obs_fl Logical. Show observed data points. Default is `FALSE`
#' @param DVID Restrict \code{SDTAB} to one observation type. Numeric values
#'   select the \code{DVID} column (default \code{1}). If \code{SDTAB} has
#'   \code{DVNAME}, a character or factor is matched to \code{DVNAME} first;
#'   otherwise a digit-only string is coerced to numeric \code{DVID}.
#' @param dv_col Character. Name of DV column in data_i. Default is `DV`
#' @param emp_perc Logical. Show empirical percentiles. Default is `TRUE`
#' @param errorbar_width Numeric. Width argument for \code{geom_errorbar}.  Default
#'   \code{0.2}.
#' @param est_covmat Data.frame. Parameter estimation covariance matrix.  The
#'   first column (\code{X1}) must list parameter names; remaining columns
#'   (named identically) form the symmetric variance–covariance matrix.
#' @param et Data.frame. Event table
#' @param eta_seq Vector of strings. Character vector of parameter names to be plotted. If `NULL`, all parameters be included. Default is `NULL`
#' @param par_seq Vector of strings. Character vector of parameter names to be plotted. If `NULL`, all parameters be included. Default is `NULL`
#' @param par_type String. A character string specifying the type of parameters used for
#'   theoretical distribution overlay (only relevant when `plot_type = 'DIST'`
#'   and `tdist = TRUE`). If 'Ind' - Individual (default), distributions are shown on the
#'   natural parameter scale assuming log-normal variability. If 'RE' - Random Effect, distributions are shown on the ETA scale,
#'   assuming a normal distribution with mean zero and covariance defined
#'   by `$OMEGAMAT`, without transformation.
#' @param excl_col Character vector. Contains column names to exclude from synthesis. Default: \code{NULL}
#' @param exclude_vars Character vector. Contains \code{VAR} levels to omit (e.g.
#'   \code{"Cc_Cmin"}).  \code{NULL} keeps all rows.
#' @param f_scales String, one of `"fixed"`, `"free"`, `"free_x"`, `"free_y"`. User can specify whether the scales (x and y axes) should be fixed across all panels (`"fixed"`), free for each panel (`"free"`), or free only in one dimension (`"free_x"` or `"free_y"`). Default is `"fixed"`
#' @param facet_i String. Column name for facet
#' @param fill_i String. Column name for fill aesthetic. Default is `NULL`
#' @param filt String. Provide a filter to apply. Default is `"T"`
#' @param fit Logical. If `TRUE`, the model fitting will be executed immediately
#'  using the specified fitter. If `FALSE`, only the fit configuration file will
#'  be generated without running the fit. Set to `FALSE` for file preparation only,
#'  or `TRUE` to run the complete fitting process. Default is `FALSE`.
#' @param fpath_i String or sg-fit object. If the string is given, the path to
#'  `.Rdata` or `.json` file with sg-fit object is expected.
#' @param free_stat String. Facet scaling option. One of `"free"`, `"free_x"`,
#'  `"free_y"`, or `"fixed"`. Default is `'free'`.
#' @param group_i String. Primary grouping variable for lines. Default is `'VAR'`.
#' @param headers List. A list specifying the data frame column headers.
#' Each element should be a list containing column information with the following structure:
#'   * `name` - string, column name in the dataset
#'   * `use` - string, column usage type. Valid values include:
#'     - "identifier" for subject ID columns
#'     - "time" for time columns
#'     - "observation" for dependent variable columns
#'     - "observationtype" for observation type identifier
#'     - "administration" for administration route
#'     - "amount" for dose amount
#'     - "eventidentifier" for event ID
#'     - "missingdependentvariable" for missing DV flag
#'     - "covariate" for covariate columns
#'   * `type` - string or NULL, data type specification. For observations use
#'    "continuous", "count/categorical" or "event", depending on the nature of
#'    observations. For covariates use "continuous" or "categorical". Can be NULL for non-covariate columns.
#' @param ncores Integer. Number of cores used for calculations. Default is 1.
#' @param id_col Character string. Specify the name of the identifier column
#'  to exclude from synthesis. Default: \code{NULL}.
#' @param indiv Logical. If `TRUE`, use individual predictions (`"IPRED"`);
#'  otherwise use population predictions (`"PRED"`). Default is `TRUE`.
#' @param inits A named vector specifying initial conditions for model variables; Default is `NULL`.
#' @param keep Vector of strings. Columns of event table to keep in the output data frame. Default is `NULL`.
#' @param lab_x String. X-axis label.
#' @param lab_y String. Y-axis label.
#' @param legend_fl Logical. Show legend. Default is `FALSE`.
#' @param levels_discrete Integer. Maximum unique values to consider a variable discrete. Default is 10.
#' @param log_x Logical. If `TRUE`, a logarithmic scale is applied to x-axis. Default is `FALSE`.
#' @param log_y Logical. If `TRUE`, a logarithmic scale is applied to y-axis. Default is `FALSE`.
#' @param log_axes Logical. If `TRUE`, a logarithmic scale is applied to all axes. Default is `FALSE`.
#' @param lty_i String. Column name for linetype aesthetic. Default is `NULL`.
#' @param max_x Numeric. X-axis maximum limit. Default is `NA`.
#' @param max_y Numeric. Y-axis maximum limit. Default is `NA`.
#' @param method One of c("liblsoda", "lsoda", "dop853", "indLin"). Method for solving ODE.
#' @param min_x Numeric. X-axis minimum limit. Default is `NA`.
#' @param min_y Numeric. Y-axis minimum limit. Default is `NA`.
#' @param model RxODE model. The model to simulate from.
#' @param n_bins Integer. Number of bins to use in the histogram. Default is 30.
#' @param n_quantiles Integer. Number of quantile groups for continuous variables in `col_i`. Default is 3.
#' @param npop Integer. Number of population replicates. Default is 1.
#' @param nsub Integer. Number of subjects sampled per population (omega/sigma matrices per ID). Default is 1.
#' @param occ A list specifying interoccasion variability properties, containing:
#'   * `init` - matrix, initial values for the interoccasion variance-covariance matrix.
#'    Same structure as `re$init` but for occasion-to-occasion variability. Use 0 for no interoccasion variability
#'   * `est` - matrix, logical matrix specifying which interoccasion variance-covariance
#'    elements to estimate. Use `TRUE` to estimate, `FALSE` to fix, `NA` for
#'    elements not applicable. Typically all elements are NA when no interoccasion variability is modeled.
#' @param omega A named matrix or vector.
#' @param opt_name String. Specify the optimizer/fitter to use for model fitting.
#'  Currently supported options:
#'   * `"Monolix"` - uses Monolix Suite for population pharmacokinetic modeling (generates .mlxtran files)
#'   * `"Simurg"` - uses SimuRg internal fitter (generates .R files with JSON control structure)
#'   Default is `"Monolix"`.
#' @param outputs Vector of strings. Names of the model variables to output.
#'  If `NULL`, all variables returned. Default is `NULL`.
#' @param path_to_save_output String. Path to save fit output files.
#'  Should be a valid directory path where the fit results and project files
#'  will be saved. If `NULL`, current working directory will be used. Default is `NULL`.
#' @param path_to_fitter String. The path to the program fitter executable.
#'  For Monolix, this should point to the monolix.bat file.
#'  If `NULL`, "C:/ProgramData/Lixoft/MonolixSuite2023R1/bin/monolix.bat" will be used as default. Default is `NULL`.
#' @param piLow Numeric. Lower prediction interval bound. Default is 0.10.
#' @param piUp Numeric. Upper prediction interval bound Default is 0.90.
#' @param plot_type Character. Type of plot to produce:
#'   * `"DIST"` (default) - histogram of individual parameters,
#'   * `"QQ"` - QQ-plot of individual parameters.
#' @param point_size Numeric. Point size for \code{geom_point}.  Default \code{2.5}.
#' @param pred.corr Logical. Apply prediction correction. Default is `FALSE`.
#' @param project_name String. The name of the Monolix project without file extension.
#'  This will be used as the base name for output files and directories.
#' @param re List. A list specifying options for random effects used in model fit, containing:
#'   * `init` - matrix, initial values for the variance-covariance matrix of random
#'    effects. Rows and columns should correspond to parameters defined in theta.
#'    Diagonal elements represent variances, off-diagonal elements represent covariances. Use 0 for no variability
#'   * `est` - matrix, logical matrix of same dimensions as `init` specifying
#'    which variance-covariance elements to estimate. Use `TRUE` to estimate,
#'    `FALSE` to fix, `NA` to not use this random effect.
#' @param ref_line_col Color for the dashed horizontal line at \code{y = 1}
#'   (no change from reference).  Default \code{"grey25"}.
#' @param rtol Numeric. A numeric relative tolerance used by the ODE solver to determine
#'  if a good solution has been achieved. This is also used in the solved linear
#'  model to check if prior doses do not add anything to the solution. Default is 1e-6.
#' @param run_id Integer. Tested model ID. Default is 1.
#' @param ruv A list specifying options for residual error used in model fit, containing:
#'   * `YNAME` - string, output name. Typically `"y1"`, `"y2"`,...
#'   * `DVID` - numeric, observation type identifier corresponding to DVID column values
#'   * `TRANS` - string, residual error distribution. Can be: `"normal"`, `"logNormal"`, `"logitNormal"`
#'   * `PRED` - string, prediction variable name from the model
#'   * `ERR` - string, error model type. Options include:
#'     - "constant" for additive error
#'     - "proportional" for proportional error
#'     - "combined1" for combined additive and proportional error
#'   * `INIT` - numeric vector, initial values for error parameters (length depends on error model)
#'   * `EST` - logical vector, whether to estimate each error parameter (same length as INIT)
#'   * `BLQM` - below limit of quantification method (can be NULL)
#'  For several observation types list of lists should be provided:
#'   one residual error (ruv) object for each observations.
#' @param sc_factor Numeric. Scaling factor for DV/PRED/IPRED values. Default is 1 (no scaling).
#' @param scale A numeric named vector. Scaling for ODE parameters of the system.
#'  The names must correspond to the parameter identifiers in the ODE specification.
#'  Each of the ODE variables will be divided by the scaling factor. Default is `NULL`.
#' @param seed Integer. Random seed for synthetic data generation reproducibility. Default is \code{123}.
#' @param shp_i String. Column name for shape aesthetic. Default is `NULL`.
#' @param sigma A named matrix representing a sigma covariance matrix or its
#'  Cholesky decomposition; Default is `NULL`.
#' @param smooth Logical. Add LOESS smooth line. Default is `TRUE`.
#' @param stimes Vector of numeric. Sampling time points. Default is `NULL`.
#' @param task_opt String. Additional task options to be passed to the fitting software.
#'  For Monolix, this can include specific task configurations or optimization settings.
#'  When `NULL`, default tasks (populationParameters, individualParameters, fim, logLikelihood)
#'  will be used.  Default is `NULL`.
#' @param tdist Logical. If `TRUE`, overlay theoretical parameter distributions
#'  based on population mean and OMEGA matrix. Default is `TRUE`.
#' @param time_col String. The column to use as a time column. Currently, can be only `TIME`. Default is `TIME`.
#' @param theor_perc Logical. Show theoretical percentiles. Default is `TRUE`.
#' @param theor_percCI Logical. Show CI around theoretical percentiles. Default is `TRUE`.
#' @param theta A named vector or data frame. Values of population parameters to simulate with. Default is `NULL`.
#' @param thetamat A named variance-covariance matrix (for parameters brought to normal distribution).
#'  Default is `NULL`.
#' @param tsld Logical. If `TRUE`, uses time since last dose instead of time from first dose. Default is `FALSE`.
#' @param quantiles A numeric vector of length 2. Lower and upper quantiles of
#'   the continuous covariate distribution to test.
#'   Default is \code{c(0.1, 0.9)}.
#' @param val_col String. Name of value column. Default is `VALUE`.
#' @param wrap_i String. Faceting formula for `facet_wrap`. Default is `NULL`.
#' @param wrap_ncol Integer. Number of columns for `facet_wrap`. Default is `NULL`.
#' @param wrap_nrow Integer. Number of rows for `facet_wrap`. Default is `NULL`.
#' @param method Character string. `"PRCC"` or `"eFAST"`.
#' @param model A model object passed to `sg_sim()`.
#' @param params Character vector of parameter names to vary.

#' @param par_bounds A tibble or data frame with columns `PAR`, `LB`, `UB`.
#' @param point_size Point size for \code{geom_point}.  Default \code{2.5}.
#' @param n_sim Integer. Number of samples (LHS size for PRCC, base frequency size for eFAST).
#' @param stimes A numeric vector of simulation times.
#' @param output A character vector of outputs to keep. Passed to `sg_sim()`.
#' @param stat_comp A character vector of summary statistics to compute.
#'        Supported internally: `"mean","median","min","max","sd","cmax","SS"`.
#' @param et An event table passed to `sg_sim()`.
#' @param theta A named numeric vector of baseline parameters. Default is `NULL`.
#'        Parameters listed in `params` are replaced by sampled values.
#' @param cov Covariates passed to `sg_sim()`. Default is `NULL`.
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
  cov,
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
  DVID = 1,
  dv_col = "DV",
  emp_perc = TRUE,
  et,
  eta_seq = NULL,
  par_seq = NULL,
  par_type = "Ind",
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
  n_sim,
  no_leg = FALSE,
  npop = 1,
  nsub = 1,
  occ,
  opt_name = "Monolix",
  omega,
  outputs = NULL,
  output = NULL,
  path_to_fitter = NULL,
  path_to_save_output = NULL,
  par_bounds,
  params =NULL,
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
  stat_comp = NULL,
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

# Subset SDTAB-like data to one endpoint: numeric DVID, or DVNAME when present.
filter_sdtab_by_DVID <- function(ds, DVID = 1) {
  if (missing(DVID) || is.null(DVID)) {
    DVID <- 1
  }
  if (!"DVID" %in% names(ds)) {
    ep_num <- suppressWarnings(as.numeric(DVID))
    if (length(ep_num) == 1L && !is.na(ep_num) && ep_num == 1) {
      return(ds)
    }
    stop("SDTAB has no DVID column; cannot filter by DVID.")
  }

  if ((is.character(DVID) || is.factor(DVID)) && "DVNAME" %in% names(ds)) {
    ep_chr <- as.character(DVID)
    dvname_vals <- unique(stats::na.omit(as.character(ds$DVNAME)))
    if (ep_chr %in% dvname_vals) {
      out <- dplyr::filter(ds, as.character(.data$DVNAME) == ep_chr)
      if (nrow(out) == 0L) {
        stop("No SDTAB rows for DVNAME = '", ep_chr, "'.")
      }
      return(out)
    }
  }

  dvid_target <- suppressWarnings(as.numeric(DVID))
  if (length(dvid_target) != 1L || is.na(dvid_target)) {
    stop("DVID must be a single DVID (numeric) or DVNAME (string present in SDTAB$DVNAME).")
  }
  out <- dplyr::filter(ds, .data$DVID == dvid_target)
  if (nrow(out) == 0L) {
    stop("No SDTAB rows for DVID = ", dvid_target, ".")
  }
  out
}

read_smrg_ctrl <- function(ctrl) {
  if (inherits(ctrl, "character")) {
    if (!file.exists(ctrl)) {
      stop("File does not exist: ", ctrl)
    }

    ext <- tools::file_ext(ctrl)

    if (tolower(ext) == "rdata") {
      result <- get(load(ctrl))
    } else if (tolower(ext) == "json") {
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("Package 'jsonlite' is required for reading JSON files.")
      }
      result <- jsonlite::fromJSON(ctrl, simplifyVector = TRUE, simplifyDataFrame = TRUE)
    } else {
      stop("Unsupported file type: ", ext, ". Supported: .RData, .json")
    }

    return(result)
  } else if (inherits(ctrl, "list")) {
    return(ctrl)
  } else if (inherits(ctrl, "data.frame")) {
    stop("ctrl cannot be a data.frame. Provide either a file path (character) or a GCO list object")
  } else {
    stop("ctrl should be either a file path (character) or a GCO list object. Got: ", class(ctrl)[1])
  }
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
gmo_converter <- function(gco_path, output_path = NULL) {
  gco <- read_smrg_ctrl(gco_path)
  theta <- gco$theta
  covs <- gco$covs
  re <- gco$re
  ruv <- gco$ruv
  model_file <- gco$model
  data_path <- gco$data
  data <- utils::read.csv(data_path, check.names = FALSE)
  model_lines <- readLines(model_file)
  input_idx <- which(grepl("^# \\[INPUT\\]", model_lines))[1]
  model_idx <- which(grepl("^# \\[MODEL\\]", model_lines))[1]
  pre_input <- model_lines[seq_len(input_idx)]
  post_input <- model_lines[seq(model_idx, length(model_lines))]
  input_lines <- character()
  if (!is.null(theta) && nrow(theta) > 0) {
    for (i in seq_len(nrow(theta))) {
      name <- theta$NAME[i]
      trans <- tolower(theta$TRANS[i])
      init <- theta$INIT[i]
      if (is.na(init)) init <- 0
      if (trans == "lognormal") {
        input_lines <- c(input_lines, paste0(name, "_pop = log(", init, ");"))
      } else if (trans == "logitnormal") {
        input_lines <- c(input_lines, paste0(name, "_pop = logit(", init, ");"))
      } else {
        input_lines <- c(input_lines, paste0(name, "_pop = ", init, ";"))
      }
    }
    input_lines <- c(input_lines, "")
  }
  cov_stat_lines <- character()
  cov_par_lines <- character()
  if (!is.null(covs) && length(covs) > 0) {
    for (j in seq_along(covs)) {
      covj <- covs[[j]]
      par <- covj$PAR
      covname <- covj$COVNAME
      func <- covj$FUNC
      cov_trans <- covj$TRANS
      init <- if (is.null(covj$INIT)) 1 else covj$INIT
      if (!is.null(func) && func == "linear" && !is.null(cov_trans) && tolower(cov_trans) == "median") {
        med_name <- paste0(covname, "_med")
        med_val <- stats::median(data[[covname]], na.rm = TRUE)
        cov_stat_lines <- c(cov_stat_lines, paste0(med_name, " = ", med_val, ";"))
        beta_name <- paste0("beta_", par, "_", covname)
        cov_par_lines <- c(cov_par_lines, paste0(beta_name, " = ", init, ";"))
      } else if (!is.null(covj$REF)) {
        # Categorical covariate: create one beta per observed non-reference category value
        ref <- covj$REF
        cov_values <- unique(stats::na.omit(data[[covname]]))
        cov_values <- sort(cov_values)
        non_ref_vals <- cov_values[cov_values != ref]
        for (val in non_ref_vals) {
          beta_name <- paste0("beta_", par, "_", val)
          cov_par_lines <- c(cov_par_lines, paste0(beta_name, " = ", init, ";"))
        }
      }
    }
  }
  if (length(cov_stat_lines) > 0) {
    input_lines <- c(input_lines, cov_stat_lines, "")
  }
  if (length(cov_par_lines) > 0) {
    input_lines <- c(input_lines, cov_par_lines, "")
  }
  omega_names <- character()
  if (!is.null(re) && !is.null(re$est)) {
    est_mat <- re$est
    if (is.matrix(est_mat)) {
      pars <- colnames(est_mat)
      for (k in seq_along(pars)) {
        par_name <- pars[k]
        est_col <- est_mat[, k]
        if (any(!is.na(est_col) & est_col)) {
          oname <- paste0("omega_", par_name)
          omega_names <- c(omega_names, par_name)
          input_lines <- c(input_lines, paste0(oname, " = 1;"))
        }
      }
    }
  }
  if (length(omega_names) > 0) {
    input_lines <- c(input_lines, "")
  }
  if (!is.null(ruv) && !is.null(ruv$ERR)) {
    ruv_list <- list(ruv)
  } else if (is.list(ruv) && length(ruv) > 0 && !is.null(ruv[[1]]$ERR)) {
    ruv_list <- ruv
  } else {
    ruv_list <- list()
  }
  err_par_lines <- character()
  for (r in seq_along(ruv_list)) {
    ruvr <- ruv_list[[r]]
    pred <- ruvr$PRED
    err_type <- ruvr$ERR
    init <- ruvr$INIT
    if (is.null(err_type) || is.null(pred) || is.null(init)) next
    if (tolower(err_type) == "combined1" && length(init) >= 2) {
      a_name <- paste0(pred, "_a")
      b_name <- paste0(pred, "_b")
      err_par_lines <- c(err_par_lines, paste0(a_name, " = ", init[1], ";"))
      err_par_lines <- c(err_par_lines, paste0(b_name, " = ", init[2], ";"))
    } else if (tolower(err_type) == "proportional" && length(init) >= 1) {
      b_name <- paste0(pred, "_b")
      err_par_lines <- c(err_par_lines, paste0(b_name, " = ", init[1], ";"))
    }
  }
  if (length(err_par_lines) > 0) {
    input_lines <- c(input_lines, err_par_lines, "")
  }
  tv_lines <- character()
  if (!is.null(theta) && nrow(theta) > 0) {
    for (i in seq_len(nrow(theta))) {
      name <- theta$NAME[i]
      trans <- tolower(theta$TRANS[i])
      if (trans == "lognormal") {
        tv_lines <- c(tv_lines, paste0(name, "_tv = exp(", name, "_pop);"))
      } else if (trans == "logitnormal") {
        tv_lines <- c(tv_lines, paste0(name, "_tv = expit(", name, "_pop);"))
      } else {
        tv_lines <- c(tv_lines, paste0(name, "_tv = ", name, "_pop;"))
      }
    }
  }
  if (length(tv_lines) > 0) {
    input_lines <- c(input_lines, tv_lines, "")
  }
  par_def_lines <- character()
  if (!is.null(theta) && nrow(theta) > 0) {
    for (i in seq_len(nrow(theta))) {
      name <- theta$NAME[i]
      par_trans <- tolower(theta$TRANS[i])
      has_omega <- name %in% omega_names
      base_expr <- paste0(name, "_tv")
      if (has_omega) {
        oname <- paste0("omega_", name)
        if (par_trans == "lognormal") {
          base_expr <- paste0(base_expr, " * exp(", oname, ")")
        } else if (par_trans == "logitnormal") {
          base_expr <- paste0(base_expr, " * expit(", oname, ")")
        } else {
          base_expr <- paste0(base_expr, " + ", oname)
        }
      }
      if (!is.null(covs) && length(covs) > 0) {
        for (j in seq_along(covs)) {
          covj <- covs[[j]]
          if (!identical(covj$PAR, name)) next
          covname <- covj$COVNAME
          func <- covj$FUNC
          cov_trans <- covj$TRANS
          if (!is.null(func) && func == "linear" && !is.null(cov_trans) && tolower(cov_trans) == "median") {
            med_name <- paste0(covname, "_med")
            beta_name <- paste0("beta_", name, "_", covname)
            base_expr <- paste0(base_expr, " * (", covname, "/", med_name, ")^", beta_name)
          } else if (!is.null(covj$REF)) {
            # Categorical covariate
            ref <- covj$REF
            cov_values <- unique(stats::na.omit(data[[covname]]))
            cov_values <- sort(cov_values)
            non_ref_vals <- cov_values[cov_values != ref]
            if (length(non_ref_vals) == 0) next
            cat_terms <- character()
            for (val in non_ref_vals) {
              beta_name <- paste0("beta_", name, "_", val)
              cat_terms <- c(cat_terms, paste0(beta_name, " * (", covname, " == ", val, ")"))
            }
            if (length(cat_terms) > 0) {
              if (par_trans == "lognormal") {
                # Multiplicative effect on typical value: tv * exp(sum(beta * indicator))
                base_expr <- paste0(base_expr, " * exp(", paste(cat_terms, collapse = " + "), ")")
              } else {
                # Additive effect on natural scale for non-lognormal
                base_expr <- paste0(base_expr, " + ", paste(cat_terms, collapse = " + "))
              }
            }
          }
        }
      }
      par_def_lines <- c(par_def_lines, paste0(name, " = ", base_expr))
    }
  }
  all_input <- c(pre_input, input_lines, "", par_def_lines, "", post_input)
  res_err_lines <- character()
  for (r in seq_along(ruv_list)) {
    ruvr <- ruv_list[[r]]
    pred <- ruvr$PRED
    err_type <- ruvr$ERR
    if (is.null(err_type) || is.null(pred)) next
    if (tolower(err_type) == "combined1") {
      res_err_lines <- c(res_err_lines, paste0(pred, "_Res_err = ", pred, " * (1 + ", pred, "_b ) + ", pred, "_a;"))
    } else if (tolower(err_type) == "proportional") {
      res_err_lines <- c(res_err_lines, paste0(pred, "_Res_err = ", pred, " * (1 + ", pred, "_b );"))
    }
  }
  if (length(res_err_lines) > 0) {
    all_input <- c(all_input, res_err_lines)
  }
  if (is_null(output_path)) {
    return(all_input)
  } else {
    writeLines(all_input, con = output_path)
    return(all_input)
  }
}

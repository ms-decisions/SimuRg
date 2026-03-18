#' Build multistart scenarios with varying initial values
#'
#' @description
#' Generates multiple scenarios with identical model structure
#' (model, RUV, RE, OCC, COVS) and different sampled initial values for
#' population parameters.
#'
#' For each start, initial values of estimated parameters are sampled
#' uniformly within specified intervals, and a corresponding `.mlxtran` file
#' is generated. A summary table describing all multistart scenarios is written
#' to disk.
#'
#' @param mod Model file path or model identifier used for scenario generation.
#' @param data Character string. Path to the dataset used in model fitting or simulation.
#' @param headers Predefined list of dataset column names (e.g., ID, TIME, DV, etc.) used by the modeling framework.
#' @param ruv Residual unexplained variability (RUV) structure.
#'   Each element contains RUV properties (e.g., `YNAME`, `ERR`) and mapping to data.
#' @param theta Tibble describing population parameter properties.
#'   Tibble should contain columns such as `NAME`, `INIT`, and `EST`.
#' @param re Random effects (RE) specification.
#'   Includes matrices `init` and `est` defining initial and estimated covariance structures.
#' @param occ Occasion effect matrices.
#'   Includes `init` and `est` matrices defining inter-occasion variability.
#' @param covs_lst Optional list describing covariate-parameter relationships.
#'   Each element should include fields such as `PAR`, `COVNAME`, and `EST` indicating inclusion in the model.
#' @param n_starts Number of multistart runs (number of initial value scenarios).
#' @param theta_intervals Optional data frame specifying sampling
#'   intervals for initial values. Must contain columns
#'   `NAME`, `lower`, and `upper`.
#'   If `NULL`, intervals are derived from `interval_factor`.
#' @param interval_factor Numeric vector of length 2 specifying
#'   multiplicative lower and upper bounds factors relative to original
#'   `INIT` values (default `c(0.2, 5)`).
#' @param vary_params List of parameters which inital values are varied.
#' Default value implies all parameters are varied. If `theta_interval` is specified,
#' list of parameters is taken from `theta_interval`
#' @param opt_name Character. Optimization engine name (e.g., `"Monolix"`, `"Simurg"`). Default `"Simurg"`.
#' @param path Character. Directory path where output files (CSV summary and model files) will be written. Default is current working directory.
#' @param project_name Character. Base name for the exported project files. Default `"my_project"`.
#' @param seed Optional random seed for reproducible sampling of
#'   initial values. If `NULL`, the global RNG state is not modified.
#'
#' @return
#' The function writes two types of outputs to the specified `path`:
#' \itemize{
#'   \item A CSV file (`scenarios_info.csv`) summarizing all generated scenarios with columns:
#'     \describe{
#'       \item{scenario_number}{Unique index of the scenario.}
#'       \item{model}{Model file or path used in the scenario.}
#'       \item{data}{Path to the dataset used.}
#'       \item{theta}{Active population parameters and initial values used.}
#'       \item{RUV}{Residual error structure(s) used.}
#'       \item{RE}{Random effect and correlation terms included.}
#'       \item{OCC}{Active occasion effects.}
#'       \item{COVS}{Included covariate-parameter relationships.}
#'     }
#'   \item Individual `.mlxtran` model files for each generated scenario,
#'   named according to the pattern `<project_name>_<i>.mlxtran`.
#' }
#'
#' @details
#' For each start:
#' \itemize{
#'   \item Population parameter initial values `INIT` for
#'   estimated parameters are sampled uniformly.
#'   \item A scenario summary is stored.
#'   \item The corresponding `.mlxtran` project files are written.
#' }
#'
#' If `seed` is provided, the function temporarily sets the
#' random number generator state and restores it upon exit.
#'
#' @examples
#' \dontrun{
#' sg_multistart(
#'   mod_lst = "model1.txt",
#'   data = "data.csv",
#'   headers = list(ID = "ID", TIME = "TIME", DV = "DV"),
#'   ruv = ruv,
#'   theta = theta,
#'   re = re,
#'   occ = occ,
#'   covs_lst = covs_list,
#'   n_start = 15,
#'   path = "results/",
#'   project_name = "test_project"
#' )
#' }
#'
#' @import dplyr
#' @importFrom readr write_csv
#' @export
#' @importFrom jsonlite fromJSON
#' @importFrom stats runif
#' @import dplyr
#' @export

sg_multistart <- function(mod, data, headers, ruv, theta, re, occ,
                          covs_lst=NULL, opt_name = "Simurg",
                          path=getwd(), project_name="multistart_project",
                          n_starts=10, theta_intervals=NULL,
                          interval_factor=c(0.2,5), vary_params=NULL,
                          seed = NULL) {

  make_theta_one <- function(theta_lst, val, trans_idx, init_idx) {
    for (i in c(1:length(trans_idx))) {
      theta_lst$TRANS[trans_idx[i]] <- val[i]
    }
    for (i in c(1:length(init_idx))) {
      theta_lst$INIT[init_idx[i]] <- val[i + length(trans_idx)]
    }
    return(theta_lst %>% mutate(TRANS = as.character(TRANS),
                                INIT = as.numeric(INIT)))
  }
  make_theta_lst <- function(theta_lst) {
    theta_lst_ <- theta_lst %>% group_by(NAME) %>%
      mutate(TRANS_ln = length(TRANS %>% unlist()),
             INIT_ln = length(INIT %>% unlist()))
    trans_idx <- which(theta_lst_$TRANS_ln > 1)
    init_idx <- which(theta_lst_$INIT_ln > 1)
    if (length(trans_idx) + length(init_idx) == 0) {
      return(list(theta_lst))
    }
    trans_val <- expand.grid(c(theta_lst$TRANS[trans_idx],
                               theta_lst$INIT[init_idx]),
                             stringsAsFactors = F)
    return(apply(trans_val, 1,function(x){make_theta_one(theta_lst, x, trans_idx, init_idx)}))
  }

  # Setup local seed and save user's global seed
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    old_seed <- .GlobalEnv$.Random.seed
    has_old_seed <- TRUE
  } else {
    has_old_seed <- FALSE
  }

  on.exit({
    if (has_old_seed) {
      .GlobalEnv$.Random.seed <- old_seed
    } else {
      if (exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }
  })

  if (!is.null(seed)) set.seed(seed)

  # Setup variation parameters for multistart
  if (is.null(vary_params)) {
    vary_params <- theta$NAME
  }
  if (is.null(theta_intervals)) {
    theta_intervals <- list()
    for (p in vary_params) {
      init_val <- theta$INIT[theta$NAME == p]
      theta_intervals[[p]] <- c(init_val * interval_factor[1], init_val * interval_factor[2])
    }
  }

  # Sample n_starts initial theta values
  theta_samples <- list()
  for (p in names(theta_intervals)) {
    bounds <- theta_intervals[[p]]
    theta_samples[[p]] <- runif(n_starts, bounds[1], bounds[2])
  }
  # Create list of theta tables with varied INIT for chosen parameters
  theta_list <- vector("list", n_starts)
  for (i in seq_len(n_starts)) {
    th <- theta
    for (p in names(theta_samples)) {
      th$INIT[th$NAME == p] <- theta_samples[[p]][i]
    }
    theta_list[[i]] <- th
  }
  # Expand if any transformation or list values
  theta_lst_exp <- lapply(theta_list, function(x) make_theta_lst(x)[[1]])

  summary_list <- vector("list", n_starts)
  for (i in seq_len(n_starts)) {
    model_path <- mod
    # Theta info string of EST parameters
    theta_true <- theta_lst_exp[[i]] %>%
      mutate(NAME_INIT = paste0(NAME,'=',INIT)) %>%
      filter(EST == TRUE) %>%
      pull(NAME_INIT) %>%
      str_c(collapse = ", ")
    # RUV info
    ruv_info <- ruv
    ruv_full <- paste0(ruv_info$YNAME, "_", ruv_info$ERR)
    # RE info
    re_mat <- re$init
    re_names <- colnames(re_mat)
    diag_idx <- which(diag(re_mat) == 1)
    re_active <- if (length(diag_idx) > 0) re_names[diag_idx] else character(0)
    corr_list <- c()
    for (d_idx in diag_idx) {
      diag_name <- re_names[d_idx]
      for (c in seq_len(ncol(re_mat))) {
        if (c != d_idx && re_mat[d_idx, c] == 1) {
          corr_list <- c(corr_list, paste0("corr_", diag_name, "_", re_names[c]))
        }
      }
      for (r in seq_len(nrow(re_mat))) {
        if (r != d_idx && re_mat[r, d_idx] == 1) {
          corr_list <- c(corr_list, paste0("corr_", diag_name, "_", re_names[r]))
        }
      }
    }
    re_full <- paste(unique(c(re_active, corr_list)), collapse = ", ")
    # OCC info
    occ_mat <- occ$init
    occ_active <- colnames(occ_mat)[colSums(occ_mat == 1, na.rm = TRUE) > 0] %>%
      str_c(collapse = ", ")
    # COVS info
    covs_active <- covs_lst %>%
      keep(~ .x$EST == TRUE) %>%
      map_chr(~ str_c(.x$PAR, "_", .x$COVNAME)) %>%
      str_c(collapse = ", ")
    summary_list[[i]] <- tibble(
      scenario_number = i,
      model = model_path,
      data = data,
      theta = theta_true,
      RUV = ruv_full,
      RE = re_full,
      OCC = occ_active,
      COVS = covs_active
    )
  }
  summary_df <- bind_rows(summary_list)
  write_csv(summary_df, paste0(path, '/multistart_info.csv'))

  for (i in seq_len(n_starts)) {
    sg_result_mod <- sg_fit(mod, data, headers, theta_lst_exp[[i]],
                            ruv, re, occ,
                            project_name = paste0(project_name,'_',i),
                            covs_lst, opt_name = "Monolix",
                            path_to_save_output = path)
    write(sg_result_mod, str_c(path, paste0(project_name,'_',i),'.mlxtran'))
  }
}


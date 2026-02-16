## Author: Ugolkov Yaroslav, Victoria Kulesh, Mikhailova Anna
## First created: 2025-10-17
## Description: Local sensitivity analysis simulations for RxODE models
## Keywords: SimuRg, model build

#' Build estimation model scenarios
#'
#' @description
#' Constructs and exports multiple model configuration scenarios by combining population parameters,
#' random effects, residual variability, occasion effects, and covariate assignments.
#' Each scenario is defined by unique combinations of these model components and exported as `.mlxtran` files.
#' Additionally, a summary CSV file describing all scenarios is generated.
#'
#' @param mod_lst List of model file paths or model identifiers used for scenario generation.
#' @param data Character string. Path to the dataset used in model fitting or simulation.
#' @param headers Predefined list of dataset column names (e.g., ID, TIME, DV, etc.) used by the modeling framework.
#' @param ruv_lst List specifying residual unexplained variability (RUV) structures.
#'   Each element contains RUV properties (e.g., `YNAME`, `ERR`) and mapping to data.
#' @param theta_lst List of tibbles describing population parameter properties.
#'   Each tibble should contain columns such as `NAME`, `INIT`, and `EST`.
#' @param re_lst List of random effects (RE) specifications.
#'   Each element includes matrices `init` and `est` defining initial and estimated covariance structures.
#' @param occ_lst List of occasion effect matrices.
#'   Each element includes `init` and `est` matrices defining inter-occasion variability.
#' @param covs_lst Optional list describing covariate-parameter relationships.
#'   Each element should include fields such as `PAR`, `COVNAME`, and `EST` indicating inclusion in the model.
#' @param task_lst Optional list defining additional tasks or configuration options for each scenario. Default is NULL.
#' @param opt_name Character. Optimization engine name (e.g., `"Monolix"`, `"Simurg"`). Default `"Simurg"`.
#' @param path Character. Directory path where output files (CSV summary and model files) will be written. Default is current working directory.
#' @param project_name Character. Base name for the exported project files. Default `"my_project"`.
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
#' This function automates the construction of multiple model configurations for model evaluation,
#' sensitivity analysis, or scenario testing. It expands parameter and variability definitions,
#' constructs all possible combinations, and exports model specifications for external fitting tools.
#'
#' @examples
#' \donttest{
#' sg_modbuild(
#'   mod_lst = list("model1.txt", "model2.txt"),
#'   data = "data.csv",
#'   headers = list(ID = "ID", TIME = "TIME", DV = "DV"),
#'   ruv_lst = ruv_lst_1,
#'   theta_lst = theta_list,
#'   re_lst = re_list,
#'   occ_lst = occ_list,
#'   covs_lst = covs_list,
#'   path = "results/",
#'   project_name = "test_project"
#' )
#' }
#'
#' @import dplyr
#' @importFrom readr write_csv
#' @export
#' @importFrom jsonlite fromJSON
#' @import dplyr
#' @export

sg_modbuild <- function(mod_lst, data, headers, ruv_lst, theta_lst, re_lst,
                        occ_lst, covs_lst=NULL, task_lst = NULL, opt_name = "Simurg",
                        path=getwd(), project_name = "my_project") {


  make_re_obj <- function(changes, est, init, variable_num, init0) {
    for (j in changes) {
      if (j == 0) {
        next
      }
      init[variable_num[j, 'row'],
           variable_num[j, 'col']]  <- init0[variable_num[j, 'row'],
                                             variable_num[j, 'col']]
      est[variable_num[j, 'row'],
          variable_num[j, 'col']]  <- T
    }
    return(list(init = init, est = est))
  }
  make_re_lst <- function(re_lst) {
    variable <- re_lst[["est"]] & !re_lst[["block"]]
    variable_num <- which(variable, arr.ind=TRUE)
    variants <- combn(c(rep(0, times = nrow(variable_num)), 1:nrow(variable_num)),
                      m = nrow(variable_num), simplify = F) %>% unique()
    init <- ifelse(!variable | is.na(variable), re_lst[["init"]], 0)
    est <- ifelse(variable, NA, !variable)
    fin_re_lst <- map(variants, function(x){
      make_re_obj(x, est, init, variable_num, re_lst[["init"]])})
    return(fin_re_lst)
  }

  make_ruv_one <- function(i, j, ruv_lst) {
    ruv_lst$EST <- ruv_lst$EST[[i]]
    ruv_lst$INIT <- ruv_lst$INIT[[i]]
    ruv_lst$ERR <- ruv_lst$ERR[[i]]
    ruv_lst$TRANS <- ruv_lst$TRANS[[j]]
    return(ruv_lst)
  }
  make_ruv_lst <- function(ruv_lst) {
    i_s <- c(1:length(ruv_lst$ERR))
    j_s <- c(1:length(ruv_lst$TRANS))
    inxs <- expand_grid(i_s, j_s)
    return(pmap(inxs, function(i_s, j_s) {make_ruv_one(i_s, j_s, ruv_lst)}))
  }
  make_ruv_all <- function(ruv_lst) {
    if (length(ruv_lst) == 1){
      return(make_ruv_lst(ruv_lst[[1]]))
    } else {
      return(map(ruv_lst, make_ruv_lst))
    }
  }

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


  theta_lst_exp <- map(theta_lst, make_theta_lst)
  re_lst_exp <- map(re_lst, make_re_lst)
  occ_lst_exp <- map(occ_lst, make_re_lst)
  ruv_lst_exp <- map(ruv_lst, function(x)  {map(x, make_ruv_lst)})

  sc_lst <- tibble()
  for (i in c(1:length(mod_lst[1]))) {
    if (length(ruv_lst_exp[[i]]) > 1) {
      ruv_lst_sev_dv <- do.call(expand.grid, ruv_lst_exp[[i]]) %>% pmap(., list) # Several DVID
    } else{
      ruv_lst_sev_dv <- ruv_lst_exp[[i]][[1]]
    }

    all_variants <- expand_grid(theta_lst_exp[[i]], re_lst_exp[[i]],
                                ruv_lst_sev_dv, occ_lst_exp[[i]])
    colnames(all_variants) <- c('theta', 're', 'ruv', 'occ')
    all_variants['mod'] <- mod_lst[[i]]
    sc_lst <- rbind(sc_lst, all_variants)
  }
  # return(sc_lst)
  n_iter <- length(sc_lst$theta)
  summary_list <- vector("list", n_iter)

  for (i in seq_len(n_iter)) {
    model_path <- sc_lst$mod[i]
    # --- theta ---
    theta_true <- sc_lst$theta[[i]] %>% mutate(NAME_INIT = paste0(NAME,'=',INIT)) %>%
      filter(EST == TRUE) %>%
      pull(NAME_INIT) %>%
      str_c(collapse = ", ")

    # --- RUV ---
    ruv_info <- sc_lst$ruv[[i]]
    ruv_full <- paste0(ruv_info$YNAME, "_", ruv_info$ERR)

    # --- RE ---
    re_mat <- sc_lst$re[[i]]$init
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

    # --- OCC ---
    occ_mat <- sc_lst$occ[[i]]$init
    occ_active <- colnames(occ_mat)[colSums(occ_mat == 1, na.rm = TRUE) > 0] %>%
      str_c(collapse = ", ")

    # --- COVS ---
    covs_active <- covs_lst %>%
      keep(~ .x$EST == TRUE) %>%
      map_chr(~ str_c(.x$PAR, "_", .x$COVNAME)) %>%
      str_c(collapse = ", ")

    # --- tibble ---
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
  write_csv(summary_df, paste0(path, '/scenarios_info.csv'))
  for (i in seq(1,length(sc_lst$mod))){
    sg_result_mod <- sg_fit(sc_lst$mod[i], data, headers, sc_lst$theta[[i]],
                            sc_lst$ruv[[i]], sc_lst$re[[i]],  sc_lst$occ[[i]],
                            project_name = paste0(project_name,'_',i),
                            covs_lst, opt_name = "Monolix",
                            path_to_save_output = path)
    write(sg_result_mod, str_c(path, paste0(project_name,'_',i),'.mlxtran'))
  }
}


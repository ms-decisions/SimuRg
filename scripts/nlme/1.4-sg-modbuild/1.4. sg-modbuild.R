## Author: Anna Mikhailova
## First created: 2025-08-14
## Description: the function to test different scenarios of models
## Keywords: simurg, model build

#####--------------- Load functions and libraries ---------------#####
library(tidyverse)
library(jsonlite)

source("V:/Collaborative_working/Simurg_R_package/Rproj_4.2.3_renv_SimurgR/functions/nlme/1.4-sg-modbuild/1.4-sg-modbuild-args.R", echo=TRUE)

re_lst <- re_lst_2[[1]]
# make_re_obj <- function(changes, est, init, variable_num, init0) {
#   for (j in changes) {
#     if (j == 0) {
#       next
#     }
#     init[variable_num[j, 'row'],
#          variable_num[j, 'col']]  <- init0[variable_num[j, 'row'],
#                                                       variable_num[j, 'col']]
#     est[variable_num[j, 'row'],
#         variable_num[j, 'col']]  <- T
#   }
#   return(list(init = init, est = est))
# }
# make_re_lst <- function(re_lst) {
#   variable <- re_lst[["est"]] & !re_lst[["block"]]
#   variable_num <- which(variable, arr.ind=TRUE)
#   variants <- combn(c(rep(0, times = nrow(variable_num)), 1:nrow(variable_num)),
#                     m = nrow(variable_num), simplify = F) %>% unique()
#   init <- ifelse(!variable | is.na(variable), re_lst[["init"]], 0)
#   est <- ifelse(variable, NA, !variable)
#   fin_re_lst <- map(variants, function(x){
#     make_re_obj(x, est, init, variable_num, re_lst[["init"]])})
#   return(fin_re_lst)
# }
#
# ruv_lst <- ruv_lst_1[[1]]
# make_ruv_one <- function(i, j, ruv_lst) {
#   ruv_lst$EST <- ruv_lst$EST[[i]]
#   ruv_lst$INIT <- ruv_lst$INIT[[i]]
#   ruv_lst$ERR <- ruv_lst$ERR[[i]]
#   ruv_lst$TRANS <- ruv_lst$TRANS[[j]]
#   return(ruv_lst)
# }
# make_ruv_lst <- function(ruv_lst) {
#   i_s <- c(1:length(ruv_lst$ERR))
#   j_s <- c(1:length(ruv_lst$TRANS))
#   inxs <- expand_grid(i_s, j_s)
#   return(pmap(inxs, function(i_s, j_s) {make_ruv_one(i_s, j_s, ruv_lst)}))
# }
#
# theta_lst <- theta_lst_2[[2]]
# make_theta_one <- function(theta_lst, val, trans_idx, init_idx) {
#   for (i in c(1:length(trans_idx))) {
#     theta_lst$TRANS[trans_idx[i]] <- val[i]
#   }
#   for (i in c(1:length(init_idx))) {
#     theta_lst$INIT[init_idx[i]] <- val[i + length(trans_idx)]
#   }
#   return(theta_lst %>% mutate(TRANS = as.character(TRANS),
#                               INIT = as.numeric(INIT)))
# }
# make_theta_lst <- function(theta_lst) {
#   theta_lst_ <- theta_lst %>% group_by(NAME) %>%
#     mutate(TRANS_ln = length(TRANS %>% unlist()),
#            INIT_ln = length(INIT %>% unlist()))
#   trans_idx <- which(theta_lst_$TRANS_ln > 1)
#   init_idx <- which(theta_lst_$INIT_ln > 1)
#   if (length(trans_idx) + length(init_idx) == 0) {
#     return(list(theta_lst))
#   }
#   trans_val <- expand.grid(c(theta_lst$TRANS[trans_idx],
#                              theta_lst$INIT[init_idx]),
#                            stringsAsFactors = F)
#   return(apply(trans_val, 1,function(x){make_theta_one(theta_lst, x, trans_idx, init_idx)}))
# }
theta_lst <- theta_lst_2
re_lst <- re_lst_2
ruv_lst <- ruv_lst_2
occ_lst <- re_lst_2
sg_modbuild <- function(mod_lst, data, headers, ruv_lst, theta_lst, re_lst,
                        occ_lst, covs_lst=NULL, task_lst = NULL, opt_name = "Simurg",path=getwd(),project_name = "my_project") {


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
  write_csv(summary_df, paste0(path, 'scenarios_info.csv'))
  for (i in seq(1,length(sc_lst$mod))){
    sg_result_mod <- sg_fit(sc_lst$mod[i], data, headers, sc_lst$theta[[i]], sc_lst$ruv[[i]], sc_lst$re[[i]],  sc_lst$occ[[i]], project_name = paste0(project_name,'_',i), covs_lst, opt_name = "Monolix")
    write(sg_result_mod, str_c(path, paste0(project_name,'_',i),'.mlxtran'))
  }
}

folder_path <- str_c(getwd(), "/scripts/nlme/1.4-sg-modbuild/")

### list of paths to structural models
mod_lst <- list(str_c(folder_path, "models/model_1c.txt"),
                str_c(folder_path, "models/model_2c.txt"))

### path to the dataset
data <- str_c(folder_path, "data/dspk-warf.csv")
re_lst_1 <- list(
  list(init = tribble(~Cl, ~Vd, ~ka, ~Vp, ~Q,
                      1, 0, 0, 0, 0,
                      0, 1, 0, 0, 0,
                      0, 0, 1, 0, 0,
                      0, 0, 0, 1, 0,
                      0, 0, 0, 0, 1) %>% as.matrix(),
       est = tribble(~Cl, ~Vd, ~ka, ~Vp, ~Q,
                     T, NA, NA, NA, NA,
                     NA, T, NA, NA, NA,
                     NA, NA, T, NA, NA,
                     NA, NA, NA, T, NA,
                     NA, NA, NA, NA, T) %>% as.matrix(),
       block =tribble(~Cl, ~Vd, ~ka, ~Vp, ~Q,
                      F, NA, NA, NA, NA,
                      NA, F, NA, NA, NA,
                      NA, NA, T, NA, NA,
                      NA, NA, NA, F, NA,
                      NA, NA, NA, NA, F) %>% as.matrix())
)
headers <- list(list(name = "ID", use = "identifier", type = NULL),
                list(name = "TIME", use = "time", type = NULL),
                list(name = "DV", use = "observation", type = "continuous"),
                list(name = "DVID", use = "observationtype", type = NULL),
                list(name = "ADM", use = "administration", type = NULL),
                list(name = "AMT", use = "amount", type = NULL),
                list(name = "EVID", use = "eventidentifier", type = NULL),
                list(name = "MDV", use = "missingdependentvariable", type = NULL),
                list(name = "AGE", use = "covariate", type = "continuous"),
                list(name = "AGE_centered", use = "covariate", type = "continuous"),
                list(name = "SEX", use = "covariate", type = "categorical"),
                list(name = "WEIGHT", use = "covariate", type = "continuous"),
                list(name = "BMI", use = "covariate", type = "continuous"))

ruv_lst_2 <- list(
  # structural model 1
  list(
    list(YNAME = "y1",
         DVID = 1,
         TRANS = "normal",
         PRED = "Cc",
         ERR = list("constant", "proportional", "combined1"), # options to test (can be length = 1, i.e., "constant" or "combined1")
         INIT = list(1, 1, c(1, 1)), # options to test (can be length = 1, i.e., 1 or c(1, 1))
         EST = list(T, T, c(T, T))) # options to test (can be length = 1, i.e., T or c(T, T))BLQM = NULL))
  ),
  # structural model 2
  list(
    list(YNAME = "y1",
         DVID = 1,
         TRANS = "normal",
         PRED = "Cc",
         ERR = list("constant", "proportional"), # options to test (can be length = 1, i.e., "constant" or "combined1")
         INIT = list(1, 1), # options to test (can be length = 1, i.e., 1 or c(1, 1))
         EST = list(T, T)) # options to test (can be length = 1, i.e., T or c(T, T))BLQM = NULL))
  )
)

theta_lst_2 <- list(
  # structural model 1
  tribble(~NAME, ~TRANS, ~INIT, ~LB, ~UB, ~EST,
          "Cl", "logNormal", 0.2, NA, NA, T,
          "Vd", "logNormal", 20, NA, NA, T,
          "ka", "logNormal", 0.2, NA, NA, T),
  # structural model 2
  tribble(~NAME, ~TRANS, ~INIT, ~LB, ~UB, ~EST,
          "Cl", "logNormal", 0.2, NA, NA, T,
          "Vd", c("Normal", "logNormal"), c(10, 20, 30), NA, NA, T,
          "ka", "logNormal", c(0.2, 0.1, 0.5), NA, NA, T, # INIT could be length > 1
          "Vp", c("Normal", "logNormal"), 10, NA, NA, T,
          "Q", "logNormal", 5, NA, NA, T))


sg_modbuild(
     mod_lst = mod_lst[1],
     data = data,
     headers = headers,
     ruv_lst = ruv_lst_2,
     theta_lst = theta_lst_2,
     re_lst = re_lst_1,
     occ_lst = re_lst_1,
     covs_lst = NULL,
     path = "V:/Collaborative_working/SimuRg_as_R_lib/SimuRg/scripts/nlme/1.4-sg-modbuild/results/",
     project_name = "test_project"
   )

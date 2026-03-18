## Author: Alina Melnikova
## First created: 2026-01-28
## Description: Simulation for covariate sensitivity analysis
## Keywords: covariates, simulations, sensitivity

###-----Functions----####

funSum_sim <- list(mean   = ~mean(., na.rm = TRUE),
                   median = ~median(., na.rm = TRUE),
                   min    = ~min(., na.rm = TRUE),
                   max    = ~max(., na.rm = TRUE),
                   sd     = ~sd(., na.rm = TRUE),
                   P025   = ~quantile(., 0.025, na.rm = TRUE),
                   P05    = ~quantile(., 0.05,  na.rm = TRUE),
                   P10    = ~quantile(., 0.10,  na.rm = TRUE),
                   P15    = ~quantile(., 0.15,  na.rm = TRUE),
                   P25    = ~quantile(., 0.25,  na.rm = TRUE),
                   P75    = ~quantile(., 0.75,  na.rm = TRUE),
                   P85    = ~quantile(., 0.85,  na.rm = TRUE),
                   P90    = ~quantile(., 0.90,  na.rm = TRUE),
                   P95    = ~quantile(., 0.95,  na.rm = TRUE),
                   P975   = ~quantile(., 0.975, na.rm = TRUE),
                   geom_mean = ~exp(mean(suppressWarnings(log(.)), na.rm = TRUE)),
                   CV     = ~sd(., na.rm = TRUE)/mean(., na.rm = TRUE)*100)

funSum_av <- list(mean   = ~mean(.),
                  median   = ~median(., na.rm = T))

funSum_exp <- list(`Cavg`   = ~mean(., na.rm = TRUE),
                   `Cmin`    = ~min(., na.rm = TRUE),
                   `Cmax`    = ~max(., na.rm = TRUE))

fun_EtCC <- function(et_base_i, cc_ds_i, cat = F){
  et_scov_i <- unique(cc_ds_i$COV) %>% map(function(n){
    cc_ds_n <- cc_ds_i %>% filter(COV == n)

    et_scov_n <- seq(nrow(cc_ds_n)) %>% map_dfr(function(m){
      row_m <- cc_ds_n %>% filter(row_number() == m)

      et_scov_m <- et_base_i %>%
        mutate_at(vars(all_of(row_m$COV)), function(k){k = row_m$COVVAL})
      if(!cat){
        et_scov_m <- et_scov_m %>% bind_cols(select(row_m, COV, PAR, BTR, KEY, COVVAL, BCOVVAL))
      } else {
        et_scov_m <- et_scov_m %>% bind_cols(select(row_m, COV, PAR, KEY, COVVAL, CATDES))
      }

      return(et_scov_m)

    }) %>% mutate(id = row_number())

    return(et_scov_n)
  })
  names(et_scov_i) <- unique(cc_ds_i$COV)
  return(et_scov_i)
}

fun_CovSens <- function(et_sim_i, cat = F, expos = F, covs_i = NULL, nsim = 100, stime_exp = NULL,
                        mod_fin_i,
                        theta_i, thetamat_i, nice_names_i,
                        var_exp = "Cc", aggr = c("min", "max", "mean")) #nsim=1000
{


  #Test
  # et_sim_i = ets_cc
  # covs_i = nice_names$COV
  # expos = T
  # stime_exp = stimes_ss
  # cat = F; nsim = 5;
  # var_exp = "Cc"
  # theta_i = par_fin_tv
  # thetamat_i = m_theta_norm_pop


  keep_i <- c("Regimen", "KEY", "COV", "COVVAL")
  if(!cat){
    keep_i <- c(keep_i, "BTR", "BCOVVAL")
  } else {
    keep_i <- c(keep_i, "CATDES")
  }

  sens_i <- et_sim_i %>% map_dfr(function(et_i){
    #Test
    #et_i <- et_sim_i[[1]]

    if(!expos){
      par_i <- unique(et_i$PAR)
      if(is.list(par_i)){ par_i <- par_i[[1]] }
      sim_i <- sg_sim(mod_fin_i, et_i, 0, outputs = par_i, theta = theta_i ,
                      thetamat = thetamat_i, covs = covs_i, npop = nsim, keep = keep_i)
    } else {
      sim_raw <- sg_sim(mod_fin_i, et_i, stime_exp, outputs = var_exp, theta = theta_i ,
                        thetamat = thetamat_i, covs = covs_i, npop = nsim,
                        keep = keep_i, ncores = max(1, parallel::detectCores()-2))

      message(sprintf("NA Cc rows: %d / %d (%.1f%%)",
                      sum(is.na(sim_raw$VALUE)), nrow(sim_raw),
                      100 * mean(is.na(sim_raw$VALUE))))

      aggr_map <- c("mean" = "Cavg", "min" = "Cmin", "max" = "Cmax")
      stopifnot("aggr must only contain 'mean', 'min', 'max'" = all(aggr %in% names(aggr_map)))
      funSum_exp_i <- funSum_exp[unname(aggr_map[aggr])]

      sim_i <- sim_raw %>% mutate(VALUE = as.numeric(VALUE))

      sim_i <- sim_i %>%
        # group_by_at(vars(all_of(c("POPN", keep_i, covs_i)))) %>% summarise_at(vars(VALUE), funSum_exp_i) %>% ungroup() %>%
        # gather("VAR", "VALUE", -all_of(c("POPN", keep_i, covs_i)))
      group_by_at(vars(all_of(c("POPN", "VAR", keep_i, covs_i)))) %>% summarise_at(vars(VALUE), funSum_exp_i) %>% ungroup() %>%
        gather("METRIC", "VALUE", -all_of(c("POPN", "VAR", keep_i, covs_i))) %>%
        mutate(VAR = paste(VAR, METRIC, sep = "_")) %>%
        select(-METRIC)

    }

    sim_i <- sim_i %>% mutate(VALUE = as.numeric(VALUE))

    sim_i_ref <- sim_i %>% filter(KEY == "REF" | KEY == 1) %>% select(POPN, Regimen, VAR, REFVAL = VALUE)
    sim_i_ch <- sim_i %>% filter(KEY != "REF" & KEY != 1) %>% left_join(sim_i_ref, by = c("POPN", "Regimen", "VAR")) %>%
      mutate(PCH = 100*(VALUE - REFVAL)/REFVAL)


    out_i <- sim_i_ch %>% group_by_at(vars(all_of(c("VAR", keep_i, covs_i)))) %>% summarise_at(vars(PCH), funSum_sim) %>% ungroup()

    # Ensure BCOVVAL is added from et_i if not already present (for continuous covariates)
    if(!cat && "BCOVVAL" %in% names(et_i)){
      if(!"BCOVVAL" %in% names(out_i)){
        bcovval_lookup <- et_i %>% select(KEY, COV, COVVAL, BCOVVAL) %>% unique()
        out_i <- out_i %>% left_join(bcovval_lookup, by = c("KEY", "COV", "COVVAL"))
      }
    }
    out_i

  })  %>% left_join(nice_names_i, by = "COV")

  if(!cat){
    sens_out <- sens_i %>% select(NICEN, VAR, KEY, mean:P975, Regimen, COVVAL:BCOVVAL) %>%
      mutate(KEY = ifelse(KEY == "LP", paste0(quantiles[[1]]*100, "th perc."), paste0(quantiles[[2]]*100, "th perc.")), LAB = str_c(NICEN, "\n", KEY, " (", round(BCOVVAL, 1), ")"))
  } else {
    sens_out <- sens_i %>% select(NICEN, VAR, CATDES, mean:P975, Regimen) %>% mutate(LAB = str_c(NICEN, " (", CATDES, ")"))
  }
  sens_out <- sens_out %>% mutate_at(vars(mean:P975), function(p){p/100+1})
  return(sens_out)
}
#####

sg_covsens_sim <- function(fpath_i = NULL, ds_parest = NULL, ds_cov = NULL,
                           mod_fin, stimes_ss, ev_t_input,
                           est_covmat,
                           Nsim = 5,
                           cont_cov_l, cat_cov_l,  quantiles = c(0.1, 0.9),
                           var_output  = "Cc", aggr = c("min", "max", "mean")){

  # --- Input validation ---
  # Data source: must provide either fpath_i alone, or both ds_parest and ds_cov
  has_fpath   <- !is.null(fpath_i)
  has_ds      <- !is.null(ds_parest) && !is.null(ds_cov)
  has_ds_part <- !is.null(ds_parest) || !is.null(ds_cov)

  if (!has_fpath && !has_ds) {
    stop(
      "No data source provided. Supply either:\n",
      "  - 'fpath_i': path to a Simurg output object, OR\n",
      "  - both 'ds_parest' (parameter estimates) and 'ds_cov' (covariate dataset)."
    )
  }
  if (has_fpath && has_ds_part) {
    stop(
      "'fpath_i' and 'ds_parest'/'ds_cov' are mutually exclusive. ",
      "Supply one data source only."
    )
  }
  if (!has_fpath && has_ds_part && !has_ds) {
    missing_ds <- if (is.null(ds_parest)) "'ds_parest'" else "'ds_cov'"
    stop(
      "Incomplete data source: ", missing_ds, " is missing. ",
      "Both 'ds_parest' and 'ds_cov' must be provided together."
    )
  }

  if (missing(mod_fin))    stop("'mod_fin' is required: provide the compiled rxode2/nlmixr model object.")
  if (missing(stimes_ss))  stop("'stimes_ss' is required: provide the steady-state simulation time points.")
  if (missing(ev_t_input)) stop("'ev_t_input' is required: provide the event table (dosing schedule).")
  if (missing(est_covmat)) stop("'est_covmat' is required: provide the parameter estimation covariance matrix.")
  if (missing(cont_cov_l)) stop("'cont_cov_l' is required: provide the continuous covariate definition list.")
  if (missing(cat_cov_l))  stop("'cat_cov_l' is required: provide the categorical covariate definition list.")

  # Warn about non-default aggregation choices to alert on typos
  valid_aggr <- c("min", "max", "mean")
  bad_aggr   <- setdiff(aggr, valid_aggr)
  if (length(bad_aggr) > 0) {
    warning("Unrecognised aggregation function(s) in 'aggr': ",
            paste(bad_aggr, collapse = ", "),
            ". Valid options are: ", paste(valid_aggr, collapse = ", "), ".")
  }
  # -------------------------

  if((!is.null(fpath_i))&(is.null(ds_parest))){
  obj_data <- read_smrg_obj(fpath_i)
  #ds_mod <-  obj_data$SDTAB
  par_sum <- obj_data$SUMTAB
  ds_catcov <- obj_data$CATAB
  ds_ccov <- obj_data$COTAB
  data_fin <- ds_ccov %>% left_join(ds_catcov, by = "ID")
  ### Population parameters
  par_fin <- par_sum %>% rename(parameter = PAR, value = VALUE)
  } else if((is.null(fpath_i))&(!is.null(ds_parest))&(!is.null(ds_cov))){
    par_fin <- ds_parest
    data_fin <- ds_cov
  }

  par_pop <- par_fin %>% filter(str_detect(parameter, "_pop")) %>% select(parameter, value) %>% deframe()
  par_fin_tv <- par_fin %>% filter(str_detect(parameter, "_pop$") | str_detect(parameter, "^beta_")) %>% select(parameter, value) %>%
    mutate(value = ifelse(str_detect(parameter, "_pop$") & !str_detect(parameter, "^beta_"), log(value), value)) %>% deframe()

  ### Reconstruct omega matrix (random effects on K_a and V_pop)
  d_omega <- par_fin %>% filter(str_detect(parameter, "omega_"))
  m_omega <- diag(d_omega$value, ncol = length(d_omega$value))
  colnames(m_omega) <- d_omega$parameter; rownames(m_omega) <- d_omega$parameter
  m_omega_full <- m_omega

  corr_set <- par_fin %>% filter(str_detect(parameter, "corr_"))
  if(nrow(corr_set) > 0){
    corr_set <- corr_set %>%
      separate(parameter, c("TYPE", "PAR1", "PAR2"), "_") %>%
      select(PAR1, PAR2, value) %>%
      mutate_at(vars(PAR1, PAR2), function(x){str_c("omega_", x)})

    omega_corr <- m_omega
    diag(omega_corr) <- rep(1, ncol(m_omega))
    for(i in 1:nrow(corr_set)){
      cs_i <- corr_set[i, ]
      omega_corr[cs_i$PAR1, cs_i$PAR2] <- cs_i$value
      omega_corr[cs_i$PAR2, cs_i$PAR1] <- cs_i$value
    }

    m_omega_full <- m_omega %*% omega_corr %*% m_omega
  }

  ### Reconstruct residual error model matrix
  d_reserr <- par_fin %>% filter(!str_detect(parameter, "_pop|omega_|corr_|beta_"))
  m_reserr <- diag(d_reserr$value, ncol = length(d_reserr$value))
  colnames(m_reserr) <- d_reserr$parameter; rownames(m_reserr) <- d_reserr$parameter

  m_theta_norm <- est_covmat %>% select_if(is.numeric) %>% as.matrix()
  colnames(m_theta_norm) <- est_covmat$X1; rownames(m_theta_norm) <- est_covmat$X1
  m_theta_norm_pop <- m_theta_norm[str_detect(rownames(m_theta_norm), "_pop|beta_"), str_detect(colnames(m_theta_norm), "_pop|beta_")]

  #####----Covariate tables----####
  # Continuous covariates table
  cont_cov_vec <- map_chr(cont_cov_l, function(x) x$NAME)

  cont_cov <- map_dfr(cont_cov_l, function(x) {
    tibble(TR = x$NAME, BTR = x$UTNAME, PAR = list(x$par_vec))
  })

  # Categorical covariates table
  cat_cov_vec <- map_chr(cat_cov_l, function(x) x$NAME)

  cat_unique <- map(cat_cov_vec, function(x){
    ds_catcov[[x]] %>% unique()
  })
  names(cat_unique) <- cat_cov_vec

  # Build cat_cov automatically from cat_cov_l, cat_cov_vec and cat_unique
  cat_cov <- map_dfr(cat_cov_vec, function(x) {
    vals <- as.character(cat_unique[[x]])
    #ref  <- cat_cov_l[[x]]$REF
    ref <- if (is.null(cat_cov_l[[x]]$REF)) {
      as.character(levels(factor(data_fin[[cat_cov_l[[x]]$NAME]]))[1])
    } else {
      cat_cov_l[[x]]$REF
    }
    tibble(
      COV    = x,
      COVVAL = vals,
      CATDES = vals,
      KEY    = as.integer(vals == ref),
      PAR    = list(cat_cov_l[[x]]$par_vec)
    )
  })

  nice_names_cont <- map_dfr(cont_cov_vec, function(x) {
    tibble(COV = x, NICEN = cont_cov_l[[x]][["NICENAME"]])
  })
  nice_names_cat <- map_dfr(cat_cov_vec, function(x) {
    tibble(COV = x, NICEN = cat_cov_l[[x]][["NICENAME"]])
  })
  nice_names <- rbind(nice_names_cont,nice_names_cat)

  ## Derive Covariate distributions
  ### Continuous covariates
  ID <- "ID"
  CCov_vec <- cont_cov$TR
  CatCov_vec <- cat_cov$COV %>% unique()

  data_fin <- data_fin %>%
    mutate(across(all_of(CatCov_vec), as.character))

  ds_cc_tr <- data_fin %>% select(all_of(c(ID)), all_of(cont_cov$TR)) %>% unique() %>%
    gather("TR", "TVALUE", -all_of(c(ID))) %>% left_join(cont_cov, by = "TR")

  ds_cc_btr <- data_fin %>% select(all_of(c(ID)), all_of(cont_cov$BTR)) %>% unique() %>%
    gather("BTR", "NVALUE", -all_of(c(ID)))
  ds_cc_btr_av <- ds_cc_btr %>% group_by(BTR) %>% summarise_at(vars(NVALUE), funSum_av)

  ds_cc <- ds_cc_tr %>%
    left_join(ds_cc_btr, by = c(ID, "BTR")) %>%
    left_join(ds_cc_btr_av, by = "BTR") %>%
    group_by(TR) %>%
    mutate(
      LP = quantile(TVALUE, quantiles[[1]]),
      UP = quantile(TVALUE, quantiles[[2]]))



  # Build per-covariate REF lookup from cont_cov_l
  ref_lookup <- map_dfr(cont_cov_l, function(cov) {
    tibble(TR = cov$NAME, REF_spec = cov$REF)
  })

  # Add REF column: use median(TVALUE) when REF_spec == "median", else use the numeric value
  ds_cc <- ds_cc %>%
    left_join(ref_lookup, by = "TR") %>%
    group_by(TR) %>%
    mutate(REF = if_else(REF_spec == "median", median(TVALUE), as.numeric(REF_spec))) %>%
    select(-REF_spec) %>%
    ungroup()


  #ds_cc_reflab <- select(ds_cc, COV = TR, median) %>% unique() %>% left_join(nice_names, by = "COV") %>% summarise(OUT = str_c(str_c(NICEN, " = ", round(median, 1)), collapse = "\n")) %>% pull(OUT)
  ds_cc_reflab <- select(ds_cc, COV = TR, median, LP, UP) %>% unique() %>% left_join(nice_names, by = "COV") %>% summarise(OUT = str_c(str_c(NICEN, " = ", round(median, 1), " (median)\n", "[",quantiles[[1]]*100,"th percentile: ", round(LP, 1), "; ",quantiles[[2]]*100,"th percentile: ", round(UP, 1), "]"), collapse = "\n")) %>% pull(OUT)

  # ds_cc_reflab <- select(ds_cc, COV = TR, median) %>% unique() %>% left_join(nice_names, by = "COV") %>% summarise(OUT = str_c(str_c(NICEN, " = ", round(median, 1)), collapse = "\n")) %>% pull(OUT) %>% str_c(., "\nResults shown for a dose calculated for the median (reference) patient")

  # Calculate low and upper quantile values for Back transformed covariates
  #Take into account that there are several continuous covariates!!!
  map_ccont <- function(target_ccont, data) {
    idx <- which.min(abs(data$TVALUE - target_ccont))
    data$NVALUE[idx]
  }
  ccont_lab_list <- list()
  for (i in c(1:length(cont_cov_l))){
    ds_cc_i <- ds_cc %>% filter(TR == cont_cov_l[[i]]$NAME)
    q_ccont <- c(unique(ds_cc_i$LP), unique(ds_cc_i$UP), unique(ds_cc_i$REF))

    ccont_labels <- sapply(
      q_ccont,
      map_ccont,
      data = ds_cc_i
    )

    names(ccont_labels) <- c("LP_BTR", "UP_BTR", "REF_BTR")
    ccont_lab_list[[cont_cov_l[[i]]$NAME]] <- ccont_labels
  }

  # Convert ccont_lab_list to a long data frame for joining: COV + KEY -> BCOVVAL
  ccont_lab_df <- map_dfr(names(ccont_lab_list), function(cov_name) {
    vals <- ccont_lab_list[[cov_name]]
    tibble(
      COV     = cov_name,
      KEY     = c("LP",       "UP",       "REF"),
      BCOVVAL = as.numeric(vals[c("LP_BTR", "UP_BTR", "REF_BTR")])
    )
  })

  cc_to_test <- ds_cc %>%
    select(COV = TR, BTR, PAR, mean, median, LP, UP, REF) %>%
    unique() %>%
    gather("KEY", "COVVAL", -COV:-median) %>%
    left_join(ccont_lab_df, by = c("COV", "KEY"))

  ### Categorical covariates
  ds_catc <- data_fin %>% select(all_of(c(ID)), all_of(cat_cov$COV)) %>% unique()

  catc_to_test <- ds_catc %>% select(-all_of(c(ID))) %>% gather("COV", "COVVAL") %>% unique() %>% left_join(cat_cov, by = c("COV", "COVVAL"))

  # Reference values for continuous covariates: log-scale (NAME) and back-transformed (UTNAME)
  cont_cov_ref <- do.call(c, unname(map(cont_cov_l, function(cov) {
    ref_row <- cc_to_test %>% filter(COV == cov$NAME, KEY == "REF")
    ut_name <- if (is.null(cov$UTNAME) || is.na(cov$UTNAME)) cov$NAME else cov$UTNAME
    setNames(
      list(ref_row %>% pull(COVVAL), ref_row %>% pull(BCOVVAL)),
      c(cov$NAME, ut_name)
    )
  })))
  cont_cov_ref <- cont_cov_ref[!duplicated(names(cont_cov_ref))]


  # Reference values for categorical covariates: use REF if defined, else first factor level
  cat_cov_ref <- setNames(
    map(cat_cov_l, function(cov) {
      if (is.null(cov$REF)) {
        as.character(levels(factor(data_fin[[cov$NAME]]))[1])
      } else {
        cov$REF
      }
    }),
    map_chr(cat_cov_l, ~ .x$NAME)
  )

  all_cov_ref <- c(cont_cov_ref, cat_cov_ref)

  # Combine with event table
  ev_t_base <- ev_t_input %>%
    mutate(!!!all_cov_ref)

  ets_cc <- fun_EtCC(ev_t_base, cc_to_test)
  ets_catc <- fun_EtCC(ev_t_base, catc_to_test, T)

   # Sensitivity of model parameters to covariate values
  out_cov_exp_sens = NULL
  t_cov_sens_par = NULL
  out_cov_par_sens = NULL

  out_cov_par_sens <- bind_rows(
    fun_CovSens(ets_cc, covs_i = nice_names$COV, nsim = Nsim,
                mod_fin_i = mod_fin,
                var_exp = var_output, aggr = aggr, nice_names_i = nice_names,
                theta_i = par_fin_tv, thetamat_i = m_theta_norm_pop) %>% mutate(Type = "Continuous"),
    fun_CovSens(ets_catc, covs_i = nice_names$COV, nsim = Nsim, cat = T,
                mod_fin_i = mod_fin,
                var_exp = var_output, aggr = aggr, nice_names_i = nice_names,
                theta_i = par_fin_tv, thetamat_i = m_theta_norm_pop) %>% mutate(Type = "Categorical")
    ) %>% mutate(LAB = fct_inorder(LAB), Type = fct_inorder(Type))

  ## Summary of sensitivity of model parameters to covariate values
  t_cov_sens_par <- out_cov_par_sens %>% mutate_at(vars(mean:P975), signif, 3) %>%
    mutate(`90%CI` = str_c(P05, ", ", P95), KEY = ifelse(is.na(KEY), "Category", KEY), BCOVVAL = as.character(BCOVVAL), BCOVVAL = ifelse(is.na(BCOVVAL), CATDES, BCOVVAL)) %>%
    select(Parameter = VAR, Covariate = NICEN, `Cov. percentile` = KEY, `Cov. value` = BCOVVAL, Mean = mean, `90%CI`)

  # Sensitivity of exposure parameters to covariate values
  out_cov_exp_sens <- bind_rows(
    fun_CovSens(ets_cc, covs_i = nice_names$COV, nsim = Nsim, expos = T, stime_exp = stimes_ss,
                mod_fin_i = mod_fin,
                var_exp = var_output, aggr = aggr, nice_names_i = nice_names,
                theta_i = par_fin_tv, thetamat_i = m_theta_norm_pop) %>% mutate(Type = "Continuous"),
    fun_CovSens(ets_catc, covs_i = nice_names$COV, nsim = Nsim, cat = T, expos = T, stime_exp = stimes_ss,
                mod_fin_i = mod_fin,
                var_exp = var_output, aggr = aggr, nice_names_i = nice_names,
                theta_i = par_fin_tv, thetamat_i = m_theta_norm_pop) %>% mutate(Type = "Categorical")
  ) %>% mutate(LAB = fct_inorder(LAB), Type = fct_inorder(Type))

  covsens_res = list(out_cov_par_sens, t_cov_sens_par, out_cov_exp_sens)

  return(covsens_res)

}

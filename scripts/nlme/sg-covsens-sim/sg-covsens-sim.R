## Author: Alina Melnikova
## First created: 2026-01-28
## Description: Simulation for covariate sensitivity analysis
## Keywords: covariates, simulations, sensitivity

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

fun_CovSens_v2 <- function(et_sim_i, cat = F, expos = F, covs_i = NULL, nsim = 1000, stime_exp = NULL, var_exp = c("Cc_mgL"), m_theta = m_theta_norm_pop_id){
  keep_i <- c("Regimen", "KEY", "COV", "COVVAL")
  if(!cat){
    keep_i <- c(keep_i, "BTR", "BCOVVAL")
  } else {
    keep_i <- c(keep_i, "CATDES")
  }

  sens_i <- et_sim_i %>% map_dfr(function(et_i){
    if(!expos){
      par_i <- unique(et_i$PAR)
      if(is.list(par_i)){ par_i <- par_i[[1]] }
      sim_i <- unique(et_i$id) %>% map_dfr(function(id_i){

        et_ii <- et_i %>% filter(id == id_i) %>% merge(m_theta) %>% mutate(id = 1:nrow(.))

        et_ii <- et_ii %>% mutate(NUM = rep(1:10, 100))
        sim_i <- unique(et_ii)$NUM %>% map_dfr(function(num){

          et_iii <- et_ii %>% filter(NUM == num)
          sim_i <- fun_ForSim(mod_fin, et_iii, 0, vars_i = par_i, theta_i = NULL, thetamat_i = NULL, cov_i = c(covs_i, colnames(m_theta)), nrep = 1, keep = keep_i)$IND %>% mutate(sim.id = id, id = id_i)

        })


      })

    } else {

      sim_i <- unique(et_i$id) %>% map_dfr(function(id_i){

        et_ii <- et_i %>% filter(id == id_i) %>% merge(m_theta) %>% mutate(id = 1:nrow(.))
        et_ii <- unique(et_ii) %>% mutate(NUM = rep(1:10, 100))
        sim_ii <- unique(et_ii$NUM) %>% map_dfr(function(num){
          print(num)
          et_iii <- et_ii %>% filter(NUM == num)

          sim_ii <- fun_ForSim(mod_fin, et_iii, stime_exp, vars_i = var_exp, theta_i = NULL, thetamat_i = NULL, cov_i = c(covs_i, colnames(m_theta)), nrep = 1, keep = keep_i, ncores = parallel::detectCores())$IND %>% mutate(sim.id = id, id = id_i)


          sim_ii_cc <- sim_ii %>%
            group_by_at(vars(all_of(c("sim.id", keep_i, covs_i)))) %>% summarise_at(vars(VALUE), funSum_exp) %>% ungroup()


          sim_ii_mod <- sim_ii_cc %>%
            gather("VAR", "VALUE", -all_of(c("sim.id", keep_i, covs_i)))

          return(sim_ii_mod)

        })



      })

    }

    #log(VALUE) ~ KEY
    #if (mean_type == "mean") {
    sim_i_ref <- sim_i %>% filter(KEY == "REF" | KEY == 1) %>% select(sim.id, VAR, REFVAL = VALUE)
    sim_i_ch <- sim_i %>% filter(KEY != "REF" & KEY != 1) %>% left_join(sim_i_ref, by = c("sim.id", "VAR")) %>%
      mutate(PCH = 100*(VALUE - REFVAL)/REFVAL#,
             #LR = log(VALUE)/log(REFVAL)
      )

    #log

    out_i <- sim_i_ch %>% group_by_at(vars(all_of(c("VAR", keep_i, covs_i)))) %>% summarise_at(vars(PCH), funSum_sim) %>% ungroup()

    #out_i_gmr <- sim_i_ch %>% group_by_at(vars(all_of(c("VAR", keep_i, covs_i)))) %>% summarise_at(vars(LR), funSum_sim) %>% ungroup() %>% select(all_of(c("VAR", keep_i, covs_i)), GMR = mean, GMR05 = P05, GMR95 = P95) #%>%
    # mutate(GMR = exp(GMR), GMR05 = exp(GMR05), GMR95 = exp(GMR95))

    #out_i_gmr2 <- sim_i_ch %>% group_by_at(vars(all_of(c("VAR", keep_i, covs_i)))) %>% summarise_at(vars(VALUE, REFVAL), funSum_sim) %>% ungroup() %>% mutate(GMR2 = VALUE_mean/REFVAL_mean, GMR205 = VALUE_P05/REFVAL_P05, GMR295 = VALUE_P95/REFVAL_P95) %>% select(all_of(c("VAR", keep_i, covs_i)), GMR2, GMR205, GMR295)

    if(!cat){
      out_i_gmr3 <- sim_i %>% mutate(KEY = factor(KEY, levels = c("REF", "LP", "UP")))
    } else {
      out_i_gmr3 <- sim_i %>% mutate(KEY = factor(KEY, levels = c(1, 0)))
    }

    out_i_gmr3_res <- unique(out_i_gmr3$VAR) %>% map_dfr(function(var){

      ds_var <- out_i_gmr3 %>% filter(VAR == var)
      mod_i <- aov(log(VALUE) ~ KEY, data = ds_var)

      out_i <- tibble(
        VALUE = coef(mod_i), SE = vcov(mod_i) %>% diag() %>% sqrt(), LCI = VALUE - SE*qnorm(0.90), UCI = VALUE + SE*qnorm(0.90)
      )%>% mutate_all(exp) %>%
        mutate(KEY = levels(fct_drop(ds_var$KEY)), VAR = var) %>%
        slice(2:nrow(.)) %>% rename(GMR = VALUE, GMR05 = LCI, GMR95 = UCI)

      if (cat) {out_i <- out_i %>% mutate(KEY = as.numeric(KEY))}

      return(out_i)

    })


    #out_i <- left_join(out_i, out_i_gmr) %>% left_join(., out_i_gmr2) %>% left_join(., out_i_gmr3_res)
    out_i <- left_join(out_i, out_i_gmr3_res)

  }) %>% left_join(nice_names, by = "COV")


  if(!cat){
    sens_out <- sens_i %>% select(NICEN, VAR, KEY, mean:GMR95, Regimen, COVVAL:BCOVVAL) %>% mutate(KEY = ifelse(KEY == "LP", "10th perc.", "90th perc."), LAB = str_c(NICEN, "\n", KEY, " (", round(BCOVVAL, 1), ")"))
  } else {
    sens_out <- sens_i %>% select(NICEN, VAR, CATDES, mean:GMR95, Regimen) %>% mutate(LAB = str_c(NICEN, " (", CATDES, ")"))
  }
  sens_out <- sens_out %>% mutate_at(vars(mean:P975), function(p){p/100+1})
  return(sens_out)
}

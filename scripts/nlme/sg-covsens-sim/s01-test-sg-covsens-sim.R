# Functions
funSum_sim <- list(mean   = ~mean(.),
                   median = ~median(.),
                   min    = ~min(.),
                   max    = ~max(.),
                   sd     = ~sd(.),
                   P025   = ~quantile(., 0.025),
                   P05    = ~quantile(., 0.05),
                   P10    = ~quantile(., 0.10),
                   P15   = ~quantile(., 0.15),
                   P25    = ~quantile(., 0.25),
                   P75    = ~quantile(., 0.75),
                   P85   = ~quantile(., 0.85),
                   P90    = ~quantile(., 0.90),
                   P95    = ~quantile(., 0.95),
                   P975   = ~quantile(., 0.975),
                   geom_mean = ~exp(mean(log(.), na.rm = T)),
                   CV     = ~sd(., na.rm = T)/mean(., na.rm = T)*100)

funSum_av <- list(mean   = ~mean(.),
                  median   = ~median(., na.rm = T))

funSum_exp <- list(`Cavg, ug/mL`   = ~mean(.),
                   `Cmin, ug/mL`    = ~min(.),
                   `Cmax, ug/mL`    = ~max(.))

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

fun_ForSim <- function(mod_i, et_i, s_times_i = NULL, vars_i = NULL, theta_i = NULL, omega_i = NULL, sigma_i = NULL, thetamat_i = NULL,
                       thlow = -Inf, nrep = 1,
                       aggr_id = F, aggr_tot = F, keep = NULL, cov_i = NULL, addcov = T, ncores = 1){

  #Test
  # mod_i = mod_fin; et_i = et_i; s_times_i = stime_exp;
  # vars_i = var_exp;
  # theta_i = par_fin_tv;
  # thetamat_i = m_theta_norm_pop; cov_i = covs_i; nrep = nsim;
  # keep = keep_i;
  # omega_i = NULL; sigma_i = NULL; thetamat_i = NULL
  # thlow = -Inf;
  # aggr_id = F;
  # aggr_tot = F;
  # addcov = T; ncores = 1


  et_i_m <- et_i %>% et()
  if(!is.null(s_times_i)){ et_i_m <- et_i_m %>% add.sampling(s_times_i) }
  if(!is.null(cov_i)){
    if(!is.null(s_times_i)){
      et_i_cov <- et_i %>% select(id, all_of(cov_i)) %>% unique()
      et_i_m <- et_i_m %>% as_tibble() %>% left_join(et_i_cov, by = "id")
    } else {
      et_i_cov <- et_i %>% select(all_of(cov_i))
      et_i_m <- et_i_m %>% bind_cols(et_i_cov)
    }
  }
  sim_i <- rxSolve(mod_i, events = et_i_m, params = theta_i, omega = omega_i, sigma = sigma_i, thetaMat = thetamat_i, nStud = nrep, thetaLower = thlow, covsInterpolation = "locf", cores = ncores, addCov = F)

  sim_i_ind <- sim_i %>% as_tibble()
  if(!"id" %in% colnames(sim_i_ind)){ sim_i_ind <- mutate(sim_i_ind, id = 1) }
  if(!"sim.id" %in% colnames(sim_i_ind)){ sim_i_ind <- mutate(sim_i_ind, sim.id = 1) }
  sim_i_ind <- sim_i_ind %>% gather("VAR", "VALUE", -id, -sim.id, -time)
  if(!is.null(vars_i)){sim_i_ind <- sim_i_ind %>% filter( VAR %in% vars_i )}

  sim_i_aggr_id <- NULL; sim_i_aggr_tot <- NULL
  #sim_i_out <- list(IND = sim_i_ind, AGGR_ID = sim_i_aggr_id, AGGR_TOT = sim_i_aggr_tot)
  sim_i_out <- list(IND = sim_i_ind, AGGR_ID = NULL, AGGR_TOT = NULL)

  if(aggr_id){
    sim_i_out$AGGR_ID <- sim_i_ind %>% group_by(id, time, VAR) %>% summarise_at(vars(VALUE), funSum_sim) %>% ungroup()
  }

  if(aggr_tot){
    sim_i_out$AGGR_TOT <- sim_i_ind %>% group_by(time, VAR) %>% summarise_at(vars(VALUE), funSum_sim) %>% ungroup()
  }

  if(!is.null(keep)){
    sim_i_out <- sim_i_out %>% map(function(x){
      if(!is.null(x) & "id" %in% colnames(x)){left_join(x, unique(select(et_i, id, all_of(keep))), by = "id")}
    })
  }

  if(addcov & !is.null(cov_i)){
    sim_i_out <- sim_i_out %>% map(function(x){
      if(!is.null(x) & "id" %in% colnames(x)){left_join(x, unique(select(et_i_cov, id, all_of(cov_i))))}
    })
  }
  return(sim_i_out)
}

fun_CovSens <- function(et_sim_i, cat = F, expos = F, covs_i = NULL, nsim = 100, stime_exp = NULL, var_exp = "Cc") #nsim=1000
  {


  #Test
  et_sim_i = ets_cc
  covs_i = nice_names$COV
  expos = T
  stime_exp = stimes_ss
  cat = F; nsim = 200;
  var_exp = "Cc"


  keep_i <- c("Regimen", "KEY", "COV", "COVVAL")
  if(!cat){
    keep_i <- c(keep_i, "BTR", "BCOVVAL")
  } else {
    keep_i <- c(keep_i, "CATDES")
  }

  sens_i <- et_sim_i %>% map_dfr(function(et_i){
    #Test
    et_i <- et_sim_i[[1]]

      if(!expos){
      par_i <- unique(et_i$PAR)
      if(is.list(par_i)){ par_i <- par_i[[1]] }
      sim_i <- fun_ForSim(mod_fin, et_i, 0, vars_i = par_i, theta_i = par_fin_tv, thetamat_i = m_theta_norm_pop, cov_i = covs_i, nrep = nsim, keep = keep_i)$IND
    } else {
      sim_i <- fun_ForSim(mod_fin, et_i, stime_exp, vars_i = var_exp, theta_i = par_fin_tv,
                          thetamat_i = m_theta_norm_pop, cov_i = covs_i, nrep = nsim,
                          keep = keep_i, ncores = max(1, parallel::detectCores()-2))$IND %>%
        group_by_at(vars(all_of(c("sim.id", keep_i, covs_i)))) %>% summarise_at(vars(VALUE), funSum_exp) %>% ungroup() %>%
        gather("VAR", "VALUE", -all_of(c("sim.id", keep_i, covs_i)))
    }

    sim_i_ref <- sim_i %>% filter(KEY == "REF" | KEY == 1) %>% select(sim.id, VAR, REFVAL = VALUE)
    sim_i_ch <- sim_i %>% filter(KEY != "REF" & KEY != 1) %>% left_join(sim_i_ref, by = c("sim.id", "VAR")) %>%
      mutate(PCH = 100*(VALUE - REFVAL)/REFVAL)

    out_i <- sim_i_ch %>% group_by_at(vars(all_of(c("VAR", keep_i, covs_i)))) %>% summarise_at(vars(PCH), funSum_sim) %>% ungroup()

    # Ensure BCOVVAL is added from et_i if not already present (for continuous covariates)
    if(!cat && "BCOVVAL" %in% names(et_i)){
      if(!"BCOVVAL" %in% names(out_i)){
        bcovval_lookup <- et_i %>% select(KEY, COV, COVVAL, BCOVVAL) %>% unique()
        out_i <- out_i %>% left_join(bcovval_lookup, by = c("KEY", "COV", "COVVAL"))
      }
    }

    })  #%>% left_join(nice_names, by = "COV")

  if(!cat){
    sens_out <- sens_i %>% select(NICEN, VAR, KEY, mean:P975, Regimen, COVVAL:BCOVVAL) %>% mutate(KEY = ifelse(KEY == "LP", "10th perc.", "90th perc."), LAB = str_c(NICEN, "\n", KEY, " (", round(BCOVVAL, 1), ")"))
  } else {
    sens_out <- sens_i %>% select(NICEN, VAR, CATDES, mean:P975, Regimen) %>% mutate(LAB = str_c(NICEN, " (", CATDES, ")"))
  }
  sens_out <- sens_out %>% mutate_at(vars(mean:P975), function(p){p/100+1})
  return(sens_out)
}




# Function parameters
quantiles <- c(0.1, 0.9)
cont_cov_l <- list(
  LG_AGE = list(
    NAME = "LG_AGE",
    UTNAME = "AGE", #(or NA or NULL if UTNAME should be = NAME),
    REF = "median", #“median” or user-defined number,
    NICENAME = "Age, years", #nice name or NULL,
    par_vec = c("CL")
  ),
  LG_WEIGHT = list(
    NAME = "LG_WEIGHT",
    UTNAME = "WEIGHT", #(or NA or NULL if UTNAME should be = NAME),
    REF = "median", #“median” or user-defined number,
    NICENAME = "Weight, kg",
    par_vec = c("Vd") #c("Vd", "CL") #nice name or NULL
  )
)
#cat_cov_vec <- c("Sex" = "SEX", "CYP2C9 genotype" = "CYP2C9")
cat_cov_l <- list(
  SEX = list(
    NAME = "SEX",
    NICENAME = "Sex", #nice name or NULL,
    REF = "0", # NULL or user-defined value,
    par_vec = c("ka")
  ),
  CYP2C9 = list(
    NAME = "CYP2C9",
    NICENAME = "CYP2C9 genotype", #nice name or NULL
    REF = "1.2", # NULL or user-defined value,
    par_vec = c("CL") #c("Vd", "CL")
  )
)


# add check of REF value
aggr = c("min", "max", "mean")

# Load simurg object with covariates
#fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK_cov.RData", package = "SimuRg")
#fpath_i <- system.file("extdata", "simurg_object", "run_4cov_smrg_results.json", package = "SimuRg")
fpath_i <- system.file("extdata", "simurg_object", "run_4cov_smrg_results.json", package = "SimuRg")

obj_data <- read_smrg_obj(fpath_i)
ds_mod <-  obj_data$SDTAB
par_sum <- obj_data$SUMTAB
ds_catcov <- obj_data$CATAB
ds_ccov <- obj_data$COTAB
data_fin <- ds_ccov %>% left_join(ds_catcov, by = "ID")
# Model
mod_fin <- RxODE({
  # Doses in mg
  # Time in hours

  ### Parameter values
  # Typical values
  ka_pop = 0.073;
  Vd_pop = 14.8;
  CL_pop = 0.347;

  # Random effects
  omega_ka = 0;
  omega_Vd = 0;
  omega_CL = 0;

  # Covariate effect
  # Continuous
  beta_CL_LG_AGE = 0.49990114;
  beta_Vd_LG_WEIGHT = 0.60529433;

  # Categorical
  beta_CL_CYP2C9_1_2 = -0.339;
  beta_CL_CYP2C9_1_3 = -0.574;
  beta_CL_CYP2C9_2_2 = -1.079;
  beta_CL_CYP2C9_2_3 = -0.745;
  beta_Cl_CYP2C9_3_3 = -2.13;

  beta_ka_SEX_1 = -0.12198035;

  # Residual error
  Cc_b = 0;

  # Transformations
  ka_tv = exp(ka_pop);
  Vd_tv = exp(Vd_pop);
  CL_tv = exp(CL_pop);

  CL_multiplier = 1.0;  # Default/reference
  ka_multiplier = 1.0;

  if (SEX == 1) {ka_multiplier = exp(beta_ka_SEX_1)}

  if (CYP2C9  == 1.2) {
    CL_multiplier = exp(beta_CL_CYP2C9_1_2);
  } else if (CYP2C9  == 1.3) {
    CL_multiplier = exp(beta_CL_CYP2C9_1_3);
  } else if (CYP2C9  == 2.2) {
    CL_multiplier = exp(beta_CL_CYP2C9_2_2);
  } else if (CYP2C9 == 2.3) {
    CL_multiplier = exp(beta_CL_CYP2C9_2_3);
  } else if (CYP2C9 == 3.3) {
    CL_multiplier = exp(beta_CL_CYP2C9_3_3);
  }

  ka = ka_tv*ka_multiplier*exp(omega_ka);
  Vd = Vd_tv*exp(beta_Vd_LG_WEIGHT * LG_WEIGHT + omega_Vd); #Vd_tv*exp(omega_Vd);
  CL = CL_tv*CL_multiplier*exp(beta_CL_LG_AGE * LG_AGE + omega_CL);


  ### Explicit functions
  Cc = Ac/Vd;

  ### Initial conditions
  Ad(0) = 0;
  Ac(0) = 0;

  ### Differential equations
  d/dt(Ad) = - ka*Ad;
  d/dt(Ac) = ka*Ad - CL*Cc;

  Cc_ResErr = Cc*(1 + Cc_b);
})


### Population parameters
par_fin <- par_sum %>% rename(parameter = PAR, value = VALUE)

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

### Mock Fisher information covariance (same parameter order as par_fin; symmetric pos-def)
pnames <- par_fin$parameter
npar   <- length(pnames)
set.seed(1)
m_cov  <- matrix(0.02, npar, npar)
diag(m_cov) <- 0.05 + runif(npar, 0, 0.05)
m_cov  <- (m_cov + t(m_cov)) / 2
est_covmat <- as_tibble(cbind(X1 = pnames, as.data.frame(m_cov)))
names(est_covmat)[-1] <- pnames

m_theta_norm <- est_covmat %>% select_if(is.numeric) %>% as.matrix()
colnames(m_theta_norm) <- est_covmat$X1; rownames(m_theta_norm) <- est_covmat$X1
m_theta_norm_pop <- m_theta_norm[str_detect(rownames(m_theta_norm), "_pop|beta_"), str_detect(colnames(m_theta_norm), "_pop|beta_")]

###
# cont_cov <- tribble(
#   ~TR,          ~BTR,      ~PAR,
#   # "WTBL",       "OWTBL",    c("Vd", "Vp", "CL", "Vmax", "Q")
#   # "LOG_CSF1",    "CSF1BL",  c("BL_CSF1", "Vd", "CL")
#   "LG_AGE",         "AGE",     c("CL"),
#   "LG_WEIGHT",     "WEIGHT",     c("Vd")
# )

cont_cov <- map_dfr(cont_cov_l, function(x) {
  tibble(TR = x$NAME, BTR = x$UTNAME, PAR = list(x$par_vec))
})

cat_cov_vec <- map_chr(cat_cov_l, function(x) x$NAME)

cat_unique <- map(cat_cov_vec, function(x){
  ds_catcov[[x]] %>% unique()
})
names(cat_unique) <- cat_cov_vec

# Build cat_cov automatically from cat_cov_l, cat_cov_vec and cat_unique
cat_cov <- map_dfr(cat_cov_vec, function(x) {
  vals <- as.character(cat_unique[[x]])
  ref  <- cat_cov_l[[x]]$REF
  tibble(
    COV    = x,
    COVVAL = vals,
    CATDES = vals,
    KEY    = as.integer(vals == ref),
    PAR    = list(cat_cov_l[[x]]$par_vec)
  )
})

# cat_cov <- tribble(
#   ~COV,         ~COVVAL,             ~CATDES,   ~KEY, ~PAR,
#   # "CMSTCAT",    "0",                 "No",      1,    c("Vd"),
#   # "CMSTCAT",    "1",                 "Yes",     0,    c("Vd"),
#   # "POP",        "cGVHD",             "cGVHD",   1,    c("Vd", "Vmax"),
#   # "POP",        "Healthy",           "Healthy", 0,    c("Vd", "Vmax"),
#   # "POP",        "Cancer",            "Cancer",  0,    c("Vd", "Vmax"),
#   # "ADACAT",     "0",                 "No",      1,    c("CL"),
#   # "ADACAT",     "1",                 "Yes",     0,    c("CL")
#   "CYP2C9",        "1.1",              "1.1",     1,    c("CL"),
#   "CYP2C9",        "1.2",              "1.2",     0,    c("CL"),
#   "CYP2C9",        "1.3",              "1.3",     0,    c("CL"),
#   "CYP2C9",        "2.2",              "2.2",     0,    c("CL"),
#   "CYP2C9",        "2.3",              "2.3",     0,    c("CL"),
#   "CYP2C9",        "3.3",              "3.3",     0,    c("CL"),
#   "SEX",           "0",                 "0",      1,    c("ka"),
#   "SEX",           "1",                 "1",      0,    c("ka")
#
# )

nice_names <- tribble(
  ~COV,        ~NICEN,
  # "WTBL",    "WTBL, kg",
  # "LOG_CSF1",  "CSF-1, ng/L",
  # "CMSTCAT",   "Corticosteroids",
  # "POP",   "Population",
  # "ADACAT",    "ADA status",
    "LG_AGE",      cont_cov_l[["LG_AGE"]][["NICENAME"]],
    "LG_WEIGHT",      cont_cov_l[["LG_WEIGHT"]][["NICENAME"]],
  "CYP2C9",      "CYP2C9 genotype",
  "SEX",         "Sex"
)

ss_cycle <- 10
fun_stimes_ss <- function(k){c(
  k*4*7*24 + c(seq(0, 23.5, 0.5), seq(24, 335, 1)),
  k*4*7*24 + 2*7*24 + c(seq(0, 23.5, 0.5), seq(24, 335, 1))
)}
stimes_ss <- fun_stimes_ss(ss_cycle)

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
ds_cc_reflab <- select(ds_cc, COV = TR, median, LP, UP) %>% unique() %>% left_join(nice_names, by = "COV") %>% summarise(OUT = str_c(str_c(NICEN, " = ", round(median, 1), " (median)\n", "[10th percentile: ", round(LP, 1), "; 90th percentile: ", round(UP, 1), "]"), collapse = "\n")) %>% pull(OUT)

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

## Event table

# et_base <- tribble(
#   ~id, ~time, ~evid, ~cmt, ~amt, ~addl, ~ii, ~IGFR, ~POPN,
#   1,   0,     1,     1,    10,   2,     24,  112,   1
# )
#
# et_tvar_cov <- bind_rows(
#   et_base,
#   tibble(id = 1, time = seq(0, 96, 0.5), evid = 0, cmt = 0, amt = 0, addl = 0, ii = 0, IGFR = 112, POPN = 1)
# ) %>% mutate(IGFR = ifelse(time > 24, 30, IGFR))

ev_t_base <- tribble(
  ~id, ~time, ~ii,  ~addl, ~dur,  ~evid, ~Regimen,        ~Dose,
  1,   0,     336,  21,    0.5,   1,     "0.3 mg/kg Q2W", 0.3
) %>% mutate(WEIGHT = unique(select(data_fin, ID, WEIGHT )) %>% pull(WEIGHT ) %>% median(),
             LG_WEIGHT = 1, LG_AGE = 1, AGE = unique(select(data_fin, ID, AGE )) %>% pull(AGE) %>% median(),
             SEX = "0", CYP2C9 = "1.1",
             amt = Dose*WEIGHT)

ets_cc <- fun_EtCC(ev_t_base, cc_to_test)
ets_catc <- fun_EtCC(ev_t_base, catc_to_test, T)

# fn_m_theta_norm_pop_id <- str_c(path_to_save, "m_theta_norm_pop_id.Rdata")
# if(!file.exists(fn_m_theta_norm_pop_id) | overwrite){
  m_theta_norm_pop_id <- MASS::mvrnorm(100, par_fin_tv, m_theta_norm_pop) %>% as.data.frame() #MASS::mvrnorm(1000, par_fin_tv, m_theta_norm_pop)
#   save(m_theta_norm_pop_id, file = fn_m_theta_norm_pop_id)
# }
# load(fn_m_theta_norm_pop_id)


# Sensitivity of model parameters to covariate values
  #fn_covsens_par <- str_c(path_to_save, "sim_covsens_par.Rdata")
  #if(!file.exists(fn_covsens_par) | overwrite){
    out_cov_par_sens <- bind_rows(
      fun_CovSens(ets_cc, covs_i = nice_names$COV) %>% mutate(Type = "Continuous"),
      fun_CovSens(ets_catc, covs_i = nice_names$COV, cat = T) %>% mutate(Type = "Categorical")
    ) %>% mutate(LAB = fct_inorder(LAB), Type = fct_inorder(Type))
   # save(out_cov_par_sens, file = fn_covsens_par)
  #}
  #load(fn_covsens_par)

## Visualization
    p_cov_sens_par <- ggplot(data = out_cov_par_sens, aes(x = LAB, y = mean, ymin = P05, ymax = P95, col = Type)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.8, ymax = 1.25, fill = "firebrick", alpha = 0.2) +
      geom_errorbar(width = 0.2) +
      geom_point(size = 2.5) +
      geom_hline(yintercept = 1, col = "grey25", lwd = 0.8, lty = "dashed") +
      geom_hline(yintercept = c(0.8, 1.25), col = "firebrick", lwd = 0.8, lty = "dotted") +
      scale_color_manual(values = MSDcol[c(1, 3, 4, 5, 6, 7)]) +
      scale_y_continuous(name = "Mean (90% CI) parameter\nchange from reference", breaks = scales::pretty_breaks(7)) +
      labs(x = NULL, caption = ds_cc_reflab) +
      facet_grid(VAR~., scales = "free") +
      coord_flip() +
      theme(legend.position = "top",
            legend.background = element_rect(fill = "white", size = 0.15, linetype = "solid", colour = "black"))

## Summary of sensitivity of model parameters to covariate values
    t_cov_sens_par <- out_cov_par_sens %>% mutate_at(vars(mean:P975), signif, 3) %>%
      mutate(`90%CI` = str_c(P05, ", ", P95), KEY = ifelse(is.na(KEY), "Category", KEY), BCOVVAL = as.character(BCOVVAL), BCOVVAL = ifelse(is.na(BCOVVAL), CATDES, BCOVVAL)) %>%
      select(Parameter = VAR, Covariate = NICEN, `Cov. percentile` = KEY, `Cov. value` = BCOVVAL, Mean = mean, `90%CI`)

# Sensitivity of exposure parameters to covariate values
    out_cov_exp_sens <- bind_rows(
      fun_CovSens(ets_cc, covs_i = nice_names$COV, expos = T, stime_exp = stimes_ss) %>% mutate(Type = "Continuous"),
      fun_CovSens(ets_catc, covs_i = nice_names$COV, cat = T, expos = T, stime_exp = stimes_ss) %>% mutate(Type = "Categorical")
    ) %>% mutate(LAB = fct_inorder(LAB), Type = fct_inorder(Type))

    et_sim_i = ets_cc
    covs_i = nice_names$COV
    expos = T
    stime_exp = stimes_ss
    cat = F; nsim = 200;
    var_exp = "Cc"
    p <- fun_CovSens(ets_cc, covs_i = nice_names$COV, expos = T, stime_exp = stimes_ss) %>% mutate(Type = "Continuous")


## Visualization
    p_cov_sens_exp <- ggplot(data = filter(out_cov_exp_sens, VAR != "Cmin, ug/mL"), aes(x = LAB, y = mean, ymin = P05, ymax = P95, col = Type)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.8, ymax = 1.25, fill = "firebrick", alpha = 0.2) +
      geom_errorbar(width = 0.2) +
      geom_point(size = 2.5) +
      geom_hline(yintercept = 1, col = "grey25", lwd = 0.8, lty = "dashed") +
      geom_hline(yintercept = c(0.8, 1.25), col = "firebrick", lwd = 0.8, lty = "dotted") +
      scale_color_manual(values = MSDcol[c(1, 3, 4, 5, 6, 7)]) +
      scale_y_continuous(name = "Mean (90% CI) parameter\nchange from reference", breaks = scales::pretty_breaks(7)) +
      labs(x = NULL, caption = ds_cc_reflab) +
      facet_grid(VAR~., scales = "free") +
      coord_flip() +
      theme(legend.position = "top",
            legend.background = element_rect(fill = "white", size = 0.15, linetype = "solid", colour = "black"))

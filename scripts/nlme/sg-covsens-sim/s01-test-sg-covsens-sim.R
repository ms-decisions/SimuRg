# Functions
funSum_av <- list(mean   = ~mean(.),
                  median   = ~median(., na.rm = T))

# Function parameters
quantiles <- c(0.1, 0.9)
cont_cov_l <- list(
  NAME = "AGE",
  UTNAME = NA, #(or NA or NULL if UTNAME should be = NAME),
  REF = "median", #“median” or user-defined number,
  NICENAME = "Age, years" #nice name or NULL
)
# add check of REF value
aggr = c("min", "max", "mean")

# Load simurg object with covariates
fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK_cov.RData", package = "SimuRg")

obj_data <- read_smrg_obj(fpath_i)
ds_mod <- obj_data$SDTAB
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
  beta_V_AGE = 0.00472;

  # Categorical
  beta_Cl_CYP2C9_gentyp_1_2 = -0.287;
  beta_Cl_CYP2C9_gentyp_1_3 = -0.566;
  beta_Cl_CYP2C9_gentyp_2_2 = -1.09;
  beta_Cl_CYP2C9_gentyp_2_3 = -0.730;
  beta_Cl_CYP2C9_gentyp_3_3 = -2.45;

  # Residual error
  Cc_b = 0;

  # Transformations
  ka_tv = exp(ka_pop);
  Vd_tv = exp(Vd_pop);
  CL_tv = exp(CL_pop);

  CL_multiplier = 1.0;  # Default/reference

  if (CYP2C9_gentyp == 1.2) {
    CL_multiplier = exp(beta_Cl_CYP2C9_gentyp_1_2);
  } else if (CYP2C9_gentyp == 1.3) {
    CL_multiplier = exp(beta_Cl_CYP2C9_gentyp_1_3);
  } else if (CYP2C9_gentyp == 2.2) {
    CL_multiplier = exp(beta_Cl_CYP2C9_gentyp_2_2);
  } else if (CYP2C9_gentyp == 2.3) {
    CL_multiplier = exp(beta_Cl_CYP2C9_gentyp_2_3);
  } else if (CYP2C9_gentyp == 3.3) {
    CL_multiplier = exp(beta_Cl_CYP2C9_gentyp_3_3);
  }

  ka = ka_tv*exp(omega_ka);
  Vd = exp(beta_V_AGE * AGE + omega_Vd); #Vd_tv*exp(omega_Vd);
  CL = CL_tv*CL_multiplier*exp(omega_CL);


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
# Helper function to create mock sg_fit object





### Population parameters
par_fin <- par_sum %>% rename(parameter = PAR, value = VAL)

# parameter <- "PAR"
# value <- "VALUE"

# par_pop <- par_fin %>% filter(str_detect(.data[[parameter]], "_pop")) %>% select(all_of(c(parameter, value))) %>% deframe()
# par_fin_tv <- par_fin %>% filter(str_detect(.data[[parameter]], "_pop$") | str_detect(.data[[parameter]], "^beta_")) %>% select(all_of(c(parameter, value))) %>%
#   mutate(value = ifelse(str_detect(.data[[parameter]], "_pop$") & !str_detect(.data[[parameter]], "^beta_"), log(value), value)) %>% deframe()
#
# par_fin_tv <- par_fin %>% filter(str_detect(.data[[parameter]], "_pop$") | str_detect(.data[[parameter]], "^beta_")) %>% select(all_of(c(parameter, value))) %>%
# mutate(!!value := ifelse(
#   str_detect(.data[[parameter]], "_pop$") & !str_detect(.data[[parameter]], "^beta_"),
#   ifelse(.data[[value]] > 0, log(.data[[value]]), NA_real_),
#   .data[[value]]
# ))

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
cont_cov <- tribble(
  ~TR,          ~BTR,      ~PAR,
  # "WTBL",       "OWTBL",    c("Vd", "Vp", "CL", "Vmax", "Q")
  # "LOG_CSF1",    "CSF1BL",  c("BL_CSF1", "Vd", "CL")
  "AGE",         "AGE",     c("Vd")
)
cat_cov <- tribble(
  ~COV,         ~COVVAL,             ~CATDES,   ~KEY, ~PAR,
  # "CMSTCAT",    "0",                 "No",      1,    c("Vd"),
  # "CMSTCAT",    "1",                 "Yes",     0,    c("Vd"),
  # "POP",        "cGVHD",             "cGVHD",   1,    c("Vd", "Vmax"),
  # "POP",        "Healthy",           "Healthy", 0,    c("Vd", "Vmax"),
  # "POP",        "Cancer",            "Cancer",  0,    c("Vd", "Vmax"),
  # "ADACAT",     "0",                 "No",      1,    c("CL"),
  # "ADACAT",     "1",                 "Yes",     0,    c("CL")
  "CYP2C9_gentyp", "1.1",              "1.1",     1,    c("CL"),
  "CYP2C9_gentyp", "1.2",              "1.2",     0,    c("CL"),
  "CYP2C9_gentyp", "1.3",              "1.3",     0,    c("CL"),
  "CYP2C9_gentyp", "2.2",              "2.2",     0,    c("CL"),
  "CYP2C9_gentyp", "2.3",              "2.3",     0,    c("CL"),
  "CYP2C9_gentyp", "3.3",              "3.3",     0,    c("CL")

)

nice_names <- tribble(
  ~COV,        ~NICEN,
  # "WTBL",    "WTBL, kg",
  # "LOG_CSF1",  "CSF-1, ng/L",
  # "CMSTCAT",   "Corticosteroids",
  # "POP",   "Population",
  # "ADACAT",    "ADA status",
    "AGE",      cont_cov_l[[1]],
  "CYP2C9_gentyp", "CYP2C9 genotype"
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
    LP = quantile(NVALUE, quantiles[[1]]),
    UP = quantile(NVALUE, quantiles[[2]]),
    REF = case_when(TR == "WTBL" ~ median(TVALUE), T ~ 0))

if (cont_cov_l$REF == "median"){
  ds_cc <- ds_cc %>% mutate(REF = median(TVALUE)) %>% ungroup()
} else {ds_cc <- ds_cc %>% mutate(REF = cont_cov_l$REF)}


#ds_cc_reflab <- select(ds_cc, COV = TR, median) %>% unique() %>% left_join(nice_names, by = "COV") %>% summarise(OUT = str_c(str_c(NICEN, " = ", round(median, 1)), collapse = "\n")) %>% pull(OUT)
ds_cc_reflab <- select(ds_cc, COV = TR, median, LP, UP) %>% unique() %>% left_join(nice_names, by = "COV") %>% summarise(OUT = str_c(str_c(NICEN, " = ", round(median, 1), " (median)\n", "[10th percentile: ", round(LP, 1), "; 90th percentile: ", round(UP, 1), "]"), collapse = "\n")) %>% pull(OUT)

# ds_cc_reflab <- select(ds_cc, COV = TR, median) %>% unique() %>% left_join(nice_names, by = "COV") %>% summarise(OUT = str_c(str_c(NICEN, " = ", round(median, 1)), collapse = "\n")) %>% pull(OUT) %>% str_c(., "\nResults shown for a dose calculated for the median (reference) patient")

# How to rewrite???
cc_to_test <- ds_cc %>% select(COV = TR, BTR, PAR, mean, median, LP, UP, REF) %>% unique() %>%
  gather("KEY", "COVVAL", -COV:-median) %>%
  mutate(BCOVVAL = case_when(COV == "NCMCT" ~ COVVAL + median, COV == "WTBL" ~ COVVAL, T ~ exp(COVVAL)*median))


### Categorical covariates
ds_catc <- data_fin %>% select(all_of(c(ID)), all_of(cat_cov$COV)) %>% unique()

catc_to_test <- ds_catc %>% select(-all_of(c(ID))) %>% gather("COV", "COVVAL") %>% unique() %>% left_join(cat_cov, by = c("COV", "COVVAL"))

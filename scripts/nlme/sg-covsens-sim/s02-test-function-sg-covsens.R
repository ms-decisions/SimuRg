source("scripts/nlme/sg-covsens-sim/sg-covsens-sim.R")

####------ Function parameters ------####
quantiles <- c(0.2, 0.8)
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
    REF = NULL, # NULL or user-defined value,
    par_vec = c("CL") #c("Vd", "CL")
  )
)

ev_t_input <- tribble(
  ~id, ~time, ~ii, ~amt,  ~addl, ~dur,  ~evid, ~Regimen,        ~Dose,
  1,   0,     336,  10,     21,    0.5,   1,     "0.3 mg/kg Q2W", 0.3
)


# add check of REF value
aggr = c("max")

# Load simurg object with covariates
#fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK_cov.RData", package = "SimuRg")
#fpath_i <- system.file("extdata", "simurg_object", "run_4cov_smrg_results.json", package = "SimuRg")
#fpath_i <- system.file("extdata", "simurg_object", "run_4cov_smrg_results.json", package = "SimuRg")
fpath_i <- system.file("scripts", "nlme", "sg-covsens-sim","run_4cov_smrg_results.json", package = "SimuRg")

obj_data <- read_smrg_obj(fpath_i)
#ds_mod <-  obj_data$SDTAB
par_sum <- obj_data$SUMTAB
par_fin <- par_sum %>% rename(parameter = PAR, value = VALUE)#Exemplar table with parameters

ds_catcov <- obj_data$CATAB
ds_ccov <- obj_data$COTAB

data_fin <- ds_ccov %>% left_join(ds_catcov, by = "ID") #Exemplar table with covariate values

### Mock Fisher information covariance (same parameter order as par_fin; symmetric pos-def)
pnames <- par_fin$parameter
npar   <- length(pnames)
set.seed(1)
m_cov  <- matrix(0.02, npar, npar)
diag(m_cov) <- 0.05 + runif(npar, 0, 0.05)
m_cov  <- (m_cov + t(m_cov)) / 2
est_covmat <- as_tibble(cbind(X1 = pnames, as.data.frame(m_cov)))
names(est_covmat)[-1] <- pnames

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

ss_cycle <- 10
fun_stimes_ss <- function(k){c(
  k*4*7*24 + c(seq(0, 23.5, 0.5), seq(24, 335, 1)),
  k*4*7*24 + 2*7*24 + c(seq(0, 23.5, 0.5), seq(24, 335, 1))
)}
stimes_ss <- fun_stimes_ss(ss_cycle)

######
#Test with GFO
output_01 <- sg_covsens_sim(fpath_i, ds_parest=NULL, ds_cov=NULL, mod_fin, stimes_ss, ev_t_input,
                         est_covmat = est_covmat,
                         Nsim=10,
                         cont_cov_l, cat_cov_l,  quantiles = c(0.2, 0.8), aggr = c("max"),
                         var_output = "Cc")

#Test with parameter and covariate datasets
output_02 <- sg_covsens_sim(fpath_i=NULL, ds_parest = par_fin, ds_cov = data_fin,
                            mod_fin, stimes_ss, ev_t_input,
                            est_covmat = est_covmat,
                            Nsim=10,
                            cont_cov_l, cat_cov_l,  quantiles = c(0.2, 0.8), aggr = c("max"),
                            var_output = "Cc")

# Load simurg object with covariates
fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK_cov.RData", package = "SimuRg")

obj_data <- read_smrg_obj(fpath_i)
ds_mod <- obj_data$SDTAB
par_sum <- obj_data$SUMTAB
ds_catcov <- obj_data$CATAB
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



# create_mock_sg_fit <- function() {
#   list(
#     SDTAB = data.frame(
#       ID = rep(1:2, each = 5),
#       TIME = c(0, 1, 2, 4, 8, 0, 1, 2, 4, 8),
#       DV = rnorm(10, mean = 10, sd = 2),
#       MDV = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#     ),
#     EVTAB = data.frame(
#       ID = c(1, 2),
#       time = c(0, 0),
#       amt = c(100, 100),
#       cmt = c(1, 1),
#       evid = c(1, 1)
#     ),
#     SUMTAB = data.frame(
#       PAR = c("ka_pop", "Cl_pop", "V_pop"),
#       VALUE = c(1.5, 2.0, 50.0),
#       TYPE = "Typical values"
#     ),
#     OMEGAMAT = matrix(c(0.1, 0, 0, 0, 0.15, 0, 0, 0, 0.2), nrow = 3, ncol = 3),
#     SIGMAMAT = matrix(c(0.05), nrow = 1, ncol = 1),
#     COTAB = NULL,
#     CATAB = NULL
#   )
# }
# EVTAB = data.frame(
#       ID = c(1, 2),
#       time = c(0, 0),
#       amt = c(100, 100),
#       cmt = c(1, 1),
#       evid = c(1, 1)
#     )
#
# ds_mod <- data.frame(
#   ID = rep(1:2, each = 5),
#   TIME = c(0, 1, 2, 4, 8, 0, 1, 2, 4, 8),
#   DV = rnorm(10, mean = 10, sd = 2),
#   MDV = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
#   WTBL = rep(c(75, 60), each = 5),
#   CRCL = rep(c(95, 72), each = 5),
#   SEX = rep(c("M", "F"), each = 5),
#   SEX_M = rep(c(1, 0), each = 5),
#   RACE = rep(c("White", "Black"), each = 5),
#   RACE_Black = rep(c(0, 1), each = 5),
#   RACE_Other = rep(c(0, 0), each = 5)
# )
#
# ## One-compartment oral model with IIV on ka and V, and covariate effects
# ## Covariates (in event table): WTBL (kg), CRCL (mL/min), SEX_M (0/1), RACE_Black (0/1), RACE_Other (0/1)
# ## Reference: WTBL = 70 kg, CRCL = 90 mL/min, SEX = F, RACE = White
# model <- rxode2::rxode2({
#   # Random effects (passed via omega matrix: omega_ka, omega_V)
#   ka = ka_pop * exp(omega_ka)
#   Cl = Cl_pop * (WTBL / 70)^beta_WTBL_Cl * (CRCL / 90)^beta_CRCL_Cl
#   V  = V_pop * (WTBL / 70)^beta_WTBL_V * exp(omega_V) *
#     (1 + beta_SEX * SEX_M) * (1 + beta_RACE_Black * RACE_Black + beta_RACE_Other * RACE_Other)
#   d/dt(Ad) = -ka * Ad
#   d/dt(Ac) = ka * Ad - Cl/V * Ac
#   Cc = Ac / V
# })
# #path_to_save <- "poppkpd/simulations/"
#
# ### Mock population parameters (covariate sensitivity: WTBL, CRCL, SEX, RACE; IIV on K_a, V_pop)
# par_fin <- tibble(
#   parameter = c(
#     "ka_pop", "Cl_pop", "V_pop",
#     "beta_WTBL_Cl", "beta_WTBL_V", "beta_CRCL_Cl",
#     "beta_SEX", "beta_RACE_Black", "beta_RACE_Other",
#     "omega_ka", "omega_V",
#     "b"
#   ),
#   value = c(
#     1.5, 2.0, 50.0,
#     0.75, 1.0, 0.5,
#     -0.15, 0.08, -0.05,
#     0.25, 0.20,
#     0.08
#   )
# )
par_fin <- par_sum

### Population parameters
parameter <- "PAR"
value <- "VALUE"
par_pop <- par_fin %>% filter(str_detect(.data[[parameter]], "_pop")) %>% select(all_of(c(parameter, value))) %>% deframe()
par_fin_tv <- par_fin %>% filter(str_detect(.data[[parameter]], "_pop$") | str_detect(.data[[parameter]], "^beta_")) %>% select(all_of(c(parameter, value))) %>%
  mutate(value = ifelse(str_detect(.data[[parameter]], "_pop$") & !str_detect(.data[[parameter]], "^beta_"), log(value), value)) %>% deframe()

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
  "WTBL",       "OWTBL",    c("Vd", "Vp", "CL", "Vmax", "Q")
  # "LOG_CSF1",    "CSF1BL",  c("BL_CSF1", "Vd", "CL")
)
cat_cov <- tribble(
  ~COV,         ~COVVAL,             ~CATDES,   ~KEY, ~PAR,
  # "CMSTCAT",    "0",                 "No",      1,    c("Vd"),
  # "CMSTCAT",    "1",                 "Yes",     0,    c("Vd"),
  # "POP",        "cGVHD",             "cGVHD",   1,    c("Vd", "Vmax"),
  # "POP",        "Healthy",           "Healthy", 0,    c("Vd", "Vmax"),
  # "POP",        "Cancer",            "Cancer",  0,    c("Vd", "Vmax"),
  "ADACAT",     "0",                 "No",      1,    c("CL"),
  "ADACAT",     "1",                 "Yes",     0,    c("CL")
)

nice_names <- tribble(
  ~COV,        ~NICEN,
  "WTBL",    "WTBL, kg",
  # "LOG_CSF1",  "CSF-1, ng/L",
  # "CMSTCAT",   "Corticosteroids",
  # "POP",   "Population",
  "ADACAT",    "ADA status",
)

ss_cycle <- 10

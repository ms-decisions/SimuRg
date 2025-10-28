## Author: Victor Sokolov
## First created: 2025-07-15
## Description:
## Keywords:


#####--------------- Load functions and libraries ---------------#####
library(tidyverse)
library(rxode2)
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

funSum_sim <- list(mean   = ~mean(., na.rm = T),
               median = ~median(., na.rm = T),
               min    = ~min(., na.rm = T),
               max    = ~max(., na.rm = T),
               sd     = ~sd(., na.rm = T),
               P025   = ~quantile(., 0.025, na.rm = T),
               P05    = ~quantile(., 0.05, na.rm = T),
               P10    = ~quantile(., 0.10, na.rm = T),
               P15   = ~quantile(., 0.15, na.rm = T),
               P25    = ~quantile(., 0.25, na.rm = T),
               P75    = ~quantile(., 0.75, na.rm = T),
               P85   = ~quantile(., 0.85, na.rm = T),
               P90    = ~quantile(., 0.90, na.rm = T),
               P95    = ~quantile(., 0.95, na.rm = T),
               P975   = ~quantile(., 0.975, na.rm = T),
               geom_mean = ~exp(mean(log(.), na.rm = T)),
               CV     = ~sd(., na.rm = T)/mean(., na.rm = T)*100)

### Simulation function
fun_ForSim <- function(mod_i, et_i, s_times_i = NULL, vars_i = NULL, theta_i = NULL, omega_i = NULL, sigma_i = NULL, thetamat_i = NULL,
                       npop = 1, nsub = 1,
                       aggr_id = F, aggr_tot = F, keep = NULL, cov_i = NULL, addcov = T, ncores = 1){
  et_i_m <- et_i %>% et()
  if(!is.null(s_times_i)){ et_i_m <- et_i_m %>% add.sampling(s_times_i) }
  if(!is.null(cov_i)){
    if(!is.null(s_times_i)){
      et_i_cov <- et_i %>% select(id, all_of(cov_i)) %>% unique()
      et_i_m <- et_i_m %>% as_tibble() %>% left_join(et_i_cov, by = "id")
    } else {
      et_i_cov <- et_i %>% select(id, time, all_of(cov_i))
      et_i_m <- et_i_m %>% bind_cols(select(et_i_cov, -id, -time))
    }
  }
  sim_i <<- rxSolve(mod_i, events = et_i_m, params = theta_i, omega = omega_i, sigma = sigma_i, thetaMat = thetamat_i, nStud = npop, nSub = nsub, covsInterpolation = "locf", cores = ncores)

  sim_i_ind <- sim_i %>% as_tibble()
  if(!"id" %in% colnames(sim_i_ind)){ sim_i_ind <- mutate(sim_i_ind, id = 1) }
  if(!"sim.id" %in% colnames(sim_i_ind)){ sim_i_ind <- mutate(sim_i_ind, sim.id = 1) }
  sim_i_ind <- sim_i_ind %>% gather("VAR", "VALUE", -id, -sim.id, -time)
  if(!is.null(vars_i)){sim_i_ind <- sim_i_ind %>% filter( VAR %in% vars_i )}

  sim_i_aggr_id <- NULL; sim_i_aggr_tot <- NULL
  sim_i_out <- list(IND = sim_i_ind, AGGR_ID = sim_i_aggr_id, AGGR_TOT = sim_i_aggr_tot)

  if(aggr_id){
    sim_i_out$AGGR_ID <- sim_i_ind %>% group_by(id, time, VAR) %>% summarise_at(vars(VALUE), funSum_sim) %>% ungroup()
  }

  if(aggr_tot){
    sim_i_out$AGGR_TOT <- sim_i_ind %>% group_by(time, VAR) %>% summarise_at(vars(VALUE), funSum_sim) %>% ungroup()
  }

  if(!is.null(keep)){
    sim_i_out <- sim_i_out %>% map(function(x){
      if(!is.null(x) & "id" %in% colnames(x)){left_join(x, select(et_i, id, all_of(keep)), by = "id")}
    })
  }

  if(addcov & !is.null(cov_i)){
    sim_i_out <- sim_i_out %>% map(function(x){
      if(!is.null(x) & "id" %in% colnames(x)){left_join(x, et_i_cov)}
    })
  }
  return(sim_i_out)
}


#####--------------- Exemplar model ---------------#####
mod_ex <- RxODE({
  # Doses in mg
  # Time in hours

  ### Parameter values
  # Typical
  POPCL = 5;
  POPVC = 180;
  POPQ = 7;
  POPVP = 52;
  POPKTR = 6;
  FBIOPAR = 1;

  # Covariate coefficients
  POPIGFRCOV = 0.8;
  POPPATCOVCLGR1 = 0.85;
  POPPATCOVCLGR2 = 0.85;
  POPPATCOVCLGR3 = 0.85;
  POPPATCOVCLGR4 = 0.85;

  # Random effects
  PPVCL = 0;
  PPVVC = 0;
  PPVKTR = 0;

  # Residual error
  RUV = 0;

  ### Covariates
  PATCOVCL = 1
  if(POPN == 1){PATCOVCL = POPPATCOVCLGR1}
  if(POPN == 2){PATCOVCL = POPPATCOVCLGR2}
  if(POPN == 3){PATCOVCL = POPPATCOVCLGR3}
  if(POPN == 4){PATCOVCL = POPPATCOVCLGR4}
  IGFRCOV = (IGFR/112)^POPIGFRCOV

  ### Parameters
  CL = POPCL * IGFRCOV * PATCOVCL * exp(PPVCL);
  VC = POPVC * exp(PPVVC);
  Q = POPQ;
  VP = POPVP;
  KTR = POPKTR * exp(PPVKTR);

  ### Explicit functions
  Cc = Ac/VC;                 # nmol/L
  Cp = Ap/VP;                 # nmol/L

  ### Initial conditions
  At1(0) = 0;         # mg
  At2(0) = 0;         # mg
  At3(0) = 0;         # mg
  At4(0) = 0;         # mg
  At5(0) = 0;         # mg
  At6(0) = 0;         # mg
  Ad(0) = 0;          # mg
  Ac(0) = 0;          # mg
  Ap(0) = 0;          # mg

  ### ODEs
  d/dt(At1) = - KTR*At1;
  d/dt(At2) = KTR*At1 - KTR*At2;
  d/dt(At3) = KTR*At2 - KTR*At3;
  d/dt(At4) = KTR*At3 - KTR*At4;
  d/dt(At5) = KTR*At4 - KTR*At5;
  d/dt(At6) = KTR*At5 - KTR*At6;
  d/dt(Ad) = KTR*At6 - KTR*Ad;
  d/dt(Ac) = KTR*Ad - CL*Cc - Q*(Cc - Cp);
  d/dt(Ap) = Q*(Cc - Cp);

  FBIO = FBIOPAR
  f(At1) = FBIO*1000000/505;
  CHECKRUV = RUV;
  Cc_ResErr = Cc + RUV*Cc;
})


#####--------------- Set up event tables ---------------#####
et_base <- tribble(
  ~id, ~time, ~evid, ~cmt, ~amt, ~addl, ~ii, ~IGFR, ~POPN,
  1,   0,     1,     1,    10,   2,     24,  112,   1
)

et_tvar_cov <- bind_rows(
  et_base,
  tibble(id = 1, time = seq(0, 96, 0.5), evid = 0, cmt = 0, amt = 0, addl = 0, ii = 0, IGFR = 112, POPN = 1)
) %>% mutate(IGFR = ifelse(time > 24, 30, IGFR))


#####--------------- Set up parameters ---------------#####
theta <- c(POPCL = 4, POPVC = 100)
omega <- matrix(c(0.2, 0.1, 0, 0.1, 0.2, 0, 0, 0, 0.2), nrow = 3, byrow = T); colnames(omega) <- c("PPVCL", "PPVVC", "PPVKTR")
sigma <- matrix(0.1); colnames(sigma) <- "RUV"
thetamat <- omega; colnames(thetamat) <- c("POPCL", "POPVC", "POPKTR")


#####--------------- Simulations ---------------#####
### Basic
sim1 <- sg_sim(model = mod_ex, et = et_base, stimes = seq(0, 168, 0.1),
               output = "Cc", covs = c("IGFR", "POPN"))
sim2 <- sg_sim(model = mod_ex, et = et_base, stimes = seq(0, 168, 0.1),
               covs = c("IGFR", "POPN"))

### BSV and RUV
sim3 <- sg_sim(model = mod_ex, et = et_base, stimes = seq(0, 168, 0.1),
               output = c("Cc", "Cc_ResErr"), covs = c("IGFR", "POPN"),
               omega = omega, sigma = sigma)
sim4 <- sg_sim(model = mod_ex, et = et_base, stimes = seq(0, 168, 0.1),
               output = c("Cc", "Cc_ResErr"), covs = c("IGFR", "POPN"),
               omega = omega, sigma = sigma,
               nsub = 10, aggr = "ID")

### Uncertainty
sim5 <- sg_sim(model = mod_ex, et = et_base, stimes = seq(0, 168, 0.1),
                   output = c("Cc", "Cc_ResErr"), covs = c("IGFR", "POPN"),
               thetamat = thetamat, npop = 10, aggr = "ID") #
sim6 <- sg_sim(model = mod_ex, et = et_base, stimes = seq(0, 168, 0.1),
               output = c("Cc", "Cc_ResErr"), covs = c("IGFR", "POPN"),
               thetamat = thetamat, nsub = 10, npop = 10, aggr = "ID") #

### BSV, RUV, Uncertainty
sim7 <- sg_sim(model = mod_ex, et = et_base, stimes = seq(0, 168, 0.1),
               output = c("Cc", "Cc_ResErr"), covs = c("IGFR", "POPN"),
                   omega = omega, sigma = sigma, thetamat = thetamat,
               nsub = 10, npop = 10, aggr = "ID")

### Time-varying covariates
sim8 <- sg_sim(model = mod_ex, et = et_tvar_cov, output = c("Cc", "IGFRCOV", "CL"),
               covs = c("IGFR", "POPN"))


#####--------------- Visualization ---------------#####
fun_Plot <- function(sim_i, aggr = F){
  if(!aggr){
    p_i <- ggplot(data = sim_i, aes(x = TIME, y = VALUE, group = interaction(ID, POPN))) +
      geom_line(alpha = 0.9, lwd = 0.8) +
      facet_wrap(~VAR, scales = "free")
  } else {
    p_i <- ggplot(data = sim_i, aes(x = TIME, y = median, ymin = P025, ymax = P975, group = ID)) +
      geom_ribbon(col = NA, alpha = 0.2) +
      geom_line(alpha = 0.9, lwd = 0.8) +
      facet_wrap(~VAR, scales = "free")
  }
  return(p_i)
}

fun_Plot(sim1)
fun_Plot(sim2)

fun_Plot(sim3)
fun_Plot(sim4, T)

# fun_Plot(sim5)
fun_Plot(sim5, T)
fun_Plot(sim6, T)

fun_Plot(sim7, T)

fun_Plot(sim8)


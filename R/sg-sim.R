## Author: Mikhailova Anna
## First created: 2025-09-11
## Description: functions for simulations form the model
## Keywords: SimuRg, simulations

#' Perform simulations from simurg object
#'
#' @inheritParams sg_dummy
#' @returns A dataset with simulation results
#' @details
#'  The dataset returned has the following columns:
#'  \describe{
#'     \item{POPN}{Optional. Population identification number}
#'     \item{ID}{Optional. Simulated subject identification number}
#'     \item{TIME}{Timepoint of the simulations}
#'     \item{VAR}{Names of the outputed variables}
#'     \item{VALUE}{Optional. Simulated value. Returned when no aggregation applied}
#'     \item{mean, median, ...}{Optional. Aggregated statistic of simulated values. Returned when aggregation is applied}
#'     \item{COV1, ...COVn}{Optional. Covariates value. Returned when `addcov` is `TRUE`}
#'     \item{KEEP1, ...KEEPn}{Optional. Columns that were specified in the `keep` argument}
#'  }
#' @examples
#' \dontrun{
#' #####--------------- Set up event tables ---------------#####
#' mod_ex <- RxODE({
#'   # Doses in mg
#'   # Time in hours
#'   ### Parameter values
#'   # Typical
#'   POPCL = 5;
#'   POPVC = 180;
#'   POPQ = 7;
#'   POPVP = 52;
#'   POPKTR = 6;
#'   FBIOPAR = 1;

#'   # Covariate coefficients
#'   POPIGFRCOV = 0.8;
#'   POPPATCOVCLGR1 = 0.85;
#'   POPPATCOVCLGR2 = 0.85;
#'   POPPATCOVCLGR3 = 0.85;
#'   POPPATCOVCLGR4 = 0.85;

#'   # Random effects
#'   PPVCL = 0;
#'   PPVVC = 0;
#'   PPVKTR = 0;

#'   # Residual error
#'   RUV = 0;

#'   ### Covariates
#'   PATCOVCL = 1
#'   if(POPN == 1){PATCOVCL = POPPATCOVCLGR1}
#'   if(POPN == 2){PATCOVCL = POPPATCOVCLGR2}
#'   if(POPN == 3){PATCOVCL = POPPATCOVCLGR3}
#'   if(POPN == 4){PATCOVCL = POPPATCOVCLGR4}
#'   IGFRCOV = (IGFR/112)^POPIGFRCOV

#'   ## Parameters
#'   CL = POPCL * IGFRCOV * PATCOVCL * exp(PPVCL);
#'   VC = POPVC * exp(PPVVC);
#'   Q = POPQ;
#'   VP = POPVP;
#'   KTR = POPKTR * exp(PPVKTR);

#'   ### Explicit functions
#'   Cc = Ac/VC;                 # nmol/L
#'   Cp = Ap/VP;                 # nmol/L

#'   ### Initial conditions
#'   At1(0) = 0;         # mg
#'   At2(0) = 0;         # mg
#'   At3(0) = 0;         # mg
#'   At4(0) = 0;         # mg
#'   At5(0) = 0;         # mg
#'   At6(0) = 0;         # mg
#'   Ad(0) = 0;          # mg
#'   Ac(0) = 0;          # mg
#'   Ap(0) = 0;          # mg

#'   ### ODEs
#'   d/dt(At1) = - KTR*At1;
#'   d/dt(At2) = KTR*At1 - KTR*At2;
#'   d/dt(At3) = KTR*At2 - KTR*At3;
#'   d/dt(At4) = KTR*At3 - KTR*At4;
#'   d/dt(At5) = KTR*At4 - KTR*At5;
#'   d/dt(At6) = KTR*At5 - KTR*At6;
#'   d/dt(Ad) = KTR*At6 - KTR*Ad;
#'   d/dt(Ac) = KTR*Ad - CL*Cc - Q*(Cc - Cp);
#'   d/dt(Ap) = Q*(Cc - Cp);

#'   FBIO = FBIOPAR
#'   f(At1) = FBIO*1000000/505;
#'   CHECKRUV = RUV;
#'   Cc_ResErr = Cc + RUV*Cc;
#' })
#'   et_base <- tribble(
#'     ~id, ~time, ~evid, ~cmt, ~amt, ~addl, ~ii, ~IGFR, ~POPN,
#'     1,   0,     1,     1,    10,   2,     24,  112,   1
#'   )
#'
#'   et_tvar_cov <- bind_rows(
#'     et_base,
#'     tibble(id = 1, time = seq(0, 96, 0.5), evid = 0, cmt = 0, amt = 0, addl = 0, ii = 0, IGFR = 112, POPN = 1)
#'   ) %>% mutate(IGFR = ifelse(time > 24, 30, IGFR))
#'
#'
#'   #####--------------- Set up parameters ---------------#####
#'   omega <- matrix(c(0.2, 0.1, 0, 0.1, 0.2, 0, 0, 0, 0.2), nrow = 3, byrow = T); colnames(omega) <- c("PPVCL", "PPVVC", "PPVKTR")
#'   sigma <- matrix(0.1); colnames(sigma) <- "RUV"
#'   thetamat <- omega; colnames(thetamat) <- c("POPCL", "POPVC", "POPKTR")
#'
#'   et <- et_base
#'   stimes <- seq(0, 168, 0.1)
#'
#'   ### Basic simulations
#'   sim1 <- sg_sim(mod_ex, et_base, stimes, outputs = "Cc", covs = c("IGFR", "POPN"))
#'   sim2 <- sg_sim(mod_ex, et_base, stimes, covs = c("IGFR", "POPN"))
#'
#'   ### BSV and RUV
#'   sim3 <- sg_sim(mod_ex, et_base, stimes, outputs = c("Cc", "Cc_ResErr"),
#'                  covs = c("IGFR", "POPN"), omega = omega, sigma = sigma)
#'   sim4 <- sg_sim(mod_ex, et_base, stimes, outputs = c("Cc", "Cc_ResErr"),
#'                  covs = c("IGFR", "POPN"), omega = omega, sigma = sigma, nsub = 10,
#'                  aggr = "ID")
#'   ### Uncertainty
#'   sim5 <- sg_sim(mod_ex, et_base, stimes, outputs = c("Cc", "Cc_ResErr"),
#'                  covs = c("IGFR", "POPN"), thetamat = thetamat, npop = 10,
#'                  aggr = "ID") #
#'   sim6 <- sg_sim(mod_ex, et_base, stimes, outputs = c("Cc", "Cc_ResErr"),
#'                  covs = c("IGFR", "POPN"), thetamat = thetamat, nsub = 10,
#'                  npop = 10, aggr = "ID")
#'
#'   ### BSV, RUV, Uncertainty
#'   sim7 <- sg_sim(mod_ex, et_base, stimes, outputs = c("Cc", "Cc_ResErr"),
#'                  covs = c("IGFR", "POPN"), omega = omega, sigma = sigma,
#'                  thetamat = thetamat, nsub = 10, npop = 10, aggr = "ID")
#'
#'   ### Time-varying covariates
#'   sim8 <- sg_sim(mod_ex, et_tvar_cov, outputs = c("Cc", "IGFRCOV", "CL"),
#'                  covs = c("IGFR", "POPN"))
#' }
#' @import rxode2
#' @importFrom purrr map
#' @import dplyr
#' @export
sg_sim <- function(model, et, stimes = NULL, outputs = NULL, theta = NULL,
                   omega = NULL, sigma = NULL, sigmaDf = NULL, sigmaLower = -Inf,
                   sigmaUpper = Inf, thetamat = NULL, covs = NULL,
                   npop = 1, nsub = 1, aggr = NULL, addcov = T, keep = NULL,
                   scale = NULL, method = "lsoda", covint = "locf", inits = NULL,
                   ncores = 1, atol = 1e-8, rtol = 1e-6, maxstep = 70000){
  et_i_m <- et %>% et()
  if(!is.null(stimes)){ et_i_m <- et_i_m %>% add.sampling(stimes) }
  if(!is.null(covs)){
    if(!is.null(stimes)){
      et_i_cov <- et %>% select(id, all_of(covs)) %>% unique()
      et_i_m <- et_i_m %>% as_tibble() %>% left_join(et_i_cov, by = "id")
    } else {
      et_i_cov <- et %>% select(id, all_of(covs))
      et_i_m <- et_i_m %>% bind_cols(et_i_cov %>% select(-id))
    }
  }
  sim_i <- rxode2::rxSolve(model, events = et_i_m, params = theta, omega = omega, inits = inits,
                   sigma = sigma, thetaMat = thetamat,nStud = npop, nSub = nsub,
                   covsInterpolation = covint, cores = ncores,
                   simVariability = T, atol = atol, rtol = rtol,
                   maxsteps = maxstep)

  sim_i_ind <- sim_i %>% as_tibble()
  if(!"id" %in% colnames(sim_i_ind)){ sim_i_ind <- mutate(sim_i_ind, id = 1) }
  if(!"sim.id" %in% colnames(sim_i_ind)){ sim_i_ind <- mutate(sim_i_ind, sim.id = 1) }
  sim_i_ind <- sim_i_ind %>% gather("VAR", "VALUE", -any_of(c("id", "sim.id", "time", "POPN")))
  if(!is.null(output)){sim_i_ind <- sim_i_ind %>% filter( VAR %in% output )}
  #
  #   sim_i_aggr_id <- NULL; sim_i_aggr_tot <- NULL
  #   sim_i_out <- list(IND = sim_i_ind, AGGR_ID = sim_i_aggr_id, AGGR_TOT = sim_i_aggr_tot)
  sim_i_out <- sim_i_ind

  if (is.null(aggr)) {
    print("No aggregation applied;")
  } else  if (aggr == "ID"){
    sim_i_out <- sim_i_ind %>% group_by(id, time, VAR) %>% summarise_at(vars(VALUE), funSum_sim) %>% ungroup()
  } else if (aggr == "NPOP"){
    sim_i_out <- sim_i_ind %>% group_by(POPN, time, VAR) %>% summarise_at(vars(VALUE), funSum_sim) %>% ungroup()
  } else if (aggr == "TIME") {
    sim_i_out <- sim_i_ind %>% group_by(time, VAR) %>% summarise_at(vars(VALUE), funSum_sim) %>% ungroup()
  }
  if(!is.null(keep) & !is.null(sim_i_out) & ("id" %in% colnames(sim_i_out))){
    suppressMessages(sim_i_out <- left_join(sim_i_out, unique(select(et, id, all_of(keep))),
                                            by = "id", relationship = "many-to-many"))
  }

  if(addcov & !is.null(covs) & !is.null(sim_i_out) & ("id" %in% colnames(sim_i_out))){
    suppressMessages(sim_i_out <- left_join(sim_i_out, unique(select(et_i_cov, id, all_of(covs))),
                                            relationship = "many-to-many"))
  }

  lookup <- c(TIME = "time", ID = "id")
  req_cols <- c("POPN", "ID", "TIME", "VAR", "VALUE", "mean", "median", "min",
                "max", "sd", "P025", "P05", "P10", "P15", "P25", "P75", "P85",
                "P90", "P95", "P975", "geom_mean", "CV", covs, keep)
  sim_i_out <- sim_i_out %>% rename(any_of(lookup)) %>% select(any_of(req_cols))
  return(sim_i_out)
}

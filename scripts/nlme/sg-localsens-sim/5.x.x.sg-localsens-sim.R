library(tidyverse)
library(rxode2)
library(SimuRg)
theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())
sg_localsens_sim <- function(model, params, stimes, output = NULL,
                             perc = 0.2, lb = NULL, ub = NULL, cov = NULL, et, theta = NULL, n_sim = 10) {

  inits <- rxInits(model)
  params_df <- data.frame(
    Parameter = names(inits),
    Value = inits,
    row.names = NULL
  )

  if(!is.null(lb) & length(lb) != length(params)) stop("The length of lb must match the number of params.")
  if(!is.null(ub) & length(ub) != length(params)) stop("The length of ub must match the number of params.")

  par_bounds <- tibble(
    PAR = as.character(params),
    LB = if(!is.null(lb)) lb else params_df$Value[match(params, params_df$Parameter)] * (1 - perc),
    UB = if(!is.null(ub)) ub else params_df$Value[match(params, params_df$Parameter)] * (1 + perc)
  )

  sens_loc <- par_bounds %>% pmap_dfr(function(PAR, LB, UB) {
    par_range <- seq(LB, UB, length.out = n_sim)

    map_dfr(par_range, function(p) {
      theta_i <- if(is.null(theta)) c() else theta
      theta_i[PAR] <- p

      sim_i <- sg_sim(
        model = model,
        stimes = stimes,
        output = output,
        theta = theta_i,
        covs = cov,
        et = et
      ) %>%
        mutate(
          PARNAME = PAR,
          PARVAL = p,
          PARVAL_NORM = (p - LB) / (UB - LB)
        )

      return(sim_i)
    })
  })

  return(sens_loc)
}




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
et_base <- tribble(
  ~id, ~time, ~evid, ~cmt, ~amt, ~addl, ~ii, ~IGFR, ~POPN,
  1,   0,     1,     1,    10,   2,     24,  112,   1
)
sim__local <- sg_localsens_sim(
     model = mod_ex,
     params = c("POPCL", "POPVC"),
     stimes = seq(0, 168, 0.1),
     output = "Cc",
     perc = 0.9,
     cov = c("IGFR", "POPN"),
     et = et_base
   )

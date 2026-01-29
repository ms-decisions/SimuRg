library(dplyr)
#library(tidyverse)
library(tidyr)
library(ggplot2)

et_base <- tribble(
      ~id, ~time, ~evid, ~cmt, ~amt, ~addl, ~ii, ~IGFR, ~POPN,
      1,   0,     1,     1,    10,   2,     24,  112,   1
    )

    et_tvar_cov <- bind_rows(
      et_base,
      tibble(id = 1, time = seq(0, 96, 0.5), evid = 0, cmt = 0, amt = 0, addl = 0, ii = 0, IGFR = 112, POPN = 1)
    ) %>% mutate(IGFR = ifelse(time > 24, 30, IGFR))

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

    #####--------------- Set up parameters ---------------#####
    omega <- matrix(c(0.2, 0.1, 0, 0.1, 0.2, 0, 0, 0, 0.2), nrow = 3, byrow = T); colnames(omega) <- c("PPVCL", "PPVVC", "PPVKTR")
    sigma <- matrix(0.1); colnames(sigma) <- "RUV"
    thetamat <- omega; colnames(thetamat) <- c("POPCL", "POPVC", "POPKTR")

    et <- et_base
    stimes <- seq(0, 168, 0.1)

    ### Basic simulations
    sim1 <- sg_sim(mod_ex, et_base, stimes, output = "Cc", covs = c("IGFR", "POPN"))
    sim2 <- sg_sim(mod_ex, et_base, stimes, covs = c("IGFR", "POPN"))

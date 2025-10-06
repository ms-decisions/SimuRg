## Author: Mikhailova Anna
## First created: 2025-08-25
## Description: functions for internal data preparation
## Keywords: SimuRg, internal data

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

fpath_i <- file.path("inst", "extdata", "simurg_object", "Warfarin_PK.RData")
obj1 <- get(load(fpath_i))
usethis::use_data(obj1, mod_ex, internal = T, overwrite = T)


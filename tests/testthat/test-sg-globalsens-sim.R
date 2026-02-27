## Author: Ugolkov Yaroslav
## First created: 2026-02-24
## Description: formal testing of sg-globalsens-sim function and its helpers
## Keywords: SimuRg, sg-modcomp, model building, comparison
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
source("inst/extdata/RxODE_model/example_rxode_model.R")
et_base <- tribble(
  ~id, ~time, ~evid, ~cmt, ~amt, ~addl, ~ii, ~IGFR, ~POPN,
  1,   0,     1,     1,    10,   2,     24,  112,   1
)
inits <- rxInits(mod_ex)
par_bounds <- tibble::tibble(
  PAR = c("POPCL", "POPVC"),
  LB  = inits[c("POPCL", "POPVC")] * (1 - 0.9),
  UB  = inits[c("POPCL", "POPVC")] * (1 + 0.9)
)
result <- sg_globalsens_sim(
  method = "PRCC",
  model = mod_ex,
  params = c("POPCL","POPVC"),
  par_bounds = par_bounds,
  n_sim = 100,
  stimes = seq(0, 168, 10),
  output = "Cc",
  cov = c("IGFR", "POPN"),
  et = et_base,
  stat_comp = c("mean")
)
test_that("sg_globalsens_sim PRCC output structure is correct", {


  # List output
  expect_type(result, "list")

  # Required list elements
  expect_true(all(c("result", "summary", "design", "bounds") %in% names(result)))

  # Result table
  expect_s3_class(result$result, "tbl_df")

  # Summary table
  expect_s3_class(result$summary, "tbl_df")

  # Design table
  expect_s3_class(result$design, "tbl_df")

  # Bounds table
  expect_s3_class(result$bounds, "tbl_df")
})
test_that("sg_globalsens_sim PRCC returns correct parameter set", {



  res_tbl <- result$result

  expect_true(all(c("VAR","STAT","PAR","estimate","p.value") %in% names(res_tbl)))

  expect_true(all(sort(unique(res_tbl$PAR)) %in%
                    c("POPCL","POPVC")))
})
result <- sg_globalsens_sim(
  method = "eFAST",
  model = mod_ex,
  params = c("POPCL","POPVC"),
  par_bounds = par_bounds,
  n_sim = 100,
  stimes = seq(0, 168, 10),
  output = "Cc",
  cov = c("IGFR", "POPN"),
  et = et_base,
  stat_comp = c("mean")
)
test_that("sg_globalsens_sim eFAST output structure is correct", {



  # List output
  expect_type(result, "list")

  # Required list elements
  expect_true(all(c("result", "summary", "design", "bounds") %in% names(result)))

  # Result table
  expect_s3_class(result$result, "tbl_df")

  # Summary table
  expect_s3_class(result$summary, "tbl_df")


  # Bounds table
  expect_s3_class(result$bounds, "tbl_df")
})
test_that("sg_globalsens_sim eFAST returns correct parameter set", {



  res_tbl <- result$result

  expect_true(all(c("VAR","STAT","PAR","TYPE","VALUE") %in% names(res_tbl)))

  expect_true(all(sort(unique(res_tbl$PAR)) %in%
                    c("POPCL","POPVC")))
})
test_that("sg_globalsens_sim handles invalid model input", {

  expect_error(
    sg_globalsens_sim(
      method = "PRCC",
      model = "invalid_model",
      params = c("POPCL","POPVC"),
      par_bounds = par_bounds,
      n_sim = 50,
      stimes = seq(0, 168, 10),
      et = et_base
    )
  )
})

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

inits <- rxInits(mod_ex)
params_ext <- c("POPCL", "POPVC", "POPQ")

par_bounds_ext <- tibble::tibble(
  PAR = params_ext,
  LB  = inits[params_ext] * 0.2,
  UB  = inits[params_ext] * 1.8
)

sim_prcc <- sg_globalsens_sim(
  method = "PRCC",
  model = mod_ex,
  params = params_ext,
  par_bounds = par_bounds_ext,
  n_sim = 120,
  stimes = seq(0, 168, 12),
  output = c("Cc", "Cp"),
  cov = c("IGFR", "POPN"),
  et = et_base,
  stat_comp = c("mean", "max", "sd")
)

sim_efast <- sg_globalsens_sim(
  method = "eFAST",
  model = mod_ex,
  params = params_ext,
  par_bounds = par_bounds_ext,
  n_sim = 120,
  stimes = seq(0, 168, 12),
  output = c("Cc", "Cp"),
  cov = c("IGFR", "POPN"),
  et = et_base,
  stat_comp = c("mean", "max", "sd")
)

test_that("PRCC: full combinatorial structure is preserved", {

  df <- sim_prcc$result

  expected_n <- length(unique(df$PAR)) *
    length(unique(df$VAR)) *
    length(unique(df$STAT))

  p <- sg_globalsens_vis(sim_prcc)

  expect_equal(nrow(p$data), expected_n)
})

test_that("eFAST: full combinatorial structure includes TYPE dimension", {

  df <- sim_efast$result

  expected_n <- length(unique(df$PAR)) *
    length(unique(df$VAR)) *
    length(unique(df$STAT)) *
    length(unique(df$TYPE))

  p <- sg_globalsens_vis(sim_efast)

  expect_equal(nrow(p$data), expected_n)
})

test_that("eFAST: each facet has all parameters and TYPE levels", {

  p <- sg_globalsens_vis(sim_efast)

  df <- p$data

  check <- df %>%
    group_by(STAT, VAR) %>%
    summarise(
      n_par = n_distinct(PAR),
      n_type = n_distinct(TYPE),
      .groups = "drop"
    )

  expect_true(all(check$n_par == length(unique(df$PAR))))
  expect_true(all(check$n_type == length(unique(df$TYPE))))
})

test_that("Filtering works correctly in high-dimensional PRCC", {

  p <- sg_globalsens_vis(
    sim_prcc,
    params = c("POPCL", "POPQ"),
    vars = "Cc",
    stats = c("mean", "max")
  )

  df <- p$data

  expect_setequal(unique(df$PAR), c("POPCL", "POPQ"))
  expect_setequal(unique(df$VAR), "Cc")
  expect_setequal(unique(df$STAT), c("mean", "max"))
})

test_that("PRCC heatmap result type is correct", {

  p <- sg_globalsens_vis(
    sim_prcc,
    type = c("heatmap"),
    params = c("POPCL", "POPQ"),
    vars = "Cc",
    stats = c("mean", "max", "sd")
  )

  expect_true(inherits(p, "ggplot"))
  #expect_snapshot(ggplot2::layer_data(p))
})

test_that("PRCC barplot result type is correct", {

  p <- sg_globalsens_vis(
    sim_prcc,
    type = c("bar"),
    params = c("POPCL", "POPQ", "POPVC"),
    vars = c("Cc","Cp"),
    stats = c("max", "sd")
  )

  expect_true(inherits(p, "ggplot"))
  #expect_snapshot(ggplot2::layer_data(p))
})

test_that("eFAST barplot result type is correct", {

  p <- sg_globalsens_vis(
    sim_efast,
    params = c("POPCL", "POPQ", "POPVC"),
    vars = "Cc",
    stats = c("mean", "max", "sd")
  )

  expect_true(inherits(p, "ggplot"))
  #expect_snapshot(ggplot2::layer_data(p))
})

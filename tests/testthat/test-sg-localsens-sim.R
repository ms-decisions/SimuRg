## Author: Ugolkov Yaroslav
## First created: 2025-31-10
## Description: formal testing of sg-localsens-sim function and its helpers
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
et_base <- tribble(
  ~id, ~time, ~evid, ~cmt, ~amt, ~addl, ~ii, ~IGFR, ~POPN,
  1,   0,     1,     1,    10,   2,     24,  112,   1
)
test_that("sg_localsens_sim works and contains all elements", {
  # Просто вызываем функцию и проверяем результат
  result <- sg_localsens_sim(
    model = mod_ex,
    params = c("POPCL", "POPVC"),
    stimes = seq(0, 168, 0.1),
    output = "Cc",
    perc = 0.9,
    cov = c("IGFR", "POPN"),
    et = et_base
  )

  # Проверяем что функция отрабатывает и возвращает tibble
  expect_s3_class(result, "tbl_df")

  # Проверяем что все основные элементы присутствуют
  required_columns <- c("TIME", "VAR", "VALUE", "PARNAME", "PARVAL", "PARVAL_NORM")
  expect_true(all(required_columns %in% names(result)))

  # Проверяем что таблица не пустая
  expect_gt(nrow(result), 0)

  # Проверяем что все параметры из списка params присутствуют
  expect_equal(sort(unique(result$PARNAME)), sort(c("POPCL", "POPVC")))

  # Проверяем что выходная переменная правильная
  expect_equal(unique(result$VAR), "Cc")
})

test_that("sg_localsens_sim fails gracefully", {
  # Ожидаем error для невалидных входных данных
  expect_error(
    sg_localsens_sim(
      model = "invalid_model",
      params = c("POPCL", "POPVC"),
      stimes = seq(0, 168, 0.1),
      et = et_base
    )
  )

  # Ожидаем error для несовпадающих длин lb/ub
  expect_error(
    sg_localsens_sim(
      model = mod_ex,
      params = c("POPCL", "POPVC"),
      stimes = seq(0, 168, 0.1),
      lb = c(1),  # Длина 1 вместо 2
      ub = c(10), # Длина 1 вместо 2
      et = et_base
    )
  )
})

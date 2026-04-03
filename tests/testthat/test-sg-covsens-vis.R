## Author: Melnikova Alina
## First created: 2026-04-03
## Description: formal testing of sg_covsens_vis function and its helpers
## Keywords: SimuRg, sg-covsens-vis

# ---- Shared fixtures -------------------------------------------------------
# Defined at file scope so all test_that() blocks can reuse them.

.cont_cov_l <- list(
  LG_AGE = list(
    NAME     = "LG_AGE",
    UTNAME   = "AGE",
    REF      = "median",
    NICENAME = "Age, years",
    par_vec  = c("CL")
  ),
  LG_WEIGHT = list(
    NAME     = "LG_WEIGHT",
    UTNAME   = "WEIGHT",
    REF      = "median",
    NICENAME = "Weight, kg",
    par_vec  = c("Vd")
  )
)

.cat_cov_l <- list(
  SEX = list(
    NAME     = "SEX",
    NICENAME = "Sex",
    REF      = "0",
    par_vec  = c("ka")
  ),
  CYP2C9 = list(
    NAME     = "CYP2C9",
    NICENAME = "CYP2C9 genotype",
    REF      = NULL,
    par_vec  = c("CL")
  )
)

.ev_t_input <- tribble(
  ~id, ~time, ~ii, ~amt, ~addl, ~dur, ~evid, ~Regimen,         ~Dose,
  1,   0,     336, 10,   21,    0.5,  1,     "0.3 mg/kg Q2W",  0.3
)

.ss_cycle <- 10
.stimes_ss <- c(
  .ss_cycle * 4 * 7 * 24 + c(seq(0, 23.5, 0.5), seq(24, 335, 1)),
  .ss_cycle * 4 * 7 * 24 + 2 * 7 * 24 + c(seq(0, 23.5, 0.5), seq(24, 335, 1))
)

.mod_fin <- RxODE({
  ka_pop = 0.073;
  Vd_pop = 14.8;
  CL_pop = 0.347;

  omega_ka = 0;
  omega_Vd = 0;
  omega_CL = 0;

  beta_CL_LG_AGE     = 0.49990114;
  beta_Vd_LG_WEIGHT  = 0.60529433;

  beta_CL_CYP2C9_1_2 = -0.339;
  beta_CL_CYP2C9_1_3 = -0.574;
  beta_CL_CYP2C9_2_2 = -1.079;
  beta_CL_CYP2C9_2_3 = -0.745;
  beta_CL_CYP2C9_3_3 = -2.13;

  beta_ka_SEX_1 = -0.12198035;

  Cc_b = 0;

  ka_tv = exp(ka_pop);
  Vd_tv = exp(Vd_pop);
  CL_tv = exp(CL_pop);

  CL_multiplier = 1.0;
  ka_multiplier = 1.0;

  if (SEX == "1") {ka_multiplier = exp(beta_ka_SEX_1)}

  if (CYP2C9 == "1") {
    CL_multiplier = exp(beta_CL_CYP2C9_1_2);
  } else if (CYP2C9 == "2") {
    CL_multiplier = exp(beta_CL_CYP2C9_1_3);
  } else if (CYP2C9 == "3") {
    CL_multiplier = exp(beta_CL_CYP2C9_2_2);
  } else if (CYP2C9 == "4") {
    CL_multiplier = exp(beta_CL_CYP2C9_2_3);
  } else if (CYP2C9 == "5") {
    CL_multiplier = exp(beta_CL_CYP2C9_3_3);
  }

  ka = ka_tv * ka_multiplier * exp(omega_ka);
  Vd = Vd_tv * exp(beta_Vd_LG_WEIGHT * LG_WEIGHT + omega_Vd);
  CL = CL_tv * CL_multiplier * exp(beta_CL_LG_AGE * LG_AGE + omega_CL);

  Cc = Ac / Vd;

  Ad(0) = 0;
  Ac(0) = 0;

  d/dt(Ad) = -ka * Ad;
  d/dt(Ac) =  ka * Ad - CL * Cc;

  Cc_ResErr = Cc * (1 + Cc_b);
})

# Mock estimation covariance matrix built from parest (symmetric, positive-definite)
.pnames <- parest$parameter
.npar   <- length(.pnames)
set.seed(1)
.m_cov  <- matrix(0.02, .npar, .npar)
diag(.m_cov) <- 0.05 + runif(.npar, 0, 0.05)
.m_cov  <- (.m_cov + t(.m_cov)) / 2
.est_covmat <- as_tibble(cbind(X1 = .pnames, as.data.frame(.m_cov)))
names(.est_covmat)[-1] <- .pnames

# ---- Cached sg_covsens_sim output (table mode) -----------------------------
.covsens_res <- sg_covsens_sim(
  fpath_i    = NULL,
  ds_parest  = parest,
  ds_covs    = ds_covval,
  model      = .mod_fin,
  stimes     = .stimes_ss,
  et         = .ev_t_input,
  est_covmat = .est_covmat,
  npop       = 10,
  cont_cov_l = .cont_cov_l,
  cat_cov_l  = .cat_cov_l,
  quantiles  = c(0.2, 0.8),
  aggr       = c("max"),
  outputs    = "Cc"
)


# ============================================================
# 1. Return type and successful build
# ============================================================

test_that("sg_covsens_vis returns a ggplot object", {
  p <- sg_covsens_vis(.covsens_res)
  expect_s3_class(p, "gg")
})

test_that("sg_covsens_vis default type builds PARSENS with ggplot_build", {
  p <- sg_covsens_vis(.covsens_res)
  pb <- ggplot2::ggplot_build(p)
  expect_s3_class(pb, "ggplot_built")
})

test_that("sg_covsens_vis type EXPSENS builds with ggplot_build", {
  p <- sg_covsens_vis(.covsens_res, type = "EXPSENS")
  pb <- ggplot2::ggplot_build(p)
  expect_s3_class(pb, "ggplot_built")
})


# ============================================================
# 2. Labels and plotting options
# ============================================================

test_that("sg_covsens_vis applies ylab and caption to plot labels", {
  ylab_txt <- "Custom y label"
  cap_txt <- "Ref: custom caption"
  p <- sg_covsens_vis(.covsens_res, ylab = ylab_txt, caption = cap_txt)
  expect_equal(p$labels$y, ylab_txt)
  expect_equal(p$labels$caption, cap_txt)
})

test_that("sg_covsens_vis accepts alternate ci_quantiles when columns exist", {
  p <- sg_covsens_vis(.covsens_res, ci_quantiles = c("P05", "P95"))
  pb <- ggplot2::ggplot_build(p)
  expect_s3_class(pb, "ggplot_built")
})


# ============================================================
# 3. exclude_vars
# ============================================================

test_that("sg_covsens_vis exclude_vars reduces facet rows to remaining VAR levels", {
  v_excl <- "CL"
  ds_kept <- dplyr::filter(.covsens_res$PARSENS, !.data$VAR %in% v_excl)
  n_var_kept <- dplyr::n_distinct(ds_kept$VAR)

  p <- sg_covsens_vis(.covsens_res, exclude_vars = v_excl)
  pb <- ggplot2::ggplot_build(p)
  layout <- pb$layout$layout

  # One row in layout per facet panel (facet_grid VAR ~ .)
  expect_equal(nrow(layout), n_var_kept)
})


# ============================================================
# 4. Input validation — error handling
# ============================================================

test_that("sg_covsens_vis errors when covsens_res lacks requested type", {
  expect_error(
    sg_covsens_vis(list(PARSENS = .covsens_res$PARSENS), type = "EXPSENS"),
    regexp = "does not contain an element named"
  )
})

test_that("sg_covsens_vis errors when ci_quantiles length is not 2", {
  expect_error(
    sg_covsens_vis(.covsens_res, ci_quantiles = "P025"),
    regexp = "length 2"
  )
})

test_that("sg_covsens_vis errors when ci_limits length is not 2", {
  expect_error(
    sg_covsens_vis(.covsens_res, ci_limits = 0.8),
    regexp = "length 2"
  )
})

test_that("sg_covsens_vis errors when ci_quantiles columns are missing", {
  expect_error(
    sg_covsens_vis(.covsens_res, ci_quantiles = c("P025", "NOT_A_COL")),
    regexp = "Column\\(s\\) not found"
  )
})

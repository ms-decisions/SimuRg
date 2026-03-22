## Author: Melnikova Alina
## First created: 2026-03-19
## Description: formal testing of sg_covsens_sim function and its helpers
## Keywords: SimuRg, sg-covsens-sim

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

# ---- Helper: run once with file mode (gfo4cov) and cache result ------------
.output_01 <- sg_covsens_sim(
  fpath_i    = gfo4cov,
  ds_parest  = NULL,
  ds_covs     = NULL,
  model      = .mod_fin,
  stimes  = .stimes_ss,
  et         = .ev_t_input,
  est_covmat = .est_covmat,
  npop       = 10,
  cont_cov_l = .cont_cov_l,
  cat_cov_l  = .cat_cov_l,
  quantiles  = c(0.2, 0.8),
  aggr       = c("max"),
  outputs    = "Cc"
)

# ---- Helper: run once with table mode (parest + ds_covval) -----------------
.output_02 <- sg_covsens_sim(
  fpath_i    = NULL,
  ds_parest  = parest,
  ds_covs     = ds_covval,
  model      = .mod_fin,
  stimes  = .stimes_ss,
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
# 1. Top-level output structure
# ============================================================

test_that("sg_covsens_sim returns a named list of length 3", {
  expect_type(.output_01, "list")
  expect_length(.output_01, 3)
  expect_named(.output_01, c("PARSENS", "SUMPARSENS", "EXPSENS"))
})

test_that("sg_covsens_sim all three list elements are data frames", {
  expect_true(is.data.frame(.output_01$PARSENS))
  expect_true(is.data.frame(.output_01$SUMPARSENS))
  expect_true(is.data.frame(.output_01$EXPSENS))
  expect_true(nrow(.output_01$PARSENS)   > 0)
  expect_true(nrow(.output_01$SUMPARSENS) > 0)
  expect_true(nrow(.output_01$EXPSENS)   > 0)
})


# ============================================================
# 2. PARSENS — parameter sensitivity
# ============================================================

test_that("PARSENS has required columns", {
  required_cols <- c("NICEN", "VAR", "LAB", "Type",
                     "mean", "median", "min", "max", "sd",
                     "P025", "P05", "P95", "P975")
  expect_true(all(required_cols %in% names(.output_01$PARSENS)))
})

test_that("PARSENS VAR contains expected model parameters", {
  vars <- unique(.output_01$PARSENS$VAR)
  expect_true("CL" %in% vars)
  expect_true("Vd" %in% vars)
  expect_true("ka" %in% vars)
})

test_that("PARSENS NICEN contains expected covariate display names", {
  nicen <- unique(.output_01$PARSENS$NICEN)
  expect_true("Age, years"        %in% nicen)
  expect_true("Weight, kg"        %in% nicen)
  expect_true("Sex"               %in% nicen)
  expect_true("CYP2C9 genotype"   %in% nicen)
})

test_that("PARSENS Type contains both Continuous and Categorical", {
  types <- unique(as.character(.output_01$PARSENS$Type))
  expect_true("Continuous"  %in% types)
  expect_true("Categorical" %in% types)
})

test_that("PARSENS LAB is a factor", {
  expect_true(is.factor(.output_01$PARSENS$LAB))
})

test_that("PARSENS Type is a factor", {
  expect_true(is.factor(.output_01$PARSENS$Type))
})

test_that("PARSENS summary statistics are numeric and non-negative", {
  stat_cols <- c("mean", "median", "P05", "P95", "P025", "P975")
  for (col in stat_cols) {
    vals <- .output_01$PARSENS[[col]]
    expect_true(is.numeric(vals),
                label = paste("PARSENS", col, "is numeric"))
    expect_false(all(is.na(vals)),
                 label = paste("PARSENS", col, "is not all NA"))
  }
})

test_that("PARSENS ratio values are centred around 1 (stored as ratio, not percent)", {
  # Values are transformed: ratio = pct_change/100 + 1, so reference = 1.
  # Non-reference rows should have values reasonably bounded away from zero.
  expect_true(all(.output_01$PARSENS$mean > 0, na.rm = TRUE))
  expect_true(all(.output_01$PARSENS$P05  > 0, na.rm = TRUE))
})

test_that("PARSENS continuous rows have KEY labels with percentile text", {
  cont_rows <- .output_01$PARSENS[.output_01$PARSENS$Type == "Continuous", ]
  key_labels <- unique(as.character(cont_rows$KEY))
  expect_true(any(grepl("perc\\.", key_labels)))
})

test_that("PARSENS BCOVVAL is present and numeric for continuous rows", {
  expect_true("BCOVVAL" %in% names(.output_01$PARSENS))
  cont_rows <- .output_01$PARSENS[.output_01$PARSENS$Type == "Continuous", ]
  expect_true(is.numeric(cont_rows$BCOVVAL) || is.character(cont_rows$BCOVVAL))
})

test_that("PARSENS CATDES is present for categorical rows", {
  expect_true("CATDES" %in% names(.output_01$PARSENS))
  cat_rows <- .output_01$PARSENS[.output_01$PARSENS$Type == "Categorical", ]
  expect_false(all(is.na(cat_rows$CATDES)))
})


# ============================================================
# 3. SUMPARSENS — compact summary table
# ============================================================

test_that("SUMPARSENS has required columns", {
  required_cols <- c("Parameter", "Covariate", "Cov. percentile",
                     "Cov. value", "Mean", "90%CI")
  expect_true(all(required_cols %in% names(.output_01$SUMPARSENS)))
})

test_that("SUMPARSENS Parameter contains expected model parameters", {
  params <- unique(.output_01$SUMPARSENS$Parameter)
  expect_true("CL" %in% params)
  expect_true("Vd" %in% params)
  expect_true("ka" %in% params)
})

test_that("SUMPARSENS 90%CI column is a character with comma-separated values", {
  ci_vals <- .output_01$SUMPARSENS$`90%CI`
  expect_true(is.character(ci_vals))
  expect_true(all(grepl(",", ci_vals)))
})

test_that("SUMPARSENS Mean is numeric", {
  expect_true(is.numeric(.output_01$SUMPARSENS$Mean))
})

test_that("SUMPARSENS row count matches PARSENS non-duplicated covariate-parameter combinations", {
  # Each unique NICEN x VAR x KEY/CATDES combination should yield one SUMPARSENS row
  expect_equal(nrow(.output_01$SUMPARSENS), nrow(.output_01$PARSENS))
})


# ============================================================
# 4. EXPSENS — exposure sensitivity
# ============================================================

test_that("EXPSENS has required columns", {
  required_cols <- c("NICEN", "VAR", "LAB", "Type",
                     "mean", "median", "P05", "P95")
  expect_true(all(required_cols %in% names(.output_01$EXPSENS)))
})

test_that("EXPSENS VAR reflects output and aggr (Cc + max -> Cc_Cmax)", {
  vars <- unique(.output_01$EXPSENS$VAR)
  expect_true("Cc_Cmax" %in% vars)
})

test_that("EXPSENS Type contains both Continuous and Categorical", {
  types <- unique(as.character(.output_01$EXPSENS$Type))
  expect_true("Continuous"  %in% types)
  expect_true("Categorical" %in% types)
})

test_that("EXPSENS summary statistics are numeric and positive", {
  expect_true(is.numeric(.output_01$EXPSENS$mean))
  expect_true(all(.output_01$EXPSENS$mean > 0, na.rm = TRUE))
})

test_that("EXPSENS LAB is a factor", {
  expect_true(is.factor(.output_01$EXPSENS$LAB))
})


# ============================================================
# 5. Table mode (ds_parest + ds_covs) produces consistent output
# ============================================================

test_that("sg_covsens_sim table mode returns a named list of length 3", {
  expect_type(.output_02, "list")
  expect_length(.output_02, 3)
  expect_named(.output_02, c("PARSENS", "SUMPARSENS", "EXPSENS"))
})

test_that("sg_covsens_sim table mode PARSENS has same columns as file mode", {
  expect_equal(sort(names(.output_02$PARSENS)), sort(names(.output_01$PARSENS)))
})

test_that("sg_covsens_sim table mode PARSENS contains same VAR set as file mode", {
  expect_equal(sort(unique(.output_02$PARSENS$VAR)),
               sort(unique(.output_01$PARSENS$VAR)))
})

test_that("sg_covsens_sim table mode EXPSENS VAR matches file mode", {
  expect_equal(sort(unique(.output_02$EXPSENS$VAR)),
               sort(unique(.output_01$EXPSENS$VAR)))
})


# ============================================================
# 6. Input validation — error handling
# ============================================================

test_that("sg_covsens_sim errors when no data source is provided", {
  expect_error(
    sg_covsens_sim(
      fpath_i    = NULL,
      ds_parest  = NULL,
      ds_covs     = NULL,
      model      = .mod_fin,
      stimes  = .stimes_ss,
      et         = .ev_t_input,
      est_covmat = .est_covmat,
      cont_cov_l = .cont_cov_l,
      cat_cov_l  = .cat_cov_l
    ),
    regexp = "No data source provided"
  )
})

test_that("sg_covsens_sim errors when fpath_i and ds_parest are both supplied", {
  expect_error(
    sg_covsens_sim(
      fpath_i    = gfo4cov,
      ds_parest  = parest,
      ds_covs     = NULL,
      model      = .mod_fin,
      stimes  = .stimes_ss,
      et         = .ev_t_input,
      est_covmat = .est_covmat,
      cont_cov_l = .cont_cov_l,
      cat_cov_l  = .cat_cov_l
    ),
    regexp = "mutually exclusive"
  )
})

test_that("sg_covsens_sim errors when ds_parest is provided without ds_covs", {
  expect_error(
    sg_covsens_sim(
      fpath_i    = NULL,
      ds_parest  = parest,
      ds_covs     = NULL,
      model      = .mod_fin,
      stimes  = .stimes_ss,
      et         = .ev_t_input,
      est_covmat = .est_covmat,
      cont_cov_l = .cont_cov_l,
      cat_cov_l  = .cat_cov_l
    ),
    regexp = "ds_covs"
  )
})

test_that("sg_covsens_sim errors when ds_covs is provided without ds_parest", {
  expect_error(
    sg_covsens_sim(
      fpath_i    = NULL,
      ds_parest  = NULL,
      ds_covs     = ds_covval,
      model      = .mod_fin,
      stimes  = .stimes_ss,
      et         = .ev_t_input,
      est_covmat = .est_covmat,
      cont_cov_l = .cont_cov_l,
      cat_cov_l  = .cat_cov_l
    ),
    regexp = "ds_parest"
  )
})


# ============================================================
# 7. aggr argument validation
# ============================================================

test_that("sg_covsens_sim warns on unrecognised aggr value", {
  expect_warning(
    sg_covsens_sim(
      fpath_i    = NULL,
      ds_parest  = parest,
      ds_covs     = ds_covval,
      model      = .mod_fin,
      stimes  = .stimes_ss,
      et         = .ev_t_input,
      est_covmat = .est_covmat,
      npop       = 5,
      cont_cov_l = .cont_cov_l,
      cat_cov_l  = .cat_cov_l,
      quantiles  = c(0.2, 0.8),
      aggr       = c("max", "invalid_metric"),
      outputs    = "Cc"
    ),
    regexp = "Unrecognised aggregation"
  )
})



# ============================================================
# 8. Multiple aggr and multiple outputs
# ============================================================

test_that("sg_covsens_sim EXPSENS VAR reflects multiple aggr metrics", {
  out <- sg_covsens_sim(
    fpath_i    = NULL,
    ds_parest  = parest,
    ds_covs     = ds_covval,
    model      = .mod_fin,
    stimes  = .stimes_ss,
    et         = .ev_t_input,
    est_covmat = .est_covmat,
    npop       = 5,
    cont_cov_l = .cont_cov_l,
    cat_cov_l  = .cat_cov_l,
    quantiles  = c(0.2, 0.8),
    aggr       = c("max", "mean"),
    outputs    = "Cc"
  )
  vars <- unique(out$EXPSENS$VAR)
  expect_true("Cc_Cmax" %in% vars)
  expect_true("Cc_Cavg" %in% vars)
})

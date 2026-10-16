## Author: Mikhailova Anna
## First created: 2025-09-05
## Description: formal testing of sg-parsum function and its helpers
## Keywords: SimuRg, sg-parsum, goodness-of-fit

#source("./scripts/nlme/sg-parsum/sg-parsum.R")

mock_parsum_obj <- function(with_cv = TRUE, with_shrinkage_sd = FALSE, drop = character()) {
  sumtab <- data.frame(
    PAR = c("CL_pop", "omega_CL", "beta_CL_WT"),
    VALUE = c(1.234, 0.3, 0.75),
    TYPE = c("Typical values", "Random effects", "Covariate coefficients"),
    SE = c(0.01234, 0.02, 0.05),
    RSE = c(12.345, 6.67, 6.67),
    LCI = c(0.2304, 0.3, 0.652),
    UCI = c(0.2696, 0.5, 0.848),
    ETAshrinkage_var = c(NA_real_, 12.3456, NA_real_),
    stringsAsFactors = FALSE
  )
  if (with_cv) {
    sumtab$CV <- c(NA_real_, 30, NA_real_)
  }
  if (with_shrinkage_sd) {
    sumtab$ETAshrinkage_sd <- c(NA_real_, 6.789, NA_real_)
    sumtab$EPSshrinkage_sd <- c(NA_real_, NA_real_, 4.321)
  }
  if (length(drop) > 0) {
    sumtab <- sumtab[, setdiff(names(sumtab), drop), drop = FALSE]
  }
  list(
    SUMTAB = sumtab,
    OFV = data.frame(LL = -100, AIC = 210, BIC = 220, BICc = 215)
  )
}

test_that("sg-parsum file load works", {
  skip_on_ci()
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  sum_tab <- sg_parsum(fpath_i)
  expect_true(inherits(sum_tab, "data.frame"))
  expect_snapshot(sum_tab)
  sum_tab <- sg_parsum(obj1)
  expect_true(inherits(sum_tab, "data.frame"))
  expect_snapshot(sum_tab)
})
test_that("sg-parsum file load works dataset with covariates", {
  skip_on_ci()
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK_cov.RData", package = "SimuRg")
  sum_tab <- sg_parsum(fpath_i)
  expect_true(inherits(sum_tab, "data.frame"))
  expect_snapshot(sum_tab)
  sum_tab <- sg_parsum(obj1)
  expect_true(inherits(sum_tab, "data.frame"))
  expect_snapshot(sum_tab)
})
test_that("sg-parsum does not work", {
  expect_error(sg_parsum(data.frame()))
})

test_that("sg-parsum rejects input that is not a fit object", {
  expect_error(
    sg_parsum(1),
    "fpath_i object should be either an sg_fit object, or a path to saved sg_fit object"
  )
})

test_that("sg-parsum errors when SUMTAB is not a data frame", {
  obj <- mock_parsum_obj()
  obj$SUMTAB <- NULL
  expect_error(sg_parsum(obj), "SUMTAB must be a data frame")
})

test_that("sg-parsum reports missing required SUMTAB columns", {
  expect_error(
    sg_parsum(mock_parsum_obj(drop = c("SE", "UCI"))),
    "SUMTAB is missing required columns: SE, UCI"
  )
})

test_that("sg-parsum omits the CV summary when SUMTAB has no CV column", {
  sum_tab <- sg_parsum(mock_parsum_obj(with_cv = FALSE), addOFV = FALSE)

  expect_named(
    sum_tab,
    c("Parameter", "Estimate", "SE", "RSE", "95% CI", "Shrinkage (var), %")
  )
  expect_false("CV % (95% CI)" %in% names(sum_tab))
  expect_equal(sum_tab$Parameter, c("CL", "omega CL", "beta CL WT"))
  expect_equal(sum_tab$Estimate, c(1.23, 0.3, 0.75))
  expect_equal(sum_tab$SE, c(0.0123, 0.02, 0.05))
  expect_equal(sum_tab$RSE, c(12.3, 6.67, 6.67))
  expect_equal(sum_tab$`95% CI`, c("0.23 - 0.27", "0.3 - 0.5", "0.652 - 0.848"))
  expect_equal(sum_tab$`Shrinkage (var), %`, c(NA_real_, 12.3, NA_real_))
})

test_that("sg-parsum adds the CV summary only for random effects", {
  sum_tab <- sg_parsum(mock_parsum_obj(with_cv = TRUE), addOFV = FALSE)

  expect_named(
    sum_tab,
    c("Parameter", "Estimate", "SE", "RSE", "95% CI", "CV % (95% CI)", "Shrinkage (var), %")
  )
  expect_equal(
    sum_tab$`CV % (95% CI)`,
    c(NA_character_, "30 (30.7, 53.3)", NA_character_)
  )
})

test_that("sg-parsum adds ETA and EPS sd shrinkage only when those columns are present", {
  sum_tab <- sg_parsum(mock_parsum_obj(with_cv = FALSE, with_shrinkage_sd = TRUE), addOFV = FALSE)

  expect_named(
    sum_tab,
    c("Parameter", "Estimate", "SE", "RSE", "95% CI", "Shrinkage (var), %", "Shrinkage (sd), %", "EPS shrinkage (sd), %")
  )
  expect_equal(sum_tab$`Shrinkage (sd), %`, c(NA_real_, 6.79, NA_real_))
  expect_equal(sum_tab$`EPS shrinkage (sd), %`, c(NA_real_, NA_real_, 4.32))

  sum_tab_eta_only <- sg_parsum(
    mock_parsum_obj(with_cv = TRUE, with_shrinkage_sd = TRUE, drop = "EPSshrinkage_sd"),
    addOFV = FALSE
  )
  expect_true("Shrinkage (sd), %" %in% names(sum_tab_eta_only))
  expect_false("EPS shrinkage (sd), %" %in% names(sum_tab_eta_only))
  expect_equal(sum_tab_eta_only$`CV % (95% CI)`, c(NA_character_, "30 (30.7, 53.3)", NA_character_))
})


test_that("sg-parsum appends OFV and AIC only when addOFV is TRUE", {
  sum_tab <- sg_parsum(mock_parsum_obj(with_cv = FALSE), addOFV = TRUE)
  expect_equal(sum_tab$OFV, rep(200, 3))
  expect_equal(sum_tab$AIC, rep(210, 3))

  sum_tab_no_ofv <- sg_parsum(mock_parsum_obj(with_cv = FALSE), addOFV = FALSE)
  expect_false(any(c("OFV", "AIC") %in% names(sum_tab_no_ofv)))
})


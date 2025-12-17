## Author: Anatoly Pokladyuk
## First created: 2025-12-17
## Description: formal testing of sg-gof-par-cov function and its helpers
## Keywords: SimuRg, sg-gof-par-cov, goodness-of-fit
test_that("sg_gof_par_cov output is correct", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  cont_cov <- tibble(
    COV = c("AGE", "WTBL"),
    COVNAME = c("Age, years", "Body weight, kg")
  )
  cat_cov <- tibble(
    COV = c("SEX", "VKORC1_gentyp"),
    COVNAME = c("Sex, M/F", "VKORC1 genotype")
  )
  x <- sg_gof_par_cov(fpath_i = fpath_i, ptype = "IndParvsCov", cont_cov = cont_cov)
  expect_true(inherits(x$ipar_vs_contcov, "ggplot"))
  expect_snapshot(ggplot2::layer_data(x$ipar_vs_contcov))
})

test_that("sg_gof_par_cov does not work", {
  expect_error(sg_gof_par_cov(data.frame()))
})

test_that("sg_gof_par_cov file load", {
  fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
  cat_cov <- tibble(
    COV = c("SEX", "VKORC1_gentyp"),
    COVNAME = c("Sex, M/F", "VKORC1 genotype")
  )
  x <- sg_gof_par_cov(fpath_i, ptype = "REvsCov", cat_cov = cat_cov)
  expect_true(inherits(x$re_vs_catcov, "ggplot"))
  expect_snapshot(ggplot2::layer_data(x$re_vs_catcov))
})

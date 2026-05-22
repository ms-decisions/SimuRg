## Author: Alina Melnikova
## First created: 2026-05-19
## Description: Stage 1 helper tests for sg-covsearch
## Keywords: SimuRg, covsearch
source("scripts/nlme/sg-covsearch/sg-covsearch.R")

test_that("get_ofv returns -2*LL and errors clearly when LL is missing", {
  gfo_ok <- list(OFV = data.frame(LL = -123.45))
  expect_equal(get_ofv(gfo_ok), 246.9, tolerance = 1e-10)

  gfo_missing_ll <- list(OFV = data.frame(AIC = 10))
  expect_error(
    get_ofv(gfo_missing_ll),
    regexp = "LL column"
  )
})


test_that("gco_to_theta_tibble preserves order and updates INIT from SUMTAB", {
  gco <- list(
    theta = data.frame(
      NAME = c("ka", "Vd", "CL"),
      TRANS = c("logNormal", "logNormal", "logNormal"),
      INIT = c(0.2, 20, 0.2),
      EST = c(TRUE, TRUE, TRUE),
      stringsAsFactors = FALSE
    )
  )
  gfo <- list(
    SUMTAB = data.frame(
      PAR = c("CL_pop", "ka_pop"),
      VALUE = c(0.31, 0.08),
      stringsAsFactors = FALSE
    )
  )

  theta_updated <- gco_to_theta_tibble(gco, gfo)

  expect_equal(theta_updated$NAME, c("ka", "Vd", "CL"))
  expect_equal(theta_updated$INIT, c(0.08, 20, 0.31), tolerance = 1e-12)
})


test_that("add_covariate builds expected structures for continuous and categorical", {
  covs <- list()

  covs <- add_covariate(
    covs_list = covs,
    param = "CL",
    cov = "AGE",
    type = "continuous"
  )
  expect_equal(covs[[1]]$PAR, "CL")
  expect_equal(covs[[1]]$COVNAME, "AGE")
  expect_equal(covs[[1]]$FUNC, "linear")
  expect_equal(covs[[1]]$TRANS, "median")
  expect_equal(covs[[1]]$INIT, 0.01)
  expect_true(covs[[1]]$EST)

  covs <- add_covariate(
    covs_list = covs,
    param = "ka",
    cov = "SEX",
    type = "categorical",
    cov_ref = "0"
  )
  expect_equal(covs[[2]]$PAR, "ka")
  expect_equal(covs[[2]]$COVNAME, "SEX")
  expect_equal(covs[[2]]$REF, "0")
  expect_equal(covs[[2]]$INIT, 0.01)
  expect_true(covs[[2]]$EST)
  expect_equal(length(covs), 2)
})


test_that("remove_covariate removes only exact match and is idempotent", {
  covs <- list(
    list(PAR = "CL", COVNAME = "AGE", FUNC = "linear", TRANS = "median", INIT = 0.01, EST = TRUE),
    list(PAR = "CL", COVNAME = "WEIGHT", FUNC = "linear", TRANS = "median", INIT = 0.01, EST = TRUE),
    list(PAR = "ka", COVNAME = "SEX", REF = "0", INIT = 0.01, EST = TRUE)
  )

  removed_once <-  remove_covariate(covs, param = "CL", cov = "AGE")
  expect_equal(length(removed_once), 2)
  expect_false(any(vapply(
    removed_once,
    function(x) identical(x$PAR, "CL") && identical(x$COVNAME, "AGE"),
    logical(1)
  )))

  removed_twice <-  remove_covariate(removed_once, param = "CL", cov = "AGE")
  expect_equal(removed_twice, removed_once)
})


test_that("Stage 1 helpers align with fitted project GCO/GFO shape", {
  gco_path <- file.path("scripts","nlme","sg-covsearch", "fitted_project", "wrfrn_pk_base_model_GCO.json")
  gfo_path <- file.path("scripts","nlme","sg-covsearch", "fitted_project", "wrfrn_pk_base_model_GFO.json")
  skip_if_not(file.exists(gco_path), "Fitted GCO JSON fixture not present.")
  skip_if_not(file.exists(gfo_path), "Fitted GFO JSON fixture not present.")

  gco <- jsonlite::fromJSON(gco_path, simplifyDataFrame = TRUE)
  gfo <- jsonlite::fromJSON(gfo_path, simplifyDataFrame = TRUE)

  expect_gt( get_ofv(gfo), 0)

  theta_updated <-  gco_to_theta_tibble(gco, gfo)
  expect_equal(theta_updated$NAME, gco$theta$NAME)
  expect_true(all(is.finite(theta_updated$INIT)))
})

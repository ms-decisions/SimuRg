## Author: Alina Melnikova
## First created: 2026-05-19
## Description: Stage 1 helper tests for sg-covsearch
## Keywords: SimuRg, covsearch

source("scripts/nlme/sg-covsearch/test-fixtures-sg-covsearch.R")
source("scripts/nlme/sg-covsearch/sg-covsearch.R")

test_that("get_ofv returns -2*LL and errors clearly when LL is missing", {
  gfo_ok <- list(OFV = data.frame(LL = -123.45))
  expect_equal( get_ofv(gfo_ok), 246.9, tolerance = 1e-10)

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

  theta_updated <-  gco_to_theta_tibble(gco, gfo)

  expect_equal(theta_updated$NAME, c("ka", "Vd", "CL"))
  expect_equal(theta_updated$INIT, c(0.08, 20, 0.31), tolerance = 1e-12)
})


test_that("add_covariate builds expected structures for continuous and categorical", {
  covs <- list()

  covs <-  add_covariate(
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

  covs <-  add_covariate(
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
  gco_path <- file.path("scripts","nlme","sg-covsearch","fitted_project", "wrfrn_pk_base_model_GCO.json")
  gfo_path <- file.path("scripts","nlme","sg-covsearch","fitted_project", "wrfrn_pk_base_model_GFO.json")
  skip_if_not(file.exists(gco_path), "Fitted GCO JSON fixture not present.")
  skip_if_not(file.exists(gfo_path), "Fitted GFO JSON fixture not present.")

  gco <- jsonlite::fromJSON(gco_path, simplifyDataFrame = TRUE)
  gfo <- jsonlite::fromJSON(gfo_path, simplifyDataFrame = TRUE)

  expect_gt( get_ofv(gfo), 0)

  theta_updated <-  gco_to_theta_tibble(gco, gfo)
  expect_equal(theta_updated$NAME, gco$theta$NAME)
  expect_true(all(is.finite(theta_updated$INIT)))
})


test_that("Stage 2 generates full candidate cross-product when test_pairs is NULL", {
  gco <- list(headers = .covsearch_headers_fixture, theta = .covsearch_theta_fixture)
  gfo <- list(COTAB = .covsearch_cotab_fixture, CATAB = .covsearch_catab_fixture)

  prep <-  stepwise_covariate_selection(
    gco = gco,
    gfo = gfo,
    covariates = c("WT", "AGE", "SEX"),
    parameters = c("CL", "V"),
    test_pairs = NULL
  )

  expect_equal(nrow(prep$candidates), 2 * 3)
  expect_true(all(c("parameter", "covariate", "type", "cov_ref", "df") %in% names(prep$candidates)))
})


test_that("Stage 2 test_pairs restricts candidates and drops invalid rows with warning", {
  gco <- list(headers = .covsearch_headers_fixture, theta = .covsearch_theta_fixture)
  gfo <- list(COTAB = .covsearch_cotab_fixture, CATAB = .covsearch_catab_fixture)

  expect_warning(

    prep <-  stepwise_covariate_selection(
      gco = gco,
      gfo = gfo,
      covariates = c("WT", "AGE", "SEX"),
      parameters = c("CL", "V"),
      test_pairs = .covsearch_test_pairs_fixture
    ),
    regexp = "dropped [0-9]+ invalid test_pairs row"
  )

  expect_equal(nrow(prep$candidates), 3)
  expect_true(all(prep$candidates$parameter %in% c("CL", "V")))
  expect_true(all(prep$candidates$covariate %in% c("WT", "AGE", "SEX")))
})


test_that("Stage 2 derives NA-safe medians for continuous covariates", {
  gco <- list(headers = .covsearch_headers_fixture, theta = .covsearch_theta_fixture)
  gfo <- list(COTAB = .covsearch_cotab_fixture, CATAB = .covsearch_catab_fixture)

  prep <-  stepwise_covariate_selection(
    gco = gco,
    gfo = gfo,
    covariates = c("WT", "AGE"),
    parameters = c("CL")
  )

  expect_equal(as.numeric(prep$cov_ref$WT), .covsearch_expected_refs_fixture$WT)
  expect_equal(as.numeric(prep$cov_ref$AGE), .covsearch_expected_refs_fixture$AGE)
})


test_that("Stage 2 categorical reference applies precedence and tie-break rules", {
  gco <- list(headers = .covsearch_headers_fixture, theta = .covsearch_theta_fixture)
  catab_tie <- data.frame(
    ID = 1:4,
    SEX = c("M", "F", "M", "F"),
    RACE = c("A", "B", "A", "C"),
    stringsAsFactors = FALSE
  )
  gfo <- list(COTAB = .covsearch_cotab_fixture, CATAB = catab_tie)
  tp <- data.frame(
    parameter = c("CL", "V"),
    covariate = c("SEX", "RACE"),
    type = c("cat", "cat"),
    reference = c("M", NA),
    center = c(NA, NA),
    stringsAsFactors = FALSE
  )

  prep <-  stepwise_covariate_selection(
    gco = gco,
    gfo = gfo,
    covariates = c("SEX", "RACE"),
    parameters = c("CL", "V"),
    test_pairs = tp
  )

  sex_row <- prep$candidates[prep$candidates$covariate == "SEX", , drop = FALSE]
  race_row <- prep$candidates[prep$candidates$covariate == "RACE", , drop = FALSE]

  expect_equal(sex_row$cov_ref[[1]], "M")
  expect_equal(race_row$cov_ref[[1]], "A")
})


test_that("Stage 2 df map uses continuous=1 and categorical k-1 excluding NA", {
  gco <- list(headers = .covsearch_headers_fixture, theta = .covsearch_theta_fixture)
  catab_with_na <- data.frame(
    ID = 1:7,
    SEX = c("M", "F", "M", "M", NA, "F", "M"),
    RACE = c("A", "B", "A", NA, "C", "B", "A"),
    stringsAsFactors = FALSE
  )
  gfo <- list(COTAB = .covsearch_cotab_fixture, CATAB = catab_with_na)

  prep <-  stepwise_covariate_selection(
    gco = gco,
    gfo = gfo,
    covariates = c("WT", "SEX", "RACE"),
    parameters = c("CL")
  )

  expect_equal(as.numeric(prep$df_map$WT), 1)
  expect_equal(as.numeric(prep$df_map$SEX), 1)
  expect_equal(as.numeric(prep$df_map$RACE), 2)
})


test_that("Stage 2 validates p_forward and p_backward defaults/range", {
  gco <- list(headers = .covsearch_headers_fixture, theta = .covsearch_theta_fixture)
  gfo <- list(COTAB = .covsearch_cotab_fixture, CATAB = .covsearch_catab_fixture)

  prep <-  stepwise_covariate_selection(gco = gco, gfo = gfo)
  expect_equal(prep$p_forward, 0.05)
  expect_equal(prep$p_backward, 0.01)

  expect_error(
     stepwise_covariate_selection(gco = gco, gfo = gfo, p_forward = 1),
    regexp = "p_forward"
  )
  expect_error(
     stepwise_covariate_selection(gco = gco, gfo = gfo, p_backward = 0),
    regexp = "p_backward"
  )
})

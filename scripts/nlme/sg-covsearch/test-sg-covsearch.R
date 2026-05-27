## Author: OpenAI Assistant
## First created: 2026-05-19
## Description: Stage 1 helper tests for sg-covsearch
## Keywords: SimuRg, covsearch

if (file.exists("scripts/nlme/sg-covsearch/test-fixtures-sg-covsearch.R")) {
  source("scripts/nlme/sg-covsearch/test-fixtures-sg-covsearch.R")
} else {
  source("test-fixtures-sg-covsearch.R")
}
if (file.exists("scripts/nlme/sg-covsearch/sg-covsearch.R")) {
  source("scripts/nlme/sg-covsearch/sg-covsearch.R")
} else {
  source("sg-covsearch.R")
}

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
  gco_path <- file.path("scripts", "nlme", "sg-covsearch", "fitted_project", "wrfrn_pk_base_model_GCO.json")
  gfo_path <- file.path("scripts", "nlme", "sg-covsearch", "fitted_project", "wrfrn_pk_base_model_GFO.json")
  if (!file.exists(gco_path)) {
    gco_path <- file.path("fitted_project", "wrfrn_pk_base_model_GCO.json")
  }
  if (!file.exists(gfo_path)) {
    gfo_path <- file.path("fitted_project", "wrfrn_pk_base_model_GFO.json")
  }
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


.stage3_mock_fit <- function(ofv_map, sumtab_map = list()) {
  force(ofv_map)
  force(sumtab_map)
  function(model, data, headers, theta, ruv, re, occ, covs, project_name,
           task_opt = NULL, opt_name = "Monolix", fit = TRUE,
           path_to_save_output = NULL, path_to_fitter = NULL) {
    ofv <- ofv_map[[project_name]]
    if (is.null(ofv)) {
      stop(sprintf("Missing mock OFV for project %s", project_name))
    }
    sumtab <- sumtab_map[[project_name]]
    if (is.null(sumtab)) {
      sumtab <- data.frame(PAR = character(0), VALUE = numeric(0), stringsAsFactors = FALSE)
    }
    list(
      OFV = data.frame(LL = -as.numeric(ofv) / 2),
      SUMTAB = sumtab,
      COTAB = data.frame(dummy = 1),
      CATAB = data.frame(dummy = 1)
    )
  }
}


.stage3_gco_fixture <- function() {
  list(
    model = "mock_model.txt",
    data = "mock_data.csv",
    headers = .covsearch_headers_fixture,
    theta = data.frame(
      NAME = c("ka", "CL", "V"),
      TRANS = c("logNormal", "logNormal", "logNormal"),
      INIT = c(0.2, 0.2, 20),
      EST = c(TRUE, TRUE, TRUE),
      stringsAsFactors = FALSE
    ),
    ruv = list(YNAME = "y1", DVID = 1, TRANS = "normal", PRED = "Cc", ERR = "proportional",
               INIT = 0.1, EST = TRUE, BLQM = NULL),
    re = list(init = diag(c(1, 1, 1)), est = diag(c(TRUE, TRUE, TRUE))),
    occ = list(init = matrix(0, nrow = 3, ncol = 3), est = matrix(NA, nrow = 3, ncol = 3)),
    covs = list()
  )
}


.stage3_gfo_fixture <- function(base_ofv = 100) {
  list(
    OFV = data.frame(LL = -base_ofv / 2),
    SUMTAB = data.frame(
      PAR = c("ka_pop", "CL_pop", "V_pop"),
      VALUE = c(0.2, 0.2, 20),
      stringsAsFactors = FALSE
    ),
    COTAB = data.frame(
      ID = 1:6,
      WT = c(70, 80, 85, 90, 75, 95),
      AGE = c(30, 40, 35, 50, 45, 55),
      stringsAsFactors = FALSE
    ),
    CATAB = data.frame(
      ID = 1:6,
      SEX = c("0", "1", "0", "0", "1", "0"),
      stringsAsFactors = FALSE
    )
  )
}


test_that("Stage 3 mocked-fit selects best significant candidate each step", {
  gco <- .stage3_gco_fixture()
  gfo <- .stage3_gfo_fixture(base_ofv = 100)
  ofv_map <- list(
    fw_s01_CL_WT = 90,   # delta 10 (best step 1)
    fw_s01_CL_AGE = 98,  # delta 2
    fw_s01_CL_SEX = 97,  # delta 3
    fw_s01_V_WT = 95,    # delta 5
    fw_s01_V_AGE = 99,   # delta 1
    fw_s01_V_SEX = 96,   # delta 4
    fw_s02_CL_AGE = 88,  # vs current 90 -> delta 2
    fw_s02_CL_SEX = 89,  # delta 1
    fw_s02_V_WT = 86,    # delta 4 (best step 2)
    fw_s02_V_AGE = 90,   # delta 0
    fw_s02_V_SEX = 89,   # delta 1
    fw_s03_CL_AGE = 85,  # vs current 86 -> delta 1
    fw_s03_CL_SEX = 85,  # delta 1
    fw_s03_V_AGE = 85,   # delta 1
    fw_s03_V_SEX = 85    # delta 1
  )

  res <- stepwise_covariate_selection(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT", "AGE", "SEX"),
    parameters = c("CL", "V"),
    p_forward = 0.05,
    fit_function = .stage3_mock_fit(ofv_map)
  )

  expect_true(res$forward_ran)
  expect_equal(res$included$parameter, c("CL", "V"))
  expect_equal(res$included$covariate, c("WT", "WT"))
  expect_equal(res$final_ofv, 86)
})


test_that("Stage 3 stops when no candidate passes qchisq threshold", {
  gco <- .stage3_gco_fixture()
  gfo <- .stage3_gfo_fixture(base_ofv = 100)
  ofv_map <- list(
    fw_s01_CL_WT = 98.5,   # delta 1.5
    fw_s01_CL_AGE = 98.9,  # delta 1.1
    fw_s01_CL_SEX = 99,    # delta 1.0
    fw_s01_V_WT = 98.6,    # delta 1.4
    fw_s01_V_AGE = 99.2,   # delta 0.8
    fw_s01_V_SEX = 98.8    # delta 1.2
  )

  res <- stepwise_covariate_selection(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT", "AGE", "SEX"),
    parameters = c("CL", "V"),
    p_forward = 0.05,
    fit_function = .stage3_mock_fit(ofv_map)
  )

  expect_equal(nrow(res$included), 0)
  expect_true(all(res$forward_history$accepted == FALSE))
})


test_that("Stage 3 removes accepted candidate from remaining set", {
  gco <- .stage3_gco_fixture()
  gfo <- .stage3_gfo_fixture(base_ofv = 100)
  ofv_map <- list(
    fw_s01_CL_WT = 90,
    fw_s01_CL_AGE = 98,
    fw_s01_CL_SEX = 97,
    fw_s01_V_WT = 95,
    fw_s01_V_AGE = 99,
    fw_s01_V_SEX = 96,
    fw_s02_CL_AGE = 85,
    fw_s02_CL_SEX = 85,
    fw_s02_V_WT = 85,
    fw_s02_V_AGE = 85,
    fw_s02_V_SEX = 85
  )

  res <- stepwise_covariate_selection(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT", "AGE", "SEX"),
    parameters = c("CL", "V"),
    p_forward = 0.05,
    fit_function = .stage3_mock_fit(ofv_map)
  )

  tested_cl_wt <- res$forward_history[
    res$forward_history$parameter == "CL" & res$forward_history$covariate == "WT",
    ,
    drop = FALSE
  ]
  expect_equal(nrow(tested_cl_wt), 1)
  expect_true(tested_cl_wt$accepted[[1]])
})


test_that("Stage 3 history stores threshold, df, and decision", {
  gco <- .stage3_gco_fixture()
  gfo <- .stage3_gfo_fixture(base_ofv = 100)
  ofv_map <- list(
    fw_s01_CL_WT = 96,
    fw_s01_CL_AGE = 99,
    fw_s01_CL_SEX = 99,
    fw_s01_V_WT = 99,
    fw_s01_V_AGE = 99,
    fw_s01_V_SEX = 99
  )

  res <- stepwise_covariate_selection(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT", "AGE", "SEX"),
    parameters = c("CL", "V"),
    p_forward = 0.05,
    fit_function = .stage3_mock_fit(ofv_map)
  )

  expect_true(all(c("threshold", "df", "decision") %in% names(res$forward_history)))
  expect_true(any(res$forward_history$decision == "accepted"))
  expect_true(all(is.finite(res$forward_history$threshold)))
})


test_that("Stage 3 updates theta INIT only from accepted step fits", {
  gco <- .stage3_gco_fixture()
  gfo <- .stage3_gfo_fixture(base_ofv = 100)
  ofv_map <- list(
    fw_s01_CL_WT = 90,
    fw_s01_CL_AGE = 99,
    fw_s01_CL_SEX = 99,
    fw_s01_V_WT = 99,
    fw_s01_V_AGE = 99,
    fw_s01_V_SEX = 99,
    fw_s02_CL_AGE = 89,
    fw_s02_CL_SEX = 89,
    fw_s02_V_WT = 89,
    fw_s02_V_AGE = 89,
    fw_s02_V_SEX = 89
  )
  sumtab_map <- list(
    fw_s01_CL_WT = data.frame(PAR = "ka_pop", VALUE = 0.5, stringsAsFactors = FALSE),
    fw_s01_CL_AGE = data.frame(PAR = "ka_pop", VALUE = 8, stringsAsFactors = FALSE),
    fw_s01_CL_SEX = data.frame(PAR = "ka_pop", VALUE = 9, stringsAsFactors = FALSE),
    fw_s01_V_WT = data.frame(PAR = "ka_pop", VALUE = 10, stringsAsFactors = FALSE),
    fw_s01_V_AGE = data.frame(PAR = "ka_pop", VALUE = 11, stringsAsFactors = FALSE),
    fw_s01_V_SEX = data.frame(PAR = "ka_pop", VALUE = 12, stringsAsFactors = FALSE)
  )

  res <- stepwise_covariate_selection(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT", "AGE", "SEX"),
    parameters = c("CL", "V"),
    p_forward = 0.05,
    update_theta_init = TRUE,
    fit_function = .stage3_mock_fit(ofv_map, sumtab_map)
  )

  ka_row <- res$final_theta[res$final_theta$NAME == "ka", , drop = FALSE]
  expect_equal(ka_row$INIT[[1]], 0.5, tolerance = 1e-12)
})




test_that("Stage 4 removes least important removable term first", {
  gco <- .stage3_gco_fixture()
  gfo <- .stage3_gfo_fixture(base_ofv = 100)
  ofv_map <- list(
    fw_s01_CL_WT = 90,
    fw_s01_CL_AGE = 95,
    fw_s02_CL_AGE = 84,
    bw_s01_CL_WT_removed = 85,
    bw_s01_CL_AGE_removed = 89,
    bw_s02_CL_AGE_removed = 93
  )

  res <- stepwise_covariate_selection(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT", "AGE"),
    parameters = c("CL"),
    p_forward = 0.05,
    p_backward = 0.01,
    run_backward = TRUE,
    fit_function = .stage4_mock_fit(ofv_map)
  )

  expect_equal(res$forward$selected$covariate, c("WT", "AGE"))
  expect_equal(res$backward$removed$covariate[[1]], "WT")
  expect_equal(res$backward$retained$covariate, "AGE")
  expect_equal(res$final_covariates, res$backward$retained)
})


test_that("Stage 4 stops without removals when all retained terms are significant", {
  gco <- .stage3_gco_fixture()
  gfo <- .stage3_gfo_fixture(base_ofv = 100)
  ofv_map <- list(
    fw_s01_CL_WT = 90,
    fw_s01_CL_AGE = 98,
    fw_s02_CL_AGE = 89,
    bw_s01_CL_WT_removed = 108
  )

  res <- stepwise_covariate_selection(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT", "AGE"),
    parameters = c("CL"),
    p_forward = 0.05,
    p_backward = 0.01,
    run_backward = TRUE,
    fit_function = .stage4_mock_fit(ofv_map)
  )

  expect_equal(nrow(res$backward$removed), 0)
  expect_equal(res$backward$retained, res$forward$selected)
  expect_equal(res$metadata$backward_steps, 0)
})


test_that("Stage 4 skips backward loop cleanly when forward selected set is empty", {
  gco <- .stage3_gco_fixture()
  gfo <- .stage3_gfo_fixture(base_ofv = 100)
  ofv_map <- list(
    fw_s01_CL_WT = 99.5,
    fw_s01_CL_AGE = 99.2
  )

  res <- stepwise_covariate_selection(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT", "AGE"),
    parameters = c("CL"),
    p_forward = 0.05,
    p_backward = 0.01,
    run_backward = TRUE,
    fit_function = .stage4_mock_fit(ofv_map)
  )

  expect_equal(nrow(res$forward$selected), 0)
  expect_equal(nrow(res$backward$history), 0)
  expect_equal(nrow(res$backward$retained), 0)
  expect_equal(res$final_covariates, res$backward$retained)
})


test_that("Stage 4 records fit failures and continues testing removable candidates", {
  gco <- .stage3_gco_fixture()
  gfo <- .stage3_gfo_fixture(base_ofv = 100)
  ofv_map <- list(
    fw_s01_CL_WT = 90,
    fw_s01_CL_AGE = 95,
    fw_s02_CL_AGE = 84,
    bw_s01_CL_AGE_removed = 86,
    bw_s02_CL_WT_removed = 92
  )

  res <- stepwise_covariate_selection(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT", "AGE"),
    parameters = c("CL"),
    p_forward = 0.05,
    p_backward = 0.01,
    run_backward = TRUE,
    fit_function = .stage4_mock_fit(
      ofv_map = ofv_map,
      fail_projects = c("bw_s01_CL_WT_removed")
    )
  )

  expect_true(any(res$backward$history$status == "fit_failed"))
  expect_true(any(grepl("Mock fit failure", res$backward$history$message, fixed = TRUE)))
  expect_true(any(res$backward$history$removed))
  expect_equal(res$backward$removed$covariate[[1]], "AGE")
})

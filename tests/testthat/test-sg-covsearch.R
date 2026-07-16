## Author: Alina Melnikova
## First created: 2026-05-19
## Description: Unit and integration tests for sg-covsearch
## Keywords: SimuRg, covsearch

#Mock objects
.covsearch_headers_fixture <- list(
  list(name = "ID", use = "identifier", type = NULL),
  list(name = "TIME", use = "time", type = NULL),
  list(name = "DV", use = "observation", type = "continuous"),
  list(name = "WT", use = "covariate", type = "continuous"),
  list(name = "AGE", use = "covariate", type = "continuous"),
  list(name = "SEX", use = "covariate", type = "categorical"),
  list(name = "RACE", use = "covariate", type = "categorical")
)

# Minimal theta fixture with at least two parameters.
.covsearch_theta_fixture <- data.frame(
  NAME = c("CL", "V"),
  TRANS = c("logNormal", "logNormal"),
  INIT = c(0.2, 20),
  EST = c(TRUE, TRUE),
  stringsAsFactors = FALSE
)

# Tiny COTAB with NA included for NA-safe median checks.
.covsearch_cotab_fixture <- data.frame(
  ID = 1:6,
  WT = c(70, 80, NA, 90, 85, 75),
  AGE = c(30, 40, NA, 50, 40, 60),
  stringsAsFactors = FALSE
)

# Tiny CATAB with class imbalance:
# - SEX: 2 levels with imbalance (expect df = 1)
# - RACE: 3 levels with imbalance (expect df = 2)
.covsearch_catab_fixture <- data.frame(
  ID = 1:8,
  SEX = c("0", "0", "0", "0", "1", "0", "0", "0"),
  RACE = c("A", "A", "B", "A", "C", "A", "B", "A"),
  stringsAsFactors = FALSE
)

# Candidate pairs fixture with valid + invalid entries.
.covsearch_test_pairs_fixture <- data.frame(
  parameter = c("CL", "CL", "V", "Q", "CL"),
  covariate = c("WT", "SEX", "AGE", "AGE", "HEIGHT"),
  type = c("cont", "cat", "cont", "cont", "cont"),
  reference = c(NA, "M", NA, NA, NA),
  center = c("median", NA, "median", "median", "median"),
  stringsAsFactors = FALSE
)

# Optional expected values for Stage 2 assertions.
.covsearch_expected_refs_fixture <- list(
  WT = 80,    # median of c(70,80,90,85,75)
  AGE = 40,   # median of c(30,40,50,40,60)
  SEX = "0",  # mode from imbalanced SEX; explicit can override to "M"
  RACE = "A"  # mode from imbalanced RACE
)

.covsearch_expected_df_fixture <- c(
  WT = 1,    # continuous
  AGE = 1,   # continuous
  SEX = 1,   # 2 levels -> k-1
  RACE = 2   # 3 levels -> k-1
)

.stage4_mock_fit <- function(ofv_map, fail_projects = character(0), sumtab_map = list()) {
  force(ofv_map)
  force(fail_projects)
  force(sumtab_map)
  function(model, data, headers, theta, ruv, re, occ, covs, project_name,
           task_opt = NULL, opt_name = "Monolix", fit = TRUE,
           path_to_save_output = NULL, path_to_fitter = NULL) {
    if (project_name %in% fail_projects) {
      stop(sprintf("Mock fit failure for project %s", project_name))
    }
    ofv <- ofv_map[[project_name]]
    if (is.null(ofv)) {
      stop(sprintf("Missing mock OFV for project %s", project_name))
    }
    sumtab <- sumtab_map[[project_name]]
    if (is.null(sumtab)) {
      sumtab <- data.frame(PAR = character(0), VALUE = numeric(0), stringsAsFactors = FALSE)
    }
    list(
      GFO = list(
        OFV = data.frame(LL = -as.numeric(ofv) / 2),
        SUMTAB = sumtab,
        COTAB = data.frame(dummy = 1),
        CATAB = data.frame(dummy = 1)
      )
    )
  }
}

#Functions

.ensure_mock_files <- function() {
  model_path <- file.path(tempdir(), "mock_model.txt")
  data_path <- file.path(tempdir(), "mock_data.csv")
  if (!file.exists(model_path)) {
    writeLines(c(
      "[LONGITUDINAL]",
      "input = {ka, CL, V}",
      "PK:"
    ), model_path)
  }
  if (!file.exists(data_path)) {
    writeLines("ID,TIME,DV,WT,AGE,SEX,RACE\n1,0,0,70,30,0,A", data_path)
  }
  list(model = normalizePath(model_path, winslash = "/", mustWork = TRUE),
       data = normalizePath(data_path, winslash = "/", mustWork = TRUE))
}

.stage2_gco_fixture <- function(project_name = NULL) {
  mock_paths <- .ensure_mock_files()
  out <- list(
    model = mock_paths$model,
    data = mock_paths$data,
    headers = .covsearch_headers_fixture,
    theta = .covsearch_theta_fixture,
    ruv = list(dummy = TRUE),
    re = list(dummy = TRUE),
    occ = list(dummy = TRUE),
    covs = list()
  )
  if (!is.null(project_name)) {
    out$project_name <- project_name
  }
  out
}

.stage2_gfo_fixture <- function(ll = -100, catab = .covsearch_catab_fixture) {
  list(
    OFV = data.frame(LL = ll),
    SUMTAB = data.frame(PAR = character(0), VALUE = numeric(0), stringsAsFactors = FALSE),
    COTAB = .covsearch_cotab_fixture,
    CATAB = catab
  )
}

.stage_mock_fit_project <- function(ofv_map, sumtab_map = list(), fail_projects = character(0)) {
  force(ofv_map)
  force(sumtab_map)
  force(fail_projects)
  function(model, data, headers, theta, ruv, re, occ, covs, project_name,
           task_opt = NULL, opt_name = "Monolix", fit = TRUE,
           path_to_save_output = NULL, path_to_fitter = NULL) {
    if (project_name %in% fail_projects) {
      stop(sprintf("Mock fit failure for project %s", project_name))
    }
    ofv <- ofv_map[[project_name]]
    if (is.null(ofv)) {
      stop(sprintf("Missing mock OFV for project %s", project_name))
    }
    sumtab <- sumtab_map[[project_name]]
    if (is.null(sumtab)) {
      sumtab <- data.frame(PAR = character(0), VALUE = numeric(0), stringsAsFactors = FALSE)
    }
    list(
      GFO = list(
        OFV = data.frame(LL = as.numeric(ofv)),
        SUMTAB = sumtab,
        COTAB = data.frame(dummy = 1),
        CATAB = data.frame(dummy = 1)
      )
    )
  }
}

.forward_history_files_exist <- function(out_dir) {
  file.exists(file.path(out_dir, "forward_history.csv")) &&
    file.exists(file.path(out_dir, "backward_history.csv"))
}

test_that("get_ofv returns LL and errors clearly when LL is missing", {
  gfo_ok <- list(OFV = data.frame(LL = -123.45))
  expect_equal(get_ofv(gfo_ok), -123.45, tolerance = 1e-10)

  expect_error(get_ofv(list()), regexp = "gfo\\$OFV is missing")
  expect_error(get_ofv(list(OFV = data.frame(AIC = 10))), regexp = "LL column")
  expect_error(get_ofv(list(OFV = data.frame(LL = NA_real_))), regexp = "missing or NA")
  expect_error(get_ofv(list(OFV = data.frame(LL = numeric(0)))), regexp = "missing or NA")

  gfo_many <- list(OFV = data.frame(LL = c(-10, -20)))
  expect_equal(get_ofv(gfo_many), -10)
})

test_that(".covsearch_sumtab_df normalizes data.frame/matrix/list records", {
  sumtab_df <- data.frame(PAR = c("CL_pop", "V_pop"), VALUE = c("0.2", "10"))
  gfo_df <- list(SUMTAB = sumtab_df)
  out_df <- .covsearch_sumtab_df(gfo_df)
  expect_equal(out_df$PAR, c("CL_pop", "V_pop"))
  expect_equal(out_df$VALUE, c(0.2, 10))

  sumtab_mx <- as.matrix(data.frame(PAR = "ka_pop", VALUE = "0.4", stringsAsFactors = FALSE))
  out_mx <- .covsearch_sumtab_df(list(SUMTAB = sumtab_mx))
  expect_equal(out_mx$PAR[[1]], "ka_pop")
  expect_equal(out_mx$VALUE[[1]], 0.4)

  sumtab_list <- list(
    list(PAR = "CL_pop", VALUE = "0.3"),
    list(PAR = "Q_pop")
  )
  out_list <- .covsearch_sumtab_df(list(SUMTAB = sumtab_list))
  expect_equal(out_list$PAR, c("CL_pop", "Q_pop"))
  expect_equal(out_list$VALUE[[1]], 0.3)
  expect_true(is.na(out_list$VALUE[[2]]))
})

test_that("gco_to_theta_tibble preserves theta order and does not update INIT", {
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

  theta_out <- gco_to_theta_tibble(gco, gfo, update_theta_init = TRUE)
  expect_equal(theta_out$NAME, c("ka", "Vd", "CL"))
  expect_equal(theta_out$INIT, c(0.2, 20, 0.2))
})

test_that("gco_to_theta_tibble handles list input and validates required fields", {
  gco_list <- list(
    theta = list(
      list(NAME = "CL", INIT = 0.2, TRANS = "logNormal"),
      list(NAME = "V", INIT = 20, TRANS = "logNormal")
    )
  )
  out <- gco_to_theta_tibble(gco_list, list())
  expect_equal(out$NAME, c("CL", "V"))

  expect_error(gco_to_theta_tibble(list(), list()), regexp = "gco\\$theta is missing")
  expect_error(
    gco_to_theta_tibble(list(theta = data.frame(NAME = "CL")), list()),
    regexp = "missing required columns: INIT"
  )
  expect_error(
    gco_to_theta_tibble(list(theta = data.frame(INIT = 1)), list()),
    regexp = "missing required columns: NAME"
  )
})

test_that("add_covariate builds structures and validates input", {
  covs <- add_covariate(
    covs_list = list(),
    param = "CL",
    cov = "AGE",
    type = "continuous"
  )
  expect_equal(covs[[1]]$PAR, "CL")
  expect_equal(covs[[1]]$COVNAME, "AGE")
  expect_equal(covs[[1]]$FUNC, "linear")
  expect_equal(covs[[1]]$TRANS, "median")
  expect_equal(covs[[1]]$INIT, 0)
  expect_true(covs[[1]]$EST)

  covs <- add_covariate(
    covs_list = covs,
    param = "ka",
    cov = "SEX",
    type = "categorical",
    cov_ref = "0"
  )
  expect_equal(covs[[2]]$REF, "0")
  expect_equal(covs[[2]]$INIT, 0)

  expect_error(add_covariate("bad", "CL", "WT", "continuous"), regexp = "must be a list")
  expect_error(add_covariate(list(), "", "WT", "continuous"), regexp = "param")
  expect_error(add_covariate(list(), "CL", "", "continuous"), regexp = "cov")
  expect_error(add_covariate(list(), "CL", "WT", "unknown"), regexp = "continuous|categorical")
  expect_error(add_covariate(list(), "CL", "SEX", "categorical"), regexp = "cov_ref")
})

test_that("remove_covariate handles null, malformed records, and duplicate pairs", {
  expect_equal(remove_covariate(NULL, "CL", "WT"), list())
  expect_error(remove_covariate("bad", "CL", "WT"), regexp = "must be a list")

  covs <- list(
    list(PAR = "CL", COVNAME = "WT"),
    list(PAR = "CL", COVNAME = "WT"),
    list(PAR = "CL", COVNAME = "AGE"),
    list(other = "keep_me")
  )
  out <- remove_covariate(covs, param = "CL", cov = "WT")
  expect_equal(length(out), 2)
  expect_equal(out[[1]]$COVNAME, "AGE")
  expect_equal(out[[2]]$other, "keep_me")

  out_again <- remove_covariate(out, param = "CL", cov = "WT")
  expect_equal(out_again, out)
})

test_that("small helpers normalize types and references consistently", {
  expect_equal(.norm_cov_type(c("cont", "continuous", "cat", "categorical", "other")),
               c("cont", "cont", "cat", "cat", NA_character_))
  expect_equal(.mode_sorted_smallest(c("B", "A", "B", "A")), "A")
  expect_true(is.na(.mode_sorted_smallest(c(NA_character_, NA_character_))))
  expect_equal(.na_safe_median(c("1", "2", NA)), 1.5)
  expect_true(is.na(.na_safe_median(c(NA_real_, NA_real_))))
  expect_equal(.covsearch_sanitize_name(" project@name! "), "project_name")
  expect_equal(.covsearch_sanitize_name("___"), "x")
})



test_that("sg_covsearch validates p-values and core arguments", {
  mock_paths <- .ensure_mock_files()
  gco <- list(
    model = mock_paths$model,
    data = mock_paths$data,
    headers = .covsearch_headers_fixture,
    theta = .covsearch_theta_fixture
  )
  gfo <- list(
    OFV = data.frame(LL = -100),
    COTAB = .covsearch_cotab_fixture,
    CATAB = .covsearch_catab_fixture
  )
  fit_stub <- function(...) stop("fit should not run")

  prep <- sg_covsearch(
    gco = gco,
    gfo = gfo,
    fit_function = fit_stub
  )
  expect_equal(prep$settings$p_forward, 0.05)
  expect_equal(prep$settings$p_backward, 0.01)
  expect_false(prep$metadata$forward_ran)

  expect_error(sg_covsearch(
    gco = gco, gfo = gfo, fit_function = fit_stub, p_forward = 1
  ), regexp = "p_forward")
  expect_error(sg_covsearch(
    gco = gco, gfo = gfo, fit_function = fit_stub, p_backward = 0
  ), regexp = "p_backward")
  expect_error(sg_covsearch(
    gco = gco, gfo = gfo, fit_function = fit_stub, run_backward = NA
  ), regexp = "run_backward")
  expect_error(sg_covsearch(
    gco = gco, gfo = gfo, fit_function = fit_stub, update_theta_init = NA
  ), regexp = "update_theta_init")
})

test_that("sg_covsearch validates covariates, parameters, and test_pairs", {
  gco <- .stage2_gco_fixture()
  gfo <- .stage2_gfo_fixture(ll = -100)

  expect_error(
    sg_covsearch(
      gco = gco, gfo = gfo, covariates = c("WT", "HEIGHT"),
      run_backward = FALSE, fit_function = .stage_mock_fit_project(list())
    ),
    regexp = "unknown covariates"
  )

  expect_error(
    sg_covsearch(
      gco = gco, gfo = gfo, parameters = c("CL", "Q"),
      run_backward = FALSE, fit_function = .stage_mock_fit_project(list())
    ),
    regexp = "unknown parameters"
  )

  expect_error(
    sg_covsearch(
      gco = gco, gfo = gfo, test_pairs = data.frame(parameter = "CL"),
      run_backward = FALSE, fit_function = .stage_mock_fit_project(list())
    ),
    regexp = "missing required columns"
  )

  bad_headers <- .covsearch_headers_fixture
  bad_headers[[4]]$type <- "binary"
  bad_gco <- gco
  bad_gco$headers <- bad_headers
  expect_error(
    sg_covsearch(
      gco = bad_gco, gfo = gfo, run_backward = FALSE, fit_function = .stage_mock_fit_project(list())
    ),
    regexp = "unsupported covariate type"
  )
})

test_that("Stage 2 generates full candidate cross-product when test_pairs is NULL", {
  gco <- .stage2_gco_fixture(project_name = "base_model")
  gfo <- .stage2_gfo_fixture(ll = -100)
  ofv_map <- list(
    base_model_fw_model_001 = -99,
    base_model_fw_model_002 = -99,
    base_model_fw_model_003 = -99,
    base_model_fw_model_004 = -99,
    base_model_fw_model_005 = -99,
    base_model_fw_model_006 = -99
  )

  out_dir <- file.path(tempdir(), paste0("covsearch_stage2_cross_", as.integer(stats::runif(1, 1, 1e9))))
  res <- sg_covsearch(
    gfo = gfo,
    gco = gco,
    output_dir = out_dir,
    covariates = c("WT", "AGE", "SEX"),
    parameters = c("CL", "V"),
    test_pairs = NULL,
    run_backward = FALSE,
    fit_function = .stage_mock_fit_project(ofv_map)
  )

  step1 <- res$forward$history[res$forward$history$step == 1, , drop = FALSE]
  expect_equal(nrow(step1), 6)
  expect_true(all(c("parameter", "covariate", "type", "df") %in% names(step1)))
  expect_equal(sort(unique(step1$parameter)), c("CL", "V"))
  expect_equal(sort(unique(step1$covariate)), c("AGE", "SEX", "WT"))
  expect_true(all(grepl("^base_model_fw_model_[0-9]{3}$", step1$project_name)))
  expect_true(.forward_history_files_exist(out_dir))
})

test_that("Stage 2 test_pairs restricts candidates and drops invalid rows with warning", {
  gco <- .stage2_gco_fixture(project_name = "base_model")
  gfo <- .stage2_gfo_fixture(ll = -100)
  ofv_map <- list(
    base_model_fw_model_001 = -99,
    base_model_fw_model_002 = -99,
    base_model_fw_model_003 = -99
  )

  expect_warning(
    res <- sg_covsearch(
      gco = gco,
      gfo = gfo,
      covariates = c("WT", "AGE", "SEX"),
      parameters = c("CL", "V"),
      test_pairs = .covsearch_test_pairs_fixture,
      run_backward = FALSE,
      fit_function = .stage_mock_fit_project(ofv_map)
    ),
    regexp = "dropped [0-9]+ invalid test_pairs row"
  )

  step1 <- res$forward$history[res$forward$history$step == 1, , drop = FALSE]
  expect_equal(nrow(step1), 3)
  expect_equal(sort(unique(step1$parameter)), c("CL", "V"))
  expect_equal(sort(unique(step1$covariate)), c("AGE", "SEX", "WT"))
})

test_that("Stage 2 categorical reference applies precedence and tie-break rules", {
  gco <- .stage2_gco_fixture(project_name = "base_model")
  catab_tie <- data.frame(
    ID = 1:4,
    SEX = c("M", "F", "M", "F"),
    RACE = c("A", "B", "A", "C"),
    stringsAsFactors = FALSE
  )
  gfo <- .stage2_gfo_fixture(ll = -100, catab = catab_tie)
  tp <- data.frame(
    parameter = c("CL", "V"),
    covariate = c("SEX", "RACE"),
    type = c("cat", "cat"),
    reference = c("M", NA),
    center = c(NA, NA),
    stringsAsFactors = FALSE
  )
  ofv_map <- list(
    base_model_fw_model_001 = -120,
    base_model_fw_model_002 = -110,
    base_model_fw_model_003 = -130
  )

  res <- sg_covsearch(
    gco = gco,
    gfo = gfo,
    covariates = c("SEX", "RACE"),
    parameters = c("CL", "V"),
    test_pairs = tp,
    run_backward = FALSE,
    fit_function = .stage_mock_fit_project(ofv_map)
  )

  covs <- res$final_gco$covs
  sex_cov <- Filter(function(x) identical(x$PAR, "CL") && identical(x$COVNAME, "SEX"), covs)
  race_cov <- Filter(function(x) identical(x$PAR, "V") && identical(x$COVNAME, "RACE"), covs)
  expect_equal(length(sex_cov), 1)
  expect_equal(length(race_cov), 1)
  expect_equal(sex_cov[[1]]$REF, "M")
  expect_equal(race_cov[[1]]$REF, "A")
})

test_that("Stage 2 df uses continuous=1 and categorical k-1 excluding NA", {
  gco <- .stage2_gco_fixture(project_name = "base_model")
  catab_with_na <- data.frame(
    ID = 1:7,
    SEX = c("M", "F", "M", "M", NA, "F", "M"),
    RACE = c("A", "B", "A", NA, "C", "B", "A"),
    stringsAsFactors = FALSE
  )
  gfo <- .stage2_gfo_fixture(ll = -100, catab = catab_with_na)
  ofv_map <- list(
    base_model_fw_model_001 = -99,
    base_model_fw_model_002 = -99,
    base_model_fw_model_003 = -99
  )

  res <- sg_covsearch(
    gco = gco,
    gfo = gfo,
    covariates = c("WT", "SEX", "RACE"),
    parameters = c("CL"),
    run_backward = FALSE,
    fit_function = .stage_mock_fit_project(ofv_map)
  )

  step1 <- res$forward$history[res$forward$history$step == 1, , drop = FALSE]
  wt_df <- unique(step1$df[step1$covariate == "WT"])
  sex_df <- unique(step1$df[step1$covariate == "SEX"])
  race_df <- unique(step1$df[step1$covariate == "RACE"])
  expect_equal(as.numeric(wt_df[[1]]), 1)
  expect_equal(as.numeric(sex_df[[1]]), 1)
  expect_equal(as.numeric(race_df[[1]]), 2)
})

test_that("Stage 3 mocked-fit selects best significant candidate each step", {
  gco <- .stage2_gco_fixture(project_name = "base_model")
  gfo <- .stage2_gfo_fixture(ll = 100)
  ofv_map <- list(
    base_model_fw_model_001 = 90,
    base_model_fw_model_002 = 95,
    base_model_fw_model_003 = 85
  )

  res <- sg_covsearch(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT"),
    parameters = c("CL", "V"),
    p_forward = 0.05,
    run_backward = FALSE,
    fit_function = .stage_mock_fit_project(ofv_map)
  )

  expect_true(res$metadata$forward_ran)
  expect_equal(res$metadata$forward_steps, 2)
  expect_equal(res$forward$selected$parameter, c("CL", "V"))
  expect_equal(res$forward$selected$covariate, c("WT", "WT"))
  expect_true(all(grepl("^base_model_fw_model_[0-9]{3}$", res$forward$history$project_name)))
})

test_that("Stage 3 stops when no candidate passes qchisq threshold", {
  gco <- .stage2_gco_fixture(project_name = "base_model")
  gfo <- .stage2_gfo_fixture(ll = 100)
  ofv_map <- list(
    base_model_fw_model_001 = 99.1,
    base_model_fw_model_002 = 99.2
  )

  res <- sg_covsearch(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT"),
    parameters = c("CL", "V"),
    p_forward = 0.05,
    run_backward = FALSE,
    fit_function = .stage_mock_fit_project(ofv_map)
  )

  expect_equal(nrow(res$forward$selected), 0)
  expect_true(all(res$forward$history$accepted == FALSE))
})

test_that("Stage 3 history stores threshold, df, and decision", {
  gco <- .stage2_gco_fixture(project_name = "base_model")
  gfo <- .stage2_gfo_fixture(ll = 100)
  ofv_map <- list(
    base_model_fw_model_001 = 95,
    base_model_fw_model_002 = 99,
    base_model_fw_model_003 = 94
  )

  res <- sg_covsearch(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT"),
    parameters = c("CL", "V"),
    p_forward = 0.05,
    run_backward = FALSE,
    fit_function = .stage_mock_fit_project(ofv_map)
  )

  expect_true(all(c("threshold", "df", "decision") %in% names(res$forward$history)))
  expect_true(any(res$forward$history$decision == "accepted"))
  expect_true(all(is.finite(res$forward$history$threshold)))
})

test_that("Stage 3 update_theta_init currently does not change theta INIT", {
  gco <- .stage2_gco_fixture(project_name = "base_model")
  gfo <- .stage2_gfo_fixture(ll = 100)
  ofv_map <- list(
    base_model_fw_model_001 = 90
  )
  sumtab_map <- list(
    base_model_fw_model_001 = data.frame(PAR = "CL_pop", VALUE = 5, stringsAsFactors = FALSE)
  )

  res <- sg_covsearch(
    gfo = gfo,
    gco = gco,
    covariates = c("WT"),
    parameters = c("CL"),
    p_forward = 0.05,
    update_theta_init = TRUE,
    run_backward = FALSE,
    fit_function = .stage_mock_fit_project(ofv_map, sumtab_map = sumtab_map)
  )

  cl_row <- res$final_gco$theta[res$final_gco$theta$NAME == "CL", , drop = FALSE]
  expect_equal(cl_row$INIT[[1]], .covsearch_theta_fixture$INIT[.covsearch_theta_fixture$NAME == "CL"])
})

test_that("Stage 4 removes least important removable term first", {
  gco <- .stage2_gco_fixture(project_name = "base_model")
  gfo <- .stage2_gfo_fixture(ll = 100)
  ofv_map <- list(
    base_model_fw_model_001 = 90,
    base_model_fw_model_002 = 95,
    base_model_fw_model_003 = 84,
    base_model_bw_model_001 = 89,
    base_model_bw_model_002 = 90,
    base_model_bw_model_003 = 97
  )

  res <- sg_covsearch(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT", "AGE"),
    parameters = c("CL"),
    p_forward = 0.05,
    p_backward = 0.01,
    run_backward = TRUE,
    fit_function = .stage_mock_fit_project(ofv_map)
  )

  expect_equal(res$forward$selected$covariate, c("WT", "AGE"))
  expect_equal(res$backward$removed$covariate[[1]], "WT")
  expect_equal(res$backward$retained$covariate, "AGE")
  expect_equal(res$final_covariates, res$backward$retained)
})

test_that("Stage 4 stops without removals when retained terms are significant", {
  gco <- .stage2_gco_fixture(project_name = "base_model")
  gfo <- .stage2_gfo_fixture(ll = 100)
  ofv_map <- list(
    base_model_fw_model_001 = 90,
    base_model_fw_model_002 = 95,
    base_model_fw_model_003 = 84,
    base_model_bw_model_001 = 93,
    base_model_bw_model_002 = 92
  )

  res <- sg_covsearch(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT", "AGE"),
    parameters = c("CL"),
    p_forward = 0.05,
    p_backward = 0.01,
    run_backward = TRUE,
    fit_function = .stage_mock_fit_project(ofv_map)
  )

  expect_equal(nrow(res$backward$removed), 0)
  expect_equal(res$backward$retained, res$forward$selected)
  expect_equal(res$metadata$backward_steps, 0)
})

test_that("Stage 4 skips backward loop cleanly when forward selected set is empty", {
  gco <- .stage2_gco_fixture(project_name = "base_model")
  gfo <- .stage2_gfo_fixture(ll = 100)
  ofv_map <- list(
    base_model_fw_model_001 = 99.5,
    base_model_fw_model_002 = 99.2
  )

  res <- sg_covsearch(
    gfo = gfo,
    gco = gco,
    output_dir = tempdir(),
    covariates = c("WT", "AGE"),
    parameters = c("CL"),
    p_forward = 0.05,
    p_backward = 0.01,
    run_backward = TRUE,
    fit_function = .stage_mock_fit_project(ofv_map)
  )

  expect_equal(nrow(res$forward$selected), 0)
  expect_equal(nrow(res$backward$history), 0)
  expect_equal(nrow(res$backward$retained), 0)
  expect_equal(res$final_covariates, res$backward$retained)
})

test_that("stepwise stops with clear error when mock fit fails", {
  gco <- .stage2_gco_fixture(project_name = "base_model")
  gfo <- .stage2_gfo_fixture(ll = 100)
  ofv_map <- list(
    base_model_fw_model_001 = 95
  )
  expect_error(
    sg_covsearch(
      gfo = gfo,
      gco = gco,
      covariates = c("WT"),
      parameters = c("CL"),
      p_forward = 0.05,
      run_backward = FALSE,
      fit_function = .stage_mock_fit_project(
        ofv_map = ofv_map,
        fail_projects = c("base_model_fw_model_001")
      )
    ),
    regexp = "fit failed in forward step"
  )
})

# ------------------------------
# Guarded fitter integration tests
# ------------------------------
test_that("sg_covsearch runs with Monolix fitter when available", {
  monolix_path <-"C:/ProgramData/Lixoft/MonolixSuite2023R1/bin/monolix.bat"
  skip_if(
    is.na(monolix_path),
    "Monolix executable not available on this machine."
  )
  #skip_if_not(exists("sg_fit", mode = "function"), "sg_fit() is not available in this session.")
  #skip_if_not(exists("sg_converter", mode = "function"), "sg_converter() is not available in this session.")


  fit_dir <- file.path( "workdir-amel", "projects", "warfarin-pk", "test-sg-covsearch")

  skip_if_not(dir.exists(fit_dir), "fitted_project fixtures are not available.")

  gco_path <- file.path(fit_dir, "wrfrn_pk_bm_ver2_GCO.json")
  gfo_path <- file.path(fit_dir, "wrfrn_pk_bm_ver2_GFO.json")
  #model_path <- file.path(fit_dir, "model-pk-1c.txt")
  #data_path <- file.path(fit_dir, "ds-warfarin-pk.csv")
  skip_if_not(file.exists(gco_path), "GCO fixture not present.")
  skip_if_not(file.exists(gfo_path), "GFO fixture not present.")
  #skip_if_not(file.exists(model_path), "Model fixture not present.")
  #skip_if_not(file.exists(data_path), "Data fixture not present.")

  #gco <- jsonlite::fromJSON(gco_path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
  #gfo <- jsonlite::fromJSON(gfo_path, simplifyVector = TRUE, simplifyDataFrame = FALSE)

  out_dir <- tempfile("covsearch-monolix-integration-")
  dir.create(out_dir, recursive = TRUE)
  on.exit(unlink(out_dir, recursive = TRUE, force = TRUE), add = TRUE)

  test_pairs <- data.frame(
    parameter = c("CL"),
    covariate = c("CYP2C9"),
    type = "cont",
    reference = NA_character_,
    center = "median",
    stringsAsFactors = FALSE
  )

  res <- sg_covsearch(
    gco = gco_path,
    gfo = gfo_path,
    #model = model_path,
    #data = data_path,
    output_dir = out_dir,
    covariates = c("CYP2C9") ,
    parameters = c("CL"),
    test_pairs = test_pairs,
    run_backward = FALSE,
    path_to_fitter = monolix_path
  )

  expect_true(nrow(res$forward$history) >= 1)
  expect_true(any(grepl("^wrfrn_pk_bm_ver2_fw_model_[0-9]{3}$", res$forward$history$project_name)))
  expect_true(all(is.finite(res$forward$history$candidate_ofv)))
  expect_true(.forward_history_files_exist(out_dir))

  project_name <- res$forward$history$project_name[[1]]
  expect_true(file.exists(file.path(out_dir, paste0(project_name, ".mlxtran"))))
  expect_true(dir.exists(file.path(out_dir, project_name)))
})



## Author: Alina Melnikova
## First created: 2026-04-22
## Description: formal testing of sg_fit file generation and validation
## Keywords: SimuRg, sg-fit

# ---- Shared fixtures -------------------------------------------------------

.fit_model <- normalizePath(
  system.file("extdata", "models", "model_PK_1c.txt", package = "SimuRg"),
  winslash = "/",
  mustWork = TRUE
)
.fit_data <- normalizePath(
  system.file("extdata", "datasets", "dspk-warf.csv", package = "SimuRg"),
  winslash = "/",
  mustWork = TRUE
)

.fit_headers <- list(
  list(name = "ID", use = "identifier", type = NULL),
  list(name = "TIME", use = "time", type = NULL),
  list(name = "DV", use = "observation", type = "continuous"),
  list(name = "DVID", use = "observationtype", type = NULL),
  list(name = "ADM", use = "administration", type = NULL),
  list(name = "AMT", use = "amount", type = NULL),
  list(name = "EVID", use = "eventidentifier", type = NULL),
  list(name = "MDV", use = "missingdependentvariable", type = NULL),
  list(name = "AGE", use = "covariate", type = "continuous")
)

.fit_theta <- data.frame(
  NAME = c("ka", "Vd", "CL"),
  TRANS = c("logNormal", "logNormal", "logNormal"),
  INIT = c(0.5, 20, 0.2),
  LB = c(NA_real_, NA_real_, NA_real_),
  UB = c(NA_real_, NA_real_, NA_real_),
  EST = c(TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

.fit_ruv <- list(
  YNAME = "y1",
  DVID = 1,
  TRANS = "normal",
  PRED = "Cc",
  ERR = "combined1",
  INIT = c(1, 0.3),
  EST = c(TRUE, FALSE),
  BLQM = NULL
)

.fit_re <- list(
  init = matrix(c(
    0.1, 0,   0,
    0,   0,   0,
    0,   0, 0.2
  ), nrow = 3, byrow = TRUE),
  est = matrix(c(
    TRUE, NA,   NA,
    NA,   NA,   NA,
    NA,   NA, TRUE
  ), nrow = 3, byrow = TRUE)
)

.fit_occ <- list(
  init = matrix(0, nrow = 3, ncol = 3),
  est = matrix(NA, nrow = 3, ncol = 3)
)

.fit_covs <- list(
  list(PAR = "CL", COVNAME = "AGE", FUNC = "linear",
       TRANS = "median", INIT = 0.01, EST = TRUE)
)

.sg_fit_args <- function(project_name, output_dir) {
  list(
    model = .fit_model,
    data = .fit_data,
    headers = .fit_headers,
    theta = .fit_theta,
    ruv = .fit_ruv,
    re = .fit_re,
    occ = .fit_occ,
    covs = .fit_covs,
    project_name = project_name,
    fit = FALSE,
    path_to_save_output = output_dir
  )
}

.find_monolix_batch <- function() {
  env_candidates <- c(
    Sys.getenv("SIMURG_MONOLIX_PATH", unset = ""),
    Sys.getenv("MONOLIX_PATH", unset = "")
  )

  lixoft_root <- "C:/ProgramData/Lixoft"
  install_dirs <- if (dir.exists(lixoft_root)) {
    list.dirs(lixoft_root, recursive = FALSE, full.names = TRUE)
  } else {
    character()
  }

  file_candidates <- c(
    env_candidates,
    file.path(install_dirs, "bin", "monolix.bat"),
    "C:/ProgramData/Lixoft/MonolixSuite2023R1/bin/monolix.bat"
  )

  file_candidates <- unique(file_candidates[nzchar(file_candidates)])
  existing <- file_candidates[file.exists(file_candidates)]

  if (length(existing) == 0) {
    return(NA_character_)
  }

  normalizePath(existing[[1]], winslash = "/", mustWork = TRUE)
}

# ============================================================
# 1. Input validation
# ============================================================

test_that("sg_fit errors when the model file is missing", {
  args <- .sg_fit_args(
    project_name = "missing-model",
    output_dir = tempdir()
  )
  args$model <- file.path(tempdir(), "does-not-exist-model.txt")

  expect_error(
    do.call(sg_fit, args),
    regexp = "model file does not exist"
  )
})

test_that("sg_fit errors when the data file is missing", {
  args <- .sg_fit_args(
    project_name = "missing-data",
    output_dir = tempdir()
  )
  args$data <- file.path(tempdir(), "does-not-exist-data.csv")

  expect_error(
    do.call(sg_fit, args),
    regexp = "data file does not exist"
  )
})

# ============================================================
# 2. Monolix project generation with fit = FALSE
# ============================================================

test_that("sg_fit writes an mlxtran project with the expected sections", {
  output_dir <- tempfile("sg-fit-write-")
  dir.create(output_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  project_name <- "unit_fit"
  result <- do.call(sg_fit, .sg_fit_args(project_name, output_dir))
  mlxtran_path <- file.path(output_dir, sprintf("%s.mlxtran", project_name))

  expect_null(result)
  expect_true(file.exists(mlxtran_path))

  mlxtran_text <- paste(readLines(mlxtran_path, warn = FALSE), collapse = "\n")

  expect_match(mlxtran_text, "<DATAFILE>", fixed = TRUE)
  expect_match(mlxtran_text, "[CONTENT]", fixed = TRUE)
  expect_match(mlxtran_text, "[INDIVIDUAL]", fixed = TRUE)
  expect_match(mlxtran_text, "[LONGITUDINAL]", fixed = TRUE)
  expect_match(mlxtran_text, "<FIT>", fixed = TRUE)
  expect_match(mlxtran_text, "<MONOLIX>", fixed = TRUE)
})

test_that("sg_fit mlxtran output reflects core headers and parameterization", {
  output_dir <- tempfile("sg-fit-content-")
  dir.create(output_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  project_name <- "content_fit"
  do.call(sg_fit, .sg_fit_args(project_name, output_dir))
  mlxtran_path <- file.path(output_dir, sprintf("%s.mlxtran", project_name))
  mlxtran_text <- paste(readLines(mlxtran_path, warn = FALSE), collapse = "\n")
  observed_dvid <- sort(unique(read.csv(.fit_data)[["DVID"]]))
  expected_yname <- paste0("'", observed_dvid, "'", collapse = ", ")
  expected_data_type <- paste0("'", observed_dvid, "'=plasma", collapse = ", ")

  expect_match(mlxtran_text, "header=\\{ID, TIME, DV, DVID", perl = TRUE)
  expect_match(mlxtran_text, "ID = \\{use=identifier\\}", perl = TRUE)
  expect_match(
    mlxtran_text,
    sprintf("DV = {use=observation, yname=%s, type=continuous}", expected_yname),
    fixed = TRUE
  )
  expect_match(
    mlxtran_text,
    sprintf("dataType = {%s}", expected_data_type),
    fixed = TRUE
  )
  expect_match(mlxtran_text, "\\[COVARIATE\\]\\ninput = \\{AGE\\}", perl = TRUE)
  expect_match(
    mlxtran_text,
    "ka = \\{distribution=logNormal, typical=ka_pop, sd=omega_ka\\}",
    perl = TRUE
  )
  expect_match(
    mlxtran_text,
    "Vd = \\{distribution=logNormal, typical=Vd_pop, no-variability\\}",
    perl = TRUE
  )
  # expect_match(
  #   mlxtran_text,
  #   "CL = \\{distribution=logNormal, typical=CL_pop, covariate=AGE, coefficient=beta_CL_AGE, sd=omega_CL\\}",
  #   perl = TRUE
  # )
  expect_match(
    mlxtran_text,
    "y1 = \\{distribution=normal, prediction=Cc, errorModel=combined1\\(Cc_a, Cc_b\\)\\}",
    perl = TRUE
  )
  expect_match(mlxtran_text, "ka_pop = \\{value=0.5, method=MLE\\}", perl = TRUE)
  expect_match(mlxtran_text, "Vd_pop = \\{value=20, method=FIXED\\}", perl = TRUE)
  expect_match(mlxtran_text, "CL_pop = \\{value=0.2, method=MLE\\}", perl = TRUE)
  expect_match(mlxtran_text, "Cc_a = \\{value=1, method=MLE\\}", perl = TRUE)
  expect_match(mlxtran_text, "Cc_b = \\{value=0.3, method=FIXED\\}", perl = TRUE)
  expect_match(mlxtran_text, "omega_ka = \\{value=0.1, method=MLE\\}", perl = TRUE)
  expect_match(mlxtran_text, "omega_CL = \\{value=0.2, method=MLE\\}", perl = TRUE)
  # expect_match(mlxtran_text, "beta_CL_AGE = \\{value=0.01, method=MLE\\}", perl = TRUE)
  expect_match(mlxtran_text, "exportpath = 'content_fit'", perl = TRUE)
})
# ============================================================
# 3. Optional Monolix integration
# ============================================================

test_that("sg_fit fit = TRUE skips cleanly when Monolix is unavailable", {
  monolix_path <- .find_monolix_batch()
  skip_if(
    is.na(monolix_path),
    "Monolix executable not available on this machine."
  )

  output_dir <- tempfile("sg-fit-monolix-")
  dir.create(output_dir, recursive = TRUE)
  on.exit(unlink(output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  project_name <- "integration_fit"
  args <- .sg_fit_args(project_name, output_dir)
  args$fit <- TRUE
  args$path_to_fitter <- monolix_path
  args$max_wait_time <- 120

  result <- do.call(sg_fit, args)
  project_dir <- file.path(output_dir, project_name)
  mlxtran_path <- file.path(output_dir, sprintf("%s.mlxtran", project_name))

  expect_true(file.exists(mlxtran_path))
  expect_true(dir.exists(project_dir))
  expect_false(is.null(result))
})


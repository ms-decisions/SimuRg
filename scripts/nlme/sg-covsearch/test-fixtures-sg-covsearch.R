## Author: Alina Melnikova
## First created: 2026-05-19
## Description: reusable fixtures for sg-covsearch tests
## Keywords: SimuRg, covsearch, fixtures

# Minimal headers fixture for covsearch preprocessing:
# - continuous covariates: WT, AGE
# - two categorical covariates: SEX (2 levels), RACE (3 levels)
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

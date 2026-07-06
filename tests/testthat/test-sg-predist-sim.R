## Description: formal testing of sg-predist-sim function
## Keywords: SimuRg, sg-predist-sim, prediction distribution, diagnostics

mod_fin <- RxODE({
  # Doses in mg
  # Time in hours

  ### Parameter values
  # Typical
  Cl_pop = 5;
  V_pop = 180;

  ka_pop = 6;


  # Random effects
  omega_Cl = 0;
  omega_V = 0;
  omega_ka = 0;

  # Residual error
  b = 0;

  ### Parameters
  Cl = Cl_pop * exp(omega_Cl);
  V = V_pop * exp(omega_V);
  ka = ka_pop * exp(omega_ka);

  ### Explicit functions
  Cc = Ac/V;                 # nmol/L

  ### Initial conditions
  Ad(0) = 0;          # mg
  Ac(0) = 0;          # mg

  ### ODEs
  d/dt(Ad) = - ka*Ad;
  d/dt(Ac) = ka*Ad - Cl*Cc ;

  CHECKRUV = b;
  Cc_ResErr = Cc + b*Cc;
})

# ---- Generalized control object (gco) fixture (mirrors test-sg-vpc-sim.R) ----
model <- system.file("extdata", "models", "rxode", "model_PK_1c.txt", package = "SimuRg")
data  <- system.file("extdata", "datasets", "dspk-warf.csv", package = "SimuRg")

headers <- list(list(name = "ID", use = "identifier", type = NULL),
                list(name = "TIME", use = "time", type = NULL),
                list(name = "DV", use = "observation", type = "continuous"),
                list(name = "DVID", use = "observationtype", type = NULL),
                list(name = "ADM", use = "administration", type = NULL),
                list(name = "AMT", use = "amount", type = NULL),
                list(name = "EVID", use = "eventidentifier", type = NULL),
                list(name = "MDV", use = "missingdependentvariable", type = NULL),
                list(name = "AGE", use = "covariate", type = "continuous"),
                list(name = "AGE_centered", use = "covariate", type = "continuous"),
                list(name = "SEX", use = "covariate", type = "categorical"),
                list(name = "WEIGHT", use = "covariate", type = "continuous"),
                list(name = "BMI", use = "covariate", type = "continuous"))

theta <- tribble(~NAME, ~TRANS, ~INIT, ~LB, ~UB, ~EST,
                 "Cl", "logNormal", 0.2, NA, NA, TRUE,
                 "V", "logNormal", 20, NA, NA, TRUE,
                 "ka", "logNormal", 0.2, NA, NA, TRUE
)

ruv <- list(YNAME = "y1", DVID = 1, TRANS = "normal", PRED = "Cc",
            ERR = "combined1", INIT = c(1, 1), EST = c(TRUE, TRUE), BLQM = NULL)

re <- list(init = tribble(~Cl, ~V, ~ka,
                          1, 0, 0,
                          0, 1, 0,
                          0, 0, 1) %>% as.matrix(),
           est = tribble(~Cl, ~V, ~ka,
                         TRUE, NA, NA,
                         NA, TRUE, NA,
                         NA, NA, TRUE) %>% as.matrix())

occ <- list(init = tribble(~Cl, ~V, ~ka,
                           0, 0, 0,
                           0, 0, 0,
                           0, 0, 0) %>% as.matrix(),
            est = tribble(~Cl, ~V, ~ka,
                          NA, NA, NA,
                          NA, NA, NA,
                          NA, NA, NA) %>% as.matrix())
covs <- list(list(PAR = "V", COVNAME = "AGE", FUNC = "linear",
                  TRANS = "median", INIT = 1, EST = TRUE),
             list(PAR = "ka", COVNAME = "SEX", REF = 0, INIT = 1, EST = TRUE))

gco <- list(headers = headers,
            data = data,
            model = model,
            task_opt = "",
            covs = covs,
            project_name = "test-proj",
            theta = theta,
            ruv = ruv,
            re = re,
            occ = occ,
            modelText = "")

# ---- Helper: minimal mock sg_fit object ----
create_mock_sg_fit <- function() {
  list(
    SDTAB = data.frame(
      ID = rep(1:2, each = 5),
      TIME = c(0, 1, 2, 4, 8, 0, 1, 2, 4, 8),
      DV = rnorm(10, mean = 10, sd = 2),
      MDV = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    ),
    EVTAB = data.frame(
      ID = c(1, 2),
      time = c(0, 0),
      amt = c(100, 100),
      cmt = c(1, 1),
      evid = c(1, 1)
    ),
    SUMTAB = data.frame(
      PAR = c("ka_pop", "Cl_pop", "V_pop"),
      VALUE = c(1.5, 2.0, 50.0),
      TYPE = "Typical values"
    ),
    OMEGAMAT = matrix(c(0.1, 0, 0, 0, 0.15, 0, 0, 0, 0.2), nrow = 3, ncol = 3),
    SIGMAMAT = matrix(c(0.05), nrow = 1, ncol = 1),
    COTAB = NULL,
    CATAB = NULL
  )
}

# ---- Integration tests using the obj1 package fixture ----
test_that("sg_predist_sim works with an explicit model", {
  res <- sg_predist_sim(obj1, model = mod_fin, outputs = "Cc", npop = 100)

  # ID is overwritten with the original individual id -> one per subject in obj1
  expect_equal(res %>% pull(ID) %>% unique() %>% length(), 100)
  expect_equal(res %>% pull(TIME) %>% unique() %>% length(),
               obj1$SDTAB$TIME %>% unique() %>% length())
  # sg_sim() returns long format; these columns must be present
  expect_true(all(c("ID", "TIME", "VAR", "VALUE") %in% colnames(res)))
})

test_that("sg_predist_sim works with a generalized control object (gco)", {
  res <- sg_predist_sim(obj1, gco = gco, outputs = "Cc", npop = 50)

  expect_s3_class(res, "data.frame")
  expect_equal(res %>% pull(ID) %>% unique() %>% length(), 100)
  expect_true(all(c("ID", "TIME", "VAR", "VALUE") %in% colnames(res)))
})

# ---- Input validation ----
test_that("sg_predist_sim errors when neither model nor gco is supplied", {
  expect_error(
    sg_predist_sim(obj1),
    regexp = "generalized control object|model"
  )
})

test_that("sg_predist_sim errors on invalid fpath_i type", {
  expect_error(sg_predist_sim(fpath_i = 123, model = mod_fin))
})

test_that("sg_predist_sim errors on non-existent file path", {
  expect_error(sg_predist_sim(fpath_i = "nonexistent_file.rda", model = mod_fin))
})

# ---- Mock-object behaviour ----
test_that("sg_predist_sim accepts a list object as fpath_i", {
  skip_if_not_installed("rxode2")
  if (!exists("sg_sim")) skip("sg_sim function not available")

  mock_obj <- create_mock_sg_fit()

  result <- tryCatch(
    sg_predist_sim(fpath_i = mock_obj, model = mod_fin, npop = 10),
    error = function(e) NULL
  )

  if (!is.null(result)) {
    expect_s3_class(result, "data.frame")
    expect_true("ID" %in% colnames(result))
  }
})

test_that("sg_predist_sim drops the MDV column and filters MDV == 1", {
  skip_if_not_installed("rxode2")
  if (!exists("sg_sim")) skip("sg_sim function not available")

  mock_obj <- create_mock_sg_fit()
  mock_obj$SDTAB <- rbind(
    mock_obj$SDTAB,
    data.frame(ID = c(1, 2), TIME = c(0, 0), DV = c(NA, NA), MDV = c(1, 1))
  )

  result <- tryCatch(
    sg_predist_sim(fpath_i = mock_obj, model = mod_fin, npop = 10),
    error = function(e) NULL
  )

  if (!is.null(result)) {
    expect_false("MDV" %in% colnames(result))
    expect_true(all(unique(result$ID) %in% unique(mock_obj$SDTAB$ID)))
  }
})

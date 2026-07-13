## Build bundled GMO (gmo_pk1c) and GSI (gsi_pk1c).
## GMO: from inst/extdata/simurg_object/gmo_pk1c_model.txt
## GSI: assembled inline (sg_sim() arguments).

root <- normalizePath(getwd(), winslash = "/")
model_path <- file.path(root, "inst", "extdata", "simurg_object", "gmo_pk1c_model.txt")

if (!file.exists(model_path)) {
  stop("Missing model source: ", model_path)
}

## --- GMO ---
model_code <- paste(readLines(model_path, warn = FALSE), collapse = "\n")

suppressPackageStartupMessages(library(rxode2))
gmo_pk1c <- eval(parse(text = paste0("rxode2::rxode2({", model_code, "})")))

## --- GSI ---
theta_par <- c(
  "ka_pop", "Vd_pop", "CL_pop", "beta_Vd_WT",
  "omega_ka", "omega_Vd", "omega_CL", "Cc_b"
)

theta <- c(
  ka_pop = -0.03180134486167579,
  Vd_pop = 1.0359883779381744,
  CL_pop = 0.24965304895423984,
  beta_Vd_WT = 0.0210163
)

thetamat <- matrix(
  c(
    0.000856896, 0.000881042, 0.000050138, -8.8635e-7, 4.62689e-13, 8.85743e-14, 1.18272e-11, 1.9793e-12,
    0.000881042, 0.337007, -0.000294842, -0.00185967, 1.8197e-10, 3.48352e-11, 4.65152e-9, 7.78437e-10,
    0.000050138, -0.000294842, 0.00435795, 0.00000246416, 2.35311e-12, 4.50466e-13, 6.01503e-11, 1.00662e-11,
    -8.8635e-7, -0.00185967, 0.00000246416, 0.0000103604, 5.59421e-15, 1.07092e-15, 1.42999e-13, 2.39311e-14,
    4.62689e-13, 1.8197e-10, 2.35311e-12, 5.59421e-15, 0.00000539959, 0.00000195658, 0.00000980014, 0.000002175,
    8.85743e-14, 3.48352e-11, 4.50466e-13, 1.07092e-15, 0.00000195658, 0.00000103366, 0.0000061597, 0.00000214883,
    1.18272e-11, 4.65152e-9, 6.01503e-11, 1.42999e-13, 0.00000980014, 0.0000061597, 0.000138024, 0.0000278457,
    1.9793e-12, 7.78437e-10, 1.00662e-11, 2.39311e-14, 0.000002175, 0.00000214883, 0.0000278457, 0.0000230985
  ),
  nrow = 8,
  ncol = 8,
  byrow = TRUE,
  dimnames = list(theta_par, theta_par)
)

et <- data.frame(
  ID = c(1, 2),
  TIME = c(0, 0),
  AMT = c(10, 10),
  ADDL = c(10, 10),
  II = c(24, 24),
  WT = c(60, 200),
  stringsAsFactors = FALSE
)

gsi_pk1c <- list(
  et = et,
  stimes = seq(0, 120, length.out = 200),
  outputs = "Cc",
  covs = "WT",
  theta = theta,
  thetamat = thetamat,
  omega = NULL,
  sigma = NULL,
  npop = 10,
  nsub = 1,
  addcov = FALSE
)

if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(root, quiet = TRUE)
  sim_check <- do.call(sg_sim, c(list(model = gmo_pk1c), gsi_pk1c))
  stopifnot(all(c("ID", "TIME", "VAR", "VALUE", "POPN") %in% names(sim_check)))
  message("Sanity check: ", nrow(sim_check), " GSO rows")
}

usethis::use_data(gmo_pk1c, gsi_pk1c, overwrite = TRUE)
message("Saved gmo_pk1c and gsi_pk1c to data/")

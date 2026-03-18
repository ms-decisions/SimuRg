## Scenarios aligned with roxygen @examples in R/sg-sim.R
## All inputs (model, theta, thetamat, omega, sigma) are explicit R objects — no file uploads.
##
## Test structure:
##   (1) Base 1-compartment PK, no covariate — sim1, sim2a–c, sim3a–b, sim4a–d, parameter override.
##   (2) 1-compartment PK with WTBL covariate on Vd — covs/keep and covariate + variability (sim5).

test_that("sg_sim runs with 1-compartment model and explicit parameters", {
  # ---- Model: 1-compartment PK, no covariate ----
  mod <- rxode2::rxode2({
    ka_pop = 0.1;
    Vd_pop = 10;
    CL_pop = 0.5;
    omega_ka = 0;
    omega_Vd = 0;
    omega_CL = 0;
    Cc_b = 0;
    ka_tv = exp(ka_pop);
    Vd_tv = exp(Vd_pop);
    CL_tv = exp(CL_pop);
    ka <- ka_tv * exp(omega_ka);
    Vd <- Vd_tv * exp(omega_Vd);
    CL <- CL_tv * exp(omega_CL);
    Cc = Ac / Vd;
    Ad(0) = 0;
    Ac(0) = 0;
    d/dt(Ad) = -ka * Ad;
    d/dt(Ac) = ka * Ad - CL * Cc;
    Cc_ResErr = Cc * (1 + Cc_b);
  })

  et_ex <- tibble::tribble(
    ~ID, ~TIME, ~EVID, ~CMT, ~AMT,
    1, 0, 1, 1, 10,
    2, 0, 1, 1, 50
  )
  stimes_ex <- seq(0, 24, 0.1)
  output_ex <- c("Ac", "Ad", "Cc", "Cc_ResErr")
  theta_ex <- c(ka_pop = log(0.1), Vd_pop = log(10), CL_pop = log(0.5))
  thetamat_ex <- matrix(c(0.05, 0.01, 0, 0.01, 0.05, 0, 0, 0, 0.05), nrow = 3, byrow = TRUE)
  rownames(thetamat_ex) <- colnames(thetamat_ex) <- c("ka_pop", "Vd_pop", "CL_pop")
  omega_ex <- matrix(c(0.2, 0.05, 0, 0.05, 0.2, 0, 0, 0, 0.2), nrow = 3, byrow = TRUE)
  rownames(omega_ex) <- colnames(omega_ex) <- c("omega_ka", "omega_Vd", "omega_CL")
  sigma_ex <- matrix(0.1)
  rownames(sigma_ex) <- colnames(sigma_ex) <- "Cc_b"

  # ---- sim1: No variability (theta only); single solve, no thetamat/omega ----
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex, outputs = output_ex)
  )

  # ---- sim2a: Population uncertainty, single scenario (ID) ----
  expect_no_error(
    sg_sim(model = mod, et = dplyr::filter(et_ex, ID == 1), stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = 1)
  )
  expect_no_error(
    sg_sim(model = mod, et = dplyr::filter(et_ex, ID == 1), stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = 10)
  )

  # ---- sim2b: Population uncertainty, multiple IDs, byID=TRUE (populations replicated per ID), shared=TRUE (same populations for each ID) ----
  expect_warning(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = nrow(et_ex), byID = TRUE, shared = TRUE)
  )
  expect_warning(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = 10, byID = TRUE, shared = TRUE)
  )

  # ---- Population uncertainty, multiple IDs, byID=TRUE shared=FALSE (separate pop per ID) ----
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = nrow(et_ex), byID = TRUE, shared = FALSE)
  )

  # ---- sim2c: Population uncertainty, byID=FALSE (one solve over full event table; npop must equal n(ID); ID = virtual subject) ----
  expect_warning(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = nrow(et_ex), byID = FALSE)
  )
  expect_warning(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = 10, byID = FALSE)
  )

  # ---- Between-subject variability (BSV), single ID ----
  expect_no_error(
    sg_sim(model = mod, et = dplyr::filter(et_ex, ID == 1), stimes = stimes_ex, theta = theta_ex,
           omega = omega_ex, nsub = 1)
  )
  expect_no_error(
    sg_sim(model = mod, et = dplyr::filter(et_ex, ID == 1), stimes = stimes_ex, theta = theta_ex,
           omega = omega_ex, nsub = 10)
  )

  # ---- BSV, multiple IDs, byID=TRUE shared=TRUE ----
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           omega = omega_ex, nsub = 10, byID = TRUE, shared = TRUE)
  )

  # ---- sim3a: BSV, multiple IDs, byID=TRUE (subjects replicated per ID), shared=FALSE (separate subjects per ID) ----
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           omega = omega_ex, nsub = 10, byID = TRUE, shared = FALSE)
  )
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           omega = omega_ex, nsub = 5, byID = TRUE, shared = FALSE)
  )

  # ---- sim3b: BSV, byID=FALSE (one solve over full event table; nsub must equal n(ID); ID = virtual subject) ----
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           omega = omega_ex, nsub = nrow(et_ex), byID = FALSE)
  )

  # ---- Both thetamat and omega, byID=FALSE (nsub = n(ID)) ----
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, omega = omega_ex, nsub = nrow(et_ex), byID = FALSE)
  )

  # ---- sim4a: BSV + uncertainty, byID=TRUE (pop/subjects replicated per ID), byPOP=TRUE (subjects replicated per population), shared=FALSE ----
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = 10, omega = omega_ex, nsub = 5,
           byID = TRUE, byPOP = TRUE, shared = FALSE)
  )

  # ---- sim4b: BSV + uncertainty, byID=TRUE, byPOP=FALSE (subjects not replicated between populations; npop = nsub) ----
  # expect_no_error(
  #   sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
  #          thetamat = thetamat_ex, npop = 5, omega = omega_ex, nsub = 5,
  #          byID = TRUE, byPOP = FALSE)
  # )

  # ---- sim4c: BSV + uncertainty, byID=FALSE (one solve; nsub/npop vs n(ID)), byPOP=TRUE, shared=FALSE ----
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = 10, omega = omega_ex, nsub = 5,
           byID = FALSE, byPOP = TRUE, shared = FALSE)
  )
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = 10, omega = omega_ex, nsub = 2,
           byID = FALSE, byPOP = TRUE, shared = FALSE)
  )
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = 2, omega = omega_ex, nsub = 10,
           byID = FALSE, byPOP = TRUE, shared = FALSE)
  )

  # ---- sim4d: BSV + uncertainty, byID=FALSE (one solve; nsub and npop equal n(ID)), byPOP=FALSE (npop = nsub) ----
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = nrow(et_ex), omega = omega_ex, nsub = nrow(et_ex),
           byID = FALSE, byPOP = FALSE)
  )
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = 10, omega = omega_ex, nsub = 5,
           byID = FALSE, byPOP = FALSE)
  )
  expect_no_error(
    sg_sim(model = mod, et = et_ex, stimes = stimes_ex, theta = theta_ex,
           thetamat = thetamat_ex, npop = 2, omega = omega_ex, nsub = 2,
           byID = FALSE, byPOP = FALSE)
  )

  # ---- Parameter override in event table (Vd_pop, omega_Vd) ----
  expect_no_error(
    sg_sim(model = mod, et = dplyr::filter(et_ex, ID == 1) %>% dplyr::mutate(Vd_pop = 2),
           stimes = stimes_ex, theta = theta_ex, thetamat = thetamat_ex, npop = 1)
  )
  expect_no_error(
    sg_sim(model = mod, et = dplyr::filter(et_ex, ID == 1) %>% dplyr::mutate(Vd_pop = 2),
           stimes = stimes_ex, theta = theta_ex, thetamat = thetamat_ex, npop = 10)
  )
  expect_no_error(
    sg_sim(model = mod, et = dplyr::filter(et_ex, ID == 1) %>% dplyr::mutate(Vd_pop = 2),
           stimes = stimes_ex, theta = theta_ex, omega = omega_ex, nsub = 1)
  )
  expect_no_error(
    sg_sim(model = mod, et = dplyr::filter(et_ex, ID == 1) %>% dplyr::mutate(Vd_pop = 2),
           stimes = stimes_ex, theta = theta_ex, omega = omega_ex, nsub = 10)
  )
  expect_no_error(
    sg_sim(model = mod, et = dplyr::filter(et_ex, ID == 1) %>%
             dplyr::mutate(Vd_pop = 2, omega_Vd = 0.5),
           stimes = stimes_ex, theta = theta_ex, omega = omega_ex, nsub = 10)
  )
})

test_that("sg_sim runs with covariate model and explicit parameters (sim5)", {
  # ---- Model: 1-compartment PK with WTBL covariate on Vd (Vd_tv = exp(Vd_pop) * (WTBL/70)^beta_WTBL_Vd_pop) ----
  mod_cov <- rxode2::rxode2({
    ka_pop = 0.1;
    Vd_pop = 10;
    CL_pop = 0.5;
    beta_WTBL_Vd_pop = 0;
    omega_ka = 0;
    omega_Vd = 0;
    omega_CL = 0;
    Cc_b = 0;
    ka_tv = exp(ka_pop);
    Vd_tv = exp(Vd_pop) * (WTBL / 70)^beta_WTBL_Vd_pop;
    CL_tv = exp(CL_pop);
    ka <- ka_tv * exp(omega_ka);
    Vd <- Vd_tv * exp(omega_Vd);
    CL <- CL_tv * exp(omega_CL);
    Cc = Ac / Vd;
    Ad(0) = 0;
    Ac(0) = 0;
    d/dt(Ad) = -ka * Ad;
    d/dt(Ac) = ka * Ad - CL * Cc;
    Cc_ResErr = Cc * (1 + Cc_b);
  })

  et_ex <- tibble::tribble(
    ~ID, ~TIME, ~EVID, ~CMT, ~AMT,
    1, 0, 1, 1, 10,
    2, 0, 1, 1, 50
  )
  stimes_ex <- seq(0, 24, 0.1)
  theta_ex <- c(ka_pop = log(0.1), Vd_pop = log(10), CL_pop = log(0.5))
  thetamat_ex <- matrix(c(0.05, 0.01, 0, 0.01, 0.05, 0, 0, 0, 0.05), nrow = 3, byrow = TRUE)
  rownames(thetamat_ex) <- colnames(thetamat_ex) <- c("ka_pop", "Vd_pop", "CL_pop")
  omega_ex <- matrix(c(0.2, 0.05, 0, 0.05, 0.2, 0, 0, 0, 0.2), nrow = 3, byrow = TRUE)
  rownames(omega_ex) <- colnames(omega_ex) <- c("omega_ka", "omega_Vd", "omega_CL")
  theta_cov <- c(theta_ex, beta_WTBL_Vd_pop = 0.5)

  # ---- sim5: Covariate (WTBL); pass covs and keep so WTBL is used and returned ----
  expect_no_error(
    sg_sim(model = mod_cov, et = dplyr::filter(et_ex, ID == 1) %>% dplyr::mutate(WTBL = 70),
           stimes = stimes_ex, theta = theta_cov, covs = c("WTBL"), keep = c("WTBL"),
           thetamat = thetamat_ex, npop = 1)
  )
  expect_no_error(
    sg_sim(model = mod_cov, et = dplyr::filter(et_ex, ID == 1) %>% dplyr::mutate(WTBL = 70),
           stimes = stimes_ex, theta = theta_cov, covs = c("WTBL"), keep = c("WTBL"),
           thetamat = thetamat_ex, npop = 5)
  )
  expect_no_error(
    sg_sim(model = mod_cov, et = dplyr::filter(et_ex, ID == 1) %>% dplyr::mutate(WTBL = 70),
           stimes = stimes_ex, theta = theta_cov, covs = c("WTBL"), keep = c("WTBL"),
           thetamat = thetamat_ex, npop = 10)
  )

  # ---- Covariate, multiple IDs with different WTBL, population uncertainty ----
  expect_warning(
    sg_sim(model = mod_cov, et = dplyr::mutate(et_ex, WTBL = c(70, 80)),
           stimes = stimes_ex, theta = theta_cov, covs = c("WTBL"), keep = c("WTBL"),
           thetamat = thetamat_ex, npop = 10, byID = TRUE, shared = TRUE)
  )

  # ---- Covariate, single ID, BSV only ----
  expect_no_error(
    sg_sim(model = mod_cov, et = dplyr::filter(et_ex, ID == 1) %>% dplyr::mutate(WTBL = 70),
           stimes = stimes_ex, theta = theta_cov, covs = c("WTBL"), keep = c("WTBL"),
           omega = omega_ex, nsub = 1)
  )
  expect_no_error(
    sg_sim(model = mod_cov, et = dplyr::filter(et_ex, ID == 1) %>% dplyr::mutate(WTBL = 70),
           stimes = stimes_ex, theta = theta_cov, covs = c("WTBL"), keep = c("WTBL"),
           omega = omega_ex, nsub = 10)
  )

  # ---- Covariate + uncertainty + BSV, byID=TRUE byPOP=TRUE shared=FALSE ----
  expect_no_error(
    sg_sim(model = mod_cov, et = dplyr::filter(et_ex, ID == 1) %>% dplyr::mutate(WTBL = 70),
           stimes = stimes_ex, theta = theta_cov, covs = c("WTBL"), keep = c("WTBL"),
           thetamat = thetamat_ex, npop = 5, omega = omega_ex, nsub = 3,
           byID = TRUE, byPOP = TRUE, shared = FALSE)
  )

  # ---- Covariate, multiple IDs, byID=TRUE shared=FALSE with omega ----
  expect_no_error(
    sg_sim(model = mod_cov, et = dplyr::mutate(et_ex, WTBL = c(70, 80)),
           stimes = stimes_ex, theta = theta_cov, covs = c("WTBL"), keep = c("WTBL"),
           omega = omega_ex, nsub = 5, byID = TRUE, shared = FALSE)
  )
})

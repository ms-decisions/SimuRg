## Author: Mikhailova Anna
## First created: 2025-09-11
## Description: functions for simulations from the model
## Keywords: SimuRg, simulations

#' Perform simulations from an rxode2 model
#'
#' Runs simulations using a [GMO] with explicit R objects for parameters (theta),
#' variance–covariance (thetamat), omega, and sigma. Argument layout is described
#' in [GSI].
#'
#' @inheritParams sg_dummy
#' @param byID logical. If \code{TRUE}, replicate populations/subjects per event-table ID.
#'   Default \code{NULL} is treated as \code{TRUE}.
#' @param byPOP logical. When both thetamat and omega are used: if \code{TRUE},
#'   replicate subjects by population. Default \code{NULL} is treated as \code{TRUE}.
#' @param shared logical. If \code{TRUE},
#'   one shared population draw for all IDs. Default \code{NULL} is treated as \code{TRUE}.
#' @param ... optional arguments passed to \code{rxode2::rxSolve} (e.g. \code{method}, \code{maxsteps}).
#' @returns [GSO]: a data frame with simulation results (long format).
#' @details
#'  Index columns \code{ID}, \code{POPN}, and \code{sim.id} are included depending on what is used:
#'  \describe{
#'     \item{No uncertainty, no variability}{Only \code{ID} (and \code{TIME}, \code{VAR}, \code{VALUE})}
#'     \item{Only uncertainty (\code{thetamat})}{ \code{POPN} and \code{ID}}
#'     \item{Only variability (\code{omega})}{\code{sim.id} and \code{ID}}
#'     \item{Uncertainty and variability}{\code{POPN}, \code{sim.id}, and \code{ID}}
#'  }
#'  Other columns:
#'  \describe{
#'     \item{TIME}{Timepoint of the simulations}
#'     \item{VAR}{Names of the output variables}
#'     \item{VALUE}{Optional. Simulated value. Returned when no aggregation applied}
#'     \item{mean, median, ...}{Optional. Aggregated statistic of simulated values. Returned when aggregation is applied}
#'     \item{COV1, ...COVn}{Optional. Covariates value. Returned when \code{addcov} is \code{TRUE}}
#'     \item{KEEP1, ...KEEPn}{Optional. Columns that were specified in the \code{keep} argument}
#'  }
#'  Event table can use either \code{id} or \code{ID}. Parameters can be fixed per ID by
#'  including parameter names (e.g. \code{*_pop}, \code{omega_*}) as columns in the event table.
#' @examples
#' \donttest{
#' library(rxode2)
#' library(dplyr)
#' library(tibble)
#'
#' # Bundled [GMO] and [GSI] (see ?gmo_pk1c, ?gsi_pk1c)
#' sim_gmo <- do.call(sg_sim, c(list(model = gmo_pk1c), gsi_pk1c))
#' head(sim_gmo)
#'
#' # Base 1-compartment PK model (no covariate)
#' mod <- rxode2::rxode2({
#'   ka_pop = 0.1;
#'   Vd_pop = 10;
#'   CL_pop = 0.5;
#'
#'   omega_ka = 0;
#'   omega_Vd = 0;
#'   omega_CL = 0;
#'
#'   Cc_b = 0;
#'   ka_tv = exp(ka_pop);
#'   Vd_tv = exp(Vd_pop);
#'   CL_tv = exp(CL_pop);
#'
#'   ka = ka_tv * exp(omega_ka);
#'   Vd = Vd_tv * exp(omega_Vd);
#'   CL = CL_tv * exp(omega_CL);
#'
#'   Cc = Ac / Vd;
#'
#'   Ad(0) = 0;
#'   Ac(0) = 0;
#'
#'   d/dt(Ad) = -ka * Ad;
#'   d/dt(Ac) = ka * Ad - CL * Cc;
#'
#'   Cc_ResErr = Cc * (1 + Cc_b);
#' })
#'
#' # Model with covariate: WTBL on Vd (Vd_tv = exp(Vd_pop) * (WTBL/70)^beta_WTBL_Vd_pop)
#'
#' mod_cov <- rxode2::rxode2({
#'   ka_pop = 0.1;
#'   Vd_pop = 10;
#'   CL_pop = 0.5;
#'
#'   beta_WTBL_Vd_pop = 0;
#'
#'   omega_ka = 0;
#'   omega_Vd = 0;
#'   omega_CL = 0;
#'
#'   Cc_b = 0;
#'
#'   ka_tv = exp(ka_pop);
#'   Vd_tv = exp(Vd_pop) * (WTBL/70)^beta_WTBL_Vd_pop;
#'   CL_tv = exp(CL_pop);
#'
#'   ka = ka_tv * exp(omega_ka);
#'   Vd = Vd_tv * exp(omega_Vd);
#'   CL = CL_tv * exp(omega_CL);
#'
#'   Cc = Ac / Vd;
#'
#'   Ad(0) = 0;
#'   Ac(0) = 0;
#'
#'   d/dt(Ad) = -ka * Ad;
#'   d/dt(Ac) = ka * Ad - CL * Cc;
#'
#'   Cc_ResErr = Cc * (1 + Cc_b);
#' })
#'
#' et_test <- tibble::tribble(
#'   ~ID, ~TIME, ~EVID, ~CMT, ~AMT,
#'   1, 0, 1, 1, 10,
#'   2, 0, 1, 1, 50
#' )
#'
#' stimes_test <- seq(0, 24, 0.1)
#' output_test <- c("Ac", "Ad", "Cc", "Cc_ResErr")
#' theta_test <- c(ka_pop = log(0.1), Vd_pop = log(10), CL_pop = log(0.5))
#'
#' thetamat_test <- matrix(c(0.05, 0.01, 0, 0.01, 0.05, 0, 0, 0, 0.05), nrow = 3, byrow = TRUE)
#' rownames(thetamat_test) <- colnames(thetamat_test) <- c("ka_pop", "Vd_pop", "CL_pop")
#'
#' omega_test <- matrix(c(0.2, 0.05, 0, 0.05, 0.2, 0, 0, 0, 0.2), nrow = 3, byrow = TRUE)
#' rownames(omega_test) <- colnames(omega_test) <- c("omega_ka", "omega_Vd", "omega_CL")
#'
#' sigma_test <- matrix(0.1); rownames(sigma_test) <- colnames(sigma_test) <- "Cc_b"
#'
#' # No variability
#' sim1 <- sg_sim(model = mod, et = et_test, stimes = stimes_test, theta = theta_test,
#'                outputs = output_test)
#'
#' # Population uncertainty, single scenario (ID)
#' sim2a <- sg_sim(model = mod, et = dplyr::filter(et_test, ID == 1),
#'                 stimes = stimes_test, theta = theta_test, thetamat = thetamat_test,
#'                 npop = 10)
#'
#' # Population uncertainty, multiple IDs, byID = TRUE (populations are replicated
#' # for each scenario (ID)),
#' # shared = TRUE (the same populations are used for each scenario (ID))
#' sim2b <- sg_sim(model = mod, et = et_test, stimes = stimes_test, theta = theta_test,
#'                thetamat = thetamat_test, npop = 10, byID = TRUE, shared = TRUE)
#'
#' # Population uncertainty, multiple IDs, byID = FALSE (one solve over full event table;
#' # npop must equal n(ID); ID - virtual subject)
#' sim2c <- sg_sim(model = mod, et = et_test, stimes = stimes_test, theta = theta_test,
#'                thetamat = thetamat_test, npop = nrow(et_test), byID = FALSE)
#'
#'
#' # Between-subject variability (BSV), multiple IDs, byID = TRUE (subjects are
#' # replicated for each scenario (ID)), shared = FALSE (separate subjects per ID)
#' sim3a <- sg_sim(model = mod, et = et_test, stimes = stimes_test, theta = theta_test,
#'                 omega = omega_test, nsub = 10, byID = TRUE, shared = FALSE)
#'
#'
#' # Between-subject variability (BSV), byID = FALSE (one solve over full event table;
#' # nsub must equal n(ID); ID - virtual subject)
#' sim3b <- sg_sim(model = mod, et = et_test, stimes = stimes_test, theta = theta_test,
#'                 omega = omega_test, nsub = nrow(et_test), byID = FALSE)
#'
#'
#' # BSV + uncertainty: byID = TRUE (populations/subjects are replicated for each
#' # scenario (ID)), byPOP = TRUE (subjects are replicated for each population),
#' # shared = FALSE (separate populations/subjects per ID)
#' sim4a <- sg_sim(model = mod, et = et_test, stimes = stimes_test, theta = theta_test,
#'                 thetamat = thetamat_test, npop = 10, omega = omega_test, nsub = 5,
#'                 byID = TRUE, byPOP = TRUE, shared = FALSE)
#'
#' # BSV + uncertainty: byID = TRUE (populations/subjects are replicated for each
#' # scenario (ID)), byPOP = FALSE (subjects are not replicated between population;
#' # npop = nsub)
#' sim4b <- sg_sim(model = mod, et = et_test, stimes = stimes_test, theta = theta_test,
#'                 thetamat = thetamat_test, npop = nrow(et_test), omega = omega_test,
#'                 nsub = nrow(et_test), byID = TRUE, byPOP = FALSE)
#'
#' # BSV + uncertainty: byID = FALSE (one solve over full event table; nsub and/or
#' # npop must equal n(ID); ID - virtual subject), byPOP = TRUE (subjects are
#' # replicated for each population), shared = FALSE (separate populations/subjects
#' # per ID)
#' sim4c <- sg_sim(model = mod, et = et_test, stimes = stimes_test, theta = theta_test,
#'                 thetamat = thetamat_test, npop = 10, omega = omega_test, nsub = 5,
#'                 byID = FALSE, byPOP = TRUE, shared = FALSE)
#'
#' # BSV + uncertainty: byID = FALSE (one solve over full event table; nsub and/or
#' # npop must equal n(ID); ID - virtual subject), byPOP = FALSE (subjects are not
#' # replicated between population; npop = nsub)
#' sim4d <- sg_sim(model = mod, et = et_test, stimes = stimes_test, theta = theta_test,
#'                 thetamat = thetamat_test, npop = nrow(et_test), omega = omega_test,
#'                 nsub = nrow(et_test), byID = FALSE, byPOP = FALSE)
#'
#' # Parameter override in event table
#' sim4 <- sg_sim(model = mod, et = dplyr::filter(et_test, ID == 1) %>%
#'                                  dplyr::mutate(Vd_pop = 2),
#'                stimes = stimes_test, theta = theta_test, thetamat = thetamat_test,
#'                npop = 1)
#'
#'
#' # Covariate (WTBL): pass covs and keep so WTBL is used and returned
#' theta_cov <- c(theta_test, beta_WTBL_Vd_pop = 0.5)
#' sim5 <- sg_sim(model = mod_cov, et = dplyr::filter(et_test, ID == 1) %>%
#'                                      dplyr::mutate(WTBL = 70),
#'                stimes = stimes_test, theta = theta_cov, covs = c("WTBL"),
#'                keep = c("WTBL"),
#'                thetamat = thetamat_test, npop = 5)
#' }
#' @import rxode2
#' @importFrom purrr map_dfr
#' @import dplyr
#' @import tidyr
#' @export
sg_sim <- function(model, et,  stimes = NULL, outputs = NULL, theta = NULL,
                   omega = NULL, sigma = NULL, thetamat = NULL, covs = NULL,
                   npop = 1, nsub = 1, aggr = NULL, addcov = TRUE, keep = NULL,
                   scale = NULL, covint = "locf", inits = NULL,
                   byID = NULL, byPOP = NULL, shared = NULL,
                   ncores = 1, atol = 1e-8, rtol = 1e-6, fpath_i = NULL, ...){
  if (!is.null(fpath_i)) {
    smrg_obj <- read_smrg_obj(fpath_i)

    if (is.null(stimes) & !is.null(smrg_obj$SDTAB)) {
      stimes <- smrg_obj$SDTAB %>% pull(TIME) %>% unique()
    }
    if (is.null(theta) & !is.null(smrg_obj$SUMTAB)) {
      theta <- smrg_obj[["SUMTAB"]][grepl("_pop", smrg_obj[["SUMTAB"]][["PAR"]]), ] %>%
        select(PAR, VALUE) %>% deframe()
    }
    if (is.null(thetamat) & !is.null(smrg_obj$COVMAT)) {
      thetamat <- smrg_obj[["COVMAT"]]
      if (nrow(thetamat) > ncol(thetamat)) {
        thetamat <- thetamat[1:ncol(thetamat), ]
      }
    }
    if (is.null(omega) & !is.null(smrg_obj$OMEGAMAT)) {
      omega <- smrg_obj[["OMEGAMAT"]]
      rownames(omega) <- colnames(omega)
    }
    if (is.null(sigma) & !is.null(smrg_obj[["SIGMAMAT"]])) {
      sigma <- smrg_obj[["SIGMAMAT"]]
      rownames(sigma) <- colnames(sigma)
    }
  }
  # ── id column: et may use "id" or "ID" ─────────────────────────────────
  et_id_col <- intersect(c("id", "ID"), colnames(et))[1]

  # ── auto-detect per-ID parameters from model + et columns ────────────
  et_standard <- c("id", "ID", "time", "TIME", "evid", "EVID", "cmt", "CMT",
                   "amt", "AMT", "addl", "ADDL", "ii", "II", "dur", "DUR",
                   "rate", "RATE", "ss", "SS", "mdv", "MDV", "tinf", "TINF",
                   "adm", "ADM")
  mod_pars    <- rxode2::rxModelVars(model)$params
  et_par_cols <- setdiff(
    intersect(mod_pars, colnames(et)),
    c(et_standard, covs, keep)
  )
  has_et_pars <- length(et_par_cols) > 0
  et_pars     <- NULL
  if (has_et_pars) {
    et_pars <- et %>%
      dplyr::select(dplyr::any_of(c("ID", "id")), dplyr::all_of(et_par_cols)) %>%
      dplyr::distinct()
    if ("ID" %in% colnames(et_pars) && !"id" %in% colnames(et_pars)) {
      et_pars <- et_pars %>% dplyr::rename(id = ID)
    }
  }

  # ── drop et-fixed parameters from variability matrices ─────────────
  if (has_et_pars) {
    if (!is.null(thetamat)) {
      drop <- intersect(et_par_cols, colnames(thetamat))
      if (length(drop) > 0) {
        idx <- !colnames(thetamat) %in% drop
        thetamat <- thetamat[idx, idx, drop = FALSE]
        if (ncol(thetamat) == 0) thetamat <- NULL
      }
    }
    if (!is.null(omega)) {
      drop <- intersect(et_par_cols, colnames(omega))
      if (length(drop) > 0) {
        idx <- !colnames(omega) %in% drop
        omega <- omega[idx, idx, drop = FALSE]
        if (ncol(omega) == 0) omega <- NULL
      }
    }
    if (!is.null(sigma)) {
      drop <- intersect(et_par_cols, colnames(sigma))
      if (length(drop) > 0) {
        idx <- !colnames(sigma) %in% drop
        sigma <- sigma[idx, idx, drop = FALSE]
        if (ncol(sigma) == 0) sigma <- NULL
      }
    }
  }

  # ── event table compilation ──────────────────────────────────────────
  et_i_m <- et %>% rxode2::et()
  if (!is.null(stimes)) et_i_m <- et_i_m %>% rxode2::add.sampling(stimes)

  et_i_cov <- NULL
  if (!is.null(covs) && !is.na(et_id_col)) {
    if (!is.null(stimes)) {
      et_i_cov <- et %>%
        dplyr::select(dplyr::all_of(et_id_col), dplyr::all_of(covs)) %>%
        dplyr::distinct()
      if (et_id_col == "ID") et_i_cov <- et_i_cov %>% dplyr::rename(id = ID)
      et_i_m   <- et_i_m %>%
        dplyr::as_tibble() %>%
        dplyr::left_join(et_i_cov, by = "id")
    } else {
      et_i_cov <- et %>%
        dplyr::select(dplyr::all_of(et_id_col), dplyr::all_of(covs))
      if (et_id_col == "ID") et_i_cov <- et_i_cov %>% dplyr::rename(id = ID)
      et_i_m   <- et_i_m %>%
        dplyr::bind_cols(et_i_cov %>% dplyr::select(-dplyr::any_of("id")))
    }
  }

  # ── internal solve helpers ───────────────────────────────────────────

  get_theta <- function(k = NULL) {
    if (!has_et_pars || is.null(k)) return(theta)
    theta_k <- et_pars %>%
      dplyr::filter(id == k) %>%
      dplyr::select(-dplyr::any_of("id")) %>%
      as.list() %>%
      unlist()
    if (!is.null(theta)) {
      theta_k <- c(theta_k, theta[!names(theta) %in% names(theta_k)])
    }
    theta_k
  }

  run_solve <- function(events, nStud, nSub, params = theta, ...) {
    rxode2::rxSolve(model, events = events, params = params, omega = omega,
                    inits = inits, sigma = sigma, thetaMat = thetamat,
                    nStud = nStud, nSub = nSub,
                    covsInterpolation = covint, cores = ncores,
                    simVariability = TRUE, atol = atol, rtol = rtol, ...)
  }

  solve_by_id <- function(nStud, nSub, post_fn = NULL, verbose = FALSE, ...) {
    et_tbl <- et_i_m %>% dplyr::as_tibble()
    id_col_et <- intersect(c("id", "ID"), colnames(et_tbl))[1]
    if (is.na(id_col_et)) {
      id_col_et <- "id"
      if (!id_col_et %in% colnames(et_tbl)) et_tbl <- dplyr::mutate(et_tbl, id = 1L)
    }
    ids <- unique(et_tbl[[id_col_et]])
    purrr::map_dfr(ids, function(k) {
      if (verbose) message("id = ", k, " (out of ", max(et_tbl[[id_col_et]]), ")")
      et_k <- et_tbl %>%
        dplyr::filter(.data[[id_col_et]] == k) %>%
        dplyr::mutate(!!id_col_et := 1L)
      theta_k <- get_theta(k)
      sim_k <- run_solve(et_k, nStud, nSub, params = theta_k, ...)
      if ("ID" %in% colnames(sim_k) && !"id" %in% colnames(sim_k)) {
        sim_k <- dplyr::rename(sim_k, id = ID)
      }
      sim_k <- dplyr::mutate(sim_k, id = k)
      if (!is.null(post_fn)) sim_k <- post_fn(sim_k, k)
      sim_k
    })
  }

  do_solve <- function(nStud, nSub, ...) {
    if (has_et_pars) solve_by_id(nStud, nSub, ...)
    else run_solve(et_i_m, nStud, nSub, ...)
  }

  # ── resolve strategy ─────────────────────────────────────────────────
  has_thetamat <- !is.null(thetamat)
  has_omega    <- !is.null(omega)
  by_id        <- is.null(byID)  || isTRUE(byID)
  by_pop       <- is.null(byPOP) || isTRUE(byPOP)
  is_shared    <- is.null(shared) || isTRUE(shared)
  n_ids        <- if (!is.na(et_id_col)) length(unique(et[[et_id_col]])) else 1L

  add_popn_ceiling <- FALSE
  sim_i <- NULL

  if (!has_thetamat && !has_omega) {
    sim_i <- do_solve(npop, nsub, ...)

  } else if (has_thetamat && !has_omega) {
    add_popn_ceiling <- TRUE
    if (by_id && is_shared) {
      sim_i <- do_solve(npop, nsub, ...)
    } else if (by_id && !is_shared) {
      sim_i <- solve_by_id(npop, nsub, ...)
    } else {
      sim_i <- run_solve(et_i_m, npop, nsub, ...)
    }

  } else if (!has_thetamat && has_omega) {
    if (by_id) {
      sim_i <- solve_by_id(npop, nsub, ...)
    } else {
      if (n_ids != nsub) {
        message("!Critical condition: nsub = n(ID); otherwise byID = TRUE")
        return(NULL)
      }
      sim_i <- do_solve(npop, nsub, ...)
    }

  } else {
    if (by_id) {
      if (by_pop) {
        add_popn_ceiling <- TRUE
        sim_i <- solve_by_id(npop, nsub, ...)
      } else {
        if (npop != nsub) {
          message("!Critical condition: nsub = npop; otherwise byPOP = TRUE")
          return(NULL)
        }
        sim_i <- do_solve(npop, nsub, ...)
      }
    } else {
      if (npop != n_ids && nsub != n_ids) {
        message("!Critical condition: nsub = n(ID) OR npop = n(ID); otherwise byID = TRUE")
        return(NULL)
      }
      if (by_pop) {
        if (is_shared) {
          sim_i <- solve_by_id(npop, nsub, ...)
        } else if (npop != n_ids) {
          sim_i <- solve_by_id(npop, 1, post_fn = function(sim_k, k) {
            sim_k %>% dplyr::mutate(sim.id = k) %>%
              dplyr::group_by(time) %>%
              dplyr::mutate(POPN = seq_len(npop)) %>%
              dplyr::ungroup()
          }, ...)
        } else {
          sim_i <- solve_by_id(1, nsub, post_fn = function(sim_k, k) {
            sim_k %>% dplyr::mutate(POPN = k) %>%
              dplyr::group_by(time) %>%
              dplyr::mutate(sim.id = seq_len(nsub)) %>%
              dplyr::ungroup()
          }, ...)
        }
      } else {
        if (npop != nsub) {
          message("!Critical condition: nsub = npop; otherwise byPOP = TRUE")
          return(NULL)
        }
        sim_i <- solve_by_id(1, 1, post_fn = function(sim_k, k) {
          sim_k %>% dplyr::mutate(POPN = k, sim.id = k)
        }, ...)
      }
    }
  }

  if (add_popn_ceiling && "sim.id" %in% colnames(sim_i)) {
    sim_i <- sim_i %>% dplyr::mutate(POPN = ceiling(sim.id / nsub))
  }

  # ── post-processing: ID / POPN / sim.id by scenario ───────────────────
  # No uncertainty & no variability → only ID
  # Only uncertainty (thetamat)     → POPN and ID
  # Only variability (omega)        → sim.id and ID
  # Uncertainty + variability      → POPN, sim.id, ID
  sim_i_ind <- sim_i %>% dplyr::as_tibble()
  if ("ID" %in% colnames(sim_i_ind) && !"id" %in% colnames(sim_i_ind)) {
    sim_i_ind <- dplyr::rename(sim_i_ind, id = ID)
  }
  if (!"id" %in% colnames(sim_i_ind)) sim_i_ind <- dplyr::mutate(sim_i_ind, id = 1)
  if (has_thetamat && !"POPN" %in% colnames(sim_i_ind)) {
    sim_i_ind <- dplyr::mutate(sim_i_ind, POPN = 1)
  }
  if (has_omega && !"sim.id" %in% colnames(sim_i_ind)) {
    sim_i_ind <- dplyr::mutate(sim_i_ind, sim.id = 1)
  }
  if (!has_thetamat && "POPN" %in% colnames(sim_i_ind)) {
    sim_i_ind <- dplyr::select(sim_i_ind, -POPN)
  }
  if (!has_omega && "sim.id" %in% colnames(sim_i_ind)) {
    sim_i_ind <- dplyr::select(sim_i_ind, -sim.id)
  }
  id_cols <- c("id", "time", if (has_thetamat) "POPN", if (has_omega) "sim.id")
  sim_i_ind <- sim_i_ind %>%
    dplyr::mutate(dplyr::across(-dplyr::any_of(id_cols), as.numeric)) %>%
    tidyr::gather("VAR", "VALUE", -dplyr::any_of(id_cols))
  if (!is.null(outputs)) sim_i_ind <- sim_i_ind %>% dplyr::filter(VAR %in% outputs)

  sim_i_out <- sim_i_ind

  if (is.null(aggr)) {
    #if (interactive()) message("No aggregation applied;")
  } else if (aggr == "ID") {
    sim_i_out <- sim_i_ind %>%
      dplyr::group_by(id, time, VAR) %>%
      dplyr::summarise_at(dplyr::vars(VALUE), funSum_sim) %>%
      dplyr::ungroup()
  } else if (aggr == "NPOP") {
    if ("POPN" %in% colnames(sim_i_ind)) {
      sim_i_out <- sim_i_ind %>%
        dplyr::group_by(POPN, time, VAR) %>%
        dplyr::summarise_at(dplyr::vars(VALUE), funSum_sim) %>%
        dplyr::ungroup()
    } else {
      sim_i_out <- sim_i_ind %>%
        dplyr::group_by(time, VAR) %>%
        dplyr::summarise_at(dplyr::vars(VALUE), funSum_sim) %>%
        dplyr::ungroup()
    }
  } else if (aggr == "TIME") {
    sim_i_out <- sim_i_ind %>%
      dplyr::group_by(time, VAR) %>%
      dplyr::summarise_at(dplyr::vars(VALUE), funSum_sim) %>%
      dplyr::ungroup()
  }

  if (!is.null(keep) && !is.null(sim_i_out) && "id" %in% colnames(sim_i_out) && !is.na(et_id_col)) {
    keep_et <- et %>%
      dplyr::select(dplyr::all_of(et_id_col), dplyr::all_of(keep)) %>%
      dplyr::distinct()
    if (et_id_col == "ID") keep_et <- keep_et %>% dplyr::rename(id = ID)
    suppressMessages(
      sim_i_out <- dplyr::left_join(sim_i_out, keep_et, by = "id", relationship = "many-to-many")
    )
  }

  if (addcov && !is.null(covs) && !is.null(sim_i_out) && "id" %in% colnames(sim_i_out) && !is.null(et_i_cov)) {
    suppressMessages(
      sim_i_out <- dplyr::left_join(
        sim_i_out,
        dplyr::distinct(dplyr::select(et_i_cov, id, dplyr::all_of(covs))),
        relationship = "many-to-many"
      )
    )
  }

  lookup <- c(TIME = "time", ID = "id")
  req_cols <- c(
    if (has_thetamat) "POPN",
    if (has_omega) "sim.id",
    "ID", "TIME", "VAR", "VALUE",
    "mean", "median", "min", "max", "sd",
    "P025", "P05", "P10", "P15", "P25", "P75", "P85", "P90", "P95", "P975",
    "geom_mean", "CV", covs, keep
  )
  sim_i_out <- sim_i_out %>%
    dplyr::rename(dplyr::any_of(lookup)) %>%
    dplyr::select(dplyr::any_of(req_cols))
  sim_i_out
}

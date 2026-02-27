#' Global sensitivity analysis via PRCC or eFAST (sg_sim-based)
#'
#' @description
#' Performs global sensitivity analysis using either PRCC (Partial Rank
#' Correlation Coefficients with LHS sampling) or eFAST (extended Fourier
#' Amplitude Sensitivity Test). For each sampled parameter set, the function
#' runs simulations via `sg_sim()`, reduces trajectories to scalar summary
#' statistics, and computes sensitivity indices for each output–statistic pair.
#'
#' @param method Character string. `"PRCC"` or `"eFAST"`.
#' @param model Model object passed to `sg_sim()`.
#' @param params Character vector of parameter names to vary.
#' @param par_bounds Tibble/data.frame with columns `PAR`, `LB`, `UB`.
#' @param n_sim Integer. Number of samples (LHS size for PRCC, base frequency size for eFAST).
#' @param stimes Numeric vector of simulation times.
#' @param output Character vector of outputs to keep. Passed to `sg_sim()`.
#' @param stat_comp Character vector of summary statistics to compute.
#'        Supported internally: `"mean","median","min","max","sd","cmax","SS"`.
#' @param et Event table passed to `sg_sim()`.
#' @param theta Named numeric vector of baseline parameters. Default NULL.
#'        Parameters listed in `params` are replaced by sampled values.
#' @param cov Covariates passed to `sg_sim()`. Default NULL.
#'
#' @return List with:
#' \describe{
#'   \item{result}{
#'     Sensitivity results.
#'
#'     If `method = "PRCC"`:
#'     Tibble with columns:
#'     \itemize{
#'       \item `PAR` — parameter
#'       \item `VAR` — model output variable
#'       \item `STAT` — summary statistic
#'       \item `estimate` — PRCC estimate (Spearman partial correlation)
#'       \item `p.value` — associated p-value
#'     }
#'
#'     If `method = "eFAST"`:
#'     Tibble with columns:
#'     \itemize{
#'       \item `PAR` — parameter
#'       \item `VAR` — model output variable
#'       \item `STAT` — summary statistic
#'       \item `TYPE` — `"First order"` or `"Total order"`
#'       \item `VALUE` — sensitivity index
#'       \item `METHOD` — `"eFAST"`
#'     }
#'   }
#'
#'   \item{summary}{
#'     Tibble of scalar summaries used for sensitivity computation:
#'     `sim.id`, `VAR`, `STAT`, `value`.
#'   }
#'
#'   \item{design}{
#'     Data.frame of sampled parameter values (real scale) with `sim.id`.
#'     For PRCC this corresponds to LHS samples.
#'     For eFAST this corresponds to the FAST design matrix.
#'   }
#'
#'   \item{bounds}{
#'     Parameter bounds used in the analysis.
#'   }
#' }
#'
#' @details
#' PRCC uses Latin Hypercube Sampling (LHS) combined with partial rank
#' correlations (Spearman-based) to quantify monotonic relationships between
#' parameters and scalar model outputs.
#'
#' eFAST performs variance-based global sensitivity analysis using spectral
#' decomposition to estimate first-order and total-order sensitivity indices.
#'
#' The function assumes that `sg_sim()` returns a data frame containing
#' at least the columns: `TIME`, `VAR`, `VALUE`, and `sim.id`.
#'
#' @import dplyr tidyr purrr lhs sensitivity ppcor
#' @export

# library(dplyr)
# library(tidyr)
# library(purrr)
# library(tibble)
library(lhs)
library(ppcor)




#### FUL FUN ####

sg_globalsens_sim <- function(method = c("PRCC", "eFAST"),
                           model,
                           params,
                           par_bounds,
                           n_sim,
                           stimes,
                           output,
                           stat_comp = c("mean", "min", "max"),
                           et,
                           theta = NULL,
                           cov = NULL) {

  method <- match.arg(method)


  req_cols <- c("PAR", "LB", "UB")
  if (!all(req_cols %in% colnames(par_bounds))) {
    stop("par_bounds must contain columns: PAR, LB, UB")
  }

  params <- as.character(params)

  bounds <- par_bounds %>%
    dplyr::filter(PAR %in% params) %>%
    dplyr::slice(match(params, PAR))

  if (nrow(bounds) != length(params))
    stop("Some parameters missing in par_bounds")


  # DESIGN GENERATION


  if (method == "PRCC") {

    X_norm <- lhs::randomLHS(n_sim, length(params)) %>% as.data.frame()
    colnames(X_norm) <- params
    X_norm$sim.id <- seq_len(n_sim)

    X_real <- purrr::imap_dfc(X_norm %>% dplyr::select(-sim.id), function(x, nm) {
      LB <- bounds$LB[match(nm, bounds$PAR)]
      UB <- bounds$UB[match(nm, bounds$PAR)]
      (UB - LB) * x + LB
    }) %>%
      dplyr::mutate(sim.id = X_norm$sim.id)

  } else {

    # eFAST design
    q <- rep("qunif", length(params))
    q.arg <- lapply(seq_len(nrow(bounds)), function(i) {
      list(min = bounds$LB[i], max = bounds$UB[i])
    })

    fast_obj <- sensitivity::fast99(
      model = NULL,
      factors = params,
      n = n_sim,
      q = q,
      q.arg = q.arg
    )

    X_real <- as.data.frame(fast_obj$X)
    X_real$sim.id <- seq_len(nrow(X_real))
    X_norm <- NULL
  }


  # SIMULATIONS


  sim_all <- purrr::map_dfr(seq_len(nrow(X_real)), function(i) {

    theta_i <- if (is.null(theta)) c() else theta
    theta_i[params] <- as.numeric(X_real[i, params])

    sg_sim(
      model = model,
      stimes = stimes,
      output = output,
      theta = theta_i,
      covs = cov,
      et = et
    ) %>%
      dplyr::mutate(sim.id = i)
  })

  if (!all(c("TIME","VAR","VALUE","sim.id") %in% colnames(sim_all)))
    stop("sg_sim() must return TIME, VAR, VALUE, sim.id")


  # SUMMARY STATISTICS

  summary_df <- sim_all %>%
    dplyr::group_by(sim.id, VAR) %>%
    dplyr::summarise(
      mean = mean(VALUE, na.rm = TRUE),
      median = median(VALUE, na.rm = TRUE),
      min = min(VALUE, na.rm = TRUE),
      max = max(VALUE, na.rm = TRUE),
      sd = sd(VALUE, na.rm = TRUE),
      cmax = max(VALUE, na.rm = TRUE),
      SS = VALUE[TIME == max(stimes)][1],
      .groups = "drop"
    ) %>%
    dplyr::select(sim.id, VAR, dplyr::all_of(stat_comp)) %>%
    tidyr::pivot_longer(
      cols = -c(sim.id, VAR),
      names_to = "STAT",
      values_to = "value"
    )


  # METHOD-SPECIFIC PART

  if (method == "PRCC") {

    prcc_df <- summary_df %>%
      dplyr::left_join(X_real, by = "sim.id")

    result <- prcc_df %>%
      dplyr::group_by(.data$VAR, .data$STAT) %>%
      dplyr::group_modify(function(dat, key) {

        dat2 <- dat %>%
          dplyr::select(.data$value, dplyr::all_of(params)) %>%
          dplyr::mutate(dplyr::across(dplyr::everything(), as.numeric)) %>%
          tidyr::drop_na()


        if (nrow(dat2) < (length(params) + 2) || stats::sd(dat2$value) == 0) {
          return(tibble::tibble(
            PAR = params,
            estimate = NA_real_,
            p.value = NA_real_
          ))
        }


        good_par <- params[sapply(dat2[, params, drop = FALSE], stats::sd) > 0]

        if (length(good_par) == 0) {
          return(tibble::tibble(
            PAR = params,
            estimate = NA_real_,
            p.value = NA_real_
          ))
        }

        # ---- pcor
        pc <- suppressWarnings(ppcor::pcor(dat2[, c("value", good_par), drop = FALSE],
                                           method = "spearman"))


        if (is.null(dimnames(pc$estimate))) {
          dimnames(pc$estimate) <- list(colnames(dat2[, c("value", good_par)]),
                                        colnames(dat2[, c("value", good_par)]))
        }
        if (is.null(dimnames(pc$p.value))) {
          dimnames(pc$p.value) <- list(colnames(dat2[, c("value", good_par)]),
                                       colnames(dat2[, c("value", good_par)]))
        }

        est <- pc$estimate["value", good_par, drop = TRUE]
        pvl <- pc$p.value["value", good_par, drop = TRUE]


        out <- tibble::tibble(
          PAR = params,
          estimate = NA_real_,
          p.value = NA_real_
        )

        out$estimate[out$PAR %in% good_par] <- as.numeric(est)
        out$p.value[out$PAR %in% good_par] <- as.numeric(pvl)

        out
      }) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(PAR = factor(.data$PAR, levels = params))

  }
  else {

    result <- purrr::map_dfr(unique(summary_df$VAR), function(v) {
      purrr::map_dfr(unique(summary_df$STAT), function(s) {

        y <- summary_df %>%
          dplyr::filter(VAR == v, STAT == s) %>%
          dplyr::pull(value)

        if (sd(y) == 0)
          return(NULL)

        fast_obj <- sensitivity::fast99(
          model = NULL,
          factors = params,
          n = n_sim,
          q = rep("qunif", length(params)),
          q.arg = lapply(seq_len(nrow(bounds)),
                         function(i) list(min = bounds$LB[i],
                                          max = bounds$UB[i]))
        )

        fast_obj$y <- y
        fast_obj <- sensitivity::tell(fast_obj)

        tibble::tibble(
          PAR = params,
          VAR = v,
          STAT = s,
          `First order` = fast_obj$D1 / fast_obj$V,
          `Total order` = 1 - fast_obj$Dt / fast_obj$V
        ) %>%
          tidyr::pivot_longer(
            cols = c(`First order`, `Total order`),
            names_to = "TYPE",
            values_to = "VALUE"
          ) %>%
          dplyr::mutate(METHOD = "eFAST")
      })
    })
  }

  return(list(
    result = result,
    summary = summary_df,
    design = X_real,
    bounds = bounds
  ))
}
mod_ex <- RxODE({
  # Doses in mg
  # Time in hours

  ### Parameter values
  # Typical
  POPCL = 5;
  POPVC = 180;
  POPQ = 7;
  POPVP = 52;
  POPKTR = 6;
  FBIOPAR = 1;

  # Covariate coefficients
  POPIGFRCOV = 0.8;
  POPPATCOVCLGR1 = 0.85;
  POPPATCOVCLGR2 = 0.85;
  POPPATCOVCLGR3 = 0.85;
  POPPATCOVCLGR4 = 0.85;

  # Random effects
  PPVCL = 0;
  PPVVC = 0;
  PPVKTR = 0;

  # Residual error
  RUV = 0;

  ### Covariates
  PATCOVCL = 1
  if(POPN == 1){PATCOVCL = POPPATCOVCLGR1}
  if(POPN == 2){PATCOVCL = POPPATCOVCLGR2}
  if(POPN == 3){PATCOVCL = POPPATCOVCLGR3}
  if(POPN == 4){PATCOVCL = POPPATCOVCLGR4}
  IGFRCOV = (IGFR/112)^POPIGFRCOV

  ### Parameters
  CL = POPCL * IGFRCOV * PATCOVCL * exp(PPVCL);
  VC = POPVC * exp(PPVVC);
  Q = POPQ;
  VP = POPVP;
  KTR = POPKTR * exp(PPVKTR);

  ### Explicit functions
  Cc = Ac/VC;                 # nmol/L
  Cp = Ap/VP;                 # nmol/L

  ### Initial conditions
  At1(0) = 0;         # mg
  At2(0) = 0;         # mg
  At3(0) = 0;         # mg
  At4(0) = 0;         # mg
  At5(0) = 0;         # mg
  At6(0) = 0;         # mg
  Ad(0) = 0;          # mg
  Ac(0) = 0;          # mg
  Ap(0) = 0;          # mg

  ### ODEs
  d/dt(At1) = - KTR*At1;
  d/dt(At2) = KTR*At1 - KTR*At2;
  d/dt(At3) = KTR*At2 - KTR*At3;
  d/dt(At4) = KTR*At3 - KTR*At4;
  d/dt(At5) = KTR*At4 - KTR*At5;
  d/dt(At6) = KTR*At5 - KTR*At6;
  d/dt(Ad) = KTR*At6 - KTR*Ad;
  d/dt(Ac) = KTR*Ad - CL*Cc - Q*(Cc - Cp);
  d/dt(Ap) = Q*(Cc - Cp);

  FBIO = FBIOPAR
  f(At1) = FBIO*1000000/505;
  CHECKRUV = RUV;
  Cc_ResErr = Cc + RUV*Cc;
})


#####--------------- Set up event tables ---------------#####
et_base <- tribble(
  ~id, ~time, ~evid, ~cmt, ~amt, ~addl, ~ii, ~IGFR, ~POPN,
  1,   0,     1,     1,    10,   2,     24,  112,   1
)
inits <- rxInits(mod_ex)
par_bounds <- tibble::tibble(
  PAR = c("POPCL", "POPVC"),
  LB  = inits[c("POPCL", "POPVC")] * (1 - 0.9),
  UB  = inits[c("POPCL", "POPVC")] * (1 + 0.9)
)
res_efast <- sg_globalsens_sim(
  method = c("eFAST"),
  model = mod_ex,
  params = c("POPCL"),
  par_bounds = par_bounds,
  n_sim = 100,
  stimes = seq(0, 168, 10),
  output = "Cc",
  cov = c("IGFR", "POPN"),
  et = et_base,
  stat_comp = c("mean")
)

res_prcc <- sg_globalsens_sim(
  method = c("PRCC"),
  model = mod_ex,
  params = c("POPCL","POPVC"),
  par_bounds = par_bounds,
  n_sim = 100,
  stimes = seq(0, 168, 10),
  output = "Cc",
  cov = c("IGFR", "POPN"),
  et = et_base,
  stat_comp = c("mean")
)

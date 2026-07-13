## Author: Alina Melnikova
## First created: 2026-01-28
## Description: Simulation for covariate sensitivity analysis
## Keywords: covariates, simulations, sensitivity

###-----Functions----####
#' @noRd
funSum_av <- list(mean   = ~mean(.),
                  median   = ~median(., na.rm = TRUE))

#' @noRd
funSum_exp <- list(`Cavg`   = ~mean(., na.rm = TRUE),
                   `Cmin`    = ~min(., na.rm = TRUE),
                   `Cmax`    = ~max(., na.rm = TRUE))

#' @noRd
calc_auc_trapz <- function(time, value){
  ok <- !(is.na(time) | is.na(value))
  t <- as.numeric(time[ok])
  y <- as.numeric(value[ok])
  if(length(t) < 2){
    return(NA_real_)
  }
  ord <- order(t)
  t <- t[ord]
  y <- y[ord]
  sum(diff(t) * (head(y, -1) + tail(y, -1)) / 2)
}

#' @noRd
fun_EtCC <- function(et_base_i, cc_ds_i, cat = FALSE){
  et_scov_i <- unique(cc_ds_i$COV) %>% map(function(n){
    cc_ds_n <- cc_ds_i %>% filter(COV == n)

    et_scov_n <- seq(nrow(cc_ds_n)) %>% map_dfr(function(m){
      row_m <- cc_ds_n %>% filter(row_number() == m)

      et_scov_m <- et_base_i %>%
        mutate_at(vars(all_of(row_m$COV)), function(k){k = row_m$COVVAL})
      if(!cat){
        et_scov_m <- et_scov_m %>% bind_cols(select(row_m, COV, PAR, BTR, KEY, COVVAL, BCOVVAL))
      } else {
        et_scov_m <- et_scov_m %>% bind_cols(select(row_m, COV, PAR, KEY, COVVAL, CATDES))
      }

      return(et_scov_m)

    }) %>% mutate(id = row_number())

    return(et_scov_n)
  })
  names(et_scov_i) <- unique(cc_ds_i$COV)
  return(et_scov_i)
}

#' @noRd
fun_CovSens <- function(et_sim_i, cat = FALSE, expos = FALSE, covs_i = NULL, nsim = 100, stime_exp = NULL,
                        mod_fin_i,
                        theta_i, omega_i, thetamat_i, nice_names_i,
                        quantiles = c(0.1,0.9),
                        var_exp = "Cc", aggr = c("min", "max", "mean", "auc")) #nsim=1000
{


  #Test
  # et_sim_i = ets_cc
  # mod_fin_i = model
  # covs_i = nice_names$COV
  # expos = T
  # stime_exp = stimes
  # cat = F; nsim = 5;
  # var_exp = "Cc"
  # theta_i = par_fin_tv
  # thetamat_i = m_theta_norm_pop
  # omega_i = m_omega_full
  # aggr = c("min", "max", "mean", "auc")




  keep_i <- c("Regimen", "KEY", "COV", "COVVAL")
  if(!cat){
    keep_i <- c(keep_i, "BTR", "BCOVVAL")
  } else {
    keep_i <- c(keep_i, "CATDES")
  }

  sens_i <- et_sim_i %>% map_dfr(function(et_i){
    #Test
    #et_i <- et_sim_i[[1]]

    if(!expos){
      par_i <- unique(et_i$PAR)
      if(is.list(par_i)){ par_i <- par_i[[1]] }
      sim_i <- sg_sim(model = mod_fin_i, et = et_i, stimes = 0, outputs = par_i, theta = theta_i, omega = omega_i,
                      thetamat = thetamat_i, covs = covs_i, npop = nsim, keep = keep_i)
    } else {
      sim_raw <- sg_sim(mod_fin_i, et_i, stimes = stime_exp, outputs = var_exp, theta = theta_i, omega = omega_i,
                        thetamat = thetamat_i, covs = covs_i, npop = nsim,
                        keep = keep_i, ncores = max(1, parallel::detectCores()-2))
      #Test diagnostic plots
      # plot_sim_raw <- ggplot2::ggplot(
      #   sim_raw,
      #   ggplot2::aes(
      #     x = TIME,
      #     y = as.numeric(VALUE),
      #     color = factor(KEY),
      #     group = interaction(KEY, POPN)
      #   )
      # ) +
      #   ggplot2::geom_line(alpha = 0.7, linewidth = 0.6) +
      #   ggplot2::facet_wrap(~VAR, scales = "free_y") +
      #   ggplot2::labs(x = "Time", y = "Value", color = "Scenario") +
      #   ggplot2::theme_bw()
      # plot_sim_raw

      # message(sprintf("NA Cc rows: %d / %d (%.1f%%)",
      #                 sum(is.na(sim_raw$VALUE)), nrow(sim_raw),
      #                 100 * mean(is.na(sim_raw$VALUE))))

      #Warning / success message
      outputs_run <- if ("VAR" %in% names(sim_raw)) unique(as.character(sim_raw$VAR)) else as.character(var_exp)
      outputs_label <- paste(outputs_run, collapse = ", ")
      na_rows <- sum(is.na(sim_raw$VALUE))
      total_rows <- nrow(sim_raw)
      na_pct <- if (total_rows > 0) 100 * na_rows / total_rows else 0
      if (na_rows > 0) {
        message(sprintf(
          "Simulation completed for output(s): %s. NA rows detected: %d / %d (%.1f%%).",
          outputs_label, na_rows, total_rows, na_pct
        ))
      }
      ###

      # aggr_map <- c("mean" = "Cavg", "min" = "Cmin", "max" = "Cmax")
      # stopifnot("aggr must only contain 'mean', 'min', 'max'" = all(aggr %in% names(aggr_map)))
      # funSum_exp_i <- funSum_exp[unname(aggr_map[aggr])]


      aggr_map <- c("mean" = "Cavg", "min" = "Cmin", "max" = "Cmax", "auc" = "Cauc")
      stopifnot("aggr must only contain 'mean', 'min', 'max', 'auc'" = all(aggr %in% names(aggr_map)))
      metric_keep <- unname(aggr_map[aggr])

      sim_i <- sim_raw %>% mutate(VALUE = as.numeric(VALUE))

      # sim_i <- sim_i %>%
      # group_by_at(vars(all_of(c("POPN", "VAR", keep_i, covs_i)))) %>% summarise_at(vars(VALUE), funSum_exp_i) %>% ungroup() %>%
      #   gather("METRIC", "VALUE", -all_of(c("POPN", "VAR", keep_i, covs_i))) %>%
      #   mutate(VAR = paste(VAR, METRIC, sep = "_")) %>%
      #   select(-METRIC)

      sim_i <- sim_i %>%
      group_by_at(vars(all_of(c("POPN", "VAR", keep_i, covs_i)))) %>%
        summarise(
          Cavg = mean(VALUE, na.rm = TRUE),
          Cmin = min(VALUE, na.rm = TRUE),
          Cmax = max(VALUE, na.rm = TRUE),
          Cauc = calc_auc_trapz(TIME, VALUE),
          .groups = "drop"
        ) %>%
        select(all_of(c("POPN", "VAR", keep_i, covs_i, metric_keep))) %>%
        gather("METRIC", "VALUE", all_of(metric_keep)) %>%
        mutate(VAR = paste(VAR, METRIC, sep = "_")) %>%
        select(-METRIC)

    }

    sim_i <- sim_i %>% mutate(VALUE = as.numeric(VALUE))

    sim_i_ref <- sim_i %>% filter(KEY == "REF" | KEY == 1) %>% select(POPN, Regimen, VAR, REFVAL = VALUE)
    sim_i_ch <- sim_i %>% filter(KEY != "REF" & KEY != 1) %>% left_join(sim_i_ref, by = c("POPN", "Regimen", "VAR")) %>%
      mutate(PCH = 100*(VALUE - REFVAL)/REFVAL)


    out_i <- sim_i_ch %>% group_by_at(vars(all_of(c("VAR", keep_i, covs_i)))) %>% summarise_at(vars(PCH), funSum_sim) %>% ungroup()

    # Ensure BCOVVAL is added from et_i if not already present (for continuous covariates)
    if(!cat && "BCOVVAL" %in% names(et_i)){
      if(!"BCOVVAL" %in% names(out_i)){
        bcovval_lookup <- et_i %>% select(KEY, COV, COVVAL, BCOVVAL) %>% unique()
        out_i <- out_i %>% left_join(bcovval_lookup, by = c("KEY", "COV", "COVVAL"))
      }
    }
    out_i

  })  %>% left_join(nice_names_i, by = "COV")

  if(!cat){
    sens_out <- sens_i %>% select(NICEN, VAR, KEY, mean:P975, Regimen, COVVAL:BCOVVAL) %>%
      mutate(KEY = ifelse(KEY == "LP", paste0(quantiles[[1]]*100, "th perc."), paste0(quantiles[[2]]*100, "th perc.")), LAB = str_c(NICEN, "\n", KEY, " (", round(BCOVVAL, 1), ")"))
  } else {
    sens_out <- sens_i %>% select(NICEN, VAR, CATDES, mean:P975, Regimen) %>% mutate(LAB = str_c(NICEN, " (", CATDES, ")"))
  }
  sens_out <- sens_out %>% mutate_at(vars(mean:P975), function(p){p/100+1})
  return(sens_out)
}

#' @noRd
normalize_cont_cov_l <- function(cont_cov_l) {
  if (length(cont_cov_l) == 0) {
    return(cont_cov_l)
  }
  purrr::imap(cont_cov_l, function(cov_def, cov_name) {
    if (!is.list(cov_def)) {
      stop("Each element of 'cont_cov_l' must be a list; invalid entry: ", cov_name, ".")
    }
    if (!"par_vec" %in% names(cov_def) || is.null(cov_def$par_vec) || length(cov_def$par_vec) == 0) {
      stop("Each continuous covariate in 'cont_cov_l' must include a non-empty 'par_vec'; missing for: ", cov_name, ".")
    }
    list(
      UTNAME = if ("UTNAME" %in% names(cov_def)) cov_def$UTNAME else NULL,
      REF = if ("REF" %in% names(cov_def)) cov_def$REF else "median",
      NICENAME = if ("NICENAME" %in% names(cov_def)) cov_def$NICENAME else NULL,
      par_vec = cov_def$par_vec
    )
  })
}

#' @noRd
normalize_cat_cov_l <- function(cat_cov_l) {
  if (length(cat_cov_l) == 0) {
    return(cat_cov_l)
  }
  purrr::imap(cat_cov_l, function(cov_def, cov_name) {
    if (!is.list(cov_def)) {
      stop("Each element of 'cat_cov_l' must be a list; invalid entry: ", cov_name, ".")
    }
    if (!"par_vec" %in% names(cov_def) || is.null(cov_def$par_vec) || length(cov_def$par_vec) == 0) {
      stop("Each categorical covariate in 'cat_cov_l' must include a non-empty 'par_vec'; missing for: ", cov_name, ".")
    }
    list(
      REF = if ("REF" %in% names(cov_def)) cov_def$REF else NULL,
      NICENAME = if ("NICENAME" %in% names(cov_def)) cov_def$NICENAME else NULL,
      par_vec = cov_def$par_vec
    )
  })
}

#' @noRd
cov_nicename <- function(cov_def, cov_name) {
  nicen <- cov_def$NICENAME
  if (is.null(nicen) || is.na(nicen) || nicen == "") cov_name else nicen
}

#####

#' Covariate sensitivity analysis via simulation
#'
#' Evaluate the impact of continuous and categorical covariates on model
#' parameters and exposure metrics.  For each covariate the function perturbs
#' its value to selected quantiles (continuous) or observed categories
#' (categorical) while holding the remaining covariates at their reference
#' values, simulates from the pharmacometric model under parameter uncertainty,
#' and summarises the resulting percent change relative to the reference.
#'
#' Two data-source modes are supported (mutually exclusive):
#' \itemize{
#'   \item \strong{File mode} — supply \code{fpath_i} (path to a SimuRg JSON /
#'     RData object that contains SUMTAB, CATAB and COTAB).
#'   \item \strong{Table mode} — supply both \code{ds_parest} (parameter
#'     estimates) and \code{ds_covs} (covariate dataset).
#' }
#'
#' @inheritParams sg_dummy
#' @param ds_covs data.frame. Subject-level covariate dataset (one row per
#'   subject) containing both continuous and categorical covariate columns.
#'   Required when \code{fpath_i} is \code{NULL}; must be provided together
#'   with \code{ds_parest}.  Default is \code{NULL}.
#' @param cat_cov_l A named list. Each element defines one categorical covariate
#'   and must itself be a list where \code{par_vec} is required; other
#'   components are optional/permissible:
#'   \describe{
#'     \item{\code{par_vec}}{Required. Character vector. Model parameter(s)
#'       affected by this covariate (e.g. \code{c("ka")}).}
#'     \item{\code{NICENAME}}{Optional. Character or \code{NULL}. Display label.}
#'     \item{\code{REF}}{Optional. Character or \code{NULL}. Reference category value.
#'       If \code{NULL}, the first factor level (alphabetically) is used.}
#'   }
#' @param cont_cov_l A named list. Each element defines one continuous covariate
#'   and must itself be a list where \code{par_vec} is required; other
#'   components are optional/permissible:
#'   \describe{
#'     \item{\code{par_vec}}{Required. Character vector. Model parameter(s)
#'       affected by this covariate (e.g. \code{c("CL")}).}
#'     \item{\code{UTNAME}}{Optional. Character or \code{NULL}. Column name of the
#'       untransformed (back-transformed) covariate
#'       (e.g. \code{"AGE"}).  If \code{NULL} or \code{NA}, defaults to
#'       \code{NAME}.}
#'     \item{\code{REF}}{Optional. Character or numeric. Reference value for the
#'       covariate.  Use \code{"median"} to derive from data, or a numeric
#'       value.}
#'     \item{\code{NICENAME}}{Optional. Character or \code{NULL}. Display label for
#'       plots and tables (e.g. \code{"Age, years"}).}
#'   }
#' @param stimes numeric vector. Sampling time points for steady-state
#'   simulation (e.g. generated by a cycling function over dosing intervals).
#' @param et data.frame. Event (dosing) table.  Must contain at least
#'   columns \code{id}, \code{time}, \code{amt}, \code{evid} and
#'   \code{Regimen}.  Additional columns such as \code{ii}, \code{addl},
#'   \code{dur} and \code{Dose} are used when present.
#' @param npop integer. Number of population-level simulation replicates drawn
#'   from the parameter uncertainty distribution.  Default is \code{5}.
#' @param outputs character vector. Name(s) of the model output variable(s)
#'   to evaluate (e.g. \code{"Cc"} or \code{c("Cc", "Effect")}).
#'   Default is \code{"Cc"}.
#' @param aggr character vector. Exposure aggregation metric(s) applied over
#'   the simulation time grid.
#'   Allowed values: \code{"min"} (Cmin), \code{"max"} (Cmax),
#'   \code{"mean"} (Cmean), \code{"auc"} (AUC).
#' @param ci integer. Confidence interval (CI) level as a percentage.
#'   Allowed values: \code{95}, \code{90}, \code{80}, \code{70}, \code{50}.
#'   Default is \code{95}.
#' @param seed integer. Seed for the random number generator.  Default is \code{NULL}.
#'
#'@returns A named list of four elements:
#' \describe{
#'   \item{\code{$PARSENS}}{data.frame.
#'     Full sensitivity results for model parameters: percent change
#'     statistics (mean, median, percentiles) for each covariate–parameter
#'     combination, with columns \code{NICEN}, \code{VAR}, \code{KEY},
#'     \code{LAB}, \code{Type} and summary statistics
#'     (\code{mean} through \code{P975}).}
#'   \item{\code{$SUMPARSENS}}{data.frame.
#'     Compact summary table with columns \code{Parameter},
#'     \code{Covariate}, \code{Cov. percentile}, \code{Cov. value},
#'     \code{Mean} and \code{90\%CI}.}
#'   \item{\code{$EXPSENS}}{data.frame.
#'     Sensitivity of exposure metrics (Cmin, Cmax, Cmean, AUC) to covariates,
#'     structured identically to \code{$PARSENS}.}
#'   \item{\code{$COVREF}}{data.frame.
#'     Reference values for each covariate used in the analysis, with columns
#'     \code{COV}, \code{NICEN}, \code{REF_VALUE} and \code{REF_SOURCE}
#'     (\code{"reference"} for user-specified values, \code{"median"} when
#'     \code{REF = "median"} was used for a continuous covariate).}
#' }
#'
#' @details
#' The analysis proceeds in two stages for each covariate type (continuous and
#' categorical):
#' \enumerate{
#'   \item \strong{Parameter sensitivity} — the covariate is set to its
#'     quantile or category value and the model is evaluated at time zero (no
#'     ODE integration) to quantify the direct effect on each parameter.
#'   \item \strong{Exposure sensitivity} — the full time-course is simulated
#'     over \code{stimes} and exposure metrics (\code{aggr}) are computed.
#' }
#'
#' Uncertainty is propagated by sampling \code{npop} parameter vectors from a
#' multivariate normal distribution parameterised by the population estimates
#' and the estimation covariance matrix (\code{est_covmat}).
#'
#' @examples
#' \donttest{
#' library(dplyr)
#' library(rxode2)
#'
#' # --- Covariate definitions (only par_vec is required) ---
#' cont_cov_l <- list(
#'   LG_AGE = list(par_vec = c("CL")),
#'   LG_WEIGHT = list(
#'     UTNAME = "WEIGHT",
#'     NICENAME = "Weight, kg",
#'     par_vec = c("Vd")
#'   )
#' )
#'
#' cat_cov_l <- list(
#'   SEX = list(par_vec = c("ka")),
#'   CYP2C9 = list(
#'     REF = "1",
#'     NICENAME = "CYP2C9 genotype",
#'     par_vec = c("CL")
#'   )
#' )
#'
#' # --- Dosing ---
#' ev_t_input <- tribble(
#'   ~id, ~time, ~ii, ~amt, ~addl, ~dur, ~evid, ~Regimen,        ~Dose,
#'   1,   0,     336, 10,   21,    0.5,  1,     "0.3 mg/kg Q2W", 0.3
#' )
#' # --- Model ---
#' model <- RxODE({
#'   # Doses in mg
#'   # Time in hours
#'
#'   ### Parameter values
#'   # Typical values
#'   ka_pop = 0.073;
#'   Vd_pop = 14.8;
#'   CL_pop = 0.347;
#'
#'   # Random effects
#'   omega_ka = 0;
#'   omega_Vd = 0;
#'   omega_CL = 0;
#'
#'   # Covariate effect
#'   # Continuous
#'   beta_CL_LG_AGE = 0.49990114;
#'   beta_Vd_LG_WEIGHT = 0.60529433;
#'
#'   # Categorical
#'   beta_CL_CYP2C9_1_2 = -0.339;
#'   beta_CL_CYP2C9_1_3 = -0.574;
#'   beta_CL_CYP2C9_2_2 = -1.079;
#'   beta_CL_CYP2C9_2_3 = -0.745;
#'   beta_CL_CYP2C9_3_3 = -2.13;
#'
#'   beta_ka_SEX_1 = -0.12198035;
#'
#'   # Residual error
#'   Cc_b = 0;
#'
#'   # Transformations
#'   ka_tv = exp(ka_pop);
#'   Vd_tv = exp(Vd_pop);
#'   CL_tv = exp(CL_pop);
#'
#'   CL_multiplier = 1.0;  # Default/reference
#'   ka_multiplier = 1.0;
#'
#'   if (SEX == "1") {ka_multiplier = exp(beta_ka_SEX_1)}
#'
#'   if (CYP2C9 == "1") {
#'     CL_multiplier = exp(beta_CL_CYP2C9_1_2);
#'   } else if (CYP2C9 == "2") {
#'     CL_multiplier = exp(beta_CL_CYP2C9_1_3);
#'   } else if (CYP2C9 == "3") {
#'     CL_multiplier = exp(beta_CL_CYP2C9_2_2);
#'   } else if (CYP2C9 == "4") {
#'     CL_multiplier = exp(beta_CL_CYP2C9_2_3);
#'   } else if (CYP2C9 == "5") {
#'     CL_multiplier = exp(beta_CL_CYP2C9_3_3);
#'   }
#'
#'   ka = ka_tv*ka_multiplier*exp(omega_ka);
#'   Vd = Vd_tv*exp(beta_Vd_LG_WEIGHT * LG_WEIGHT + omega_Vd); #Vd_tv*exp(omega_Vd);
#'   CL = CL_tv*CL_multiplier*exp(beta_CL_LG_AGE * LG_AGE + omega_CL);
#'
#'
#'   ### Explicit functions
#'   Cc = Ac/Vd;
#'
#'   ### Initial conditions
#'   Ad(0) = 0;
#'   Ac(0) = 0;
#'
#'   ### Differential equations
#'   d/dt(Ad) = - ka*Ad;
#'   d/dt(Ac) = ka*Ad - CL*Cc;
#'
#'   Cc_ResErr = Cc*(1 + Cc_b);
#' })
#'
#' # --- Estimation covariance (mock) ---
#' pnames     <- parest$parameter
#' npar       <- length(pnames)
#' set.seed(1)
#' m_cov      <- matrix(0.02, npar, npar)
#' diag(m_cov) <- 0.05 + runif(npar, 0, 0.05)
#' m_cov      <- (m_cov + t(m_cov)) / 2
#' est_covmat <- as_tibble(cbind(X1 = pnames, as.data.frame(m_cov)))
#' names(est_covmat)[-1] <- pnames
#'
#' # --- Simulation times (steady-state cycle 10) ---
#' ss_cycle <- 10
#' stimes_ss <- c(
#'   ss_cycle * 4 * 7 * 24 + c(seq(0, 23.5, 0.5), seq(24, 335, 1)),
#'   ss_cycle * 4 * 7 * 24 + 2 * 7 * 24 + c(seq(0, 23.5, 0.5), seq(24, 335, 1))
#' )
#'
#' # --- Run ---
#' result <- sg_covsens_sim(
#'   fpath_i=NULL, ds_parest = parest, ds_covs = ds_covval,
#'   model = model, stimes = stimes_ss, et = ev_t_input,
#'   est_covmat = est_covmat, npop = 10,
#'   cont_cov_l = cont_cov_l, cat_cov_l = cat_cov_l,
#'   quantiles = c(0.2, 0.8), aggr = c("max"),
#'   outputs = "Cc"
#' )
#' print(result[["PARSENS"]])
#' print(result[["SUMPARSENS"]])
#' print(result[["EXPSENS"]])
#'
#' }
#'
#' @seealso \code{\link{sg_sim}}, \code{\link{read_smrg_obj}}
#'
#' @import dplyr
#' @importFrom purrr map map_dfr map_chr
#' @importFrom tidyr gather separate
#' @importFrom stringr str_detect str_c
#' @importFrom forcats fct_inorder
#' @export
sg_covsens_sim <- function(fpath_i = NULL, ds_parest = NULL, ds_covs = NULL,
                           model, stimes, et,
                           est_covmat = NULL,
                           npop = 5,
                           cont_cov_l, cat_cov_l,  quantiles = c(0.1, 0.9),
                           outputs  = "Cc", aggr = c("min", "max", "mean", "auc"),
                           seed = NULL, ci = 95){

  # --- Input validation ---
  # Data source: must provide either fpath_i alone, or both ds_parest,  ds_covs and est_covmat
  has_fpath   <- !is.null(fpath_i)
  has_ds      <- !is.null(ds_parest) && !is.null(ds_covs) && !is.null(est_covmat)
  has_ds_part <- !is.null(ds_parest) || !is.null(ds_covs) || !is.null(est_covmat)

  if (!has_fpath && !has_ds) {
    stop(
      "No data source provided. Supply either:\n",
      "  - 'fpath_i': path to a Simurg output object, OR\n",
      "  - all of the following: 'ds_parest' (parameter estimates), 'ds_covs' (covariate dataset) and 'est_covmat' (parameter estimation covariance matrix)."
    )
  }
  if (has_fpath && has_ds_part) {
    stop(
      "'fpath_i' and 'ds_parest'/'ds_covs'/'est_covmat' are mutually exclusive. ",
      "Supply one data source only."
    )
  }
  if (!has_fpath && has_ds_part && !has_ds) {
    missing_ds <- c()
    if (is.null(ds_parest)) missing_ds <- c(missing_ds, "'ds_parest'")
    if (is.null(ds_covs)) missing_ds <- c(missing_ds, "'ds_covs'")
    if (is.null(est_covmat)) missing_ds <- c(missing_ds, "'est_covmat'")
    stop(
      "Incomplete data source: missing ", paste(missing_ds, collapse = ", "), ". ",
      "When 'fpath_i' is NULL, provide all three: 'ds_parest', 'ds_covs' and 'est_covmat'."
    )
  }

  if (missing(model))    stop("'model' is required: provide the compiled rxode2/nlmixr model object.")
  if (missing(stimes))  stop("'stimes' is required: provide the steady-state simulation time points.")
  if (missing(et)) stop("'et' is required: provide the event table (dosing schedule).")
  if (missing(cont_cov_l)) stop("'cont_cov_l' is required: provide the continuous covariate definition list.")
  if (missing(cat_cov_l))  stop("'cat_cov_l' is required: provide the categorical covariate definition list.")

  cont_cov_l <- normalize_cont_cov_l(cont_cov_l)
  cat_cov_l <- normalize_cat_cov_l(cat_cov_l)

  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(as.integer(seed))
  }
  # Warn about non-default aggregation choices to alert on typos
  valid_aggr <- c("min", "max", "mean", "auc")
  bad_aggr   <- setdiff(aggr, valid_aggr)
  if (length(bad_aggr) > 0) {
    warning("Unrecognised aggregation function(s) in 'aggr': ",
            paste(bad_aggr, collapse = ", "),
            ". Valid options are: ", paste(valid_aggr, collapse = ", "), ".")
    aggr <- valid_aggr
  }
  # Ci-quantiles tibble
  ci <- as.numeric(ci)
  ci_allowed <- c(95, 90, 80, 70, 50)
  if (!ci %in% ci_allowed) {
    warning("'ci' must be one of: ", paste(ci_allowed, collapse = ", "), ". Using 95.")
    ci <- 95
  }
  # Tibble for CI - quantile correspondence (symmetric central intervals)
  #quantiles <- c("P025", "P05", "P10", "P20", "P30", "P50", "P70", "P80", "P90", "P95")
  ci_table <- tibble(
    CI   = c(95, 90, 80, 70, 50),
    LOW  = c("P025", "P05", "P10", "P15", "P25"),
    HIGH = c("P975", "P95", "P90", "P85", "P75")
  )
  ci_bounds <- ci_table %>% filter(CI == ci)
  ci_low  <- ci_bounds$LOW
  ci_high <- ci_bounds$HIGH
  ci_col  <- paste0(ci, "%CI")

  # -------------------------
  ds_catcov <- NULL
  if((!is.null(fpath_i))&(is.null(ds_parest))){
  obj_data <- read_smrg_obj(fpath_i)
  #ds_mod <-  obj_data$SDTAB
  par_sum <- obj_data$SUMTAB
  ds_catcov <- obj_data$CATAB
  ds_ccov <- obj_data$COTAB
  ds_covmat <- obj_data$COVMAT
  par_names <- par_sum %>% filter(EST=="ESTIMATED") %>% select(PAR) %>% pull()
  names(ds_covmat)[-1] <- par_names #Rename columns to parameter names
  ds_covmat <- as_tibble(ds_covmat) %>%
    mutate(X1 = par_names) %>%
    select(X1, everything())
  est_covmat <- ds_covmat
  
  data_fin <- ds_ccov %>% left_join(ds_catcov, by = "ID")
  ### Population parameters
  par_fin <- par_sum %>% rename(parameter = PAR, value = VALUE)
  } else if((is.null(fpath_i))&(!is.null(ds_parest))&(!is.null(ds_covs))){
    par_fin <- ds_parest
    data_fin <- ds_covs
  }

  par_pop <- par_fin %>% filter(str_detect(parameter, "_pop")) %>% select(parameter, value) %>% deframe()
  par_fin_tv <- par_fin %>% filter(str_detect(parameter, "_pop$") | str_detect(parameter, "^beta_")) %>% select(parameter, value) %>%
  #mutate(value = ifelse(str_detect(parameter, "_pop$") & !str_detect(parameter, "^beta_"), log(value), value)) %>% deframe()
  mutate(value = {
    cond <- str_detect(parameter, "_pop$") & !str_detect(parameter, "^beta_")
    ifelse(cond, log(ifelse(cond, value, 1)), value)
  }) %>% deframe()

  ### Reconstruct omega matrix (random effects on K_a and V_pop)
  d_omega <- par_fin %>% filter(str_detect(parameter, "omega_"))
  m_omega <- diag(d_omega$value, ncol = length(d_omega$value))
  colnames(m_omega) <- d_omega$parameter; rownames(m_omega) <- d_omega$parameter
  m_omega_full <- m_omega

  corr_set <- par_fin %>% filter(str_detect(parameter, "corr_"))
  if(nrow(corr_set) > 0){
    corr_set <- corr_set %>%
      separate(parameter, c("TYPE", "PAR1", "PAR2"), "_") %>%
      select(PAR1, PAR2, value) %>%
      mutate_at(vars(PAR1, PAR2), function(x){str_c("omega_", x)})

    omega_corr <- m_omega
    diag(omega_corr) <- rep(1, ncol(m_omega))
    for(i in 1:nrow(corr_set)){
      cs_i <- corr_set[i, ]
      omega_corr[cs_i$PAR1, cs_i$PAR2] <- cs_i$value
      omega_corr[cs_i$PAR2, cs_i$PAR1] <- cs_i$value
    }

    m_omega_full <- m_omega %*% omega_corr %*% m_omega
  }

  # ### Reconstruct residual error model matrix
  # d_reserr <- par_fin %>% filter(!str_detect(parameter, "_pop|omega_|corr_|beta_"))
  # m_reserr <- diag(d_reserr$value, ncol = length(d_reserr$value))
  # colnames(m_reserr) <- d_reserr$parameter; rownames(m_reserr) <- d_reserr$parameter

  m_theta_norm <- est_covmat %>% select_if(is.numeric) %>% as.matrix()
  colnames(m_theta_norm) <- est_covmat$X1; rownames(m_theta_norm) <- est_covmat$X1
  m_theta_norm_pop <- m_theta_norm[str_detect(rownames(m_theta_norm), "_pop|beta_"), str_detect(colnames(m_theta_norm), "_pop|beta_")]

  #####----Covariate tables----####
  # Validate and filter cont_cov_l against data_fin columns
  .missing_cont_name <- setdiff(names(cont_cov_l), names(data_fin))
  if (length(.missing_cont_name) > 0) {
    warning(
      "The following continuous covariate names from 'cont_cov_l' were not found ",
      "as columns in the covariate data and will be skipped: ",
      paste(.missing_cont_name, collapse = ", "), "."
    )
    cont_cov_l <- cont_cov_l[setdiff(names(cont_cov_l), .missing_cont_name)]
    if (length(cont_cov_l) == 0) {
      stop(
        "All continuous covariates in 'cont_cov_l' were removed because their named columns ",
        "are absent from the covariate data. Check 'cont_cov_l' and the dataset."
      )
    }
  }

  .missing_cont_utname <- unlist(Filter(Negate(is.null), map(names(cont_cov_l), function(cov_name) {
    ut <- cont_cov_l[[cov_name]]$UTNAME
    if (!is.null(ut) && !isTRUE(is.na(ut)) && ut != cov_name) ut else NULL
  }))) %>% setdiff(names(data_fin))
  if (length(.missing_cont_utname) > 0) {
    warning(
      "The following continuous covariate UTNAME(s) from 'cont_cov_l' were not found ",
      "as columns in the covariate data; affected covariates will be skipped: ",
      paste(.missing_cont_utname, collapse = ", "), "."
    )
    cont_cov_l <- cont_cov_l[purrr::map_lgl(names(cont_cov_l), function(cov_name) {
      ut <- cont_cov_l[[cov_name]]$UTNAME
      if (!is.null(ut) && !isTRUE(is.na(ut)) && ut != cov_name) !(ut %in% .missing_cont_utname) else TRUE
    })]
    if (length(cont_cov_l) == 0) {
      stop(
        "All continuous covariates in 'cont_cov_l' were removed because their UTNAME columns ",
        "are absent from the covariate data. Check 'cont_cov_l' and the dataset."
      )
    }
  }

  # Validate and filter cat_cov_l against data_fin columns
  .missing_cat_name <- setdiff(names(cat_cov_l), names(data_fin))
  if (length(.missing_cat_name) > 0) {
    warning(
      "The following categorical covariate names from 'cat_cov_l' were not found ",
      "as columns in the covariate data and will be skipped: ",
      paste(.missing_cat_name, collapse = ", "), "."
    )
    cat_cov_l <- cat_cov_l[setdiff(names(cat_cov_l), .missing_cat_name)]
    if (length(cat_cov_l) == 0) {
      stop(
        "All categorical covariates in 'cat_cov_l' were removed because their named columns ",
        "are absent from the covariate data. Check 'cat_cov_l' and the dataset."
      )
    }
  }

  # Continuous covariates table
  cont_cov_vec <- names(cont_cov_l)

  cont_cov <- map_dfr(cont_cov_vec, function(cov_name) {
    cov <- cont_cov_l[[cov_name]]
    btr <- if (is.null(cov$UTNAME) || isTRUE(is.na(cov$UTNAME))) cov_name else cov$UTNAME
    tibble(TR = cov_name, BTR = btr, PAR = list(cov$par_vec))
  })

  # Categorical covariates table
  cat_cov_vec <- names(cat_cov_l)

  if (is.null(ds_catcov)){
    ds_catcov <- ds_covs %>% select(all_of(c("ID",cat_cov_vec)))
  }

  cat_unique <- map(cat_cov_vec, function(x){
    ds_catcov[[x]] %>% unique()
  })
  names(cat_unique) <- cat_cov_vec

  # Build cat_cov automatically from cat_cov_l, cat_cov_vec and cat_unique
  cat_cov <- map_dfr(cat_cov_vec, function(x) {
    vals <- as.character(cat_unique[[x]])
    #ref  <- cat_cov_l[[x]]$REF
    ref <- if (is.null(cat_cov_l[[x]]$REF)) {
      as.character(levels(factor(data_fin[[x]]))[1])
    } else {
      cat_cov_l[[x]]$REF
    }
    tibble(
      COV    = x,
      COVVAL = vals,
      CATDES = vals,
      KEY    = as.integer(vals == ref),
      PAR    = list(cat_cov_l[[x]]$par_vec)
    )
  })

  nice_names_cont <- map_dfr(cont_cov_vec, function(x) {
    tibble(COV = x, NICEN = cov_nicename(cont_cov_l[[x]], x))
  })
  nice_names_cat <- map_dfr(cat_cov_vec, function(x) {
    tibble(COV = x, NICEN = cov_nicename(cat_cov_l[[x]], x))
  })
  nice_names <- rbind(nice_names_cont,nice_names_cat)

  ## Derive Covariate distributions
  ### Continuous covariates
  ID <- "ID"
  CCov_vec <- cont_cov$TR
  CatCov_vec <- cat_cov$COV %>% unique()

  data_fin <- data_fin %>%
    mutate(across(all_of(CatCov_vec), as.character))

  ds_cc_tr <- data_fin %>% select(all_of(c(ID)), all_of(cont_cov$TR)) %>% unique() %>%
    gather("TR", "TVALUE", -all_of(c(ID))) %>% left_join(cont_cov, by = "TR")

  ds_cc_btr <- data_fin %>% select(all_of(c(ID)), all_of(cont_cov$BTR)) %>% unique() %>%
    gather("BTR", "NVALUE", -all_of(c(ID)))
  ds_cc_btr_av <- ds_cc_btr %>% group_by(BTR) %>% summarise_at(vars(NVALUE), funSum_av)

  ds_cc <- ds_cc_tr %>%
    left_join(ds_cc_btr, by = c(ID, "BTR")) %>%
    left_join(ds_cc_btr_av, by = "BTR") %>%
    group_by(TR) %>%
    mutate(
      LP = quantile(TVALUE, quantiles[[1]]),
      UP = quantile(TVALUE, quantiles[[2]]))



  # Build per-covariate REF lookup from cont_cov_l
  ref_lookup <- map_dfr(cont_cov_vec, function(cov_name) {
    tibble(TR = cov_name, REF_spec = cont_cov_l[[cov_name]]$REF)
  })

  # Add REF column: use median(TVALUE) when REF_spec == "median", else use the numeric value
  ds_cc <- ds_cc %>%
    left_join(ref_lookup, by = "TR") %>%
    group_by(TR) %>%
    #mutate(REF = if_else(REF_spec == "median", median(TVALUE), as.numeric(REF_spec))) %>%
    mutate(REF = if_else(REF_spec == "median", median(TVALUE), suppressWarnings(as.numeric(REF_spec)))) %>%
    select(-REF_spec) %>%
    ungroup()


  #ds_cc_reflab <- select(ds_cc, COV = TR, median) %>% unique() %>% left_join(nice_names, by = "COV") %>% summarise(OUT = str_c(str_c(NICEN, " = ", round(median, 1)), collapse = "\n")) %>% pull(OUT)
  ds_cc_reflab <- select(ds_cc, COV = TR, median, LP, UP) %>% unique() %>% left_join(nice_names, by = "COV") %>% summarise(OUT = str_c(str_c(NICEN, " = ", round(median, 1), " (median)\n", "[",quantiles[[1]]*100,"th percentile: ", round(LP, 1), "; ",quantiles[[2]]*100,"th percentile: ", round(UP, 1), "]"), collapse = "\n")) %>% pull(OUT)

  # ds_cc_reflab <- select(ds_cc, COV = TR, median) %>% unique() %>% left_join(nice_names, by = "COV") %>% summarise(OUT = str_c(str_c(NICEN, " = ", round(median, 1)), collapse = "\n")) %>% pull(OUT) %>% str_c(., "\nResults shown for a dose calculated for the median (reference) patient")

  # Calculate low and upper quantile values for Back transformed covariates
  #Take into account that there are several continuous covariates!!!
  map_ccont <- function(target_ccont, data) {
    idx <- which.min(abs(data$TVALUE - target_ccont))
    data$NVALUE[idx]
  }
  ccont_lab_list <- list()
  for (i in seq_along(cont_cov_l)){
    cov_name_i <- names(cont_cov_l)[i]
    ds_cc_i <- ds_cc %>% filter(TR == cov_name_i)
    q_ccont <- c(unique(ds_cc_i$LP), unique(ds_cc_i$UP), unique(ds_cc_i$REF))

    ccont_labels <- sapply(
      q_ccont,
      map_ccont,
      data = ds_cc_i
    )

    names(ccont_labels) <- c("LP_BTR", "UP_BTR", "REF_BTR")
    ccont_lab_list[[cov_name_i]] <- ccont_labels
  }

  # Convert ccont_lab_list to a long data frame for joining: COV + KEY -> BCOVVAL
  ccont_lab_df <- map_dfr(names(ccont_lab_list), function(cov_name) {
    vals <- ccont_lab_list[[cov_name]]
    tibble(
      COV     = cov_name,
      KEY     = c("LP",       "UP",       "REF"),
      BCOVVAL = as.numeric(vals[c("LP_BTR", "UP_BTR", "REF_BTR")])
    )
  })

  cc_to_test <- ds_cc %>%
    select(COV = TR, BTR, PAR, mean, median, LP, UP, REF) %>%
    unique() %>%
    mutate(across(c(LP, UP, REF), as.numeric)) %>%
    gather("KEY", "COVVAL", -COV:-median) %>%
    left_join(ccont_lab_df, by = c("COV", "KEY"))

  ### Categorical covariates
  ds_catc <- data_fin %>% select(all_of(c(ID)), all_of(cat_cov$COV)) %>% unique()

  catc_to_test <- ds_catc %>% select(-all_of(c(ID))) %>% gather("COV", "COVVAL") %>% unique() %>% left_join(cat_cov, by = c("COV", "COVVAL"))

  # Reference values for continuous covariates: model scale (list key) and back-transformed (UTNAME)
  cont_cov_ref <- do.call(c, unname(map(names(cont_cov_l), function(cov_name) {
    cov <- cont_cov_l[[cov_name]]
    ref_row <- cc_to_test %>% filter(COV == cov_name, KEY == "REF")
    ut_name <- if (is.null(cov$UTNAME) || is.na(cov$UTNAME)) cov_name else cov$UTNAME
    setNames(
      list(ref_row %>% pull(COVVAL), ref_row %>% pull(BCOVVAL)),
      c(cov_name, ut_name)
    )
  })))
  cont_cov_ref <- cont_cov_ref[!duplicated(names(cont_cov_ref))]


  # Reference values for categorical covariates: use REF if defined, else first factor level
  cat_cov_ref <- setNames(
    map(cat_cov_vec, function(cov_name) {
      cov <- cat_cov_l[[cov_name]]
      if (is.null(cov$REF)) {
        as.character(levels(factor(data_fin[[cov_name]]))[1])
      } else {
        cov$REF
      }
    }),
    cat_cov_vec
  )

  all_cov_ref <- c(cont_cov_ref, cat_cov_ref)

  ### Create cov ref table for output
  cont_cov_ref_tbl <- map_dfr(cont_cov_vec, function(cov_name) {
    cov_def <- cont_cov_l[[cov_name]]
    nicen <- cov_nicename(cov_def, cov_name)

    ref_row <- cc_to_test %>%
      filter(COV == cov_name, KEY == "REF") %>%
      slice(1)

    ref_value <- ref_row$BCOVVAL
    if (length(ref_value) == 0 || is.na(ref_value)) {
      ref_value <- ref_row$COVVAL
    }

    if (length(ref_value) == 0 || is.na(ref_value)) {
      ref_value <- NA_character_
    } else {
      ref_value <- as.character(ref_value)
    }

    ref_source <- "reference"
    if (is.character(cov_def$REF) && tolower(cov_def$REF) == "median") {
      ref_source <- "median"
    }

    tibble(
      COV = cov_name,
      NICEN = nicen,
      REF_VALUE = ref_value,
      REF_SOURCE = ref_source
    )
  })

  cat_cov_ref_tbl <- map_dfr(cat_cov_vec, function(cov_name) {
    cov_def <- cat_cov_l[[cov_name]]
    nicen <- cov_nicename(cov_def, cov_name)

    tibble(
      COV = cov_name,
      NICEN = nicen,
      REF_VALUE = as.character(cat_cov_ref[[cov_name]]),
      REF_SOURCE = "reference"
    )
  })

  cov_ref_tbl <- bind_rows(cont_cov_ref_tbl, cat_cov_ref_tbl)

  # Combine with event table
  ev_t_base <- et %>%
    mutate(!!!all_cov_ref)

  ets_cc <- fun_EtCC(ev_t_base, cc_to_test)
  ets_catc <- fun_EtCC(ev_t_base, catc_to_test, T)

   # Sensitivity of model parameters to covariate values
  out_cov_exp_sens = NULL
  t_cov_sens_par = NULL
  out_cov_par_sens = NULL

  out_cov_par_sens <- bind_rows(
    fun_CovSens(et_sim_i = ets_cc, covs_i = nice_names$COV, nsim = npop,
                mod_fin_i = model,
                theta_i = par_fin_tv, omega_i = m_omega_full,
                thetamat_i = m_theta_norm_pop,
                nice_names_i = nice_names,quantiles = quantiles,
                var_exp = outputs, aggr = aggr
                ) %>% mutate(Type = "Continuous"),
    fun_CovSens(et_sim_i = ets_catc, cat = TRUE, covs_i = nice_names$COV, nsim = npop,
                mod_fin_i = model,
                theta_i = par_fin_tv, omega_i = m_omega_full,
                thetamat_i = m_theta_norm_pop,
                nice_names_i = nice_names,quantiles = quantiles,
                var_exp = outputs, aggr = aggr
                ) %>% mutate(Type = "Categorical")
    ) %>% mutate(LAB = fct_inorder(LAB), Type = fct_inorder(Type))

  ## Summary of sensitivity of model parameters to covariate values
  t_cov_sens_par <- out_cov_par_sens %>% mutate_at(vars(mean:P975), signif, 3) %>%
    mutate(!!ci_col := str_c(.data[[ci_low]], ", ", .data[[ci_high]]),
           KEY = ifelse(is.na(KEY), "Category", KEY),
           BCOVVAL = as.character(BCOVVAL),
           BCOVVAL = ifelse(is.na(BCOVVAL), CATDES, BCOVVAL)) %>%
    select(Parameter = VAR, Covariate = NICEN, `Cov. percentile` = KEY, `Cov. value` = BCOVVAL, Mean = mean, all_of(ci_col))
    # mutate(`90%CI` = str_c(P05, ", ", P95), KEY = ifelse(is.na(KEY), "Category", KEY), BCOVVAL = as.character(BCOVVAL), BCOVVAL = ifelse(is.na(BCOVVAL), CATDES, BCOVVAL)) %>%
    # select(Parameter = VAR, Covariate = NICEN, `Cov. percentile` = KEY, `Cov. value` = BCOVVAL, Mean = mean, `90%CI`)

  # Sensitivity of exposure parameters to covariate values
  out_cov_exp_sens <- bind_rows(
    fun_CovSens(ets_cc, expos = TRUE, covs_i = nice_names$COV, nsim = npop,
                stime_exp = stimes,
                mod_fin_i = model,
                theta_i = par_fin_tv, omega_i = m_omega_full,
                thetamat_i = m_theta_norm_pop,
                nice_names_i = nice_names, quantiles = quantiles,
                var_exp = outputs, aggr = aggr
                ) %>% mutate(Type = "Continuous"),
    fun_CovSens(ets_catc, cat = TRUE, expos = TRUE, covs_i = nice_names$COV, nsim = npop,
                stime_exp = stimes,
                mod_fin_i = model,
                theta_i = par_fin_tv, omega_i = m_omega_full,
                thetamat_i = m_theta_norm_pop,
                nice_names_i = nice_names, quantiles = quantiles,
                var_exp = outputs, aggr = aggr
                ) %>% mutate(Type = "Categorical")
  ) %>% mutate(LAB = fct_inorder(LAB), Type = fct_inorder(Type))

  covsens_res = list(PARSENS = out_cov_par_sens,
                     SUMPARSENS = t_cov_sens_par,
                     EXPSENS = out_cov_exp_sens,
                     COVREF = cov_ref_tbl)

  return(covsens_res)

}

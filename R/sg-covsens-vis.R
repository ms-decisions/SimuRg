#' Visualisation of covariate sensitivity analysis results
#'
#' Draw a forest-style graphic (points with uncertainty intervals) from the
#' sensitivity tables produced by \code{\link{sg_covsens_sim}}.  Each facet row
#' corresponds to a model output or exposure metric (\code{VAR}); each point is a
#' covariate scenario (\code{LAB}), coloured by covariate type (\code{Type}).
#'
#' Two views are available, matching the named elements of the simulation
#' output:
#' \itemize{
#'   \item \strong{\code{PARSENS}} — sensitivity of individual parameters
#'     (simulation at time zero, no ODE time course).
#'   \item \strong{\code{EXPSENS}} — sensitivity of exposure summaries
#'     (e.g. Cmin, Cmax, Cavg) after full simulation over \code{stimes} in
#'     \code{sg_covsens_sim}.
#' }
#' @inheritParams sg_dummy
#' @param plot_type Character scalar: \code{"PARSENS"} (default) or \code{"EXPSENS"}.
#' @param palette Colors for the \code{Type} scale (continuous vs
#'   categorical covariates, etc.).  Recycled if there are more levels than
#'   colors.  Default \code{MSDcol[c(1, 3, 4, 5, 6, 7)]}.
#' @param lab_y Axis label for the numeric scale (this becomes the horizontal
#'   axis after \code{coord_flip()}).  Default mentions a 95\% interval; change
#'   if you use different \code{ci_quantiles}.
#' @param cap Optional figure caption, passed to
#'   \code{ggplot2::labs(caption = ...)} (e.g. text describing reference
#'   covariate values).  Default \code{NULL}.
#'
#' @details
#' Values on the y-axis are ratios relative to the reference scenario: \code{1}
#' means no change.  In \code{sg_covsens_sim}, percent change relative to
#' reference is transformed to this scale before tabulation.  The shaded region
#' between \code{ci_limits} highlights a target interval; points whose
#' intervals lie largely inside can be read as scenarios consistent with that
#' criterion, subject to study-specific rules.
#'
#' The plot layers are drawn in order: reference band, error bars, points,
#' reference line at 1, dotted lines at \code{ci_limits}, then faceting by
#' \code{VAR} and flipped coordinates so labels read along the vertical axis.
#'
#' @returns A \code{ggplot2} object (inactive until printed or saved).  You can
#'   add further layers or themes with the usual \pkg{ggplot2} API.
#'
#' @examples
#' \dontrun{
#' # Typical workflow: run the simulation (see examples in ?sg_covsens_sim),
#' # then visualise parameter and exposure sensitivity.
#'
#' cont_cov_l <- list(
#' LG_AGE = list(NAME = "LG_AGE", UTNAME = "AGE",
#'               REF = "median", NICENAME = "Age, years",
#'               par_vec = c("CL")),
#' LG_WEIGHT = list(NAME = "LG_WEIGHT", UTNAME = "WEIGHT",
#'                  REF = "median", NICENAME = "Weight, kg",
#'                  par_vec = c("Vd"))
#' )
#'
#' cat_cov_l <- list(
#'   SEX = list(NAME = "SEX", NICENAME = "Sex",
#'              REF = "0", par_vec = c("ka")),
#'   CYP2C9 = list(NAME = "CYP2C9", NICENAME = "CYP2C9 genotype",
#'                 REF = NULL, par_vec = c("CL"))
#' )
#'
#' # --- Dosing ---
#' ev_t_input <- tribble(
#'  ~id, ~time, ~ii, ~amt, ~addl, ~dur, ~evid, ~Regimen,        ~Dose,
#'   1,   0,     336, 10,   21,    0.5,  1,     "0.3 mg/kg Q2W", 0.3
#' )
#' # --- Model ---
#' model <- RxODE({
#'  # Doses in mg
#'   # Time in hours
#'
#'  ### Parameter values
#'   # Typical values
#'  ka_pop = 0.073;
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
#'  # Transformations
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
#' result <- sg_covsens_sim(
#'   fpath_i = NULL, ds_parest = parest, ds_covs = ds_covval,
#'   model = model, stimes = stimes_ss, et = ev_t_input,
#'   est_covmat = est_covmat, npop = 10,
#'   cont_cov_l = cont_cov_l, cat_cov_l = cat_cov_l,
#'   quantiles = c(0.1, 0.9), aggr = c("min", "max", "mean"),
#'   outputs = "Cc"
#' )
#'
#' p_par <- sg_covsens_vis(result, plot_type = "PARSENS")
#' p_exp <- sg_covsens_vis(result, plot_type = "EXPSENS")
#' print(p_par)
#' print(p_exp)
#'
#' # Alternate interval columns (must exist in the sensitivity tables)
#' sg_covsens_vis(result, ci_quantiles = c("P05", "P95"))
#'
#' # Drop selected metrics from the exposure panel
#' sg_covsens_vis(result, plot_type = "EXPSENS", exclude_vars = c("Cc_Cmin"))
#' }
#'
#' @seealso \code{\link{sg_covsens_sim}}
#'
#' @export
sg_covsens_vis <- function(
    covsens_res,
    plot_type      = c("PARSENS", "EXPSENS"),
    exclude_vars   = NULL,
    ci_quantiles   = c("P025", "P975"),
    ci_limits      = c(0.8, 1.25),
    ci_band_alpha  = 0.2,
    ci_band_col    = "firebrick",
    ref_line_col   = "grey25",
    palette    = MSDcol[c(1, 3, 4, 5, 6, 7)],
    point_size     = 2.5,
    errorbar_width = 0.2,
    lab_y           = "Mean (95% CI)\nchange from reference",
    cap        = NULL
) {
  plot_type <- match.arg(plot_type)

  if (!plot_type %in% names(covsens_res)) {
    stop("'covsens_res' does not contain an element named '", plot_type, "'.")
  }
  if (length(ci_quantiles) != 2) {
    stop("'ci_quantiles' must be a character vector of length 2.")
  }
  if (length(ci_limits) != 2) {
    stop("'ci_limits' must be a numeric vector of length 2.")
  }

  ds <- covsens_res[[plot_type]]

  missing_q <- setdiff(ci_quantiles, colnames(ds))
  if (length(missing_q) > 0) {
    stop("Column(s) not found in data: ", paste(missing_q, collapse = ", "),
         ". Adjust 'ci_quantiles'.")
  }

  if (!is.null(exclude_vars)) {
    ds <- dplyr::filter(ds, !VAR %in% exclude_vars)
  }

  p <- ggplot2::ggplot(
    data = ds,
    ggplot2::aes(
      x    = LAB,
      y    = mean,
      ymin = .data[[ci_quantiles[1]]],
      ymax = .data[[ci_quantiles[2]]],
      col  = Type
    )
  ) +
    ggplot2::annotate(
      "rect",
      xmin  = -Inf, xmax = Inf,
      ymin  = ci_limits[1], ymax = ci_limits[2],
      fill  = ci_band_col,
      alpha = ci_band_alpha
    ) +
    ggplot2::geom_errorbar(width = errorbar_width) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::geom_hline(
      yintercept = 1,
      col = ref_line_col, lwd = 0.8, lty = "dashed"
    ) +
    ggplot2::geom_hline(
      yintercept = ci_limits,
      col = ci_band_col, lwd = 0.8, lty = "dotted"
    ) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(7)) +
    ggplot2::labs(x = NULL, y = lab_y, caption = cap) +
    ggplot2::facet_grid(VAR ~ ., scales = "free") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.minor  = ggplot2::element_blank(),
      legend.position   = "top",
      legend.background = ggplot2::element_rect(
        fill      = "white",
        linewidth = 0.15,
        linetype  = "solid",
        colour    = "black"
      )
    )

  p
}

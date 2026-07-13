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
#'   axis after \code{coord_flip()}).  Default mentions a 95% interval; change
#'   if you use different \code{ci_quantiles}.
#' @param cap Optional figure caption, passed to
#'   \code{ggplot2::labs(caption = ...)} (e.g. text describing reference
#'   covariate values).  Default \code{NULL}.
#' @param var_nice_names Named character vector. Display labels for exposure
#'   facet rows when \code{plot_type = "EXPSENS"}. Names must match the unique
#'   values in \code{VAR} (after \code{exclude_vars}); length equals the number
#'   of those unique values. Default \code{NULL} (facets use \code{VAR}).
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
#' \donttest{
#' library(tibble)
#' library(dplyr)
#' library(rxode2)
#' # Typical workflow: run the simulation (see examples in ?sg_covsens_sim),
#' # then visualise parameter and exposure sensitivity.
#'
#' cont_cov_l <- list(
#'   LG_AGE = list(NAME = "LG_AGE", UTNAME = "AGE",
#'                 REF = "median", NICENAME = "Age, years",
#'                 par_vec = c("CL")),
#'   LG_WEIGHT = list(NAME = "LG_WEIGHT", UTNAME = "WEIGHT",
#'                 REF = "median", NICENAME = "Weight, kg",
#'                 par_vec = c("Vd"))
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
#' p_par
#' p_exp
#'
#' # Alternate interval columns (must exist in the sensitivity tables)
#' p <- sg_covsens_vis(result, ci_quantiles = c("P05", "P95"))
#' p
#'
#' # Drop selected metrics from the exposure panel
#' p <- sg_covsens_vis(result, plot_type = "EXPSENS", exclude_vars = c("Cc_Cmin"))
#' p
#' }
#'
#' @seealso \code{\link{sg_covsens_sim}}
#'
#' @export
sg_covsens_vis <- function(
    covsens_res,
    plot_type      = c("PARSENS", "EXPSENS"),
    exclude_vars   = NULL,
    #ci_quantiles   = c("P025", "P975"),
    ci = 95,
    ci_limits      = c(0.8, 1.25),
    ci_band_alpha  = 0.2,
    ci_band_col    = "firebrick",
    ref_line_col   = "grey25",
    palette    = MSDcol[c(1, 3, 4, 5, 6, 7)],
    point_size     = 2.5,
    errorbar_width = 0.2,
    lab_y           = "standard",
    cap        = "standard",
    var_nice_names = NULL
) {
  plot_type <- match.arg(plot_type)

  if (!plot_type %in% names(covsens_res)) {
    stop("'covsens_res' does not contain an element named '", plot_type, "'.")
  }
  # if (length(ci_quantiles) != 2) {
  #   stop("'ci_quantiles' must be a character vector of length 2.")
  # }
  if (length(ci_limits) != 2) {
    stop("'ci_limits' must be a numeric vector of length 2.")
  }

  ds <- covsens_res[[plot_type]]



  if (!is.null(exclude_vars)) {
    ds <- dplyr::filter(ds, !VAR %in% exclude_vars)
  }



  if (!is.null(var_nice_names)) {
    use_var_nice_names <- TRUE
    if (!is.character(var_nice_names)) {
      warning("'var_nice_names' must be a named character vector. Ignoring 'var_nice_names'.")
      use_var_nice_names <- FALSE
    }
    if (use_var_nice_names && (is.null(names(var_nice_names)) || any(names(var_nice_names) == ""))) {
      warning("'var_nice_names' must be a named character vector; names must match unique 'VAR' values. Ignoring 'var_nice_names'.")
      use_var_nice_names <- FALSE
    }
    var_unique <- unique(as.character(ds$VAR))
    if (use_var_nice_names && length(var_nice_names) != length(var_unique)) {
      warning(
        "'var_nice_names' must have length ", length(var_unique),
        " (one label per unique 'VAR' value). Ignoring 'var_nice_names'."
      )
      use_var_nice_names <- FALSE
    }
    if (use_var_nice_names && !setequal(names(var_nice_names), var_unique)) {
      warning(
        "'var_nice_names' names must match the unique values in 'VAR': ",
        paste(var_unique, collapse = ", "),
        ". Ignoring 'var_nice_names'."
      )
      use_var_nice_names <- FALSE
    }
    if (use_var_nice_names) {

      # Allow duplicated display labels by de-duplicating factor levels.
      nice_levels <- unique(unname(var_nice_names[as.character(var_unique)]))
      ds <- ds %>%
        mutate(
          VAR_NICE = factor(
            unname(var_nice_names[as.character(VAR)]),
            levels = nice_levels
          )
        ) %>%
        relocate(VAR_NICE, .after = VAR)
    }
  }
 ###
 # Ci-quantiles tibble
  ci <- as.numeric(ci)
  ci_allowed <- c(95, 90, 80, 70, 50)
  if (!ci %in% ci_allowed) {
    warning("'ci' must be one of: ", paste(ci_allowed, collapse = ", "), ". Using 95.")
    ci <- 95
  }
  # Tibble for CI - quantile correspondence (symmetric central intervals)
  ci_table <- tibble(
    CI   = c(95, 90, 80, 70, 50),
    LOW  = c("P025", "P05", "P10", "P15", "P25"),
    HIGH = c("P975", "P95", "P90", "P85", "P75")
  )
  ci_bounds <- ci_table %>% filter(CI == ci)
  ci_low  <- ci_bounds$LOW
  ci_high <- ci_bounds$HIGH
  ci_col  <- paste0(ci, "%CI")
  ci_quantiles <- c(ci_low, ci_high)
  ###
  if (lab_y == "standard"){lab_y <- paste0("Mean (", ci, "% CI)\nchange from reference")}
  ### Caption writer
  if (!is.null(exclude_vars)) {
    ds <- dplyr::filter(ds, !VAR %in% exclude_vars)
  }
  .format_ref_value <- function(x) {
    if (length(x) == 0 || is.null(x) || is.na(x)) return(NA_character_)
    x_num <- suppressWarnings(as.numeric(x))
    if (!is.na(x_num)) {
      return(format(signif(x_num, 4), trim = TRUE, scientific = FALSE))
    }
    as.character(x)
  }

  .build_covref_caption <- function(x) {
    if (!"SUMPARSENS" %in% names(x)) return(NULL)
    sum_ds <- x$SUMPARSENS

    cov_col <- dplyr::case_when(
      "Covariate" %in% colnames(sum_ds) ~ "Covariate",
      "NICEN" %in% colnames(sum_ds) ~ "NICEN",
      TRUE ~ NA_character_
    )
    perc_col <- dplyr::case_when(
      "Cov. percentile" %in% colnames(sum_ds) ~ "Cov. percentile",
      "Cov.percentile" %in% colnames(sum_ds) ~ "Cov.percentile",
      TRUE ~ NA_character_
    )
    val_col <- dplyr::case_when(
      "Cov. value" %in% colnames(sum_ds) ~ "Cov. value",
      "Cov.value" %in% colnames(sum_ds) ~ "Cov.value",
      TRUE ~ NA_character_
    )

    if (any(is.na(c(cov_col, perc_col, val_col)))) return(NULL)

    sum_cov <- sum_ds %>%
      dplyr::select(COVARIATE = dplyr::all_of(cov_col),
                    PERCENTILE = dplyr::all_of(perc_col),
                    VALUE = dplyr::all_of(val_col)) %>%
      dplyr::filter(!is.na(COVARIATE), !is.na(PERCENTILE), !is.na(VALUE)) %>%
      dplyr::mutate(
        PERCENTILE = as.character(PERCENTILE),
        VALUE = as.character(VALUE)
      ) %>%
      dplyr::distinct()

    if (nrow(sum_cov) == 0) return(NULL)

    ref_cov <- NULL
    if ("COVREF" %in% names(x)) {
      ref_cov <- x$COVREF
      required_ref_cols <- c("NICEN", "REF_VALUE", "REF_SOURCE")
      if (!all(required_ref_cols %in% colnames(ref_cov))) {
        ref_cov <- NULL
      } else {
        ref_cov <- ref_cov %>%
          dplyr::select(COVARIATE = NICEN, REF_VALUE, REF_SOURCE) %>%
          dplyr::filter(!is.na(COVARIATE), !is.na(REF_VALUE)) %>%
          dplyr::distinct(COVARIATE, .keep_all = TRUE)
      }
    }

    cov_order <- unique(sum_cov$COVARIATE)

    ref_line_for <- function(cov_name) {
      if (is.null(ref_cov)) return(NULL)
      r <- ref_cov %>% dplyr::filter(COVARIATE == cov_name)
      if (nrow(r) == 0) return(NULL)
      ref_value <- .format_ref_value(r$REF_VALUE[[1]])
      ref_source <- r$REF_SOURCE[[1]]
      if (is.null(ref_source) || is.na(ref_source) || ref_source == "") {
        ref_source <- "reference"
      }
      paste0(cov_name, " = ", ref_value, " (", ref_source, ")")
    }

    .format_percentile_label <- function(x) {
      x_chr <- as.character(x)
      x_chr <- gsub("perc\\.", "percentile", x_chr, fixed = FALSE)
      x_chr <- gsub("\\s+", " ", x_chr)
      trimws(x_chr)
    }

    perc_line_for <- function(cov_name) {
      cdat <- sum_cov %>% dplyr::filter(COVARIATE == cov_name)
      if (nrow(cdat) == 0) return(NULL)

      cdat <- cdat %>%
        dplyr::mutate(
          PCT_NUM = suppressWarnings(as.numeric(sub("^([0-9.]+).*", "\\1", PERCENTILE)))
        ) %>%
        dplyr::filter(!is.na(PCT_NUM)) %>%
        dplyr::arrange(PCT_NUM) %>%
        dplyr::distinct(PCT_NUM, .keep_all = TRUE)

      if (nrow(cdat) == 0) return(NULL)

      if (nrow(cdat) == 1) {
        paste0("[",
               .format_percentile_label(cdat$PERCENTILE[[1]]), ": ",
               .format_ref_value(cdat$VALUE[[1]]), "]")
      } else {
        low <- cdat[1, , drop = FALSE]
        high <- cdat[nrow(cdat), , drop = FALSE]
        paste0("[",
               .format_percentile_label(low$PERCENTILE[[1]]), ": ", .format_ref_value(low$VALUE[[1]]), "; ",
               .format_percentile_label(high$PERCENTILE[[1]]), ": ", .format_ref_value(high$VALUE[[1]]), "]")
      }
    }

    cov_blocks <- lapply(cov_order, function(cov_name) {
      parts <- c(ref_line_for(cov_name), perc_line_for(cov_name))
      parts <- parts[!is.na(parts) & parts != ""]
      if (length(parts) == 0) return(NULL)
      paste(parts, collapse = "\n")
    })
    cov_blocks <- unlist(cov_blocks, use.names = FALSE)
    cov_blocks <- cov_blocks[!is.na(cov_blocks) & cov_blocks != ""]
    if (length(cov_blocks) == 0) return(NULL)
    paste(cov_blocks, collapse = "\n\n")
  }

  if (cap=="standard") {
    cap <- .build_covref_caption(covsens_res)
  }

  #facet_layer <- if (plot_type == "EXPSENS" && !is.null(var_nice_names)) {
  facet_layer <- if (plot_type == "EXPSENS" && "VAR_NICE" %in% names(ds)) {
    ggplot2::facet_grid(VAR_NICE ~ ., scales = "free")
  } else {
    ggplot2::facet_grid(VAR ~ ., scales = "free")
  }

  ###

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
    #ggplot2::facet_grid(VAR ~ ., scales = "free") +
    facet_layer +
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

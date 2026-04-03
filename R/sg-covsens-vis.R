#' Visualise covariate sensitivity analysis results
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
#'
#' @param covsens_res Named list as returned by \code{sg_covsens_sim()}.  Must
#'   contain the element selected by \code{type} (\code{PARSENS} and/or
#'   \code{EXPSENS} data.frames with columns \code{LAB}, \code{VAR},
#'   \code{mean}, \code{Type}, and the interval columns named by
#'   \code{ci_quantiles}).
#' @param type Character scalar: \code{"PARSENS"} (default) or \code{"EXPSENS"}.
#' @param exclude_vars Character vector of \code{VAR} levels to omit (e.g.
#'   \code{"Cc_Cmin"}).  \code{NULL} keeps all rows.
#' @param ci_quantiles Character vector of length 2: names of the lower and
#'   upper uncertainty columns in the sensitivity table, in that order.
#'   Defaults \code{c("P025", "P975")} to match the default percentiles in
#'   \code{sg_covsens_sim}.  Use other names (e.g. \code{c("P05", "P95")}) if
#'   you changed \code{quantiles} in the simulation and the columns exist.
#' @param ci_limits Numeric vector of length 2: lower and upper bounds of the
#'   shaded acceptance band and of the dotted horizontal guides.  Default
#'   \code{c(0.8, 1.25)} is a common bioequivalence-style window on the ratio
#'   scale.
#' @param ci_band_alpha Numeric in \eqn{[0,1]}: transparency of the shaded band.
#'   Default \code{0.2}.
#' @param ci_band_col Colour for the band fill and dotted limit lines.
#'   Default \code{"firebrick"}.
#' @param ref_line_col Colour for the dashed horizontal line at \code{y = 1}
#'   (no change from reference).  Default \code{"grey25"}.
#' @param col_palette Colours for the \code{Type} scale (continuous vs
#'   categorical covariates, etc.).  Recycled if there are more levels than
#'   colours.  Default \code{MSDcol[c(1, 3, 4, 5, 6, 7)]}.
#' @param point_size Point size for \code{geom_point}.  Default \code{2.5}.
#' @param errorbar_width Width argument for \code{geom_errorbar}.  Default
#'   \code{0.2}.
#' @param ylab Axis label for the numeric scale (this becomes the horizontal
#'   axis after \code{coord_flip()}).  Default mentions a 95\% interval; change
#'   if you use different \code{ci_quantiles}.
#' @param caption Optional figure caption, passed to
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
#' result <- sg_covsens_sim(
#'   fpath_i = NULL, ds_parest = parest, ds_covs = ds_covval,
#'   model = model, stimes = stimes_ss, et = ev_t_input,
#'   est_covmat = est_covmat, npop = 10,
#'   cont_cov_l = cont_cov_l, cat_cov_l = cat_cov_l,
#'   quantiles = c(0.1, 0.9), aggr = c("min", "max", "mean"),
#'   outputs = "Cc"
#' )
#'
#' p_par <- sg_covsens_vis(result, type = "PARSENS")
#' p_exp <- sg_covsens_vis(result, type = "EXPSENS")
#' print(p_par)
#' print(p_exp)
#'
#' # Alternate interval columns (must exist in the sensitivity tables)
#' sg_covsens_vis(result, ci_quantiles = c("P05", "P95"))
#'
#' # Drop selected metrics from the exposure panel
#' sg_covsens_vis(result, type = "EXPSENS", exclude_vars = c("Cc_Cmax"))
#' }
#'
#' @seealso \code{\link{sg_covsens_sim}}
#'
#' @export
sg_covsens_vis <- function(
    covsens_res,
    type           = c("PARSENS", "EXPSENS"),
    exclude_vars   = NULL,
    ci_quantiles   = c("P025", "P975"),
    ci_limits      = c(0.8, 1.25),
    ci_band_alpha  = 0.2,
    ci_band_col    = "firebrick",
    ref_line_col   = "grey25",
    col_palette    = MSDcol[c(1, 3, 4, 5, 6, 7)],
    point_size     = 2.5,
    errorbar_width = 0.2,
    ylab           = "Mean (95% CI)\nchange from reference",
    caption        = NULL
) {
  type <- match.arg(type)

  if (!type %in% names(covsens_res)) {
    stop("'covsens_res' does not contain an element named '", type, "'.")
  }
  if (length(ci_quantiles) != 2) {
    stop("'ci_quantiles' must be a character vector of length 2.")
  }
  if (length(ci_limits) != 2) {
    stop("'ci_limits' must be a numeric vector of length 2.")
  }

  ds <- covsens_res[[type]]

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
    ggplot2::scale_color_manual(values = col_palette) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(7)) +
    ggplot2::labs(x = NULL, y = ylab, caption = caption) +
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

#' Visualise covariate sensitivity analysis results
#'
#' Produces a forest-style dot-and-errorbar plot of the covariate sensitivity
#' results returned by \code{\link{sg_covsens_sim}}.  Two output types are
#' supported: \code{"PARSENS"} (sensitivity of model parameters) and
#' \code{"EXPSENS"} (sensitivity of simulated exposure metrics).
#'
#' @param covsens_res Named list returned by \code{sg_covsens_sim()}, which
#'   must contain elements \code{PARSENS} and \code{EXPSENS}.
#' @param type Character scalar; which sensitivity output to plot.
#'   One of \code{"PARSENS"} (default) or \code{"EXPSENS"}.
#' @param exclude_vars Character vector of \code{VAR} values to suppress from
#'   the plot (e.g. \code{"Cmin, ug/mL"}).  Default \code{NULL} (show all).
#' @param ci_quantiles Character vector of length 2 giving the names of the
#'   lower and upper confidence-limit columns in the data, in that order.
#'   Default \code{c("P025", "P975")} (i.e. 95 \% CI).
#' @param ci_limits Numeric vector of length 2 defining the lower and upper
#'   bounds of the reference (acceptance) band drawn as a shaded rectangle and
#'   dotted lines.  Default \code{c(0.8, 1.25)}.
#' @param ci_band_alpha Numeric scalar; alpha transparency of the filled
#'   reference band.  Default \code{0.2}.
#' @param ci_band_col Colour of the reference-band fill and dotted boundary
#'   lines.  Default \code{"firebrick"}.
#' @param ref_line_col Colour of the horizontal reference line at \code{y = 1}.
#'   Default \code{"grey25"}.
#' @param col_palette Character vector of colours used for the \code{Type}
#'   aesthetic (Continuous, Categorical, …).  The vector is recycled if there
#'   are more levels than colors.  Default uses the six-colour MSD palette
#'   subset \code{MSDcol[c(1, 3, 4, 5, 6, 7)]}.
#' @param point_size Numeric scalar; size of the point geom.  Default \code{2.5}.
#' @param errorbar_width Numeric scalar; width of the error-bar caps.
#'   Default \code{0.2}.
#' @param ylab Character string for the y-axis (horizontal after
#'   \code{coord_flip}) label.
#'   Default \code{"Mean (95\% CI)\nchange from reference"}.
#' @param caption Character string passed to \code{labs(caption = ...)}.
#'   Typically the reference-value label produced by \code{sg_covsens_sim()}.
#'   Default \code{NULL}.
#'   @param panel_height Numeric scalar; fixed height in cm applied to every
#'   facet panel.  Useful when many covariate rows and multiple \code{VAR}
#'   facets make panels too narrow to read.  Requires the \pkg{ggh4x} package.
#'   Default \code{NULL} (let ggplot2 determine panel height automatically).
#'
#' @return A \code{ggplot} object.
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
    caption        = NULL,
    panel_height   = NULL
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
    ggplot2::scale_y_continuous(
      name   = ylab,
      breaks = scales::pretty_breaks(7)
    ) +
    ggplot2::labs(x = NULL, caption = caption) +
    ggplot2::facet_grid(VAR ~ ., scales = "free") +
    ggplot2::coord_flip() +
    ggplot2::theme(
      legend.position   = "top",
      legend.background = ggplot2::element_rect(
        fill     = "white",
        linewidth = 0.15,
        linetype  = "solid",
        colour    = "black"
      )
    )

  if (!is.null(panel_height)) {
    if (requireNamespace("ggh4x", quietly = TRUE)) {
      p <- p + ggh4x::force_panelsizes(rows = grid::unit(panel_height, "cm"))
    } else {
      warning("Package 'ggh4x' is required for 'panel_height'. ",
              "Install it with: install.packages('ggh4x')")
    }
  }

  p
}

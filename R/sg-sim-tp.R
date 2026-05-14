## Author: Victor Sokolov
## First created: 2025-11-13
## Description: Time-profile plotting with optional summary bands
## Keywords: SimuRg, time profiles, diagnostics


#' Time-profile plotting with optional summary bands
#'
#' Creates a ggplot time-profile plot from longitudinal data with optional
#' summarization (mean/median/geometric mean), variance bands (SD/SE/IQR),
#' percentile ribbons, faceting, and log scaling.
#'
#' @inheritParams sg_dummy
#' @param bands_i vector of characters. Character vector of length 2 with column names for custom ymin/ymax ribbons. Default is `NULL`
#' @param cent_i string. Central tendency measure. One of `'mean'`, `'median'`, `'geom_mean'`, or `NULL` for raw values. Default is `NULL`
#' @param vrns_i string. Variance band type. One of `'SD'`, `'SE'`, `'IQR'`, or `NULL`. Default is `NULL`
#' @param lperc_i numeric. Lower percentile for percentile ribbons (must be provided together with `uperc_i`). Default is `NULL`
#' @param uperc_i numeric. Upper percentile for percentile ribbons (must be provided together with `lperc_i`). Default is `NULL`
#' @param add_points numeric. Size of points to add to the plot. If > 0, points will be added. Default is `0`
#' @param grid_i string. Faceting formula for `facet_grid` (e.g., `'~VAR'` or `'A~B'`). Default is `NULL`
#'
#' @return A ggplot object with lines and optional ribbons/points/facets.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom rlang sym syms expr
#'
#' @examples
#' make_extended_mock_data <- function() {
#'  data.frame(
#'    TIME = rep(1:4, times = 6),
#'    VALUE = rnorm(24, mean = 10, sd = 2),
#'    VAR = rep(rep(c("A", "B"), each = 4), times = 3),
#'    Regimen = rep(c("R1", "R2", "R3"), each = 8),
#'    CAT1 = rep(c("Cat1", "Cat2"), each = 12)
#'  )
#' }
#' ds_sim <- make_extended_mock_data()
#' p <- sg_sim_tp(ds_i = ds_sim, group_i = 'VAR', col_i = 'VAR', fill_i = 'VAR',
#'                wrap_i = '~VAR', wrap_ncol = 2)
#' @export
sg_sim_tp <- function(
  ds_i,
  time_col = 'TIME',
  val_col = 'VALUE',
  group_i = 'VAR',
  bands_i = NULL,
  cent_i = NULL,
  vrns_i = NULL,
  lperc_i = NULL,
  uperc_i = NULL,
  add_points = 0,
  col_i = NULL,
  fill_i = NULL,
  lty_i = NULL,
  shp_i = NULL,
  grid_i = NULL,
  wrap_i = NULL,
  free_stat = 'free',
  wrap_ncol = NULL,
  wrap_nrow = NULL,
  min_x = NA,
  max_x = NA,
  min_y = NA,
  max_y = NA,
  log_y = FALSE,
  log_x = FALSE
) {

  # Determine all grouping columns
  group_cols <- unique(c(group_i, col_i, fill_i, lty_i, shp_i))
  group_cols <- group_cols[!is.null(group_cols)]

  # Check if summarization is needed
  needs_summary <- !is.null(cent_i) || !is.null(vrns_i) || !is.null(lperc_i) || !is.null(uperc_i)

  if (needs_summary) {
    # Prepare grouping for summarization
    group_vars <- c(time_col, group_cols)

    # Summarize data
    ds_plot <- ds_i %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(
        n = n(),
        mean_val = if (!is.null(cent_i) && cent_i == "mean") mean(!!sym(val_col), na.rm = TRUE) else NA_real_,
        median_val = if (!is.null(cent_i) && cent_i == "median") median(!!sym(val_col), na.rm = TRUE) else NA_real_,
        geom_mean_val = if (!is.null(cent_i) && cent_i == "geom_mean") exp(mean(log(!!sym(val_col)), na.rm = TRUE)) else NA_real_,
        sd_val = if (!is.null(vrns_i) && vrns_i == "SD") sd(!!sym(val_col), na.rm = TRUE) else NA_real_,
        se_val = if (!is.null(vrns_i) && vrns_i == "SE") sd(!!sym(val_col), na.rm = TRUE) / sqrt(n()) else NA_real_,
        q25 = if (!is.null(vrns_i) && vrns_i == "IQR") quantile(!!sym(val_col), 0.25, na.rm = TRUE) else NA_real_,
        q75 = if (!is.null(vrns_i) && vrns_i == "IQR") quantile(!!sym(val_col), 0.75, na.rm = TRUE) else NA_real_,
        lperc = if (!is.null(lperc_i)) quantile(!!sym(val_col), lperc_i, na.rm = TRUE) else NA_real_,
        uperc = if (!is.null(uperc_i)) quantile(!!sym(val_col), uperc_i, na.rm = TRUE) else NA_real_,
        .groups = 'drop'
      )

    # Calculate central tendency value
    if (!is.null(cent_i)) {
      if (cent_i == "mean") {
        ds_plot <- ds_plot %>% mutate(y_central = mean_val)
      } else if (cent_i == "median") {
        ds_plot <- ds_plot %>% mutate(y_central = median_val)
      } else if (cent_i == "geom_mean") {
        ds_plot <- ds_plot %>% mutate(y_central = geom_mean_val)
      }
    } else {
      # If no central tendency, use original values (shouldn't happen in summary mode)
      ds_plot <- ds_plot %>% mutate(y_central = mean_val)
    }

    # Calculate variance bands if specified
    if (!is.null(vrns_i)) {
      if (vrns_i == "SD") {
        ds_plot <- ds_plot %>%
          mutate(
            y_lower = y_central - sd_val,
            y_upper = y_central + sd_val
          )
      } else if (vrns_i == "SE") {
        ds_plot <- ds_plot %>%
          mutate(
            y_lower = y_central - se_val,
            y_upper = y_central + se_val
          )
      } else if (vrns_i == "IQR") {
        ds_plot <- ds_plot %>%
          mutate(
            y_lower = q25,
            y_upper = q75
          )
      }
    }

    # Add percentile bands if specified
    if (!is.null(lperc_i) && !is.null(uperc_i)) {
      ds_plot <- ds_plot %>%
        mutate(
          y_lower_perc = lperc,
          y_upper_perc = uperc
        )
    }

  } else {
    # No summarization needed - use raw data
    ds_plot <- ds_i %>%
      rename(y_central = !!sym(val_col))
  }

  # Initialize ggplot
  p <- ggplot(ds_plot, aes(x = !!sym(time_col), y = y_central))

  # Add ribbon for bands_i if specified
  if (!is.null(bands_i) && length(bands_i) == 2) {
    ribbon_mapping <- aes(ymin = !!sym(bands_i[1]), ymax = !!sym(bands_i[2]))
    if (!is.null(fill_i)) {
      ribbon_mapping$fill <- sym(fill_i)
    }
    p <- p + geom_ribbon(ribbon_mapping, alpha = 0.3)
  }

  # Add ribbon for percentile bands
  if (!is.null(lperc_i) && !is.null(uperc_i) && needs_summary) {
    ribbon_mapping <- aes(ymin = y_lower_perc, ymax = y_upper_perc)
    if (!is.null(fill_i)) {
      ribbon_mapping$fill <- sym(fill_i)
    }
    p <- p + geom_ribbon(ribbon_mapping, alpha = 0.3)
  }

  # Add ribbon for variance bands
  if (!is.null(vrns_i) && needs_summary) {
    ribbon_mapping <- aes(ymin = y_lower, ymax = y_upper)
    if (!is.null(fill_i)) {
      ribbon_mapping$fill <- sym(fill_i)
    }
    p <- p + geom_ribbon(ribbon_mapping, alpha = 0.3)
  }

  # Build aesthetic mappings for line
  line_mapping <- aes()
  if (!is.null(col_i)) {
    line_mapping$colour <- sym(col_i)
  }
  if (!is.null(lty_i)) {
    line_mapping$linetype <- sym(lty_i)
  }
  if (length(group_cols) > 0) {
    line_mapping$group <- expr(interaction(!!!syms(group_cols)))
  }

  # Add line layer
  p <- p + geom_line(line_mapping, linewidth = 0.8)

  # Add points if requested
  if (add_points > 0) {
    point_mapping <- aes()
    if (!is.null(col_i)) {
      point_mapping$colour <- sym(col_i)
    }
    if (!is.null(shp_i)) {
      point_mapping$shape <- sym(shp_i)
    }
    p <- p + geom_point(point_mapping, size = add_points)
  }

  # Add faceting
  if (!is.null(grid_i)) {
    p <- p + facet_grid(as.formula(grid_i), scales = free_stat)
  } else if (!is.null(wrap_i)) {
    p <- p + facet_wrap(as.formula(wrap_i), scales = free_stat,
                        ncol = wrap_ncol, nrow = wrap_nrow)
  }

  # Apply log scales
  if (!is.na(log_y) && log_y) {
    p <- p + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10", scales::math_format(10^.x)))
  }

  if (!is.na(log_x) && log_x) {
    p <- p + scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10", scales::math_format(10^.x)))
  }

  # Apply coordinate limits
  xlim_vec <- c(
    if (!is.na(min_x)) min_x else NA,
    if (!is.na(max_x)) max_x else NA
  )

  ylim_vec <- c(
    if (!is.na(min_y)) min_y else NA,
    if (!is.na(max_y)) max_y else NA
  )

  if (!all(is.na(xlim_vec)) || !all(is.na(ylim_vec))) {
    p <- p + coord_cartesian(
      xlim = if (!all(is.na(xlim_vec))) xlim_vec else NULL,
      ylim = if (!all(is.na(ylim_vec))) ylim_vec else NULL
    )
  }

  return(p)
}



## Author: Victor Sokolov
## First created: 2025-11-13
## Description:
## Keywords:

#####--------------- Load functions and libraries ---------------#####
# library(tidyverse)
# theme_set(theme_bw())
# theme_update(panel.grid.minor = element_blank())
#
#
# #####--------------- Read the data ---------------#####
# ds_sim <- read_csv("sim_timeprof_test.csv", col_types = cols())
#####--------------- Function body ---------------#####
#' Time-profile plotting with optional summary bands
#'
#' Creates a ggplot time-profile plot from longitudinal data with optional
#' summarization (mean/median/geometric mean), variance bands (SD/SE/IQR),
#' percentile ribbons, faceting, and log scaling.
#'
#' @param ds_i Data frame with source data.
#' @param x_col Name of x-axis (time) column (default: 'TIME').
#' @param y_col Name of y-axis (value) column (default: 'VALUE').
#' @param group_i Primary grouping variable for lines (default: 'VAR').
#' @param bands_i Character vector of length 2 with column names for custom ymin/ymax ribbons.
#' @param cent_i Central tendency: one of 'mean', 'median', 'geom_mean', or NULL for raw values.
#' @param vrns_i Variance band type: one of 'SD', 'SE', 'IQR', or NULL.
#' @param lperc_i,uperc_i Lower/upper percentiles for percentile ribbons (both must be provided).
#' @param add_points Numeric > 0 to add points (size).
#' @param col_i,fill_i,lty_i,shp_i Optional aesthetics (colour, fill, linetype, shape).
#' @param grid_i,wrap_i Faceting formulas for grid or wrap;
#' @param free_stat Facet scaling ('free', 'free_x', etc.).
#' @param wrap_ncol,wrap_nrow Number of columns/rows for facet_wrap.
#' @param x_min,x_max,y_min,y_max Axis limits (use NA to skip).
#' @param log_y,log_x Logical; apply log10 scale to y or x.
#'
#' @return A ggplot object with lines and optional ribbons/points/facets.
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom rlang sym syms expr
#'
#' @examples
#' \dontrun{
#' p <- sg_sim_tp(ds_i = ds_sim, group_i = 'VAR', col_i = 'VAR')
#' }

#####--------------- Function body ---------------#####
sg_sim_tp <- function(
  ds_i,
  x_col = 'TIME',
  y_col = 'VALUE',
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
  x_min = NA,
  x_max = NA,
  y_min = NA,
  y_max = NA,
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
    group_vars <- c(x_col, group_cols)

    # Summarize data
    ds_plot <- ds_i %>%
      group_by(across(all_of(group_vars))) %>%
      summarise(
        n = n(),
        mean_val = if (!is.null(cent_i) && cent_i == "mean") mean(!!sym(y_col), na.rm = TRUE) else NA_real_,
        median_val = if (!is.null(cent_i) && cent_i == "median") median(!!sym(y_col), na.rm = TRUE) else NA_real_,
        geom_mean_val = if (!is.null(cent_i) && cent_i == "geom_mean") exp(mean(log(!!sym(y_col)), na.rm = TRUE)) else NA_real_,
        sd_val = if (!is.null(vrns_i) && vrns_i == "SD") sd(!!sym(y_col), na.rm = TRUE) else NA_real_,
        se_val = if (!is.null(vrns_i) && vrns_i == "SE") sd(!!sym(y_col), na.rm = TRUE) / sqrt(n()) else NA_real_,
        q25 = if (!is.null(vrns_i) && vrns_i == "IQR") quantile(!!sym(y_col), 0.25, na.rm = TRUE) else NA_real_,
        q75 = if (!is.null(vrns_i) && vrns_i == "IQR") quantile(!!sym(y_col), 0.75, na.rm = TRUE) else NA_real_,
        lperc = if (!is.null(lperc_i)) quantile(!!sym(y_col), lperc_i, na.rm = TRUE) else NA_real_,
        uperc = if (!is.null(uperc_i)) quantile(!!sym(y_col), uperc_i, na.rm = TRUE) else NA_real_,
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
      rename(y_central = !!sym(y_col))
  }

  # Initialize ggplot
  p <- ggplot(ds_plot, aes(x = !!sym(x_col), y = y_central))

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
    p <- p + scale_y_log10()
  }

  if (!is.na(log_x) && log_x) {
    p <- p + scale_x_log10()
  }

  # Apply coordinate limits
  xlim_vec <- c(
    if (!is.na(x_min)) x_min else NA,
    if (!is.na(x_max)) x_max else NA
  )

  ylim_vec <- c(
    if (!is.na(y_min)) y_min else NA,
    if (!is.na(y_max)) y_max else NA
  )

  if (!all(is.na(xlim_vec)) || !all(is.na(ylim_vec))) {
    p <- p + coord_cartesian(
      xlim = if (!all(is.na(xlim_vec))) xlim_vec else NULL,
      ylim = if (!all(is.na(ylim_vec))) ylim_vec else NULL
    )
  }

  return(p)
}


#####--------------- Example run ---------------#####

# # Example 1: Simple plot with raw data grouped by VAR
# p1 <- sg_sim_tp(ds_i = ds_sim, group_i = 'VAR', col_i = 'VAR')
#
# # Example 2: Plot with mean and SD bands
# p2 <- sg_sim_tp(
#   ds_i = ds_sim,
#   cent_i = 'mean',
#   vrns_i = 'SD',
#   col_i = 'VAR',
#   fill_i = 'VAR'
# )
#
# # Example 3: Plot with median and 5th-95th percentile bands
# p3 <- sg_sim_tp(
#   ds_i = ds_sim,
#   cent_i = 'median',
#   lperc_i = 0.05,
#   uperc_i = 0.95,
#   col_i = 'VAR',
#   fill_i = 'VAR',
#   wrap_i = '~VAR',
#   wrap_ncol = 2
# )
#
# # Example 4: Plot grouped by Dose with points and faceting
# p4 <- sg_sim_tp(
#   ds_i = ds_sim,
#   cent_i = 'mean',
#   vrns_i = 'SE',
#   col_i = 'Regimen',
#   fill_i = 'Regimen',
#   add_points = 2,
#   grid_i = '~VAR',
#   log_y = TRUE
# )
#
# # Example 5: Plot with custom axis limits and IQR
# p5 <- sg_sim_tp(
#   ds_i = ds_sim,
#   cent_i = 'median',
#   vrns_i = 'IQR',
#   col_i = 'CAT1',
#   fill_i = 'CAT1',
#   x_min = 0,
#   x_max = 100,
#   y_min = 0
# )

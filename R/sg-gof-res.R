## Author: Ugolkov Yaroslav, Anatoly Pokladyuk
## First created: 2025-09-05
## Description: sg-gof-res and its helper functions
## Keywords: SimuRg, sg-gof-res, goodness-of-fit
#' Create residual diagnostic plots (WRES, IWRES, CWRES, etc.)
#'
#' @description
#' Generates residual diagnostic plots either versus time/predictions
#' or as residual distributions. Supports faceting, covariate coloring,
#' quantile binning for continuous covariates, and optional smoothing.
#' @inheritParams sg_dummy
#' @param res_type Character. Name of residual column to plot (e.g., `"WRES"`, `"IWRES"`, `"CWRES"`)
#' @param vs_time Logical. If `TRUE`, plot residuals vs TIME, else vs predictions (IPRED/PRED)
#' @param addline Logical. Add connecting lines per subject. Default `TRUE`
#' @param dens Logical. If `TRUE`, plot histogram/density of residuals instead of scatter
#' @param levels_discrete Integer. Threshold for treating covariates as discrete. Default 10
#'
#' @return A `ggplot2` object
#'
#' @examples
#' \donttest{
#' # Basic residuals vs time
#' sg_gof_res(
#'   fpath_i = paste0("functions/nlme/3.4-sg-gof-res/simurg-object/Warfarin_PK.RData"),
#'   res_type = "IWRES",
#'   vs_time = TRUE
#' )
#' # IWRES dist
#' sg_gof_res(
#'   fpath_i = paste0("functions/nlme/3.4-sg-gof-res/simurg-object/Warfarin_PK.RData"),
#'   res_type = "IWRES",
#'   cov_cols = "SEX",
#'   dens = T,
#'   lab_x = "IWRES",
#'   facet_i = "SEX",
#'   addline = T) + labs(subtitle = "Individual")
#' }
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom scales pretty_breaks pretty_breaks trans_format math_format
#' @export
sg_gof_res <- function(
    fpath_i, res_type, cov_cols = NULL, indiv = TRUE, vs_time = TRUE,
    addline = TRUE, alpha_i = 0.5, smooth = TRUE, log_x = FALSE,
    sc_factor = 1, abreaks = scales::pretty_breaks(7),
    lab_y = NULL, lab_x = NULL, col_i = NULL, col_lab = NULL,
    facet_i = NULL, f_scales = "fixed",
    dens = FALSE, n_bins = 50,
    ymin = NA, ymax = NA, xmin = NA, xmax = NA, no_leg = FALSE,
    n_quantiles = 3, levels_discrete = 10
) {

  # --- helper functions ---
  is_discrete <- function(x, max_levels = levels_discrete) {
    n_unique <- length(unique(na.omit(x)))
    n_unique <= max_levels
  }
  continuous_to_categories <- function(x, n_quant = 3) {
    cut(
      x,
      breaks = quantile(x, probs = seq(0, 1, length.out = n_quant + 1), na.rm = TRUE),
      include.lowest = TRUE,
      labels = paste0("Q", 1:n_quant)
    )
  }

  # --- data prep ---
  X <- if (vs_time) "TIME" else if (indiv) "IPRED" else "PRED"
  smrg_obj <- read_smrg_obj(fpath_i)

  ds_i <- smrg_obj$SDTAB %>%
    filter(MDV != 1) %>%
    mutate_at(vars(TIME, IPRED, PRED, DV), function(s) s / sc_factor) %>%
    rename_at(vars(one_of(res_type)), ~"Y") %>%
    rename_at(vars(one_of(X)), ~"X")

  # --- covariates ---
  if (!is.null(cov_cols)) {
    suppressMessages({
      ds_covs <- left_join(smrg_obj$COTAB, smrg_obj$CATAB)
    })
    ds_covs_i <- ds_covs %>% select(ID, one_of(cov_cols)) %>% distinct()
    stopifnot(n_distinct(ds_covs_i$ID) == nrow(ds_covs_i))
    suppressMessages({
      ds_i <- ds_i %>% left_join(ds_covs_i, by = "ID")
    })

    if (!is_discrete(ds_i[[cov_cols]])) {
      ds_i[[cov_cols]] <- continuous_to_categories(ds_i[[cov_cols]], n_quant = n_quantiles)
    } else {
      ds_i[[cov_cols]] <- factor(ds_i[[cov_cols]])
    }
  }

  # --- build plot ---
  if (!dens) {
    # scatter vs time/pred
    p_char <- list(
      labs(y = lab_y, x = lab_x),
      geom_hline(yintercept = c(-2, 2), col = "firebrick", size = 0.5, linetype = "dashed"),
      geom_hline(yintercept = 0, col = "black", size = 0.7, linetype = "dashed"),
      scale_y_continuous(breaks = scales::pretty_breaks(7), limits = c(ymin, ymax)),
      scale_color_manual(values = rep(MSDcol, 20)),
      theme(legend.justification = c("left", "center"),
            legend.box.just = "left",
            legend.background = element_rect(fill = "white", size = 0.15, linetype = "solid", colour = "black"),
            legend.key.size = unit(0.38, "cm"),
            legend.title = element_text(size = 8),
            legend.text = element_text(size = 8),
            plot.title = element_text(size = 12))
    )
    if (log_x) {
      p_char <- c(p_char, scale_x_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ))
    } else {
      p_char <- c(p_char, scale_x_continuous(breaks = abreaks))
    }

    p_Res <- ggplot(ds_i, aes(x = X, y = Y)) + p_char

    if (!is.null(col_i)) {
      p_Res <- p_Res +
        geom_point(aes(col = !!sym(col_i)), size = 1.5, alpha = alpha_i) +
        labs(col = col_lab) +
        guides(color = guide_legend(override.aes = list(alpha = 1)))

      if (addline) p_Res <- p_Res + geom_line(aes(col = !!sym(col_i), group = ID), linewidth = 0.4, alpha = alpha_i)
      if (smooth)  p_Res <- p_Res + geom_smooth(aes(col = !!sym(col_i), group = !!sym(col_i)), formula = "y ~ x", method = "loess", linewidth = 1.2, se = FALSE)
    } else {
      p_Res <- p_Res +
        geom_point(col = MSDcol[1], size = 1.5, alpha = alpha_i)

      if (addline) p_Res <- p_Res + geom_line(aes(group = ID), col = MSDcol[1], linewidth = 0.4, alpha = alpha_i)
      if (smooth)  p_Res <- p_Res + geom_smooth(formula = "y ~ x", method = "loess", linewidth = 1.2, se = FALSE, col = MSDcol[3])
    }

  } else {
    # histogram/density of residuals
    p_char <- list(
      scale_y_continuous(name = "Density", breaks = scales::pretty_breaks(7), expand = c(0, 0), limits = c(0, ymax)),
      scale_x_continuous(name = lab_x, breaks = scales::pretty_breaks(7), expand = c(0, 0)),
      scale_fill_manual(values = rep(MSDcol, 20)),
      coord_cartesian(xlim = c(xmin, xmax)),
      theme(plot.title = element_text(size = 12))
    )

    if (!is.null(col_i)) {
      p_Res <- ggplot(ds_i, aes(x = Y, y = ..density..)) +
        geom_histogram(aes(fill = !!sym(col_i)), bins = n_bins, col = "grey25", alpha = alpha_i, position = "identity") +
        labs(fill = col_lab) +
        guides(fill = guide_legend(override.aes = list(alpha = 1))) +
        p_char
    } else {
      p_Res <- ggplot(ds_i, aes(x = Y, y = ..density..)) +
        geom_histogram(bins = n_bins, col = "grey25", fill = MSDcol[2]) +
        p_char
    }
    if (addline) {
      p_Res <- p_Res + annotate("line", x = seq(-4, 4, 0.01), y = dnorm(seq(-4, 4, 0.01)), size = 0.8, linetype = "dashed")
    }
  }

  # --- facets & legend ---
  if (!is.null(facet_i)) {
    p_Res <- p_Res + facet_wrap(as.formula(paste0("~", facet_i)), scales = f_scales)
  }
  if (no_leg) {
    p_Res <- p_Res + theme(legend.position = "none")
  }

  return(p_Res)
}

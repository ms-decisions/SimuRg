## Author: Ugolkov Yaroslav, Anatoly Pokladyuk
## First created: 2025-09-05
## Description: sg-gof-res and its helper functions
## Keywords: SimuRg, sg-gof-res, goodness-of-fit
#' Create residual diagnostic plots
#'
#' @description
#' Generates residual diagnostic plots versus time/predictions.
#' Supports faceting, covariate coloring, quantile binning for
#' continuous covariates, and optional smoothing.
#' @inheritParams sg_dummy
#' @param vs_time Logical. If `TRUE`, plot residuals vs TIME; otherwise, plot
#'  residuals vs predictions (IPRED/PRED).
#' @param weighted Logical. If `TRUE`, use weighted residuals; otherwise, use
#'  RES/IRES
#'
#' @return A `ggplot2` object
#'
#' @examples
#' # Basic example with mock data
#' set.seed(123)  # For reproducibility
#' n_subjects <- 50
#' mock_obj <- list(
#'   SDTAB = do.call(rbind, lapply(1:n_subjects, function(id) {
#'     n_obs <- 6
#'     times <- sort(runif(n_obs, min = 0, max = 24))  # random times between 0 and 24h
#'     data.frame(
#'       ID    = id,
#'       TIME  = times,
#'       DV    = rnorm(n_obs, mean = 10, sd = 2),
#'       PRED  = rnorm(n_obs, mean = 10, sd = 1.5),
#'       IPRED = rnorm(n_obs, mean = 10, sd = 1.2),
#'       IWRES = rnorm(n_obs, mean = 0, sd = 0.8),
#'       IRES  = rnorm(n_obs, mean = 0, sd = 1.2),
#'       MDV   = 0
#'     )
#'   })),
#'   COTAB = data.frame(
#'     ID  = 1:n_subjects,
#'     AGE = sample(20:80, n_subjects, replace = TRUE)
#'   ),
#'   CATAB = data.frame(
#'     ID   = 1:n_subjects,
#'     RACE = sample(c("Hispanic", "Asian", "Caucasian"), n_subjects, replace = TRUE)
#'   )
#' )
#'
#' # Basic plot: individual weighted residuals vs TIME (weighted = TRUE, vs_time = TRUE)
#' p <- sg_gof_res(mock_obj, smooth = FALSE)
#'
#' # With covariates and faceting (use RACE as facet and AGE as color)
#' p <- sg_gof_res(
#'   mock_obj,
#'   smooth = TRUE,
#'   cov_cols = c("RACE","AGE"),
#'   col_i = "RACE",
#'   facet_i = "AGE",
#'   indiv = TRUE,
#'   weighted = TRUE,
#'   vs_time = FALSE,
#'   legend_fl = TRUE
#' )
#' p
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom scales pretty_breaks pretty_breaks trans_format math_format
#' @export
sg_gof_res <- function(
    fpath_i, DVID = 1, cov_cols = NULL, indiv = TRUE, vs_time = TRUE, weighted = TRUE,
    addline = TRUE, alpha_i = 0.5, smooth = TRUE, log_x = FALSE,
    abreaks = scales::pretty_breaks(7), lab_y = NULL, lab_x = NULL,
    col_i = NULL, col_lab = NULL, facet_i = NULL, f_scales = "fixed",
    n_bins = 50, min_y = NA, max_y = NA, min_x = NA, max_x = NA, legend_fl = FALSE,
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

  if (is.null(smrg_obj$SDTAB)) {
    stop("sg_fit object must contain SDTAB component")
  }

  sdtab <- smrg_obj$SDTAB

  if (is.data.frame(sdtab) && nrow(sdtab) == 0) {
    stop("SDTAB is empty (no rows)")
  }
  if (is.list(sdtab) && length(sdtab) == 0) {
    stop("SDTAB is empty (no elements)")
  }
  if (is.data.frame(sdtab)) {
    ds_i <- sdtab
  } else if (is.list(sdtab)) {
    ds_i <- as.data.frame(do.call(rbind, sdtab))
  } else {
    stop("SDTAB must be a data frame or a list of data frames")
  }

  ds_i <- filter_sdtab_by_DVID(ds_i, DVID)

  if (indiv) {
    if (weighted) {
      # Weighted individual residuals
      if ("IWRES" %in% colnames(ds_i)) {
        res_type <- "IWRES"
      } else {
        stop("IWRES column not found in data for individual weighted residuals.")
      }
    } else {
      # Unweighted individual residuals
      if ("IRES" %in% colnames(ds_i)) {
        res_type <- "IRES"
      } else {
        stop("IRES column not found in data for individual residuals.")
      }
    }
  } else {
    # Population plots
    if (weighted) {
      # Weighted population residuals: CWRES preferred, then WRES
      if ("CWRES" %in% colnames(ds_i)) {
        res_type <- "CWRES"
      } else if ("WRES" %in% colnames(ds_i)) {
        res_type <- "WRES"
      } else {
        stop("No weighted population residual column (CWRES or WRES) found in data.")
      }
    } else {
      # Unweighted population residuals
      if ("RES" %in% colnames(ds_i)) {
        res_type <- "RES"
      } else {
        stop("RES column not found in data for population residuals.")
      }
    }
  }

  if ("MDV" %in% colnames(ds_i)) {
    ds_i <- ds_i %>% filter(.data$MDV != 1)
  }
  ds_i <- ds_i %>%
    rename_at(vars(any_of(res_type)), ~"Y") %>%
    rename_at(vars(any_of(X)), ~"X")

  # --- covariates ---
  if (!is.null(cov_cols)) {
    suppressMessages({
      ds_covs <- left_join(smrg_obj$COTAB, smrg_obj$CATAB)
    })
    ds_covs_i <- ds_covs %>% select(ID, one_of(cov_cols)) %>% unique()
    stopifnot(n_distinct(ds_covs_i$ID) == nrow(ds_covs_i))
    suppressMessages({
      ds_i <- ds_i %>% left_join(ds_covs_i, by = "ID")
    })

    for (col in cov_cols) {
      if (!is_discrete(ds_i[[col]])) {
        ds_i[[col]] <- continuous_to_categories(ds_i[[col]], n_quant = n_quantiles)
        #message(paste(col, "was converted into", n_quantiles, "quantile-based categories."))
      } else {
        ds_i[[col]] <- factor(ds_i[[col]])
      }
    }
  }

  # Compute y-axis limits (pass NA for auto limits)
  if (!all(is.na(c(min_y, max_y)))) {
    y_limits <- c(min_y, max_y)
    y_scale <- scale_y_continuous(breaks = abreaks, limits = y_limits)
  } else {
    y_scale <- scale_y_continuous(breaks = abreaks)
  }

  if (!all(is.na(c(min_x, max_x)))) {
    x_limits <- c(min_x, max_x)
  } else {
    x_limits <- NULL
  }

  y_label <- if (is.null(lab_y)) res_type else lab_y
  x_label <- if (is.null(lab_x)) X else lab_x

  # scatter vs time/pred
  p_char <- c(list(
    labs(y = y_label, x = x_label),
    #geom_hline(yintercept = c(-2, 2), col = "firebrick", linewidth = 0.5, linetype = "dashed"),
    geom_hline(yintercept = 0, col = "black", linewidth = 0.7, linetype = "dashed"),
    y_scale,
#    scale_color_manual(values = rep(MSDcol, 20)),
    theme(legend.justification = c("left", "center"),
          legend.box.just = "left",
          legend.background = element_rect(fill = "white", linewidth = 0.15,
                                           linetype = "solid", colour = "black"),
          legend.key.size = unit(0.38, "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          plot.title = element_text(size = 12),
          panel.grid.minor = element_blank())
    ),
  if (weighted) list(geom_hline(yintercept = c(-2, 2),
                              col = "firebrick", linewidth = 0.5, linetype = "dashed"))
  )
  if (log_x) {
    p_char <- c(p_char, scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = x_limits
    ))
  } else {
    p_char <- c(p_char, scale_x_continuous(breaks = abreaks, limits = x_limits))
  }
    p_Res <- ggplot(ds_i, aes(x = X, y = Y)) + p_char + theme_bw()

  if (!is.null(col_i)) {
    if (!is.factor(ds_i[[col_i]])) {
      ds_i[[col_i]] <- factor(ds_i[[col_i]])
    }
    col_levels <- levels(ds_i[[col_i]])
    n_cols <- length(col_levels)

    col_values <- rep(MSDcol, length.out = n_cols)

    col_labels <- paste0(col_i, ": ", col_levels)

    # Add the custom color scale with formatted labels
    p_Res <- p_Res +
      scale_color_manual(values = col_values, labels = col_labels)
  }

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

  # --- facets & legend ---
  if (!is.null(facet_i)) {
    p_Res <- p_Res + facet_wrap(as.formula(paste0("~", facet_i)), scales = f_scales, labeller = label_both)
  }
  if (!legend_fl) {
    p_Res <- p_Res + theme(legend.position = "none")
  }
  return(p_Res)
}

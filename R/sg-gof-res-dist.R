## Author: Ugolkov Yaroslav
## First created: 2025-09-05
## Description: sg-gof-res-dist and its helper functions
## Keywords: SimuRg, sg-gof-res-dist, goodness-of-fit
#'
#' Plot distribution or QQ-plot of residuals from NONMEM/Simurg run
#'@description
#' This function visualizes residuals from [GFO] `$SDTAB`. The residuals can be plotted either as
#' histograms with overlaid normal density curves, or as QQ-plots to assess
#' normality. Pass `cov_cols` and `facet_i` to draw one panel per covariate
#' level (for example a population column in `CATAB`).
#'
#' @inheritParams sg_dummy
#' @param res_type Character vector. One or several types of residuals to plot.
#'   Values must correspond to column name(s) in `SDTAB`, for example
#'   `"RES"`, `"WRES"`, `"CWRES"`, `"IWRES"`, `"IRES"`. If multiple residual
#'   types are provided, plots will be generated for each of them.
#' @param ndist Logical. If `TRUE`, overlays the corresponding normal density
#'   curve on the histogram. Default value is `TRUE`.
#' @param plot_type Character. Type of plot to produce:
#'   * `"DIST"` - histogram of residuals, with a normal density overlay when `ndist` is `TRUE`,
#'   * `"QQ"` - QQ-plot of residuals.
#'   Default value is `"DIST"`.
#' @param f_scales String, one of `"fixed"`, `"free"`, `"free_x"`, `"free_y"`.
#'   Scales of the facet panels. Default value is `"free"` when `facet_i` is
#'   `NULL`, and `"fixed"` when `facet_i` is set, so densities can be compared
#'   across covariate levels.
#' @param lab_x String. X-axis label. Default value is `NULL`, which uses the
#'   residual name when `res_type` has one element and `"Residuals"` otherwise.
#'   Ignored when `plot_type` is `"QQ"`.
#' @param n_quantiles Integer. Number of quantile groups for a continuous
#'   covariate in `cov_cols` or `facet_i`. Facet labels show the group and its
#'   range. A square bracket means that endpoint is included, and a round
#'   bracket means it is excluded. The lowest group includes its left
#'   endpoint, for example `Q1 [20, 36.67]`. Each following group excludes its
#'   left endpoint and includes its right endpoint, for example
#'   `Q2 (36.67, 53.33]`, so a value on a shared boundary belongs to the
#'   higher group. Default value is `3`.
#'
#' @return A `ggplot` object.
#' @examples
#' \donttest{
#' # Plot the distribution of individual weighted residuals (IWRES)
#' # (default `DVID = 1`; pass `DVID` or a `DVNAME` string to select another endpoint)
#' fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData",
#'                         package = "SimuRg")
#' p <- sg_gof_res_dist(fpath_i = fpath_i, res_type = "IWRES")
#' p
#'
#' # One panel per level of a categorical covariate
#' p <- sg_gof_res_dist(
#'   fpath_i = fpath_i,
#'   res_type = "IWRES",
#'   cov_cols = "SEX",
#'   facet_i = "SEX"
#' )
#' p
#'
#' # QQ-plots for the same default endpoint
#' p <- sg_gof_res_dist(fpath_i = fpath_i, res_type = "RES", plot_type = "QQ")
#' p
#'
#' # Several residual columns at once (still one `DVID` unless you change it)
#' p <- sg_gof_res_dist(fpath_i = fpath_i, res_type = c("IWRES", "IRES"))
#' p
#' }
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @importFrom scales pretty_breaks
#' @export


sg_gof_res_dist <- function(fpath_i, res_type = "RES", DVID = 1, n_bins = 30,
                            ndist = TRUE, plot_type = "DIST",
                            cov_cols = NULL, facet_i = NULL, f_scales = NULL,
                            lab_x = NULL, n_quantiles = 3, levels_discrete = 10) {
  is_discrete <- function(x, max_levels = levels_discrete) {
    if (is.character(x) || is.factor(x) || is.logical(x)) {
      return(TRUE)
    }
    length(unique(stats::na.omit(x))) <= max_levels
  }

  continuous_to_categories <- function(x, n_quant = 3) {
    brks <- unname(stats::quantile(
      x,
      probs = seq(0, 1, length.out = n_quant + 1),
      na.rm = TRUE
    ))
    lo <- formatC(brks[seq_len(n_quant)], format = "g")
    hi <- formatC(brks[2:(n_quant + 1)], format = "g")
    open <- c("[", rep("(", n_quant - 1L))
    labels <- paste0("Q", seq_len(n_quant), " ", open, lo, ", ", hi, "]")
    cut(
      x,
      breaks = brks,
      include.lowest = TRUE,
      labels = labels
    )
  }

  as_table <- function(x, name) {
    if (is.null(x)) {
      return(NULL)
    }
    if (is.data.frame(x)) {
      return(x)
    }
    if (is.list(x) && length(x) > 0) {
      return(as.data.frame(do.call(rbind, x)))
    }
    stop(name, " must be a data frame or a list of data frames")
  }

  loaded <- read_smrg_obj(fpath_i)
  if (is.null(loaded$SDTAB) && !is.null(loaded$GFO)) {
    smrg_obj <- loaded$GFO
  } else {
    smrg_obj <- loaded
  }

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
    sdtab_i <- sdtab
  } else if (is.list(sdtab)) {
    sdtab_i <- as.data.frame(do.call(rbind, sdtab))
  } else {
    stop("SDTAB must be a data frame or a list of data frames")
  }

  sdtab_i <- filter_sdtab_by_DVID(sdtab_i, DVID)

  required_cols <- unique(c("DVID", res_type))
  missing_cols <- setdiff(required_cols, colnames(sdtab_i))
  if (length(missing_cols) > 0) {
    stop("SDTAB is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  sdtab_i <- sdtab_i %>%
    mutate(
      DVID = if ("DVID" %in% colnames(.)) as.numeric(DVID) else NA_real_
    )

  for (res_col in res_type) {
    if (res_col %in% colnames(sdtab_i)) {
      sdtab_i <- sdtab_i %>% mutate(across(all_of(res_col), as.numeric))
    }
  }

  facet_cols <- if (is.null(facet_i)) character() else as.character(facet_i)
  cov_requested <- if (is.null(cov_cols)) character() else as.character(cov_cols)
  need_cols <- unique(c(cov_requested, facet_cols))
  missing_in_sdtab <- setdiff(need_cols, colnames(sdtab_i))

  if (length(missing_in_sdtab) > 0) {
    if (!"ID" %in% colnames(sdtab_i)) {
      stop(
        "SDTAB must contain ID to join covariates: ",
        paste(missing_in_sdtab, collapse = ", ")
      )
    }
    cotab <- as_table(smrg_obj$COTAB, "COTAB")
    catab <- as_table(smrg_obj$CATAB, "CATAB")
    if (is.null(cotab) && is.null(catab)) {
      stop(
        "Covariate columns not found in SDTAB, and COTAB/CATAB are missing: ",
        paste(missing_in_sdtab, collapse = ", ")
      )
    }
    if (is.null(cotab)) {
      ds_covs <- catab
    } else if (is.null(catab)) {
      ds_covs <- cotab
    } else {
      suppressMessages({
        ds_covs <- left_join(cotab, catab)
      })
    }
    still_missing <- setdiff(missing_in_sdtab, colnames(ds_covs))
    if (length(still_missing) > 0) {
      stop(
        "Covariate columns not found in SDTAB, COTAB, or CATAB: ",
        paste(still_missing, collapse = ", ")
      )
    }
    ds_covs_i <- ds_covs %>%
      select(ID, all_of(missing_in_sdtab)) %>%
      distinct()
    if (n_distinct(ds_covs_i$ID) != nrow(ds_covs_i)) {
      stop("Covariate tables must contain one row per ID")
    }
    suppressMessages({
      sdtab_i <- sdtab_i %>% left_join(ds_covs_i, by = "ID")
    })
  }

  for (col in need_cols) {
    if (!is_discrete(sdtab_i[[col]])) {
      sdtab_i[[col]] <- continuous_to_categories(sdtab_i[[col]], n_quant = n_quantiles)
    } else {
      sdtab_i[[col]] <- factor(
        sdtab_i[[col]],
        levels = unique(stats::na.omit(sdtab_i[[col]]))
      )
    }
  }

  if (is.null(f_scales)) {
    f_scales <- if (length(facet_cols) == 0) "free" else "fixed"
  }
  f_scales <- match.arg(f_scales, c("fixed", "free", "free_x", "free_y"))

  if (length(facet_cols) == 0) {
    panel_formula <- DVID ~ residual_type
  } else {
    panel_vars <- facet_cols
    if (length(res_type) > 1) {
      panel_vars <- c(panel_vars, "residual_type")
    }
    panel_formula <- stats::as.formula(
      paste("~", paste(panel_vars, collapse = " + "))
    )
  }

  group_cols <- unique(c("DVID", facet_cols))
  res_for_plot2 <- sdtab_i %>%
    select(all_of(unique(c(group_cols, res_type)))) %>%
    pivot_longer(
      cols = all_of(res_type),
      names_to = "residual_type",
      values_to = "value"
    )

  x_name <- if (!is.null(lab_x)) {
    lab_x
  } else if (length(res_type) == 1L) {
    res_type
  } else {
    "Residuals"
  }

  p_i <- ggplot(data = res_for_plot2, aes(x = value, y = after_stat(density))) +
    geom_histogram(bins = n_bins, col = "grey25", fill = MSDcol[2]) +
    facet_wrap(panel_formula, scales = f_scales) +
    scale_y_continuous(name = "Density", breaks = scales::pretty_breaks(7), expand = c(0, 0), limits = c(0, NA)) +
    scale_x_continuous(name = x_name, breaks = scales::pretty_breaks(7), expand = c(0, 0)) +
    theme(panel.grid.minor = element_blank()) +
    theme_bw()

  if (ndist) {
    norm_group <- unique(c(group_cols, "residual_type"))
    norm_params <- res_for_plot2 %>%
      group_by(across(all_of(norm_group))) %>%
      summarise(
        mean_val = mean(value, na.rm = TRUE),
        sd_val = sd(value, na.rm = TRUE),
        min_val = min(value, na.rm = TRUE),
        max_val = max(value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      filter(is.finite(sd_val), sd_val > 0, is.finite(min_val), is.finite(max_val))

    if (nrow(norm_params) > 0) {
      norm_curves <- norm_params %>%
        group_by(across(all_of(norm_group))) %>%
        summarise(
          x = list(seq(min_val - 0.5 * sd_val, max_val + 0.5 * sd_val, length.out = 1000)),
          y = list(dnorm(x[[1]], mean = mean_val, sd = sd_val)),
          .groups = "drop"
        ) %>%
        unnest(c(x, y))

      p_i <- p_i + geom_line(
        data = norm_curves,
        aes(x = x, y = y),
        color = "black",
        linewidth = 0.8,
        linetype = "dashed"
      )
    }
  }

  if (plot_type == "QQ") {
    p_i <- ggplot(res_for_plot2, aes(sample = value)) +
      stat_qq(size = 1.75, color = MSDcol[1], alpha = 0.8) +
      stat_qq_line(col = "firebrick") +
      labs(x = "Theoretical quantiles", y = "Sample quantiles") +
      geom_abline(linewidth = 0.5, col = "black", linetype = "dashed") +
      scale_x_continuous(breaks = scales::pretty_breaks(7)) +
      scale_y_continuous(breaks = scales::pretty_breaks(7)) +
      theme(panel.grid.minor = element_blank()) +
      theme_bw() +
      facet_wrap(panel_formula, scales = f_scales)
  }

  p_i
}

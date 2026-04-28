## Author: Mikhailova Anna
## First created: 2025-08-25
## Description: sg-gof-tp and its helper functions
## Keywords: SimuRg, sg-gof-tp, goodness-of-fit

#' Plot time profiles of the fitted data
#'
#' @inheritParams sg_dummy
#' @param pop Logical. If `TRUE` (default), adds population predictions (PRED) as dashed lines and points.
#' @param cap String. Plot caption. Default is "empty circles - observed data solid lines with point - individual predictions dashed grey lines with point - population predictions"
#' @param lab_x String. X-ax label. Default is "Time since first dose, h"
#' @param lab_y String. Y-ax label. Default is "Plasma concentration, mmol/L"
#' @param sort_by Character vector of column names to sort ID facets by.
#'        Use "DOSE" to sort by last dose. Other names (e.g. covariates) must exist in SimuRg object.
#' @param desc Logical. If TRUE, sort in descending order for all columns in `sort_by`.
#' @returns A list of plots with predicted time profiles, faceted by id
#' @examples
#' \donttest{
#' fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData",
#'                         package = "SimuRg")
#' plot_list <- sg_gof_tp(fpath_i)
#' print(plot_list[[1]])
#' }
#' @import dplyr
#' @import ggplot2
#' @importFrom scales pretty_breaks trans_format
#' @importFrom forcats fct_inorder
#' @importFrom zoo na.locf
#' @importFrom purrr map
#' @import stringr
#' @export
sg_gof_tp <- function(fpath_i, pop = T, filt = "T", cov_cols = NULL, col_i = NULL,
                      DVID = 1, tsld = F, f_scales = "free", sort_by = NULL,
                      desc = F, log_y = F, lab_x = "Time since first dose, h",
                      lab_y = "Plasma concentration, mmol/L",
                      cap = str_c("empty circles - observed data\n",
                                  "solid lines with point - ",
                                  "individual predictions\n",
                                  "dashed grey lines with ",
                                  "point - population predictions"),
                      n_quantiles = 3, levels_discrete = 10){

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

  if (!is_null(read_smrg_obj(fpath_i)$GFO))
    smrg_obj <- read_smrg_obj(fpath_i)$GFO
  else
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

  ds_tp_pre <- ds_i
  if ("MDV" %in% colnames(ds_i)) {
    ds_tp_pre <- ds_tp_pre %>% filter(.data$MDV != 1)
  }

  continuous_factors <- list()

  if (!is.null(cov_cols)) {
    suppressMessages({
      ds_covs <- left_join(smrg_obj$COTAB, smrg_obj$CATAB)
    })
    ds_covs_i <- ds_covs %>% select(ID, one_of(cov_cols)) %>% unique()
    stopifnot(n_distinct(ds_covs_i$ID) == nrow(ds_covs_i))
    suppressMessages({
      ds_tp_pre <- ds_tp_pre %>% left_join(ds_covs_i, by = "ID")
    })

    for (col in cov_cols) {
      if (!is_discrete(ds_tp_pre[[col]])) {
        cat_col <- paste0(col, "_cat")
        ds_tp_pre[[cat_col]] <- continuous_to_categories(ds_tp_pre[[col]], n_quant = n_quantiles)
        #message(paste(col, "was converted into", n_quantiles, "quantile-based categories."))
        continuous_factors[[col]] <- cat_col
      } else {
        ds_tp_pre[[col]] <- factor(ds_tp_pre[[col]])
      }
    }
  }

  if (!is.null(sort_by) && "DOSE" %in% sort_by) {
    if (is.null(smrg_obj$EVTAB)) {
      stop("EVTAB not found in smrg_obj; cannot sort by DOSE.")
    }
    evtab <- smrg_obj$EVTAB %>%
      filter(EVID == 1) %>%
      group_by(ID) %>%
      arrange(TIME) %>%
      summarise(DOSE = first(AMT), .groups = "drop")
    ds_tp_pre <- ds_tp_pre %>% left_join(evtab, by = "ID")
  }

  ds_tp_pre <- ds_tp_pre %>% filter(eval(rlang::parse_expr(filt)))

  if (!is.null(col_i)) {
    actual_col <- if (col_i %in% names(continuous_factors)) continuous_factors[[col_i]] else col_i
    if (!is.factor(ds_tp_pre[[actual_col]])) {
      ds_tp_pre[[actual_col]] <- factor(ds_tp_pre[[actual_col]])
    }
    col_levels <- levels(ds_tp_pre[[actual_col]])
    n_cols <- length(col_levels)
    col_values <- rep(MSDcol, length.out = n_cols)   # or use a better color palette
    col_labels <- paste0(col_i, ": ", col_levels)
    color_column <- actual_col
  }

  if (!is.null(sort_by) && length(sort_by) > 0) {
    # Check that all sort_by columns exist
    missing <- setdiff(sort_by, colnames(ds_tp_pre))
    if (length(missing) > 0) {
      stop("sort_by column(s) not found: ", paste(missing, collapse = ", "))
    }

    # Obtain desired ID order
    if (desc) {
      id_order <- ds_tp_pre %>%
        distinct(ID, !!!syms(sort_by)) %>%
        arrange(across(all_of(sort_by), desc)) %>%
        pull(ID)
    } else {
      id_order <- ds_tp_pre %>%
        distinct(ID, !!!syms(sort_by)) %>%
        arrange(across(all_of(sort_by))) %>%
        pull(ID)
    }

    # Reorder ID factor in the whole dataset
    ds_tp_pre <- ds_tp_pre %>%
      mutate(ID = factor(ID, levels = id_order))

    # Build facet labels: "ID <id>; col1=val1; col2=val2"
    label_lookup <- ds_tp_pre %>%
      distinct(ID, !!!syms(sort_by)) %>%
      rowwise() %>%
      mutate(
        label_part = paste0(sort_by, "=", c_across(all_of(sort_by)), collapse = "; "),
        facet_label = paste0("ID ", ID, "; ", label_part)
      ) %>%
      ungroup() %>%
      select(ID, facet_label) %>%
      mutate(facet_label = factor(facet_label, levels = unique(facet_label[order(ID)])))

    ds_tp_pre <- ds_tp_pre %>%
      left_join(label_lookup, by = "ID")
  } else {
    ds_tp_pre <- ds_tp_pre %>%
      mutate(ID = fct_inorder(as.character(ID))) %>%
      mutate(facet_label = as.character(ID)) %>%
      mutate(facet_label = factor(facet_label, levels = levels(ID)))
  }

  if ((tsld) & !("ATSLD" %in% colnames(ds_i))) {
    stop("No column specified for time since last dose")
  } else if (tsld) {
    ds_tp_pre <- ds_tp_pre %>%
      mutate(TSLD = ATSLD) %>% group_by(ID) %>% # TSLD calculations should be implemented
      mutate(OCCMAN = ifelse(lag(TSLD) > TSLD | is.na(lag(TSLD)),
                             row_number(), NA), OCCMAN = zoo::na.locf(OCCMAN)) %>% ungroup() %>%
      mutate(DAYMAN = ceiling(TIME/24), Period = case_when(DAYMAN == 1 ~ "Day 1",
                                                           DAYMAN == 2 ~ "Day 2",
                                                           DAYMAN == 3 ~ "Day 3",
                                                           DAYMAN > 3 ~ "Day 3+"))
  }

  ds_tp_filt <- ds_tp_pre %>%
    arrange(ID) %>%
    mutate(NTILE = ceiling(as.numeric(ID) / 20))
    # arrange(DVID) %>%
    # mutate(NTILE = ceiling(as.numeric(fct_inorder(as.character(ID)))/20),
    #        ID = fct_inorder(as.character(ID)))

  TVAR <- "TIME";
  if(tsld){ TVAR <- "TSLD";}

  p_i <- unique(ds_tp_filt$NTILE) %>% map(function(n){
    set_n <- ds_tp_filt %>% filter(NTILE == n)

    p_n <- ggplot(data = set_n, aes(x = !!sym(TVAR), y = DV)) + # , group = OCCMAN)) +
      geom_point(size = 1.5, shape = 21, stroke = 0.8, fill = "white", col = "grey15")

    if(pop) {   # add PRED layers only when pop is TRUE
      p_n <- p_n +
        geom_line(aes(y = PRED), lty = "dashed", col = "grey75") +
        geom_point(aes(y = PRED), size = 1, col = "grey75")
    }

    if(!is.null(col_i)) {
      p_n <- p_n +
        geom_line(aes(y = IPRED, color = .data[[color_column]])) +
        geom_point(aes(y = IPRED, color = .data[[color_column]]), size = 1)
    } else {
      p_n <- p_n +
        geom_line(aes(y = IPRED)) +
        geom_point(aes(y = IPRED), size = 1)
    }

    # Rest of the layers (scales, facets, theme)
    p_n <- p_n +
      scale_x_continuous(name = lab_x, breaks = scales::pretty_breaks(7)) +
      scale_y_continuous(name = lab_y, breaks = scales::pretty_breaks(7))

    # Add manual color scale if needed
    if(!is.null(col_i)) {
      # col_values and col_labels should have been prepared earlier (as in your code)
      p_n <- p_n + scale_color_manual(
        name = col_i,
        values = col_values,
        labels = col_levels,
        drop = FALSE
      )
      # Show legend, hide grid minor
      p_n <- p_n + theme(legend.position = "right", panel.grid.minor = element_blank())
    } else {
      # Original scale_color_manual? You had scale_color_manual(values = MSDcol) but no mapping – remove it
      # Actually, better to keep legend none and no color scale
      p_n <- p_n + theme(legend.position = "none", panel.grid.minor = element_blank())
    }

    p_n <- p_n +
      facet_wrap(~facet_label, scales = f_scales) +
      theme_bw() +
      labs(subtitle = str_c(" (part ", n, " out of ", max(ds_tp_filt$NTILE), ")"),
           caption = cap)

    # p_n <- p_n +
    #   geom_line(aes(y = IPRED)) +
    #   geom_point(aes(y = IPRED), size = 1) + #, shape = Period
    #   scale_x_continuous(name = lab_x, breaks = scales::pretty_breaks(7)) +
    #   scale_y_continuous(name = lab_y, breaks = scales::pretty_breaks(7)) +
    #   scale_color_manual(values = MSDcol) +
    #   facet_wrap(~ID, scales = f_scales) +
    #   theme_bw() +
    #   theme(legend.position = "none", panel.grid.minor = element_blank()) +
    #   labs(subtitle = str_c( " (part ", n, " out of ", max(ds_tp_filt$NTILE), ")"),
    #        caption = cap)
    if(log_y){
      p_n <- p_n + scale_y_log10(
        name = lab_y,
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      )
    }
    return(p_n)
  })
  return(p_i)
}

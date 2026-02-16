## Author: Ugolkov Yaroslav
## First created: 2025-09-05
## Description: sg-gof-obpr and its helper functions
## Keywords: SimuRg, sg-gof-obpr, goodness-of-fit
#'
#'Observed vs predicted plot
#'
#'
#' @description
#' Function generates observed versus predicted (OBS vs PRED/IPRED) scatter plots,
#' a fundamental goodness-of-fit diagnostic tool in pharmacometric modeling.
#' This visualization assesses model adequacy by comparing observed clinical measurements against model-predicted values,
#' enabling identification of systematic bias, heteroscedasticity, and model misspecification patterns.
#' Function contains options for faceting, coloring by covariates, and trend lines.
#'
#' @inheritParams sg_dummy
#' @param lab_x X-axis label. Default "Model-predicted values"
#' @param lab_y Y-axis label. Default "Observed values"
#'
#' @return A ggplot2 object
#'

#' @examples
#' \donttest{
#' # Basic example with mock data
#' set.seed(123)  # For reproducibility
#' mock_obj <- list(
#'   SDTAB = data.frame(
#'     ID = rep(1:3, each = 5),
#'     TIME = rep(c(0, 1, 2, 4, 8), 3),
#'     DV = rnorm(15, mean = 10, sd = 2),
#'     PRED = rnorm(15, mean = 10, sd = 1.5),
#'     IPRED = rnorm(15, mean = 10, sd = 1.2),
#'     MDV = rep(0, 15)
#'   ),
#'   COTAB = data.frame(ID = 1:3, AGE = c(30, 45, 60)),
#'   CATAB = data.frame(ID = 1:3, RACE = c("Hispanic", "Hispanic", "Asian"))
#' )
#'
#' # Basic plot
#' p <- sg_gof_obpr(mock_obj)
#'
#' # With covariates and faceting
#' p <- sg_gof_obpr(
#'   mock_obj,
#'   cov_cols = "RACE",
#'   col_i = "RACE",
#'   facet_i = "RACE"
#' )
#' }
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom jsonlite fromJSON
#' @importFrom scales pretty_breaks trans_breaks trans_format math_format number_format
#' @export

sg_gof_obpr <- function(
    fpath_i, cov_cols = NULL, indiv = T, addline = T, alpha_i = 0.5,
    smooth = T, log_axes = F, sc_factor = 1, abreaks = scales::pretty_breaks(7),
    xlab = "Model-predicted values", ylab = "Observed values", col_i = NULL, col_lab = NULL,
    facet_i = NULL, f_scales = "fixed",
    no_leg = F, n_quantiles = 3, levels_discrete = 10
){

  #browser()
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

  X <- ifelse(indiv, "IPRED", "PRED")
  smrg_obj <- read_smrg_obj(fpath_i) #get(load(fpath_i))

  # Validate that smrg_obj has SDTAB
  if (is.null(smrg_obj$SDTAB)) {
    stop("sg_fit object must contain SDTAB component")
  }

  sdtab <- smrg_obj$SDTAB

  # Check if SDTAB is empty
  if (is.data.frame(sdtab) && nrow(sdtab) == 0) {
    stop("SDTAB is empty (no rows)")
  }
  if (is.list(sdtab) && length(sdtab) == 0) {
    stop("SDTAB is empty (no elements)")
  }

  # ds_i <- as.data.frame(do.call(rbind, sdtab)) %>%
  #   mutate(IPRED = as.numeric(IPRED), PRED = as.numeric(PRED),
  #          TIME = as.numeric(TIME), WRES = as.numeric(WRES),
  #          IWRES = as.numeric(IWRES), DV = as.numeric(DV),
  #          MDV = as.numeric(MDV)) %>%
  #   #read_table(fpath_i, skip = 1, col_names = T, col_types = cols()) %>%
  #   filter(MDV != 1) %>%
  #   mutate_at(vars(IPRED, PRED, DV), function(s){s/sc_factor}) %>%
  #   rename_at(vars(one_of(X)), function(n){n = "X"})



  # Handle SDTAB as either a data frame or a list of data frames
  if (is.data.frame(sdtab)) {
    ds_i <- sdtab
  } else if (is.list(sdtab)) {
    ds_i <- as.data.frame(do.call(rbind, sdtab))
  } else {
    stop("SDTAB must be a data frame or a list of data frames")
  }

  # Check for required columns and convert to numeric
  required_cols <- c("IPRED", "PRED", "TIME", "DV", "MDV")
  missing_cols <- setdiff(required_cols, colnames(ds_i))
  if (length(missing_cols) > 0) {
    stop("SDTAB is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Convert columns to numeric (only if they exist)
  ds_i <- ds_i %>%
    mutate(
      IPRED = if ("IPRED" %in% colnames(.)) as.numeric(IPRED) else NA_real_,
      PRED = if ("PRED" %in% colnames(.)) as.numeric(PRED) else NA_real_,
      TIME = as.numeric(TIME),
      DV = as.numeric(DV),
      MDV = as.numeric(MDV)
    )

  # Add optional columns if they exist
  if ("WRES" %in% colnames(ds_i)) {
    ds_i <- ds_i %>% mutate(WRES = as.numeric(WRES))
  }
  if ("IWRES" %in% colnames(ds_i)) {
    ds_i <- ds_i %>% mutate(IWRES = as.numeric(IWRES))
  }

  ds_i <- ds_i %>%
    filter(MDV != 1)

  # Apply scaling factor to IPRED, PRED, and DV if they exist
  scale_cols <- intersect(c("IPRED", "PRED", "DV"), colnames(ds_i))
  if (length(scale_cols) > 0) {
    ds_i <- ds_i %>%
      mutate(across(all_of(scale_cols), ~ .x / sc_factor))
  }

  # Rename the selected prediction column to "X"
  if (X %in% colnames(ds_i)) {
    ds_i <- ds_i %>% rename(X = all_of(X))
  } else {
    stop("Column '", X, "' not found in SDTAB. Available columns: ", paste(colnames(ds_i), collapse = ", "))
  }

  if(!is.null(cov_cols)){
    suppressMessages({
      ds_covs <- left_join(smrg_obj$COTAB,smrg_obj$CATAB)})
    ds_covs_i <- ds_covs %>% select(ID, one_of(cov_cols)) %>% unique()
    stopifnot(n_distinct(ds_covs_i$ID) == nrow(ds_covs_i))
    suppressMessages({
      ds_i <- ds_i %>% left_join(ds_covs_i, by = "ID")})
    if (!is_discrete(ds_i[[cov_cols]])) {
      ds_i[[cov_cols]] <- continuous_to_categories(ds_i[[cov_cols]], n_quant = n_quantiles)
      #message(paste(cov_cols, "was converted into", n_quantiles, "quantile-based categories."))
    } else {
      ds_i[[cov_cols]] <- factor(ds_i[[cov_cols]])
    }
  }

  lim_obpr <- c(min(ds_i$X, ds_i$DV), max(ds_i$X, ds_i$DV))

  p_char <- list(
    labs(y = ylab, x = xlab),
    geom_abline(linewidth = 0.5, col = "black", linetype = "dashed"),
    scale_color_manual(values = rep(MSDcol[c(2, 3, 4, 5, 7, 9, 6, 8, 1)], 30)),
    theme(legend.justification = c("left", "center"),
          legend.box.just = "left",
          legend.background = element_rect(fill = "white", linewidth = 0.15,
                                           linetype = "solid", colour = "black"),
          legend.key.size = unit(0.38, "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          plot.title = element_text(size = 12),
          panel.grid.minor = element_blank())
  )

  if(log_axes){
    p_char <- c(p_char,
                scale_x_log10(
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))),
                scale_y_log10(
                  breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))))
  } else {
    p_char <- c(p_char,
                scale_x_continuous(limits = lim_obpr, breaks = abreaks,
                                   labels = scales::number_format()),
                scale_y_continuous(limits = lim_obpr, breaks = abreaks,
                                   labels = scales::number_format()))
  }

  p_ObPr <- ggplot(data = ds_i, aes(x = X, y = DV)) + p_char + theme_bw()

  if(!is.null(col_i)){
    p_ObPr <- p_ObPr +
      geom_point(aes(col = !!sym(col_i)), size = 1.5, alpha = alpha_i) +
      guides(color = guide_legend(override.aes = list(alpha = 1))) +
      labs(col = col_lab)

    if(addline){ p_ObPr <- p_ObPr + geom_line(aes(col = !!sym(col_i), group = ID), lwd = 0.4, alpha = alpha_i) }
    if(smooth){ p_ObPr <- p_ObPr + geom_smooth(aes(col = !!sym(col_i), group = !!sym(col_i)),
                                               formula = "y ~ x", method = "loess",
                                               linewidth = 1.2, se = F) }

  } else {
    p_ObPr <- p_ObPr +
      geom_point(size = 1.5, alpha = alpha_i) +
      guides(color = guide_legend(override.aes = list(alpha = 1))) +
      labs(col = col_lab)

    if(addline){ p_ObPr <- p_ObPr + geom_line(aes(group = ID), col = MSDcol[1], lwd = 0.4, alpha = alpha_i) }
    if(smooth){ p_ObPr <- p_ObPr + geom_smooth(formula = "y ~ x", method = "loess",
                                               linewidth = 1.2, se = F, col = MSDcol[3]) }
  }

  if(!is.null(facet_i)){
    p_ObPr <- p_ObPr + facet_wrap(as.formula(str_c("~", facet_i)), scales = f_scales)
  }

  if(no_leg){
    p_ObPr <- p_ObPr + theme(legend.position = "none")
  }

  return(p_ObPr)
}


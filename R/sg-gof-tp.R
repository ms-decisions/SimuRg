## Author: Mikhailova Anna
## First created: 2025-08-25
## Description: sg-gof-tp and its helper functions
## Keywords: SimuRg, sg-gof-tp, foodness-of-fit

#' Plot time profiles of the fitted data
#'
#' @inheritParams sg_dummy
#' @param cap string. Plot caption. Default is "empty circles - observed data\n solid lines with point - individual predictions\n dashed grey lines with point - population predictions"
#' @param lab_x string. X-ax label. Default is "Time since first dose, h"
#' @param lab_y string. Y-ax label. Default is "Plasma concentration, mmol/L"
#' @returns A list of plots with predicted time profiles, faceted by id
#' @examples
#' \dontrun{
#' fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg"))
#' plot_list <- sg_gof_tp(fpath_i)
#' print(plot_list[[1]])
#' }
#' @import dplyr
#' @import ggplot2
#' @importFrom scales pretty_breaks trans_format
#' @importFrom forcats fct_inorder
#' @importFrom purrr map
#' @import stringr
#' @export
sg_gof_tp <- function(fpath_i, filt = "T",
                      tsld = F, f_scales = "free",
                      log_y = F, lab_x = "Time since first dose, h",
                      lab_y = "Plasma concentration, mmol/L",
                      cap = "empty circles - observed data\nsolid lines with point - individual predictions\ndashed grey lines with point - population predictions"){
  if (inherits(fpath_i, "character") ) {
    if (file.exists(fpath_i)) {
      obj <- get(load(fpath_i))
    } else {
      stop("File specified by fpath_i does not exist")
    }

  } else if (inherits(fpath_i, "list")) {
    obj <- fpath_i
  } else {
    stop("fpath_i object should be either an sg_fit object, or a path to saved sg_fit object")
  }
  ds_tp_pre <- obj1$SDTAB %>%
    filter(MDV != 1) %>% filter(eval(rlang::parse_expr(filt)))

  if ((tsld) & !("ATSLD" %in% colnames(obj1$SDTAB))) {
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
    arrange(DVID) %>%
    mutate(NTILE = ceiling(as.numeric(fct_inorder(as.character(ID)))/20),
           ID = fct_inorder(as.character(ID)))

  TVAR <- "TIME";
  if(tsld){ TVAR <- "TSLD";}

  p_i <- unique(ds_tp_filt$NTILE) %>% map(function(n){
    set_n <- ds_tp_filt %>% filter(NTILE == n)

    p_n <- ggplot(data = set_n, aes(x = !!sym(TVAR), y = DV)) + # , group = OCCMAN)) +
      geom_point(size = 1.5, shape = 21, stroke = 0.8, fill = "white", col = "grey15") +
      geom_line(aes(y = PRED), lty = "dashed", col = "grey75") +
      geom_point(aes(y = PRED), size = 1, col = "grey75") + #, shape = Period
      geom_line(aes(y = IPRED)) +
      geom_point(aes(y = IPRED), size = 1) + #, shape = Period
      scale_x_continuous(name = lab_x, breaks = scales::pretty_breaks(7)) +
      scale_y_continuous(name = lab_y, breaks = scales::pretty_breaks(7)) +
      scale_color_manual(values = MSDcol) +
      facet_wrap(~ID, scales = f_scales) +
      theme(legend.position = "none") +
      labs(subtitle = str_c( " (part ", n, " out of ", max(ds_tp_filt$NTILE), ")"),
           caption = cap)
    if(log_y){
      p_n <- p_n + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x)))
    }
    return(p_n)
  })
  return(p_i)
}

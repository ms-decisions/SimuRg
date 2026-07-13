## Author: Melnikova Alina,  Mikhailova Anna
## First created: 2025-09-04
## Description: sg-parsum and its helper functions
## Keywords: SimuRg, sg-parsum, goodness-of-fit

#' Extract parameter summary (theta, eta, SE, RSE, CI, CV, shrinkage)
#'
#' @inheritParams sg_dummy
#' @returns Table derived from [GFO] `$SUMTAB` (and `$OFV` when `addOFV = TRUE`).
#' @examples
#' fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
#' sum_tab <- sg_parsum(fpath_i)
#' print(sum_tab)
#' @import dplyr
#' @import stringr
#' @export
sg_parsum <- function(fpath_i, addOFV = TRUE){
  if (inherits(fpath_i, "character")) {
    obj <- get(load(fpath_i))
  } else if (inherits(fpath_i, "list")) {
    obj <- fpath_i
  } else {
    stop("fpath_i object should be either an sg_fit object, or a path to saved sg_fit object")
  }
  obj_partab_i <- obj$SUMTAB
  par_i <- obj_partab_i %>% select(Parameter = PAR, TYPE, Estimate = VALUE, SE, RSE, CV, LCI, UCI, Shrinkage_var = ETAshrinkage_var) %>%
    mutate(`95% CI` = paste0(signif(LCI,3)," - ", signif(UCI,3)), Estimate = signif(Estimate,3), SE = signif(SE,3), RSE = signif(RSE,3), `Shrinkage (var), %` = signif(Shrinkage_var,3))
  par_out <- par_i %>% mutate( `CV % (95% CI)` = ifelse(TYPE == "Random effects", str_c(round(CV, 3), " (", signif(100*sqrt(exp(LCI^2) - 1), 3), ", ", signif(100*sqrt(exp(UCI^2) - 1), 3), ")"), NA)) %>%
    #select(-UCI, -LCI, -CV, -TYPE)
    select(Parameter, Estimate, SE, RSE, `95% CI`, `CV % (95% CI)`, `Shrinkage (var), %`)

  par_out <- par_out %>%
    mutate(Parameter = str_remove_all(Parameter, "_pop")) %>% mutate(Parameter = str_replace_all(Parameter, "_", " "))

  if (addOFV) {
    ofv_lst <- obj$OFV

    par_out <- par_out %>% mutate(OFV = -2*ofv_lst[["LL"]], AIC = ofv_lst[["AIC"]])
  }

  return(par_out)
}

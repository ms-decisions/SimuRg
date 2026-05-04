## Author: Mikhailova Anna
## First created: 2025-10-27
## Description: functions for vpc calculations
## Keywords: SimuRg, vpc, diagnostics

#' Perform simulations for VPC plot
#'
#' @inheritParams sg_dummy
#' @param fpath_i Either a character string specifying the file path to a saved `sg_fit` object (R data file), or a list object containing the `sg_fit` results directly. The object must contain: `SDTAB`, `EVTAB`, `SUMTAB`, `OMEGAMAT`, `SIGMAMAT`, and optionally `COTAB` and `CATAB`
#' @param gco Generalized control object. Either an R list or path to Rdata/json file. Either `gco` or `model` file should be specified for simulations model specification. Defauls is `NULL`
#' @param npop Integer specifying the number of virtual subjects to simulate per original individual. Higher values provide more robust percentile estimates but increase computation time. Default is `100`
#' @returns A dataset with simulation results
#' @examples
#' \donttest{
#' library(rxode2)
#' library(tibble)
#' fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
#' # Simulate VPC from with generalized model object
#' mod <- rxode2::rxode2({
#'   ka_pop = 0.1;
#'   Vd_pop = 10;
#'   Cl_pop = 0.5;
#'
#'   omega_ka = 0;
#'   omega_Vd = 0;
#'   omega_Cl = 0;
#'
#'   Cc_b = 0;
#'   ka_tv = exp(ka_pop);
#'   Vd_tv = exp(Vd_pop);
#'   Cl_tv = exp(Cl_pop);
#'
#'   ka = ka_tv * exp(omega_ka);
#'   Vd = Vd_tv * exp(omega_Vd);
#'   Cl = Cl_tv * exp(omega_Cl);
#'
#'   Cc = Ac / Vd;
#'
#'   Ad(0) = 0;
#'   Ac(0) = 0;
#'
#'   d/dt(Ad) = -ka * Ad;
#'   d/dt(Ac) = ka * Ad - Cl * Cc;
#'
#'   Cc_ResErr = Cc * (1 + Cc_b);
#' })
#'
#' sg_vpc_sim(fpath_i, model=mod, outputs = "Cc_ResErr")
#' # Make VPC plot with the generalized control object
#' model <- system.file("extdata", "models", "rxode", "model_PK_1c.txt", package = "SimuRg")
#' data  <- system.file("extdata", "datasets", "dspk-warf.csv", package = "SimuRg")
#' headers <- list(list(name = "ID", use = "identifier", type = NULL),
#'              list(name = "TIME", use = "time", type = NULL),
#'                           list(name = "DV", use = "observation", type = "continuous"),
#'                           list(name = "DVID", use = "observationtype", type = NULL),
#'                           list(name = "ADM", use = "administration", type = NULL),
#'                           list(name = "AMT", use = "amount", type = NULL),
#'                           list(name = "EVID", use = "eventidentifier", type = NULL),
#'                           list(name = "MDV", use = "missingdependentvariable", type = NULL),
#'                           list(name = "AGE", use = "covariate", type = "continuous"),
#'                           list(name = "AGE_centered", use = "covariate", type = "continuous"),
#'                           list(name = "SEX", use = "covariate", type = "categorical"),
#'                           list(name = "WEIGHT", use = "covariate", type = "continuous"),
#'                           list(name = "BMI", use = "covariate", type = "continuous"))
#'
#' theta <- tribble(~NAME, ~TRANS, ~INIT, ~LB, ~UB, ~EST,
#'                   "Cl", "logNormal", 0.2, NA, NA, TRUE,
#'                 "V", "logNormal", 20, NA, NA, TRUE,
#'                 "ka", "logNormal", 0.2, NA, NA, TRUE
#')
#'
#' ruv <- list(YNAME = "y1", DVID = 1, TRANS = "normal", PRED = "Cc",
#'            ERR = "combined1", INIT = c(1, 1), EST = c(TRUE, TRUE), BLQM = NULL)
#'
#' re <- list(init = tribble(~Cl, ~V, ~ka,
#'                             1, 0, 0,
#'                             0, 1, 0,
#'                             0, 0, 1) %>% as.matrix(),
#'             est = tribble(~Cl, ~V, ~ka,
#'                          TRUE, NA, NA,
#'                            NA, TRUE, NA,
#'                           NA, NA, TRUE) %>% as.matrix())
#' occ <- list(init = tribble(~Cl, ~V, ~ka,
#'                         0, 0, 0,
#'                           0, 0, 0,
#'                           0, 0, 0) %>% as.matrix(),
#'            est = tribble(~Cl, ~V, ~ka,
#'                          NA, NA, NA,
#'                          NA, NA, NA,
#'                          NA, NA, NA) %>% as.matrix())
#' covs <- list(list(PAR = "V", COVNAME = "AGE", FUNC = "linear",
#'                  TRANS = "median", INIT = 1, EST = TRUE),
#'             list(PAR = "ka", COVNAME = "SEX", REF = 0, INIT = 1, EST = TRUE))
#'
#' gco <-list(headers = headers, data = data, model = model, task_opt = "",
#'           covs = covs, project_name = "test-proj", theta = theta,
#'           ruv = ruv, re = re, occ = occ, modelText = "")
#'
#' res <- sg_vpc_sim(fpath_i, gco = gco, output = "Cc")
#' }
#' @import rxode2
#' @importFrom purrr map_dfr
#' @import dplyr
#' @importFrom stringr str_remove
#' @importFrom rlang .data
#' @export
sg_vpc_sim <- function(fpath_i, gco = NULL, model = NULL, time_col = "TIME", outputs = NULL, npop = 100){
  obj <- read_smrg_obj(fpath_i)
  if (is_null(model) & is_null(gco)) {
    stop("Specify either a generalized control object (gco) or model to simulate from")
  } else if (is_null(model) & !is_null(gco)) {
    model <- rxode2::rxode2(gmo_converter(gco, output_path = NULL))
  } else if (!is_null(model) & !is_null(gco)) {
    message("Both gco and model specified. Model from gco is used for simulations")
    model <- rxode2::rxode2(gmo_converter(gco))
  }
  data_fin.noex <-  obj$SDTAB %>% filter(MDV != 1) %>% select(-MDV)
  data_fin.noex$time  <- data_fin.noex[[time_col]]
  ev_tab <- obj$EVTAB

  if (!is.null(obj$COTAB)) ev_tab <-  merge(ev_tab, obj$COTAB, by = "ID", all.x = T)
  if (!is.null(obj$CATAB)) ev_tab <-  merge(ev_tab, obj$CATAB, by = "ID", all.x = T)
  if (!is_empty(obj$REGTAB)) ev_tab <-  merge(ev_tab, obj$REGTAB, by = c("ID", time_col), all.x = T)
  covs_i <- c(colnames(obj$COTAB), colnames(obj$CATAB))
  covs_i <- covs_i[covs_i != "ID"]
  id_seq <- unique(data_fin.noex$ID)
  par_fin_tv <- obj$SUMTAB %>% filter(TYPE == "Typical values") %>%
    select(PAR, VALUE) %>% #mutate(PAR = str_remove(PAR, "_pop")) %>%
    deframe()
  sim_vpc_full <- id_seq %>% map_dfr(function(id_seq.i){
    data_fin.noex.i <- data_fin.noex %>% filter(ID == id_seq.i) %>% pull(time)
    ev_tab.i <- ev_tab %>% filter(ID == id_seq.i) %>%
      rename(DEFID = ID) %>% mutate(id = 1)
    sim.i <- sg_sim(model = model, et = ev_tab.i, stimes = data_fin.noex.i,
                    outputs = outputs, theta = par_fin_tv, omega = obj$OMEGAMAT,
                    sigma = obj$SIGMAMAT, covs = covs_i, nsub = npop, byID = TRUE,
                    addcov = F, ncores = parallel::detectCores()-1) %>%
      mutate(ID = id_seq.i)
    return(sim.i)
  })
  return(sim_vpc_full)
}

## Author: Mikhailova Anna
## First created: 2025-10-27
## Description: functions for vpc calculations
## Keywords: SimuRg, vpc, diagnostics

#' Perform simulations for VPC plot
#'
#' @inheritParams sg_dummy
#' @returns A dataset with simulation results
#' @examples
#' @import rxode2
#' @importFrom purrr map_dfr
#' @import dplyr
#' @export
sg_vpc_sim <- function(fpath_i, model, time_col = "TIME", output = NULL, npop = 100){
  if (inherits(fpath_i, "character")) {
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
  data_fin.noex <-  obj$SDTAB %>% filter(MDV != 1) %>% select(-MDV)
  data_fin.noex$time  <- data_fin.noex[[time_col]]
  ev_tab <- obj$EVTAB
  if (!is.null(obj$COTAB)) ev_tab <-  merge(ev_tab, obj$COTAB, by = "ID", all.x = T)
  if (!is.null(obj$CATAB)) ev_tab <-  merge(ev_tab, obj$CATAB, by = "ID", all.x = T)
  covs_i <- c(colnames(obj$COTAB), colnames(obj$CATAB))
  covs_i <- covs_i[covs_i != "ID"]
  id_seq <- unique(data_fin.noex$ID)
  par_fin_tv <- obj$SUMTAB %>% filter(TYPE == "Typical values") %>%
    select(PAR, VALUE) %>% mutate(PAR = str_remove(PAR, "_pop")) %>% deframe()
  sim_vpc_full <- id_seq %>% map_dfr(function(id_seq.i){
    data_fin.noex.i <- data_fin.noex %>% filter(ID == id_seq.i) %>% pull(time)
    ev_tab.i <- ev_tab %>% filter(ID == id_seq.i) %>%
      rename(DEFID = ID) %>% mutate(id = 1)
    sim.i <- sg_sim(model = model, et = ev_tab.i, stimes = data_fin.noex.i,
                    output = output, theta = par_fin_tv, omega = obj$OMEGAMAT,
                    sigma = obj$SIGMAMAT, covs = covs_i, npop = npop,
                    addcov = F, ncores = parallel::detectCores()-1) %>%
      mutate(ID = id_seq.i)
    return(sim.i)
  })
  return(sim_vpc_full)
}

library(tidyverse)
load("data-raw/example_code/simurg_object/Warfarin_PK.RData")
sg_vpc_sim <- function(fpath_i, mod_fin, time_col = "TIME", output = NULL, nrep = 100){
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
    # message(str_c(which(id_seq == id_seq.i), " out of ", length(id_seq)))
    data_fin.noex.i <- data_fin.noex %>% filter(ID == id_seq.i) %>% pull(time)
    ev_tab.i <- ev_tab %>% filter(ID == id_seq.i) %>%
      rename(DEFID = ID) %>% mutate(id = 1)
    sim.i <- sg_sim(model = mod_fin, et = ev_tab.i, stimes = data_fin.noex.i,
                    output = output, theta = par_fin_tv, omega = obj$OMEGAMAT,
                    sigma = obj$SIGMAMAT, covs = covs_i, npop = nrep,
                    addcov = F, ncores = parallel::detectCores()-1) %>%
      mutate(ID = unique(data_fin.noex.i$DEFID)) %>% select(-id)
    return(sim.i)
  })
  return(sim_vpc_full)
}

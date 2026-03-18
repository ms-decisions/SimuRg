library(tidyverse)
library(rxode2)
mod_fin <- RxODE({
  # Doses in mg
  # Time in hours

  ### Parameter values
  # Typical
  Cl_pop = 5;
  V_pop = 180;

  ka_pop = 6;


  # Random effects
  omega_Cl = 0;
  omega_V = 0;
  omega_ka = 0;

  # Residual error
  b = 0;

  ### Parameters
  Cl = Cl_pop * exp(omega_Cl);
  V = V_pop * exp(omega_V);
  ka = ka_pop * exp(omega_ka);

  ### Explicit functions
  Cc = Ac/V;                 # nmol/L

  ### Initial conditions
  Ad(0) = 0;          # mg
  Ac(0) = 0;          # mg

  ### ODEs
  d/dt(Ad) = - ka*Ad;
  d/dt(Ac) = ka*Ad - Cl*Cc ;

  CHECKRUV = b;
  Cc_ResErr = Cc + b*Cc;
})

load("inst/extdata/simurg_object/Warfarin_PK.RData")
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
    select(PAR, VALUE)  %>% deframe()#%>% mutate(PAR = str_remove(PAR, "_pop"))
  sim_vpc_full <- id_seq %>% map_dfr(function(id_seq.i){
    data_fin.noex.i <- data_fin.noex %>% filter(ID == id_seq.i) %>% pull(time)
    ev_tab.i <- ev_tab %>% filter(ID == id_seq.i) %>%
      rename(DEFID = ID) %>% mutate(id = 1)
    omegamat <- obj$OMEGAMAT
    rownames(omegamat) <- colnames(omegamat)
    sim.i <- sg_sim(model = mod_fin, et = ev_tab.i, stimes = data_fin.noex.i,
                    output = output, theta = par_fin_tv, omega = omegamat,
                    sigma = obj$SIGMAMAT, covs = covs_i, nsub = nrep,
                    addcov = F, ncores = parallel::detectCores()-1) %>%
      mutate(ID = id_seq.i)
    return(sim.i)
  })
  return(sim_vpc_full)
}
res <- sg_vpc_sim(obj1, mod_fin, output = "Cc_ResErr")


sg_vpc_sim <- function(fpath_i, mod_fin, time_col = "TIME", output = NULL, grps_vpc = NULL, tsld = F, pc = T,
                       ds_obs_i = data_fin.vpc, ds_pred_i = sim_vpc_full, ds_pred_tv = sim_vpc_tv,
                       n_bins = 15, lq_i = 0.1, uq_i = 0.9,
                       predp = F, obsdat = F, binlim = F, logy = F, logx = F, xlab = NULL, ylab = NULL,
                       sc_factor = 1, bincent = "BINM", fcsc = "free", xlim = NULL, ylim = NULL){
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
  data_fin.noex <-  obj$SDTAB %>% filter(MDV != 1) %>% rename(time = ATSFD) %>% select(-MDV)
  ev_tab <- obj$EVTAB
  id_seq <- unique(data_fin.noex$ID)
  sim_vpc_full <- id_seq %>% map_dfr(function(id_seq.i){
    message(str_c(which(id_seq == id_seq.i), " out of ", length(id_seq)))
    data_fin.noex.i <- data_fin.noex %>% filter(ID == id_seq.i) %>% pull(time)
    ev_tab.i <- ev_tab %>% fiter(ID == id_seq.i) %>%
      rename(DEFID = ID) %>% mutate(id = 1)
    sim.i <- sg_sim(model = mod_fin, et = ev_tab_i, stimes = data_fin.noex.i,
                    output = output, theta_i = par_fin_tv, omega_i = m_omega_full,
                    sigma_i = m_reserr, cov_i = covs_i, nrep = nrep_i,
                    addcov = F, ncores = parallel::detectCores()-1) %>%
      mutate(ID = unique(data_fin.noex.i$DEFID)) %>% select(-id)
    return(sim.i)
  })
  return(sim_vpc_full)
}

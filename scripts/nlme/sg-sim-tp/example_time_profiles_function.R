fun_TProf <- function(
    path_i, yname, ylabel, xlabel = "Time, hours",
    log_axes = F, sc_factor = 1, sc_times = 1, y_lim = NULL,
    id_char = ds_id_char, id_dose = NULL, nt = 30, f_nrow = 7, spec_IDs = NULL, tsld = F, ds_i = ds_poppkpd, route_i = NULL
){
  path_IndFits <- "ChartsData/IndividualFits/"
  flist_IndFits <- list.files(path = str_c(path_i, "/", path_IndFits)) %>% grep(yname, ., value = TRUE)
  print(flist_IndFits)
  
  ds_fits <- read_csv(str_c(path_i, "/", path_IndFits, "/", yname, "_fits.txt"), col_types = cols()) %>% 
    select(ID, TIME = time, Population = popPred, Individual = indivPredMode) %>% 
    gather("STATUS", "PRED", -ID, -TIME) %>% 
    mutate(STATUS = fct_relevel(STATUS, "Population", "Individual"),
           # IDmod = str_remove(ID, "SNDX-6352-"),
           PRED = PRED*sc_factor, TIME = TIME*sc_times)
  ds_obs <- read_csv(str_c(path_i, "/", path_IndFits, "/", yname, "_observations.txt"), col_types = cols()) %>% 
    select(ID, TIME = time, DV = !!yname, censored) %>% 
    mutate(# IDmod = str_remove(ID, "SNDX-6352-"), 
           DV = DV*sc_factor, TIME = TIME*sc_times)
  
  if(tsld & !is.null(ds_i)){
    # ds_obs <- ds_obs %>% mutate(TIMER = round(TIME, 1)) %>% 
    #   left_join(ds_i %>% 
    #               filter(DVIDN == parse_number(yname)) %>% 
    #               mutate(TIMER = round(AFRLT, 1), TSLD = APRLT*sc_times) %>% 
    #               select(ID = USUBJID, TIMER, TSLD) %>% unique(), by = c("ID", "TIMER")) %>% 
    #   select(-TIMER)
    doses_i <- ds_i %>% filter(EVID == 1) %>% select(ID = USUBJID, TIME = AFRLT) %>% mutate(TIME = TIME*sc_times)
    ds_fits <- bind_rows(mutate(ds_fits, type = "pk"), mutate(doses_i, type = "dose")) %>%
      group_by(ID) %>% arrange(TIME, .by_group = T) %>%
      mutate(
        LASTDT = if_else(type == "dose", TIME, NA_real_),
        OCC = cumsum(type == "dose")
      ) %>%
      fill(LASTDT, OCC, .direction = "down") %>%
      filter(type == "pk") %>%
      mutate(TSLD = TIME - LASTDT) %>%
      ungroup() %>% select(-type, -LASTDT)
    
    ds_obs <- bind_rows(mutate(ds_obs, type = "pk"), mutate(doses_i, type = "dose")) %>%
      group_by(ID) %>% arrange(TIME, .by_group = T) %>%
      mutate(
        LASTDT = if_else(type == "dose", TIME, NA_real_),
        OCC = cumsum(type == "dose")
      ) %>%
      fill(LASTDT, OCC, .direction = "down") %>%
      filter(type == "pk") %>%
      mutate(TSLD = TIME - LASTDT) %>%
      ungroup() %>% select(-type, -LASTDT)
  }
  
  if(!is.null(spec_IDs)){
    ds_fits <- ds_fits %>% filter(ID %in% spec_IDs)
    ds_obs <- ds_obs %>% filter(ID %in% spec_IDs)
  }
  if(!is.null(id_char)){
    ds_fits <- ds_fits %>% left_join(ds_id_char, by = "ID") %>% mutate(IDmod = str_c(str_remove(ID, "SNDX-6352-|INCA034176-|INCA34176-"), "\n", Dose, " mg/kg ", Route, " ", Regimen), IDmod = fct_inorder(IDmod))
    ds_obs <- ds_obs %>% left_join(ds_id_char, by = "ID") %>% mutate(IDmod = str_c(str_remove(ID, "SNDX-6352-|INCA034176-|INCA34176-"), "\n", Dose, " mg/kg ", Route, " ", Regimen), IDmod = fct_inorder(IDmod))
    if(!is.null(route_i)){
      ds_fits <- ds_fits %>% filter(Route == route_i)
      ds_obs <- ds_obs %>% filter(Route == route_i)
    }
  }
  
  id_tile <- ds_obs %>% select(ID) %>% unique() %>% mutate(SET = ntile(ID, ceiling(n_distinct(ds_obs$ID)/nt)))
  ds_fits <- ds_fits %>% left_join(id_tile, by = "ID")
  ds_obs <- ds_obs %>% left_join(id_tile, by = "ID")
  if(!is.null(id_dose)){
    id_dose <- id_dose %>% left_join(unique(select(ds_fits, ID, IDmod, SET)), by = "ID") %>% 
      mutate(TIME = TIME*sc_times)
  }
  
  p_TimeProf <- unique(ds_obs$SET) %>% map(function(p){
    ds_fits_p <- ds_fits %>% filter(SET == p)
    ds_obs_p <- ds_obs %>% filter(SET == p)
    
    if(tsld){
      ds_fits_p["TIME"] <- ds_fits_p["TSLD"]
      ds_obs_p["TIME"] <- ds_obs_p["TSLD"]
    } else {
      ds_fits_p <- ds_fits_p %>% mutate(OCC = 1)
      ds_obs_p <- ds_obs_p %>% mutate(OCC = 1)
    }
    
    p_p <- ggplot() +
      geom_line(data = ds_fits_p, aes(x = TIME, y = PRED, col = STATUS, group = interaction(ID, OCC, STATUS)), size = 0.6)
    if(!is.null(id_dose)){
      id_dose_p <- id_dose %>% filter(SET == p)
      p_p <- p_p +
        geom_vline(data = id_dose_p, aes(xintercept = TIME), col = "grey25", lty = "longdash", size = 0.4, alpha = 0.5)
    }
    if(nrow(filter(ds_obs_p, censored == 1)) > 0){
      if(log_axes){
        p_p <- p_p +
          geom_point(data = filter(ds_obs_p, censored == 1), aes(x = TIME, y = 0.15), col = "firebrick", shape = 15, size = 1.5)
      } else {
        p_p <- p_p +
          geom_point(data = filter(ds_obs_p, censored == 1), aes(x = TIME, y = -Inf), col = "firebrick", shape = 15, size = 1.5)
      }
      
    }
    p_p <- p_p +
      geom_point(data = filter(ds_obs_p, censored != 1), aes(x = TIME, y = DV), col = "black", 
                 shape = 21, fill = "grey85", size = 1.25, stroke = 0.8) +
      facet_wrap(~IDmod, scales = "free", nrow = f_nrow) +
      scale_color_manual(name = "Prediction", values = MSDcol[c(1, 4)]) +
      scale_x_continuous(breaks = scales::pretty_breaks(7)) +
      coord_cartesian(ylim = y_lim) +
      labs(x = xlabel, y = ylabel)
    if(log_axes){
      p_p <- p_p +
        scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                      labels = scales::trans_format("log10", scales::math_format(10^.x)))
    } else {
      p_p <- p_p +
        scale_y_continuous(breaks = scales::pretty_breaks(7))
    }
    return(p_p)
  })
  return(p_TimeProf)
}
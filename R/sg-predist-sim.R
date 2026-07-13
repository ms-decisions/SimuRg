#' Perform simulations for prediction distribution plots
#'
#' @inheritParams sg_dummy
#' @param gco Generalized control object (see [GMO]). Either an R list or path to
#'   Rdata/json file. Either `gco` or `model` should be specified for the
#'   simulation model specification. Default is `NULL`
#' @returns [GSO]: a data frame with simulation results.
#' @seealso [GFO], [GMO], [GSO], [sg_sim()], [sg_vpc_sim()]
#' @examples
#' \donttest{
#' library(rxode2)
#' fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
#' mod_fin <- rxode2({
#'   # Differential equations
#'   d/dt(Ad) = -ka * Ad
#'   d/dt(Ac) = ka * Ad - Cl/V * Ac
#'   # Concentration calculations
#'   Cc = Ac / V
#' })
#' sg_predist_sim(fpath_i, model = mod_fin, outputs = "Cc")
#' }
#' @import rxode2
#' @importFrom purrr map_dfr
#' @import dplyr
#' @importFrom stringr str_remove
#' @importFrom rlang .data is_null is_empty
#' @export
sg_predist_sim <- function(fpath_i, gco = NULL, model = NULL, time_col = "TIME", outputs = NULL, npop = 500){
  obj <- read_smrg_obj(fpath_i)
  if (is_null(model) & is_null(gco)) {
    stop("Specify either a generalized control object (gco) or model to simulate from")
  } else if (is_null(model) & !is_null(gco)) {
    model <- rxode2::rxode2(gmo_converter(gco, output_path = NULL))
  } else if (!is_null(model) & !is_null(gco)) {
    message("Both gco and model specified. Model from gco is used for simulations")
    model <- rxode2::rxode2(gmo_converter(gco))
  }

  data_fin.noex <- obj$SDTAB %>% filter(MDV != 1) %>% select(-MDV)
  data_fin.noex$time <- data_fin.noex[[time_col]]
  ev_tab <- obj$EVTAB

  if (!is.null(obj$COTAB)) ev_tab <- merge(ev_tab, obj$COTAB, by = "ID", all.x = TRUE)
  if (!is.null(obj$CATAB)) ev_tab <- merge(ev_tab, obj$CATAB, by = "ID", all.x = TRUE)
  if (!is_empty(obj$REGTAB)) ev_tab <- merge(ev_tab, obj$REGTAB, by = c("ID", time_col), all.x = TRUE)

  covs_i <- c(colnames(obj$COTAB), colnames(obj$CATAB))
  covs_i <- covs_i[covs_i != "ID"]
  id_seq <- unique(data_fin.noex$ID)

  par_fin_tv <- obj$SUMTAB %>%
    filter(TYPE == "Typical values") %>%
    select(PAR, VALUE) %>%
    #mutate(PAR = str_remove(PAR, "_pop")) %>%
    deframe()

  sim_predist_full <- id_seq %>% map_dfr(function(id_seq.i){
    data_fin.noex.i <- data_fin.noex %>%
      filter(ID == id_seq.i) %>%
      pull(time)

    ev_tab.i <- ev_tab %>%
      filter(ID == id_seq.i) %>%
      rename(DEFID = ID) %>%
      mutate(id = 1)

    sim.i <- sg_sim(model = model,
                    et = ev_tab.i,
                    stimes = data_fin.noex.i,
                    outputs = outputs,
                    theta = par_fin_tv,
                    omega = obj$OMEGAMAT,
                    sigma = NULL,
                    covs = covs_i,
                    nsub = npop,
                    byID = TRUE,
                    addcov = FALSE,
                    ncores = parallel::detectCores()-1) %>%
      mutate(ID = id_seq.i)
    return(sim.i)
  })

  return(sim_predist_full)
}

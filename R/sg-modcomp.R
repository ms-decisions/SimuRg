## Author: Melnikova Alina,  Mikhailova Anna
## First created: 2025-09-05
## Description: sg-modcomp and its helper functions
## Keywords: SimuRg, sg-modcomp, model building, comparison

#' Returns dataframe with summary statistics for model comparison
#'
#' @inheritParams sg_dummy
#' @returns Data frame with model comparison metrics derived from [GFO] `$OFV` and
#'   related fit outputs.
#' @examples
#' fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
#' sum_tab <- sg_modcomp(fpath_i, 3)
#' print(sum_tab)
#' @import dplyr
#' @import stringr
#' @import stringr
#' @export
sg_modcomp <- function(fpath_i, run_id = 1){

  # Load the object from the specified path
  if (inherits(fpath_i, "character")) {
    obj <- get(load(fpath_i))
  } else if (inherits(fpath_i, "list")) {
    obj <- fpath_i
  } else {
    stop("fpath_i object should be either an sg_fit object, or a path to saved sg_fit object")
  }


  # Initialize empty dataframe for model comparisons
  modcomp <- data.frame()

  # Extract Project name
  project_name <- obj$PROJNAME

  # Extract parameters from parsum
  parameters <- obj$SUMTAB %>% filter(TYPE=="Typical values") %>% select(PAR) %>% pull()
  parameters <- parameters %>% str_remove("_pop")

  # Extract ETAs - get column names starting with "eta_" from PATAB
  eta_cols <- grep("^eta_", names(obj$PATAB), value = TRUE)
  eta_names <- gsub("^eta_", "", eta_cols)  # Remove "eta_" prefix

  # Extract Covariates - process SUMTAB$PAR for beta_ terms

  covar_vec <- c()
  cont_cov <- names(obj$CATAB)
  cat_cov <- names(obj$COTAB)

  if ("SUMTAB" %in% names(obj)) {
    beta_terms <- grep("^beta_", obj$SUMTAB$PAR, value = TRUE)

    cov_pairs <- gsub("^beta_", "", beta_terms)
    for (term in cov_pairs) {
      # Split beta_[covariate1]_[covariate2] into components
      cont_cov <- names(obj$CATAB)
      cat_cov <- names(obj$COTAB)

      for (par in parameters){
        if (str_detect(term,fixed(par))){
          term1 <- str_remove(term, fixed(par))
          par_i <- par
          break} else {term1 <- term}
      }
      for (cont in cont_cov){
        if (str_detect(term1,fixed(cont))){
          term2 <- str_remove(term1, fixed(cont))
          cov_i <- cont
          break} else {term2 <- term1}
      }

      for (cat in cat_cov){
        if (str_detect(term2,fixed(cat))){
          term3 <- str_remove(term2, fixed(cat))
          cov_i <- cat
          break}
        else {term3 <- term2}
      }

      covar_vec <- c(covar_vec, str_c( cov_i,par_i, sep = "~"))
    }
  }

  # Count identical elements and rename if n>1
  if (length(covar_vec) > 0) {
    # Count occurrences of each element
    covar_counts <- table(covar_vec)

    # Create new vector with renamed elements
    covar_vec_renamed <- covar_vec

    for (i in seq_along(covar_vec)) {
      element <- covar_vec[i]
      count <- covar_counts[element]

      if (count > 1) {
        # Extract cov_i and par_i from the element (format: "cov_i~par_i")
        parts <- str_split(element, "~", simplify = TRUE)
        if (length(parts) == 2) {
          cov_i <- parts[1]
          par_i <- parts[2]
          # Rename with count: "cov_i(n levels)~par_i"
          covar_vec_renamed[i] <- str_c(cov_i, "(", count, " levels)", "~", par_i)
        }
      }
    }

    # Replace original vector with renamed version
    covar_vec <- covar_vec_renamed
  }

  # Remove identical elements (duplicates) from covar_vec
  if (length(covar_vec) > 0) {
    covar_vec <- unique(covar_vec)
  }

  # Create dataframe with covariate information
  modcomp <- data.frame(
    Run_ID = run_id,
    `Project_name` = project_name,
    `ETAs` = paste(eta_names, collapse = ", "),
    `Covariates` = paste(unique(covar_vec), collapse = ", "),
    stringsAsFactors = FALSE
  )

  # Create dataframe with parameter Identifiability inforamtion
  ds_out_i <- sg_parsum(fpath_i, addOFV = TRUE)

  ds_idetn_i <- ds_out_i %>% filter(RSE > 50)
  ident_stat <- NULL
  if(nrow(ds_idetn_i) > 0){
    ident_stat <- str_c(str_remove(ds_idetn_i$Parameter, "[,].*?$"), collapse = "; ")
  }

  ds_shr_i <- ds_out_i %>% filter(`Shrinkage (var), %` > 50)
  shr_stat <- NULL
  if(nrow(ds_shr_i) > 0){
    shr_stat <- str_c(str_remove(ds_shr_i$Parameter, "[,].*?$"), collapse = "; ")
  }


  ds_out <- modcomp %>% mutate(Identifiability = ifelse(nrow(ds_idetn_i) == 0, "Yes", "No"),
                               Shrinkage = ifelse(nrow(ds_shr_i) == 0, "Yes", "No"),
                               `Pars. with RSE > 50%` = ident_stat, `Pars. with shr. > 50%` = shr_stat,
                               `Max RSE, %` = max(ds_out_i$RSE, na.rm = TRUE),
                               `Max shrinkage, %` = max(ds_out_i$`Shrinkage (var), %`, na.rm = TRUE),
                               OFV = unique(ds_out_i$OFV), AIC = unique(ds_out_i$AIC)) %>%
    rename(`Run ID` = Run_ID, `Project name` = Project_name)
  return(ds_out)
}


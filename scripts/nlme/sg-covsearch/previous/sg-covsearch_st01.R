## Author: Alina Melnikova
## First created: 2026-05-19
## Description: utility helpers for staged covariate search
## Keywords: SimuRg, covsearch

#' Calculate objective function value from fit output
#'
#' Converts log-likelihood (`LL`) from `gfo$OFV` into objective function value:
#' `OFV = -2 * LL`.
#'
#' @param gfo A generalized fit output object containing `OFV`.
#'
#' @return Numeric scalar OFV.
#' @keywords internal
get_ofv <- function(gfo) {
  if (is.null(gfo$OFV)) {
    stop("get_ofv: gfo$OFV is missing.")
  }

  ofv_tab <- as.data.frame(gfo$OFV, stringsAsFactors = FALSE)
  if (!"LL" %in% names(ofv_tab)) {
    stop("get_ofv: gfo$OFV must contain an LL column.")
  }
  if (nrow(ofv_tab) < 1 || is.na(ofv_tab$LL[[1]])) {
    stop("get_ofv: gfo$OFV$LL is missing or NA.")
  }

  -2 * as.numeric(ofv_tab$LL[[1]])
}


#' Build theta tibble updated with fitted typical values
#'
#' Preserves original parameter row order from `gco$theta`. When present,
#' values from `gfo$SUMTAB` are mapped by `PAR = paste0(NAME, "_pop")` and
#' overwrite `INIT`.
#'
#' @param gco A generalized control object containing `theta`.
#' @param gfo A generalized fit output object containing `SUMTAB`.
#'
#' @return Tibble with theta columns and updated `INIT`.
#' @keywords internal
gco_to_theta_tibble <- function(gco, gfo) {
  if (is.null(gco$theta)) {
    stop("gco_to_theta_tibble: gco$theta is missing.")
  }

  theta_tb <- tibble::as_tibble(gco$theta)
  required_cols <- c("NAME", "INIT")
  missing_cols <- setdiff(required_cols, names(theta_tb))
  if (length(missing_cols) > 0) {
    stop(
      sprintf(
        "gco_to_theta_tibble: gco$theta missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      )
    )
  }

  if (is.null(gfo$SUMTAB)) {
    return(theta_tb)
  }

  sumtab <- as.data.frame(gfo$SUMTAB, stringsAsFactors = FALSE)
  if (!all(c("PAR", "VALUE") %in% names(sumtab))) {
    return(theta_tb)
  }

  map_par <- paste0(theta_tb$NAME, "_pop")
  idx <- match(map_par, sumtab$PAR)
  matched <- !is.na(idx)
  if (any(matched)) {
    theta_tb$INIT[matched] <- sumtab$VALUE[idx[matched]]
  }

  theta_tb
}


#' Append one covariate relationship definition
#'
#' @param covs_list Existing list of covariate definitions, or `NULL`.
#' @param param Parameter name to modify.
#' @param cov Covariate name.
#' @param type Covariate type: `"continuous"` or `"categorical"`.
#' @param cov_ref Reference category value for categorical covariates.
#'
#' @return Updated list with appended covariate record.
#' @keywords internal
add_covariate <- function(covs_list, param, cov, type, cov_ref = NULL) {
  if (is.null(covs_list)) {
    covs_list <- list()
  }
  if (!is.list(covs_list)) {
    stop("add_covariate: covs_list must be a list or NULL.")
  }
  if (!is.character(param) || length(param) != 1 || !nzchar(param)) {
    stop("add_covariate: param must be a non-empty string.")
  }
  if (!is.character(cov) || length(cov) != 1 || !nzchar(cov)) {
    stop("add_covariate: cov must be a non-empty string.")
  }
  if (!is.character(type) || length(type) != 1 || !nzchar(type)) {
    stop("add_covariate: type must be a non-empty string.")
  }

  type_norm <- tolower(type)
  if (!type_norm %in% c("continuous", "categorical")) {
    stop("add_covariate: type must be either 'continuous' or 'categorical'.")
  }

  if (identical(type_norm, "continuous")) {
    new_cov <- list(
      PAR = param,
      COVNAME = cov,
      FUNC = "linear",
      TRANS = "median",
      INIT = 0.01,
      EST = TRUE
    )
  } else {
    if (is.null(cov_ref) || (length(cov_ref) == 1 && is.na(cov_ref))) {
      stop("add_covariate: cov_ref is required for categorical covariates.")
    }
    new_cov <- list(
      PAR = param,
      COVNAME = cov,
      REF = cov_ref,
      INIT = 0.01,
      EST = TRUE
    )
  }

  c(covs_list, list(new_cov))
}


#' Remove an exact parameter-covariate pair from covariate list
#'
#' @param covs_list Existing list of covariate definitions, or `NULL`.
#' @param param Parameter name.
#' @param cov Covariate name.
#'
#' @return Updated list with matching `(PAR, COVNAME)` entry removed.
#' @keywords internal
remove_covariate <- function(covs_list, param, cov) {
  if (is.null(covs_list)) {
    return(list())
  }
  if (!is.list(covs_list)) {
    stop("remove_covariate: covs_list must be a list or NULL.")
  }

  keep_idx <- vapply(
    covs_list,
    FUN.VALUE = logical(1),
    FUN = function(x) {
      has_pair <- is.list(x) && !is.null(x$PAR) && !is.null(x$COVNAME)
      if (!has_pair) {
        return(TRUE)
      }
      !(identical(x$PAR, param) && identical(x$COVNAME, cov))
    }
  )

  covs_list[keep_idx]
}

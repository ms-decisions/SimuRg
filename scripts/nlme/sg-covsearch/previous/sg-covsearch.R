## Author: Alina Melnikova
## First created: 2026-05-19
## Description: utility helpers for staged covariate search
## Keywords: SimuRg, covsearch

source("scripts/nlme/sg-covsearch/test-fixtures-sg-covsearch.R")

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


# .as_covsearch_df <- function(x, name) {
#   if (is.null(x)) {
#     stop(sprintf("stepwise_covariate_selection: %s is missing.", name))
#   }
#   as.data.frame(x, stringsAsFactors = FALSE)
# }

.as_covsearch_df <- function(x, name) {
  if (is.null(x)) {
    stop(sprintf("stepwise_covariate_selection: %s is missing.", name))
  }
  if (is.data.frame(x)) {
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }

  is_list_of_lists <- is.list(x) && length(x) > 0 &&
    all(vapply(x, is.list, logical(1)))
  if (is_list_of_lists) {
    col_names <- unique(unlist(lapply(x, names), use.names = FALSE))
    rows <- lapply(x, function(rec) {
      row <- setNames(vector("list", length(col_names)), col_names)
      for (nm in col_names) {
        val <- rec[[nm]]
        row[[nm]] <- if (is.null(val) || length(val) == 0) NA else val[[1]]
      }
      as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
    })
    return(do.call(rbind, rows))
  }

  as.data.frame(x, stringsAsFactors = FALSE)
}


.norm_cov_type <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- ifelse(x %in% c("cont", "continuous"), "cont",
                ifelse(x %in% c("cat", "categorical"), "cat", NA_character_))
  out
}


.mode_sorted_smallest <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  freq <- table(x)
  top <- names(freq)[freq == max(freq)]
  sort(top)[1]
}


.na_safe_median <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) {
    return(NA_real_)
  }
  stats::median(x, na.rm = TRUE)
}


#' Prepare validated candidate pairs and references for stepwise covsearch
#'
#' Stage 2 implementation:
#' - validates user inputs (`covariates`, `parameters`, `test_pairs`, p-values)
#' - derives covariate references from `COTAB` (continuous median, NA-safe) and
#'   `CATAB` (categorical mode)
#' - applies optional `test_pairs` filtering with warn+drop for invalid rows
#' - computes candidate-specific degrees of freedom
#'
#' @param gco Generalized control object containing at least `headers` and `theta`.
#' @param gfo Generalized fit output object containing at least `COTAB` and `CATAB`.
#' @param covariates Optional character vector of covariate names to consider.
#' @param parameters Optional character vector of parameter names to consider.
#' @param test_pairs Optional data.frame with columns:
#'   `parameter`, `covariate`, `type`, `reference`, `center`.
#' @param p_forward Numeric in (0,1), default 0.05.
#' @param p_backward Numeric in (0,1), default 0.01.
#'
#' @return A list with validated candidate table and preparation metadata.
#' @keywords internal
stepwise_covariate_selection <- function(gco, gfo,
                                         covariates = NULL,
                                         parameters = NULL,
                                         test_pairs = NULL,
                                         p_forward = 0.05,
                                         p_backward = 0.01) {
  if (!is.numeric(p_forward) || length(p_forward) != 1 || is.na(p_forward) ||
      p_forward <= 0 || p_forward >= 1) {
    stop("stepwise_covariate_selection: p_forward must be numeric in (0,1).")
  }
  if (!is.numeric(p_backward) || length(p_backward) != 1 || is.na(p_backward) ||
      p_backward <= 0 || p_backward >= 1) {
    stop("stepwise_covariate_selection: p_backward must be numeric in (0,1).")
  }

  #For test
  #gco <- list(headers = .covsearch_headers_fixture, theta = .covsearch_theta_fixture)
  #gfo <- list(COTAB = .covsearch_cotab_fixture, CATAB = .covsearch_catab_fixture)

  ###


  headers_df <- .as_covsearch_df(gco$headers, "gco$headers")
  theta_df <- .as_covsearch_df(gco$theta, "gco$theta")
  cotab_df <- .as_covsearch_df(gfo$COTAB, "gfo$COTAB")
  catab_df <- .as_covsearch_df(gfo$CATAB, "gfo$CATAB")

  if (!all(c("name", "use", "type") %in% names(headers_df))) {
    stop("stepwise_covariate_selection: gco$headers must contain name/use/type.")
  }
  if (!"NAME" %in% names(theta_df)) {
    stop("stepwise_covariate_selection: gco$theta must contain NAME.")
  }

  cov_headers <- headers_df[headers_df$use == "covariate", c("name", "type"), drop = FALSE]
  if (nrow(cov_headers) == 0) {
    stop("stepwise_covariate_selection: no covariates found in gco$headers.")
  }
  cov_headers$type <- .norm_cov_type(cov_headers$type)
  if (any(is.na(cov_headers$type))) {
    bad <- unique(cov_headers$name[is.na(cov_headers$type)])
    stop(
      sprintf(
        "stepwise_covariate_selection: unsupported covariate type for: %s",
        paste(bad, collapse = ", ")
      )
    )
  }

  valid_covariates <- as.character(cov_headers$name)
  valid_parameters <- as.character(theta_df$NAME)

  if (is.null(covariates)) {
    covariates <- valid_covariates
  } else {
    if (!is.character(covariates)) {
      stop("stepwise_covariate_selection: covariates must be a character vector.")
    }
    invalid_cov <- setdiff(unique(covariates), valid_covariates)
    if (length(invalid_cov) > 0) {
      stop(sprintf(
        "stepwise_covariate_selection: unknown covariates: %s",
        paste(invalid_cov, collapse = ", ")
      ))
    }
    covariates <- unique(covariates)
  }

  if (is.null(parameters)) {
    parameters <- valid_parameters
  } else {
    if (!is.character(parameters)) {
      stop("stepwise_covariate_selection: parameters must be a character vector.")
    }
    invalid_par <- setdiff(unique(parameters), valid_parameters)
    if (length(invalid_par) > 0) {
      stop(sprintf(
        "stepwise_covariate_selection: unknown parameters: %s",
        paste(invalid_par, collapse = ", ")
      ))
    }
    parameters <- unique(parameters)
  }

  cov_type_map <- setNames(cov_headers$type, cov_headers$name)
  cont_covs <- names(cov_type_map)[cov_type_map == "cont"]
  cat_covs <- names(cov_type_map)[cov_type_map == "cat"]

  cov_ref_cont <- setNames(vector("list", length(cont_covs)), cont_covs)
  for (cv in cont_covs) {
    cov_ref_cont[[cv]] <- .na_safe_median(cotab_df[[cv]])
  }

  cov_ref_cat <- setNames(vector("list", length(cat_covs)), cat_covs)
  for (cv in cat_covs) {
    cov_ref_cat[[cv]] <- .mode_sorted_smallest(catab_df[[cv]])
  }

  df_map <- setNames(numeric(length(cov_type_map)), names(cov_type_map))
  for (cv in names(cov_type_map)) {
    if (identical(cov_type_map[[cv]], "cont")) {
      df_map[[cv]] <- 1
    } else {
      non_na_vals <- unique(as.character(catab_df[[cv]][!is.na(catab_df[[cv]])]))
      df_map[[cv]] <- max(length(non_na_vals) - 1, 0)
    }
  }

  if (is.null(test_pairs)) {
    candidates <- expand.grid(
      parameter = parameters,
      covariate = covariates,
      stringsAsFactors = FALSE
    )
    candidates$type <- cov_type_map[candidates$covariate]
    candidates$reference <- NA_character_
    candidates$center <- NA_character_
  } else {
    if (!is.data.frame(test_pairs)) {
      stop("stepwise_covariate_selection: test_pairs must be a data.frame or NULL.")
    }
    needed_cols <- c("parameter", "covariate", "type")
    missing_cols <- setdiff(needed_cols, names(test_pairs))
    if (length(missing_cols) > 0) {
      stop(sprintf(
        "stepwise_covariate_selection: test_pairs missing required columns: %s",
        paste(missing_cols, collapse = ", ")
      ))
    }

    candidates <- as.data.frame(test_pairs, stringsAsFactors = FALSE)
    if (!"reference" %in% names(candidates)) candidates$reference <- NA_character_
    if (!"center" %in% names(candidates)) candidates$center <- NA_character_

    candidates$type <- .norm_cov_type(candidates$type)
    header_type <- cov_type_map[candidates$covariate]
    has_type <- !is.na(candidates$type)
    has_known_cov <- candidates$covariate %in% covariates
    has_known_par <- candidates$parameter %in% parameters
    type_match <- has_known_cov & has_type & (candidates$type == header_type)

    valid_row <- has_known_par & has_known_cov & has_type & type_match
    if (any(!valid_row)) {
      bad_n <- sum(!valid_row)
      warning(sprintf(
        "stepwise_covariate_selection: dropped %d invalid test_pairs row(s).",
        bad_n
      ))
    }
    candidates <- candidates[valid_row, c("parameter", "covariate", "type", "reference", "center"), drop = FALSE]
  }

  if (nrow(candidates) == 0) {
    stop("stepwise_covariate_selection: no valid candidate pairs remain.")
  }

  candidates$type <- .norm_cov_type(candidates$type)
  candidates$center <- as.character(candidates$center)
  candidates$reference <- as.character(candidates$reference)
  need_median_center <- candidates$type == "cont" &
    (is.na(candidates$center) | !nzchar(candidates$center))
  candidates$center[need_median_center] <- "median"

  candidates$cov_ref <- vapply(
    seq_len(nrow(candidates)),
    FUN.VALUE = character(1),
    FUN = function(i) {
      cv <- candidates$covariate[i]
      if (identical(candidates$type[i], "cont")) {
        as.character(cov_ref_cont[[cv]])
      } else {
        user_ref <- candidates$reference[i]
        if (!is.na(user_ref) && nzchar(user_ref)) {
          user_ref
        } else {
          as.character(cov_ref_cat[[cv]])
        }
      }
    }
  )
  candidates$df <- unname(df_map[candidates$covariate])

  cov_ref <- as.list(c(cov_ref_cont, cov_ref_cat))
  for (i in seq_len(nrow(candidates))) {
    if (identical(candidates$type[i], "cat")) {
      user_ref <- candidates$reference[i]
      if (!is.na(user_ref) && nzchar(user_ref)) {
        cov_ref[[candidates$covariate[i]]] <- user_ref
      }
    }
  }

  list(
    candidates = candidates,
    cov_ref = cov_ref,
    df_map = as.list(df_map),
    parameters = parameters,
    covariates = covariates,
    p_forward = p_forward,
    p_backward = p_backward
  )
}

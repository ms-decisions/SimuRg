## Author: Alina Melnikova
## First created: 2026-05-19
## Description: covariate search tool
## Keywords: SimuRg, covsearch

#' Calculate objective function value from fit output
#'
#' Use OFV = -2*LL from gfo$OFV
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

  as.numeric(ofv_tab$LL[[1]])
}


#' Normalize SUMTAB to a proper data frame
#'
#' Attempts to coerce `gfo$SUMTAB` into a tabular structure even when JSON was
#' loaded as nested row-lists. When available, this uses `read_smrg_obj()` first
#' to reuse common SimuRg table-normalization logic.
#'
#' @param gfo A generalized fit output object.
#'
#' @return Data frame with at least standard SUMTAB columns when available.
#' @keywords internal
.covsearch_sumtab_df <- function(gfo) {
  if (is.null(gfo) || is.null(gfo$SUMTAB)) {
    return(NULL)
  }

  # Reuse shared reader/normalizer when available in the current session.
  if (exists("read_smrg_obj", mode = "function")) {
    gfo <- tryCatch(read_smrg_obj(gfo), error = function(e) gfo)
  }

  sumtab_raw <- gfo$SUMTAB
  sumtab <- NULL
  sumtab_01 <- NULL

  if (is.data.frame(sumtab_raw)) {
    sumtab <- sumtab_raw
  } else if (is.matrix(sumtab_raw)) {
    sumtab <- as.data.frame(sumtab_raw, stringsAsFactors = FALSE)
  } else if (is.list(sumtab_raw) && length(sumtab_raw) > 0) {
    # Typical case for JSON read with simplifyDataFrame = FALSE:
    # list of row-records; convert each row then bind.
    row_dfs <- lapply(sumtab_raw, function(row) {
      df <- as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
      for (nm in c("PAR", "VALUE")) {
        if (!nm %in% names(df)) {
          df[[nm]] <- NA
        }
      }
      df[, c("PAR", "VALUE"), drop = FALSE]
    })
    for (i in seq_along(row_dfs)) {
      if (i==1) {
        sumtab_01 <- row_dfs[[i]]
      } else {
        sumtab_01 <- rbind(sumtab_01, row_dfs[[i]])
      }
    }
    sumtab <- sumtab_01
  } else {
    sumtab <- tryCatch(
      as.data.frame(sumtab_raw, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  }

  if (is.null(sumtab)) {
    return(NULL)
  }

  # Ensure standard Monolix SUMTAB columns exist.
  # standard_cols <- c("PAR", "VALUE", "TYPE", "DISTRIBUTION", "EST", "SE", "RSE", "LCI", "UCI")
  # for (nm in setdiff(standard_cols, names(sumtab))) {
  #   sumtab[[nm]] <- NA
  # }
  # sumtab <- sumtab[, c(standard_cols, setdiff(names(sumtab), standard_cols)), drop = FALSE]

  # Keep these two parse-friendly for matching and numeric INIT overwrite.
  sumtab$PAR <- as.character(sumtab$PAR)
  sumtab$VALUE <- suppressWarnings(as.numeric(sumtab$VALUE))

  sumtab
}


#' Build theta tibble updated with fitted typical values
#'
#' Preserves original parameter row order from `gco$theta`. When present,
#' values from `gfo$SUMTAB` are mapped by `PAR = paste0(NAME, "_pop")` and
#' overwrite `INIT`.
#'
#' @param gco A generalized control object containing `theta`.
#' @param gfo A generalized fit output object containing `SUMTAB`.
#' @param update_theta_init Logical; when `TRUE`, update INIT from `gfo$SUMTAB`.
#'
#' @return Tibble with theta columns and updated `INIT`.
#' @keywords internal
gco_to_theta_tibble <- function(gco, gfo, update_theta_init = FALSE) {
  if (is.null(gco$theta)) {
    stop("gco_to_theta_tibble: gco$theta is missing.")
  }

  theta_raw <- gco$theta
  if (is.data.frame(theta_raw) || inherits(theta_raw, "tbl_df")) {
    # Preserve GCO values exactly (including INIT) when already tabular.
    theta_tb <- tibble::as_tibble(theta_raw)
  } else if (is.list(theta_raw) && length(theta_raw) > 0) {
    # Support list-of-records shape from some JSON loaders.
    theta_tb <- do.call(
      rbind,
      lapply(theta_raw, function(row) {
        as.data.frame(row, stringsAsFactors = FALSE, check.names = FALSE)
      })
    )
    theta_tb <- tibble::as_tibble(theta_tb)
  } else {
    stop("gco_to_theta_tibble: gco$theta must be a data.frame/tibble or non-empty list.")
  }

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

  # if (!isTRUE(update_theta_init)) {
  #   return(theta_tb)
  # } 

  # sumtab <- .covsearch_sumtab_df(gfo)
  # if (is.null(sumtab)) {
  #   return(theta_tb)
  # }
  # if (!all(c("PAR", "VALUE") %in% names(sumtab))) {
  #   return(theta_tb)
  # }

  # map_par <- paste0(theta_tb$NAME, "_pop")
  # idx <- match(map_par, sumtab$PAR)
  # matched <- !is.na(idx)
  # if (any(matched)) {
  #   theta_tb$INIT[matched] <- sumtab$VALUE[idx[matched]]
  # }

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
      INIT = 0,
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


.covsearch_null_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}


.covsearch_sanitize_name <- function(x) {
  out <- gsub("[^[:alnum:]_]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  if (!nzchar(out)) "x" else out
}


.covsearch_existing_covs <- function(x) {
  if (is.null(x)) {
    return(list())
  }
  if (!is.list(x) || length(x) == 0) {
    return(list())
  }
  all_records <- all(vapply(x, is.list, logical(1)))
  if (!all_records) {
    return(list())
  }
  x
}


.covsearch_json_to_object <- function(x, arg_name) {
  if (is.character(x) && length(x) == 1 && nzchar(x) && file.exists(x)) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop(sprintf(
        "stepwise_covariate_selection: package 'jsonlite' is required to read %s JSON path.",
        arg_name
      ))
    }
    parsed <- jsonlite::fromJSON(x, simplifyVector = FALSE)
    attr(parsed, "source_json_path") <- normalizePath(x, winslash = "/", mustWork = TRUE)
    return(parsed)
  }
  x
}


.covsearch_unwrap_scalar <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.atomic(x) && length(x) == 1) {
    return(x)
  }
  if (is.list(x) && length(x) == 0) {
    return(NULL)
  }
  if (is.list(x) && length(x) == 1 && !is.null(names(x)) && nzchar(names(x)[1])) {
    return(x[[1]])
  }
  x
}


.covsearch_to_header_records <- function(headers) {
  headers_df <- .as_covsearch_df(headers, "gco$headers")
  req <- c("name", "use")
  miss <- setdiff(req, names(headers_df))
  if (length(miss) > 0) {
    stop(sprintf(
      "stepwise_covariate_selection: gco$headers missing required columns: %s",
      paste(miss, collapse = ", ")
    ))
  }
  if (!"type" %in% names(headers_df)) {
    headers_df$type <- NA_character_
  }

  out <- lapply(seq_len(nrow(headers_df)), function(i) {
    type_i <- .covsearch_unwrap_scalar(headers_df$type[[i]])
    if (length(type_i) != 1 || is.na(type_i) || identical(type_i, "") || identical(type_i, "{}")) {
      type_i <- NULL
    } else {
      type_i <- as.character(type_i)
    }
    list(
      name = as.character(headers_df$name[[i]]),
      use = as.character(headers_df$use[[i]]),
      type = type_i
    )
  })
  out
}


.covsearch_to_matrix <- function(x, storage_mode = c("double", "logical"), field_name = "matrix") {
  storage_mode <- match.arg(storage_mode)
  if (is.null(x)) {
    return(NULL)
  }
  mat <- NULL
  if (is.matrix(x)) {
    mat <- x
  } else if (is.data.frame(x)) {
    mat <- as.matrix(x)
  } else if (is.list(x) && length(x) > 0) {
    rows <- lapply(x, function(row) {
      if (is.null(row)) return(NA)
      if (is.atomic(row)) return(unname(row))
      if (is.list(row)) return(unname(unlist(row, recursive = FALSE, use.names = FALSE)))
      row
    })
    mat <- tryCatch(do.call(rbind, rows), error = function(e) NULL)
  }
  if (is.null(mat)) {
    stop(sprintf("stepwise_covariate_selection: cannot coerce %s to matrix.", field_name))
  }

  if (identical(storage_mode, "double")) {
    suppressWarnings(storage.mode(mat) <- "double")
  } else {
    mat[] <- lapply(as.data.frame(mat, stringsAsFactors = FALSE), as.logical)
    mat <- as.matrix(mat)
    storage.mode(mat) <- "logical"
  }
  mat
}


.covsearch_sumtab_value <- function(gfo, par_name) {
  sumtab <- .covsearch_sumtab_df(gfo)
  if (is.null(sumtab) || nrow(sumtab) == 0 || !"PAR" %in% names(sumtab) || !"VALUE" %in% names(sumtab)) {
    return(NA_real_)
  }
  idx <- which(as.character(sumtab$PAR) == par_name)
  if (length(idx) < 1) {
    return(NA_real_)
  }
  suppressWarnings(as.numeric(sumtab$VALUE[[idx[1]]]))
}


.covsearch_normalize_ruv_entry <- function(ruv_entry, gfo) {
  if (is.null(ruv_entry$YNAME)) {
    stop("stepwise_covariate_selection: each ruv entry must include YNAME.")
  }
  if (is.null(ruv_entry$DVID)) {
    stop("stepwise_covariate_selection: each ruv entry must include DVID.")
  }
  if (is.null(ruv_entry$PRED)) {
    stop("stepwise_covariate_selection: each ruv entry must include PRED.")
  }
  if (is.null(ruv_entry$ERR)) {
    stop("stepwise_covariate_selection: each ruv entry must include ERR.")
  }
  if (is.null(ruv_entry$TRANS)) {
    stop("stepwise_covariate_selection: each ruv entry must include TRANS.")
  }

  pred <- as.character(.covsearch_unwrap_scalar(ruv_entry$PRED))
  err <- tolower(as.character(.covsearch_unwrap_scalar(ruv_entry$ERR)))

  expected_n <- if (grepl("combined", err, ignore.case = TRUE)) 2L else 1L
  init_vec <- suppressWarnings(as.numeric(unlist(ruv_entry$INIT, use.names = FALSE)))
  est_vec <- as.logical(unlist(ruv_entry$EST, use.names = FALSE))

  if (length(init_vec) < expected_n) {
    init_vec <- c(init_vec, rep(NA_real_, expected_n - length(init_vec)))
  } else if (length(init_vec) > expected_n) {
    init_vec <- init_vec[seq_len(expected_n)]
  }
  if (length(est_vec) < expected_n) {
    est_vec <- c(est_vec, rep(TRUE, expected_n - length(est_vec)))
  } else if (length(est_vec) > expected_n) {
    est_vec <- est_vec[seq_len(expected_n)]
  }

  if (grepl("constant", err, ignore.case = TRUE)) {
    if (is.na(init_vec[[1]])) init_vec[[1]] <- .covsearch_sumtab_value(gfo, paste0(pred, "_a"))
    if (is.na(init_vec[[1]])) init_vec[[1]] <- 1
  } else if (grepl("proportional", err, ignore.case = TRUE)) {
    if (is.na(init_vec[[1]])) init_vec[[1]] <- .covsearch_sumtab_value(gfo, paste0(pred, "_b"))
    if (is.na(init_vec[[1]])) init_vec[[1]] <- 0.1
  } else if (grepl("combined", err, ignore.case = TRUE)) {
    if (is.na(init_vec[[1]])) init_vec[[1]] <- .covsearch_sumtab_value(gfo, paste0(pred, "_a"))
    if (is.na(init_vec[[2]])) init_vec[[2]] <- .covsearch_sumtab_value(gfo, paste0(pred, "_b"))
    if (is.na(init_vec[[1]])) init_vec[[1]] <- 1
    if (is.na(init_vec[[2]])) init_vec[[2]] <- 0.1
  }

  list(
    YNAME = as.character(.covsearch_unwrap_scalar(ruv_entry$YNAME)),
    DVID = suppressWarnings(as.numeric(.covsearch_unwrap_scalar(ruv_entry$DVID))),
    TRANS = as.character(.covsearch_unwrap_scalar(ruv_entry$TRANS)),
    PRED = pred,
    ERR = as.character(.covsearch_unwrap_scalar(ruv_entry$ERR)),
    INIT = init_vec,
    EST = est_vec,
    BLQM = if (!is.null(ruv_entry$BLQM)) .covsearch_unwrap_scalar(ruv_entry$BLQM) else NULL
  )
}


.covsearch_resolve_path <- function(path_value, base_dir) {
  if (is.null(path_value) || !is.character(path_value) || length(path_value) != 1 || !nzchar(path_value)) {
    return(path_value)
  }
  candidate <- path_value
  if (!is.null(base_dir) && !grepl("^[A-Za-z]:[/\\\\]|^/", candidate)) {
    candidate <- file.path(base_dir, candidate)
  }
  candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
  candidate
}


.covsearch_preprocess_inputs <- function(gfo, gco) {
  gfo <- .covsearch_json_to_object(gfo, "gfo")
  gco <- .covsearch_json_to_object(gco, "gco")

  gco_src <- attr(gco, "source_json_path")
  gco_base <- if (is.null(gco_src)) NULL else dirname(gco_src)

  gco$headers <- .covsearch_to_header_records(gco$headers)

  theta_df <- .as_covsearch_df(gco$theta, "gco$theta")
  if (!"NAME" %in% names(theta_df)) {
    stop("stepwise_covariate_selection: gco$theta must contain NAME.")
  }
  if (!"INIT" %in% names(theta_df)) {
    stop("stepwise_covariate_selection: gco$theta must contain INIT.")
  }
  if (!"EST" %in% names(theta_df)) {
    theta_df$EST <- TRUE
  }
  if (!"TRANS" %in% names(theta_df)) {
    theta_df$TRANS <- "logNormal"
  }
  theta_df$NAME <- as.character(theta_df$NAME)
  theta_df$INIT <- suppressWarnings(as.numeric(theta_df$INIT))
  theta_df$EST <- as.logical(theta_df$EST)
  theta_df$TRANS <- as.character(theta_df$TRANS)
  gco$theta <- theta_df

  if (!is.null(gco$re)) {
    gco$re$init <- .covsearch_to_matrix(gco$re$init, storage_mode = "double", field_name = "gco$re$init")
    gco$re$est <- .covsearch_to_matrix(gco$re$est, storage_mode = "logical", field_name = "gco$re$est")
  }
  if (!is.null(gco$occ)) {
    gco$occ$init <- .covsearch_to_matrix(gco$occ$init, storage_mode = "double", field_name = "gco$occ$init")
    gco$occ$est <- .covsearch_to_matrix(gco$occ$est, storage_mode = "logical", field_name = "gco$occ$est")
  }

  if (!is.null(gco$ruv)) {
    ruv_in <- gco$ruv
    ruv_list <- if (!is.null(ruv_in$YNAME)) list(ruv_in) else ruv_in
    if (!is.list(ruv_list) || length(ruv_list) == 0) {
      stop("stepwise_covariate_selection: gco$ruv must be a list with at least one entry.")
    }
    gco$ruv <- lapply(ruv_list, .covsearch_normalize_ruv_entry, gfo = gfo)
    if (length(gco$ruv) == 1) {
      gco$ruv <- gco$ruv[[1]]
    }
  }

  gco$model <- .covsearch_resolve_path(gco$model, gco_base)
  gco$data <- .covsearch_resolve_path(gco$data, gco_base)
  gco$path_to_fitter <- .covsearch_resolve_path(gco$path_to_fitter, gco_base)

  list(gfo = gfo, gco = gco)
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
#' @param gfo Generalized fit output object containing at least `COTAB` and `CATAB`.
#' @param gco Generalized control object containing at least `headers` and `theta`.
#' @param output_dir Path where fit projects are written.
#' @param covariates Optional character vector of covariate names to consider.
#' @param parameters Optional character vector of parameter names to consider.
#' @param test_pairs Optional data.frame with columns:
#'   `parameter`, `covariate`, `type`, `reference`, `center`.
#' @param p_forward Numeric in (0,1), default 0.05.
#' @param p_backward Numeric in (0,1), default 0.01.
#' @param fit_function Function used to run model fits (default: `sg_fit`).
#'   Must return a fit-like object consumable by `get_ofv`.
#' @param update_theta_init Logical; when `TRUE`, refreshes theta INIT values from
#'   accepted fit only (never from rejected candidates).
#' @param run_backward Logical; when `TRUE`, run Stage 4 backward elimination after
#'   forward inclusion converges.
#' @param update_theta_init_backward Logical; when `TRUE`, refresh theta INIT only
#'   after accepted backward removals.
#' @param path_to_fitter Optional path to fitter executable.
#'
#' @return A list with final model state, forward/backward summaries, runtime
#'   settings, and execution metadata.
#' @keywords internal
stepwise_covariate_selection <- function(gfo, gco, output_dir = tempdir(),
                                         covariates = NULL,
                                         parameters = NULL,
                                         test_pairs = NULL,
                                         p_forward = 0.05,
                                         p_backward = 0.01,
                                         fit_function = sg_fit,
                                         update_theta_init = FALSE,
                                         run_backward = TRUE,
                                         update_theta_init_backward = FALSE,
                                         path_to_fitter = NULL) {
  prepped <- .covsearch_preprocess_inputs(gfo = gfo, gco = gco)
  gfo <- prepped$gfo
  gco <- prepped$gco

  if (!is.character(output_dir) || length(output_dir) != 1 || !nzchar(output_dir)) {
    stop("stepwise_covariate_selection: output_dir must be a non-empty string.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

  if (!is.function(fit_function)) {
    stop("stepwise_covariate_selection: fit_function must be a function.")
  }
  if (!is.logical(update_theta_init) || length(update_theta_init) != 1 || is.na(update_theta_init)) {
    stop("stepwise_covariate_selection: update_theta_init must be TRUE or FALSE.")
  }
  if (!is.logical(run_backward) || length(run_backward) != 1 || is.na(run_backward)) {
    stop("stepwise_covariate_selection: run_backward must be TRUE or FALSE.")
  }
  if (!is.logical(update_theta_init_backward) ||
      length(update_theta_init_backward) != 1 ||
      is.na(update_theta_init_backward)) {
    stop("stepwise_covariate_selection: update_theta_init_backward must be TRUE or FALSE.")
  }

  if (!is.numeric(p_forward) || length(p_forward) != 1 || is.na(p_forward) ||
      p_forward <= 0 || p_forward >= 1) {
    stop("stepwise_covariate_selection: p_forward must be numeric in (0,1).")
  }
  if (!is.numeric(p_backward) || length(p_backward) != 1 || is.na(p_backward) ||
      p_backward <= 0 || p_backward >= 1) {
    stop("stepwise_covariate_selection: p_backward must be numeric in (0,1).")
  }

  #### For testing ####
  # gco_path <- file.path("scripts", "nlme", "sg-covsearch", "fitted_project", "wrfrn_pk_base_model_GCO.json")
  # gfo_path <- file.path("scripts", "nlme", "sg-covsearch", "fitted_project", "wrfrn_pk_base_model_GFO.json")
  # gco <- jsonlite::fromJSON(gco_path, simplifyVector = TRUE, simplifyDataFrame = TRUE)
  # gfo <- jsonlite::fromJSON(gfo_path, simplifyVector = TRUE, simplifyDataFrame = TRUE)
  # covariates = c("CLCR", "SEX", "VKORC1")
  # parameters = c("CL", "Vd")
  # test_pairs = NULL
  # data_path <- file.path("scripts", "nlme", "sg-covsearch", "fitted_project", "ds-warfarin-pk.csv")
  # data <- read.csv(data_path)
  #####################


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
    #  candidates <- expand.grid(
    #   parameter = parameters,
    #   covariate = covariates,
    #   stringsAsFactors = FALSE
    # )[ , c("parameter", "covariate") ]
    candidates <- do.call(
  rbind,
  lapply(parameters, function(p) {
    data.frame(
      parameter = p,
      covariate = covariates,
      stringsAsFactors = FALSE
    )
  })
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

  base_ofv <- get_ofv(gfo)
  forward_history <- data.frame(
    step = integer(0),
    parameter = character(0),
    covariate = character(0),
    type = character(0),
    df = numeric(0),
    current_ofv = numeric(0),
    candidate_ofv = numeric(0),
    delta_ofv = numeric(0),
    threshold = numeric(0),
    significant = logical(0),
    accepted = logical(0),
    decision = character(0),
    project_name = character(0),
    model_name = character(0),
    candidate_name = character(0),
    status = character(0),
    message = character(0),
    stringsAsFactors = FALSE
  )
  included <- data.frame(
    step = integer(0),
    parameter = character(0),
    covariate = character(0),
    type = character(0),
    df = numeric(0),
    delta_ofv = numeric(0),
    threshold = numeric(0),
    project_name = character(0),
    stringsAsFactors = FALSE
  )
  backward_history <- data.frame(
    step = integer(0),
    parameter = character(0),
    covariate = character(0),
    type = character(0),
    df = numeric(0),
    current_ofv = numeric(0),
    removed_ofv = numeric(0),
    delta_ofv = numeric(0),
    threshold = numeric(0),
    significant = logical(0),
    removed = logical(0),
    decision = character(0),
    project_name = character(0),
    model_name = character(0),
    candidate_name = character(0),
    status = character(0),
    message = character(0),
    stringsAsFactors = FALSE
  )
  backward_removed <- data.frame(
    step = integer(0),
    parameter = character(0),
    covariate = character(0),
    type = character(0),
    df = numeric(0),
    delta_ofv = numeric(0),
    threshold = numeric(0),
    project_name = character(0),
    stringsAsFactors = FALSE
  )

  settings <- list(
    output_dir = output_dir,
    covariates = covariates,
    parameters = parameters,
    p_forward = p_forward,
    p_backward = p_backward,
    update_theta_init = update_theta_init,
    run_backward = run_backward,
    update_theta_init_backward = update_theta_init_backward
  )

  needs_fit_fields <- c("model", "data", "headers", "theta", "ruv", "re", "occ")
  can_run_forward <- all(needs_fit_fields %in% names(gco))
  if (!can_run_forward) {
    forward_selected <- included
    final_covariates <- forward_selected
    final_gco <- gco
    return(list(
      final_gco = final_gco,
      final_covariates = final_covariates,
      forward = list(
        selected = forward_selected,
        history = forward_history
      ),
      backward = list(
        removed = backward_removed,
        retained = final_covariates,
        history = backward_history
      ),
      settings = settings,
      metadata = list(
        forward_ran = FALSE,
        backward_ran = FALSE,
        forward_steps = 0L,
        backward_steps = 0L
      )
    ))
  }

  current_covs <- .covsearch_existing_covs(gco$covs)
  remaining <- candidates
  current_ofv <- base_ofv
  current_gfo <- gfo
  current_theta <- gco_to_theta_tibble(gco, gfo, update_theta_init = update_theta_init)
  base_project_name <- .covsearch_sanitize_name(
    .covsearch_null_coalesce(gco$project_name, "base_model")
  )
  if (!nzchar(base_project_name)) {
    base_project_name <- "base_model"
  }
  step_id <- 1L

  while (nrow(remaining) > 0) {
    tested_rows <- list()
    best_idx <- NA_integer_
    best_delta <- -Inf
    best_fit <- NULL

    for (i in seq_len(nrow(remaining))) {
      cand <- remaining[i, , drop = FALSE]
      is_cont <- identical(cand$type[[1]], "cont")
      cov_type <- if (is_cont) "continuous" else "categorical"
      model_name <- paste0(base_project_name, "_model", sprintf("%02d", step_id))
      candidate_covs <- add_covariate(
        covs_list = current_covs,
        param = cand$parameter[[1]],
        cov = cand$covariate[[1]],
        type = cov_type,
        cov_ref = cand$cov_ref[[1]]
      )

      proj_name <- paste0(
        model_name,
        "_cand", sprintf("%03d", i)
      )

      task_opt_fast_fit <- paste(
  "populationParameters()",
  "individualParameters(method = conditionalMean)",
  "fim(method = StochasticApproximation)",
  "logLikelihood(method = ImportanceSampling)",
  "plotResult(run = false, method = none)",
  sep = "\n"
)
      fit_args <- list(
        model = gco$model,
        data = gco$data,
        headers = gco$headers,
        theta = current_theta,
        ruv = gco$ruv,
        re = gco$re,
        occ = gco$occ,
        covs = candidate_covs,
        project_name = proj_name,
        task_opt =  task_opt_fast_fit, #.covsearch_null_coalesce(gco$task_opt, NULL),
        opt_name = .covsearch_null_coalesce(gco$opt_name, "Monolix"),
        fit = TRUE,
        path_to_save_output = output_dir,
        path_to_fitter = .covsearch_null_coalesce(path_to_fitter, gco$path_to_fitter)
      )



      fit_res <- tryCatch(
        do.call(fit_function, fit_args),
        error = function(e) e
      )

      if (inherits(fit_res, "error")) {
        row_i <- data.frame(
          step = step_id,
          parameter = cand$parameter[[1]],
          covariate = cand$covariate[[1]],
          type = cand$type[[1]],
          df = as.numeric(cand$df[[1]]),
          current_ofv = current_ofv,
          candidate_ofv = NA_real_,
          delta_ofv = NA_real_,
          threshold = stats::qchisq(1 - p_forward, df = as.numeric(cand$df[[1]])),
          significant = FALSE,
          accepted = FALSE,
          decision = "failed",
          project_name = proj_name,
          model_name = model_name,
          candidate_name = proj_name,
          status = "fit_failed",
          message = conditionMessage(fit_res),
          stringsAsFactors = FALSE
        )
        tested_rows[[length(tested_rows) + 1L]] <- row_i
        next
      }

      cand_ofv <- tryCatch(get_ofv(fit_res$GFO), error = function(e) NA_real_)
      thr <- stats::qchisq(1 - p_forward, df = as.numeric(cand$df[[1]]))
      delta <- current_ofv - cand_ofv
      significant <- is.finite(delta) && !is.na(thr) && delta > thr

      row_i <- data.frame(
        step = step_id,
        parameter = cand$parameter[[1]],
        covariate = cand$covariate[[1]],
        type = cand$type[[1]],
        df = as.numeric(cand$df[[1]]),
        current_ofv = current_ofv,
        candidate_ofv = cand_ofv,
        delta_ofv = delta,
        threshold = thr,
        significant = significant,
        accepted = FALSE,
        decision = "rejected",
        project_name = proj_name,
        model_name = model_name,
        candidate_name = proj_name,
        status = "ok",
        message = "",
        stringsAsFactors = FALSE
      )
      tested_rows[[length(tested_rows) + 1L]] <- row_i

      if (isTRUE(significant) && delta > best_delta) {
        best_delta <- delta
        best_idx <- i
        best_fit <- fit_res
      }
    }

    step_df <- do.call(rbind, tested_rows)
    if (!is.na(best_idx)) {
      accepted_mask <- step_df$parameter == remaining$parameter[[best_idx]] &
        step_df$covariate == remaining$covariate[[best_idx]] &
        step_df$status == "ok"
      step_df$accepted[accepted_mask] <- TRUE
      step_df$decision[accepted_mask] <- "accepted"
    }
    forward_history <- rbind(forward_history, step_df)

    if (is.na(best_idx)) {
      break
    }

    accepted <- remaining[best_idx, , drop = FALSE]
    accepted_is_cont <- identical(accepted$type[[1]], "cont")
    accepted_type <- if (accepted_is_cont) "continuous" else "categorical"
    current_covs <- add_covariate(
      covs_list = current_covs,
      param = accepted$parameter[[1]],
      cov = accepted$covariate[[1]],
      type = accepted_type,
      cov_ref = accepted$cov_ref[[1]]
    )
    current_gfo <- best_fit$GFO
    current_ofv <- get_ofv(current_gfo)

    if (isTRUE(update_theta_init)) {
      current_theta <- gco_to_theta_tibble(
        list(theta = current_theta),
        current_gfo,
        update_theta_init = update_theta_init
      )
    }

    included <- rbind(
      included,
      data.frame(
        step = step_id,
        parameter = accepted$parameter[[1]],
        covariate = accepted$covariate[[1]],
        type = accepted$type[[1]],
        df = as.numeric(accepted$df[[1]]),
        delta_ofv = best_delta,
        threshold = stats::qchisq(1 - p_forward, df = as.numeric(accepted$df[[1]])),
        project_name = paste0(base_project_name, "_model", sprintf("%02d", step_id)),
        stringsAsFactors = FALSE
      )
    )

    keep <- !(remaining$parameter == accepted$parameter[[1]] &
      remaining$covariate == accepted$covariate[[1]])
    remaining <- remaining[keep, , drop = FALSE]
    step_id <- step_id + 1L
  }

  forward_steps <- nrow(included)
  backward_ran <- isTRUE(run_backward)

  retained <- included
  if (backward_ran && nrow(retained) > 0) {
    bw_step <- 1L
    repeat {
      if (nrow(retained) == 0) {
        break
      }

      tested_rows <- list()
      removable_term_idx <- NA_integer_
      removable_delta <- Inf
      removable_fit <- NULL

      for (i in seq_len(nrow(retained))) {
        term <- retained[i, , drop = FALSE]
        target_model_idx <- max(nrow(retained) - 1L, 0L)
        model_name <- paste0(base_project_name, "_model", sprintf("%02d", target_model_idx))
        proj_name <- paste0(
          model_name,
          "_bw", sprintf("%02d", bw_step),
          "_cand", sprintf("%03d", i)
        )
        candidate_covs <- remove_covariate(
          covs_list = current_covs,
          param = term$parameter[[1]],
          cov = term$covariate[[1]]
        )

        fit_args <- list(
          model = gco$model,
          data = gco$data,
          headers = gco$headers,
          theta = current_theta,
          ruv = gco$ruv,
          re = gco$re,
          occ = gco$occ,
          covs = candidate_covs,
          project_name = proj_name,
          task_opt = .covsearch_null_coalesce(gco$task_opt, NULL),
          opt_name = .covsearch_null_coalesce(gco$opt_name, "Monolix"),
          fit = TRUE,
          path_to_save_output = output_dir,
          path_to_fitter = .covsearch_null_coalesce(path_to_fitter, gco$path_to_fitter)
        )
        fit_res <- tryCatch(
          do.call(fit_function, fit_args),
          error = function(e) e
        )

        thr <- stats::qchisq(1 - p_backward, df = as.numeric(term$df[[1]]))
        if (inherits(fit_res, "error")) {
          tested_rows[[length(tested_rows) + 1L]] <- data.frame(
            step = bw_step,
            parameter = term$parameter[[1]],
            covariate = term$covariate[[1]],
            type = term$type[[1]],
            df = as.numeric(term$df[[1]]),
            current_ofv = current_ofv,
            removed_ofv = NA_real_,
            delta_ofv = NA_real_,
            threshold = thr,
            significant = NA,
            removed = FALSE,
            decision = "failed",
            project_name = proj_name,
            model_name = model_name,
            candidate_name = proj_name,
            status = "fit_failed",
            message = conditionMessage(fit_res),
            stringsAsFactors = FALSE
          )
          next
        }

        removed_ofv <- tryCatch(get_ofv(fit_res$GFO), error = function(e) NA_real_)
        delta <- removed_ofv - current_ofv
        significant <- is.finite(delta) && !is.na(thr) && delta >= thr
        is_removable <- is.finite(delta) && !is.na(thr) && delta < thr

        tested_rows[[length(tested_rows) + 1L]] <- data.frame(
          step = bw_step,
          parameter = term$parameter[[1]],
          covariate = term$covariate[[1]],
          type = term$type[[1]],
          df = as.numeric(term$df[[1]]),
          current_ofv = current_ofv,
          removed_ofv = removed_ofv,
          delta_ofv = delta,
          threshold = thr,
          significant = significant,
          removed = FALSE,
          decision = if (is_removable) "candidate_remove" else "retain",
          project_name = proj_name,
          model_name = model_name,
          candidate_name = proj_name,
          status = "ok",
          message = "",
          stringsAsFactors = FALSE
        )

        if (is_removable && delta < removable_delta) {
          removable_delta <- delta
          removable_term_idx <- i
          removable_fit <- fit_res
        }
      }

      step_df <- do.call(rbind, tested_rows)
      if (!is.na(removable_term_idx)) {
        removed_mask <- step_df$parameter == retained$parameter[[removable_term_idx]] &
          step_df$covariate == retained$covariate[[removable_term_idx]] &
          step_df$status == "ok"
        step_df$removed[removed_mask] <- TRUE
        step_df$decision[removed_mask] <- "removed"
      }
      backward_history <- rbind(backward_history, step_df)

      if (is.na(removable_term_idx)) {
        break
      }

      removed_term <- retained[removable_term_idx, , drop = FALSE]
      target_model_idx <- max(nrow(retained) - 1L, 0L)
      removed_project <- paste0(
        base_project_name, "_model", sprintf("%02d", target_model_idx)
      )
      removed_thr <- stats::qchisq(1 - p_backward, df = as.numeric(removed_term$df[[1]]))
      backward_removed <- rbind(
        backward_removed,
        data.frame(
          step = bw_step,
          parameter = removed_term$parameter[[1]],
          covariate = removed_term$covariate[[1]],
          type = removed_term$type[[1]],
          df = as.numeric(removed_term$df[[1]]),
          delta_ofv = removable_delta,
          threshold = removed_thr,
          project_name = removed_project,
          stringsAsFactors = FALSE
        )
      )

      current_covs <- remove_covariate(
        covs_list = current_covs,
        param = removed_term$parameter[[1]],
        cov = removed_term$covariate[[1]]
      )
      current_gfo <- removable_fit$GFO
      current_ofv <- get_ofv(current_gfo)
      if (isTRUE(update_theta_init_backward)) {
        current_theta <- gco_to_theta_tibble(
          list(theta = current_theta),
          current_gfo,
          update_theta_init = update_theta_init_backward
        )
      }

      retained <- retained[-removable_term_idx, , drop = FALSE]
      bw_step <- bw_step + 1L
    }
  }

  backward_steps <- nrow(backward_removed)
  final_covariates <- if (backward_ran) retained else included
  final_covs <- current_covs
  final_gco <- gco
  final_gco$covs <- final_covs
  final_gco$theta <- current_theta

  list(
    final_gco = final_gco,
    final_covariates = final_covariates,
    forward = list(
      selected = included,
      history = forward_history
    ),
    backward = list(
      removed = backward_removed,
      retained = retained,
      history = backward_history
    ),
    settings = settings,
    metadata = list(
      forward_ran = TRUE,
      backward_ran = backward_ran,
      forward_steps = forward_steps,
      backward_steps = backward_steps
    )
  )
}

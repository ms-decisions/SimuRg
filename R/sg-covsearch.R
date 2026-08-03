## Author: Alina Melnikova
## First created: 2026-05-19
## Description: covariate search tool
## Keywords: SimuRg, covsearch

#' Calculate objective function value from fit output
#'
#' Use OFV = -2*LL from gfo$OFV
#'
#' @param gfo [GFO] containing `OFV`.
#'
#' @return Numeric scalar OFV.
#' @noRd
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
#' @param gfo [GFO].
#'
#' @return Data frame with at least standard SUMTAB columns when available.
#' @noRd
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
#' @param gco [GCO] containing `theta`.
#' @param gfo [GFO] containing `SUMTAB`.
#' @param update_theta_init Logical; when `TRUE`, update INIT from `gfo$SUMTAB`.
#'
#' @return Tibble with theta columns and updated `INIT`.
#' @noRd
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
#' @noRd
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
      INIT = 0,
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
#' @noRd
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

#' Coerce covariate search input to a data frame
#'
#' @param x Object. Input object expected to be coercible to a `data.frame`.
#' @param name Character. Label used in error messages when `x` is missing.
#'
#' @return Data.frame. A normalized `data.frame` representation of `x`.
#' @noRd
.as_covsearch_df <- function(x, name) {
  if (is.null(x)) {
    stop(sprintf("sg_covsearch: %s is missing.", name))
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

#' Normalize covariate type labels
#'
#' @param x Character. Covariate type values to normalize.
#'
#' @return Character. Normalized covariate types (`cont`, `cat`, or `NA`).
#' @noRd
.norm_cov_type <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- ifelse(x %in% c("cont", "continuous"), "cont",
                ifelse(x %in% c("cat", "categorical"), "cat", NA_character_))
  out
}

#' Compute most frequent value with deterministic tie-break
#'
#' @param x Vector. Values used to compute the mode.
#'
#' @return Character. Most frequent non-missing value, choosing the
#'   lexicographically smallest value on ties.
#' @noRd
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

#' Compute median while preserving all-missing input
#'
#' @param x Vector. Numeric-like values to summarize.
#'
#' @return Numeric. Median of non-missing values, or `NA_real_` when all values
#'   are missing.
#' @noRd
.na_safe_median <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) {
    return(NA_real_)
  }
  stats::median(x, na.rm = TRUE)
}

#' Return fallback value when input is `NULL`
#'
#' @param x Object. Primary value.
#' @param y Object. Fallback value returned when `x` is `NULL`.
#'
#' @return Object. `x` when non-`NULL`, otherwise `y`.
#' @noRd
.covsearch_null_coalesce <- function(x, y) {
  if (is.null(x)) y else x
}

#' Sanitize strings for safe identifier usage
#'
#' @param x Character. Value to sanitize.
#'
#' @return Character. Sanitized identifier containing only alphanumeric
#'   characters and underscores.
#' @noRd
.covsearch_sanitize_name <- function(x) {
  out <- gsub("[^[:alnum:]_]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  if (!nzchar(out)) "x" else out
}

#' Validate and normalize existing covariate records
#'
#' @param x List. Candidate covariate records.
#'
#' @return List. Input list when it is a non-empty list of lists; otherwise an
#'   empty list.
#' @noRd
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



#####

#' Stepwise covariate modeling (SCM) search
#'
#' Runs forward inclusion and optional backward elimination of
#' parameter–covariate relationships, selecting terms by likelihood-ratio tests
#' against `p_forward` and `p_backward`.
#'
#' @param gfo [GFO] or character. Baseline fit object, or path to a GFO
#'   `.json`/`.RData`, containing at least `COTAB` and `CATAB`.
#' @param gco [GCO] or character. Control object (list), or path to a GCO
#'   `.json`/`.RData`, containing at least `headers` and `theta`.
#' @param output_dir Character. Directory where fit projects and search history
#'   files are written. Default value is `tempdir()`.
#' @param covariates Character vector or `NULL`. Covariate names to consider.
#'   When `NULL`, all covariates from `gco$headers` are used. Default value is
#'   `NULL`.
#' @param parameters Character vector or `NULL`. Parameter names to consider.
#'   When `NULL`, all parameters from `gco$theta` are used. Default value is
#'   `NULL`.
#' @param test_pairs Data.frame or `NULL`. Candidate pairs with columns
#'   `parameter`, `covariate`, `type`, `reference`, and `center`. When `NULL`,
#'   all parameter–covariate combinations are generated. Default value is `NULL`.
#' @param p_forward Numeric in (0,1). Significance level for forward inclusion.
#'   Default value is `0.05`.
#' @param p_backward Numeric in (0,1). Significance level for backward
#'   elimination. Default value is `0.01`.
#' @param fit_function Function. Fitting function that returns a fit-like object
#'   consumable by `get_ofv`. Default value is `sg_fit`.
#' @param update_theta_init Logical. If `TRUE`, refreshes theta `INIT` values from
#'   accepted forward fits only (never from rejected candidates). Default value
#'   is `FALSE`.
#' @param run_backward Logical. If `TRUE`, runs Stage 4 backward elimination after
#'   forward inclusion converges. Default value is `TRUE`.
#' @param update_theta_init_backward Logical. If `TRUE`, refreshes theta `INIT`
#'   only after accepted backward removals. Default value is `FALSE`.
#' @param path_to_fitter Character or `NULL`. Path to the fitter executable.
#'   When `NULL`, `gco$path_to_fitter` is used if present. Default value is
#'   `NULL`.
#'
#' @export
#' @examples
#' model_path <- tempfile(fileext = ".txt")
#' data_path <- tempfile(fileext = ".csv")
#' writeLines(c("[LONGITUDINAL]", "input = {CL, V}", "PK:"), model_path)
#' writeLines(c("ID,TIME,DV,WT", "1,0,0,70", "2,0,0,80"), data_path)
#'
#' gco <- list(
#'   model = model_path,
#'   data = data_path,
#'   headers = list(
#'     list(name = "ID", use = "identifier", type = NULL),
#'     list(name = "TIME", use = "time", type = NULL),
#'     list(name = "DV", use = "observation", type = "continuous"),
#'     list(name = "WT", use = "covariate", type = "continuous")
#'   ),
#'   theta = data.frame(
#'     NAME = c("CL", "V"),
#'     TRANS = c("logNormal", "logNormal"),
#'     INIT = c(0.2, 20),
#'     EST = c(TRUE, TRUE),
#'     stringsAsFactors = FALSE
#'   ),
#'   ruv = list(dummy = TRUE),
#'   re = list(dummy = TRUE),
#'   occ = list(dummy = TRUE),
#'   covs = list(),
#'   project_name = "base_model"
#' )
#'
#' gfo <- list(
#'   OFV = data.frame(LL = 100),
#'   SUMTAB = data.frame(PAR = character(0), VALUE = numeric(0), stringsAsFactors = FALSE),
#'   COTAB = data.frame(ID = 1:4, WT = c(70, 80, 90, 75), stringsAsFactors = FALSE),
#'   CATAB = data.frame(ID = 1:4, stringsAsFactors = FALSE)
#' )
#'
#' mock_fit <- function(model, data, headers, theta, ruv, re, occ, covs, project_name,
#'                      task_opt = NULL, opt_name = "Monolix", fit = TRUE,
#'                      path_to_save_output = NULL, path_to_fitter = NULL) {
#'   ofv <- if (grepl("_001$", project_name)) 95 else 99
#'   list(
#'     GFO = list(
#'       OFV = data.frame(LL = ofv),
#'       SUMTAB = data.frame(PAR = character(0), VALUE = numeric(0), stringsAsFactors = FALSE),
#'       COTAB = data.frame(dummy = 1),
#'       CATAB = data.frame(dummy = 1)
#'     )
#'   )
#' }
#'
#' result <- sg_covsearch(
#'   gfo = gfo,
#'   gco = gco,
#'   output_dir = tempfile("covsearch-example-"),
#'   covariates = "WT",
#'   parameters = "CL",
#'   run_backward = FALSE,
#'   fit_function = mock_fit
#' )
#'
#' result$forward$selected[, c("parameter", "covariate")]
#'
#' @return A list with `final_gco`, `final_covariates`, forward/backward
#'   summaries, runtime `settings`, and execution `metadata`.

sg_covsearch<- function(gfo, gco, output_dir = NULL,
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


  if (!is.function(fit_function)) {
    stop("sg_covsearch: fit_function must be a function.")
  }
  if (!is.logical(update_theta_init) || length(update_theta_init) != 1 || is.na(update_theta_init)) {
    stop("sg_covsearch: update_theta_init must be TRUE or FALSE.")
  }
  if (!is.logical(run_backward) || length(run_backward) != 1 || is.na(run_backward)) {
    stop("sg_covsearch: run_backward must be TRUE or FALSE.")
  }
  if (!is.logical(update_theta_init_backward) ||
      length(update_theta_init_backward) != 1 ||
      is.na(update_theta_init_backward)) {
    stop("sg_covsearch: update_theta_init_backward must be TRUE or FALSE.")
  }

  if (!is.numeric(p_forward) || length(p_forward) != 1 || is.na(p_forward) ||
      p_forward <= 0 || p_forward >= 1) {
    stop("sg_covsearch: p_forward must be numeric in (0,1).")
  }
  if (!is.numeric(p_backward) || length(p_backward) != 1 || is.na(p_backward) ||
      p_backward <= 0 || p_backward >= 1) {
    stop("sg_covsearch: p_backward must be numeric in (0,1).")
  }

  gco_dir <- NULL
  if (is.character(gco)) {
    if (length(gco) != 1 || !nzchar(gco)) {
      stop("sg_covsearch: gco path must be a non-empty string.")
    }
    if (!file.exists(gco)) {
      stop("sg_covsearch: gco file does not exist: ", gco)
    }
    gco_path <- normalizePath(gco, winslash = "/", mustWork = TRUE)
    gco_dir <- dirname(gco_path)
    gco <- read_smrg_ctrl(gco_path)
    #gco <- jsonlite::fromJSON(gco_path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
  } else if (!is.list(gco)) {
    stop("sg_covsearch: gco must be a list or path to a GCO file.")
  }

  if (is.null(output_dir) || (length(output_dir) == 1 && is.na(output_dir))) {
    if (is.list(gco)) {
      stop("sg_covsearch: output_dir must be provided when gco is a list.")
    }
    if (is.null(gco_dir) || !nzchar(gco_dir)) {
      stop("sg_covsearch: output_dir is NULL/NA, but gco_dir is unavailable. Pass output_dir explicitly or provide gco as a file path.")
    }
    output_dir <- gco_dir
  }
  if (!is.character(output_dir) || length(output_dir) != 1 || is.na(output_dir) || !nzchar(output_dir)) {
    stop("sg_covsearch: output_dir must be a non-empty string, NULL, or NA.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

  .is_abs_path <- function(x) {
    grepl("^(?:[A-Za-z]:[/\\\\]|/|\\\\\\\\)", x)
  }

  .resolve_gco_ref_path <- function(x, field_name) {
    if (is.null(x) || !is.character(x) || length(x) != 1 || !nzchar(x)) {
      stop("sg_covsearch: gco$", field_name, " must be a non-empty string.")
    }
    if (!is.null(gco_dir) && !.is_abs_path(x)) {
      x <- file.path(gco_dir, x)
    }
    normalizePath(x, winslash = "/", mustWork = TRUE)
  }

  if (is.character(gfo)) {
    if (length(gfo) != 1 || !nzchar(gfo)) {
      stop("sg_covsearch: gfo path must be a non-empty string.")
    }
    if (!file.exists(gfo)) {
      stop("sg_covsearch: gfo file does not exist: ", gfo)
    }
    gfo_path <- normalizePath(gfo, winslash = "/", mustWork = TRUE)
    gfo_dir <- dirname(gfo_path)
    gfo <- read_smrg_obj(gfo_path)
    #gfo <- jsonlite::fromJSON(gfo_path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
  } else if (!is.list(gfo)) {
    stop("sg_covsearch: gfo must be a list or path to a GFO file.")
  }

  headers_df <- .as_covsearch_df(gco$headers, "gco$headers")
  theta_df <- .as_covsearch_df(gco$theta, "gco$theta")
  cotab_df <- .as_covsearch_df(gfo$COTAB, "gfo$COTAB")
  catab_df <- .as_covsearch_df(gfo$CATAB, "gfo$CATAB")

  model <- .resolve_gco_ref_path(gco$model, "model")
  data <- .resolve_gco_ref_path(gco$data, "data")

  if (!all(c("name", "use", "type") %in% names(headers_df))) {
    stop("sg_covsearch: gco$headers must contain name/use/type.")
  }
  if (!"NAME" %in% names(theta_df)) {
    stop("sg_covsearch: gco$theta must contain NAME.")
  }

  cov_headers <- headers_df[headers_df$use == "covariate", c("name", "type"), drop = FALSE]
  if (nrow(cov_headers) == 0) {
    stop("sg_covsearch: no covariates found in gco$headers.")
  }
  cov_headers$type <- .norm_cov_type(cov_headers$type)
  if (any(is.na(cov_headers$type))) {
    bad <- unique(cov_headers$name[is.na(cov_headers$type)])
    stop(
      sprintf(
        "sg_covsearch: unsupported covariate type for: %s",
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
      stop("sg_covsearch: covariates must be a character vector.")
    }
    invalid_cov <- setdiff(unique(covariates), valid_covariates)
    if (length(invalid_cov) > 0) {
      stop(sprintf(
        "sg_covsearch: unknown covariates: %s",
        paste(invalid_cov, collapse = ", ")
      ))
    }
    covariates <- unique(covariates)
  }

  if (is.null(parameters)) {
    parameters <- valid_parameters
  } else {
    if (!is.character(parameters)) {
      stop("sg_covsearch: parameters must be a character vector.")
    }
    invalid_par <- setdiff(unique(parameters), valid_parameters)
    if (length(invalid_par) > 0) {
      stop(sprintf(
        "sg_covsearch: unknown parameters: %s",
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
      stop("sg_covsearch: test_pairs must be a data.frame or NULL.")
    }
    needed_cols <- c("parameter", "covariate", "type")
    missing_cols <- setdiff(needed_cols, names(test_pairs))
    if (length(missing_cols) > 0) {
      stop(sprintf(
        "sg_covsearch: test_pairs missing required columns: %s",
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
        "sg_covsearch: dropped %d invalid test_pairs row(s).",
        bad_n
      ))
    }
    candidates <- candidates[valid_row, c("parameter", "covariate", "type", "reference", "center"), drop = FALSE]
  }

  if (nrow(candidates) == 0) {
    stop("sg_covsearch: no valid candidate pairs remain.")
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
    fh_csv <- file.path(output_dir, "forward_history.csv")
    bh_csv <- file.path(output_dir, "backward_history.csv")
    utils::write.csv(as.data.frame(forward_history), fh_csv, row.names = FALSE, na = "")
    utils::write.csv(as.data.frame(backward_history), bh_csv, row.names = FALSE, na = "")
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
  fw_checked_idx <- 0L
  bw_checked_idx <- 0L

  while (nrow(remaining) > 0) {
    tested_rows <- list()
    best_idx <- NA_integer_
    best_delta <- -Inf
    best_fit <- NULL

    for (i in seq_len(nrow(remaining))) {
      fw_checked_idx <- fw_checked_idx + 1L
      cand <- remaining[i, , drop = FALSE]
      is_cont <- identical(cand$type[[1]], "cont")
      cov_type <- if (is_cont) "continuous" else "categorical"
      proj_name <- paste0(base_project_name, "_fw_model_", sprintf("%03d", fw_checked_idx))
      candidate_covs <- add_covariate(
        covs_list = current_covs,
        param = cand$parameter[[1]],
        cov = cand$covariate[[1]],
        type = cov_type,
        cov_ref = cand$cov_ref[[1]]
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
        model =  model,
        data = data,
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
        error = function(e) {
          stop(
            sprintf(
              "sg_covsearch: fit failed in forward step %d for parameter '%s' and covariate '%s' (project '%s'): %s",
              step_id,
              cand$parameter[[1]],
              cand$covariate[[1]],
              proj_name,
              conditionMessage(e)
            ),
            call. = FALSE
          )
        }
      )

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
        step_df$covariate == remaining$covariate[[best_idx]]
      step_df$accepted[accepted_mask] <- TRUE
      step_df$decision[accepted_mask] <- "accepted"
    }
    forward_history <- rbind(forward_history, step_df)

    if (is.na(best_idx)) {
      break
    }

    accepted_project <- step_df$project_name[
      step_df$parameter == remaining$parameter[[best_idx]] &
        step_df$covariate == remaining$covariate[[best_idx]]
    ]
    accepted_project <- as.character(accepted_project)
    accepted_project <- accepted_project[!is.na(accepted_project) & nzchar(accepted_project)]
    if (length(accepted_project) == 0) {
      accepted_project <- paste0(base_project_name, "_fw_model_", sprintf("%03d", fw_checked_idx))
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
        project_name = accepted_project,
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
        bw_checked_idx <- bw_checked_idx + 1L
        term <- retained[i, , drop = FALSE]
        proj_name <- paste0(base_project_name, "_bw_model_", sprintf("%03d", bw_checked_idx))
        candidate_covs <- remove_covariate(
          covs_list = current_covs,
          param = term$parameter[[1]],
          cov = term$covariate[[1]]
        )

        fit_args <- list(
          model = model,
          data = data,
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
          error = function(e) {
            stop(
              sprintf(
                "sg_covsearch: fit failed in backward step %d for parameter '%s' and covariate '%s' (project '%s'): %s",
                bw_step,
                term$parameter[[1]],
                term$covariate[[1]],
                proj_name,
                conditionMessage(e)
              ),
              call. = FALSE
            )
          }
        )

        thr <- stats::qchisq(1 - p_backward, df = as.numeric(term$df[[1]]))

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
          step_df$covariate == retained$covariate[[removable_term_idx]]
        step_df$removed[removed_mask] <- TRUE
        step_df$decision[removed_mask] <- "removed"
      }
      backward_history <- rbind(backward_history, step_df)

      if (is.na(removable_term_idx)) {
        break
      }

      removed_term <- retained[removable_term_idx, , drop = FALSE]
      removed_project <- step_df$project_name[
        step_df$parameter == retained$parameter[[removable_term_idx]] &
          step_df$covariate == retained$covariate[[removable_term_idx]]
      ]
      removed_project <- as.character(removed_project)
      removed_project <- removed_project[!is.na(removed_project) & nzchar(removed_project)]
      if (length(removed_project) == 0) {
        removed_project <- paste0(base_project_name, "_bw_model_", sprintf("%03d", bw_checked_idx))
      } else {
        removed_project <- removed_project[[1]]
      }
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
  fh_csv <- file.path(output_dir, "forward_history.csv")
  bh_csv <- file.path(output_dir, "backward_history.csv")
  utils::write.csv(as.data.frame(forward_history), fh_csv, row.names = FALSE, na = "")
  utils::write.csv(as.data.frame(backward_history), bh_csv, row.names = FALSE, na = "")

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

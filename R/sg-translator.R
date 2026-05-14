## Author: Victor Sokolov
## First created: 2026-03-05
## Description: Translate structural models between rxode2 and MLXTRAN syntax
## Keywords: SimuRg, monolix, rxode2, translator, converter

#' Translate structural models between rxode2 and MLXTRAN syntax
#'
#' @description
#' Converts structural pharmacometric model files between rxode2 (R package)
#' and MLXTRAN (Monolix) syntax. Supports ODE-based models and Monolix pkmodel macros.
#'
#' @param input_path Character string. Path to the input structural model file (.txt).
#' @param to Character string. Target format: \code{"mlxtran"} or \code{"rxode"}.
#' @param output_path Character string. Path to write the translated model file (.txt).
#' @param dm_list List of dosing macros with \code{cmt} and \code{adm} vectors.
#'   For example, \code{list(cmt = 2, adm = 1)} or \code{list(cmt = c(2, 1), adm = c(1, 2))}.
#'   Used to generate \code{iv()} dosing macros in MLXTRAN output.
#' @param regressors Character vector of regressor variable names.
#'   For example, \code{c("WTBL", "ADAN")} generates \code{WTBL = \{use = regressor\}} in MLXTRAN.
#' @param output_vars Character vector of output variable names for the MLXTRAN \code{OUTPUT:} section.
#'   For example, \code{c("Cc_mgL", "CSF1_ngmL")} generates \code{output = \{Cc_mgL, CSF1_ngmL\}}.
#'   Required when \code{to = "mlxtran"}; ignored when \code{to = "rxode"}.
#' @param macros Logical. If \code{TRUE}, attempt to convert rxode2 model to pkmodel macro
#'   format in MLXTRAN (only for simple 1- or 2-compartment oral PK models). Default \code{TRUE}.
#' @param stiff Logical. If \code{TRUE}, add \code{odeType = stiff} to the MLXTRAN EQUATION section.
#'   Default \code{TRUE}.
#'
#' @return Invisibly returns the translated model as a character string.
#'   The translated model is also written to \code{output_path}.
#'
#' @details
#' The function detects the source format automatically from file content.
#'
#' \strong{rxode2 format} uses \code{# [INPUT]}, \code{# [MODEL]}, \code{# [OUTPUT]} section markers,
#' \code{d/dt(X)} ODE notation, \code{f(X)} for bioavailability, and \code{alag(X)} for lag time.
#'
#' \strong{MLXTRAN format} uses \code{[LONGITUDINAL]}, \code{PK:}/\code{EQUATION:}/\code{OUTPUT:} sections,
#' \code{ddt_X} ODE notation, \code{compartment()}/\code{iv()} dosing macros, and optionally
#' \code{pkmodel()} macros.
#'
#' When converting rxode2 to MLXTRAN with \code{macros = TRUE}, the function checks if the model
#' structure matches a simple \code{pkmodel()} pattern (1- or 2-compartment with first-order oral
#' absorption). If not, it falls back to full ODE format with \code{compartment()}/\code{iv()} macros.
#'
#' MLXTRAN models using \code{absorption()}/\code{elimination()}/\code{peripheral()} macros (type 2)
#' are expanded to explicit ODEs when converted to rxode2.
#'
#' @examples
#' \donttest{
#' # Convert rxode2 model to MLXTRAN ODE format
#' rxode_model_hv_iv <- system.file("extdata", "models", "rxode",
#'                                  "model_PK_hv_iv.txt", package = "SimuRg")
#' rxode_model_1c <- system.file("extdata", "models", "rxode",
#'                               "model_PK_1c.txt", package = "SimuRg")
#' monolix_model_hv_iv <- system.file("extdata", "models", "monolix",
#'                                    "model_PK_hv_iv.txt", package = "SimuRg")
#' monolix_model_1c <- system.file("extdata", "models", "monolix",
#'                                 "model_PK_1c.txt", package = "SimuRg")
#'
#' path_to_save <- tempdir()
#'
#' sg_translator(rxode_model_hv_iv, to = "mlxtran",
#'               output_path = normalizePath(file.path(path_to_save,
#'                                                     "model_PK_hv_iv_mlx.txt"),
#'                                           mustWork = FALSE),
#'               output_vars = "Cc_mgL",
#'               dm_list = list(cmt = 1, adm = 1), stiff = TRUE)
#'
#' # Convert rxode2 model to MLXTRAN with pkmodel macros
#' sg_translator(rxode_model_1c, to = "mlxtran",
#'               output_path = normalizePath(file.path(path_to_save,
#'                                                     "model_PK_1c_mlx.txt"),
#'                                           mustWork = FALSE),
#'               output_vars = "Cc_nM", macros = TRUE)
#'
#' # Convert MLXTRAN model to rxode2
#' sg_translator(monolix_model_hv_iv, to = "rxode",
#'               output_path = normalizePath(file.path(path_to_save,
#'                                                     "model_PK_hv_iv_rx.txt"),
#'                                            mustWork = FALSE))
#' }
#'
#' @export
sg_translator <- function(input_path, to, output_path, dm_list = NULL,
                          regressors = NULL, output_vars = NULL,
                          macros = TRUE, stiff = TRUE) {

  input_path <- normalizePath(input_path, mustWork = FALSE)
  if(!file.exists(input_path)) {
    stop("Input model file does not exist. Check file existance or try to use absolute path")
  }
  raw <- readLines(input_path, warn = FALSE)
  lines <- stringr::str_trim(raw, side = "right")
  is_mlxtran <- any(stringr::str_detect(lines, "^\\s*\\[LONGITUDINAL\\]"))

  if (to == "mlxtran" && !is_mlxtran && is.null(output_vars)) {
    stop("'output_vars' is required when converting to MLXTRAN (e.g., output_vars = c(\"Cc_mgL\"))")
  }

  if (to == "rxode" && is_mlxtran) {
    result <- .tr_mlxtran_to_rxode(lines)
  } else if (to == "mlxtran" && !is_mlxtran) {
    result <- .tr_rxode_to_mlxtran(lines, dm_list, regressors, output_vars, macros, stiff)
  } else {
    stop("Source format already matches target or format is unrecognized")
  }

  output <- paste(result, collapse = "\n")
  writeLines(output, output_path)
  invisible(output)
}


# =============================================================================
# COMMENT AND LINE CONVERSION HELPERS
# =============================================================================

.tr_comment_mlx_to_rx <- function(line) {
  line %>%
    stringr::str_replace(";;;", "###") %>%
    stringr::str_replace(";;", "##") %>%
    stringr::str_replace(";", "#")
}

.tr_comment_rx_to_mlx <- function(line) {
  line %>%
    stringr::str_replace("###", ";;;") %>%
    stringr::str_replace("##", ";;") %>%
    stringr::str_replace("#", ";")
}

.tr_mlx_code_to_rx <- function(line, add_semicolon = TRUE) {
  trimmed <- stringr::str_trim(line)
  if (trimmed == "") return("")

  sc_pos <- stringr::str_locate(trimmed, ";")[1, "start"]

  if (!is.na(sc_pos) && sc_pos == 1) {
    return(.tr_comment_mlx_to_rx(trimmed))
  }

  if (!is.na(sc_pos)) {
    code_part <- stringr::str_sub(trimmed, 1, sc_pos - 1) %>% stringr::str_trim("right")
    comment_part <- stringr::str_sub(trimmed, sc_pos) %>% .tr_comment_mlx_to_rx()
    sep <- if (add_semicolon) ";\t\t\t\t" else "\t\t\t\t"
    return(paste0(code_part, sep, comment_part))
  }

  if (add_semicolon) {
    return(paste0(trimmed, ";"))
  }
  trimmed
}

.tr_rx_code_to_mlx <- function(line) {
  trimmed <- stringr::str_trim(line)
  if (trimmed == "") return("")

  if (stringr::str_detect(trimmed, "^#")) {
    return(.tr_comment_rx_to_mlx(trimmed))
  }

  hash_pos <- stringr::str_locate(trimmed, "#")[1, "start"]

  if (!is.na(hash_pos) && hash_pos > 1) {
    code_part <- stringr::str_sub(trimmed, 1, hash_pos - 1) %>%
      stringr::str_trim("right") %>%
      stringr::str_remove(";$")
    comment_part <- stringr::str_sub(trimmed, hash_pos) %>%
      .tr_comment_rx_to_mlx()
    return(paste0(code_part, "\t\t\t\t", comment_part))
  }

  stringr::str_remove(trimmed, ";\\s*$")
}

.tr_extract_macro_arg <- function(macro_str, arg_name) {
  pattern <- paste0("(?:,|\\()\\s*", arg_name, "\\s*=\\s*")
  m <- regexpr(pattern, macro_str, perl = TRUE)
  if (m == -1) return(NA_character_)

  start <- m + attr(m, "match.length")
  depth <- 0
  chars <- strsplit(macro_str, "")[[1]]
  end <- start

  for (j in start:length(chars)) {
    ch <- chars[j]
    if (ch == "(") depth <- depth + 1
    if (ch == ")") {
      if (depth == 0) { end <- j - 1; break }
      depth <- depth - 1
    }
    if (ch == "," && depth == 0) { end <- j - 1; break }
    end <- j
  }

  stringr::str_trim(substr(macro_str, start, end))
}

.tr_extract_rhs <- function(line_trimmed) {
  full_rhs <- stringr::str_match(line_trimmed, "=\\s*(.*)")[, 2]
  hash_pos <- stringr::str_locate(full_rhs, "#")[1, "start"]
  sc_pos <- stringr::str_locate(full_rhs, ";")[1, "start"]

  if (!is.na(hash_pos) && hash_pos > 1) {
    rhs <- stringr::str_sub(full_rhs, 1, hash_pos - 1) %>%
      stringr::str_trim() %>% stringr::str_remove(";$") %>% stringr::str_trim()
    comment <- stringr::str_sub(full_rhs, hash_pos)
  } else if (!is.na(sc_pos)) {
    rhs <- stringr::str_sub(full_rhs, 1, sc_pos - 1) %>% stringr::str_trim()
    remaining <- stringr::str_sub(full_rhs, sc_pos + 1) %>% stringr::str_trim()
    comment <- if (nchar(remaining) > 0) paste0("; ", remaining) else ""
  } else {
    rhs <- stringr::str_trim(full_rhs)
    comment <- ""
  }

  list(rhs = rhs, comment = comment)
}


# =============================================================================
# MLXTRAN → RXODE2 CONVERSION
# =============================================================================

.tr_mlxtran_to_rxode <- function(lines) {
  n <- length(lines)

  idx_desc <- which(stringr::str_detect(lines, "^\\s*DESCRIPTION:"))[1]
  idx_long <- which(stringr::str_detect(lines, "^\\s*\\[LONGITUDINAL\\]"))[1]
  idx_pk   <- which(stringr::str_detect(lines, "^\\s*PK:"))[1]
  idx_eq   <- which(stringr::str_detect(lines, "^\\s*EQUATION:"))[1]
  idx_out  <- which(stringr::str_detect(lines, "^\\s*OUTPUT:"))[1]

  # --- Description ---
  desc_rx <- character(0)
  if (!is.na(idx_desc)) {
    desc_end <- ifelse(!is.na(idx_long), idx_long - 1, n)
    if (idx_desc + 1 <= desc_end) {
      desc_raw <- lines[(idx_desc + 1):desc_end]
      desc_raw <- desc_raw[stringr::str_trim(desc_raw) != ""]
      desc_rx <- vapply(desc_raw, .tr_comment_mlx_to_rx, character(1), USE.NAMES = FALSE)
    }
  }

  # --- Input parameters and regressors ---
  input_idx <- which(stringr::str_detect(lines, "^\\s*input\\s*="))[1]
  input_str <- stringr::str_extract(lines[input_idx], "\\{[^}]+\\}") %>%
    stringr::str_remove_all("[{}]") %>%
    stringr::str_split(",") %>% .[[1]] %>% stringr::str_trim()

  reg_idx <- which(stringr::str_detect(lines, "=\\s*\\{\\s*use\\s*=\\s*regressor"))
  reg_vars <- if (length(reg_idx) > 0) {
    stringr::str_extract(lines[reg_idx], "^\\s*\\w+") %>% stringr::str_trim()
  } else {
    character(0)
  }
  model_params <- setdiff(input_str, reg_vars)

  # --- Parse PK section ---
  compartments <- list()
  iv_macros <- list()
  pk_equations <- character(0)
  pkmodel_line <- NULL
  pkmodel_output_var <- NULL
  absorption_macros <- list()
  peripheral_macros <- list()
  elimination_macros <- list()

  if (!is.na(idx_pk)) {
    pk_end <- min(c(idx_eq, idx_out, n + 1)[!is.na(c(idx_eq, idx_out, n + 1))]) - 1
    pk_raw <- lines[(idx_pk + 1):pk_end]

    for (pl in pk_raw) {
      pl_t <- stringr::str_trim(pl)
      if (pl_t == "" || stringr::str_detect(pl_t, "^;")) next

      if (stringr::str_detect(pl_t, "^compartment\\s*\\(")) {
        cmt <- as.integer(stringr::str_extract(pl_t, "(?<=cmt\\s{0,5}=\\s{0,5})\\d+"))
        amt <- stringr::str_match(pl_t, "amount\\s*=\\s*(\\w+)")[, 2]
        vol <- stringr::str_match(pl_t, "volume\\s*=\\s*(\\w+)")[, 2]
        conc <- stringr::str_match(pl_t, "concentration\\s*=\\s*(\\w+)")[, 2]
        compartments[[length(compartments) + 1]] <- list(
          cmt = cmt, amount = amt, volume = vol, concentration = conc
        )

      } else if (stringr::str_detect(pl_t, "^iv\\s*\\(")) {
        macro_part <- stringr::str_remove(pl_t, ";.*$")
        cmt <- as.integer(stringr::str_extract(macro_part, "(?<=cmt\\s{0,5}=\\s{0,5})\\d+"))
        adm <- as.integer(stringr::str_extract(macro_part, "(?<=adm\\s{0,5}=\\s{0,5})\\d+"))
        p_val <- .tr_extract_macro_arg(macro_part, "p")
        tlag_val <- .tr_extract_macro_arg(macro_part, "Tlag")
        inline_comment <- stringr::str_extract(pl_t, ";.*$")
        iv_macros[[length(iv_macros) + 1]] <- list(
          cmt = cmt, adm = adm, p = p_val, Tlag = tlag_val, comment = inline_comment
        )

      } else if (stringr::str_detect(pl_t, "pkmodel\\s*\\(")) {
        pkmodel_line <- pl_t
        pkmodel_output_var <- stringr::str_extract(pl_t, "^\\w+")

      } else if (stringr::str_detect(pl_t, "^absorption\\s*\\(")) {
        cmt <- as.integer(stringr::str_extract(pl_t, "(?<=cmt\\s{0,5}=\\s{0,5})\\d+"))
        tk0 <- .tr_extract_macro_arg(pl_t, "Tk0")
        ka_val <- .tr_extract_macro_arg(pl_t, "ka")
        tlag_val <- .tr_extract_macro_arg(pl_t, "Tlag")
        p_val <- .tr_extract_macro_arg(pl_t, "p")
        absorption_macros[[length(absorption_macros) + 1]] <- list(
          cmt = cmt,
          Tk0 = if (!is.na(tk0)) tk0 else NULL,
          ka  = if (!is.na(ka_val)) ka_val else NULL,
          Tlag = if (!is.na(tlag_val)) tlag_val else NULL,
          p   = if (!is.na(p_val)) p_val else NULL
        )

      } else if (stringr::str_detect(pl_t, "^peripheral\\s*\\(")) {
        args <- stringr::str_extract(pl_t, "\\([^)]+\\)") %>%
          stringr::str_remove_all("[()]") %>%
          stringr::str_split(",") %>% .[[1]] %>% stringr::str_trim()
        peripheral_macros[[length(peripheral_macros) + 1]] <- list(k12 = args[1], k21 = args[2])

      } else if (stringr::str_detect(pl_t, "^elimination\\s*\\(")) {
        cmt <- as.integer(stringr::str_extract(pl_t, "(?<=cmt\\s{0,5}=\\s{0,5})\\d+"))
        inner <- stringr::str_extract(pl_t, "\\([^)]+\\)") %>% stringr::str_remove_all("[()]")
        parts <- stringr::str_split(inner, ",")[[1]] %>% stringr::str_trim()
        k_val <- parts[!stringr::str_detect(parts, "cmt")]
        if (length(k_val) > 0) {
          k_val <- stringr::str_remove(k_val, "^k\\s*=\\s*") %>% stringr::str_trim()
        }
        elimination_macros[[length(elimination_macros) + 1]] <- list(cmt = cmt, k = k_val)

      } else {
        pk_equations <- c(pk_equations, pl_t)
      }
    }
  }

  # --- Parse EQUATION section ---
  eq_raw <- character(0)
  if (!is.na(idx_eq)) {
    eq_end <- ifelse(!is.na(idx_out), idx_out - 1, n)
    eq_raw <- lines[(idx_eq + 1):eq_end]
  }

  # --- Parse OUTPUT ---
  output_line <- NULL
  if (!is.na(idx_out)) {
    out_lines <- lines[(idx_out + 1):n]
    out_match <- out_lines[stringr::str_detect(out_lines, "output\\s*=")]
    if (length(out_match) > 0) {
      braced <- stringr::str_extract(out_match[1], "\\{[^}]+\\}")
      if (!is.na(braced)) {
        ov <- stringr::str_remove_all(braced, "[{}]") %>% stringr::str_trim()
      } else {
        ov <- stringr::str_match(out_match[1], "output\\s*=\\s*(.+)")[, 2] %>%
          stringr::str_trim()
      }
      output_line <- paste0("# output = {", ov, "}")
    }
  }

  # --- Build rxode2 output ---
  result <- character(0)

  if (length(desc_rx) > 0) result <- c(result, desc_rx, "")

  result <- c(result, "# [INPUT]")
  for (p in model_params) result <- c(result, paste0(p, " = 1;"))
  result <- c(result, "")

  result <- c(result, "# [MODEL]")

  if (!is.null(pkmodel_line)) {
    result <- c(result, .tr_expand_pkmodel(pkmodel_line, pkmodel_output_var,
                                           pk_equations, eq_raw))
  } else if (length(absorption_macros) > 0) {
    result <- c(result, .tr_expand_type2(
      absorption_macros, peripheral_macros, elimination_macros,
      pk_equations, eq_raw, compartments
    ))
  } else {
    for (eq in pk_equations) result <- c(result, .tr_mlx_code_to_rx(eq))
    result <- c(result, .tr_convert_eq_section(eq_raw))

    f_alag <- character(0)
    for (iv in iv_macros) {
      comp_idx <- which(vapply(compartments, function(x) x$cmt, integer(1)) == iv$cmt)
      if (length(comp_idx) == 0) next
      state_var <- compartments[[comp_idx[1]]]$amount

      if (!is.na(iv$p)) {
        cmt_str <- ""
        if (!is.na(iv$comment)) {
          cmt_str <- paste0(" ", .tr_comment_mlx_to_rx(iv$comment))
        }
        f_alag <- c(f_alag, paste0("f(", state_var, ") = ", iv$p, ";", cmt_str))
      }
      if (!is.na(iv$Tlag)) {
        f_alag <- c(f_alag, paste0("alag(", state_var, ") = ", iv$Tlag, ";"))
      }
    }
    if (length(f_alag) > 0) result <- c(result, "", f_alag)
  }

  if (!is.null(output_line)) {
    result <- c(result, "", "# [OUTPUT]", output_line)
  }

  result
}


# =============================================================================
# EQUATION SECTION CONVERTER (MLXTRAN → RXODE2)
# =============================================================================

.tr_convert_eq_section <- function(eq_lines) {
  result <- character(0)
  i <- 1
  n <- length(eq_lines)

  ddt_states <- stringr::str_extract(
    eq_lines[stringr::str_detect(eq_lines, "^\\s*ddt_\\w+")],
    "(?<=ddt_)\\w+"
  )

  while (i <= n) {
    trimmed <- stringr::str_trim(eq_lines[i])

    if (trimmed == "") {
      result <- c(result, "")
      i <- i + 1
      next
    }

    if (stringr::str_detect(trimmed, "^\\s*odeType\\s*=")) {
      i <- i + 1
      next
    }

    if (stringr::str_detect(trimmed, "^if\\s+") &&
        !stringr::str_detect(trimmed, "^if\\s*\\(")) {
      block <- .tr_parse_mlx_if_block(eq_lines, i)
      result <- c(result, block$rx_lines)
      i <- block$next_i
      next
    }

    if (stringr::str_detect(trimmed, "^ddt_\\w+\\s*=")) {
      state <- stringr::str_extract(trimmed, "(?<=^ddt_)\\w+")
      rhs_raw <- stringr::str_match(trimmed, "=\\s*(.*)")[, 2]
      rhs_converted <- .tr_mlx_code_to_rx(paste0("X = ", rhs_raw)) %>%
        stringr::str_remove("^X = ")
      result <- c(result, paste0("d/dt(", state, ") = ", rhs_converted))
      i <- i + 1
      next
    }

    if (stringr::str_detect(trimmed, "^\\w+_0\\s*=")) {
      var_name <- stringr::str_extract(trimmed, "^\\w+(?=_0\\s*=)")
      if (var_name %in% ddt_states) {
        rhs_raw <- stringr::str_match(trimmed, "_0\\s*=\\s*(.*)")[, 2]
        rhs_converted <- .tr_mlx_code_to_rx(paste0("X = ", rhs_raw)) %>%
          stringr::str_remove("^X = ")
        result <- c(result, paste0(var_name, "(0) = ", rhs_converted))
        i <- i + 1
        next
      }
    }

    result <- c(result, .tr_mlx_code_to_rx(trimmed))
    i <- i + 1
  }

  result
}


# =============================================================================
# IF/ELSE BLOCK CONVERTERS
# =============================================================================

.tr_parse_mlx_if_block <- function(lines, start_i) {
  n <- length(lines)
  condition <- stringr::str_trim(lines[start_i]) %>% stringr::str_remove("^if\\s+")

  i <- start_i + 1
  if_body <- character(0)
  else_body <- character(0)
  has_else <- FALSE

  while (i <= n) {
    lt <- stringr::str_trim(lines[i])
    if (lt == "else") { has_else <- TRUE; i <- i + 1; next }
    if (lt == "end")  { i <- i + 1; break }
    if (!has_else) {
      if_body <- c(if_body, lt)
    } else {
      else_body <- c(else_body, lt)
    }
    i <- i + 1
  }

  rx <- character(0)
  if (length(if_body) == 1 && length(else_body) <= 1) {
    rx <- c(rx, paste0("if(", condition, ")"))
    rx <- c(rx, paste0("{ ", .tr_mlx_code_to_rx(if_body[1], add_semicolon = FALSE), " }"))
    if (has_else && length(else_body) == 1) {
      rx <- c(rx, "else")
      rx <- c(rx, paste0("{ ", .tr_mlx_code_to_rx(else_body[1], add_semicolon = FALSE), " };"))
    }
  } else {
    rx <- c(rx, paste0("if(", condition, ") {"))
    for (bl in if_body) rx <- c(rx, paste0("  ", .tr_mlx_code_to_rx(bl, add_semicolon = FALSE)))
    rx <- c(rx, "}")
    if (has_else) {
      rx <- c(rx, "else {")
      for (bl in else_body) rx <- c(rx, paste0("  ", .tr_mlx_code_to_rx(bl, add_semicolon = FALSE)))
      rx <- c(rx, "};")
    }
  }

  list(rx_lines = rx, next_i = i)
}

.tr_convert_inline_if_rx_to_mlx <- function(line) {
  var_name <- stringr::str_extract(line, "^\\w+")
  cond <- stringr::str_match(line, "if\\s*\\((.+?)\\)\\s*\\{")[, 2]
  if_expr <- stringr::str_match(line, "\\{\\s*([^}]+?)\\s*\\}")[, 2] %>%
    stringr::str_trim()
  else_expr <- stringr::str_match(line, "else\\s*\\{\\s*([^}]+?)\\s*\\}")[, 2] %>%
    stringr::str_trim()

  c(
    paste0("if ", cond),
    paste0("  ", var_name, " = ", if_expr),
    "else",
    paste0("  ", var_name, " = ", else_expr),
    "end"
  )
}


# =============================================================================
# PKMODEL MACRO EXPANSION (MLXTRAN → RXODE2)
# =============================================================================

.tr_expand_pkmodel <- function(pkmodel_line, output_var, pk_equations, eq_raw) {
  args_str <- stringr::str_extract(pkmodel_line, "(?<=pkmodel\\s{0,3}\\().*(?=\\))") %>%
    stringr::str_trim()
  all_args <- stringr::str_split(args_str, ",")[[1]] %>% stringr::str_trim()

  named_args <- all_args[stringr::str_detect(all_args, "=")]
  positional <- all_args[!stringr::str_detect(all_args, "=")]

  p_val <- NULL
  for (na in named_args) {
    if (stringr::str_detect(na, "^p\\s*=")) {
      p_val <- stringr::str_remove(na, "^p\\s*=\\s*") %>% stringr::str_trim()
    }
  }

  n_pos <- length(positional)
  has_tlag <- n_pos %in% c(4, 6)
  is_2c <- n_pos >= 5

  tlag_var <- NULL
  if (is_2c) {
    if (has_tlag) {
      tlag_var <- positional[1]; ka <- positional[2]; v <- positional[3]
      cl <- positional[4]; k12 <- positional[5]; k21 <- positional[6]
    } else {
      ka <- positional[1]; v <- positional[2]; cl <- positional[3]
      k12 <- positional[4]; k21 <- positional[5]
    }
  } else {
    if (has_tlag) {
      tlag_var <- positional[1]; ka <- positional[2]; v <- positional[3]; cl <- positional[4]
    } else {
      ka <- positional[1]; v <- positional[2]; cl <- positional[3]
    }
  }

  result <- character(0)

  for (eq in pk_equations) result <- c(result, .tr_mlx_code_to_rx(eq))

  eq_converted <- .tr_convert_eq_section(eq_raw)
  if (length(eq_converted) > 0) {
    eq_converted <- eq_converted[stringr::str_trim(eq_converted) != ""]
    if (length(eq_converted) > 0) result <- c(result, eq_converted)
  }

  result <- c(result, paste0(output_var, " = Ac/", v, ";"))

  result <- c(result, "", "### Initial conditions")
  result <- c(result, "Ad(0) = 0;")
  result <- c(result, "Ac(0) = 0;")
  if (is_2c) result <- c(result, "Ap(0) = 0;")

  result <- c(result, "", "### Differential equations")
  result <- c(result, paste0("d/dt(Ad) = -", ka, "*Ad;"))
  if (is_2c) {
    result <- c(result, paste0("d/dt(Ac) = ", ka, "*Ad - ", cl, "*Ac/", v,
                               " - ", k12, "*Ac + ", k21, "*Ap;"))
    result <- c(result, paste0("d/dt(Ap) = ", k12, "*Ac - ", k21, "*Ap;"))
  } else {
    result <- c(result, paste0("d/dt(Ac) = ", ka, "*Ad - ", cl, "*Ac/", v, ";"))
  }

  if (!is.null(p_val)) result <- c(result, "", paste0("f(Ad) = ", p_val, ";"))
  if (!is.null(tlag_var)) result <- c(result, paste0("alag(Ad) = ", tlag_var, ";"))

  result
}


# =============================================================================
# MACROS TYPE 2 EXPANSION (MLXTRAN → RXODE2)
# =============================================================================

.tr_expand_type2 <- function(absorption_macros, peripheral_macros, elimination_macros,
                             pk_equations, eq_raw, compartments) {
  result <- character(0)

  for (eq in pk_equations) result <- c(result, .tr_mlx_code_to_rx(eq))

  if (length(compartments) > 0) {
    comp1 <- compartments[[1]]
    if (!is.na(comp1$concentration) && !is.na(comp1$volume)) {
      result <- c(result, paste0(comp1$concentration, " = Ac/", comp1$volume, ";"))
    }
  }

  eq_converted <- .tr_convert_eq_section(eq_raw)
  if (length(eq_converted) > 0) result <- c(result, eq_converted)

  depot_names <- character(0)
  central_inflow <- character(0)
  f_lines_out <- character(0)
  alag_lines_out <- character(0)
  zero_order_defs <- character(0)
  depot_odes <- character(0)

  for (j in seq_along(absorption_macros)) {
    am <- absorption_macros[[j]]

    if (!is.null(am$Tk0)) {
      dn <- paste0("Ad", j - 1)
      depot_names <- c(depot_names, dn)
      rate_var <- paste0("ka", j - 1)
      dose_expr <- if (!is.null(am$p)) paste0("amtDose*", am$p) else "amtDose"
      zero_order_defs <- c(zero_order_defs,
        paste0(rate_var, " = if(", dn, " > 0){(", dose_expr, ")/", am$Tk0, "} else {0};"))
      central_inflow <- c(central_inflow, rate_var)
      depot_odes <- c(depot_odes, paste0("d/dt(", dn, ") = -", rate_var, ";"))
      if (!is.null(am$p)) f_lines_out <- c(f_lines_out, paste0("f(", dn, ") = ", am$p, ";"))
      if (!is.null(am$Tlag)) alag_lines_out <- c(alag_lines_out, paste0("alag(", dn, ") = ", am$Tlag, ";"))

    } else if (!is.null(am$ka)) {
      dn <- "Ad"
      depot_names <- c(depot_names, dn)
      central_inflow <- c(central_inflow, paste0(am$ka, "*", dn))
      depot_odes <- c(depot_odes, paste0("d/dt(", dn, ") = -", am$ka, "*", dn, ";"))
      if (!is.null(am$p)) f_lines_out <- c(f_lines_out, paste0("f(", dn, ") = ", am$p, ";"))
      if (!is.null(am$Tlag)) alag_lines_out <- c(alag_lines_out, paste0("alag(", dn, ") = ", am$Tlag, ";"))
    }
  }

  elim_term <- ""
  if (length(elimination_macros) > 0) {
    elim_term <- paste0(" - ", elimination_macros[[1]]$k, "*Ac")
  }

  periph_central <- ""
  periph_ode <- character(0)
  if (length(peripheral_macros) > 0) {
    pm <- peripheral_macros[[1]]
    periph_central <- paste0(" - ", pm$k12, "*Ac + ", pm$k21, "*Ap")
    periph_ode <- paste0("d/dt(Ap) = ", pm$k12, "*Ac - ", pm$k21, "*Ap;")
  }

  if (length(zero_order_defs) > 0) result <- c(result, zero_order_defs)

  result <- c(result, "", "### Initial conditions")
  for (dn in depot_names) result <- c(result, paste0(dn, "(0) = 0;"))
  result <- c(result, "Ac(0) = 0;")
  if (length(peripheral_macros) > 0) result <- c(result, "Ap(0) = 0;")

  result <- c(result, "", "### Differential equations")
  result <- c(result, depot_odes)

  inflow_str <- paste(central_inflow, collapse = " + ")
  result <- c(result, paste0("d/dt(Ac) = ", inflow_str, elim_term, periph_central, ";"))
  if (length(periph_ode) > 0) result <- c(result, periph_ode)

  if (length(f_lines_out) > 0 || length(alag_lines_out) > 0) {
    result <- c(result, "", f_lines_out, alag_lines_out)
  }

  result
}


# =============================================================================
# RXODE2 → MLXTRAN CONVERSION
# =============================================================================

.tr_rxode_to_mlxtran <- function(lines, dm_list, regressors, output_vars, macros, stiff) {
  n <- length(lines)

  idx_input  <- which(stringr::str_detect(lines, "^\\s*#?\\s*\\[INPUT\\]"))[1]
  idx_model  <- which(stringr::str_detect(lines, "^\\s*#?\\s*\\[MODEL\\]"))[1]
  idx_output <- which(stringr::str_detect(lines, "^\\s*#?\\s*\\[OUTPUT\\]"))[1]

  # --- Description ---
  desc_mlx <- character(0)
  if (!is.na(idx_input) && idx_input > 1) {
    desc_raw <- lines[1:(idx_input - 1)]
    desc_raw <- desc_raw[stringr::str_trim(desc_raw) != ""]
    if (length(desc_raw) > 0) {
      desc_mlx <- vapply(desc_raw, .tr_comment_rx_to_mlx, character(1), USE.NAMES = FALSE)
    }
  }

  # --- Input parameters ---
  model_params <- character(0)
  if (!is.na(idx_input) && !is.na(idx_model)) {
    input_lines <- lines[(idx_input + 1):(idx_model - 1)]
    for (il in input_lines) {
      il_t <- stringr::str_trim(il)
      if (il_t == "" || stringr::str_detect(il_t, "^#")) next
      pname <- stringr::str_extract(il_t, "^\\w+")
      if (!is.na(pname)) model_params <- c(model_params, pname)
    }
  }

  # --- Model section ---
  model_end <- ifelse(!is.na(idx_output), idx_output - 1, n)
  model_section <- lines[(idx_model + 1):model_end]

  # First pass: extract structure
  ode_states <- character(0)
  f_map <- list()
  alag_map <- list()
  f_comments <- list()

  for (line in model_section) {
    trimmed <- stringr::str_trim(line)

    if (stringr::str_detect(trimmed, "^d/dt\\(\\w+\\)\\s*=")) {
      state <- stringr::str_extract(trimmed, "(?<=d/dt\\()\\w+")
      ode_states <- c(ode_states, state)
    }

    if (stringr::str_detect(trimmed, "^f\\(\\w+\\)\\s*=")) {
      state <- stringr::str_extract(trimmed, "(?<=f\\()\\w+")
      parsed <- .tr_extract_rhs(trimmed)
      f_map[[state]] <- parsed$rhs
      f_comments[[state]] <- parsed$comment
    }

    if (stringr::str_detect(trimmed, "^alag\\(\\w+\\)\\s*=")) {
      state <- stringr::str_extract(trimmed, "(?<=alag\\()\\w+")
      parsed <- .tr_extract_rhs(trimmed)
      alag_map[[state]] <- parsed$rhs
    }
  }

  # Check pkmodel eligibility
  use_pkmodel <- macros && length(ode_states) %in% c(2, 3) &&
    .tr_can_use_pkmodel(ode_states, model_section)

  # Second pass: convert model body lines for EQUATION section
  eq_body <- character(0)
  i <- 1
  n_model <- length(model_section)

  while (i <= n_model) {
    trimmed <- stringr::str_trim(model_section[i])

    if (trimmed == "") { eq_body <- c(eq_body, ""); i <- i + 1; next }

    # Skip f() and alag() (handled via iv macros)
    if (stringr::str_detect(trimmed, "^f\\(\\w+\\)\\s*=") ||
        stringr::str_detect(trimmed, "^alag\\(\\w+\\)\\s*=")) {
      i <- i + 1
      next
    }

    # d/dt(X) = expr; → ddt_X = expr
    if (stringr::str_detect(trimmed, "^d/dt\\(\\w+\\)\\s*=")) {
      state <- stringr::str_extract(trimmed, "(?<=d/dt\\()\\w+")
      rhs_raw <- stringr::str_match(trimmed, "=\\s*(.*)")[, 2]
      rhs_mlx <- .tr_rx_code_to_mlx(paste0("X = ", rhs_raw)) %>%
        stringr::str_remove("^X = ")
      eq_body <- c(eq_body, paste0("ddt_", state, " = ", rhs_mlx))
      i <- i + 1
      next
    }

    # X(0) = expr; → X_0 = expr
    if (stringr::str_detect(trimmed, "^\\w+\\(0\\)\\s*=")) {
      state <- stringr::str_extract(trimmed, "^\\w+(?=\\(0\\))")
      rhs_raw <- stringr::str_match(trimmed, "=\\s*(.*)")[, 2]
      rhs_mlx <- .tr_rx_code_to_mlx(paste0("X = ", rhs_raw)) %>%
        stringr::str_remove("^X = ")
      eq_body <- c(eq_body, paste0(state, "_0 = ", rhs_mlx))
      i <- i + 1
      next
    }

    # Inline if: var = if(cond){expr} else {expr};
    if (stringr::str_detect(trimmed, "^\\w+\\s*=\\s*if\\s*\\(")) {
      eq_body <- c(eq_body, .tr_convert_inline_if_rx_to_mlx(trimmed))
      i <- i + 1
      next
    }

    # Multi-line if(condition) block
    if (stringr::str_detect(trimmed, "^if\\s*\\(") &&
        !stringr::str_detect(trimmed, "^\\w+\\s*=\\s*if")) {
      block <- .tr_parse_rx_if_block(model_section, i)
      eq_body <- c(eq_body, block$mlx_lines)
      i <- block$next_i
      next
    }

    # Regular line
    eq_body <- c(eq_body, .tr_rx_code_to_mlx(model_section[i]))
    i <- i + 1
  }

  # output_vars is passed from the caller and validated there

  # --- Build MLXTRAN output ---
  result <- character(0)

  if (length(desc_mlx) > 0) result <- c(result, "DESCRIPTION:", desc_mlx, "")

  result <- c(result, "[LONGITUDINAL]")
  all_input <- model_params
  if (!is.null(regressors)) all_input <- c(all_input, regressors)
  result <- c(result, paste0("input = {", paste(all_input, collapse = ", "), "}"))
  if (!is.null(regressors)) {
    for (r in regressors) result <- c(result, paste0(r, " = {use = regressor}"))
  }
  result <- c(result, "")

  if (use_pkmodel) {
    result <- c(result, "PK:")
    result <- c(result, .tr_build_pkmodel(ode_states, model_section, f_map,
                                          alag_map, eq_body))
  } else {
    result <- c(result, "PK:")
    for (j in seq_along(ode_states)) {
      result <- c(result, paste0("compartment(cmt = ", j, ", amount = ",
                                 ode_states[j], ")"))
    }

    if (!is.null(dm_list)) {
      result <- c(result, "")
      for (j in seq_along(dm_list$cmt)) {
        cmt_num <- dm_list$cmt[j]
        adm_num <- dm_list$adm[j]
        state <- ode_states[cmt_num]

        iv_str <- paste0("iv(cmt = ", cmt_num, ", adm = ", adm_num)
        if (!is.null(f_map[[state]])) iv_str <- paste0(iv_str, ", p = ", f_map[[state]])
        if (!is.null(alag_map[[state]])) iv_str <- paste0(iv_str, ", Tlag = ",
                                                          alag_map[[state]])
        iv_str <- paste0(iv_str, ")")

        if (!is.null(f_comments[[state]]) && f_comments[[state]] != "") {
          iv_str <- paste0(iv_str, "\t\t", .tr_comment_rx_to_mlx(f_comments[[state]]))
        }
        result <- c(result, iv_str)
      }
    }

    result <- c(result, "", "")
    result <- c(result, "EQUATION:")
    if (stiff) result <- c(result, "odeType = stiff")
    result <- c(result, eq_body)
  }

  result <- c(result, "", "OUTPUT:")
  result <- c(result, paste0("output = {", paste(output_vars, collapse = ", "), "}"))

  result
}


# =============================================================================
# MULTI-LINE IF/ELSE PARSER (RXODE2 → MLXTRAN)
# =============================================================================

.tr_parse_rx_if_block <- function(lines, start_i) {
  n <- length(lines)
  trimmed <- stringr::str_trim(lines[start_i])
  condition <- stringr::str_match(trimmed, "^if\\s*\\((.+?)\\)")[, 2]

  i <- start_i + 1
  if_body <- character(0)
  else_body <- character(0)
  has_else <- FALSE

  while (i <= n) {
    lt <- stringr::str_trim(lines[i])

    if (lt == "else" || stringr::str_detect(lt, "^\\}\\s*else")) {
      has_else <- TRUE
      i <- i + 1
      next
    }

    content <- lt %>%
      stringr::str_remove("^\\{\\s*") %>%
      stringr::str_remove("\\s*\\}\\s*;?$") %>%
      stringr::str_trim()

    has_closing_brace <- stringr::str_detect(lt, "\\}\\s*;?\\s*$")

    if (content != "" && content != "{" && content != "}") {
      if (!has_else) {
        if_body <- c(if_body, content)
      } else {
        else_body <- c(else_body, content)
      }
    }

    i <- i + 1
    if (has_closing_brace && has_else) break
    if (has_closing_brace && !has_else) {
      if (i <= n && stringr::str_detect(stringr::str_trim(lines[i]), "^else")) next
      break
    }
  }

  mlx <- c(paste0("if ", condition))
  for (bl in if_body) mlx <- c(mlx, paste0("  ", .tr_rx_code_to_mlx(bl)))
  if (has_else) {
    mlx <- c(mlx, "else")
    for (bl in else_body) mlx <- c(mlx, paste0("  ", .tr_rx_code_to_mlx(bl)))
  }
  mlx <- c(mlx, "end")

  list(mlx_lines = mlx, next_i = i)
}


# =============================================================================
# PKMODEL ELIGIBILITY CHECK (RXODE2 → MLXTRAN MACROS)
# =============================================================================

.tr_can_use_pkmodel <- function(ode_states, model_lines) {
  n_states <- length(ode_states)
  if (!(n_states %in% c(2, 3))) return(FALSE)

  ddt_rhs <- character(0)
  for (line in model_lines) {
    trimmed <- stringr::str_trim(line)
    if (stringr::str_detect(trimmed, "^d/dt\\(")) {
      rhs <- stringr::str_match(trimmed, "=\\s*(.*)")[, 2] %>%
        stringr::str_remove("#.*") %>% stringr::str_remove(";") %>%
        stringr::str_trim()
      ddt_rhs <- c(ddt_rhs, rhs)
    }
  }
  if (length(ddt_rhs) != n_states) return(FALSE)

  depot <- ode_states[1]
  if (!stringr::str_detect(ddt_rhs[1], paste0("^-\\s*\\w+\\s*\\*\\s*", depot, "$"))) return(FALSE)

  central <- ode_states[2]
  if (!stringr::str_detect(ddt_rhs[2], paste0("\\w+\\s*\\*\\s*", depot))) return(FALSE)
  if (!stringr::str_detect(ddt_rhs[2], paste0("\\w+\\s*\\*\\s*", central, "\\s*/\\s*\\w+"))) return(FALSE)

  if (n_states == 3) {
    periph <- ode_states[3]
    if (!stringr::str_detect(ddt_rhs[3], paste0("\\w+\\s*\\*\\s*", central))) return(FALSE)
    if (!stringr::str_detect(ddt_rhs[3], paste0("\\w+\\s*\\*\\s*", periph))) return(FALSE)
  }

  TRUE
}


# =============================================================================
# PKMODEL BUILDER (RXODE2 → MLXTRAN MACROS)
# =============================================================================

.tr_build_pkmodel <- function(ode_states, model_section, f_map, alag_map, eq_body) {
  depot <- ode_states[1]
  central <- ode_states[2]

  ka <- NULL; v <- NULL; cl <- NULL; k12 <- NULL; k21 <- NULL

  for (line in model_section) {
    trimmed <- stringr::str_trim(line)
    if (stringr::str_detect(trimmed, paste0("^d/dt\\(", depot, "\\)"))) {
      rhs <- stringr::str_match(trimmed, "=\\s*(.*)")[, 2] %>%
        stringr::str_remove("#.*") %>% stringr::str_remove(";") %>% stringr::str_trim()
      ka <- stringr::str_extract(rhs, "(?<=^-\\s{0,5})\\w+")
      if (is.na(ka)) ka <- stringr::str_extract(rhs, "(?<=^-)\\w+")
    }
    if (stringr::str_detect(trimmed, paste0("^d/dt\\(", central, "\\)"))) {
      rhs <- stringr::str_match(trimmed, "=\\s*(.*)")[, 2] %>%
        stringr::str_remove("#.*") %>% stringr::str_remove(";") %>% stringr::str_trim()
      cl_v_match <- stringr::str_match(rhs, "(\\w+)\\s*\\*\\s*\\w+\\s*/\\s*(\\w+)")
      cl <- cl_v_match[, 2]
      v <- cl_v_match[, 3]
    }
    if (length(ode_states) == 3) {
      periph <- ode_states[3]
      if (stringr::str_detect(trimmed, paste0("^d/dt\\(", periph, "\\)"))) {
        rhs <- stringr::str_match(trimmed, "=\\s*(.*)")[, 2] %>%
          stringr::str_remove("#.*") %>% stringr::str_remove(";") %>% stringr::str_trim()
        pm <- stringr::str_match(rhs, "(\\w+)\\s*\\*\\s*\\w+\\s*-\\s*(\\w+)")
        k12 <- pm[, 2]
        k21 <- pm[, 3]
      }
    }
  }

  tlag_val <- alag_map[[depot]]
  p_val <- f_map[[depot]]

  conc_var <- "Cc"
  if (!is.null(v)) {
    for (line in model_section) {
      trimmed <- stringr::str_trim(line)
      if (stringr::str_detect(trimmed, paste0("^\\w+\\s*=\\s*", central, "\\s*/\\s*", v))) {
        conc_var <- stringr::str_extract(trimmed, "^\\w+")
        break
      }
    }
  }

  result <- character(0)

  # Include non-special equation body lines, skip ODE/init headers and content
  skip_next_empty <- FALSE
  for (idx in seq_along(eq_body)) {
    el_t <- stringr::str_trim(eq_body[idx])

    if (stringr::str_detect(el_t, "^ddt_")) { skip_next_empty <- TRUE; next }
    if (stringr::str_detect(el_t, "^\\w+_0\\s*=")) { skip_next_empty <- TRUE; next }
    if (!is.null(v) && stringr::str_detect(el_t, paste0("^", conc_var, "\\s*=\\s*", central, "\\s*/\\s*", v))) next

    if (stringr::str_detect(el_t, "^;;;\\s*(Initial conditions|Differential equations|ODEs)")) next

    if (el_t == "" && skip_next_empty) { skip_next_empty <- FALSE; next }
    skip_next_empty <- FALSE

    result <- c(result, eq_body[idx])
  }

  # Build pkmodel() call
  pkm_args <- character(0)
  if (!is.null(tlag_val)) pkm_args <- c(pkm_args, tlag_val)
  pkm_args <- c(pkm_args, ka, v, cl)
  if (!is.null(k12)) pkm_args <- c(pkm_args, k12, k21)

  call_str <- paste(pkm_args, collapse = ", ")
  if (!is.null(p_val)) call_str <- paste0(call_str, ", p = ", p_val)

  result <- c(result, paste0(conc_var, " = pkmodel(", call_str, ")"))

  result
}

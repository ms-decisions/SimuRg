gmo_converter <- function(gco_path, output_path=NULL) {
  gco <- read_smrg_obj(gco_path)
  theta <- gco$theta
  covs <- gco$covs
  re <- gco$re
  ruv <- gco$ruv
  model_file <- gco$model
  data_path <- gco$data
  data <- utils::read.csv(data_path, check.names = FALSE)
  model_lines <- readLines(model_file)
  input_idx <- which(grepl("^# \\[INPUT\\]", model_lines))[1]
  model_idx <- which(grepl("^# \\[MODEL\\]", model_lines))[1]
  pre_input <- model_lines[seq_len(input_idx)]
  post_input <- model_lines[seq(model_idx, length(model_lines))]
  input_lines <- character()
  if (!is.null(theta) && nrow(theta) > 0) {
    for (i in seq_len(nrow(theta))) {
      name <- theta$NAME[i]
      trans <- tolower(theta$TRANS[i])
      init <- theta$INIT[i]
      if (is.na(init)) init <- 0
      if (trans == "lognormal") {
        input_lines <- c(input_lines, paste0(name, "_pop = log(", init, ");"))
      } else if (trans == "logitnormal") {
        input_lines <- c(input_lines, paste0(name, "_pop = logit(", init, ");"))
      } else {
        input_lines <- c(input_lines, paste0(name, "_pop = ", init, ";"))
      }
    }
    input_lines <- c(input_lines, "")
  }
  cov_stat_lines <- character()
  cov_par_lines <- character()
  if (!is.null(covs) && length(covs) > 0) {
    for (j in seq_along(covs)) {
      covj <- covs[[j]]
      par <- covj$PAR
      covname <- covj$COVNAME
      func <- covj$FUNC
      cov_trans <- covj$TRANS
      init <- if (is.null(covj$INIT)) 1 else covj$INIT
      if (!is.null(func) && func == "linear" && !is.null(cov_trans) && tolower(cov_trans) == "median") {
        med_name <- paste0(covname, "_med")
        med_val <- stats::median(data[[covname]], na.rm = TRUE)
        cov_stat_lines <- c(cov_stat_lines, paste0(med_name, " = ", med_val, ";"))
        beta_name <- paste0("beta_", par, "_", covname)
        cov_par_lines <- c(cov_par_lines, paste0(beta_name, " = ", init, ";"))
      } else if (!is.null(covj$REF)) {
        # Categorical covariate: create one beta per observed non-reference category value
        ref <- covj$REF
        cov_values <- unique(stats::na.omit(data[[covname]]))
        cov_values <- sort(cov_values)
        non_ref_vals <- cov_values[cov_values != ref]
        for (val in non_ref_vals) {
          beta_name <- paste0("beta_", par, "_", val)
          cov_par_lines <- c(cov_par_lines, paste0(beta_name, " = ", init, ";"))
        }
      }
    }
  }
  if (length(cov_stat_lines) > 0) {
    input_lines <- c(input_lines, cov_stat_lines, "")
  }
  if (length(cov_par_lines) > 0) {
    input_lines <- c(input_lines, cov_par_lines, "")
  }
  omega_names <- character()
  if (!is.null(re) && !is.null(re$est)) {
    est_mat <- re$est
    if (is.matrix(est_mat)) {
      pars <- colnames(est_mat)
      for (k in seq_along(pars)) {
        par_name <- pars[k]
        est_col <- est_mat[, k]
        if (any(!is.na(est_col) & est_col)) {
          oname <- paste0("omega_", par_name)
          omega_names <- c(omega_names, par_name)
          input_lines <- c(input_lines, paste0(oname, " = 1;"))
        }
      }
    }
  }
  if (length(omega_names) > 0) {
    input_lines <- c(input_lines, "")
  }
  if (!is.null(ruv) && !is.null(ruv$ERR)) {
    ruv_list <- list(ruv)
  } else if (is.list(ruv) && length(ruv) > 0 && !is.null(ruv[[1]]$ERR)) {
    ruv_list <- ruv
  } else {
    ruv_list <- list()
  }
  err_par_lines <- character()
  for (r in seq_along(ruv_list)) {
    ruvr <- ruv_list[[r]]
    pred <- ruvr$PRED
    err_type <- ruvr$ERR
    init <- ruvr$INIT
    if (is.null(err_type) || is.null(pred) || is.null(init)) next
    if (tolower(err_type) == "combined1" && length(init) >= 2) {
      a_name <- paste0(pred, "_a")
      b_name <- paste0(pred, "_b")
      err_par_lines <- c(err_par_lines, paste0(a_name, " = ", init[1], ";"))
      err_par_lines <- c(err_par_lines, paste0(b_name, " = ", init[2], ";"))
    } else if (tolower(err_type) == "proportional" && length(init) >= 1) {
      b_name <- paste0(pred, "_b")
      err_par_lines <- c(err_par_lines, paste0(b_name, " = ", init[1], ";"))
    }
  }
  if (length(err_par_lines) > 0) {
    input_lines <- c(input_lines, err_par_lines, "")
  }
  tv_lines <- character()
  if (!is.null(theta) && nrow(theta) > 0) {
    for (i in seq_len(nrow(theta))) {
      name <- theta$NAME[i]
      trans <- tolower(theta$TRANS[i])
      if (trans == "lognormal") {
        tv_lines <- c(tv_lines, paste0(name, "_tv = exp(", name, "_pop);"))
      } else if (trans == "logitnormal") {
        tv_lines <- c(tv_lines, paste0(name, "_tv = expit(", name, "_pop);"))
      } else {
        tv_lines <- c(tv_lines, paste0(name, "_tv = ", name, "_pop;"))
      }
    }
  }
  if (length(tv_lines) > 0) {
    input_lines <- c(input_lines, tv_lines, "")
  }
  par_def_lines <- character()
  if (!is.null(theta) && nrow(theta) > 0) {
    for (i in seq_len(nrow(theta))) {
      name <- theta$NAME[i]
      par_trans <- tolower(theta$TRANS[i])
      has_omega <- name %in% omega_names
      base_expr <- paste0(name, "_tv")
      if (has_omega) {
        oname <- paste0("omega_", name)
        if (par_trans == "lognormal") {
          base_expr <- paste0(base_expr, " * exp(", oname, ")")
        } else if (par_trans == "logitnormal") {
          base_expr <- paste0(base_expr, " * expit(", oname, ")")
        } else {
          base_expr <- paste0(base_expr, " + ", oname)
        }
      }
      if (!is.null(covs) && length(covs) > 0) {
        for (j in seq_along(covs)) {
          covj <- covs[[j]]
          if (!identical(covj$PAR, name)) next
          covname <- covj$COVNAME
          func <- covj$FUNC
          cov_trans <- covj$TRANS
          if (!is.null(func) && func == "linear" && !is.null(cov_trans) && tolower(cov_trans) == "median") {
            med_name <- paste0(covname, "_med")
            beta_name <- paste0("beta_", name, "_", covname)
            base_expr <- paste0(base_expr, " * (", covname, "/", med_name, ")^", beta_name)
          } else if (!is.null(covj$REF)) {
            # Categorical covariate
            ref <- covj$REF
            cov_values <- unique(stats::na.omit(data[[covname]]))
            cov_values <- sort(cov_values)
            non_ref_vals <- cov_values[cov_values != ref]
            if (length(non_ref_vals) == 0) next
            cat_terms <- character()
            for (val in non_ref_vals) {
              beta_name <- paste0("beta_", name, "_", val)
              cat_terms <- c(cat_terms, paste0(beta_name, " * (", covname, " == ", val, ")"))
            }
            if (length(cat_terms) > 0) {
              if (par_trans == "lognormal") {
                # Multiplicative effect on typical value: tv * exp(sum(beta * indicator))
                base_expr <- paste0(base_expr, " * exp(", paste(cat_terms, collapse = " + "), ")")
              } else {
                # Additive effect on natural scale for non-lognormal
                base_expr <- paste0(base_expr, " + ", paste(cat_terms, collapse = " + "))
              }
            }
          }
        }
      }
      par_def_lines <- c(par_def_lines, paste0(name, " = ", base_expr))
    }
  }
  all_input <- c(pre_input, input_lines, "", par_def_lines, "", post_input)
  res_err_lines <- character()
  for (r in seq_along(ruv_list)) {
    ruvr <- ruv_list[[r]]
    pred <- ruvr$PRED
    err_type <- ruvr$ERR
    if (is.null(err_type) || is.null(pred)) next
    if (tolower(err_type) == "combined1") {
      res_err_lines <- c(res_err_lines, paste0(pred, "_Res_err = ", pred, " * (1 + ", pred, "_b ) + ", pred, "_a;"))
    } else if (tolower(err_type) == "proportional") {
      res_err_lines <- c(res_err_lines, paste0(pred, "_Res_err = ", pred, " * (1 + ", pred, "_b );"))
    }
  }
  if (length(res_err_lines) > 0) {
    all_input <- c(all_input, res_err_lines)
  }
  writeLines(all_input, con = output_path)
}


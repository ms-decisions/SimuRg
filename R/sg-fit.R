## Author: Mikhailova Anna, Kulesh Victoria
## First created: 2025-11-05
## Description: function for fit
## Keywords: SimuRg, fit

#' Run fit with monolix/simurg/nonmem fitter
#'
#' @inheritParams sg_dummy
#' @param model string. Path to a txt file file with model structure in Monolix syntax. This should be a valid file path pointing to a model file that defines the pharmacokinetic/pharmacodynamic model structure
#' @param theta tibble or data.frame. Should contain the following columns:
#'   * `NAME` - string, parameter name. Should match names specified in the model file
#'   * `TRANS` - string, parameter transformation type. Can be: `"normal"`, `"logNormal"`, `"logitNormal"`
#'   * `INIT` - numeric, initial value for the parameter
#'   * `LB` - numeric, lower bound for logit transformation, `NA` for other transformations
#'   * `UB` - numeric, upper bound for logit transformation,  `NA` for other transformation
#'   * `EST` - logical, whether to estimate this parameter.
#' @param max_wait_time numeric. Maximum time in seconds to wait for fit results to complete. Default is 3600 seconds (1 hour). Set to `Inf` for no timeout.
#' @returns if option `fit = TRUE`, generalized simurg output object is returned. Otherwise, the file for fit is written and no output is returned
#' @examples
#' \donttest{
#' library(tibble)
#' library(dplyr)
#' library(stringr)
#' # Specify structural model path. For fit, should be in Monolix syntax
#' model <- system.file("extdata", "models", "model_PK_1c.txt", package = "SimuRg")
#' # Specify data path. Should be ADPPK-like format
#' data  <- system.file("extdata", "datasets", "dspk-warf.csv", package = "SimuRg")
#'
#' # Specify headers. Should be a list of lists with the following elements:
#' # name - string, name of the column
#' # use - string, use of the column from Monolix documentation
#' # type - string, type of the column. for use = "covariate", should be
#' # "continuous" or "categorical", for other uses should be NULL
#'  headers <- list(list(name = "ID", use = "identifier", type = NULL),
#'                  list(name = "TIME", use = "time", type = NULL),
#'                  list(name = "DV", use = "observation", type = "continuous"),
#'                  list(name = "DVID", use = "observationtype", type = NULL),
#'                  list(name = "ADM", use = "administration", type = NULL),
#'                  list(name = "AMT", use = "amount", type = NULL),
#'                  list(name = "EVID", use = "eventidentifier", type = NULL),
#'                  list(name = "MDV", use = "missingdependentvariable", type = NULL),
#'                  list(name = "AGE", use = "covariate", type = "continuous"),
#'                  list(name = "AGE_centered", use = "covariate", type = "continuous"),
#'                  list(name = "SEX", use = "covariate", type = "categorical"),
#'                  list(name = "WEIGHT", use = "covariate", type = "continuous"),
#'                  list(name = "BMI", use = "covariate", type = "continuous"))
#'
#' # Dataset with the parameters properties. Should be a tibble with the following columns:
#' # NAME - string, name of the parameter
#' # TRANS - string, distribution of the parameter. Should be one of the following:
#' # "normal", "logNormal", "logitNormal"
#' # INIT - numeric, initial value of the parameter or its fixed value
#' # LB - numeric, lower bound of the parameter for logit transformation
#' # UB - numeric, upper bound of the parameter for logit transformation
#' # EST - logical, estimation status of the parameter. For estimation, should
#' # be TRUE, for fixed, should be FALSE
#'  theta <- tribble(~NAME, ~TRANS, ~INIT, ~LB, ~UB, ~EST,
#'                   "Cl", "logNormal", 0.2, NA, NA, TRUE,
#'                   "V", "logNormal", 20, NA, NA, TRUE,
#'                   "ka", "logNormal", 0.2, NA, NA, TRUE
#'  )
#'  # Examples of random effect model specification
#'  # Examples of random effect model specification for fit
#'  # Single observation (legacy format):
#'  # YNAME - string, name of the observation, ususally y1, y2, ...
#'  # DVID - numeric, observation type identifier corresponding to DVID column values
#'  # TRANS - string, residual error distribution. Can be: "normal", "logNormal",
#'  # "logitNormal"
#'  # PRED - string, prediction variable name from the model
#'  # ERR - string, error model type. Options include: "constant" for additive
#'  # error, "proportional" for proportional error, "combined1" for combined
#'  # additive and proportional error
#'  # INIT - numeric vector, initial values for error parameters
#'  # (length depends on error model)
#'  # EST - logical vector, whether to estimate each error parameter
#'  # (same length as INIT)
#'  # BLQM - below limit of quantification method (can be NULL)
#'  # Single observation (legacy format):
#'  ruv <- list(YNAME = "y1", DVID = 1, TRANS = "normal", PRED = "Cc",
#'              ERR = "combined1", INIT = c(1, 1), EST = c(TRUE, TRUE), BLQM = NULL)
#'
#'  # Multiple observations (recommended format):
#' ruv <- list(
#'   list(YNAME = "y1", DVID = 1, TRANS = "normal", PRED = "Cc",
#'        ERR = "combined1", INIT = c(1, 1), EST = c(TRUE, TRUE), BLQM = NULL),
#'   list(YNAME = "y2", DVID = 2, TRANS = "normal", PRED = "EFF",
#'       ERR = "proportional", INIT = c(0.1), EST = c(TRUE), BLQM = NULL)
#' )
#'
#'  # Example of random effects (RE) specification.
#'  #
#'  # init matrix: provides the initial values for the random effects.
#'  # est matrix: controls how each random effect is handled:
#'  #   TRUE  - the random effect will be estimated,
#'  #   FALSE - the random effect will be fixed at its initial value,
#'  #   NA    - no random effect will be applied.
#'  #
#'  # To fit a model without random effects, set all entries in the
#'  # est matrix to NA.
#'  #
#'  # The same logic applies to the between-occasion variability
#'  # matrix (occ).
#'
#'  re <- list(init = tribble(~Cl, ~V, ~ka,
#'                            1, 0, 0,
#'                            0, 0, 0,
#'                            0, 0, 1) %>% as.matrix(),
#'             est = tribble(~Cl, ~V, ~ka,
#'                           TRUE, NA, NA,
#'                           NA, NA, NA,
#'                           NA, NA, TRUE) %>% as.matrix())
#' # Example of between-occasion variability (BOV) specification. The structure
#' # is the same as for RE, but for BOV
#'
#'  occ <- list(init = tribble(~Cl, ~V, ~ka,
#'                             0, 0, 0,
#'                             0, 0, 0,
#'                             0, 0, 0) %>% as.matrix(),
#'              est = tribble(~Cl, ~V, ~ka,
#'                            NA, NA, NA,
#'                            NA, NA, NA,
#'                            NA, NA, NA) %>% as.matrix())
#' # Example of covariate specification. Should be a list of lists with the following elements:
#' # PAR - string, name of the parameter to which the covariate is applied
#' # COVNAME - string, name of the covariate
#' # FUNC - string, function to apply to the covariate. Should be "linear" for
#' # continuous covariates, "categorical" for categorical covariates
#' # TRANS - string, transformation of the covariate. Should be "median" for
#' # continuous covariates, "reference" for categorical covariates
#' # INIT - numeric, initial value of the covariate
#' # EST - logical, estimation status of the covariate. For estimation, should
#' #be TRUE, for fixed, should be FALSE
#'  covs <- list(list(PAR = "V", COVNAME = "AGE", FUNC = "linear",
#'                    TRANS = "median", INIT = 1, EST = TRUE),
#'               list(PAR = "ka", COVNAME = "SEX", REF = 0, INIT = 1, EST = TRUE))
#'  output_path <- str_c(tempdir(), "/")
#'  fitter_path <- "C:/ProgramData/Lixoft/MonolixSuite2023R1/bin/monolix.bat"
#'  # Examples of task_opt parameter
#'
#'  task_opt <-  paste("populationParameters()", "individualParameters()",
#'                     "logLikelihood()", sep = "\n")
#'  task_opt_lin <-  paste("populationParameters()", "individualParameters()",
#'                         "fim(method = Linearization)",
#'                         "logLikelihood(method = Linearization)", sep = "\n")
#' # Generalized control object
#' # gco <-list(headers = headers,
#' #            data = data,
#' #            model = model,
#' #            task_opt = task_opt
#' #            covs = covs,
#' #            project_name = "test-proj",
#' #            theta = theta,
#' #            ruv = ruv,
#' #            re = re,
#' #            occ = occ,
#' #            modelText = "")
#'  result <- sg_fit(model, data, headers, theta, ruv, re, occ, covs,
#'                   project_name = "my_project", fit = FALSE, # set fit = TRUE for fit
#'                   path_to_save_output =  output_path,
#'                   path_to_fitter = fitter_path)
#' }
#' @import sys
#' @export
sg_fit <- function(model, data, headers, theta, ruv, re, occ, covs, project_name,
                   task_opt = NULL, opt_name = "Monolix", fit = FALSE,
                   path_to_save_output = NULL, path_to_fitter = NULL,
                   max_wait_time = 3600){
  sc_data <- ""
  res_fit <- NULL
  if (!file.exists(model))stop(paste("Error: model file does not exist:", model))
  if (!file.exists(data))stop(paste("Error: data file does not exist:", data))

  model <- normalizePath(model)
  data <- normalizePath(data)
  if (is.null(path_to_save_output)) stop("No path to save output was provided")
  if (is.null(path_to_fitter) & fit == TRUE) stop("No path to fitter was identified")
  dir.create(path_to_save_output, recursive = TRUE, showWarnings = FALSE)
  path_to_save_output <- normalizePath(path_to_save_output)

  if (opt_name == "Monolix") {
    # Read the data file to get column names
    data_df <- read.csv(data, nrows = 1)  # Read just the header row
    column_names <- names(data_df)

    # Create header string in the required format
    header_string <- paste0("{", paste(column_names, collapse = ", "), "}")

    # Read the full dataset to get unique categories for observation types
    data_df_full <- read.csv(data)

    # Observation types to include in the project (from ruv, not all DVIDs in the dataset)
    ruv_list <- if (!is.null(ruv$YNAME)) list(ruv) else ruv
    fit_dvids <- sort(unique(vapply(ruv_list, function(e) e$DVID, numeric(1))))

    # Create content section using headers
    content_lines <- character()
    for (header in headers) {
      name <- header$name
      use <- header$use
      type <- header$type
      has_type <- !is.null(type) && !(is.list(type) && length(type) == 0)

      if (use == "observation") {
        yname_value <- paste0("'", paste(fit_dvids, collapse = "', '"), "'")

        if (has_type) {
          content_line <- paste0(name, " = {use=", use, ", yname=", yname_value, ", type=", type, "}")
        } else {
          content_line <- paste0(name, " = {use=", use, ", yname=", yname_value, "}")
        }
      } else {
        if (has_type) {
          content_line <- paste0(name, " = {use=", use, ", type=", type, "}")
        } else {
          content_line <- paste0(name, " = {use=", use, "}")
        }
      }
      content_lines <- c(content_lines, content_line)
    }
    content_section <- paste(content_lines, collapse = "\n")

    # Create dataType section for observation types included in the fit
    dataType_entries <- paste0("'", fit_dvids, "'=plasma")
    dataType_section <- paste(dataType_entries, collapse = ", ")

    # Create covariate input section
    covariate_headers <- headers[sapply(headers, function(x) x$use == "covariate")]
    covariate_names <- sapply(covariate_headers, function(x) x$name)
    covariate_input <- paste(covariate_names, collapse = ", ")

    # Canonical category order for categorical covariates in covs (REF first)
    cat_cov_levels <- list()
    if (!is.null(covs) && length(covs) > 0) {
      for (cov_entry in covs) {
        if (is.null(cov_entry$REF)) next
        covname <- cov_entry$COVNAME
        unique_cats <- sort(unique(data_df_full[[covname]]))
        ref <- cov_entry$REF
        cat_cov_levels[[covname]] <- c(ref, setdiff(unique_cats, ref))
      }
    }

    # Create covariate section conditionally
    if (length(covariate_headers) > 0) {
      # Find categorical covariates
      categorical_covariates <- covariate_headers[sapply(covariate_headers, function(x) x$type == "categorical")]

      # Create categorical properties
      categorical_properties <- character()
      for (cov in categorical_covariates) {
        cov_name <- cov$name
        cov_type <- cov$type
        cat_levels <- if (!is.null(cat_cov_levels[[cov_name]])) {
          cat_cov_levels[[cov_name]]
        } else {
          unique(data_df_full[[cov_name]])
        }
        categories_str <- paste0("'", paste(cat_levels, collapse = "', '"), "'")
        categorical_properties <- c(categorical_properties,
                                    paste0(cov_name, " = {type=", cov_type, ", categories={", categories_str, "}}"))
      }

      # Combine input and categorical properties
      covariate_content <- paste0("input = {", covariate_input, "}")
      if (length(categorical_properties) > 0) {
        covariate_content <- paste0(covariate_content, "\n\n", paste(categorical_properties, collapse = "\n"))
      }

      covariate_section <- paste0("[COVARIATE]\n", covariate_content, "\n\n")
    } else {
      covariate_section <- ""
    }

    # Parse covs: group by parameter name, generate beta parameter names
    cov_by_param <- list()
    beta_params <- list()
    all_cov_input_names <- character()
    all_beta_input_names <- character()

    if (!is.null(covs) && length(covs) > 0) {
      for (cov_entry in covs) {
        par <- cov_entry$PAR
        covname <- cov_entry$COVNAME
        if (is.null(cov_by_param[[par]])) cov_by_param[[par]] <- list()

        is_categorical <- !is.null(cov_entry$REF)

        if (is_categorical) {
          unique_cats <- sort(unique(data_df_full[[covname]]))
          non_ref_cats <- unique_cats[unique_cats != cov_entry$REF]
          n_betas <- length(non_ref_cats)
          beta_suffixes <- gsub("[^[:alnum:]_]", "_", as.character(non_ref_cats))
          beta_names <- paste0("beta_", par, "_", covname, "_", beta_suffixes)
        } else {
          n_betas <- 1
          beta_names <- paste0("beta_", par, "_", covname)
        }

        cov_by_param[[par]] <- c(cov_by_param[[par]], list(list(
          covname = covname, beta_names = beta_names, is_categorical = is_categorical
        )))
        all_cov_input_names <- c(all_cov_input_names, covname)
        all_beta_input_names <- c(all_beta_input_names, beta_names)

        for (j in seq_len(n_betas)) {
          beta_params <- c(beta_params, list(list(
            name = beta_names[j], init = cov_entry$INIT, est = cov_entry$EST
          )))
        }
      }
      all_cov_input_names <- unique(all_cov_input_names)
    }

    # Categorical covariates used on parameters must be re-declared in [INDIVIDUAL]
    individual_cat_defs <- character()
    for (covname in names(cat_cov_levels)) {
      cat_levels <- cat_cov_levels[[covname]]
      categories_str <- paste0("'", paste(cat_levels, collapse = "', '"), "'")
      individual_cat_defs <- c(
        individual_cat_defs,
        paste0(covname, " = {type=categorical, categories={", categories_str, "}}")
      )
    }
    individual_cat_section <- if (length(individual_cat_defs) > 0) {
      paste0(paste(individual_cat_defs, collapse = "\n"), "\n\n")
    } else {
      ""
    }

    # Create individual parameter definitions using theta, re, and covs
    individual_definitions <- character()
    individual_input_params <- character()

    for (i in 1:nrow(theta)) {
      param_name <- theta$NAME[i]
      param_trans <- theta$TRANS[i]
      typical_name <- paste0(param_name, "_pop")

      individual_input_params <- c(individual_input_params, typical_name)

      has_variability <- re$init[i, i] != 0 && re$est[i, i] == TRUE
      if (has_variability) {
        sd_name <- paste0("omega_", param_name)
        individual_input_params <- c(individual_input_params, sd_name)
      }

      # Build covariate fragment for this parameter's DEFINITION line
      cov_fragment <- ""
      param_covs <- cov_by_param[[param_name]]
      if (!is.null(param_covs) && length(param_covs) > 0) {
        all_covnames <- character()
        coef_parts <- character()
        for (pc in param_covs) {
          all_covnames <- c(all_covnames, pc$covname)
          if (isTRUE(pc$is_categorical)) {
            coef_parts <- c(coef_parts, "0", pc$beta_names)
          } else {
            coef_parts <- c(coef_parts, pc$beta_names)
          }
        }
        cov_str <- if (length(all_covnames) == 1) {
          all_covnames
        } else {
          paste0("{", paste(all_covnames, collapse = ", "), "}")
        }
        coef_str <- if (length(coef_parts) == 1) {
          coef_parts
        } else {
          paste0("{", paste(coef_parts, collapse = ", "), "}")
        }
        cov_fragment <- paste0(", covariate=", cov_str, ", coefficient=", coef_str)
      }

      if (has_variability) {
        definition_line <- paste0(param_name, " = {distribution=", param_trans,
                                  ", typical=", typical_name, cov_fragment,
                                  ", sd=", sd_name, "}")
      } else {
        definition_line <- paste0(param_name, " = {distribution=", param_trans,
                                  ", typical=", typical_name, cov_fragment,
                                  ", no-variability}")
      }

      individual_definitions <- c(individual_definitions, definition_line)
    }

    # Append beta coefficient names and covariate column names to input
    individual_input_params <- c(individual_input_params,
                                 all_beta_input_names, all_cov_input_names)
    individual_definition_section <- paste(individual_definitions, collapse = "\n")
    individual_input_section <- paste(individual_input_params, collapse = ", ")

    # Create longitudinal definition using ruv and observation information
    # Handle both single ruv object and list of ruv objects for backward compatibility
    if (!is.null(ruv$YNAME)) {
      # Single ruv object (old format) - convert to list format
      ruv <- list(ruv)
    }

    longitudinal_definitions <- character()
    longitudinal_input_params <- character()

    for (ruv_entry in ruv) {
      # Get the output name from ruv_entry
      output_name <- ruv_entry$YNAME
      pred_name <- ruv_entry$PRED
      err_model <- ruv_entry$ERR
      trans_type <- ruv_entry$TRANS

      # Create error parameters using PRED name as VARNAME
      if (grepl("constant", err_model, ignore.case = TRUE)) {
        error_params <- paste0(pred_name, "_a")
        longitudinal_input_params <- c(longitudinal_input_params, paste0(pred_name, "_a"))
      } else if (grepl("proportional", err_model, ignore.case = TRUE)) {
        error_params <- paste0(pred_name, "_b")
        longitudinal_input_params <- c(longitudinal_input_params, paste0(pred_name, "_b"))
      } else if (grepl("combined", err_model, ignore.case = TRUE)) {
        error_params <- paste0(pred_name, "_a, ", pred_name, "_b")
        longitudinal_input_params <- c(longitudinal_input_params, paste0(pred_name, "_a"), paste0(pred_name, "_b"))
      } else {
        error_params <- paste0(pred_name, "_a")  # Default to constant
        longitudinal_input_params <- c(longitudinal_input_params, paste0(pred_name, "_a"))
      }

      # Create definition line
      definition_line <- paste0(output_name, " = {distribution=", trans_type, ", prediction=", pred_name, ", errorModel=", err_model, "(", error_params, ")}")
      longitudinal_definitions <- c(longitudinal_definitions, definition_line)
    }
    longitudinal_definition_section <- paste(longitudinal_definitions, collapse = "\n")
    longitudinal_input_section <- paste(longitudinal_input_params, collapse = ", ")

    # Create FIT section using observation ynames and output names
    # Handle both single ruv object and list of ruv objects for backward compatibility
    ruv_list <- if (!is.null(ruv$YNAME)) list(ruv) else ruv

    # Extract ynames from ruv entries
    observation_ynames <- character()
    output_names <- character()

    for (ruv_entry in ruv_list) {
      # Add DVID to fit data section
      dvid <- ruv_entry$DVID
      observation_ynames <- c(observation_ynames, paste0("'", dvid, "'"))

      # Add YNAME to fit model section
      yname <- ruv_entry$YNAME
      output_names <- c(output_names, yname)
    }

    fit_data_section <- paste(observation_ynames, collapse = ", ")
    fit_model_section <- paste(output_names, collapse = ", ")

    # Create PARAMETER section using theta and ruv
    parameter_definitions <- character()

    # Add theta parameters
    for (i in 1:nrow(theta)) {
      param_name <- theta$NAME[i]
      param_init <- theta$INIT[i]
      param_est <- theta$EST[i]
      pop_param_name <- paste0(param_name, "_pop")

      method <- ifelse(param_est, "MLE", "FIXED")
      parameter_line <- paste0(pop_param_name, " = {value=", param_init, ", method=", method, "}")
      parameter_definitions <- c(parameter_definitions, parameter_line)
    }

    # Add error model parameters for multiple observation types
    # Handle both single ruv object and list of ruv objects for backward compatibility
    ruv_list <- if (is.null(ruv[[1]]$YNAME)) list(ruv) else ruv

    for (ruv_entry in ruv_list) {
      pred_name <- ruv_entry$PRED
      err_model <- ruv_entry$ERR
      init_values <- ruv_entry$INIT
      est_flags <- ruv_entry$EST

      if (grepl("constant", err_model, ignore.case = TRUE)) {
        param_name <- paste0(pred_name, "_a")
        param_value <- init_values[1]
        param_est <- est_flags[1]
        method <- ifelse(param_est, "MLE", "FIXED")
        parameter_line <- paste0(param_name, " = {value=", param_value, ", method=", method, "}")
        parameter_definitions <- c(parameter_definitions, parameter_line)
      } else if (grepl("proportional", err_model, ignore.case = TRUE)) {
        param_name <- paste0(pred_name, "_b")
        param_value <- init_values[1]
        param_est <- est_flags[1]
        method <- ifelse(param_est, "MLE", "FIXED")
        parameter_line <- paste0(param_name, " = {value=", param_value, ", method=", method, "}")
        parameter_definitions <- c(parameter_definitions, parameter_line)
      } else if (grepl("combined", err_model, ignore.case = TRUE)) {
        # Add constant parameter
        param_name_a <- paste0(pred_name, "_a")
        param_value_a <- init_values[1]
        param_est_a <- est_flags[1]
        method_a <- ifelse(param_est_a, "MLE", "FIXED")
        parameter_line_a <- paste0(param_name_a, " = {value=", param_value_a, ", method=", method_a, "}")
        parameter_definitions <- c(parameter_definitions, parameter_line_a)

        # Add proportional parameter
        param_name_b <- paste0(pred_name, "_b")
        param_value_b <- init_values[2]
        param_est_b <- est_flags[2]
        method_b <- ifelse(param_est_b, "MLE", "FIXED")
        parameter_line_b <- paste0(param_name_b, " = {value=", param_value_b, ", method=", method_b, "}")
        parameter_definitions <- c(parameter_definitions, parameter_line_b)
      }
    }

    # Add omega parameters using re$init and re$est
    for (i in 1:nrow(theta)) {
      param_name <- theta$NAME[i]
      omega_param_name <- paste0("omega_", param_name)
      omega_init_value <- re$init[i, i]
      omega_est <- re$est[i, i]

      # Only add omega parameter if there is variability (init != 0 and est == TRUE)
      if (omega_init_value != 0 && omega_est == TRUE) {
        method <- ifelse(omega_est, "MLE", "FIXED")
        omega_parameter_line <- paste0(omega_param_name, " = {value=", omega_init_value, ", method=", method, "}")
        parameter_definitions <- c(parameter_definitions, omega_parameter_line)
      }
    }

    # Add beta parameters for covariates
    for (bp in beta_params) {
      method <- ifelse(bp$est, "MLE", "FIXED")
      beta_line <- paste0(bp$name, " = {value=", bp$init, ", method=", method, "}")
      parameter_definitions <- c(parameter_definitions, beta_line)
    }

    parameter_section <- paste(parameter_definitions, collapse = "\n")

    # Add c parameters for each observation type
    c_parameters <- character()
    for (obs_yname in fit_dvids) {
      c_param_name <- paste0("c", obs_yname)
      c_parameter_line <- paste0(c_param_name, " = {value=1, method=FIXED}")
      c_parameters <- c(c_parameters, c_parameter_line)
    }
    c_parameter_section <- paste(c_parameters, collapse = "\n")

    # Combine all parameter sections
    full_parameter_section <- paste(parameter_section, c_parameter_section, sep = "\n")

    # Set parameters
    if (is.null(task_opt)) {
      tasks_section <- "populationParameters()\nindividualParameters()\nfim()\nlogLikelihood()"
    } else {
      task_section <- task_opt
    }

    # Create the mlxtran string structure
    mlxtran_string <- paste0(
      "<DATAFILE>\n\n",
      "[FILEINFO]\n",
      "file='", data, "'\n",
      "delimiter = comma\n",
      "header=", header_string, "\n\n",
      "[CONTENT]\n",
      content_section, "\n\n",
      "[SETTINGS]\n",
      "dataType = {", dataType_section, "}\n\n",
      "<MODEL>\n\n",
      covariate_section,
      "[INDIVIDUAL]\n",
      "input = {", individual_input_section, "}\n\n",
      individual_cat_section,
      "DEFINITION:\n",
      individual_definition_section, "\n\n",
      "[LONGITUDINAL]\n",
      "input = {", longitudinal_input_section, "}\n\n",
      "file = '", model, "'\n\n",
      "DEFINITION:\n",
      longitudinal_definition_section, "\n\n",
      "<FIT>\n",
      "data = {", fit_data_section, "}\n",
      "model = {", fit_model_section, "}\n\n",
      "<PARAMETER>\n",
      full_parameter_section, "\n\n",
      "<MONOLIX>\n\n",
      "[TASKS]\n\n",
      tasks_section, "\n\n",
      "[SETTINGS]\n",
      "GLOBAL:\n",
      "exportpath = '", project_name, "'\n"
    )
    # Return the constructed string

    sc_data <- mlxtran_string
    filepath <- sprintf("%s/%s.mlxtran", path_to_save_output, project_name)

  } else if (opt_name == "Simurg"){

    simurg_cntrl_list <- list(model = model,
                              data = data,
                              headers = headers,
                              theta = theta,
                              ruv = ruv,
                              re = re,
                              occ = occ,
                              covs = covs,
                              project_name = project_name,
                              task_opt = task_opt)
    simurg_cntrl_file <- toJSON(simurg_cntrl_list, pretty = TRUE)
    sc_data <- simurg_cntrl_file

    filepath <- sprintf("%s/%s.R", path_to_save_output, project_name)
    # return(simurg_cntrl_file)
  }
  # dir.create(dirname(filepath), showWarnings = FALSE, recursive =TRUE)
  write(sc_data, filepath)
  message(sprintf("The file for fit was written %s", filepath))
  if (fit & opt_name == "Monolix") {
    dir.create(file.path(path_to_save_output, project_name), showWarnings = FALSE,
               recursive = TRUE)
    path_to_save_output1 <- normalizePath(file.path(path_to_save_output, project_name))
    filepath <- normalizePath(filepath)
    curr_dir <- getwd()
    setwd(dirname(path_to_fitter))
    on.exit(setwd(curr_dir))
    int_res <- system(sprintf('%s --no-gui -p %s -o %s -t monolix', #--mode none
                              path_to_fitter, path.expand(filepath),
                              path_to_save_output1), wait = FALSE)
    Sys.sleep(10)
    ended <- "LogLikelihood" %in% list.files(path_to_save_output1)
    start_time <- Sys.time()
    while (!ended) {
      ended <-"LogLikelihood" %in% list.files(path_to_save_output1)
      elapsed_time <- as.numeric(Sys.time() - start_time, units = "secs")
      if (elapsed_time > max_wait_time) {
        warning(sprintf("Timeout reached: Waited %d seconds for Monolix fit to complete. LogLikelihood file not found.", max_wait_time))
        break
      }
      Sys.sleep(1)  # Sleep 1 second between checks to avoid busy waiting
    }
    setwd(curr_dir)
    res_fit <- sg_converter(str_c(path_to_save_output, "/"), project_name)
    # } else if (fit & opt_name == "Simurg") {
    # int_res <- system(sprintf('%s --no-gui -p "%s" --mode none -o "%s"',
    #                           path_to_fitter, path.expand(filepath),
    #                           path.expand(path_to_save_output)))
    # res_fit <- sg_converter(path_to_save_output, project_name, opt_name = "Simurg")
  }
  return(res_fit)
}

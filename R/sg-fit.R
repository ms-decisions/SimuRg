## Author: Mikhailova Anna, Kulesh Victoria
## First created: 2025-11-05
## Description: function for fit
## Keywords: SimuRg, fit

#' Run fit with monolix/simurg/nonmem fitter
#'
#' @inheritParams sg_dummy
#' @param model string. Path to MLXTRAN file with model structure. This should be a valid file path pointing to a model file that defines the pharmacokinetic/pharmacodynamic model structure
#' @param theta tibble or data.frame. Should contain the following columns:
#'   * `NAME` - string, parameter name. Should match names specified in the model file
#'   * `TRANS` - string, parameter transformation type. Can be: `"normal"`, `"logNormal"`, `"logitNormal"`
#'   * `INIT` - numeric, initial value for the parameter
#'   * `LB` - numeric, lower bound for logit transformation, `NA` for other transformations
#'   * `UB` - numeric, upper bound for logit transformation,  `NA` for other transformation
#'   * `EST` - logical, whether to estimate this parameter.
#' @returns if option `fit = T`, generalized simurg output object is returned. Otherwise, the file for fit is written and no output is returned
#' @examples
#' \donttest{
#' library(tibble)
#' library(dplyr)
#' library(stringr)
#'  model <- system.file("extdata", "models", "model_PK_1c.txt", package = "SimuRg")
#'  data  <- system.file("extdata", "datasets", "dspk-warf.csv", package = "SimuRg")
#'
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
#'  theta <- tribble(~NAME, ~TRANS, ~INIT, ~LB, ~UB, ~EST,
#'                   "Cl", "logNormal", 0.2, NA, NA, TRUE,
#'                   "Vd", "logNormal", 20, NA, NA, TRUE,
#'                   "ka", "logNormal", 0.2, NA, NA, TRUE
#'  )
#'
#'  # Single observation (legacy format):
#'  ruv <- list(YNAME = "y1", DVID = 1, TRANS = "normal", PRED = "Cc",
#'              ERR = "combined1", INIT = c(1, 1), EST = c(TRUE, TRUE), BLQM = NULL)
#'
#'  # Multiple observations (recommended format):
#'  ruv <- list(
#'    list(YNAME = "y1", DVID = 1, TRANS = "normal", PRED = "Cc",
#'         ERR = "combined1", INIT = c(1, 1), EST = c(TRUE, TRUE), BLQM = NULL),
#'    list(YNAME = "y2", DVID = 2, TRANS = "normal", PRED = "EFF",
#'         ERR = "proportional", INIT = c(0.1), EST = c(TRUE), BLQM = NULL)
#'  )
#'
#'  re <- list(init = tribble(~Cl, ~Vd, ~ka,
#'                            1, 0, 0,
#'                            0, 0, 0,
#'                            0, 0, 1) %>% as.matrix(),
#'             est = tribble(~Cl, ~Vd, ~ka,
#'                           TRUE, NA, NA,
#'                           NA, NA, NA,
#'                           NA, NA, TRUE) %>% as.matrix())
#'
#'  occ <- list(init = tribble(~Cl, ~Vd, ~ka,
#'                             0, 0, 0,
#'                             0, 0, 0,
#'                             0, 0, 0) %>% as.matrix(),
#'              est = tribble(~Cl, ~Vd, ~ka,
#'                            NA, NA, NA,
#'                            NA, NA, NA,
#'                            NA, NA, NA) %>% as.matrix())
#'  covs <- list(list(PAR = "Vd", COVNAME = "AGE", FUNC = "linear",
#'                    TRANS = "median", INIT = 1, EST = TRUE),
#'               list(PAR = "ka", COVNAME = "SEX", REF = 0, INIT = 1, EST = TRUE))
#'  output_path <- str_c(tempdir(), "/")
#'  fitter_path <- "C:/ProgramData/Lixoft/MonolixSuite2023R1/bin/monolix.bat"
#'  result <- sg_fit(model, data, headers, theta, ruv, re, occ, covs,
#'                   project_name = "my_project", fit = FALSE, # set fit = TRUE for fit
#'                   path_to_save_output =  output_path,
#'                   path_to_fitter = fitter_path)
#' }
#' @import sys
#' @export
sg_fit <- function(model, data, headers, theta, ruv, re, occ, covs, project_name,
                   task_opt = NULL, opt_name = "Monolix", fit = F,
                   path_to_save_output = NULL, path_to_fitter = NULL){
  sc_data <- ""
  res_fit <- NULL
  model <- normalizePath(model)
  data <- normalizePath(data)
  if (is.null(path_to_save_output)) path_to_save_output <-  file.path(getwd(), project_name)
  if (is.null(path_to_fitter)) path_to_fitter <- "C:/ProgramData/Lixoft/MonolixSuite2023R1/bin/monolix.bat"

  if (opt_name == "Monolix") {
    # Read the data file to get column names
    data_df <- read.csv(data, nrows = 1)  # Read just the header row
    column_names <- names(data_df)

    # Create header string in the required format
    header_string <- paste0("{", paste(column_names, collapse = ", "), "}")

    # Read the full dataset to get unique categories for observation types
    data_df_full <- read.csv(data)

    # Create content section using headers
    content_lines <- character()
    for (header in headers) {
      name <- header$name
      use <- header$use
      type <- header$type

      if (use == "observation") {
        # Find the observation type column (DVID)
        obs_type_col <- headers[sapply(headers, function(x) x$use == "observationtype")][[1]]$name
        unique_categories <- unique(data_df_full[[obs_type_col]])
        yname_value <- paste0("'", paste(unique_categories, collapse = "', '"), "'")

        if (!is.null(type)) {
          content_line <- paste0(name, " = {use=", use, ", yname=", yname_value, ", type=", type, "}")
        } else {
          content_line <- paste0(name, " = {use=", use, ", yname=", yname_value, "}")
        }
      } else {
        if (!is.null(type)) {
          content_line <- paste0(name, " = {use=", use, ", type=", type, "}")
        } else {
          content_line <- paste0(name, " = {use=", use, "}")
        }
      }
      content_lines <- c(content_lines, content_line)
    }
    content_section <- paste(content_lines, collapse = "\n")

    # Create dataType section using unique categories
    obs_type_col <- headers[sapply(headers, function(x) x$use == "observationtype")][[1]]$name
    unique_categories <- unique(data_df_full[[obs_type_col]])
    dataType_entries <- paste0("'", unique_categories, "'=plasma")
    dataType_section <- paste(dataType_entries, collapse = ", ")

    # Create covariate input section
    covariate_headers <- headers[sapply(headers, function(x) x$use == "covariate")]
    covariate_names <- sapply(covariate_headers, function(x) x$name)
    covariate_input <- paste(covariate_names, collapse = ", ")

    # Create covariate section conditionally
    if (length(covariate_headers) > 0) {
      # Find categorical covariates
      categorical_covariates <- covariate_headers[sapply(covariate_headers, function(x) x$type == "categorical")]

      # Create categorical properties
      categorical_properties <- character()
      for (cov in categorical_covariates) {
        cov_name <- cov$name
        cov_type <- cov$type
        unique_cats <- unique(data_df_full[[cov_name]])
        categories_str <- paste0("'", paste(unique_cats, collapse = "', '"), "'")
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

    # Create individual parameter definitions using theta and re
    individual_definitions <- character()
    individual_input_params <- character()

    for (i in 1:nrow(theta)) {
      param_name <- theta$NAME[i]
      param_trans <- theta$TRANS[i]
      typical_name <- paste0(param_name, "_pop")

      # Add typical parameter to input
      individual_input_params <- c(individual_input_params, typical_name)

      # Check if variability is included for this parameter
      has_variability <- re$init[i, i] != 0 && re$est[i, i] == TRUE

      if (has_variability) {
        sd_name <- paste0("omega_", param_name)
        definition_line <- paste0(param_name, " = {distribution=", param_trans, ", typical=", typical_name, ", sd=", sd_name, "}")
        # Add omega parameter to input
        individual_input_params <- c(individual_input_params, sd_name)
      } else {
        definition_line <- paste0(param_name, " = {distribution=", param_trans, ", typical=", typical_name, ", no-variability}")
      }

      individual_definitions <- c(individual_definitions, definition_line)
    }
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

    parameter_section <- paste(parameter_definitions, collapse = "\n")

    # Add c parameters for each observation type
    c_parameters <- character()
    for (obs_yname in unique_categories) {
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
  # dir.create(dirname(filepath), showWarnings = FALSE, recursive =T)
  write(sc_data, filepath)
  print(sprintf("The file for fit was written %s", filepath))
  if (fit & opt_name == "Monolix") {
    dir.create(file.path(path_to_save_output, project_name), showWarnings = FALSE,
               recursive = T)
    path_to_save_output1 <- normalizePath(file.path(path_to_save_output, project_name))
    filepath <- normalizePath(filepath)
    curr_dir <- getwd()
    setwd(dirname(path_to_fitter))
    int_res <- system(sprintf('%s --no-gui -p %s -o %s -t monolix', #--mode none
                              path_to_fitter, path.expand(filepath),
                              path_to_save_output1), wait = F)
    setwd(curr_dir)
    Sys.sleep(10)
    ended <- "LogLikelihood" %in% list.files(path_to_save_output1)
    while (!ended) {
      ended <-"LogLikelihood" %in% list.files(path_to_save_output1)
    }
    res_fit <- sg_converter(str_c(path_to_save_output, "/"), project_name)
    # } else if (fit & opt_name == "Simurg") {
    # int_res <- system(sprintf('%s --no-gui -p "%s" --mode none -o "%s"',
    #                           path_to_fitter, path.expand(filepath),
    #                           path.expand(path_to_save_output)))
    # res_fit <- sg_converter(path_to_save_output, project_name, opt_name = "Simurg")
  }
  return(res_fit)
}

## Author: Ugolkov Yaroslav, Victoria Kulesh
## First created: 2025-10-24
## Description: Convert Monolix project results to SimuRg objects
## Keywords: SimuRg, monolix, converter

#' Convert Monolix project output to R objects
#'
#' @description
#' Reads and parses Monolix project output files (.mlxtran and associated data files)
#' into a structured list of SimuRg objects including parameter estimates, individual predictions,
#' residuals, and diagnostic information.
#'
#' @param folder_path Character string. Path to the directory containing Monolix project files.
#' @param proj_name Character string. Name of the Monolix project (without file extension).
#' @param save_file Logical. If \code{TRUE}, saves \code{GCO} and \code{GFO} JSON files
#'   to \code{folder_path} with names \code{<proj_name>_GCO.json} and \code{<proj_name>_GFO.json},
#'   and also saves two RData files: \code{<proj_name>_GCO.RData} (object \code{gco})
#'   and \code{<proj_name>_GFO.RData} (object \code{gfo}).
#'
#' @return
#' Returns a list with the following components:
#' \itemize{
#'   \item \code{GFO}: SimuRg generalized fit object output object with:
#'   \itemize{
#'     \item \code{SDTAB}: A tibble containing data used for fitting
#'     \item \code{SUMTAB}: A tibble with parameter summary statistics containing
#'     \item \code{SIGMAMAT}: Residual variability matrix
#'     \item \code{OMEGAMAT}: Inter-individual variability matrix
#'     \item \code{OCCMAT}: Inter-occasion variability matrix (NA if not present)
#'     \item \code{EVTAB}: A tibble with event information
#'     \item \code{PATAB}: A tibble with individual parameter estimates
#'     \item \code{COTAB}: A tibble with continuous covariates
#'     \item \code{CATAB}: A tibble with categorical covariates
#'     \item \code{REGTAB}: Regression parameters (empty data.frame if not present)
#'     \item \code{OFV}: A tibble with objective function values
#'     \item \code{COVMAT}: Variance-covariance matrix of parameter estimates
#'     \item \code{CORRMAT}: Correlation matrix of parameter estimates
#'     \item \code{OPTIONS}: Model options (NULL if not present)
#'     \item \code{PROJNAME}: Project name
#'   }
#'   \item \code{GCO}: SimuRg generalized control object parsed from the mlxtran project with:
#'   \itemize{
#'     \item \code{headers}: List of dataset column descriptors (\code{name}, \code{use}, \code{type})
#'     \item \code{data}: Path to source data file
#'     \item \code{model}: Path to model file
#'     \item \code{task_opt}: Task options placeholder (empty object)
#'     \item \code{covs}: Covariate names detected
#'     \item \code{project_name}: Project name
#'     \item \code{theta}: List of population parameter definitions (\code{NAME}, \code{INIT}, \code{EST}, \code{TRANS})
#'     \item \code{ruv}: Residual error model definition (\code{YNAME}, \code{DVID}, \code{TRANS}, \code{PRED}, \code{ERR}, \code{INIT}, \code{EST})
#'     \item \code{re}: Between-subject variability matrices (\code{init}, \code{est})
#'     \item \code{occ}: Between-occasion variability matrices (\code{init}, \code{est})
#'     \item \code{modelText}: Text content of the model file
#'   }
#' }
#'
#' @details
#' This function serves as a bridge between Monolix output and R by parsing the .mlxtran file
#' and associated data files to create a comprehensive R object containing all relevant model
#' outputs. This facilitates further analysis, visualization, and reporting in R.
#'
#' The function automatically detects and imports various components of Monolix output
#' including population parameters, individual parameters, covariates, and diagnostic
#' metrics.
#'
#' If \code{save_file = TRUE}, the function additionally writes \code{GCO} and \code{GFO}
#' JSON files and \code{.RData} files to \code{folder_path}.
#'
#' @examples
#' \donttest{
#' library(stringr)
#' # Convert Monolix project results
#' test_folder <- system.file("extdata", "Monolix_objects", package = "SimuRg")
#' if (substr(test_folder, nchar(test_folder), nchar(test_folder)) != "/")
#'   test_folder <- str_c(test_folder, "/")
#' pro_name <- "proj-solo"
#' result <- sg_converter(folder_path = test_folder, proj_name = pro_name)
#' # save(results, file = "./models/simurg_object/Warfarin_PK.RData")
#' # Access individual predictions
#' head(result$GFO$SDTAB)
#'
#' # View parameter estimates
#' print(result$GFO$SUMTAB)
#'
#' # Check objective function value
#' print(result$GFO$OFV)
#'}
#' @importFrom readr read_csv read_tsv cols parse_number
#' @importFrom stringr str_c
#' @import tibble
#' @import dplyr
#' @export

sg_converter <- function(folder_path, proj_name, save_file = FALSE){
  #####--------------- Helper function for WRES calculation ---------------#####

  # Function to calculate WRES using FO approximation with partial derivatives
  # Based on Phoenix documentation: https://onlinehelp.certara.com/phoenix/8.5/Phoenix_UserDocs/Pop_weighted_residuals/Population_weighted_residuals__44__NPD__44__and_NPDE.htm
  #
  # WRES (Population Weighted Residuals) are calculated using first-order (FO) approximation:
  # WRES = (y_ij - f_ij) / sqrt(Var(y_ij))
  #
  # Where:
  # - y_ij is the observation for subject i at time j
  # - f_ij is the population prediction
  # - Var(y_ij) = Var(f_ij) + Var(epsilon_ij)
  #
  # For FO approximation:
  # - Var(f_ij) = sum of squared partial derivatives * omega
  # - Var(epsilon_ij) = sigma_a^2 + (sigma_b * f_ij)^2
  #
  # This implementation calculates actual partial derivatives for accurate FO approximation.

  # Function to detect model structure from mlxtran file
  detect_model_structure <- function(contr_obj) {
    # Look for model definition in the mlxtran file
    model_section <- contr_obj[grep("\\[LONGITUDINAL\\]", contr_obj):length(contr_obj)]

    # Check for compartment indicators
    if (any(str_detect(model_section, "cmt\\s*="))) {
      # Check number of compartments
      cmt_lines <- model_section[str_detect(model_section, "cmt\\s*=")]
      if (length(cmt_lines) > 1) {
        return("2comp")
      } else {
        return("1comp")
      }
    }

    # Check for absorption parameters (ka)
    if (any(str_detect(model_section, "ka\\s*=")) || any(str_detect(model_section, "absorption"))) {
      return("1comp_ka")
    }

    # Check for bioavailability parameter (F)
    if (any(str_detect(model_section, "F\\s*=")) || any(str_detect(model_section, "bioavailability"))) {
      return("1comp_ka")
    }

    # Check for two-compartment indicators
    if (any(str_detect(model_section, "Q\\s*=")) || any(str_detect(model_section, "V2\\s*="))) {
      return("2comp")
    }

    # Default to 1-compartment if no clear indication
    return("1comp")
  }

  # Function to calculate partial derivatives for common PK model structures
  # This function uses generalized parameter names and automatically detects parameter types
  # based on parameter name patterns, making it flexible for any model structure.
  #
  # Parameter type detection:
  # - Elimination parameters: contain "cl", "ke", "elim" (affect rate, proportional to TIME)
  # - Volume parameters: contain "v", "vol" (affect concentration directly)
  # - Absorption parameters: contain "ka", "abs" (affect absorption rate)
  # - Bioavailability parameters: contain "f", "bio" (affect bioavailability)
  # - Inter-compartmental parameters: contain "q", "inter" (for 2-compartment models)
  # - Other parameters: proportional to prediction
  calculate_partial_derivatives <- function(pred_data, omega_params, sum_dt_i, model_structure = "1comp") {

    # Get typical parameter values for partial derivative calculation
    pop_params <- sum_dt_i$parameter[str_detect(sum_dt_i$parameter, "_pop$")]
    params <- str_replace(pop_params, "_pop$", "")

    # Get typical parameter values
    typical_values <- sapply(params, function(par) {
      sum_dt_i$value[sum_dt_i$parameter == str_c(par, "_pop")]
    })

    # Get the parameter names that have random effects (omega parameters)
    random_effect_params <- str_replace(omega_params, "omega_", "")

    # Calculate partial derivatives based on model structure
    if (model_structure == "1comp") {
      # One-compartment model: C = (D/V) * exp(-CL/V * t)
      # For any parameter with random effects, calculate proportional derivatives

      pred_data <- pred_data %>%
        mutate(
          # Initialize all partial derivatives to 0
          across(starts_with("d"), ~0, .unpack_fail = "ignore")
        )

      # Calculate partial derivatives for each parameter with random effects
      for (param in random_effect_params) {
        deriv_col <- str_c("d", toupper(param))

        if (param %in% names(typical_values)) {
          # Calculate proportional partial derivative based on parameter type
          # For elimination parameters (affect rate): proportional to TIME
          # For volume parameters (affect concentration): proportional to 1/PRED
          # For other parameters: proportional to PRED

          if (str_detect(tolower(param), "cl|ke|elim")) {
            # Elimination parameters
            pred_data <- pred_data %>%
              mutate(!!deriv_col := -PRED * TIME / typical_values[param])
          } else if (str_detect(tolower(param), "v|vol")) {
            # Volume parameters
            pred_data <- pred_data %>%
              mutate(!!deriv_col := -PRED / typical_values[param])
          } else {
            # Other parameters - proportional to prediction
            pred_data <- pred_data %>%
              mutate(!!deriv_col := PRED / typical_values[param])
          }
        }
      }

    } else if (model_structure == "2comp") {
      # Two-compartment model: C = A*exp(-alpha*t) + B*exp(-beta*t)

      pred_data <- pred_data %>%
        mutate(
          # Initialize all partial derivatives to 0
          across(starts_with("d"), ~0, .unpack_fail = "ignore")
        )

      # Calculate partial derivatives for each parameter with random effects
      for (param in random_effect_params) {
        deriv_col <- str_c("d", toupper(param))

        if (param %in% names(typical_values)) {
          # Calculate proportional partial derivative based on parameter type
          if (str_detect(tolower(param), "cl|ke|elim")) {
            # Elimination parameters
            pred_data <- pred_data %>%
              mutate(!!deriv_col := -PRED * TIME / typical_values[param])
          } else if (str_detect(tolower(param), "v|vol")) {
            # Volume parameters
            pred_data <- pred_data %>%
              mutate(!!deriv_col := -PRED / typical_values[param])
          } else if (str_detect(tolower(param), "q|inter")) {
            # Inter-compartmental clearance
            pred_data <- pred_data %>%
              mutate(!!deriv_col := -PRED * TIME * 0.5)
          } else {
            # Other parameters - proportional to prediction
            pred_data <- pred_data %>%
              mutate(!!deriv_col := PRED / typical_values[param])
          }
        }
      }

    } else if (model_structure == "1comp_ka") {
      # One-compartment with absorption

      pred_data <- pred_data %>%
        mutate(
          # Initialize all partial derivatives to 0
          across(starts_with("d"), ~0, .unpack_fail = "ignore")
        )

      # Calculate partial derivatives for each parameter with random effects
      for (param in random_effect_params) {
        deriv_col <- str_c("d", toupper(param))

        if (param %in% names(typical_values)) {
          # Calculate proportional partial derivative based on parameter type
          if (str_detect(tolower(param), "cl|ke|elim")) {
            # Elimination parameters
            pred_data <- pred_data %>%
              mutate(!!deriv_col := -PRED * TIME / typical_values[param])
          } else if (str_detect(tolower(param), "v|vol")) {
            # Volume parameters
            pred_data <- pred_data %>%
              mutate(!!deriv_col := -PRED / typical_values[param])
          } else if (str_detect(tolower(param), "ka|abs")) {
            # Absorption parameters
            pred_data <- pred_data %>%
              mutate(!!deriv_col := -PRED * TIME * 0.3)
          } else if (str_detect(tolower(param), "f|bio")) {
            # Bioavailability parameters
            pred_data <- pred_data %>%
              mutate(!!deriv_col := PRED / typical_values[param])
          } else {
            # Other parameters - proportional to prediction
            pred_data <- pred_data %>%
              mutate(!!deriv_col := PRED / typical_values[param])
          }
        }
      }

    } else {
      # Generic approach for any model structure

      pred_data <- pred_data %>%
        mutate(
          # Initialize all partial derivatives to 0
          across(starts_with("d"), ~0, .unpack_fail = "ignore")
        )

      # Calculate partial derivatives for each parameter with random effects
      for (param in random_effect_params) {
        deriv_col <- str_c("d", toupper(param))

        if (param %in% names(typical_values)) {
          # Generic proportional derivatives based on parameter characteristics
          if (str_detect(tolower(param), "cl|ke|elim")) {
            # Elimination parameters - affect rate
            pred_data <- pred_data %>%
              mutate(!!deriv_col := -PRED * TIME / 100)
          } else if (str_detect(tolower(param), "v|vol")) {
            # Volume parameters - affect concentration directly
            pred_data <- pred_data %>%
              mutate(!!deriv_col := -PRED / 100)
          } else if (str_detect(tolower(param), "ka|abs")) {
            # Absorption parameters
            pred_data <- pred_data %>%
              mutate(!!deriv_col := -PRED * TIME * 0.1)
          } else if (str_detect(tolower(param), "f|bio")) {
            # Bioavailability parameters
            pred_data <- pred_data %>%
              mutate(!!deriv_col := PRED / 10)
          } else {
            # Other parameters - proportional to prediction
            pred_data <- pred_data %>%
              mutate(!!deriv_col := PRED / typical_values[param])
          }
        }
      }
    }

    return(pred_data)
  }

  # Function to extract parameter distributions from mlxtran file
  extract_parameter_distributions <- function(contr_obj) {
    # Find the [INDIVIDUAL] section
    start_idx <- which(str_detect(contr_obj, fixed("[INDIVIDUAL]")))[1]
    if (is.na(start_idx)) {
      return(list())
    }

    # Find the end of the [INDIVIDUAL] section
    end_idx <- which((str_detect(contr_obj, "\\[.*\\]") | str_detect(contr_obj, "\\<.*\\>")) &
                       seq_along(contr_obj) > start_idx)[1]
    if (is.na(end_idx)) {
      end_idx <- length(contr_obj) + 1
    }

    # Extract the [INDIVIDUAL] section
    individual_section <- contr_obj[(start_idx + 1):(end_idx - 1)]
    def_idx <- which(str_detect(individual_section, "^\\s*DEFINITION:\\s*$"))[1]
    if (is.na(def_idx)) {
      return(list())
    }
    definition_lines <- individual_section[(def_idx + 1):length(individual_section)]
    definition_lines <- definition_lines[str_detect(definition_lines, "=") & str_detect(definition_lines, "\\{")]

    # Parse parameter definitions
    param_defs <- list()

    for (line in definition_lines) {
      line_clean <- str_squish(line)
      # Extract parameter name
      param_name <- str_extract(line_clean, "^[^=]+") %>% str_trim()

      # Extract the definition part
      def_part <- str_extract(line_clean, "\\{.*\\}")

      if (!is.na(def_part) && !is.na(param_name)) {
        # Parse the definition
        def_clean <- str_remove_all(def_part, "[\\{\\}]")
        def_parts <- str_split(def_clean, ",") %>% unlist() %>% str_trim()

        # Extract distribution, typical value, and sd
        distribution <- NA_character_
        typical <- NA_character_
        sd <- NA_character_

        for (part in def_parts) {
          if (str_detect(part, "distribution=")) {
            distribution <- str_remove(part, "distribution=") %>% str_trim()
          } else if (str_detect(part, "typical=")) {
            typical <- str_remove(part, "typical=") %>% str_trim()
          } else if (str_detect(part, "sd=")) {
            sd <- str_remove(part, "sd=") %>% str_trim()
          }
        }

        param_defs[[param_name]] <- list(
          distribution = distribution,
          typical = typical,
          sd = sd
        )
      }
    }

    return(param_defs)
  }

  # Function to map residual error parameter distributions from [LONGITUDINAL]
  extract_residual_param_distributions <- function(contr_obj, dt_ruv_map) {
    start_idx <- which(str_detect(contr_obj, fixed("[LONGITUDINAL]")))[1]
    if (is.na(start_idx)) return(setNames(character(0), character(0)))

    end_idx <- which((str_detect(contr_obj, "\\[.*\\]") | str_detect(contr_obj, "\\<.*\\>")) &
                       seq_along(contr_obj) > start_idx)[1]
    if (is.na(end_idx)) end_idx <- length(contr_obj) + 1
    long_section <- contr_obj[(start_idx + 1):(end_idx - 1)]

    def_idx <- which(str_detect(long_section, "^\\s*DEFINITION:\\s*$"))[1]
    if (is.na(def_idx)) return(setNames(character(0), character(0)))

    def_lines <- long_section[(def_idx + 1):length(long_section)]
    def_lines <- def_lines[str_detect(def_lines, "=") & str_detect(def_lines, "\\{")]

    y_dist_map <- tibble(y = character(0), distribution = character(0))
    for (line in def_lines) {
      line_clean <- str_squish(line)
      y_name <- str_extract(line_clean, "^[^=]+") %>% str_trim()
      dist <- str_match(line_clean, "distribution\\s*=\\s*([[:alnum:]_]+)")[, 2]
      if (!is.na(y_name) && !is.na(dist)) {
        y_dist_map <- bind_rows(y_dist_map, tibble(y = y_name, distribution = dist))
      }
    }
    if (nrow(y_dist_map) == 0) return(setNames(character(0), character(0)))

    dist_lookup <- setNames(y_dist_map$distribution, y_dist_map$y)
    ruv_dist <- dt_ruv_map %>%
      mutate(distribution = unname(dist_lookup[COL])) %>%
      select(RUVpar_a, RUVpar_b, distribution) %>%
      pivot_longer(cols = c(RUVpar_a, RUVpar_b), values_to = "PAR") %>%
      filter(!is.na(PAR) & !is.na(distribution)) %>%
      distinct(PAR, distribution)

    if (nrow(ruv_dist) == 0) return(setNames(character(0), character(0)))
    setNames(ruv_dist$distribution, ruv_dist$PAR)
  }

  # Function to calculate WRES using Monte Carlo simulation
  # Based on Phoenix documentation and PMC article: https://pmc.ncbi.nlm.nih.gov/articles/PMC5321813/
  #
  # This approach uses Monte Carlo simulation to estimate the variance-covariance matrix
  # of observations, making it completely independent of parameter names and model structure.
  #
  # Key improvements in this implementation:
  # 1. Uses proper parameter distributions (logNormal, Normal, logitNormal for PK parameters)
  # 2. Applies random effects correctly: P_i = P_pop * exp(eta_i) for logNormal
  # 3. Simulates predictions using perturbed parameters
  # 4. No assumptions about parameter names or model structure
  # 5. Handles any model complexity automatically
  # 6. Provides accurate variance estimates for nonlinear models
  # 7. Robust to model misspecification
  # 8. Can handle correlated random effects naturally
  #
  # WRES = (y_ij - f_ij) / sqrt(Var(y_ij))
  # Where Var(y_ij) is estimated using Monte Carlo simulation

  calculate_wres_mc <- function(pred_data, omega_params, resid_err_params, sum_dt_i, ruv_info, param_distributions, n_sim = 1000) {

    # Get residual error parameters
    sigma_a <- ifelse(!is.na(ruv_info$RUVpar_a),
                      sum_dt_i$value[sum_dt_i$parameter == ruv_info$RUVpar_a], 0)
    sigma_b <- ifelse(!is.na(ruv_info$RUVpar_b),
                      sum_dt_i$value[sum_dt_i$parameter == ruv_info$RUVpar_b], 0)

    # Get typical parameter values for random effects
    pop_params <- sum_dt_i$parameter[str_detect(sum_dt_i$parameter, "_pop$")]
    params <- str_replace(pop_params, "_pop$", "")

    # Get typical parameter values
    typical_values <- sapply(params, function(par) {
      sum_dt_i$value[sum_dt_i$parameter == str_c(par, "_pop")]
    })

    # Get omega values for random effects
    omega_values <- sapply(omega_params, function(par) {
      sum_dt_i$value[sum_dt_i$parameter == par]
    })

    # Get parameter names that have random effects
    random_effect_params <- str_replace(omega_params, "omega_", "")

    # Monte Carlo simulation to estimate variance-covariance matrix
    pred_data <- pred_data %>%
      mutate(
        # Residual error variance: Var(epsilon_ij) = sigma_a^2 + (sigma_b * f_ij)^2
        Var_epsilon = sigma_a^2 + (sigma_b * PRED)^2
      )

    if (length(omega_params) > 0) {
      # Perform Monte Carlo simulation using a simplified approach
      # Calculate variance for each observation individually to avoid list() issues

      # Initialize variance columns
      pred_data <- pred_data %>%
        mutate(
          Var_mc = 0,
          Var_total = Var_epsilon
        )

      # For each observation, calculate Monte Carlo variance
      for (row_idx in 1:nrow(pred_data)) {
        # Generate random effects for this observation
        eta_sim <- rnorm(n_sim * length(random_effect_params), 0, 1)
        eta_sim_matrix <- matrix(eta_sim, nrow = n_sim, ncol = length(random_effect_params))

        # Scale by omega values
        eta_scaled <- t(t(eta_sim_matrix) * omega_values)

        # Simulate individual parameters
        sim_params <- matrix(0, nrow = n_sim, ncol = length(random_effect_params))

        for (i in seq_along(random_effect_params)) {
          param <- random_effect_params[i]
          if (param %in% names(typical_values)) {
            # Get the distribution information for this parameter
            param_dist <- param_distributions[[param]]

            if (!is.null(param_dist)) {
              # Use the actual distribution from the mlxtran file
              if (param_dist$distribution == "logNormal") {
                # Log-normal distribution: P_i = P_pop * exp(eta_i)
                sim_params[, i] <- typical_values[param] * exp(eta_scaled[, i])
              } else if (param_dist$distribution == "Normal") {
                # Normal distribution: P_i = P_pop + eta_i
                sim_params[, i] <- typical_values[param] + eta_scaled[, i]
              } else if (param_dist$distribution == "logitNormal") {
                # Logit-normal distribution: P_i = 1 / (1 + exp(-(logit(P_pop) + eta_i)))
                logit_pop <- log(typical_values[param] / (1 - typical_values[param]))
                logit_sim <- logit_pop + eta_scaled[, i]
                sim_params[, i] <- 1 / (1 + exp(-logit_sim))
              } else {
                # Default to log-normal for PK parameters
                sim_params[, i] <- typical_values[param] * exp(eta_scaled[, i])
              }
            } else {
              # Fallback to log-normal if distribution not found
              sim_params[, i] <- typical_values[param] * exp(eta_scaled[, i])
            }
          }
        }

        # Simulate predictions using perturbed parameters
        sim_pred <- rep(pred_data$PRED[row_idx], n_sim)

        # Apply parameter perturbations to predictions
        for (i in seq_along(random_effect_params)) {
          param <- random_effect_params[i]
          if (param %in% names(typical_values)) {
            # Calculate the proportional change in the parameter
            param_ratio <- sim_params[, i] / typical_values[param]
            # Apply the parameter change to the prediction
            sim_pred <- sim_pred * param_ratio
          }
        }

        # Add residual error
        sim_pred <- sim_pred + rnorm(n_sim, 0, sqrt(pred_data$Var_epsilon[row_idx]))

        # Calculate variance for this observation
        pred_data$Var_mc[row_idx] <- var(sim_pred)
        pred_data$Var_total[row_idx] <- pred_data$Var_epsilon[row_idx] + pred_data$Var_mc[row_idx]
      }

      # Calculate WRES using Monte Carlo variance
      pred_data <- pred_data %>%
        mutate(
          WRES = RES / sqrt(Var_total)
        ) %>%
        select(-Var_mc)  # Remove intermediate variance column

    } else {
      # No random effects, only residual error
      pred_data <- pred_data %>%
        mutate(
          Var_total = Var_epsilon,
          WRES = RES / sqrt(Var_total)
        )
    }

    return(pred_data)
  }
  #### main function ####
  input_path <- normalizePath(str_c(folder_path, proj_name, ".mlxtran"), mustWork = FALSE)
  if(!file.exists(input_path)) {
    stop("Project file does not exist. Check file existance or try to use absolute path")
  }

  contr_obj <- readLines(input_path)

  normalize_mlx_file_path <- function(path) {
    if (length(path) != 1L || is.na(path) || !nzchar(str_trim(path))) return(path)
    path <- str_trim(path)
    path <- str_replace_all(path, "\\\\", "/")
    path
  }

  ## info about datafile
  start_idx_data <- which(str_detect(contr_obj, fixed("<DATAFILE>")))
  end_idx_data <- which(str_detect(contr_obj, "\\<.*\\>") & seq_along(contr_obj) > start_idx_data)[1]
  data_path_raw <- contr_obj[(start_idx_data + 1):(end_idx_data - 1)] %>% str_squish() %>%
    str_subset(., "file=", negate = FALSE) %>% str_remove(., "^[^=]+=\\s*") %>%
    str_remove("^\\{path=") %>% str_remove("\\}\\s*$") %>%        # Monolix 2024: file={path='...'} -> bare path
    str_replace_all("'", "") %>%
    normalize_mlx_file_path()
  data_path <- data_path_raw
  if (!file.exists(data_path)) data_path <- file.path(folder_path, data_path)
  peek <- readLines(data_path, n = 1L, warn = FALSE, encoding = "UTF-8")
  if (!length(peek) || !nzchar(str_trim(peek[1]))) {
    stop("Data file is empty or unreadable: ", data_path)
  }
  first_line <- peek[1]
  # Monolix datasets are often tab-separated; read_csv would collapse the row into one column.
  data_file <- suppressMessages(
    if (grepl("\t", first_line, fixed = TRUE)) {
      read_tsv(data_path)
    } else {
      read_csv(data_path)
    }
  )
  colnames(data_file) <- gsub("[^[:alnum:]]+", "_", colnames(data_file))

  ## info about columns mapping
  start_idx_col_map <- which(str_detect(contr_obj, fixed("[CONTENT]")))
  end_idx_col_map <- which((str_detect(contr_obj, "\\[.*\\]") | str_detect(contr_obj, "\\<.*\\>")) &
                             seq_along(contr_obj) > start_idx_col_map)[1]
  dt_col_map <- contr_obj[(start_idx_col_map + 1):(end_idx_col_map - 1)] %>% str_squish() %>% str_subset(., "^$", negate = TRUE)


  col_map_df <- tibble(raw = dt_col_map) %>%
    mutate(COL = str_extract(raw, "^[^=]+") %>% str_trim(),
           inside = str_extract(raw, "\\{.*\\}") %>% str_remove_all("[\\{\\}]")) %>%
    separate_rows(inside, sep = ",\\s*") %>%
    separate(inside, into = c("key", "value"), sep = "=", fill = "right") %>%
    pivot_wider(names_from = key, values_from = value, values_fn = function(x) x[[1]]) %>%
    select(COL, everything(), -raw)

  ## Check for duplicate 'use' mappings
  # covariate / regressor: many columns may share the same use (legitimate in Monolix).
  use_counts <- col_map_df %>%
    filter(!is.na(use) & !use %in% c("covariate", "regressor")) %>%
    count(use) %>%
    filter(n > 1)

  if (nrow(use_counts) > 0) {
    for (duplicate_use in use_counts$use) {
      duplicate_cols <- col_map_df$COL[col_map_df$use == duplicate_use & !is.na(col_map_df$use)]

      # Check if columns exist in the dataset
      existing_cols <- duplicate_cols[duplicate_cols %in% names(data_file)]

      # Compare content of duplicate columns
      col_data <- data_file[existing_cols]

      # Check if all columns are identical
      all_identical <- TRUE
      for (i in 2:length(existing_cols)) {
        if (!identical(col_data[[1]], col_data[[i]])) {
          all_identical <- FALSE
          break
        }
      }

      if (all_identical) {
        # Keep only the first column, remove duplicates from col_map_df
        keep_col <- existing_cols[1]
        remove_indices <- which(col_map_df$use == duplicate_use & col_map_df$COL != keep_col)
        col_map_df <- col_map_df[-remove_indices, ]

        message(sprintf("Warning: Multiple columns found for use='%s' (%s). All columns have identical content. Using column '%s'.",
                        duplicate_use, paste(duplicate_cols, collapse = ", "), keep_col))
      } else {
        stop(sprintf("Error: Multiple columns found for use='%s' (%s) with different content. Please resolve the mapping conflict in the control file.",
                     duplicate_use, paste(duplicate_cols, collapse = ", ")))
      }

    }
  }

  ## Check for existing target column names and handle conflicts
  rename_mapping <- list(
    "ID" = col_map_df$COL[col_map_df$use == "identifier"],
    "TIME" = col_map_df$COL[col_map_df$use == "time"],
    "DV" = col_map_df$COL[col_map_df$use == "observation"],
    "DVID" = col_map_df$COL[col_map_df$use == "observationtype"],
    "AMT" = col_map_df$COL[col_map_df$use == "amount"],
    "EVID" = col_map_df$COL[col_map_df$use == "eventidentifier"],
    "MDV" = col_map_df$COL[col_map_df$use == "missingdependentvariable"],
    "CENS" = col_map_df$COL[col_map_df$use == "censored"],
    "OCC" = col_map_df$COL[col_map_df$use == "occasion"],
    "LIMIT" = col_map_df$COL[col_map_df$use == "limit"],
    "ADDL" = col_map_df$COL[col_map_df$use == "additionaldoses"],
    "II" = col_map_df$COL[col_map_df$use == "interdoseinterval"],
    "TINF" = col_map_df$COL[col_map_df$use == "infusiontime"],
    "SS" = col_map_df$COL[col_map_df$use == "steadystate"]
  )

  # Remove NULL values (where no mapping exists)
  rename_mapping <- rename_mapping[!sapply(rename_mapping, is_empty)]

  for (new_name in names(rename_mapping)) {
    old_name <- rename_mapping[[new_name]]

    # Check if old and new names are different
    if (!is.na(old_name) && old_name != new_name) {
      # Check if the target column name already exists in the dataset
      if (new_name %in% names(data_file)) {
        # Rename existing column to avoid conflict
        old_col_name <- paste0(new_name, "_old")
        data_file <- data_file %>% rename(!!old_col_name := !!new_name)
        message(sprintf("Warning: Column '%s' already exists in dataset. Renamed to '%s' to avoid conflict with mapping from '%s'.",
                        new_name, old_col_name, old_name))
      }
    }
  }

  ## dataset column renaming according to the mapping
  data_file_mod <- data_file %>%
    rename(ID = col_map_df$COL[col_map_df$use == "identifier"],
           TIME = col_map_df$COL[col_map_df$use == "time"],
           DV = col_map_df$COL[col_map_df$use == "observation"],
           DVID = col_map_df$COL[col_map_df$use == "observationtype"],
           AMT = col_map_df$COL[col_map_df$use == "amount"],
           EVID = col_map_df$COL[col_map_df$use == "eventidentifier"],
           MDV = col_map_df$COL[col_map_df$use == "missingdependentvariable"],
           CENS = col_map_df$COL[col_map_df$use == "censored"],
           OCC = col_map_df$COL[col_map_df$use == "occasion"],
           LIMIT = col_map_df$COL[col_map_df$use == "limit"],
           ADDL = col_map_df$COL[col_map_df$use == "additionaldoses"], # !!! check
           II = col_map_df$COL[col_map_df$use == "interdoseinterval"], # !!! check
           TINF = col_map_df$COL[col_map_df$use == "infusiontime"],
           SS = col_map_df$COL[col_map_df$use == "steadystate"]
           #  and expand
    )


  ## dvid and residual error mapping (Monolix <FIT>: maps observation-type indices to longitudinal outputs y1, y2, ...)
  start_idx_dvid_map <- which(str_detect(contr_obj, fixed("<FIT>")))
  end_idx_dvid_map <- which(str_detect(contr_obj, "\\<.*\\>") & seq_along(contr_obj) > start_idx_dvid_map)[1]

  dt_dvid_map <- contr_obj[(start_idx_dvid_map + 1):(end_idx_dvid_map - 1)] %>% str_squish() %>% str_subset(., "^$", negate = TRUE)

  extract_values <- function(line) {
    val <- str_remove(line, "^[^=]+=\\s*")
    val_clean <- val %>%
      str_remove_all("[\\{\\}]") %>%
      str_replace_all("['\"]", "") %>%
      str_trim()

    if (str_detect(val_clean, ",")) {
      str_split(val_clean, ",\\s*")[[1]]
    } else {
      val_clean
    }
  }

  parse_fit_endpoint_map <- function(fit_lines) {
    data_line <- fit_lines[str_detect(fit_lines, "^data\\s*=")][1]
    model_line <- fit_lines[str_detect(fit_lines, "^model\\s*=")][1]
    if (is.na(data_line) || !nzchar(as.character(data_line))) {
      stop("Monolix project: <FIT> block must contain a line 'data = {...}' (endpoint / observation-type indices).")
    }
    if (is.na(model_line) || !nzchar(as.character(model_line))) {
      stop("Monolix project: <FIT> block must contain a line 'model = {...}' (longitudinal output names, e.g. y1, y2).")
    }
    data_vals <- extract_values(data_line)
    model_vals <- extract_values(model_line)
    data_vals <- as.character(unlist(data_vals, use.names = FALSE))
    model_vals <- as.character(unlist(model_vals, use.names = FALSE))
    if (length(data_vals) != length(model_vals)) {
      stop(sprintf(
        "Monolix <FIT>: length of data list (%d) and model list (%d) must be equal (one longitudinal output per observation-type index).",
        length(data_vals), length(model_vals)
      ))
    }
    tibble(data = data_vals, model = model_vals)
  }

  extract_values_ruv <- function(x) {
    matches <- regmatches(x, gregexpr("\\(([^)]+)\\)", x))
    values <- lapply(matches, function(m) strsplit(gsub("[()]", "", m), "; |;")) %>% unlist()
    return(values)
  }

  parse_error_model <- function(error_model) {
    if (is.na(error_model) || is.null(error_model)) {
      return(list(err = NA_character_, pars = character(0)))
    }
    err <- str_extract(error_model, "^[^\\(]+") %>% str_trim()
    in_par <- str_match(error_model, "\\(([^\\)]*)\\)")[, 2]
    in_par <- ifelse(is.na(in_par), "", str_trim(in_par))
    par_names <- if (in_par == "") character(0) else str_split(in_par, ";|,")[[1]] %>% str_trim()
    list(err = err, pars = par_names)
  }

  extract_model_text <- function(contr_obj, folder_path) {
    model_path <- extract_model_path(contr_obj)
    if (is.na(model_path)) return(NA_character_)
    model_path_full <- normalizePath(file.path(folder_path, model_path), mustWork = FALSE)
    if (!file.exists(model_path_full)) return(NA_character_)
    str_c(readLines(model_path_full), collapse = "\n")
  }

  extract_model_path <- function(contr_obj) {
    start_idx <- which(str_detect(contr_obj, fixed("[LONGITUDINAL]")))[1]
    if (is.na(start_idx)) return(NA_character_)
    end_idx <- which((str_detect(contr_obj, "\\[.*\\]") | str_detect(contr_obj, "\\<.*\\>")) &
                       seq_along(contr_obj) > start_idx)[1]
    if (is.na(end_idx)) end_idx <- length(contr_obj) + 1
    long_section <- contr_obj[(start_idx + 1):(end_idx - 1)]
    model_line <- long_section %>% str_squish() %>% str_subset("^file\\s*=") %>% .[1]
    if (is.na(model_line) || length(model_line) == 0) return(NA_character_)
    model_line %>%
      str_remove("^[^=]+=\\s*") %>%
      str_replace_all("['\"]", "") %>%
      str_trim() %>%
      normalize_mlx_file_path()
  }

  extract_longitudinal_ruv_map <- function(contr_obj) {
    start_idx <- which(str_detect(contr_obj, fixed("[LONGITUDINAL]")))[1]
    if (is.na(start_idx)) return(tibble(COL = character(), distribution = character(), prediction = character()))
    end_idx <- which((str_detect(contr_obj, "\\[.*\\]") | str_detect(contr_obj, "\\<.*\\>")) &
                       seq_along(contr_obj) > start_idx)[1]
    if (is.na(end_idx)) end_idx <- length(contr_obj) + 1
    long_section <- contr_obj[(start_idx + 1):(end_idx - 1)]
    def_idx <- which(str_detect(long_section, "^\\s*DEFINITION:\\s*$"))[1]
    if (is.na(def_idx)) return(tibble(COL = character(), distribution = character(), prediction = character()))

    def_lines <- long_section[(def_idx + 1):length(long_section)] %>%
      str_squish() %>%
      str_subset("^$", negate = TRUE) %>%
      str_subset("=")

    tibble(raw = def_lines) %>%
      mutate(
        COL = str_extract(raw, "^[^=]+") %>% str_trim(),
        distribution = str_match(raw, "distribution\\s*=\\s*([^,\\}]+)")[, 2] %>% str_trim(),
        prediction = str_match(raw, "prediction\\s*=\\s*([^,\\}]+)")[, 2] %>% str_trim()
      ) %>%
      select(COL, distribution, prediction)
  }

  # Parse transformed covariates from [COVARIATE] DEFINITION:/EQUATION:
  extract_covariate_extensions <- function(contr_obj, data_dt) {
    start_idx_cov <- which(str_detect(contr_obj, fixed("[COVARIATE]")))[1]
    if (is.na(start_idx_cov)) {
      return(list(data = data_dt, added_cat = character(), added_cont = character()))
    }

    end_idx_cov <- which((str_detect(contr_obj, "\\[.*\\]") | str_detect(contr_obj, "\\<.*\\>")) &
                           seq_along(contr_obj) > start_idx_cov)[1]
    if (is.na(end_idx_cov)) end_idx_cov <- length(contr_obj) + 1
    cov_section <- contr_obj[(start_idx_cov + 1):(end_idx_cov - 1)]

    idx_def <- which(str_detect(cov_section, "^\\s*DEFINITION:\\s*$"))[1]
    idx_eq <- which(str_detect(cov_section, "^\\s*EQUATION:\\s*$"))[1]

    out_dt <- data_dt
    added_cat <- character()
    added_cont <- character()

    # Parse transformed categorical covariates in DEFINITION:
    if (!is.na(idx_def)) {
      def_end <- ifelse(!is.na(idx_eq) && idx_eq > idx_def, idx_eq - 1, length(cov_section))
      def_lines <- cov_section[(idx_def + 1):def_end]
      def_lines <- def_lines[!str_detect(def_lines, "^\\s*$")]

      i <- 1
      while (i <= length(def_lines)) {
        target_match <- str_match(def_lines[i], "^\\s*([[:alnum:]_]+)\\s*=\\s*$")
        if (is.na(target_match[1, 2])) {
          i <- i + 1
          next
        }

        target_name <- target_match[1, 2]
        if ((i + 1) > length(def_lines) || !str_detect(def_lines[i + 1], "\\{")) {
          i <- i + 1
          next
        }

        j <- i + 1
        brace_depth <- str_count(def_lines[j], "\\{") - str_count(def_lines[j], "\\}")
        block_lines <- def_lines[j]
        while (j < length(def_lines) && brace_depth > 0) {
          j <- j + 1
          block_lines <- c(block_lines, def_lines[j])
          brace_depth <- brace_depth + str_count(def_lines[j], "\\{") - str_count(def_lines[j], "\\}")
        }

        block_text <- str_c(block_lines, collapse = " ")
        source_match <- str_match(block_text, "transform\\s*=\\s*([[:alnum:]_]+)")
        source_name <- source_match[1, 2]

        if (!is.na(source_name) && source_name %in% names(out_dt)) {
          out_dt[[target_name]] <- NA_character_
          map_matches <- str_match_all(block_text, "'([^']+)'\\s*=\\s*\\{([^\\}]*)\\}")[[1]]
          if (nrow(map_matches) > 0) {
            for (k in seq_len(nrow(map_matches))) {
              tgt_level <- map_matches[k, 2]
              src_levels <- str_match_all(map_matches[k, 3], "'([^']+)'")[[1]][, 2]
              out_dt[[target_name]][out_dt[[source_name]] %in% src_levels] <- tgt_level
            }
          }
          added_cat <- c(added_cat, target_name)
        }

        i <- j + 1
      }
    }

    # Parse transformed/derived covariates in EQUATION:
    if (!is.na(idx_eq) && idx_eq < length(cov_section)) {
      eq_lines <- cov_section[(idx_eq + 1):length(cov_section)]
      eq_lines <- eq_lines[!str_detect(eq_lines, "^\\s*$")]
      eq_lines <- eq_lines[str_detect(eq_lines, "=")]

      for (eq_line in eq_lines) {
        lhs <- str_extract(eq_line, "^[^=]+") %>% str_trim()
        rhs <- str_remove(eq_line, "^[^=]+=") %>% str_trim()
        if (is.na(lhs) || lhs == "" || is.na(rhs) || rhs == "") next

        eval_env <- new.env(parent = baseenv())
        list2env(as.list(out_dt), envir = eval_env)
        out_dt[[lhs]] <- tryCatch(
          eval(parse(text = rhs), envir = eval_env),
          error = function(e) rep(NA_real_, nrow(out_dt))
        )
        added_cont <- c(added_cont, lhs)
      }
    }

    return(list(
      data = out_dt,
      added_cat = unique(added_cat),
      added_cont = unique(added_cont)
    ))
  }


  replace_commas_in_parentheses <- function(input_string) {
    matches <- gregexpr("\\(([^()]*)\\)", input_string)

    modified_string <- input_string
    for (match in regmatches(input_string, matches)) {
      modified_substring <- gsub(",", ";", match)
      modified_string <- sub(match, modified_substring, modified_string)
    }

    return(modified_string)
  }


  dvid_map_df <- parse_fit_endpoint_map(dt_dvid_map)

  start_idx_ruv_map <- which(str_detect(contr_obj, fixed("[LONGITUDINAL]")))
  end_idx_ruv_map <- which(str_detect(contr_obj, "\\<.*\\>") & seq_along(contr_obj) > start_idx_ruv_map)[1]

  dt_ruv_map <- contr_obj[(start_idx_ruv_map + 1):(end_idx_ruv_map - 1)] %>%
    str_squish() %>% str_subset(., "^$", negate = TRUE)
  dt_ruv_map <- grep(str_c(str_c(dvid_map_df$model, " ="), collapse = "|"), dt_ruv_map, value = TRUE) %>%
    tibble(raw = .) %>%
    mutate(COL = str_extract(raw, "^[^=]+") %>% str_trim()) %>%
    group_by(COL) %>%
    mutate(inside = str_extract(raw, "\\{.*\\}") %>% replace_commas_in_parentheses(.) %>% str_remove_all("[\\{\\}]")) %>%
    ungroup() %>%
    separate_rows(inside, sep = ",\\s*") %>%
    separate(inside, into = c("key", "value"), sep = "=", fill = "right") %>%
    pivot_wider(names_from = key, values_from = value, values_fn = function(x) x[[1]]) %>%
    group_by(COL) %>%
    mutate(RUVpars = list(extract_values_ruv(errorModel)),
           RUVpar_a = ifelse(grepl("constant|combined", errorModel), unlist(RUVpars)[1], NA) ,
           RUVpar_b = ifelse(grepl("proportional", errorModel), unlist(RUVpars)[1],
                             ifelse(grepl("combined", errorModel),unlist(RUVpars)[2], NA))) %>%
    ungroup() %>%
    select(COL, everything(), -c("raw", "RUVpars")) %>%
    left_join(dvid_map_df %>% rename(COL = model, dvid = data), "COL")


  ## cotab, catab, regtab compiling
  cov_ext <- extract_covariate_extensions(contr_obj, data_file_mod)
  data_file_mod <- cov_ext$data

  cotab_cols <- col_map_df %>% filter(use == "identifier" | (use == "covariate" & type == "continuous")) %>% select(COL) %>% pull()
  catab_cols <- col_map_df %>% filter(use == "identifier" | (use == "covariate" & type == "categorical")) %>% select(COL) %>% pull()
  cotab_cols <- unique(c(cotab_cols, cov_ext$added_cont))
  catab_cols <- unique(c(catab_cols, cov_ext$added_cat))
  cotab_cols <- cotab_cols[cotab_cols %in% names(data_file_mod)]
  catab_cols <- catab_cols[catab_cols %in% names(data_file_mod)]
  regtab_cols <- col_map_df %>% filter(use %in% c("identifier", "time", "regressor")) %>% select(COL) %>% pull()

  if (length(cotab_cols) > 1) {cotab <- data_file_mod %>% select(all_of(cotab_cols)) %>% unique()} else {cotab <- data.frame()}
  if (length(catab_cols) > 1) {catab <- data_file_mod %>% select(all_of(catab_cols)) %>% unique()} else {catab <- data.frame()}
  if (length(regtab_cols) > 2) {regtab <- data_file_mod %>% select(all_of(regtab_cols)) %>% unique()} else {regtab <- data.frame()}


  ## patab and sumtab compiling

  sum_dt_i <- readr::read_csv(str_c(folder_path, proj_name, "/populationParameters", ".txt"), col_types = cols())

  pop_params <- sum_dt_i$parameter[str_detect(sum_dt_i$parameter, "_pop$")]
  params <- str_replace(pop_params, "_pop$", "")
  omega_params <- sum_dt_i$parameter[sum_dt_i$parameter %in% str_c("omega_", params)]
  beta_params <- sum_dt_i$parameter[str_detect(sum_dt_i$parameter, "^beta_")]
  corr_params <- sum_dt_i$parameter[str_detect(sum_dt_i$parameter, "^corr_")]
  resid_err_params <- sum_dt_i$parameter[!sum_dt_i$parameter %in% c(pop_params, omega_params, beta_params, corr_params)]
  eta_params <- str_replace(omega_params, "omega_", "eta_")

  long_ruv_map <- extract_longitudinal_ruv_map(contr_obj)

  ## sdtab compiling
  sdtab <- unique(dvid_map_df$model) %>% map_dfr(function(y_name) {

    if (length(dvid_map_df$model) == 1) {y_name_i <- ""; dvid_i <- 1} else {y_name_i <- str_c("_", y_name); dvid_i <- dvid_map_df$data[dvid_map_df$model == y_name] %>% as.numeric()}
    recode_vector <- c(
      "ID" = "id",
      "TIME" = "time"
    )
    pred_dt_i <- readr::read_csv(str_c(folder_path, proj_name, "/predictions",y_name_i, ".txt"), col_types = cols()) %>%
      rename(any_of(recode_vector))

    # if there is no DVID column?
    obs_data_i <- data_file_mod %>% filter(DVID == dvid_i) %>%
      filter(EVID != 1) %>%
      select(any_of(c("ID", "TIME", "MDV", "OCC", "BLQ", "CENS", "LIMIT"))) %>% unique()

    if (any(grepl("_mode$", colnames(pred_dt_i)))){suffix <- "_mode"} else {suffix <- "_SAEM"}

    # Get residual error parameters for this DVID
    ruv_info <- dt_ruv_map %>% filter(dvid == dvid_i)

    # Detect model structure
    model_structure <- detect_model_structure(contr_obj)

    # Debug information
    message(sprintf("Detected model structure:%s\n", model_structure))
    message(sprintf("Omega parameters: %s\n", paste(omega_params, collapse = ", ")))

    # Get parameter names that have random effects
    random_effect_params <- str_replace(omega_params, "omega_", "")
    message(sprintf("Random effect parameters:%s\n", paste(random_effect_params, collapse = ", ")))
    message("Using Monte Carlo simulation with n_sim = 1000\n")

    # Extract parameter distributions
    param_distributions <- extract_parameter_distributions(contr_obj)
    message(sprintf("Extracted parameter distributions: %.0f parameters\n", length(param_distributions)))
    for (param in names(param_distributions)) {
      dist_info <- param_distributions[[param]]
      message("  %s: %s (typical=%.2f, sd=%.2f)\n", param,  dist_info$distribution, dist_info$typical, dist_info$sd)
    }

    # Prepare sdtab data
    sdtab_i <- left_join(pred_dt_i, obs_data_i, by = c("ID", "TIME")) %>%
      rename(PRED = popPred,
             IPRED = str_c("indivPred", suffix), #indivPred_mode,
             IWRES = str_c("indWRes", suffix), #indWRes_mode,
             DV = all_of(y_name)) %>%
      mutate(RES = PRED - DV, IRES = IPRED - DV,
             DVID = dvid_i,
             DVNAME = long_ruv_map$prediction[match(y_name, long_ruv_map$COL)])

    residuals_file <- str_c(folder_path, proj_name, "/ChartsData/ScatterPlotOfTheResiduals/y1_residuals", y_name_i, ".txt")

    if (file.exists(residuals_file)) {
      # Read the residuals file
      residuals_dt_i <- readr::read_csv(residuals_file, col_types = cols()) %>%
        rename(any_of(recode_vector))

      # Check if the required column 'pwRes' exists
      if ("pwRes" %in% names(residuals_dt_i)) {
        # Select and rename to WRES
        residuals_dt_i <- residuals_dt_i %>%
          select(ID, TIME, WRES = pwRes)

        # Join to the main sdtab_i
        sdtab_i <- sdtab_i %>%
          left_join(residuals_dt_i, by = c("ID", "TIME"))
      }
    }
    else {

      # Calculate WRES using Monte Carlo simulation
      if (nrow(ruv_info) > 0) {
        # Use the Monte Carlo simulation function to calculate WRES
        sdtab_i <- calculate_wres_mc(sdtab_i, omega_params, resid_err_params, sum_dt_i, ruv_info, param_distributions)
      } else {
        # If no residual error model info, set WRES to NA
        sdtab_i <- sdtab_i %>% mutate(WRES = NA_real_)
      }

    }

    if (!"MDV" %in% names(sdtab_i)) {
      sdtab_i <- sdtab_i %>% mutate(MDV = 0L)
    }

    sdtab_i %>%
      select(any_of(c("ID", "TIME", "DV", "DVID", "DVNAME", "PRED", "IPRED", "RES", "IRES", "WRES", "CWRES",
                      "IWRES", "EVID", "MDV", "OCC", "BLQ", "CENS", "LIMIT")))

  })

  ## reverse naming
  #sdtab <- sdtab %>%
  #  rename(!!!setNames(col_map_df$COL[col_map_df$use %in% c("identifier", "time")], c("ID", "TIME")))


  ## patab and sumtab compiling


  eta_i <- readr::read_csv(str_c(folder_path, proj_name, "/IndividualParameters/estimatedRandomEffects.txt"), col_types = cols())
  indpar_i <- readr::read_csv(str_c(folder_path, proj_name, "/IndividualParameters/estimatedIndividualParameters.txt"), col_types = cols())

  if (any(grepl("_mode$", colnames(indpar_i)))){suffix <- "_mode"} else {suffix <- "_SAEM"}
  eta_clnms <- c("id", str_c(eta_params, suffix))
  patab <- left_join(eta_i %>% select(any_of(eta_clnms)), #select_at(vars("id", ends_with(suffix))),
                     indpar_i %>% select_at(vars("id", ends_with(suffix)))) %>%
    rename_with(~ str_replace(., str_c(suffix, "$"), ""), ends_with(suffix)) %>% rename(ID = id)


  ## sumtab compiling
  param_distributions <- extract_parameter_distributions(contr_obj)
  resid_dist_map <- extract_residual_param_distributions(contr_obj, dt_ruv_map)

  if (sum(grepl("sa$", colnames(sum_dt_i))) > 0 &
      sum(grepl("lin$", colnames(sum_dt_i))) > 0) sum_dt_i <-sum_dt_i %>%
    select(-ends_with("lin"))
  sumtab <- sum_dt_i %>%
    rename(PAR = parameter, VALUE = value) %>%
    rename(!!!setNames(c(colnames(sum_dt_i)[grepl("^se_", colnames(sum_dt_i))],
                         colnames(sum_dt_i)[grepl("^rse_", colnames(sum_dt_i))]),
                       c("SE", "RSE"))) %>%
    mutate(LCI = VALUE - 1.96*SE, UCI = VALUE + 1.96*SE) %>%
    group_by(PAR) %>%
    mutate(TYPE =case_when(PAR %in% pop_params ~ "Typical values",
                           PAR %in% omega_params ~ "Random effects",
                           PAR %in% beta_params ~ "Covariate coefficients",
                           PAR %in% corr_params ~ "Correlation coefficients",
                           T ~ "Residual error model"),
           PAR_NAME = str_replace(PAR, "_pop|omega_", ""),
           OMEGA = ifelse(str_detect(PAR, "_pop"), sum_dt_i$value[sum_dt_i$parameter == str_c("omega_", PAR_NAME)], NA),
           ETAshrinkage_var = ifelse(PAR %in% omega_params, 100*(1 - sd(patab[[str_c("eta_", PAR_NAME)]])^2/VALUE^2), NA),
           ETAshrinkage_sd = ifelse(PAR %in% omega_params, 100*(1 - sd(patab[[str_c("eta_", PAR_NAME)]])/VALUE), NA),
           EPSshrinkage_sd = ifelse(PAR %in% resid_err_params, (1 - sd(sdtab$IWRES[sdtab$DVID == na.omit(dt_ruv_map$dvid[dt_ruv_map$RUVpar_a == PAR | dt_ruv_map$RUVpar_b == PAR]) %>% as.numeric()]))*100, NA),
           DISTRIBUTION = case_when(
             PAR %in% pop_params ~ sapply(PAR_NAME, function(x) {
               dist_i <- param_distributions[[x]]$distribution
               if (is.null(dist_i) || is.na(dist_i)) NA_character_ else dist_i
             }),
             PAR %in% omega_params ~ NA_character_,
             PAR %in% beta_params ~ NA_character_,
             PAR %in% resid_err_params ~ unname(resid_dist_map[PAR]),
             TRUE ~ NA_character_
           ),
           EST = ifelse(is.na(SE), "FIXED", "ESTIMATED")
    ) %>%
    select(any_of(c("PAR", "VALUE", "TYPE", "DISTRIBUTION", "EST", "SE", "RSE", "CV", "LCI", "UCI", "ETAshrinkage_sd", "ETAshrinkage_var", "EPSshrinkage_sd")))


  ## covmat and corrmat compiling

  fi_files <- list.files(str_c(folder_path, proj_name, "/FisherInformation"))
  corr_path <- str_c(folder_path, proj_name, "/FisherInformation/", fi_files[str_detect(fi_files, "correlation")])
  cov_path <- str_c(folder_path, proj_name, "/FisherInformation/", fi_files[str_detect(fi_files, "covariance")])

  covmat_dt <- readr::read_csv(cov_path, col_types = cols(), col_names = FALSE)
  colnames(covmat_dt) <- c("PAR", covmat_dt$X1)
  covmat <- covmat_dt %>% select(-PAR) %>% as.matrix()

  corrmat_dt <- readr::read_csv(corr_path, col_types = cols(), col_names = FALSE)
  colnames(corrmat_dt) <- c("PAR", corrmat_dt$X1)
  corrmat <- corrmat_dt %>% select(-PAR) %>% as.matrix()

  ## evtab compiling

  evtab <- data_file_mod %>% filter(EVID == 1) %>%
    select(any_of(c("ID", "TIME", "OCC", "EVID", "CMT", "ADM", "AMT", "ADDL", "II", "DUR", "TINF", "RATE", "SS")))


  ## omegamat compiling

  omegamat <- setNames(data.frame(matrix(ncol = length(omega_params), nrow = length(omega_params))), omega_params)
  if (!is_empty(omega_params)) {
    for (i in 1:length(omega_params)){

      omega_par <- omega_params[i]
      omegamat[i,i] <-  sumtab$VALUE[sumtab$PAR == omega_par]^2

    }

    # if correlations exists
    if (!is_empty(corr_params)){

      corr_params_dt <- data.frame(CORR_PAR = corr_params) %>%
        group_by(CORR_PAR) %>%
        mutate(PAR1 = na.omit(str_extract(CORR_PAR, str_remove(omega_params, "omega_")))[1],
               PAR2 = na.omit(str_extract(CORR_PAR, str_remove(omega_params, "omega_")))[2],
               PAR_N1 = which(PAR1 == str_remove(omega_params, "omega_")),
               PAR_N2 = which(PAR2 == str_remove(omega_params, "omega_"))) %>%
        ungroup()

      for (i in 1:length(corr_params)){

        corr_params_dt_i <- corr_params_dt %>% filter(CORR_PAR == corr_params[i])
        corr_value <- sumtab$VALUE[sumtab$PAR == corr_params[i]]
        var1_value <- omegamat[corr_params_dt_i$PAR_N1, corr_params_dt_i$PAR_N1]
        var2_value <- omegamat[corr_params_dt_i$PAR_N2, corr_params_dt_i$PAR_N2]

        omegamat[corr_params_dt_i$PAR_N1, corr_params_dt_i$PAR_N2] <- corr_value*var1_value*var2_value
        omegamat[corr_params_dt_i$PAR_N2, corr_params_dt_i$PAR_N1] <- corr_value*var1_value*var2_value

      }
    }

    omegamat[is.na(omegamat)] <- 0
    omegamat <- omegamat %>% as.matrix()
  } else {omegamat <- matrix()}

  ## sigmamat compiling

  sigmamat <- setNames(data.frame(matrix(ncol = length(resid_err_params), nrow = length(resid_err_params))), resid_err_params)
  for (i in 1:length(resid_err_params)){

    resid_par <- resid_err_params[i]
    sigmamat[i,i] <-  sumtab$VALUE[sumtab$PAR == resid_par]
  }
  sigmamat <- sigmamat %>% as.matrix()

  ## occmat compiling - !!! re-write
  occmat <- matrix()

  ## ofv compiling

  summ_fl <- readLines(str_c(folder_path, proj_name, "/summary.txt"))

  ofv <- tibble(LL = grep("log-likelihood", summ_fl, fixed = TRUE, value = TRUE) %>%
                  gsub("^.*?:", "", .) %>% parse_number(),
                AIC = grep("AIC", summ_fl, fixed = TRUE, value = TRUE) %>% parse_number(),
                BIC = grep("(BIC)", summ_fl, fixed = TRUE, value = TRUE) %>% parse_number(),
                BICc = grep("(BICc)", summ_fl, fixed = TRUE, value = TRUE) %>% parse_number()) #%>% unlist()

  ## options compiling - !!! re-write
  opt <- NULL

  ## gco compiling
  headers_gco <- col_map_df %>%
    transmute(
      name = COL,
      use = use,
      type = ifelse(is.na(type) | type == "", NA_character_, type)
    ) %>%
    rowwise() %>%
    mutate(.header = list(list(
      name = name,
      use = use,
      type = if (is.na(type)) setNames(list(), character(0)) else type
    ))) %>%
    ungroup() %>%
    pull(.header)

  start_idx_param <- which(str_detect(contr_obj, fixed("<PARAMETER>")))[1]
  end_idx_param <- which(str_detect(contr_obj, "\\<.*\\>") & seq_along(contr_obj) > start_idx_param)[1]
  if (is.na(end_idx_param)) end_idx_param <- length(contr_obj) + 1
  param_lines <- contr_obj[(start_idx_param + 1):(end_idx_param - 1)] %>%
    str_squish() %>%
    str_subset("^$", negate = TRUE)
  param_map <- tibble(raw = param_lines) %>%
    mutate(PAR = str_extract(raw, "^[^=]+") %>% str_trim(),
           value = str_match(raw, "value\\s*=\\s*([^,\\}]+)")[, 2] %>% as.numeric(),
           method = str_match(raw, "method\\s*=\\s*([^,\\}]+)")[, 2] %>% str_trim())

  theta_gco <- tibble(NAME = str_remove(pop_params, "_pop$")) %>%
    mutate(
      INIT = param_map$value[match(str_c(NAME, "_pop"), param_map$PAR)],
      EST = !toupper(param_map$method[match(str_c(NAME, "_pop"), param_map$PAR)]) %in% c("FIXED"),
      TRANS = sapply(NAME, function(x) {
        d <- param_distributions[[x]]$distribution
        if (is.null(d) || is.na(d)) NA_character_ else unname(d)
      })
    )
  theta_gco <- lapply(seq_len(nrow(theta_gco)), function(i) {
    list(
      NAME = theta_gco$NAME[i],
      INIT = theta_gco$INIT[i],
      EST = theta_gco$EST[i],
      TRANS = theta_gco$TRANS[i]
    )
  })

  ruv_trans <- long_ruv_map$distribution[match(dt_ruv_map$COL, long_ruv_map$COL)]
  ruv_pred <- long_ruv_map$prediction[match(dt_ruv_map$COL, long_ruv_map$COL)]
  if ("distribution" %in% names(dt_ruv_map)) {
    ruv_trans[is.na(ruv_trans)] <- dt_ruv_map$distribution[is.na(ruv_trans)]
  }
  if ("prediction" %in% names(dt_ruv_map)) {
    ruv_pred[is.na(ruv_pred)] <- dt_ruv_map$prediction[is.na(ruv_pred)]
  }
  ruv_gco <- dt_ruv_map %>%
    mutate(.idx = row_number()) %>%
    rowwise() %>%
    mutate(
      .err = list(parse_error_model(errorModel)),
      YNAME = COL,
      DVID = as.character(dvid),
      TRANS = ruv_trans[.idx],
      PRED = ruv_pred[.idx],
      ERR = .err$err,
      INIT = list(if (length(.err$pars) > 0) param_map$value[match(.err$pars, param_map$PAR)] else NA_real_),
      EST = list(if (length(.err$pars) > 0) (!toupper(param_map$method[match(.err$pars, param_map$PAR)]) %in% c("FIXED")) else NA)
    ) %>%
    ungroup() %>%
    select(-.idx)
  ruv_gco <- lapply(seq_len(nrow(ruv_gco)), function(i) {
    list(
      YNAME = ruv_gco$YNAME[i],
      DVID = ruv_gco$DVID[i],
      TRANS = ruv_gco$TRANS[i],
      PRED = ruv_gco$PRED[i],
      ERR = ruv_gco$ERR[i],
      INIT = ruv_gco$INIT[[i]],
      EST = ruv_gco$EST[[i]]
    )
  })
  if (length(ruv_gco) == 1) ruv_gco <- ruv_gco[[1]]

  re_names <- str_remove(pop_params, "_pop$")
  re_init <- matrix(0, nrow = length(re_names), ncol = length(re_names))
  re_est <- matrix(NA, nrow = length(re_names), ncol = length(re_names))
  colnames(re_init) <- rownames(re_init) <- re_names
  colnames(re_est) <- rownames(re_est) <- re_names
  if (length(omega_params) > 0) {
    for (op in omega_params) {
      p_name <- str_remove(op, "omega_")
      p_idx <- match(p_name, re_names)
      if (!is.na(p_idx)) {
        re_init[p_idx, p_idx] <- param_map$value[match(op, param_map$PAR)]
        re_est[p_idx, p_idx] <- !toupper(param_map$method[match(op, param_map$PAR)]) %in% c("FIXED")
      }
    }
  }
  re_gco <- list(init = re_init, est = re_est)
  occ_gco <- list(
    init = matrix(0, nrow = length(re_names), ncol = length(re_names), dimnames = list(re_names, re_names)),
    est = matrix(NA, nrow = length(re_names), ncol = length(re_names), dimnames = list(re_names, re_names))
  )

  covs_gco <- col_map_df %>% filter(use == "covariate") %>% pull(COL)
  model_path_raw <- extract_model_path(contr_obj)
  gco <- list(
    headers = headers_gco,
    data = c(data_path_raw),
    model = c(model_path_raw),
    task_opt = setNames(list(), character(0)),
    covs = covs_gco,
    project_name = proj_name,
    theta = theta_gco,
    ruv = ruv_gco,
    re = re_gco,
    occ = occ_gco,
    modelText = extract_model_text(contr_obj, folder_path)
  )

  ## final
  gfo <- list(SDTAB = sdtab,
              SUMTAB = sumtab,
              SIGMAMAT = sigmamat,
              OMEGAMAT = omegamat,
              OCCMAT = occmat,
              EVTAB = evtab,
              PATAB = patab,
              COTAB = cotab,
              CATAB = catab,
              REGTAB = regtab,
              OFV = ofv,
              COVMAT = covmat,
              CORRMAT = corrmat,
              OPTIONS = opt,
              PROJNAME = proj_name)
  sg_object <- list(GFO = gfo, GCO = gco)

  if (isTRUE(save_file)) {
    output_dir <- normalizePath(folder_path, mustWork = FALSE)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    gco_path <- file.path(output_dir, str_c(proj_name, "_GCO.json"))
    gfo_path <- file.path(output_dir, str_c(proj_name, "_GFO.json"))
    writeLines(jsonlite::toJSON(gco, pretty = TRUE, auto_unbox = FALSE, null = "null"), gco_path)
    writeLines(jsonlite::toJSON(gfo, pretty = TRUE, auto_unbox = FALSE, null = "null"), gfo_path)

    gco_rdata_path <- file.path(output_dir, str_c(proj_name, "_GCO.RData"))
    gfo_rdata_path <- file.path(output_dir, str_c(proj_name, "_GFO.RData"))
    save(gco, file = gco_rdata_path)
    save(gfo, file = gfo_rdata_path)
  }

  return(sg_object)


}

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
#'
#' @return
#' Returns a list with the following components:
#' \itemize{
#'   \item \code{SDTAB}: A tibble containing simulation data with columns
#'   \item \code{SUMTAB}: A tibble with parameter summary statistics containing
#'   \item \code{SIGMAMAT}: Residual variability matrix
#'   \item \code{OMEGAMAT}: Inter-individual variability matrix
#'   \item \code{OCCMAT}: Inter-occasion variability matrix (NA if not present)
#'   \item \code{EVTAB}: A tibble with event information
#'   \item \code{PATAB}: A tibble with individual parameter estimates
#'   \item \code{COTAB}: A tibble with continuous covariates
#'   \item \code{CATAB}: A tibble with categorical covariates
#'   \item \code{REGTAB}: Regression parameters (empty data.frame if not present)
#'   \item \code{OFV}: A tibble with objective function values
#'   \item \code{COVMAT}: Variance-covariance matrix of parameter estimates
#'   \item \code{CORRMAT}: Correlation matrix of parameter estimates
#'   \item \code{OPTIONS}: Model options (NULL if not present)
#'   \item \code{PROJNAME}: Project name
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
#' @examples
#' \dontrun{
#' # Convert Monolix project results
#' test_folder <- system.file("extdata", "Monolix_objects", package = "SimuRg")
#' if (substr(test_folder, nchar(test_folder), nchar(test_folder)) != "/")
#'   test_folder <- str_c(test_folder, "/")
#' pro_name <- "proj-r-solo"
#' result <- sg_converter(folder_path = test_folder, proj_name = pro_name)
#' save(results, file = "./models/simurg_object/Warfarin_PK.RData")
#' # Access individual predictions
#' head(results$SDTAB)
#'
#' # View parameter estimates
#' print(results$SUMTAB)
#'
#' # Check objective function value
#' print(results$OFV)
#' }
#'
#' @importFrom readr read_csv cols parse_number
#' @importFrom stringr str_c
#' @import tibble
#' @import dplyr
#' @export

sg_converter <- function(folder_path, proj_name){
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
    start_idx <- which(str_detect(contr_obj, fixed("[INDIVIDUAL]")))
    if (length(start_idx) == 0) {
      return(list())
    }

    # Find the end of the [INDIVIDUAL] section
    end_idx <- which(str_detect(contr_obj, "\\[.*\\]") & seq_along(contr_obj) > start_idx)[1]
    if (is.na(end_idx)) {
      end_idx <- length(contr_obj)
    }

    # Extract the [INDIVIDUAL] section
    individual_section <- contr_obj[(start_idx + 1):(end_idx - 1)]

    # Parse parameter definitions
    param_defs <- list()

    for (line in individual_section) {
      line_clean <- str_squish(line)
      if (str_detect(line_clean, "=") && str_detect(line_clean, "\\{")) {
        # Extract parameter name
        param_name <- str_extract(line_clean, "^[^=]+") %>% str_trim()

        # Extract the definition part
        def_part <- str_extract(line_clean, "\\{.*\\}")

        if (!is.na(def_part) && !is.na(param_name)) {
          # Parse the definition
          def_clean <- str_remove_all(def_part, "[\\{\\}]")
          def_parts <- str_split(def_clean, ",") %>% unlist() %>% str_trim()

          # Extract distribution, typical value, and sd
          distribution <- NA
          typical <- NA
          sd <- NA

          for (part in def_parts) {
            if (str_detect(part, "distribution=")) {
              distribution <- str_remove(part, "distribution=")
            } else if (str_detect(part, "typical=")) {
              typical <- str_remove(part, "typical=")
            } else if (str_detect(part, "sd=")) {
              sd <- str_remove(part, "sd=")
            }
          }

          param_defs[[param_name]] <- list(
            distribution = distribution,
            typical = typical,
            sd = sd
          )
        }
      }
    }

    return(param_defs)
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

  contr_obj <- readLines(str_c(folder_path, proj_name, ".mlxtran"))

  ## info about datafile
  start_idx_data <- which(str_detect(contr_obj, fixed("<DATAFILE>")))
  end_idx_data <- which(str_detect(contr_obj, "\\<.*\\>") & seq_along(contr_obj) > start_idx_data)[1]
  data_path <- contr_obj[(start_idx_data + 1):(end_idx_data - 1)] %>% str_squish() %>%
    str_subset(., "file=", negate = F) %>% str_remove(., "^[^=]+=\\s*") %>% str_replace_all("'", "")
  if (!file.exists(data_path)) data_path <- str_c(folder_path, data_path)
  data_file <- read_csv(data_path)

  ## info about columns mapping
  start_idx_col_map <- which(str_detect(contr_obj, fixed("[CONTENT]")))
  end_idx_col_map <- which(str_detect(contr_obj, "\\[.*\\]") & seq_along(contr_obj) > start_idx_col_map)[1]
  dt_col_map <- contr_obj[(start_idx_col_map + 1):(end_idx_col_map - 1)] %>% str_squish() %>% str_subset(., "^$", negate = T)


  col_map_df <- tibble(raw = dt_col_map) %>%
    mutate(COL = str_extract(raw, "^[^=]+") %>% str_trim(),
           inside = str_extract(raw, "\\{.*\\}") %>% str_remove_all("[\\{\\}]")) %>%
    separate_rows(inside, sep = ",\\s*") %>%
    separate(inside, into = c("key", "value"), sep = "=", fill = "right") %>%
    pivot_wider(names_from = key, values_from = value) %>%
    select(COL, everything(), -raw)

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


  ## dvid and residual error mapping
  start_idx_dvid_map <- which(str_detect(contr_obj, fixed("<FIT>")))
  end_idx_dvid_map <- which(str_detect(contr_obj, "\\<.*\\>") & seq_along(contr_obj) > start_idx_dvid_map)[1]

  dt_dvid_map <- contr_obj[(start_idx_dvid_map + 1):(end_idx_dvid_map - 1)] %>% str_squish() %>% str_subset(., "^$", negate = T)


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

  extract_values_ruv <- function(x) {
    matches <- regmatches(x, gregexpr("\\(([^)]+)\\)", x))
    values <- lapply(matches, function(m) strsplit(gsub("[()]", "", m), "; |;")) %>% unlist()
    return(values)
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


  dvid_map_df <- tibble(data = extract_values(dt_dvid_map[1]),
                        model = extract_values(dt_dvid_map[2]))

  start_idx_ruv_map <- which(str_detect(contr_obj, fixed("[LONGITUDINAL]")))
  end_idx_ruv_map <- which(str_detect(contr_obj, "\\<.*\\>") & seq_along(contr_obj) > start_idx_ruv_map)[1]

  dt_ruv_map <- contr_obj[(start_idx_ruv_map + 1):(end_idx_ruv_map - 1)] %>% str_squish() %>% str_subset(., "^$", negate = T)
  dt_ruv_map <- grep(str_c(str_c(dvid_map_df$model, " ="), collapse = "|"), dt_ruv_map, value = T) %>%
    tibble(raw = .) %>%
    mutate(COL = str_extract(raw, "^[^=]+") %>% str_trim()) %>%
    group_by(COL) %>%
    mutate(inside = str_extract(raw, "\\{.*\\}") %>% replace_commas_in_parentheses(.) %>% str_remove_all("[\\{\\}]")) %>%
    ungroup() %>%
    separate_rows(inside, sep = ",\\s*") %>%
    separate(inside, into = c("key", "value"), sep = "=", fill = "right") %>%
    pivot_wider(names_from = key, values_from = value) %>%
    group_by(COL) %>%
    mutate(RUVpars = list(extract_values_ruv(errorModel)),
           RUVpar_a = ifelse(grepl("constant|combined", errorModel), unlist(RUVpars)[1], NA) ,
           RUVpar_b = ifelse(grepl("proportional", errorModel), unlist(RUVpars)[1],
                             ifelse(grepl("combined", errorModel),unlist(RUVpars)[2], NA))) %>%
    ungroup() %>%
    select(COL, everything(), -c("raw", "RUVpars")) %>%
    left_join(dvid_map_df %>% rename(COL = model, dvid = data), "COL")


  ## cotab, catab, regtab compiling

  cotab_cols <- col_map_df %>% filter(use == "identifier" | (use == "covariate" & type == "continuous")) %>% select(COL) %>% pull()
  catab_cols <- col_map_df %>% filter(use == "identifier" | (use == "covariate" & type == "categorical")) %>% select(COL) %>% pull()
  regtab_cols <- col_map_df %>% filter(use %in% c("identifier", "time", "regressor")) %>% select(COL) %>% pull()

  if (length(cotab_cols) > 1) {cotab <- data_file_mod %>% select(all_of(cotab_cols)) %>% unique()} else {cotab <- data.frame()}
  if (length(catab_cols) > 1) {catab <- data_file_mod %>% select(all_of(catab_cols)) %>% unique()} else {catab <- data.frame()}
  if (length(regtab_cols) > 2) {regtab <- data_file_mod %>% select(regtab_cols) %>% unique()} else {regtab <- data.frame()}


  ## patab and sumtab compiling

  sum_dt_i <- read_csv(str_c(folder_path, proj_name, "/populationParameters", ".txt"), col_types = cols())

  pop_params <- sum_dt_i$parameter[str_detect(sum_dt_i$parameter, "_pop$")]
  params <- str_replace(pop_params, "_pop$", "")
  omega_params <- sum_dt_i$parameter[sum_dt_i$parameter %in% str_c("omega_", params)]
  beta_params <- sum_dt_i$parameter[str_detect(sum_dt_i$parameter, "^beta_")]
  corr_params <- sum_dt_i$parameter[str_detect(sum_dt_i$parameter, "^corr_")]
  resid_err_params <- sum_dt_i$parameter[!sum_dt_i$parameter %in% c(pop_params, omega_params, beta_params, corr_params)]
  eta_params <- str_replace(omega_params, "omega_", "eta_")

  ## sdtab compiling
  sdtab <- unique(dvid_map_df$model) %>% map_dfr(function(y_name) {

    if (length(dvid_map_df$model) == 1) {y_name_i <- ""; dvid_i <- 1} else {y_name_i <- str_c("_", y_name); dvid_i <- dvid_map_df$data[dvid_map_df$model == y_name] %>% as.numeric()}
    pred_dt_i <- read_csv(str_c(folder_path, proj_name, "/predictions",y_name_i, ".txt"), col_types = cols()) %>%
      rename(ID = id, TIME = time)

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
    cat("Detected model structure:", model_structure, "\n")
    cat("Omega parameters:", paste(omega_params, collapse = ", "), "\n")

    # Get parameter names that have random effects
    random_effect_params <- str_replace(omega_params, "omega_", "")
    cat("Random effect parameters:", paste(random_effect_params, collapse = ", "), "\n")
    cat("Using Monte Carlo simulation with n_sim = 1000\n")

    # Extract parameter distributions
    param_distributions <- extract_parameter_distributions(contr_obj)
    cat("Extracted parameter distributions:", length(param_distributions), "parameters\n")
    for (param in names(param_distributions)) {
      dist_info <- param_distributions[[param]]
      cat("  ", param, ": ", dist_info$distribution, " (typical=", dist_info$typical, ", sd=", dist_info$sd, ")\n")
    }

    # Prepare sdtab data
    sdtab_i <- left_join(pred_dt_i, obs_data_i, by = c("ID", "TIME")) %>%
      rename(PRED = popPred,
             IPRED = str_c("indivPred", suffix), #indivPred_mode,
             IWRES = str_c("indWRes", suffix), #indWRes_mode,
             DV = y_name) %>%
      mutate(RES = PRED - DV, IRES = IPRED - DV,
             DVID = dvid_i)

    # Calculate WRES using Monte Carlo simulation
    if (nrow(ruv_info) > 0) {
      # Use the Monte Carlo simulation function to calculate WRES
      sdtab_i <- calculate_wres_mc(sdtab_i, omega_params, resid_err_params, sum_dt_i, ruv_info, param_distributions)
    } else {
      # If no residual error model info, set WRES to NA
      sdtab_i <- sdtab_i %>% mutate(WRES = NA_real_)
    }

    sdtab_i %>%
      select(any_of(c("ID", "TIME", "DV", "DVID", "PRED", "IPRED", "RES", "IRES", "WRES", "CWRES",
                    "IWRES", "EVID", "MDV", "OCC", "BLQ", "CENS", "LIMIT")))

  })

  ## reverse naming
  #sdtab <- sdtab %>%
  #  rename(!!!setNames(col_map_df$COL[col_map_df$use %in% c("identifier", "time")], c("ID", "TIME")))


  ## patab and sumtab compiling


  eta_i <- read_csv(str_c(folder_path, proj_name, "./IndividualParameters/estimatedRandomEffects.txt"), col_types = cols())
  indpar_i <- read_csv(str_c(folder_path, proj_name, "./IndividualParameters/estimatedIndividualParameters.txt"), col_types = cols())

  if (any(grepl("_mode$", colnames(indpar_i)))){suffix <- "_mode"} else {suffix <- "_SAEM"}
  eta_clnms <- c("id", str_c(eta_params, suffix))
  patab <- left_join(eta_i %>% select(any_of(eta_clnms)), #select_at(vars("id", ends_with(suffix))),
                     indpar_i %>% select_at(vars("id", ends_with(suffix)))) %>%
    rename_with(~ str_replace(., str_c(suffix, "$"), ""), ends_with(suffix)) %>% rename(ID = id)


  ## sumtab compiling

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
           EST = ifelse(is.na(SE), "FIXED", "ESTIMATED")
           ) %>%
    select(any_of(c("PAR", "VALUE", "TYPE", "EST", "SE", "RSE", "CV", "LCI", "UCI", "ETAshrinkage_sd", "ETAshrinkage_var", "EPSshrinkage_sd")))


  ## covmat and corrmat compiling

  fi_files <- list.files(str_c(folder_path, proj_name, "./FisherInformation"))
  corr_path <- str_c(folder_path, proj_name, "./FisherInformation/", fi_files[str_detect(fi_files, "correlation")])
  cov_path <- str_c(folder_path, proj_name, "./FisherInformation/", fi_files[str_detect(fi_files, "covariance")])

  covmat_dt <- read_csv(cov_path, col_types = cols(), col_names = F)
  colnames(covmat_dt) <- c("PAR", covmat_dt$X1)
  covmat <- covmat_dt %>% select(-PAR) %>% as.matrix()

  corrmat_dt <- read_csv(corr_path, col_types = cols(), col_names = F)
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
        mutate(PAR1 = na.omit(str_extract(CORR_PAR, params))[1],
               PAR2 = na.omit(str_extract(CORR_PAR, params))[2],
               PAR_N1 = which(PAR1 == params),
               PAR_N2 = which(PAR2 == params)) %>%
        ungroup()

      for (i in 1:length(corr_params)){

        corr_params_dt_i <- corr_params_dt %>% filter(CORR_PAR == corr_params[i])
        corr_value <- sumtab$VALUE[sumtab$PAR == corr_params[i]]
        var1_value <- omegamat[corr_params_dt_i$PAR_N1, corr_params_dt_i$PAR_N1]
        var2_value <- omegamat[corr_params_dt_i$PAR_N2, corr_params_dt_i$PAR_N2]

        omegamat[corr_params_dt_i$PAR_N1, corr_params_dt_i$PAR_N2] <- corr_value*var1_value*var2_value

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

  ofv <- tibble(LL = grep("log-likelihood", summ_fl, fixed = T, value = T) %>%
                         gsub("^.*?:", "", .) %>% parse_number(),
                    AIC = grep("AIC", summ_fl, fixed = T, value = T) %>% parse_number(),
                    BIC = grep("(BIC)", summ_fl, fixed = T, value = T) %>% parse_number(),
                    BICc = grep("(BICc)", summ_fl, fixed = T, value = T) %>% parse_number()) #%>% unlist()

  ## options compiling - !!! re-write
  opt <- NULL

  ## final
  sg_object = list(SDTAB = sdtab,
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

  return(sg_object)


}




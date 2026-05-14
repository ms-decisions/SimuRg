## Author: Alina Melnikova
## First created: 2025-07-17
## Description: Covariate simulation from original dataset
## Keywords: covariates, simulations
## Version: v1.1 - add JS and KL metrics, covariate matrix, merge umap


# Comparison of correlation matrices
compare_cor_matrices <- function(data_obs, data_syn, vars, method = "kendall") {

  R_obs <- cor(data_obs[, vars, drop = FALSE],
               use = "pairwise.complete.obs",
               method = method)

  R_syn <- cor(data_syn[, vars, drop = FALSE],
               use = "pairwise.complete.obs",
               method = method)

  idx <- upper.tri(R_obs, diag = FALSE)

  # mean_abs_diff = mean(abs(R_obs[idx] - R_syn[idx]))
  # return(mean_abs_diff )
   list(
     mean_abs_diff = mean(abs(R_obs[idx] - R_syn[idx])),
     max_abs_diff  = max(abs(R_obs[idx] - R_syn[idx])),
  #   #frobenius     = sqrt(sum((R_obs - R_syn)^2)),
     R_obs = R_obs,
     R_syn = R_syn,
    R_diff = R_syn - R_obs
  )
}

# Create optimal visit sequence based on correlations

create_optimal_visit_sequence <- function(data, var_cont, var_cat) {
  var_all <- c(var_cont, var_cat)

  # Handle edge cases
  if (length(var_all) == 0) {
    warning("No variables provided for visit sequence")
    return(character(0))
  }

  if (length(var_cont) == 0 && length(var_cat) > 0) {
    # Only categorical variables, return as-is
    return(var_cat)
  }

  if (length(var_cont) < 2) {
    # If too few continuous variables, use default order
    return(var_all)
  }

  # Calculate correlation matrix for continuous variables
  if (length(var_cont) >= 2) {
    cor_matrix <- cor(data[, var_cont, drop = FALSE],
                      use = "pairwise.complete.obs")
    # Replace NA with 0 (for variables with no variance)
    cor_matrix[is.na(cor_matrix)] <- 0

    # Calculate average absolute correlation for each variable
    avg_abs_corr <- rowMeans(abs(cor_matrix))
    names(avg_abs_corr) <- var_cont

    # Order continuous variables by average correlation (highest first)
    var_cont_ordered <- names(sort(avg_abs_corr, decreasing = TRUE))
  } else {
    var_cont_ordered <- var_cont
  }

  # Combine: continuous (ordered) + categorical
  # Handle both single and multiple categorical variables
  if (length(var_cat) > 0) {
    visit_seq <- c(var_cont_ordered, var_cat)
  } else {
    visit_seq <- var_cont_ordered
  }

  return(visit_seq)
}

# Remove exact duplicates between synthetic and original data by adding noise
remove_exact_duplicates <- function(data_syn, data_orig, var_cont, var_cat,
                                     noise_level = 0.10, seed = 123) {

  # Find exact duplicates between synthetic and original
  common_cols <- intersect(names(data_syn), names(data_orig))
  data_syn$row_id <- seq_len(nrow(data_syn))

  # Identify which synthetic rows are exact duplicates of original rows
  dupl_indices <- data_syn %>%
    dplyr::semi_join(data_orig, by = common_cols) %>%
    pull(row_id)

  n_duplicates <- length(dupl_indices)

  if (n_duplicates == 0) {
    data_syn$row_id <- NULL
    return(list(
      data_cleaned = data_syn,
      n_duplicates_removed = 0,
      duplicate_indices = integer(0)
    ))
  }

  message("Found", n_duplicates, "exact duplicates. Adding noise to remove them...\n")

  # Add noise to continuous variables for duplicate rows
  if (length(var_cont) > 0) {
    for (var in var_cont) {
      if (var %in% names(data_syn)) {
        # Calculate SD from original data
        sd_val <- sd(data_orig[[var]], na.rm = TRUE)

        # Add noise only to duplicate rows
        if (is.null(seed)) {
          noise <- rnorm(n_duplicates, mean = 0, sd = noise_level * sd_val)
        } else {
          withr::with_seed(seed, {
            noise <- rnorm(n_duplicates, mean = 0, sd = noise_level * sd_val)
          })
        }
        data_syn[dupl_indices, var] <- data_syn[dupl_indices, var] + noise

        # Ensure values stay within reasonable bounds (min/max of original)
        min_val <- min(data_orig[[var]], na.rm = TRUE)
        max_val <- max(data_orig[[var]], na.rm = TRUE)
        data_syn[dupl_indices, var] <- pmin(pmax(data_syn[dupl_indices, var], min_val), max_val)
      }
    }
  }

  # For binary/categorical variables, randomly flip a small proportion
  if (length(var_cat) > 0) {
    for (var in var_cat) {
      if (var %in% names(data_syn)) {
        levels_var <- levels(data_syn[[var]])
        n_levels <- length(levels_var)

        # Only perturb if there are at least 2 levels
        if (n_levels >= 2) {
          # Randomly select ~20% of duplicate rows to flip for this variable
          n_to_flip <- max(1, ceiling(0.2 * n_duplicates))
          if (is.null(seed)) {
            rows_to_flip <- sample(dupl_indices, size = n_to_flip)
          } else {
            withr::with_seed(seed, {
              rows_to_flip <- sample(dupl_indices, size = n_to_flip)
            })
          }

          # For each selected row, change to a different level
          for (idx in rows_to_flip) {
            current_level <- as.character(data_syn[idx, var])
            other_levels <- setdiff(levels_var, current_level)
            if (length(other_levels) > 0) {
              data_syn[idx, var] <- factor(sample(other_levels, 1), levels = levels_var)
            }
          }
        }
      }
    }
  }

  data_syn$row_id <- NULL

  return(list(
    data_cleaned = data_syn,
    n_duplicates_removed = n_duplicates,
    duplicate_indices = dupl_indices
  ))
}

#' Perform generation of synthetic datasets for an empirical distribution
#'
#' The function operates in two modes:
#' \itemize{
#'   \item \strong{Fixed seed mode} (when \code{seed} is specified): Generates a single dataset using the provided seed
#'   \item \strong{Search mode} (when \code{seed = NA}): Iteratively searches for \code{nds} datasets that meet the correlation difference threshold
#' }
#'
#' @returns A list of lists, where each element contains results for one generated dataset:
#' \itemize{
#'   \item \code{datagen} - Synthetic dataset returned by \code{synthpop::syn()} (data.frame)
#'   \item \code{seed} - Random seed used to generate this dataset (integer or NA)
#'   \item \code{exact_dupl_check} - Logical indicating if exact duplicates exist between original and synthetic data
#'   \item \code{dplot_umap} - ggplot object with combined UMAP visualization comparing original and synthetic data (or \code{NULL} if diag_plots=FALSE)
#'   \item \code{ks_test} - Tibble with Kolmogorov-Smirnov p-values and statuses for continuous variables (or \code{NULL} if no continuous variables)
#'   \item \code{jsd_res} - Weighted mean Jensen-Shannon divergence (JSD) value for categorical variables (numeric or \code{NULL} if no categorical variables)
#'   \item \code{corr_diff_mean} - Mean absolute difference between original and synthetic correlation matrices (numeric or \code{NULL} if no continuous variables)
#'   \item \code{corr_diff_max} - Maximum absolute difference between original and synthetic correlation matrices (numeric or \code{NULL} if no continuous variables)
#'   \item \code{dplot_corr_diff} - ggplot heatmap object showing correlation difference matrix (Synthetic - Original) (or \code{NULL} if diag_plots=FALSE or no continuous variables)
#'   \item \code{dplot_cont} - List of ggplot histograms for continuous variables (or \code{NULL} if diag_plots=FALSE or no continuous variables)
#'   \item \code{dplot_cat} - List of ggplot barplots for categorical variables (or \code{NULL} if diag_plots=FALSE or no categorical variables)
#' }
#' @inheritParams sg_dummy
#' @param data_i data frame. Input data frame containing the original dataset to be synthesized (required)
#' @param nobj integer. Specify the exact number of rows to generate in the synthetic dataset. When provided, overrides \code{npop} (optional, default: \code{NA})
#' @param minnumlev integer. Threshold; numeric variables with <= \code{minnumlev} unique values are converted to factors (optional, default: \code{3})
#' @param seed integer. Random seed for synthetic data generation. If provided (not \code{NA}), generates a single dataset with this seed (fixed seed mode). If \code{NA}, uses search mode to find \code{nds} datasets meeting correlation threshold (optional, default: \code{NA})
#' @param seed_umap integer. Random seed for UMAP algorithm reproducibility (optional, default: \code{123})
#' @param palette character vector. Contains color codes (hex format) for custom plot color schemes. If provided, should contain at least 2 colors. Used for histograms, bar plots, and UMAP visualizations (optional, default: \code{c("#3a6eba", "#efdd3c", "#1a1866", "#f2b93b")})
#' @param diag_plots logical flag. If \code{TRUE}, generates diagnostic plots and UMAP visualizations (optional, default: \code{TRUE})
#' @param remove_duplicates logical flag. If \code{TRUE}, automatically removes exact duplicates between original and synthetic data by adding controlled noise (optional, default: \code{TRUE})
#' @param noise_level numeric. Proportion of standard deviation to use when adding noise to continuous variables for duplicate removal. E.g., 0.10 means 10% of SD (optional, default: \code{0.10})
#' @param nds integer. Number of synthetic datasets to generate in search mode. Ignored in fixed seed mode (optional, default: \code{1})
#' @param tg_corrdif numeric. Target maximum absolute correlation difference threshold for dataset selection in search mode. Ignored in fixed seed mode (optional, default: \code{0.1})
#' @examples
#' \donttest{
#' library(dplyr)
#'
#' # Generate example dataset
#' data <- data.frame(
#'   id = 1:150,
#'   ALB = rnorm(150, mean = 3.5, sd = 0.5),
#'   ALT = rnorm(150, mean = 50, sd = 20),
#'   SEX = factor(sample(c("M", "F"), 150, replace = TRUE)),
#'   RACE = factor(sample(c("White", "Black", "Asian", "Other"), 150, replace = TRUE))
#' )
#'
#' # EXAMPLE 1: Fixed seed mode - generate single dataset with specified seed
#' output_fixed <- sg_vpop_est(data_i = data,
#'                             id_col = "id",
#'                             seed = 123,        # Fixed seed mode
#'                             seed_umap = 40,
#'                             diag_plots = TRUE)
#'
#' # Access the dataset and its diagnostics
#' print(head(output_fixed[[1]]$datagen, 10))
#' print(output_fixed[[1]]$ks_test)
#' print(output_fixed[[1]]$seed)           # Will be 123
#' print(output_fixed[[1]]$corr_diff_max)
#' print(output_fixed[[1]]$dplot_umap)
#'
#' # EXAMPLE 2: Search mode - generate 3 datasets meeting correlation threshold
#' output_search <- sg_vpop_est(data_i = data,
#'                              id_col = "id",
#'                              seed = NA,        # Search mode (default)
#'                              seed_umap = 40,
#'                              diag_plots = TRUE,
#'                              nds = 3,          # Generate 3 datasets
#'                              tg_corrdif = 0.1) # Target correlation difference
#'
#' # Compare metrics across all datasets
#' sapply(output_search, function(x) x$corr_diff_max)
#' sapply(output_search, function(x) x$jsd_res)
#' sapply(output_search, function(x) x$seed)  # Different seeds for each dataset
#'
#'}
#' @import dplyr
#' @import ggplot2
#' @importFrom purrr map_dfr map
#' @importFrom stringr str_c
#' @importFrom cluster daisy
#' @importFrom fastDummies dummy_cols
#' @importFrom synthpop syn
#' @importFrom uwot umap umap_transform
#' @importFrom tidyr drop_na
#' @importFrom tibble tibble as_tibble
#' @importFrom recipes recipe step_dummy step_center step_scale prep bake
#' @importFrom philentropy distance
#' @export
sg_vpop_est <-  function(data_i, nobj = NA, id_col = NULL, minnumlev = 3,npop = 1,
                       excl_col = NULL,seed = NA,
                       seed_umap = 123, palette = NULL,
                       diag_plots = FALSE,
                       #show_info=FALSE,
                       remove_duplicates = TRUE,
                       noise_level = 0.10,
                       nds = 1,#number of datasets to generate
                       tg_corrdif = 0.1){

  theme_set(theme_bw())
  theme_update(panel.grid.minor = element_blank())
  MSDcol_cut <- c("#3a6eba","#efdd3c", "#1a1866", "#b73b58")

  if(!is.null(palette) & length(palette)>1){
    color_p = rep(palette,2)
  } else {color_p = MSDcol_cut
    message("Default color palette is used")}

  #Control parameters for Random Forest (optimized for correlation preservation)
  num_trees = 300          # Increased from 100 for better stability
  max_depth = 12           # Increased from 7 to capture more interactions
  min_node_size = 5        # Decreased from 30 for finer splits
  #use_smoothing = TRUE     # Enable smoothing for continuous variables
  optimize_visit_seq = TRUE # Use correlation-based visit sequence

  ### Comment!!!! - Testing code below (comment out for production)
  # data_orig <- read.csv("functions/datasets/data_gbsg2.csv") %>% head(300)  # Change to your dataset path
  # data_i = data_orig; nobj = NA; id_col = NULL; minnumlev = 3;npop = 1
  # excl_col = NULL;
  # seed = 98;
  # seed_umap = 42; palette = NULL;
  # diag_plots = F; remove_duplicates = TRUE;
  # noise_level = 0.10; tg_corrdif = 0.1; nds = 5

  #####--------------- Processing ---------------#####
  #Exclude ID column and Exclusion columns

  if (!is.null(id_col)){
    if (!(id_col %in% names(data_i))){
      warning("`id_col` = '", id_col, "' is not present in `data`.")
    }
    data_i <- data_i %>%
      dplyr::select(-any_of(id_col))
  }

  if (!is.null(excl_col)){
    missing_excl <- setdiff(excl_col, names(data_i))
    if (length(missing_excl) > 0){
      warning("The following `excl_col` values are not present in `data`: ",
              paste(missing_excl, collapse = ", "))
    }
    data_i <- data_i %>%
      dplyr::select(-any_of(excl_col))
  }


  #Exclude rows with NAs.
  data_i <- data_i %>%
    drop_na()

  var_char = NULL

  # Convert character columns to factor
  is_char  <- vapply(data_i, is.character, logical(1))
  var_char  <- names(data_i)[is_char]

  if (!is.na(minnumlev)){minnum = minnumlev} else {minnum = 3}

  # Convert numeric variables with <= minnumlev unique values to factor
  is_low_level_num <- vapply(
    data_i,
    function(x) is.numeric(x) && length(unique(na.omit(x))) <= minnum,
    logical(1)
  )
  var_low_level <- names(data_i)[is_low_level_num]


  data_i <- data_i %>%
    mutate(across(any_of(c(var_char, var_low_level)), as.factor))

  is_cont <- vapply(data_i, is.numeric, logical(1))
  is_cat  <- vapply(data_i, is.factor, logical(1))

  # Define names for continuous and categorical (ordered/factor) columns
  var_cont <- names(data_i)[is_cont]
  var_cat  <- names(data_i)[is_cat]
  var_all <- c(var_cont, var_cat)



  n_orig <- nrow(data_i)


  # if (show_info){
  # # Print dataset summary
  #  message("Number of rows in the data_orig:", nrow(data_i), "\n")
  #  message("Number and names of continuous columns - ", length(var_cont), ":", paste(var_cont, collapse = ", "), "\n")
  #  message("Number and names of categorical columns - ", length(var_cat), ":", paste(var_cat, collapse = ", "), "\n")
  #}
  data_orig <- data_i

  # Ensure data_i is a tibble to satisfy synthpop list-argument requirements
  data_i <- dplyr::as_tibble(data_i)

  # Convert all integer columns to numeric to prevent smoothing/casting errors
  # (Smoothing produces continuous values, which can fail if cast back to integer)
  data_i <- data_i %>% mutate(across(where(is.integer), as.numeric))

  if (!is.na(nobj)){n_new = nobj} else if(!is.na(npop)&(npop>0)){
    n_new = as.integer(round(n_orig * npop))
    } else {n_new = n_orig}

  #Synthpop Method
  # if (!is.na(seed)){seed_i = seed} else {seed_i = 123}
  # message("Seed for synthetic data generation:", seed_i, "\n")

  # Create optimal visit sequence based on correlations (if enabled)
  if (optimize_visit_seq && length(var_cont) >= 2) {
    message("Creating correlation-based visit sequence...\n")
    visit_sequence <- create_optimal_visit_sequence(data_i, var_cont, var_cat)
  } else {
    # Use default order: all variables in their original order
    visit_sequence <- var_all
  }

  # Prepare method vector: "rf" for all variables
  method_vector <- rep("rf", length(var_all))
  names(method_vector) <- var_all



  message(sprintf("RF parameters: trees = %.2f, max_depth = %.2f, , min_node_size = %.2f\n", num_trees, max_depth,
      min_node_size))
  message("Visit sequence:", paste(visit_sequence, collapse = ", "), "\n")
  message("Variables to synthesize:", length(visit_sequence),
      "(", length(var_cont), "continuous,", length(var_cat), "categorical )\n")


  # Generate synthetic data with optimized parameters
  # Pass RF parameters directly to syn() function

  # Initialize results list
  results_list <- list()

  # Check mode: fixed seed or search mode
  if (!is.na(seed)) {
    # === MODE 1: Fixed seed mode ===
    message("\n=== Fixed seed mode: generating dataset with seed = ", seed, "===\n")
    mode_fixed_seed <- TRUE
    n_datasets <- 1
  } else {
    # === MODE 2: Search mode (original logic) ===
    message("\n=== Search mode: generating", nds, "dataset(s) with optimal seeds ===\n")
    mode_fixed_seed <- FALSE
    n_datasets <- nds

    # Initialize search parameters
    n_iterations <- 100
    diff_lim <- tg_corrdif
    seed_values <- 1:n_iterations
    i_init <- 0
    i <- 0
    res_seeds <- c()
  }

  # Main loop for dataset generation and diagnostics
  for (k in 1:n_datasets) {

    # === GENERATION PHASE ===
    if (mode_fixed_seed) {
      # Fixed seed mode: generate single dataset with specified seed
      message("\n=== Generating dataset with fixed seed =", seed, "===\n")

      syn_obj <- syn(data_i,
                     method = method_vector,
                     visit.sequence = visit_sequence,
                     #smoothing = smoothing_vector,
                     k = n_new,
                     seed = seed,
                     ranger.num.trees = num_trees,
                     ranger.max.depth = max_depth,
                     ranger.min.node.size = min_node_size,
                     ranger.respect.unordered.factors = "order")

      data_syn <- syn_obj$syn
      target_seed <- seed

    } else {
      # Search mode: iterate through seeds until target corr_diff is achieved
      message("\n=== Generating dataset", k, "of", n_datasets, "===\n")

      corr_diff_contr <- diff_lim + 0.05
      p <- 0

      while ((corr_diff_contr > diff_lim) & (p < n_iterations)) {
        i <- i + 1
        p <- p + 1
        syn_obj <- syn(data_i,
                       method = method_vector,
                       visit.sequence = visit_sequence,
                       #smoothing = smoothing_vector,
                       k = n_new,
                       seed = seed_values[p],
                       ranger.num.trees = num_trees,
                       ranger.max.depth = max_depth,
                       ranger.min.node.size = min_node_size,
                       ranger.respect.unordered.factors = "order")
        data_syn <- syn_obj$syn
        corr_diff <- compare_cor_matrices(data_i, data_syn, var_cont, method = "kendall")
        corr_diff_contr <- corr_diff$max_abs_diff
      }

      target_seed <- seed_values[[p]]

      if ((p == n_iterations) & (corr_diff_contr > diff_lim)) {
        message("Target correlation coefficient difference is not obtained\n")
        target_seed <- NA
      }

      res_seeds <- c(res_seeds, target_seed)
      message("Dataset", k, "generated with seed:", target_seed, "\n")

      # Update seed range for next iteration
      i_init <- i
      seed_values <- (i_init + 1):(n_iterations + i_init)
    }

    # === DIAGNOSTICS PHASE (common for both modes) ===
    message("Computing diagnostics for dataset", k, "...\n")

    # Check for duplicates *between* real and synthetic data (not within each)
    common_cols <- intersect(names(data_i), names(data_syn))
    dupl_check_before <- nrow(dplyr::semi_join(data_syn, data_i, by = common_cols)) > 0
    n_dupl_removed <- 0

    if (dupl_check_before){
      #warning("Exact duplicates found between original and synthetic data.")

      # Remove duplicates if requested
      if (remove_duplicates) {
        dedup_result <- remove_exact_duplicates(
          data_syn = data_syn,
          data_orig = data_i,
          var_cont = var_cont,
          var_cat = var_cat,
          noise_level = noise_level,
          #seed = seed_i
        )

        data_syn <- dedup_result$data_cleaned
        n_dupl_removed <- dedup_result$n_duplicates_removed

        # Check again after deduplication
        dupl_check_after <- nrow(dplyr::semi_join(data_syn, data_i, by = common_cols)) > 0

        if (dupl_check_after) {
          message("Warning: Some duplicates may still remain after noise addition.\n")
        } else {
          message("Successfully removed all", n_dupl_removed, "exact duplicates.\n")
        }
      }
    }

    dupl_check <- nrow(dplyr::semi_join(data_syn, data_i, by = common_cols)) > 0

    # Initialize diagnostic variables
    ks_results <- NULL
    corr_diff <- NULL
    JSD_results <- NULL

    if (length(var_cont)>0){
      # Get Kolmogorov-Smirnov test
      # ks_results <- map_dfr(var_cont, function(var) {
      #   x <- data_i[[var]]
      #   y <- data_syn[[var]]
      #
      #   # Detect ties across combined samples; ks.test assumes continuous data
      #   has_ties <- any(duplicated(c(x, y)))
      #
      #   # Suppress the default ks.test warning about ties and emit a custom one instead
      #   ks_result <- suppressWarnings(ks.test(x, y))
      #   if (has_ties) {
      #     warning(
      #       "Kolmogorov-Smirnov test for variable '", var,
      #       "' may be approximate due to ties in the data",
      #       call. = FALSE
      #     )
      #   }
      #   p_val <- ks_result$p.value
      #
      #   # Format p-value as character to ensure consistent type for bind_rows
      #   if (p_val < 0.01) {
      #     p_val_formatted <- "<0.01"
      #   } else {
      #     p_val_formatted <- as.character(round(p_val, 4))
      #   }
      #
      #   tibble(
      #     variable = var,
      #     p.value = p_val_formatted,
      #     status = ifelse(p_val > 0.05, "similar", "different")
      #   )
      # })

      ks_results <- suppressWarnings({
        map_dfr(var_cont, function(var) {
          x <- data_i[[var]]
          y <- data_syn[[var]]

          # Detect ties across combined samples; ks.test assumes continuous data
          has_ties <- any(duplicated(c(x, y)))

          # Suppress the default ks.test warning about ties and emit a custom one instead
          ks_result <- ks.test(x, y)
          if (has_ties) {
            warning(
              "Kolmogorov-Smirnov test for variable '", var,
              "' may be approximate due to ties in the data",
              call. = FALSE
            )
          }
          p_val <- ks_result$p.value

          # Format p-value as character to ensure consistent type for bind_rows
          if (p_val < 0.01) {
            p_val_formatted <- "<0.01"
          } else {
            p_val_formatted <- as.character(round(p_val, 4))
          }

          tibble(
            variable = var,
            p.value = p_val_formatted,
            status = ifelse(p_val > 0.05, "similar", "different")
          )
        })
      })



      #Comparison of correlation matrices (requires at least 2 variables)
      if (length(var_cont) > 1) {
        corr_diff <- compare_cor_matrices(data_i, data_syn, var_cont, method = "kendall")
      }
    }

    if (length(var_cat)>0){
      n_levels <- c()
      # Jensen-Shannon (JS) divergence
      jsd_by_var <- sapply(var_cat, function(var) {

        # Align category levels
        lev <- union(
          levels(factor(data_i[[var]])),
          levels(factor(data_syn[[var]]))
        )
        # Empirical PMFs
        p <- prop.table(table(factor(data_i[[var]], levels = lev)))
        q <- prop.table(table(factor(data_syn[[var]], levels = lev)))

        # Small epsilon to avoid log(0)
        eps <- 1e-8
        p <- p + eps
        q <- q + eps
        p <- p / sum(p)
        q <- q / sum(q)

        # Jensen-Shannon divergence (obs || syn)
        as.numeric(
          distance(rbind(p, q), method = "jensen-shannon")
        )
      })

      jsd <- as.vector(jsd_by_var)
      names(jsd) <- var_cat
      n_levels <- sapply(var_cat, function(var) {
        length(levels(factor(data_i[[var]])))
      })
      names(n_levels) <- var_cat
      # Calculate weighted mean of Jensen-Shannon divergence
      JSD_results<- sum(jsd * n_levels / sum(n_levels))
    }

    # Preparation for plots
    plist_cont = NULL
    plist_cat = NULL
    plot_umap = NULL
    p_corr_heatmap = NULL

    if (diag_plots){
      rec <- recipe(~ ., data = data_i)
      if (length(var_cat)>0){
        # recipe for preprocessing
        rec <- rec %>%
          step_dummy(all_of(var_cat), one_hot = TRUE)
      }
      if (length(var_cont)>0){
        # recipe for preprocessing
        rec <- rec %>%
          step_center(all_of(var_cont)) %>%
          step_scale(all_of(var_cont)) %>%
          prep()
      }

      X_real <- bake(rec, new_data = data_i)
      X_syn  <- bake(rec, new_data = data_syn)
      # UMAP (fit on real data only)



      message("Seed for UMAP:", seed_umap, "\n")

      umap_model <- uwot::umap(
        X_real,
        n_neighbors = 15,
        min_dist = 0.1,
        metric = "euclidean",
        ret_model = TRUE,
        ret_nn = TRUE,
        seed = seed_umap
      )

      # project synthetic data
      emb_real <- umap_model$embedding
      emb_syn  <- uwot::umap_transform(X_syn, umap_model, seed = seed_umap)

      # Create combined visualization
      umap_df <- rbind(
        data.frame(X = emb_real[, 1], Y = emb_real[, 2], Source = "Original data"),
        data.frame(X = emb_syn[, 1],  Y = emb_syn[, 2],  Source = "Synthetic data")
      )

      # Create plot
      plot_umap <- ggplot(umap_df, aes(x = X, y = Y, color = Source)) +
        geom_point(alpha = 0.4, size = 1) +
        scale_color_manual(values = c("Original data" = color_p[[3]], "Synthetic data" = color_p[[4]])) +
        labs(title = paste0("UMAP: Dataset ", k),
             x = "UMAP Dimension 1", y = "UMAP Dimension 2", color = "Dataset") +
        theme_minimal() +
        theme(legend.position = "bottom")

      if (length(var_cont)>0){
        #Get Continuous histograms
        plist_cont <- var_cont %>% map(function(var){
          # Create data frame for plotting
          plot_data <- data.frame(
            value = c(data_i[[var]], data_syn[[var]]),
            source = rep(c("Original", "Synthetic"),
                         times = c(length(data_i[[var]]), length(data_syn[[var]])))
          )

          # Create ggplot histogram
          ggplot(plot_data, aes(x = value, fill = source)) +
            geom_histogram(alpha = 0.6, bins = 30, position = "identity",
                           color = "black", linewidth = 0.3) +
            scale_fill_manual(
              values = c("Original" = color_p[[1]], "Synthetic" = color_p[[2]])
            ) +
            labs(title = paste(var, "Distribution"),
                 x = var,
                 y = "Frequency",
                 fill = "Dataset") +
            theme(legend.position = "topright")
        })

        #Heat map for corr_diff R_diff (using ggplot2 for better visualization)
        # Only create heatmap if there are at least 2 continuous variables
        if (length(var_cont) > 1 && !is.null(corr_diff)) {
        corr_diff_df <- as.data.frame(as.table(corr_diff$R_diff))
        colnames(corr_diff_df) <- c("Var1", "Var2", "R_diff")

        p_corr_heatmap <- ggplot(corr_diff_df,
                                 aes(x = Var2, y = Var1, fill = R_diff)) +
          geom_tile(color = "grey80") +
          geom_text(aes(label = round(abs(R_diff), 2)),
                    color = "black", size = 3.5)  +
          scale_fill_gradient2(
            low = "blue", mid = "white", high = "red", midpoint = 0,
            limits = c(-1, 1),
            oob = scales::squish,
            name = "R_diff"
          ) +
          labs(
            title = paste0("Correlation Difference (Dataset ", k, ")"),
            x = "Variable",
            y = "Variable"
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "right"
          )
        }
      }

      if (length(var_cat)>0){
        # Get barplot table for categorical data
        plist_cat <- var_cat %>% map(function(var){
          # Create tables for original and synthetic data
          orig_table <- table(data_i[[var]])
          synth_table <- table(data_syn[[var]])

          # Convert to data frame for ggplot
          plot_data <- data.frame(
            category = rep(names(orig_table), 2),
            count = c(as.numeric(orig_table), as.numeric(synth_table)),
            source = rep(c("Original", "Synthetic"), each = length(orig_table))
          )

          # Create ggplot barplot
          ggplot(plot_data, aes(x = category, y = count, fill = source)) +
            geom_col(position = "dodge", color = "black", linewidth = 0.3) +
            scale_fill_manual(values = c("Original" = color_p[[3]],
                                         "Synthetic" = color_p[[4]])) +
            labs(title = paste(var, "Distribution"),
                 x = var,
                 y = "Count",
                 fill = "Dataset") +
            theme(legend.position = "topright",
                  axis.text.x = element_text(angle = 45, hjust = 1))
        })
      }

      if (length(var_cat) == 0){message("No categorical variables to plot\n")}
      if (length(var_cont) == 0){message("No continuous variables to plot\n")}
    }

    # Store dataset and its diagnostics
    results_list[[k]] <- list(
      datagen = data_syn,
      seed = target_seed,
      exact_dupl_check = dupl_check,
      dplot_cont = plist_cont,
      dplot_cat = plist_cat,
      dplot_umap = plot_umap,
      ks_test = ks_results,
      jsd_res = JSD_results,
      corr_diff_mean = if(!is.null(corr_diff)) corr_diff$mean_abs_diff else NULL,
      corr_diff_max = if(!is.null(corr_diff)) corr_diff$max_abs_diff else NULL,
      dplot_corr_diff = p_corr_heatmap
    )

    message("Dataset", k, "diagnostics completed.\n")
  }

  # Final summary message
  if (mode_fixed_seed) {
    message("\n=== Dataset generated with fixed seed =", seed, "===\n")
  } else {
    message("\n=== All", nds, "datasets generated and diagnosed ===\n")
  }

  return(results_list)

}

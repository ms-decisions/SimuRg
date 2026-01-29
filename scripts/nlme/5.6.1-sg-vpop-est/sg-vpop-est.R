## Author: Alina Melnikova
## First created: 2025-07-17
## Description: Covariate simulation from original dataset
## Keywords: covariates, simulations
## Version: v1.1 - add JS and KL metrics, covariate matrix, merge umap

#' Perform generation of synthetic dataset for an empirical distribution
#'
#' @returns A list containing:
#' \itemize{
#'   \item \code{datagen} - Synthetic dataset returned by \code{synthpop::syn()} (data.frame)
#'   \item \code{exact_dupl_check} - Logical indicating if exact duplicates exist between original and synthetic data
#'   \item \code{dplot_umap} - ggplot object with combined UMAP visualization comparing original and synthetic data (or \code{NULL} if none)
#'   \item \code{ks_test} - Tibble with Kolmogorov–Smirnov p-values and statuses for continuous variables
#'   \item \code{jsd_res} - Weighted mean Jensen–Shannon divergence (JSD) value for categorical variables (numeric)
#'   \item \code{corr_diff_mean} - Mean absolute difference between original and synthetic correlation matrices (numeric)
#'   \item \code{corr_diff_max} - Maximum absolute difference between original and synthetic correlation matrices (numeric)
#'   \item \code{dplot_corr_diff} - ggplot heatmap object showing correlation difference matrix (Synthetic - Original) (or \code{NULL} if none)
#' }
#' @inheritParams sg_dummy
#' @param data_i data frame. Input data frame containing the original dataset to be synthesized (required)
#' @param nobj integer. Specify the exact number of rows to generate in the synthetic dataset. When provided, overrides \code{npop} (optional, default: \code{NA})
#' @param minnumlev integer. Threshold; numeric variables with ≤ \code{minnumlev} unique values are converted to factors (optional, default: \code{3})
#' @param seed_umap integer. Random seed for UMAP algorithm reproducibility (optional, default: \code{123})
#' @param palette character vector. Contains color codes (hex format) for custom plot color schemes. If provided, should contain at least 2 colors. Used for histograms, bar plots, and UMAP visualizations (optional, default: \code{c("#3a6eba", "#efdd3c", "#1a1866", "#f2b93b")})
#' @param diag_plots logical flag. If \code{TRUE}, generates diagnostic plots and UMAP visualizations (optional, default: \code{TRUE})
#' @param remove_duplicates logical flag. If \code{TRUE}, automatically removes exact duplicates between original and synthetic data by adding controlled noise (optional, default: \code{TRUE})
#' @param noise_level numeric. Proportion of standard deviation to use when adding noise to continuous variables for duplicate removal. E.g., 0.10 means 10% of SD (optional, default: \code{0.10})
#' @examples
#' \dontrun{
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
#' output <- sg_vpop_est(data_i = data,#original data
#'                      id_col = "id", #name of ID column
#'                       seed = 123,#provide reproducibility
#'                       seed_umap = 40,#provide reproducibility for umap
#'                       diag_plots = T)
#' print(head(output$datagen,10))
#' print(output$ks_test)
#' print(output$dplot_umap[[1]])
#' print(output$dplot_umap)
#'
#'}
#' @import dplyr
#' @import ggplot2
#' @importFrom purrr map_dfr map
#' @importFrom stringr str_c
#' @importFrom cluster daisy
#' @importFrom fastDummies dummy_cols
#' @importFrom synthpop syn
#' @importFrom umap umap umap.defaults
#' @importFrom tidyr drop_na
#' @export
#Comparison of correlation matrices
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

#' Create optimal visit sequence based on correlations
#'
#' Orders variables so that highly correlated variables are synthesized in a way
#' that preserves their relationships. Variables with highest average correlations
#' are synthesized first.
#'
#' @param data Data frame with variables to order
#' @param var_cont Names of continuous variables
#' @param var_cat Names of categorical variables
#' @return Character vector of variable names in optimal order
create_optimal_visit_sequence <- function(data, var_cont, var_cat) {
  var_all <- c(var_cont, var_cat)

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
  visit_seq <- c(var_cont_ordered, var_cat)

  return(visit_seq)
}

#' Remove exact duplicates between synthetic and original data by adding noise
#'
#' @param data_syn Synthetic dataset
#' @param data_orig Original dataset
#' @param var_cont Names of continuous variables
#' @param var_cat Names of categorical variables
#' @param noise_level Proportion of SD to use for noise (e.g., 0.10 = 10%)
#' @param seed Random seed for reproducibility
#' @return List with cleaned synthetic data and number of duplicates removed
remove_exact_duplicates <- function(data_syn, data_orig, var_cont, var_cat,
                                     noise_level = 0.10, seed = 123) {

  set.seed(seed)

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

  cat("Found", n_duplicates, "exact duplicates. Adding noise to remove them...\n")

  # Add noise to continuous variables for duplicate rows
  if (length(var_cont) > 0) {
    for (var in var_cont) {
      if (var %in% names(data_syn)) {
        # Calculate SD from original data
        sd_val <- sd(data_orig[[var]], na.rm = TRUE)

        # Add noise only to duplicate rows
        noise <- rnorm(n_duplicates, mean = 0, sd = noise_level * sd_val)
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
          rows_to_flip <- sample(dupl_indices, size = n_to_flip)

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
sg_vpop_est <-  function(data_i, nobj = NA, id_col = NULL, minnumlev = 3,npop = 1,
                       excl_col = NULL,seed = NA, seed_umap = NA, palette = NULL,
                       diag_plots = F, show_info=F, remove_duplicates = TRUE,
                       noise_level = 0.10){

  theme_set(theme_bw())
  theme_update(panel.grid.minor = element_blank())
  MSDcol_cut <- c("#3a6eba","#efdd3c", "#1a1866", "#b73b58")

  if(!is.null(palette) & length(palette)>1){
    color_p = rep(palette,2)
  } else {color_p = MSDcol_cut
    cat("Default color palette is used")}

  #Control parameters for Random Forest (optimized for correlation preservation)
  num_trees = 300          # Increased from 100 for better stability
  max_depth = 12           # Increased from 7 to capture more interactions
  min_node_size = 5        # Decreased from 30 for finer splits
  use_smoothing = TRUE     # Enable smoothing for continuous variables
  optimize_visit_seq = TRUE # Use correlation-based visit sequence

  ### Comment!!!!
  # data_orig <- read.csv("datasets/data_gbsg2.csv")  # Change to your dataset path
  # data_i = data_orig; nobj = NA; id_col = NULL; minnumlev = 3;npop = 1
  # excl_col = NULL;seed = 42; seed_umap = 42; palette = NULL;
  # diag_plots = T; show_info=F; remove_duplicates = TRUE;
  # noise_level = 0.10

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
  #var_cont <- var_cont[var_cont != id_col]
  var_all <- c(var_cont, var_cat)



  n_orig <- nrow(data_i)


  if (show_info){
  # Print dataset summary
   cat("Number of rows in the data_orig:", nrow(data_i), "\n")
   cat("Number and names of continuous columns - ", length(var_cont), ":", paste(var_cont, collapse = ", "), "\n")
   cat("Number and names of categorical columns - ", length(var_cat), ":", paste(var_cat, collapse = ", "), "\n")
  }
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
  if (!is.na(seed)){seed_i = seed} else {seed_i = 123}
  cat("Seed for synthetic data generation:", seed_i, "\n")

  # Create optimal visit sequence based on correlations (if enabled)
  visit_sequence <- if (optimize_visit_seq && length(var_cont) >= 2) {
    cat("Creating correlation-based visit sequence...\n")
    create_optimal_visit_sequence(data_i, var_cont, var_cat)
  } else {
    NULL  # Use default order
  }

  # Prepare method vector: "rf" for all variables
  method_vector <- rep("rf", length(var_all))
  names(method_vector) <- var_all

  # Prepare smoothing vector for continuous variables (if enabled)
  smoothing_vector <- if (use_smoothing && length(var_cont) > 0) {
    # Create a named list as required by synthpop
    sm_list <- as.list(rep("density", length(var_cont)))
    names(sm_list) <- var_cont
    sm_list
  } else {
    NULL
  }

  cat("RF parameters: trees =", num_trees, ", max_depth =", max_depth,
      ", min_node_size =", min_node_size, "\n")
  if (!is.null(visit_sequence)) {
    cat("Using optimized visit sequence with", length(visit_sequence), "variables\n")
  }
  if (!is.null(smoothing_vector)) {
    cat("Smoothing enabled for", length(smoothing_vector), "continuous variables\n")
  }

  # Generate synthetic data with optimized parameters
  # Pass RF parameters directly to syn() function
  syn_obj <- syn(data_i,
                 method = method_vector,
                 visit.sequence = visit_sequence,
                 smoothing = smoothing_vector,
                 k = n_new,
                 seed = seed_i,
                 ranger.num.trees = num_trees,
                 ranger.max.depth = max_depth,
                 ranger.min.node.size = min_node_size,
                 ranger.respect.unordered.factors = "order")
  data_syn <- syn_obj$syn  # Extract synthetic data

  # Check for duplicates *between* real and synthetic data (not within each)
  common_cols <- intersect(names(data_i), names(data_syn))
  dupl_check_before <- nrow(dplyr::semi_join(data_syn, data_i, by = common_cols)) > 0
  n_dupl_removed <- 0

  if (dupl_check_before){
    warning("Exact duplicates found between original and synthetic data.")

    # Remove duplicates if requested
    if (remove_duplicates) {
      dedup_result <- remove_exact_duplicates(
        data_syn = data_syn,
        data_orig = data_i,
        var_cont = var_cont,
        var_cat = var_cat,
        noise_level = noise_level,
        seed = seed_i
      )

      data_syn <- dedup_result$data_cleaned
      n_dupl_removed <- dedup_result$n_duplicates_removed

      # Check again after deduplication
      dupl_check_after <- nrow(dplyr::semi_join(data_syn, data_i, by = common_cols)) > 0

      if (dupl_check_after) {
        cat("Warning: Some duplicates may still remain after noise addition.\n")
      } else {
        cat("Successfully removed all", n_dupl_removed, "exact duplicates.\n")
      }
    }
  }

  dupl_check <- nrow(dplyr::semi_join(data_syn, data_i, by = common_cols)) > 0

  if (length(var_cont)>0){
  # Get Kolmogorov-Smirnov test
  ks_results <- map_dfr(var_cont, function(var) {
    x <- data_i[[var]]
    y <- data_syn[[var]]

    # Detect ties across combined samples; ks.test assumes continuous data
    has_ties <- any(duplicated(c(x, y)))

    # Suppress the default ks.test warning about ties and emit a custom one instead
    ks_result <- suppressWarnings(ks.test(x, y))
    if (has_ties) {
      warning(
        "Kolmogorov–Smirnov test for variable '", var,
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

  #Comparison of correlation matrices
  corr_diff <- compare_cor_matrices(data_i, data_syn, var_cont, method = "kendall")

  }
  if (length(var_cat)>0){


    n_levels <- c()
    # Jensen–Shannon (JS) divergence
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
    print(jsd_by_var)

    jsd <- as.vector(jsd_by_var)
    names(jsd) <- var_cat
    n_levels <- sapply(var_cat, function(var) {
      length(levels(factor(data_i[[var]])))
    })
    names(n_levels) <- var_cat
    # Calculate weighted mean of Jensen-Shannon divergence
    JSD_results<- sum(jsd * n_levels / sum(n_levels))
  }
  #Preparation for UMAP

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
    if (!is.na(seed_umap)){
      seed_ii = seed_umap
      set.seed(seed_ii)
    } else {seed_ii = seed_i}

    cat("Seed for synthetic data generation:", seed_ii, "\n")

      umap_model <- umap(
      X_real,
      n_neighbors = 15,
      min_dist = 0.1,
      metric = "euclidean",
      ret_model = TRUE
    )

    # project synthetic data
    emb_real <- umap_model$embedding
    emb_syn  <- umap_transform(X_syn, umap_model)

    # Create combined visualization
    umap_df <- rbind(
      data.frame(X = emb_real[, 1], Y = emb_real[, 2], Source = "Original data"),
      data.frame(X = emb_syn[, 1],  Y = emb_syn[, 2],  Source = "Synthetic data")
    )

    # Create plot
    plot_umap <- ggplot(umap_df, aes(x = X, y = Y, color = Source)) +
      geom_point(alpha = 0.4, size = 1) +
      scale_color_manual(values = c("Original data" = color_p[[3]], "Synthetic data" = color_p[[4]])) +
      labs(title = "UMAP: Combined Data Comparison", x = "UMAP Dimension 1", y = "UMAP Dimension 2", color = "Dataset") +
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
        #theme_minimal() +
        theme(legend.position = "topright")
    })

    #Heat map for corr_diff R_diff (using ggplot2 for better visualization)
    corr_diff_df <- as.data.frame(as.table(corr_diff$R_diff))
    colnames(corr_diff_df) <- c("Var1", "Var2", "R_diff")

    p_corr_heatmap <- ggplot(corr_diff_df,
                             aes(x = Var2, y = Var1, fill = R_diff)) +
            geom_tile(color = "grey80") +
      geom_text(aes(label = round(abs(R_diff), 2)),
                color = "black", size = 3.5)  +
      #geom_text(aes(label = "TEST"), color = "black", size = 10) +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red", midpoint = 0,
        limits = c(-1, 1),                # fix color scale to absolute range
        oob = scales::squish,             # clamp any rounding noise to [-1, 1]
        name = "R_diff"
      ) +
      labs(
        title = "Correlation Difference (Synthetic - Original)",
        x = "Variable",
        y = "Variable"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )
    #ptable_hist_cont <- do.call("arrangeGrob", c(plist_cont, ncol = 2))
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
    if (length(var_cat) == 0){cat("No categorical variables to plot")}
    if (length(var_cont) == 0){cat("No continuous variables to plot")}
  }


  return(list(datagen = data_syn,
              exact_dupl_check = dupl_check,
              #n_duplicates_removed = n_dupl_removed,
              dplot_cont =plist_cont,
              dplot_cat = plist_cat,
              dplot_umap = plot_umap,
              ks_test = ks_results,
              #kld_res = KLD_results,
              jsd_res = JSD_results,
              corr_diff_mean = corr_diff$mean_abs_diff,
              corr_diff_max = corr_diff$max_abs_diff,
              dplot_corr_diff = p_corr_heatmap
  )
  )

}

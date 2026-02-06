## Author: Alina Melnikova
## First created: 2025-07-17
## Description: Covariate simulation from original dataset
## Keywords: covariates, simulations

#' Perform generation of synthetic dataset for an empirical distribution
#'
#' @returns A list containing:
#' \itemize{
#'   \item \code{datagen} - Synthetic dataset returned by \code{synthpop::syn()} (data.frame)
#'   \item \code{dplot_cont} - List of ggplot histogram objects, one per continuous variable (or \code{NULL} if none)
#'   \item \code{dplot_cat} - List of ggplot bar-plot objects, one per categorical variable (or \code{NULL} if none)
#'   \item \code{dplot_umap_cont} - ggplot object with the continuous UMAP comparison (or \code{NULL})
#'   \item \code{dplot_umap_cat} - ggplot object with the categorical UMAP comparison (or \code{NULL})
#'   \item \code{ks_test} - Tibble with Kolmogorov–Smirnov p-values and statuses for continuous variables
#' }
#' @inheritParams sg_dummy
#' @param data_i data frame. Input data frame containing the original dataset to be synthesized (required)
#' @param nobj integer. Specify the exact number of rows to generate in the synthetic dataset. When provided, overrides \code{npop} (optional, default: \code{NA})
#' @param minnumlev integer. Threshold; numeric variables with ≤ \code{minnumlev} unique values are converted to factors (optional, default: \code{3})
#' @param seed_umap integer. Random seed for UMAP algorithm reproducibility (optional, default: \code{123})
#' @param palette character vector. Contains color codes (hex format) for custom plot color schemes. If provided, should contain at least 2 colors. Used for histograms, bar plots, and UMAP visualizations (optional, default: \code{c("#3a6eba", "#efdd3c", "#1a1866", "#f2b93b")})
#' @param diag_plots logical flag. If \code{TRUE}, generates diagnostic plots and UMAP visualizations (optional, default: \code{TRUE})
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
#' print(output$dplot_umap_cont[[1]])
#' print(output$dplot_umap_cat[[1]])
#' print(output$dplot_umap_cont)
#' print(output$dplot_umap_cat)
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
sg_vpop_est <-  function(data_i, nobj = NA, id_col = NULL, minnumlev = 3,npop = 1,
                       excl_col = NULL,seed = NA, seed_umap = NA, palette = c("#3a6eba","#efdd3c", "#1a1866", "#f2b93b"),
                       diag_plots = T){

  theme_set(theme_bw())
  theme_update(panel.grid.minor = element_blank())
  MSDcol_cut <- c("#3a6eba","#efdd3c", "#1a1866", "#f2b93b")

  if(!is.null(palette) & length(palette)>1){
    color_p = rep(palette,2)
  } else {color_p = MSDcol_cut
    cat("Default color palette is used")}
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



  # Print dataset summary
   cat("Number of rows in the data_orig:", nrow(data_i), "\n")
   cat("Number and names of continuous columns - ", length(var_cont), ":", paste(var_cont, collapse = ", "), "\n")
   cat("Number and names of categorical columns - ", length(var_cat), ":", paste(var_cat, collapse = ", "), "\n")

  data_orig <- data_i

  if (!is.na(nobj)){n_new = nobj} else if(!is.na(npop)&(npop>0)){
    n_new = as.integer(round(n_orig * npop))
    } else {n_new = n_orig}

  #Synthpop Method
  if (!is.na(seed)){seed_i = seed} else {seed_i = 123}
  print(seed_i)
  syn_obj <- syn(data_i, method = "rf", k = n_new, seed = seed_i)
  data_syn <- syn_obj$syn  # Extract synthetic data



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

    p_val_formatted <- ifelse(p_val < 0.01, "<0.01", p_val)
    tibble(
      variable = var,
      p.value = p_val_formatted,
      status = ifelse(p_val > 0.05, "similar", "different")
    )
  })
  # ks_results <- map_dfr(var_cont, function(var) {
  #   ks_result <- ks.test(data_i[[var]], data_syn[[var]])
  #   p_val <- ks_result$p.value
  #
  #   p_val_formatted <- ifelse(p_val < 0.01, "<0.01", p_val)
  #   tibble(
  #     variable = var,
  #     p.value = p_val_formatted,
  #     status = ifelse(p_val > 0.05, "similar", "different")
  #   )
  # })

  plist_cont = NULL
  plist_cat = NULL
  plot_umap_cont = NULL
  plot_umap_cat = NULL

  if (diag_plots){
    if (length(var_cont)>0){

    # Get UMAP for continuous data
    combined_data <- rbind(
      cbind(data_i[,var_cont], Source = "Original data"),
      cbind(data_syn[,var_cont], Source = "Synthetic data")
    )


    data_for_umap <- combined_data[, -which(colnames(combined_data) == "Source")]

    # Remove duplicate rows to satisfy UMAP requirements
    dup_mask_umap <- duplicated(data_for_umap)
    if (any(dup_mask_umap)) {
      data_for_umap <- data_for_umap[!dup_mask_umap, , drop = FALSE]
      combined_data_umap <- combined_data[!dup_mask_umap, , drop = FALSE]
    } else {
      combined_data_umap <- combined_data
    }

    # Run UMAP on combined data
    if (!is.na(seed_umap)){seed_ii = seed_umap} else {seed_ii = seed_i}

    umap_combined <- umap(data_for_umap, n_neighbors = 15, min_dist = 0.1, random_state = seed_ii)

    #
    # Create combined visualization
    umap_df <- data.frame(
      X = umap_combined$layout[,1],
      Y = umap_combined$layout[,2],
      Source = combined_data_umap$Source
    )

    plot_umap_cont <- ggplot(umap_df, aes(x = X, y = Y, color = Source)) +
      geom_point(alpha = 0.4, size = 1) +
      scale_color_manual(values = c("Original data" = color_p[[3]], "Synthetic data" = color_p[[4]])) +
      labs(
        title = str_c("UMAP: Continuous Data Comparison. Seed_umap = ", seed_ii),
        subtitle = paste0("Original: N=", sum(umap_df$Source == "Original data"),
                          ", Synthetic: N=", sum(umap_df$Source == "Synthetic data")),
        x = "UMAP Dimension 1",
        y = "UMAP Dimension 2",
        color = "Dataset"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      )

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
        # scale_fill_manual(
        #   values = setNames(
        #     color_p[c(1, 2, 3)],
        #     c("Original", "Synthetic", "Difference")
        #   )
        # ) +
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

    #ptable_hist_cont <- do.call("arrangeGrob", c(plist_cont, ncol = 2))
    }



    #Umap for categorical data
    if (length(var_cat)>0){
    cat_orig <- data_i[, var_cat, drop = FALSE]
    cat_syn <- data_syn[, var_cat, drop = FALSE]

    # Combine for processing
    combined_cat <- rbind(
      cbind(cat_orig, Source = "Original data"),
      cbind(cat_syn, Source = "Synthetic data")
    )

    # One-hot encode categorical data for UMAP
    # Remove the Source column temporarily
    combined_cat_no_source <- combined_cat[, -which(colnames(combined_cat) == "Source")]

    # Convert to dummy variables using fastDummies
    cat_encoded <- dummy_cols(combined_cat_no_source,
                              remove_first_dummy = FALSE,
                              remove_selected_columns = TRUE)


    # Run UMAP
    umap_result <- umap(cat_encoded, n_neighbors = 15, min_dist = 0.1, random_state = seed_ii)
    # umap_result <- umap(cat_encoded, config = umap_config)

    # Create data frame for visualization
    umap_cat_df <- data.frame(
      X = umap_result$layout[, 1],
      Y = umap_result$layout[, 2],
      Source = combined_cat$Source
    )

    # Create plot
    plot_umap_cat <-  ggplot(umap_cat_df, aes(x = X, y = Y, color = Source)) +
      geom_point(alpha = 0.4, size = 1) +
      scale_color_manual(values = c("Original data" = color_p[[3]],
                                    "Synthetic data" = color_p[[4]])) +
      labs(
        title = str_c("UMAP: Categorical Data Comparison. Seed_umap = ", seed_ii),
        subtitle = paste0("Original: N=", sum(umap_cat_df$Source == "Original data"),
                          ", Synthetic: N=", sum(umap_cat_df$Source == "Synthetic data"),
                          " | ", length(var_cat), " categorical variables"),
        x = "UMAP Dimension 1",
        y = "UMAP Dimension 2",
        color = "Dataset"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      )

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


    #ptable_bp_cat <- do.call("arrangeGrob", c(plist_cat, ncol = 3))
}
    if (length(var_cat) == 0){cat("No categorical variables to plot")}
    if (length(var_cont) == 0){cat("No continuous variables to plot")}
  }


  return(list(datagen = data_syn,
              dplot_cont =plist_cont,
              dplot_cat =plist_cat,
              dplot_umap_cont = plot_umap_cont,
              dplot_umap_cat = plot_umap_cat,
              ks_test = ks_results
  )
  )

}

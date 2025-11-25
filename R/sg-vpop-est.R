## Author: Alina Melnikova
## First created: 2025-07-17
## Description: Covariate simulation from original dataset
## Keywords: covariates, simulations

#' Perform generation of synthetic dataset
#'
#' @inheritParams sg_dummy
#' @returns A dataset with simulation results, 4 diagnostic plots and results of Kolmogorov-Smirnov test
#' @param data - input dataset
#' @param idcol - string name of ID column
#' @param minnumlev - maximum number of numeric variable levels that will be regarded as factor. Set to 3 by default
#' @param seed - seed number for data synthesis. Set to 123 by default
#' @param seed_umap - seed number for umap plots. Set to 123 by default
#' @param palette - vector with user color palette
#' @param diag_plots - logical, TRUE by default. Set TRUE to generate diagnostic plots, set FALSE otherwise.
#' @examples
#' \dontrun{
#' library(tidyverse)
#' library(readr)
#' library(dplyr)
#' library(grid)
#'
#' load("data/data_pbc.rda")
#' sd_umap = 40
#' output <- sg_vpop_est(data = data,#original data
#'                      idcol = "id", #name of ID column
#'                       seed = 123,#provide reproducibility
#'                       seed_umap = sd_umap,#provide reproducibility for umap
#'                       diag_plots = T
#' print(head(output$datagen,10))
#' print(output$ks_test)
#' grid.draw(output$dplot_cont)
#' grid.draw(output$dplot_cat)
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
#' @importFrom gridExtra arrangeGrob
#' @importFrom synthpop syn
#' @importFrom umap umap umap.defaults
#' @importFrom tidyr drop_na
#' @export


sg_vpop_est = function(data,
                       nobj = NA,
                       idcol = NULL,
                       minnumlev = 3,
                       expfctr = 1,
                       exclcol = NULL,
                       seed = NA,
                       seed_umap = NA,
                       palette = NULL,
                       diag_plots = T){


  theme_set(theme_bw())
  theme_update(panel.grid.minor = element_blank())
  MSDcol <- c("#1a1866", "#f2b93b", "#b73b58", "#a2d620", "#5839bb", "#9c4ec7",
                  "#3a6eba","#efdd3c", "#69686d")

  if(!is.null(palette)){color_p = palette} else {color_p = MSDcol}
  #####--------------- Processing ---------------#####

  if (!is.null(idcol)){data <- data %>%
    dplyr::select(-all_of(idcol))}

  if (!is.null(exclcol)){
    data <- data %>%
      dplyr::select(-all_of(exclcol))
  }


  #Exclude rows with NAs.
  data <- data %>%
    drop_na()

  var_char = NULL

  # Convert character columns to factor
  is_char  <- vapply(data, is.character, logical(1))
  var_char  <- names(data)[is_char]

  if (!is.na(minnumlev)){minnum = minnumlev} else {minnum = 3}

  # Convert numeric variables with <= minnumlev unique values to factor
  is_low_level_num <- vapply(
    data,
    function(x) is.numeric(x) && length(unique(na.omit(x))) <= minnum,
    logical(1)
  )
  var_low_level <- names(data)[is_low_level_num]


  data <- data %>%
    mutate(across(any_of(c(var_char, var_low_level)), as.factor))

  is_cont <- vapply(data, is.numeric, logical(1))
  is_cat  <- vapply(data, is.factor, logical(1))

  # Define names for continuous and categorical (ordered/factor) columns
  var_cont <- names(data)[is_cont]
  var_cat  <- names(data)[is_cat]
  #var_cont <- var_cont[var_cont != idcol]
  var_all <- c(var_cont, var_cat)



  n_orig <- nrow(data)



  # Print dataset summary
   cat("Number of rows in the data_orig:", nrow(data), "\n")
   cat("Number and names of continuous columns - ", length(var_cont), ":", paste(var_cont, collapse = ", "), "\n")
   cat("Number and names of categorical columns - ", length(var_cat), ":", paste(var_cat, collapse = ", "), "\n")

  data_orig <- data

  if (!is.na(nobj)){n_new = nobj} else if(!is.na(expfctr)&(expfctr>0)){
    n_new = as.integer(round(n_orig * expfctr))
    } else {n_new = n_orig}

  #Synthpop Method
  if (!is.na(seed)){seed_i = seed} else {seed_i = 123}
  print(seed_i)
  syn_obj <- syn(data, method = "rf", k = n_new, seed = seed_i)
  data_syn <- syn_obj$syn  # Extract synthetic data



  # Get Kolmogorov-Smirnov test

  ks_results <- map_dfr(var_cont, function(var) {
    ks_result <- ks.test(data[[var]], data_syn[[var]])
    p_val <- ks_result$p.value

    p_val_formatted <- ifelse(p_val < 0.01, "<0.01", p_val)
    tibble(
      variable = var,
      p.value = p_val_formatted,
      status = ifelse(p_val > 0.05, "similar", "different")
    )
  })

  ptable_hist_cont = NULL
  ptable_bp_cat  = NULL
  plot_umap_cont = NULL
  plot_umap_cat = NULL

  if (diag_plots){
    if (length(var_cont)>0){
    # Get UMAP
    combined_data <- rbind(
      cbind(data[,var_cont], Source = "Original data"),
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
      scale_color_manual(values = c("Original data" = color_p[[1]], "Synthetic data" = color_p[[2]])) +
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
        value = c(data[[var]], data_syn[[var]]),
        source = rep(c("Original", "Synthetic"),
                     times = c(length(data[[var]]), length(data_syn[[var]])))
      )

      # Create ggplot histogram
      ggplot(plot_data, aes(x = value, fill = source)) +
        geom_histogram(alpha = 0.6, bins = 30, position = "identity",
                       color = "black", linewidth = 0.3) +
        scale_fill_manual(
          values = setNames(
            color_p[c(7, 8, 4)],
            c("Original", "Synthetic", "Difference")
          )
        ) +
        labs(title = paste(var, "Distribution"),
             x = var,
             y = "Frequency",
             fill = "Dataset") +
        #theme_minimal() +
        theme(legend.position = "topright")
    })

    ptable_hist_cont <- do.call("arrangeGrob", c(plist_cont, ncol = 2))
    }



    #Umap for categorical data
    if (length(var_cat)>0){
    cat_orig <- data[, var_cat, drop = FALSE]
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
      scale_color_manual(values = c("Original data" = color_p[[1]],
                                    "Synthetic data" = color_p[[2]])) +
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
      orig_table <- table(data[[var]])
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
        scale_fill_manual(values = c("Original" = color_p[[1]],
                                     "Synthetic" = color_p[[2]])) +
        labs(title = paste(var, "Distribution"),
             x = var,
             y = "Count",
             fill = "Dataset") +
        theme(legend.position = "topright",
              axis.text.x = element_text(angle = 45, hjust = 1))
    })


    ptable_bp_cat <- do.call("arrangeGrob", c(plist_cat, ncol = 3))
}
    if (length(var_cat) == 0){cat("No categorical variables to plot")}
    if (length(var_cont) == 0){cat("No continuous variables to plot")}
  }


  return(list(datagen = data_syn,
              dplot_cont = ptable_hist_cont,
              dplot_cat  = ptable_bp_cat,
              dplot_umap_cont = plot_umap_cont,
              dplot_umap_cat = plot_umap_cat,
              ks_test = ks_results
  )
  )

}

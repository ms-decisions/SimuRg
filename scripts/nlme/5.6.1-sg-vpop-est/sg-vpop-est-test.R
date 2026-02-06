## Author: Alina Melnikova
## First created: 2025-07-17
## Description: Covariate simulation from original dataset
## Keywords: covariates, simulations


#library(tidyverse)
library(dplyr)
library(cluster)      # For distance calculations (Gower distance)
library(fastDummies)
library(synthpop)
library(umap)

sg_vpop_est = function(data,
                       nobj = NA,
                       idcol = NULL,
                       expfctr = 1,
                       exclcol = NULL,
                       seed = NA,
                       seed_umap = NULL,
                       diag_plots = T){


theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())
MSDcol <- rep(c("#1a1866", "#f2b93b", "#b73b58", "#a2d620", "#5839bb", "#9c4ec7",
                "#3a6eba","#efdd3c", "#69686d"), 5)


#####--------------- Processing ---------------#####
# data <- read_csv("functions/datasets/data_aids.csv")
# nobj = NA
# idcol = "id"
# expfctr = 1
# exclcol = NULL
# seed = NA
# seed_umap = c(40,41,42)
# diag_plots = F

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
var_binary = NULL
var_cat3 = NULL

# Convert character columns to factor
is_char  <- vapply(data, is.character, logical(1))
var_char  <- names(data)[is_char]

# Also convert binary numeric variables to factor
is_binary <- vapply(data, function(x) length(unique(x)) == 2, logical(1))
var_binary <- names(data)[is_binary & !is_char]

is_cat3 <- vapply(data, function(x) length(unique(x)) == 3, logical(1))
var_cat3 <- names(data)[is_cat3 & !is_char]


  data <- data %>%
    mutate(across(any_of(c(var_char, var_binary, var_cat3)), as.factor))


# if (!is.null(catcol)){
#   data <- data %>%
#     mutate(across(all_of(c(var_char, catcol)), as.factor))
# } else {
#   data <- data %>%
#     mutate(across(all_of(var_char), as.factor))
# }

n_orig <- nrow(data)

is_cont <- vapply(data, is.numeric, logical(1))
is_cat  <- vapply(data, is.factor, logical(1))

# Define names for continuous and categorical (ordered/factor) columns
var_cont <- names(data)[is_cont]
var_cat  <- names(data)[is_cat]
#var_cont <- var_cont[var_cont != idcol]
var_all <- c(var_cont, var_cat)

# Print dataset summary
# cat("Number of rows in the data_orig:", nrow(data), "\n")
# cat("Number and names of continuous columns - ", length(var_cont), ":", paste(var_cont, collapse = ", "), "\n")
# cat("Number and names of categorical columns - ", length(var_cat), ":", paste(var_cat, collapse = ", "), "\n")

data_orig <- data

if (!is.na(nobj)){n_new = nobj} else {n_new = n_orig*expfctr}


#Synthpop Method
if (!is.na(seed)){seed_i = seed} else {seed_i = 123}
print(seed_i)

syn_obj <- syn(data, method = "rf", k = n_new, seed = seed_i)
data_syn <- syn_obj$syn  # Extract synthetic data


# Get Kolmogorov-Smirnov test

ks_results <- map_dfr(var_cont, function(var) {
  ks_result <- ks.test(data[[var]], data_syn[[var]])
  p_val <- ks_result$p.value
  #p_val_formatted <- ifelse(p_val < 0.01, "<0.01", as.character(round(p_val, 2)))
  p_val_formatted <- ifelse(p_val < 0.01, "<0.01", p_val)
  tibble(
    variable = var,
    p.value = p_val_formatted,
    status = ifelse(p_val > 0.05, "similar", "different")
  )
})

ptable_hist_cont = NULL
ptable_bp_cat  = NULL
plot_table_umap_cont = NULL
plot_table_umap_cat = NULL

if (diag_plots){
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
# if (!is.na(seed_umap)){seed_ii = seed_umap} else {seed_ii = seed_i}
#
#   umap_combined <- umap(data_for_umap, n_neighbors = 15, min_dist = 0.1, random_state = seed_ii)

if (!is.null(seed_umap)){seed_vec = seed_umap} else {seed_vec = c(40,41,42) }

  plot_umap_cont <- seed_vec %>% map(function(x){
    set.seed(x)
    umap_config <- umap.defaults
    umap_config$n_neighbors <- 15
    umap_config$min_dist <- 0.1
    umap_combined <- umap(data_for_umap, config = umap_config)

  # Create combined visualization
  umap_df <- data.frame(
    X = umap_combined$layout[,1],
    Y = umap_combined$layout[,2],
    Source = combined_data_umap$Source
  )

  plot_umap_cont <- ggplot(umap_df, aes(x = X, y = Y, color = Source)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Original data" = "#1a1866", "Synthetic data" = "#f2b93b")) +
    labs(
      title = str_c("UMAP: Continuous Data Comparison. Seed_umap = ", x),
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
})
plot_table_umap_cont <- do.call("arrangeGrob", c(plot_umap_cont, ncol = 3))
# grid.draw(plot_table_umap_cont)
#print(plot_umap_cont)

#Umap for categorical data
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

# Remove duplicate rows if any
# dup_mask_cat <- duplicated(cat_encoded)
# if (any(dup_mask_cat)) {
#   cat_encoded <- cat_encoded[!dup_mask_cat, , drop = FALSE]
#   combined_cat_filtered <- combined_cat[!dup_mask_cat, , drop = FALSE]
# } else {
#   combined_cat_filtered <- combined_cat
# }

# Run UMAP on categorical data with multiple seeds
plot_umap_cat <- seed_vec %>% map(function(x) {
  set.seed(x)

  # Configure UMAP parameters
  umap_config <- umap.defaults
  umap_config$n_neighbors <- 15
  umap_config$min_dist <- 0.1
  umap_config$metric <- "euclidean"  # Good for one-hot encoded data

  # Run UMAP
  umap_result <- umap(cat_encoded, config = umap_config)

  # Create data frame for visualization
  umap_cat_df <- data.frame(
    X = umap_result$layout[, 1],
    Y = umap_result$layout[, 2],
    Source = combined_cat$Source
  )

  # Create plot
  ggplot(umap_cat_df, aes(x = X, y = Y, color = Source)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_color_manual(values = c("Original data" = "#1a1866",
                                  "Synthetic data" = "#f2b93b")) +
    labs(
      title = str_c("UMAP: Categorical Data Comparison. Seed umap = ", x),
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
})

# Arrange plots in grid
plot_table_umap_cat <- do.call("arrangeGrob", c(plot_umap_cat, ncol = 3))
#grid.draw(plot_table_umap_cat)


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
    scale_fill_manual(values = c("Original" = rgb(1,0,0,0.4),
                                 "Synthetic" = rgb(0,0,1,0.4))) +
    labs(title = paste(var, "Distribution"),
         x = var,
         y = "Frequency",
         fill = "Dataset") +
    #theme_minimal() +
    theme(legend.position = "topright")
})

ptable_hist_cont <- do.call("arrangeGrob", c(plist_cont, ncol = 2))
#grid.draw(ptable_hist_cont )


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
    scale_fill_manual(values = c("Original" = MSDcol[[1]],
                                 "Synthetic" = MSDcol[[2]])) +
    labs(title = paste(var, "Distribution"),
         x = var,
         y = "Count",
         fill = "Dataset") +
    theme(legend.position = "topright",
          axis.text.x = element_text(angle = 45, hjust = 1))
})


ptable_bp_cat <- do.call("arrangeGrob", c(plist_cat, ncol = 3))
#grid.draw(ptable_bp_cat)
}


return(list(datagen = data_syn,
            dplot_cont = ptable_hist_cont,
            dplot_cat  = ptable_bp_cat,
            #dplot_umap_cont = plot_umap_cont,
            dplot_umap_cont = plot_table_umap_cont,
            dplot_umap_cat = plot_table_umap_cat,
            ks_test = ks_results
            )
       )

}

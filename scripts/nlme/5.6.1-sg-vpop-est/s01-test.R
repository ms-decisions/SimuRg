rm(list=ls())
#library(tidyverse)
#library(readr)
library(dplyr)
#library(gridExtra)
library(grid)
#library(cowplot)

theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

source("scripts/nlme/5.6.1-sg-vpop-est/sg-vpop-est.R")

# data <- readRDS("data/data_pbc.rda")
# load("data/data_pbc.rda")
# #data <- read_csv("scripts/nlme/5.6.1-sg-vpop-est/datasets/data_pbc.csv")

# Generate example dataset
data <- data.frame(
  id = 1:150,
  ALB = rnorm(150, mean = 3.5, sd = 0.5),
  ALT = rnorm(150, mean = 50, sd = 20),
  SEX = factor(sample(c("M", "F"), 150, replace = TRUE)),
  RACE = factor(sample(c("White", "Black", "Asian", "Other"), 150, replace = TRUE))
)

# Get syntesized dataset. To get another syntesized dataset - change seed
output <- sg_vpop_est(data_i = data,#original data
                      id_col = "id", #name of ID column
                      seed = 123,#provide reproducibility
                      seed_umap = 40,#provide reproducibility for umap
                      diag_plots = T
)
View(output$datagen)

print(output$dplot_cont)
print(output$dplot_cat)
print(output$dplot_umap_cont)
print(output$dplot_umap_cat)
print(output$ks_test)

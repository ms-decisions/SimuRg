rm(list=ls())
#library(tidyverse)
library(readr)
library(dplyr)
#library(gridExtra)
library(grid)
#library(cowplot)

theme_set(theme_bw())
theme_update(panel.grid.minor = element_blank())

source("scripts/nlme/5.6.1-sg-vpop-est/sg-vpop-est.R")

data <- readRDS("data/data_pbc.rda")
load("data/data_pbc.rda")
#data <- read_csv("scripts/nlme/5.6.1-sg-vpop-est/datasets/data_pbc.csv")

# Get syntesized dataset. To get another syntesized dataset - change seed
sd_umap = c(40,41)
output <- sg_vpop_est(data = data,#original data
                      idcol = "id", #name of ID column
                      seed = 123,#provide reproducibility
                      seed_umap = sd_umap,#provide reproducibility for umap
                      diag_plots = T
)
View(output$datagen)

grid.draw(output$dplot_cont)
grid.draw(output$dplot_cat)
grid.draw(output$dplot_umap_cont)
grid.draw(output$dplot_umap_cat)
print(output$ks_test)

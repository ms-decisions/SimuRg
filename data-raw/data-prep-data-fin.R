## Author: Alina Melnikova
## First created: 2026-03-20
## Description: functions for external data preparation
## Keywords: SimuRg, external data

ds_covval <- read.csv(file.path("inst", "extdata",  "datasets", "data_fin_i.csv"))

usethis::use_data(ds_covval, overwrite = TRUE)

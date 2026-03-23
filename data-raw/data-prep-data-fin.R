## Author: Alina Melnikova
## First created: 2026-03-20
## Description: functions for external data preparation
## Keywords: SimuRg, external data

estcovmat <- read.csv(file.path("inst", "extdata",  "datasets", "data_fin_i.csv"))

usethis::use_data(estcovmat, overwrite = TRUE)

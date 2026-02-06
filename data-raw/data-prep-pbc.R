## Author: Alina Melnikova
## First created: 2025-11-24
## Description: functions for external data preparation - pbc dataset
## Keywords: SimuRg, external data

data_pbc<- read.csv(file.path("inst", "extdata",  "datasets", "data_pbc.csv"))
usethis::use_data(data_pbc, overwrite = TRUE)

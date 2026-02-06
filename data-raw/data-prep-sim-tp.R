## Author: Melnikova Alina
## First created: 2025-12-14
## Description: functions for external data preparation
## Keywords: SimuRg, external data

sim_timeprof <- read.csv(file.path("inst", "extdata",  "datasets", "sim_timeprof_test.csv"))
usethis::use_data(sim_timeprof, overwrite = TRUE)

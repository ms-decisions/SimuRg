## Author: Alina Melnikova
## First created: 2026-03-18
## Description: functions for external data preparation
## Keywords: SimuRg, external data

#gfo4cov <- read.csv(file.path("inst", "extdata",  "datasets", "run_4cov_smrg_results_mod.json"))
fpath <- system.file("extdata",  "simurg_object", "run_4cov_smrg_results_mod.json", package = "SimuRg")
gfo4cov <- read_smrg_obj(fpath)
usethis::use_data(gfo4cov, overwrite = TRUE, internal = TRUE)

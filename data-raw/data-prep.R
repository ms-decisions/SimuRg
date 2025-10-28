## Author: Mikhailova Anna
## First created: 2025-08-25
## Description: functions for external data preparation
## Keywords: SimuRg, external data

warfarin <- read.csv(file.path("inst", "extdata",  "datasets", "dspk-warf.csv"))
usethis::use_data(warfarin, overwrite = TRUE)

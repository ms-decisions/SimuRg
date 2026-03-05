## Author: Mikhailova Anna
## First created: 2025-08-25
## Description: functions for external data preparation
## Keywords: SimuRg, external data

warfarin <- read.csv(file.path("inst", "extdata",  "datasets", "dspk-warf.csv")) %>%
  select(-c(AGE_centered, G1_1, G1_2, G1_3, G2_2, G2_3, G3_3, GG, AG, AA))
usethis::use_data(warfarin, overwrite = TRUE)

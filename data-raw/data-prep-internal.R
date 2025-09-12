## Author: Mikhailova Anna
## First created: 2025-08-25
## Description: functions for internal data preparation
## Keywords: SimuRg, internal data

fpath_i <- file.path("inst", "extdata", "simurg_object", "Warfarin_PK.RData")
obj1 <- get(load(fpath_i))
usethis::use_data(obj1, internal = T, overwrite = T)

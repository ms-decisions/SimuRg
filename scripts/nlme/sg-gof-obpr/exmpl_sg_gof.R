library(tidyverse)
# library(rxode2)
# library(SimuRg)
devtools::load_all()

smrg <- read_smrg_obj(fpath)
#### AV tests ####
fpath <- "../simurg_r/scripts/nlme/sg-gof-obpr/1cmt-RE-Vd-CL-prop-FEMALE-on-Vd-CRCL-on-CL.RData"
# Basic plot
p0 <- sg_gof_obpr(fpath_i = fpath)
p0
# With categorical covariate and faceting
p1 <- sg_gof_obpr(
  fpath_i = fpath,
  cov_cols = "FEMALE",
  col_i = "FEMALE",
  col_lab = "Female",
  facet_i = "FEMALE",
  addline = F,
  f_scales = "free_x",
  lab_y = "X-axis",
  lab_x = "Y-axis",
  #no_leg = T
  #abreaks = seq(0, 70, 5)
  #smooth = F,
  #log_axes = T
  #alpha_i = 0.1
)

p1
# With several covariates and faceting
p3 <- sg_gof_obpr(
  fpath_i = fpath,
  cov_cols = c("FEMALE", "SLE"),
  #col_i = "FEMALE",
  #col_lab = "Female",
  facet_i = c("SLE", "FEMALE"),
  addline = F,
  f_scales = "free",
  #no_leg = T
  #abreaks = seq(0, 70, 5)
  #smooth = F,
  #log_axes = T
  #alpha_i = 0.1
)
p3
# With continuous covariate and faceting
p4 <- sg_gof_obpr(
  fpath_i = fpath,
  cov_cols = c("WT"),
  #col_i = "FEMALE",
  #col_lab = "Female",
  facet_i = "WT",
  addline = F,
  f_scales = "free",
  no_leg = T
  #abreaks = seq(0, 70, 5)
  #smooth = F,
  #log_axes = T
  #alpha_i = 0.1
)
p4
# With continuous covariate and faceting
p5 <- sg_gof_obpr(
  fpath_i = fpath,
  cov_cols = c("FEMALE","WT"),
  col_i = "WT",
  col_lab = "Weight",
  facet_i = c("FEMALE","WT"),
  addline = F,
  f_scales = "free",
  n_quantiles = 4
  #no_leg = T
  #abreaks = seq(0, 70, 5)
  #smooth = F,
  #log_axes = T
  #alpha_i = 0.1
)
p5



devtools::load_all()

fpath <- "scripts/nlme/sg-gof-obpr/1cmt-RE-Vd-CL-prop-FEMALE-on-Vd-CRCL-on-CL.RData"
if (!file.exists(fpath)) {
  fpath <- file.path("..", "sg-gof-obpr", "1cmt-RE-Vd-CL-prop-FEMALE-on-Vd-CRCL-on-CL.RData")
}

# Basic residual plots
p0 <- sg_gof_res(fpath_i = fpath, smooth = FALSE)
p0

# Versus predictions, population WRES (this fit has WRES, not CWRES)
p1 <- sg_gof_res(
  fpath_i = fpath,
  indiv = FALSE,
  vs_time = FALSE,
  weighted = TRUE,
  addline = FALSE
)
p1

# sc_factor on time and on predictions
p_sc_time <- sg_gof_res(fpath_i = fpath, vs_time = TRUE, sc_factor = 2, addline = FALSE)
p_sc_pred <- sg_gof_res(
  fpath_i = fpath,
  vs_time = FALSE,
  sc_factor = 10,
  addline = FALSE
)
p_sc_time
p_sc_pred

# exp_dv: IPRED/DV in this file are already on the linear scale, so exponentiate
# a log-transformed copy to check the back-transform.
smrg <- read_smrg_obj(fpath)
smrg$SDTAB$IPRED <- log(pmax(smrg$SDTAB$IPRED, 1e-8))
smrg$SDTAB$PRED <- log(pmax(smrg$SDTAB$PRED, 1e-8))
smrg$SDTAB$DV <- log(pmax(smrg$SDTAB$DV, 1e-8))
p_exp <- sg_gof_res(
  smrg,
  vs_time = FALSE,
  exp_dv = TRUE,
  sc_factor = 1,
  addline = FALSE,
  smooth = FALSE
)
p_exp

# Continuous covariate: legend labels include quantile bounds (Q1, lo-hi)
p_q <- sg_gof_res(
  fpath_i = fpath,
  cov_cols = c("WT", "FEMALE"),
  col_i = "WT",
  col_lab = "Weight",
  facet_i = "FEMALE",
  legend_fl = TRUE,
  n_quantiles = 4,
  addline = FALSE
)
p_q

# Log x labels: ordinary numbers, or powers of ten
p_log_num <- sg_gof_res(
  fpath_i = fpath,
  vs_time = FALSE,
  log_x = TRUE,
  log_labs = "number",
  addline = FALSE
)
p_log_math <- sg_gof_res(
  fpath_i = fpath,
  vs_time = FALSE,
  log_x = TRUE,
  log_labs = "math",
  addline = FALSE
)
p_log_num
p_log_math

# CWRES is preferred over WRES when both columns exist
smrg_cw <- read_smrg_obj(fpath)
smrg_cw$SDTAB$CWRES <- smrg_cw$SDTAB$WRES
p_cwres <- sg_gof_res(
  smrg_cw,
  indiv = FALSE,
  weighted = TRUE,
  vs_time = TRUE,
  addline = FALSE,
  smooth = FALSE
)
p_cwres

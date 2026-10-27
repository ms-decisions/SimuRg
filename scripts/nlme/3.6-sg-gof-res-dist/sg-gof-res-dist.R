devtools::load_all()

fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")

p_iwres <- sg_gof_res_dist(fpath_i, res_type = "IWRES")
p_iwres

p_sex <- sg_gof_res_dist(
  fpath_i,
  res_type = "IWRES",
  cov_cols = "SEX",
  facet_i = "SEX"
)
p_sex

p_age <- sg_gof_res_dist(
  fpath_i,
  res_type = "IWRES",
  cov_cols = "AGE",
  facet_i = "AGE",
  n_quantiles = 3
)
p_age

p_qq <- sg_gof_res_dist(fpath_i, res_type = "RES", plot_type = "QQ")
p_qq

# Warfarin SDTAB has WRES, not CWRES. CWRES is used the same way when the column exists.
p_wres <- sg_gof_res_dist(
  fpath_i,
  res_type = "WRES",
  cov_cols = c("SEX", "AGE"),
  facet_i = "SEX"
)
p_wres

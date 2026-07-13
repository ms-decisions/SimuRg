## Author: Ugolkov Yaroslav
## First created: 2026-08-07
## Description: Perform simulations for prediction distribution plots
## Keywords: SimuRg, monolix, converter


devtools::load_all()
devtools::document()

library(rxode2)
fpath_i <- system.file("extdata", "simurg_object", "Warfarin_PK.RData", package = "SimuRg")
mod_fin <- rxode2({
  # Differential equations
  ka = ka_pop*exp(omega_ka)
  V = V_pop
  Cl = Cl_pop*exp(omega_Cl)
  d/dt(Ad) = -ka * Ad
  d/dt(Ac) = ka * Ad - Cl/V * Ac
  # Concentration calculations
  Cc = Ac / V
})
obj1$SUMTAB$PAR
predist_sim <- sg_predist_sim(fpath_i, model =mod_fin, outputs = "Cc")
sg_predist_vis(fpath_i,predist_sim)
?sg_predist_vis
?sg_predist_sim

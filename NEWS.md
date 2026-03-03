# SimuRg 0.1.4(2026-02-24)
## Documentation
* Examples in documentation fix
* Warfarin, sim_timeprof and data_pbc datasets documentation addition

## Bug fixes
* `sg_fit()` now waits for `summary.txt` filecreation, before calling `sg_converter()` 

# SimuRg 0.1.3 (2026-02-17)

## NEW FEATURES
* 8 Goodness of Fit functions: `sg_parsum()`, `sg_modcomp()`,
`sg_gof_obpr()`, `sg_gof_res()`, `sg_gof_par_dist()`, `sg_gof_res_dist()`,
`sg_gof_par_cov()`, `sg_gof_tp()`
* Function to perform simulations (rxode2-based solver): `sg_sim()`, `sg_sim_tp()`
* Functions for simulation-based analysis: `sg_vpc_sim()`, `sg_vpc_vis()`, `sg_preddist_sim()`
* Model fit and model building functions(Monolix fitter): `sg_fit()`, `sg_modbuild()`
* Function to convert projects across different systems. In the current version, 
conversion from Monolix(2023) to GCO (Simurg) was realised: `sg_converter()`
* Function to generate virtual populations based on distributions, or a dataset: `sg_vpop_est()`


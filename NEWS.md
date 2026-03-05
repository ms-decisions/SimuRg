# SimuRg 0.1.5(2026-03-03)
## New features
* `sg_globalsens_sim()` function for global sensetivity analysis was added to the package
*  Added possibility use `sg_fit()` with several observations

## Bug fixes
* `sg_fit()` error was fixed. The function will work with different residual error types.
  Now it will wait for `summary.txt` file appearance while fitting with Monolix
* Usage of `f_scales`, `lab_x`, `lab_y`, `facet_i`, `cov_cols`arguments were fixed for 
  `sg_gof_obpr()` function.
* Visualization in the `sg_gof_par_cov()` function was improved

## Documentation
* `sg_fit()` documentation was improved due to the commits of the users

# SimuRg 0.1.4(2026-02-24)

## Bug fixes
* `Bugs in examples of documentation were fixed

## Documentation
* Examples in documentation fix
* Warfarin, sim_timeprof and data_pbc datasets documentation addition


# SimuRg 0.1.3 (2026-02-17)

## New features
* 8 Goodness of Fit functions: `sg_parsum()`, `sg_modcomp()`,
`sg_gof_obpr()`, `sg_gof_res()`, `sg_gof_par_dist()`, `sg_gof_res_dist()`,
`sg_gof_par_cov()`, `sg_gof_tp()`
* Function to perform simulations (rxode2-based solver): `sg_sim()`, `sg_sim_tp()`
* Functions for simulation-based analysis: `sg_vpc_sim()`, `sg_vpc_vis()`, `sg_preddist_sim()`
* Model fit and model building functions(Monolix fitter): `sg_fit()`, `sg_modbuild()`
* Function to convert projects across different systems. In the current version, 
conversion from Monolix(2023) to GCO (Simurg) was realised: `sg_converter()`
* Function to generate virtual populations based on distributions, or a dataset: `sg_vpop_est()`


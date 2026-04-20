# SimuRg 0.1.9 (2026-04-09)

## New features:
* `sg_covsens_vis()` function was added to the project with tests and documentation
* `sg_converter()` now works with two DVIDs

## Bug fixes:
* `sg_converter()` and `sg_vpc_sim()` : fixed bugs reported in [SD-1360]
* `sg_vpc_vis()`: fixed bugs reported in [SD-1361]

## Documentation:
* `sg_vpc_vis()` documentation was fixed with the bug [SD-1361]

# SimuRg 0.1.8 (2026-04-03)

## New features 
* `sg_converter()` now returns two objects: GFO - generalized fit object with fit results and GCO - generalized control object, that contain information about the fit options and see the 

## Bug fixes 
* `sg_gof_res()` was fixed with the [SD-1391] support request
* `sg_covsens_sim()` was tested on additional projects

## Documentation 
* all functions arguments descriptions now starts with a capital letter

# SimuRg 0.1.7 (2026-03-23)

## New features
* `sg_globalsens_vis()` function was added to the package


## Bug fixes 
* `sg_converter()` error with correlations between random effects [SD-1340] has been solved
* `sg_converter()` error duplicated DVID columns was solved [SD-1354]. Now the package compares the column
  names in the initial dataset and the ones specified in headers to resolve possible conflicts.
* `sg_converter()` error in dataset column names resolved [SD-1355]. Now the columns of the 
  initial dataset are renamed after loading: all special signs in the column names 
  are replaced with `_`

## Documentation
* examples in `sg_covsens_sim()` was updated

# SimuRg 0.1.6 (2026-03-19)

## New features
* `sg_multistart()` function for model assessment was added to the package
* `sg_covsensns_sim()` function was added for covariate sensitivity
* `sg_sim()` function was extended to work with multiple simulations scenarios
* `sg_translator()` function for model syntax translation from rxode to mlxtran and vice-versa was added
* `sg_converter()` function now works with Monolix 2023 and 2024 versions
* Example of GFO was updated in the package

## Bug fixes 
* `sg_pgof_par_dist()` errors with correlation plots was fixed

## Documentation
* package description is now available after calling `help(SimuRg)`
* examples in `sg_vpc_vis()` was updated

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


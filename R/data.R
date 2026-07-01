## Author: Mikhailova Anna
## First created: 2025-08-25
## Description: functions for common internal data structures generation
## Keywords: SimuRg, internal data

MSDcol <- c("#1a1866", "#f2b93b", "#b73b58", "#a2d620", "#14D98E", "#9c4ec7",
            "#3a6eba", "#efdd3c", "#69686d", '#844538', '#D91477', '#F3A9FF')

funSum_sim <- list(mean   = ~mean(., na.rm = TRUE),
                   median = ~median(., na.rm = TRUE),
                   min    = ~min(., na.rm = TRUE),
                   max    = ~max(., na.rm = TRUE),
                   sd     = ~sd(., na.rm = TRUE),
                   P025   = ~quantile(., 0.025, na.rm = TRUE),
                   P05    = ~quantile(., 0.05, na.rm = TRUE),
                   P10    = ~quantile(., 0.10, na.rm = TRUE),
                   P15   = ~quantile(., 0.15, na.rm = TRUE),
                   P25    = ~quantile(., 0.25, na.rm = TRUE),
                   P75    = ~quantile(., 0.75, na.rm = TRUE),
                   P85   = ~quantile(., 0.85, na.rm = TRUE),
                   P90    = ~quantile(., 0.90, na.rm = TRUE),
                   P95    = ~quantile(., 0.95, na.rm = TRUE),
                   P975   = ~quantile(., 0.975, na.rm = TRUE),
                   geom_mean = ~exp(mean(log(.[. > 0]), na.rm = TRUE)),
                   CV     = ~sd(., na.rm = TRUE)/mean(., na.rm = TRUE)*100)

#' Generalized Control Object (GCO)
#'
#' Generalized control object (GCO) stores the model setup, dataset mapping,
#' and estimation options used for fitting or downstream simulation configuration.
#'
#' Create a GCO with [sg_converter()] from a Monolix project, or with [sg_fit()].
#'
#' @format ## `GCO`
#' A named list with the following components:
#' \describe{
#'   \item{`headers`}{List of lists. One element per dataset column with fields:
#'     `name` (character, column name),
#'     `use` (character, Monolix column role such as `"identifier"`, `"time"`,
#'     `"observation"`, `"covariate"`, etc.),
#'     `type` (character or empty; for covariates: `"continuous"` or `"categorical"`).}
#'   \item{`data`}{Character. Path to the source dataset file.}
#'   \item{`model`}{Character. Path to the structural model file (Monolix syntax).}
#'   \item{`task_opt`}{List. Task options passed to the fitter; empty list when not set.}
#'   \item{`covs`}{Character vector. Names of covariate columns detected in the dataset
#'     headers.}
#'   \item{`project_name`}{Character. Project or run name.}
#'   \item{`theta`}{List of lists or data frame. Population (typical) parameter
#'     definitions. Each record contains:
#'     `NAME` (character, parameter name without `_pop` suffix),
#'     `TRANS` (character, `"normal"`, `"logNormal"`, or `"logitNormal"`),
#'     `INIT` (numeric, initial or fixed value on the natural scale),
#'     `EST` (logical, whether the parameter is estimated). When passed to
#'     [sg_fit()], a tibble with optional `LB` and `UB` for logit bounds is also
#'     accepted.}
#'   \item{`ruv`}{List, or list of lists for multiple observation types. Residual
#'     error model definition with fields:
#'     `YNAME` (character, longitudinal output name, e.g. `"y1"`),
#'     `DVID` (numeric, observation-type identifier),
#'     `TRANS` (character, residual distribution),
#'     `PRED` (character, prediction variable name),
#'     `ERR` (character, error model type such as `"constant"`, `"proportional"`,
#'     `"combined1"`),
#'     `INIT` (numeric vector of initial residual-error parameter values),
#'     `EST` (logical vector, estimation flags),
#'     `BLQM` (optional BLQ handling method; `NULL` when not used).}
#'   \item{`re`}{List with `init` and `est` numeric matrices. Between-subject
#'     variability: diagonal elements are initial variances of `omega_*` parameters,
#'     off-diagonal elements are covariances; `est` uses `TRUE`/`FALSE`/`NA` in
#'     the same layout as [sg_fit()] input.}
#'   \item{`occ`}{List with `init` and `est` matrices. Between-occasion variability
#'     in the same layout as `re`; all-zero/`NA` when occasion variability is not
#'     modeled.}
#'   \item{`modelText`}{Character. Full text of the structural model file.}
#' }
#' @examples
#' \donttest{
#' # Bundled GCO (Simurg fit configuration, JSON format)
#' gco_path <- system.file("extdata", "simurg_object", "sg_object.json",
#'                         package = "SimuRg")
#' gco <- read_smrg_ctrl(gco_path)
#' names(gco)
#' gco$project_name
#' gco$theta
#' gco$ruv
#' gco$re$init
#' }
#' @seealso [GFO], [GMO], [sg_fit()], [sg_converter()], [read_smrg_ctrl()]
#' @name GCO
NULL

#' Generalized Fit Object (GFO)
#'
#' Generalized fit object (GFO) stores population and individual estimation results,
#' diagnostic tables, and variance–covariance matrices produced by model fitting.
#'
#' Create a GFO with [sg_converter()] from a Monolix project, or with [sg_fit()]
#' when `fit = TRUE`.
#'
#' @format ## `GFO`
#' A named list with the following components:
#' \describe{
#'   \item{`SDTAB`}{Data frame. Standard diagnostic table of observations and
#'     model predictions. Typical columns:
#'     `ID`, `TIME`, `DV` (observed), `DVID`, `DVNAME`,
#'     `PRED` (population prediction), `IPRED` (individual prediction),
#'     `RES`, `IRES`, `WRES`, `IWRES` (residuals),
#'     `MDV`, and optionally `OCC`, `BLQ`, `CENS`, `LIMIT`.}
#'   \item{`SUMTAB`}{Data frame. Parameter summary table. Typical columns:
#'     `PAR` (parameter name), `VALUE` (estimate),
#'     `TYPE` (`"Typical values"`, `"Random effects"`, `"Covariate coefficients"`,
#'     `"Correlation coefficients"`, or `"Residual error model"`),
#'     `DISTRIBUTION`, `EST` (`"ESTIMATED"` or `"FIXED"`),
#'     `SE`, `RSE`, `CV`, `LCI`, `UCI`,
#'     `ETAshrinkage_var`, `ETAshrinkage_sd`, `EPSshrinkage_sd`.}
#'   \item{`SIGMAMAT`}{Numeric matrix. Residual error variance–covariance matrix
#'     (diagonal elements are residual variance parameters).}
#'   \item{`OMEGAMAT`}{Numeric matrix. Inter-individual variability
#'     variance–covariance matrix on the `omega_*` scale.}
#'   \item{`OCCMAT`}{Numeric matrix. Inter-occasion variability matrix; empty when
#'     occasion variability is not present.}
#'   \item{`EVTAB`}{Data frame. Dosing/event table extracted from the dataset.
#'     Typical columns: `ID`, `TIME`, `EVID`, `CMT`, `ADM`, `AMT`, and optionally
#'     `OCC`, `ADDL`, `II`, `DUR`, `TINF`, `RATE`, `SS`.}
#'   \item{`PATAB`}{Data frame. Individual random effects and post hoc parameter
#'     estimates. Columns `ID`, `eta_*` (individual ETAs), and individual PK/PD
#'     parameter columns (e.g. `ka`, `Vd`, `CL`).}
#'   \item{`COTAB`}{Data frame. Continuous covariates: `ID` plus continuous
#'     covariate columns (including derived columns such as log-transformed
#'     covariates). Empty data frame when none are present.}
#'   \item{`CATAB`}{Data frame. Categorical covariates: `ID` plus categorical
#'     covariate columns. Empty data frame when none are present.}
#'   \item{`REGTAB`}{Data frame. Time-varying regressors (`ID`, `TIME`, regressor
#'     columns). Empty data frame when none are present.}
#'   \item{`OFV`}{Data frame with one row. Model fit criteria parsed from Monolix
#'     `summary.txt`: `LL` (minus 2 times log-likelihood, Monolix OFV), `AIC`,
#'     `BIC`, `BICc`.}
#'   \item{`COVMAT`}{Numeric matrix. Variance–covariance matrix of population
#'     parameter estimates.}
#'   \item{`CORRMAT`}{Numeric matrix. Correlation matrix of population parameter
#'     estimates.}
#'   \item{`OPTIONS`}{List or `NULL`. Additional model or task options when
#'     available.}
#'   \item{`PROJNAME`}{Character. Project or run name.}
#' }
#' @examples
#' \donttest{
#' gfo <- read_smrg_obj(gfo4cov)
#' names(gfo)
#' head(gfo$SDTAB)
#' head(gfo$SUMTAB)
#' head(gfo$PATAB)
#' gfo$OFV
#' gfo$OMEGAMAT
#' }
#' @seealso [GCO], [sg_fit()], [sg_converter()], [read_smrg_obj()], [gfo4cov]
#' @name GFO
NULL

#' Generalized Model Object (GMO)
#'
#' Generalized model object (GMO) is an RxODE2 model that encodes the full
#' pharmacometric model in simulation syntax: structural parameters, covariate
#' effects, inter-individual variability, and residual error components.
#'
#' Build a GMO with `rxode2::rxode2()` or use the bundled [gmo_pk1c] example.
#'
#' @format ## `GMO`
#' An RxODE2 model object (`rxode2` class). Typical contents:
#' \describe{
#'   \item{Typical values}{Fixed effects on the log scale (`ka_pop`, `Vd_pop`,
#'     `CL_pop`, …) with transformation lines such as `ka_tv = exp(ka_pop)`.}
#'   \item{Random effects}{Initial values for ETAs (`omega_ka`, `omega_Vd`,
#'     `omega_CL`, …) and individual parameters such as
#'     `ka = ka_tv * exp(omega_ka)`.}
#'   \item{Covariate effects}{Covariate coefficients (e.g. `beta_Vd_WT`) applied
#'     to typical values, e.g. `Vd = Vd_tv * exp(omega_Vd) * exp(beta_Vd_WT * WT)`.}
#'   \item{Structural model}{Initial conditions and ODE equations (`Ad`, `Ac`,
#'     `Cc`, …).}
#'   \item{Residual error}{Error model parameters (e.g. `Cc_b`) and observation
#'     outputs such as `Cc_ResErr = Cc * (1 + Cc_b)`.}
#' }
#' @examples
#' \donttest{
#' class(gmo_pk1c)
#' gmo_pk1c$model
#' gmo_pk1c$params
#' }
#' @seealso [GCO], [GSI], [gmo_pk1c], [sg_sim()]
#' @name GMO
NULL

#' One-compartment PK generalized model object (GMO)
#'
#' Bundled [GMO]: one-compartment oral PK model with inter-individual variability
#' on `ka`, `Vd`, and `CL`, linear covariate effect of `WT` on `Vd`
#' (`beta_Vd_WT`), and proportional residual error on `Cc_ResErr`.
#'
#' @format ## `gmo_pk1c`
#' An `rxode2` model object with population parameters (`*_pop`), random-effect
#' placeholders (`omega_*`), covariate effect `beta_Vd_WT`, structural ODEs
#' (`Ad`, `Ac`, `Cc`), and residual error output (`Cc_ResErr`).
"gmo_pk1c"

#' Generalized simulations input (GSI) for one-compartment PK example
#'
#' Bundled [GSI] matching the standard simulation scenario: two dosing regimens
#' with different `WT`, population parameter uncertainty (`thetamat`, `npop = 10`),
#' and output `Cc` only.
#'
#' @format ## `gsi_pk1c`
#' A named list ready for [sg_sim()], with components:
#' \describe{
#'   \item{`et`}{Data frame. Event table (`ID`, `TIME`, `AMT`, `ADDL`, `II`, `WT`).}
#'   \item{`stimes`}{Numeric vector. Output time grid (0–120 h).}
#'   \item{`outputs`}{Character vector. Model outputs to retain (`"Cc"`).}
#'   \item{`covs`}{Character vector. Covariate columns joined from `et` (`"WT"`).}
#'   \item{`theta`}{Named numeric vector. Population parameters (`ka_pop`, `Vd_pop`,
#'     `CL_pop`, `beta_Vd_WT`).}
#'   \item{`thetamat`}{Numeric matrix. Parameter estimation covariance (8 x 8).}
#'   \item{`omega`, `sigma`}{`NULL` when inter-individual and residual variability
#'     are not resampled.}
#'   \item{`npop`, `nsub`}{Integers. Numbers of population (`10`) and subject (`1`) replicates.}
#'   \item{`addcov`}{Logical. `FALSE`; covariates are already in `et`.}
#' }
"gsi_pk1c"

#' Generalized Simulations Input (GSI)
#'
#' Generalized simulations input (GSI) is the set of arguments passed to [sg_sim()]
#' to define a simulation scenario: model, event table, parameter values, and
#' optional uncertainty or variability matrices. See bundled [gsi_pk1c] for a
#' complete example.
#'
#' @format ## `GSI`
#' A named list (not a formal S3 class) with components matching [sg_sim()]
#' arguments. When serialized to JSON, the same information may use `output`
#' instead of `outputs` and store `model` as a character string; in R, pass an
#' [GMO] object to [sg_sim()].
#' \describe{
#'   \item{`model`}{[GMO]. RxODE2 model (`gmo_pk1c` in bundled examples).}
#'   \item{`et`}{Data frame. Event table; columns such as `ID`, `TIME`, `AMT`,
#'     `ADDL`, `II`, and covariates (e.g. `WT`).}
#'   \item{`stimes`}{Numeric vector. Output sampling time grid.}
#'   \item{`outputs`}{Character vector. Model output names to retain (e.g. `"Cc"`).}
#'   \item{`covs`}{Character vector. Covariate column names taken from `et`.}
#'   \item{`theta`}{Named numeric vector. Population parameters on the model scale
#'     (`ka_pop`, `Vd_pop`, `CL_pop`, covariate coefficients, …).}
#'   \item{`thetamat`}{Numeric matrix. Parameter-estimation covariance for
#'     population uncertainty (`npop` > 1).}
#'   \item{`omega`, `sigma`}{Matrix or `NULL`. Inter-individual and residual
#'     variability; empty when not resampled.}
#'   \item{`npop`, `nsub`}{Integers. Population and subject replicate counts.}
#'   \item{`byID`, `byPOP`, `shared`, `aggr`, `addcov`, `keep`, …}{Optional
#'     [sg_sim()] controls.}
#' }
#' @examples
#' \donttest{
#' names(gsi_pk1c)
#' head(gsi_pk1c$et)
#' gsi_pk1c$theta
#' gsi_pk1c$npop
#' # Pass to sg_sim() together with the GMO:
#' do.call(sg_sim, c(list(model = gmo_pk1c), gsi_pk1c))
#' }
#' @seealso [GSO], [GMO], [gmo_pk1c], [gsi_pk1c], [sg_sim()]
#' @name GSI
NULL

#' Generalized Simulations Output (GSO)
#'
#' Generalized simulations output (GSO) is the long-format data frame returned
#' by [sg_sim()] and related simulation helpers. Run the bundled [gsi_pk1c]
#' scenario with [gmo_pk1c] to reproduce the standard example output.
#'
#' @format ## `GSO`
#' A data frame in long format with one row per `ID`, replicate, time point, and
#' output variable. With population uncertainty (`thetamat`, `npop` > 1) typical
#' columns are:
#' \describe{
#'   \item{`ID`}{Dosing scenario identifier from the event table.}
#'   \item{`POPN`}{Population replicate index (present when `thetamat` is used).}
#'   \item{`sim.id`}{Subject replicate index (present when `omega` is used).}
#'   \item{`TIME`}{Simulation time (h).}
#'   \item{`VAR`}{Output variable name (e.g. `"Cc"`).}
#'   \item{`VALUE`}{Simulated value.}
#' }
#' When `aggr` is set, aggregation statistic columns (`mean`, `median`, `P05`,
#' …) replace or supplement `VALUE`; see the Details section of [sg_sim()].
#' @examples
#' \donttest{
#' gso <- do.call(sg_sim, c(list(model = gmo_pk1c), gsi_pk1c))
#' names(gso)
#' head(gso)
#' subset(gso, VAR == "Cc" & ID == 1 & POPN == 1)
#' }
#' @seealso [GSI], [GMO], [gmo_pk1c], [gsi_pk1c], [sg_sim()]
#' @name GSO
NULL

#' Warfarin PKPD dataset
#'
#' Generated dataset with warfarin PK/PD measurements and covariates for 100 patients.
#'
#' @format ## `warfarin`
#' A data frame with 1700 rows and 17 columns:
#' \describe{
#'   \item{ID}{identifier of the individual}
#'   \item{TIME}{time of the dose or observation record (nominal)}
#'   \item{DV}{records the measurement data}
#'   \item{DVID}{identifier for the observation type (to distinguish different types of observations, e.g PK and PD)}
#'   \item{DVNAME}{name for the observation type (to distinguish different types of observations, e.g PK and PD)}
#'   \item{CMT}{compartment of the event}
#'   \item{ADM}{identifier for the type of dose}
#'   \item{AMT}{dose amount}
#'   \item{EVID}{identifier to indicate if the line is a dose-line or a response-line}
#'   \item{MDV}{identifier to ignore the observation information of that line}
#'   \item{AGE}{age of the individual}
#'   \item{SEX}{sex of the individual}
#'   \item{WEIGHT}{weight of the individual}
#'   \item{BMI}{BMI of the individual}
#'   \item{CLCR}{creatinine clearance of the individual}
#'   \item{CYP2C9_gentyp}{CYP2C9 genotype of the individual}
#'   \item{VKORC1_gentyp}{VKORC1 genotype of the individual}
#' }
"warfarin"

#' Primary Biliary Cirrhosis (PBC) dataset
#'
#' Clinical trial dataset with Primary Biliary Cirrhosis patient data including survival outcomes and covariates.
#'
#' @format ## `data_pbc`
#' A data frame with 280 rows and 15 columns:
#' \describe{
#'   \item{id}{identifier of the individual}
#'   \item{years}{time in years (survival time)}
#'   \item{status2}{survival status (0 = censored, 1 = event)}
#'   \item{drug}{treatment group ("D-penicil" or "placebo")}
#'   \item{age}{age of the individual in years}
#'   \item{sex}{sex of the individual ("male" or "female")}
#'   \item{ascites}{presence of ascites ("Yes" or "No")}
#'   \item{hepatomegaly}{presence of hepatomegaly ("Yes" or "No")}
#'   \item{spiders}{presence of spider angiomas ("Yes" or "No")}
#'   \item{edema}{edema status ("No edema", "edema no diuretics", or "edema despite diuretics")}
#'   \item{serBilir}{serum bilirubin level (mg/dl)}
#'   \item{serChol}{serum cholesterol level (mg/dl)}
#'   \item{albumin}{albumin level (g/dl)}
#'   \item{platelets}{platelet count (per cubic mm/1000)}
#' }
"data_pbc"

#' Warfarin population PK parameter estimates
#'
#' Final parameter estimates from a warfarin population pharmacokinetic model,
#' including typical values, covariate effects, random effects, and residual
#' error parameters with associated uncertainty measures.
#'
#' @format ## `parest`
#' A data frame with 15 rows and 12 columns:
#' \describe{
#'   \item{parameter}{parameter name as used in the model}
#'   \item{value}{final parameter estimate}
#'   \item{TYPE}{parameter category:`"Typical values"`, `"Covariate effects"`, `"Random effects"`, or `"Residual error model"`}
#'   \item{EST}{estimation status (`"ESTIMATED"`)}
#'   \item{SE}{standard error of the estimate}
#'   \item{RSE}{relative standard error (%)}
#'   \item{LCI}{lower bound of the 95 % confidence interval}
#'   \item{UCI}{upper bound of the 95 % confidence interval}
#'   \item{ETAshrinkage_var}{ETA shrinkage on the variance scale (%); `NA` for non-random-effect parameters}
#'   \item{ETAshrinkage_sd}{ETA shrinkage on the SD scale (%); `NA` for non-random-effect parameters}
#'   \item{EPSshrinkage_sd}{EPS shrinkage on the SD scale (%); `NA` for non-residual-error parameters}
#' }
#'@source Derived from \code{par_fin_i.csv} — warfarin PopPK model estimation output.
"parest"

#' Warfarin covariate values dataset
#'
#' Individual covariate values for 100 patients used in the warfarin
#' covariate sensitivity simulation, including demographic and pharmacogenomic variables.
#'
#' @format ## `ds_covval`
#' A data frame with 100 rows and 9 columns:
#' \describe{
#'   \item{ID}{individual patient identifier}
#'   \item{AGE}{age in years}
#'   \item{WEIGHT}{body weight in kg}
#'   \item{BMI}{body mass index (kg/m\ifelse{html}{\out{<sup>2</sup>}}{\eqn{^2}})}
#'   \item{LG_WEIGHT}{log-transformed body weight (log10 scale), used as continuous covariate in the PK model}
#'   \item{LG_AGE}{log-transformed age (log10 scale), used as continuous covariate in the PK model}
#'   \item{SEX}{sex of the individual (0 = female, 1 = male)}
#'   \item{CYP2C9}{CYP2C9 genotype encoded as integer: 0 = \emph{*1/*1} (reference), 1 = \emph{*1/*2}, 2 = \emph{*1/*3}, 3 = \emph{*2/*2}, 4 = \emph{*2/*3}, 5 = \emph{*3/*3}}
#'   \item{VKORC1}{VKORC1 genotype (\code{"GG"}, \code{"AG"}, or \code{"AA"})}
#' }
#'@source Derived from \code{data_fin_i.csv} — individual covariate dataset for the warfarin PopPK covariate sensitivity analysis.
"ds_covval"

#' Warfarin PopPK model fitting output (4-covariate model)
#'
#' A list containing the full output of a warfarin population pharmacokinetic
#' model fit with four covariates (SEX, LG_WEIGHT, CYP2C9, VKORC1), including
#' parameter estimates, model diagnostics, individual predictions, and covariate data.
#'
#' @format ## `gfo4cov`
#' A named list with 12 elements:
#' \describe{
#'   \item{PROJNAME}{character string identifying the project run (\code{"run_4cov"})}
#'   \item{SUMTAB}{data frame with 15 rows and 11 columns containing the parameter
#'     summary table. Columns: \code{PAR} (parameter name), \code{VALUE} (estimate),
#'     \code{TYPE} (category: \code{"Typical values"}, \code{"Covariate effects"},
#'     \code{"Random effects"}, or \code{"Residual error model"}), \code{EST}
#'     (estimation status), \code{SE} (standard error), \code{RSE} (relative SE, \%),
#'     \code{LCI} and \code{UCI} (95\% confidence interval bounds),
#'     \code{ETAshrinkage_var}, \code{ETAshrinkage_sd}, \code{EPSshrinkage_sd} (shrinkage metrics).}
#'   \item{SDTAB}{data frame with 1600 rows and 11 columns — the standard table of
#'     observations and model predictions. Columns: \code{ID}, \code{TIME}, \code{DV}
#'     (observed), \code{DVID}, \code{PRED} (population prediction), \code{IPRED}
#'     (individual prediction), \code{RES} (residual), \code{IRES} (individual residual),
#'     \code{WRES} (weighted residual), \code{IWRES} (individual weighted residual), \code{MDV}.}
#'   \item{PATAB}{data frame with 100 rows and 7 columns — individual parameter table.
#'     Columns: \code{ID}, \code{eta_ka}, \code{eta_Vd}, \code{eta_CL}
#'     (individual ETA random effects), \code{ka}, \code{Vd}, \code{CL}
#'     (individual post-hoc PK parameter estimates).}
#'   \item{COVMAT}{15 \eqn{\times} 15 numeric matrix — covariance matrix of the
#'     population parameter estimates.}
#'   \item{CORRMAT}{15 \eqn{\times} 15 numeric matrix — correlation matrix of the
#'     population parameter estimates.}
#'   \item{OFV}{data frame with 1 row and 4 columns — model fit criteria:
#'     \code{LL} (\code{-2 x log-likelihood}, i.e. Monolix OFV), \code{AIC}, \code{BIC}, \code{BICc}.}
#'   \item{OMEGAMAT}{3 \eqn{\times} 3 numeric matrix — OMEGA variance-covariance
#'     matrix of the inter-individual random effects (\code{eta_ka}, \code{eta_Vd}, \code{eta_CL}).}
#'   \item{SIGMAMAT}{1 \eqn{\times} 1 numeric matrix — SIGMA residual error
#'     variance matrix.}
#'   \item{EVTAB}{data frame with 100 rows and 6 columns — dosing/event table.
#'     Columns: \code{ID}, \code{TIME}, \code{EVID}, \code{AMT}, \code{ADM}, \code{CMT}.}
#'   \item{COTAB}{data frame with 100 rows and 6 columns — continuous covariate table.
#'     Columns: \code{ID}, \code{AGE}, \code{WEIGHT}, \code{BMI},
#'     \code{LG_WEIGHT} (log10-transformed weight), \code{LG_AGE} (log10-transformed age).}
#'   \item{CATAB}{data frame with 100 rows and 4 columns — categorical covariate table.
#'     Columns: \code{ID}, \code{SEX} (0 = female, 1 = male),
#'     \code{CYP2C9} (genotype integer code, 0–5),
#'     \code{VKORC1} (genotype: \code{"GG"}, \code{"AG"}, or \code{"AA"}).}
#' }
#' @source Generated from the warfarin PopPK model fitting run \code{"run_4cov"} used
#' in the covariate sensitivity analysis workflow.
"gfo4cov"

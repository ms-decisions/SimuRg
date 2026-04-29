## Author: Mikhailova Anna
## First created: 2025-08-25
## Description: functions for common internal data structures generation
## Keywords: SimuRg, internal data

MSDcol <- c("#1a1866", "#f2b93b", "#b73b58", "#a2d620", "#14D98E", "#9c4ec7",
            "#3a6eba", "#efdd3c", "#69686d", '#844538', '#D91477', '#F3A9FF')

funSum_sim <- list(mean   = ~mean(., na.rm = T),
                   median = ~median(., na.rm = T),
                   min    = ~min(., na.rm = T),
                   max    = ~max(., na.rm = T),
                   sd     = ~sd(., na.rm = T),
                   P025   = ~quantile(., 0.025, na.rm = T),
                   P05    = ~quantile(., 0.05, na.rm = T),
                   P10    = ~quantile(., 0.10, na.rm = T),
                   P15   = ~quantile(., 0.15, na.rm = T),
                   P25    = ~quantile(., 0.25, na.rm = T),
                   P75    = ~quantile(., 0.75, na.rm = T),
                   P85   = ~quantile(., 0.85, na.rm = T),
                   P90    = ~quantile(., 0.90, na.rm = T),
                   P95    = ~quantile(., 0.95, na.rm = T),
                   P975   = ~quantile(., 0.975, na.rm = T),
                   geom_mean = ~exp(mean(log(.[. > 0]), na.rm = T)),
                   CV     = ~sd(., na.rm = T)/mean(., na.rm = T)*100)

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
#'     \code{LL} (log-likelihood), \code{AIC}, \code{BIC}, \code{BICc}.}
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

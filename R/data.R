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

#' Simulated time profile dataset
#'
#' Simulated pharmacokinetic time profile data with concentration measurements and covariates.
#'
#' @format ## `sim_timeprof`
#' A data frame with 2000 rows and 9 columns:
#' \describe{
#'   \item{ID}{identifier of the individual}
#'   \item{TIME}{time points}
#'   \item{VAR}{variable name (type of measurement)}
#'   \item{VALUE}{concentration value}
#'   \item{Regimen}{dosing regimen description}
#'   \item{Dose}{dose amount}
#'   \item{WTBL}{weight baseline}
#'   \item{CAT1}{categorical covariate}
#' }
"sim_timeprof"

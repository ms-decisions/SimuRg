## Author: Ugolkov Yaroslav, Victoria Kulesh
## First created: 2025-10-24
## Description: Convert Monolix project results to SimuRg objects
## Keywords: SimuRg, monolix, converter

#' Convert Monolix project output to R objects
#'
#' @description
#' Reads and parses Monolix project output files (.mlxtran and associated data files)
#' into a structured list of SimuRg objects including parameter estimates, individual predictions,
#' residuals, and diagnostic information.
#'
#' @param folder_path Character string. Path to the directory containing Monolix project files.
#' @param proj_name Character string. Name of the Monolix project (without file extension).
#' @param save_file Logical. If \code{TRUE}, saves \code{GCO} and \code{GFO} JSON files
#'   to \code{folder_path} with names \code{<proj_name>_GCO.json} and \code{<proj_name>_GFO.json},
#'   and also saves two RData files: \code{<proj_name>_GCO.RData} (object \code{gco})
#'   and \code{<proj_name>_GFO.RData} (object \code{gfo}).
#'
#' @return
#' Returns a list with the following components:
#' \itemize{
#'   \item \code{GFO}: SimuRg generalized fit object output object with:
#'   \itemize{
#'     \item \code{SDTAB}: A tibble containing simulation data with columns
#'     \item \code{SUMTAB}: A tibble with parameter summary statistics containing
#'     \item \code{SIGMAMAT}: Residual variability matrix
#'     \item \code{OMEGAMAT}: Inter-individual variability matrix
#'     \item \code{OCCMAT}: Inter-occasion variability matrix (NA if not present)
#'     \item \code{EVTAB}: A tibble with event information
#'     \item \code{PATAB}: A tibble with individual parameter estimates
#'     \item \code{COTAB}: A tibble with continuous covariates
#'     \item \code{CATAB}: A tibble with categorical covariates
#'     \item \code{REGTAB}: Regression parameters (empty data.frame if not present)
#'     \item \code{OFV}: A tibble with objective function values
#'     \item \code{COVMAT}: Variance-covariance matrix of parameter estimates
#'     \item \code{CORRMAT}: Correlation matrix of parameter estimates
#'     \item \code{OPTIONS}: Model options (NULL if not present)
#'     \item \code{PROJNAME}: Project name
#'   }
#'   \item \code{GCO}: SimuRg generalized control object parsed from the mlxtran project with:
#'   \itemize{
#'     \item \code{headers}: List of dataset column descriptors (\code{name}, \code{use}, \code{type})
#'     \item \code{data}: Path to source data file
#'     \item \code{model}: Path to model file
#'     \item \code{task_opt}: Task options placeholder (empty object)
#'     \item \code{covs}: Covariate names detected
#'     \item \code{project_name}: Project name
#'     \item \code{theta}: List of population parameter definitions (\code{NAME}, \code{INIT}, \code{EST}, \code{TRANS})
#'     \item \code{ruv}: Residual error model definition (\code{YNAME}, \code{DVID}, \code{TRANS}, \code{PRED}, \code{ERR}, \code{INIT}, \code{EST})
#'     \item \code{re}: Between-subject variability matrices (\code{init}, \code{est})
#'     \item \code{occ}: Between-occasion variability matrices (\code{init}, \code{est})
#'     \item \code{modelText}: Text content of the model file
#'   }
#' }
#'
#' @details
#' This function serves as a bridge between Monolix output and R by parsing the .mlxtran file
#' and associated data files to create a comprehensive R object containing all relevant model
#' outputs. This facilitates further analysis, visualization, and reporting in R.
#'
#' The function automatically detects and imports various components of Monolix output
#' including population parameters, individual parameters, covariates, and diagnostic
#' metrics.
#'
#' If \code{save_file = TRUE}, the function additionally writes \code{GCO} and \code{GFO}
#' JSON files and \code{.RData} files to \code{folder_path}.
#'
#' @examples
#' \donttest{
#' library(stringr)
#' # Convert Monolix project results
#' test_folder <- system.file("extdata", "Monolix_objects", package = "SimuRg")
#' if (substr(test_folder, nchar(test_folder), nchar(test_folder)) != "/")
#'   test_folder <- str_c(test_folder, "/")
#' pro_name <- "proj-solo"
#' result <- sg_converter(folder_path = test_folder, proj_name = pro_name)
#' # save(results, file = "./models/simurg_object/Warfarin_PK.RData")
#' # Access individual predictions
#' head(result$GFO$SDTAB)
#'
#' # View parameter estimates
#' print(result$GFO$SUMTAB)
#'
#' # Check objective function value
#' print(result$GFO$OFV)
#'}
#' @importFrom readr read_csv read_tsv cols parse_number
#' @importFrom stringr str_c
#' @import tibble
#' @import dplyr
#' @export


devtools::load_all()
devtools::document()

project_name <- "1cmt-RE-Vd-CL-prop-FEMALE-on-Vd-CRCL-on-CL"
folder_path <- "./scripts/nlme/2.1-sg-converter/monolix-2023/fenoprofen-pk/Monolix/"
folder_path_IFn <- "./scripts/nlme/2.1-sg-converter/monolix-2023/IFN_full/"
project_name_IFN <- "ifn_full_Vmax_bcell_2"
result <- sg_converter(folder_path = folder_path, proj_name = project_name)

result$GFO$SUMTAB
result$GFO$COTAB
result$GFO
result$GCO$modelText
result$GCO$ruv$INIT  # uppercase INIT (not init); expected 0.3 for proportional(b)
class(result$GCO$ruv$DVID)

result_ifn <- sg_converter(folder_path = folder_path_IFn, proj_name = project_name_IFN, save_file = F)

result_ifn$GFO$COTAB
result_ifn$GFO$CATAB
result_ifn$GFO$SUMTAB
result_ifn$GCO
result_ifn$GFO$PROJNAME

read_json("./scripts/nlme/2.1-sg-converter/derived-data/aaa_test.json")

#### 2 outputs ####


folder_path_MBMA <- "./scripts/nlme/2.1-sg-converter/monolix-2023/2output/"
project_name_MBMA <- "ifn_full_Vmax_bcell_outputs2"
result_MBMA <- sg_converter(
  folder_path = folder_path_MBMA,
  proj_name = project_name_MBMA,
  save_file = FALSE
)
result_MBMA$GFO$SUMTAB
result_MBMA$GFO$SDTAB
result_MBMA$GFO$SDTAB$DVNAME
result_MBMA$GCO$ruv[[1]]
result_MBMA$GCO$modelText
result_MBMA$GCO$ruv[[1]]$INIT
class(result_MBMA$GCO$ruv[[1]]$DVID)


sg_gof_tp(result_MBMA$GFO, DVID = 1)
sg_gof_tp?
help(sg_gof_tp)


sg_gof_res_dist(result_MBMA$GFO, DVID = 0)
sg_gof_obpr(result_MBMA$GFO, DVID = 2)


#### initial value for error ####

result <- sg_converter(folder_path = "./scripts/nlme/2.1-sg-converter/derived-data/", proj_name = "wrfrn_pk_base_model_02")

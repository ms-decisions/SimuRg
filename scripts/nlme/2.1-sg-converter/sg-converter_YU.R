## Author: Ugolkov Yaroslav, Victoria Kulesh
## First created: 2025-10-24
## Description: Convert Monolix project results to SimuRg objects
## Keywords: SimuRg, monolix, converter


devtools::load_all()
# devtools::document()

project_name <- "1cmt-RE-Vd-CL-prop-FEMALE-on-Vd-CRCL-on-CL"
folder_path <- "./scripts/nlme/2.1-sg-converter/monolix-2023/fenoprofen-pk/Monolix/"
folder_path_IFn <- "./scripts/nlme/2.1-sg-converter/monolix-2023/IFN_full/"
project_name_IFN <- "ifn_full_Vmax_bcell_2"
result <- sg_converter(folder_path = folder_path, proj_name = project_name)

result$GFO$SUMTAB
result$GFO$COTAB
result$GFO$COVMAT
result$GCO$modelText
result$GCO$ruv$INIT  # uppercase INIT (not init); expected 0.3 for proportional(b)
class(result$GCO$ruv$DVID)

result_ifn <- sg_converter(folder_path = folder_path_IFn, proj_name = project_name_IFN, save_file = T)

result_ifn$GFO$COTAB
result_ifn$GFO$COVMAT
result_ifn$GFO$SUMTAB
result_ifn$GCO
result_ifn$GFO$PROJNAME

read_json("./scripts/nlme/2.1-sg-converter/monolix-2023/IFN_full/ifn_full_Vmax_bcell_2_GFO.json")

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

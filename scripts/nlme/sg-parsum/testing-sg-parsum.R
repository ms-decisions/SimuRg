devtools::load_all()

project_name <- "1cmt-RE-ka-Vd-CL-prop"
folder_path <- "./scripts/nlme/sg-parsum/Monolix/"

tf <- sg_parsum(fpath_i = str_c(folder_path, project_name, "_GFO.RData"), addOFV = F)
tf

tt<- sg_parsum(fpath_i = str_c(folder_path, project_name, "_GFO.RData"), addOFV = T)
tt

## Author: Victoria Kulesh
## First created: 12-08-2025
## Description: arguments for sg-modbuild functions
## Keywords: sg-modbuild


#####--------------- Load functions and libraries ---------------#####
library(tidyverse)


#####--------------- Input structure ---------------#####
folder_path <- str_c(getwd(), "functions/nlme/1.4-sg-modbuild/")

### list of paths to structural models
mod_lst <- list(str_c(folder_path, "models/model_1c.txt"),
                str_c(folder_path, "models/model_2c.txt"))

### path to the dataset
data <- str_c(folder_path, "data/dspk-warf.csv")

### predefined list of names in the dataset
headers <- list(list(name = "ID", use = "identifier", type = NULL),
                list(name = "TIME", use = "time", type = NULL),
                list(name = "DV", use = "observation", type = "continuous"),
                list(name = "DVID", use = "observationtype", type = NULL),
                list(name = "ADM", use = "administration", type = NULL),
                list(name = "AMT", use = "amount", type = NULL),
                list(name = "EVID", use = "eventidentifier", type = NULL),
                list(name = "MDV", use = "missingdependentvariable", type = NULL),
                list(name = "AGE", use = "covariate", type = "continuous"),
                list(name = "AGE_centered", use = "covariate", type = "continuous"),
                list(name = "SEX", use = "covariate", type = "categorical"),
                list(name = "WEIGHT", use = "covariate", type = "continuous"),
                list(name = "BMI", use = "covariate", type = "continuous"))

### list for residual values options

# option 1: 1 list element (for all structural models)
ruv_lst_1 <- list(list(YNAME = "y1", 
                 DVID = 1, 
                 TRANS = "normal", 
                 PRED = "Cc", 
                 ERR = list("constant", "proportional", "combined1"), # options to test (can be length = 1, i.e., "constant" or "combined1") 
                 INIT = list(1, 1, c(1, 1)), # options to test (can be length = 1, i.e., 1 or c(1, 1))
                 EST = list(T, T, c(T, T)), # options to test (can be length = 1, i.e., T or c(T, T))
                 BLQM = NULL))
# if you have DVID > 1 - add new list wit the same structure

# option 2: list of n list elements (n - number of structural models)
ruv_lst_2 <- list(
  # structural model 1
         list(
           list(YNAME = "y1", 
                DVID = 1, 
                TRANS = "normal", 
                PRED = "Cc", 
                ERR = list("constant", "proportional", "combined1"), # options to test (can be length = 1, i.e., "constant" or "combined1") 
                INIT = list(1, 1, c(1, 1)), # options to test (can be length = 1, i.e., 1 or c(1, 1))
                EST = list(T, T, c(T, T))) # options to test (can be length = 1, i.e., T or c(T, T))BLQM = NULL))
              ),
  # structural model 2
         list(
           list(YNAME = "y1", 
                DVID = 1, 
                TRANS = "normal", 
                PRED = "Cc", 
                ERR = list("constant", "proportional"), # options to test (can be length = 1, i.e., "constant" or "combined1") 
                INIT = list(1, 1), # options to test (can be length = 1, i.e., 1 or c(1, 1))
                EST = list(T, T)) # options to test (can be length = 1, i.e., T or c(T, T))BLQM = NULL))
         )
        )


### list for theta values options
#! INIT, TRANS, LB, UB could be length > 1


# option 1: 1 list element (for all structural models - all parameters from each model structure in 1 tibble)
theta_lst_1 <- list(tribble(~NAME, ~TRANS, ~INIT, ~LB, ~UB, ~EST,
                      "Cl", "logNormal", 0.2, NA, NA, T,
                      "Vd", "logNormal", 20, NA, NA, T, 
                      "ka", "logNormal", c(0.2, 0.1, 0.5), NA, NA, T,
                      "Vp", "logNormal", 10, NA, NA, T,  
                      "Q", "logNormal", 5, NA, NA, T ))


# option 2:  list of n list elements (n - number of structural models)
theta_lst_2 <- list(
  # structural model 1
  tribble(~NAME, ~TRANS, ~INIT, ~LB, ~UB, ~EST,
                      "Cl", "logNormal", 0.2, NA, NA, T,
                      "Vd", "logNormal", 20, NA, NA, T, 
                      "ka", "logNormal", 0.2, NA, NA, T),
  # structural model 2
              tribble(~NAME, ~TRANS, ~INIT, ~LB, ~UB, ~EST, 
                      "Cl", "logNormal", 0.2, NA, NA, T,
                      "Vd", c("Normal", "logNormal"), c(10, 20, 30), NA, NA, T, 
                      "ka", "logNormal", c(0.2, 0.1, 0.5), NA, NA, T, # INIT could be length > 1
                      "Vp", c("Normal", "logNormal"), 10, NA, NA, T, 
                      "Q", "logNormal", 5, NA, NA, T))


# option 3:  option for multistart 
#if VAR = True - initial value for the parameter is obtained from the unifrom distribution with the ranges:
# option 3.1:  +- % in LEVEL from initial values set in INIT
# option 3.2:  ranges in INIT 
theta_lst_3 <- list(tribble(~NAME, ~TRANS, ~INIT, ~LB, ~UB, ~EST, ~VAR, ~LEVEL,
                      "Cl", "logNormal", 0.2, NA, NA, T, T, 25, # option 3.1
                      "Vd", "logNormal", 20, NA, NA, T,  NA, NA,
                      "ka", "logNormal", 0.2, NA, NA, T,  T, NA, # option 3.2
                      "Vp", "logNormal", 10, NA, NA, T,  NA, NA,
                      "Q", "logNormal", 5, NA, NA, T,  NA, NA))

            

### list for random effects options
# init - matrix with initial values for omegas (diagonal elements) of correlations (non-diagonal elements)
# est - matrix with identificator, whether the RE/Corr parameter need to be estimated (T - estimated, F - fixed to value in init matrix, NA - not set)
# block - matrix with identificator, whether the RE/Corr parameter need to be blocked (reserved for all options) in model building process (T - blocked RE/Corr parameter (reserved in all scenarious), F - engage in scenarios crossing, NA - not set)
# block matrix is not mandatory (if "block" is omitted  - all RE/Corr parameters are engaged in scenarious crossing)


# option 1: 1 list element (for all structural models - all parameters from each model structure in 1 matrix)
re_lst_1 <- list(
  list(init = tribble(~Cl, ~Vd, ~ka, ~Vp, ~Q,
                      1, 0, 0, 0, 0,
                      0, 1, 0, 0, 0,
                      0, 0, 1, 0, 0,
                      0, 0, 0, 1, 0,
                      0, 0, 0, 0, 1) %>% as.matrix(),
       est = tribble(~Cl, ~Vd, ~ka, ~Vp, ~Q,
                     T, NA, NA, NA, NA,
                     NA, T, NA, NA, NA,
                     NA, NA, T, NA, NA,
                     NA, NA, NA, T, NA,
                     NA, NA, NA, NA, T) %>% as.matrix(),
       block =tribble(~Cl, ~Vd, ~ka, ~Vp, ~Q,
                      F, NA, NA, NA, NA,
                      NA, F, NA, NA, NA,
                      NA, NA, T, NA, NA,
                      NA, NA, NA, F, NA,
                      NA, NA, NA, NA, F) %>% as.matrix())
)



# option 2:  list of n list elements (n - number of structural models)
re_lst_2 <- list(
  # structural model 1 
  # this option corresponds to 4 scenarious:
  # 1 - omega_Vd (blocked parameter)
  # 2 - omega_Cl and omega_Vd (blocked parameter)
  # 3 - omega_ka and omega_Vd (blocked parameter)
  # 4 - omega_Cl, omega_ka and omega_Vd (blocked parameter)

  list(init = tribble(~Cl, ~Vd, ~ka,
                      1, 0, 0, 
                      0, 1, 0,
                      0, 0, 1) %>% as.matrix(),
       est = tribble(~Cl, ~Vd, ~ka,
                     T, NA, NA, 
                     NA, T, NA, 
                     NA, NA, T) %>% as.matrix(),
       block =tribble(~Cl, ~Vd, ~ka,
                      F, NA, NA, 
                      NA, T, NA, 
                      NA, NA, F) %>% as.matrix()),
  
  # structural model 2
  # this option corresponds to 4 scenarious:
  # 1 - no RE
  # 2-7 - 1RE
  # 8-18 - 2RE
  # 19-29 - 3RE
  # 30-34 - 4RE
  # 35 - 5RE (all)
  
  list(init = tribble(~Cl, ~Vd, ~ka, ~Vp, ~Q,
                      1, 0, 0, 0, 0,
                      0, 1, 0, 0, 0,
                      0, 0, 1, 0, 0,
                      0, 0, 0, 1, 0,
                      0, 0, 0, 0, 1) %>% as.matrix(),
       est = tribble(~Cl, ~Vd, ~ka, ~Vp, ~Q,
                     T, NA, NA, NA, NA,
                     NA, T, NA, NA, NA,
                     NA, NA, T, NA, NA,
                     NA, NA, NA, T, NA,
                     NA, NA, NA, NA, T) %>% as.matrix(),
       block =tribble(~Cl, ~Vd, ~ka, ~Vp, ~Q,
                      F, NA, NA, NA, NA,
                      NA, F, NA, NA, NA,
                      NA, NA, F, NA, NA,
                      NA, NA, NA, F, NA,
                      NA, NA, NA, NA, F) %>% as.matrix())
)

# occ_lst - the came logic as for re_lst
# covs_lst - need to compile




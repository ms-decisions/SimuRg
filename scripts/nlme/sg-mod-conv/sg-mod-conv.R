# Specify structural model path. For fit, should be in Monolix syntax
model <-  "scripts/nlme/sg-mod-conv/model_PK_1c.txt"
# Specify data path. Should be ADPPK-like format
data  <- system.file("extdata", "datasets", "dspk-warf.csv", package = "SimuRg")

# Specify headers. Should be a list of lists with the following elements:
# name - string, name of the column
# use - string, use of the column from Monolix documentation
# type - string, type of the column. for use = "covariate", should be "continuous" or "categorical", for other uses should be NULL

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

# Dataset with the parameters properties. Should be a tibble with the following columns:
# NAME - string, name of the parameter
# TRANS - string, distribution of the parameter. Should be one of the following: "normal", "logNormal", "logitNormal"
# INIT - numeric, initial value of the parameter or its fixed value
# LB - numeric, lower bound of the parameter for logit transformation
# UB - numeric, upper bound of the parameter for logit transformation
# EST - logical, estimation status of the parameter. For estimation, should be TRUE, for fixed, should be FALSE

 theta <- tribble(~NAME, ~TRANS, ~INIT, ~LB, ~UB, ~EST,
                  "Cl", "logNormal", 0.2, NA, NA, TRUE,
                  "V", "logNormal", 20, NA, NA, TRUE,
                  "ka", "logNormal", 0.2, NA, NA, TRUE
 )
 # Examples of random effect model specification for fit
 # Single observation (legacy format):
 # YNAME - string, name of the observation, ususally y1, y2, ...
 # DVID - numeric, observation type identifier corresponding to DVID column values
 # TRANS - string, residual error distribution. Can be: "normal", "logNormal", "logitNormal"
 # PRED - string, prediction variable name from the model
 # ERR - string, error model type. Options include: "constant" for additive error, "proportional" for proportional error, "combined1" for combined additive and proportional error
 # INIT - numeric vector, initial values for error parameters (length depends on error model)
 # EST - logical vector, whether to estimate each error parameter (same length as INIT)
 # BLQM - below limit of quantification method (can be NULL)
 ruv <- list(YNAME = "y1", DVID = 1, TRANS = "normal", PRED = "Cc",
             ERR = "combined1", INIT = c(1, 1), EST = c(TRUE, TRUE), BLQM = NULL)

 # Example of random effects (RE) specification.
 #
 # init matrix: provides the initial values for the random effects.
 # est matrix: controls how each random effect is handled:
 #   TRUE  - the random effect will be estimated,
 #   FALSE - the random effect will be fixed at its initial value,
 #   NA    - no random effect will be applied.
 #
 # To fit a model without random effects, set all entries in the
 # est matrix to NA.
 #
 # The same logic applies to the between-occasion variability
 # matrix (occ).

 re <- list(init = tribble(~Cl, ~V, ~ka,
                           1, 0, 0,
                           0, 1, 0,
                           0, 0, 1) %>% as.matrix(),
            est = tribble(~Cl, ~V, ~ka,
                          TRUE, NA, NA,
                          NA, TRUE, NA,
                          NA, NA, TRUE) %>% as.matrix())
# Example of between-occasion variability (BOV) specification. The structure is the same as for RE, but for BOV
 occ <- list(init = tribble(~Cl, ~V, ~ka,
                            0, 0, 0,
                            0, 0, 0,
                            0, 0, 0) %>% as.matrix(),
             est = tribble(~Cl, ~V, ~ka,
                           NA, NA, NA,
                           NA, NA, NA,
                           NA, NA, NA) %>% as.matrix())

# Example of covariate specification. Should be a list of lists with the following elements:
# PAR - string, name of the parameter to which the covariate is applied
# COVNAME - string, name of the covariate
# FUNC - string, function to apply to the covariate. Should be "linear" for continuous covariates, "categorical" for categorical covariates
# TRANS - string, transformation of the covariate. Should be "median" for continuous covariates, "reference" for categorical covariates
# INIT - numeric, initial value of the covariate
# EST - logical, estimation status of the covariate. For estimation, should be TRUE, for fixed, should be FALSE
 covs <- list(list(PAR = "V", COVNAME = "AGE", FUNC = "linear",
                   TRANS = "median", INIT = 1, EST = TRUE),
              list(PAR = "ka", COVNAME = "SEX", REF = 0, INIT = 1, EST = TRUE))
 output_path <- "scripts/nlme/sg-mod-conv/"
 gco <-list(headers = headers,
            data = data,
            model = model,
            task_opt = "",
            covs = covs,
            project_name = "test-proj",
            theta = theta,
            ruv = ruv,
            re = re,
            occ = occ,
            modelText = "")
 save(gco, file = "scripts/nlme/sg-mod-conv/gco_example.Rdata")

 sg_gmo_build(gco = gco, output_path = "scripts/nlme/sg-mod-conv/gmo6.R")


gco <-list(headers = headers,
           data = data,
           model = model,
           task_opt = "",
           covs = covs,
           project_name = "test-proj",
           theta = theta,
           ruv = ruv,
           re = list(init = tribble(~Cl, ~V, ~ka,
                                    1, 0, 0,
                                    0, 0, 0,
                                    0, 0, 1) %>% as.matrix(),
                     est = tribble(~Cl, ~V, ~ka,
                                   TRUE, NA, NA,
                                   NA, NA, NA,
                                   NA, NA, TRUE) %>% as.matrix()),
           occ = occ,
           modelText = "")

sg_gmo_build(gco = gco, output_path = "scripts/nlme/sg-mod-conv/gmo5.R")
gco <-list(headers = headers,
           data = data,
           model = model,
           task_opt = "",
           covs = covs,
           project_name = "test-proj",
           theta = tribble(~NAME, ~TRANS, ~INIT, ~LB, ~UB, ~EST,
                           "Cl", "Normal", 0.2, NA, NA, TRUE,
                           "V", "logNormal", 10, NA, NA, TRUE,
                           "ka", "Normal", 0.2, NA, NA, TRUE
           ),
           ruv = ruv,
           re = re,
           occ = occ,
           modelText = "")

sg_gmo_build(gco = gco, output_path = "scripts/nlme/sg-mod-conv/gmo7.R")

gco <-list(headers = headers,
           data = data,
           model = model,
           task_opt = "",
           covs = covs,
           project_name = "test-proj",
           theta = theta,
           ruv = list(
             list(YNAME = "y1", DVID = 1, TRANS = "normal", PRED = "Cc",
                  ERR = "combined1", INIT = c(1, 1), EST = c(TRUE, TRUE), BLQM = NULL),
             list(YNAME = "y2", DVID = 2, TRANS = "normal", PRED = "Cc_nM",
                  ERR = "proportional", INIT = c(0.1), EST = c(TRUE), BLQM = NULL)
           ),
           re = re,
           occ = occ,
           modelText = "")

sg_mod_conv(gco = gco, output_path = "scripts/nlme/sg-mod-conv/gmo3.R")



 # Multiple observations (recommended format):

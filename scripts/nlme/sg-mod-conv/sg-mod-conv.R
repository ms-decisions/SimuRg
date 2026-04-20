model <-  "scripts/nlme/sg-mod-conv/model_PK_1c.txt"
data  <- system.file("extdata", "datasets", "dspk-warf.csv", package = "SimuRg")

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

 theta <- tribble(~NAME, ~TRANS, ~INIT, ~LB, ~UB, ~EST,
                  "Cl", "logNormal", 0.2, NA, NA, TRUE,
                  "V", "logNormal", 20, NA, NA, TRUE,
                  "ka", "logNormal", 0.2, NA, NA, TRUE
 )
 # Examples of random effect model specification
 # Single observation (legacy format):
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

 occ <- list(init = tribble(~Cl, ~V, ~ka,
                            0, 0, 0,
                            0, 0, 0,
                            0, 0, 0) %>% as.matrix(),
             est = tribble(~Cl, ~V, ~ka,
                           NA, NA, NA,
                           NA, NA, NA,
                           NA, NA, NA) %>% as.matrix())
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

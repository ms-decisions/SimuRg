test_that("sg_multistart creates expected number of files", {
  # str_c("V:/Collaborative_working/SimuRg_as_R_lib/SimuRg/scripts/nlme/1.4-sg-modbuild/")
  folder_path <- system.file("extdata", package = "SimuRg")

  mod <- paste(folder_path, "/models/model_PK_1c.txt", sep = "/")
  ### path to the dataset
  data <- paste(folder_path, "datasets", "dspk-warf.csv", sep = "/")

  re <- list(init = tribble(~Cl, ~Vd, ~ka, ~Vp, ~Q,
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
                   "Cl", "logNormal", 0.2, NA, NA, T,
                   "Vd", "logNormal", 20, NA, NA, T,
                   "ka", "logNormal", 0.2, NA, NA, T)

  theta_intervals <- list(
    Cl = c(0.25*theta$INIT[theta$NAME == "Cl"], 4*theta$INIT[theta$NAME == "Cl"]),
    #Vd = c(0.5*theta$INIT[theta$NAME == "Vd"], 2*theta$INIT[theta$NAME == "Vd"]),
    ka   = c(1.25*theta$INIT[theta$NAME == "ka"], 1.5*theta$INIT[theta$NAME == "ka"])
  )

  path <- tempdir()# tempdir() #system.file("extdata", package = "SimuRg")

  ruv <- list(YNAME = "y1",
              DVID = 1,
              TRANS = "normal",
              PRED = "Cc",
              ERR = "proportional",
              INIT = 1,
              EST = T)

  n_starts <- 20

  sg_multistart(
    mod = mod,
    data = data,
    headers = headers,
    ruv = ruv,
    theta = theta,
    re = re,
    occ = re,
    n_starts = n_starts,
    theta_intervals = theta_intervals,
    covs_lst = NULL,
    path = paste0(path, "\\"),
    project_name = "multistart_test_project",
    seed = 126
  )

  # Проверяем что создался CSV файл с информацией о сценариях
  expect_true(file.exists(paste0(path, "/multistart_info.csv")))

  # Проверяем наличие первого файла
  mlxtran_files <- list.files(path, pattern = "\\.mlxtran$")
  expect_true("multistart_test_project_1.mlxtran" %in% mlxtran_files)

  # Проверяем что CSV файл не пустой
  csv_content <- read_csv(paste0(path, "/multistart_info.csv"))
  # expect_gt(nrow(csv_content), 0)
  expect_equal(nrow(csv_content), n_starts)
  # Clear tempdir
  clr_files <- list.files(path, full.names = TRUE)
  unlink(clr_files, recursive = TRUE, force = TRUE)


})

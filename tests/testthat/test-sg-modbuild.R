## Author: Ugolkov Yaroslav
## First created: 2025-31-10
## Description: formal testing of sg-modbuild function
## Keywords: SimuRg, sg-modbuild, model building

test_that("sg_modbuild creates expected number of files", {

  folder_path <- system.file("extdata", package = "SimuRg")

  ### list of paths to structural models
  mod_lst <- list(paste(folder_path, "/models/model_PK_1c.txt", sep = "/"),
                  paste(folder_path, "/models/model_PK_2c.txt", sep = "/"))

  ### path to the dataset
  data <- paste(folder_path, "datasets", "dspk-warf.csv", sep = "/")
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

  path <- tempdir() #system.file("extdata", package = "SimuRg")
  sg_modbuild(
    mod_lst = mod_lst[1],
    data = data,
    headers = headers,
    ruv_lst = ruv_lst_2,
    theta_lst = theta_lst_2,
    re_lst = re_lst_1,
    occ_lst = re_lst_1,
    covs_lst = NULL,
    path = paste0(path, "\\"),
    project_name = "tests_test_project"
  )

  # Проверяем что создался CSV файл с информацией о сценариях
  expect_true(file.exists(paste0(path, "/scenarios_info.csv")))

  # Проверяем что создались mlxtran файлы
  # В данном случае ожидаем 2 файла (2 варианта RUV: constant + proportional)
  mlxtran_files <- list.files(path, pattern = "\\.mlxtran$")
  expect_length(mlxtran_files, 768)

  # Проверяем имена файлов
  expect_true(all(c("tests_test_project_1.mlxtran", "tests_test_project_768.mlxtran") %in% mlxtran_files))

  # Проверяем что CSV файл не пустой
  csv_content <- read_csv(paste0(path, "/scenarios_info.csv"))
  # expect_gt(nrow(csv_content), 0)
  expect_equal(nrow(csv_content), 768) # 2 сценария
  clr_files <- list.files(path, full.names = TRUE)
  unlink(clr_files, recursive = TRUE, force = TRUE)


})

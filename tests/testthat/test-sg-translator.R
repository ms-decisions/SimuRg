## Author: Anna Mikhailova
## First created: 2026-05-04
## Description: Basic testing of sg-converter function
## Keywords: SimuRg, sg-converter

test_that("sg-translator from monolix to rxode2 works", {
  path_to_save <- tempdir()
  monolix_model_hv_iv <- system.file("extdata", "models", "monolix",
                                     "model_PK_hv_iv.txt", package = "SimuRg")
  monolix_model_1c <- system.file("extdata", "models", "monolix",
                                  "model_PK_1c.txt", package = "SimuRg")



  # Convert MLXTRAN model to rxode2
  sg_translator(monolix_model_hv_iv, to = "rxode",
                output_path = normalizePath(file.path(path_to_save,
                                                      "model_PK_hv_iv_rx.txt"),
                                             mustWork = FALSE))
  sg_translator(monolix_model_1c, to = "rxode",
                output_path = normalizePath(file.path(path_to_save,
                                                      "model_PK_1c_rx.txt"),
                                            mustWork = FALSE))
  expect_true(file.exists(file.path(path_to_save, "model_PK_hv_iv_rx.txt")))
  expect_true(file.exists(file.path(path_to_save, "model_PK_1c_rx.txt")))
  clr_files <- list.files(path_to_save, full.names = TRUE)
  unlink(clr_files, recursive = TRUE, force = TRUE)
})

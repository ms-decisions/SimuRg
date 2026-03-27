## Author: Yaroslav Ugolkov?
## First created: 2025-10-17
## Description: Basic testing of sg-converter function
## Keywords: SimuRg, sg-converter

test_that("sg-converter works and contains all elements", {
  # skip_if(T)
  test_folder <- system.file("extdata", "Monolix_objects", package = "SimuRg")
  if (substr(test_folder, nchar(test_folder), nchar(test_folder)) != "/")
    test_folder <- str_c(test_folder, "/")
  pro_name <- "proj-solo"

  # Simple function call
  result <- sg_converter(folder_path = test_folder, proj_name = pro_name)
  expect_type(result, "list")

  required_elements <- c("SDTAB", "SUMTAB", "SIGMAMAT", "OMEGAMAT", "OCCMAT",
                         "EVTAB", "PATAB", "COTAB", "CATAB", "REGTAB",
                         "OFV", "COVMAT", "CORRMAT", "OPTIONS", "PROJNAME")
  expect_true(all(required_elements %in% names(result$GFO)))

  expect_gt(nrow(result$GFO$SDTAB), 0)
  expect_gt(nrow(result$GFO$SUMTAB), 0)
  expect_gt(nrow(result$GFO$PATAB), 0)
})

test_that("sg-converter fails gracefully", {
  expect_error(
    suppressWarnings(sg_converter("invalid_path", "invalid_project")),
    "Project file does not exist. Check file existance or try to use absolute path"
  )
})


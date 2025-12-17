## Author: [Ваше имя]
## First created: 2025-10-17
## Description: Formal testing of sg-converter function
## Keywords: SimuRg, sg-converter, monolix

## Author: [Ваше имя]
## First created: 2025-10-17
## Description: Basic testing of sg-converter function
## Keywords: SimuRg, sg-converter
test_that("sg-converter works and contains all elements", {
  test_folder <- system.file("extdata", "Monolix_objects", package = "SimuRg")
  if (substr(test_folder, nchar(test_folder), nchar(test_folder)) != "/") test_folder <- str_c(test_folder, "/")
  pro_name <- "proj-r-solo"

  # Просто вызываем функцию и проверяем результат
  result <- sg_converter(folder_path = test_folder, proj_name = pro_name)

  # Проверяем что функция отрабатывает и возвращает список
  expect_type(result, "list")

  # Проверяем что все основные элементы присутствуют
  required_elements <- c("SDTAB", "SUMTAB", "SIGMAMAT", "OMEGAMAT", "OCCMAT",
                         "EVTAB", "PATAB", "COTAB", "CATAB", "REGTAB",
                         "OFV", "COVMAT", "CORRMAT", "OPTIONS", "PROJNAME")
  expect_true(all(required_elements %in% names(result)))

  # Проверяем что ключевые таблицы не пустые
  expect_gt(nrow(result$SDTAB), 0)
  expect_gt(nrow(result$SUMTAB), 0)
  expect_gt(nrow(result$PATAB), 0)
})

test_that("sg-converter fails gracefully", {
  # Ожидаем error для несуществующих файлов (подавляем warning)
  expect_error(
    suppressWarnings(sg_converter("invalid_path", "invalid_project")),
    "cannot open"
  )
})

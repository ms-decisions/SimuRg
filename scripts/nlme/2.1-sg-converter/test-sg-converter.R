## Author: Yaroslav Ugolkov
## First created: 2025-10-17
## Description: Unit tests for sg-converter requirements/risks
## Keywords: SimuRg, sg-converter

source("./scripts/nlme/2.1-sg-converter/sg-converter_AMel.R")

fixture_copy <- function() {
  src_extdata <- system.file("extdata", package = "SimuRg")
  dst_root <- tempfile("sg-converter-fixture-")
  dir.create(dst_root, recursive = TRUE, showWarnings = FALSE)
  copied <- file.copy(src_extdata, dst_root, recursive = TRUE)
  expect_true(copied)

  folder <- file.path(dst_root, "extdata", "Monolix_objects")
  list(
    root = file.path(dst_root, "extdata"),
    folder = folder,
    folder_path = paste0(normalizePath(folder, winslash = "/"), "/"),
    proj = "proj-solo"
  )
}

project_files <- function(fx) {
  list(
    mlx = file.path(fx$folder, paste0(fx$proj, ".mlxtran")),
    csv = file.path(fx$root, "datasets", "dspk-warf.csv"),
    pop = file.path(fx$folder, fx$proj, "populationParameters.txt"),
    fi_dir = file.path(fx$folder, fx$proj, "FisherInformation")
  )
}

replace_first_matching_line <- function(lines, pattern, replacement) {
  idx <- which(grepl(pattern, lines, perl = TRUE))[1]
  if (is.na(idx)) stop("Pattern not found: ", pattern)
  lines[idx] <- replacement
  lines
}

insert_line_after <- function(lines, pattern, line_to_insert) {
  idx <- which(grepl(pattern, lines, perl = TRUE))[1]
  if (is.na(idx)) stop("Anchor pattern not found: ", pattern)
  append(lines, values = line_to_insert, after = idx)
}

update_mlx <- function(fx, edit_fun) {
  p <- project_files(fx)$mlx
  lines <- readLines(p, warn = FALSE)
  writeLines(edit_fun(lines), p)
}

run_converter <- function(fx) {
  suppressMessages(
    sg_converter(folder_path = fx$folder_path, proj_name = fx$proj)
  )
}

test_that("1) missing mlxtran project file errors clearly", {
  expect_error(
    suppressWarnings(sg_converter("invalid_path", "invalid_project")),
    "Project file does not exist. Check file existance or try to use absolute path"
  )
})

test_that("2) DATAFILE path parsing supports absolute and Monolix 2024 syntax", {
  # Absolute path case
  fx_abs <- fixture_copy()
  files_abs <- project_files(fx_abs)
  abs_data_path <- normalizePath(files_abs$csv, winslash = "/")
  update_mlx(fx_abs, function(lines) {
    replace_first_matching_line(lines, "^file\\s*=", sprintf("file='%s'", abs_data_path))
  })
  res_abs <- run_converter(fx_abs)
  expect_type(res_abs, "list")
  expect_gt(nrow(res_abs$GFO$SDTAB), 0)

  # Monolix 2024 file={path='...'} case
  fx_2024 <- fixture_copy()
  update_mlx(fx_2024, function(lines) {
    replace_first_matching_line(lines, "^file\\s*=", "file={path='../datasets/dspk-warf.csv'}")
  })
  res_2024 <- run_converter(fx_2024)
  expect_type(res_2024, "list")
  expect_gt(nrow(res_2024$GFO$SDTAB), 0)
})

test_that("3) both CSV and TSV datasets are parsed", {
  fx <- fixture_copy()
  files <- project_files(fx)

  dt <- read.csv(files$csv, check.names = FALSE)
  tsv_path <- file.path(fx$root, "datasets", "dspk-warf.tsv")
  write.table(dt, tsv_path, sep = "\t", row.names = FALSE, quote = FALSE)

  update_mlx(fx, function(lines) {
    replace_first_matching_line(lines, "^file\\s*=", "file='../datasets/dspk-warf.tsv'")
  })

  result <- run_converter(fx)
  expect_gt(nrow(result$GFO$SDTAB), 0)
  expect_true(all(c("ID", "TIME", "DV", "DVID") %in% names(result$GFO$SDTAB)))
})

test_that("4) required mapping columns are renamed without losing target-name conflicts", {
  fx <- fixture_copy()
  files <- project_files(fx)
  dt <- read.csv(files$csv, check.names = FALSE)

  # Create remapped source columns and keep conflicting target names with bogus values.
  dt$SUBJ <- dt$ID
  dt$TM <- dt$TIME
  dt$OBS <- dt$DV
  dt$TYPE <- dt$DVID
  dt$DOSE <- dt$AMT
  dt$EVT <- dt$EVID
  dt$ID <- -999
  dt$TIME <- -999
  dt$DV <- -999
  dt$DVID <- -999
  dt$AMT <- -999
  dt$EVID <- -999
  write.csv(dt, files$csv, row.names = FALSE)

  update_mlx(fx, function(lines) {
    lines <- replace_first_matching_line(lines, "^ID\\s*=\\s*\\{use=identifier\\}", "SUBJ = {use=identifier}")
    lines <- replace_first_matching_line(lines, "^TIME\\s*=\\s*\\{use=time\\}", "TM = {use=time}")
    lines <- replace_first_matching_line(lines, "^DV\\s*=\\s*\\{use=observation, yname='1', type=continuous\\}", "OBS = {use=observation, yname='1', type=continuous}")
    lines <- replace_first_matching_line(lines, "^DVID\\s*=\\s*\\{use=observationtype\\}", "TYPE = {use=observationtype}")
    lines <- replace_first_matching_line(lines, "^AMT\\s*=\\s*\\{use=amount\\}", "DOSE = {use=amount}")
    lines <- replace_first_matching_line(lines, "^EVID\\s*=\\s*\\{use=eventidentifier\\}", "EVT = {use=eventidentifier}")
    lines
  })

  result <- run_converter(fx)
  # SDTAB carries observation-level fields; dose/event fields are checked in EVTAB.
  expect_true(all(c("ID", "TIME", "DV", "DVID") %in% names(result$GFO$SDTAB)))
  expect_true(all(c("ID", "TIME", "AMT", "EVID") %in% names(result$GFO$EVTAB)))
  expect_true(all(result$GFO$SDTAB$ID != -999))
  expect_true(all(result$GFO$SDTAB$TIME != -999))
  expect_true(all(result$GFO$SDTAB$DV != -999))
})

test_that("5) duplicate non-covariate `use` mappings pass if identical and fail if different", {
  # Pass case: duplicated time mapping with identical values
  fx_ok <- fixture_copy()
  files_ok <- project_files(fx_ok)
  dt_ok <- read.csv(files_ok$csv, check.names = FALSE)
  dt_ok$TIME_DUP <- dt_ok$TIME
  write.csv(dt_ok, files_ok$csv, row.names = FALSE)

  update_mlx(fx_ok, function(lines) {
    lines <- insert_line_after(lines, "^TIME\\s*=\\s*\\{use=time\\}", "TIME_DUP = {use=time}")
    lines
  })
  expect_no_error(run_converter(fx_ok))

  # Fail case: duplicated time mapping with different values
  fx_bad <- fixture_copy()
  files_bad <- project_files(fx_bad)
  dt_bad <- read.csv(files_bad$csv, check.names = FALSE)
  dt_bad$TIME_DUP <- dt_bad$TIME + 0.5
  write.csv(dt_bad, files_bad$csv, row.names = FALSE)

  update_mlx(fx_bad, function(lines) {
    lines <- insert_line_after(lines, "^TIME\\s*=\\s*\\{use=time\\}", "TIME_DUP = {use=time}")
    lines
  })
  expect_error(
    run_converter(fx_bad),
    "Multiple columns found for use='time'.*different content"
  )
})

test_that("6) FIT endpoint map enforces matching data/model lengths", {
  fx <- fixture_copy()
  update_mlx(fx, function(lines) {
    replace_first_matching_line(lines, "^model\\s*=\\s*\\{y1\\}", "model = {y1, y2}")
  })
  expect_error(
    run_converter(fx),
    "length of data list.*model list.*must be equal"
  )
})

test_that("7) residual error mapping covers proportional, constant, and combined", {
  # proportional (fixture default)
  fx_prop <- fixture_copy()
  res_prop <- run_converter(fx_prop)
  expect_equal(res_prop$GCO$ruv$ERR, "proportional")
  expect_length(res_prop$GCO$ruv$INIT, 1)

  # constant
  fx_const <- fixture_copy()
  update_mlx(fx_const, function(lines) {
    replace_first_matching_line(
      lines,
      "^y1\\s*=\\s*\\{distribution=normal, prediction=Cc, errorModel=proportional\\(Cc_b\\)\\}",
      "y1 = {distribution=normal, prediction=Cc, errorModel=constant(Cc_b)}"
    )
  })
  res_const <- run_converter(fx_const)
  expect_equal(res_const$GCO$ruv$ERR, "constant")
  expect_length(res_const$GCO$ruv$INIT, 1)

  # combined
  fx_comb <- fixture_copy()
  files_comb <- project_files(fx_comb)
  update_mlx(fx_comb, function(lines) {
    lines <- replace_first_matching_line(
      lines,
      "^y1\\s*=\\s*\\{distribution=normal, prediction=Cc, errorModel=proportional\\(Cc_b\\)\\}",
      "y1 = {distribution=normal, prediction=Cc, errorModel=combined(Cc_a, Cc_b)}"
    )
    lines <- insert_line_after(lines, "^Cc_b\\s*=\\s*\\{value=1, method=MLE\\}", "Cc_a = {value=0.5, method=MLE}")
    lines
  })
  pop_dt <- read.csv(files_comb$pop, check.names = FALSE)
  pop_dt <- rbind(
    pop_dt,
    data.frame(
      parameter = "Cc_a",
      value = 0.05,
      CV = NA_real_,
      se_lin = NA_real_,
      rse_lin = NA_real_,
      se_sa = NA_real_,
      rse_sa = NA_real_,
      stringsAsFactors = FALSE
    )
  )
  write.csv(pop_dt, files_comb$pop, row.names = FALSE)

  res_comb <- run_converter(fx_comb)
  expect_equal(res_comb$GCO$ruv$ERR, "combined")
  expect_length(res_comb$GCO$ruv$INIT, 2)
})

test_that("8) SDTAB joins by ID/TIME and fills MDV=0 when MDV source is absent", {
  fx <- fixture_copy()
  files <- project_files(fx)

  dt <- read.csv(files$csv, check.names = FALSE)
  dt$MDV <- NULL
  write.csv(dt, files$csv, row.names = FALSE)

  update_mlx(fx, function(lines) {
    lines[!grepl("^MDV\\s*=\\s*\\{use=missingdependentvariable\\}", lines, perl = TRUE)]
  })

  result <- run_converter(fx)
  expect_true("MDV" %in% names(result$GFO$SDTAB))
  expect_true(all(result$GFO$SDTAB$MDV == 0))
  expect_false(any(is.na(result$GFO$SDTAB$PRED)))
  expect_false(any(is.na(result$GFO$SDTAB$IPRED)))
})

test_that("9) Monte Carlo WRES is reproducible with an explicit seed", {
  fx <- fixture_copy()
  set.seed(2026)
  res_1 <- run_converter(fx)
  set.seed(2026)
  res_2 <- run_converter(fx)

  expect_equal(res_1$GFO$SDTAB$WRES, res_2$GFO$SDTAB$WRES, tolerance = 1e-12)
})

test_that("10) GFO/GCO preserve required shape and optional placeholders", {
  fx <- fixture_copy()
  files <- project_files(fx)
  if (dir.exists(files$fi_dir)) unlink(files$fi_dir, recursive = TRUE)

  result <- run_converter(fx)
  expect_type(result, "list")
  expect_true(all(c("GFO", "GCO") %in% names(result)))

  required_gfo <- c(
    "SDTAB", "SUMTAB", "SIGMAMAT", "OMEGAMAT", "OCCMAT", "EVTAB", "PATAB",
    "COTAB", "CATAB", "REGTAB", "OFV", "COVMAT", "CORRMAT", "OPTIONS", "PROJNAME"
  )
  expect_true(all(required_gfo %in% names(result$GFO)))
  expect_true(all(c("headers", "data", "model", "theta", "ruv", "re", "occ") %in% names(result$GCO)))

  expect_true(is.matrix(result$GFO$COVMAT))
  expect_true(is.matrix(result$GFO$CORRMAT))
  expect_true(is.matrix(result$GFO$OCCMAT))
  # Empty placeholders are created as matrix(), which is a 1x1 matrix with NA.
  expect_lte(length(result$GFO$COVMAT), 1)
  expect_lte(length(result$GFO$CORRMAT), 1)
  expect_lte(length(result$GFO$OCCMAT), 1)
  expect_null(result$GFO$OPTIONS)
})


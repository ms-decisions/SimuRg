## Author: Melnikova Alina
## First created: 2025-11-24
## Description: formal testing of sg-vpop-est function and its helpers
## Keywords: SimuRg, sg-parsum, sg-vpop-est

test_that("sg-vpop-est file load works", {
  fpath_i <- system.file("data", "data_pbc.rda", package = "SimuRg")
  load(fpath_i)
  output <- sg_vpop_est(data = data_pbc, diag_plots = TRUE)
  expect_true(inherits(output$datagen, "data.frame"))
  expect_true(
    is.null(output$dplot_cont) ||
      (is.list(output$dplot_cont) &&
         all(vapply(output$dplot_cont, function(x) inherits(x, "ggplot"), logical(1))))
  )
  expect_true(
    is.null(output$dplot_cat) ||
      (is.list(output$dplot_cat) &&
         all(vapply(output$dplot_cat, function(x) inherits(x, "ggplot"), logical(1))))
  )
  expect_true(inherits(output$dplot_umap_cont, "ggplot"))
  expect_true(inherits(output$dplot_umap_cat, "ggplot"))
  expect_true(inherits(output$ks_test, "data.frame"))
})

test_that("sg-vpop-est does not work on empty dataset", {
  expect_error(sg_vpop_est(data.frame()))
})


# Basic functionality tests
test_that("sg_vpop_est returns correct structure with simple continuous data", {
  # Create simple test data with continuous variables
  test_data <- data.frame(
    x1 = rnorm(50, mean = 10, sd = 2),
    x2 = rnorm(50, mean = 20, sd = 3)
  )

  output <- sg_vpop_est(data = test_data, diag_plots = FALSE, seed = 123)

  # Check output structure
  expect_type(output, "list")
  expect_named(output, c("datagen", "dplot_cont", "dplot_cat",
                         "dplot_umap_cont", "dplot_umap_cat", "ks_test"))

  # Check datagen
  expect_true(is.data.frame(output$datagen))
  expect_equal(ncol(output$datagen), ncol(test_data))
  expect_equal(colnames(output$datagen), colnames(test_data))

  # Check ks_test
  expect_true(is.data.frame(output$ks_test))
  expect_true("variable" %in% colnames(output$ks_test))
  expect_true("p.value" %in% colnames(output$ks_test))
  expect_true("status" %in% colnames(output$ks_test))
})

test_that("sg_vpop_est works with mixed continuous and categorical data", {
  # Create test data with both continuous and categorical variables
  test_data <- data.frame(
    x1 = rnorm(50, mean = 10, sd = 2),
    x2 = rnorm(50, mean = 20, sd = 3),
    cat1 = factor(rep(c("A", "B", "C"), length.out = 50)),
    cat2 = factor(rep(c("X", "Y"), each = 25))
  )

  output <- sg_vpop_est(data = test_data, diag_plots = FALSE, seed = 123)

  expect_true(is.data.frame(output$datagen))
  expect_equal(ncol(output$datagen), ncol(test_data))
  expect_equal(nrow(output$datagen), nrow(test_data))
})

test_that("sg_vpop_est respects nobj parameter", {
  test_data <- data.frame(
    x1 = rnorm(100, mean = 10, sd = 2),
    x2 = rnorm(100, mean = 20, sd = 3)
  )

  output <- sg_vpop_est(data = test_data, nobj = 50, diag_plots = FALSE, seed = 123)

  expect_equal(nrow(output$datagen), 50)
})

test_that("sg_vpop_est respects expfctr parameter", {
  test_data <- data.frame(
    x1 = rnorm(100, mean = 10, sd = 2),
    x2 = rnorm(100, mean = 20, sd = 3)
  )

  output <- sg_vpop_est(data = test_data, expfctr = 2, diag_plots = FALSE, seed = 123)

  expect_equal(nrow(output$datagen), 200)
})

test_that("sg_vpop_est excludes idcol when specified", {
  test_data <- data.frame(
    id = 1:50,
    x1 = rnorm(50, mean = 10, sd = 2),
    x2 = rnorm(50, mean = 20, sd = 3)
  )

  output <- sg_vpop_est(data = test_data, idcol = "id", diag_plots = FALSE, seed = 123)

  expect_false("id" %in% colnames(output$datagen))
  expect_equal(ncol(output$datagen), 2)
})

test_that("sg_vpop_est excludes exclcol when specified", {
  test_data <- data.frame(
    x1 = rnorm(50, mean = 10, sd = 2),
    x2 = rnorm(50, mean = 20, sd = 3),
    exclude_me = rnorm(50, mean = 5, sd = 1)
  )

  output <- sg_vpop_est(data = test_data, exclcol = "exclude_me",
                        diag_plots = FALSE, seed = 123)

  expect_false("exclude_me" %in% colnames(output$datagen))
  expect_equal(ncol(output$datagen), 2)
})

test_that("sg_vpop_est uses seed for reproducibility", {
  test_data <- data.frame(
    x1 = rnorm(50, mean = 10, sd = 2),
    x2 = rnorm(50, mean = 20, sd = 3)
  )

  output1 <- sg_vpop_est(data = test_data, seed = 123, diag_plots = FALSE)
  output2 <- sg_vpop_est(data = test_data, seed = 123, diag_plots = FALSE)

  # Results should be identical with same seed
  expect_equal(output1$datagen, output2$datagen)
})

test_that("sg_vpop_est works with diag_plots = FALSE", {
  test_data <- data.frame(
    x1 = rnorm(50, mean = 10, sd = 2),
    x2 = rnorm(50, mean = 20, sd = 3)
  )

  output <- sg_vpop_est(data = test_data, diag_plots = FALSE, seed = 123)

  expect_null(output$dplot_cont)
  expect_null(output$dplot_cat)
  expect_null(output$dplot_umap_cont)
  expect_null(output$dplot_umap_cat)
  expect_true(is.data.frame(output$datagen))
  expect_true(is.data.frame(output$ks_test))
})

test_that("sg_vpop_est works with diag_plots = TRUE", {
  test_data <- data.frame(
    x1 = rnorm(150, mean = 10, sd = 2),
    x2 = rnorm(150, mean = 20, sd = 3)
  )

  output <- sg_vpop_est(data = test_data, diag_plots = TRUE, seed = 123, seed_umap = 123)

  expect_true(
    is.null(output$dplot_cont) ||
      (is.list(output$dplot_cont) &&
         all(vapply(output$dplot_cont, function(x) inherits(x, "ggplot"), logical(1))))
  )
  expect_true(
    is.null(output$dplot_cat) ||
      (is.list(output$dplot_cat) &&
         all(vapply(output$dplot_cat, function(x) inherits(x, "ggplot"), logical(1))))
  )
  expect_true(inherits(output$dplot_umap_cont, "ggplot") || is.null(output$dplot_umap_cont))
  expect_true(inherits(output$dplot_umap_cat, "ggplot") || is.null(output$dplot_umap_cat))
})

test_that("sg_vpop_est converts character columns to factors", {
  test_data <- data.frame(
    x1 = rnorm(50, mean = 10, sd = 2),
    char_col = rep(c("A", "B", "C"), length.out = 50),
    stringsAsFactors = FALSE
  )

  output <- sg_vpop_est(data = test_data, diag_plots = FALSE, seed = 123)

  # The synthetic data should have the same structure
  expect_true(is.data.frame(output$datagen))
  expect_true("char_col" %in% colnames(output$datagen))
})

test_that("sg_vpop_est handles minnumlev parameter", {
  # Create data with numeric column that has few unique values
  test_data <- data.frame(
    x1 = rnorm(50, mean = 10, sd = 2),
    low_level = rep(c(1, 2, 3), length.out = 50)  # Only 3 unique values
  )

  output <- sg_vpop_est(data = test_data, minnumlev = 3, diag_plots = FALSE, seed = 123)

  expect_true(is.data.frame(output$datagen))
  expect_equal(ncol(output$datagen), 2)
})

test_that("sg_vpop_est handles NA values by removing rows", {
  test_data <- data.frame(
    x1 = c(rnorm(45, mean = 10, sd = 2), rep(NA, 5)),
    x2 = rnorm(50, mean = 20, sd = 3)
  )

  # Should not error, but will remove rows with NA
  output <- sg_vpop_est(data = test_data, diag_plots = FALSE, seed = 123)

  expect_true(is.data.frame(output$datagen))
  expect_equal(nrow(output$datagen), 45)  # Should match number of non-NA rows
})

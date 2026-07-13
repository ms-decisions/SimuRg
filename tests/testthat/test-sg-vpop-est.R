## Author: Melnikova Alina
## First created: 2025-11-24
## Description: formal testing of sg-vpop-est function and its helpers
## Keywords: SimuRg, sg-vpop-est



test_that("sg-vpop-est file load works", {

  # fpath_i <- system.file("data", "data_pbc.rda", package = "SimuRg")
  # load(fpath_i)
  output <- sg_vpop_est(data_i = data_pbc, diag_plots = TRUE, id_col = "id",
                        excl_col = "years", seed = 123)

  # Output is now a list of lists
  #expect_type(output, "list")
  #expect_length(output, 1)  # Default generates 1 dataset

  # Check first dataset structure
  result <- output
  expect_true(inherits(result$datagen, "data.frame"))
  expect_true(is.numeric(result$seed) || is.na(result$seed))
  expect_true(is.logical(result$exact_dupl_check))

  # Check plots
  expect_true(
    is.null(result$dplot_cont) ||
      (is.list(result$dplot_cont) &&
         all(vapply(result$dplot_cont, function(x) inherits(x, "ggplot"), logical(1))))
  )
  expect_true(
    is.null(result$dplot_cat) ||
      (is.list(result$dplot_cat) &&
         all(vapply(result$dplot_cat, function(x) inherits(x, "ggplot"), logical(1))))
  )
  expect_true(inherits(result$dplot_umap, "ggplot"))
  expect_true(
    is.null(result$ks_test) || inherits(result$ks_test, "data.frame")
  )
  expect_true(
    is.null(result$dplot_corr_diff) || inherits(result$dplot_corr_diff, "ggplot")
  )
})

test_that("sg-vpop-est does not work on empty dataset", {
  expect_error(sg_vpop_est(data.frame()))
})


# Basic functionality tests
test_that("sg_vpop_est returns correct structure with simple continuous data", {
  # Create larger test data to avoid duplicate warnings
  #skip_if(T)
  set.seed(42)
  # test_data <- data.frame(
  #   x1 = rnorm(200, mean = 10, sd = 2),
  #   x2 = rnorm(200, mean = 20, sd = 3)
  # )
  # fpath_i <- system.file("data", "data_pbc.rda", package = "SimuRg")
  # load(fpath_i)
  test_data <- data_pbc %>% select(years, age)

  output <-
    sg_vpop_est(data_i = test_data, diag_plots = FALSE, seed = 123, remove_duplicates = TRUE)


  # Check output structure - now a list of lists
  # expect_type(output, "list")
  # expect_length(output, 1)

  result <- output
  expect_named(result, c("datagen", "seed", "exact_dupl_check", "dplot_cont",
                         "dplot_cat", "dplot_umap", "ks_test", "jsd_res",
                         "corr_diff_mean", "corr_diff_max", "dplot_corr_diff"))

  # Check datagen
  expect_true(is.data.frame(result$datagen))
  expect_equal(ncol(result$datagen), ncol(test_data))
  expect_equal(sort(colnames(result$datagen)), sort(colnames(test_data)))

  # Check seed
  expect_equal(result$seed, 123)

  # Check exact_dupl_check
  expect_true(is.logical(result$exact_dupl_check))

  # Check ks_test
  expect_true(is.data.frame(result$ks_test))
  expect_true("variable" %in% colnames(result$ks_test))
  expect_true("p.value" %in% colnames(result$ks_test))
  expect_true("status" %in% colnames(result$ks_test))

  # Check correlation diff metrics
  expect_true(is.numeric(result$corr_diff_mean))
  expect_true(is.numeric(result$corr_diff_max))

  # JSD should be NULL for continuous-only data
  expect_null(result$jsd_res)
})

test_that("sg_vpop_est works with mixed continuous and categorical data", {
  # Create larger test data with both continuous and categorical variables
  #skip_if(T)
  set.seed(43)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3),
    cat1 = factor(sample(c("A", "B", "C", "D"), 200, replace = TRUE)),
    cat2 = factor(sample(c("X", "Y", "Z"), 200, replace = TRUE))
  )

  output <-
    sg_vpop_est(data_i = test_data, diag_plots = FALSE, seed = 123, remove_duplicates = TRUE)


  result <- output
  expect_true(is.data.frame(result$datagen))
  expect_equal(ncol(result$datagen), ncol(test_data))

  # Check that JSD is computed for categorical variables as a tibble/data frame
  expect_true(tibble::is_tibble(result$jsd_res) || is.data.frame(result$jsd_res))
  expect_true(all(c("variable", "jsd", "n_levels") %in% colnames(result$jsd_res)))
  expect_true(all(result$jsd_res$variable %in% c("cat1", "cat2")))
  expect_true(all(is.finite(result$jsd_res$jsd)))
  expect_true(all(result$jsd_res$n_levels >= 1))
  expect_true(is.numeric(attr(result$jsd_res, "jsd_weighted_mean")))

  # Check that correlation diff is computed for continuous variables
  expect_true(is.numeric(result$corr_diff_mean))
  expect_true(is.numeric(result$corr_diff_max))
})

test_that("sg_vpop_est respects nobj parameter", {
  #skip_if(T)
  set.seed(44)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3)
  )

  output <-
    sg_vpop_est(data_i = test_data, nobj = 100, diag_plots = FALSE, seed = 123, remove_duplicates = TRUE)


  result <- output
  expect_equal(nrow(result$datagen), 100)
})

test_that("sg_vpop_est respects npop parameter", {
  #skip_if(T)
  set.seed(45)
  test_data <- data.frame(
    x1 = rnorm(150, mean = 10, sd = 2),
    x2 = rnorm(150, mean = 20, sd = 3)
  )

  output <-
    sg_vpop_est(data_i = test_data, npop = 2, diag_plots = FALSE, seed = 123, remove_duplicates = TRUE)


  result <- output
  expect_equal(nrow(result$datagen), 300)
})

test_that("sg_vpop_est excludes idcol when specified", {
  #skip_if(T)
  set.seed(46)
  test_data <- data.frame(
    id = 1:200,
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3)
  )

  output <-
    sg_vpop_est(data_i = test_data, id_col = "id", diag_plots = FALSE, seed = 123, remove_duplicates = TRUE)


  result <- output
  expect_false("id" %in% colnames(result$datagen))
  expect_equal(ncol(result$datagen), 2)
})

test_that("sg_vpop_est excludes exclcol when specified", {
  #skip_if(T)
  set.seed(47)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3),
    exclude_me = rnorm(200, mean = 5, sd = 1)
  )

  output <-
    sg_vpop_est(data_i = test_data, excl_col = "exclude_me",
                diag_plots = FALSE, seed = 123, remove_duplicates = TRUE)


  result <- output
  expect_false("exclude_me" %in% colnames(result$datagen))
  expect_equal(ncol(result$datagen), 2)
})

test_that("sg_vpop_est uses seed for reproducibility", {
  #skip_if(T)
  set.seed(48)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3)
  )

  output1 <-
    sg_vpop_est(data_i = test_data, seed = 123, diag_plots = FALSE, remove_duplicates = TRUE)

  output2 <-
    sg_vpop_est(data_i = test_data, seed = 123, diag_plots = FALSE, remove_duplicates = TRUE)


  # Results should be identical with same seed
  expect_equal(output1$datagen, output2$datagen)
  expect_equal(output1$seed, 123)
  expect_equal(output2$seed, 123)
})

test_that("sg_vpop_est works with diag_plots = FALSE", {
  #skip_if(T)
  set.seed(49)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3)
  )

  output <-
    sg_vpop_est(data_i = test_data, diag_plots = FALSE, seed = 123, remove_duplicates = TRUE)


  result <- output
  expect_null(result$dplot_cont)
  expect_null(result$dplot_cat)
  expect_null(result$dplot_umap)
  expect_null(result$dplot_corr_diff)
  expect_true(is.data.frame(result$datagen))
  expect_true(is.data.frame(result$ks_test))
})

test_that("sg_vpop_est works with diag_plots = TRUE", {
  #skip_if(T)
  test_data <- data.frame(
    x1 = rnorm(300, mean = 10, sd = 2),
    x2 = rnorm(300, mean = 20, sd = 3)
  )

  output <- sg_vpop_est(data_i = test_data, diag_plots = TRUE, seed = 123, seed_umap = 123)

  result <- output
  expect_true(
    is.null(result$dplot_cont) ||
      (is.list(result$dplot_cont) &&
         all(vapply(result$dplot_cont, function(x) inherits(x, "ggplot"), logical(1))))
  )
  expect_true(
    is.null(result$dplot_cat) ||
      (is.list(result$dplot_cat) &&
         all(vapply(result$dplot_cat, function(x) inherits(x, "ggplot"), logical(1))))
  )
  expect_true(inherits(result$dplot_umap, "ggplot"))
  expect_true(inherits(result$dplot_corr_diff, "ggplot") || is.null(result$dplot_corr_diff))
})

test_that("sg_vpop_est converts character columns to factors", {
  #skip_if(T)
  set.seed(50)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    char_col = sample(c("A", "B", "C", "D", "E"), 200, replace = TRUE),
    stringsAsFactors = FALSE
  )

  output <-
    sg_vpop_est(data_i = test_data, diag_plots = FALSE, seed = 123, remove_duplicates = TRUE)


  result <- output
  # The synthetic data should have the same structure
  expect_true(is.data.frame(result$datagen))
  expect_true("char_col" %in% colnames(result$datagen))

  # Check that JSD is computed for categorical variable
  expect_true(tibble::is_tibble(result$jsd_res) || is.data.frame(result$jsd_res))
  expect_true(all(c("variable", "jsd", "n_levels") %in% colnames(result$jsd_res)))
})

test_that("sg_vpop_est handles minnumlev parameter", {
  #skip_if(T)
  # Create data with numeric column that has few unique values
  set.seed(51)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    low_level = sample(c(1, 2, 3, 4), 200, replace = TRUE)  # Only 4 unique values
  )

  output <-
    sg_vpop_est(data_i = test_data, minnumlev = 4, diag_plots = FALSE, seed = 123, remove_duplicates = TRUE)


  result <- output
  expect_true(is.data.frame(result$datagen))
  expect_equal(ncol(result$datagen), 2)

  # With minnumlev = 4, low_level should be treated as categorical
  # Check that JSD is computed
  expect_true(tibble::is_tibble(result$jsd_res) || is.data.frame(result$jsd_res))
  expect_true(all(c("variable", "jsd", "n_levels") %in% colnames(result$jsd_res)))
})

test_that("sg_vpop_est handles NA values by removing rows", {
  #skip_if(T)
  set.seed(52)
  test_data <- data.frame(
    x1 = c(rnorm(195, mean = 10, sd = 2), rep(NA, 5)),
    x2 = rnorm(200, mean = 20, sd = 3)
  )

  # Should not error, but will remove rows with NA
  output <-
    sg_vpop_est(data_i = test_data, diag_plots = FALSE, seed = 123, remove_duplicates = TRUE)


  result <- output
  expect_true(is.data.frame(result$datagen))
  expect_equal(nrow(result$datagen), 195)  # Should match number of non-NA rows
})

# New tests for updated features

test_that("sg_vpop_est fixed seed mode generates single dataset", {
  #skip_if(T)
  set.seed(53)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3)
  )

  output <-
    sg_vpop_est(data_i = test_data, seed = 456, diag_plots = FALSE, remove_duplicates = TRUE)


  # Should generate exactly 1 dataset in fixed seed mode

  expect_equal(output$seed, 456)
})

test_that("sg_vpop_est search mode generates multiple datasets", {
  #skip_if(T)
  set.seed(54)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3),
    x3 = rnorm(200, mean = 5, sd = 1)
  )

  # Search mode with nds = 3
  output <-
    sg_vpop_est(data_i = test_data, seed = NA, nds = 3,
                tg_corrdif = 0.15, diag_plots = FALSE, remove_duplicates = TRUE)


  # Should generate 3 datasets
  expect_length(output, 3)

  # Each dataset should have different seeds
  seeds <- sapply(output, function(x) x$seed)
  expect_true(all(!is.na(seeds)))

  # Check that all datasets meet structure requirements
  for (i in 1:3) {
    expect_true(is.data.frame(output[[i]]$datagen))
    expect_true(is.numeric(output[[i]]$corr_diff_max))
  }
})

test_that("sg_vpop_est computes correlation difference metrics", {
  #skip_if(T)
  set.seed(55)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3),
    x3 = rnorm(200, mean = 5, sd = 1)
  )

  output <-
    sg_vpop_est(data_i = test_data, seed = 789, diag_plots = FALSE, remove_duplicates = TRUE)


  result <- output

  # Check correlation metrics
  expect_true(is.numeric(result$corr_diff_mean))
  expect_true(is.numeric(result$corr_diff_max))
  expect_true(result$corr_diff_mean >= 0)
  expect_true(result$corr_diff_max >= result$corr_diff_mean)
})

test_that("sg_vpop_est computes JSD for categorical variables", {
  #skip_if(T)
  set.seed(56)
  test_data <- data.frame(
    cat1 = factor(sample(c("A", "B", "C", "D"), 200, replace = TRUE)),
    cat2 = factor(sample(c("X", "Y", "Z"), 200, replace = TRUE)),
    x1 = rnorm(200, mean = 10, sd = 2)
  )

  output <-
    sg_vpop_est(data_i = test_data, seed = 321, diag_plots = FALSE, remove_duplicates = TRUE)


  result <- output

  # Check JSD metric
  expect_true(tibble::is_tibble(result$jsd_res) || is.data.frame(result$jsd_res))
  expect_true(all(c("variable", "jsd", "n_levels") %in% colnames(result$jsd_res)))
  expect_true(all(result$jsd_res$jsd >= 0))
  expect_true(all(result$jsd_res$jsd <= 1))
})

test_that("sg_vpop_est checks for exact duplicates", {
  #skip_if(T)
  set.seed(57)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3)
  )

  output <-
    sg_vpop_est(data_i = test_data, seed = 111,
                diag_plots = FALSE, remove_duplicates = TRUE)


  result <- output

  # Should have checked for duplicates
  expect_true(is.logical(result$exact_dupl_check))
})

test_that("sg_vpop_est generates correct plots with diag_plots = TRUE", {
  #skip_if(T)
  set.seed(58)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3),
    cat1 = factor(sample(c("A", "B", "C"), 200, replace = TRUE))
  )

  output <-
    sg_vpop_est(data_i = test_data, seed = 555,
                diag_plots = TRUE, seed_umap = 42, remove_duplicates = TRUE)


  result <- output

  # Check continuous plots
  expect_true(is.list(result$dplot_cont))
  expect_length(result$dplot_cont, 2)  # 2 continuous variables
  expect_true(all(vapply(result$dplot_cont,
                         function(x) inherits(x, "ggplot"), logical(1))))

  # Check categorical plots
  expect_true(is.list(result$dplot_cat))
  expect_length(result$dplot_cat, 1)  # 1 categorical variable
  expect_true(all(vapply(result$dplot_cat,
                         function(x) inherits(x, "ggplot"), logical(1))))

  # Check UMAP plot
  expect_true(inherits(result$dplot_umap, "ggplot"))

  # Check correlation difference heatmap
  expect_true(inherits(result$dplot_corr_diff, "ggplot"))
})

test_that("sg_vpop_est respects tg_corrdif in search mode", {
  #skip_if(T)
  set.seed(59)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3),
    x3 = rnorm(200, mean = 5, sd = 1)
  )

  # Use a stricter correlation threshold
  output <-
    sg_vpop_est(data_i = test_data, seed = NA, nds = 2,
                tg_corrdif = 0.05, diag_plots = FALSE, remove_duplicates = TRUE)


  # Check that generated datasets meet (or try to meet) the threshold
  for (i in 1:length(output)) {
    result <- output[[i]]
    # If seed is not NA, it means it found a dataset
    if (!is.na(result$seed)) {
      expect_true(is.numeric(result$corr_diff_max))
    }
  }
})

test_that("sg_vpop_est handles categorical-only data", {
  #skip_if(T)
  set.seed(60)
  test_data <- data.frame(
    cat1 = factor(sample(c("A", "B", "C", "D"), 200, replace = TRUE)),
    cat2 = factor(sample(c("X", "Y", "Z"), 200, replace = TRUE)),
    cat3 = factor(sample(c("Low", "Med", "High"), 200, replace = TRUE))
  )

  output <- sg_vpop_est(data_i = test_data, seed = 999, diag_plots = FALSE, remove_duplicates = TRUE)


  result <- output

  # Check that it works with only categorical data
  expect_true(is.data.frame(result$datagen))
  expect_equal(ncol(result$datagen), 3)

  # JSD should be computed
  expect_true(tibble::is_tibble(result$jsd_res) || is.data.frame(result$jsd_res))
  expect_true(all(c("variable", "jsd", "n_levels") %in% colnames(result$jsd_res)))

  # Correlation metrics should be NULL (no continuous variables)
  expect_null(result$corr_diff_mean)
  expect_null(result$corr_diff_max)
  expect_null(result$ks_test)
})

test_that("sg_vpop_est noise_level parameter affects duplicate removal", {
  #skip_if(T)
  set.seed(61)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3)
  )

  # Test with different noise levels
  output1 <-
    sg_vpop_est(data_i = test_data, seed = 123,
                remove_duplicates = TRUE, noise_level = 0.05,
                diag_plots = FALSE)


  output2 <-
    sg_vpop_est(data_i = test_data, seed = 123,
                remove_duplicates = TRUE, noise_level = 0.20,
                diag_plots = FALSE)


  # Both should complete without error
  expect_true(is.data.frame(output1$datagen))
  expect_true(is.data.frame(output2$datagen))
})

test_that("compare_cor_matrices returns expected structure and metrics", {
  data_obs <- data.frame(
    a = c(1, 2, 3, 4, 5),
    b = c(2, 4, 5, 7, 9),
    c = c(5, 4, 4, 2, 1)
  )
  data_syn <- data.frame(
    a = c(1, 2, 3, 4, 5),
    b = c(2, 3, 5, 8, 9),
    c = c(4, 4, 3, 2, 1)
  )
  vars <- c("a", "b", "c")

  res <- compare_cor_matrices(data_obs, data_syn, vars = vars, method = "pearson")

  expect_type(res, "list")
  expect_true(all(c("mean_abs_diff", "max_abs_diff", "R_obs", "R_syn", "R_diff") %in% names(res)))
  expect_equal(dim(res$R_obs), c(3, 3))
  expect_equal(dim(res$R_syn), c(3, 3))
  expect_equal(res$R_diff, res$R_syn - res$R_obs)

  idx <- upper.tri(res$R_obs, diag = FALSE)
  expected_diff <- abs(res$R_obs[idx] - res$R_syn[idx])
  expect_equal(res$mean_abs_diff, mean(expected_diff))
  expect_equal(res$max_abs_diff, max(expected_diff))
})

test_that("create_optimal_visit_sequence handles edge cases and ordering", {
  x <- data.frame(
    c1 = c(1, 2, 3, 4, 5, 6),
    c2 = c(1, 2, 2, 4, 5, 6),
    c3 = c(6, 1, 5, 2, 4, 3),
    cat1 = factor(c("A", "B", "A", "B", "A", "B"))
  )

  expect_warning(
    create_optimal_visit_sequence(x, var_cont = character(0), var_cat = character(0)),
    "No variables provided"
  )

  expect_equal(
    create_optimal_visit_sequence(x, var_cont = character(0), var_cat = "cat1"),
    "cat1"
  )

  expect_equal(
    create_optimal_visit_sequence(x, var_cont = "c1", var_cat = "cat1"),
    c("c1", "cat1")
  )

  seq_all <- create_optimal_visit_sequence(
    x,
    var_cont = c("c1", "c2", "c3"),
    var_cat = "cat1"
  )
  expect_equal(tail(seq_all, 1), "cat1")
  expect_equal(head(seq_all, 3), c("c2", "c1", "c3"))
})

test_that("remove_exact_duplicates reports and perturbs matching rows", {
  data_orig <- data.frame(
    x = c(0, 1, 2, 3, 4, 5),
    y = c(10, 12, 14, 16, 18, 20),
    cat1 = factor(c("A", "A", "B", "B", "C", "C"))
  )
  data_syn <- data.frame(
    x = c(0, 1, 9),
    y = c(10, 12, 99),
    cat1 = factor(c("A", "A", "C"), levels = levels(data_orig$cat1))
  )

  res <- remove_exact_duplicates(
    data_syn = data_syn,
    data_orig = data_orig,
    var_cont = c("x", "y"),
    var_cat = "cat1",
    noise_level = 0.2,
    seed = 777
  )

  expect_true(is.data.frame(res$data_cleaned))
  expect_equal(res$n_duplicates_removed, 2)
  expect_equal(sort(res$duplicate_indices), c(1, 2))

  common_cols <- intersect(names(res$data_cleaned), names(data_orig))
  expect_equal(nrow(dplyr::semi_join(res$data_cleaned, data_orig, by = common_cols)), 0)



  dup_idx <- res$duplicate_indices
  expect_true(all(res$data_cleaned$x[dup_idx] >= min(data_orig$x) & res$data_cleaned$x[dup_idx] <= max(data_orig$x)))
  expect_true(all(res$data_cleaned$y[dup_idx] >= min(data_orig$y) & res$data_cleaned$y[dup_idx] <= max(data_orig$y)))
  expect_equal(levels(res$data_cleaned$cat1), levels(data_syn$cat1))
})

test_that("sg_vpop_est return shape differs by fixed-seed and search mode", {
  set.seed(62)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3)
  )

  fixed_output <- sg_vpop_est(
    data_i = test_data,
    seed = 321,
    diag_plots = FALSE,
    remove_duplicates = TRUE
  )
  expect_true(is.data.frame(fixed_output$datagen))
  expect_false(is.list(fixed_output[[1]]) && "datagen" %in% names(fixed_output[[1]]))

  search_output <- sg_vpop_est(
    data_i = test_data,
    seed = NA,
    nds = 2,
    tg_corrdif = 0.2,
    diag_plots = FALSE,
    remove_duplicates = TRUE
  )
  expect_length(search_output, 2)
  expect_true(all(vapply(search_output, function(x) is.list(x) && is.data.frame(x$datagen), logical(1))))
})

test_that("sg_vpop_est enforces tg_corrdif when search seed is found", {
  set.seed(63)
  test_data <- data.frame(
    x1 = rnorm(200, mean = 10, sd = 2),
    x2 = rnorm(200, mean = 20, sd = 3),
    x3 = rnorm(200, mean = 5, sd = 1)
  )

  output <- sg_vpop_est(
    data_i = test_data,
    seed = NA,
    nds = 2,
    tg_corrdif = 0.1,
    diag_plots = FALSE,
    remove_duplicates = TRUE
  )

  for (i in seq_along(output)) {
    if (!is.na(output[[i]]$seed)) {
      expect_lte(output[[i]]$corr_diff_max, 0.1)
    }
  }
})

test_that("sg_vpop_est warns for missing id_col and excl_col values", {
  set.seed(64)
  test_data <- data.frame(
    x1 = rnorm(150, mean = 10, sd = 2),
    x2 = rnorm(150, mean = 20, sd = 3)
  )

  expect_warning(
    sg_vpop_est(
      data_i = test_data,
      id_col = "missing_id",
      diag_plots = FALSE,
      seed = 123,
      remove_duplicates = TRUE
    ),
    "id_col"
  )

  expect_warning(
    sg_vpop_est(
      data_i = test_data,
      excl_col = c("missing_excluded"),
      diag_plots = FALSE,
      seed = 123,
      remove_duplicates = TRUE
    ),
    "excl_col"
  )
})

test_that("sg_vpop_est handles integer columns without casting failures", {
  set.seed(67)
  test_data <- data.frame(
    x_int = sample(1:8, 200, replace = TRUE),
    x_num = rnorm(200, mean = 30, sd = 5)
  )

  output <- sg_vpop_est(
    data_i = test_data,
    seed = 123,
    diag_plots = FALSE,
    remove_duplicates = TRUE
  )

  expect_true(is.data.frame(output$datagen))
  expect_true(all(c("x_int", "x_num") %in% names(output$datagen)))
  expect_true(is.numeric(output$datagen$x_int))
})

test_that("sg_vpop_est categorical-only with diag_plots TRUE returns only categorical plots", {
  set.seed(69)
  test_data <- data.frame(
    cat1 = factor(sample(c("A", "B", "C"), 220, replace = TRUE)),
    cat2 = factor(sample(c("X", "Y"), 220, replace = TRUE))
  )

  output <- sg_vpop_est(
    data_i = test_data,
    seed = 123,
    diag_plots = TRUE,
    seed_umap = 42,
    remove_duplicates = TRUE
  )

  expect_null(output$dplot_cont)
  expect_true(is.list(output$dplot_cat))
  expect_null(output$dplot_corr_diff)
  expect_true(inherits(output$dplot_umap, "ggplot"))
})

test_that("remove_exact_duplicates noise_level increases perturbation magnitude", {
  data_orig <- data.frame(
    x = seq(1, 30, by = 1),
    y = seq(50, 79, by = 1),
    cat1 = factor(rep(c("A", "B", "C"), length.out = 30))
  )
  data_syn <- data.frame(
    x = c(5, 6, 7, 400),
    y = c(54, 55, 56, 900),
    cat1 = factor(c("B", "C", "A", "C"), levels = levels(data_orig$cat1))
  )

  low_noise <- remove_exact_duplicates(
    data_syn = data_syn,
    data_orig = data_orig,
    var_cont = c("x", "y"),
    var_cat = "cat1",
    noise_level = 0.01,
    seed = 123
  )

  high_noise <- remove_exact_duplicates(
    data_syn = data_syn,
    data_orig = data_orig,
    var_cont = c("x", "y"),
    var_cat = "cat1",
    noise_level = 0.50,
    seed = 123
  )

  dupl_idx <- low_noise$duplicate_indices
  low_shift <- mean(abs(low_noise$data_cleaned$x[dupl_idx] - data_syn$x[dupl_idx]))
  high_shift <- mean(abs(high_noise$data_cleaned$x[dupl_idx] - data_syn$x[dupl_idx]))

  expect_gt(high_shift, low_shift)
})

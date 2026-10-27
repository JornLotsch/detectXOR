test_that("the designed XOR pair in XOR_data is the only detection", {
  set.seed(1)
  result <- run_xor(XOR_data)
  expect_equal(detected_pairs(result), "Variable_C||Variable_D")
  expect_equal(result$results_df["Variable_C||Variable_D", "xor_pattern"], "default")
  expect_equal(rownames(result$results_df)[1], "Variable_A||Variable_B")
})

test_that("a class column not named 'class' is handled", {
  data <- XOR_data
  names(data)[names(data) == "class"] <- "Group"
  set.seed(1)
  expect_equal(detected_pairs(run_xor(data, class_col = "Group")), "Variable_C||Variable_D")
})

test_that("detection does not depend on class coding", {
  swapped <- XOR_data
  swapped$class <- 3 - swapped$class
  labelled <- XOR_data
  labelled$class <- ifelse(labelled$class == 1, "b", "a")
  set.seed(1)
  original <- run_xor(XOR_data)$results_df
  set.seed(1)
  expect_equal(run_xor(swapped)$results_df$xor_shape_detected, original$xor_shape_detected)
  set.seed(1)
  result_labelled <- run_xor(labelled)$results_df
  expect_equal(result_labelled$xor_shape_detected, original$xor_shape_detected)
  expect_equal(result_labelled$xor_purity, original$xor_purity)
})

test_that("extreme_handling accepts 'trim' and its synonym 'remove'", {
  set.seed(1)
  expect_equal(detected_pairs(run_xor(XOR_data, extreme_handling = "trim")), "Variable_C||Variable_D")
  set.seed(1)
  expect_equal(detected_pairs(run_xor(XOR_data, extreme_handling = "remove")), "Variable_C||Variable_D")
  expect_error(run_xor(XOR_data, extreme_handling = "drop"))
})

test_that("rotated XOR patterns are detected as rotated", {
  set.seed(2)
  result <- run_xor(simulate_xor(rotated = TRUE))
  expect_equal(detected_pairs(result), "X||Y")
  expect_equal(result$results_df["X||Y", "xor_pattern"], "rotated")
})

test_that("binary XOR data are detected", {
  set.seed(3)
  n <- 400
  data <- data.frame(A = stats::rbinom(n, 1, 0.5), B = stats::rbinom(n, 1, 0.5), C = stats::rbinom(n, 1, 0.5))
  data$class <- ifelse(xor(data$A, data$B), 2, 1)
  result <- run_xor(data, extreme_handling = "none", compute_axes_parallel_significance = TRUE)
  expect_equal(detected_pairs(result), "A||B")
  expect_true(is.na(result$results_df["A||B", "axis_tests_passed"]))
})

test_that("incomplete XOR patterns are rejected", {
  set.seed(4)
  result <- run_xor(utils::head(XOR_data, 250))
  expect_length(detected_pairs(result), 0)
  expect_false(result$results_df["Variable_C||Variable_D", "coverage_ok"])
})

test_that("detection of a moderately noisy XOR does not deteriorate with sample size", {
  for (n in c(400, 2000)) {
    set.seed(5)
    expect_equal(detected_pairs(run_xor(simulate_xor(n = n, label_noise = 0.05))), "X||Y")
  }
})

test_that("no XOR is detected in data without structure", {
  set.seed(6)
  data <- data.frame(matrix(stats::rnorm(400 * 6), 400), class = rep(1:2, each = 200))
  expect_length(detected_pairs(run_xor(data)), 0)
})

test_that("the Wilcoxon criterion enters the decision when requested", {
  set.seed(7)
  result <- run_xor(XOR_data, compute_axes_parallel_significance = TRUE)
  expect_true(result$results_df["Variable_C||Variable_D", "axis_tests_passed"])
  expect_equal(detected_pairs(result), "Variable_C||Variable_D")

  # A failed Wilcoxon criterion overrides the tile and tau criteria
  pair_list <- result$pair_list
  pair_list[["Variable_C||Variable_D"]]$results$tile_wilcox_significance$pass_groups <- FALSE
  df <- create_results_df(pair_list, result$settings)
  expect_false(df["Variable_C||Variable_D", "xor_shape_detected"])
})

test_that("the tau criterion enters the decision when requested", {
  set.seed(8)
  result <- run_xor(XOR_data)
  pair_list <- result$pair_list
  pair_list[["Variable_C||Variable_D"]]$results$classwise_tau$default$is_strong_reverse_correlated <- FALSE
  df <- create_results_df(pair_list, result$settings)
  expect_false(df["Variable_C||Variable_D", "xor_shape_detected"])
})

test_that("results are reproducible with set.seed", {
  set.seed(9)
  first <- run_xor(XOR_data)$results_df
  set.seed(9)
  second <- run_xor(XOR_data)$results_df
  expect_identical(first, second)
})

test_that("p-value adjustment methods are available", {
  set.seed(10)
  maxT <- run_xor(XOR_data)$results_df
  set.seed(10)
  none <- run_xor(XOR_data, p_adjust_method = "none")$results_df
  expect_equal(none$p_adjusted, none$p_value)
  expect_true(all(maxT$p_adjusted >= maxT$p_value))
})

test_that("deprecated abs_diff_threshold gives a warning", {
  set.seed(11)
  expect_warning(run_xor(XOR_data, abs_diff_threshold = 20), "deprecated")
})

test_that("parallel and sequential processing give identical results", {
  skip_on_cran()
  skip_if_not_installed("pbmcapply")
  skip_on_os("windows")
  set.seed(12)
  sequential <- run_xor(XOR_data, max_cores = 1)$results_df
  set.seed(12)
  parallel <- suppressWarnings(run_xor(XOR_data, max_cores = 2))$results_df
  expect_identical(sequential, parallel)
})

test_that("tile boundaries can be changed", {
  set.seed(13)
  result <- run_xor(XOR_data, split_quantiles = c(0.25, 0.75))
  expect_equal(detected_pairs(result), "Variable_C||Variable_D")
  expect_equal(result$settings$split_quantiles, c(0.25, 0.75))
  expect_error(run_xor(XOR_data, split_quantiles = c(0.7, 0.3)))
})

test_that("binary XOR with unbalanced value frequencies is detected", {
  set.seed(14)
  n <- 400
  # Dominant model: value 1 in 75% of cases, so the (0, 0) cell is rare
  data <- data.frame(A = stats::rbinom(n, 1, 0.75), B = stats::rbinom(n, 1, 0.75), C = stats::rbinom(n, 1, 0.75))
  data$class <- ifelse(xor(data$A, data$B), 2, 1)
  expect_equal(detected_pairs(run_xor(data, extreme_handling = "none")), "A||B")
})

test_that("class balancing restores detection of rotated XOR under imbalance", {
  set.seed(15)
  n <- 400
  cls <- sample(c(rep(1, 100), rep(2, 300)))
  a <- stats::rbinom(n, 1, 0.5)
  b <- ifelse(cls == 1, a, 1 - a)
  x <- 3 * a + stats::rnorm(n)
  y <- 3 * b + stats::rnorm(n)
  data <- data.frame(X = x + y, Y = x - y, Z = stats::rnorm(n), class = cls)
  set.seed(16)
  result <- run_xor(data, balance_classes = TRUE)
  expect_equal(detected_pairs(result), "X||Y")
  expect_true(result$settings$balance_classes)
})

test_that("purity of a perfect XOR table is 1 and invariant to label assignment", {
  tiles <- factor(c("Lower Left", "Upper Right", "Upper Left", "Lower Right"), levels = TILE_LEVELS)
  class2 <- matrix(c(0, 0, 1, 1))
  counts <- tile_class_counts(tiles, class2)
  expect_equal(compute_xor_balanced_accuracy(counts, "default"), 1)
  expect_equal(compute_xor_balanced_accuracy(tile_class_counts(tiles, 1 - class2), "default"), 0)
  expect_equal(compute_max_xor_purity(tile_class_counts(tiles, 1 - class2)), 1)
})

test_that("purity is 0.5 when a class is absent from the critical tiles", {
  tiles <- factor(c("Lower Left", "Upper Right", "Center", "Center"), levels = TILE_LEVELS)
  counts <- tile_class_counts(tiles, matrix(c(0, 0, 1, 1)))
  expect_equal(compute_xor_balanced_accuracy(counts, "default"), 0.5)
})

test_that("the purity lower bound is below the estimate and increases with n", {
  small <- purity_lower_bound(10, 10, 10, 10)
  large <- purity_lower_bound(100, 100, 100, 100)
  expect_lt(small, 1)
  expect_lt(small, large)
})

test_that("binary variables are split into low and high tiles", {
  data <- data.frame(a = c(0, 0, 1, 1), b = c(0, 1, 0, 1))
  tiles <- assign_3x3_tiles(data, "a", "b")
  expect_equal(as.character(tiles), c("Lower Left", "Upper Left", "Lower Right", "Upper Right"))
})

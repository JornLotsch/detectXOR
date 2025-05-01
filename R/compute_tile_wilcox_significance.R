#' Compute group-wise Wilcoxon significance for tile groups
#'
#' @param orig_pair_list List. Output from compute_tile_patterns_for_pairs.
#' @param class_col Character. Name of the class column (default: "class").
#' @param p_threshold Numeric. Significance threshold (default: 0.05).
#' @param split_method Character. Tile assignment method ("quantile"/"range") (default: "quantile").
#' @return List with results$tile_wilcox_significance added for each pair.
#' @export
compute_tile_wilcox_significance_for_pairs <- function(
    orig_pair_list,
    class_col = "class",
    p_threshold = 0.05,
    split_method = "quantile"  # Added for future use
) {
  for (nm in names(orig_pair_list)) {
    pair <- orig_pair_list[[nm]]
    df <- pair$data

    var_names <- setdiff(names(df), class_col)
    if (length(var_names) != 2) stop("Data must contain exactly 2 variables + class column.")

    tile_pattern <- pair$results$tile_pattern
    if (is.null(tile_pattern)) stop("Run compute_tile_patterns_for_pairs first.")

    tile_wilcox_result <- check_tile_groupwise_wilcox_significance(
      data = df,
      var1_name = var_names[1],
      var2_name = var_names[2],
      class_name = class_col,
      chosen_tiles = colnames(tile_pattern$observed_table),
      scenario = tile_pattern$xor_pattern,
      p_threshold = p_threshold,
      split_method = split_method  # Forwarded (even if unused)
    )

    orig_pair_list[[nm]]$results$tile_wilcox_significance <- tile_wilcox_result
  }
  orig_pair_list
}

#' Check group-wise Wilcoxon significance for tile groups (for a single pair)
#' Uses scenario and chosen tiles as determined by tile_pattern.
#' @keywords internal
check_tile_groupwise_wilcox_significance <- function(
    data, var1_name, var2_name, class_name, chosen_tiles, scenario,
    p_threshold = 0.05,
    split_method = "quantile"  # Added for consistency
) {
  # Helper to assign tiles
  assign_3x3_tiles <- function(data, var1_name, var2_name, split_method = "quantile") {
    if (split_method == "quantile") {
      x_limits <- quantile(data[[var1_name]], probs = c(1 / 3, 2 / 3), na.rm = TRUE)
      y_limits <- quantile(data[[var2_name]], probs = c(1 / 3, 2 / 3), na.rm = TRUE)
    } else {
      x_limits <- split_range_three(data[[var1_name]])
      y_limits <- split_range_three(data[[var2_name]])
    }
    x_group <- cut(data[[var1_name]], c(-Inf, x_limits, Inf), labels = c("Left", "Mid", "Right"), include.lowest = TRUE)
    y_group <- cut(data[[var2_name]], c(-Inf, y_limits, Inf), labels = c("Lower", "Mid", "Upper"), include.lowest = TRUE)
    tile <- mapply(function(x, y) {
      switch(paste(x, y),
             "Left Upper" = "Upper Left",
             "Right Upper" = "Upper Right",
             "Left Lower" = "Lower Left",
             "Right Lower" = "Lower Right",
             "Mid Upper" = "Upper Mid",
             "Left Mid" = "Mid Left",
             "Right Mid" = "Mid Right",
             "Mid Lower" = "Lower Mid",
             "Mid Mid" = "Center",
             NA)
    }, x_group, y_group)
    levels_order <- c(
      "Upper Left", "Upper Mid", "Upper Right",
      "Mid Left", "Center", "Mid Right",
      "Lower Left", "Lower Mid", "Lower Right"
    )
    factor(tile, levels = levels_order)
  }
  split_range_three <- function(x) {
    min_val <- min(x, na.rm = TRUE)
    max_val <- max(x, na.rm = TRUE)
    break1 <- min_val + (max_val - min_val) / 3
    break2 <- min_val + 2 * (max_val - min_val) / 3
    return(c(break1, break2))
  }

  # Assign tiles (once, for all data)
  tiles_population <- assign_3x3_tiles(data, var1_name, var2_name)
  # Subset data to chosen tiles
  data_tiles <- data[tiles_population %in% chosen_tiles, ]
  data_tiles$tile <- factor(tiles_population[tiles_population %in% chosen_tiles], levels = chosen_tiles)
  data_tiles$class <- as.factor(data_tiles[[class_name]])

  # Helper to run group-wise combined Wilcoxon test
  run_combined_wilcox <- function(df, tiles, var1, var2, class_col) {
    idx <- df$tile %in% tiles
    values <- c(df[idx, var1], df[idx, var2])
    classes <- c(df[idx, class_col], df[idx, class_col])
    if (length(unique(classes)) == 2 && length(values) > 1) {
      test <- tryCatch(
        wilcox.test(values ~ classes),
        error = function(e) NULL
      )
      if (!is.null(test)) {
        means <- tapply(values, classes, mean, na.rm = TRUE)
        medians <- tapply(values, classes, median, na.rm = TRUE)
        direction <- if (!is.null(means)) names(means)[which.max(means)] else NA
        return(list(
          p_value = test$p.value,
          statistic = test$statistic,
          means = means,
          medians = medians,
          direction = direction
        ))
      }
    }
    list(p_value = NA, statistic = NA, means = NA, medians = NA, direction = NA)
  }

  # Prepare groupings and possibly rotate data if scenario is "rotated"
  if (scenario == "default") {
    left_tiles <- c("Lower Left", "Upper Left")
    right_tiles <- c("Lower Right", "Upper Right")
    upper_tiles <- c("Upper Left", "Upper Right")
    lower_tiles <- c("Lower Left", "Lower Right")
    dq <- data_tiles
    v1 <- var1_name
    v2 <- var2_name
  } else {
    dq <- data_tiles
    dq$rot_plus <- dq[[var1_name]] + dq[[var2_name]]
    dq$rot_minus <- dq[[var1_name]] - dq[[var2_name]]
    dq$tile_rot <- NA
    dq$tile_rot[dq$tile == "Mid Left"] <- "Upper Left"
    dq$tile_rot[dq$tile == "Upper Mid"] <- "Upper Right"
    dq$tile_rot[dq$tile == "Mid Right"] <- "Lower Right"
    dq$tile_rot[dq$tile == "Lower Mid"] <- "Lower Left"
    dq$tile_rot <- factor(dq$tile_rot, levels = c("Lower Left", "Lower Right", "Upper Left", "Upper Right"))
    left_tiles <- c("Lower Left", "Upper Left")
    right_tiles <- c("Lower Right", "Upper Right")
    upper_tiles <- c("Upper Left", "Upper Right")
    lower_tiles <- c("Lower Left", "Lower Right")
    dq$tile <- dq$tile_rot
    v1 <- "rot_plus"
    v2 <- "rot_minus"
  }

  test_left <- run_combined_wilcox(dq, left_tiles, v1, v2, class_name)
  test_right <- run_combined_wilcox(dq, right_tiles, v1, v2, class_name)
  test_upper <- run_combined_wilcox(dq, upper_tiles, v1, v2, class_name)
  test_lower <- run_combined_wilcox(dq, lower_tiles, v1, v2, class_name)

  summary_combined <- function(test, group) {
    if (!is.na(test$p_value)) {
      paste0(
        group, ": class 1 mean = ", signif(test$means[1], 3),
        ", class 2 mean = ", signif(test$means[2], 3),
        ", p = ", signif(test$p_value, 3),
        ", higher = class ", test$direction
      )
    } else {
      paste0(group, ": insufficient data or only one class present")
    }
  }
  summary_left <- summary_combined(test_left, "Left")
  summary_right <- summary_combined(test_right, "Right")
  summary_upper <- summary_combined(test_upper, "Upper")
  summary_lower <- summary_combined(test_lower, "Lower")

  pvals_groups <- c(test_left$p_value, test_right$p_value, test_upper$p_value, test_lower$p_value)
  pass_groups <- length(pvals_groups) == 4 && all(!is.na(pvals_groups)) && all(pvals_groups < p_threshold)
  tile_group_message <- if (pass_groups) {
    "All four group-wise combined tests are significant: XOR pattern is statistically supported."
  } else {
    "Not all four group-wise combined tests are significant or present: XOR pattern is not statistically supported."
  }

  list(
    pass_groups = pass_groups,
    summary_left = summary_left,
    summary_right = summary_right,
    summary_upper = summary_upper,
    summary_lower = summary_lower,
    tile_group_message = tile_group_message,
    pvals_groups = pvals_groups,
    scenario = scenario,
    chosen_tiles = chosen_tiles
  )
}

# pairwise_dfs <- compute_tile_wilcox_significance_for_pairs(pairwise_dfs)

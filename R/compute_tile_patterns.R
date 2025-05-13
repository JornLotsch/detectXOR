#' Compute tile pattern analysis for all pairwise datasets
#'
#' This function iterates over a list of pairwise datasets (as produced by create_pairwise_datasets)
#' and fills the results$tile_pattern entry for each pair using the provided tile pattern analysis.
#'
#' Assign 3x3 tiles to data points based on two variables
#' @keywords internal
assign_3x3_tiles <- function(data, var1_name, var2_name, split_method = "quantile") {
  if (split_method == "quantile") {
    xq <- quantile(data[[var1_name]], probs = c(1/3, 2/3), na.rm = TRUE)
    yq <- quantile(data[[var2_name]], probs = c(1/3, 2/3), na.rm = TRUE)
    
    # Check if quantile limits are unique for x
    if (length(unique(xq)) < 2) {
      x_limits <- split_range_three(data[[var1_name]])
    } else {
      x_limits <- xq
    }
    
    # Check if quantile limits are unique for y
    if (length(unique(yq)) < 2) {
      y_limits <- split_range_three(data[[var2_name]])
    } else {
      y_limits <- yq
    }
    
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

#' Split a numeric vector into three equal-width bins
#' @keywords internal
split_range_three <- function(x) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  break1 <- min_val + (max_val - min_val) / 3
  break2 <- min_val + 2 * (max_val - min_val) / 3
  return(c(break1, break2))
}

#' Chi-squared test and difference check for observed vs expected counts
#' @keywords internal
check_chisquared_difference <- function(
    observed_counts, expected_counts,
    abs_diff_threshold = 20,
    filter_expected = 0,
    p_threshold = 0.05,
    monte_carlo_B = 0
) {
  observed_counts <- as.matrix(observed_counts)
  expected_counts <- as.matrix(expected_counts)
  expected_counts[expected_counts == 0] <- 0.5
  observed_vector <- as.vector(t(observed_counts))
  expected_vector <- as.vector(t(expected_counts))
  mask <- expected_vector > filter_expected
  observed_vector_filt <- observed_vector[mask]
  expected_vector_filt <- expected_vector[mask]
  abs_difference_sum <- sum(abs(observed_vector - expected_vector))
  # Standard chi-square test
  if (length(observed_vector_filt) > 1) {
    chisq_test_result <- suppressWarnings(
      chisq.test(
        x = observed_vector_filt,
        p = expected_vector_filt / sum(expected_vector_filt),
        rescale.p = TRUE
      )
    )
    p_value <- chisq_test_result$p.value
    statistic <- chisq_test_result$statistic
    df <- chisq_test_result$parameter
    # Monte Carlo version
    chisq_mc_result <- suppressWarnings(
      chisq.test(
        x = observed_vector_filt,
        p = expected_vector_filt / sum(expected_vector_filt),
        rescale.p = TRUE,
        simulate.p.value = TRUE,
        B = monte_carlo_B
      )
    )
    p_value_mc <- chisq_mc_result$p.value
    statistic_mc <- chisq_mc_result$statistic
  } else {
    p_value <- NA
    statistic <- NA
    df <- NA
    p_value_mc <- NA
    statistic_mc <- NA
  }
  differs <- (!is.na(p_value) && p_value < p_threshold) || abs_difference_sum > abs_diff_threshold
  differs_mc <- (!is.na(p_value_mc) && p_value_mc < p_threshold) || abs_difference_sum > abs_diff_threshold
  list(
    differs = differs,
    p_value = p_value,
    statistic = statistic,
    df = df,
    differs_mc = differs_mc,
    p_value_mc = p_value_mc,
    statistic_mc = statistic_mc,
    abs_difference_sum = abs_difference_sum,
    observed = observed_counts,
    expected = expected_counts,
    message = if (differs | differs_mc) {
      "Observed counts differ from expectations (statistically or practically)."
    } else {
      "No significant or practical difference detected."
    }
  )
}

#' Main tile pattern analysis for a single pair
#' @keywords internal
analyze_tile_dominance <- function(
    data, var1_name, var2_name, class_name,
    p_threshold = 0.05,
    abs_diff_threshold = 20,
    split_method = "quantile"
) {
  TILES_DEFAULT <- c("Lower Left", "Lower Right", "Upper Left", "Upper Right")
  TILES_ROTATED <- c("Mid Left", "Upper Mid", "Lower Mid", "Mid Right")
  IDEAL_PATTERN <- c(1, 2, 2, 1)
  # Assign tiles and create observed table
  tiles_population <- assign_3x3_tiles(data, var1_name, var2_name)
  observed_table <- table(data[[class_name]], tiles_population)
  # Helper to compute dominance
  compute_dominance <- function(obs_table, tiles, ideal_pattern) {
    obs_sub <- obs_table[, tiles, drop = FALSE]
    ideal_cells <- mapply(function(class, col) obs_sub[class, col], ideal_pattern, 1:4)
    dominance <- sum(ideal_cells) / sum(obs_sub)
    dominance
  }
  # Helper to build expected table for scenario
  create_expected_table_2x4 <- function(obs_table, tiles, ideal_pattern, strict = FALSE) {
    obs_sub <- obs_table[, tiles, drop = FALSE]
    expected_table <- matrix(0, nrow = 2, ncol = 4, dimnames = list(rownames(obs_table), tiles))
    if (strict) {
      total <- sum(obs_sub)
      for (i in seq_along(ideal_pattern)) {
        expected_class <- ideal_pattern[i]
        expected_table[expected_class, i] <- total / 4
      }
    } else {
      tile_totals <- colSums(obs_sub)
      for (i in seq_along(ideal_pattern)) {
        expected_class <- ideal_pattern[i]
        expected_table[expected_class, i] <- tile_totals[i]
      }
    }
    expected_table
  }
  # Helper to build expected table for scenario (neutral)
  create_expected_table_2x4_neutral <- function(obs_table, tiles, strict = FALSE) {
    obs_sub <- obs_table[, tiles, drop = FALSE]
    n_classes <- rowSums(obs_sub)
    rel_classes <- n_classes / sum(n_classes)
    expected_table <- matrix(0, nrow = 2, ncol = 4, dimnames = list(rownames(obs_table), tiles))
    if (strict) {
      expected_table[1,] <- 1/4 * n_classes[1]
      expected_table[2,] <- 1/4 * n_classes[2]
    } else {
      tile_totals <- colSums(obs_sub)
      expected_table[1,] <- round(tile_totals * rel_classes[1])
      expected_table[2,] <- round(tile_totals * rel_classes[2])
    }
    expected_table
  }
  # Compute dominance for both scenarios
  dominance_default <- compute_dominance(observed_table, TILES_DEFAULT, IDEAL_PATTERN)
  dominance_rotated <- compute_dominance(observed_table, TILES_ROTATED, IDEAL_PATTERN)
  if (is.na(dominance_default) && is.na(dominance_rotated)) {
    chosen_tiles <- TILES_DEFAULT
    scenario <- "default"
  } else if (is.na(dominance_rotated) || dominance_default >= dominance_rotated) {
    chosen_tiles <- TILES_DEFAULT
    scenario <- "default"
  } else {
    chosen_tiles <- TILES_ROTATED
    scenario <- "rotated"
  }
  # Prepare observed table for the dominant scenario
  observed_sub <- observed_table[, chosen_tiles, drop = FALSE]
  tile_totals <- colSums(observed_sub)
  # Automatically choose strict or not
  if (all(tile_totals > 0)) {
    strict_mode <- FALSE
    strict_note <- NULL
  } else {
    strict_mode <- TRUE
    strict_note <- "Some tiles empty: using strict (ideal) expected table."
  }
  expected_table <- create_expected_table_2x4(
    observed_table, chosen_tiles, IDEAL_PATTERN, strict = strict_mode
  )
  diff_result <- check_chisquared_difference(
    observed_sub, expected_table,
    abs_diff_threshold = abs_diff_threshold,
    p_threshold = p_threshold
  )
  # Neutral expected table
  expected_table_neutral <- create_expected_table_2x4_neutral(
    observed_table, chosen_tiles, strict = strict_mode
  )
  diff_result_neutral <- check_chisquared_difference(
    observed_sub, expected_table_neutral,
    abs_diff_threshold = abs_diff_threshold,
    p_threshold = p_threshold
  )
  # --- 2/9 CORNER/ROTATED TILE CONDITION ---
  tiles_population_thirds <- assign_3x3_tiles(data, var1_name, var2_name, split_method = "quantile")
  observed_table_thirds <- table(data[[class_name]], tiles_population_thirds)
  percent_in_corners <- function(tab, tiles) {
    apply(tab, 1, function(row) {
      sum(row[tiles]) / sum(row)
    })
  }
  corner_percentages_corners <-
    percent_in_corners(observed_table, c("Upper Left", "Lower Left", "Upper Right", "Lower Right"))
  corner_percentages_mids <-
    percent_in_corners(observed_table, c("Mid Left", "Mid Right", "Upper Mid", "Lower Mid"))
  corner_condition1_met <- FALSE
  corner_condition2_met <- FALSE
  if (scenario == "default") {
    if (all(corner_percentages_corners >= 2 / 9)) corner_condition1_met <- TRUE
    if (all(corner_percentages_mids < 1 / 9)) corner_condition2_met <- TRUE
    corner_message <- paste0(
      "Proportion of class 1 in target tiles: ", signif(corner_percentages_corners[1], 3), "; ",
      "class 2: ", signif(corner_percentages_corners[2], 3), "."
    )
  } else {
    if (all(corner_percentages_mids >= 2 / 9)) corner_condition1_met <- TRUE
    if (all(corner_percentages_corners < 1 / 9)) corner_condition2_met <- TRUE
    corner_message <- paste0(
      "Proportion of class 1 in target tiles: ", signif(corner_percentages_mids[1], 3), "; ",
      "class 2: ", signif(corner_percentages_mids[2], 3), "."
    )
  }
  # Compose interpretation
  interpretation <- paste0(
    strict_note, " ",
    "Standard p = ", signif(diff_result$p_value, 3),
    ", Monte Carlo p = ", signif(diff_result$p_value_mc, 3), ". ",
    diff_result$message, "\n",
    corner_message, "\n"
  )
  xor_shape_detected <- (!(diff_result$differs | diff_result$differs_mc)) &&
    diff_result_neutral$differs
  list(
    observed_table = observed_sub,
    expected_table = expected_table,
    chi_squared_test = list(
      statistic = diff_result$statistic,
      df = diff_result$df,
      p_value = diff_result$p_value,
      statistic_mc = diff_result$statistic_mc,
      p_value_mc = diff_result$p_value_mc
    ),
    abs_difference_sum = diff_result$abs_difference_sum,
    interpretation = interpretation,
    xor_shape_detected = xor_shape_detected,
    xor_pattern = scenario,
    dominance_default = dominance_default,
    dominance_rotated = dominance_rotated,
    strict_mode = strict_mode
  )
}


#' Compute tile pattern analysis for all pairwise datasets (modularized)
#'
#' @param orig_pair_list List. Output from create_pairwise_datasets.
#' @param class_col Character. Name of the class column (default: "class").
#' @param p_threshold Numeric. Significance threshold for p-values (default: 0.05).
#' @param abs_diff_threshold Numeric. Absolute difference threshold (default: 20).
#' @param split_method Character. Method for splitting tiles ("quantile" or "range") (default: "quantile").
#' @return The input list, but with results$tile_pattern filled for each pair.
#' @importFrom stats quantile chisq.test wilcox.test
#' @export
compute_tile_patterns_for_pairs <- function(
    orig_pair_list,
    class_col = "class",
    p_threshold = 0.05,
    abs_diff_threshold = 20,
    split_method = "quantile"
) {
  for (nm in names(orig_pair_list)) {
    pair <- orig_pair_list[[nm]]
    df <- pair$data
    var_names <- setdiff(names(df), class_col)
    if (length(var_names) != 2) stop("Each pairwise data frame must have exactly 2 variables + class column.")

    tile_pattern_result <- analyze_tile_dominance(
      data = df,
      var1_name = var_names[1],
      var2_name = var_names[2],
      class_name = class_col,
      p_threshold = p_threshold,
      abs_diff_threshold = abs_diff_threshold,
      split_method = split_method  # Critical addition
    )

    orig_pair_list[[nm]]$results$tile_pattern <- tile_pattern_result
  }
  orig_pair_list
}




# Suppose you have already created orig_pair_list with create_pairwise_datasets()
# pairwise_dfs <- compute_tile_patterns_for_pairs(pairwise_dfs)

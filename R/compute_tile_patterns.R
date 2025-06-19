#' Compute tile pattern analysis for all pairwise datasets
#'
#' This function iterates over a list of pairwise datasets (as produced by create_pairwise_datasets)
#' and fills the results$tile_pattern entry for each pair using the provided tile pattern analysis.
#'
#' Assign 3x3 tiles to data points based on two variables
#' @keywords internal
assign_3x3_tiles <- function(data, var1_name, var2_name, split_method = "quantile") {
  # Input validation
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("'data' must be a non-empty data frame", call. = FALSE)
  }
  if (!all(c(var1_name, var2_name) %in% names(data))) {
    stop("Variables '", var1_name, "' and '", var2_name, "' must exist in data", call. = FALSE)
  }
  if (!split_method %in% c("quantile", "range")) {
    stop("'split_method' must be either 'quantile' or 'range'", call. = FALSE)
  }

  # Get variable vectors
  var1_values <- data[[var1_name]]
  var2_values <- data[[var2_name]]

  # Check for sufficient variation
  if (length(unique(var1_values[!is.na(var1_values)])) < 2) {
    warning("Variable '", var1_name, "' has insufficient variation for tile assignment", call. = FALSE)
  }
  if (length(unique(var2_values[!is.na(var2_values)])) < 2) {
    warning("Variable '", var2_name, "' has insufficient variation for tile assignment", call. = FALSE)
  }

  # Determine split points
  if (split_method == "quantile") {
    x_limits <- compute_quantile_splits(var1_values)
    y_limits <- compute_quantile_splits(var2_values)
  } else {
    x_limits <- split_range_into_thirds(var1_values)
    y_limits <- split_range_into_thirds(var2_values)
  }

  # Internal function for flexible group assignment
  assign_flex_groups <- function(values, limits, labels) {
    breaks <- c(-Inf, limits, Inf)
    groups <- cut(values, breaks = breaks, labels = labels, include.lowest = TRUE)
    if (!any(values > limits[1] & values < limits[2]) &&
        any(values == limits[2]) &&
        !any(values > limits[2])) {
      groups[values == limits[2]] <- labels[length(labels)]
      groups <- factor(groups, levels = labels)
    }
    return(groups)
  }

  # Create group assignments
  x_group <- assign_flex_groups(var1_values, x_limits, c("Left", "Mid", "Right"))
  y_group <- assign_flex_groups(var2_values, y_limits, c("Lower", "Mid", "Upper"))

  # Map to tile names
  tile_names <- mapply(map_coordinates_to_tile, x_group, y_group, USE.NAMES = FALSE)

  # Create factor with proper ordering
  tile_levels <- c(
    "Upper Left", "Upper Mid", "Upper Right",
    "Mid Left", "Center", "Mid Right",
    "Lower Left", "Lower Mid", "Lower Right"
  )

  factor(tile_names, levels = tile_levels)
}


#' Compute robust quantile splits with fallback to range splits
#' @keywords internal
compute_quantile_splits <- function(x) {
  # Try quantile method first
  quantile_splits <- quantile(x, probs = c(1/3, 2/3), na.rm = TRUE)

  # Check if quantiles provide sufficient separation
  if (length(unique(quantile_splits)) < 2) {
    # Fall back to range-based splitting
    return(split_range_into_thirds(x))
  }

  return(quantile_splits)
}

#' Split a numeric vector into three equal-width bins
#' @keywords internal
split_range_into_thirds <- function(x) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)

  if (is.infinite(min_val) || is.infinite(max_val) || min_val == max_val) {
    # Handle edge cases
    return(c(min_val, min_val))
  }

  range_size <- max_val - min_val
  break1 <- min_val + range_size / 3
  break2 <- min_val + 2 * range_size / 3

  c(break1, break2)
}

#' Map x,y coordinates to tile names
#' @keywords internal
map_coordinates_to_tile <- function(x, y) {
  if (is.na(x) || is.na(y)) return(NA)

  coordinate_key <- paste(x, y)

  tile_mapping <- c(
    "Left Upper" = "Upper Left",
    "Right Upper" = "Upper Right",
    "Left Lower" = "Lower Left",
    "Right Lower" = "Lower Right",
    "Mid Upper" = "Upper Mid",
    "Left Mid" = "Mid Left",
    "Right Mid" = "Mid Right",
    "Mid Lower" = "Lower Mid",
    "Mid Mid" = "Center"
  )

  result <- tile_mapping[coordinate_key]
  if (is.na(result)) {
    warning("Unknown coordinate combination: ", coordinate_key, call. = FALSE)
  }

  return(result)
}

#' Chi-squared test and difference check for observed vs expected counts
#' @keywords internal
check_chisquared_difference <- function(
  observed_counts, expected_counts,
  abs_diff_threshold = 20,
  filter_expected = 0,
  p_threshold = 0.05,
  monte_carlo_B = 2000
) {
  # Input validation and conversion
  observed_counts <- as.matrix(observed_counts)
  expected_counts <- as.matrix(expected_counts)

  if (!identical(dim(observed_counts), dim(expected_counts))) {
    stop("Observed and expected count matrices must have the same dimensions", call. = FALSE)
  }

  # Handle zero expected counts
  expected_counts[expected_counts == 0] <- 0.5

  # Convert to vectors for analysis
  observed_vector <- as.vector(t(observed_counts))
  expected_vector <- as.vector(t(expected_counts))

  # Apply filtering if specified
  if (filter_expected > 0) {
    mask <- expected_vector > filter_expected
    observed_vector_filt <- observed_vector[mask]
    expected_vector_filt <- expected_vector[mask]
  } else {
    observed_vector_filt <- observed_vector
    expected_vector_filt <- expected_vector
  }

  # Calculate absolute difference
  abs_difference_sum <- sum(abs(observed_vector - expected_vector))

  # Initialize result variables
  p_value <- NA
  statistic <- NA
  df <- NA
  p_value_mc <- NA
  statistic_mc <- NA

  # Perform chi-square tests if we have sufficient data
  if (length(observed_vector_filt) > 1 && sum(expected_vector_filt) > 0) {

    # Standard chi-square test
    tryCatch({
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
    }, error = function(e) {
      warning("Standard chi-square test failed: ", e$message, call. = FALSE)
    })

    # Monte Carlo chi-square test
    if (monte_carlo_B > 0) {
      tryCatch({
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
      }, error = function(e) {
        warning("Monte Carlo chi-square test failed: ", e$message, call. = FALSE)
      })
    }
  }

  # Determine significance
  differs <- (!is.na(p_value) && p_value < p_threshold) ||
    abs_difference_sum > abs_diff_threshold

  differs_mc <- (!is.na(p_value_mc) && p_value_mc < p_threshold) ||
    abs_difference_sum > abs_diff_threshold

  # Generate message
  message <- if (differs || differs_mc) {
    "Observed counts differ from expectations (statistically or practically)."
  } else {
    "No significant or practical difference detected."
  }

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
    message = message
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
  # Input validation
  validate_tile_analysis_inputs(data, var1_name, var2_name, class_name, p_threshold, abs_diff_threshold)

  # Define tile sets and ideal pattern
  TILES_DEFAULT <- c("Lower Left", "Lower Right", "Upper Left", "Upper Right")
  TILES_ROTATED <- c("Mid Left", "Upper Mid", "Lower Mid", "Mid Right")
  IDEAL_PATTERN <- c(1, 2, 2, 1)

  # Assign tiles and create observed table
  tiles_population <- assign_3x3_tiles(data, var1_name, var2_name, split_method)
  observed_table <- table(data[[class_name]], tiles_population)

  # Ensure we have the expected number of classes
  if (nrow(observed_table) != 2) {
    stop("Analysis requires exactly 2 classes, found ", nrow(observed_table), call. = FALSE)
  }

  # Compute dominance for both scenarios
  dominance_results <- compute_scenario_dominance(observed_table, TILES_DEFAULT, TILES_ROTATED, IDEAL_PATTERN)

  # Select best scenario
  scenario_selection <- select_best_scenario(dominance_results)
  chosen_tiles <- scenario_selection$tiles
  scenario <- scenario_selection$scenario

  # Analyze chosen scenario
  analysis_result <- analyze_chosen_scenario(
    observed_table, chosen_tiles, IDEAL_PATTERN, scenario,
    p_threshold, abs_diff_threshold
  )

  # Add scenario information
  analysis_result$xor_pattern <- scenario
  analysis_result$dominance_default <- dominance_results$default
  analysis_result$dominance_rotated <- dominance_results$rotated

  return(analysis_result)
}

#' Validate inputs for tile analysis
#' @keywords internal
validate_tile_analysis_inputs <- function(data, var1_name, var2_name, class_name, p_threshold, abs_diff_threshold) {
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("'data' must be a non-empty data frame", call. = FALSE)
  }

  required_vars <- c(var1_name, var2_name, class_name)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing variables in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  if (!is.numeric(p_threshold) || p_threshold <= 0 || p_threshold >= 1) {
    stop("'p_threshold' must be between 0 and 1", call. = FALSE)
  }

  if (!is.numeric(abs_diff_threshold) || abs_diff_threshold < 0) {
    stop("'abs_diff_threshold' must be non-negative", call. = FALSE)
  }
}

#' Compute dominance scores for different scenarios
#' @keywords internal
compute_scenario_dominance <- function(observed_table, tiles_default, tiles_rotated, ideal_pattern) {
  compute_dominance <- function(obs_table, tiles, ideal_pattern) {
    if (!all(tiles %in% colnames(obs_table))) {
      return(NA)
    }

    obs_sub <- obs_table[, tiles, drop = FALSE]
    if (sum(obs_sub) == 0) {
      return(NA)
    }

    ideal_cells <- mapply(function(class, col) {
      if (class <= nrow(obs_sub) && col <= ncol(obs_sub)) {
        obs_sub[class, col]
      } else {
        0
      }
    }, ideal_pattern, seq_along(ideal_pattern))

    sum(ideal_cells) / sum(obs_sub)
  }

  list(
    default = compute_dominance(observed_table, tiles_default, ideal_pattern),
    rotated = compute_dominance(observed_table, tiles_rotated, ideal_pattern)
  )
}

#' Select the best scenario based on dominance scores
#' @keywords internal
select_best_scenario <- function(dominance_results) {
  TILES_DEFAULT <- c("Lower Left", "Lower Right", "Upper Left", "Upper Right")
  TILES_ROTATED <- c("Mid Left", "Upper Mid", "Lower Mid", "Mid Right")

  dom_default <- dominance_results$default
  dom_rotated <- dominance_results$rotated

  if (is.na(dom_default) && is.na(dom_rotated)) {
    # Both scenarios invalid, default to default
    return(list(tiles = TILES_DEFAULT, scenario = "default"))
  } else if (is.na(dom_rotated) || (!is.na(dom_default) && dom_default >= dom_rotated)) {
    return(list(tiles = TILES_DEFAULT, scenario = "default"))
  } else {
    return(list(tiles = TILES_ROTATED, scenario = "rotated"))
  }
}

#' Analyze the chosen scenario in detail
#' @keywords internal
analyze_chosen_scenario <- function(observed_table, chosen_tiles, ideal_pattern, scenario, p_threshold, abs_diff_threshold) {
  # Extract observed data for chosen tiles
  observed_sub <- observed_table[, chosen_tiles, drop = FALSE]
  tile_totals <- colSums(observed_sub)

  # Determine if we need strict mode
  strict_mode <- any(tile_totals == 0)
  strict_note <- if (strict_mode) "Some tiles empty: using strict (ideal) expected table." else ""

  # Create expected tables
  expected_table <- create_expected_table_2x4(observed_table, chosen_tiles, ideal_pattern, strict_mode)
  expected_table_neutral <- create_expected_table_2x4_neutral(observed_table, chosen_tiles, strict_mode)

  # Perform statistical tests
  diff_result <- check_chisquared_difference(
    observed_sub, expected_table,
    abs_diff_threshold = abs_diff_threshold,
    p_threshold = p_threshold
  )

  diff_result_neutral <- check_chisquared_difference(
    observed_sub, expected_table_neutral,
    abs_diff_threshold = abs_diff_threshold,
    p_threshold = p_threshold
  )

  # Compute corner conditions
  corner_analysis <- analyze_corner_conditions(observed_table, scenario)

  # Create interpretation
  interpretation <- create_interpretation_message(
    strict_note, diff_result, corner_analysis$message
  )

  # Determine XOR detection
  xor_shape_detected <- (!(diff_result$differs || diff_result$differs_mc)) && diff_result_neutral$differs

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
    strict_mode = strict_mode
  )
}

#' Create expected table for XOR pattern
#' @keywords internal
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

#' Create neutral expected table
#' @keywords internal
create_expected_table_2x4_neutral <- function(obs_table, tiles, strict = FALSE) {
  obs_sub <- obs_table[, tiles, drop = FALSE]
  n_classes <- rowSums(obs_sub)

  if (sum(n_classes) == 0) {
    return(matrix(0, nrow = 2, ncol = 4, dimnames = list(rownames(obs_table), tiles)))
  }

  rel_classes <- n_classes / sum(n_classes)
  expected_table <- matrix(0, nrow = 2, ncol = 4, dimnames = list(rownames(obs_table), tiles))

  if (strict) {
    expected_table[1, ] <- n_classes[1] / 4
    expected_table[2, ] <- n_classes[2] / 4
  } else {
    tile_totals <- colSums(obs_sub)
    expected_table[1, ] <- tile_totals * rel_classes[1]
    expected_table[2, ] <- tile_totals * rel_classes[2]
  }

  expected_table
}

#' Analyze corner/edge conditions for different scenarios
#' @keywords internal
analyze_corner_conditions <- function(observed_table, scenario) {
  percent_in_tiles <- function(tab, tiles) {
    if (!all(tiles %in% colnames(tab))) {
      return(rep(NA, nrow(tab)))
    }
    apply(tab, 1, function(row) {
      sum(row[tiles]) / sum(row)
    })
  }

  corner_percentages <- percent_in_tiles(observed_table, c("Upper Left", "Lower Left", "Upper Right", "Lower Right"))
  mid_percentages <- percent_in_tiles(observed_table, c("Mid Left", "Mid Right", "Upper Mid", "Lower Mid"))

  if (scenario == "default") {
    target_percentages <- corner_percentages
    message <- sprintf(
      "Proportion of class 1 in target tiles: %s; class 2: %s.",
      signif(target_percentages[1], 3),
      signif(target_percentages[2], 3)
    )
  } else {
    target_percentages <- mid_percentages
    message <- sprintf(
      "Proportion of class 1 in target tiles: %s; class 2: %s.",
      signif(target_percentages[1], 3),
      signif(target_percentages[2], 3)
    )
  }

  list(
    target_percentages = target_percentages,
    message = message
  )
}

#' Create interpretation message
#' @keywords internal
create_interpretation_message <- function(strict_note, diff_result, corner_message) {
  paste0(
    if (nchar(strict_note) > 0) paste0(strict_note, " ") else "",
    "Standard p = ", signif(diff_result$p_value, 3),
    ", Monte Carlo p = ", signif(diff_result$p_value_mc, 3), ". ",
    diff_result$message, "\n",
    corner_message, "\n"
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
#' @importFrom stats quantile chisq.test
#' @export
compute_tile_patterns_for_pairs <- function(
  orig_pair_list,
  class_col = "class",
  p_threshold = 0.05,
  abs_diff_threshold = 20,
  split_method = "quantile"
) {
  # Input validation
  if (!is.list(orig_pair_list) || length(orig_pair_list) == 0) {
    stop("'orig_pair_list' must be a non-empty list", call. = FALSE)
  }

  if (!is.character(class_col) || length(class_col) != 1) {
    stop("'class_col' must be a single character string", call. = FALSE)
  }

  if (!is.numeric(p_threshold) || p_threshold <= 0 || p_threshold >= 1) {
    stop("'p_threshold' must be between 0 and 1", call. = FALSE)
  }

  if (!is.numeric(abs_diff_threshold) || abs_diff_threshold < 0) {
    stop("'abs_diff_threshold' must be non-negative", call. = FALSE)
  }

  if (!split_method %in% c("quantile", "range")) {
    stop("'split_method' must be either 'quantile' or 'range'", call. = FALSE)
  }

  # Process each pair
  for (nm in names(orig_pair_list)) {
    tryCatch({
      pair <- orig_pair_list[[nm]]

      # Validate pair structure
      if (is.null(pair$data) || !is.data.frame(pair$data)) {
        warning("Pair '", nm, "' missing valid data, skipping", call. = FALSE)
        next
      }

      df <- pair$data
      var_names <- setdiff(names(df), class_col)

      if (length(var_names) != 2) {
        warning("Pair '", nm, "' must contain exactly 2 variables + class column, skipping", call. = FALSE)
        next
      }

      if (!class_col %in% names(df)) {
        warning("Class column '", class_col, "' not found in pair '", nm, "', skipping", call. = FALSE)
        next
      }

      # Perform tile pattern analysis
      tile_pattern_result <- analyze_tile_dominance(
        data = df,
        var1_name = var_names[1],
        var2_name = var_names[2],
        class_name = class_col,
        p_threshold = p_threshold,
        abs_diff_threshold = abs_diff_threshold,
        split_method = split_method
      )

      # Store results
      if (is.null(orig_pair_list[[nm]]$results)) {
        orig_pair_list[[nm]]$results <- list()
      }
      orig_pair_list[[nm]]$results$tile_pattern <- tile_pattern_result

    }, error = function(e) {
      warning("Error processing pair '", nm, "': ", e$message, call. = FALSE)
    })
  }

  orig_pair_list
}

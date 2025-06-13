#' Compute group-wise Wilcoxon significance for tile groups
#'
#' @param orig_pair_list List. Output from compute_tile_patterns_for_pairs.
#' @param class_col Character. Name of the class column (default: "class").
#' @param p_threshold Numeric. Significance threshold (default: 0.05).
#' @param split_method Character. Tile assignment method ("quantile" or "range") (default: "quantile").
#' @return List with results$tile_wilcox_significance added for each pair.
#' @export
compute_tile_wilcox_significance_for_pairs <- function(
  orig_pair_list,
  class_col = "class",
  p_threshold = 0.05,
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
    stop("'p_threshold' must be a numeric value between 0 and 1", call. = FALSE)
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

      tile_pattern <- pair$results$tile_pattern
      if (is.null(tile_pattern)) {
        warning("Pair '", nm, "' missing tile_pattern results, skipping", call. = FALSE)
        next
      }

      # Compute Wilcoxon significance
      tile_wilcox_result <- check_tile_groupwise_wilcox_significance(
        data = df,
        var1_name = var_names[1],
        var2_name = var_names[2],
        class_name = class_col,
        chosen_tiles = colnames(tile_pattern$observed_table),
        scenario = tile_pattern$xor_pattern,
        p_threshold = p_threshold,
        split_method = split_method
      )

      orig_pair_list[[nm]]$results$tile_wilcox_significance <- tile_wilcox_result

    }, error = function(e) {
      warning("Error processing pair '", nm, "': ", e$message, call. = FALSE)
    })
  }

  orig_pair_list
}

#' Check group-wise Wilcoxon significance for tile groups (for a single pair)
#' Uses scenario and chosen tiles as determined by tile_pattern.
#' @keywords internal
check_tile_groupwise_wilcox_significance <- function(
  data, var1_name, var2_name, class_name, chosen_tiles, scenario,
  p_threshold = 0.05,
  split_method = "quantile"
) {
  # Input validation
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("'data' must be a non-empty data frame", call. = FALSE)
  }

  required_vars <- c(var1_name, var2_name, class_name)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing variables in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  # Check for sufficient data
  if (length(unique(data[[class_name]])) < 2) {
    return(create_empty_wilcox_result("Less than 2 classes present", scenario, chosen_tiles))
  }

  # Assign tiles to all data points
  tiles_population <- assign_tiles_to_data(data, var1_name, var2_name, split_method)

  # Filter to chosen tiles only
  data_tiles <- filter_to_chosen_tiles(data, tiles_population, chosen_tiles, class_name)

  if (nrow(data_tiles) == 0) {
    return(create_empty_wilcox_result("No data in chosen tiles", scenario, chosen_tiles))
  }

  # Apply rotation transformation if needed
  processed_data <- apply_scenario_transformation(data_tiles, var1_name, var2_name, scenario)

  # Run Wilcoxon tests for all four groups
  test_results <- run_all_group_tests(processed_data, class_name)

  # Generate summary and results
  create_wilcox_result(test_results, p_threshold, scenario, chosen_tiles)
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Assign tiles to data points using 3x3 grid
assign_tiles_to_data <- function(data, var1_name, var2_name, split_method) {
  if (split_method == "quantile") {
    x_limits <- quantile(data[[var1_name]], probs = c(1/3, 2/3), na.rm = TRUE)
    y_limits <- quantile(data[[var2_name]], probs = c(1/3, 2/3), na.rm = TRUE)
  } else {
    x_limits <- split_range_into_thirds(data[[var1_name]])
    y_limits <- split_range_into_thirds(data[[var2_name]])
  }

  x_group <- cut(data[[var1_name]],
                 breaks = c(-Inf, x_limits, Inf),
                 labels = c("Left", "Mid", "Right"),
                 include.lowest = TRUE)

  y_group <- cut(data[[var2_name]],
                 breaks = c(-Inf, y_limits, Inf),
                 labels = c("Lower", "Mid", "Upper"),
                 include.lowest = TRUE)

  # Create tile names from x and y groups
  tile_names <- mapply(combine_to_tile_name, x_group, y_group, USE.NAMES = FALSE)

  # Define proper factor levels for tiles
  tile_levels <- c(
    "Upper Left", "Upper Mid", "Upper Right",
    "Mid Left", "Center", "Mid Right",
    "Lower Left", "Lower Mid", "Lower Right"
  )

  factor(tile_names, levels = tile_levels)
}

#' Split numeric range into three equal parts
split_range_into_thirds <- function(x) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  range_size <- max_val - min_val

  c(min_val + range_size/3, min_val + 2*range_size/3)
}

#' Combine x and y group labels into tile names
combine_to_tile_name <- function(x, y) {
  tile_map <- c(
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

  key <- paste(x, y)
  tile_map[key]
}

#' Filter data to chosen tiles and prepare for analysis
filter_to_chosen_tiles <- function(data, tiles_population, chosen_tiles, class_name) {
  # Keep only observations in chosen tiles
  in_chosen_tiles <- tiles_population %in% chosen_tiles
  data_subset <- data[in_chosen_tiles, , drop = FALSE]

  # Add tile and class factor columns
  data_subset$tile <- factor(tiles_population[in_chosen_tiles], levels = chosen_tiles)
  data_subset$class <- as.factor(data_subset[[class_name]])

  data_subset
}

#' Apply rotation transformation for rotated XOR scenarios
apply_scenario_transformation <- function(data_tiles, var1_name, var2_name, scenario) {
  if (scenario == "default") {
    # No transformation needed
    processed_data <- data_tiles
    processed_data$v1 <- var1_name
    processed_data$v2 <- var2_name

  } else if (scenario == "rotated") {
    # Apply rotation transformation
    processed_data <- data_tiles
    processed_data$rot_plus <- data_tiles[[var1_name]] + data_tiles[[var2_name]]
    processed_data$rot_minus <- data_tiles[[var1_name]] - data_tiles[[var2_name]]

    # Map original tiles to rotated coordinate system
    tile_rotation_map <- c(
      "Mid Left" = "Upper Left",
      "Upper Mid" = "Upper Right",
      "Mid Right" = "Lower Right",
      "Lower Mid" = "Lower Left"
    )

    processed_data$tile_rot <- tile_rotation_map[as.character(processed_data$tile)]
    processed_data$tile_rot <- factor(processed_data$tile_rot,
                                      levels = c("Lower Left", "Lower Right", "Upper Left", "Upper Right"))

    # Update tile assignment and variable names
    processed_data$tile <- processed_data$tile_rot
    processed_data$v1 <- "rot_plus"
    processed_data$v2 <- "rot_minus"

  } else {
    stop("Unknown scenario: ", scenario, call. = FALSE)
  }

  processed_data
}

#' Run Wilcoxon tests for all four tile groups
run_all_group_tests <- function(processed_data, class_name) {
  # Define tile groupings
  group_definitions <- list(
    left = c("Lower Left", "Upper Left"),
    right = c("Lower Right", "Upper Right"),
    upper = c("Upper Left", "Upper Right"),
    lower = c("Lower Left", "Lower Right")
  )

  # Run tests for each group
  test_results <- lapply(group_definitions, function(tiles) {
    run_combined_wilcox_test(processed_data, tiles,
                             processed_data$v1[1], processed_data$v2[1], class_name)
  })

  names(test_results) <- names(group_definitions)
  test_results
}

#' Run combined Wilcoxon test for a specific tile group
run_combined_wilcox_test <- function(df, tiles, var1, var2, class_col) {
  # Select observations in specified tiles
  in_tiles <- df$tile %in% tiles

  if (sum(in_tiles) == 0) {
    return(create_empty_test_result())
  }

  # Combine values from both variables
  values <- c(df[in_tiles, var1], df[in_tiles, var2])
  classes <- c(df[in_tiles, class_col], df[in_tiles, class_col])

  # Check if we have enough data for testing
  if (length(unique(classes)) < 2 || length(values) <= 1) {
    return(create_empty_test_result())
  }

  # Perform Wilcoxon test
  test_result <- tryCatch({
    wilcox.test(values ~ classes)
  }, error = function(e) {
    NULL
  })

  if (is.null(test_result)) {
    return(create_empty_test_result())
  }

  # Calculate summary statistics
  means <- tapply(values, classes, mean, na.rm = TRUE)
  medians <- tapply(values, classes, median, na.rm = TRUE)
  direction <- if (length(means) == 2) names(means)[which.max(means)] else NA

  list(
    p_value = test_result$p.value,
    statistic = test_result$statistic,
    means = means,
    medians = medians,
    direction = direction,
    valid = TRUE
  )
}

#' Create empty test result for cases with insufficient data
create_empty_test_result <- function() {
  list(
    p_value = NA,
    statistic = NA,
    means = NA,
    medians = NA,
    direction = NA,
    valid = FALSE
  )
}

#' Generate summary string for test results
generate_test_summary <- function(test, group_name) {
  if (test$valid && !is.na(test$p_value)) {
    sprintf(
      "%s: class 1 mean = %s, class 2 mean = %s, p = %s, higher = class %s",
      group_name,
      signif(test$means[1], 3),
      signif(test$means[2], 3),
      signif(test$p_value, 3),
      test$direction
    )
  } else {
    sprintf("%s: insufficient data or only one class present", group_name)
  }
}

#' Create final Wilcoxon test result
create_wilcox_result <- function(test_results, p_threshold, scenario, chosen_tiles) {
  # Generate summaries
  summaries <- mapply(generate_test_summary, test_results,
                      c("Left", "Right", "Upper", "Lower"),
                      SIMPLIFY = FALSE)

  # Extract p-values
  p_values <- sapply(test_results, function(x) x$p_value)

  # Check if all tests pass significance threshold
  valid_p_values <- p_values[!is.na(p_values)]
  pass_groups <- length(valid_p_values) == 4 && all(valid_p_values < p_threshold)

  # Generate overall message
  tile_group_message <- if (pass_groups) {
    "All four group-wise combined tests are significant: XOR pattern is statistically supported."
  } else {
    "Not all four group-wise combined tests are significant or present: XOR pattern is not statistically supported."
  }

  list(
    pass_groups = pass_groups,
    summary_left = summaries[[1]],
    summary_right = summaries[[2]],
    summary_upper = summaries[[3]],
    summary_lower = summaries[[4]],
    tile_group_message = tile_group_message,
    pvals_groups = p_values,
    scenario = scenario,
    chosen_tiles = chosen_tiles
  )
}

#' Create empty result for cases where analysis cannot proceed
create_empty_wilcox_result <- function(reason, scenario, chosen_tiles) {
  list(
    pass_groups = FALSE,
    summary_left = paste("Left:", reason),
    summary_right = paste("Right:", reason),
    summary_upper = paste("Upper:", reason),
    summary_lower = paste("Lower:", reason),
    tile_group_message = paste("Analysis failed:", reason),
    pvals_groups = rep(NA, 4),
    scenario = scenario,
    chosen_tiles = chosen_tiles
  )
}
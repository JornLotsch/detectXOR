
#' @title Generate Spaghetti Plots for XOR-Detected Pairs
#' @description Creates connected line plots for variable pairs flagged with XOR patterns.
#'   Each line connects the same observation across the two variables, showing
#'   how values change between variables for each class.
#' 
#' @param results Either a data frame from \code{detect_xor()$results_df} or 
#'   the full list returned by \code{detect_xor()}.
#' @param data Original dataset containing variables and classes used in XOR detection.
#' @param class_col Character string specifying the class column name.
#' @param scale_data Logical indicating whether to scale variables to unit variance. 
#'   Should match the scaling used in \code{detect_xor()}. Default: TRUE.
#' 
#' @return A ggplot2 object. Returns an empty plot with message if no XOR patterns detected.
#' 
#' @details The function creates "spaghetti plots" where each line represents one observation,
#'   connecting its values across the two variables in each detected XOR pair. 
#'   Different colors represent different classes, helping visualize the XOR pattern.
#'   
#'   Facet headers are colored based on chi-square p-values, and rotation is applied
#'   automatically for pairs detected as "rotated" XOR patterns.
#' 
#' @examples
#' \dontrun{
#' # Basic usage
#' data(iris)
#' iris$class <- iris$Species
#' results <- detect_xor(iris, class_col = "class")
#' plot <- generate_spaghetti_plot_from_results(results, iris, "class")
#' print(plot)
#' 
#' # Using just the results dataframe
#' plot <- generate_spaghetti_plot_from_results(results$results_df, iris, "class")
#' }
#' 
#' @seealso \code{\link{detect_xor}}, \code{\link{generate_xy_plot_from_results}}
#' @export
generate_spaghetti_plot_from_results <- function(results, data, class_col, scale_data = TRUE) {

  # Extract results dataframe and validate inputs
  results_info <- extract_and_validate_results(results, data, class_col)
  results_df <- results_info$results_df
  results_name <- results_info$results_name

  # Get XOR-detected variable pairs
  variable_pairs <- get_xor_variable_pairs(results_df)

  if (nrow(variable_pairs) == 0) {
    return(create_empty_plot("No XOR patterns detected\nTry adjusting detection parameters",
                             "Spaghetti Plot - XOR Pattern Detection Results"))
  }

  # Process data for all variable pairs
  long_data <- process_pairs_for_spaghetti(variable_pairs, data, results_df, class_col, scale_data)

  # Create facet annotations
  facet_colors <- create_facet_annotations(results_df, variable_pairs)

  # Build and return plot
  create_spaghetti_plot(long_data, facet_colors, results_name)
}

#' @title Generate XY Scatter Plots for XOR-Detected Pairs
#' @description Creates scatterplots with decision boundaries for detected XOR patterns.
#'   Shows the spatial distribution of classes in the variable space with
#'   quantile or range-based reference lines.
#' 
#' @param results Either a data frame from \code{detect_xor()$results_df} or 
#'   the full list returned by \code{detect_xor()}.
#' @param data Original dataset containing variables and classes used in XOR detection.
#' @param class_col Character string specifying the class column name.
#' @param scale_data Logical indicating whether to scale variables to unit variance. 
#'   Should match the scaling used in \code{detect_xor()}. Default: TRUE.
#' @param quantile_lines Numeric vector of quantiles for reference lines when 
#'   line_method = "quantile". Default: c(1/3, 2/3).
#' @param line_method Character string specifying boundary calculation method:
#'   "quantile" or "range". Default: "quantile".
#' 
#' @return A ggplot2 object. Returns an empty plot with message if no XOR patterns detected.
#' 
#' @details Creates scatter plots for each XOR-detected variable pair with reference lines
#'   that help visualize the decision boundaries. The plots show how different classes
#'   are distributed in the 2D variable space, with the XOR pattern visible as
#'   diagonal separation of classes.
#' 
#' @examples
#' \dontrun{
#' # Basic usage
#' data(iris)
#' iris$class <- iris$Species
#' results <- detect_xor(iris, class_col = "class")
#' plot <- generate_xy_plot_from_results(results, iris, "class")
#' print(plot)
#' 
#' # With custom quantile lines
#' plot <- generate_xy_plot_from_results(results, iris, "class", 
#'                                      quantile_lines = c(0.25, 0.75))
#' 
#' # Using range-based boundaries
#' plot <- generate_xy_plot_from_results(results, iris, "class", 
#'                                      line_method = "range")
#' }
#' 
#' @seealso \code{\link{detect_xor}}, \code{\link{generate_spaghetti_plot_from_results}}
#' @export
generate_xy_plot_from_results <- function(results, data, class_col, scale_data = TRUE,
                                          quantile_lines = c(1/3, 2/3), line_method = "quantile") {

  # Input validation
  if (!line_method %in% c("quantile", "range")) {
    stop("'line_method' must be either 'quantile' or 'range'", call. = FALSE)
  }

  if (!is.numeric(quantile_lines) || any(quantile_lines < 0) || any(quantile_lines > 1)) {
    stop("'quantile_lines' must be numeric values between 0 and 1", call. = FALSE)
  }

  # Extract results dataframe and validate inputs
  results_info <- extract_and_validate_results(results, data, class_col)
  results_df <- results_info$results_df
  results_name <- results_info$results_name

  # Get XOR-detected variable pairs
  variable_pairs <- get_xor_variable_pairs(results_df)

  if (nrow(variable_pairs) == 0) {
    return(create_empty_plot("No XOR patterns detected\nTry adjusting detection parameters",
                             "XY Scatter Plot - XOR Pattern Detection Results"))
  }

  # Process data for all variable pairs
  combined_data <- process_pairs_for_xy(variable_pairs, data, results_df, class_col, scale_data)

  # Calculate boundaries
  boundaries <- calculate_boundaries(combined_data, quantile_lines, line_method)

  # Create facet annotations
  facet_colors <- create_facet_annotations(results_df, variable_pairs)

  # Build and return plot
  create_xy_plot(combined_data, boundaries, facet_colors, results_name)
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Extract results dataframe and validate all inputs
extract_and_validate_results <- function(results, data, class_col) {
  # Handle both list and dataframe inputs
  if (is.list(results) && "results_df" %in% names(results)) {
    results_df <- results$results_df
    results_name <- deparse(substitute(results))
  } else if (is.data.frame(results)) {
    results_df <- results
    results_name <- deparse(substitute(results))
  } else {
    stop("'results' must be either a data frame or a list containing 'results_df'", call. = FALSE)
  }

  # Validate inputs
  validate_visualization_inputs(results_df, data, class_col)

  list(results_df = results_df, results_name = results_name)
}

#' Comprehensive input validation for visualization functions
validate_visualization_inputs <- function(results_df, data, class_col) {
  # Check required columns in results
  required_cols <- c("var1", "var2", "xor_shape_detected")
  missing_cols <- setdiff(required_cols, names(results_df))
  if (length(missing_cols) > 0) {
    stop("Results missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  # Check data frame
  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame", call. = FALSE)
  }

  if (nrow(data) == 0) {
    stop("'data' cannot be empty", call. = FALSE)
  }

  # Check class column
  if (!is.character(class_col) || length(class_col) != 1) {
    stop("'class_col' must be a single character string", call. = FALSE)
  }

  if (!class_col %in% names(data)) {
    stop("Column '", class_col, "' not found in data", call. = FALSE)
  }

  # Check if variables in results exist in data
  all_vars <- unique(c(results_df$var1, results_df$var2))
  missing_vars <- setdiff(all_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  # Check for numeric variables
  non_numeric <- all_vars[!sapply(data[all_vars], is.numeric)]
  if (length(non_numeric) > 0) {
    warning("Non-numeric variables detected: ", paste(non_numeric, collapse = ", "),
            ". They will be converted to numeric.", call. = FALSE)
  }
}

#' Get variable pairs that have XOR patterns detected
get_xor_variable_pairs <- function(results_df) {
  xor_pairs <- results_df[results_df$xor_shape_detected == TRUE, c("var1", "var2"), drop = FALSE]

  if (nrow(xor_pairs) == 0) {
    message("No XOR patterns detected in results")
  } else {
    message("Found ", nrow(xor_pairs), " XOR-detected variable pairs")
  }

  return(xor_pairs)
}

#' Create empty plot for cases with no XOR patterns
create_empty_plot <- function(message_text, title_text) {
  ggplot() +
    geom_text(aes(x = 0.5, y = 0.5, label = message_text),
              size = 5, hjust = 0.5, vjust = 0.5) +
    theme_void() +
    labs(title = title_text,
         subtitle = "No patterns found")
}

#' Apply XOR transformations (rotation and scaling) to pair data
apply_xor_transformations <- function(pair_data, rotate_condition, scale_data, class_col) {
  # Apply rotation if needed
  if (rotate_condition) {
    x_orig <- pair_data[, 2]
    y_orig <- pair_data[, 3]
    pair_data[, 2] <- x_orig + y_orig
    pair_data[, 3] <- x_orig - y_orig
  }

  # Apply scaling if requested
  if (scale_data) {
    cols_to_scale <- setdiff(names(pair_data), class_col)
    numeric_cols <- sapply(pair_data[cols_to_scale], is.numeric)

    if (any(numeric_cols)) {
      pair_data[cols_to_scale][numeric_cols] <- scale(pair_data[cols_to_scale][numeric_cols])
    }
  }

  return(pair_data)
}

#' Process variable pairs for spaghetti plot
process_pairs_for_spaghetti <- function(variable_pairs, data, results_df, class_col, scale_data) {
  message("Processing ", nrow(variable_pairs), " variable pairs for spaghetti plot...")

  # Process each pair
  processed_pairs <- lapply(seq_len(nrow(variable_pairs)), function(i) {
    pair <- variable_pairs[i, ]

    # Find corresponding row in results
    results_row_idx <- which(results_df$var1 == pair$var1 & results_df$var2 == pair$var2)

    # Check rotation condition
    rotate_condition <- ("xor_pattern" %in% names(results_df) &&
      length(results_row_idx) > 0 &&
      !is.na(results_df$xor_pattern[results_row_idx[1]]) &&
      results_df$xor_pattern[results_row_idx[1]] == "rotated")

    # Extract pair data
    pair_data <- data[, c(class_col, pair$var1, pair$var2), drop = FALSE]

    # Apply transformations
    pair_data <- apply_xor_transformations(pair_data, rotate_condition, scale_data, class_col)

    # Rename columns for consistency
    colnames(pair_data) <- c("class", "var1", "var2")

    # Add identifiers
    pair_data$var1_var2 <- paste(pair$var1, pair$var2, sep = "||")
    if (rotate_condition) {
      pair_data$var1_var2 <- paste(pair_data$var1_var2, "rotated", sep = "||")
    }
    pair_data$ID <- seq_len(nrow(pair_data))

    return(pair_data)
  })

  # Combine all pairs
  long_data <- do.call(rbind, processed_pairs)

  # Reshape for spaghetti plot
  long_data <- reshape2::melt(long_data, id.vars = c("ID", "class", "var1_var2"))

  return(long_data)
}

#' Process variable pairs for XY scatter plot
process_pairs_for_xy <- function(variable_pairs, data, results_df, class_col, scale_data) {
  message("Processing ", nrow(variable_pairs), " variable pairs for XY plot...")

  # Process each pair
  processed_pairs <- lapply(seq_len(nrow(variable_pairs)), function(i) {
    pair <- variable_pairs[i, ]

    # Find corresponding row in results
    results_row_idx <- which(results_df$var1 == pair$var1 & results_df$var2 == pair$var2)

    # Check rotation condition
    rotate_condition <- ("xor_pattern" %in% names(results_df) &&
      length(results_row_idx) > 0 &&
      !is.na(results_df$xor_pattern[results_row_idx[1]]) &&
      results_df$xor_pattern[results_row_idx[1]] == "rotated")

    # Extract pair data
    pair_data <- data[, c(class_col, pair$var1, pair$var2), drop = FALSE]

    # Apply transformations
    pair_data <- apply_xor_transformations(pair_data, rotate_condition, scale_data, class_col)

    # Rename columns for XY plot
    colnames(pair_data) <- c("class", "x", "y")

    # Add identifier
    pair_data$var1_var2 <- paste(pair$var1, pair$var2, sep = "||")
    if (rotate_condition) {
      pair_data$var1_var2 <- paste(pair_data$var1_var2, "rotated", sep = "||")
    }

    return(pair_data)
  })

  # Combine all pairs
  do.call(rbind, processed_pairs)
}

#' Calculate boundaries for XY plot
calculate_boundaries <- function(combined_data, quantile_lines, line_method) {
  if (line_method == "quantile") {
    boundaries <- combined_data %>%
      dplyr::group_by(var1_var2) %>%
      dplyr::summarise(
        x_q1 = quantile(x, quantile_lines[1], na.rm = TRUE),
        x_q2 = quantile(x, quantile_lines[2], na.rm = TRUE),
        y_q1 = quantile(y, quantile_lines[1], na.rm = TRUE),
        y_q2 = quantile(y, quantile_lines[2], na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    boundaries <- combined_data %>%
      dplyr::group_by(var1_var2) %>%
      dplyr::summarise(
        x_q1 = min(x, na.rm = TRUE) + (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/3,
        x_q2 = min(x, na.rm = TRUE) + 2*(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))/3,
        y_q1 = min(y, na.rm = TRUE) + (max(y, na.rm = TRUE) - min(y, na.rm = TRUE))/3,
        y_q2 = min(y, na.rm = TRUE) + 2*(max(y, na.rm = TRUE) - min(y, na.rm = TRUE))/3,
        .groups = "drop"
      )
  }

  return(boundaries)
}

#' Create facet annotations with chi-square p-values
create_facet_annotations <- function(results_df, variable_pairs) {
  tryCatch({
    # Filter for XOR-detected pairs and create annotations
    facet_colors <- results_df %>%
      dplyr::filter(xor_shape_detected == TRUE) %>%
      dplyr::mutate(
        pair_key = paste(var1, var2, sep = "||"),
        pair_key = ifelse("xor_pattern" %in% names(results_df) &
                            !is.na(xor_pattern) & xor_pattern == "rotated",
                          paste0(pair_key, "||rotated"), pair_key),
        chi_sq_label = dplyr::case_when(
          is.na(chi_sq_p_value) ~ "Chi sq: N/A",
          round(chi_sq_p_value, 3) == 0 ~ paste("Chi sq =", formatC(chi_sq_p_value, format = "e", digits = 3)),
          TRUE ~ paste("Chi sq =", round(chi_sq_p_value, 3))
        )
      ) %>%
      dplyr::select(pair_key, chi_sq_label, chi_sq_p_value)

    # Rename for consistency with plot data
    names(facet_colors)[names(facet_colors) == "pair_key"] <- "var1_var2"

    return(facet_colors)
  }, error = function(e) {
    warning("Could not create facet annotations: ", e$message, call. = FALSE)
    return(data.frame(var1_var2 = character(0), chi_sq_label = character(0),
                      chi_sq_p_value = numeric(0)))
  })
}

#' Create color mapping for facet strips
create_strip_colors <- function(facet_colors) {
  if (nrow(facet_colors) == 0) {
    return(NULL)
  }

  unique_facet_values <- unique(facet_colors$var1_var2)
  score_colors <- grDevices::colorRampPalette(c("gold", "cornsilk"))(length(unique_facet_values))
  names(score_colors) <- unique_facet_values

  return(score_colors)
}

#' Create the actual spaghetti plot
create_spaghetti_plot <- function(long_data, facet_colors, results_name) {
  # Create strip colors
  score_colors <- create_strip_colors(facet_colors)

  # Define custom strip styling
  if (!is.null(score_colors)) {
    strip <- ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(fill = score_colors)
    )
  } else {
    strip <- ggh4x::strip_themed()
  }

  # Build plot
  p <- ggplot(long_data, aes(x = variable, y = value, color = as.factor(class), group = ID)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_line(linewidth = 0.3, alpha = 0.8) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::theme_light(base_size = 14) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(colour = "black"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggthemes::scale_color_colorblind() +
    ggplot2::labs(
      title = "Spaghetti Plot of Variable Pairs",
      subtitle = paste("Derived from results:", results_name),
      x = "Variables",
      y = "Values",
      color = "Class",
      caption = "Facet by variable pairs: var1_var2"
    ) +
    ggh4x::facet_wrap2(~var1_var2, scales = "free_y", strip = strip)

  # Add chi-square annotations if available
  if (nrow(facet_colors) > 0 && "chi_sq_label" %in% names(facet_colors)) {
    p <- p +
      ggplot2::geom_text(
        data = facet_colors,
        aes(x = -Inf, y = Inf, label = chi_sq_label),
        inherit.aes = FALSE,
        hjust = -0.1, vjust = 1.5,
        size = 4,
        color = "black"
      )
  }

  return(p)
}

#' Create the actual XY scatter plot
create_xy_plot <- function(combined_data, boundaries, facet_colors, results_name) {
  # Create strip colors
  score_colors <- create_strip_colors(facet_colors)

  # Define custom strip styling
  if (!is.null(score_colors)) {
    strip <- ggh4x::strip_themed(
      background_x = ggh4x::elem_list_rect(fill = score_colors)
    )
  } else {
    strip <- ggh4x::strip_themed()
  }

  # Build plot
  p <- ggplot(combined_data, aes(x = x, y = y, color = as.factor(class))) +
    geom_point(alpha = 0.7, size = 2) +
    ggplot2::geom_vline(data = boundaries, aes(xintercept = x_q1), linetype = "dashed", color = "grey40") +
    ggplot2::geom_vline(data = boundaries, aes(xintercept = x_q2), linetype = "dashed", color = "grey40") +
    ggplot2::geom_hline(data = boundaries, aes(yintercept = y_q1), linetype = "dashed", color = "grey40") +
    ggplot2::geom_hline(data = boundaries, aes(yintercept = y_q2), linetype = "dashed", color = "grey40") +
    ggh4x::facet_wrap2(~var1_var2, scales = "free", strip = strip) +
    ggthemes::scale_color_colorblind() +
    ggplot2::theme_light(base_size = 14) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(color = "black"),
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "XY Scatter Plots of Variable Pairs",
      subtitle = paste("Derived from results:", results_name),
      x = "Variable 1",
      y = "Variable 2",
      color = "Class",
      caption = "Facet by variable pairs: var1_var2"
    )

  return(p)
}
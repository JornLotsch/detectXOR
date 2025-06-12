#' @title Generate Spaghetti Plots for XOR-Detected Pairs
#' @description Creates connected line plots for variable pairs flagged with XOR patterns.
#' @param results Either a data frame from \code{detect_xor$results_df} or the full list returned by \code{detect_xor}.
#' @param data Original dataset containing variables and classes.
#' @param class_col Character specifying the class column name.
#' @param scale_data Logical indicating whether to scale variables. Default: TRUE.
#' @return ggplot object (does not save files automatically).
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_line theme_light labs theme element_rect element_text expansion
#' @importFrom ggh4x facet_wrap2 strip_themed elem_list_rect
#' @importFrom dplyr filter mutate
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom ggthemes scale_color_colorblind
#' @importFrom glue glue
#' @importFrom utils globalVariables
generate_spaghetti_plot_from_results <- function(results, data, class_col, scale_data = TRUE) {
  # Handle both list and dataframe inputs
  if (is.list(results) && "results_df" %in% names(results)) {
    results_df <- results$results_df
    results_name <- deparse(substitute(results))
  } else if (is.data.frame(results)) {
    results_df <- results
    results_name <- deparse(substitute(results))
  } else {
    stop("Results must be either a data frame or a list containing 'results_df'")
  }

  # Validate inputs
  if (!all(c("var1", "var2", "xor_shape_detected") %in% colnames(results_df))) {
    stop("Results must contain var1, var2, and xor_shape_detected columns")
  }

  variable_pairs <- results_df[results_df$xor_shape_detected == TRUE, c("var1", "var2")]

  if (nrow(variable_pairs) == 0) {
    empty_plot <- ggplot() +
      geom_text(aes(x = 0.5, y = 0.5, label = "No XOR patterns detected"), size = 5) +
      theme_void()
    return(empty_plot)
  }

  # Data processing
  long_data <- do.call(rbind, lapply(seq_len(nrow(variable_pairs)), function(i) {
    rotate_condition <- ("xor_pattern" %in% names(results_df) &&
      results_df$xor_pattern[which(rownames(results_df) == rownames(variable_pairs)[i])] == "rotated")

    pair <- variable_pairs[i, ]
    pair_data <- data[, c(class_col, pair$var1, pair$var2)]

    if (rotate_condition) {
      pair_data_original <- pair_data
      pair_data[, 2] <- pair_data_original[, 2] + pair_data_original[, 3]
      pair_data[, 3] <- pair_data_original[, 2] - pair_data_original[, 3]
    }

    if (scale_data) pair_data[-which(names(pair_data) == class_col)] <- scale(pair_data[-which(names(pair_data) == class_col)])

    colnames(pair_data) <- c("class", "var1", "var2")
    pair_data$var1_var2 <- paste(pair$var1, pair$var2, sep = "||")
    if (rotate_condition) pair_data$var1_var2 <- paste(pair_data$var1_var2, "rotated", sep = "||")
    pair_data$ID <- seq_len(nrow(pair_data))
    return(pair_data)
  }))

  # Reshape the data into long format for plotting
  long_data <- melt(long_data, id.vars = c("ID", "class", "var1_var2"))

  # Create annotation data
  facet_colors <- results_df %>%
    tibble::rownames_to_column(var = "var1_var2") %>%
    filter(xor_shape_detected == TRUE) %>%
    mutate(
      chi_sq_p_value = chi_sq_p_value,
      var1_var2 = ifelse(xor_pattern == "rotated", paste0(var1_var2, "||rotated"), var1_var2)
    ) %>%
    select(var1_var2, chi_sq_p_value)

  # Create color mapping
  unique_facet_values <- unique(facet_colors$var1_var2)
  score_colors <- colorRampPalette(c("gold", "cornsilk"))(length(unique_facet_values))
  names(score_colors) <- unique_facet_values

  # Define custom strip styling for facet headers
  strip <- strip_themed(
    background_x = elem_list_rect(fill = score_colors)
  )

  # Build plot
  p <- ggplot(long_data, aes(x = variable, y = value, color = as.factor(class), group = ID)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_line(linewidth = 0.3, alpha = 0.8) +
    scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    theme_light(base_size = 14) +
    theme(
      strip.text = element_text(colour = "black"),
      legend.position = "bottom"
    ) +
    scale_color_colorblind() +
    labs(
      title = "Spaghetti Plot of Variable Pairs",
      subtitle = paste("Derived from results:", results_name),
      x = "Variables",
      y = "Values",
      color = "Class",
      caption = "Facet by variable pairs: var1_var2"
    ) +
    facet_wrap2(~var1_var2, scales = "free_y", strip = strip) +
    geom_text(
      data = facet_colors,
      aes(x = -Inf, y = Inf, label = paste("Chi square =", ifelse(round(chi_sq_p_value, 3) == 0,
                                                                  formatC(chi_sq_p_value, format = "e", digits = 3),
                                                                  round(chi_sq_p_value, 3))
      )),
      inherit.aes = FALSE,
      hjust = -0.1, vjust = 1.5,
      size = 4,
      color = "black"
    )

  return(p)
}

#' @title Generate XY Scatter Plots for XOR-Detected Pairs
#' @description Creates scatterplots with decision boundaries for detected XOR patterns.
#' @param results Either a data frame from \code{detect_xor$results_df} or the full list returned by \code{detect_xor}.
#' @inheritParams generate_spaghetti_plot_from_results
#' @param quantile_lines Numeric vector of quantiles for reference lines. Default: c(1/3, 2/3).
#' @param line_method Method for boundary calculation ("quantile" or "range"). Default: "quantile".
#' @return ggplot object (does not save files automatically).
#' @export
#' @importFrom ggplot2 geom_vline geom_hline theme element_rect element_text
#' @importFrom dplyr group_by summarise
#' @importFrom stats quantile
generate_xy_plot_from_results <- function(results, data, class_col,
                                          scale_data = TRUE,
                                          quantile_lines = c(1/3, 2/3),
                                          line_method = "quantile") {

  # Handle both list and dataframe inputs
  if (is.list(results) && "results_df" %in% names(results)) {
    results_df <- results$results_df
    results_name <- deparse(substitute(results))
  } else if (is.data.frame(results)) {
    results_df <- results
    results_name <- deparse(substitute(results))
  } else {
    stop("Results must be either a data frame or a list containing 'results_df'")
  }

  # Input validation
  if (!line_method %in% c("quantile", "range")) {
    stop("line_method must be either 'quantile' or 'range'")
  }

  variable_pairs <- results_df[results_df$xor_shape_detected == TRUE, c("var1", "var2")]

  if (nrow(variable_pairs) == 0) {
    empty_plot <- ggplot() +
      geom_text(aes(x = 0.5, y = 0.5, label = "No XOR patterns detected"), size = 5) +
      theme_void()
    return(empty_plot)
  }

  # Data processing - similar structure to working version
  combined_data <- do.call(rbind, lapply(seq_len(nrow(variable_pairs)), function(i) {
    pair <- variable_pairs[i, ]
    pair_data <- data[, c(class_col, pair$var1, pair$var2)]

    rotate_condition <- ("xor_pattern" %in% names(results_df) &&
      results_df$xor_pattern[i] == "rotated")

    if (rotate_condition) {
      pair_data_original <- pair_data
      pair_data[, 2] <- pair_data_original[, 2] + pair_data_original[, 3]
      pair_data[, 3] <- pair_data_original[, 2] - pair_data_original[, 3]
    }

    if (scale_data) {
      cols_to_scale <- setdiff(names(pair_data), class_col)
      pair_data[cols_to_scale] <- scale(pair_data[cols_to_scale])
    }

    colnames(pair_data) <- c("class", "x", "y")
    pair_data$var1_var2 <- paste(pair$var1, pair$var2, sep = "||")
    if (rotate_condition) {
      pair_data$var1_var2 <- paste(pair_data$var1_var2, "rotated", sep = "||")
    }
    return(pair_data)
  }))

  # Boundary calculation
  if (line_method == "quantile") {
    boundaries <- combined_data %>%
      group_by(var1_var2) %>%
      summarise(
        x_q1 = quantile(x, quantile_lines[1], na.rm = TRUE),
        x_q2 = quantile(x, quantile_lines[2], na.rm = TRUE),
        y_q1 = quantile(y, quantile_lines[1], na.rm = TRUE),
        y_q2 = quantile(y, quantile_lines[2], na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    boundaries <- combined_data %>%
      group_by(var1_var2) %>%
      summarise(
        x_q1 = min(x) + (max(x) - min(x))/3,
        x_q2 = min(x) + 2*(max(x) - min(x))/3,
        y_q1 = min(y) + (max(y) - min(y))/3,
        y_q2 = min(y) + 2*(max(y) - min(y))/3,
        .groups = "drop"
      )
  }

  # Create annotation data
  facet_colors <- results_df %>%
    tibble::rownames_to_column(var = "var1_var2") %>%
    filter(xor_shape_detected == TRUE) %>%
    mutate(
      chi_sq_p_value = chi_sq_p_value,
      var1_var2 = ifelse(xor_pattern == "rotated", paste0(var1_var2, "||rotated"), var1_var2)
    ) %>%
    select(var1_var2, chi_sq_p_value)

  # Create color mapping
  unique_facet_values <- unique(facet_colors$var1_var2)
  score_colors <- colorRampPalette(c("gold", "cornsilk"))(length(unique_facet_values))
  names(score_colors) <- unique_facet_values

  # Define custom strip styling for facet headers
  strip <- strip_themed(
    background_x = elem_list_rect(fill = score_colors)
  )

  # Build plot
  p <- ggplot(combined_data, aes(x = x, y = y, color = as.factor(class))) +
    geom_point(alpha = 0.7, size = 2) +
    geom_vline(data = boundaries, aes(xintercept = x_q1), linetype = "dashed", color = "grey40") +
    geom_vline(data = boundaries, aes(xintercept = x_q2), linetype = "dashed", color = "grey40") +
    geom_hline(data = boundaries, aes(yintercept = y_q1), linetype = "dashed", color = "grey40") +
    geom_hline(data = boundaries, aes(yintercept = y_q2), linetype = "dashed", color = "grey40") +
    facet_wrap2(~var1_var2, scales = "free", strip = strip) +
    scale_color_colorblind() +
    theme_light(base_size = 14) +
    theme(
      strip.text = element_text(color = "black"),
      legend.position = "bottom"
    ) +
    labs(
      title = "XY Scatter Plots of Variable Pairs",
      subtitle = paste("Derived from results:", results_name),
      x = "Variable 1",
      y = "Variable 2",
      color = "Class",
      caption = "Facet by variable pairs: var1_var2"
    )

  return(p)
}

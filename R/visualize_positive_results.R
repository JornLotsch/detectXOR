#' @title Generate Spaghetti Plots for XOR-Detected Pairs
#' @description Creates connected line plots for variable pairs flagged with XOR patterns.
#' @param results Data frame from \code{detect_xor} containing detected pairs.
#' @param data Original dataset containing variables and classes.
#' @param class_col Character specifying the class column name.
#' @param output_dir Directory path for saving plots. Default: current directory.
#' @param scale_data Logical indicating whether to scale variables. Default: TRUE.
#' @return ggplot object and saves SVG file.
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_line theme_light labs theme element_rect element_text
#' @importFrom ggh4x facet_wrap2 strip_themed
#' @importFrom dplyr filter mutate
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom scales colorblind_pal
#' @importFrom glue glue
#' @importFrom utils globalVariables
generate_spaghetti_plot_from_results <- function(results, data, class_col, output_dir = ".", scale_data = TRUE) {
  # Validate inputs
  if (!all(c("var1", "var2", "xor_shape_detected") %in% colnames(results))) {
    stop("Results must contain var1, var2, and xor_shape_detected columns")
  }

  # Create output directory if needed
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  results_name <- deparse(substitute(results))
  variable_pairs <- results[results$xor_shape_detected, c("var1", "var2")]

  if (nrow(variable_pairs) == 0) {
    empty_plot <- ggplot() +
      geom_text(aes(0.5, 0.5, label = "No XOR patterns detected"), size = 5) +
      theme_void()

    ggsave(
      file.path(output_dir, glue("{results_name}_spaghetti_plot.svg")),
      empty_plot,
      width = 5, height = 5
    )
    return(invisible(empty_plot))
  }

  # Data processing
  long_data <- do.call(rbind, lapply(seq_len(nrow(variable_pairs)), function(i) {
    pair <- variable_pairs[i, ]
    pair_data <- data[, c(class_col, pair$var1, pair$var2)]

    if ("xor_pattern" %in% colnames(results) && results$xor_pattern[i] == "rotated") {
      pair_data_original <- pair_data
      pair_data[, 2] <- pair_data_original[, 2] + pair_data_original[, 3]
      pair_data[, 3] <- pair_data_original[, 2] - pair_data_original[, 3]
    }

    if (scale_data) pair_data[-1] <- scale(pair_data[-1])

    colnames(pair_data) <- c("class", "var1", "var2")
    pair_data$pair_id <- glue("{pair$var1}_{pair$var2}")
    if ("xor_pattern" %in% colnames(results) && results$xor_pattern[i] == "rotated") {
      pair_data$pair_id <- glue("{pair_data$pair_id}_rotated")
    }
    pair_data$ID <- seq_len(nrow(pair_data))
    melt(pair_data, id.vars = c("ID", "class", "pair_id"))
  }))

  # Create annotation data
  annotation_data <- results %>%
    tibble::rownames_to_column("pair_id") %>%
    filter(xor_shape_detected) %>%
    mutate(pair_id = ifelse(xor_pattern == "rotated", glue("{pair_id}_rotated"), pair_id))

  # Create color mapping
  unique_pairs <- unique(annotation_data$pair_id)
  color_palette <- colorRampPalette(c("gold", "cornsilk"))(length(unique_pairs))
  names(color_palette) <- unique_pairs

  # Build plot
  p <- ggplot(long_data, aes(variable, value, color = class, group = ID)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_line(linewidth = 0.3, alpha = 0.8) +
    facet_wrap2(~pair_id, scales = "free_y",
                strip = strip_themed(background_x = lapply(long_data$pair_id, function(x) {
                  elem_list_rect(fill = color_palette[x])
                }))) +
    geom_text(
      data = annotation_data,
      aes(x = -Inf, y = Inf, label = glue("\\u00cf\\u0087\\u00c2\\u00b2={signif(chi_sq_p_value,3)}") ),
      inherit.aes = FALSE, hjust = -0.1, vjust = 1.5, size = 3
    ) +
    labs(title = "XOR Pattern Spaghetti Plots", x = "Variables", y = "Values") +
    theme_light() +
    theme(legend.position = "bottom",
          strip.text = element_text(color = "black"),
          panel.border = element_rect(color = "grey80", fill = NA))

  # Save with dynamic sizing
  n_facets <- length(unique(long_data$pair_id))
  ggsave(
    file.path(output_dir, glue("{results_name}_spaghetti_plot.svg")),
    plot = p,
    width = max(7, 3 * ceiling(sqrt(n_facets))),
    height = max(5, 3 * ceiling(n_facets/ceiling(sqrt(n_facets)))),
    limitsize = FALSE
  )

  return(invisible(p))
}

#' @title Generate XY Scatter Plots for XOR-Detected Pairs
#' @description Creates scatterplots with decision boundaries for detected XOR patterns.
#' @inheritParams generate_spaghetti_plot_from_results
#' @param quantile_lines Numeric vector of quantiles for reference lines. Default: c(1/3, 2/3).
#' @param line_method Method for boundary calculation ("quantile" or "range"). Default: "quantile".
#' @return ggplot object and saves SVG file.
#' @export
#' @importFrom ggplot2 geom_vline geom_hline theme element_rect element_text
#' @importFrom dplyr group_by summarise
#' @importFrom stats quantile
generate_xy_plot_from_results <- function(results, data, class_col, output_dir = ".",
                                          scale_data = TRUE,
                                          quantile_lines = c(1/3, 2/3),
                                          line_method = "quantile") {

  # Input validation
  if (!line_method %in% c("quantile", "range")) {
    stop("line_method must be either 'quantile' or 'range'")
  }

  results_name <- deparse(substitute(results))
  variable_pairs <- results[results$xor_shape_detected, c("var1", "var2")]

  if (nrow(variable_pairs) == 0) {
    empty_plot <- ggplot() +
      geom_text(aes(0.5, 0.5, label = "No XOR patterns detected"), size = 5) +
      theme_void()

    ggsave(
      file.path(output_dir, glue("{results_name}_xy_plot.svg")),
      empty_plot,
      width = 5, height = 5
    )
    return(invisible(empty_plot))
  }

  # Data processing
  plot_data <- do.call(rbind, lapply(seq_len(nrow(variable_pairs)), function(i) {
    pair <- variable_pairs[i, ]
    pair_data <- data[, c(class_col, pair$var1, pair$var2)]

    if ("xor_pattern" %in% colnames(results) && results$xor_pattern[i] == "rotated") {
      pair_data_original <- pair_data
      pair_data[, 2] <- pair_data_original[, 2] + pair_data_original[, 3]
      pair_data[, 3] <- pair_data_original[, 2] - pair_data_original[, 3]
    }

    if (scale_data) pair_data[-1] <- scale(pair_data[-1])

    colnames(pair_data) <- c("class", "x", "y")
    pair_data$pair_id <- glue("{pair$var1}_{pair$var2}")
    if ("xor_pattern" %in% colnames(results) && results$xor_pattern[i] == "rotated") {
      pair_data$pair_id <- glue("{pair_data$pair_id}_rotated")
    }
    pair_data
  }))

  # Boundary calculation
  boundaries <- if (line_method == "quantile") {
    plot_data %>%
      group_by(pair_id) %>%
      summarise(
        x_q1 = quantile(x, quantile_lines[1], na.rm = TRUE),
        x_q2 = quantile(x, quantile_lines[2], na.rm = TRUE),
        y_q1 = quantile(y, quantile_lines[1], na.rm = TRUE),
        y_q2 = quantile(y, quantile_lines[2], na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    plot_data %>%
      group_by(pair_id) %>%
      summarise(
        x_q1 = min(x) + (max(x) - min(x))/3,
        x_q2 = min(x) + 2*(max(x) - min(x))/3,
        y_q1 = min(y) + (max(y) - min(y))/3,
        y_q2 = min(y) + 2*(max(y) - min(y))/3,
        .groups = "drop"
      )
  }

  # Create annotation data
  annotation_data <- results %>%
    tibble::rownames_to_column("pair_id") %>%
    filter(xor_shape_detected) %>%
    mutate(pair_id = ifelse(xor_pattern == "rotated", glue("{pair_id}_rotated"), pair_id))

  # Create color mapping
  unique_pairs <- unique(annotation_data$pair_id)
  color_palette <- colorRampPalette(c("gold", "cornsilk"))(length(unique_pairs))
  names(color_palette) <- unique_pairs

  # Build plot
  p <- ggplot(plot_data, aes(x, y, color = class)) +
    geom_point(alpha = 0.7) +
    geom_vline(data = boundaries, aes(xintercept = x_q1), linetype = 2) +
    geom_vline(data = boundaries, aes(xintercept = x_q2), linetype = 2) +
    geom_hline(data = boundaries, aes(yintercept = y_q1), linetype = 2) +
    geom_hline(data = boundaries, aes(yintercept = y_q2), linetype = 2) +
    facet_wrap2(~pair_id, scales = "free",
                strip = strip_themed(background_x = lapply(plot_data$pair_id, function(x) {
                  elem_list_rect(fill = color_palette[x])
                }))) +
    geom_text(
      data = annotation_data,
      aes(x = -Inf, y = Inf, label = glue("\\u00cf\\u0087\\u00c2\\u00b2={signif(chi_sq_p_value,3)}") ),
      inherit.aes = FALSE, hjust = -0.1, vjust = 1.5, size = 3
    ) +
    labs(title = "XOR Pattern Scatterplots", x = "Variable 1", y = "Variable 2") +
    theme_light() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(color = "black"),
      panel.border = element_rect(color = "grey80", fill = NA)
    )

  # Save with dynamic sizing
  n_facets <- length(unique(plot_data$pair_id))
  ggsave(
    file.path(output_dir, glue("{results_name}_xy_plot.svg")),
    plot = p,
    width = max(7, 3 * ceiling(sqrt(n_facets))),
    height = max(5, 3 * ceiling(n_facets/ceiling(sqrt(n_facets)))),
    limitsize = FALSE
  )

  return(invisible(p))
}

# Suppress package check notes
utils::globalVariables(c(
  "pair_id", "variable", "value", "x", "y",
  "x_q1", "x_q2", "y_q1", "y_q2", "chi_sq_p_value"
))


############################################## Plot data ##################################

# Required Libraries ---------------------------------------------------------
REQUIRED_PACKAGES <- c(
  "ggplot2", "grid", "gridExtra", "GGally", "ggthemes",
  "reshape2", "cowplot", "ggh4x", "ggvoronoi", "parallel", "ComplexHeatmap",
  "dplyr", "smotefamily", "patchwork", "viridis", "scales", "tidyr"
)

# Install and load required packages
install_and_load_packages <- function(packages) {
  suppressMessages({
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        tryCatch({
          install.packages(pkg)
        }, error = function(e) {
          message("Error installing package: ", pkg)
        })
      }
      # Try to load the package
      tryCatch({
        library(pkg, character.only = TRUE)
      }, error = function(e) {
        message("Error loading package: ", pkg)
      })
    }
  })
}

# Install and load packages
install_and_load_packages(REQUIRED_PACKAGES)

base_dir <- "/home/joern/.Datenplatte/Joerns Dateien/Aktuell/FeatureSelectionInformationTheory/08AnalyseProgramme/R/"
setwd(base_dir)

# Helper function
split_range_three <- function(x) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  break1 <- min_val + (max_val - min_val) / 3
  break2 <- min_val + 2 * (max_val - min_val) / 3
  c(break1, break2)
}

# XOR Matrix Plot
XOR_matrix_plot <- function(data, Cls, data_name = "dataset", true_pairs_list = NULL, true_pairs_color = "lightyellow") {
  
  # Initialize true_pairs as an empty list if it doesn't exist
  if (is.null(true_pairs_list)) {
    true_pairs_list <- list()  # Empty list if no pairs are defined
  }
  
  # Create logical matrix (handles empty lists gracefully)
  create_logical_matrix <- function(data, true_pairs_list) {
    vars <- colnames(data)
    n <- length(vars)
    logical_matrix <- matrix(FALSE, nrow = n, ncol = n, dimnames = list(vars, vars))
    
    if (length(true_pairs_list) > 0) {  # Skip processing if list is empty
      for (pair in true_pairs_list) {
        var1 <- pair[1]
        var2 <- pair[2]
        idx1 <- which(vars == var1)
        idx2 <- which(vars == var2)
        logical_matrix[idx1, idx2] <- TRUE
        logical_matrix[idx2, idx1] <- TRUE
      }
    }
    
    return(logical_matrix)
  }
  
  logical_matrix <- create_logical_matrix(data, true_pairs_list)
  
  # 1. Create base plot 
  plot <- GGally::ggpairs(
    data = as.data.frame(data),
    aes(color = factor(Cls)),
    upper = list(continuous = GGally::wrap("cor", method = "kendall", size = 4)),
    lower = list(continuous = "points"),
    diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.3)),
    progress = FALSE
  ) +
    labs(title = paste("Scatter Plot Matrix of", data_name)) +
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind() +
    theme(
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black"),
      panel.border = element_rect(color = "grey70", fill = NA)
    )
  
  # 2. Modify all panels
  for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
      if (i == j) {
        # Diagonal panels: Transparent background
        plot[i, j] <- plot[i, j] +
          theme(
            panel.background = element_rect(fill = NA, color = NA),
            panel.border = element_rect(color = "grey70", fill = NA)
          )
      } else if (i > j) {
        # Lower triangle panels
        if (logical_matrix[j, i]) {
          # TRUE cells: Pink background
          plot[i, j] <- plot[i, j] +
            theme(
              panel.background = element_rect(fill = true_pairs_color, color = NA),
              panel.border = element_rect(color = "grey70", fill = NA)
            )
        } else {
          # FALSE cells: Transparent background
          plot[i, j] <- plot[i, j] +
            theme(
              panel.background = element_rect(fill = NA, color = NA),
              panel.border = element_rect(color = "grey70", fill = NA)
            )
        }
      }
    }
  }
  plot
}

# XOR Density Plot
XOR_density_plot <- function(data, Cls, fill = FALSE ) {
  # Reshape data for long format
  if (is.null(rownames(data))) rownames(data) <- seq_len(nrow(data)) # Assign row indices
  data_long <- reshape2::melt(cbind.data.frame(ID = rownames(data), Cls = factor(Cls), data))

  dens_alpha <- ifelse(fill, 0.1, 0)
  ggplot(data_long, aes(x = value, color = as.factor(Cls), fill = as.factor(Cls))) +
    geom_density(alpha = dens_alpha) +
    facet_wrap(variable ~ ., nrow = 1) +
    theme_light() +
    guides(color = "none", fill = "none") +
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind() +
    labs(title = "Variable Density") +
    theme(
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black")
    )
}

# XOR Raw Data Plot
XOR_raw_plot <- function(data, Cls) {
  # Reshape data to long format
  if (is.null(rownames(data))) rownames(data) <- seq_len(nrow(data)) # Assign row indices
  data_long <- reshape2::melt(cbind.data.frame(ID = rownames(data), Cls = factor(Cls), data))

  ggplot(data_long, aes(x = variable, y = value, color = as.factor(Cls), fill = as.factor(Cls))) +
    geom_violin(alpha = 0.1, width = 0.4, position = position_dodge(width = 0.8)) +
    geom_boxplot(alpha = 0.3, width = 0.2, position = position_dodge(width = 0.8), outlier.shape = NA) +
    geom_jitter(alpha = 0.6, size = 0.3, position = position_dodge(width = 0.8)) +
    ggpubr::stat_compare_means(vjust = .2) +
    facet_wrap(variable ~ ., nrow = 1, scales = "free") +
    theme_light() +
    guides(color = "none") +
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind() +
    labs(title = "Raw Data per Variable", fill = "Class") +
    theme(
      legend.position = "bottom", legend.direction = "horizontal",
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black")
    )
}

# XOR Spaghetti Plot
XOR_spaghetti_plot <- function(data, Cls) {
  # Reshape data for long format
  if (is.null(rownames(data))) rownames(data) <- seq_len(nrow(data)) # Assign row indices
  data_long <- reshape2::melt(cbind.data.frame(ID = rownames(data), Cls = factor(Cls), data))

  ggplot(data_long, aes(x = variable, y = value, color = Cls, group = ID)) +
    geom_point() +
    geom_line(linewidth = .1) +
    scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +
    theme_light() +
    guides(color = "none") +
    ggthemes::scale_color_colorblind() +
    labs(title = "Raw Data Connected") +
    theme(axis.text.x = element_text(angle = 90, hjust = 0.5))
}

# Wrapper Function: generate_XOR_plots
generate_XOR_plots <- function(data, Cls, data_name = "dataset", true_pairs_list = NULL, true_pairs_color = "lightyellow", plot_dir = ".") {
  # Dynamically adjust plot width and height for ggpairs plot
  num_vars <- ncol(data) -1
  cell_size <- 3 # The size of each cell (in inches), adjust as needed
  plot_width <- num_vars * cell_size
  plot_height <- num_vars * cell_size

  # Exploration Plot 1 uses 2/3 of the height of the full `ggpairs` plot.
  exploration_plot1_height <- 2 / 3 * plot_height

  # Create the plots
  matrix_plot <- XOR_matrix_plot(data, Cls, data_name, true_pairs_list, true_pairs_color) # Matrix plot (ggpairs)
  density_plot <- XOR_density_plot(data, Cls) # Density plot
  raw_plot <- XOR_raw_plot(data, Cls) # Raw data plot
  spaghetti_plot <- XOR_spaghetti_plot(data, Cls) # Spaghetti plot

  # Combined Plot 1: Raw, Spaghetti, Density
  XOR_exploration_plot1 <- cowplot::plot_grid(
    raw_plot,
    spaghetti_plot,
    density_plot,
    labels = "AUTO",
    ncol = 1,
    align = "h",
    axis = "lr"
  )
  # print(XOR_exploration_plot1)

  # Combined Plot 2: GGPairs + Spaghetti Plot
  XOR_exploration_plot2 <- cowplot::plot_grid(
    raw_plot,
    grid::grid.grabExpr(print(matrix_plot)), # Matrix plot as one element
    spaghetti_plot,
    labels = "AUTO",
    ncol = 1,
    rel_heights = c(3, 7, 3), # Larger height for the matrix plot
    align = "h",
    axis = "lr"
  )
  # print(XOR_exploration_plot2)

  # Save the plots
  filename1 <- paste0(data_name, "_exploration_plot1.svg")
  filename2 <- paste0(data_name, "_exploration_plot2.svg")
  filename_matrix <- paste0(data_name, "_matrix_plot.svg") # Stand-alone matrix plot

  # Saving Combined Plot 1 with 2/3 height of the full plot
  ggsave(file.path(plot_dir, filename1),
         plot = XOR_exploration_plot1,
         width = plot_width,
         height = 1.75 * exploration_plot1_height  ,
         limitsize = FALSE)

  # Saving Combined Plot 2 with the same dimensions as the GGPairs plot
  ggsave(file.path(plot_dir, filename2),
         plot = XOR_exploration_plot2,
         width = plot_width,
         height = 1.75 * plot_height,
         limitsize = FALSE)

  # Saving the matrix plot (ggpairs) with quadratic cells
  ggsave(file.path(plot_dir, filename_matrix),
         plot = matrix_plot,
         width = plot_width,
         height = plot_height,
         limitsize = FALSE)

  # Return all plot objects as a list
  return(list(
    XOR_matrix_plot = matrix_plot,
    XOR_density_plot = density_plot,
    XOR_raw_plot = raw_plot,
    XOR_spaghetti_plot = spaghetti_plot,
    XOR_exploration_plot1 = XOR_exploration_plot1,
    XOR_exploration_plot2 = XOR_exploration_plot2,
    saved_files = list(filename1, filename2, filename_matrix) # Include saved file paths
  ))
}


# Function Definition
generate_spaghetti_plot_from_results <- function(results, data, class_col, output_dir = ".", scale_data = TRUE) {
  # Automatically derive the name of the results data frame
  results_name <- deparse(substitute(results))

  # Step 1: Filter pairs where xor_shape_detected == TRUE
  variable_pairs <- results[results$xor_shape_detected == TRUE, c("var1", "var2")]


  # Handle case where no variable pairs satisfy xor_shape_detected == TRUE
  if (nrow(variable_pairs) == 0) {
    empty_plot <- ggplot() +
      geom_text(aes(x = 0.5, y = 0.5, label = "No variable pairs with xor_shape_detected == TRUE found.")) +
      theme_void() +
      labs(title = "Empty Plot: No Data to Display")

    # Save the empty plot as an SVG
    ggsave(
      filename = file.path(output_dir, paste0(results_name, "_spaghetti_plot.svg")),
      plot = empty_plot,
      device = "svg",
      width = 5, # Default width
      height = 5, # Default height
      limitsize = FALSE
    )

    return(empty_plot) # Return the empty plot
  }


  # Step 2: Assemble a data frame for plotting
  long_data <- do.call(rbind, lapply(seq_len(nrow(variable_pairs)), function(i) {
    rotate_condition <- ("xor_pattern" %in% names(results) && 
                           results$xor_pattern[which(rownames(results) == rownames(variable_pairs)[i])] == "rotated")

    pair <- variable_pairs[i,]
    pair_data <- data[, c(class_col, pair$var1, pair$var2)]
    if (rotate_condition) {
      pair_data_original <- pair_data <- data[, c(class_col, pair$var1, pair$var2)]
      pair_data[, 2] <- pair_data_original[, 2] + pair_data_original[, 3]
      pair_data[, 3] <- pair_data_original[, 2] - pair_data_original[, 3]
    }

    # Function to apply scaling and ranking locally
    transform_data <- function(data, scale_data, class_col) {
      if (scale_data) {
        data[-which(names(data) == class_col)] <- scale(data[-which(names(data) == class_col)])
      }
      return(data)
    }
    pair_data <- transform_data(pair_data, scale_data, class_col)

    colnames(pair_data) <- c("class", "var1", "var2")
    pair_data$var1_var2 <- paste(pair$var1, pair$var2, sep = "_")
    if (rotate_condition) pair_data$var1_var2 <- paste(pair_data$var1_var2, "rotated", sep = "_")
    pair_data$ID <- seq_len(nrow(pair_data))
    return(pair_data)
  }))

  # Reshape the data into long format for plotting
  long_data <- melt(long_data, id.vars = c("ID", "class", "var1_var2"))

  # Create a data frame with facet-specific coordinates for ablines in the plot
  facet_colors <- results %>%
    tibble::rownames_to_column(var = "var1_var2") %>% # Convert row names to a column
    filter(xor_shape_detected == TRUE) %>%          # Filter rows where xor_shape_detected is TRUE
    mutate(
      chi_sq_p_value = chi_sq_p_value
    ) %>%
    select(var1_var2, chi_sq_p_value, xor_pattern) # Keep only relevant columns
  
  # Modify the var1_var2 column based on xor_shape_detected_rotated
  facet_colors <- facet_colors %>%
    mutate(var1_var2 = ifelse(xor_pattern == "rotated", paste0(var1_var2, "_rotated"), var1_var2))
  
  # Generate colors mapping from XOR score (e.g., "gold" to "cornsilk")
  unique_facet_values <- unique(facet_colors$var1_var2)
  score_colors <- colorRampPalette(c("gold", "cornsilk"))(length(unique_facet_values))
  names(score_colors) <- unique_facet_values
  
  # Define custom strip styling for facet headers
  strip <- strip_themed(
    background_x = elem_list_rect(fill = score_colors)
  )

  # Spaghetti plot with annotations
  spaghetti_plot <- ggplot(long_data, aes(x = variable, y = value, color = as.factor(class), group = ID)) +
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
    facet_wrap2(~var1_var2, scales = "free_y", strip = strip)  +
    
    # Add text annotations for chi_sq_p_value
    geom_text(
      data = facet_colors,
      aes(x = -Inf, y = Inf, label = paste("χ² =", ifelse(round(chi_sq_p_value, 3) == 0, 
                                                          formatC(chi_sq_p_value, format = "e", digits = 3), 
                                                          round(chi_sq_p_value, 3))
      )),
      inherit.aes = FALSE,
      hjust = -0.1, vjust = 1.5,
      size = 4,
      color = "black"
    )
  
  # Dynamically adapt plot size
  num_facets <- length(unique(long_data$var1_var2))
  num_columns <- ceiling(sqrt(num_facets))
  num_rows <- ceiling(num_facets / num_columns)
  height_per_facet <- 4 
  width_per_facet <- 4  
  plot_height <- max(5, height_per_facet * num_rows)  
  plot_width <- max(7, width_per_facet * num_columns)  
  
  # Save the spaghetti plot as an SVG
  ggsave(
    filename = file.path(output_dir, paste0(results_name, "_spaghetti_plot.svg")),
    plot = spaghetti_plot,
    device = "svg",
    width = plot_width,
    height = 1.2*plot_height,
    limitsize = FALSE
  )
  
  # Return the final plot
  return(spaghetti_plot)
  
}

generate_xy_plot_from_results <- function(results, data, class_col, output_dir = ".", 
                                          scale_data = TRUE, 
                                          quantile_lines = c(1/3, 2/3),
                                          line_method = "quantile") {
  # Automatically derive the name of the results data frame
  results_name <- deparse(substitute(results))
  
  # Step 1: Filter pairs where xor_shape_detected == TRUE
  variable_pairs <- results[results$xor_shape_detected == TRUE, c("var1", "var2")]
  
  # Handle case where no variable pairs satisfy xor_shape_detected == TRUE
  if (nrow(variable_pairs) == 0) {
    empty_plot <- ggplot() +
      geom_text(aes(x = 0.5, y = 0.5, label = "No variable pairs with xor_shape_detected == TRUE found.")) +
      theme_void() +
      labs(title = "Empty Plot: No Data to Display")
    
    # Save the empty plot as an SVG
    ggsave(
      filename = file.path(output_dir, paste0(results_name, "_xy_plot.svg")),
      plot = empty_plot,
      device = "svg",
      width = 5,
      height = 5,
      limitsize = FALSE
    )
    
    return(empty_plot)
  }
  
  # Step 2: Assemble a data frame for plotting
  combined_data <- do.call(rbind, lapply(seq_len(nrow(variable_pairs)), function(i) {
    pair <- variable_pairs[i, ]
    pair_data <- data[, c(class_col, pair$var1, pair$var2)]
    
    # Check if rotation is needed and handle NA values explicitly
    rotate_condition <- ("xor_pattern" %in% names(results) && 
                           results$xor_pattern[i] == "rotated" )
    
    if (rotate_condition) {
      pair_data_original <- pair_data
      pair_data[, 2] <- pair_data_original[, 2] + pair_data_original[, 3]
      pair_data[, 3] <- pair_data_original[, 2] - pair_data_original[, 3]
    }
    
    # Scale data if required
    if (scale_data) {
      cols_to_scale <- setdiff(names(pair_data), class_col)
      pair_data[cols_to_scale] <- scale(pair_data[cols_to_scale])
    }
    
    colnames(pair_data) <- c("class", "x", "y")
    pair_data$var1_var2 <- paste(pair$var1, pair$var2, sep = "_")
    if (rotate_condition) {
      pair_data$var1_var2 <- paste(pair_data$var1_var2, "rotated", sep = "_")
    }
    return(pair_data)
  }))
  
  # Step 3: Add additional columns to combined data
  combined_data$temp_var1_var2 <- gsub("_rotated", "", combined_data$var1_var2)
  
  # Helper function to add columns (chi_sq_p_value, border_x, border_y)
  add_column_from_results <- function(temp_var1_var2, var1_var2, col_original, col_rotated) {
    sapply(seq_along(temp_var1_var2), function(i) {
      x <- temp_var1_var2[i]
      if (grepl("_rotated", var1_var2[i])) {
        if (x %in% rownames(results)) results[x, col_rotated] else NA
      } else {
        if (x %in% rownames(results)) results[x, col_original] else NA
      }
    })
  }

  # Clean up temporary column
  combined_data$temp_var1_var2 <- NULL
  
  # Create a data frame with facet-specific coordinates for ablines in the plot
  facet_colors <- results %>%
    tibble::rownames_to_column(var = "var1_var2") %>% # Convert row names to a column
    filter(xor_shape_detected == TRUE) %>%          # Filter rows where xor_shape_detected is TRUE
    mutate(
      chi_sq_p_value = chi_sq_p_value
    ) %>%
    select(var1_var2, chi_sq_p_value, xor_pattern) # Keep only relevant columns
  
  # Modify the var1_var2 column based on xor_shape_detected_rotated
  facet_colors <- facet_colors %>%
    mutate(var1_var2 = ifelse(xor_pattern == "rotated", paste0(var1_var2, "_rotated"), var1_var2))
  
  # Generate colors mapping from XOR score (e.g., "gold" to "cornsilk")
  unique_facet_values <- unique(facet_colors$var1_var2)
  score_colors <- colorRampPalette(c("gold", "cornsilk"))(length(unique_facet_values))
  names(score_colors) <- unique_facet_values
  
  
  # Define custom strip styling for facet headers
  strip <- strip_themed(
    background_x = elem_list_rect(fill = score_colors)
  )

  plot_with_quantile_lines <- function(
    data,
    x_var = "x",
    y_var = "y",
    facet_var = "var1_var2",
    class_var = "class",
    quantiles = quantile_lines,
    ...
  ) {
    # 1. Calculate quantiles per facet
    if (!is.null(quantile_lines)) {
      if (line_method == "quantile") {
        quantile_lines <- data %>%
          group_by(.data[[facet_var]]) %>%
          summarize(
            !!!setNames(
              lapply(quantiles, function(q) expr(quantile(!!sym(x_var), !!q))),
              paste0("x_q", seq_along(quantiles))
            ),
            !!!setNames(
              lapply(quantiles, function(q) expr(quantile(!!sym(y_var), !!q))),
              paste0("y_q", seq_along(quantiles))
            ),
            .groups = "drop"
          ) 
      } else {
        quantile_lines <- data %>%
          group_by(.data[[facet_var]]) %>%
          summarize(
            x_q1 = split_range_three(.data[[x_var]])[1],
            x_q2 = split_range_three(.data[[x_var]])[2],
            y_q1 = split_range_three(.data[[y_var]])[1],
            y_q2 = split_range_three(.data[[y_var]])[2],
            .groups = "drop"
          )
      }
      
    
    # 2. Reshape for plotting
    vlines <- quantile_lines %>%
      select(all_of(facet_var), starts_with("x_q")) %>%
      pivot_longer(
        cols = starts_with("x_q"),
        names_to = "which",
        values_to = "xintercept"
      )
    
    hlines <- quantile_lines %>%
      select(all_of(facet_var), starts_with("y_q")) %>%
      pivot_longer(
        cols = starts_with("y_q"),
        names_to = "which",
        values_to = "yintercept"
      )
    }
    # 3. Plot
    p <- ggplot(data, aes_string(x = x_var, y = y_var, color = paste0("as.factor(", class_var, ")"))) +
      geom_point(alpha = 0.7, size = 2) +
      facet_wrap2(~ var1_var2, scales = "free", strip = strip)   +
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
    if (!is.null(quantile_lines)) {
      p <- p +
      geom_vline(
        data = vlines,
        aes(xintercept = xintercept),
        linetype = "dashed",
        color = "grey40"
      ) +
        geom_hline(
          data = hlines,
          aes(yintercept = yintercept),
          linetype = "dashed",
          color = "grey40"
        ) }
    
    
    return(p)
  }
  
  # Step 4: Create the scatter plot
  xy_plot <- plot_with_quantile_lines(combined_data, quantile_lines = quantile_lines) 
  
  # Dynamically adapt plot size
  num_facets <- length(unique(combined_data$var1_var2))
  num_columns <- ceiling(sqrt(num_facets))
  num_rows <- ceiling(num_facets / num_columns)
  
  plot_height <- max(5, num_rows * 4)
  plot_width <- max(7, num_columns * 4)
  
  # Save the final plot as an SVG
  ggsave(
    filename = file.path(output_dir, paste0(results_name, "_xy_plot.svg")),
    plot = xy_plot,
    device = "svg",
    width = plot_width,
    height = plot_height * 1.2,
    limitsize = FALSE
  )
  
  # Return the final plot
  return(xy_plot)
}

# Function to generate and plot heatmap
generate_condition_heatmap <- function(data, clustering = TRUE, clustering_method = "ward.D2") {
  # Check if required columns exist in the data frame
  required_columns <- c("var1", "var2", "xor_shape_detected")
  missing_columns <- setdiff(required_columns, colnames(data))
  if (length(missing_columns) > 0) {
    stop(paste("The following required columns are missing in the data frame:",
               paste(missing_columns, collapse = ", ")))
  }

  # Load necessary libraries
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    install.packages("ComplexHeatmap")
  }
  if (!requireNamespace("grid", quietly = TRUE)) {
    install.packages("grid")
  }
  if (!requireNamespace("cowplot", quietly = TRUE)) {
    install.packages("cowplot")
  }
  if (!requireNamespace("ggthemes", quietly = TRUE)) {
    install.packages("ggthemes")
  }
  library(ComplexHeatmap)
  library(grid)
  library(cowplot) # For combining and saving plots
  library(ggthemes)

  # Extract the name of the data frame as a string
  data_name <- deparse(substitute(data))

  # Create the var1 by var2 xor_shape_detected matrix
  heatmap_data <- xtabs(xor_shape_detected ~ var1 + var2, data = data)

  # Convert to a matrix
  heatmap_matrix <- as.matrix(heatmap_data)

  # Ensure proper clustering with ward.D2
  row_dend <- hclust(dist(heatmap_matrix), method = clustering_method)
  col_dend <- hclust(dist(t(heatmap_matrix)), method = clustering_method)

  # Define colors for binary values: FALSE -> light gray, TRUE -> ggthemes color scheme
  binary_colors <- c("FALSE" = "lightgray", "TRUE" = ggthemes::colorblind_pal()(8)[2])

  # Ensure `xor_shape_detected` is treated as binary (logical or factor)
  heatmap_matrix <- ifelse(heatmap_matrix > 0, "TRUE", "FALSE")

  # Dynamically calculate label font size based on matrix dimensions
  n_rows <- nrow(heatmap_matrix)
  n_cols <- ncol(heatmap_matrix)

  # Use proportional font sizes (adjust values as needed)
  row_font_size <- max(6, min(12, 120 / n_rows)) # Between 6 and 12, scales with rows
  col_font_size <- max(6, min(12, 120 / n_cols)) # Between 6 and 12, scales with columns

  # Dynamically calculate the plot size based on matrix dimensions
  base_size <- 6 # Base size for smaller matrices
  size_factor <- 0.1 # Smaller scaling factor for better compactness
  plot_width <- base_size + size_factor * n_cols
  plot_height <- base_size + size_factor * n_rows

  # Create the heatmap using ComplexHeatmap
  heatmap_obj <- Heatmap(
    heatmap_matrix,
    name = "Conditions Met", # Legend title
    row_dend_width = unit(4, "cm"), # Set row dendrogram width
    column_dend_height = unit(4, "cm"), # Set column dendrogram height
    cluster_rows = row_dend, # Use ward.D2 clustering for rows
    cluster_columns = col_dend, # Use ward.D2 clustering for columns
    show_row_dend = clustering, # Show row dendrogram only if clustering is enabled
    show_column_dend = clustering, # Show column dendrogram only if clustering is enabled
    col = binary_colors, # Use binary colors
    rect_gp = gpar(col = "white", lwd = 1), # Add white margins between cells
  # Dynamically adjust row and column label sizes
    row_names_gp = gpar(fontsize = row_font_size),
    column_names_gp = gpar(fontsize = col_font_size),
    heatmap_legend_param = list(
      title = "Conditions Met",
      at = c("FALSE", "TRUE"), # Values in legend
      labels = c("False", "True"), # Legend labels
      legend_direction = "horizontal", # Make the legend horizontal
      title_position = "topcenter", # Center the title above the legend
      nrow = 1 # Make the entire legend layout horizontal
    ),
    column_title = data_name # Add data frame name as heatmap title
  )
  
  # Capture the heatmap as a grob object for compatibility with cowplot
  heatmap_plot <- grid::grid.grabExpr(draw(
    heatmap_obj,
    heatmap_legend_side = "bottom", # Move the legend to the bottom
    annotation_legend_side = "bottom" # Ensure annotation legend is also at the bottom
  ))

  # Use cowplot::plot_grid to wrap the heatmap plot
  final_plot <- cowplot::plot_grid(heatmap_plot)
  # print(final_plot)

  # File names for saving
  svg_output_filename <- paste0(data_name, "_results_heatmap.svg")
  png_output_filename <- paste0(data_name, "_results_heatmap.png")

  # Save the plot as SVG using ggsave
  ggsave(
    filename = svg_output_filename,
    plot = final_plot,
    width = plot_width,
    height = plot_height,
    units = "in",
    limitsize = FALSE
  )

  # Save the plot as PNG using ggsave
  ggsave(
    filename = png_output_filename,
    plot = final_plot,
    width = plot_width,
    height = plot_height,
    units = "in",
    dpi = 300, # High resolution for PNG
    limitsize = FALSE
  )

  # Message to confirm saving
  message("Heatmap saved as: ", svg_output_filename)
  message("Heatmap saved as: ", png_output_filename)

  return(final_plot) # Return the final cowplot object
}


# Optional Plot Function for Specific Experiment
generate_champagne_glass__plot <- function(params_data, file_output, title_prefix) {
  plot <- ggplot(
    params_data,
    aes(x = m1 - m2, y = S, color = xor_shape_detected, fill = xor_shape_detected)
  ) +
    geom_point(shape = 21, size = 6, alpha = 0.5) +
    ggvoronoi::geom_voronoi(alpha = 0.3) +
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind() +
    theme_light() +
    labs(title = paste0(title_prefix, " XOR Detections")) +
    theme(
      legend.position.inside = TRUE, legend.position = c(0.2, 0.8),
      legend.background = element_rect(fill = ggplot2::alpha("white", 0.5)),
      strip.background = element_rect(fill = "cornsilk"),
      strip.text = element_text(colour = "black")
    )

  # Save Plot
  ggsave(file_output, plot = plot, width = 16, height = 8)
  print(paste("Plot saved as", file_output))

  return(plot)
}

print("Script XOR_plot_functions_5 loaded")

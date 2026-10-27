# Plot functions for the descriptive figures of the manuscript (Figures 1-4),
# taken from the plot functions used for the original submission.

library(ggplot2)
library(dplyr)

# Scatter plot matrix with Kendall's tau; true XOR pairs highlighted
XOR_matrix_plot <- function(data, Cls, data_name = "dataset", true_pairs_list = NULL, true_pairs_color = "lightyellow") {
  if (is.null(true_pairs_list)) true_pairs_list <- list()
  vars <- colnames(data)
  logical_matrix <- matrix(FALSE, nrow = length(vars), ncol = length(vars), dimnames = list(vars, vars))
  for (pair in true_pairs_list) {
    logical_matrix[pair[1], pair[2]] <- TRUE
    logical_matrix[pair[2], pair[1]] <- TRUE
  }

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

  for (i in seq_len(ncol(data))) {
    for (j in seq_len(ncol(data))) {
      if (i == j) {
        plot[i, j] <- plot[i, j] + theme(panel.background = element_rect(fill = NA, color = NA),
                                         panel.border = element_rect(color = "grey70", fill = NA))
      } else if (i > j) {
        fill <- if (logical_matrix[j, i]) true_pairs_color else NA
        plot[i, j] <- plot[i, j] + theme(panel.background = element_rect(fill = fill, color = NA),
                                         panel.border = element_rect(color = "grey70", fill = NA))
      }
    }
  }
  plot
}

# Density per variable
XOR_density_plot <- function(data, Cls, fill = FALSE) {
  data_long <- reshape2::melt(cbind.data.frame(ID = seq_len(nrow(data)), Cls = factor(Cls), data), id.vars = c("ID", "Cls"))
  ggplot(data_long, aes(x = value, color = as.factor(Cls), fill = as.factor(Cls))) +
    geom_density(alpha = ifelse(fill, 0.1, 0)) +
    facet_wrap(variable ~ ., nrow = 1) +
    theme_light() +
    guides(color = "none", fill = "none") +
    ggthemes::scale_color_colorblind() +
    ggthemes::scale_fill_colorblind() +
    labs(title = "Variable Density") +
    theme(strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "black"))
}

# Violin, box and jitter plots per variable with class comparisons
XOR_raw_plot <- function(data, Cls) {
  data_long <- reshape2::melt(cbind.data.frame(ID = seq_len(nrow(data)), Cls = factor(Cls), data), id.vars = c("ID", "Cls"))
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
    theme(legend.position = "bottom", legend.direction = "horizontal",
          strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "black"))
}

# Stacked bar plots for categorical (0/1/2) variables with chi-square tests
XOR_raw_cat_plot <- function(data, Cls) {
  data_long <- reshape2::melt(cbind.data.frame(ID = seq_len(nrow(data)), Cls = factor(Cls), data), id.vars = c("ID", "Cls"))
  data_agg <- data_long %>%
    filter(value %in% c(0, 1, 2)) %>%
    group_by(variable, Cls, value) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(value = factor(value, levels = c(0, 1, 2)))
  stat_df <- data_long %>%
    filter(value %in% c(0, 1, 2)) %>%
    group_by(variable) %>%
    summarise(pval = chisq.test(table(Cls, value))$p.value)
  ggplot(data_agg, aes(x = Cls, y = count, fill = value)) +
    geom_col(position = "stack", width = 0.7) +
    facet_wrap(~variable, nrow = 1, scales = "free") +
    theme_light() +
    ggthemes::scale_fill_colorblind() +
    labs(title = "Raw Data per Variable", fill = "Value") +
    theme(legend.position = "bottom", legend.direction = "horizontal",
          strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "black")) +
    geom_text(data = stat_df, aes(x = 1.5, y = Inf, label = paste0("Chi-sq p = ", signif(pval, 3))),
              vjust = 1.5, hjust = 0.5, size = 3, inherit.aes = FALSE)
}

# Cases connected across variables
XOR_spaghetti_plot <- function(data, Cls) {
  data_long <- reshape2::melt(cbind.data.frame(ID = seq_len(nrow(data)), Cls = factor(Cls), data), id.vars = c("ID", "Cls"))
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

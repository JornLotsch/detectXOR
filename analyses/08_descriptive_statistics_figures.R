# Conventional statistics and descriptive figures of the introductory data sets
# (manuscript Table 1 and Figures 1-4).
# The ESOM/U-matrix panels of Figure 2 are drawn from the trained maps in
# data/esom; the ESOMs were trained interactively with Umatrix::iEsomTrain
# (50 x 80 planar grid, 20 epochs).
# Output: results/table1_*.csv, figures/Figure1-4.{svg,eps}

source("00_setup.R")
source("descriptive_plot_functions.R")

datasets <- load_example_datasets()
intro <- datasets[c("XOR_data", "XOR_data_nonsig")]

save_figure <- function(plot, name, width, height) {
  ggsave(file.path(FIGURES_DIR, paste0(name, ".svg")), plot, width = width, height = height, limitsize = FALSE)
  ggsave(file.path(FIGURES_DIR, paste0(name, ".eps")), plot, width = width, height = height,
         device = grDevices::cairo_ps, limitsize = FALSE)
}

# Table 1: conventional tests
for (name in names(intro)) {
  data <- intro[[name]]
  vars <- setdiff(names(data), "class")
  tests <- do.call(rbind, lapply(vars, function(v) {
    w <- stats::wilcox.test(data[[v]] ~ factor(data$class))
    t <- stats::t.test(data[[v]] ~ factor(data$class))
    data.frame(variable = v, wilcoxon_W = unname(w$statistic), wilcoxon_p = w$p.value,
               t = unname(t$statistic), t_p = t$p.value)
  }))
  logistic <- summary(stats::glm(factor(class) ~ ., data = data, family = stats::binomial))$coefficients
  long <- reshape2::melt(cbind(data, rowID = seq_len(nrow(data))), id.vars = c("class", "rowID"),
                         variable.name = "variable", value.name = "value")
  long$class <- factor(long$class)
  anova <- suppressMessages(afex::aov_car(value ~ class * variable + Error(rowID / variable), data = long, factorize = FALSE))
  anova_table <- afex::nice(anova, correction = "GG", es = "none")
  write.csv(tests, file.path(RESULTS_DIR, paste0("table1_tests_", name, ".csv")), row.names = FALSE)
  write.csv(data.frame(term = rownames(logistic), logistic, check.names = FALSE),
            file.path(RESULTS_DIR, paste0("table1_logistic_", name, ".csv")), row.names = FALSE)
  write.csv(as.data.frame(anova_table), file.path(RESULTS_DIR, paste0("table1_rmANOVA_", name, ".csv")), row.names = FALSE)
}

# Figure 1: raw data and scatter plot matrices
figure1 <- cowplot::plot_grid(
  XOR_raw_plot(intro$XOR_data[, -1], intro$XOR_data$class),
  XOR_raw_plot(intro$XOR_data_nonsig[, -1], intro$XOR_data_nonsig$class),
  grid::grid.grabExpr(print(XOR_matrix_plot(intro$XOR_data[, -1], factor(intro$XOR_data$class), "XOR_data",
                                            true_pairs_list = list(c("Variable_C", "Variable_D"))))),
  grid::grid.grabExpr(print(XOR_matrix_plot(intro$XOR_data_nonsig[, -1], factor(intro$XOR_data_nonsig$class), "XOR_data_nonsig",
                                            true_pairs_list = list(c("Variable_C", "Variable_D"))))),
  ncol = 2, rel_heights = c(1, 3), align = "v", axis = "lr", labels = "AUTO"
)
save_figure(figure1, "Figure1", width = 24, height = 15)

# Figure 2: projections
projection_plot <- function(coords, class, title) {
  df <- data.frame(x = coords[, 1], y = coords[, 2], group = factor(class))
  p <- ggplot(df) +
    geom_point(aes(x = x, y = y, color = group, shape = group)) +
    stat_ellipse(aes(x = x, y = y, color = group, group = group), level = 0.95, linetype = "dashed", linewidth = 0.8) +
    labs(title = title, x = "Dim1", y = "Dim2", color = "Class") +
    guides(shape = "none") +
    theme_light() +
    ggthemes::scale_color_colorblind() +
    theme(legend.position = "inside", legend.position.inside = c(0.1, 0.1),
          legend.background = element_rect(fill = ggplot2::alpha("white", 0.2))) +
    geom_vline(xintercept = 0, color = "grey20", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "grey20", linetype = "dashed")
  grid::grid.grabExpr(print(ggExtra::ggMarginal(p, type = "densigram", groupColour = TRUE, groupFill = TRUE)))
}
panels <- list()
for (name in names(intro)) {
  data <- intro[[name]]
  x <- data[, -1]
  pca <- mixOmics::pca(X = x, scale = TRUE)$variates$X
  plsda_fit <- mixOmics::plsda(X = x, Y = factor(data$class), scale = TRUE)
  plsda <- stats::predict(plsda_fit, newdata = x)$variates
  set.seed(42)
  umap_layout <- umap::umap(scale(x))$layout
  esom <- file.path(DATA_DIR, "esom", paste0("Umx_", name))
  umx_plot <- Umatrix::plotMatrix(
    Matrix = dbt.DataIO::ReadUMX(FileName = paste0(esom, ".umx")),
    BestMatches = dbt.DataIO::ReadBM(FileName = paste0(esom, ".bm"))$BestMatches,
    Cls = data$class, ClsColors = ggthemes::colorblind_pal()(8), Toroid = FALSE,
    BmSize = 9, RemoveOcean = TRUE, ColorStyle = "Umatrix"
  ) + theme_light() + guides(color = "none", fill = "none") + labs(title = paste0(name, ": U-matrix projection"))
  panels[[name]] <- list(
    PCA = projection_plot(pca, data$class, paste0(name, ": PCA projection")),
    PLSDA = projection_plot(plsda, data$class, paste0(name, ": PLSDA projection")),
    UMAP = projection_plot(umap_layout, data$class, paste0(name, ": UMAP projection")),
    Umatrix = umx_plot
  )
}
figure2 <- cowplot::plot_grid(
  plotlist = unlist(lapply(c("PCA", "PLSDA", "UMAP", "Umatrix"), function(p) list(panels$XOR_data[[p]], panels$XOR_data_nonsig[[p]])), recursive = FALSE),
  labels = "AUTO", align = "hv", axis = "tblr", ncol = 2, rel_heights = c(3, 3, 3, 2)
)
save_figure(figure2, "Figure2", width = 12, height = 24)

# Figure 3: XOR constellation (spaghetti plot, default and rotated XOR)
xor_cd <- data.frame(class = intro$XOR_data$class, scale(intro$XOR_data[, c("Variable_C", "Variable_D")]))
xor_hi <- data.frame(class = datasets$XOR_data_plus_rotated$class, scale(datasets$XOR_data_plus_rotated[, c("Variable_H", "Variable_I")]))
tile_rects <- function(xb, yb, tiles) {
  coords <- list(
    corners = list(c(-Inf, xb[1], yb[2], Inf), c(xb[2], Inf, -Inf, yb[1]), c(-Inf, xb[1], -Inf, yb[1]), c(xb[2], Inf, yb[2], Inf)),
    midfields = list(c(-Inf, xb[1], yb[1], yb[2]), c(xb[2], Inf, yb[1], yb[2]), c(xb[1], xb[2], yb[2], Inf), c(xb[1], xb[2], -Inf, yb[1]))
  )[[tiles]]
  lapply(coords, function(r) annotate("rect", xmin = r[1], xmax = r[2], ymin = r[3], ymax = r[4], fill = "lightyellow", alpha = 0.3))
}
xb <- stats::quantile(xor_cd$Variable_C, c(1/3, 2/3))
yb <- stats::quantile(xor_cd$Variable_D, c(1/3, 2/3))
dot_default <- ggplot(xor_cd, aes(x = Variable_C, y = Variable_D, color = factor(class))) +
  tile_rects(xb, yb, "corners") +
  geom_point() +
  geom_vline(xintercept = xb, linetype = "dashed") +
  geom_hline(yintercept = yb, linetype = "dashed") +
  geom_vline(xintercept = stats::median(xor_cd$Variable_C), color = "salmon", linewidth = 1, linetype = "dashed") +
  geom_hline(yintercept = stats::median(xor_cd$Variable_D), color = "salmon", linewidth = 1, linetype = "dashed") +
  theme_light() + ggthemes::scale_color_colorblind() +
  labs(title = "XOR dotplot", color = "Class", x = "Variable 1", y = "Variable 2") +
  theme(legend.position = "inside", legend.position.inside = c(0.5, 0.5),
        legend.background = element_rect(fill = ggplot2::alpha("white", 0.5)))
xb <- stats::quantile(xor_hi$Variable_H, c(1/3, 2/3))
yb <- stats::quantile(xor_hi$Variable_I, c(1/3, 2/3))
xm <- stats::median(xor_hi$Variable_H)
ym <- stats::median(xor_hi$Variable_I)
dot_rotated <- ggplot(xor_hi, aes(x = Variable_H, y = Variable_I, color = factor(class))) +
  tile_rects(xb, yb, "midfields") +
  geom_point(show.legend = FALSE) +
  geom_vline(xintercept = xb, linetype = "dashed") +
  geom_hline(yintercept = yb, linetype = "dashed") +
  geom_abline(intercept = ym - xm, slope = 1, color = "salmon", linewidth = 1, linetype = "dashed") +
  geom_abline(intercept = ym + xm, slope = -1, color = "salmon", linewidth = 1, linetype = "dashed") +
  theme_light() + ggthemes::scale_color_colorblind() +
  labs(title = "XOR rotated dotplot", x = "Variable 1'", y = "Variable 2'")
figure3 <- cowplot::plot_grid(
  XOR_spaghetti_plot(xor_cd[, -1], xor_cd$class) + labs(title = "XOR variable raw data connected"),
  dot_default, dot_rotated,
  labels = "AUTO", nrow = 1, align = "v", axis = "tb"
)
save_figure(figure3, "Figure3", width = 13, height = 5)

# Figure 4: simulated genetic data
genetic <- datasets$XOR_genetic_data
set.seed(42)
figure4 <- cowplot::plot_grid(
  XOR_raw_cat_plot(genetic[, -1], genetic$class),
  grid::grid.grabExpr(print(XOR_matrix_plot(apply(genetic[, -1], 2, function(x) jitter(x, amount = 0.15)),
                                            factor(genetic$class), "XOR_genetic_data",
                                            true_pairs_list = list(c("Variable_C", "Variable_D"))))),
  ncol = 1, rel_heights = c(1, 3), align = "v", axis = "lr", labels = "AUTO"
)
save_figure(figure4, "Figure4", width = 12, height = 15)

write_session_info("08_descriptive_statistics_figures")

# XOR detection in the example data sets (manuscript Sections 3.4.1-3.4.5 and 4.2,
# Table 3, Figures 5-8).
# Output: results/detection_*.csv, results/detection_summary.csv,
#         results/sliding_means_summary.csv, figures/Figure5-8.{svg,eps}

source("00_setup.R")
library(ggplot2)
library(circlize)

datasets <- load_example_datasets()

run_detection <- function(data, ...) {
  set.seed(42)
  suppressMessages(detect_xor(data, class_col = "class", n_perm = N_PERM, max_cores = N_CORES, ...))
}

# Down-sample class 2 to a fraction of its size (as in the original analysis)
create_imbalanced_data <- function(data, class_column, sampling_fraction) {
  classes <- unique(data[[class_column]])
  subset_class1 <- data[data[[class_column]] == classes[1], ]
  subset_class2 <- data[data[[class_column]] == classes[2], ]
  subset_class2 <- subset_class2[sample(seq_len(nrow(subset_class2)),
                                        size = ceiling(sampling_fraction * nrow(subset_class2))), ]
  imbalanced <- rbind(subset_class1, subset_class2)
  imbalanced[sample(seq_len(nrow(imbalanced))), ]
}

# Data sets for the chord diagram, in display order
screened <- list(
  XOR_data = datasets$XOR_data,
  XOR_data_nonsig = datasets$XOR_data_nonsig
)
for (fraction in c(0.2, 0.4, 0.6, 0.8)) {
  set.seed(42)
  screened[[paste0("XOR_data_plus_rotated_", 100 * fraction)]] <-
    create_imbalanced_data(datasets$XOR_data_plus_rotated, "class", fraction)
}
screened$XOR_data_plus_rotated_100 <- datasets$XOR_data_plus_rotated
screened$pseudoXOR_data_first250 <- datasets$pseudoXOR_data_first250
screened$XOR_genetic_data <- datasets$XOR_genetic_data

# XOR pairs by design, and pairs that form XOR patterns as a consequence of the
# construction: C and F, and D and G, share the same mode assignment per case,
# and the latent variables of the rotated pair H = Hx + Ix, I = Hx - Ix have the
# mode assignment of D (Hx) and C (Ix), so that pairs of C or D with H or I form
# sheared XOR patterns (e.g., C and H correspond to a and a + b)
pair_key <- function(v1, v2) paste(pmin(v1, v2), pmax(v1, v2), sep = "||")
designed_pairs <- list(
  XOR_data = "Variable_C||Variable_D",
  XOR_data_nonsig = "Variable_C||Variable_D",
  plus_rotated = c("Variable_C||Variable_D", "Variable_F||Variable_G", "Variable_H||Variable_I"),
  pseudoXOR_data_first250 = character(0),
  XOR_genetic_data = "Variable_C||Variable_D"
)
implied_pairs <- list(plus_rotated = c("Variable_C||Variable_G", "Variable_D||Variable_F",
                                        "Variable_C||Variable_H", "Variable_C||Variable_I",
                                        "Variable_D||Variable_H", "Variable_D||Variable_I"))
design_group <- function(name) if (grepl("plus_rotated", name)) "plus_rotated" else name

# Run the screening
results <- lapply(names(screened), function(name) {
  message("Screening ", name)
  res <- run_detection(screened[[name]])$results_df
  write.csv(res, file.path(RESULTS_DIR, paste0("detection_", name, ".csv")), row.names = FALSE)
  key <- pair_key(res$var1, res$var2)
  designed <- key %in% designed_pairs[[design_group(name)]]
  implied <- key %in% implied_pairs[[design_group(name)]]
  res$type <- ifelse(designed & res$xor_shape_detected, "Correct, as designed",
              ifelse(designed, "Missing",
              ifelse(implied & res$xor_shape_detected, "Correct, implied by design",
              ifelse(implied, "Missing",
              ifelse(res$xor_shape_detected, "False positive", "No XOR present or detected")))))
  cbind(dataset = name, n_class1 = sum(screened[[name]]$class == 1),
        n_class2 = sum(screened[[name]]$class == 2), res)
})
results <- do.call(rbind, results)
write.csv(results, file.path(RESULTS_DIR, "detection_example_datasets.csv"), row.names = FALSE)

# Numerical summary per data set
summary_table <- do.call(rbind, lapply(split(results, factor(results$dataset, levels = names(screened))), function(d) {
  data.frame(
    dataset = d$dataset[1], n_class1 = d$n_class1[1], n_class2 = d$n_class2[1], pairs = nrow(d),
    expected_xor = sum(d$type %in% c("Correct, as designed", "Correct, implied by design", "Missing")),
    detected_expected = sum(d$type %in% c("Correct, as designed", "Correct, implied by design")),
    missed = sum(d$type == "Missing"),
    false_positive = sum(d$type == "False positive"),
    detected_pairs = paste(sub("Variable_", "", d$var1[d$xor_shape_detected]),
                           sub("Variable_", "", d$var2[d$xor_shape_detected]), sep = "-", collapse = ", ")
  )
}))
write.csv(summary_table, file.path(RESULTS_DIR, "detection_summary.csv"), row.names = FALSE)
print(summary_table)

# Figure 5: chord diagram of the detection results
type_colors <- c("Correct, as designed" = "#009E73", "Correct, implied by design" = "#CC79A7",
                 "Missing" = "#0072B2", "False positive" = "#D55E00",
                 "No XOR present or detected" = "cornsilk")
links <- data.frame(
  from = paste0(results$dataset, "_", results$var1),
  to = paste0(results$dataset, "_", results$var2),
  weight = 10,
  col = type_colors[results$type]
)
all_vars <- unique(c(links$from, links$to))
sector_dataset <- sub("_Variable_.*", "", all_vars)
dataset_colors <- setNames(
  c("chartreuse", "chartreuse3", colorRampPalette(c("dodgerblue", "dodgerblue4"))(5), "salmon", "gold"),
  names(screened)
)

draw_chord_diagram <- function() {
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE,
             canvas.xlim = c(-1.3, 1.3), canvas.ylim = c(-1.3, 1.3))
  par(mar = rep(0, 4))
  chordDiagram(
    x = links[, c("from", "to", "weight")],
    group = setNames(sector_dataset, all_vars),
    grid.col = setNames(dataset_colors[sector_dataset], all_vars),
    col = links$col,
    transparency = 0.1,
    annotationTrack = "grid",
    annotationTrackHeight = c(0.05, 0.1),
    link.sort = TRUE,
    link.largest.ontop = TRUE
  )
  circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    circos.text(x = mean(get.cell.meta.data("xlim")), y = 2.5,
                labels = sub(".*_", "", get.cell.meta.data("sector.index")),
                facing = "bending", cex = 0.7)
  })
  for (ds in names(screened)) {
    ds_vars <- all_vars[sector_dataset == ds]
    mid <- ds_vars[ceiling(length(ds_vars) / 2)]
    circos.text(x = mean(get.cell.meta.data("xlim", sector.index = mid)), y = 6,
                labels = ds, sector.index = mid, facing = "bending", niceFacing = TRUE,
                cex = 0.6, col = dataset_colors[ds])
  }
  legend("topleft", legend = names(type_colors), col = type_colors, lwd = 3, cex = 0.5,
         title = "XOR pair type", box.lty = 0)
  circos.clear()
}
svglite::svglite(file.path(FIGURES_DIR, "Figure5.svg"), width = 7, height = 7)
draw_chord_diagram()
dev.off()
grDevices::cairo_ps(file.path(FIGURES_DIR, "Figure5.eps"), width = 7, height = 7)
draw_chord_diagram()
dev.off()

save_figure <- function(plot, name, width, height) {
  ggsave(file.path(FIGURES_DIR, paste0(name, ".svg")), plot, width = width, height = height)
  ggsave(file.path(FIGURES_DIR, paste0(name, ".eps")), plot, width = width, height = height, device = grDevices::cairo_ps)
}

# Figure 6: detected pairs in XOR_data_plus_rotated
res_rotated <- results[results$dataset == "XOR_data_plus_rotated_100", ]
figure6 <- cowplot::plot_grid(
  generate_spaghetti_plot_from_results(res_rotated, datasets$XOR_data_plus_rotated, "class"),
  generate_xy_plot_from_results(res_rotated, datasets$XOR_data_plus_rotated, "class"),
  ncol = 1, labels = "AUTO"
)
save_figure(figure6, "Figure6", width = 12, height = 14)

# Figure 7: sliding means data, all pairs screened
message("Screening XOR_sliding_means_data (", choose(200, 2), " pairs)")
res_slide <- run_detection(datasets$XOR_sliding_means_data)$results_df
write.csv(res_slide, file.path(RESULTS_DIR, "detection_XOR_sliding_means_data.csv"), row.names = FALSE)
params <- merge(datasets$XOR_sliding_means_params,
                res_slide[, c("var1", "var2", "xor_shape_detected", "xor_purity", "purity_lower", "p_adjusted")],
                by = c("var1", "var2"), all.x = TRUE)
# All pairs of an A and a B variable share the case-wise mode assignment and are
# therefore XOR patterns by construction; pairs of two A or two B variables are not
designed_keys <- pair_key(params$var1, params$var2)
non_designed <- res_slide[!(pair_key(res_slide$var1, res_slide$var2) %in% designed_keys), ]
non_designed$kind <- ifelse(substr(non_designed$var1, 1, 1) != substr(non_designed$var2, 1, 1),
                            "non-designed A-B pairs (XOR by construction)", "A-A and B-B pairs (no XOR)")
write.csv(params, file.path(RESULTS_DIR, "detection_sliding_means_designed_pairs.csv"), row.names = FALSE)

params$sd_level <- cut(params$S, breaks = c(0, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 6),
                       labels = seq(0.5, 5, by = 0.5))
params$mean_difference <- abs(params$m1 - params$m2)
slide_summary <- do.call(rbind, lapply(split(params, params$sd_level), function(d) {
  ci <- wilson_ci(sum(d$xor_shape_detected), nrow(d))
  data.frame(sd = d$sd_level[1], pairs = nrow(d), detected = sum(d$xor_shape_detected),
             rate = mean(d$xor_shape_detected), lower = ci[1], upper = ci[2],
             median_purity = stats::median(d$xor_purity))
}))
slide_summary <- rbind(slide_summary, do.call(rbind, lapply(split(non_designed, non_designed$kind), function(d) {
  ci <- wilson_ci(sum(d$xor_shape_detected), nrow(d))
  data.frame(sd = d$kind[1], pairs = nrow(d), detected = sum(d$xor_shape_detected),
             rate = mean(d$xor_shape_detected), lower = ci[1], upper = ci[2],
             median_purity = stats::median(d$xor_purity))
})))
write.csv(slide_summary, file.path(RESULTS_DIR, "sliding_means_summary.csv"), row.names = FALSE)
print(slide_summary)

params$XOR <- factor(ifelse(params$xor_shape_detected, "detected", "not detected"), levels = c("detected", "not detected"))
figure7 <- ggplot(params, aes(x = m1 - m2, y = S)) +
  ggforce::geom_voronoi_tile(aes(fill = XOR, group = -1L), alpha = 0.3, bound = c(-8, 8, 0, 5.6), colour = "grey50",
                             linewidth = 0.5) +
  geom_point(aes(color = XOR, fill = XOR), shape = 21, size = 4, alpha = 0.6) +
  scale_color_manual(values = c("detected" = "#E69F00", "not detected" = "#000000")) +
  scale_fill_manual(values = c("detected" = "#E69F00", "not detected" = "#000000")) +
  theme_light() +
  labs(x = expression(m[1] - m[2]), y = "Standard deviation S", color = "XOR", fill = "XOR",
       title = "XOR detection in sliding means data") +
  theme(legend.position = "bottom")
save_figure(figure7, "Figure7", width = 8, height = 8)

# Figure 8: incomplete XOR pattern (no detection); variables C and D shown
res_pseudo <- results[results$dataset == "pseudoXOR_data_first250", ]
show_pair <- res_pseudo[res_pseudo$var1 == "Variable_C" & res_pseudo$var2 == "Variable_D", ]
show_pair$xor_shape_detected <- TRUE # flag for display only
figure8 <- cowplot::plot_grid(
  generate_spaghetti_plot_from_results(show_pair, datasets$pseudoXOR_data_first250, "class"),
  generate_xy_plot_from_results(show_pair, datasets$pseudoXOR_data_first250, "class"),
  ncol = 2, labels = "AUTO"
)
save_figure(figure8, "Figure8", width = 12, height = 6)

write_session_info("02_detection_example_datasets")

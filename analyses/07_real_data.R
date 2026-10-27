# Application to real data (manuscript Sections 3.4.9 and 4.6).
#  1. Serum lipid mediators in multiple sclerosis (102 patients, 301 controls,
#     44 lipid mediators; Lötsch et al., Sci Rep 2018). The data are not part of
#     the repository; set the environment variable XOR_LIPIDOMICS_DIR to the
#     directory containing the .lrn and .cls files to include this analysis.
#  2. Pima Indians Diabetes data (R package MASS, 532 complete cases).
# Output: results/real_data_*.csv, results/real_data_summary.csv,
#         figures/Figure12_*.{svg,eps} for data sets with detected XOR pairs

source("00_setup.R")
library(ggplot2)

N_PERM_REAL <- 9999

real_data <- list()

# 1. Lipidomics data
lipid_dir <- Sys.getenv("XOR_LIPIDOMICS_DIR", "/home/joern/Dokumente/BiomarkerLipide")
lrn_file <- file.path(lipid_dir, "09Originale", "DatenAgeSexCorrectedImputed.lrn")
cls_file <- file.path(lipid_dir, "04ESOMUmatrix", "DatenAgeSexCorrectedImputedAllUmxTransformed.cls")
if (file.exists(lrn_file) && file.exists(cls_file)) {
  lrn <- utils::read.table(lrn_file, header = FALSE, comment.char = "%", sep = "\t")
  header <- strsplit(grep("^% Key", readLines(lrn_file), value = TRUE), "\t")[[1]][-1]
  names(lrn) <- c("Key", header)
  cls <- utils::read.table(cls_file, header = FALSE, comment.char = "%", skip = 2, col.names = c("Key", "class"))
  lipids <- merge(cls, lrn, by = "Key")
  # Class 1 = multiple sclerosis, class 2 = healthy controls
  real_data$MS_lipidomics <- lipids[, setdiff(names(lipids), "Key")]
} else {
  message("Lipidomics data not found; skipping. Set XOR_LIPIDOMICS_DIR to include them.")
}

# 2. Pima Indians Diabetes data (complete cases, R package MASS)
pima <- rbind(MASS::Pima.tr, MASS::Pima.te)
real_data$Pima_diabetes <- data.frame(class = ifelse(pima$type == "Yes", 2, 1), pima[, setdiff(names(pima), "type")])

summary_rows <- list()
for (name in names(real_data)) {
  data <- real_data[[name]]
  message("Screening ", name, ": ", nrow(data), " cases, ", ncol(data) - 1, " variables")
  set.seed(42)
  res <- suppressMessages(detect_xor(data, class_col = "class", n_perm = N_PERM_REAL, max_cores = N_CORES))$results_df
  write.csv(res, file.path(RESULTS_DIR, paste0("real_data_", name, ".csv")), row.names = FALSE)

  # Univariate class differences for comparison (Wilcoxon tests, Benjamini-Hochberg)
  vars <- setdiff(names(data), "class")
  p_univariate <- stats::p.adjust(vapply(vars, function(v) stats::wilcox.test(data[[v]] ~ data$class)$p.value, numeric(1)), "BH")

  summary_rows[[name]] <- data.frame(
    dataset = name, n_class1 = sum(data$class == 1), n_class2 = sum(data$class == 2),
    variables = length(vars), pairs = nrow(res),
    univariate_significant = sum(p_univariate < 0.05),
    purity_significant = sum(res$p_adjusted < 0.05),
    tile_criteria_met = sum(res$p_adjusted < 0.05 & res$purity_lower >= 0.8 & res$coverage_ok),
    xor_detected = sum(res$xor_shape_detected),
    max_purity = max(res$xor_purity),
    detected_pairs = paste(res$var1[res$xor_shape_detected], res$var2[res$xor_shape_detected], sep = "-", collapse = ", ")
  )

  if (any(res$xor_shape_detected)) {
    figure <- cowplot::plot_grid(
      generate_spaghetti_plot_from_results(res, data, "class"),
      generate_xy_plot_from_results(res, data, "class"),
      ncol = 1, labels = "AUTO"
    )
    ggsave(file.path(FIGURES_DIR, paste0("Figure12_", name, ".svg")), figure, width = 12, height = 12)
    ggsave(file.path(FIGURES_DIR, paste0("Figure12_", name, ".eps")), figure, width = 12, height = 12, device = grDevices::cairo_ps)
  }
}
summary_table <- do.call(rbind, summary_rows)
write.csv(summary_table, file.path(RESULTS_DIR, "real_data_summary.csv"), row.names = FALSE)
print(summary_table)

# Most XOR-like pairs per data set, for the manuscript text
top <- do.call(rbind, lapply(names(real_data), function(name) {
  res <- utils::read.csv(file.path(RESULTS_DIR, paste0("real_data_", name, ".csv")))
  head(cbind(dataset = name, res[order(-res$xor_purity), c("var1", "var2", "xor_shape_detected", "candidate_pattern",
                                                            "xor_purity", "purity_lower", "coverage_ok", "p_adjusted",
                                                            "tau_class1", "tau_class2")]), 5)
}))
write.csv(top, file.path(RESULTS_DIR, "real_data_top_pairs.csv"), row.names = FALSE)
print(top, digits = 3)

write_session_info("07_real_data")

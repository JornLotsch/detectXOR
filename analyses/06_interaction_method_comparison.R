# Comparison of detectXOR with established interaction-detection approaches
# (manuscript Sections 3.4.8 and 4.5). Each simulated data set contains only the
# target pair (X, Y); a method "flags" the pair at alpha = 0.05.
# Methods:
#   detectXOR            complete XOR decision rule
#   detectXOR_purity     permutation test of the XOR purity only
#   logistic_interaction likelihood-ratio test of the X:Y term in logistic regression
#   interaction_info     interaction information I(X;Y;C) on tertile-discretized
#                        variables, permutation test (synergy > 0)
#   MDR                  multifactor dimensionality reduction on the 3 x 3 tertile
#                        cells, balanced accuracy, permutation test
# Output: results/method_comparison.csv, figures/Figure11.{svg,eps}

source("00_setup.R")
source("simulation_functions.R")
library(ggplot2)

COMPARISON_SAMPLE_SIZES <- c(100, 400)
COMPARISON_DISTRIBUTIONS <- c("normal", "binary")
REPLICATES <- 500

# Class counts per 3 x 3 cell for the observed (column 1) and permuted labelings
cell_counts <- function(tiles, cls, n_perm) {
  is2 <- as.numeric(cls == 2)
  z <- cbind(is2, vapply(seq_len(n_perm), function(b) sample(is2), numeric(length(is2))))
  onehot <- outer(as.integer(tiles), seq_len(9), "==") * 1
  n2 <- crossprod(onehot, z)
  list(n1 = colSums(onehot) - n2, n2 = n2)
}

# Column and row index of the 9 tiles (row-wise from upper left)
TILE_X <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)
TILE_Y <- c(3, 3, 3, 2, 2, 2, 1, 1, 1)

joint_entropy <- function(counts_list) {
  total <- Reduce(`+`, lapply(counts_list, colSums))
  -Reduce(`+`, lapply(counts_list, function(m) {
    p <- sweep(m, 2, total, "/")
    colSums(ifelse(p > 0, p * log(p), 0))
  }))
}

interaction_information <- function(cc) {
  collapse <- function(m, idx) rowsum(m, idx)
  # Only the terms that change under label permutation are needed for the test;
  # the full value is I(XY;C) - I(X;C) - I(Y;C) = -H(XY,C) + H(X,C) + H(Y,C) + const
  h_xyc <- joint_entropy(list(cc$n1, cc$n2))
  h_xc <- joint_entropy(list(collapse(cc$n1, TILE_X), collapse(cc$n2, TILE_X)))
  h_yc <- joint_entropy(list(collapse(cc$n1, TILE_Y), collapse(cc$n2, TILE_Y)))
  -h_xyc + h_xc + h_yc
}

mdr_balanced_accuracy <- function(cc) {
  N1 <- colSums(cc$n1)
  N2 <- colSums(cc$n2)
  high <- sweep(cc$n2, 2, N1, "*") > sweep(cc$n1, 2, N2, "*")
  0.5 * (colSums(cc$n2 * high) / N2 + colSums(cc$n1 * !high) / N1)
}

permutation_p <- function(stat) (1 + sum(stat[-1] >= stat[1] - 1e-12)) / length(stat)

compare_methods <- function(data) {
  res <- suppressWarnings(suppressMessages(detect_xor(data, n_perm = N_PERM, extreme_handling = "none")))$results_df
  d <- data
  d$X <- as.numeric(scale(d$X))
  d$Y <- as.numeric(scale(d$Y))
  d$y01 <- as.numeric(d$class == 2)
  fit0 <- suppressWarnings(stats::glm(y01 ~ X + Y, family = stats::binomial, data = d))
  fit1 <- suppressWarnings(stats::glm(y01 ~ X * Y, family = stats::binomial, data = d))
  p_lr <- stats::pchisq(2 * (as.numeric(stats::logLik(fit1)) - as.numeric(stats::logLik(fit0))), df = 1, lower.tail = FALSE)
  tiles <- detectXOR:::assign_3x3_tiles(data, "X", "Y")
  cc <- cell_counts(tiles, data$class, N_PERM)
  c(detectXOR = res$xor_shape_detected[1],
    detectXOR_purity = res$p_value[1] < 0.05,
    logistic_interaction = p_lr < 0.05,
    interaction_info = permutation_p(interaction_information(cc)) < 0.05,
    MDR = permutation_p(mdr_balanced_accuracy(cc)) < 0.05)
}

conditions <- expand.grid(scenario = SCENARIOS, n = COMPARISON_SAMPLE_SIZES,
                          distribution = COMPARISON_DISTRIBUTIONS, stringsAsFactors = FALSE)
jobs <- expand.grid(condition = seq_len(nrow(conditions)), replicate = seq_len(REPLICATES))
message(nrow(jobs), " data sets")
out <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  cond <- conditions[jobs$condition[j], ]
  set.seed(4e6 + j)
  data <- simulate_dataset(cond$scenario, cond$n, cond$distribution, n_noise = 0)
  compare_methods(data)
}, mc.cores = N_CORES)
failed <- vapply(out, inherits, logical(1), what = "try-error")
if (any(failed)) stop("Method comparison failed: ", as.character(out[[which(failed)[1]]]))
raw <- cbind(conditions[jobs$condition, ], as.data.frame(do.call(rbind, out)))

methods <- c("detectXOR", "detectXOR_purity", "logistic_interaction", "interaction_info", "MDR")
summary_table <- do.call(rbind, lapply(methods, function(m) {
  agg <- stats::aggregate(raw[[m]], by = raw[c("scenario", "n", "distribution")],
                          FUN = function(v) c(x = sum(v), n = length(v)))
  ci <- wilson_ci(agg$x[, "x"], agg$x[, "n"])
  data.frame(method = m, agg[c("scenario", "n", "distribution")], x = agg$x[, "x"], replicates = agg$x[, "n"],
             rate = agg$x[, "x"] / agg$x[, "n"], lower = ci[, "lower"], upper = ci[, "upper"])
}))
summary_table$xor_structure <- summary_table$scenario %in% XOR_SCENARIOS
write.csv(summary_table, file.path(RESULTS_DIR, "method_comparison.csv"), row.names = FALSE)

# Figure 11
scenario_labels <- c(
  null = "Null", xor = "XOR", xor_rotated = "Rotated XOR", xor_noise_10 = "XOR, 10% label noise",
  xor_noise_20 = "XOR, 20% label noise", xor_noise_30 = "XOR, 30% label noise",
  incomplete_xor = "Incomplete XOR", and_interaction = "AND interaction",
  nonmonotone_single = "Non-monotone single variable", linear_diagonal = "Linear diagonal",
  main_effect = "Main effect"
)
method_labels <- c(detectXOR = "detectXOR (complete rule)", detectXOR_purity = "detectXOR purity test only",
                   logistic_interaction = "Logistic regression, interaction term",
                   interaction_info = "Interaction information", MDR = "MDR")
plot_data <- summary_table
plot_data$scenario <- factor(scenario_labels[plot_data$scenario], levels = rev(scenario_labels))
plot_data$method <- factor(method_labels[plot_data$method], levels = method_labels)
plot_data$condition <- paste0(plot_data$distribution, ", n = ", plot_data$n)
figure11 <- ggplot(plot_data, aes(x = rate, y = scenario, color = method)) +
  geom_pointrange(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.7), size = 0.25) +
  facet_wrap(~condition, nrow = 1) +
  scale_color_manual(values = c("#E69F00", "#F0E442", "#0072B2", "#009E73", "#CC79A7")) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_light() +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "cornsilk"),
        strip.text = element_text(colour = "black")) +
  guides(color = guide_legend(nrow = 2)) +
  labs(x = "Rate at which the variable pair is flagged (alpha = 0.05)", y = NULL, color = "Method")
ggsave(file.path(FIGURES_DIR, "Figure11.svg"), figure11, width = 14, height = 7)
ggsave(file.path(FIGURES_DIR, "Figure11.eps"), figure11, width = 14, height = 7, device = grDevices::cairo_ps)

write_session_info("06_interaction_method_comparison")

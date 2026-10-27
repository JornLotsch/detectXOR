# Simulation study of the XOR detection rule (manuscript Sections 3.4.6 and 4.3).
# Study A: family-wise false-positive rate under the null hypothesis.
# Study B: detection rates of the target pair for XOR and non-XOR structures,
#          and false detections among the noise pairs; for imbalanced classes
#          also with class balancing by down-sampling (balance_classes = TRUE).
# Output: results/simulation_null.csv, results/simulation_detection.csv,
#         results/simulation_raw_*.rds, figures/Figure9.{svg,eps}

source("00_setup.R")
source("simulation_functions.R")
library(ggplot2)

SAMPLE_SIZES <- c(50, 100, 200, 400, 1000)

run_replicate <- function(condition, seed) {
  set.seed(seed)
  data <- simulate_dataset(condition$scenario, condition$n, condition$distribution, condition$prop_class1)
  res <- suppressWarnings(suppressMessages(
    detect_xor(data, n_perm = N_PERM, balance_classes = condition$balance_classes,
               extreme_handling = if (condition$distribution == "binary") "none" else "winsorize")
  ))$results_df
  summarize_detection(res)
}

run_conditions <- function(conditions, replicates, seed_offset) {
  jobs <- expand.grid(condition = seq_len(nrow(conditions)), replicate = seq_len(replicates))
  out <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
    run_replicate(conditions[jobs$condition[j], ], seed = seed_offset + j)
  }, mc.cores = N_CORES)
  failed <- vapply(out, inherits, logical(1), what = "try-error")
  if (any(failed)) stop("Simulation failed: ", as.character(out[[which(failed)[1]]]))
  cbind(conditions[jobs$condition, ], replicate = jobs$replicate, as.data.frame(do.call(rbind, out)))
}

aggregate_rates <- function(raw, outcome) {
  keys <- setdiff(names(raw), c("replicate", "target_detected", "target_significant",
                                "other_detected", "any_detected", "any_significant"))
  agg <- stats::aggregate(raw[[outcome]], by = raw[keys], FUN = function(v) c(x = sum(v), n = length(v)))
  agg <- cbind(agg[keys], x = agg$x[, "x"], replicates = agg$x[, "n"])
  ci <- wilson_ci(agg$x, agg$replicates)
  cbind(agg, rate = agg$x / agg$replicates, lower = ci[, "lower"], upper = ci[, "upper"])
}

# Study A: null hypothesis, 10 variables (45 pairs) unrelated to the class
null_conditions <- expand.grid(scenario = "null", n = SAMPLE_SIZES, distribution = DISTRIBUTIONS,
                               prop_class1 = 0.5, balance_classes = FALSE, stringsAsFactors = FALSE)
message("Study A: ", nrow(null_conditions), " conditions x 500 replicates")
raw_null <- run_conditions(null_conditions, replicates = 500, seed_offset = 1e6)
saveRDS(raw_null, file.path(RESULTS_DIR, "simulation_raw_null.rds"))
null_summary <- merge(
  aggregate_rates(raw_null, "any_detected"),
  aggregate_rates(raw_null, "any_significant"),
  by = c("scenario", "n", "distribution", "prop_class1", "balance_classes"), suffixes = c("_xor_rule", "_permutation_test")
)
write.csv(null_summary, file.path(RESULTS_DIR, "simulation_null.csv"), row.names = FALSE)
print(null_summary[, c("n", "distribution", "rate_xor_rule", "rate_permutation_test", "upper_permutation_test")])

# Study B: target pair with XOR or non-XOR structure plus 8 noise variables
detection_conditions <- rbind(
  expand.grid(scenario = setdiff(SCENARIOS, "null"), n = SAMPLE_SIZES, distribution = DISTRIBUTIONS,
              prop_class1 = c(0.5, 0.25), balance_classes = FALSE, stringsAsFactors = FALSE),
  expand.grid(scenario = setdiff(SCENARIOS, "null"), n = SAMPLE_SIZES, distribution = DISTRIBUTIONS,
              prop_class1 = 0.25, balance_classes = TRUE, stringsAsFactors = FALSE)
)
message("Study B: ", nrow(detection_conditions), " conditions x 200 replicates")
raw_detection <- run_conditions(detection_conditions, replicates = 200, seed_offset = 2e6)
saveRDS(raw_detection, file.path(RESULTS_DIR, "simulation_raw_detection.rds"))
detection_summary <- merge(
  aggregate_rates(raw_detection, "target_detected"),
  aggregate_rates(raw_detection, "other_detected"),
  by = c("scenario", "n", "distribution", "prop_class1", "balance_classes"), suffixes = c("_target", "_noise_pairs")
)
detection_summary$xor_structure <- detection_summary$scenario %in% XOR_SCENARIOS
write.csv(detection_summary, file.path(RESULTS_DIR, "simulation_detection.csv"), row.names = FALSE)

# Figure 9: detection rates of the target pair
scenario_labels <- c(
  xor = "XOR", xor_rotated = "Rotated XOR", xor_noise_10 = "XOR, 10% label noise",
  xor_noise_20 = "XOR, 20% label noise", xor_noise_30 = "XOR, 30% label noise",
  incomplete_xor = "Incomplete XOR", and_interaction = "AND interaction",
  nonmonotone_single = "Non-monotone single variable", linear_diagonal = "Linear diagonal",
  main_effect = "Main effect"
)
plot_data <- detection_summary
plot_data$scenario <- factor(scenario_labels[plot_data$scenario], levels = scenario_labels)
plot_data$classes <- factor(ifelse(plot_data$prop_class1 == 0.5, "balanced (1:1)",
                                   ifelse(plot_data$balance_classes, "imbalanced (1:3), down-sampled", "imbalanced (1:3)")),
                            levels = c("balanced (1:1)", "imbalanced (1:3)", "imbalanced (1:3), down-sampled"))
figure9 <- ggplot(plot_data, aes(x = n, y = rate_target, color = distribution, linetype = classes)) +
  geom_ribbon(aes(ymin = lower_target, ymax = upper_target, fill = distribution), alpha = 0.08, color = NA) +
  geom_line() +
  geom_point(size = 1) +
  scale_x_log10(breaks = SAMPLE_SIZES) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~scenario, ncol = 5) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme_light() +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "cornsilk"),
        strip.text = element_text(colour = "black")) +
  labs(x = "Sample size n", y = "Rate of XOR detection of the target pair",
       color = "Distribution", fill = "Distribution", linetype = "Classes")
ggsave(file.path(FIGURES_DIR, "Figure9.svg"), figure9, width = 14, height = 7)
ggsave(file.path(FIGURES_DIR, "Figure9.eps"), figure9, width = 14, height = 7, device = grDevices::cairo_ps)

write_session_info("04_simulation_study")

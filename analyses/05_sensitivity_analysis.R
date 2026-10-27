# Sensitivity of the detection results to the parameters of the decision rule
# (manuscript Sections 3.4.7 and 4.4). One parameter is varied at a time while
# the others are kept at their defaults; all settings are applied to the same
# simulated data sets (common random numbers).
# Output: results/sensitivity_analysis.csv, figures/Figure10.{svg,eps}

source("00_setup.R")
source("simulation_functions.R")
library(ggplot2)

SENSITIVITY_SCENARIOS <- c("null", "xor", "xor_rotated", "xor_noise_20", "incomplete_xor",
                           "and_interaction", "nonmonotone_single")
N_SAMPLE <- 400
REPLICATES <- 200

defaults <- list(split_quantiles = c(1/3, 2/3), check_tau = TRUE, tau_threshold = 0.3,
                 purity_threshold = 0.8, min_tile_share = 0.05)
settings <- c(
  lapply(list(c(0.25, 0.75), c(1/3, 2/3), c(0.4, 0.6)), function(v)
    list(parameter = "Tile boundaries (quantiles)", value = paste(round(v, 2), collapse = "/"),
         args = utils::modifyList(defaults, list(split_quantiles = v)))),
  list(list(parameter = "Kendall's tau threshold", value = "no tau criterion",
            args = utils::modifyList(defaults, list(check_tau = FALSE)))),
  lapply(c(0.1, 0.2, 0.3, 0.4, 0.5), function(v)
    list(parameter = "Kendall's tau threshold", value = as.character(v),
         args = utils::modifyList(defaults, list(tau_threshold = v)))),
  lapply(c(0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95), function(v)
    list(parameter = "Minimum purity s0", value = as.character(v),
         args = utils::modifyList(defaults, list(purity_threshold = v)))),
  lapply(c(0, 0.025, 0.05, 0.1, 0.2, 0.3), function(v)
    list(parameter = "Minimum tile share rho", value = as.character(v),
         args = utils::modifyList(defaults, list(min_tile_share = v))))
)

jobs <- expand.grid(scenario = SENSITIVITY_SCENARIOS, replicate = seq_len(REPLICATES), stringsAsFactors = FALSE)
message(nrow(jobs), " data sets x ", length(settings), " settings")
out <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  # Same seed for a given scenario and replicate in all settings
  seed <- 3e6 + 1000 * match(jobs$scenario[j], SENSITIVITY_SCENARIOS) + jobs$replicate[j]
  set.seed(seed)
  data <- simulate_dataset(jobs$scenario[j], N_SAMPLE, "normal")
  do.call(rbind, lapply(settings, function(s) {
    set.seed(seed)
    res <- suppressWarnings(suppressMessages(do.call(detect_xor, c(list(data = data, n_perm = N_PERM), s$args))))$results_df
    det <- summarize_detection(res)
    data.frame(scenario = jobs$scenario[j], replicate = jobs$replicate[j], parameter = s$parameter,
               value = s$value, target_detected = det[["target_detected"]], any_detected = det[["any_detected"]])
  }))
}, mc.cores = N_CORES)
failed <- vapply(out, inherits, logical(1), what = "try-error")
if (any(failed)) stop("Sensitivity analysis failed: ", as.character(out[[which(failed)[1]]]))
raw <- do.call(rbind, out)

# Outcome: detection of the target pair; for null data, any detection among the 45 pairs
raw$outcome <- ifelse(raw$scenario == "null", raw$any_detected, raw$target_detected)
agg <- stats::aggregate(outcome ~ scenario + parameter + value, data = raw, FUN = function(v) c(x = sum(v), n = length(v)))
agg <- cbind(agg[c("scenario", "parameter", "value")], x = agg$outcome[, "x"], n = agg$outcome[, "n"])
ci <- wilson_ci(agg$x, agg$n)
agg <- cbind(agg, rate = agg$x / agg$n, lower = ci[, "lower"], upper = ci[, "upper"])
write.csv(agg, file.path(RESULTS_DIR, "sensitivity_analysis.csv"), row.names = FALSE)

# Figure 10
scenario_labels <- c(null = "Null (any of 45 pairs)", xor = "XOR", xor_rotated = "Rotated XOR",
                     xor_noise_20 = "XOR, 20% label noise", incomplete_xor = "Incomplete XOR",
                     and_interaction = "AND interaction", nonmonotone_single = "Non-monotone single variable")
# Unique axis positions per parameter; default values are labeled as such
default_labels <- c("Tile boundaries (quantiles)" = "0.33/0.67", "Kendall's tau threshold" = "0.3",
                    "Minimum purity s0" = "0.8", "Minimum tile share rho" = "0.05")
setting_keys <- vapply(settings, function(s) paste(s$parameter, s$value, sep = "|"), character(1))
agg$x <- factor(paste(agg$parameter, agg$value, sep = "|"), levels = setting_keys)
agg$scenario <- factor(scenario_labels[agg$scenario], levels = scenario_labels)
axis_label <- function(keys) {
  parameter <- sub("\\|.*", "", keys)
  value <- sub(".*\\|", "", keys)
  ifelse(value == default_labels[parameter], paste0(value, " (default)"), value)
}
figure10 <- ggplot(agg, aes(x = x, y = rate, color = scenario, group = scenario)) +
  geom_line() +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.2) +
  facet_wrap(~parameter, scales = "free_x", nrow = 1) +
  scale_x_discrete(labels = axis_label) +
  scale_color_manual(values = c("#000000", "#009E73", "#56B4E9", "#0072B2", "#E69F00", "#D55E00", "#CC79A7")) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_light() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "cornsilk"), strip.text = element_text(colour = "black")) +
  labs(x = "Parameter value", y = "Rate of XOR detection", color = "Scenario")
ggsave(file.path(FIGURES_DIR, "Figure10.svg"), figure10, width = 14, height = 6)
ggsave(file.path(FIGURES_DIR, "Figure10.eps"), figure10, width = 14, height = 6, device = grDevices::cairo_ps)

write_session_info("05_sensitivity_analysis")

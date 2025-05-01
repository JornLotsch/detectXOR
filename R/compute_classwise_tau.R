#' Compute classwise Kendall tau correlations for all pairwise datasets
#'
#' This function iterates over a list of pairwise datasets (as produced by create_pairwise_datasets)
#' and fills the results$classwise_tau entry for each pair using classwise Kendall tau analysis.
#'
#' @param orig_pair_list List. Output from create_pairwise_datasets.
#' @param class_col Character. Name of the class column (default: "class").
#' @param tau_threshold Numeric. Minimum absolute tau for "strong" correlation (default: 0.3).
#' @return The input list, but with results$classwise_tau filled for each pair.
#' @importFrom stats cor.test
#' @export
compute_classwise_tau_for_pairs <- function(
    orig_pair_list,
    class_col = "class",
    tau_threshold = 0.3
) {
  for (nm in names(orig_pair_list)) {
    pair <- orig_pair_list[[nm]]
    df <- pair$data
    var_names <- setdiff(names(df), class_col)
    if (length(var_names) != 2) stop("Each pairwise data frame must have exactly 2 variables + class column.")
    classwise_tau_result <- check_kendall_correlation_pair(
      data = df,
      var1_name = var_names[1],
      var2_name = var_names[2],
      class_col = class_col,
      tau_threshold = tau_threshold
    )
    orig_pair_list[[nm]]$results$classwise_tau <- classwise_tau_result
  }
  orig_pair_list
}

#' Main function for XOR Cartesian reversed tau correlation analysis
#' @keywords internal
check_kendall_correlation_pair <- function(data, var1_name, var2_name, class_col,
                                           tau_threshold = 0.3) {
  stopifnot(all(c(var1_name, var2_name, class_col) %in% colnames(data)))
  cls <- unique(data[[class_col]])
  stopifnot(length(cls) == 2)

  test_correlation <- function(df) {
    class1_data <- df[df$class == cls[1],]
    class2_data <- df[df$class == cls[2],]
    tau1_test <- cor.test(class1_data$var1, class1_data$var2, method = "kendall")
    tau2_test <- cor.test(class2_data$var1, class2_data$var2, method = "kendall")
    tau1 <- as.numeric(tau1_test$estimate)
    tau2 <- as.numeric(tau2_test$estimate)
    is_reverse <- all(abs(tau1) >= tau_threshold, abs(tau2) >= tau_threshold, sign(tau1) != sign(tau2))
    list(
      tau_class1 = tau1,
      tau_class2 = tau2,
      is_strong_reverse_correlated = is_reverse,
      tau_class1_p_value = tau1_test$p.value,
      tau_class2_p_value = tau2_test$p.value
    )
  }

  orig_data <- data.frame(
    class = data[[class_col]],
    var1 = data[[var1_name]],
    var2 = data[[var2_name]]
  )
  rotated_data <- data.frame(
    class = data[[class_col]],
    var1 = data[[var1_name]] + data[[var2_name]],
    var2 = data[[var1_name]] - data[[var2_name]]
  )
  classwise_tau_orig <- test_correlation(orig_data)
  classwise_tau_rotated <- test_correlation(rotated_data)
  list(
    default = classwise_tau_orig,
    rotated = classwise_tau_rotated
  )
}


# pairwise_dfs <- compute_classwise_tau_for_pairs(pairwise_dfs)

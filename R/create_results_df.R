#' Combine the XOR criteria of all pairs into the final decision and summary table
#' @param orig_pair_list List of analyzed pairs.
#' @param settings List of detect_xor parameter values.
#' @return Data frame with one row per variable pair.
#' @noRd
create_results_df <- function(orig_pair_list, settings) {
  if (is.null(orig_pair_list) || length(orig_pair_list) == 0) return(NULL)

  # Only keep pairs with tile_pattern (minimal requirement)
  valid_pairs <- Filter(function(x) {
    is.list(x) && !is.null(x$results$tile_pattern)
  }, orig_pair_list)

  if (length(valid_pairs) == 0) return(NULL)

  results_df <- do.call(rbind, lapply(valid_pairs, function(pair) {
    tile_pat <- pair$results$tile_pattern
    var_names <- setdiff(names(pair$data), settings$class_col)

    # Class-wise tau in the coordinates of the chosen pattern
    taus <- pair$results$classwise_tau[[tile_pat$xor_pattern]]
    tau_ok <- if (is.null(taus)) NA else taus$is_strong_reverse_correlated

    # Axis-wise Wilcoxon tests (NA if not computed or not applicable)
    wilcox <- pair$results$tile_wilcox_significance
    axis_ok <- if (is.null(wilcox)) NA else wilcox$pass_groups

    # Collect unmet criteria
    reasons <- c(
      if (tile_pat$p_adjusted >= settings$p_threshold) "purity not above chance",
      if (tile_pat$purity_lower < settings$purity_threshold) "purity below threshold",
      if (!tile_pat$coverage_ok) "incomplete pattern (expected tile empty or sparse)",
      if (settings$check_tau && tile_pat$candidate && !isTRUE(tau_ok)) "class-wise tau criterion not met",
      if (settings$compute_axes_parallel_significance && tile_pat$candidate && identical(axis_ok, FALSE))
        "group-wise Wilcoxon tests not all significant"
    )
    detected <- length(reasons) == 0

    interpretation <- paste0(
      tile_pat$xor_pattern, " pattern: purity = ", signif(tile_pat$purity, 3),
      " (lower bound ", signif(tile_pat$purity_lower, 3), "), p = ", signif(tile_pat$p_value, 3),
      ", adjusted p = ", signif(tile_pat$p_adjusted, 3), ". ",
      if (detected) "XOR pattern detected." else paste0("No XOR: ", paste(reasons, collapse = "; "), ".")
    )

    data.frame(
      var1 = var_names[1],
      var2 = var_names[2],
      xor_shape_detected = detected,
      xor_pattern = if (detected) tile_pat$xor_pattern else "none",
      candidate_pattern = tile_pat$xor_pattern,
      xor_purity = tile_pat$purity,
      purity_lower = tile_pat$purity_lower,
      coverage_ok = tile_pat$coverage_ok,
      p_value = tile_pat$p_value,
      p_adjusted = tile_pat$p_adjusted,
      tau_class1 = if (is.null(taus)) NA else taus$tau_class1,
      tau_class1_p_value = if (is.null(taus)) NA else taus$tau_class1_p_value,
      tau_class2 = if (is.null(taus)) NA else taus$tau_class2,
      tau_class2_p_value = if (is.null(taus)) NA else taus$tau_class2_p_value,
      is_strong_reverse_correlated = tau_ok,
      axis_tests_passed = axis_ok,
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  }))

  results_df
}

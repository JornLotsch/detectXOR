create_results_df_from_orig_pair_list <- function(orig_pair_list) {
  if (is.null(orig_pair_list) || length(orig_pair_list) == 0) return(NULL)

  # Only keep pairs with tile_pattern (minimal requirement)
  valid_pairs <- Filter(function(x) {
    is.list(x) && !is.null(x$results$tile_pattern)
  }, orig_pair_list)

  if (length(valid_pairs) == 0) return(NULL)

  results_df <- do.call(rbind, lapply(valid_pairs, function(pair) {
    tile_pat <- pair$results$tile_pattern
    taus <- pair$results$classwise_tau

    var_names <- setdiff(names(pair$data), "class")

    # Defensive: If classwise_tau is missing or incomplete, fill with NA
    pattern_exists <- !is.null(taus) && !is.null(tile_pat$xor_pattern) &&
      !is.null(taus[[tile_pat$xor_pattern]])
    is_strong_reverse <- if (pattern_exists) taus[[tile_pat$xor_pattern]]$is_strong_reverse_correlated else NA

    xor_override <- pattern_exists && !is_strong_reverse

    base_interpretation <- tile_pat$interpretation
    should_add_warning <- grepl("No significant or practical difference detected", base_interpretation) &&
      !is.na(is_strong_reverse) &&
      !is_strong_reverse

    if (should_add_warning) {
      base_interpretation <- paste0(
        base_interpretation,
        " However, class-wise correlation does not meet XOR expectations (XOR pattern not detected)."
      )
    }

    final_xor_detected <- if (should_add_warning) FALSE else tile_pat$xor_shape_detected

    data.frame(
      var1 = var_names[1],
      var2 = var_names[2],
      xor_shape_detected = final_xor_detected,
      xor_pattern = if (!is.null(final_xor_detected) && final_xor_detected) tile_pat$xor_pattern else "none",
      chi_sq = tile_pat$chi_squared_test$statistic,
      chi_sq_p_value = tile_pat$chi_squared_test$p_value,
      tau_class1 = if (pattern_exists) taus[[tile_pat$xor_pattern]]$tau_class1 else NA,
      tau_class1_p_value = if (pattern_exists) taus[[tile_pat$xor_pattern]]$tau_class1_p_value else NA,
      tau_class2 = if (pattern_exists) taus[[tile_pat$xor_pattern]]$tau_class2 else NA,
      tau_class2_p_value = if (pattern_exists) taus[[tile_pat$xor_pattern]]$tau_class2_p_value else NA,
      is_strong_reverse_correlated = if (pattern_exists) is_strong_reverse else NA,
      interpretation = base_interpretation,
      stringsAsFactors = FALSE
    )
  }))

  results_df
}


# results_df <- create_results_df_from_orig_pair_list(pairwise_dfs)


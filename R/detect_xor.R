#' @title Detect XOR-like patterns in pairwise variable relationships
#' @description Orchestrates the full XOR detection pipeline: tile pattern analysis, classwise Kendall tau correlation, and group-wise Wilcoxon significance tests.
#' @param data Data frame containing variables and class column.
#' @param class_col Character. Name of the class column (default: "class").
#' @param check_tau Logical. Whether to compute classwise Kendall tau correlations (default: TRUE).
#' @param compute_axes_parallel_significance Logical. Whether to compute group-wise Wilcoxon tests (default: TRUE).
#' @param p_threshold Numeric. Significance threshold for statistical tests (default: 0.05).
#' @param tau_threshold Numeric. Minimum absolute tau for "strong" correlation (default: 0.3).
#' @param abs_diff_threshold Numeric. Minimum absolute difference for practical significance (default: 20).
#' @param split_method Character. Method for splitting tiles ("quantile" or "range") (default: "quantile").
#' @return Data frame summarizing XOR detection results for all variable pairs.
#' @importFrom stats wilcox.test
#' @importFrom stats cor.test
#' @importFrom stats chisq.test
#' @export
#' @examples
#' \dontrun{
#' results <- detect_xor(
#'   data = my_data,
#'   class_col = "class",
#'   p_threshold = 0.05,
#'   tau_threshold = 0.3
#' )
#' }
detect_xor <- function(
  data,
  class_col = "class",
  check_tau = TRUE,
  compute_axes_parallel_significance = TRUE,
  p_threshold = 0.05,
  tau_threshold = 0.3,
  abs_diff_threshold = 20,
  split_method = "quantile",
  max_cores = NULL,
  extreme_handling = "winsorize",
  winsor_limits = c(0.05, 0.95),
  scale_data = TRUE,
  use_complete = TRUE
) {
  # Determine number of cores to use
  n_cores <- determine_n_cores(max_cores)
  use_parallel <- n_cores > 1

  # Check required packages for parallel processing
  if (use_parallel && !check_parallel_packages()) {
    warning("Parallel processing packages not available. Using sequential processing.")
    use_parallel <- FALSE
    n_cores <- 1
  }

  message("\nStep 1: Creating pairwise datasets...")
  orig_pair_list <- create_pairwise_datasets(
    data,
    class_col = class_col,
    extreme_handling = extreme_handling,
    winsor_limits = winsor_limits,
    scale_data = scale_data,
    use_complete = use_complete
  )

  message("\nStep 2: Computing tile patterns...")
  if (use_parallel) {
    chunks <- create_balanced_chunks(length(orig_pair_list), n_cores)

    processed_chunks <- safe_parallel_execution(
      chunks,
      function(chunk_indices) {
        chunk_data <- orig_pair_list[chunk_indices]
        compute_tile_patterns_for_pairs(
          chunk_data,
          class_col = class_col,
          p_threshold = p_threshold,
          abs_diff_threshold = abs_diff_threshold,
          split_method = split_method
        )
      },
      n_cores = n_cores
    )

    orig_pair_list <- processed_chunks
    names(orig_pair_list) <- names(create_pairwise_datasets(
      data,
      class_col = class_col,
      extreme_handling = extreme_handling,
      winsor_limits = winsor_limits,
      scale_data = scale_data,
      use_complete = use_complete
    ))
  } else {
    orig_pair_list <- compute_tile_patterns_for_pairs(
      orig_pair_list,
      class_col = class_col,
      p_threshold = p_threshold,
      abs_diff_threshold = abs_diff_threshold,
      split_method = split_method
    )
  }

  # Validate Step 2 output
  if (is.null(orig_pair_list[[1]]$results$tile_pattern)) {
    stop("Tile pattern computation failed. Check compute_tile_patterns_for_pairs()")
  }

  # Step 3: Classwise tau with safe indexing
  if (check_tau) {
    xor_pairs <- names(orig_pair_list)[vapply(orig_pair_list, function(x) {
      isTRUE(x$results$tile_pattern$xor_shape_detected)
    }, logical(1))]

    if (length(xor_pairs) > 0) {
      message("\nStep 3: Computing classwise tau for ", length(xor_pairs), " XOR pairs...")
      if (use_parallel) {
        chunks <- create_balanced_chunks(length(xor_pairs), n_cores)

        processed_tau_chunks <- safe_parallel_execution(
          chunks,
          function(chunk_indices) {
            chunk_pairs <- xor_pairs[chunk_indices]
            chunk_data <- orig_pair_list[chunk_pairs]
            compute_classwise_tau_for_pairs(
              chunk_data,
              class_col = class_col,
              tau_threshold = tau_threshold
            )
          },
          n_cores = n_cores
        )

        orig_pair_list[xor_pairs] <- processed_tau_chunks
      } else {
        orig_pair_list[xor_pairs] <- compute_classwise_tau_for_pairs(
          orig_pair_list[xor_pairs],
          class_col = class_col,
          tau_threshold = tau_threshold
        )
      }
    }
  }

  # Step 4: Wilcoxon tests
  if (compute_axes_parallel_significance) {
    message("\nStep 4: Computing Wilcoxon significance tests...")
    if (use_parallel) {
      chunks <- create_balanced_chunks(length(orig_pair_list), n_cores)

      processed_wilcox_chunks <- safe_parallel_execution(
        chunks,
        function(chunk_indices) {
          chunk_data <- orig_pair_list[chunk_indices]
          compute_tile_wilcox_significance_for_pairs(
            chunk_data,
            class_col = class_col,
            p_threshold = p_threshold,
            split_method = split_method
          )
        },
        n_cores = n_cores
      )

      orig_pair_list <- processed_wilcox_chunks
    } else {
      orig_pair_list <- compute_tile_wilcox_significance_for_pairs(
        orig_pair_list,
        class_col = class_col,
        p_threshold = p_threshold,
        split_method = split_method
      )
    }
  }

  message("\nFinalizing results...")
  if (!exists("create_results_df_from_orig_pair_list")) {
    stop("Missing results formatting function")
  }

  # Return both the detailed pair list and the results dataframe
  list(
    results_df = create_results_df_from_orig_pair_list(orig_pair_list),
    pair_list = orig_pair_list
  )
}
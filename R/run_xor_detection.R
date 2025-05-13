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
#' results <- run_xor_detection(
#'   data = my_data,
#'   class_col = "class",
#'   p_threshold = 0.05,
#'   tau_threshold = 0.3
#' )
#' }
# run_xor_detection <- function(
#     data,
#     class_col = "class",
#     check_tau = TRUE,
#     compute_axes_parallel_significance = TRUE,
#     p_threshold = 0.05,
#     tau_threshold = 0.3,
#     abs_diff_threshold = 20,
#     split_method = "quantile"
# ) {
#   # Step 1: Create pairwise datasets
#   orig_pair_list <- create_pairwise_datasets(data, class_col = class_col)
#
#   # Step 2: Compute tile patterns with validation
#   orig_pair_list <- compute_tile_patterns_for_pairs(
#     orig_pair_list,
#     class_col = class_col,
#     p_threshold = p_threshold,
#     abs_diff_threshold = abs_diff_threshold,
#     split_method = split_method
#   )
#
#   # Validate Step 2 output
#   if (is.null(orig_pair_list[[1]]$results$tile_pattern)) {
#     stop("Tile pattern computation failed. Check compute_tile_patterns_for_pairs()")
#   }
#
#   # Step 3: Classwise tau with safe indexing
#   if (check_tau) {
#     xor_pairs <- names(orig_pair_list)[vapply(orig_pair_list, function(x) {
#       isTRUE(x$results$tile_pattern$xor_shape_detected)
#     }, logical(1))]
#
#     if (length(xor_pairs) > 0) {
#       orig_pair_list[xor_pairs] <- compute_classwise_tau_for_pairs(
#         orig_pair_list[xor_pairs],
#         class_col = class_col,
#         tau_threshold = tau_threshold
#       )
#     }
#   }
#
#   # Step 4: Wilcoxon tests
#   if (compute_axes_parallel_significance) {
#     orig_pair_list <- compute_tile_wilcox_significance_for_pairs(
#       orig_pair_list,
#       class_col = class_col,
#       p_threshold = p_threshold,
#       split_method = split_method
#     )
#   }
#
#   # Final output validation
#   if (!exists("create_results_df_from_orig_pair_list")) {
#     stop("Missing results formatting function")
#   }
#   create_results_df_from_orig_pair_list(orig_pair_list)
# }
#
# # Note: The following helper functions must be in your package's namespace:
# # - create_pairwise_datasets()
# # - compute_tile_patterns_for_pairs()
# # - compute_classwise_tau_for_pairs()
# # - compute_tile_wilcox_significance_for_pairs()
# # - create_results_df_from_orig_pair_list()
#
# run_xor_detection_parallel <- function(
#   data,
#   class_col = "class",
#   check_tau = TRUE,
#   compute_axes_parallel_significance = TRUE,
#   p_threshold = 0.05,
#   tau_threshold = 0.3,
#   abs_diff_threshold = 20,
#   split_method = "quantile",
#   n_cores = parallel::detectCores() - 1
# ) {
#   # Required packages
#   require(future)
#   require(future.apply)
#   require(pbmcapply)
#
#   message("\nStep 1: Creating pairwise datasets...")
#   orig_pair_list <- create_pairwise_datasets(data, class_col = class_col)
#
#   is_windows <- .Platform$OS.type == "windows"
#   chunk_size <- ceiling(length(orig_pair_list) / n_cores)
#   chunks <- split(seq_along(orig_pair_list),
#                   ceiling(seq_along(orig_pair_list) / chunk_size))
#
#   message("\nStep 2: Computing tile patterns...")
#   if (is_windows) {
#     future::plan(future::multisession, workers = n_cores)
#     processed_chunks <- future_lapply(chunks, function(chunk_indices) {
#       chunk_data <- orig_pair_list[chunk_indices]
#       compute_tile_patterns_for_pairs(
#         chunk_data,
#         class_col = class_col,
#         p_threshold = p_threshold,
#         abs_diff_threshold = abs_diff_threshold,
#         split_method = split_method
#       )
#     })
#   } else {
#     processed_chunks <- pbmclapply(chunks, function(chunk_indices) {
#       chunk_data <- orig_pair_list[chunk_indices]
#       compute_tile_patterns_for_pairs(
#         chunk_data,
#         class_col = class_col,
#         p_threshold = p_threshold,
#         abs_diff_threshold = abs_diff_threshold,
#         split_method = split_method
#       )
#     }, mc.cores = n_cores)
#   }
#
#   orig_pair_list <- unlist(processed_chunks, recursive = FALSE)
#   names(orig_pair_list) <- names(create_pairwise_datasets(data, class_col = class_col))
#
#   if (is.null(orig_pair_list[[1]]$results$tile_pattern)) {
#     stop("Tile pattern computation failed. Check compute_tile_patterns_for_pairs()")
#   }
#
#   if (check_tau) {
#     xor_pairs <- names(orig_pair_list)[vapply(orig_pair_list, function(x) {
#       isTRUE(x$results$tile_pattern$xor_shape_detected)
#     }, logical(1))]
#
#     if (length(xor_pairs) > 0) {
#       message("\nStep 3: Computing classwise tau for ", length(xor_pairs), " XOR pairs...")
#       if (is_windows) {
#         chunks <- split(xor_pairs,
#                         ceiling(seq_along(xor_pairs) / chunk_size))
#
#         processed_tau_chunks <- future_lapply(chunks, function(chunk_pairs) {
#           chunk_data <- orig_pair_list[chunk_pairs]
#           compute_classwise_tau_for_pairs(
#             chunk_data,
#             class_col = class_col,
#             tau_threshold = tau_threshold
#           )
#         })
#         tau_results <- unlist(processed_tau_chunks, recursive = FALSE)
#
#       } else {
#         chunks <- split(xor_pairs,
#                         ceiling(seq_along(xor_pairs) / chunk_size))
#
#         processed_tau_chunks <- pbmclapply(chunks, function(chunk_pairs) {
#           chunk_data <- orig_pair_list[chunk_pairs]
#           compute_classwise_tau_for_pairs(
#             chunk_data,
#             class_col = class_col,
#             tau_threshold = tau_threshold
#           )
#         }, mc.cores = n_cores)
#         tau_results <- unlist(processed_tau_chunks, recursive = FALSE)
#       }
#       orig_pair_list[xor_pairs] <- tau_results
#     }
#   }
#
#   if (compute_axes_parallel_significance) {
#     message("\nStep 4: Computing Wilcoxon significance tests...")
#     chunks <- split(seq_along(orig_pair_list),
#                     ceiling(seq_along(orig_pair_list) / chunk_size))
#
#     if (is_windows) {
#       processed_wilcox_chunks <- future_lapply(chunks, function(chunk_indices) {
#         chunk_data <- orig_pair_list[chunk_indices]
#         compute_tile_wilcox_significance_for_pairs(
#           chunk_data,
#           class_col = class_col,
#           p_threshold = p_threshold,
#           split_method = split_method
#         )
#       })
#     } else {
#       processed_wilcox_chunks <- pbmclapply(chunks, function(chunk_indices) {
#         chunk_data <- orig_pair_list[chunk_indices]
#         compute_tile_wilcox_significance_for_pairs(
#           chunk_data,
#           class_col = class_col,
#           p_threshold = p_threshold,
#           split_method = split_method
#         )
#       }, mc.cores = n_cores)
#     }
#
#     orig_pair_list <- unlist(processed_wilcox_chunks, recursive = FALSE)
#   }
#
#   if (is_windows) {
#     future::plan(future::sequential)
#   }
#
#   message("\nFinalizing results...")
#   if (!exists("create_results_df_from_orig_pair_list")) {
#     stop("Missing results formatting function")
#   }
#   create_results_df_from_orig_pair_list(orig_pair_list)
# }


# run_xor_detection <- function(
#   data,
#   class_col = "class",
#   check_tau = TRUE,
#   compute_axes_parallel_significance = TRUE,
#   p_threshold = 0.05,
#   tau_threshold = 0.3,
#   abs_diff_threshold = 20,
#   split_method = "quantile",
#   max_cores = NULL
# ) {
#   # Determine number of cores to use
#   n_cores <- determine_n_cores(max_cores)
#   use_parallel <- n_cores > 1

#   # Required packages for parallel processing
#   if (use_parallel) {
#     require(future)
#     require(future.apply)
#     require(pbmcapply)
#   }

#   message("\nStep 1: Creating pairwise datasets...")
#   orig_pair_list <- create_pairwise_datasets(data, class_col = class_col)

#   message("\nStep 2: Computing tile patterns...")
#   if (use_parallel) {
#     is_windows <- .Platform$OS.type == "windows"
#     chunk_size <- ceiling(length(orig_pair_list) / n_cores)
#     chunks <- split(seq_along(orig_pair_list),
#                     ceiling(seq_along(orig_pair_list) / chunk_size))

#     if (is_windows) {
#       future::plan(future::multisession, workers = n_cores)
#       processed_chunks <- future_lapply(chunks, function(chunk_indices) {
#         chunk_data <- orig_pair_list[chunk_indices]
#         compute_tile_patterns_for_pairs(
#           chunk_data,
#           class_col = class_col,
#           p_threshold = p_threshold,
#           abs_diff_threshold = abs_diff_threshold,
#           split_method = split_method
#         )
#       })
#     } else {
#       processed_chunks <- pbmclapply(chunks, function(chunk_indices) {
#         chunk_data <- orig_pair_list[chunk_indices]
#         compute_tile_patterns_for_pairs(
#           chunk_data,
#           class_col = class_col,
#           p_threshold = p_threshold,
#           abs_diff_threshold = abs_diff_threshold,
#           split_method = split_method
#         )
#       }, mc.cores = n_cores)
#     }
#     orig_pair_list <- unlist(processed_chunks, recursive = FALSE)
#     names(orig_pair_list) <- names(create_pairwise_datasets(data, class_col = class_col))
#   } else {
#     orig_pair_list <- compute_tile_patterns_for_pairs(
#       orig_pair_list,
#       class_col = class_col,
#       p_threshold = p_threshold,
#       abs_diff_threshold = abs_diff_threshold,
#       split_method = split_method
#     )
#   }

#   # Validate Step 2 output
#   if (is.null(orig_pair_list[[1]]$results$tile_pattern)) {
#     stop("Tile pattern computation failed. Check compute_tile_patterns_for_pairs()")
#   }

#   # Step 3: Classwise tau with safe indexing
#   if (check_tau) {
#     xor_pairs <- names(orig_pair_list)[vapply(orig_pair_list, function(x) {
#       isTRUE(x$results$tile_pattern$xor_shape_detected)
#     }, logical(1))]

#     if (length(xor_pairs) > 0) {
#       message("\nStep 3: Computing classwise tau for ", length(xor_pairs), " XOR pairs...")
#       if (use_parallel) {
#         chunks <- split(xor_pairs,
#                         ceiling(seq_along(xor_pairs) / chunk_size))

#         if (is_windows) {
#           processed_tau_chunks <- future_lapply(chunks, function(chunk_pairs) {
#             chunk_data <- orig_pair_list[chunk_pairs]
#             compute_classwise_tau_for_pairs(
#               chunk_data,
#               class_col = class_col,
#               tau_threshold = tau_threshold
#             )
#           })
#         } else {
#           processed_tau_chunks <- pbmclapply(chunks, function(chunk_pairs) {
#             chunk_data <- orig_pair_list[chunk_pairs]
#             compute_classwise_tau_for_pairs(
#               chunk_data,
#               class_col = class_col,
#               tau_threshold = tau_threshold
#             )
#           }, mc.cores = n_cores)
#         }
#         tau_results <- unlist(processed_tau_chunks, recursive = FALSE)
#         orig_pair_list[xor_pairs] <- tau_results
#       } else {
#         orig_pair_list[xor_pairs] <- compute_classwise_tau_for_pairs(
#           orig_pair_list[xor_pairs],
#           class_col = class_col,
#           tau_threshold = tau_threshold
#         )
#       }
#     }
#   }

#   # Step 4: Wilcoxon tests
#   if (compute_axes_parallel_significance) {
#     message("\nStep 4: Computing Wilcoxon significance tests...")
#     if (use_parallel) {
#       chunks <- split(seq_along(orig_pair_list),
#                       ceiling(seq_along(orig_pair_list) / chunk_size))

#       if (is_windows) {
#         processed_wilcox_chunks <- future_lapply(chunks, function(chunk_indices) {
#           chunk_data <- orig_pair_list[chunk_indices]
#           compute_tile_wilcox_significance_for_pairs(
#             chunk_data,
#             class_col = class_col,
#             p_threshold = p_threshold,
#             split_method = split_method
#           )
#         })
#       } else {
#         processed_wilcox_chunks <- pbmclapply(chunks, function(chunk_indices) {
#           chunk_data <- orig_pair_list[chunk_indices]
#           compute_tile_wilcox_significance_for_pairs(
#             chunk_data,
#             class_col = class_col,
#             p_threshold = p_threshold,
#             split_method = split_method
#           )
#         }, mc.cores = n_cores)
#       }
#       orig_pair_list <- unlist(processed_wilcox_chunks, recursive = FALSE)
#     } else {
#       orig_pair_list <- compute_tile_wilcox_significance_for_pairs(
#         orig_pair_list,
#         class_col = class_col,
#         p_threshold = p_threshold,
#         split_method = split_method
#       )
#     }
#   }

#   if (use_parallel && is_windows) {
#     future::plan(future::sequential)
#   }

#   message("\nFinalizing results...")
#   if (!exists("create_results_df_from_orig_pair_list")) {
#     stop("Missing results formatting function")
#   }

#   # Return both the detailed pair list and the results dataframe
#   list(
#     results_df = create_results_df_from_orig_pair_list(orig_pair_list),
#     pair_list = orig_pair_list
#   )
# }


run_xor_detection <- function(
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

  # Required packages for parallel processing
  if (use_parallel) {
    requireNamespace("future")
    requireNamespace("future.apply")
    requireNamespace("pbmcapply")
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
    is_windows <- .Platform$OS.type == "windows"
    chunk_size <- ceiling(length(orig_pair_list) / n_cores)
    chunks <- split(seq_along(orig_pair_list),
                    ceiling(seq_along(orig_pair_list) / chunk_size))

    if (is_windows) {
      future::plan(future::multisession, workers = n_cores)
      processed_chunks <- future.apply::future_lapply(chunks, function(chunk_indices) {
        chunk_data <- orig_pair_list[chunk_indices]
        compute_tile_patterns_for_pairs(
          chunk_data,
          class_col = class_col,
          p_threshold = p_threshold,
          abs_diff_threshold = abs_diff_threshold,
          split_method = split_method
        )
      })
    } else {
      processed_chunks <- pbmcapply::pbmclapply(chunks, function(chunk_indices) {
        chunk_data <- orig_pair_list[chunk_indices]
        compute_tile_patterns_for_pairs(
          chunk_data,
          class_col = class_col,
          p_threshold = p_threshold,
          abs_diff_threshold = abs_diff_threshold,
          split_method = split_method
        )
      }, mc.cores = n_cores)
    }
    orig_pair_list <- unlist(processed_chunks, recursive = FALSE)
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
        chunks <- split(xor_pairs,
                        ceiling(seq_along(xor_pairs) / chunk_size))

        if (is_windows) {
          processed_tau_chunks <- future.apply::future_lapply(chunks, function(chunk_pairs) {
            chunk_data <- orig_pair_list[chunk_pairs]
            compute_classwise_tau_for_pairs(
              chunk_data,
              class_col = class_col,
              tau_threshold = tau_threshold
            )
          })
        } else {
          processed_tau_chunks <- pbmcapply::pbmclapply(chunks, function(chunk_pairs) {
            chunk_data <- orig_pair_list[chunk_pairs]
            compute_classwise_tau_for_pairs(
              chunk_data,
              class_col = class_col,
              tau_threshold = tau_threshold
            )
          }, mc.cores = n_cores)
        }
        tau_results <- unlist(processed_tau_chunks, recursive = FALSE)
        orig_pair_list[xor_pairs] <- tau_results
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
      chunks <- split(seq_along(orig_pair_list),
                      ceiling(seq_along(orig_pair_list) / chunk_size))

      if (is_windows) {
        processed_wilcox_chunks <- future.apply::future_lapply(chunks, function(chunk_indices) {
          chunk_data <- orig_pair_list[chunk_indices]
          compute_tile_wilcox_significance_for_pairs(
            chunk_data,
            class_col = class_col,
            p_threshold = p_threshold,
            split_method = split_method
          )
        })
      } else {
        processed_wilcox_chunks <- pbmcapply::pbmclapply(chunks, function(chunk_indices) {
          chunk_data <- orig_pair_list[chunk_indices]
          compute_tile_wilcox_significance_for_pairs(
            chunk_data,
            class_col = class_col,
            p_threshold = p_threshold,
            split_method = split_method
          )
        }, mc.cores = n_cores)
      }
      orig_pair_list <- unlist(processed_wilcox_chunks, recursive = FALSE)
    } else {
      orig_pair_list <- compute_tile_wilcox_significance_for_pairs(
        orig_pair_list,
        class_col = class_col,
        p_threshold = p_threshold,
        split_method = split_method
      )
    }
  }

  if (use_parallel && is_windows) {
    future::plan(future::sequential)
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


# Helper function to determine number of cores
determine_n_cores <- function(max_cores = NULL) {
  if (!is.null(max_cores)) {
    return(max(1, min(parallel::detectCores() - 1, max_cores)))
  }
  max(1, parallel::detectCores() - 1)
}

# Helper function for progress bar in sequential processing
lapply_with_bar <- function(X, FUN, ...) {
  pb <- txtProgressBar(min = 0, max = length(X), style = 3)
  result <- vector("list", length(X))
  for (i in seq_along(X)) {
    result[[i]] <- FUN(X[[i]], ...)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  result
}

# results <- run_xor_detection(data = XOR_data, class_col = "class")
# results <- run_xor_detection_parallel(data = XOR_data, class_col = "class")

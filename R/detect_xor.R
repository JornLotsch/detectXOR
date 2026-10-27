#' @title Detect XOR-like patterns in pairwise variable relationships
#' @description Screens all pairs of numeric variables in a two-class data set for
#'   XOR-like patterns. Each pair is partitioned into a 3 x 3 tertile grid, and
#'   the class composition of the four critical tiles of the default (corner
#'   tiles) and rotated (edge-center tiles) XOR manifestation is evaluated.
#' @param data Data frame containing variables and class column.
#' @param class_col Character. Name of the class column (default: "class").
#' @param check_tau Logical. Whether to require opposite class-wise Kendall tau
#'   correlations for a positive XOR call (default: TRUE).
#' @param compute_axes_parallel_significance Logical. Whether to compute the four
#'   group-wise Wilcoxon tests along the axes of the critical tiles and require
#'   all of them to be significant for a positive XOR call (default: FALSE). The
#'   tests are not applied to pairs of binary variables.
#' @param p_threshold Numeric. Significance level for the adjusted permutation
#'   p-value, the one-sided lower confidence bound of the XOR purity, and the
#'   Wilcoxon tests (default: 0.05).
#' @param tau_threshold Numeric. Minimum absolute class-wise Kendall tau (default: 0.3).
#' @param purity_threshold Numeric in (0.5, 1). Minimum XOR purity that the lower
#'   confidence bound must reach (default: 0.8).
#' @param min_tile_share Numeric in [0, 0.5). Minimum share of a class' points in
#'   the critical tiles that each of its two expected tiles must contain
#'   (default: 0.05). Rejects incomplete XOR patterns with an empty or nearly
#'   empty quadrant, while tolerating unequal tile populations, e.g., for binary
#'   variables with unbalanced value frequencies.
#' @param n_perm Integer. Number of class-label permutations for the null
#'   distribution (default: 999).
#' @param p_adjust_method Character. Multiple-testing adjustment across all
#'   variable pairs: "maxT" (Westfall-Young single-step max-T, controls the
#'   family-wise error rate; default) or any method of \code{stats::p.adjust}
#'   (e.g., "BH", "holm", "none").
#' @param split_method Character. Method for splitting tiles ("quantile" or "range") (default: "quantile").
#' @param split_quantiles Numeric vector of length 2. Lower and upper quantile of the
#'   tile boundaries for \code{split_method = "quantile"} (default: tertiles, c(1/3, 2/3)).
#' @param max_cores Integer. Maximum number of cores to use for parallel processing (default: 1).
#' @param extreme_handling Character. Method for handling extreme values: "winsorize",
#'   "trim" (remove rows with extreme values; "remove" is accepted as a synonym),
#'   or "none" (default: "winsorize").
#' @param winsor_limits Numeric vector of length 2. Lower and upper quantiles for
#'   winsorization or trimming (default: c(0.05, 0.95)).
#' @param scale_data Logical. Whether to scale variables to unit variance (default: TRUE).
#' @param use_complete Logical. Whether to use only complete cases (default: TRUE).
#' @param balance_classes Logical. Whether to randomly down-sample the larger class
#'   to the size of the smaller class before screening (default: FALSE). Recommended
#'   for marked class imbalance, under which the quantile-based tile boundaries can
#'   misrepresent rotated XOR patterns.
#' @param abs_diff_threshold Deprecated and ignored (used by detectXOR < 0.2.0).
#' @return List containing:
#'   \itemize{
#'     \item results_df: Data frame summarizing XOR detection results for all variable pairs
#'     \item pair_list: List of detailed results for each variable pair
#'     \item settings: List of the parameter values used
#'   }
#' @details For a given XOR manifestation, the XOR rule assigns one class to two
#'   diagonally opposite critical tiles and the other class to the remaining two.
#'   The XOR purity is the balanced accuracy of this rule among the observations
#'   in the critical tiles, maximized over both label assignments (so that the
#'   result does not depend on class coding) and over the two manifestations.
#'   It is 0.5 for no XOR structure and 1 for a perfect XOR pattern.
#'
#'   A pair is called XOR-positive if all of the following hold:
#'   \enumerate{
#'     \item the purity exceeds chance: the permutation p-value, adjusted for the
#'       number of variable pairs, is below \code{p_threshold}. All pairs share
#'       the same label permutations, which allows max-T adjustment;
#'     \item the purity is large: its one-sided \code{1 - p_threshold} lower
#'       confidence bound is at least \code{purity_threshold};
#'     \item the pattern is complete: each class populates both of its expected
#'       tiles with at least \code{min_tile_share} of its critical-tile points;
#'     \item if \code{check_tau}, the class-wise Kendall tau correlations have
#'       opposite signs and absolute values of at least \code{tau_threshold};
#'     \item if \code{compute_axes_parallel_significance}, all four group-wise
#'       Wilcoxon tests are significant (continuous variables only).
#'   }
#'   Criteria 4 and 5 are evaluated only for pairs meeting criteria 1 to 3.
#'   Permutations are drawn from the session's random number generator; call
#'   \code{set.seed()} beforehand for reproducible p-values. Parallel processing
#'   is used when \code{max_cores > 1} and the required packages are available.
#'
#' @examples
#' data(XOR_data)
#' set.seed(42)
#' result <- detect_xor(XOR_data, class_col = "class", n_perm = 199)
#' result$results_df[, c("var1", "var2", "xor_shape_detected", "xor_pattern", "xor_purity")]
#' @seealso \code{\link{generate_spaghetti_plot_from_results}}, \code{\link{generate_xy_plot_from_results}}
#' @author Jorn Lotsch, Alfred Ultsch
#' @importFrom stats wilcox.test cor.test p.adjust
#' @export
detect_xor <- function(
  data,
  class_col = "class",
  check_tau = TRUE,
  compute_axes_parallel_significance = FALSE,
  p_threshold = 0.05,
  tau_threshold = 0.3,
  purity_threshold = 0.8,
  min_tile_share = 0.05,
  n_perm = 999,
  p_adjust_method = "maxT",
  split_method = "quantile",
  split_quantiles = c(1/3, 2/3),
  max_cores = 1,
  extreme_handling = "winsorize",
  winsor_limits = c(0.05, 0.95),
  scale_data = TRUE,
  use_complete = TRUE,
  balance_classes = FALSE,
  abs_diff_threshold = NULL
) {
  # Input validation
  if (!is.null(abs_diff_threshold)) {
    warning("'abs_diff_threshold' is deprecated and ignored since detectXOR 0.2.0.", call. = FALSE)
  }
  if (!is.numeric(p_threshold) || p_threshold <= 0 || p_threshold >= 1) {
    stop("'p_threshold' must be between 0 and 1", call. = FALSE)
  }
  if (!is.numeric(purity_threshold) || purity_threshold <= 0.5 || purity_threshold >= 1) {
    stop("'purity_threshold' must be between 0.5 and 1", call. = FALSE)
  }
  if (!is.numeric(min_tile_share) || min_tile_share < 0 || min_tile_share >= 0.5) {
    stop("'min_tile_share' must be in [0, 0.5)", call. = FALSE)
  }
  if (!is.numeric(n_perm) || length(n_perm) != 1 || n_perm < 1) {
    stop("'n_perm' must be a positive integer", call. = FALSE)
  }
  n_perm <- as.integer(n_perm)
  if (1 / (n_perm + 1) >= p_threshold) {
    warning("'n_perm' is too small to reach p-values below 'p_threshold'.", call. = FALSE)
  }
  p_adjust_method <- match.arg(p_adjust_method, c("maxT", stats::p.adjust.methods))
  split_method <- match.arg(split_method, c("quantile", "range"))
  if (!is.numeric(split_quantiles) || length(split_quantiles) != 2 ||
      split_quantiles[1] <= 0 || split_quantiles[2] >= 1 || split_quantiles[1] >= split_quantiles[2]) {
    stop("'split_quantiles' must be two increasing values in (0, 1)", call. = FALSE)
  }
  extreme_handling <- match.arg(extreme_handling, c("winsorize", "trim", "remove", "none"))
  if (extreme_handling == "remove") extreme_handling <- "trim"

  # Determine number of cores to use
  n_cores <- determine_n_cores(max_cores)
  if (n_cores > 1 && !check_parallel_packages()) {
    warning("Parallel processing packages not available. Using sequential processing.")
    n_cores <- 1
  }

  settings <- list(
    class_col = class_col, check_tau = check_tau,
    compute_axes_parallel_significance = compute_axes_parallel_significance,
    p_threshold = p_threshold, tau_threshold = tau_threshold,
    purity_threshold = purity_threshold, min_tile_share = min_tile_share,
    n_perm = n_perm, p_adjust_method = p_adjust_method, split_method = split_method,
    split_quantiles = split_quantiles,
    extreme_handling = extreme_handling, winsor_limits = winsor_limits,
    scale_data = scale_data, use_complete = use_complete, balance_classes = balance_classes
  )

  # Optional class balancing by random down-sampling of the larger class
  if (balance_classes) {
    if (!class_col %in% names(data)) stop("Class column '", class_col, "' not found in data.", call. = FALSE)
    rows_by_class <- split(seq_len(nrow(data)), data[[class_col]], drop = TRUE)
    n_min <- min(lengths(rows_by_class))
    keep <- unlist(lapply(rows_by_class, function(i) i[sample.int(length(i), n_min)]), use.names = FALSE)
    data <- data[sort(keep), , drop = FALSE]
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
  if (length(orig_pair_list) == 0) {
    return(list(results_df = NULL, pair_list = orig_pair_list, settings = settings))
  }

  # Class labels and label permutations, shared by all pairs
  class_vec <- as.character(orig_pair_list[[1]]$data[[class_col]])
  class_levels <- sort(unique(class_vec))
  perm_class2 <- vapply(
    seq_len(n_perm),
    function(b) as.numeric(sample(class_vec) == class_levels[2]),
    numeric(length(class_vec))
  )
  perm_class2 <- matrix(perm_class2, nrow = length(class_vec))

  message("\nStep 2: Computing tile patterns and permutation tests...")
  chunks <- create_balanced_chunks(length(orig_pair_list), n_cores)
  chunk_results <- safe_parallel_execution(
    chunks,
    function(chunk_indices) {
      compute_tile_patterns_for_pairs(
        orig_pair_list[chunk_indices],
        class_col = class_col,
        class_levels = class_levels,
        perm_class2 = perm_class2,
        p_threshold = p_threshold,
        split_method = split_method,
        split_quantiles = split_quantiles,
        min_tile_share = min_tile_share
      )
    },
    n_cores = n_cores,
    flatten = FALSE
  )
  orig_pair_list <- do.call(c, unname(lapply(chunk_results, `[[`, "pair_list")))
  perm_max <- do.call(pmax, lapply(chunk_results, `[[`, "perm_max"))

  # Multiple-testing adjustment of the permutation p-values
  has_tiles <- vapply(orig_pair_list, function(x) !is.null(x$results$tile_pattern), logical(1))
  if (!any(has_tiles)) {
    stop("Tile pattern computation failed for all variable pairs.", call. = FALSE)
  }
  purity <- vapply(orig_pair_list[has_tiles], function(x) x$results$tile_pattern$purity, numeric(1))
  p_raw <- vapply(orig_pair_list[has_tiles], function(x) x$results$tile_pattern$p_value, numeric(1))
  p_adj <- if (p_adjust_method == "maxT") {
    vapply(purity, function(s) (1 + sum(perm_max >= s - 1e-12)) / (n_perm + 1), numeric(1))
  } else {
    stats::p.adjust(p_raw, method = p_adjust_method)
  }
  for (i in seq_along(p_adj)) {
    nm <- names(p_adj)[i]
    tp <- orig_pair_list[[nm]]$results$tile_pattern
    tp$p_adjusted <- unname(p_adj[i])
    tp$candidate <- tp$p_adjusted < p_threshold && tp$purity_lower >= purity_threshold && tp$coverage_ok
    orig_pair_list[[nm]]$results$tile_pattern <- tp
  }

  candidates <- names(orig_pair_list)[vapply(orig_pair_list, function(x) {
    isTRUE(x$results$tile_pattern$candidate)
  }, logical(1))]

  # Step 3: Classwise tau for candidate pairs
  if (check_tau && length(candidates) > 0) {
    message("\nStep 3: Computing classwise tau for ", length(candidates), " candidate pairs...")
    orig_pair_list[candidates] <- compute_classwise_tau_for_pairs(
      orig_pair_list[candidates],
      class_col = class_col,
      tau_threshold = tau_threshold
    )
  }

  # Step 4: Wilcoxon tests for candidate pairs
  if (compute_axes_parallel_significance && length(candidates) > 0) {
    message("\nStep 4: Computing Wilcoxon significance tests for ", length(candidates), " candidate pairs...")
    orig_pair_list[candidates] <- compute_tile_wilcox_significance_for_pairs(
      orig_pair_list[candidates],
      class_col = class_col,
      p_threshold = p_threshold,
      split_method = split_method,
      split_quantiles = split_quantiles
    )
  }

  message("\nFinalizing results...")
  list(
    results_df = create_results_df(orig_pair_list, settings),
    pair_list = orig_pair_list,
    settings = settings
  )
}

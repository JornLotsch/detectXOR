#' Tile levels of the 3x3 grid (row-wise from upper left)
#' @noRd
TILE_LEVELS <- c(
  "Upper Left", "Upper Mid", "Upper Right",
  "Mid Left", "Center", "Mid Right",
  "Lower Left", "Lower Mid", "Lower Right"
)

#' Critical tiles of the two XOR manifestations
#'
#' Within each pattern, the first two tiles share one XOR parity (expected to be
#' populated by one class) and the last two tiles share the opposite parity
#' (expected to be populated by the other class).
#' @noRd
XOR_PATTERNS <- list(
  default = c("Lower Left", "Upper Right", "Upper Left", "Lower Right"),
  rotated = c("Mid Left", "Mid Right", "Upper Mid", "Lower Mid")
)

#' Assign 3x3 tiles to data points based on two variables
#' @noRd
assign_3x3_tiles <- function(data, var1_name, var2_name, split_method = "quantile",
                             split_quantiles = c(1/3, 2/3)) {
  # Input validation
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("'data' must be a non-empty data frame", call. = FALSE)
  }
  if (!all(c(var1_name, var2_name) %in% names(data))) {
    stop("Variables '", var1_name, "' and '", var2_name, "' must exist in data", call. = FALSE)
  }
  if (!split_method %in% c("quantile", "range")) {
    stop("'split_method' must be either 'quantile' or 'range'", call. = FALSE)
  }

  # Get variable vectors
  var1_values <- data[[var1_name]]
  var2_values <- data[[var2_name]]

  # Check for sufficient variation
  if (length(unique(var1_values[!is.na(var1_values)])) < 2) {
    warning("Variable '", var1_name, "' has insufficient variation for tile assignment", call. = FALSE)
  }
  if (length(unique(var2_values[!is.na(var2_values)])) < 2) {
    warning("Variable '", var2_name, "' has insufficient variation for tile assignment", call. = FALSE)
  }

  # Determine split points
  if (split_method == "quantile") {
    x_limits <- compute_quantile_splits(var1_values, split_quantiles)
    y_limits <- compute_quantile_splits(var2_values, split_quantiles)
  } else {
    x_limits <- split_range_into_thirds(var1_values)
    y_limits <- split_range_into_thirds(var2_values)
  }

  # Internal function for flexible group assignment: if no values fall strictly
  # between the limits (e.g., binary data), values at the upper limit are "high"
  assign_flex_groups <- function(values, limits, labels) {
    breaks <- c(-Inf, limits, Inf)
    groups <- cut(values, breaks = breaks, labels = labels, include.lowest = TRUE)
    if (!any(values > limits[1] & values < limits[2], na.rm = TRUE) &&
        any(values == limits[2], na.rm = TRUE) &&
        !any(values > limits[2], na.rm = TRUE)) {
      groups[which(values == limits[2])] <- labels[length(labels)]
      groups <- factor(groups, levels = labels)
    }
    return(groups)
  }

  # Create group assignments
  x_group <- assign_flex_groups(var1_values, x_limits, c("Left", "Mid", "Right"))
  y_group <- assign_flex_groups(var2_values, y_limits, c("Lower", "Mid", "Upper"))

  factor(map_coordinates_to_tile(x_group, y_group), levels = TILE_LEVELS)
}


#' Compute robust quantile splits with fallback to range splits
#' @noRd
compute_quantile_splits <- function(x, probs = c(1/3, 2/3)) {
  # Try quantile method first
  quantile_splits <- quantile(x, probs = probs, na.rm = TRUE)

  # Check if quantiles provide sufficient separation
  if (length(unique(quantile_splits)) < 2) {
    # Fall back to range-based splitting
    return(split_range_into_thirds(x))
  }

  return(quantile_splits)
}

#' Split a numeric vector into three equal-width bins
#' @noRd
split_range_into_thirds <- function(x) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)

  if (is.infinite(min_val) || is.infinite(max_val) || min_val == max_val) {
    # Handle edge cases
    return(c(min_val, min_val))
  }

  range_size <- max_val - min_val
  break1 <- min_val + range_size / 3
  break2 <- min_val + 2 * range_size / 3

  c(break1, break2)
}

#' Map x,y group labels to tile names (vectorized)
#' @noRd
map_coordinates_to_tile <- function(x, y) {
  tile_mapping <- c(
    "Left Upper" = "Upper Left",
    "Right Upper" = "Upper Right",
    "Left Lower" = "Lower Left",
    "Right Lower" = "Lower Right",
    "Mid Upper" = "Upper Mid",
    "Left Mid" = "Mid Left",
    "Right Mid" = "Mid Right",
    "Mid Lower" = "Lower Mid",
    "Mid Mid" = "Center"
  )

  result <- unname(tile_mapping[paste(x, y)])
  result[is.na(x) | is.na(y)] <- NA
  result
}

#' Count class members per tile for one or many class labelings
#'
#' @param tiles Factor of tile assignments with levels \code{TILE_LEVELS}.
#' @param class2 Numeric 0/1 matrix (n x B); 1 marks the second class. Each
#'   column is one labeling (observed or permuted).
#' @return List with 9 x B matrices \code{n1} and \code{n2} of class counts per tile.
#' @noRd
tile_class_counts <- function(tiles, class2) {
  ok <- !is.na(tiles)
  onehot <- outer(as.integer(tiles[ok]), seq_along(TILE_LEVELS), "==") * 1
  n2 <- crossprod(onehot, class2[ok, , drop = FALSE])
  n1 <- colSums(onehot) - n2
  rownames(n1) <- rownames(n2) <- TILE_LEVELS
  list(n1 = n1, n2 = n2)
}

#' Balanced accuracy of the XOR rule within the critical tiles of a pattern
#'
#' The rule assigns class 1 to the first parity (first two critical tiles) and
#' class 2 to the second parity. Values below 0.5 indicate the reversed label
#' assignment; the label-invariant XOR purity is \code{0.5 + abs(ba - 0.5)}.
#' @param counts Output of \code{tile_class_counts}.
#' @param pattern Character, "default" or "rotated".
#' @return Numeric vector of length B.
#' @noRd
compute_xor_balanced_accuracy <- function(counts, pattern) {
  tiles <- XOR_PATTERNS[[pattern]]
  n1 <- counts$n1[tiles, , drop = FALSE]
  n2 <- counts$n2[tiles, , drop = FALSE]
  m1 <- colSums(n1)
  m2 <- colSums(n2)
  ba <- 0.5 * (colSums(n1[1:2, , drop = FALSE]) / m1 + colSums(n2[3:4, , drop = FALSE]) / m2)
  ba[m1 == 0 | m2 == 0] <- 0.5
  ba
}

#' Label-invariant XOR purity, maximized over both XOR manifestations
#' @noRd
compute_max_xor_purity <- function(counts) {
  purities <- lapply(names(XOR_PATTERNS), function(p) {
    0.5 + abs(compute_xor_balanced_accuracy(counts, p) - 0.5)
  })
  do.call(pmax, purities)
}

#' One-sided lower confidence bound of the XOR purity
#'
#' Purity is the mean of two binomial proportions (class 1 and class 2 points in
#' their expected tiles). The bound uses the Agresti-Caffo adjustment (one
#' pseudo-success and one pseudo-failure per class) and a normal approximation.
#' @noRd
purity_lower_bound <- function(x1, m1, x2, m2, alpha = 0.05) {
  p1 <- (x1 + 1) / (m1 + 2)
  p2 <- (x2 + 1) / (m2 + 2)
  se <- 0.5 * sqrt(p1 * (1 - p1) / (m1 + 2) + p2 * (1 - p2) / (m2 + 2))
  0.5 * (p1 + p2) - stats::qnorm(1 - alpha) * se
}

#' Main tile pattern analysis for a single pair
#'
#' @param data Data frame with the two variables and the class column.
#' @param var1_name,var2_name Variable names.
#' @param class_name Name of the class column.
#' @param class_levels Character vector of the two class labels.
#' @param perm_class2 Numeric 0/1 matrix (n x B) of permuted class-2 indicators,
#'   shared by all pairs.
#' @param p_threshold Alpha for the lower confidence bound of the purity.
#' @param split_method "quantile" or "range".
#' @param split_quantiles Lower and upper quantile of the tile boundaries.
#' @param min_tile_share Minimum share of a class' critical-tile points that each
#'   of its two expected tiles must hold.
#' @return List with tile pattern results, including the permutation null
#'   distribution \code{perm_stat} of the maximized purity.
#' @noRd
analyze_tile_dominance <- function(
  data, var1_name, var2_name, class_name, class_levels, perm_class2,
  p_threshold = 0.05,
  split_method = "quantile",
  split_quantiles = c(1/3, 2/3),
  min_tile_share = 0.05
) {
  # Input validation
  validate_tile_analysis_inputs(data, var1_name, var2_name, class_name, p_threshold)

  # Assign tiles and count classes per tile
  tiles <- assign_3x3_tiles(data, var1_name, var2_name, split_method, split_quantiles)
  is_class2 <- as.character(data[[class_name]]) == class_levels[2]
  observed_counts <- tile_class_counts(tiles, matrix(as.numeric(is_class2)))

  # Observed balanced accuracy and purity for both manifestations
  ba <- vapply(names(XOR_PATTERNS), function(p) compute_xor_balanced_accuracy(observed_counts, p), numeric(1))
  purity <- 0.5 + abs(ba - 0.5)
  scenario <- names(XOR_PATTERNS)[which.max(purity)]

  # Permutation null distribution of the maximized purity
  perm_stat <- compute_max_xor_purity(tile_class_counts(tiles, perm_class2))
  p_value <- (1 + sum(perm_stat >= max(purity) - 1e-12)) / (length(perm_stat) + 1)

  # Observed table of the chosen pattern, oriented so that class 1 is expected
  # in the first two columns and class 2 in the last two columns
  chosen_tiles <- XOR_PATTERNS[[scenario]]
  if (ba[scenario] < 0.5) chosen_tiles <- chosen_tiles[c(3, 4, 1, 2)]
  observed_table <- rbind(observed_counts$n1[chosen_tiles, 1], observed_counts$n2[chosen_tiles, 1])
  dimnames(observed_table) <- list(class_levels, chosen_tiles)

  x1 <- sum(observed_table[1, 1:2])
  m1 <- sum(observed_table[1, ])
  x2 <- sum(observed_table[2, 3:4])
  m2 <- sum(observed_table[2, ])
  purity_lower <- purity_lower_bound(x1, m1, x2, m2, alpha = p_threshold)

  # Coverage: each class must populate both of its expected tiles
  tile_shares <- c(
    observed_table[1, 1] / x1, observed_table[1, 2] / x1,
    observed_table[2, 3] / x2, observed_table[2, 4] / x2
  )
  names(tile_shares) <- chosen_tiles
  coverage_ok <- x1 > 0 && x2 > 0 && all(tile_shares >= min_tile_share)

  list(
    observed_table = observed_table,
    xor_pattern = scenario,
    purity = unname(purity[scenario]),
    purity_default = unname(purity["default"]),
    purity_rotated = unname(purity["rotated"]),
    purity_lower = purity_lower,
    tile_shares = tile_shares,
    coverage_ok = coverage_ok,
    p_value = p_value,
    perm_stat = perm_stat
  )
}

#' Validate inputs for tile analysis
#' @noRd
validate_tile_analysis_inputs <- function(data, var1_name, var2_name, class_name, p_threshold) {
  if (!is.data.frame(data) || nrow(data) == 0) {
    stop("'data' must be a non-empty data frame", call. = FALSE)
  }

  required_vars <- c(var1_name, var2_name, class_name)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Missing variables in data: ", paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  if (!is.numeric(p_threshold) || p_threshold <= 0 || p_threshold >= 1) {
    stop("'p_threshold' must be between 0 and 1", call. = FALSE)
  }
}

#' Compute tile pattern analysis for all pairwise datasets
#'
#' @param orig_pair_list List. Output from create_pairwise_datasets.
#' @param class_col Character. Name of the class column.
#' @param class_levels Character vector of the two class labels.
#' @param perm_class2 Numeric 0/1 matrix (n x B) of permuted class-2 indicators.
#' @param p_threshold Numeric. Alpha for the lower confidence bound of the purity.
#' @param split_method Character. Method for splitting tiles ("quantile" or "range").
#' @param split_quantiles Numeric. Lower and upper quantile of the tile boundaries.
#' @param min_tile_share Numeric. Minimum share per expected tile.
#' @return List with \code{pair_list} (the input list with results$tile_pattern
#'   filled for each pair) and \code{perm_max}, the maximum of the permuted
#'   purities over all pairs for each permutation (for max-T adjustment).
#' @noRd
compute_tile_patterns_for_pairs <- function(
  orig_pair_list,
  class_col,
  class_levels,
  perm_class2,
  p_threshold = 0.05,
  split_method = "quantile",
  split_quantiles = c(1/3, 2/3),
  min_tile_share = 0.05
) {
  perm_max <- rep(0.5, ncol(perm_class2))

  for (nm in names(orig_pair_list)) {
    df <- orig_pair_list[[nm]]$data
    var_names <- setdiff(names(df), class_col)

    tile_pattern_result <- tryCatch(
      analyze_tile_dominance(
        data = df,
        var1_name = var_names[1],
        var2_name = var_names[2],
        class_name = class_col,
        class_levels = class_levels,
        perm_class2 = perm_class2,
        p_threshold = p_threshold,
        split_method = split_method,
        split_quantiles = split_quantiles,
        min_tile_share = min_tile_share
      ),
      error = function(e) {
        warning("Error processing pair '", nm, "': ", e$message, call. = FALSE)
        NULL
      }
    )
    if (is.null(tile_pattern_result)) next

    # Keep only the running maximum of the permutation statistics
    perm_max <- pmax(perm_max, tile_pattern_result$perm_stat)
    tile_pattern_result$perm_stat <- NULL
    orig_pair_list[[nm]]$results$tile_pattern <- tile_pattern_result
  }

  list(pair_list = orig_pair_list, perm_max = perm_max)
}

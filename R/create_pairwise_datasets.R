#' Create Pairwise Datasets for Two-Class Data (with Structure and Checks)
#'
#' This function preprocesses a data frame (optionally winsorizing, trimming, and/or scaling numeric variables)
#' and returns a named list of lists, each with two elements:
#' \itemize{
#'   \item \code{data}: a data frame containing the class column and a pair of original variables.
#'   \item \code{results}: set to \code{NULL} (for future use).
#' }
#' Each data frame corresponds to a unique pair of numeric variables, and contains three columns:
#' the class column and the two variables.
#'
#' @param data A data frame containing the data.
#' @param class_col Character string. The name of the class column (must have exactly two unique values).
#' @param extreme_handling Character. How to handle extreme values in numeric columns.
#'   One of \code{"none"} (default, do nothing), \code{"winsorize"} (replace extremes with quantile values),
#'   or \code{"trim"} (remove rows with extreme values).
#' @param winsor_limits Numeric vector of length 2. The lower and upper quantile limits for winsorizing or trimming (default: c(0.05, 0.95)).
#' @param scale_data Logical. Whether to scale numeric columns to mean 0, sd 1 (default: TRUE).
#' @param use_complete Logical. Whether to use only rows with complete data (default: TRUE).
#'
#' @return A named list. Each element is a list with two elements:
#'   \itemize{
#'     \item \code{data}: a data frame with the class column and two numeric variables.
#'     \item \code{results}: \code{NULL} (for later use).
#'   }
#'   The names of the list elements are of the form "var1||var2".
#'
#' @details
#' Only the first two unique values in the class column are used. Non-numeric columns (other than the class column)
#' are not allowed for variable pairs. If \code{extreme_handling = "trim"}, rows with any numeric variable
#' outside the specified quantiles are removed before further processing.
#' Variables with no variance (containing only one unique value) are excluded from creating pairs.
#'
#' @importFrom DescTools Winsorize
#'
#' @examples
#' # pairwise_dfs <- create_pairwise_datasets(
#' #   data = XOR_data,
#' #   class_col = "class"
#' # )
#' # str(pairwise_dfs[[1]])
#'
create_pairwise_datasets <- function(
    data,
    class_col,
    extreme_handling = c("none", "winsorize", "trim"),
    winsor_limits = c(0.05, 0.95),
    scale_data = TRUE,
    use_complete = TRUE
) {
  extreme_handling <- match.arg(extreme_handling)

  # --- Checks ---
  if (!is.data.frame(data)) stop("Input 'data' must be a data frame.")
  if (!class_col %in% names(data)) stop(sprintf("Class column '%s' not found in data.", class_col))
  class_vals <- unique(na.omit(data[[class_col]]))
  if (length(class_vals) != 2) stop(sprintf("Class column '%s' must have exactly 2 unique values (after removing NAs).", class_col))

  vars <- setdiff(names(data), class_col)
  if (length(vars) < 2) stop("At least two variable columns (excluding the class column) are required.")
  is_num <- sapply(data[vars], is.numeric)
  if (!all(is_num)) {
    non_num_vars <- vars[!is_num]
    stop(sprintf("All variable columns must be numeric. Non-numeric columns: %s", paste(non_num_vars, collapse = ", ")))
  }

  # --- Check for variables with no variance ---
  var_unique_counts <- sapply(data[vars], function(x) length(unique(na.omit(x))))
  no_variance_vars <- vars[var_unique_counts <= 1]
  if (length(no_variance_vars) > 0) {
    message(sprintf("Removing %d variable(s) with no variance: %s",
                   length(no_variance_vars),
                   paste(no_variance_vars, collapse = ", ")))

    # Remove variables with no variance
    vars <- setdiff(vars, no_variance_vars)

    # Check if we still have enough variables
    if (length(vars) < 2) {
      stop("Not enough variables with variance to create pairs (at least 2 needed).")
    }
  }

  # Helper: Filter complete cases
  filter_complete_cases <- function(df, use_complete = TRUE) {
    if (use_complete) {
      df <- df[complete.cases(df), ]
    }
    return(df)
  }

  # Helper: Winsorize numeric columns
  winsorize_data <- function(df, winsor_limits) {
    as.data.frame(lapply(df, function(x) DescTools::Winsorize(x, val = quantile(x, probs = winsor_limits))))
  }

  # Helper: Trim (cut extremes) from numeric columns
  trim_data <- function(df, winsor_limits) {
    is_num <- sapply(df, is.numeric)
    num_df <- df[, is_num, drop = FALSE]
    bounds <- lapply(num_df, function(x) quantile(x, probs = winsor_limits))
    mask <- rep(TRUE, nrow(df))
    for (v in names(bounds)) {
      mask <- mask & (df[[v]] >= bounds[[v]][1]) & (df[[v]] <= bounds[[v]][2])
    }
    df[mask, , drop = FALSE]
  }

  # Helper: Preprocess data (winsorize/trim/none, scale)
  preprocess_data <- function(df, class_col, extreme_handling, winsor_limits, scale_data) {
    class_vec <- df[[class_col]]
    num_df <- df[, !(names(df) %in% class_col), drop = FALSE]
    is_num <- sapply(num_df, is.numeric)
    num_df_num <- num_df[, is_num, drop = FALSE]
    # Handle extremes
    if (extreme_handling == "winsorize") {
      num_df_num <- winsorize_data(num_df_num, winsor_limits)
    }
    # Scaling
    if (scale_data) {
      num_df_num <- as.data.frame(scale(num_df_num))
    }
    # Recombine with non-numeric columns, if any
    non_num_df <- num_df[, !is_num, drop = FALSE]
    num_df_final <- cbind(num_df_num, non_num_df)
    # Reattach class column as the last column
    num_df_final[[class_col]] <- class_vec
    return(num_df_final)
  }

  # Helper: Create original pairs (now returns a list of lists: list(data=..., results=NULL))
  create_orig_pair_dfs <- function(df, class_col, vars_with_variance) {
    # Use only variables with variance
    vars <- intersect(names(df), vars_with_variance)
    if (length(vars) < 2) return(list())
    pairs <- t(combn(vars, 2))
    pair_list <- list()
    for (i in 1:nrow(pairs)) {
      v1 <- pairs[i, 1]
      v2 <- pairs[i, 2]
      name <- paste(v1, v2, sep = "||")
      pair_data <- data.frame(
        class = df[[class_col]],
        var1 = df[[v1]],
        var2 = df[[v2]]
      )
      names(pair_data) <- c("class", v1, v2)
      pair_list[[name]] <- list(
        data = pair_data,
        results = list(tile_pattern = NULL, tile_wilcox_significance = NULL, classwise_tau = NULL)
      )
    }
    return(pair_list)
  }


  # ---- MAIN FUNCTION BODY ----
  # Step 1: Only keep two classes
  data <- data[data[[class_col]] %in% class_vals, ]

  # Step 2: Handle extremes if trimming (cutting extremes)
  if (extreme_handling == "trim") {
    is_num <- sapply(data, is.numeric)
    num_df <- data[, is_num, drop = FALSE]
    bounds <- lapply(num_df, function(x) quantile(x, probs = winsor_limits))
    mask <- rep(TRUE, nrow(data))
    for (v in names(bounds)) {
      mask <- mask & (data[[v]] >= bounds[[v]][1]) & (data[[v]] <= bounds[[v]][2])
    }
    data <- data[mask, , drop = FALSE]
  }

  # Step 3: Filter for complete cases if requested
  data <- filter_complete_cases(data, use_complete = use_complete)

  # Step 4: Preprocess (winsorize or none, scale)
  data_prep <- preprocess_data(
    data, class_col,
    extreme_handling = extreme_handling,
    winsor_limits = winsor_limits,
    scale_data = scale_data
  )

  # Step 5: Generate all pairwise data frames (as list of lists)
  # Pass the list of variables with variance to the pair creation function
  orig_pair_list <- create_orig_pair_dfs(data_prep, class_col, vars)

  # Return message if no pairs were created
  if (length(orig_pair_list) == 0) {
    message("No variable pairs were created. Check if you have at least two numeric variables with variance.")
  }

  return(orig_pair_list)
}

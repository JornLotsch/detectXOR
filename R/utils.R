#' @title Determine optimal number of cores for parallel processing
#' @description Calculates the number of CPU cores to use for parallel operations,
#' leaving one core free for system operations and respecting user-specified limits.
#' @param max_cores Integer or NULL. Maximum number of cores to use. If NULL,
#' uses all available cores minus one (default: NULL).
#' @return Integer. Number of cores to use for parallel processing (minimum 1).
#' @details This function ensures at least one core remains available for system
#' operations by using \code{parallel::detectCores() - 1}. If \code{max_cores} is
#' specified, it returns the minimum of the detected cores and the user limit.
#' @keywords internal
#' @importFrom parallel detectCores
#' @examples
#' \dontrun{
#' # Use all available cores minus one
#' n_cores <- determine_n_cores()
#'
#' # Limit to maximum 4 cores
#' n_cores <- determine_n_cores(max_cores = 4)
#' }
determine_n_cores <- function(max_cores = NULL) {
  if (!is.null(max_cores)) {
    return(max(1, min(parallel::detectCores() - 1, max_cores)))
  }
  max(1, parallel::detectCores() - 1)
}

#' @title Apply function with text progress bar
#' @description Sequential version of lapply that displays a text progress bar
#' to monitor computation progress, useful for long-running operations.
#' @param X List or vector. Input data to iterate over.
#' @param FUN Function. Function to apply to each element of X.
#' @param ... Additional arguments passed to FUN.
#' @return List. Results of applying FUN to each element of X.
#' @details Creates a text progress bar using \code{txtProgressBar} and updates
#' it after processing each element. The progress bar is automatically closed
#' upon completion. This function is particularly useful when parallel processing
#' is not available or desired.
#' @keywords internal
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @examples
#' \dontrun{
#' # Apply function with progress bar
#' results <- lapply_with_bar(1:100, function(x) {
#'   Sys.sleep(0.1)  # Simulate computation
#'   x^2
#' })
#' }
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

#' @title Check if parallel processing packages are available
#' @description Verifies that required packages for parallel processing are
#' installed and can be loaded, with informative error messages if missing.
#' @param packages Character vector. Names of packages to check (default:
#' c("future", "future.apply", "pbmcapply")).
#' @return Logical. TRUE if all packages are available, FALSE otherwise.
#' @details This function uses \code{requireNamespace} to check package
#' availability without loading them. It's designed to be used before
#' attempting parallel operations to provide clear error messages.
#' @keywords internal
#' @examples
#' \dontrun{
#' # Check if parallel packages are available
#' if (check_parallel_packages()) {
#'   # Proceed with parallel processing
#' } else {
#'   # Fall back to sequential processing
#' }
#' }
check_parallel_packages <- function(packages = c("future", "future.apply", "pbmcapply")) {
  missing_packages <- packages[!vapply(packages, requireNamespace,
                                       logical(1), quietly = TRUE)]

  if (length(missing_packages) > 0) {
    warning("Missing packages for parallel processing: ",
            paste(missing_packages, collapse = ", "),
            ". Install with: install.packages(c(",
            paste0("'", missing_packages, "'", collapse = ", "), "))")
    return(FALSE)
  }
  return(TRUE)
}

#' @title Create balanced chunks for parallel processing
#' @description Splits a vector of indices into approximately equal-sized chunks
#' for efficient parallel processing across multiple cores.
#' @param n_items Integer. Total number of items to split.
#' @param n_cores Integer. Number of chunks/cores to create.
#' @return List. Each element contains indices for one chunk.
#' @details This function ensures balanced workload distribution by creating
#' chunks of similar sizes. The last chunk may contain fewer items if
#' \code{n_items} is not evenly divisible by \code{n_cores}.
#' @keywords internal
#' @examples
#' \dontrun{
#' # Split 100 items across 4 cores
#' chunks <- create_balanced_chunks(100, 4)
#' # Each chunk will have ~25 items
#' }
create_balanced_chunks <- function(n_items, n_cores) {
  chunk_size <- ceiling(n_items / n_cores)
  split(seq_len(n_items), ceiling(seq_len(n_items) / chunk_size))
}

#' @title Safe parallel execution wrapper
#' @description Executes a function in parallel with proper error handling and
#' fallback to sequential processing if parallel execution fails.
#' @param chunks List. Data chunks for parallel processing.
#' @param FUN Function. Function to apply to each chunk.
#' @param n_cores Integer. Number of cores to use.
#' @param use_windows_future Logical. Whether to use future package on Windows
#' (default: TRUE).
#' @param ... Additional arguments passed to FUN.
#' @return List. Combined results from all chunks.
#' @details This function automatically detects the operating system and uses
#' appropriate parallel processing methods. On Windows, it uses the future
#' package; on Unix-like systems, it uses pbmcapply. If parallel processing
#' fails, it falls back to sequential execution.
#' @keywords internal
#' @examples
#' \dontrun{
#' chunks <- create_balanced_chunks(100, 4)
#' results <- safe_parallel_execution(
#'   chunks,
#'   function(chunk) lapply(chunk, sqrt),
#'   n_cores = 4
#' )
#' }
safe_parallel_execution <- function(chunks, FUN, n_cores,
                                    use_windows_future = TRUE, ...) {
  is_windows <- .Platform$OS.type == "windows"

  tryCatch({
    if (is_windows && use_windows_future) {
      if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
        stop("Windows parallel processing requires 'future' and 'future.apply' packages")
      }

      future::plan(future::multisession, workers = n_cores)
      on.exit(future::plan(future::sequential), add = TRUE)

      results <- future.apply::future_lapply(chunks, FUN, ...)
    } else {
      if (!requireNamespace("pbmcapply", quietly = TRUE)) {
        stop("Unix parallel processing requires 'pbmcapply' package")
      }

      results <- pbmcapply::pbmclapply(chunks, FUN, mc.cores = n_cores, ...)
    }

    return(unlist(results, recursive = FALSE))

  }, error = function(e) {
    warning("Parallel processing failed: ", e$message,
            ". Falling back to sequential processing.")

    # Flatten chunks and apply function sequentially
    all_items <- unlist(chunks, recursive = FALSE)
    return(lapply_with_bar(all_items, function(item) FUN(list(item), ...)[[1]]))
  })
}
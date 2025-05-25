\name{detectXOR}
\alias{detectXOR}
\title{Detect XOR Patterns in Classification Data}
\description{
Identifies XOR-shaped relationships between variables using statistical tests and pattern detection.
}
\usage{
detectXOR(
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
)
}
\arguments{
  \item{data}{Data frame containing features and class column}
  \item{class_col}{Name of class column (default: "class")}
  \item{check_tau}{Logical - compute classwise tau coefficients (default: TRUE)}
  \item{compute_axes_parallel_significance}{Logical - compute Wilcoxon tests (default: TRUE)}
  \item{p_threshold}{Significance threshold (default: 0.05)}
  \item{tau_threshold}{Tau coefficient threshold (default: 0.3)}
  \item{abs_diff_threshold}{Absolute difference threshold for patterns (default: 20)}
  \item{split_method}{Method for splitting data ("quantile" or "range") (default: "quantile")}
  \item{max_cores}{Maximum cores for parallel processing (default: NULL = automatic)}
  \item{extreme_handling}{Method for handling extreme values; options include "winsorize" or "none" (default: "winsorize")}
  \item{winsor_limits}{Numeric vector of length 2 specifying lower and upper quantiles for winsorization (default: c(0.05, 0.95))}
  \item{scale_data}{Logical; whether to scale/standardize the data before analysis (default: TRUE)}
  \item{use_complete}{Logical; whether to use only complete cases (default: TRUE)}
}
\value{
List containing:
  \item{results_df}{Data frame with detection results}
  \item{pair_list}{Full analysis details for each pair}
}
\examples{
\dontrun{
data <- data.frame(class = rep(letters[1:2], each = 50),
                   x = c(rnorm(50), rnorm(50, 2)),
                   y = c(rnorm(50), rnorm(50, 2)))
results <- detectXOR(data, class_col = "class")
}
}
\seealso{
\code{\link{generate_spaghetti_plot_from_results}} for visualization of results
}
\keyword{XOR}
\keyword{pattern}
\keyword{detection}

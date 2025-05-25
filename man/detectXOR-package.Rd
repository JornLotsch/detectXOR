\name{detectXOR-package}
\alias{detectXOR-package}
\docType{package}
\title{XOR Pattern Detection and Visualization}
\description{
Provides statistical detection and advanced visualization of XOR-shaped relationships in classification data. Implements parallel processing for large datasets and produces publication-quality ggplot2 visualizations.
}
\details{
\strong{Core Features:}
\enumerate{
  \item Statistical detection using chi-square tests and Kendall's tau
  \item Cross-platform parallel processing (future + pbmcapply)
  \item Interactive spaghetti plots and decision boundary visualizations
  \item Automatic pattern rotation detection
}

\strong{Main Functions:}
\itemize{
  \item \code{\link{detectXOR}}: Core detection algorithm
  \item \code{\link{generate_spaghetti_plot_from_results}}: Interactive line plots
  \item \code{\link{generate_xy_plot_from_results}}: Decision boundary plots
}
}
\author{
Jorn Lotsch <j.lotsch@em.uni-frankfurt.de>
}
\references{
Methodological foundations:
\itemize{
  \item Non-linear pattern detection in machine learning
  \item Statistical dependency measures (Kendall's tau)
  \item Parallel computing frameworks in R
}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/JornLotsch/detectXOR}
  \item Report bugs at \url{https://github.com/JornLotsch/detectXOR/issues}
}

Related packages:
\itemize{
  \item \code{\link[DescTools]{DescTools}} for statistical tests
  \item \code{\link[ggh4x]{ggh4x}} for advanced plotting
  \item \code{\link[future]{future}} for parallel processing
}
}
\examples{
\dontrun{
# Basic workflow with included dataset
data(XOR_data)

# Detect XOR patterns
results <- detectXOR(XOR_data, class_col = "class")

# Generate visualizations
generate_spaghetti_plot_from_results(
  results$results_df,
  XOR_data,
  class_col = "class"
)

generate_xy_plot_from_results(
  results$results_df,
  XOR_data,
  class_col = "class"
)
}
}
\keyword{package}
\keyword{XOR}
\keyword{machine learning}
\keyword{pattern detection}

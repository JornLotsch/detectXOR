\name{detectXOR-package}
\alias{detectXOR-package}
\docType{package}
\title{XOR Pattern Detection and Visualization}
\description{
Provides tools for detecting XOR-like patterns in variable pairs in 
    two-class data sets. Includes visualizations for pattern exploration and reporting
    capabilities with both text and HTML output formats.
}
\details{
\strong{Core Features:}
\enumerate{
  \item Statistical detection using chi-square tests and Kendall's tau
  \item Spaghetti plots and xy plot for pattern visualization
}
\strong{Main Functions:}
\itemize{
  \item \code{\link{detect_xor}}: Core detection algorithm
  \item \code{\link{generate_spaghetti_plot_from_results}}: Line plots
  \item \code{\link{generate_xy_plot_from_results}}: Plot for pattern visualization
}
}
\author{
Jorn Lotsch <j.lotsch@em.uni-frankfurt.de>
}
\references{
Methodological foundations:
\itemize{
  \item Pattern detection in machine learning
  \item Statistical dependency measures (Kendall's tau)
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
\donttest{
# Basic workflow with included dataset
data(XOR_data)

# Detect XOR patterns
results <- detect_xor(XOR_data, class_col = "class")

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

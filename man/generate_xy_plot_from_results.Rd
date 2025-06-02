\name{generate_xy_plot_from_results}
\alias{generate_xy_plot_from_results}
\title{Generate XOR Scatter Plots}
\description{
Creates scatterplots with decision boundaries for XOR patterns.
}
\usage{
generate_xy_plot_from_results(
  results,
  data,
  class_col,
  output_dir = ".",
  scale_data = TRUE,
  quantile_lines = c(1/3, 2/3),
  line_method = "quantile"
)
}
\arguments{
\item{results}{Detection results from detect_xor}
\item{data}{Original dataset}
\item{class_col}{Name of class column}
\item{output_dir}{Output directory path (default: ".")}
\item{scale_data}{Logical - scale variables before plotting (default: TRUE)}
\item{quantile_lines}{Quantiles for reference lines (default: c(1/3, 2/3))}
\item{line_method}{Boundary method ("quantile" or "range") (default: "quantile")}
}
\value{
Invisible ggplot object
}
\examples{
\dontrun{
generate_xy_plot_from_results(
  results = xor_results,
  data = my_data,
  class_col = "class"
)
}
}
\seealso{
\code{\link{detect_xor}} for detection function
}
\keyword{XOR}
\keyword{visualization}

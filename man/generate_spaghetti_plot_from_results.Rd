\name{generate_spaghetti_plot_from_results}
\alias{generate_spaghetti_plot_from_results}
\title{Generate XOR Spaghetti Plots}
\description{
Creates connected line plots for variables showing XOR patterns.
}
\usage{
generate_spaghetti_plot_from_results(
  results,
  data,
  class_col,
  output_dir = ".",
  scale_data = TRUE
)
}
\arguments{
\item{results}{Detection results from detect_xor}
\item{data}{Original dataset}
\item{class_col}{Name of class column}
\item{output_dir}{Output directory path (default: ".")}
\item{scale_data}{Logical - scale variables before plotting (default: TRUE)}
}
\value{
Invisible ggplot object
}
\examples{
\dontrun{
generate_spaghetti_plot_from_results(
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

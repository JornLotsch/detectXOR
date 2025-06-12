\name{generate_xy_plot_from_results}
\alias{generate_xy_plot_from_results}
\title{Generate XOR Scatter Plots}
\description{
Creates scatterplots with decision boundaries for variable pairs showing XOR patterns.
}
\usage{
generate_xy_plot_from_results(
  results,
  data,
  class_col,
  scale_data = TRUE,
  quantile_lines = c(1/3, 2/3),
  line_method = "quantile"
)
}
\arguments{
\item{results}{Either a data frame from \code{detect_xor()$results_df} or the full list object returned by \code{detect_xor()}}
\item{data}{Original dataset containing variables and classes}
\item{class_col}{Character string specifying the name of the class column}
\item{scale_data}{Logical indicating whether to scale variables before plotting (default: TRUE)}
\item{quantile_lines}{Numeric vector of length 2 specifying quantiles for reference lines (default: c(1/3, 2/3))}
\item{line_method}{Character string specifying the boundary calculation method, either "quantile" or "range" (default: "quantile")}
}
\value{
Returns a ggplot object. No files are saved automatically.
}
\details{
This function creates scatter plots for variable pairs that have been flagged as showing XOR patterns by \code{detect_xor()}. The plots include dashed reference lines that help visualize the decision boundaries used in XOR pattern detection.

The function automatically handles both original and rotated XOR patterns, applying the appropriate coordinate transformation when necessary. Variable pairs are separated using "||" as the delimiter in plot labels.

The \code{line_method} parameter controls how reference lines are calculated:
\itemize{
  \item "quantile": Lines are placed at the specified quantiles of the data distribution
  \item "range": Lines divide the data range into three equal parts
}

If no XOR patterns are detected, an empty plot with an appropriate message is returned.

To save the plot, use \code{ggplot2::ggsave()} or other standard R plotting save methods.
}
\examples{
\donttest{
# Using full results object (recommended)
data(XOR_data)
results <- detect_xor(data = XOR_data, class_col = "class")
xy_plot <- generate_xy_plot_from_results(
  results = results,
  data = XOR_data,
  class_col = "class"
)

# Display the plot
print(xy_plot)

# Using different boundary method
xy_plot_range <- generate_xy_plot_from_results(
  results = results,
  data = XOR_data,
  class_col = "class",
  line_method = "range"
)

# Save the plot if needed
# ggplot2::ggsave("my_xy_plot.png", xy_plot)

# Using extracted results_df (also works)
xy_plot_df <- generate_xy_plot_from_results(
  results = results$results_df,
  data = XOR_data,
  class_col = "class"
)
}
}
\seealso{
\code{\link{detect_xor}} for XOR pattern detection, \code{\link{generate_spaghetti_plot_from_results}} for spaghetti plots
}
\keyword{XOR}
\keyword{visualization}
\keyword{scatter}
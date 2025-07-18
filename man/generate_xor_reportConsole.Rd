% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_xor_reportConsole.R
\name{generate_xor_reportConsole}
\alias{generate_xor_reportConsole}
\title{Generate XOR Detection Report (Console-friendly)}
\description{
Creates a report with formatted table and plots for XOR pattern detection results.
}
\usage{
generate_xor_reportConsole(
  results,
  data,
  class_col,
  scale_data = TRUE,
  show_plots = TRUE,
  quantile_lines = c(1/3, 2/3),
  line_method = "quantile"
)
}
\arguments{
\item{results}{Either a data frame from \code{detect_xor$results_df} or the full list returned by \code{detect_xor}.}

\item{data}{Original dataset containing variables and classes.}

\item{class_col}{Character specifying the class column name.}

\item{scale_data}{Logical indicating whether to scale variables in plots. Default: TRUE.}

\item{show_plots}{Logical indicating whether to display plots. Default: TRUE.}

\item{quantile_lines}{Numeric vector of quantiles for reference lines in XY plots. Default: c(1/3, 2/3).}

\item{line_method}{Method for boundary calculation ("quantile" or "range"). Default: "quantile".}
}
\value{
Invisibly returns a list containing the formatted table and plots (if generated).
}
\seealso{
\code{\link{detect_xor}} for XOR pattern detection,
\code{\link{generate_xor_reportHTML}} for HTML report generation
}
\name{generate_spaghetti_plot_from_results}
 \alias{generate_spaghetti_plot_from_results}
 \title{Generate XOR Spaghetti Plots}
 \description{
 Creates connected line plots for variable pairs showing XOR patterns.
 }
 \usage{
 generate_spaghetti_plot_from_results(
   results,
   data,
   class_col,
   scale_data = TRUE
 )
 }
 \arguments{
 \item{results}{Either a data frame from \code{detect_xor()$results_df} or the full list object returned by \code{detect_xor()}}
 \item{data}{Original dataset containing variables and classes}
 \item{class_col}{Character string specifying the name of the class column}
 \item{scale_data}{Logical indicating whether to scale variables before plotting (default: TRUE)}
 }
 \value{
 Returns a ggplot object. No files are saved automatically.
 }
 \details{
 This function creates spaghetti plots (connected line plots) for variable pairs that have been flagged as showing XOR patterns by \code{detect_xor()}. The function automatically handles both original and rotated XOR patterns, applying the appropriate coordinate transformation when necessary.

 The function accepts either the full results object returned by \code{detect_xor()} or just the \code{results_df} component extracted from it. Variable pairs are separated using "||" as the delimiter in plot labels.

 If no XOR patterns are detected, an empty plot with an appropriate message is returned.

 To save the plot, use \code{ggplot2::ggsave()} or other standard R plotting save methods.
 }
 \examples{
 \donttest{
 # Using full results object (recommended)
 data(XOR_data)
 results <- detect_xor(data = XOR_data, class_col = "class")
 spaghetti_plot <- generate_spaghetti_plot_from_results(
   results = results,
   data = XOR_data,
   class_col = "class"
 )

 # Display the plot
 print(spaghetti_plot)

 # Save the plot if needed
 # ggplot2::ggsave("my_spaghetti_plot.png", spaghetti_plot)

 # Using extracted results_df (also works)
 xy_plot <- generate_spaghetti_plot_from_results(
   results = results$results_df,
   data = XOR_data,
   class_col = "class"
 )
 }
 }
 \seealso{
 \code{\link{detect_xor}} for XOR pattern detection, \code{\link{generate_xy_plot_from_results}} for scatter plots
 }
 \keyword{XOR}
 \keyword{visualization}
 \keyword{spaghetti}
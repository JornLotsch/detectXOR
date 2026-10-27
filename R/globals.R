# Declare global variables for non-standard evaluation
# Suppresses R CMD check notes about variables used in data manipulation and plotting expressions (e.g., with dplyr or ggplot2), as well as imported functions referenced without explicit namespace.
# This is necessary for variables and functions such as \code{ID}, \code{xor_shape_detected}, \code{xor_pattern}, \code{future.apply}, and \code{pbmcapply} that are used in pipelines or evaluated indirectly.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "ID", "xor_shape_detected", "xor_pattern",
    "pair_id", "variable", "value", "x", "y",
    "x_q1", "x_q2", "y_q1", "y_q2", "xor_purity", "p_adjusted",
    "var1_var2", "var1", "var2", "pair_key", "stat_label"
  ))
}
# Simulate a two-class data set with an XOR pair (X, Y) and a noise variable Z
simulate_xor <- function(n = 400, label_noise = 0, rotated = FALSE) {
  a <- stats::rbinom(n, 1, 0.5)
  b <- stats::rbinom(n, 1, 0.5)
  cls <- ifelse(xor(a, b), 2, 1)
  flip <- stats::runif(n) < label_noise
  cls[flip] <- 3 - cls[flip]
  x <- stats::rnorm(n, ifelse(a == 1, 10, 3))
  y <- stats::rnorm(n, ifelse(b == 1, 10, 3))
  if (rotated) {
    xy <- cbind(x + y, x - y)
    x <- xy[, 1]
    y <- xy[, 2]
  }
  data.frame(X = x, Y = y, Z = stats::rnorm(n), class = cls)
}

run_xor <- function(data, ...) {
  suppressMessages(detect_xor(data, n_perm = 199, ...))
}

detected_pairs <- function(result) {
  df <- result$results_df
  paste(df$var1, df$var2, sep = "||")[df$xor_shape_detected]
}

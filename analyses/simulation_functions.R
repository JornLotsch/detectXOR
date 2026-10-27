# Data generators for the simulation studies (scripts 04-06).
#
# Each data set contains a target pair (X, Y) with a given structure and,
# optionally, noise variables that are unrelated to the class. Variables are
# built from latent binary states (low/high) plus an error term, so that XOR
# and non-XOR structures have comparable marginal distributions.

SCENARIOS <- c(
  "null", "xor", "xor_rotated", "xor_noise_10", "xor_noise_20", "xor_noise_30",
  "incomplete_xor", "and_interaction", "nonmonotone_single", "linear_diagonal", "main_effect"
)
XOR_SCENARIOS <- c("xor", "xor_rotated", "xor_noise_10", "xor_noise_20", "xor_noise_30")
DISTRIBUTIONS <- c("normal", "uniform", "skewed", "binary")

# Error term with mean 0 and standard deviation 1
draw_error <- function(n, distribution) {
  switch(distribution,
    normal = stats::rnorm(n),
    uniform = stats::runif(n, -sqrt(3), sqrt(3)),
    skewed = {
      e <- stats::rlnorm(n, sdlog = 0.75)
      (e - exp(0.75^2 / 2)) / sqrt((exp(0.75^2) - 1) * exp(0.75^2))
    },
    binary = rep(0, n)
  )
}

# Variable from latent binary states: separation between modes in SD units
make_variable <- function(state, distribution, separation) {
  if (distribution == "binary") return(as.numeric(state))
  separation * state + draw_error(length(state), distribution)
}

# Class labels with a given proportion of class 1
make_classes <- function(n, prop_class1) {
  n1 <- round(n * prop_class1)
  sample(c(rep(1, n1), rep(2, n - n1)))
}

# Latent states (a, b) for an XOR pattern: class 1 in (0,0)/(1,1), class 2 in (0,1)/(1,0)
xor_states <- function(cls) {
  a <- stats::rbinom(length(cls), 1, 0.5)
  b <- ifelse(cls == 1, a, 1 - a)
  list(a = a, b = b)
}

#' Simulate one data set
#' @param scenario One of SCENARIOS.
#' @param n Sample size.
#' @param distribution One of DISTRIBUTIONS.
#' @param prop_class1 Proportion of class 1.
#' @param n_noise Number of additional noise variables.
#' @param separation Distance between the latent modes in SD units.
simulate_dataset <- function(scenario, n, distribution = "normal", prop_class1 = 0.5,
                             n_noise = 8, separation = 3) {
  cls <- make_classes(n, prop_class1)
  noise_rate <- switch(scenario, xor_noise_10 = 0.1, xor_noise_20 = 0.2, xor_noise_30 = 0.3, 0)

  if (scenario %in% c("xor", "xor_rotated", "xor_noise_10", "xor_noise_20", "xor_noise_30")) {
    s <- xor_states(cls)
    x <- make_variable(s$a, distribution, separation)
    y <- make_variable(s$b, distribution, separation)
    if (scenario == "xor_rotated") {
      xy <- cbind(x + y, x - y)
      x <- xy[, 1]
      y <- xy[, 2]
    }
    # Label noise: a share of the cases receives the wrong class
    flip <- stats::runif(n) < noise_rate
    cls[flip] <- 3 - cls[flip]
  } else if (scenario == "incomplete_xor") {
    # XOR pattern with one empty quadrant: class 2 only in (0, 1)
    s <- xor_states(cls)
    s$a[cls == 2] <- 0
    s$b[cls == 2] <- 1
    x <- make_variable(s$a, distribution, separation)
    y <- make_variable(s$b, distribution, separation)
  } else if (scenario == "and_interaction") {
    # Interaction that is not XOR: class 2 if both states are high
    a <- stats::rbinom(n, 1, 0.5)
    b <- stats::rbinom(n, 1, 0.5)
    cls <- ifelse(a == 1 & b == 1, 2, 1)
    x <- make_variable(a, distribution, separation)
    y <- make_variable(b, distribution, separation)
  } else if (scenario == "nonmonotone_single") {
    # Class depends non-monotonically on X alone (class 2 in the middle);
    # Y has the same bimodal marginal distribution but is unrelated to the class
    a <- stats::rbinom(n, 1, 0.5)
    b <- stats::rbinom(n, 1, 0.5)
    cls <- ifelse(a != b, 2, 1)
    x <- make_variable(a, distribution, separation) + make_variable(b, distribution, separation)
    y <- make_variable(stats::rbinom(n, 1, 0.5), distribution, separation)
  } else if (scenario == "linear_diagonal") {
    # Classes separated along the diagonal x + y
    a <- stats::rbinom(n, 1, 0.5)
    b <- stats::rbinom(n, 1, 0.5)
    x <- make_variable(a, distribution, separation)
    y <- make_variable(b, distribution, separation)
    score <- x + y + stats::rnorm(n, sd = 0.5)
    cls <- ifelse(score > stats::median(score), 2, 1)
  } else if (scenario == "main_effect") {
    # Class difference in X only
    x <- make_variable(cls == 2, distribution, separation)
    y <- make_variable(stats::rbinom(n, 1, 0.5), distribution, separation)
  } else if (scenario == "null") {
    x <- make_variable(stats::rbinom(n, 1, 0.5), distribution, separation)
    y <- make_variable(stats::rbinom(n, 1, 0.5), distribution, separation)
  } else {
    stop("Unknown scenario: ", scenario)
  }

  data <- data.frame(class = cls, X = x, Y = y)
  for (j in seq_len(n_noise)) {
    data[[paste0("N", j)]] <- make_variable(stats::rbinom(n, 1, 0.5), distribution, separation)
  }
  data
}

# Summarize a detect_xor result for the target pair X-Y and the noise pairs
summarize_detection <- function(result_df) {
  target <- result_df$var1 == "X" & result_df$var2 == "Y"
  c(
    target_detected = any(result_df$xor_shape_detected[target]),
    target_significant = any(result_df$p_adjusted[target] < 0.05),
    other_detected = any(result_df$xor_shape_detected[!target]),
    any_detected = any(result_df$xor_shape_detected),
    any_significant = any(result_df$p_adjusted < 0.05)
  )
}

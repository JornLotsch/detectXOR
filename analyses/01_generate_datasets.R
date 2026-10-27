# Generation of the example data sets of the manuscript.
# The code reproduces the data generation used for the original submission
# (same random seeds and sampling order); plotting code was removed.
# Output: data/example_datasets.rds

source("00_setup.R")

datasets <- list()

# Common parameters
n <- 400 # Total number of samples
n1 <- 200 # Samples in class 1
n2 <- n - n1 # Samples in class 2
class <- c(rep(1, n1), rep(2, n2)) # Class labels

span <- 5 # Range for uniform sampling
A1 <- 1
B1 <- A1 + span # Range for class 1
A2 <- B1 + 1
B2 <- A2 + span # Range for class 2

mean_XOR1 <- 3
mean_XOR2 <- mean_XOR1 + 7
sd_XOR <- 1
half_n1 <- n1 %/% 2
half_n2 <- n1 - half_n1

# Generate the introductory data set; delta_A and delta_B are the class
# differences of the non-XOR variables A and B
generate_intro_data <- function(delta_A, delta_B) {
  # Variable A: Two normal distributions
  mean_A1 <- 3
  sd_A1 <- 1
  Variable_A <- c(rnorm(n1, mean = mean_A1, sd = sd_A1), rnorm(n2, mean = mean_A1 + delta_A, sd = sd_A1))

  # Variable B: Two different normal distributions
  mean_B1 <- 1
  sd_B1 <- 3
  Variable_B <- c(rnorm(n1, mean = mean_B1, sd = sd_B1), rnorm(n2, mean = mean_B1 + delta_B, sd = sd_B1))

  # Variables C and D (XOR problem simulation)
  dist_X1_group1 <- c(rnorm(half_n1, mean = mean_XOR1, sd = sd_XOR), rnorm(half_n1, mean = mean_XOR2, sd = sd_XOR))
  dist_Y1_group1 <- c(rnorm(half_n1, mean = mean_XOR1, sd = sd_XOR), rnorm(half_n1, mean = mean_XOR2, sd = sd_XOR))
  dist_X2_group2 <- c(rnorm(half_n2, mean = mean_XOR2, sd = sd_XOR), rnorm(half_n2, mean = mean_XOR1, sd = sd_XOR))
  dist_Y2_group2 <- c(rnorm(half_n2, mean = mean_XOR1, sd = sd_XOR), rnorm(half_n2, mean = mean_XOR2, sd = sd_XOR))

  # Variable E: Random uniform distribution
  Variable_E <- runif(n, min = 1, max = 10)

  data.frame(
    class = class,
    Variable_A = Variable_A,
    Variable_B = Variable_B,
    Variable_C = c(dist_X1_group1, dist_X2_group2),
    Variable_D = c(dist_Y1_group1, dist_Y2_group2),
    Variable_E = Variable_E
  )
}

# XOR_data: variables A and B with narrowly significant class differences
set.seed(123)
datasets$XOR_data <- generate_intro_data(delta_A = 0.25, delta_B = -0.7)

# XOR_data_nonsig: no significant class differences in any variable
set.seed(123)
datasets$XOR_data_nonsig <- generate_intro_data(delta_A = 0, delta_B = 0)

# Additional XOR variables F and G: uniform sampling in range mean +/- 1 SD
set.seed(123)
dist_X1_uniform_group1 <- c(runif(half_n1, min = mean_XOR1 - sd_XOR, max = mean_XOR1 + sd_XOR),
                            runif(half_n1, min = mean_XOR2 - sd_XOR, max = mean_XOR2 + sd_XOR))
dist_Y1_uniform_group1 <- c(runif(half_n1, min = mean_XOR1 - sd_XOR, max = mean_XOR1 + sd_XOR),
                            runif(half_n1, min = mean_XOR2 - sd_XOR, max = mean_XOR2 + sd_XOR))
dist_X2_uniform_group2 <- c(runif(half_n2, min = mean_XOR2 - sd_XOR, max = mean_XOR2 + sd_XOR),
                            runif(half_n2, min = mean_XOR1 - sd_XOR, max = mean_XOR1 + sd_XOR))
dist_Y2_uniform_group2 <- c(runif(half_n2, min = mean_XOR1 - sd_XOR, max = mean_XOR1 + sd_XOR),
                            runif(half_n2, min = mean_XOR2 - sd_XOR, max = mean_XOR2 + sd_XOR))
XOR_data_extended <- cbind.data.frame(
  datasets$XOR_data,
  Variable_F = c(dist_X1_uniform_group1, dist_X2_uniform_group2),
  Variable_G = c(dist_Y1_uniform_group1, dist_Y2_uniform_group2)
)

# Rotated XOR variables H and I
set.seed(123)
H1 <- runif(n1 / 2, min = A1, max = B1)
H2 <- runif(n2 / 2, min = A2, max = B2)
I1 <- runif(n1 / 2, min = A1, max = B1)
I2 <- runif(n2 / 2, min = A2, max = B2)
Ha <- c(H1, H2)
Ia <- c(I1, I2)
H1 <- runif(n1 / 2, min = A1, max = B1)
H2 <- runif(n2 / 2, min = A2, max = B2)
I1 <- runif(n1 / 2, min = A1, max = B1)
I2 <- runif(n2 / 2, min = A2, max = B2)
Hb <- c(H1, H2)
Ib <- c(I2, I1) # reversed I variables for XOR structure
Hx <- c(Ha, Hb)
Ix <- c(Ia, Ib)
datasets$XOR_data_plus_rotated <- cbind.data.frame(XOR_data_extended, Variable_H = Hx + Ix, Variable_I = Hx - Ix)

# Incomplete XOR pattern: first 250 cases of XOR_data
datasets$pseudoXOR_data_first250 <- head(datasets$XOR_data, 250)

# Simulated genetic data (dominant hereditary model)
set.seed(123)
n_total <- 1000
pA1 <- 0.55; pA2 <- 0.45
pB1 <- 0.55; pB2 <- 0.45
pC <- 0.5; pD <- 0.5
pE <- 0.5
genotype <- function(n, p) {
  q <- 1 - p
  sample(0:2, n, replace = TRUE, prob = c(p^2, 2 * p * q, q^2))
}
to_dominant <- function(g) ifelse(g >= 1, 1, 0)
A_all <- genotype(n_total, pA1)
B_all <- genotype(n_total, pB1)
C_all <- genotype(n_total, pC)
D_all <- genotype(n_total, pD)
E_all <- genotype(n_total, pE)
class_all <- ifelse(to_dominant(C_all) == to_dominant(D_all), 1, 2)
set.seed(5678)
idx1 <- sample(which(class_all == 1), n1)
idx2 <- sample(which(class_all == 2), n2)
selected_idx <- c(idx1, idx2)
A_selected <- c(genotype(n1, pA1), genotype(n2, pA2))
B_selected <- c(genotype(n1, pB1), genotype(n2, pB2))
datasets$XOR_genetic_data <- data.frame(
  class = c(rep(1, n1), rep(2, n2)),
  Variable_A = to_dominant(A_selected),
  Variable_B = to_dominant(B_selected),
  Variable_C = to_dominant(C_all[selected_idx]),
  Variable_D = to_dominant(D_all[selected_idx]),
  Variable_E = to_dominant(E_all[selected_idx])
)

# Sliding means data: 100 XOR pairs with shifting means and increasing SD
set.seed(123)
XOR_sliding_means_data <- data.frame(class = factor(c(rep(1, n1), rep(2, n1))))
XOR_sliding_means_params <- data.frame()
m1_values <- seq(3, 10, length.out = 10)
m2_values <- seq(10, 3, length.out = 10)
std_devs <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
var_index <- 1
for (sd_orig in std_devs) {
  for (i in seq_along(m1_values)) {
    m1 <- jitter(m1_values[i], factor = 1)
    m2 <- jitter(m2_values[i], factor = 1)
    sd_current <- jitter(sd_orig, factor = 1)
    nn1 <- n1 %/% 2
    A_class1 <- c(rnorm(nn1, mean = m1, sd = sd_current), rnorm(nn1, mean = m2, sd = sd_current))
    B_class1 <- c(rnorm(nn1, mean = m1, sd = sd_current), rnorm(nn1, mean = m2, sd = sd_current))
    A_class2 <- c(rnorm(nn1, mean = m2, sd = sd_current), rnorm(nn1, mean = m1, sd = sd_current))
    B_class2 <- c(rnorm(nn1, mean = m1, sd = sd_current), rnorm(nn1, mean = m2, sd = sd_current))
    XOR_sliding_means_data[[paste0("A", var_index)]] <- c(A_class1, A_class2)
    XOR_sliding_means_data[[paste0("B", var_index)]] <- c(B_class1, B_class2)
    XOR_sliding_means_params <- rbind(XOR_sliding_means_params, data.frame(
      var1 = paste0("A", var_index), var2 = paste0("B", var_index),
      m1 = m1, m2 = m2, S = sd_current
    ))
    var_index <- var_index + 1
  }
}
datasets$XOR_sliding_means_data <- XOR_sliding_means_data
datasets$XOR_sliding_means_params <- XOR_sliding_means_params

# The introductory data set must be identical to the one shipped with the package
data("XOR_data", package = "detectXOR", envir = environment())
stopifnot(isTRUE(all.equal(datasets$XOR_data, XOR_data, check.attributes = FALSE)))

saveRDS(datasets, file.path(DATA_DIR, "example_datasets.rds"))
message("Saved ", length(datasets), " objects to ", file.path(DATA_DIR, "example_datasets.rds"))
write_session_info("01_generate_datasets")

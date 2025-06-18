# Load libraries
library(ggplot2)
library(cowplot)
library(dplyr)
library(ggthemes)
library(grid)

# Plot functions (custom for XOR_data and XOR_data_nonsig)
source("/home/joern/.Datenplatte/Joerns Dateien/Aktuell/FeatureSelectionInformationTheory/08AnalyseProgramme/R/XOR_plot_functions5.R")

# Seed
seed <- 123

# PARAMETERS
n <- 400
n1 <- 200
n2 <- n - n1
class <- c(rep(1, n1), rep(2, n2))
span <- 5
A1 <- 1; B1 <- A1 + span
A2 <- B1 + 1; B2 <- A2 + span

##############################################
# 1. XOR dataset
##############################################

set.seed(seed)
delta_A <- 0.25; mean_A1 <- 3; sd_A1 <- 1
mean_A2 <- mean_A1 + delta_A
Variable_A <- c(rnorm(n1, mean = mean_A1, sd = sd_A1), rnorm(n2, mean = mean_A2, sd = sd_A1))

delta_B <- -0.7; mean_B1 <- 1; sd_B1 <- 3
mean_B2 <- mean_B1 + delta_B
Variable_B <- c(rnorm(n1, mean = mean_B1, sd = sd_B1), rnorm(n2, mean = mean_B2, sd = sd_B1))

mean_XOR1 <- 3; mean_XOR2 <- mean_XOR1 + 7; sd_XOR <- 1
half_n1 <- n1 %/% 2; half_n2 <- n1 - half_n1
dist_X1_group1 <- c(rnorm(half_n1, mean = mean_XOR1, sd = sd_XOR), rnorm(half_n1, mean = mean_XOR2, sd = sd_XOR))
dist_Y1_group1 <- c(rnorm(half_n1, mean = mean_XOR1, sd = sd_XOR), rnorm(half_n1, mean = mean_XOR2, sd = sd_XOR))
dist_X2_group2 <- c(rnorm(half_n2, mean = mean_XOR2, sd = sd_XOR), rnorm(half_n2, mean = mean_XOR1, sd = sd_XOR))
dist_Y2_group2 <- c(rnorm(half_n2, mean = mean_XOR1, sd = sd_XOR), rnorm(half_n2, mean = mean_XOR2, sd = sd_XOR))
Variable_C <- c(dist_X1_group1, dist_X2_group2)
Variable_D <- c(dist_Y1_group1, dist_Y2_group2)
Variable_E <- runif(n, min = 1, max = 10)

XOR_data <- data.frame(
  class = class,
  Variable_A = Variable_A,
  Variable_B = Variable_B,
  Variable_C = Variable_C,
  Variable_D = Variable_D,
  Variable_E = Variable_E
)

# Custom plot for XOR_data
initial_check_plot_XOR_data <-
  cowplot::plot_grid(
    XOR_raw_plot(data = XOR_data[, -1], Cls = as.factor(XOR_data$class)),
    grid::grid.grabExpr(print(
      XOR_matrix_plot(XOR_data[, -1], Cls = as.factor(XOR_data$class), data_name = "actual supervised processed XOR_data")
    )),
    ncol = 1,
    rel_heights = c(1, 3),
    align = "v", axis = "lr",
    labels = "AUTO"
  )
print(initial_check_plot_XOR_data)

##############################################
# 2. XOR nonsignificant dataset
##############################################

set.seed(seed)
delta_A <- 0; mean_A1 <- 3; sd_A1 <- 1
mean_A2 <- mean_A1 + delta_A
Variable_A <- c(rnorm(n1, mean = mean_A1, sd = sd_A1), rnorm(n2, mean = mean_A2, sd = sd_A1))

delta_B <- 0; mean_B1 <- 1; sd_B1 <- 3
mean_B2 <- mean_B1 + delta_B
Variable_B <- c(rnorm(n1, mean = mean_B1, sd = sd_B1), rnorm(n2, mean = mean_B2, sd = sd_B1))

dist_X1_group1 <- c(rnorm(half_n1, mean = mean_XOR1, sd = sd_XOR), rnorm(half_n1, mean = mean_XOR2, sd = sd_XOR))
dist_Y1_group1 <- c(rnorm(half_n1, mean = mean_XOR1, sd = sd_XOR), rnorm(half_n1, mean = mean_XOR2, sd = sd_XOR))
dist_X2_group2 <- c(rnorm(half_n2, mean = mean_XOR2, sd = sd_XOR), rnorm(half_n2, mean = mean_XOR1, sd = sd_XOR))
dist_Y2_group2 <- c(rnorm(half_n2, mean = mean_XOR1, sd = sd_XOR), rnorm(half_n2, mean = mean_XOR2, sd = sd_XOR))
Variable_C <- c(dist_X1_group1, dist_X2_group2)
Variable_D <- c(dist_Y1_group1, dist_Y2_group2)
Variable_E <- runif(n, min = 1, max = 10)

XOR_data_nonsig <- data.frame(
  class = class,
  Variable_A = Variable_A,
  Variable_B = Variable_B,
  Variable_C = Variable_C,
  Variable_D = Variable_D,
  Variable_E = Variable_E
)

# Custom plot for XOR_data_nonsig
initial_check_plot_XOR_data_nonsig <-
  cowplot::plot_grid(
    XOR_raw_plot(data = XOR_data_nonsig[, -1], Cls = XOR_data_nonsig$class),
    grid::grid.grabExpr(print(
      XOR_matrix_plot(XOR_data_nonsig[, -1], Cls = as.factor(XOR_data_nonsig$class), data_name = "actual supervised processed XOR_data_nonsig")
    )),
    ncol = 1,
    rel_heights = c(1, 3),
    align = "v", axis = "lr",
    labels = "AUTO"
  )
print(initial_check_plot_XOR_data_nonsig)

##############################################
# 3. Rotated XOR dataset
##############################################

set.seed(seed)
# If XOR_data_extended is not defined, use XOR_data as base
if (!exists("XOR_data_extended")) XOR_data_extended <- XOR_data

# Generate rotated variables
set.seed(123)
H1 <- runif(n1 / 2, min = A1, max = B1)
H2 <- runif(n2 / 2, min = A2, max = B2)
I1 <- runif(n1 / 2, min = A1, max = B1)
I2 <- runif(n2 / 2, min = A2, max = B2)
Ha <- c(H1, H2); Ia <- c(I1, I2)
H1 <- runif(n1 / 2, min = A1, max = B1)
H2 <- runif(n2 / 2, min = A2, max = B2)
I1 <- runif(n1 / 2, min = A1, max = B1)
I2 <- runif(n2 / 2, min = A2, max = B2)
Hb <- c(H1, H2); Ib <- c(I2, I1)
Hx <- c(Ha, Hb); Ix <- c(Ia, Ib)
H <- Hx + Ix; I <- Hx - Ix

XOR_data_plus_rotated <- cbind.data.frame(XOR_data_extended, Variable_H = H, Variable_I = I)

initial_check_plot_XOR_data_plus_rotated <-
  cowplot::plot_grid(
    XOR_raw_plot(data = XOR_data_plus_rotated[, -1], Cls = as.factor(XOR_data_plus_rotated$class)),
    grid::grid.grabExpr(print(
      XOR_matrix_plot(XOR_data_plus_rotated[, -1], Cls = as.factor(XOR_data_plus_rotated$class), data_name = "actual supervised processed XOR_data_plus_rotated")
    )),
    ncol = 1,
    rel_heights = c(1, 3),
    align = "v", axis = "lr",
    labels = "AUTO"
  )
print(initial_check_plot_XOR_data_plus_rotated)

##############################################
# 4. XOR-like non-modal (uniform) dataset
##############################################

set.seed(seed)
A <- runif(n, min = -10, max = 10)
B <- runif(n, min = -10, max = 10)
class <- ifelse((A > 0 & B <= 0) | (A <= 0 & B > 0), 2, 1)
C <- rnorm(n, mean = 5, sd = 2)
D <- runif(n, min = 0, max = 20)
E <- rnorm(n, mean = -3, sd = 1)
XOR_uniform_data <- data.frame(class = class, Variable_A = A, Variable_B = B, Variable_C = C, Variable_D = D, Variable_E = E)

initial_check_plot_XOR_uniform_data <-
  cowplot::plot_grid(
    XOR_raw_plot(data = XOR_uniform_data[, -1], Cls = as.factor(XOR_uniform_data$class)),
    grid::grid.grabExpr(print(
      XOR_matrix_plot(XOR_uniform_data[, -1], Cls = as.factor(XOR_uniform_data$class), data_name = "actual supervised processed XOR_uniform_data")
    )),
    ncol = 1,
    rel_heights = c(1, 3),
    align = "v", axis = "lr",
    labels = "AUTO"
  )
print(initial_check_plot_XOR_uniform_data)

##############################################
# 4. XOR dataset only 1st 250 cases
##############################################

pseudoXOR_data_first250 <- head(XOR_data, 250)

##############################################
# 6. XOR linear shapes
##############################################

set.seed(123)

# Set parameters
n <- 400
n1 <- 200
nn1 <- n1 %/% 2
jitter_amount <- 0.2

# Helper to generate a block of data
make_block <- function(A, B) {
  data.frame(A = A + runif(length(A), -jitter_amount, jitter_amount),
             B = B + runif(length(B), -jitter_amount, jitter_amount))
}

# Generate all blocks
block1 <- make_block(
  c(seq(10, 1, length.out = nn1), seq(1, 10, length.out = nn1), seq(1, 10, length.out = nn1), seq(10, 1, length.out = nn1)),
  c(seq(1, 10, length.out = nn1), seq(10, 1, length.out = nn1), seq(1, 10, length.out = nn1), seq(10, 1, length.out = nn1))
)
block2 <- make_block(
  c(seq(5, 0, length.out = nn1), seq(10, 5, length.out = nn1), seq(10, 5, length.out = nn1), seq(5, 0, length.out = nn1)),
  c(seq(0, 5, length.out = nn1), seq(5, 10, length.out = nn1), seq(0, 5, length.out = nn1), seq(5, 10, length.out = nn1))
)
block3 <- make_block(
  c(seq(5, 0, length.out = nn1), seq(0, 5, length.out = nn1), seq(10, 5, length.out = nn1), seq(10, 5, length.out = nn1)),
  c(seq(0, 5, length.out = nn1), seq(10, 5, length.out = nn1), seq(0, 5, length.out = nn1), seq(5, 10, length.out = nn1))
)
block4 <- make_block(
  c(seq(5, 0, length.out = nn1), seq(10, 5, length.out = nn1), seq(5, 10, length.out = nn1), seq(0, 5, length.out = nn1)),
  c(seq(0, 5, length.out = nn1), seq(5, 10, length.out = nn1), seq(0, 5, length.out = nn1), seq(5, 10, length.out = nn1))
)
block5 <- make_block(
  c(seq(0, 5, length.out = nn1), seq(10, 5, length.out = nn1), seq(5, 10, length.out = nn1), seq(0, 5, length.out = nn1)),
  c(seq(0, 5, length.out = nn1), seq(5, 10, length.out = nn1), seq(5, 0, length.out = nn1), seq(5, 10, length.out = nn1))
)
block6 <- make_block(
  c(seq(0, 5, length.out = nn1), seq(5, 10, length.out = nn1), seq(10, 5, length.out = nn1), seq(0, 5, length.out = nn1)),
  c(seq(0, 5, length.out = nn1), seq(5, 0, length.out = nn1), seq(5, 10, length.out = nn1), seq(5, 10, length.out = nn1))
)
block7 <- make_block(
  c(seq(0, 5, length.out = nn1), seq(5, 10, length.out = nn1), seq(10, 5, length.out = nn1), seq(10, 5, length.out = nn1)),
  c(seq(0, 5, length.out = nn1), seq(5, 10, length.out = nn1), seq(0, 5, length.out = nn1), seq(5, 10, length.out = nn1))
)
block8 <- make_block(
  c(seq(0, 5, length.out = nn1), seq(5, 10, length.out = nn1), seq(10, 5, length.out = nn1), seq(10, 5, length.out = nn1)),
  c(seq(0, 5, length.out = nn1), seq(5, 10, length.out = nn1), seq(5, 10, length.out = nn1), seq(5, 10, length.out = nn1))
)

# Combine blocks into a single data frame
XOR_data_linear <- data.frame(
  class = as.factor(rep(1:2, each = n/2)),
  Variable_A = block1$A,
  Variable_B = block1$B,
  Variable_C = block2$A,
  Variable_D = block2$B,
  Variable_E = block3$A,
  Variable_F = block3$B,
  Variable_G = block4$A,
  Variable_H = block4$B,
  Variable_I = block5$A,
  Variable_J = block5$B,
  Variable_K = block6$A,
  Variable_L = block6$B,
  Variable_M = block7$A,
  Variable_N = block7$B,
  Variable_O = block8$A,
  Variable_P = block8$B
)

# Standard plot for XOR_data_linear
initial_check_plot_XOR_data_linear <-
  cowplot::plot_grid(
    XOR_raw_plot(data = XOR_data_linear[, -1], Cls = as.factor(XOR_data_linear$class)),
    grid::grid.grabExpr(print(
      XOR_matrix_plot(XOR_data_linear[, -1], Cls = as.factor(XOR_data_linear$class), data_name = "actual supervised processed XOR_data_linear")
    )),
    ncol = 1,
    rel_heights = c(1, 3),
    align = "v", axis = "lr",
    labels = "AUTO"
  )
print(initial_check_plot_XOR_data_linear)


##############################################
# 7. XOR pseudo genetic dataset
##############################################
set.seed(123)

n_total <- 1000  # Start with more cases to ensure enough per class
n1 <- 200
n2 <- 200

# Allele frequencies (choose as desired, but here we use p = 0.5 for C and D for simplicity)
pA1 <- 0.55; pA2 <- 0.45  # Slightly different for A
pB1 <- 0.55; pB2 <- 0.45  # Slightly different for B
pC <- 0.5; pD <- 0.5      # Same for C and D (could be different if you wish)
pE <- 0.5                 # Same for E

# Function to generate genotype (0,1,2) from allele frequency p under HWE
genotype <- function(n, p) {
  q <- 1 - p
  sample(0:2, n, replace = TRUE, prob = c(p^2, 2*p*q, q^2))
}

# Function to transform genotype to binary using dominant model
to_dominant <- function(g) ifelse(g >= 1, 1, 0)

# Generate all variables for all cases (genotype, not binary)
A_all <- genotype(n_total, pA1)  # Will be overwritten for class 2
B_all <- genotype(n_total, pB1)  # Will be overwritten for class 2
C_all <- genotype(n_total, pC)
D_all <- genotype(n_total, pD)
E_all <- genotype(n_total, pE)

# Transform C and D to binary (dominant model) for class assignment
C_bin_all <- to_dominant(C_all)
D_bin_all <- to_dominant(D_all)

# Assign class by XOR on (C_bin, D_bin)
class_all <- ifelse(C_bin_all == D_bin_all, 1, 2)

# Sample exactly n1 and n2 cases from each class
idx1_all <- which(class_all == 1)
idx2_all <- which(class_all == 2)

if (length(idx1_all) < n1 || length(idx2_all) < n2) {
  stop("Not enough cases in one of the classes. Increase n_total or adjust allele frequencies.")
}

set.seed(5678)  # For reproducibility in sampling
idx1 <- sample(idx1_all, n1)
idx2 <- sample(idx2_all, n2)
selected_idx <- c(idx1, idx2)

# Generate A and B with class-specific allele frequencies for selected cases
A_selected <- c(genotype(n1, pA1), genotype(n2, pA2))
B_selected <- c(genotype(n1, pB1), genotype(n2, pB2))

# Subset C, D, E to the selected cases (still as genotypes)
C_selected <- C_all[selected_idx]
D_selected <- D_all[selected_idx]
E_selected <- E_all[selected_idx]

# Transform all to binary (dominant model)
A_bin <- to_dominant(A_selected)
B_bin <- to_dominant(B_selected)
C_bin <- to_dominant(C_selected)
D_bin <- to_dominant(D_selected)
E_bin <- to_dominant(E_selected)

# Construct final data frame
XOR_genetic_data <- data.frame(
  class = c(rep(1, n1), rep(2, n2)),
  Variable_A = A_bin,
  Variable_B = B_bin,
  Variable_C = C_bin,
  Variable_D = D_bin,
  Variable_E = E_bin
)

# Genotype version (before dominant transformation)
XOR_genetic_data_genotype <- data.frame(
  class = c(rep(1, n1), rep(2, n2)),
  Variable_A = A_selected,
  Variable_B = B_selected,
  Variable_C = C_selected,
  Variable_D = D_selected,
  Variable_E = E_selected
)

# --- Verification ---
cat("Number of cases per variable (binary):", nrow(XOR_genetic_data), "\n")
cat("Number of cases per variable (genotype):", nrow(XOR_genetic_data_genotype), "\n")

cat("\nClass distribution:\n")
print(table(XOR_genetic_data$class))

cat("\nChi-square test for Variable_C (binary) vs class:\n")
print(chisq.test(table(XOR_genetic_data$class, XOR_genetic_data$Variable_C)))

cat("\nChi-square test for Variable_D (binary) vs class:\n")
print(chisq.test(table(XOR_genetic_data$class, XOR_genetic_data$Variable_D)))

cat("\nCheck for missing values:\n")
print(sum(is.na(XOR_genetic_data)))
print(sum(is.na(XOR_genetic_data_genotype)))

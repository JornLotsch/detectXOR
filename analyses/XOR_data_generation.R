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
# 5. XOR linear shapes 
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


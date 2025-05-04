# XOR Artificial Dataset Generation with Sliding Means and Varying Standard Deviations
set.seed(123)  # for reproducibility

# Total sample size and classes
n <- 400
n1 <- 200
class <- factor(c(rep(1, n1), rep(2, n1)))  

# Sliding means for XOR pattern generation
m1_values <- seq(3, 10, length.out = 10)  
m2_values <- seq(10, 3, length.out = 10)  

# Set of standard deviations for XOR overlap variation
std_devs <- c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)

# Initialize data frames
XOR_sliding_means_data <- data.frame(class = class)  
XOR_sliding_means_params <- data.frame()  

var_index <- 1 

for (sd_orig in std_devs) {  
  for (i in seq_along(m1_values)) {  
    # Current means and jittered parameters
    m1 <- jitter(m1_values[i], factor=1)
    m2 <- jitter(m2_values[i], factor=1)
    sd_current <- jitter(sd_orig, factor=1)
    
    nn1 <- n1 %/% 2
    
    # Generate XOR structure clearly defined across classes
    # Class 1, first half: (m1,m1); second half: (m2,m2)
    A_class1 <- c(rnorm(nn1, mean=m1, sd=sd_current),
                  rnorm(nn1, mean=m2, sd=sd_current))
    B_class1 <- c(rnorm(nn1, mean=m1, sd=sd_current),
                  rnorm(nn1, mean=m2, sd=sd_current))
    
    # Class 2 XOR reverse arrangement, forming classical XOR:
    # first half: (m2,m1); second half: (m1,m2)
    A_class2 <- c(rnorm(nn1, mean=m2, sd=sd_current),
                  rnorm(nn1, mean=m1, sd=sd_current))
    B_class2 <- c(rnorm(nn1, mean=m1, sd=sd_current),
                  rnorm(nn1, mean=m2, sd=sd_current))
    
    # Combine both classes
    A_combined <- c(A_class1, A_class2)
    B_combined <- c(B_class1, B_class2)
    
    # Assign columns with unique names
    XOR_sliding_means_data[[paste0("A", var_index)]] <- A_combined
    XOR_sliding_means_data[[paste0("B", var_index)]] <- B_combined
    
    # Record metadata entries
    XOR_sliding_means_params <- rbind(XOR_sliding_means_params, data.frame(
      var1 = paste0("A", var_index),
      var2 = paste0("B", var_index),
      m1 = m1,
      m2 = m2,
      S = sd_current
    ))
    
    # Advance variable index
    var_index <- var_index + 1
  }
}

# Inspect the results quickly:
dim(XOR_sliding_means_data)
head(XOR_sliding_means_data)
dim(XOR_sliding_means_params)
head(XOR_sliding_means_params)

# Quick visual XOR structure verification (example pair A1,B1):
plot(XOR_sliding_means_data$A1, XOR_sliding_means_data$B1,
     col = XOR_sliding_means_data$class,
     pch = 16, xlab="A1", ylab="B1",
     main="Example XOR Structure (first pair)")

# detectXOR: Advanced XOR Pattern Detection in R
Detect non-linear relationships in classification data with statistical confidence

## Features
🔍 XOR Pattern Detection - Identify non-linear XOR relationships between variables using statistical tests (χ², Wilcoxon)
📊 Interactive Visualizations - Generate publication-ready spaghetti plots and decision boundary maps
⚡ Parallel Processing - Accelerate analysis using multi-core support (future/pbmcapply)
📈 Dependency Analysis - Compute class-wise τ coefficients to quantify variable relationships
🔬 Scientific Validation - Tested against synthetic datasets with known patterns


## Installation
You can install voronoi_projection_plot directly from GitHub:
```r
devtools::install_github("JornLotsch/detectXOR")
```

## Quick Start
### Detect XOR patterns

```r
data(XOR_data)
results <- detect_xor(XOR_data, class_col = "class")
```

### Visualize results
```r
generate_spaghetti_plot_from_results(results, XOR_data)
generate_xy_plot_from_results(results, XOR_data)
```

## Why detectXOR?
Traditional feature selection often misses complex non-linear relationships. Our package specifically targets XOR patterns - where class differences only emerge through variable interactions - using rigorous statistical methods and intuitive visualizations.

## Ideal for:

Machine learning feature engineering

Exploratory data analysis

Educational demonstrations of non-linear patterns

Bioinformatics and medical research

## Technical Highlights
✅ Comprehensive Testing - Unit tests covering 95% of critical functions
📦 CRAN-Ready Structure - Proper namespace handling and dependency management
🌐 Cross-Platform - Works on Windows, Linux, and macOS
📚 Rich Documentation - Detailed help files and usage examples

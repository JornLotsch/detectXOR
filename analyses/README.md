# Analyses for the detectXOR manuscript

R scripts that reproduce all data sets, tables and figures of the manuscript
"Generalized XOR patterns reveal accurate predictors without statistical
significance" (Lötsch, Kringel, Ultsch). The folder is not part of the R package
build.

## Requirements

- R >= 4.5 and detectXOR 0.2.0 (`install.packages("detectXOR")`)
- CRAN packages: ggplot2, dplyr, cowplot, circlize, svglite, ggforce, ggthemes,
  GGally, ggpubr, ggExtra, reshape2, randomForest, e1071, pROC, afex, mixOmics
  (Bioconductor), umap, Umatrix, dbt.DataIO, MASS

## Usage

Run the scripts from this directory in numerical order, e.g.

```sh
cd analyses
Rscript 01_generate_datasets.R
Rscript 02_detection_example_datasets.R
```

Script 01 must be run first; the others are independent of each other. The
number of cores used for the simulations is set by the environment variable
`XOR_CORES` (default: all cores minus one). All random processes are seeded, and
results do not depend on the number of cores.

| Script | Content | Manuscript |
| --- | --- | --- |
| `00_setup.R` | Common settings, package version check, helper functions | |
| `01_generate_datasets.R` | Example data sets (introductory, rotated/uniform, incomplete, genetic, sliding means) | Sections 2, 3.4 |
| `02_detection_example_datasets.R` | XOR screening of the example data sets, class imbalance, sliding means | Table 3, Figures 5-8 |
| `03_ml_benchmark.R` | Random forest, SVM and logistic regression (with and without interaction) in repeated stratified cross-validation | Table 2 |
| `04_simulation_study.R` | Family-wise false-positive rate under the null hypothesis; detection rates for XOR and non-XOR structures across sample sizes, distributions and class ratios | Figure 9 |
| `05_sensitivity_analysis.R` | Sensitivity to tile boundaries, Kendall's tau threshold, minimum purity and minimum tile share | Figure 10 |
| `06_interaction_method_comparison.R` | Comparison with logistic regression interaction test, interaction information and MDR | Figure 11 |
| `07_real_data.R` | Lipid mediators in multiple sclerosis and Pima Indians Diabetes data | Section 4.6 |
| `08_descriptive_statistics_figures.R` | Conventional tests and descriptive figures | Table 1, Figures 1-4 |
| `simulation_functions.R` | Data generators for scripts 04-06 | |
| `descriptive_plot_functions.R` | Plot functions for Figures 1-4 | |

The lipidomics data of script 07 are not included in the repository; set the
environment variable `XOR_LIPIDOMICS_DIR` to their location to include this part.
The trained emergent self-organizing maps for Figure 2 are provided in
`data/esom` (trained interactively with `Umatrix::iEsomTrain`).

Results are written to `results/` (CSV tables, session information) and figures
to `figures/` (SVG and EPS).

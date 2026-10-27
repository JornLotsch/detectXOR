# Common setup for the analysis scripts of the detectXOR manuscript.
# All scripts are run from the 'analyses' directory, e.g.
#   cd analyses; Rscript 02_detection_example_datasets.R

if (!file.exists("00_setup.R")) {
  stop("Please run the analysis scripts from the 'analyses' directory.")
}

# Package version used for all analyses of the manuscript
REQUIRED_DETECTXOR_VERSION <- "0.2.0"

suppressPackageStartupMessages(library(detectXOR))
if (packageVersion("detectXOR") < REQUIRED_DETECTXOR_VERSION) {
  stop("detectXOR >= ", REQUIRED_DETECTXOR_VERSION, " is required, found ", packageVersion("detectXOR"))
}

# Output directories
DATA_DIR <- "data"
RESULTS_DIR <- "results"
FIGURES_DIR <- "figures"
for (d in c(DATA_DIR, RESULTS_DIR, FIGURES_DIR)) dir.create(d, showWarnings = FALSE)

# Number of cores for simulations; override with the environment variable XOR_CORES
N_CORES <- as.integer(Sys.getenv("XOR_CORES", max(1, parallel::detectCores() - 1)))

# Number of permutations used in all detect_xor() calls
N_PERM <- 999

# Load the example data sets created by 01_generate_datasets.R
load_example_datasets <- function() {
  file <- file.path(DATA_DIR, "example_datasets.rds")
  if (!file.exists(file)) stop("Run 01_generate_datasets.R first.")
  readRDS(file)
}

# Wilson score confidence interval for a proportion
wilson_ci <- function(x, n, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- x / n
  centre <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  cbind(lower = pmax(0, centre - half), upper = pmin(1, centre + half))
}

# Write session information next to the results of a script
write_session_info <- function(script_name) {
  writeLines(capture.output(sessionInfo()),
             file.path(RESULTS_DIR, paste0("sessionInfo_", script_name, ".txt")))
}

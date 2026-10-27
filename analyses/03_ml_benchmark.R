# Classification performance of variable pairs in the introductory data sets
# (manuscript Section 2.1.2, Table 2).
#
# Design: repeated stratified 5-fold cross-validation (20 repetitions, i.e.,
# 100 test folds). In each fold, models are trained and tuned on the training
# folds only and evaluated on the held-out test fold. Class 2 is the positive
# class. Performance is summarized by the median and the 2.5% and 97.5%
# percentiles across the 100 test folds.
#
# Classifiers: random forest (mtry chosen by out-of-bag error), support vector
# machine with radial kernel (cost and gamma tuned by 5-fold cross-validation
# within the training folds), logistic regression with main effects only, and
# logistic regression including the interaction term.
# Output: results/ml_benchmark_folds.csv, results/ml_benchmark_table2.csv

source("00_setup.R")
suppressPackageStartupMessages({
  library(randomForest)
  library(e1071)
})

datasets <- load_example_datasets()
FEATURE_SETS <- list("A, B" = c("Variable_A", "Variable_B"), "C, D" = c("Variable_C", "Variable_D"))
N_REPETITIONS <- 20
N_FOLDS <- 5
SEED <- 42

# Stratified fold assignment
stratified_folds <- function(cls, k) {
  folds <- integer(length(cls))
  for (cl in unique(cls)) {
    idx <- which(cls == cl)
    folds[idx] <- sample(rep_len(seq_len(k), length(idx)))
  }
  folds
}

metrics <- function(truth, prob_class2) {
  pred <- ifelse(prob_class2 > 0.5, 2, 1)
  sensitivity <- mean(pred[truth == 2] == 2)
  specificity <- mean(pred[truth == 1] == 1)
  # Fixed direction: higher predicted probability of class 2 indicates class 2
  auc <- as.numeric(pROC::auc(pROC::roc(truth, prob_class2, levels = c(1, 2), direction = "<", quiet = TRUE)))
  c(balanced_accuracy = (sensitivity + specificity) / 2, sensitivity = sensitivity,
    specificity = specificity, roc_auc = auc)
}

fit_and_predict <- function(train, test) {
  train$class <- factor(train$class, levels = c(1, 2))
  features <- setdiff(names(train), "class")
  formula_main <- stats::as.formula(paste("class ~", paste(features, collapse = " + ")))
  formula_interaction <- stats::as.formula(paste("class ~", paste(features, collapse = " * ")))

  # Random forest: mtry by out-of-bag error
  rf_fits <- lapply(seq_along(features), function(m) randomForest(formula_main, data = train, ntree = 500, mtry = m))
  rf <- rf_fits[[which.min(vapply(rf_fits, function(f) f$err.rate[f$ntree, "OOB"], numeric(1)))]]

  # SVM: grid search by 5-fold cross-validation within the training data
  tuned <- e1071::tune.svm(formula_main, data = train, kernel = "radial",
                           cost = 2^(-2:6), gamma = 2^(-4:2),
                           tunecontrol = e1071::tune.control(sampling = "cross", cross = 5))
  svm_fit <- e1071::svm(formula_main, data = train, kernel = "radial", probability = TRUE,
                        cost = tuned$best.parameters$cost, gamma = tuned$best.parameters$gamma)

  logit_main <- suppressWarnings(stats::glm(formula_main, data = train, family = stats::binomial))
  logit_interaction <- suppressWarnings(stats::glm(formula_interaction, data = train, family = stats::binomial))

  svm_prob <- attr(stats::predict(svm_fit, test, probability = TRUE), "probabilities")[, "2"]
  list(
    "Random forest" = stats::predict(rf, test, type = "prob")[, "2"],
    "SVM" = svm_prob,
    "Logistic regression, main effects" = stats::predict(logit_main, test, type = "response"),
    "Logistic regression with interaction" = stats::predict(logit_interaction, test, type = "response")
  )
}

jobs <- expand.grid(dataset = c("XOR_data", "XOR_data_nonsig"), feature_set = names(FEATURE_SETS),
                    repetition = seq_len(N_REPETITIONS), stringsAsFactors = FALSE)
out <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  job <- jobs[j, ]
  data <- datasets[[job$dataset]][, c("class", FEATURE_SETS[[job$feature_set]])]
  set.seed(SEED + job$repetition)
  folds <- stratified_folds(data$class, N_FOLDS)
  do.call(rbind, lapply(seq_len(N_FOLDS), function(f) {
    probs <- fit_and_predict(data[folds != f, ], data[folds == f, ])
    do.call(rbind, lapply(names(probs), function(m) {
      data.frame(dataset = job$dataset, feature_set = job$feature_set, repetition = job$repetition,
                 fold = f, classifier = m, t(metrics(data$class[folds == f], probs[[m]])))
    }))
  }))
}, mc.cores = min(N_CORES, nrow(jobs)))
failed <- vapply(out, inherits, logical(1), what = "try-error")
if (any(failed)) stop("Benchmark failed: ", as.character(out[[which(failed)[1]]]))
folds <- do.call(rbind, out)
write.csv(folds, file.path(RESULTS_DIR, "ml_benchmark_folds.csv"), row.names = FALSE)

# Table 2: median (2.5%-97.5% percentiles) across the 100 test folds
metric_names <- c("balanced_accuracy", "sensitivity", "specificity", "roc_auc")
table2 <- do.call(rbind, lapply(split(folds, list(folds$dataset, folds$feature_set, folds$classifier), drop = TRUE), function(d) {
  do.call(rbind, lapply(metric_names, function(m) {
    q <- stats::quantile(d[[m]], c(0.5, 0.025, 0.975))
    data.frame(dataset = d$dataset[1], feature_set = d$feature_set[1], classifier = d$classifier[1], metric = m,
               median = q[1], lower = q[2], upper = q[3], excludes_0.5 = q[2] > 0.5 | q[3] < 0.5,
               formatted = sprintf("%.3f (%.3f--%.3f)", q[1], q[2], q[3]))
  }))
}))
table2 <- table2[order(table2$feature_set, table2$classifier, match(table2$metric, metric_names), table2$dataset), ]
write.csv(table2, file.path(RESULTS_DIR, "ml_benchmark_table2.csv"), row.names = FALSE)
print(stats::reshape(table2[, c("dataset", "feature_set", "classifier", "metric", "formatted")],
                     idvar = c("feature_set", "classifier", "metric"), timevar = "dataset", direction = "wide"),
      row.names = FALSE)

write_session_info("03_ml_benchmark")

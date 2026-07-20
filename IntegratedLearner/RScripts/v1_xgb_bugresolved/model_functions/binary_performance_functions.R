#This file is specifically for evaluating performance metrics for continuous data

# Train metrics collections
kappa_train_results <- list(
  kappa_fold_metabolomics, kappa_fold_mss, kappa_fold_concat,
  kappa_averaged_stacked, kappa_weighted_nnls, kappa_sparse_nnls, kappa_pls
)

auroc_train_results <- list(
  auroc_fold_metabolomics, auroc_fold_mss, auroc_fold_concat,
  auroc_averaged_stacked, auroc_weighted_nnls, auroc_sparse_nnls, auroc_pls
)

accuracy_train_results <- list(
  accuracy_fold_metabolomics, accuracy_fold_mss, accuracy_fold_concat,
  accuracy_averaged_stacked, accuracy_weighted_nnls, accuracy_sparse_nnls, accuracy_pls
)

auprc_train_results <- list(
  auprc_fold_metabolomics, auprc_fold_mss, auprc_fold_concat,
  auprc_averaged_stacked, auprc_weighted_nnls, auprc_sparse_nnls, auprc_pls
)

sensitivity_train_results <- list(
  sensitivity_fold_metabolomics, sensitivity_fold_mss, sensitivity_fold_concat,
  sensitivity_averaged_stacked, sensitivity_weighted_nnls, sensitivity_sparse_nnls, sensitivity_pls
)

specificity_train_results <- list(
  specificity_fold_metabolomics, specificity_fold_mss, specificity_fold_concat,
  specificity_averaged_stacked, specificity_weighted_nnls, specificity_sparse_nnls, specificity_pls
)

precision_train_results <- list(
  precision_fold_metabolomics, precision_fold_mss, precision_fold_concat,
  precision_averaged_stacked, precision_weighted_nnls, precision_sparse_nnls, precision_pls
)

f1_train_results <- list(
  f1_fold_metabolomics, f1_fold_mss, f1_fold_concat,
  f1_averaged_stacked, f1_weighted_nnls, f1_sparse_nnls, f1_pls
)

balanced_accuracy_train_results <- list(
  balanced_accuracy_fold_metabolomics, balanced_accuracy_fold_mss, balanced_accuracy_fold_concat,
  balanced_accuracy_averaged_stacked, balanced_accuracy_weighted_nnls,
  balanced_accuracy_sparse_nnls, balanced_accuracy_pls
)

# Test metrics collections
kappa_test_results <- list(
  test_kappa_fold_metabolomics, test_kappa_fold_mss, test_kappa_fold_concat,
  test_kappa_averaged_stacked, test_kappa_weighted_nnls,
  test_kappa_sparse_nnls, test_kappa_pls
)

auroc_test_results <- list(
  test_auroc_fold_metabolomics, test_auroc_fold_mss, test_auroc_fold_concat,
  test_auroc_averaged_stacked, test_auroc_weighted_nnls,
  test_auroc_sparse_nnls, test_auroc_pls
)

accuracy_test_results <- list(
  test_accuracy_fold_metabolomics, test_accuracy_fold_mss, test_accuracy_fold_concat,
  test_accuracy_averaged_stacked, test_accuracy_weighted_nnls,
  test_accuracy_sparse_nnls, test_accuracy_pls
)

auprc_test_results <- list(
  test_auprc_fold_metabolomics, test_auprc_fold_mss, test_auprc_fold_concat,
  test_auprc_averaged_stacked, test_auprc_weighted_nnls,
  test_auprc_sparse_nnls, test_auprc_pls
)

sensitivity_test_results <- list(
  test_sensitivity_fold_metabolomics, test_sensitivity_fold_mss, test_sensitivity_fold_concat,
  test_sensitivity_averaged_stacked, test_sensitivity_weighted_nnls,
  test_sensitivity_sparse_nnls, test_sensitivity_pls
)

specificity_test_results <- list(
  test_specificity_fold_metabolomics, test_specificity_fold_mss, test_specificity_fold_concat,
  test_specificity_averaged_stacked, test_specificity_weighted_nnls,
  test_specificity_sparse_nnls, test_specificity_pls
)

precision_test_results <- list(
  test_precision_fold_metabolomics, test_precision_fold_mss, test_precision_fold_concat,
  test_precision_averaged_stacked, test_precision_weighted_nnls,
  test_precision_sparse_nnls, test_precision_pls
)

f1_test_results <- list(
  test_f1_fold_metabolomics, test_f1_fold_mss, test_f1_fold_concat,
  test_f1_averaged_stacked, test_f1_weighted_nnls,
  test_f1_sparse_nnls, test_f1_pls
)

balanced_accuracy_test_results <- list(
  test_balanced_accuracy_fold_metabolomics, test_balanced_accuracy_fold_mss,
  test_balanced_accuracy_fold_concat, test_balanced_accuracy_averaged_stacked,
  test_balanced_accuracy_weighted_nnls, test_balanced_accuracy_sparse_nnls,
  test_balanced_accuracy_pls
)


# Define the row labels
row_labels <- c(
  "Metabolomics", "MSS", "Concatenated", "Averaged Stacked",
  "Weighted NNLS", "Lasso Stacked", "PLS"
)

# Summarize training metrics
kappa_train_df             <- stats_results(kappa_train_results, row_labels, "Kappa Train")
auroc_train_df             <- stats_results(auroc_train_results, row_labels, "AUROC Train")
accuracy_train_df          <- stats_results(accuracy_train_results, row_labels, "Accuracy Train")
auprc_train_df             <- stats_results(auprc_train_results, row_labels, "AUPRC Train")
sensitivity_train_df       <- stats_results(sensitivity_train_results, row_labels, "Sensitivity Train")
specificity_train_df       <- stats_results(specificity_train_results, row_labels, "Specificity Train")
precision_train_df         <- stats_results(precision_train_results, row_labels, "Precision Train")
f1_train_df                <- stats_results(f1_train_results, row_labels, "F1 Train")
balanced_accuracy_train_df <- stats_results(balanced_accuracy_train_results, row_labels, "Balanced Accuracy Train")

# Summarize test metrics
kappa_test_df             <- stats_results(kappa_test_results, row_labels, "Kappa Test")
auroc_test_df             <- stats_results(auroc_test_results, row_labels, "AUROC Test")
accuracy_test_df          <- stats_results(accuracy_test_results, row_labels, "Accuracy Test")
auprc_test_df             <- stats_results(auprc_test_results, row_labels, "AUPRC Test")
sensitivity_test_df       <- stats_results(sensitivity_test_results, row_labels, "Sensitivity Test")
specificity_test_df       <- stats_results(specificity_test_results, row_labels, "Specificity Test")
precision_test_df         <- stats_results(precision_test_results, row_labels, "Precision Test")
f1_test_df                <- stats_results(f1_test_results, row_labels, "F1 Test")
balanced_accuracy_test_df <- stats_results(balanced_accuracy_test_results, row_labels, "Balanced Accuracy Test")

all_results_df <- reduce(
  list(
    kappa_train_df, auroc_train_df, accuracy_train_df,
    auprc_train_df, sensitivity_train_df, specificity_train_df,
    precision_train_df, f1_train_df, balanced_accuracy_train_df,
    kappa_test_df, auroc_test_df, accuracy_test_df,
    auprc_test_df, sensitivity_test_df, specificity_test_df,
    precision_test_df, f1_test_df, balanced_accuracy_test_df
  ),
  function(x, y) full_join(x, y, by = "Category")
)

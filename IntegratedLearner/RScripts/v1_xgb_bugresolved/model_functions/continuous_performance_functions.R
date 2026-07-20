# Combine training lists for each metric
rmse_train_results      <- list(rmse_fold_metabolomics, rmse_fold_mss, rmse_fold_concat,
                                rmse_averaged_stacked, rmse_weighted_nnls, rmse_sparse_nnls, rmse_pls)
r2_train_results        <- list(r2_fold_metabolomics, r2_fold_mss, r2_fold_concat,
                                r2_averaged_stacked, r2_weighted_nnls, r2_sparse_nnls, r2_pls)
mae_train_results       <- list(mae_fold_metabolomics, mae_fold_mss, mae_fold_concat,
                                mae_averaged_stacked, mae_weighted_nnls, mae_sparse_nnls, mae_pls)
mape_train_results      <- list(mape_fold_metabolomics, mape_fold_mss, mape_fold_concat,
                                mape_averaged_stacked, mape_weighted_nnls, mape_sparse_nnls, mape_pls)
median_abs_error_train_results <- list(median_abs_error_fold_metabolomics,
                                       median_abs_error_fold_mss,
                                       median_abs_error_fold_concat,
                                       median_abs_error_averaged_stacked,
                                       median_abs_error_weighted_nnls,
                                       median_abs_error_sparse_nnls, 
                                       median_abs_error_pls)

# Combine test lists for each metric
rmse_test_results  <- list(test_rmse_fold_metabolomics, test_rmse_fold_mss, test_rmse_fold_concat,
                           test_rmse_averaged_stacked, test_rmse_weighted_nnls,
                           test_rmse_sparse_nnls, test_rmse_pls)
r2_test_results    <- list(test_r2_fold_metabolomics, test_r2_fold_mss, test_r2_fold_concat,
                           test_r2_averaged_stacked, test_r2_weighted_nnls,
                           test_r2_sparse_nnls, test_r2_pls)
mae_test_results   <- list(test_mae_fold_metabolomics, test_mae_fold_mss, test_mae_fold_concat,
                           test_mae_averaged_stacked, test_mae_weighted_nnls,
                           test_mae_sparse_nnls, test_mae_pls)
mape_test_results  <- list(test_mape_fold_metabolomics, test_mape_fold_mss, test_mape_fold_concat,
                           test_mape_averaged_stacked, test_mape_weighted_nnls,
                           test_mape_sparse_nnls, test_mape_pls)
median_abs_error_test_results <- list(test_median_abs_error_fold_metabolomics,
                                      test_median_abs_error_fold_mss,
                                      test_median_abs_error_fold_concat,
                                      test_median_abs_error_averaged_stacked,
                                      test_median_abs_error_weighted_nnls,
                                      test_median_abs_error_sparse_nnls,
                                      test_median_abs_error_pls)

# Define the row labels
row_labels <- c("Metabolomics", "MSS", "Concatenated", "Averaged Stacked",
                "Weighted NNLS", "Lasso Stacked", "PLS")

# Summarize training metrics
rmse_train_df            <- stats_results(rmse_train_results, row_labels, "RMSE Train")
r2_train_df              <- stats_results(r2_train_results, row_labels, "R2 Train")
mae_train_df             <- stats_results(mae_train_results, row_labels, "MAE Train")
mape_train_df            <- stats_results(mape_train_results, row_labels, "MAPE Train")
median_abs_error_train_df<- stats_results(median_abs_error_train_results, row_labels, "Median Abs Error Train")

# Summarize test metrics
rmse_test_df             <- stats_results(rmse_test_results, row_labels, "RMSE Test")
r2_test_df               <- stats_results(r2_test_results, row_labels, "R2 Test")
mae_test_df              <- stats_results(mae_test_results, row_labels, "MAE Test")
mape_test_df             <- stats_results(mape_test_results, row_labels, "MAPE Test")
median_abs_error_test_df <- stats_results(median_abs_error_test_results, row_labels, "Median Abs Error Test")

all_results_df <- reduce(
  list(
    rmse_train_df, r2_train_df, mae_train_df,
    mape_train_df, median_abs_error_train_df,
    rmse_test_df, r2_test_df, mae_test_df,
    mape_test_df, median_abs_error_test_df
  ),
  function(x, y) {
    full_join(x, y, by = "Category")
  }
)

# MicrobiomeMetabolomeIntegration2026


## DataPreprocessing

This subdirectory contains R scripts used to prepare the Erawijantari, Franzosa, Wang, and Yachida datasets for downstream modeling. Each script combines genus- and species-level taxonomic profiles, metabolomics data, and selected metadata by sample. Metabolite and taxonomic features are labeled with `m__` and `t__` prefixes, respectively. Study-specific datasets and outcome comparisons are then exported as CSV files.

Before running a script, update `setwd()` to the corresponding study-data directory containing `genera.tsv`, `species.tsv`, `mtb.tsv`, `mtb.map.tsv`, and `metadata.tsv`.

The processed datasets can be found at: https://github.com/borenstein-lab/microbiome-metabolome-curated-data/wiki/Data-overview#datasets-included

## `IntegratedLearner/RScripts/v2_Rdata_implementation`

This subdirectory contains the primary modeling workflow and its supporting learner and performance functions. The pipeline expects metabolite features to begin with `m__` and taxonomic features with `t__`.

### Main pipeline

- **`automated_pipeline_v2.R`** runs the complete integrative machine-learning analysis. It:
  - reads command-line or configuration-file parameters;
  - performs metabolomics and taxonomic transformation and feature reduction;
  - creates 20 training/testing splits and repeated cross-validation folds;
  - trains metabolomics-only, taxonomy-only, and concatenated models;
  - evaluates unscaled, scaled, and block-weighted concatenation;
  - integrates metabolomics and taxonomic predictions using averaging, non-negative least squares, sparse regression, and partial least squares; and
  - exports performance summaries, feature importance, selected-feature stability, split diagnostics, run documentation, and `.RData` files.

Preprocessing and feature selection are performed within each training fold to prevent information from validation or test samples from influencing model development. Completed-run markers allow interrupted analyses to resume without repeating successful runs.

### Supporting scripts

- **`model_functions/enet_function.R`** fits and predicts binary or continuous outcomes using Elastic Net.
- **`model_functions/rf_function.R`** fits and predicts outcomes using Random Forest.
- **`model_functions/xgboost_function.R`** fits XGBoost models with CPU or GPU support and validates predictor names to prevent outcome leakage.
- **`model_functions/binary_functions.R`** calculates accuracy, AUROC, AUPRC, kappa, sensitivity, specificity, precision, F1 score, and balanced accuracy.
- **`model_functions/continuous_functions.R`** calculates \(R^2\), RMSE, MAE, MAPE, and median absolute error.
- **`model_functions/binary_performance_functions.R`** and **`continuous_performance_functions.R`** organize stored performance values into summary tables.

## `ModelPerformanceExtraction`

This subdirectory contains scripts that aggregate and visualize model performance across seeds and analytical configurations. Before running either script, update `input_base` and `sub_dir` to identify the directory containing the pipeline performance CSV files.

- **`model_performance_binary.R`** summarizes binary-outcome performance using accuracy, kappa, AUROC, AUPRC, sensitivity, specificity, precision, F1 score, and balanced accuracy.
- **`model_performance_continuous.R`** summarizes continuous-outcome performance using \(R^2\), RMSE, MAE, MAPE, and median absolute error.

Both scripts parse the learner, transformation, feature-reduction, and seed information from the pipeline filenames. They also check for missing runs, reshape results into a long format, calculate across-seed means and 95% confidence intervals, and compare model and integration categories.

Results are written to a study-specific `results` directory and include:

- CSV summary tables and missing-run reports;
- ranked dot-and-whisker performance plots;
- heatmaps comparing pipeline configurations and integration categories;
- summaries of how frequently each category achieved the best performance; and
- an HTML report combining the principal tables and figures.

## `ModelPerformanceSummaries`

This subdirectory contains scripts that summarize model performance across datasets and generate the analyses used for Figures 2 and 3. These scripts use the tables produced by the `ModelPerformanceExtraction` workflow. File paths and dataset lists must be updated before running the scripts.

### Continuous outcomes: Figure 2

- **`Figure2_perform_CIs.R`** reads each dataset’s `all_metrics_long.csv` file and calculates the mean and 95% confidence interval across seeds for RMSE and MAE. It generates separate training- and test-performance plots for each dataset, faceted by model or integration category.
- **`Figure2_summary_trends_v2.R`** compares RMSE and MAE results across continuous-outcome datasets. It defines a competitive set containing configurations within 2% of the lowest error and summarizes how frequently each learner, transformation, dimensionality, and integration strategy appears in this set.

### Binary outcomes: Figure 3

- **`Figure3_perform_CIs.R`** calculates across-seed means and 95% confidence intervals for binary classification metrics, including accuracy, AUROC, AUPRC, sensitivity, specificity, precision, F1 score, and balanced accuracy. Separate training- and test-performance plots are generated for each dataset.
- **`Figure3_summary_trends_v2.R`** compares binary-outcome configurations across datasets using AUROC and AUPRC. It defines the competitive set as configurations within 2% of the highest performance and summarizes the analytical choices represented among these configurations.

The summary-trend scripts provide both raw counts, with exact ties sharing credit, and dataset-equal counts, where each dataset contributes equal total weight. Outputs include PDF summary panels, competitive-set tile plots, CSV frequency tables, and an Excel workbook containing training- and testing-performance matrices.

## `FeatureSelectionStabilityCode`

This subdirectory contains the workflow used to evaluate whether different learners, preprocessing choices, and integration strategies prioritize the same biological features. The analysis ranks features using their importance and recurrence across resampled datasets, selects the top 20 features from each analytical condition, and calculates the percentage shared between pairs of conditions.

### Feature-importance preprocessing

- **`PreprocessingFeatures_1_condensed.R`** combines validation-fold feature importance across cross-validation folds and 20 train/test splits. It calculates each feature’s mean importance, percentile rank, frequency across seeds, and frequency across cross-validation fits. Elastic Net coefficients are converted to absolute values before summarization.
- **`PreprocessingFeatures_2_condensed.R`** calculates stability scores by multiplying each feature’s importance percentile rank by its frequency across seeds or cross-validation fits. These scores are used to rank features for the overlap analyses.

### Pairwise feature-overlap generation

- **`Concat_heatmap_generation_3_condensed_diagonal100.R`** selects the top 20 features for each learner and preprocessing configuration and calculates pairwise overlap within metabolomics, taxonomy, and the three concatenation strategies. Each view contains 24 analytical conditions.
- **`SingleOmics_heatmap_generation_3_condensed_diagonal100_revised.R`** performs the corresponding analysis for metabolomics-only and taxonomy-only models. Preprocessing choices from the other omics type are ignored, producing 12 relevant conditions per view.

Feature overlap is reported as:

\[
\text{Top-20 overlap (\%)} =
\frac{\text{Number of shared top-20 features}}{20}\times 100.
\]

Diagonal cells represent comparisons of a condition with itself and are displayed as 100%.

### Cross-dataset summaries

- **`Concat_binary_summaries_4_diagonal100.R`** aggregates the 24-by-24 overlap matrices across binary-outcome datasets.
- **`Concat_continuous_summaries_5_diagonal100.R`** aggregates the 24-by-24 overlap matrices across continuous-outcome datasets.
- **`SingleOmics_binary_summaries_4_diagonal100_revised.R`** aggregates the 12-by-12 metabolomics and taxonomy overlap matrices across binary-outcome datasets.
- **`SingleOmics_continuous_summaries_5_diagonal100_revised.R`** aggregates the corresponding single-omics matrices across continuous-outcome datasets.

Outputs include top-feature tables, pairwise-overlap CSV files, dataset-level and cross-dataset summary statistics, diagnostic tables, and lower-triangle PDF heatmaps. Update the root directory and dataset lists in each script before running the workflow.




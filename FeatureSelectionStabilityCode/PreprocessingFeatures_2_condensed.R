library(tidyverse)

# Add the stable-score columns required by heatmap_generation_3.R to the
# per-combination seed summaries created by PreprocessingFeatures_1_condensed.R.

# ---- User settings ----

root_dir <- "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/FeatureImportanceLogs/Validation_Feature_Importance"

seed_summary_dir <- file.path(
  root_dir,
  "_feature_importance_seed_summaries"
)

stable_score_dir <- file.path(
  root_dir,
  "_feature_importance_seed_summaries_with_stable_scores"
)

# ---- Helper function ----

add_heatmap_scores <- function(dat, path) {
  required_cols <- c(
    "CombinationID",
    "View",
    "Feature",
    "mean_importance_percentile_rank",
    "seed_presence_frequency"
  )

  missing_cols <- setdiff(required_cols, names(dat))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required column(s) in ", path, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # CV-based stability is retained because heatmap_generation_3.R uses
  # stable_score_cv as a tie-breaker and in its recurrence summaries.
  if (!"overall_cv_fit_presence_frequency" %in% names(dat)) {
    dat <- dat %>%
      mutate(overall_cv_fit_presence_frequency = NA_real_)
  }

  dat %>%
    mutate(
      Feature = str_trim(as.character(Feature)),
      stable_score_seed =
        mean_importance_percentile_rank * seed_presence_frequency,
      stable_score_cv =
        mean_importance_percentile_rank *
        overall_cv_fit_presence_frequency
    ) %>%
    select(
      any_of(c("DatasetID", "CombinationID", "View", "ModelType")),
      Feature,
      mean_importance_percentile_rank,
      seed_presence_frequency,
      stable_score_seed,
      stable_score_cv
    )
}

# ---- Create the files read by heatmap_generation_3.R ----

dir.create(stable_score_dir, recursive = TRUE, showWarnings = FALSE)

seed_summary_files <- list.files(
  path = seed_summary_dir,
  pattern = "__seed_importance_summary\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(seed_summary_files) == 0L) {
  stop("No per-combination seed-summary files found under:\n", seed_summary_dir)
}

walk(seed_summary_files, function(path) {
  dat <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  heatmap_dat <- add_heatmap_scores(dat, path)

  dataset_id <- basename(dirname(path))
  dataset_output_dir <- file.path(stable_score_dir, dataset_id)

  dir.create(dataset_output_dir, recursive = TRUE, showWarnings = FALSE)

  output_file <- file.path(
    dataset_output_dir,
    basename(path) %>%
      str_replace(
        "__seed_importance_summary\\.csv$",
        "__seed_importance_summary_with_stable_scores.csv"
      )
  )

  write_csv(heatmap_dat, output_file, na = "")
})

message("Done.")
message(
  length(seed_summary_files),
  " stable-score files written under: ",
  stable_score_dir
)

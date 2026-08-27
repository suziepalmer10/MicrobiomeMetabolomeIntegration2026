library(tidyverse)

# Convert validation-fold feature importance directly into the seed-summary
# files read by PreprocessingFeatures_2.R.

# ---- User settings ----

root_dir <- "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/FeatureImportanceLogs/Validation_Feature_Importance"

seed_summary_dir <- file.path(
  root_dir,
  "_feature_importance_seed_summaries"
)

expected_n_runs <- 20L
expected_cv_fits_per_run <- 15L
strict_run_check <- TRUE

run_file_re <- "_run([0-9]+)_validation_fold_feature_importance\\.csv$"

feature_id_cols <- c(
  "View",
  "ModelType",
  "Split",
  "Feature",
  "ImportanceType"
)

# ---- Helper functions ----

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

percent_rank_safe <- function(x) {
  n_observed <- sum(!is.na(x))

  if (n_observed == 0L) {
    rep(NA_real_, length(x))
  } else if (n_observed == 1L) {
    if_else(is.na(x), NA_real_, 1)
  } else {
    percent_rank(x)
  }
}

collapse_one_run <- function(path, run_num, combination_id) {
  dat <- read_csv(path, show_col_types = FALSE, progress = FALSE)

  required_cols <- c(
    "Seed",
    "View",
    "Feature",
    "Importance"
  )

  missing_cols <- setdiff(required_cols, names(dat))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required column(s) in ", path, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  seed_values <- unique(na.omit(dat$Seed))

  if (length(seed_values) != 1L) {
    stop("Expected exactly one Seed value in: ", path)
  }

  importance_text <- as.character(dat$Importance)
  importance_numeric <- suppressWarnings(as.numeric(importance_text))

  invalid_importance <-
    !is.na(dat$Importance) &
    str_trim(importance_text) != "" &
    is.na(importance_numeric)

  if (any(invalid_importance)) {
    stop("Importance contains non-numeric values in: ", path)
  }

  ids <- intersect(feature_id_cols, names(dat))

  dat <- dat %>%
    mutate(
      Feature = str_trim(as.character(Feature)),
      .importance = importance_numeric
    )

  if (str_detect(str_to_lower(combination_id), "model-enet|elastic")) {
    dat <- dat %>%
      mutate(.importance = abs(.importance))
  }

  has_cv_ids <- all(c("RepeatID", "Fold") %in% names(dat))

  if (has_cv_ids) {
    dat <- dat %>%
      unite(
        ".cv_fit_id",
        RepeatID,
        Fold,
        sep = "__",
        remove = FALSE
      )

    n_cv_fits_possible <- n_distinct(dat$.cv_fit_id)

    if (n_cv_fits_possible != expected_cv_fits_per_run) {
      warning(
        "Expected ", expected_cv_fits_per_run,
        " CV fits but found ", n_cv_fits_possible,
        " in: ", path
      )
    }
  } else {
    dat <- dat %>%
      mutate(.cv_fit_id = NA_character_)

    n_cv_fits_possible <- NA_integer_

    warning(
      "RepeatID and Fold were not both found; CV-fit frequency will be NA in: ",
      path
    )
  }

  dat %>%
    group_by(across(all_of(ids))) %>%
    summarise(
      # Conditional mean among CV fits where the feature was reported.
      .run_mean_importance = mean_or_na(.importance),
      .n_cv_fits_present = if (has_cv_ids) {
        n_distinct(.cv_fit_id[!is.na(.importance)])
      } else {
        NA_integer_
      },
      .groups = "drop"
    ) %>%
    mutate(
      .run_num = as.integer(run_num),
      .seed = as.character(seed_values[[1]]),
      .n_cv_fits_possible = as.integer(n_cv_fits_possible)
    )
}

summarize_one_combination <- function(file_tbl) {
  file_tbl <- arrange(file_tbl, run_num)

  dataset_id <- unique(file_tbl$dataset_id)
  combination_id <- unique(file_tbl$combination_id)

  duplicate_runs <- file_tbl %>%
    count(run_num) %>%
    filter(n > 1L)

  if (nrow(duplicate_runs) > 0L) {
    stop(
      "Duplicate run files for dataset ", dataset_id,
      ", combination ", combination_id, ": ",
      paste(duplicate_runs$run_num, collapse = ", ")
    )
  }

  observed_runs <- sort(file_tbl$run_num)
  expected_runs <- seq_len(expected_n_runs)

  if (strict_run_check && !identical(observed_runs, expected_runs)) {
    stop(
      "Run-number check failed for ", combination_id,
      ". Missing: ",
      paste(setdiff(expected_runs, observed_runs), collapse = ", "),
      "; unexpected: ",
      paste(setdiff(observed_runs, expected_runs), collapse = ", ")
    )
  }

  run_dat <- map_dfr(seq_len(nrow(file_tbl)), function(i) {
    collapse_one_run(
      path = file_tbl$path[[i]],
      run_num = file_tbl$run_num[[i]],
      combination_id = combination_id
    )
  })

  run_metadata <- run_dat %>%
    distinct(.run_num, .seed, .n_cv_fits_possible)

  if (nrow(run_metadata) != nrow(file_tbl)) {
    stop("Inconsistent run metadata for: ", combination_id)
  }

  if (n_distinct(run_metadata$.seed) != nrow(run_metadata)) {
    warning(
      "Seed values are reused for ", combination_id,
      "; run files will remain the seed-frequency denominator."
    )
  }

  n_seeds_total <- nrow(file_tbl)

  n_cv_fits_total <- if (
    all(!is.na(run_metadata$.n_cv_fits_possible))
  ) {
    sum(run_metadata$.n_cv_fits_possible)
  } else {
    NA_integer_
  }

  ids <- intersect(feature_id_cols, names(run_dat))

  seed_summary <- run_dat %>%
    group_by(across(all_of(ids))) %>%
    summarise(
      # Conditional mean among runs/seeds where the feature was reported.
      mean_importance = mean_or_na(.run_mean_importance),
      n_seeds_present = n_distinct(
        .run_num[!is.na(.run_mean_importance)]
      ),
      n_cv_fits_with_importance = if (
        all(is.na(.n_cv_fits_present))
      ) {
        NA_integer_
      } else {
        sum(.n_cv_fits_present, na.rm = TRUE)
      },
      .groups = "drop"
    ) %>%
    mutate(
      DatasetID = .env$dataset_id,
      CombinationID = .env$combination_id,
      n_seeds_total = .env$n_seeds_total,
      seed_presence_frequency = n_seeds_present / n_seeds_total,
      n_cv_fits_total = .env$n_cv_fits_total,
      overall_cv_fit_presence_frequency = if_else(
        !is.na(n_cv_fits_total) & n_cv_fits_total > 0L,
        as.numeric(n_cv_fits_with_importance) / n_cv_fits_total,
        NA_real_
      ),
      .before = 1
    )

  rank_groups <- intersect(
    c(
      "DatasetID",
      "CombinationID",
      "View",
      "ModelType",
      "Split",
      "ImportanceType"
    ),
    names(seed_summary)
  )

  seed_summary <- seed_summary %>%
    group_by(across(all_of(rank_groups))) %>%
    mutate(
      mean_importance_percentile_rank = percent_rank_safe(mean_importance)
    ) %>%
    ungroup() %>%
    select(
      DatasetID,
      CombinationID,
      all_of(ids),
      mean_importance,
      mean_importance_percentile_rank,
      n_seeds_present,
      n_seeds_total,
      seed_presence_frequency,
      n_cv_fits_with_importance,
      n_cv_fits_total,
      overall_cv_fit_presence_frequency
    )

  output_dir <- file.path(seed_summary_dir, dataset_id)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  output_file <- file.path(
    output_dir,
    paste0(combination_id, "__seed_importance_summary.csv")
  )

  write_csv(seed_summary, output_file, na = "")

  tibble(
    DatasetID = dataset_id,
    CombinationID = combination_id,
    InputRunFiles = nrow(file_tbl),
    FeaturesWritten = nrow(seed_summary),
    OutputFile = output_file
  )
}

# ---- Create the summaries used by PreprocessingFeatures_2.R ----

dir.create(seed_summary_dir, recursive = TRUE, showWarnings = FALSE)

dataset_dirs <- list.dirs(
  root_dir,
  recursive = FALSE,
  full.names = TRUE
)

# Exclude the directory where this script writes its results.
dataset_dirs <- dataset_dirs[
  basename(dataset_dirs) != "_feature_importance_seed_summaries"
]

files <- unlist(
  map(
    dataset_dirs,
    ~ list.files(
      path = .x,
      pattern = run_file_re,
      recursive = FALSE,
      full.names = TRUE
    )
  ),
  use.names = FALSE
)

if (length(files) == 0L) {
  stop("No validation-fold feature-importance files found under: ", root_dir)
}

file_index <- tibble(
  path = files,
  file = basename(files),
  dataset_id = basename(dirname(files))
) %>%
  mutate(
    run_num = as.integer(str_match(file, run_file_re)[, 2]),
    combination_id = str_remove(file, run_file_re)
  ) %>%
  filter(!is.na(run_num))

summary_index <- file_index %>%
  group_by(dataset_id, combination_id) %>%
  group_split() %>%
  map_dfr(summarize_one_combination)

print(summary_index)
message("Seed summaries written under: ", seed_summary_dir)

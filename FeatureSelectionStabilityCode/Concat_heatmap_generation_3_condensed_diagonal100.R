library(tidyverse)

# Create the per-dataset pairwise overlap files used by
# binary_summaries_4.R and continuous_summaries_5.R.
# Also export one lower-triangle top-N overlap matrix for each of the five views,
# including 100% diagonal self-comparison cells.

# ================================
# User settings
# ================================

summary_root_dir <- "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/FeatureImportanceLogs/Validation_Feature_Importance/_feature_importance_seed_summaries_with_stable_scores"

top_n_features <- 20L
top_n_label <- paste0("top", top_n_features)
maximum_overlap_percent <- 100

ranking_metric <- "stable_score_seed"
require_positive_stable_score <- FALSE

expected_views <- c(
  "metabolomics",
  "mss",
  "concat_none",
  "concat_scaling",
  "concat_scaling_weights"
)

# Each view is expected to contain the 24 model/preprocessing conditions used
# by binary_summaries_4.R and continuous_summaries_5.R.
expected_conditions_per_view <- 24L

preferred_model_order <- c("enet", "rf", "xgb")
preferred_metab_transform_order <- c("log2", "none")
preferred_taxa_transform_order <- c("clr", "none")
preferred_metab_reduction_order <- c("limma", "none")
preferred_taxa_reduction_order <- c("wilcox", "none")

seed_summary_pattern <- "__seed_importance_summary_with_stable_scores\\.csv$"


# ================================
# Helper functions
# ================================

clean_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9_+-]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_remove("^_") %>%
    str_remove("_$")
}

validate_required_columns <- function(dat, required_cols, file_label = "data frame") {
  missing_cols <- setdiff(required_cols, names(dat))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required column(s) in ", file_label, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  invisible(TRUE)
}

safe_factor <- function(x, preferred_levels) {
  observed <- unique(as.character(x))
  extra <- setdiff(sort(observed), preferred_levels)

  factor(
    as.character(x),
    levels = c(preferred_levels, extra)
  )
}

parse_combination_metadata <- function(dat) {
  combo_parts <- str_match(
    dat$CombinationID,
    "^(.+?)__model-([^_]+)__metabT-([^_]+)__taxaT-([^_]+)__metabR-([^_]+)__taxaR-([^_]+)$"
  )

  invalid_rows <- is.na(combo_parts[, 1])

  if (any(invalid_rows)) {
    invalid_ids <- unique(dat$CombinationID[invalid_rows])

    stop(
      "Could not parse the following CombinationID value(s):\n",
      paste(head(invalid_ids, 10L), collapse = "\n")
    )
  }

  dat %>%
    mutate(
      StudyOutcome = combo_parts[, 2],
      ModelParsed = combo_parts[, 3],
      MetabTransform = combo_parts[, 4],
      TaxaTransform = combo_parts[, 5],
      MetabReduction = combo_parts[, 6],
      TaxaReduction = combo_parts[, 7],
      CategoryID = paste(CombinationID, View, sep = "__VIEW-"),
      CategoryLabel = paste0(
        ModelParsed,
        "\n",
        "mT=", MetabTransform,
        "; tT=", TaxaTransform,
        "\n",
        "mR=", MetabReduction,
        "; tR=", TaxaReduction
      )
    )
}

make_category_info <- function(dat) {
  category_info <- dat %>%
    distinct(
      CategoryID,
      CombinationID,
      View,
      ModelParsed,
      MetabTransform,
      TaxaTransform,
      MetabReduction,
      TaxaReduction,
      CategoryLabel
    )

  conflicting_categories <- category_info %>%
    count(CategoryID, name = "n_metadata_rows") %>%
    filter(n_metadata_rows > 1L)

  if (nrow(conflicting_categories) > 0L) {
    stop("At least one CategoryID maps to conflicting metadata.")
  }

  category_info %>%
    mutate(
      View = safe_factor(View, expected_views),
      ModelParsed = safe_factor(ModelParsed, preferred_model_order),
      MetabTransform = safe_factor(
        MetabTransform,
        preferred_metab_transform_order
      ),
      TaxaTransform = safe_factor(
        TaxaTransform,
        preferred_taxa_transform_order
      ),
      MetabReduction = safe_factor(
        MetabReduction,
        preferred_metab_reduction_order
      ),
      TaxaReduction = safe_factor(
        TaxaReduction,
        preferred_taxa_reduction_order
      )
    ) %>%
    arrange(
      View,
      ModelParsed,
      MetabTransform,
      TaxaTransform,
      MetabReduction,
      TaxaReduction,
      CategoryID
    ) %>%
    mutate(CategoryOrder = row_number())
}

select_top_features <- function(feature_dat, category_info) {
  # Retain only the strongest row if a feature appears more than once within
  # the same CombinationID x View category.
  feature_dat <- feature_dat %>%
    arrange(
      CategoryID,
      Feature,
      desc(.data[[ranking_metric]]),
      desc(seed_presence_frequency),
      desc(stable_score_cv),
      desc(mean_importance_percentile_rank)
    ) %>%
    group_by(CategoryID, Feature) %>%
    slice(1L) %>%
    ungroup()

  if (require_positive_stable_score) {
    eligible_dat <- feature_dat %>%
      filter(
        !is.na(.data[[ranking_metric]]),
        .data[[ranking_metric]] > 0
      )
  } else {
    eligible_dat <- feature_dat %>%
      filter(!is.na(.data[[ranking_metric]]))
  }

  category_status <- category_info %>%
    left_join(
      eligible_dat %>%
        count(CategoryID, name = "n_features_eligible_for_topn"),
      by = "CategoryID"
    ) %>%
    mutate(
      n_features_eligible_for_topn = replace_na(
        n_features_eligible_for_topn,
        0L
      ),
      valid_topn_category =
        n_features_eligible_for_topn >= top_n_features
    )

  top_feature_table <- eligible_dat %>%
    arrange(
      CategoryID,
      desc(.data[[ranking_metric]]),
      desc(seed_presence_frequency),
      desc(stable_score_cv),
      desc(mean_importance_percentile_rank),
      Feature
    ) %>%
    group_by(CategoryID) %>%
    slice_head(n = top_n_features) %>%
    ungroup()

  list(
    top_feature_table = top_feature_table,
    category_status = category_status
  )
}

make_within_view_pairwise_overlap <- function(
    top_feature_table,
    category_status
) {
  top_feature_sets <- top_feature_table %>%
    group_by(CategoryID) %>%
    summarise(
      top_features = list(unique(as.character(Feature))),
      n_top_features_selected = n_distinct(Feature),
      .groups = "drop"
    )

  category_sets <- category_status %>%
    left_join(top_feature_sets, by = "CategoryID") %>%
    mutate(
      top_features = map(
        top_features,
        function(x) {
          if (is.null(x) || length(x) == 0L || all(is.na(x))) {
            character(0)
          } else {
            unique(as.character(x))
          }
        }
      ),
      n_top_features_selected = replace_na(
        n_top_features_selected,
        0L
      ),
      valid_topn_category =
        valid_topn_category &
        n_top_features_selected >= top_n_features
    )

  pairwise_overlap <- map_dfr(
    expected_views,
    function(this_view) {
      this_view_categories <- category_sets %>%
        filter(as.character(View) == this_view)

      category_ids <- this_view_categories$CategoryID

      feature_lookup <- set_names(
        this_view_categories$top_features,
        category_ids
      )

      valid_lookup <- set_names(
        this_view_categories$valid_topn_category,
        category_ids
      )

      crossing(
        CategoryID_A = category_ids,
        CategoryID_B = category_ids
      ) %>%
        mutate(
          valid_comparison =
            valid_lookup[CategoryID_A] &
            valid_lookup[CategoryID_B],
          n_shared_topn = map2_int(
            CategoryID_A,
            CategoryID_B,
            function(category_a, category_b) {
              if (
                !isTRUE(valid_lookup[[category_a]]) ||
                !isTRUE(valid_lookup[[category_b]])
              ) {
                return(NA_integer_)
              }

              length(
                intersect(
                  feature_lookup[[category_a]],
                  feature_lookup[[category_b]]
                )
              )
            }
          ),
          topn_overlap_percent = case_when(
            valid_comparison ~
              pmin(
                maximum_overlap_percent,
                maximum_overlap_percent * n_shared_topn / top_n_features
              ),
            TRUE ~ NA_real_
          )
        )
    }
  )

  metadata_a <- category_status %>%
    select(
      CategoryID,
      CategoryOrder,
      View,
      ModelParsed,
      MetabTransform,
      TaxaTransform,
      MetabReduction,
      TaxaReduction
    ) %>%
    rename_with(~ paste0(.x, "_A"), -CategoryID) %>%
    rename(CategoryID_A = CategoryID)

  metadata_b <- category_status %>%
    select(
      CategoryID,
      CategoryOrder,
      View,
      ModelParsed,
      MetabTransform,
      TaxaTransform,
      MetabReduction,
      TaxaReduction
    ) %>%
    rename_with(~ paste0(.x, "_B"), -CategoryID) %>%
    rename(CategoryID_B = CategoryID)

  pairwise_overlap %>%
    left_join(metadata_a, by = "CategoryID_A") %>%
    left_join(metadata_b, by = "CategoryID_B") %>%
    select(
      CategoryID_A,
      CategoryID_B,
      CategoryOrder_A,
      CategoryOrder_B,
      View_A,
      View_B,
      ModelParsed_A,
      ModelParsed_B,
      MetabTransform_A,
      MetabTransform_B,
      TaxaTransform_A,
      TaxaTransform_B,
      MetabReduction_A,
      MetabReduction_B,
      TaxaReduction_A,
      TaxaReduction_B,
      valid_comparison,
      topn_overlap_percent
    ) %>%
    arrange(CategoryOrder_A, CategoryOrder_B)
}

make_view_overlap_matrix <- function(
    pairwise_overlap,
    category_info,
    view_name
) {
  category_meta <- category_info %>%
    filter(as.character(View) == view_name) %>%
    arrange(CategoryOrder) %>%
    mutate(
      MatrixOrder = row_number(),
      MatrixBaseLabel = str_replace_all(CategoryLabel, "\\n", " | ")
    )

  if (nrow(category_meta) == 0L) {
    stop("No categories found for View = ", view_name)
  }

  n_digits <- nchar(as.character(nrow(category_meta)))

  category_meta <- category_meta %>%
    mutate(
      MatrixLabel = sprintf(
        paste0("%0", n_digits, "d | %s"),
        MatrixOrder,
        MatrixBaseLabel
      )
    )

  matrix_metadata_a <- category_meta %>%
    select(
      CategoryID_A = CategoryID,
      MatrixOrder_A = MatrixOrder,
      MatrixLabel_A = MatrixLabel
    )

  matrix_metadata_b <- category_meta %>%
    select(
      CategoryID_B = CategoryID,
      MatrixOrder_B = MatrixOrder,
      MatrixLabel_B = MatrixLabel
    )

  matrix_values <- pairwise_overlap %>%
    filter(
      as.character(View_A) == view_name,
      as.character(View_B) == view_name
    ) %>%
    inner_join(matrix_metadata_a, by = "CategoryID_A") %>%
    inner_join(matrix_metadata_b, by = "CategoryID_B") %>%
    transmute(
      RowLabel = MatrixLabel_A,
      ColLabel = MatrixLabel_B,
      MatrixValue = case_when(
        MatrixOrder_A == MatrixOrder_B ~ maximum_overlap_percent,
        MatrixOrder_A > MatrixOrder_B & valid_comparison ~
          pmin(topn_overlap_percent, maximum_overlap_percent),
        TRUE ~ NA_real_
      )
    )

  expand_grid(
    RowLabel = category_meta$MatrixLabel,
    ColLabel = category_meta$MatrixLabel
  ) %>%
    left_join(
      matrix_values,
      by = c("RowLabel", "ColLabel")
    ) %>%
    mutate(
      RowLabel = factor(
        RowLabel,
        levels = category_meta$MatrixLabel
      ),
      ColLabel = factor(
        ColLabel,
        levels = category_meta$MatrixLabel
      )
    ) %>%
    arrange(RowLabel, ColLabel) %>%
    pivot_wider(
      names_from = ColLabel,
      values_from = MatrixValue
    ) %>%
    mutate(
      PipelineCondition = as.character(RowLabel),
      .before = 1
    ) %>%
    select(
      PipelineCondition,
      all_of(category_meta$MatrixLabel)
    )
}

remove_stale_core_outputs <- function(output_dir) {
  # Remove only the pairwise input and base view matrices generated by this
  # condensed script. Other historical files are left unchanged.
  pairwise_files <- list.files(
    path = output_dir,
    pattern = paste0("_pairwise_", top_n_label, "_overlap\\.csv$"),
    full.names = TRUE
  )

  view_matrix_files <- list.files(
    path = output_dir,
    pattern = paste0(
      "_",
      top_n_label,
      "_overlap_VIEW_.+_[0-9]+x[0-9]+\\.csv$"
    ),
    full.names = TRUE
  )

  stale_files <- unique(c(pairwise_files, view_matrix_files))

  if (length(stale_files) > 0L) {
    unlink(stale_files)
  }
}


# ================================
# Process one dataset directory
# ================================

process_one_dataset_dir <- function(dataset_dir) {
  message("\nProcessing: ", basename(dataset_dir))

  seed_summary_files <- list.files(
    path = dataset_dir,
    pattern = seed_summary_pattern,
    full.names = TRUE
  )

  if (length(seed_summary_files) == 0L) {
    stop("No stable-score seed-summary files found in:\n", dataset_dir)
  }

  all_features <- map_dfr(
    seed_summary_files,
    function(path) {
      read_csv(
        path,
        show_col_types = FALSE,
        progress = FALSE
      )
    }
  )

  required_cols <- c(
    "CombinationID",
    "View",
    "Feature",
    "stable_score_seed",
    "seed_presence_frequency",
    "mean_importance_percentile_rank"
  )

  validate_required_columns(
    all_features,
    required_cols,
    basename(dataset_dir)
  )

  if (!"stable_score_cv" %in% names(all_features)) {
    all_features <- all_features %>%
      mutate(stable_score_cv = NA_real_)
  }

  all_features <- all_features %>%
    mutate(
      View = as.character(View),
      Feature = str_trim(as.character(Feature))
    )

  if (any(is.na(all_features$Feature) | all_features$Feature == "")) {
    stop("Missing or empty Feature values were found in ", dataset_dir)
  }

  all_features <- all_features %>%
    parse_combination_metadata()

  study_outcomes <- unique(all_features$StudyOutcome)

  if (length(study_outcomes) != 1L) {
    stop(
      "Expected exactly one study/outcome identifier in ",
      dataset_dir,
      ", but found: ",
      paste(study_outcomes, collapse = ", ")
    )
  }

  observed_views <- sort(unique(all_features$View))
  missing_views <- setdiff(expected_views, observed_views)
  unexpected_views <- setdiff(observed_views, expected_views)

  if (length(missing_views) > 0L || length(unexpected_views) > 0L) {
    stop(
      "The dataset must contain exactly these five views: ",
      paste(expected_views, collapse = ", "),
      ". Missing: ",
      ifelse(length(missing_views) == 0L, "none", paste(missing_views, collapse = ", ")),
      ". Unexpected: ",
      ifelse(length(unexpected_views) == 0L, "none", paste(unexpected_views, collapse = ", "))
    )
  }

  category_info <- make_category_info(all_features)

  view_condition_counts <- category_info %>%
    mutate(View = as.character(View)) %>%
    count(View, name = "n_conditions")

  invalid_view_counts <- view_condition_counts %>%
    filter(n_conditions != expected_conditions_per_view)

  if (nrow(invalid_view_counts) > 0L) {
    count_message <- invalid_view_counts %>%
      transmute(
        text = paste0(View, "=", n_conditions)
      ) %>%
      pull(text) %>%
      paste(collapse = ", ")

    stop(
      "Expected ", expected_conditions_per_view,
      " conditions in each view for ", basename(dataset_dir),
      ", but found: ", count_message
    )
  }

  top_results <- select_top_features(
    feature_dat = all_features,
    category_info = category_info
  )

  top_feature_table <- top_results$top_feature_table
  category_status <- top_results$category_status

  invalid_categories <- category_status %>%
    filter(!valid_topn_category)

  if (nrow(invalid_categories) > 0L) {
    warning(
      nrow(invalid_categories),
      " categories in ",
      basename(dataset_dir),
      " had fewer than ",
      top_n_features,
      " eligible features. Their overlap cells will be blank."
    )
  }

  pairwise_overlap <- make_within_view_pairwise_overlap(
    top_feature_table = top_feature_table,
    category_status = category_status
  )

  downstream_required_cols <- c(
    "CategoryID_A",
    "CategoryID_B",
    "CategoryOrder_A",
    "CategoryOrder_B",
    "View_A",
    "View_B",
    "ModelParsed_A",
    "ModelParsed_B",
    "MetabTransform_A",
    "MetabTransform_B",
    "TaxaTransform_A",
    "TaxaTransform_B",
    "MetabReduction_A",
    "MetabReduction_B",
    "TaxaReduction_A",
    "TaxaReduction_B",
    "valid_comparison",
    "topn_overlap_percent"
  )

  validate_required_columns(
    pairwise_overlap,
    downstream_required_cols,
    "pairwise overlap output"
  )

  dataset_file_prefix <- clean_filename(study_outcomes[[1]])

  output_dir <- file.path(
    dataset_dir,
    paste0("_", top_n_label, "_stable_score_seed_overlap")
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  remove_stale_core_outputs(output_dir)

  pairwise_file <- file.path(
    output_dir,
    paste0(
      dataset_file_prefix,
      "_pairwise_",
      top_n_label,
      "_overlap.csv"
    )
  )

  write_csv(pairwise_overlap, pairwise_file, na = "")

  view_matrix_files <- map_chr(
    expected_views,
    function(this_view) {
      view_matrix <- make_view_overlap_matrix(
        pairwise_overlap = pairwise_overlap,
        category_info = category_info,
        view_name = this_view
      )

      n_conditions <- nrow(view_matrix)
      dimension_label <- paste0(n_conditions, "x", n_conditions)

      matrix_file <- file.path(
        output_dir,
        paste0(
          dataset_file_prefix,
          "_",
          top_n_label,
          "_overlap_VIEW_",
          clean_filename(this_view),
          "_",
          dimension_label,
          ".csv"
        )
      )

      write_csv(view_matrix, matrix_file, na = "")
      matrix_file
    }
  )

  if (length(view_matrix_files) != length(expected_views)) {
    stop("Did not create all five view-specific overlap matrices.")
  }

  message("Pairwise input: ", basename(pairwise_file))
  message("View matrices written: ", length(view_matrix_files))

  tibble(
    DatasetDirectory = basename(dataset_dir),
    DatasetFilePrefix = dataset_file_prefix,
    InputFiles = length(seed_summary_files),
    Categories = nrow(category_info),
    InvalidCategories = nrow(invalid_categories),
    PairwiseFile = pairwise_file,
    ViewMatrixFiles = paste(view_matrix_files, collapse = "; "),
    Status = "OK",
    ErrorMessage = NA_character_
  )
}


# ================================
# Batch process all dataset directories
# ================================

dataset_dirs <- list.dirs(
  path = summary_root_dir,
  recursive = FALSE,
  full.names = TRUE
)

dataset_dirs <- dataset_dirs[
  map_lgl(
    dataset_dirs,
    function(path) {
      length(
        list.files(
          path = path,
          pattern = seed_summary_pattern,
          full.names = TRUE
        )
      ) > 0L
    }
  )
]

if (length(dataset_dirs) == 0L) {
  stop(
    "No dataset directories containing stable-score seed-summary files were found under:\n",
    summary_root_dir
  )
}

batch_results <- map_dfr(
  dataset_dirs,
  function(dataset_dir) {
    tryCatch(
      process_one_dataset_dir(dataset_dir),
      error = function(e) {
        tibble(
          DatasetDirectory = basename(dataset_dir),
          DatasetFilePrefix = NA_character_,
          InputFiles = NA_integer_,
          Categories = NA_integer_,
          InvalidCategories = NA_integer_,
          PairwiseFile = NA_character_,
          ViewMatrixFiles = NA_character_,
          Status = "ERROR",
          ErrorMessage = conditionMessage(e)
        )
      }
    )
  }
)

print(batch_results, width = Inf)

if (any(batch_results$Status == "ERROR")) {
  stop("One or more dataset directories could not be processed. See the printed results.")
}

message("\nDone.")
message(
  "Each dataset now contains one pairwise overlap input file and five view-specific ",
  top_n_label,
  " overlap matrices."
)

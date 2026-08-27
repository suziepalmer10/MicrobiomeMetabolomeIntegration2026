library(tidyverse)

# Create per-dataset pairwise feature-overlap files for the two single-omics
# views. Metabolomics conditions are defined only by the metabolomics
# transformation and reduction settings; taxa conditions (input View = "mss")
# are defined only by the taxa transformation and reduction settings.
#
# Therefore, each view contains 12 conditions:
#   3 models x 4 relevant preprocessing combinations.
#
# The exported lower-triangle matrices include 100% diagonal self-comparisons.

# ================================
# User settings
# ================================

summary_root_dir <- "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/FeatureImportanceLogs/Validation_Feature_Importance/_feature_importance_seed_summaries_with_stable_scores"

top_n_features <- 20L
top_n_label <- paste0("top", top_n_features)
maximum_overlap_percent <- 100

ranking_metric <- "stable_score_seed"
require_positive_stable_score <- FALSE

# The source feature-importance files use View = "mss" for the taxa view.
expected_views <- c("metabolomics", "mss")
view_display_labels <- c(
  "metabolomics" = "Metabolomics",
  "mss" = "Taxa"
)

expected_conditions_per_view <- 12L
not_applicable_label <- "not_applicable"

preferred_model_order <- c("enet", "rf", "xgb")

# Exact single-omics preprocessing order requested for the heatmaps:
# metabolomics: (log2, none), (log2, limma), (none, none), (none, limma)
# taxa:         (clr, none),  (clr, wilcox), (none, none), (none, wilcox)
preferred_metab_transform_order <- c("log2", "none")
preferred_metab_reduction_order <- c("none", "limma")
preferred_taxa_transform_order <- c("clr", "none")
preferred_taxa_reduction_order <- c("none", "wilcox")

seed_summary_pattern <- "__seed_importance_summary_with_stable_scores\\.csv$"

single_omics_overlap_dir_name <- paste0(
  "_",
  top_n_label,
  "_single_omics_stable_score_seed_overlap"
)


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

make_single_omics_condition_id <- function(
    view,
    model,
    metab_t,
    taxa_t,
    metab_r,
    taxa_r
) {
  view <- as.character(view)

  case_when(
    view == "metabolomics" ~ paste0(
      "model-", model,
      "__metabT-", metab_t,
      "__metabR-", metab_r
    ),
    view == "mss" ~ paste0(
      "model-", model,
      "__taxaT-", taxa_t,
      "__taxaR-", taxa_r
    ),
    TRUE ~ NA_character_
  )
}

make_single_omics_condition_label <- function(
    view,
    model,
    metab_t,
    taxa_t,
    metab_r,
    taxa_r
) {
  view <- as.character(view)

  case_when(
    view == "metabolomics" ~ paste0(
      model,
      "\n(", metab_t, ", ", metab_r, ")"
    ),
    view == "mss" ~ paste0(
      model,
      "\n(", taxa_t, ", ", taxa_r, ")"
    ),
    TRUE ~ NA_character_
  )
}

parse_combination_metadata <- function(dat) {
  combo_parts <- str_match(
    as.character(dat$CombinationID),
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
      SourceCombinationID = as.character(CombinationID),
      StudyOutcome = combo_parts[, 2],
      ModelParsed = combo_parts[, 3],
      SourceMetabTransform = combo_parts[, 4],
      SourceTaxaTransform = combo_parts[, 5],
      SourceMetabReduction = combo_parts[, 6],
      SourceTaxaReduction = combo_parts[, 7],
      View = as.character(View)
    ) %>%
    filter(View %in% expected_views) %>%
    mutate(
      # Irrelevant-omics settings are replaced with a common value so they
      # cannot create additional single-omics conditions.
      MetabTransform = if_else(
        View == "metabolomics",
        SourceMetabTransform,
        not_applicable_label
      ),
      TaxaTransform = if_else(
        View == "mss",
        SourceTaxaTransform,
        not_applicable_label
      ),
      MetabReduction = if_else(
        View == "metabolomics",
        SourceMetabReduction,
        not_applicable_label
      ),
      TaxaReduction = if_else(
        View == "mss",
        SourceTaxaReduction,
        not_applicable_label
      ),
      SingleOmicsConditionID = make_single_omics_condition_id(
        View,
        ModelParsed,
        MetabTransform,
        TaxaTransform,
        MetabReduction,
        TaxaReduction
      ),
      SingleOmicsCombinationID = paste0(
        StudyOutcome,
        "__",
        SingleOmicsConditionID
      ),
      CategoryID = paste(
        SingleOmicsCombinationID,
        View,
        sep = "__VIEW-"
      ),
      CategoryLabel = make_single_omics_condition_label(
        View,
        ModelParsed,
        MetabTransform,
        TaxaTransform,
        MetabReduction,
        TaxaReduction
      )
    )
}

make_category_info <- function(dat) {
  category_info <- dat %>%
    distinct(
      CategoryID,
      SingleOmicsCombinationID,
      SingleOmicsConditionID,
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
    stop("At least one single-omics CategoryID maps to conflicting metadata.")
  }

  category_info %>%
    mutate(
      View = safe_factor(View, expected_views),
      ModelParsed = safe_factor(ModelParsed, preferred_model_order),
      MetabTransform = safe_factor(
        MetabTransform,
        c(preferred_metab_transform_order, not_applicable_label)
      ),
      TaxaTransform = safe_factor(
        TaxaTransform,
        c(preferred_taxa_transform_order, not_applicable_label)
      ),
      MetabReduction = safe_factor(
        MetabReduction,
        c(preferred_metab_reduction_order, not_applicable_label)
      ),
      TaxaReduction = safe_factor(
        TaxaReduction,
        c(preferred_taxa_reduction_order, not_applicable_label)
      ),
      RelevantTransformOrder = case_when(
        as.character(View) == "metabolomics" ~ match(
          as.character(MetabTransform),
          preferred_metab_transform_order
        ),
        as.character(View) == "mss" ~ match(
          as.character(TaxaTransform),
          preferred_taxa_transform_order
        ),
        TRUE ~ NA_integer_
      ),
      RelevantReductionOrder = case_when(
        as.character(View) == "metabolomics" ~ match(
          as.character(MetabReduction),
          preferred_metab_reduction_order
        ),
        as.character(View) == "mss" ~ match(
          as.character(TaxaReduction),
          preferred_taxa_reduction_order
        ),
        TRUE ~ NA_integer_
      )
    ) %>%
    arrange(
      View,
      ModelParsed,
      RelevantTransformOrder,
      RelevantReductionOrder,
      CategoryID
    ) %>%
    group_by(View) %>%
    mutate(CategoryOrder = row_number()) %>%
    ungroup()
}

make_expected_condition_grid <- function() {
  metabolomics_grid <- expand_grid(
    ModelParsed = preferred_model_order,
    MetabTransform = preferred_metab_transform_order,
    MetabReduction = preferred_metab_reduction_order
  ) %>%
    mutate(
      View = "metabolomics",
      TaxaTransform = not_applicable_label,
      TaxaReduction = not_applicable_label
    )

  taxa_grid <- expand_grid(
    ModelParsed = preferred_model_order,
    TaxaTransform = preferred_taxa_transform_order,
    TaxaReduction = preferred_taxa_reduction_order
  ) %>%
    mutate(
      View = "mss",
      MetabTransform = not_applicable_label,
      MetabReduction = not_applicable_label
    )

  bind_rows(metabolomics_grid, taxa_grid) %>%
    mutate(
      SingleOmicsConditionID = make_single_omics_condition_id(
        View,
        ModelParsed,
        MetabTransform,
        TaxaTransform,
        MetabReduction,
        TaxaReduction
      )
    ) %>%
    select(View, SingleOmicsConditionID)
}

validate_expected_conditions <- function(category_info, dataset_label) {
  expected_conditions <- make_expected_condition_grid()

  observed_conditions <- category_info %>%
    transmute(
      View = as.character(View),
      SingleOmicsConditionID = as.character(SingleOmicsConditionID)
    ) %>%
    distinct()

  missing_conditions <- anti_join(
    expected_conditions,
    observed_conditions,
    by = c("View", "SingleOmicsConditionID")
  )

  unexpected_conditions <- anti_join(
    observed_conditions,
    expected_conditions,
    by = c("View", "SingleOmicsConditionID")
  )

  condition_counts <- observed_conditions %>%
    count(View, name = "n_conditions")

  invalid_counts <- condition_counts %>%
    filter(n_conditions != expected_conditions_per_view)

  if (
    nrow(missing_conditions) > 0L ||
      nrow(unexpected_conditions) > 0L ||
      nrow(invalid_counts) > 0L
  ) {
    stop(
      "The single-omics condition grid was not the expected 12 conditions per view in ",
      dataset_label,
      ".\nMissing conditions: ",
      ifelse(
        nrow(missing_conditions) == 0L,
        "none",
        paste(
          paste0(
            missing_conditions$View,
            ":",
            missing_conditions$SingleOmicsConditionID
          ),
          collapse = "; "
        )
      ),
      "\nUnexpected conditions: ",
      ifelse(
        nrow(unexpected_conditions) == 0L,
        "none",
        paste(
          paste0(
            unexpected_conditions$View,
            ":",
            unexpected_conditions$SingleOmicsConditionID
          ),
          collapse = "; "
        )
      )
    )
  }

  invisible(condition_counts)
}

choose_representative_source_combinations <- function(feature_dat) {
  # A single-omics model may be repeated under multiple settings for the other
  # omics block. Select one source CombinationID per canonical condition so
  # those repeated copies are not counted as distinct models. When possible,
  # prefer "none" for the irrelevant transformation and reduction settings.
  source_map <- feature_dat %>%
    distinct(
      CategoryID,
      SingleOmicsCombinationID,
      SingleOmicsConditionID,
      View,
      ModelParsed,
      MetabTransform,
      TaxaTransform,
      MetabReduction,
      TaxaReduction,
      SourceCombinationID,
      SourceMetabTransform,
      SourceTaxaTransform,
      SourceMetabReduction,
      SourceTaxaReduction
    ) %>%
    mutate(
      IrrelevantTransformPenalty = case_when(
        View == "metabolomics" ~ as.integer(SourceTaxaTransform != "none"),
        View == "mss" ~ as.integer(SourceMetabTransform != "none"),
        TRUE ~ 1L
      ),
      IrrelevantReductionPenalty = case_when(
        View == "metabolomics" ~ as.integer(SourceTaxaReduction != "none"),
        View == "mss" ~ as.integer(SourceMetabReduction != "none"),
        TRUE ~ 1L
      )
    ) %>%
    arrange(
      CategoryID,
      IrrelevantTransformPenalty,
      IrrelevantReductionPenalty,
      SourceCombinationID
    ) %>%
    group_by(CategoryID) %>%
    mutate(
      RepresentativeSourceRank = row_number(),
      IsRepresentativeSource = RepresentativeSourceRank == 1L
    ) %>%
    ungroup()

  representative_sources <- source_map %>%
    filter(IsRepresentativeSource) %>%
    select(CategoryID, SourceCombinationID)

  list(
    source_map = source_map,
    representative_sources = representative_sources
  )
}

make_redundant_source_topn_consistency <- function(feature_dat) {
  # Verify that source combinations differing only in settings from the unused
  # omics block produce the same top-N feature set. The heatmap still
  # uses one deterministic representative source combination, but any
  # unexpected differences are recorded explicitly.
  source_feature_dat <- feature_dat %>%
    arrange(
      CategoryID,
      SourceCombinationID,
      Feature,
      desc(.data[[ranking_metric]]),
      desc(seed_presence_frequency),
      desc(stable_score_cv),
      desc(mean_importance_percentile_rank)
    ) %>%
    group_by(CategoryID, SourceCombinationID, Feature) %>%
    slice(1L) %>%
    ungroup()

  if (require_positive_stable_score) {
    eligible_source_features <- source_feature_dat %>%
      filter(
        !is.na(.data[[ranking_metric]]),
        .data[[ranking_metric]] > 0
      )
  } else {
    eligible_source_features <- source_feature_dat %>%
      filter(!is.na(.data[[ranking_metric]]))
  }

  source_topn_signatures <- eligible_source_features %>%
    arrange(
      CategoryID,
      SourceCombinationID,
      desc(.data[[ranking_metric]]),
      desc(seed_presence_frequency),
      desc(stable_score_cv),
      desc(mean_importance_percentile_rank),
      Feature
    ) %>%
    group_by(CategoryID, SourceCombinationID) %>%
    slice_head(n = top_n_features) %>%
    summarise(
      n_top_features = n_distinct(Feature),
      topn_signature = paste(
        sort(unique(as.character(Feature))),
        collapse = " || "
      ),
      .groups = "drop"
    )

  source_grid <- source_feature_dat %>%
    distinct(CategoryID, SourceCombinationID) %>%
    left_join(
      source_topn_signatures,
      by = c("CategoryID", "SourceCombinationID")
    ) %>%
    mutate(
      n_top_features = replace_na(n_top_features, 0L),
      valid_source_topn = n_top_features >= top_n_features
    )

  category_metadata <- feature_dat %>%
    distinct(
      CategoryID,
      View,
      ModelParsed,
      MetabTransform,
      TaxaTransform,
      MetabReduction,
      TaxaReduction,
      CategoryLabel
    )

  source_grid %>%
    group_by(CategoryID) %>%
    summarise(
      n_source_combinations = n_distinct(SourceCombinationID),
      n_valid_source_topn_sets = sum(valid_source_topn),
      n_unique_valid_topn_sets = n_distinct(
        topn_signature[valid_source_topn],
        na.rm = TRUE
      ),
      all_redundant_topn_sets_identical =
        n_source_combinations >= 1L &
        n_valid_source_topn_sets == n_source_combinations &
        n_unique_valid_topn_sets == 1L,
      .groups = "drop"
    ) %>%
    left_join(category_metadata, by = "CategoryID") %>%
    mutate(View = safe_factor(View, expected_views)) %>%
    arrange(View, ModelParsed, CategoryID)
}

select_top_features <- function(feature_dat, category_info) {
  # Retain only the strongest row if a feature appears more than once within
  # the selected representative CategoryID.
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
    arrange(View_A, CategoryOrder_A, CategoryOrder_B)
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

remove_stale_single_omics_outputs <- function(output_dir) {
  stale_files <- list.files(
    path = output_dir,
    pattern = paste0(
      "(_single_omics_pairwise_",
      top_n_label,
      "_overlap\\.csv$)|(_single_omics_",
      top_n_label,
      "_overlap_VIEW_.+_[0-9]+x[0-9]+\\.csv$)"
    ),
    full.names = TRUE
  )

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

  observed_source_views <- sort(unique(all_features$View))
  missing_single_omics_views <- setdiff(expected_views, observed_source_views)

  if (length(missing_single_omics_views) > 0L) {
    stop(
      "Missing required single-omics view(s) in ",
      basename(dataset_dir),
      ": ",
      paste(missing_single_omics_views, collapse = ", ")
    )
  }

  # parse_combination_metadata() also drops concatenated views.
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

  category_info <- make_category_info(all_features)
  condition_counts <- validate_expected_conditions(
    category_info,
    basename(dataset_dir)
  )

  redundant_source_consistency <-
    make_redundant_source_topn_consistency(all_features)

  inconsistent_redundant_categories <- redundant_source_consistency %>%
    filter(!all_redundant_topn_sets_identical)

  if (nrow(inconsistent_redundant_categories) > 0L) {
    warning(
      nrow(inconsistent_redundant_categories),
      " single-omics categories in ",
      basename(dataset_dir),
      " did not have identical ",
      top_n_label,
      " feature sets across the available redundant source combinations. ",
      "The representative-source map and consistency CSV identify these cases."
    )
  }

  source_selection <- choose_representative_source_combinations(all_features)
  source_map <- source_selection$source_map
  representative_sources <- source_selection$representative_sources

  representative_features <- all_features %>%
    semi_join(
      representative_sources,
      by = c("CategoryID", "SourceCombinationID")
    )

  top_results <- select_top_features(
    feature_dat = representative_features,
    category_info = category_info
  )

  top_feature_table <- top_results$top_feature_table
  category_status <- top_results$category_status

  invalid_categories <- category_status %>%
    filter(!valid_topn_category)

  if (nrow(invalid_categories) > 0L) {
    warning(
      nrow(invalid_categories),
      " single-omics categories in ",
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
    "single-omics pairwise overlap output"
  )

  dataset_file_prefix <- clean_filename(study_outcomes[[1]])

  output_dir <- file.path(
    dataset_dir,
    single_omics_overlap_dir_name
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  remove_stale_single_omics_outputs(output_dir)

  pairwise_file <- file.path(
    output_dir,
    paste0(
      dataset_file_prefix,
      "_single_omics_pairwise_",
      top_n_label,
      "_overlap.csv"
    )
  )

  condition_info_file <- file.path(
    output_dir,
    paste0(
      dataset_file_prefix,
      "_single_omics_condition_info.csv"
    )
  )

  source_map_file <- file.path(
    output_dir,
    paste0(
      dataset_file_prefix,
      "_single_omics_representative_source_map.csv"
    )
  )

  consistency_file <- file.path(
    output_dir,
    paste0(
      dataset_file_prefix,
      "_single_omics_redundant_source_topn_consistency.csv"
    )
  )

  top_feature_file <- file.path(
    output_dir,
    paste0(
      dataset_file_prefix,
      "_single_omics_",
      top_n_label,
      "_features.csv"
    )
  )

  write_csv(pairwise_overlap, pairwise_file, na = "")
  write_csv(category_info, condition_info_file, na = "")
  write_csv(source_map, source_map_file, na = "")
  write_csv(redundant_source_consistency, consistency_file, na = "")
  write_csv(top_feature_table, top_feature_file, na = "")

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
      output_view_name <- clean_filename(view_display_labels[[this_view]])

      matrix_file <- file.path(
        output_dir,
        paste0(
          dataset_file_prefix,
          "_single_omics_",
          top_n_label,
          "_overlap_VIEW_",
          output_view_name,
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
    stop("Did not create both single-omics view-specific overlap matrices.")
  }

  message("Single-omics pairwise input: ", basename(pairwise_file))
  message("Single-omics view matrices written: ", length(view_matrix_files))
  print(condition_counts)

  tibble(
    DatasetDirectory = basename(dataset_dir),
    DatasetFilePrefix = dataset_file_prefix,
    InputFiles = length(seed_summary_files),
    SingleOmicsCategories = nrow(category_info),
    RepresentativeSourceModels = nrow(representative_sources),
    InconsistentRedundantCategories =
      nrow(inconsistent_redundant_categories),
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
          SingleOmicsCategories = NA_integer_,
          RepresentativeSourceModels = NA_integer_,
          InconsistentRedundantCategories = NA_integer_,
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
  "Each dataset now contains one single-omics pairwise overlap file and two ",
  top_n_label,
  " overlap matrices: a 12 x 12 metabolomics matrix and a 12 x 12 taxa matrix."
)

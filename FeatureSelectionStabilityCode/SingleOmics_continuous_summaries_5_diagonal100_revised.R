library(tidyverse)

# Aggregate the per-dataset single-omics pairwise overlap files across
# continuous endpoints. Each heatmap contains 12 conditions per view:
# three models x four relevant preprocessing combinations.
#
# Metabolomics conditions ignore taxa preprocessing settings.
# Taxa conditions (input View = "mss") ignore metabolomics preprocessing settings.
# Diagonal self-comparisons are displayed as 100% squares.

# ================================
# User settings
# ================================

summary_root_dir <- "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/FeatureImportanceLogs/Validation_Feature_Importance/_feature_importance_seed_summaries_with_stable_scores"

top_n_features <- 20L

top_n_label <- paste0("top", top_n_features)
top_n_title <- paste0("Top-", top_n_features)
top_n_lower <- paste0("top-", top_n_features)

triangle_to_show <- "lower_left"
show_diagonal <- TRUE
heatmap_max_percent <- 100

endpoint_type <- "continuous"
endpoint_title <- "Continuous"

dataset_subdirs <- c(
  "Era_Cholesterol",
  "Era_Glucose",
  "Franzosa_CD_Fp",
  "Franzosa_UC_Fp",
  "Franzosa_IBD_Fp",
  "Wang_creatinine",
  "Wang_eGFR",
  "Wang_urea"
)

# The source pairwise files use View = "mss" for the taxa view.
expected_views <- c("metabolomics", "mss")
view_display_labels <- c(
  "metabolomics" = "Metabolomics",
  "mss" = "Taxa"
)
view_output_names <- c(
  "metabolomics" = "metabolomics",
  "mss" = "taxa"
)

view_heatmap_colors <- c(
  "metabolomics" = "#9D1B28",
  "mss" = "#0152A1"
)

not_applicable_label <- "not_applicable"
expected_conditions_per_view <- 12L

preferred_model_order <- c("enet", "rf", "xgb")

# Exact single-omics preprocessing order requested for the heatmaps:
# metabolomics: (log2, none), (log2, limma), (none, none), (none, limma)
# taxa:         (clr, none),  (clr, wilcox), (none, none), (none, wilcox)
preferred_metab_transform_order <- c("log2", "none")
preferred_metab_reduction_order <- c("none", "limma")
preferred_taxa_transform_order <- c("clr", "none")
preferred_taxa_reduction_order <- c("none", "wilcox")

single_omics_overlap_dir_name <- paste0(
  "_",
  top_n_label,
  "_single_omics_stable_score_seed_overlap"
)

pairwise_pattern <- paste0(
  "_single_omics_pairwise_",
  top_n_label,
  "_overlap\\.csv$"
)

output_dir <- file.path(
  summary_root_dir,
  paste0(
    "_",
    endpoint_type,
    "_",
    top_n_label,
    "_aggregate_single_omics_12x12_condition_heatmaps"
  )
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1.0 means a heatmap cell is colored only if all loaded datasets have a
# valid value. Lower values permit cells supported by a subset of datasets.
min_valid_dataset_fraction <- 1.0
show_numbers_in_heatmaps <- TRUE

heatmap_width <- 10
heatmap_height <- 8.5


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

pretty_dataset_label <- function(x) {
  x %>%
    str_replace_all("_", " ")
}

get_view_display_label <- function(view_name) {
  view_name <- as.character(view_name)

  if (view_name %in% names(view_display_labels)) {
    view_display_labels[[view_name]]
  } else {
    view_name
  }
}

get_view_output_name <- function(view_name) {
  view_name <- as.character(view_name)

  if (view_name %in% names(view_output_names)) {
    view_output_names[[view_name]]
  } else {
    clean_filename(view_name)
  }
}

get_view_heatmap_color <- function(view_name) {
  view_name <- as.character(view_name)

  if (view_name %in% names(view_heatmap_colors)) {
    view_heatmap_colors[[view_name]]
  } else {
    warning(
      "No heatmap color specified for View = ",
      view_name,
      ". Using default blue."
    )
    "#2166AC"
  }
}

safe_mean <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(!is.na(x)) < 2L) return(NA_real_)
  sd(x, na.rm = TRUE)
}

safe_iqr <- function(x) {
  if (sum(!is.na(x)) < 2L) return(NA_real_)
  IQR(x, na.rm = TRUE)
}

safe_min <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (length(x) == 0L || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

validate_required_columns <- function(dat, required_cols, file_label = "data frame") {
  missing_cols <- setdiff(required_cols, names(dat))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required column(s) in ",
      file_label,
      ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  invisible(TRUE)
}

get_single_file <- function(path, pattern) {
  files <- list.files(
    path = path,
    pattern = pattern,
    full.names = TRUE
  )

  if (length(files) == 0L) {
    return(NA_character_)
  }

  if (length(files) > 1L) {
    warning(
      "More than one matching file found in:\n",
      path,
      "\nUsing the first file:\n",
      files[[1]]
    )
  }

  files[[1]]
}

standardize_pairwise_overlap <- function(dat, file_label) {
  overlap_col <- case_when(
    "topn_overlap_percent" %in% names(dat) ~ "topn_overlap_percent",
    paste0(top_n_label, "_overlap_percent") %in% names(dat) ~
      paste0(top_n_label, "_overlap_percent"),
    "top20_overlap_percent" %in% names(dat) ~ "top20_overlap_percent",
    "top10_overlap_percent" %in% names(dat) ~ "top10_overlap_percent",
    TRUE ~ NA_character_
  )

  if (is.na(overlap_col)) {
    stop("Could not find an overlap-percent column in: ", file_label)
  }

  required_cols <- c(
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
    "valid_comparison"
  )

  validate_required_columns(dat, required_cols, file_label)

  dat %>%
    mutate(
      overlap_percent = as.numeric(.data[[overlap_col]]),
      valid_comparison =
        as.character(valid_comparison) %in% c("TRUE", "true", "1")
    )
}

make_condition_id <- function(
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

make_condition_label <- function(
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

canonical_metab_transform <- function(view, x) {
  if_else(
    as.character(view) == "metabolomics",
    as.character(x),
    not_applicable_label
  )
}

canonical_taxa_transform <- function(view, x) {
  if_else(
    as.character(view) == "mss",
    as.character(x),
    not_applicable_label
  )
}

canonical_metab_reduction <- function(view, x) {
  if_else(
    as.character(view) == "metabolomics",
    as.character(x),
    not_applicable_label
  )
}

canonical_taxa_reduction <- function(view, x) {
  if_else(
    as.character(view) == "mss",
    as.character(x),
    not_applicable_label
  )
}

build_expected_condition_info <- function() {
  metabolomics_grid <- expand_grid(
    ModelParsed = preferred_model_order,
    MetabTransform = preferred_metab_transform_order,
    MetabReduction = preferred_metab_reduction_order
  ) %>%
    mutate(
      View = "metabolomics",
      TaxaTransform = not_applicable_label,
      TaxaReduction = not_applicable_label,
      RelevantTransformOrder = match(
        MetabTransform,
        preferred_metab_transform_order
      ),
      RelevantReductionOrder = match(
        MetabReduction,
        preferred_metab_reduction_order
      )
    )

  taxa_grid <- expand_grid(
    ModelParsed = preferred_model_order,
    TaxaTransform = preferred_taxa_transform_order,
    TaxaReduction = preferred_taxa_reduction_order
  ) %>%
    mutate(
      View = "mss",
      MetabTransform = not_applicable_label,
      MetabReduction = not_applicable_label,
      RelevantTransformOrder = match(
        TaxaTransform,
        preferred_taxa_transform_order
      ),
      RelevantReductionOrder = match(
        TaxaReduction,
        preferred_taxa_reduction_order
      )
    )

  bind_rows(metabolomics_grid, taxa_grid) %>%
    mutate(
      ViewOrder = match(View, expected_views),
      ModelOrder = match(ModelParsed, preferred_model_order),
      ConditionID = make_condition_id(
        View,
        ModelParsed,
        MetabTransform,
        TaxaTransform,
        MetabReduction,
        TaxaReduction
      ),
      ConditionLabel = make_condition_label(
        View,
        ModelParsed,
        MetabTransform,
        TaxaTransform,
        MetabReduction,
        TaxaReduction
      )
    ) %>%
    arrange(
      ViewOrder,
      ModelOrder,
      RelevantTransformOrder,
      RelevantReductionOrder
    ) %>%
    group_by(View) %>%
    mutate(ConditionOrder = row_number()) %>%
    ungroup() %>%
    select(
      View,
      ConditionID,
      ConditionLabel,
      ConditionOrder,
      ModelParsed,
      MetabTransform,
      TaxaTransform,
      MetabReduction,
      TaxaReduction
    )
}

extract_observed_condition_info <- function(pairwise_dat) {
  condition_info_a <- pairwise_dat %>%
    transmute(
      View = as.character(View_A),
      ConditionID = make_condition_id(
        View_A,
        ModelParsed_A,
        canonical_metab_transform(View_A, MetabTransform_A),
        canonical_taxa_transform(View_A, TaxaTransform_A),
        canonical_metab_reduction(View_A, MetabReduction_A),
        canonical_taxa_reduction(View_A, TaxaReduction_A)
      )
    )

  condition_info_b <- pairwise_dat %>%
    transmute(
      View = as.character(View_B),
      ConditionID = make_condition_id(
        View_B,
        ModelParsed_B,
        canonical_metab_transform(View_B, MetabTransform_B),
        canonical_taxa_transform(View_B, TaxaTransform_B),
        canonical_metab_reduction(View_B, MetabReduction_B),
        canonical_taxa_reduction(View_B, TaxaReduction_B)
      )
    )

  bind_rows(condition_info_a, condition_info_b) %>%
    filter(View %in% expected_views, !is.na(ConditionID)) %>%
    distinct()
}

validate_observed_condition_grid <- function(pairwise_dat, expected_info) {
  observed_info <- extract_observed_condition_info(pairwise_dat)
  expected_ids <- expected_info %>%
    select(View, ConditionID) %>%
    distinct()

  missing_conditions <- anti_join(
    expected_ids,
    observed_info,
    by = c("View", "ConditionID")
  )

  unexpected_conditions <- anti_join(
    observed_info,
    expected_ids,
    by = c("View", "ConditionID")
  )

  if (nrow(missing_conditions) > 0L || nrow(unexpected_conditions) > 0L) {
    stop(
      "The loaded pairwise files do not contain the expected single-omics ",
      "12-condition grid for both views.\nMissing conditions: ",
      ifelse(
        nrow(missing_conditions) == 0L,
        "none",
        paste(
          paste0(missing_conditions$View, ":", missing_conditions$ConditionID),
          collapse = "; "
        )
      ),
      "\nUnexpected conditions: ",
      ifelse(
        nrow(unexpected_conditions) == 0L,
        "none",
        paste(
          paste0(unexpected_conditions$View, ":", unexpected_conditions$ConditionID),
          collapse = "; "
        )
      )
    )
  }

  invisible(observed_info)
}


# ================================
# Locate and load pairwise overlap files
# ================================

dataset_map <- tibble(
  RequestedSubdir = dataset_subdirs,
  DatasetDir = file.path(summary_root_dir, dataset_subdirs),
  Found = dir.exists(file.path(summary_root_dir, dataset_subdirs))
)

write_csv(
  dataset_map,
  file.path(
    output_dir,
    paste0(endpoint_type, "_", top_n_label, "_single_omics_dataset_input_resolution.csv")
  ),
  na = ""
)

missing_dirs <- dataset_map %>%
  filter(!Found)

if (nrow(missing_dirs) > 0L) {
  warning(
    "The following ",
    endpoint_type,
    " dataset directories were not found:\n",
    paste(missing_dirs$RequestedSubdir, collapse = "\n")
  )
}

loaded_pairwise <- list()
load_status <- list()

for (i in seq_len(nrow(dataset_map))) {
  this_subdir <- dataset_map$RequestedSubdir[i]
  this_dir <- dataset_map$DatasetDir[i]
  this_found <- dataset_map$Found[i]

  dataset_label <- pretty_dataset_label(this_subdir)
  dataset_file_prefix <- clean_filename(this_subdir)

  if (!this_found) {
    load_status[[length(load_status) + 1L]] <- tibble(
      RequestedSubdir = this_subdir,
      DatasetLabel = dataset_label,
      Status = "MISSING_DATASET_DIR",
      PairwiseFile = NA_character_,
      Message = "Dataset directory was not found."
    )
    next
  }

  overlap_dir <- file.path(
    this_dir,
    single_omics_overlap_dir_name
  )

  if (!dir.exists(overlap_dir)) {
    load_status[[length(load_status) + 1L]] <- tibble(
      RequestedSubdir = this_subdir,
      DatasetLabel = dataset_label,
      Status = "MISSING_SINGLE_OMICS_OVERLAP_DIR",
      PairwiseFile = NA_character_,
      Message = paste0("Missing overlap directory: ", overlap_dir)
    )
    next
  }

  pairwise_file <- get_single_file(overlap_dir, pairwise_pattern)

  if (is.na(pairwise_file)) {
    load_status[[length(load_status) + 1L]] <- tibble(
      RequestedSubdir = this_subdir,
      DatasetLabel = dataset_label,
      Status = "MISSING_SINGLE_OMICS_PAIRWISE_FILE",
      PairwiseFile = NA_character_,
      Message = paste0("Could not find pairwise overlap file in: ", overlap_dir)
    )
    next
  }

  pairwise_dat <- read_csv(
    pairwise_file,
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    standardize_pairwise_overlap(pairwise_file) %>%
    filter(
      as.character(View_A) %in% expected_views,
      as.character(View_B) %in% expected_views
    ) %>%
    mutate(
      RequestedSubdir = this_subdir,
      DatasetLabel = dataset_label,
      DatasetFilePrefix = dataset_file_prefix,
      EndpointType = endpoint_type,
      PairwiseFile = pairwise_file,
      .before = 1
    )

  loaded_pairwise[[length(loaded_pairwise) + 1L]] <- pairwise_dat

  load_status[[length(load_status) + 1L]] <- tibble(
    RequestedSubdir = this_subdir,
    DatasetLabel = dataset_label,
    Status = "OK",
    PairwiseFile = pairwise_file,
    Message = NA_character_
  )
}

load_status_tbl <- bind_rows(load_status)

write_csv(
  load_status_tbl,
  file.path(
    output_dir,
    paste0(endpoint_type, "_", top_n_label, "_single_omics_load_status.csv")
  ),
  na = ""
)

if (length(loaded_pairwise) == 0L) {
  stop(
    "No ",
    endpoint_type,
    " single-omics pairwise overlap files were loaded. Check the load-status CSV."
  )
}

all_pairwise <- bind_rows(loaded_pairwise)

n_loaded_datasets <- all_pairwise %>%
  distinct(RequestedSubdir) %>%
  nrow()

min_valid_datasets_per_cell <- ceiling(
  min_valid_dataset_fraction * n_loaded_datasets
)

message("Loaded ", endpoint_type, " datasets: ", n_loaded_datasets)
message(
  "Minimum valid datasets per colored summary cell: ",
  min_valid_datasets_per_cell
)

write_csv(
  all_pairwise,
  file.path(
    output_dir,
    paste0(endpoint_type, "_", top_n_label, "_single_omics_all_pairwise_overlap_loaded.csv")
  ),
  na = ""
)


# ================================
# Define and validate the 12-condition grid for each view
# ================================

condition_info <- build_expected_condition_info()
validate_observed_condition_grid(all_pairwise, condition_info)

condition_counts_by_view <- condition_info %>%
  count(View, name = "n_conditions")

if (any(condition_counts_by_view$n_conditions != expected_conditions_per_view)) {
  stop("Internal error: expected exactly 12 conditions per single-omics view.")
}

write_csv(
  condition_info,
  file.path(
    output_dir,
    paste0(endpoint_type, "_", top_n_label, "_single_omics_condition_info.csv")
  ),
  na = ""
)

write_csv(
  condition_counts_by_view,
  file.path(
    output_dir,
    paste0(endpoint_type, "_", top_n_label, "_single_omics_condition_counts_by_view.csv")
  ),
  na = ""
)

print(condition_counts_by_view)


# ================================
# Create dataset-level condition-pair values
# ================================

dataset_condition_pair_values <- all_pairwise %>%
  filter(
    as.character(View_A) == as.character(View_B),
    as.character(View_A) %in% expected_views
  ) %>%
  mutate(
    View = as.character(View_A),

    CanonicalMetabTransform_A = canonical_metab_transform(
      View_A,
      MetabTransform_A
    ),
    CanonicalTaxaTransform_A = canonical_taxa_transform(
      View_A,
      TaxaTransform_A
    ),
    CanonicalMetabReduction_A = canonical_metab_reduction(
      View_A,
      MetabReduction_A
    ),
    CanonicalTaxaReduction_A = canonical_taxa_reduction(
      View_A,
      TaxaReduction_A
    ),

    CanonicalMetabTransform_B = canonical_metab_transform(
      View_B,
      MetabTransform_B
    ),
    CanonicalTaxaTransform_B = canonical_taxa_transform(
      View_B,
      TaxaTransform_B
    ),
    CanonicalMetabReduction_B = canonical_metab_reduction(
      View_B,
      MetabReduction_B
    ),
    CanonicalTaxaReduction_B = canonical_taxa_reduction(
      View_B,
      TaxaReduction_B
    ),

    ConditionID_A = make_condition_id(
      View_A,
      ModelParsed_A,
      CanonicalMetabTransform_A,
      CanonicalTaxaTransform_A,
      CanonicalMetabReduction_A,
      CanonicalTaxaReduction_A
    ),

    ConditionID_B = make_condition_id(
      View_B,
      ModelParsed_B,
      CanonicalMetabTransform_B,
      CanonicalTaxaTransform_B,
      CanonicalMetabReduction_B,
      CanonicalTaxaReduction_B
    ),

    is_self_condition = ConditionID_A == ConditionID_B,

    valid_cell_comparison =
      is_self_condition |
      (valid_comparison & !is.na(overlap_percent)),

    cell_overlap_percent = case_when(
      is_self_condition ~ heatmap_max_percent,
      valid_comparison & !is.na(overlap_percent) ~
        pmin(overlap_percent, heatmap_max_percent),
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(ConditionID_A), !is.na(ConditionID_B)) %>%
  group_by(
    RequestedSubdir,
    DatasetLabel,
    EndpointType,
    View,
    ConditionID_A,
    ConditionID_B
  ) %>%
  summarise(
    dataset_overlap_percent = safe_mean(
      cell_overlap_percent[valid_cell_comparison]
    ),
    valid_dataset_cell = any(valid_cell_comparison, na.rm = TRUE),
    n_source_rows_for_cell = n(),
    min_source_overlap_percent = safe_min(
      cell_overlap_percent[valid_cell_comparison]
    ),
    max_source_overlap_percent = safe_max(
      cell_overlap_percent[valid_cell_comparison]
    ),
    source_overlap_range = case_when(
      is.na(min_source_overlap_percent) | is.na(max_source_overlap_percent) ~
        NA_real_,
      TRUE ~ max_source_overlap_percent - min_source_overlap_percent
    ),
    .groups = "drop"
  )

write_csv(
  dataset_condition_pair_values,
  file.path(
    output_dir,
    paste0(endpoint_type, "_", top_n_label, "_single_omics_dataset_condition_pair_values.csv")
  ),
  na = ""
)


# ================================
# Aggregate condition-pair values across datasets
# ================================

condition_grid <- condition_info %>%
  select(
    View,
    ConditionID_A = ConditionID
  ) %>%
  inner_join(
    condition_info %>%
      select(
        View,
        ConditionID_B = ConditionID
      ),
    by = "View",
    relationship = "many-to-many"
  )

condition_grid_counts <- condition_grid %>%
  count(View, name = "n_condition_pairs")

write_csv(
  condition_grid_counts,
  file.path(
    output_dir,
    paste0(endpoint_type, "_", top_n_label, "_single_omics_condition_grid_counts_by_view.csv")
  ),
  na = ""
)

aggregate_condition_pair_summary <- condition_grid %>%
  left_join(
    dataset_condition_pair_values,
    by = c("View", "ConditionID_A", "ConditionID_B")
  ) %>%
  group_by(
    View,
    ConditionID_A,
    ConditionID_B
  ) %>%
  summarise(
    n_possible_datasets = n_loaded_datasets,
    n_valid_datasets = sum(!is.na(dataset_overlap_percent), na.rm = TRUE),
    mean_topn_overlap_across_datasets = safe_mean(dataset_overlap_percent),
    median_topn_overlap_across_datasets = safe_median(dataset_overlap_percent),
    sd_topn_overlap_across_datasets = safe_sd(dataset_overlap_percent),
    iqr_topn_overlap_across_datasets = safe_iqr(dataset_overlap_percent),
    .groups = "drop"
  ) %>%
  mutate(
    is_self_condition = ConditionID_A == ConditionID_B,

    valid_summary_cell =
      is_self_condition |
      n_valid_datasets >= min_valid_datasets_per_cell,

    mean_plot_value = case_when(
      is_self_condition ~ heatmap_max_percent,
      valid_summary_cell ~
        pmin(mean_topn_overlap_across_datasets, heatmap_max_percent),
      TRUE ~ NA_real_
    ),

    median_plot_value = case_when(
      is_self_condition ~ heatmap_max_percent,
      valid_summary_cell ~
        pmin(median_topn_overlap_across_datasets, heatmap_max_percent),
      TRUE ~ NA_real_
    )
  ) %>%
  left_join(
    condition_info %>%
      select(
        View,
        ConditionID_A = ConditionID,
        ConditionLabel_A = ConditionLabel,
        ConditionOrder_A = ConditionOrder,
        ModelParsed_A = ModelParsed,
        MetabTransform_A = MetabTransform,
        TaxaTransform_A = TaxaTransform,
        MetabReduction_A = MetabReduction,
        TaxaReduction_A = TaxaReduction
      ),
    by = c("View", "ConditionID_A")
  ) %>%
  left_join(
    condition_info %>%
      select(
        View,
        ConditionID_B = ConditionID,
        ConditionLabel_B = ConditionLabel,
        ConditionOrder_B = ConditionOrder,
        ModelParsed_B = ModelParsed,
        MetabTransform_B = MetabTransform,
        TaxaTransform_B = TaxaTransform,
        MetabReduction_B = MetabReduction,
        TaxaReduction_B = TaxaReduction
      ),
    by = c("View", "ConditionID_B")
  ) %>%
  arrange(
    match(View, expected_views),
    ConditionOrder_A,
    ConditionOrder_B
  )

write_csv(
  aggregate_condition_pair_summary,
  file.path(
    output_dir,
    paste0(
      endpoint_type,
      "_",
      top_n_label,
      "_single_omics_aggregate_condition_pair_summary_across_datasets.csv"
    )
  ),
  na = ""
)

invalid_summary_cells <- aggregate_condition_pair_summary %>%
  filter(!valid_summary_cell) %>%
  arrange(
    match(View, expected_views),
    ConditionOrder_A,
    ConditionOrder_B
  )

write_csv(
  invalid_summary_cells,
  file.path(
    output_dir,
    paste0(
      endpoint_type,
      "_",
      top_n_label,
      "_single_omics_invalid_summary_cells.csv"
    )
  ),
  na = ""
)


# ================================
# Plot aggregate 12 x 12 heatmaps
# ================================

plot_aggregate_condition_heatmap <- function(
    aggregate_tbl,
    view_name,
    metric = c("median", "mean"),
    show_numbers = TRUE
) {
  metric <- match.arg(metric)

  metric_col <- case_when(
    metric == "median" ~ "median_plot_value",
    metric == "mean" ~ "mean_plot_value"
  )

  metric_title <- case_when(
    metric == "median" ~ "Median",
    metric == "mean" ~ "Mean"
  )

  view_high_color <- get_view_heatmap_color(view_name)
  view_display_label <- get_view_display_label(view_name)

  condition_order <- condition_info %>%
    filter(View == view_name) %>%
    arrange(ConditionOrder) %>%
    pull(ConditionID)

  label_lookup <- condition_info %>%
    filter(View == view_name) %>%
    arrange(ConditionOrder) %>%
    select(ConditionID, ConditionLabel) %>%
    deframe()

  plot_tbl <- aggregate_tbl %>%
    filter(View == view_name) %>%
    mutate(
      keep_triangle = case_when(
        triangle_to_show == "lower_left" ~ ConditionOrder_A > ConditionOrder_B,
        triangle_to_show == "upper_right" ~ ConditionOrder_A < ConditionOrder_B,
        TRUE ~ TRUE
      ),

      keep_diagonal = show_diagonal & ConditionOrder_A == ConditionOrder_B,
      keep_cell = keep_triangle | keep_diagonal
    ) %>%
    filter(keep_cell) %>%
    mutate(
      ConditionID_A = factor(ConditionID_A, levels = rev(condition_order)),
      ConditionID_B = factor(ConditionID_B, levels = condition_order),

      plot_value = case_when(
        ConditionOrder_A == ConditionOrder_B ~ heatmap_max_percent,
        TRUE ~ pmin(.data[[metric_col]], heatmap_max_percent)
      ),

      overlap_label = case_when(
        is.na(plot_value) ~ "",
        TRUE ~ sprintf("%.0f", plot_value)
      ),

      overlap_text_color = case_when(
        is.na(plot_value) ~ "black",
        plot_value >= 55 ~ "white",
        TRUE ~ "black"
      )
    )

  p <- ggplot(
    plot_tbl,
    aes(
      x = ConditionID_B,
      y = ConditionID_A,
      fill = plot_value
    )
  ) +
    geom_tile(
      color = "white",
      linewidth = 0.25
    ) +
    scale_fill_gradient(
      name = paste0(
        metric_title,
        " ",
        top_n_title,
        "\noverlap (%)"
      ),
      low = "white",
      high = view_high_color,
      limits = c(0, heatmap_max_percent),
      breaks = seq(0, heatmap_max_percent, by = 25),
      na.value = "grey80"
    ) +
    coord_fixed() +
    scale_x_discrete(
      labels = label_lookup,
      drop = FALSE
    ) +
    scale_y_discrete(
      labels = label_lookup,
      drop = FALSE
    ) +
    labs(
      title = paste0(
        endpoint_title,
        " endpoints: ",
        metric_title,
        " ",
        top_n_lower,
        " feature-overlap summary, ",
        view_display_label
      ),
      subtitle = paste0(
        "Conditions are model + relevant preprocessing pair only. ",
        "The ",
        str_replace_all(triangle_to_show, "_", "-"),
        " triangle and the ",
        heatmap_max_percent,
        "% diagonal are shown."
      ),
      x = "Model and preprocessing pair",
      y = "Model and preprocessing pair",
      caption = paste0(
        "Gray cells indicate fewer than ",
        min_valid_datasets_per_cell,
        " valid ",
        endpoint_type,
        " datasets for that condition pair. Blank cells are mirrored ",
        "comparisons that were intentionally removed."
      )
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 8
      ),
      axis.text.y = element_text(size = 8),
      legend.position = "right"
    )

  if (show_numbers) {
    p <- p +
      geom_text(
        aes(
          label = overlap_label,
          color = overlap_text_color
        ),
        size = 2.7
      ) +
      scale_color_identity()
  }

  p
}

view_names <- expected_views
median_heatmaps <- list()
mean_heatmaps <- list()

for (this_view in view_names) {
  message(
    "Plotting aggregate median and mean heatmaps for ",
    get_view_display_label(this_view),
    "."
  )

  median_heatmaps[[this_view]] <- plot_aggregate_condition_heatmap(
    aggregate_tbl = aggregate_condition_pair_summary,
    view_name = this_view,
    metric = "median",
    show_numbers = show_numbers_in_heatmaps
  )

  mean_heatmaps[[this_view]] <- plot_aggregate_condition_heatmap(
    aggregate_tbl = aggregate_condition_pair_summary,
    view_name = this_view,
    metric = "mean",
    show_numbers = show_numbers_in_heatmaps
  )

  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "figure_",
        endpoint_type,
        "_",
        top_n_label,
        "_MEDIAN_aggregate_single_omics_12x12_VIEW_",
        get_view_output_name(this_view),
        ".pdf"
      )
    ),
    plot = median_heatmaps[[this_view]],
    width = heatmap_width,
    height = heatmap_height
  )

  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "figure_",
        endpoint_type,
        "_",
        top_n_label,
        "_MEAN_aggregate_single_omics_12x12_VIEW_",
        get_view_output_name(this_view),
        ".pdf"
      )
    ),
    plot = mean_heatmaps[[this_view]],
    width = heatmap_width,
    height = heatmap_height
  )
}


# ================================
# Multi-page PDFs
# ================================

median_multipage_file <- file.path(
  output_dir,
  paste0(
    "figure_",
    endpoint_type,
    "_",
    top_n_label,
    "_MEDIAN_aggregate_single_omics_12x12_all_views_multipage.pdf"
  )
)

pdf(
  file = median_multipage_file,
  width = heatmap_width,
  height = heatmap_height
)

for (this_view in view_names) {
  print(median_heatmaps[[this_view]])
}

dev.off()

mean_multipage_file <- file.path(
  output_dir,
  paste0(
    "figure_",
    endpoint_type,
    "_",
    top_n_label,
    "_MEAN_aggregate_single_omics_12x12_all_views_multipage.pdf"
  )
)

pdf(
  file = mean_multipage_file,
  width = heatmap_width,
  height = heatmap_height
)

for (this_view in view_names) {
  print(mean_heatmaps[[this_view]])
}

dev.off()


# ================================
# Compact diagnostics
# ================================

summary_cell_diagnostics <- aggregate_condition_pair_summary %>%
  group_by(View) %>%
  summarise(
    n_cells_total = n(),
    n_self_cells = sum(is_self_condition, na.rm = TRUE),
    n_valid_summary_cells = sum(valid_summary_cell, na.rm = TRUE),
    n_invalid_nonself_cells =
      sum(!valid_summary_cell & !is_self_condition, na.rm = TRUE),
    median_n_valid_datasets =
      safe_median(n_valid_datasets[!is_self_condition]),
    min_n_valid_datasets =
      safe_min(n_valid_datasets[!is_self_condition]),
    max_n_valid_datasets =
      safe_max(n_valid_datasets[!is_self_condition]),
    .groups = "drop"
  ) %>%
  mutate(
    ViewDisplayLabel = map_chr(View, get_view_display_label),
    .after = View
  )

write_csv(
  summary_cell_diagnostics,
  file.path(
    output_dir,
    paste0(
      endpoint_type,
      "_",
      top_n_label,
      "_single_omics_summary_cell_diagnostics_by_view.csv"
    )
  ),
  na = ""
)

print(summary_cell_diagnostics)


# ================================
# Final messages
# ================================

message("\nDone.")
message(
  endpoint_title,
  " single-omics aggregate 12 x 12 heatmap output directory: ",
  output_dir
)
message("Main files created:")
message(
  "  - figure_",
  endpoint_type,
  "_",
  top_n_label,
  "_MEDIAN_aggregate_single_omics_12x12_all_views_multipage.pdf"
)
message(
  "  - figure_",
  endpoint_type,
  "_",
  top_n_label,
  "_MEAN_aggregate_single_omics_12x12_all_views_multipage.pdf"
)
message("  - Individual metabolomics and taxa median/mean PDF heatmaps")
message(
  "  - ",
  endpoint_type,
  "_",
  top_n_label,
  "_single_omics_aggregate_condition_pair_summary_across_datasets.csv"
)
message(
  "  - ",
  endpoint_type,
  "_",
  top_n_label,
  "_single_omics_dataset_condition_pair_values.csv"
)

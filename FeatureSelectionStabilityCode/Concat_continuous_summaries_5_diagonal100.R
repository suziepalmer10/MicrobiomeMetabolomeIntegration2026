library(tidyverse)

# ================================
# User settings
# ================================

summary_root_dir <- "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/FeatureImportanceLogs/Validation_Feature_Importance/_feature_importance_seed_summaries_with_stable_scores"

top_n_features <- 20L

top_n_label <- paste0("top", top_n_features)
top_n_title <- paste0("Top-", top_n_features)
top_n_lower <- paste0("top-", top_n_features)

triangle_to_show <- "lower_left"

# Display self-comparisons as diagonal squares at the maximum overlap value.
show_diagonal <- TRUE
heatmap_max_percent <- 100

continuous_subdirs <- c(
  "Era_Cholesterol",
  "Era_Glucose",
  "Franzosa_CD_Fp",
  "Franzosa_UC_Fp",
  "Franzosa_IBD_Fp",
  "Wang_creatinine",
  "Wang_eGFR",
  "Wang_urea"
)

# View-specific heatmap colors
view_heatmap_colors <- c(
  "metabolomics" = "#9D1B28",
  "mss" = "#0152A1",
  "concat_none" = "#BDBADC",
  "concat_scaling" = "#E4A0F7",
  "concat_scaling_weights" = "#551053"
)

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

continuous_output_dir <- file.path(
  summary_root_dir,
  paste0("_continuous_", top_n_label, "_aggregate_24x24_condition_heatmaps")
)

dir.create(continuous_output_dir, recursive = TRUE, showWarnings = FALSE)

# Main choice:
# 1.0 means a heatmap cell is colored only if all loaded continuous datasets have a valid value.
# 0.8 means at least 80% of loaded continuous datasets must have a valid value.
# 0.5 means at least half of loaded continuous datasets must have a valid value.
min_valid_dataset_fraction <- 1.0

# Diagonal self-comparisons are displayed at 100%, rather than gray or blank.

show_numbers_in_heatmaps <- TRUE

preferred_view_order <- c(
  "metabolomics",
  "mss",
  "concat_none",
  "concat_scaling",
  "concat_scaling_weights"
)

preferred_model_order <- c("enet", "rf", "xgb")
preferred_metab_transform_order <- c("log2", "none")
preferred_taxa_transform_order <- c("clr", "none")
preferred_metab_reduction_order <- c("limma", "none")
preferred_taxa_reduction_order <- c("wilcox", "none")

pairwise_pattern <- paste0("_pairwise_", top_n_label, "_overlap\\.csv$")

heatmap_width <- 13
heatmap_height <- 11


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

safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_median <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  median(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  if (sum(!is.na(x)) < 2) return(NA_real_)
  sd(x, na.rm = TRUE)
}

safe_iqr <- function(x) {
  if (sum(!is.na(x)) < 2) return(NA_real_)
  IQR(x, na.rm = TRUE)
}

safe_min <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

validate_required_columns <- function(dat, required_cols, file_label = "data frame") {
  missing_cols <- setdiff(required_cols, names(dat))
  
  if (length(missing_cols) > 0) {
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
  
  if (length(files) == 0) {
    return(NA_character_)
  }
  
  if (length(files) > 1) {
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
    paste0(top_n_label, "_overlap_percent") %in% names(dat) ~ paste0(top_n_label, "_overlap_percent"),
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
      valid_comparison = as.character(valid_comparison) %in% c("TRUE", "true", "1")
    )
}

make_condition_id <- function(model, metab_t, taxa_t, metab_r, taxa_r) {
  paste0(
    "model-", model,
    "__metabT-", metab_t,
    "__taxaT-", taxa_t,
    "__metabR-", metab_r,
    "__taxaR-", taxa_r
  )
}

make_condition_label <- function(model, metab_t, taxa_t, metab_r, taxa_r) {
  paste0(
    model,
    "\n",
    "mT=", metab_t,
    "; tT=", taxa_t,
    "\n",
    "mR=", metab_r,
    "; tR=", taxa_r
  )
}

safe_factor <- function(x, preferred_levels) {
  observed <- unique(as.character(x))
  extra <- setdiff(sort(observed), preferred_levels)
  
  factor(
    as.character(x),
    levels = c(preferred_levels, extra)
  )
}


# ================================
# Locate and load continuous pairwise overlap files
# ================================

continuous_dataset_map <- tibble(
  RequestedSubdir = continuous_subdirs,
  DatasetDir = file.path(summary_root_dir, continuous_subdirs),
  Found = dir.exists(file.path(summary_root_dir, continuous_subdirs))
)

write_csv(
  continuous_dataset_map,
  file.path(continuous_output_dir, paste0("continuous_", top_n_label, "_dataset_input_resolution.csv")),
  na = ""
)

missing_dirs <- continuous_dataset_map %>%
  filter(!Found)

if (nrow(missing_dirs) > 0) {
  warning(
    "The following continuous dataset directories were not found:\n",
    paste(missing_dirs$RequestedSubdir, collapse = "\n")
  )
}

loaded_pairwise <- list()
load_status <- list()

for (i in seq_len(nrow(continuous_dataset_map))) {
  this_subdir <- continuous_dataset_map$RequestedSubdir[i]
  this_dir <- continuous_dataset_map$DatasetDir[i]
  this_found <- continuous_dataset_map$Found[i]
  
  dataset_label <- pretty_dataset_label(this_subdir)
  dataset_file_prefix <- clean_filename(this_subdir)
  
  if (!this_found) {
    load_status[[length(load_status) + 1]] <- tibble(
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
    paste0("_", top_n_label, "_stable_score_seed_overlap")
  )
  
  if (!dir.exists(overlap_dir)) {
    load_status[[length(load_status) + 1]] <- tibble(
      RequestedSubdir = this_subdir,
      DatasetLabel = dataset_label,
      Status = "MISSING_OVERLAP_DIR",
      PairwiseFile = NA_character_,
      Message = paste0("Missing overlap directory: ", overlap_dir)
    )
    next
  }
  
  pairwise_file <- get_single_file(overlap_dir, pairwise_pattern)
  
  if (is.na(pairwise_file)) {
    load_status[[length(load_status) + 1]] <- tibble(
      RequestedSubdir = this_subdir,
      DatasetLabel = dataset_label,
      Status = "MISSING_PAIRWISE_FILE",
      PairwiseFile = NA_character_,
      Message = paste0("Could not find pairwise overlap file in: ", overlap_dir)
    )
    next
  }
  
  pairwise_dat <- read_csv(pairwise_file, show_col_types = FALSE, progress = FALSE) %>%
    standardize_pairwise_overlap(pairwise_file) %>%
    mutate(
      RequestedSubdir = this_subdir,
      DatasetLabel = dataset_label,
      DatasetFilePrefix = dataset_file_prefix,
      EndpointType = "continuous",
      PairwiseFile = pairwise_file,
      .before = 1
    )
  
  loaded_pairwise[[length(loaded_pairwise) + 1]] <- pairwise_dat
  
  load_status[[length(load_status) + 1]] <- tibble(
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
  file.path(continuous_output_dir, paste0("continuous_", top_n_label, "_load_status.csv")),
  na = ""
)

if (length(loaded_pairwise) == 0) {
  stop("No continuous pairwise overlap files were loaded. Check the load-status CSV.")
}

all_pairwise <- bind_rows(loaded_pairwise)

n_loaded_continuous_datasets <- all_pairwise %>%
  distinct(RequestedSubdir) %>%
  nrow()

min_valid_datasets_per_cell <- ceiling(
  min_valid_dataset_fraction * n_loaded_continuous_datasets
)

message("Loaded continuous datasets: ", n_loaded_continuous_datasets)
message("Minimum valid datasets per colored summary cell: ", min_valid_datasets_per_cell)

write_csv(
  all_pairwise,
  file.path(continuous_output_dir, paste0("continuous_", top_n_label, "_all_pairwise_overlap_loaded.csv")),
  na = ""
)


# ================================
# Create condition metadata
# ================================

condition_info_A <- all_pairwise %>%
  transmute(
    View = as.character(View_A),
    ConditionID = make_condition_id(
      ModelParsed_A,
      MetabTransform_A,
      TaxaTransform_A,
      MetabReduction_A,
      TaxaReduction_A
    ),
    ConditionLabel = make_condition_label(
      ModelParsed_A,
      MetabTransform_A,
      TaxaTransform_A,
      MetabReduction_A,
      TaxaReduction_A
    ),
    ModelParsed = as.character(ModelParsed_A),
    MetabTransform = as.character(MetabTransform_A),
    TaxaTransform = as.character(TaxaTransform_A),
    MetabReduction = as.character(MetabReduction_A),
    TaxaReduction = as.character(TaxaReduction_A)
  )

condition_info_B <- all_pairwise %>%
  transmute(
    View = as.character(View_B),
    ConditionID = make_condition_id(
      ModelParsed_B,
      MetabTransform_B,
      TaxaTransform_B,
      MetabReduction_B,
      TaxaReduction_B
    ),
    ConditionLabel = make_condition_label(
      ModelParsed_B,
      MetabTransform_B,
      TaxaTransform_B,
      MetabReduction_B,
      TaxaReduction_B
    ),
    ModelParsed = as.character(ModelParsed_B),
    MetabTransform = as.character(MetabTransform_B),
    TaxaTransform = as.character(TaxaTransform_B),
    MetabReduction = as.character(MetabReduction_B),
    TaxaReduction = as.character(TaxaReduction_B)
  )

condition_info <- bind_rows(condition_info_A, condition_info_B) %>%
  distinct() %>%
  mutate(
    View = safe_factor(View, preferred_view_order),
    ModelParsed = safe_factor(ModelParsed, preferred_model_order),
    MetabTransform = safe_factor(MetabTransform, preferred_metab_transform_order),
    TaxaTransform = safe_factor(TaxaTransform, preferred_taxa_transform_order),
    MetabReduction = safe_factor(MetabReduction, preferred_metab_reduction_order),
    TaxaReduction = safe_factor(TaxaReduction, preferred_taxa_reduction_order)
  ) %>%
  arrange(
    View,
    ModelParsed,
    MetabTransform,
    TaxaTransform,
    MetabReduction,
    TaxaReduction,
    ConditionID
  ) %>%
  group_by(View) %>%
  mutate(
    ConditionOrder = row_number()
  ) %>%
  ungroup()

write_csv(
  condition_info,
  file.path(continuous_output_dir, paste0("continuous_", top_n_label, "_condition_info.csv")),
  na = ""
)

condition_counts_by_view <- condition_info %>%
  count(View, name = "n_conditions")

write_csv(
  condition_counts_by_view,
  file.path(continuous_output_dir, paste0("continuous_", top_n_label, "_condition_counts_by_view.csv")),
  na = ""
)

print(condition_counts_by_view)


# ================================
# Create dataset-level condition-pair values
# ================================

dataset_condition_pair_values <- all_pairwise %>%
  filter(
    as.character(View_A) == as.character(View_B)
  ) %>%
  mutate(
    View = as.character(View_A),
    
    ConditionID_A = make_condition_id(
      ModelParsed_A,
      MetabTransform_A,
      TaxaTransform_A,
      MetabReduction_A,
      TaxaReduction_A
    ),
    
    ConditionID_B = make_condition_id(
      ModelParsed_B,
      MetabTransform_B,
      TaxaTransform_B,
      MetabReduction_B,
      TaxaReduction_B
    ),
    
    is_self_condition = ConditionID_A == ConditionID_B,
    
    # Self-comparisons are defined as complete overlap. Non-self values are
    # retained only when the underlying pairwise comparison is valid.
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
    n_raw_rows_for_cell = n(),
    .groups = "drop"
  )

write_csv(
  dataset_condition_pair_values,
  file.path(continuous_output_dir, paste0("continuous_", top_n_label, "_dataset_condition_pair_values.csv")),
  na = ""
)


# ================================
# Aggregate condition-pair values across continuous datasets
# ================================

condition_info_unique <- condition_info %>%
  distinct(
    View,
    ConditionID,
    .keep_all = TRUE
  )

condition_grid <- condition_info_unique %>%
  select(
    View,
    ConditionID_A = ConditionID
  ) %>%
  inner_join(
    condition_info_unique %>%
      select(
        View,
        ConditionID_B = ConditionID
      ),
    by = "View",
    relationship = "many-to-many"
  )

condition_grid_counts <- condition_grid %>%
  count(View, name = "n_condition_pairs")

print(condition_grid_counts)

write_csv(
  condition_grid_counts,
  file.path(continuous_output_dir, paste0("continuous_", top_n_label, "_condition_grid_counts_by_view.csv")),
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
    n_possible_continuous_datasets = n_loaded_continuous_datasets,
    
    n_valid_continuous_datasets =
      sum(!is.na(dataset_overlap_percent), na.rm = TRUE),
    
    mean_topn_overlap_across_datasets =
      safe_mean(dataset_overlap_percent),
    
    median_topn_overlap_across_datasets =
      safe_median(dataset_overlap_percent),
    
    sd_topn_overlap_across_datasets =
      safe_sd(dataset_overlap_percent),
    
    iqr_topn_overlap_across_datasets =
      safe_iqr(dataset_overlap_percent),
    
    .groups = "drop"
  ) %>%
  mutate(
    is_self_condition = ConditionID_A == ConditionID_B,
    
    valid_summary_cell =
      is_self_condition |
      n_valid_continuous_datasets >= min_valid_datasets_per_cell,
    
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
    View,
    ConditionOrder_A,
    ConditionOrder_B
  )

write_csv(
  aggregate_condition_pair_summary,
  file.path(
    continuous_output_dir,
    paste0("continuous_", top_n_label, "_aggregate_condition_pair_summary_across_datasets.csv")
  ),
  na = ""
)

invalid_summary_cells <- aggregate_condition_pair_summary %>%
  filter(!valid_summary_cell) %>%
  arrange(View, ConditionOrder_A, ConditionOrder_B)

write_csv(
  invalid_summary_cells,
  file.path(
    continuous_output_dir,
    paste0("continuous_", top_n_label, "_invalid_or_self_summary_cells.csv")
  ),
  na = ""
)


# ================================
# Plot aggregate 24 x 24 heatmaps
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
  
  # View-specific high-gradient color
  view_high_color <- get_view_heatmap_color(view_name)
  
  condition_order <- condition_info %>%
    filter(as.character(View) == view_name) %>%
    arrange(ConditionOrder) %>%
    pull(ConditionID)
  
  label_lookup <- condition_info %>%
    filter(as.character(View) == view_name) %>%
    arrange(ConditionOrder) %>%
    select(ConditionID, ConditionLabel) %>%
    deframe()
  
  plot_tbl <- aggregate_tbl %>%
    filter(as.character(View) == view_name) %>%
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
      
      # White text improves readability for darker high-value tiles.
      overlap_text_color = case_when(
        is.na(plot_value) ~ "black",
        view_name %in% c("metabolomics", "mss", "taxa", "concat_scaling_weights") &
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
      linewidth = 0.18
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
        "Continuous endpoints: ",
        metric_title,
        " ",
        top_n_lower,
        " feature-overlap summary, View = ",
        view_name
      ),
      subtitle = paste0(
        "Only the ",
        str_replace_all(triangle_to_show, "_", "-"),
        " triangle and the ",
        heatmap_max_percent,
        "% diagonal are shown because the condition-pair matrix is symmetric."
      ),
      x = "Pipeline condition",
      y = "Pipeline condition",
      caption = paste0(
        "Gray cells indicate fewer than ",
        min_valid_datasets_per_cell,
        " valid continuous datasets for that condition pair. ",
        "Blank cells are mirrored comparisons that were intentionally removed."
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
        size = 7
      ),
      axis.text.y = element_text(size = 7),
      legend.position = "right"
    )
  
  if (show_numbers) {
    p <- p +
      geom_text(
        aes(
          label = overlap_label,
          color = overlap_text_color
        ),
        size = 2.2
      ) +
      scale_color_identity()
  }
  
  p
}


view_names <- condition_info %>%
  mutate(View = as.character(View)) %>%
  distinct(View) %>%
  mutate(
    View = safe_factor(View, preferred_view_order)
  ) %>%
  arrange(View) %>%
  pull(View) %>%
  as.character()

median_heatmaps <- list()
mean_heatmaps <- list()

for (this_view in view_names) {
  message("Plotting aggregate median and mean heatmaps for View = ", this_view)
  
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
      continuous_output_dir,
      paste0(
        "figure_continuous_",
        top_n_label,
        "_MEDIAN_aggregate_24x24_VIEW_",
        clean_filename(this_view),
        ".pdf"
      )
    ),
    plot = median_heatmaps[[this_view]],
    width = heatmap_width,
    height = heatmap_height
  )
  
  ggsave(
    filename = file.path(
      continuous_output_dir,
      paste0(
        "figure_continuous_",
        top_n_label,
        "_MEAN_aggregate_24x24_VIEW_",
        clean_filename(this_view),
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
  continuous_output_dir,
  paste0("figure_continuous_", top_n_label, "_MEDIAN_aggregate_24x24_all_views_multipage.pdf")
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
  continuous_output_dir,
  paste0("figure_continuous_", top_n_label, "_MEAN_aggregate_24x24_all_views_multipage.pdf")
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
    
    median_n_valid_continuous_datasets =
      safe_median(n_valid_continuous_datasets[!is_self_condition]),
    
    min_n_valid_continuous_datasets =
      safe_min(n_valid_continuous_datasets[!is_self_condition]),
    
    max_n_valid_continuous_datasets =
      safe_max(n_valid_continuous_datasets[!is_self_condition]),
    
    .groups = "drop"
  )

write_csv(
  summary_cell_diagnostics,
  file.path(
    continuous_output_dir,
    paste0("continuous_", top_n_label, "_summary_cell_diagnostics_by_view.csv")
  ),
  na = ""
)

print(summary_cell_diagnostics)


# ================================
# Final messages
# ================================

message("\nDone.")
message("Continuous aggregate 24 x 24 heatmap output directory: ", continuous_output_dir)
message("Main files created:")
message("  - figure_continuous_", top_n_label, "_MEDIAN_aggregate_24x24_all_views_multipage.pdf")
message("  - figure_continuous_", top_n_label, "_MEAN_aggregate_24x24_all_views_multipage.pdf")
message("  - Individual median and mean PDF heatmaps for each View")
message("  - continuous_", top_n_label, "_aggregate_condition_pair_summary_across_datasets.csv")
message("  - continuous_", top_n_label, "_dataset_condition_pair_values.csv")
message("  - continuous_", top_n_label, "_summary_cell_diagnostics_by_view.csv")
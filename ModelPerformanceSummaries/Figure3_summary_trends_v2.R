# ============================================================
# COMPETITIVE-SET SUMMARY USING WITHIN 2% OF BEST
# WITH EXACT-TIE-ADJUSTED RAW COUNTS
# ------------------------------------------------------------
# Main idea:
#   For each dataset and metric, keep all TEST configurations
#   whose mean_across_seeds is within 2% of the best value.
#
#   For metrics where higher is better (e.g., AUROC, AUPRC):
#
#     best_value = max(mean_across_seeds)
#     keep if mean_across_seeds >= best_value * 0.98
#
# Counting modes:
#
#   raw:
#     each qualifying configuration counts 1,
#     BUT exact ties at the same performance value
#     share credit 1 / n_ties
#
#   dataset_equal:
#     each dataset contributes total credit 1
#     within each metric, split evenly across
#     all qualifying configurations
#
# Notes:
#   - AUROC summarizes overall discriminative ability
#   - AUPRC emphasizes recovery of the positive class
#     and is more sensitive to class imbalance
# ============================================================


library(tidyverse)
library(stringr)
library(forcats)
library(patchwork)
library(scales)

# ------------------------------------------------------------
# 1) Input files
# ------------------------------------------------------------

setwd('/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/Performance_CSVs/Binary')
file_map <- c(
  "Yachida_control_vs_stage3_4" = "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/Performance_CSVs/Binary/Yachida_control_vs_stage3_4/results/tables/Yachida_control_vs_stage3_4_heatmap_summary.csv",
  "Yachida_control_vs_stage1_2" = "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/Performance_CSVs/Binary/Yachida_control_vs_stage1_2/results/tables/Yachida_control_vs_stage1_2_heatmap_summary.csv",
  "Yachida_control_vs_MP_stage0" = "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/Performance_CSVs/Binary/Yachida_control_vs_MP_stage0/results/tables/Yachida_control_vs_MP_stage0_heatmap_summary.csv",
  "Yachida_control_vs_CRC" = "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/Performance_CSVs/Binary/Yachida_control_vs_CRC/results/tables/Yachida_control_vs_CRC_heatmap_summary.csv",
  "Wang_convskidneyfailure" = "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/Performance_CSVs/Binary/Wang_convskidneyfailure/results/tables/Wang_convskidneyfailure_heatmap_summary.csv",
  "Franz_con_vs_UC" = "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/Performance_CSVs/Binary/Franzosa_con_vs_UC/results/tables/Franzosa_con_vs_UC_heatmap_summary.csv",
  "Franz_con_vs_IBD" = "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/Performance_CSVs/Binary/Franzosa_con_vs_IBD/results/tables/Franzosa_con_vs_IBD_heatmap_summary.csv",
  "Franz_con_vs_CD" = "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/Performance_CSVs/Binary/Franzosa_con_vs_CD/results/tables/Franzosa_con_vs_CD_heatmap_summary.csv",
  "Erawijantari_convsgas" = "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/Performance_CSVs/Binary/Erawijantari_convsgas/results/tables/Erawijantari_convsgas_heatmap_summary.csv"
)

dataset_pretty <- c(
  "Yachida_control_vs_stage3_4" = "CRC_Stage3_4",
  "Yachida_control_vs_stage1_2" = "CRC_Stage1_2",
  "Yachida_control_vs_MP_stage0" = "CRC_MP_stage0",
  "Yachida_control_vs_CRC" = "CRC",
  "Wang_convskidneyfailure" = 'kidneyfailure',
  "Franz_con_vs_UC" = "UC",
  "Franz_con_vs_IBD" = "IBD",
  "Franz_con_vs_CD" = "CD",
  "Erawijantari_convsgas" = 'gastrectomy'
)

# ------------------------------------------------------------
# 2) Helpers
# ------------------------------------------------------------
extract_field <- function(x, key) {
  out <- str_match(x, paste0(key, "=([^|]+)"))[, 2]
  str_trim(out)
}

recode_learner <- function(x) {
  recode(
    x,
    "rf" = "Random Forest",
    "enet" = "Elastic Net",
    "xgb" = "XGBoost",
    "xgboost" = "XGBoost",
    .default = x
  )
}

recode_category_display <- function(x) {
  recode(
    x,
    "concat_none" = "Concatenation: none",
    "concat_scaling" = "Concatenation: scaling",
    "concat_scaling_weights" = "Concatenation: weighted",
    "stacked_average" = "Stacking: average",
    "stacked_averaging" = "Stacking: average",
    "stacked_pls" = "Stacking: PLS",
    "stacked_nnls" = "Stacking: NNLS",
    "stacked_sparse" = "Stacking: sparse",
    "mss" = "MSS only",
    "metabolomics" = "Metabolomics only",
    .default = x
  )
}

integration_family_fn <- function(category) {
  case_when(
    str_detect(category, "^concat")  ~ "Concatenation",
    str_detect(category, "^stacked") ~ "Stacking",
    TRUE ~ "Single-omic"
  )
}

dimensionality_label <- function(metab_r, taxa_r) {
  case_when(
    metab_r == "limma" & taxa_r == "wilcox" ~ "Feature reduced",
    metab_r == "none"  & taxa_r == "none"   ~ "Full dimensional",
    TRUE ~ "Other"
  )
}

metab_label <- function(x) {
  case_when(
    x == "log2" ~ "Metab: log2",
    x == "none" ~ "Metab: none",
    TRUE ~ "Metab: other"
  )
}

taxa_label <- function(x) {
  case_when(
    x == "clr" ~ "MSS: CLR",
    x == "none" ~ "MSS: none",
    TRUE ~ "MSS: other"
  )
}

# ------------------------------------------------------------
# 3) Read and standardize one file
# ------------------------------------------------------------
read_perf_file <- function(path, dataset_name) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      dataset = dataset_name,
      dataset_label = dataset_pretty[[dataset_name]],
      learner_raw = extract_field(pipeline_id, "model"),
      metab_scaling_raw = extract_field(pipeline_id, "metabT"),
      taxa_scaling_raw = extract_field(pipeline_id, "taxaT"),
      metab_reduction = extract_field(pipeline_id, "metabR"),
      taxa_reduction = extract_field(pipeline_id, "taxaR")
    ) %>%
    mutate(
      learner = recode_learner(learner_raw),
      category_display = recode_category_display(category),
      integration_family = integration_family_fn(category),
      dimensionality = dimensionality_label(metab_reduction, taxa_reduction),
      metab_scaling = metab_label(metab_scaling_raw),
      taxa_scaling = taxa_label(taxa_scaling_raw)
    ) %>%
    mutate(
      config_display = paste(
        learner,
        category_display,
        dimensionality,
        metab_scaling,
        taxa_scaling,
        sep = " | "
      )
    )
}

all_dat <- purrr::imap_dfr(file_map, read_perf_file)



# ------------------------------------------------------------
# 4) Keep test auroc and auprc only
# ------------------------------------------------------------
test_dat <- all_dat %>%
  filter(split == "Test", metric %in% c("AUROC", "AUPRC")) %>%
  mutate(
    dataset_label = factor(dataset_label, levels = unname(dataset_pretty))
  )

# ------------------------------------------------------------
# 5) Competitive set: within 2% of best
# ------------------------------------------------------------
threshold_pct <- 0.02

competitive_set <- test_dat %>%
  group_by(dataset, dataset_label, metric) %>%
  mutate(
    best_value = max(mean_across_seeds, na.rm = TRUE),
    threshold_value = best_value * (1 - threshold_pct),
    delta_from_best = best_value - mean_across_seeds,
    pct_from_best = (best_value - mean_across_seeds) / best_value
  ) %>%
  filter(mean_across_seeds >= threshold_value) %>%
  ungroup() %>%
  group_by(dataset, metric, mean_across_seeds) %>%
  mutate(
    n_exact_ties = n(),
    raw_tie_credit = 1 / n_exact_ties
  ) %>%
  ungroup() %>%
  group_by(dataset, metric) %>%
  mutate(
    n_in_competitive_set = n(),
    dataset_equal_credit = 1 / n_in_competitive_set
  ) %>%
  ungroup()

# ------------------------------------------------------------
# 6) Count helper
# modes:
#   raw              = exact ties share credit 1/n_ties
#   dataset_equal    = each dataset contributes total credit 1
# ------------------------------------------------------------
count_variant <- function(data, var, mode = c("raw", "dataset_equal")) {
  mode <- match.arg(mode)
  var_sym <- rlang::ensym(var)
  
  if (mode == "raw") {
    data %>% count(variant = !!var_sym, wt = raw_tie_credit, name = "n")
  } else {
    data %>% count(variant = !!var_sym, wt = dataset_equal_credit, name = "n")
  }
}

make_summary_counts <- function(df, metric_name, mode = c("raw", "dataset_equal")) {
  mode <- match.arg(mode)
  dat <- df %>% filter(metric == metric_name)
  
  list(
    learner = count_variant(dat, learner, mode),
    dimensionality = count_variant(dat, dimensionality, mode),
    metabolomics = count_variant(dat, metab_scaling, mode),
    mss = count_variant(dat, taxa_scaling, mode),
    integration = count_variant(dat, integration_family, mode),
    single_omic_split = count_variant(
      dat %>% filter(integration_family == "Single-omic"),
      category_display,
      mode
    ),
    concat = count_variant(dat %>% filter(integration_family == "Concatenation"), category_display, mode),
    stacking = count_variant(dat %>% filter(integration_family == "Stacking"), category_display, mode)
  )
}

# ------------------------------------------------------------
# 7) Horizontal bar plot
# ------------------------------------------------------------
make_hbar <- function(count_df, title_text, order_levels = NULL, fill_color = "#4C78A8") {
  count_df <- count_df %>%
    filter(!is.na(variant), n > 0)
  
  if (nrow(count_df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = paste(title_text, "\nNo entries"), size = 5) +
        theme_void()
    )
  }
  
  if (!is.null(order_levels)) {
    count_df <- count_df %>%
      mutate(variant = factor(variant, levels = order_levels)) %>%
      filter(!is.na(variant))
  } else {
    count_df <- count_df %>%
      arrange(n, variant) %>%
      mutate(variant = factor(variant, levels = variant))
  }
  
  total_n <- sum(count_df$n)
  
  count_df <- count_df %>%
    mutate(
      pct = n / total_n,
      label = paste0(round(n, 2), " (", percent(pct, accuracy = 1), ")")
    )
  
  ggplot(count_df, aes(x = n, y = variant)) +
    geom_col(fill = fill_color, width = 0.7) +
    geom_text(aes(label = label), hjust = -0.1, size = 3.5) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
    labs(
      title = title_text,
      subtitle = paste0("Total credit = ", round(total_n, 2)),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 10)
    )
}

# ------------------------------------------------------------
# 8) Summary panel using 2% competitive set
# ------------------------------------------------------------
make_bar_panel <- function(df, metric_name, mode = c("raw", "dataset_equal")) {
  mode <- match.arg(mode)
  counts <- make_summary_counts(df, metric_name, mode = mode)
  
  title_suffix <- case_when(
    mode == "raw" ~ " (within 2% of best; exact ties share credit)",
    mode == "dataset_equal" ~ " (within 2% of best; dataset-equal credit)"
  )
  
  p1 <- make_hbar(
    counts$learner,
    paste0("Learner", title_suffix),
    order_levels = c("Random Forest", "Elastic Net", "XGBoost"),
    fill_color = "#1b9e77"
  )
  
  p2 <- make_hbar(
    counts$dimensionality,
    paste0("Dimensionality", title_suffix),
    order_levels = c("Full dimensional", "Feature reduced", "Other"),
    fill_color = "#d95f02"
  )
  
  p3 <- make_hbar(
    counts$metabolomics,
    paste0("Metabolomics scaling", title_suffix),
    order_levels = c("Metab: log2", "Metab: none", "Metab: other"),
    fill_color = "#7570b3"
  )
  
  p4 <- make_hbar(
    counts$mss,
    paste0("MSS scaling", title_suffix),
    order_levels = c("MSS: CLR", "MSS: none", "MSS: other"),
    fill_color = "#e7298a"
  )
  
  p5 <- make_hbar(
    counts$integration,
    paste0("Integration family", title_suffix),
    order_levels = c("Concatenation", "Stacking", "Single-omic"),
    fill_color = "#66a61e"
  )
  
  p6 <- make_hbar(
    counts$single_omic_split,
    paste0("Single-omic split", title_suffix),
    order_levels = c("Metabolomics only", "MSS only"),
    fill_color = "#1f78b4"
  )
  
  p7 <- make_hbar(
    counts$concat,
    paste0("Concatenation method", title_suffix),
    order_levels = c("Concatenation: none", "Concatenation: scaling", "Concatenation: weighted"),
    fill_color = "#e6ab02"
  )
  
  p8 <- make_hbar(
    counts$stacking,
    paste0("Stacking method", title_suffix),
    order_levels = c("Stacking: average", "Stacking: PLS", "Stacking: NNLS", "Stacking: sparse"),
    fill_color = "#a6761d"
  )
  
  
  ((p1 | p2 | p3) / (p4 | p5| p6) / (p7 | p8)) +
    plot_annotation(
      title = paste0(metric_name, " competitive-set summary"),
      subtitle = "All test configurations within 2% of the best value in each dataset",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11)
      )
    )
}

# ------------------------------------------------------------
# 9) Optional tile plot for competitive set membership
# ------------------------------------------------------------
make_competitive_tile_plot <- function(df, metric_name) {
  row_levels <- df %>%
    filter(metric == metric_name) %>%
    distinct(config_display, learner, category_display, dimensionality, metab_scaling, taxa_scaling) %>%
    arrange(learner, category_display, dimensionality, metab_scaling, taxa_scaling, config_display) %>%
    pull(config_display)
  
  plot_dat <- df %>%
    filter(metric == metric_name) %>%
    mutate(config_display = factor(config_display, levels = rev(unique(row_levels)))) %>%
    distinct(dataset_label, config_display)
  
  ggplot(plot_dat, aes(x = dataset_label, y = config_display)) +
    geom_tile(fill = "#4C78A8", color = "white", linewidth = 0.5, width = 0.95, height = 0.9) +
    labs(
      title = paste0(metric_name, ": configurations within 2% of best"),
      subtitle = "Filled cells indicate inclusion in the competitive set",
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold", angle = 35, hjust = 1),
      axis.text.y = element_text(size = 8)
    )
}

# ------------------------------------------------------------
# 10) Build plots
# ------------------------------------------------------------
p_auroc_bar_raw <- make_bar_panel(competitive_set, "AUROC", mode = "raw")
p_auprc_bar_raw  <- make_bar_panel(competitive_set, "AUPRC", mode = "raw")

p_auroc_bar_equal <- make_bar_panel(competitive_set, "AUROC", mode = "dataset_equal")
p_auprc_bar_equal  <- make_bar_panel(competitive_set, "AUPRC", mode = "dataset_equal")

p_auroc_tile <- make_competitive_tile_plot(competitive_set, "AUROC")
p_auprc_tile  <- make_competitive_tile_plot(competitive_set, "AUPRC")

print(p_auroc_bar_raw)
print(p_auprc_bar_raw)
print(p_auroc_bar_equal)
print(p_auprc_bar_equal)
print(p_auroc_tile)
print(p_auprc_tile)

# ------------------------------------------------------------
# 11) Save plots
# ------------------------------------------------------------
#ggsave("auroc_bar_summary_within_2pct_raw.png", p_auroc_bar_raw, width = 15, height = 12, dpi = 300)
ggsave("auroc_bar_summary_within_2pct_raw.pdf", p_auroc_bar_raw, width = 15, height = 12)

#ggsave("auprc_bar_summary_within_2pct_raw.png", p_auprc_bar_raw, width = 15, height = 12, dpi = 300)
ggsave("auprc_bar_summary_within_2pct_raw.pdf", p_auprc_bar_raw, width = 15, height = 12)

#ggsave("auroc_bar_summary_within_2pct_dataset_equal.png", p_auroc_bar_equal, width = 15, height = 12, dpi = 300)
ggsave("auroc_bar_summary_within_2pct_dataset_equal.pdf", p_auroc_bar_equal, width = 15, height = 12)

#ggsave("auprc_bar_summary_within_2pct_dataset_equal.png", p_auprc_bar_equal, width = 15, height = 12, dpi = 300)
ggsave("auprc_bar_summary_within_2pct_dataset_equal.pdf", p_auprc_bar_equal, width = 15, height = 12)

#ggsave("auroc_competitive_set_tile_plot.png", p_auroc_tile, width = 14, height = 9, dpi = 300)
ggsave("auroc_competitive_set_tile_plot.pdf", p_auroc_tile, width = 14, height = 9)

#ggsave("auprc_competitive_set_tile_plot.png", p_auprc_tile, width = 14, height = 9, dpi = 300)
ggsave("auprc_competitive_set_tile_plot.pdf", p_auprc_tile, width = 14, height = 9)

# ------------------------------------------------------------
# 12) Save summary tables
# ------------------------------------------------------------
competitive_table <- competitive_set %>%
  select(
    dataset_label, metric, mean_across_seeds,
    best_value, threshold_value, delta_from_best, pct_from_best,
    n_exact_ties, raw_tie_credit,
    n_in_competitive_set, dataset_equal_credit,
    learner, category_display, integration_family,
    dimensionality, metab_scaling, taxa_scaling,
    config_display, pipeline_id
  ) %>%
  arrange(metric, dataset_label, mean_across_seeds)

readr::write_csv(competitive_table, "competitive_set_within_2pct_summary.csv")

make_frequency_table <- function(df, metric_name, mode = c("raw", "dataset_equal")) {
  mode <- match.arg(mode)
  counts <- make_summary_counts(df, metric_name, mode = mode)
  
  bind_rows(
    counts$learner %>% mutate(group = "Learner"),
    counts$dimensionality %>% mutate(group = "Dimensionality"),
    counts$metabolomics %>% mutate(group = "Metabolomics scaling"),
    counts$mss %>% mutate(group = "MSS scaling"),
    counts$integration %>% mutate(group = "Integration family"),
    counts$single_omic_split %>% mutate(group = "Single-omic split"),
    counts$concat %>% mutate(group = "Concatenation method"),
    counts$stacking %>% mutate(group = "Stacking method")
  ) %>%
    mutate(metric = metric_name, mode = mode) %>%
    select(metric, mode, group, variant, n)
}

freq_table <- bind_rows(
  make_frequency_table(competitive_set, "AUROC", mode = "raw"),
  make_frequency_table(competitive_set, "AUPRC", mode = "raw"),
  make_frequency_table(competitive_set, "AUROC", mode = "dataset_equal"),
  make_frequency_table(competitive_set, "AUPRC", mode = "dataset_equal")
)

readr::write_csv(freq_table, "competitive_set_within_2pct_frequency_summary.csv")



##########################


library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)

heatmap_value_col <- "mean_across_seeds"

category_levels <- c(
  "concat_none",
  "concat_scaling",
  "concat_scaling_weights",
  "metabolomics",
  "mss",
  "stacked_average",
  "stacked_nnls",
  "stacked_pls",
  "stacked_sparse"
)

make_heatmap_matrix <- function(df, dataset_name, metric_name, split_name) {
  df %>%
    filter(dataset == dataset_name, metric == metric_name, split == split_name) %>%
    mutate(category = factor(category, levels = category_levels)) %>%
    arrange(pipeline_id, category) %>%
    select(
      row_label = pipeline_id,
      col_label = category,
      value = all_of(heatmap_value_col)
    ) %>%
    distinct() %>%
    pivot_wider(
      names_from = col_label,
      values_from = value
    ) %>%
    arrange(row_label)
}

sanitize_sheet_name <- function(x) {
  x %>%
    str_replace_all("[\\\\/\\?\\*\\[\\]:]", "_") %>%
    str_squish()
}

metric_short <- c(
  AUROC = "AUROC",
  AUPRC = "AUPRC",
  ACCURACY = "ACC",
  BALANCED_ACCURACY = "BACC",
  F1 = "F1",
  PRECISION = "PREC",
  SENSITIVITY = "SENS"
)

make_sheet_name <- function(ds_label, metric_name, split_name, existing_names = character()) {
  metric_abbrev <- if (metric_name %in% names(metric_short)) metric_short[[metric_name]] else metric_name
  
  base <- paste(ds_label, metric_abbrev, split_name, sep = "_")
  base <- sanitize_sheet_name(base)
  
  # First truncate to Excel's 31-character limit
  base <- substr(base, 1, 31)
  
  # Ensure uniqueness even after truncation
  candidate <- base
  i <- 1
  while (tolower(candidate) %in% tolower(existing_names)) {
    suffix <- paste0("_", i)
    candidate <- substr(base, 1, 31 - nchar(suffix))
    candidate <- paste0(candidate, suffix)
    i <- i + 1
  }
  
  candidate
}

write_heatmap_sheet <- function(wb, sheet_name, title_text, subtitle_text, mat_tbl) {
  addWorksheet(wb, sheet_name)
  
  writeData(wb, sheet_name, title_text, startRow = 1, startCol = 1)
  mergeCells(wb, sheet_name, cols = 1:(ncol(mat_tbl) + 1), rows = 1)
  
  writeData(wb, sheet_name, subtitle_text, startRow = 2, startCol = 1)
  mergeCells(wb, sheet_name, cols = 1:(ncol(mat_tbl) + 1), rows = 2)
  
  writeData(wb, sheet_name, mat_tbl, startRow = 4, startCol = 1, colNames = TRUE)
  
  title_style <- createStyle(
    textDecoration = "bold",
    fontSize = 14,
    halign = "center"
  )
  
  subtitle_style <- createStyle(
    textDecoration = "italic",
    fontSize = 11,
    halign = "center"
  )
  
  header_style <- createStyle(
    textDecoration = "bold",
    fgFill = "#F2F2F2",
    halign = "center",
    valign = "center",
    border = "Bottom"
  )
  
  row_style <- createStyle(
    halign = "left",
    valign = "center"
  )
  
  num_style <- createStyle(
    numFmt = "0.0000",
    halign = "center",
    valign = "center"
  )
  
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet_name, subtitle_style, rows = 2, cols = 1, gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet_name, header_style, rows = 4, cols = 1:ncol(mat_tbl), gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet_name, row_style, rows = 5:(nrow(mat_tbl) + 4), cols = 1, gridExpand = TRUE, stack = TRUE)
  
  if (ncol(mat_tbl) > 1) {
    addStyle(
      wb, sheet_name, num_style,
      rows = 5:(nrow(mat_tbl) + 4),
      cols = 2:ncol(mat_tbl),
      gridExpand = TRUE,
      stack = TRUE
    )
  }
  
  freezePane(wb, sheet_name, firstActiveRow = 5, firstActiveCol = 2)
  
  setColWidths(wb, sheet_name, cols = 1, widths = 60)
  if (ncol(mat_tbl) > 1) {
    setColWidths(wb, sheet_name, cols = 2:ncol(mat_tbl), widths = 14)
  }
  
  setRowHeights(wb, sheet_name, rows = 1, heights = 22)
  setRowHeights(wb, sheet_name, rows = 2, heights = 18)
  setRowHeights(wb, sheet_name, rows = 4, heights = 20)
}

heatmap_wb <- createWorkbook()

for (ds in names(dataset_pretty)) {
  ds_label <- dataset_pretty[[ds]]
  
  for (metric_name in c("AUROC", "AUPRC", "ACCURACY", "BALANCED_ACCURACY", "F1", "PRECISION", "SENSITIVITY")) {
    for (split_name in c("Train", "Test")) {
      
      mat_tbl <- make_heatmap_matrix(
        df = all_dat,
        dataset_name = ds,
        metric_name = metric_name,
        split_name = split_name
      )
      
      sheet_name <- make_sheet_name(
        ds_label = ds_label,
        metric_name = metric_name,
        split_name = split_name,
        existing_names = names(heatmap_wb)
      )
      
      title_text <- paste0(ds_label, " - ", split_name, " ", metric_name, " heatmap")
      subtitle_text <- "Cell values are mean across seeds"
      
      write_heatmap_sheet(
        wb = heatmap_wb,
        sheet_name = sheet_name,
        title_text = title_text,
        subtitle_text = subtitle_text,
        mat_tbl = mat_tbl
      )
    }
  }
}

saveWorkbook(
  heatmap_wb,
  "heatmap_training_testing_matrices.xlsx",
  overwrite = TRUE
)
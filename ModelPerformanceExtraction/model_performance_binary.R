# Binary outcome ML pipeline analysis
# ------------------------------------------------------------
# Purpose:
#   Read all seed-level CSV files for a binary classification analysis,
#   parse the pipeline settings from the filenames, reshape the metric tables,
#   summarize performance across seeds, and generate:
#     1) ranked dot-whisker plots
#     2) heatmaps across pipeline/category combinations
#     3) winner-category frequency summaries and plots
#     4) HTML summary report for easy viewing
#
# Notes:
#   - Test SD/CI inside a single file may be NA because each seed may have one
#     test set, but across-seed mean/CI are still computed here.
#   - Heatmaps print values to 4 decimal places.
#   - Output filenames are prefixed with the outcome label.
#
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(fs)
  library(janitor)
  library(stringr)
  library(forcats)
  library(scales)
  library(patchwork)
  library(RColorBrewer)
})

# -----------------------------
# User settings
# -----------------------------
input_base <- "/endosome/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code/A_Integrative_Pipeline_Scripts/results/performance_csv"
sub_dir    <- "Erawijantari_Glucose"

input_dir  <- file.path(input_base, sub_dir)
output_dir <- file.path(input_dir, "results")

dir_create(output_dir)
dir_create(file.path(output_dir, "plots"))
dir_create(file.path(output_dir, "tables"))

outcome_pattern <- paste0(
  "^", sub_dir,
  "__model-(enet|rf|xgb)__metabT-(log2|none)__taxaT-(clr|none)__metabR-(limma|none)__taxaR-(wilcox|none)_run\\d{2}\\.csv$"
)

# Primary category used for pipeline-level ranking
# Options:
#   "best_per_seed"
#   "all_categories"
#   anything else -> fixed category via primary_fixed_category
primary_category_mode <- "best_per_seed"
primary_fixed_category <- "concat_none"

# Metric used to choose the best category per seed/file
primary_rank_metric <- "AUROC"

# Binary classification metrics
main_metrics <- c(
  "ACCURACY", "KAPPA", "AUROC", "AUPRC",
  "SENSITIVITY", "SPECIFICITY", "PRECISION",
  "F1", "BALANCED_ACCURACY"
)

# All are higher-is-better for this binary performance table
metric_direction <- c(
  ACCURACY = "higher",
  KAPPA = "higher",
  AUROC = "higher",
  AUPRC = "higher",
  SENSITIVITY = "higher",
  SPECIFICITY = "higher",
  PRECISION = "higher",
  F1 = "higher",
  BALANCED_ACCURACY = "higher"
)

# Metrics to display in the 8 dot-whisker plots (2 splits x 4 metrics)
plot_metrics <- c("AUROC", "AUPRC", "F1", "BALANCED_ACCURACY")

# Metrics to display as heatmaps
heatmap_metrics <- c(
  "ACCURACY", "KAPPA", "AUROC", "AUPRC",
  "SENSITIVITY", "SPECIFICITY", "PRECISION",
  "F1", "BALANCED_ACCURACY"
)

# Metrics for "winner category frequency" summaries
winner_metrics <- c("AUROC", "AUPRC", "F1", "BALANCED_ACCURACY")

# -----------------------------
# Helper functions
# -----------------------------
parse_filename <- function(path) {
  nm <- basename(path)
  
  m <- stringr::str_match(
    nm,
    "^(.*?)__model-([^_]+)__metabT-([^_]+)__taxaT-([^_]+)__metabR-([^_]+)__taxaR-([^_]+)_run(\\d{2})\\.csv$"
  )
  
  if (is.na(m[1, 1])) {
    stop(paste("Filename did not match expected pattern:", nm))
  }
  
  tibble(
    file = path,
    file_name = nm,
    outcome = m[1, 2],
    model = m[1, 3],
    metab_transform = m[1, 4],
    taxa_transform = m[1, 5],
    metab_reduction = m[1, 6],
    taxa_reduction = m[1, 7],
    seed = as.integer(m[1, 8])
  ) %>%
    mutate(
      metabolite_log2 = if_else(metab_transform == "log2", "log2", "none"),
      taxa_clr = if_else(taxa_transform == "clr", "clr", "none"),
      feature_reduction = if_else(
        metab_reduction != "none" | taxa_reduction != "none",
        "reduced",
        "full"
      ),
      reduction_signature = paste0(
        "metabR=", metab_reduction,
        "; taxaR=", taxa_reduction
      ),
      pipeline_id = paste(
        paste0("model=", model),
        paste0("metabT=", metab_transform),
        paste0("taxaT=", taxa_transform),
        paste0("metabR=", metab_reduction),
        paste0("taxaR=", taxa_reduction),
        sep = " | "
      )
    )
}
read_one_file <- function(path) {
  meta <- parse_filename(path)
  
  dat <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  
  dat_clean <- dat %>%
    clean_names() %>%
    rename(category = 1)
  
  # Robust parsing after clean_names():
  # e.g. "accuracy_train_mean", "balanced_accuracy_test_ci_lower"
  dat_long <- dat_clean %>%
    pivot_longer(
      cols = -category,
      names_to = "metric_col",
      values_to = "value"
    ) %>%
    mutate(
      metric_col = str_replace_all(metric_col, "__+", "_"),
      metric_col = str_remove(metric_col, "^x_")
    ) %>%
    extract(
      metric_col,
      into = c("metric", "split", "stat"),
      regex = "^(.*)_(train|test)_(mean|sd|n|ci_lower|ci_upper)$",
      remove = TRUE
    ) %>%
    filter(!is.na(metric), !is.na(split), !is.na(stat)) %>%
    pivot_wider(names_from = stat, values_from = value) %>%
    mutate(
      metric = str_to_upper(metric),
      split = str_to_title(split)
    )
  
  bind_cols(dat_long, meta[rep(1, nrow(dat_long)), ]) %>%
    relocate(
      file, file_name, outcome, seed, model, metabolite_log2, taxa_clr,
      feature_reduction, metab_transform, taxa_transform,
      metab_reduction, taxa_reduction, reduction_signature,
      pipeline_id, category, split, metric
    )
}

mean_se_ci <- function(x, conf = 0.95) {
  x <- x[is.finite(x)]
  n <- length(x)
  
  if (n == 0) {
    return(tibble(
      n = 0, mean = NA_real_, sd = NA_real_, se = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_
    ))
  }
  
  m <- mean(x)
  s <- sd(x)
  se <- s / sqrt(n)
  tcrit <- if (n > 1) qt((1 + conf) / 2, df = n - 1) else NA_real_
  ci_low <- if (n > 1) m - tcrit * se else NA_real_
  ci_high <- if (n > 1) m + tcrit * se else NA_real_
  
  tibble(n = n, mean = m, sd = s, se = se, ci_low = ci_low, ci_high = ci_high)
}

pick_primary_rows <- function(long_df,
                              mode = "best_per_seed",
                              fixed_category = "concat_none",
                              rank_metric = "AUROC",
                              metric_direction) {
  
  rank_df <- long_df %>%
    filter(split == "Test", metric == rank_metric) %>%
    select(file_name, pipeline_id, seed, category, mean)
  
  if (mode == "best_per_seed") {
    if (metric_direction[[rank_metric]] == "higher") {
      chosen <- rank_df %>%
        group_by(file_name) %>%
        slice_max(order_by = mean, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(file_name, chosen_category = category)
    } else {
      chosen <- rank_df %>%
        group_by(file_name) %>%
        slice_min(order_by = mean, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(file_name, chosen_category = category)
    }
  } else if (mode == "all_categories") {
    return(long_df)
  } else {
    chosen <- rank_df %>%
      filter(category == fixed_category) %>%
      distinct(file_name, chosen_category = category)
  }
  
  long_df %>%
    inner_join(chosen, by = "file_name") %>%
    filter(category == chosen_category) %>%
    select(-chosen_category)
}

safe_outcome_label <- function(x) {
  x %>%
    unique() %>%
    .[!is.na(.)] %>%
    .[1] %>%
    str_replace_all("[^A-Za-z0-9]+", "_")
}

metric_pretty <- function(x) {
  case_when(
    x == "BALANCED_ACCURACY" ~ "Balanced Accuracy",
    x == "AUROC" ~ "AUROC",
    x == "AUPRC" ~ "AUPRC",
    x == "ACCURACY" ~ "Accuracy",
    x == "KAPPA" ~ "Kappa",
    x == "SENSITIVITY" ~ "Sensitivity",
    x == "SPECIFICITY" ~ "Specificity",
    x == "PRECISION" ~ "Precision",
    x == "F1" ~ "F1",
    TRUE ~ x
  )
}

df_to_html_table <- function(df, digits = 4) {
  esc <- function(x) {
    x <- as.character(x)
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }
  
  df_fmt <- df %>%
    mutate(across(where(is.numeric), ~ ifelse(is.na(.x), "", format(round(.x, digits), nsmall = digits))))
  
  header <- paste0("<tr>", paste0("<th>", esc(names(df_fmt)), "</th>", collapse = ""), "</tr>")
  rows <- apply(df_fmt, 1, function(row) {
    paste0("<tr>", paste0("<td>", esc(row), "</td>", collapse = ""), "</tr>")
  })
  
  paste0(
    "<table border='1' cellspacing='0' cellpadding='5' style='border-collapse:collapse; margin-bottom:20px;'>",
    header,
    paste(rows, collapse = "\n"),
    "</table>"
  )
}

# -----------------------------
# Discover files and report completeness
# -----------------------------
all_files_found <- list.files(
  path = input_dir,
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

files <- all_files_found[
  grepl(outcome_pattern, basename(all_files_found))
]

message("Total CSVs found recursively: ", length(all_files_found))
message("Matched binary seed-result CSVs: ", length(files))

if (length(files) == 0) {
  stop("No files matched outcome_pattern. Check input_dir and filename pattern.")
}

bad_files <- files[!grepl(outcome_pattern, basename(files))]
if (length(bad_files) > 0) {
  print(basename(bad_files))
  stop("Some files do not match the expected naming convention.")
}

file_index <- purrr::map_dfr(files, parse_filename)

expected_grid <- tidyr::expand_grid(
  model = c("enet", "rf", "xgb"),
  metab_transform = c("log2", "none"),
  taxa_transform = c("clr", "none"),
  metab_reduction = c("limma", "none"),
  taxa_reduction = c("wilcox", "none"),
  seed = 1:20
) %>%
  mutate(
    valid_reduction_combo = case_when(
      metab_reduction == "limma" & taxa_reduction == "wilcox" ~ TRUE,
      metab_reduction == "none" & taxa_reduction == "none" ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(valid_reduction_combo) %>%
  select(-valid_reduction_combo) %>%
  mutate(expected_key = paste(model, metab_transform, taxa_transform, metab_reduction, taxa_reduction, seed, sep = "__"))

observed_grid <- file_index %>%
  mutate(observed_key = paste(model, metab_transform, taxa_transform, metab_reduction, taxa_reduction, seed, sep = "__"))

missing_runs <- expected_grid %>%
  anti_join(observed_grid, by = c(
    "model", "metab_transform", "taxa_transform",
    "metab_reduction", "taxa_reduction", "seed"
  ))

write_csv(file_index, file.path(output_dir, "tables", "file_index.csv"))
write_csv(missing_runs, file.path(output_dir, "tables", "missing_expected_runs.csv"))

message(glue("Found {nrow(file_index)} files."))
message(glue("Expected 480 total files for the full design; missing {nrow(missing_runs)} runs based on the naming scheme."))

# -----------------------------
# Read and tidy all files
# -----------------------------
all_long <- map_dfr(files, read_one_file) %>%
  filter(metric %in% main_metrics) %>%
  mutate(
    category = factor(category),
    model = factor(model, levels = c("enet", "rf", "xgb")),
    metabolite_log2 = factor(metabolite_log2, levels = c("none", "log2")),
    taxa_clr = factor(taxa_clr, levels = c("none", "clr")),
    feature_reduction = factor(feature_reduction, levels = c("full", "reduced"))
  )

write_csv(all_long, file.path(output_dir, "tables", "all_metrics_long.csv"))

outcome_label <- safe_outcome_label(all_long$outcome)
message("Outcome label for output files: ", outcome_label)

# -----------------------------
# Choose rows to plot
# -----------------------------
plot_long <- pick_primary_rows(
  all_long,
  mode = primary_category_mode,
  fixed_category = primary_fixed_category,
  rank_metric = primary_rank_metric,
  metric_direction = metric_direction
)

# -----------------------------
# Across-seed summaries for the dot-whisker plots
# -----------------------------
plot_summary <- plot_long %>%
  filter(
    split %in% c("Train", "Test"),
    metric %in% plot_metrics
  ) %>%
  group_by(pipeline_id, split, metric) %>%
  group_modify(~ mean_se_ci(.x$mean)) %>%
  ungroup() %>%
  mutate(direction = unname(metric_direction[metric]))

write_csv(
  plot_summary,
  file.path(output_dir, "tables", paste0(outcome_label, "_train_test_metric_summary_across_seeds.csv"))
)

make_metric_plot <- function(df, metric_name, split_name) {
  d <- df %>%
    filter(metric == metric_name, split == split_name) %>%
    mutate(
      pipeline_id = if (metric_direction[[metric_name]] == "higher") {
        fct_reorder(pipeline_id, mean, .desc = FALSE)
      } else {
        fct_reorder(pipeline_id, mean, .desc = TRUE)
      }
    )
  
  y_lab <- paste0(metric_pretty(metric_name), " (mean across seeds, 95% CI)")
  
  ggplot(d, aes(x = pipeline_id, y = mean)) +
    geom_errorbar(
      aes(ymin = ci_low, ymax = ci_high),
      width = 0.15,
      linewidth = 0.5,
      na.rm = TRUE
    ) +
    geom_point(size = 2.6) +
    coord_flip() +
    labs(
      title = paste(outcome_label, "-", split_name, metric_pretty(metric_name)),
      subtitle = "Points = mean across seeds; whiskers = 95% CI",
      x = NULL,
      y = y_lab
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )
}

plot_list <- crossing(
  split = c("Train", "Test"),
  metric = plot_metrics
) %>%
  mutate(
    plot_obj = purrr::map2(metric, split, ~ make_metric_plot(plot_summary, .x, .y))
  )

pwalk(
  list(plot_list$plot_obj, plot_list$metric, plot_list$split),
  function(p, metric_name, split_name) {
    file_stub <- paste0(
      outcome_label, "_performance_",
      tolower(split_name), "_",
      tolower(gsub("[^A-Za-z0-9]+", "_", metric_name))
    )
    
    ggsave(
      filename = file.path(output_dir, "plots", paste0(file_stub, ".png")),
      plot = p,
      width = 9,
      height = 7,
      dpi = 300
    )
    
    ggsave(
      filename = file.path(output_dir, "plots", paste0(file_stub, ".pdf")),
      plot = p,
      width = 9,
      height = 7
    )
  }
)

train_panel <- (
  plot_list$plot_obj[[which(plot_list$split == "Train" & plot_list$metric == "AUROC")]] +
    plot_list$plot_obj[[which(plot_list$split == "Train" & plot_list$metric == "AUPRC")]]
) / (
  plot_list$plot_obj[[which(plot_list$split == "Train" & plot_list$metric == "F1")]] +
    plot_list$plot_obj[[which(plot_list$split == "Train" & plot_list$metric == "BALANCED_ACCURACY")]]
)

test_panel <- (
  plot_list$plot_obj[[which(plot_list$split == "Test" & plot_list$metric == "AUROC")]] +
    plot_list$plot_obj[[which(plot_list$split == "Test" & plot_list$metric == "AUPRC")]]
) / (
  plot_list$plot_obj[[which(plot_list$split == "Test" & plot_list$metric == "F1")]] +
    plot_list$plot_obj[[which(plot_list$split == "Test" & plot_list$metric == "BALANCED_ACCURACY")]]
)

ggsave(
  file.path(output_dir, "plots", paste0(outcome_label, "_performance_train_panel.png")),
  train_panel,
  width = 14,
  height = 12,
  dpi = 300
)

ggsave(
  file.path(output_dir, "plots", paste0(outcome_label, "_performance_test_panel.png")),
  test_panel,
  width = 14,
  height = 12,
  dpi = 300
)

# -----------------------------
# Heatmaps with values in each cell (4 decimals)
# -----------------------------
heatmap_summary <- all_long %>%
  filter(metric %in% heatmap_metrics) %>%
  group_by(pipeline_id, category, split, metric) %>%
  summarise(
    mean_across_seeds = mean(mean, na.rm = TRUE),
    sd_across_seeds = sd(mean, na.rm = TRUE),
    n_seeds = sum(is.finite(mean)),
    .groups = "drop"
  ) %>%
  mutate(metric_pretty = metric_pretty(metric))

write_csv(
  heatmap_summary,
  file.path(output_dir, "tables", paste0(outcome_label, "_heatmap_summary.csv"))
)

make_perf_heatmap <- function(df, split_name, metric_name) {
  d <- df %>%
    filter(split == split_name, metric == metric_name) %>%
    mutate(
      pipeline_id = if (metric_direction[[metric_name]] == "higher") {
        fct_reorder(pipeline_id, mean_across_seeds, .desc = FALSE)
      } else {
        fct_reorder(pipeline_id, mean_across_seeds, .desc = TRUE)
      }
    )
  
  ggplot(d, aes(x = category, y = pipeline_id, fill = mean_across_seeds)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(aes(label = sprintf("%.4f", mean_across_seeds)), size = 2.6) +
    scale_fill_gradient(
      low = "white",
      high = "red",
      na.value = "grey90"
    ) +
    labs(
      title = paste(outcome_label, "-", split_name, metric_pretty(metric_name), "heatmap"),
      subtitle = "Cell values are mean across seeds",
      x = "Category",
      y = "Pipeline",
      fill = "Mean"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 7)
    )
}

heatmap_grid <- crossing(
  split = c("Train", "Test"),
  metric = heatmap_metrics
) %>%
  mutate(
    plot_obj = purrr::map2(split, metric, ~ make_perf_heatmap(heatmap_summary, .x, .y))
  )

pwalk(
  list(heatmap_grid$plot_obj, heatmap_grid$split, heatmap_grid$metric),
  function(p, split_name, metric_name) {
    stub <- paste0(
      outcome_label, "_heatmap_",
      tolower(split_name), "_",
      tolower(gsub("[^A-Za-z0-9]+", "_", metric_name))
    )
    
    ggsave(
      file.path(output_dir, "plots", paste0(stub, ".png")),
      p,
      width = 12,
      height = 8,
      dpi = 300
    )
    
    ggsave(
      file.path(output_dir, "plots", paste0(stub, ".pdf")),
      p,
      width = 12,
      height = 8
    )
  }
)

# -----------------------------
# Winner category frequency for higher-is-better metrics
# -----------------------------
winner_category_by_seed <- all_long %>%
  filter(metric %in% winner_metrics) %>%
  group_by(file_name, pipeline_id, seed, split, metric) %>%
  slice_max(order_by = mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(metric_pretty = metric_pretty(metric))

write_csv(
  winner_category_by_seed,
  file.path(output_dir, "tables", paste0(outcome_label, "_winner_category_by_seed.csv"))
)

winner_category_frequency <- winner_category_by_seed %>%
  count(split, metric, metric_pretty, category, name = "wins") %>%
  group_by(split, metric, metric_pretty) %>%
  mutate(
    total = sum(wins),
    prop = wins / total
  ) %>%
  ungroup()

write_csv(
  winner_category_frequency,
  file.path(output_dir, "tables", paste0(outcome_label, "_winner_category_frequency.csv"))
)

make_winner_category_plot <- function(df, split_name, metric_name_pretty) {
  d <- df %>%
    filter(split == split_name, metric_pretty == metric_name_pretty) %>%
    mutate(category = fct_reorder(category, wins, .desc = FALSE))
  
  ggplot(d, aes(x = category, y = wins)) +
    geom_col(width = 0.75) +
    geom_text(aes(label = wins), hjust = -0.15, size = 3.5) +
    coord_flip() +
    labs(
      title = paste(outcome_label, "-", split_name, metric_name_pretty, "winner category frequency"),
      subtitle = "Count of pipeline-seed combinations where each category had the best value",
      x = NULL,
      y = "Number of wins"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold")
    ) +
    expand_limits(y = max(d$wins, na.rm = TRUE) * 1.08)
}

winner_plot_grid <- tidyr::crossing(
  split = c("Train", "Test"),
  metric_pretty = c("AUROC", "AUPRC", "F1", "Balanced Accuracy")
) %>%
  mutate(
    plot_obj = purrr::map2(split, metric_pretty,
                           ~ make_winner_category_plot(winner_category_frequency, .x, .y))
  )

pwalk(
  list(winner_plot_grid$plot_obj, winner_plot_grid$split, winner_plot_grid$metric_pretty),
  function(p, split_name, metric_name_pretty) {
    stub <- paste0(
      outcome_label, "_winner_category_frequency_",
      tolower(split_name), "_",
      tolower(gsub("[^A-Za-z0-9]+", "_", metric_name_pretty))
    )
    
    ggsave(
      file.path(output_dir, "plots", paste0(stub, ".png")),
      p,
      width = 8,
      height = 5.5,
      dpi = 300
    )
    
    ggsave(
      file.path(output_dir, "plots", paste0(stub, ".pdf")),
      p,
      width = 8,
      height = 5.5
    )
  }
)

train_winner_panel <- (
  winner_plot_grid$plot_obj[[which(winner_plot_grid$split == "Train" & winner_plot_grid$metric_pretty == "AUROC")]] /
    winner_plot_grid$plot_obj[[which(winner_plot_grid$split == "Train" & winner_plot_grid$metric_pretty == "AUPRC")]] /
    winner_plot_grid$plot_obj[[which(winner_plot_grid$split == "Train" & winner_plot_grid$metric_pretty == "F1")]] /
    winner_plot_grid$plot_obj[[which(winner_plot_grid$split == "Train" & winner_plot_grid$metric_pretty == "Balanced Accuracy")]]
)

test_winner_panel <- (
  winner_plot_grid$plot_obj[[which(winner_plot_grid$split == "Test" & winner_plot_grid$metric_pretty == "AUROC")]] /
    winner_plot_grid$plot_obj[[which(winner_plot_grid$split == "Test" & winner_plot_grid$metric_pretty == "AUPRC")]] /
    winner_plot_grid$plot_obj[[which(winner_plot_grid$split == "Test" & winner_plot_grid$metric_pretty == "F1")]] /
    winner_plot_grid$plot_obj[[which(winner_plot_grid$split == "Test" & winner_plot_grid$metric_pretty == "Balanced Accuracy")]]
)

ggsave(
  file.path(output_dir, "plots", paste0(outcome_label, "_winner_category_frequency_train_panel.png")),
  train_winner_panel,
  width = 10,
  height = 17,
  dpi = 300
)

ggsave(
  file.path(output_dir, "plots", paste0(outcome_label, "_winner_category_frequency_test_panel.png")),
  test_winner_panel,
  width = 10,
  height = 17,
  dpi = 300
)

# -----------------------------
# HTML summary report
# -----------------------------
html_file <- file.path(output_dir, paste0(outcome_label, "_summary_report.html"))

summary_table_html <- df_to_html_table(
  plot_summary %>%
    arrange(split, metric, desc(mean)) %>%
    mutate(metric = metric_pretty(metric)),
  digits = 4
)

winner_table_html <- df_to_html_table(
  winner_category_frequency %>%
    arrange(split, metric_pretty, desc(wins)),
  digits = 4
)

heatmap_table_html <- df_to_html_table(
  heatmap_summary %>%
    select(split, metric_pretty, pipeline_id, category, mean_across_seeds, sd_across_seeds, n_seeds) %>%
    arrange(split, metric_pretty, pipeline_id, category),
  digits = 4
)

heatmap_img_lines <- unlist(lapply(c("Train", "Test"), function(sp) {
  unlist(lapply(heatmap_metrics, function(mt) {
    stub <- paste0(
      outcome_label, "_heatmap_",
      tolower(sp), "_",
      tolower(gsub("[^A-Za-z0-9]+", "_", mt))
    )
    paste0("<h3>", sp, " ", metric_pretty(mt), "</h3><img src='plots/", stub, ".png'/>")
  }))
}))

html_lines <- c(
  "<!DOCTYPE html>",
  "<html>",
  "<head>",
  paste0("<title>", outcome_label, " summary report</title>"),
  "<meta charset='utf-8'/>",
  "<style>",
  "body { font-family: Arial, sans-serif; margin: 30px; }",
  "h1, h2, h3 { color: #222; }",
  "img { max-width: 100%; height: auto; border: 1px solid #ddd; margin-bottom: 20px; }",
  "table { font-size: 12px; }",
  "th { background: #f2f2f2; }",
  ".section { margin-bottom: 40px; }",
  "</style>",
  "</head>",
  "<body>",
  paste0("<h1>", outcome_label, " summary report</h1>"),
  
  "<div class='section'>",
  "<h2>Performance panels</h2>",
  paste0("<h3>Train</h3><img src='plots/", outcome_label, "_performance_train_panel.png'/>"),
  paste0("<h3>Test</h3><img src='plots/", outcome_label, "_performance_test_panel.png'/>"),
  "</div>",
  
  "<div class='section'>",
  "<h2>Winner category frequency panels</h2>",
  paste0("<h3>Train</h3><img src='plots/", outcome_label, "_winner_category_frequency_train_panel.png'/>"),
  paste0("<h3>Test</h3><img src='plots/", outcome_label, "_winner_category_frequency_test_panel.png'/>"),
  "</div>",
  
  "<div class='section'>",
  "<h2>Heatmaps</h2>",
  heatmap_img_lines,
  "</div>",
  
  "<div class='section'>",
  "<h2>Across-seed summary table</h2>",
  summary_table_html,
  "</div>",
  
  "<div class='section'>",
  "<h2>Winner category frequency table</h2>",
  winner_table_html,
  "</div>",
  
  "<div class='section'>",
  "<h2>Heatmap summary table</h2>",
  heatmap_table_html,
  "</div>",
  
  "</body>",
  "</html>"
)

writeLines(html_lines, con = html_file)

message("HTML summary report written to: ", html_file)
message("Done.")
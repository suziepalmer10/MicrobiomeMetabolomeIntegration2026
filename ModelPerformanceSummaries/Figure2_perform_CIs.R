library(tidyverse)

base_dir <- "/Users/suzettepalmer/Library/CloudStorage/Box-Box/Personal - Suzette Palmer/IntegratedLearner_Revisions/Results/Performance_CSVs"

datasets <- c(
  "Erawijantari_Cholesterol",
  "Erawijantari_Glucose",
  "Franzosa_CD_Fp",
  "Franzosa_IBD_Fp",
  "Franzosa_UC_Fp",
  "Wang_Creatinine",
  "Wang_eGFR",
  "Wang_Urea"
)

pdf_device <- function(filename, width, height, ...) {
  grDevices::pdf(filename, width = width, height = height, useDingbats = FALSE)
}

safe_ci <- function(x) {
  x <- x[!is.na(x)]
  n <- length(x)
  
  if (n == 0) {
    return(tibble(
      mean_value = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      n = 0L
    ))
  }
  
  m <- mean(x)
  
  if (n == 1) {
    return(tibble(
      mean_value = m,
      se = NA_real_,
      ci_lower = m,
      ci_upper = m,
      n = 1L
    ))
  }
  
  s <- sd(x)
  se <- s / sqrt(n)
  tcrit <- qt(0.975, df = n - 1)
  
  tibble(
    mean_value = m,
    se = se,
    ci_lower = m - tcrit * se,
    ci_upper = m + tcrit * se,
    n = n
  )
}

make_plot <- function(summary_df, split_value, metric_value, title_text) {
  plot_df <- summary_df %>%
    filter(split == split_value, metric == metric_value, n > 0) %>%
    mutate(
      config = factor(config, levels = unique(config))
    )
  
  ggplot(plot_df, aes(x = config, y = mean_value)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15) +
    facet_wrap(~ category, scales = "free_y", drop = TRUE) +
    labs(
      x = "Configuration",
      y = "Mean across seeds (95% CI)",
      title = title_text
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold")
    )
}

for (dataset in datasets) {
  csv_path <- file.path(base_dir, dataset, "results", "tables", "all_metrics_long.csv")
  
  if (!file.exists(csv_path)) {
    message("Skipping missing file: ", csv_path)
    next
  }
  
  df <- read_csv(csv_path, show_col_types = FALSE)
  
  summary_df <- df %>%
    filter(
      metric %in% c("RMSE", "MAE"),
      split %in% c("Train", "Test")
    ) %>%
    mutate(
      config = paste(
        model,
        paste0("metT-", metab_transform),
        paste0("taxT-", taxa_transform),
        paste0("metR-", metab_reduction),
        paste0("taxR-", taxa_reduction),
        sep = "_"
      )
    ) %>%
    group_by(category, split, metric, config) %>%
    summarise(
      safe_ci(mean),
      .groups = "drop"
    )
  
  out_dir <- file.path(base_dir, dataset, "results")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (split_value in c("Train", "Test")) {
    for (metric_value in c("RMSE", "MAE")) {
      plot_title <- paste0(dataset, " - ", split_value, " ", metric_value)
      out_file <- file.path(out_dir, paste0(dataset, "_", split_value, "_", metric_value, ".pdf"))
      
      p <- make_plot(summary_df, split_value, metric_value, plot_title)
      
      ggsave(out_file, p, width = 16, height = 8, device = pdf_device)
      message("Wrote: ", out_file)
    }
  }
}
# ---------------------------
# Project root (shared folder)
# ---------------------------
PROJECT_ROOT <- "/work/PCDC/shared/IntegratedLearner_202412/A_Integrative_Pipeline_Code"
SCRIPTS_DIR  <- file.path(PROJECT_ROOT, "A_Integrative_Pipeline_Scripts")
CONFIG_DIR   <- file.path(PROJECT_ROOT, "C_Configuration_Files")

# Config file to run
cfg <- file.path(CONFIG_DIR, "Erawijantari", "Erawijantari_con_vs_gas.txt")

# SECTION 1: Load required packages
pkgs_needed_cran <- c(
  "caret", "kernlab", "xgboost", "tidyverse", "MLmetrics",
  "nnls", "argparse", "glmnet", "pls", "pROC", "irr", "randomForest",
  "readr", "PRROC"
)

# Install/load CRAN packages
for (pkg in pkgs_needed_cran) {
  if (!requireNamespace(pkg, quietly=TRUE)) {
    stop("Missing package: ", pkg, " (install it before running on the cluster).")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# limma is Bioconductor
if (!requireNamespace("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("limma", ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages(library(limma))

set.seed(123)

# SECTION 2: Parse CLI args or source config in interactive mode
parser <- ArgumentParser(description = "Process integrated learner arguments")
parser$add_argument("--model_to_run", type = "character")
parser$add_argument("--file_path", type = "character")
parser$add_argument("--input_file", type = "character")
parser$add_argument("--study_name", type = "character")
parser$add_argument("--type_of_analysis", type = "character")
parser$add_argument("--response_variable", type = "character")
parser$add_argument("--stratify_variable", type = "character", default = "")
parser$add_argument("--training_proportion", type = "numeric", default = 0.8)
parser$add_argument("--num_repeats", type = "numeric", default = 3)
parser$add_argument("--num_folds", type = "numeric", default = 5)
parser$add_argument("--metab_transform", type = "character", default = "log2")
parser$add_argument("--metab_reduction", type = "character", default = "none")
parser$add_argument("--taxa_sum_threshold", type = "numeric", default = 1e-4)
parser$add_argument("--taxa_reduction", type = "character", default = "none")
parser$add_argument("--corr_cutoff", type = "numeric", default = 0.90)
parser$add_argument("--pseudocount", type = "numeric", default = 1e-9)
parser$add_argument("--taxa_transform", type = "character", default = "none",
                    help = "Taxa transform: none|clr (centered log-ratio)")
parser$add_argument("--use_gpu", type = "character", default = "FALSE",
                    help = "Whether to use GPU for xgboost: TRUE|FALSE")
parser$add_argument("--gpu_device", type = "character", default = "cuda",
                    help = "XGBoost device string, e.g. cuda or cuda:0")
parser$add_argument("--nthread", type = "integer", default = 16,
                    help = "Number of CPU threads for model fitting")

# ---- Helpers to match bash naming exactly ----
sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

infer_model_key <- function(path) {
  fname <- tolower(basename(path))
  if (grepl("enet", fname)) return("enet")
  if (grepl("rf", fname)) return("rf")
  if (grepl("xgboost|xgb", fname)) return("xgb")
  stop("Could not infer model key from model_to_run: ", fname)
}

build_study_name <- function(base, model_key, metab_t, taxa_t, metab_r, taxa_r) {
  paste0(
    sanitize(base),
    "__model-", sanitize(model_key),
    "__metabT-", sanitize(metab_t),
    "__taxaT-",  sanitize(taxa_t),
    "__metabR-", sanitize(metab_r),
    "__taxaR-",  sanitize(taxa_r)
  )
}

if (interactive()) {
  if (!file.exists(cfg)) stop("Config file not found: ", cfg)
  source(cfg)
  file_path <- PROJECT_ROOT
  
  model_to_run <- file.path(SCRIPTS_DIR, "model_functions", "enet_function.R")
  metab_transform <- "none"
  metab_reduction <- "limma"
  taxa_reduction  <- "wilcox"
  taxa_transform  <- "none"
  corr_cutoff <- 0.98
  taxa_sum_threshold <- 1e-4
  pseudocount <- 1e-9
  use_gpu <- TRUE
  gpu_device <- "cuda"
  nthread <- 16
  
  base <- if (exists("study_name_base") && nzchar(study_name_base)) {
    study_name_base
  } else if (exists("study_name") && nzchar(study_name)) {
    study_name
  } else {
    stop("Config must define study_name_base or study_name")
  }
  
  model_key <- infer_model_key(model_to_run)
  study_name <- build_study_name(
    base      = base,
    model_key = model_key,
    metab_t   = metab_transform,
    taxa_t    = taxa_transform,
    metab_r   = metab_reduction,
    taxa_r    = taxa_reduction
  )
  
  parsed_args <- list(
    model_to_run = model_to_run,
    file_path = file_path,
    input_file = input_file,
    study_name = study_name,
    type_of_analysis = type_of_analysis,
    response_variable = response_variable,
    stratify_variable = stratify_variable,
    training_proportion = training_proportion,
    num_repeats = num_repeats,
    num_folds = num_folds,
    metab_transform = metab_transform,
    metab_reduction = metab_reduction,
    taxa_reduction = taxa_reduction,
    taxa_transform = taxa_transform,
    corr_cutoff = corr_cutoff,
    taxa_sum_threshold = taxa_sum_threshold,
    pseudocount = pseudocount,
    use_gpu = use_gpu,
    gpu_device = gpu_device,
    nthread = nthread
  )
} else {
  parsed_args <- as.list(parser$parse_args())
}

as_logical_flag <- function(x) {
  if (is.logical(x)) return(x)
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop("Could not parse logical flag from: ", x)
}

parsed_args$use_gpu <- as_logical_flag(parsed_args$use_gpu)

# keep the rest of the script using args
args <- parsed_args

# SECTION 3: Resolve paths + basic checks
is_abs <- function(p) grepl("^(/|[A-Za-z]:\\\\|[A-Za-z]:/|~)", p)

expand_path <- function(p) {
  if (is.null(p) || is.na(p) || p == "") return(p)
  path.expand(p)
}

args$file_path    <- expand_path(args$file_path)
args$model_to_run <- expand_path(args$model_to_run)
args$input_file   <- expand_path(args$input_file)



# Global results directory (HPC storage location)
# Global results directory (within the project tree)
results_base_dir <- file.path(args$file_path, "A_Integrative_Pipeline_Scripts", "results")
results_base_dir <- path.expand(results_base_dir)


################################
dir.create(results_base_dir, recursive = TRUE, showWarnings = FALSE)

message("All results will be written to: ", results_base_dir)


# ---------------------------
# Checkpoint directory for completed seeds/runs
# ---------------------------
completed_dir <- file.path(results_base_dir, "CompletedSeeds")
dir.create(completed_dir, recursive = TRUE, showWarnings = FALSE)

# Helper: given base args$study_name and run_id, compute the study_name used inside run_one
study_name_for_run <- function(base_study_name, run_id) {
  paste0(base_study_name, "_run", sprintf("%02d", run_id))
}

# Helper: path to done marker for a run (atomic file)
done_marker_for_run <- function(study_name_run) {
  file.path(completed_dir, paste0(study_name_run, ".done"))
}

# Helper: path to saved rdata (where run_one writes RData)
rdata_for_run <- function(study_name_run) {
  file.path(results_base_dir, "RDataFiles", paste0(study_name_run, ".RData"))
}

if (is.null(args$file_path) || !dir.exists(args$file_path)) {
  stop("file_path does not exist: ", args$file_path)
}

model_path <- if (is_abs(args$model_to_run)) {
  args$model_to_run
} else {
  file.path(args$file_path, args$model_to_run)
}

input_path <- if (is_abs(args$input_file)) {
  args$input_file
} else {
  # try relative to scripts dir first (most common)
  file.path(args$file_path, "A_Integrative_Pipeline_Scripts", args$input_file)
}

if (!file.exists(model_path)) stop("Model file not found: ", model_path)
if (!file.exists(input_path)) stop("Input file not found: ", input_path)


# Learner loader (robust, prevents function overwrites)

infer_model_type_from_path <- function(path) {
  fname <- tolower(basename(path))
  if (grepl("enet", fname)) return("enet")
  if (grepl("rf", fname)) return("rf")
  if (grepl("xgboost|xgb", fname)) return("xgboost")
  stop("Could not infer model type from model_to_run: ", fname)
}

load_learner <- function(model_path, model_type) {
  env <- new.env(parent = baseenv())
  sys.source(model_path, envir = env)
  
  required <- c("fit", "predict")
  missing <- required[!vapply(required, function(nm) exists(nm, envir = env, inherits = FALSE), logical(1))]
  if (length(missing) > 0) {
    stop("Model script must define functions: ", paste(missing, collapse = ", "),
         ". File: ", model_path)
  }
  
  # Optional: allow model scripts to supply importance directly
  if (!exists("get_importance", envir = env, inherits = FALSE)) {
    env$get_importance <- NULL
  }
  
  # Helpful metadata
  env$model_type <- model_type
  env
}


# SECTION 4: Load model + data
model_type <- infer_model_type_from_path(args$model_to_run)
message("Detected base learner: ", model_type)

model_env <- load_learner(model_path, model_type)
if (!file.exists(model_path)) stop("Model file not found: ", model_path)
if (!file.exists(input_path)) stop("Input file not found: ", input_path)

data <- readr::read_csv(input_path, show_col_types = FALSE)

# Ensure stable, unique column names everywhere (prevents limma/df mismatch)
names(data) <- make.unique(names(data))

if (!(args$response_variable %in% names(data))) {
  stop("response_variable not found in data: ", args$response_variable)
}

# Variables for downstream code (optional convenience)
type_of_analysis <- args$type_of_analysis
response_variable <- args$response_variable
stratify_variable <- args$stratify_variable
training_proportion <- args$training_proportion
num_repeats <- args$num_repeats
num_folds <- args$num_folds
file_path <- args$file_path
model_to_run <- args$model_to_run
input_file <- args$input_file

run_one <- function(args, run_id, seed, return_big = FALSE) {
  set.seed(seed)
  
  # Make a unique study_name per run so files don't overwrite
  args$study_name <- paste0(args$study_name, "_run", sprintf("%02d", run_id))
  study_name <- args$study_name
  while (sink.number() > 0) sink()
  
  output_file <- file.path(
    results_base_dir,
    "AnalysisRunDocumentation",
    paste0(study_name, ".txt")
  )
  
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(output_file)) {
    file.remove(output_file)
    message("Existing file removed: ", output_file)
  }
  sink(output_file)
  options(warn = 1)
  
  # file metadata
  print(paste('The model ran for this analysis is: ', model_to_run))
  print(paste('The filepath for this analysis is: ', file_path))
  print(paste('The input file for this analysis is: ', input_file))
  print(paste('The study name for this analysis is: ', study_name))
  print(paste('The type of performance metrics for this analysis is: ', type_of_analysis))
  print(paste('The name of the response variable for this anlaysis is: ', response_variable))
  print(paste('The name of the stratify variable, if used, for this analysis is: ', stratify_variable))
  print(paste('The proportion of data used for training is: ', training_proportion))
  print(paste('The number of repeats used for Repeated V-fold cross validation is: ', num_repeats))
  print(paste('The number of folds used for Repeated V-fold cross validation is: ', num_folds))
  message("Metabolomics transform: ", args$metab_transform, " (pseudocount=", args$pseudocount, ")")
  message("Taxa transform: ", args$taxa_transform, " (pseudocount=", args$pseudocount, ")")
  message("use_gpu: ", args$use_gpu)
  message("gpu_device: ", args$gpu_device)
  message("nthread: ", args$nthread)
  
  start_time <- Sys.time()
  
  # SECTION 2: LOAD IN DATA AND CHECK DATA TYPE AND LENGTH
  get_columns_by_prefix <- function(prefix) {
    grep(paste0('^', prefix), colnames(data), value = TRUE)
  }
  
  summarise_missing <- function(cols, label) {
    print(paste('The number of', label, 'are:', length(cols)))
    n_na <- sum(sapply(data[cols], function(x) any(is.na(x))))
    print(paste(label, 'number of columns with NA:', n_na))
  }
  
  convert_and_check_numeric <- function(cols) {
    for (col in cols) {
      if (!is.numeric(data[[col]])) {
        data[[col]] <<- as.numeric(data[[col]])
        if (any(is.na(data[[col]]))) {
          warning(paste('Some values in column', col,
                        'could not be converted to numeric and have been set to NA.'))
        }
      }
    }
    if (!all(sapply(data[cols], is.numeric))) {
      stop('Not all specified columns are numeric. Please check the data.')
    }
  }
  
  # features for metabolites and taxa
  m_columns <- get_columns_by_prefix('m__')
  t_columns <- get_columns_by_prefix('t__')
  
  # Apply global taxa sum threshold across ALL samples
  taxa_sum_threshold <- args$taxa_sum_threshold
  if (length(t_columns) > 0) {
    taxa_sums_all <- colSums(data[, t_columns, drop = FALSE], na.rm = TRUE)
    keep_taxa_global <- names(taxa_sums_all)[taxa_sums_all >= taxa_sum_threshold]
    drop_taxa_global <- setdiff(t_columns, keep_taxa_global)
    message(sprintf("Global taxa filter: dropping %d taxa with column sum < %g (keeping %d).",
                    length(drop_taxa_global), taxa_sum_threshold, length(keep_taxa_global)))
    t_columns <- keep_taxa_global
  } else {
    message("No taxa columns (t__) found in data.")
  }
  
  # Helper to unwrap models that return list(model=caret_train_obj)
  unwrap_fit <- function(fit_obj) {
    if (is.list(fit_obj) && "model" %in% names(fit_obj)) return(fit_obj$model)
    fit_obj
  }
  
  extract_fold_importance <- function(fit_obj, model_type, view_name,
                                      run_id, seed, repeat_id, fold,
                                      split = c("validation", "test")) {
    split <- match.arg(split)
    
    fit <- unwrap_fit(fit_obj)
    
    meta_cols <- data.frame(
      View = view_name,
      ModelType = model_type,
      Split = split,
      RunID = run_id,
      Seed = seed,
      RepeatID = repeat_id,
      Fold = fold,
      stringsAsFactors = FALSE
    )
    
    # ---------------- ENET / caret glmnet ----------------
    if (model_type == "enet") {
      bt <- fit$bestTune
      lam <- bt$lambda
      alpha <- bt$alpha
      
      cm <- as.matrix(glmnet::coef.glmnet(fit$finalModel, s = lam))
      df <- data.frame(
        Feature = rownames(cm),
        Coef = as.numeric(cm[, 1]),
        stringsAsFactors = FALSE
      )
      df <- df[df$Feature != "(Intercept)", , drop = FALSE]
      df$Importance <- abs(df$Coef)
      df$ImportanceType <- "abs(coef)"
      df$Alpha <- alpha
      df$Lambda <- lam
      
      out <- cbind(meta_cols, df)
      out <- out[order(-out$Importance), , drop = FALSE]
      rownames(out) <- NULL
      return(out)
    }
    
    # ---------------- RF / caret rf ----------------
    if (model_type == "rf") {
      vi <- tryCatch(caret::varImp(fit, scale = FALSE), error = function(e) NULL)
      if (is.null(vi)) {
        return(cbind(meta_cols, data.frame(
          Feature=character(), Importance=numeric(), ImportanceType=character(),
          stringsAsFactors=FALSE
        )))
      }
      
      imp <- vi$importance
      imp$Feature <- rownames(imp)
      
      if ("Overall" %in% colnames(imp)) {
        imp$Importance <- imp$Overall
        itype <- "varImp(Overall)"
      } else {
        cols <- setdiff(colnames(imp), "Feature")
        imp$Importance <- rowMeans(imp[, cols, drop = FALSE], na.rm = TRUE)
        itype <- "varImp(mean)"
      }
      
      df <- data.frame(
        Feature = imp$Feature,
        Importance = as.numeric(imp$Importance),
        ImportanceType = itype,
        stringsAsFactors = FALSE
      )
      
      out <- cbind(meta_cols, df)
      out <- out[order(-out$Importance), , drop = FALSE]
      rownames(out) <- NULL
      return(out)
    }
    
    # ---------------- XGBOOST / custom xgb wrapper ----------------
    if (model_type == "xgboost") {
      booster <- fit$finalModel
      
      if (is.null(booster)) {
        stop("XGBoost fit object does not contain finalModel.")
      }
      
      if (is.null(fit$feature_names)) {
        stop(
          "XGBoost fit object does not contain fit$feature_names. ",
          "Cannot safely extract XGBoost feature importance."
        )
      }
      
      # Use feature names captured before .outcome was added to trainingData
      feat_names <- as.character(fit$feature_names)
      
      if (length(feat_names) == 0) {
        stop("fit$feature_names is empty for XGBoost model.")
      }
      
      if (anyDuplicated(feat_names)) {
        dup_names <- unique(feat_names[duplicated(feat_names)])
        
        stop(
          "Duplicate XGBoost feature names detected: ",
          paste(dup_names, collapse = ", ")
        )
      }
      
      if (".outcome" %in% feat_names) {
        stop(
          ".outcome is present in fit$feature_names. ",
          "This indicates true outcome bleed-over into the XGBoost predictor matrix."
        )
      }
      
      bad_feature_names <- feat_names[
        !stringr::str_detect(
          stringr::str_trim(as.character(feat_names)),
          "^(m__|t__)"
        )
      ]
      
      if (length(bad_feature_names) > 0) {
        print(unique(bad_feature_names))
        
        stop(
          "XGBoost feature names contain non-metabolite/non-taxon predictors. ",
          "Expected all predictors to start with m__ or t__."
        )
      }
      
      xi <- tryCatch(
        xgboost::xgb.importance(
          feature_names = feat_names,
          model = booster
        ),
        error = function(e) {
          stop(
            "xgb.importance failed for XGBoost model: ",
            conditionMessage(e)
          )
        }
      )
      
      if (is.null(xi) || nrow(xi) == 0) {
        return(cbind(meta_cols, data.frame(
          Feature = character(),
          Importance = numeric(),
          ImportanceType = character(),
          stringsAsFactors = FALSE
        )))
      }
      
      if (any(xi$Feature == ".outcome")) {
        print(xi[xi$Feature == ".outcome", , drop = FALSE])
        
        stop(
          "xgb.importance returned .outcome as a feature even though fit$feature_names was used. ",
          "This should not happen; inspect feature-name mapping before using this output."
        )
      }
      
      bad_importance_features <- xi$Feature[
        !stringr::str_detect(
          stringr::str_trim(as.character(xi$Feature)),
          "^(m__|t__)"
        )
      ]
      
      if (length(bad_importance_features) > 0) {
        print(unique(bad_importance_features))
        
        stop(
          "xgb.importance returned non-metabolite/non-taxon features. ",
          "Expected all returned XGBoost features to start with m__ or t__."
        )
      }
      
      df <- data.frame(
        Feature = xi$Feature,
        Importance = as.numeric(xi$Gain),
        ImportanceType = "xgb(Gain)",
        Gain = as.numeric(xi$Gain),
        Cover = as.numeric(xi$Cover),
        Frequency = as.numeric(xi$Frequency),
        stringsAsFactors = FALSE
      )
      
      out <- cbind(meta_cols, df)
      out <- out[order(-out$Importance), , drop = FALSE]
      rownames(out) <- NULL
      return(out)
    }
    
    stop("Unknown model_type: ", model_type)
  }
  
  check_keep_present <- function(keep_vec, train_df, val_df, label="") {
    miss_train <- setdiff(keep_vec, names(train_df))
    miss_val   <- setdiff(keep_vec, names(val_df))
    if (length(miss_train) || length(miss_val)) {
      stop(
        label, " keep contains columns not present.\n",
        "Missing in train: ", paste(head(miss_train, 20), collapse=", "), "\n",
        "Missing in val: ", paste(head(miss_val, 20), collapse=", ")
      )
    }
  }
  
  # apply log2 scaling for metabolomics data
  apply_metab_transform <- function(df, metab_cols, transform = "log2", pseudocount = 1e-9) {
    if (length(metab_cols) == 0) return(df)
    if (transform == "none") {
      # still coerce to numeric safely, but preserve names
      df[, metab_cols] <- lapply(df[, metab_cols, drop = FALSE], as.numeric)
      return(df)
    }
    if (transform == "log2") {
      df[, metab_cols] <- lapply(df[, metab_cols, drop = FALSE], function(x) log2(as.numeric(x) + pseudocount))
      return(df)
    }
    stop("Unknown metab_transform: ", transform)
  }
  
  # NZV + correlation filter (robust)
  fit_nzv_corr <- function(train_df, feature_cols, corr_cutoff = 0.90) {
    if (length(feature_cols) == 0) return(character(0))
    X <- train_df[, feature_cols, drop = FALSE]
    na_cols <- sapply(X, function(x) all(is.na(x)))
    if (any(na_cols)) {
      X <- X[, !na_cols, drop = FALSE]
      feature_cols <- colnames(X)
    }
    if (length(feature_cols) == 0) return(character(0))
    
    nzv_idx <- caret::nearZeroVar(X)
    keep_cols <- feature_cols
    if (length(nzv_idx) > 0) keep_cols <- setdiff(feature_cols, feature_cols[nzv_idx])
    
    if (length(keep_cols) >= 2) {
      vars <- apply(train_df[, keep_cols, drop = FALSE], 2, function(x) var(as.numeric(x), na.rm = TRUE))
      keep_cols <- keep_cols[vars > 0 & !is.na(vars)]
    }
    
    if (length(keep_cols) >= 2) {
      Xk <- as.data.frame(lapply(train_df[, keep_cols, drop = FALSE], function(x) as.numeric(x)))
      cmat <- suppressWarnings(cor(Xk, use = "pairwise.complete.obs"))
      if (all(is.na(cmat))) {
        warning("Correlation matrix all NA; skipping correlation filtering.")
      } else {
        bad <- caret::findCorrelation(cmat, cutoff = corr_cutoff, names = TRUE, exact = TRUE)
        if (length(bad) > 0) keep_cols <- setdiff(keep_cols, bad)
      }
    }
    
    if (length(keep_cols) == 0) {
      keep_cols <- feature_cols
      warning("NZV+corr removed all features; returning original feature set as fallback.")
    }
    keep_cols
  }
  
  apply_taxa_transform <- function(df, taxa_cols, transform = "none",
                                   pseudocount = 1e-6,
                                   negative_policy = c("error", "warn_clamp0", "warn_drop_rows")) {
    if (transform == "none" || length(taxa_cols) == 0) return(df)
    negative_policy <- match.arg(negative_policy)
    
    if (transform != "clr") stop("Unknown taxa_transform: ", transform)
    
    X <- as.data.frame(lapply(df[, taxa_cols, drop = FALSE], function(x) as.numeric(x)))
    
    if (anyNA(X)) stop("CLR: NA detected in taxa columns; remove/impute before CLR.")
    
    # Guard: negatives
    if (any(X < 0)) {
      neg_count <- sum(X < 0)
      msg <- paste0("CLR: found ", neg_count, " negative taxa values. Negatives should not be present for CLR.")
      if (negative_policy == "error") stop(msg)
      
      if (negative_policy == "warn_clamp0") {
        warning(msg, " Clamping negatives to 0.")
        X[X < 0] <- 0
      }
      
      if (negative_policy == "warn_drop_rows") {
        warning(msg, " Dropping rows containing any negative taxa value.")
        bad_rows <- apply(X, 1, function(r) any(r < 0))
        df <- df[!bad_rows, , drop = FALSE]
        X  <- X[!bad_rows, , drop = FALSE]
      }
    }
    
    # Drop all-zero rows in taxa block (CLR undefined)
    row_sum <- rowSums(as.matrix(X))
    bad_rows <- which(row_sum <= 0)
    if (length(bad_rows) > 0) {
      warning("CLR: ", length(bad_rows),
              " rows have zero total abundance in taxa block. Dropping these rows.")
      df <- df[-bad_rows, , drop = FALSE]
      X  <- X[-bad_rows, , drop = FALSE]
    }
    
    # Replace zeros with pseudocount
    if (any(X == 0)) {
      message("CLR: replacing zeros with pseudocount = ", pseudocount)
      X[X == 0] <- pseudocount
    }
    
    # CLR transform
    logX <- log(as.matrix(X))
    row_mean <- rowMeans(logX)
    clr_mat <- sweep(logX, 1, row_mean, "-")
    
    # Sanity check
    rs <- rowSums(clr_mat)
    max_abs <- max(abs(rs))
    if (!is.finite(max_abs)) stop("CLR produced non-finite row sums (Inf/NaN).")
    if (max_abs > 1e-4) {
      warning(sprintf("CLR rowSum check large: max |rowSum| = %.3e (nrow=%d). Continuing.",
                      max_abs, nrow(clr_mat)))
      worst <- head(order(abs(rs), decreasing = TRUE), 5)
      message("Worst CLR rowSums: ", paste(sprintf("%.3e", rs[worst]), collapse = ", "))
    }
    
    # Write back + return
    df[, taxa_cols] <- clr_mat
    return(df)
  }
  
  zscore_fit <- function(df, feature_cols) {
    mu  <- sapply(df[, feature_cols, drop = FALSE], function(x) mean(as.numeric(x), na.rm = TRUE))
    sdv <- sapply(df[, feature_cols, drop = FALSE], function(x) sd(as.numeric(x), na.rm = TRUE))
    sdv[is.na(sdv) | sdv == 0] <- 1
    list(mu = mu, sd = sdv)
  }
  
  zscore_apply <- function(df, feature_cols, fit) {
    X <- as.data.frame(lapply(df[, feature_cols, drop = FALSE], as.numeric))
    X <- sweep(X, 2, fit$mu, "-")
    X <- sweep(X, 2, fit$sd, "/")
    df[, feature_cols] <- X
    df
  }
  
  apply_block_weights_inplace <- function(df, metab_cols, taxa_cols) {
    pm <- length(metab_cols)
    pt <- length(taxa_cols)
    if (pm > 0) df[, metab_cols] <- df[, metab_cols, drop = FALSE] * (1 / sqrt(pm))
    if (pt > 0) df[, taxa_cols]  <- df[, taxa_cols,  drop = FALSE] * (1 / sqrt(pt))
    df
  }
  
  apply_concat_balance_variants <- function(train_concat_red, val_concat_red,
                                            response_variable, stratify_variable,
                                            kept_metab_cols, kept_taxa_cols) {
    # predictor columns (exclude outcome + optional stratify)
    drop_cols <- c(response_variable)
    if (!is.null(stratify_variable) && stratify_variable != "") drop_cols <- c(drop_cols, stratify_variable)
    concat_pred_cols <- setdiff(names(train_concat_red), drop_cols)
    
    # blocks (only those present)
    metab_block <- intersect(kept_metab_cols, concat_pred_cols)
    taxa_block  <- intersect(kept_taxa_cols,  concat_pred_cols)
    
    # 1) none
    out_none_train <- train_concat_red
    out_none_val   <- val_concat_red
    
    # 2) scaling (fold-wise z-score on predictors)
    out_scaling_train <- train_concat_red
    out_scaling_val   <- val_concat_red
    if (length(concat_pred_cols) > 0) {
      zfit <- zscore_fit(out_scaling_train, concat_pred_cols)
      out_scaling_train <- zscore_apply(out_scaling_train, concat_pred_cols, zfit)
      out_scaling_val   <- zscore_apply(out_scaling_val,   concat_pred_cols, zfit)
    }
    
    # 3) scaling + weights
    out_sw_train <- out_scaling_train
    out_sw_val   <- out_scaling_val
    out_sw_train <- apply_block_weights_inplace(out_sw_train, metab_block, taxa_block)
    out_sw_val   <- apply_block_weights_inplace(out_sw_val,   metab_block, taxa_block)
    
    list(
      concat_none = list(train = out_none_train, val = out_none_val),
      concat_scaling = list(train = out_scaling_train, val = out_scaling_val),
      concat_scaling_weights = list(train = out_sw_train, val = out_sw_val)
    )
  }
  
  
  make_limma_X <- function(df, feature_cols) {
    if (length(feature_cols) == 0) return(NULL)
    
    M <- as.matrix(df[, feature_cols, drop = FALSE])
    storage.mode(M) <- "double"     # safely coerce to numeric
    X <- t(M)                       # features x samples for limma
    X
  }
  
  # Effect-size helpers for coefficient stability
  # Direction is deterministic: level2 (sorted) minus level1 (sorted)
  compute_limma_effects <- function(train_df, stratify_variable, feature_cols) {
    if (length(feature_cols) == 0) return(data.frame())
    grp_raw <- train_df[[stratify_variable]]
    grp_levels <- sort(unique(as.character(grp_raw)))
    if (length(grp_levels) != 2) return(data.frame())
    grp <- factor(as.character(grp_raw), levels = grp_levels)
    
    X <- make_limma_X(train_df, feature_cols)
    if (is.null(X) || nrow(X) == 0 || ncol(X) == 0) return(data.frame())
    
    design <- model.matrix(~ grp)
    fit <- limma::eBayes(limma::lmFit(X, design))
    tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
    if (nrow(tt) == 0) return(data.frame())
    
    data.frame(
      Feature = rownames(tt),
      Effect = as.numeric(tt$logFC),
      PValue = as.numeric(tt$P.Value),
      AdjP   = as.numeric(tt$adj.P.Val),
      Level1 = grp_levels[1],
      Level2 = grp_levels[2],
      stringsAsFactors = FALSE
    )
  }
  
  compute_wilcox_effects <- function(train_df, stratify_variable, feature_cols) {
    if (length(feature_cols) == 0) return(data.frame())
    grp_raw <- train_df[[stratify_variable]]
    grp_levels <- sort(unique(as.character(grp_raw)))
    if (length(grp_levels) != 2) return(data.frame())
    grp <- factor(as.character(grp_raw), levels = grp_levels)
    
    pvals <- rep(NA_real_, length(feature_cols))
    effs  <- rep(NA_real_, length(feature_cols))
    names(pvals) <- feature_cols
    names(effs)  <- feature_cols
    
    for (col in feature_cols) {
      x <- as.numeric(train_df[[col]])
      g1 <- x[grp == grp_levels[1]]
      g2 <- x[grp == grp_levels[2]]
      
      # degenerate
      if (length(unique(c(g1, g2))) < 2) {
        pvals[col] <- 1
        effs[col] <- 0
        next
      }
      
      effs[col] <- stats::median(g2, na.rm = TRUE) - stats::median(g1, na.rm = TRUE)
      
      p <- tryCatch({
        suppressWarnings(wilcox.test(g1, g2, exact = FALSE)$p.value)
      }, error = function(e) 1)
      if (!is.finite(p)) p <- 1
      pvals[col] <- p
    }
    
    qvals <- p.adjust(pvals, method = "fdr")
    
    data.frame(
      Feature = feature_cols,
      Effect = as.numeric(effs[feature_cols]),
      PValue = as.numeric(pvals[feature_cols]),
      AdjP = as.numeric(qvals[feature_cols]),
      Level1 = grp_levels[1],
      Level2 = grp_levels[2],
      stringsAsFactors = FALSE
    )
  }
  
  ensure_min_keep <- function(keep, ranked, min_keep_p) {
    # keep: current keep vector (may be empty)
    # ranked: vector of all features sorted best->worst (e.g., by p)
    # returns at least min_keep_p if possible
    if (length(ranked) == 0) return(character(0))
    if (length(keep) >= min_keep_p) return(unique(keep))
    # pad by taking from ranked until reach min_keep_p
    padded <- unique(c(keep, ranked))
    padded[seq_len(min(min_keep_p, length(padded)))]
  }
  
  fit_limma_select <- function(train_df, stratify_variable, feature_cols,
                               top_k = NULL, fdr_q = 0.05,
                               fallback_top_k = 200, min_keep_p = 10) {
    
    if (is.null(stratify_variable) || stratify_variable == "")
      stop("Limma reduction requires a non-null stratify_variable.")
    if (!(stratify_variable %in% colnames(train_df)))
      stop("stratify_variable not found in training data.")
    
    grp_raw <- train_df[[stratify_variable]]
    grp_levels <- sort(unique(as.character(grp_raw)))
    if (length(grp_levels) != 2)
      stop("stratify_variable must be binary (exactly 2 levels) for limma.")
    grp <- factor(as.character(grp_raw), levels = grp_levels)
    
    if (length(feature_cols) == 0) return(character(0))
    
    X <- make_limma_X(train_df, feature_cols)
    if (is.null(X) || nrow(X) == 0 || ncol(X) == 0) return(character(0))
    
    design <- model.matrix(~ grp)
    fit <- limma::eBayes(limma::lmFit(X, design))
    tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
    
    ranked <- rownames(tt)
    cat("Missing from train_df:", sum(!ranked %in% names(train_df)), "\n")
    print(head(setdiff(ranked, names(train_df)), 20))
    
    if (nrow(tt) == 0) {
      warning("Limma returned zero rows — check input data. Returning original features.")
      return(feature_cols)
    }
    
    ranked <- rownames(tt)
    
    if (!is.null(top_k)) {
      keep <- ranked[seq_len(min(top_k, length(ranked)))]
    } else {
      keep <- ranked[tt$adj.P.Val <= fdr_q]
      if (length(keep) == 0) {
        message("No features passed FDR in limma. Using fallback_top_k = ", fallback_top_k)
        keep <- ranked[seq_len(min(fallback_top_k, length(ranked)))]
      }
    }
    
    if (length(keep) < min_keep_p) {
      message("Limma selected ", length(keep), " < ", min_keep_p,
              "; padding to ", min_keep_p, " using smallest raw P values.")
      keep <- ensure_min_keep(keep, ranked, min_keep_p)
    }
    
    if (length(keep) == 0) keep <- feature_cols
    keep
  }
  
  fit_wilcox_select <- function(train_df,
                                stratify_variable,
                                feature_cols,
                                top_k = 200,
                                fdr_q = NULL,
                                fallback_top_k = 200,
                                min_keep_p = 10) {
    
    if (is.null(stratify_variable) || stratify_variable == "")
      stop("Wilcoxon reduction requires a non-null stratify_variable.")
    if (!(stratify_variable %in% colnames(train_df)))
      stop("stratify_variable not found in training data.")
    
    grp_raw <- train_df[[stratify_variable]]
    grp_levels <- sort(unique(as.character(grp_raw)))
    if (length(grp_levels) != 2)
      stop("stratify_variable must be binary (exactly 2 levels) for wilcox.")
    grp <- factor(as.character(grp_raw), levels = grp_levels)
    
    if (length(feature_cols) == 0) return(character(0))
    
    pvals <- sapply(feature_cols, function(col) {
      x <- as.numeric(train_df[[col]])
      g1 <- x[grp == grp_levels[1]]
      g2 <- x[grp == grp_levels[2]]
      if (length(unique(c(g1, g2))) < 2) return(1)
      p <- tryCatch(suppressWarnings(wilcox.test(g1, g2, exact = FALSE)$p.value),
                    error = function(e) 1)
      if (!is.finite(p)) p <- 1
      p
    })
    
    ranked <- names(sort(pvals, decreasing = FALSE))  # best -> worst
    
    # Primary selection
    if (!is.null(fdr_q)) {
      qvals <- p.adjust(pvals, method = "fdr")
      keep <- names(qvals)[qvals <= fdr_q]
      if (length(keep) == 0) {
        message("No features passed FDR in wilcox; using fallback_top_k = ", fallback_top_k)
        keep <- ranked[seq_len(min(fallback_top_k, length(ranked)))]
      }
    } else {
      keep <- ranked[seq_len(min(top_k, length(ranked)))]
    }
    
    # Enforce minimum keep
    if (length(keep) < min_keep_p) {
      message("Wilcox selected ", length(keep), " < ", min_keep_p,
              "; padding to ", min_keep_p, " using smallest raw p-values.")
      keep <- ensure_min_keep(keep, ranked, min_keep_p)
    }
    
    if (length(keep) == 0) {
      warning("Wilcox produced 0 features; returning original feature set as fallback.")
      keep <- feature_cols
    }
    
    keep
  }
  
  # apply_view_reduction with soft fallback if stratify missing/not binary
  # NEW: optional return_effects=TRUE to collect coefficient/effect sizes
  apply_view_reduction <- function(train_df, val_df,
                                   feature_cols,
                                   response_variable,
                                   method,
                                   response_type,
                                   stratify_variable = NULL,
                                   corr_cutoff = 0.90,
                                   selector_top_k = 200,
                                   selector_fdr_q = 0.05,
                                   min_keep_p = 10,
                                   return_effects = FALSE,
                                   view_name = NA_character_) {
    
    original_method <- method
    stratify_valid <- TRUE
    effects_df <- data.frame()
    
    if (method %in% c("limma", "wilcox")) {
      if (is.null(stratify_variable) || stratify_variable == "" || !(stratify_variable %in% colnames(train_df))) {
        message("Stratify variable missing for method '", method, "'. Falling back to 'none'.")
        stratify_valid <- FALSE
      } else {
        grp <- as.factor(train_df[[stratify_variable]])
        if (length(unique(as.character(grp))) != 2) {
          message("Stratify variable is not binary for method '", method, "'. Falling back to 'none'.")
          stratify_valid <- FALSE
        }
      }
      if (!stratify_valid) method <- "none"
    }
    
    if (method == "none") {
      keep <- feature_cols
      
    } else if (method == "nzv_corr") {
      keep <- fit_nzv_corr(train_df, feature_cols, corr_cutoff = corr_cutoff)
      
    } else if (method == "limma") {
      # collect full effects table if requested
      if (return_effects) {
        effects_df <- compute_limma_effects(train_df, stratify_variable, feature_cols)
        if (nrow(effects_df) == 0) message("Limma effects table empty (view=", view_name, ").")
      }
      keep <- fit_limma_select(
        train_df = train_df,
        stratify_variable = stratify_variable,
        feature_cols = feature_cols,
        top_k = selector_top_k,
        fdr_q = selector_fdr_q,
        min_keep_p = min_keep_p
      )
      
    } else if (method == "wilcox") {
      if (return_effects) {
        effects_df <- compute_wilcox_effects(train_df, stratify_variable, feature_cols)
        if (nrow(effects_df) == 0) message("Wilcox effects table empty (view=", view_name, ").")
      }
      keep <- fit_wilcox_select(
        train_df = train_df,
        stratify_variable = stratify_variable,
        feature_cols = feature_cols,
        top_k = selector_top_k,
        fdr_q = selector_fdr_q,
        min_keep_p = min_keep_p
      )
      
    } else {
      stop("Unknown reduction method: ", original_method)
    }
    
    if (length(keep) == 0) {
      message("Reduction method '", original_method, "' produced 0 features. Falling back to 'none'.")
      keep <- feature_cols
    }
    
    train_red <- train_df %>% dplyr::select(all_of(response_variable), all_of(keep))
    val_red   <- val_df   %>% dplyr::select(all_of(response_variable), all_of(keep))
    
    # If we gathered effects: subset to kept features and annotate
    if (return_effects && nrow(effects_df) > 0) {
      effects_df <- effects_df[effects_df$Feature %in% keep, , drop = FALSE]
      effects_df$View <- view_name
      effects_df$Method <- method
    }
    
    list(train = train_red, val = val_red, keep = keep, effects = effects_df, method_used = method)
  }
  
  # Build concat automatically (UNION of kept features)
  build_concat_auto <- function(train_data, val_data, response_variable, kept_metab, kept_taxa) {
    keep_feats <- union(kept_metab, kept_taxa)
    if (length(keep_feats) == 0) {
      keep_feats <- setdiff(names(train_data), response_variable)
      warning("Concat union empty -> using all features from concat as fallback.")
    }
    keep_cols <- c(response_variable, keep_feats)
    train_concat_red <- train_data %>% dplyr::select(all_of(keep_cols))
    val_concat_red   <- val_data   %>% dplyr::select(all_of(keep_cols))
    list(train = train_concat_red, val = val_concat_red, keep = keep_feats)
  }
  
  # ---------------------------
  # NA handling (row-wise removal on critical columns)
  # ---------------------------
  all_columns_to_check <- c(response_variable, m_columns, t_columns)
  if (!is.null(stratify_variable) && stratify_variable != "") {
    all_columns_to_check <- unique(c(all_columns_to_check, stratify_variable))
  }
  missing_check_cols <- setdiff(all_columns_to_check, colnames(data))
  if (length(missing_check_cols) > 0) {
    warning("Some check columns are not present in data and will be ignored for NA removal: ",
            paste(missing_check_cols, collapse = ", "))
    all_columns_to_check <- setdiff(all_columns_to_check, missing_check_cols)
  }
  
  if (length(all_columns_to_check) > 0) {
    rows_with_na <- apply(data[, all_columns_to_check, drop = FALSE], 1, function(r) any(is.na(r)))
    n_na_rows <- sum(rows_with_na)
    if (n_na_rows > 0) {
      message(sprintf("Removing %d rows with NA in critical columns (response/metab/taxa/stratify).", n_na_rows))
      data <- data[!rows_with_na, , drop = FALSE]
    } else {
      message("No rows with NA detected in critical columns.")
    }
  } else {
    message("No columns available for NA-check; skipping NA removal.")
  }
  
  convert_and_check_numeric(c(m_columns, t_columns))
  print('All metabolite and taxa columns are numeric.')
  
  # ---------------------------
  # Helper: infer positive class (case-insensitive)
  # Prefers Control/Healthy if present; otherwise falls back to 2nd level.
  # ---------------------------
  infer_positive_class <- function(y_factor,
                                   preferred = c("Control", "Healthy"),
                                   synonyms = c("control", "healthy", "normal")) {
    if (!is.factor(y_factor)) y_factor <- as.factor(y_factor)
    levs <- levels(y_factor)
    levs_lc <- tolower(levs)
    
    preferred_lc <- tolower(preferred)
    synonyms_lc  <- tolower(synonyms)
    
    # exact match to preferred
    hit <- which(levs_lc %in% preferred_lc)
    if (length(hit) >= 1) return(levs[hit[1]])
    
    # exact match to synonyms
    hit <- which(levs_lc %in% synonyms_lc)
    if (length(hit) >= 1) return(levs[hit[1]])
    
    # partial match (e.g., "Control_Group", "Healthy_Subjects")
    pat <- paste0("(", paste(c(preferred_lc, synonyms_lc), collapse="|"), ")")
    hit <- which(grepl(pat, levs_lc))
    if (length(hit) >= 1) return(levs[hit[1]])
    
    # fallback: old convention
    levs[2]
  }
  
  # Metrics and response variable type
  if (type_of_analysis == 'binary') {
    source(file.path(file_path, 'A_Integrative_Pipeline_Scripts/model_functions/binary_functions.R'))
    metrics <- list(
      accuracy = accuracy_calculation,
      kappa    = kappa_calculation,
      auroc    = auroc_calculation,
      auprc    = auprc_calculation,
      sensitivity = sensitivity_calculation,
      specificity = specificity_calculation,
      precision = precision_calculation,
      f1 = f1_calculation,
      balanced_accuracy = balanced_accuracy_calculation
    )
    print('For binary table this is the counts recorded:')
    print(table(data[[response_variable]]))
    print('Binary functions loaded.')
    
    if (!is.factor(data[[response_variable]])) {
      data[[response_variable]] <- as.factor(data[[response_variable]])
      print(paste('The response variable', response_variable, 'has been converted to a factor.'))
    }
    
    if (length(levels(data[[response_variable]])) != 2) {
      stop("Binary response must have exactly 2 levels.")
    }
    
    # Enforce deterministic class order (important for ROC and truth=1 mapping)
    data[[response_variable]] <- factor(
      data[[response_variable]],
      levels = sort(levels(data[[response_variable]]))
    )
    
    print(paste("Binary levels set to:", paste(levels(data[[response_variable]]), collapse = ", ")))
    # Auto-detect positive class (case-insensitive): prefers Control/Healthy
    args$positive_class <- infer_positive_class(data[[response_variable]])
    print(paste("Auto positive_class selected as:", args$positive_class))
    
  } else if (type_of_analysis == 'continuous') {
    source(file.path(file_path, 'A_Integrative_Pipeline_Scripts/model_functions/continuous_functions.R'))
    metrics <- list(
      rmse = rmse_calculation,
      r2   = r2_calculation,
      mae  = mae_calculation,
      mape = mape_calculation,
      median_abs_error = median_absolute_error_calculation
    )
    print('Continuous functions loaded.')
    
    if (!is.numeric(data[[response_variable]])) {
      data[[response_variable]] <- as.numeric(data[[response_variable]])
      print(paste('The response variable', response_variable, 'has been converted to numeric.'))
    }
    
  } else {
    stop("Invalid type_of_analysis. Please enter 'binary' or 'continuous'.")
  }
  

  
  # ---------------------------
  # SECTION 3 - OOF STACKING + Progress/Timing/Debug
  # ---------------------------
  if (type_of_analysis == "binary") {
    if (!is.factor(data[[response_variable]])) {
      if (all(data[[response_variable]] %in% c(0, 1))) {
        data[[response_variable]] <- factor(data[[response_variable]],
                                            levels = c(0, 1),
                                            labels = c("Class0", "Class1"))
      } else {
        data[[response_variable]] <- factor(data[[response_variable]])
        if (length(levels(data[[response_variable]])) != 2) stop("Binary response must have exactly 2 levels.")
      }
    } else {
      if (length(levels(data[[response_variable]])) != 2) stop("Binary response must have exactly 2 levels.")
    }
  } else {
    data[[response_variable]] <- as.numeric(data[[response_variable]])
  }
  
  # 1) Split train and test
  if (!is.null(stratify_variable) && stratify_variable != "" && stratify_variable %in% names(data)) {
    parts <- createDataPartition(data[[stratify_variable]], p = training_proportion, list = FALSE)
  } else {
    parts <- createDataPartition(data[[response_variable]], p = training_proportion, list = FALSE)
  }
  train <- data[parts, , drop = FALSE]
  test  <- data[-parts, , drop = FALSE]
  
  # Subset the testing data
  testset_concatenation <- test %>% dplyr::select(all_of(response_variable), all_of(m_columns), all_of(t_columns))
  testset_metabolomics  <- test %>% dplyr::select(all_of(response_variable), all_of(m_columns))
  testset_mss           <- test %>% dplyr::select(all_of(response_variable), all_of(t_columns))
  
  # 2) Prepare OOF storage aligned to rows of 'train'
  n_train <- nrow(train)
  oof_metab_sum   <- rep(0, n_train)
  oof_mss_sum     <- rep(0, n_train)
  oof_concat_none_sum   <- rep(0, n_train)
  oof_concat_scaling_sum <- rep(0, n_train)
  oof_concat_sw_sum      <- rep(0, n_train)  # scaling+weights
  oof_count       <- rep(0, n_train)
  
  train_truth <- if (type_of_analysis == "binary") {
    as.numeric(train[[response_variable]] == args$positive_class)
  } else {
    as.numeric(train[[response_variable]])
  }
  
  metab_model_folds  <- list()
  mss_model_folds    <- list()
  concat_none_model_folds <- list()
  concat_scaling_model_folds <- list()
  concat_sw_model_folds <- list()
  
  # Debug/progress controls
  verbose_progress <- TRUE
  save_oof_debug   <- TRUE
  save_fold_index  <- TRUE
  progress_every_fold <- TRUE
  
  t0_all <- Sys.time()
  fold_timing <- data.frame(
    repeat_id = integer(),
    fold = integer(),
    seconds = numeric(),
    n_train_fold = integer(),
    n_val_fold = integer(),
    stringsAsFactors = FALSE
  )
  
  oof_debug <- list()
  fold_indices_by_repeat <- vector("list", length = num_repeats)
  .elapsed_sec <- function(t_start) as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  
  model_categories <- c("metabolomics", "mss", "concat_none", "concat_scaling", "concat_scaling_weights")
  fold_metric_store <- lapply(model_categories, function(.) {
    as.list(setNames(vector("list", length(metrics)), names(metrics)))
  })
  names(fold_metric_store) <- model_categories
  for (cat in names(fold_metric_store)) {
    for (m in names(metrics)) fold_metric_store[[cat]][[m]] <- numeric(0)
  }
  
  # ---------------------------
  # per-fold feature count logging
  # ---------------------------
  fold_feature_log <- data.frame(
    run_id = integer(),
    seed = integer(),
    repeat_id = integer(),
    fold = integer(),
    metab_method = character(),
    taxa_method = character(),
    n_metab_keep = integer(),
    n_taxa_keep = integer(),
    n_concat_keep = integer(),
    stringsAsFactors = FALSE
  )
  
  fold_metrics_raw <- data.frame(
    run_id = integer(),
    seed = integer(),
    repeat_id = integer(),
    fold = integer(),
    model = character(),
    metric = character(),
    value = numeric(),
    stringsAsFactors = FALSE
  )
  #store fold importance
  fold_importance_all <- list()
  # 3) Repeated CV loop - produce OOF predictions (NO test predictions here)
  set.seed(1234)
  for (repeat_ in 1:num_repeats) {
    t0_rep <- Sys.time()
    if (verbose_progress) message(sprintf("Repeat %d/%d started...", repeat_, num_repeats))
    
    fold_y <- if (!is.null(stratify_variable) && stratify_variable != "" && stratify_variable %in% names(train)) {
      train[[stratify_variable]]
    } else {
      train[[response_variable]]
    }
    
    set.seed(1234 + repeat_)
    indices <- createFolds(fold_y, k = num_folds, list = TRUE)
    if (save_fold_index) fold_indices_by_repeat[[repeat_]] <- indices
    
    for (fold in seq_along(indices)) {
      t0_fold <- Sys.time()
      if (verbose_progress && progress_every_fold) {
        message(sprintf("  Fold %d/%d (repeat %d/%d) ...", fold, num_folds, repeat_, num_repeats))
      }
      
      train_index      <- unlist(indices[-fold])
      validation_index <- unlist(indices[fold])
      
      train_data      <- train[train_index, , drop = FALSE]
      validation_data <- train[validation_index, , drop = FALSE]
      
      # helper: columns we want to carry along (NOT predictors)
      extra_cols <- c(response_variable)
      if (!is.null(stratify_variable) && stratify_variable != "" && stratify_variable %in% names(train_data)) {
        extra_cols <- c(extra_cols, stratify_variable)
      }
      
      # Subset per view (carry stratify column through)
      train_concat <- train_data %>% dplyr::select(all_of(extra_cols), all_of(m_columns), all_of(t_columns))
      train_metab  <- train_data %>% dplyr::select(all_of(extra_cols), all_of(m_columns))
      train_mss    <- train_data %>% dplyr::select(all_of(extra_cols), all_of(t_columns))
      
      val_concat <- validation_data %>% dplyr::select(all_of(extra_cols), all_of(m_columns), all_of(t_columns))
      val_metab  <- validation_data %>% dplyr::select(all_of(extra_cols), all_of(m_columns))
      val_mss    <- validation_data %>% dplyr::select(all_of(extra_cols), all_of(t_columns))
      
      
      # Determine pseudocount to use for taxa CLR (prefer taxa_pseudocount; fallback to metab pseudocount)
      taxa_pc <- args$pseudocount
      
      # Apply metabolomics transforms (existing)
      train_metab <- apply_metab_transform(train_metab, m_columns, transform = args$metab_transform, pseudocount = args$pseudocount)
      val_metab   <- apply_metab_transform(val_metab,   m_columns, transform = args$metab_transform, pseudocount = args$pseudocount)
      train_concat <- apply_metab_transform(train_concat, m_columns, transform = args$metab_transform, pseudocount = args$pseudocount)
      val_concat   <- apply_metab_transform(val_concat,   m_columns, transform = args$metab_transform, pseudocount = args$pseudocount)
      
      summarize_taxa_sign <- function(df, taxa_cols, label = "taxa") {
        if (length(taxa_cols) == 0) return(invisible(NULL))
        X <- as.matrix(df[, taxa_cols, drop = FALSE])
        message(label, ": negatives=", sum(X < 0, na.rm = TRUE),
                " zeros=", sum(X == 0, na.rm = TRUE),
                " min=", min(X, na.rm = TRUE))
      }
      
      summarize_taxa_sign(train_mss, t_columns, "train_mss pre-CLR")
      summarize_taxa_sign(val_mss,   t_columns, "val_mss pre-CLR")
      
      # --- NEW: apply taxa transform (CLR) if requested ---
      if (!is.null(args$taxa_transform) && tolower(args$taxa_transform) == "clr") {
        train_mss <- apply_taxa_transform(train_mss, t_columns, transform = "clr", pseudocount = taxa_pc, negative_policy = "error")
        val_mss   <- apply_taxa_transform(val_mss,   t_columns, transform = "clr", pseudocount = taxa_pc, negative_policy = "error")
        # also apply to the concat frames so concat features receive the same transform
        train_concat <- apply_taxa_transform(train_concat, t_columns, transform = "clr", pseudocount = taxa_pc, negative_policy = "error")
        val_concat   <- apply_taxa_transform(val_concat,   t_columns, transform = "clr", pseudocount = taxa_pc, negative_policy = "error")
        message(sprintf("Applied CLR taxa transform (pseudocount=%g) for repeat=%d fold=%d", taxa_pc, repeat_, fold))
      } else {
        # leave taxa as-is
      }
      
      
      
      metab_red <- apply_view_reduction(
        train_df = train_metab, val_df = val_metab,
        feature_cols = m_columns,
        response_variable = response_variable,
        method = args$metab_reduction,
        response_type = type_of_analysis,
        stratify_variable = stratify_variable,
        corr_cutoff = args$corr_cutoff,
        selector_top_k = 200,
        selector_fdr_q = 0.05,
        return_effects = FALSE,
        view_name = "metabolomics"
      )
      message("[LIMMA CHECK] method_used=", metab_red$method_used,
              " keep=", length(metab_red$keep))
      train_metab_red <- metab_red$train
      val_metab_red   <- metab_red$val
      kept_metab_cols <- metab_red$keep
      check_keep_present(kept_metab_cols, train_metab, val_metab, label="metabolomics")
      
      taxa_red <- apply_view_reduction(
        train_df = train_mss, val_df = val_mss,
        feature_cols = t_columns,
        response_variable = response_variable,
        method = args$taxa_reduction,
        response_type = type_of_analysis,
        stratify_variable = stratify_variable,
        corr_cutoff = args$corr_cutoff,
        selector_top_k = 200,
        selector_fdr_q = 0.05,
        return_effects = FALSE,
        view_name = "taxa"
      )
      message("[TAXA CHECK] method_used=", taxa_red$method_used,
              " keep=", length(taxa_red$keep))
      train_mss_red <- taxa_red$train
      val_mss_red   <- taxa_red$val
      kept_taxa_cols <- taxa_red$keep
      check_keep_present(kept_taxa_cols, train_mss, val_mss, label="taxa")
      
      
      
      # ---- DEBUG: print reduction outcomes per fold ----
      message(sprintf(
        "[REDUCTION DEBUG] run=%d seed=%d repeat=%d fold=%d | metab_method=%s kept=%d | taxa_method=%s kept=%d",
        run_id, seed, repeat_, fold,
        as.character(args$metab_reduction), length(kept_metab_cols),
        as.character(args$taxa_reduction),  length(kept_taxa_cols)
      ))
      
      # Optional: show first few feature names when tiny (useful for diagnosing n=0/1)
      if (length(kept_metab_cols) < 5) {
        message("[REDUCTION DEBUG] kept_metab_cols: ", paste(kept_metab_cols, collapse = ", "))
      }
      if (length(kept_taxa_cols) < 5) {
        message("[REDUCTION DEBUG] kept_taxa_cols: ", paste(kept_taxa_cols, collapse = ", "))
      }
      
      # Optional: hard warning when this will break glmnet
      if (length(kept_metab_cols) < 2) {
        warning(sprintf("[REDUCTION WARNING] metab kept <2 features (kept=%d). glmnet may fail.", length(kept_metab_cols)))
      }
      if (length(kept_taxa_cols) < 2) {
        warning(sprintf("[REDUCTION WARNING] taxa kept <2 features (kept=%d). glmnet may fail.", length(kept_taxa_cols)))
      }
      
      
      
      concat_red <- build_concat_auto(
        train_data = train_concat,
        val_data   = val_concat,
        response_variable = response_variable,
        kept_metab = kept_metab_cols,
        kept_taxa  = kept_taxa_cols
      )
      train_concat_red <- concat_red$train
      val_concat_red   <- concat_red$val
      
      concat_variants <- apply_concat_balance_variants(
        train_concat_red, val_concat_red,
        response_variable = response_variable,
        stratify_variable = stratify_variable,
        kept_metab_cols = kept_metab_cols,
        kept_taxa_cols  = kept_taxa_cols
      )
      
      train_concat_none <- concat_variants$concat_none$train
      val_concat_none   <- concat_variants$concat_none$val
      
      train_concat_scaling <- concat_variants$concat_scaling$train
      val_concat_scaling   <- concat_variants$concat_scaling$val
      
      train_concat_sw <- concat_variants$concat_scaling_weights$train
      val_concat_sw   <- concat_variants$concat_scaling_weights$val

      
      
      # ---------------------------
      # NEW: log per-fold keep counts
      # ---------------------------
      fold_feature_log <- rbind(
        fold_feature_log,
        data.frame(
          run_id = run_id,
          seed = seed,
          repeat_id = repeat_,
          fold = fold,
          metab_method = as.character(args$metab_reduction),
          taxa_method = as.character(args$taxa_reduction),
          n_metab_keep = length(kept_metab_cols),
          n_taxa_keep = length(kept_taxa_cols),
          n_concat_keep = ncol(train_concat_red) - 1L,
          stringsAsFactors = FALSE
        )
      )
      balance_label <- if (!is.null(args$concat_balance)) as.character(args$concat_balance) else "NA"
      
      message(sprintf(
        "Feature keep counts | repeat=%d fold=%d | balance=%s | metab(%s)=%d taxa(%s)=%d concat=%d",
        repeat_, fold,
        balance_label,
        as.character(args$metab_reduction), length(kept_metab_cols),
        as.character(args$taxa_reduction), length(kept_taxa_cols),
        ncol(train_concat_red) - 1L
      ))
      
      # Fit base models on fold-train only (exclude stratify_variable from predictors)
      fit_metab  <- model_env$fit(
        trainData = train_metab_red,
        response_variable = response_variable,
        response_type = type_of_analysis,
        stratify_variable = stratify_variable,
        seed = 100 + repeat_ * 1000 + fold,
        positive_class = args$positive_class,
        use_gpu = args$use_gpu,
        gpu_device = args$gpu_device,
        nthread = args$nthread
      )
      
      fit_mss <- model_env$fit(
        trainData = train_mss_red,
        response_variable = response_variable,
        response_type = type_of_analysis,
        stratify_variable = stratify_variable,
        seed = 200 + repeat_ * 1000 + fold,
        positive_class = args$positive_class,
        use_gpu = args$use_gpu,
        gpu_device = args$gpu_device,
        nthread = args$nthread
      )
      
      fit_concat_none <- model_env$fit(
        trainData = train_concat_none,
        response_variable = response_variable,
        response_type = type_of_analysis,
        stratify_variable = stratify_variable,
        seed = 300 + repeat_ * 1000 + fold + 1,
        positive_class = args$positive_class,
        use_gpu = args$use_gpu,
        gpu_device = args$gpu_device,
        nthread = args$nthread
        )
      
      fit_concat_scaling <- model_env$fit(
        trainData = train_concat_scaling,
        response_variable = response_variable,
        response_type = type_of_analysis,
        stratify_variable = stratify_variable,
        seed = 300 + repeat_ * 1000 + fold + 2,
        positive_class = args$positive_class,
        use_gpu = args$use_gpu,
        gpu_device = args$gpu_device,
        nthread = args$nthread
      )
      
      fit_concat_sw <- model_env$fit(
        trainData = train_concat_sw,
        response_variable = response_variable,
        response_type = type_of_analysis,
        stratify_variable = stratify_variable,
        seed = 300 + repeat_ * 1000 + fold + 3,
        positive_class = args$positive_class,
        use_gpu = args$use_gpu,
        gpu_device = args$gpu_device,
        nthread = args$nthread
      )
      
      # export per-fold feature importance (validation)
      fold_importance_all[[length(fold_importance_all) + 1]] <-
        extract_fold_importance(fit_metab, model_type, "metabolomics",
                                run_id, seed, repeat_, fold, split = "validation")
      
      fold_importance_all[[length(fold_importance_all) + 1]] <-
        extract_fold_importance(fit_mss, model_type, "mss",
                                run_id, seed, repeat_, fold, split = "validation")
      
      fold_importance_all[[length(fold_importance_all) + 1]] <-
        extract_fold_importance(fit_concat_none, model_type, "concat_none",
                                run_id, seed, repeat_, fold, split = "validation")
      
      fold_importance_all[[length(fold_importance_all) + 1]] <-
        extract_fold_importance(fit_concat_scaling, model_type, "concat_scaling",
                                run_id, seed, repeat_, fold, split = "validation")
      
      fold_importance_all[[length(fold_importance_all) + 1]] <-
        extract_fold_importance(fit_concat_sw, model_type, "concat_scaling_weights",
                                run_id, seed, repeat_, fold, split = "validation")
      
      
      metab_model_folds[[length(metab_model_folds) + 1]]   <- fit_metab
      mss_model_folds[[length(mss_model_folds) + 1]]       <- fit_mss
      concat_none_model_folds[[length(concat_none_model_folds) + 1]] <- fit_concat_none
      concat_scaling_model_folds[[length(concat_scaling_model_folds) + 1]] <- fit_concat_scaling
      concat_sw_model_folds[[length(concat_sw_model_folds) + 1]] <- fit_concat_sw
      
      message("DEBUG fit class: ", paste(class(fit_metab), collapse = ", "))
      if (is.list(fit_metab)) message("DEBUG fit has names: ", paste(names(fit_metab), collapse = ", "))
      
      # Predict only on the validation fold (OOF) (exclude stratify_variable from predictors)
      pred_metab <- model_env$predict(
        fit = fit_metab,
        newData = val_metab_red,
        response_variable = response_variable,
        response_type = type_of_analysis,
        stratify_variable = stratify_variable,
        positive_class = args$positive_class
      )
      
      pred_mss <- model_env$predict(
        fit = fit_mss,
        newData = val_mss_red,
        response_variable = response_variable,
        response_type = type_of_analysis,
        stratify_variable = stratify_variable,
        positive_class = args$positive_class
      )
      
      pred_concat_none <- model_env$predict(
        fit = fit_concat_none,
        newData = val_concat_none,
        response_variable = response_variable,
        response_type = type_of_analysis,
        stratify_variable = stratify_variable,
        positive_class = args$positive_class
      )
      
      pred_concat_scaling <- model_env$predict(
        fit = fit_concat_scaling,
        newData = val_concat_scaling,
        response_variable = response_variable,
        response_type = type_of_analysis,
        stratify_variable = stratify_variable,
        positive_class = args$positive_class
      )
      
      pred_concat_sw <- model_env$predict(
        fit = fit_concat_sw,
        newData = val_concat_sw,
        response_variable = response_variable,
        response_type = type_of_analysis,
        stratify_variable = stratify_variable,
        positive_class = args$positive_class
      )
      
      oof_metab_sum[validation_index]  <- oof_metab_sum[validation_index]  + pred_metab$pred
      oof_mss_sum[validation_index]    <- oof_mss_sum[validation_index]    + pred_mss$pred
      oof_concat_none_sum[validation_index]    <- oof_concat_none_sum[validation_index]    + pred_concat_none$pred
      oof_concat_scaling_sum[validation_index] <- oof_concat_scaling_sum[validation_index] + pred_concat_scaling$pred
      oof_concat_sw_sum[validation_index]      <- oof_concat_sw_sum[validation_index]      + pred_concat_sw$pred
      oof_count[validation_index]      <- oof_count[validation_index] + 1
      
      for (metric_name in names(metrics)) {
        fn <- metrics[[metric_name]]
        
        val_metab <- fn(pred_metab$truth, pred_metab$pred)
        val_mss   <- fn(pred_mss$truth, pred_mss$pred)
        val_cn    <- fn(pred_concat_none$truth, pred_concat_none$pred)
        val_cs    <- fn(pred_concat_scaling$truth, pred_concat_scaling$pred)
        val_csw   <- fn(pred_concat_sw$truth, pred_concat_sw$pred)
        
        # Store in aggregated object (existing behavior)
        fold_metric_store[["metabolomics"]][[metric_name]] <- c(
          fold_metric_store[["metabolomics"]][[metric_name]], val_metab
        )
        fold_metric_store[["mss"]][[metric_name]] <- c(
          fold_metric_store[["mss"]][[metric_name]], val_mss
        )
        fold_metric_store[["concat_none"]][[metric_name]] <- c(
          fold_metric_store[["concat_none"]][[metric_name]], val_cn
        )
        fold_metric_store[["concat_scaling"]][[metric_name]] <- c(
          fold_metric_store[["concat_scaling"]][[metric_name]], val_cs
        )
        fold_metric_store[["concat_scaling_weights"]][[metric_name]] <- c(
          fold_metric_store[["concat_scaling_weights"]][[metric_name]], val_csw
        )
        
        # Store raw fold metrics
        fold_metrics_raw <- rbind(
          fold_metrics_raw,
          data.frame(run_id=run_id, seed=seed, repeat_id=repeat_, fold=fold,
                     model="metabolomics", metric=metric_name, value=val_metab),
          data.frame(run_id=run_id, seed=seed, repeat_id=repeat_, fold=fold,
                     model="mss", metric=metric_name, value=val_mss),
          data.frame(run_id=run_id, seed=seed, repeat_id=repeat_, fold=fold,
                     model="concat_none", metric=metric_name, value=val_cn),
          data.frame(run_id=run_id, seed=seed, repeat_id=repeat_, fold=fold,
                     model="concat_scaling", metric=metric_name, value=val_cs),
          data.frame(run_id=run_id, seed=seed, repeat_id=repeat_, fold=fold,
                     model="concat_scaling_weights", metric=metric_name, value=val_csw)
        )
      }
      
      if (save_oof_debug) {
        oof_debug[[length(oof_debug) + 1]] <- list(
          repeat_id = repeat_,
          fold = fold,
          validation_index = validation_index,
          truth = pred_metab$truth,
          pred_metab = pred_metab$pred,
          pred_mss = pred_mss$pred,
          pred_concat_none = pred_concat_none$pred, 
          pred_concat_sw = pred_concat_sw$pred, 
          pred_concat_scaling = pred_concat_scaling$pred 
          
        )
      }
      
      fold_timing <- rbind(
        fold_timing,
        data.frame(
          repeat_id = repeat_,
          fold = fold,
          seconds = .elapsed_sec(t0_fold),
          n_train_fold = length(train_index),
          n_val_fold = length(validation_index),
          stringsAsFactors = FALSE
        )
      )
      
      if (verbose_progress && progress_every_fold) {
        message(sprintf("    done in %.2fs", tail(fold_timing$seconds, 1)))
      }
    } # end fold
    
    if (verbose_progress) {
      message(sprintf("Repeat %d/%d finished in %.2fs", repeat_, num_repeats, .elapsed_sec(t0_rep)))
    }
  } # end repeat
  
  # write validation-fold importance to CSV
  imp_dir <- file.path(results_base_dir, "FeatureImportanceLogs")
  dir.create(imp_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (length(fold_importance_all) > 0) {
    val_imp_df <- do.call(rbind, fold_importance_all)
    val_imp_file <- file.path(imp_dir, paste0(study_name, "_validation_fold_feature_importance.csv"))
    write.csv(val_imp_df, val_imp_file, row.names = FALSE)
    message("Saved fold feature importance to: ", val_imp_file)
  } else {
    message("No fold feature importance collected.")
  }
  
  perf_dir <- file.path(results_base_dir, "FoldLevelPerformance")
  dir.create(perf_dir, recursive = TRUE, showWarnings = FALSE)
  
  fold_perf_file <- file.path(perf_dir, paste0(study_name, "_fold_level_performance.csv"))
  write.csv(fold_metrics_raw, fold_perf_file, row.names = FALSE)
  
  message("Saved fold-level performance to: ", fold_perf_file)
  
  repeat_summary <- fold_metrics_raw %>%
    dplyr::group_by(run_id, seed, repeat_id, model, metric) %>%
    dplyr::summarise(
      mean_value = mean(value, na.rm=TRUE),
      sd_value   = sd(value, na.rm=TRUE),
      .groups="drop"
    )
  
  repeat_perf_file <- file.path(perf_dir, paste0(study_name, "_repeat_level_performance.csv"))
  write.csv(repeat_summary, repeat_perf_file, row.names = FALSE)
  
  message("Saved repeat-level performance to: ", repeat_perf_file)
  
  if (any(oof_count == 0)) stop("Some training rows never received OOF predictions. Check fold creation.")
  oof_metab  <- oof_metab_sum  / oof_count
  oof_mss    <- oof_mss_sum    / oof_count
  oof_concat_none    <- oof_concat_none_sum    / oof_count
  oof_concat_scaling <- oof_concat_scaling_sum / oof_count
  oof_concat_sw      <- oof_concat_sw_sum      / oof_count
  
  if (verbose_progress) {
    message(sprintf("All repeats/folds complete. Total CV time: %.2fs", .elapsed_sec(t0_all)))
    message(sprintf("Average fold time: %.2fs (min %.2fs, max %.2fs)",
                    mean(fold_timing$seconds),
                    min(fold_timing$seconds),
                    max(fold_timing$seconds)))
  }
  
  meta_train <- data.frame(
    Metabolomics = oof_metab,
    MSS          = oof_mss,
    Concat_None  = oof_concat_none,
    Concat_Scaling = oof_concat_scaling,
    Concat_Scaling_Weights = oof_concat_sw,
    GroundTruth  = train_truth
  )
  
  # ---------------------------
  # Save per-fold feature count log for this run
  # ---------------------------
  foldlog_dir <- file.path(results_base_dir, "FeatureCountLogs")
  dir.create(foldlog_dir, recursive = TRUE, showWarnings = FALSE)
  foldlog_file <- file.path(foldlog_dir, paste0(study_name, "_fold_feature_counts.csv"))
  write.csv(fold_feature_log, foldlog_file, row.names = FALSE)
  message("Per-fold feature count log saved to: ", foldlog_file)
  
  # ---------------------------
  # SECTION 4 - META-LEARNER & FINAL TEST EVALUATION
  # ---------------------------
  if (type_of_analysis == "binary") family <- "binomial" else family <- "gaussian"
  
  meta_train_stack <- data.frame(
    Metabolomics = meta_train$Metabolomics,
    MSS          = meta_train$MSS,
    GroundTruth  = meta_train$GroundTruth
  )
  
  nnls_fit <- nnls::nnls(as.matrix(meta_train_stack[, c("Metabolomics", "MSS")]), meta_train_stack$GroundTruth)
  nnls_weights <- nnls_fit$x
  meta_train_stack$nnls_weighted <- meta_train_stack$Metabolomics * nnls_weights[1] + meta_train_stack$MSS * nnls_weights[2]
  meta_train_stack$average <- 0.5 * meta_train_stack$Metabolomics + 0.5 * meta_train_stack$MSS
  
  alpha_values  <- seq(0, 1, by = 0.1)
  lambda_values <- 10^seq(2, -2, by = -0.1)
  best_alpha <- NULL; best_lambda <- NULL; best_cv_error <- Inf
  for (alpha in alpha_values) {
    cv_fit <- glmnet::cv.glmnet(
      x = as.matrix(meta_train_stack[, c("Metabolomics", "MSS")]),
      y = meta_train_stack$GroundTruth,
      alpha = alpha,
      nfolds = 5,
      family = family,
      lambda = lambda_values,
      standardize = TRUE
    )
    cv_error <- min(cv_fit$cvm)
    if (cv_error < best_cv_error) {
      best_cv_error <- cv_error
      best_alpha    <- alpha
      best_lambda   <- cv_fit$lambda.min
    }
  }
  
  final_glmnet <- glmnet::glmnet(
    x = as.matrix(meta_train_stack[, c("Metabolomics", "MSS")]),
    y = meta_train_stack$GroundTruth,
    alpha = best_alpha,
    lambda = best_lambda,
    family = family,
    standardize = TRUE
  )
  sparse_nnls_weights <- as.numeric(coef(final_glmnet))
  meta_train_stack$sparse_weighted <- as.numeric(predict(final_glmnet, newx = as.matrix(meta_train_stack[, c("Metabolomics", "MSS")]), s = best_lambda))
  
  pls_cv <- pls::plsr(GroundTruth ~ Metabolomics + MSS, data = meta_train_stack, validation = "CV")
  best_ncomp <- which.min(pls_cv$validation$PRESS)
  pls_model <- pls::plsr(GroundTruth ~ Metabolomics + MSS, data = meta_train_stack, ncomp = best_ncomp)
  meta_train_stack$pls_pred <- drop(as.data.frame(predict(pls_model, meta_train_stack[, c("Metabolomics", "MSS")], ncomp = best_ncomp))[, 1])
  
  # --- Refit base models ONCE on full train (Metab, MSS, Concat separately) ---
  extra_cols_full <- c(response_variable)
  if (!is.null(stratify_variable) && stratify_variable != "" && stratify_variable %in% names(train)) {
    extra_cols_full <- c(extra_cols_full, stratify_variable)
  }
  
  train_metab_full  <- train %>% dplyr::select(all_of(extra_cols_full), all_of(m_columns))
  test_metab_full   <- test  %>% dplyr::select(all_of(extra_cols_full), all_of(m_columns))
  
  train_mss_full    <- train %>% dplyr::select(all_of(extra_cols_full), all_of(t_columns))
  test_mss_full     <- test  %>% dplyr::select(all_of(extra_cols_full), all_of(t_columns))
  
  train_concat_full <- train %>% dplyr::select(all_of(extra_cols_full), all_of(m_columns), all_of(t_columns))
  test_concat_full  <- test  %>% dplyr::select(all_of(extra_cols_full), all_of(m_columns), all_of(t_columns))
  
  # choose taxa pseudocount
  taxa_pc <- args$pseudocount
  
  # apply metab transform (already present)
  train_metab_full <- apply_metab_transform(train_metab_full, m_columns, args$metab_transform, args$pseudocount)
  test_metab_full  <- apply_metab_transform(test_metab_full,  m_columns, args$metab_transform, args$pseudocount)
  
  train_concat_full <- apply_metab_transform(train_concat_full, m_columns, args$metab_transform, args$pseudocount)
  test_concat_full  <- apply_metab_transform(test_concat_full,  m_columns, args$metab_transform, args$pseudocount)
  
  # --- NEW: taxa CLR if requested ---
  if (!is.null(args$taxa_transform) && tolower(args$taxa_transform) == "clr") {
    train_mss_full <- apply_taxa_transform(train_mss_full, t_columns, transform = "clr", pseudocount = taxa_pc, negative_policy = "error")
    test_mss_full  <- apply_taxa_transform(test_mss_full,  t_columns, transform = "clr", pseudocount = taxa_pc, negative_policy = "error")
    train_concat_full <- apply_taxa_transform(train_concat_full, t_columns, transform = "clr", pseudocount = taxa_pc, negative_policy = "error")
    test_concat_full  <- apply_taxa_transform(test_concat_full,  t_columns, transform = "clr", pseudocount = taxa_pc, negative_policy = "error")
    message("Applied CLR taxa transform to full train/test (pseudocount=", taxa_pc, ").")
  }
  
  # NEW: collect effects tables here for coefficient stability
  metab_full_red <- apply_view_reduction(
    train_df = train_metab_full, val_df = test_metab_full,
    feature_cols = m_columns,
    response_variable = response_variable,
    method = args$metab_reduction,
    response_type = type_of_analysis,
    stratify_variable = stratify_variable,
    corr_cutoff = args$corr_cutoff,
    selector_top_k = 200,
    selector_fdr_q = 0.05,
    return_effects = TRUE,
    view_name = "metabolomics"
  )
  train_metab_final <- metab_full_red$train
  test_metab_final  <- metab_full_red$val
  kept_metab_final  <- metab_full_red$keep
  metab_effects_kept <- metab_full_red$effects
  
  taxa_full_red <- apply_view_reduction(
    train_df = train_mss_full, val_df = test_mss_full,
    feature_cols = t_columns,
    response_variable = response_variable,
    method = args$taxa_reduction,
    response_type = type_of_analysis,
    stratify_variable = stratify_variable,
    corr_cutoff = args$corr_cutoff,
    selector_top_k = 200,
    selector_fdr_q = 0.05,
    return_effects = TRUE,
    view_name = "taxa"
  )
  train_mss_final <- taxa_full_red$train
  test_mss_final  <- taxa_full_red$val
  kept_taxa_final <- taxa_full_red$keep
  taxa_effects_kept <- taxa_full_red$effects
  
  concat_full_red <- build_concat_auto(
    train_data = train_concat_full,
    val_data   = test_concat_full,
    response_variable = response_variable,
    kept_metab = kept_metab_final,
    kept_taxa  = kept_taxa_final
  )
  train_concat_final <- concat_full_red$train
  test_concat_final  <- concat_full_red$val
  
  final_concat_variants <- apply_concat_balance_variants(
    train_concat_final, test_concat_final,
    response_variable = response_variable,
    stratify_variable = stratify_variable,
    kept_metab_cols = kept_metab_final,
    kept_taxa_cols  = kept_taxa_final
  )
  
  train_concat_none_final <- final_concat_variants$concat_none$train
  test_concat_none_final  <- final_concat_variants$concat_none$val
  
  train_concat_scaling_final <- final_concat_variants$concat_scaling$train
  test_concat_scaling_final  <- final_concat_variants$concat_scaling$val
  
  train_concat_sw_final <- final_concat_variants$concat_scaling_weights$train
  test_concat_sw_final  <- final_concat_variants$concat_scaling_weights$val
  
  final_metab_fit  <- model_env$fit(train_metab_final,  response_variable, type_of_analysis, seed = 1000, positive_class = args$positive_class, use_gpu = args$use_gpu,
                                    gpu_device = args$gpu_device,
                                    nthread = args$nthread)
  final_mss_fit    <- model_env$fit(train_mss_final,    response_variable, type_of_analysis, seed = 2000, positive_class = args$positive_class, use_gpu = args$use_gpu,
                                    gpu_device = args$gpu_device,
                                    nthread = args$nthread)
  final_concat_none_fit <- model_env$fit(train_concat_none_final, response_variable, type_of_analysis, seed = 3001, positive_class = args$positive_class, use_gpu = args$use_gpu,
                                         gpu_device = args$gpu_device,
                                         nthread = args$nthread)
  final_concat_scaling_fit <- model_env$fit(train_concat_scaling_final, response_variable, type_of_analysis, seed = 3002, positive_class = args$positive_class, use_gpu = args$use_gpu,
                                            gpu_device = args$gpu_device,
                                            nthread = args$nthread)
  final_concat_sw_fit <- model_env$fit(train_concat_sw_final, response_variable, type_of_analysis,  seed = 3003, positive_class = args$positive_class, use_gpu = args$use_gpu,
                                       gpu_device = args$gpu_device,
                                       nthread = args$nthread)
  
  # write test importance from final refit models
  test_imp_list <- list(
    extract_fold_importance(final_metab_fit, model_type, "metabolomics", run_id, seed, NA_integer_, NA_integer_, split="test"),
    extract_fold_importance(final_mss_fit, model_type, "mss", run_id, seed, NA_integer_, NA_integer_, split="test"),
    extract_fold_importance(final_concat_none_fit, model_type, "concat_none", run_id, seed, NA_integer_, NA_integer_, split="test"),
    extract_fold_importance(final_concat_scaling_fit, model_type, "concat_scaling", run_id, seed, NA_integer_, NA_integer_, split="test"),
    extract_fold_importance(final_concat_sw_fit, model_type, "concat_scaling_weights", run_id, seed, NA_integer_, NA_integer_, split="test")
  )
  test_imp_df <- do.call(rbind, test_imp_list)
  
  imp_dir <- file.path(results_base_dir,  "FeatureImportanceLogs")
  dir.create(imp_dir, recursive = TRUE, showWarnings = FALSE)
  
  test_imp_file <- file.path(imp_dir, paste0(study_name, "_test_feature_importance.csv"))
  write.csv(test_imp_df, test_imp_file, row.names = FALSE)
  message("Saved test feature importance to: ", test_imp_file)
  
  test_pred_metab  <- model_env$predict(final_metab_fit,  test_metab_final,  response_variable, type_of_analysis, positive_class = args$positive_class)
  test_pred_mss    <- model_env$predict(final_mss_fit,    test_mss_final,    response_variable, type_of_analysis, positive_class = args$positive_class)
  test_pred_concat_none <- model_env$predict(final_concat_none_fit, test_concat_none_final, response_variable, type_of_analysis, positive_class = args$positive_class )
  test_pred_concat_scaling <- model_env$predict(final_concat_scaling_fit, test_concat_scaling_final, response_variable, type_of_analysis, positive_class = args$positive_class)
  test_pred_concat_sw <- model_env$predict(final_concat_sw_fit, test_concat_sw_final, response_variable, type_of_analysis, positive_class = args$positive_class)
  
  test_concat_none_pred <- test_pred_concat_none$pred
  test_concat_scaling_pred <- test_pred_concat_scaling$pred
  test_concat_sw_pred <- test_pred_concat_sw$pred
  
  
  test_truth <- test_pred_metab$truth
  
  
  
  stack_test <- data.frame(
    Metabolomics = test_pred_metab$pred,
    MSS          = test_pred_mss$pred
  )
  
  test_nnls_pred   <- as.numeric(as.matrix(stack_test[, c("Metabolomics", "MSS")]) %*% nnls_weights)
  test_average     <- 0.5 * stack_test$Metabolomics + 0.5 * stack_test$MSS
  test_sparse_pred <- as.numeric(predict(final_glmnet, newx = as.matrix(stack_test[, c("Metabolomics", "MSS")]), s = best_lambda))
  test_pls_pred <- drop(as.data.frame(predict(pls_model, stack_test[, c("Metabolomics", "MSS")], ncomp = best_ncomp))[, 1])
  
  final_test_metrics <- list(
    stacked_average = sapply(metrics, function(fn) fn(test_truth, test_average)),
    stacked_nnls    = sapply(metrics, function(fn) fn(test_truth, test_nnls_pred)),
    stacked_sparse  = sapply(metrics, function(fn) fn(test_truth, test_sparse_pred)),
    stacked_pls     = sapply(metrics, function(fn) fn(test_truth, test_pls_pred)),
    
    concat_none           = sapply(metrics, function(fn) fn(test_truth, test_concat_none_pred)),
    concat_scaling        = sapply(metrics, function(fn) fn(test_truth, test_concat_scaling_pred)),
    concat_scaling_weights= sapply(metrics, function(fn) fn(test_truth, test_concat_sw_pred)),
    
    metabolomics    = sapply(metrics, function(fn) fn(test_truth, test_pred_metab$pred)),
    mss             = sapply(metrics, function(fn) fn(test_truth, test_pred_mss$pred))
  )
  
  final_results <- list(
    meta_train = meta_train,
    meta_train_stack = meta_train_stack,
    meta_fit = list(
      nnls_weights = nnls_weights,
      best_alpha = best_alpha,
      best_lambda = best_lambda,
      sparse_nnls_weights = sparse_nnls_weights,
      best_ncomp = best_ncomp
    ),
    test_preds = list(
      metabolomics = test_pred_metab$pred,
      mss = test_pred_mss$pred,
      concat_none = test_concat_none_pred,
      concat_scaling = test_concat_scaling_pred,
      concat_scaling_weights = test_concat_sw_pred,
      stacked_average = test_average,
      stacked_nnls = test_nnls_pred,
      stacked_sparse = test_sparse_pred,
      stacked_pls = test_pls_pred
    ),
    test_truth = test_truth,
    final_test_metrics = final_test_metrics,
    kept_features = list(
      metabolomics = kept_metab_final,
      taxa = kept_taxa_final,
      concat = concat_full_red$keep
    ),
    kept_effects = list(
      metabolomics = metab_effects_kept,
      taxa = taxa_effects_kept
    ),
    fold_feature_log = fold_feature_log,
    run_id = run_id,
    seed = seed
  )
  

  # BUILD PERFORMANCE DATAFRAMES & SAVE RESULTS 
  stats_results <- function(lists, row_labels, metric, conf_level = 0.95) {
    values <- lapply(lists, function(x) {
      x <- unlist(x)
      x <- suppressWarnings(as.numeric(x))
      x[is.finite(x)]
    })
    means <- sapply(values, function(x) mean(x, na.rm = TRUE))
    ns <- sapply(values, function(x) sum(is.finite(x)))
    sds <- sapply(values, function(x) {
      x <- x[is.finite(x)]
      if (length(x) < 2) return(NA_real_)
      sd(x)
    })
    ses <- rep(NA_real_, length(ns))
    ok <- ns >= 2
    ses[ok] <- sds[ok] / sqrt(ns[ok])
    alpha <- 1 - conf_level
    t_crit <- rep(NA_real_, length(ns))
    t_crit[ok] <- qt(1 - alpha/2, df = ns[ok] - 1)
    lower_ci <- rep(NA_real_, length(ns))
    upper_ci <- rep(NA_real_, length(ns))
    lower_ci[ok] <- means[ok] - t_crit[ok] * ses[ok]
    upper_ci[ok] <- means[ok] + t_crit[ok] * ses[ok]
    result_df <- data.frame(
      Category = row_labels,
      Mean = means,
      SD = sds,
      N = ns,
      CI_Lower = lower_ci,
      CI_Upper = upper_ci
    )
    colnames(result_df)[2:6] <- paste(metric, c("_Mean", "_SD", "_N", "_CI_Lower", "_CI_Upper"), sep = " ")
    return(result_df)
  }
  
  make_lists_for_metric <- function(store, metric_name, categories) {
    out <- lapply(categories, function(cat) store[[cat]][[metric_name]])
    names(out) <- categories
    out
  }
  
  categories_train <- c(
    "metabolomics", "mss",
    "concat_none", "concat_scaling", "concat_scaling_weights",
    "stacked_average", "stacked_nnls", "stacked_sparse", "stacked_pls"
  )
  
  train_metric_dfs <- lapply(names(metrics), function(metric_name) {
    lists <- make_lists_for_metric(fold_metric_store, metric_name, categories_train)
    stats_results(lists, row_labels = categories_train, metric = metric_name)
  })
  
  train_df <- Reduce(function(x, y) dplyr::full_join(x, y, by = "Category"), train_metric_dfs)
  
  rename_train_cols <- function(df) {
    nm <- names(df)
    for (metric_name in names(metrics)) {
      MET <- toupper(metric_name)
      nm <- gsub(paste0("^", metric_name, " _Mean$"), paste0(MET, " Train _Mean"), nm)
      nm <- gsub(paste0("^", metric_name, " _SD$"),   paste0(MET, " Train _SD"), nm)
      nm <- gsub(paste0("^", metric_name, " _N$"),    paste0(MET, " Train _N"), nm)
      nm <- gsub(paste0("^", metric_name, " _CI_Lower$"), paste0(MET, " Train _CI_Lower"), nm)
      nm <- gsub(paste0("^", metric_name, " _CI_Upper$"), paste0(MET, " Train _CI_Upper"), nm)
    }
    names(df) <- nm
    df
  }
  train_df <- rename_train_cols(train_df)
  
  test_categories <- names(final_test_metrics)
  test_metric_dfs <- lapply(names(metrics), function(metric_name) {
    lists <- lapply(test_categories, function(cat) final_test_metrics[[cat]][[metric_name]])
    names(lists) <- test_categories
    stats_results(lists, row_labels = test_categories, metric = metric_name)
  })
  
  test_df <- Reduce(function(x, y) dplyr::full_join(x, y, by = "Category"), test_metric_dfs)
  
  rename_test_cols <- function(df) {
    nm <- names(df)
    for (metric_name in names(metrics)) {
      MET <- toupper(metric_name)
      nm <- gsub(paste0("^", metric_name, " _Mean$"), paste0(MET, " Test _Mean"), nm)
      nm <- gsub(paste0("^", metric_name, " _SD$"),   paste0(MET, " Test _SD"), nm)
      nm <- gsub(paste0("^", metric_name, " _N$"),    paste0(MET, " Test _N"), nm)
      nm <- gsub(paste0("^", metric_name, " _CI_Lower$"), paste0(MET, " Test _CI_Lower"), nm)
      nm <- gsub(paste0("^", metric_name, " _CI_Upper$"), paste0(MET, " Test _CI_Upper"), nm)
    }
    names(df) <- nm
    df
  }
  test_df <- rename_test_cols(test_df)
  
  all_results_df <- dplyr::full_join(train_df, test_df, by = "Category")
  
  dir.create(results_base_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------
  # SAVE PERFORMANCE CSV FIRST
  # ---------------------------
  performance_dir <- file.path(results_base_dir, "performance_csv")
  dir.create(performance_dir, recursive = TRUE, showWarnings = FALSE)
  performance_filename <- file.path(performance_dir, paste0(study_name, ".csv"))
  tmp_csv <- paste0(performance_filename, ".tmp")
  
  write.csv(all_results_df, file = tmp_csv, row.names = FALSE)
  if (!file.rename(tmp_csv, performance_filename)) {
    if (file.exists(tmp_csv)) unlink(tmp_csv)
    stop("Failed to rename CSV temp file to final output: ", tmp_csv)
  }
  message("Performance CSV saved to: ", performance_filename)
  
  # ---------------------------
  # WRITE DONE MARKER BEFORE RDATA
  # ---------------------------
  done_marker <- done_marker_for_run(study_name)
  tmp_done <- paste0(done_marker, ".tmp")
  ok_dm <- tryCatch({
    writeLines(as.character(Sys.time()), tmp_done)
    file.rename(tmp_done, done_marker)
    TRUE
  }, error = function(e) {
    if (file.exists(tmp_done)) unlink(tmp_done)
    warning("Failed to write done marker for run: ", study_name, " (", conditionMessage(e), ")")
    FALSE
  })
  if (isTRUE(ok_dm)) {
    message("Wrote done marker: ", done_marker)
  }
  
  # ---------------------------
  # SAVE RDATA LAST
  # If this fails, warn only.
  # ---------------------------
  rdata_dir <- file.path(results_base_dir, "RDataFiles")
  dir.create(rdata_dir, recursive = TRUE, showWarnings = FALSE)
  rdata_filename <- file.path(rdata_dir, paste0(study_name, ".RData"))
  tmp_rdata <- paste0(rdata_filename, ".tmp")
  
  testfile <- file.path(rdata_dir, paste0("write_test_", Sys.getpid(), ".tmp"))
  ok <- tryCatch({
    writeLines("test", testfile)
    unlink(testfile)
    TRUE
  }, error = function(e) e)
  
  if (!isTRUE(ok)) {
    warning("Results directory not writable for RData: ", conditionMessage(ok))
    rdata_filename <- NA_character_
  } else {
    save_ok <- tryCatch({
      save(all_results_df, final_results, fold_metric_store, meta_train, meta_train_stack,
           file = tmp_rdata, compress = "xz")
      if (!file.rename(tmp_rdata, rdata_filename)) {
        stop("Failed to rename RData temp file to final output: ", tmp_rdata)
      }
      TRUE
    }, error = function(e) e)
    
    if (!isTRUE(save_ok)) {
      if (file.exists(tmp_rdata)) unlink(tmp_rdata)
      warning(sprintf(
        "Failed to save RData [%s]. Object size: %.2f MB. Error: %s",
        tmp_rdata,
        as.numeric(object.size(final_results)) / 1024^2,
        conditionMessage(save_ok)
      ))
      rdata_filename <- NA_character_
    } else {
      message("RData file saved to: ", rdata_filename)
    }
  }
  
  
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  hours <- floor(elapsed_time / 3600); minutes <- floor((elapsed_time %% 3600) / 60)
  cat(sprintf("Elapsed time: %d hours and %d minutes or %f seconds\n", hours, minutes, elapsed_time))
  print(warnings())
  sink()
  
  # ---------------------------
  # Return SMALL summary (do not keep large model objects in memory)
  # ---------------------------
  small <- list(
    run_id = run_id,
    seed   = seed,
    study_name = study_name,
    
    # paths to outputs for debugging / provenance
    rdata_file = rdata_filename,
    performance_csv = performance_filename,
    fold_perf_file = fold_perf_file,
    repeat_perf_file = repeat_perf_file,
    val_imp_file = if (exists("val_imp_file")) val_imp_file else NA_character_,
    test_imp_file = if (exists("test_imp_file")) test_imp_file else NA_character_,
    
    # items your stability code expects
    kept_features = final_results$kept_features,
    kept_effects  = final_results$kept_effects,
    fold_feature_log = final_results$fold_feature_log
  )
  
  if (isTRUE(return_big)) {
    small$final_results <- final_results
  }
  
  # aggressively free memory inside this run
  rm(final_results)
  gc()
  
  return(small)
  } # end run_one # end run_one

# ---------------------------
# Diagnostics across seeds (unchanged logic) and run loop
# ---------------------------
#seeds <- 1001:1020
seeds <- 1001:1020
run_summaries <- vector("list", length(seeds))

get_test_indices_for_seed <- function(seed) {
  set.seed(seed)
  if (!is.null(stratify_variable) && stratify_variable != "" && stratify_variable %in% names(data)) {
    parts <- caret::createDataPartition(data[[stratify_variable]], p = training_proportion, list = FALSE)
  } else {
    parts <- caret::createDataPartition(data[[response_variable]], p = training_proportion, list = FALSE)
  }
  test_idx <- setdiff(seq_len(nrow(data)), parts)
  return(test_idx)
}

test_indices_list <- lapply(seeds, get_test_indices_for_seed)
names(test_indices_list) <- as.character(seeds)

overlap_matrix <- matrix(NA_real_, length(seeds), length(seeds),
                         dimnames = list(as.character(seeds), as.character(seeds)))
jaccard_matrix <- matrix(NA_real_, length(seeds), length(seeds),
                         dimnames = list(as.character(seeds), as.character(seeds)))

for (i in seq_along(seeds)) {
  A <- test_indices_list[[i]]
  for (j in seq_along(seeds)) {
    B <- test_indices_list[[j]]
    overlap_matrix[i, j] <- length(intersect(A, B)) / length(A)
    jaccard_matrix[i, j] <- length(intersect(A, B)) / length(union(A, B))
  }
}

#results_base_dir <- file.path(args$file_path, "A_Integrative_Pipeline_Scripts", "results")
dir.create(results_base_dir, recursive = TRUE, showWarnings = FALSE)
diag_dir <- file.path(results_base_dir, "SplitDiagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(overlap_matrix, file = file.path(diag_dir, paste0(args$study_name, "_test_overlap_matrix.csv")))
write.csv(jaccard_matrix, file = file.path(diag_dir, paste0(args$study_name, "_test_jaccard_matrix.csv")))

png(filename = file.path(diag_dir, paste0(args$study_name, "_test_overlap_heatmap.png")), width = 1200, height = 1000, res = 150)
heatmap(overlap_matrix, symm = TRUE, main = "Test-set overlap (|A∩B| / |A|)")
dev.off()
png(filename = file.path(diag_dir, paste0(args$study_name, "_test_jaccard_heatmap.png")), width = 1200, height = 1000, res = 150)
heatmap(jaccard_matrix, symm = TRUE, main = "Test-set Jaccard similarity (|A∩B| / |A∪B|)")
dev.off()

off_diag_overlap <- overlap_matrix[lower.tri(overlap_matrix)]
message(sprintf("Test overlap across seeds (%d splits): mean=%.3f, min=%.3f, max=%.3f",
                length(seeds), mean(off_diag_overlap), min(off_diag_overlap), max(off_diag_overlap)))
off_diag_jacc <- jaccard_matrix[lower.tri(jaccard_matrix)]
message(sprintf("Test Jaccard across seeds (%d splits): mean=%.3f, min=%.3f, max=%.3f",
                length(seeds), mean(off_diag_jacc), min(off_diag_jacc), max(off_diag_jacc)))

for (i in seq_along(seeds)) {
  run_id <- i
  seed_val <- seeds[i]
  # compute expected study_name for this run (same logic as run_one uses)
  study_name_run <- study_name_for_run(args$study_name, run_id)
  done_marker <- done_marker_for_run(study_name_run)
  rdata_file  <- rdata_for_run(study_name_run)
  
  if (file.exists(done_marker)) {
    message(sprintf("Skipping run %d (seed=%d): done marker exists: %s", run_id, seed_val, done_marker))
    # Try to recover a small summary into run_summaries (kept features, fold logs)
    s <- tryCatch({
      if (file.exists(rdata_file)) {
        # load only into a temporary environment to avoid polluting workspace & not keep giant objects
        env_tmp <- new.env()
        load(rdata_file, envir = env_tmp)
        # heuristics: expected objects might be 'final_results' or 'all_results_df' etc.
        small_res <- list(run_id = run_id, seed = seed_val)
        # attempt to pull kept_features and fold_feature_log from final_results if present
        if (exists("final_results", envir = env_tmp)) {
          fr <- get("final_results", envir = env_tmp)
          if (is.list(fr) && !is.null(fr$kept_features)) small_res$kept_features <- fr$kept_features
          if (is.list(fr) && !is.null(fr$fold_feature_log)) small_res$fold_feature_log <- fr$fold_feature_log
        }
        # fallback: if variables stored directly
        if (exists("kept_features", envir = env_tmp) && is.null(small_res$kept_features)) {
          small_res$kept_features <- get("kept_features", envir = env_tmp)
        }
        if (exists("fold_feature_log", envir = env_tmp) && is.null(small_res$fold_feature_log)) {
          small_res$fold_feature_log <- get("fold_feature_log", envir = env_tmp)
        }
        # clean up the temporary environment
        rm(env_tmp)
        gc()
        small_res
      } else {
        list(run_id = run_id, seed = seed_val, note = "done marker exists but rdata missing")
      }
    }, error = function(e) {
      warning("Failed to load rdata for skipped run ", study_name_run, ": ", conditionMessage(e))
      list(run_id = run_id, seed = seed_val, note = paste("load_error:", conditionMessage(e)))
    })
    run_summaries[[i]] <- s
    next
  }
  
  # not done: execute run
  message(sprintf("Starting run %d (seed=%d) — no done marker found.", run_id, seed_val))
  rs <- tryCatch(
    run_one(args = args, run_id = run_id, seed = seed_val, return_big = FALSE),
    error = function(e) {
      warning(sprintf("run_one failed for run %d (seed=%d): %s", run_id, seed_val, conditionMessage(e)))
      NULL
    }
  )
  
  # store summary returned (should be a small summary if you implemented that earlier)
  run_summaries[[i]] <- rs
  
  # free some memory after each run
  gc()
}

# ---------------------------
# NEW: Stability reports across the 20 splits (feature frequency + coefficient stability)
# ---------------------------

base_study_name <- args$study_name

kept_metab_list <- list()
kept_taxa_list  <- list()
kept_concat_list <- list()

effects_metab_list <- list()
effects_taxa_list  <- list()

fold_logs_all <- list()

for (i in seq_along(run_summaries)) {
  rs <- run_summaries[[i]]
  if (is.null(rs)) next
  
  if (!is.null(rs$kept_features)) {
    kept_metab_list[[as.character(seeds[i])]]  <- rs$kept_features$metabolomics
    kept_taxa_list[[as.character(seeds[i])]]   <- rs$kept_features$taxa
    kept_concat_list[[as.character(seeds[i])]] <- rs$kept_features$concat
  }
  
  if (!is.null(rs$kept_effects)) {
    # annotate with seed/run for traceability
    if (!is.null(rs$kept_effects$metabolomics) && nrow(rs$kept_effects$metabolomics) > 0) {
      dfm <- rs$kept_effects$metabolomics
      dfm$Seed <- seeds[i]
      dfm$RunID <- i
      effects_metab_list[[as.character(seeds[i])]] <- dfm
    }
    if (!is.null(rs$kept_effects$taxa) && nrow(rs$kept_effects$taxa) > 0) {
      dft <- rs$kept_effects$taxa
      dft$Seed <- seeds[i]
      dft$RunID <- i
      effects_taxa_list[[as.character(seeds[i])]] <- dft
    }
  }
  
  if (!is.null(rs$fold_feature_log) && nrow(rs$fold_feature_log) > 0) {
    fold_logs_all[[length(fold_logs_all) + 1]] <- rs$fold_feature_log
  }
}

n_runs_ok <- length(kept_metab_list)
message("Stability report: number of successful runs included = ", n_runs_ok)

build_freq_df <- function(kept_list, view_name) {
  if (length(kept_list) == 0) {
    return(data.frame(
      View = view_name,
      Feature = character(),
      Count_Selected = integer(),
      Proportion_Selected = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  all_feats <- sort(unique(unlist(kept_list)))
  counts <- sapply(all_feats, function(f) sum(sapply(kept_list, function(x) f %in% x)))
  
  df <- data.frame(
    View = view_name,
    Feature = names(counts),
    Count_Selected = as.integer(counts),
    Proportion_Selected = as.numeric(counts) / length(kept_list),
    stringsAsFactors = FALSE
  )
  df <- df[order(df$Count_Selected, decreasing = TRUE), , drop = FALSE]
  rownames(df) <- NULL
  df
}

metab_stab  <- build_freq_df(kept_metab_list,  "metabolomics")
taxa_stab   <- build_freq_df(kept_taxa_list,   "taxa")
concat_stab <- build_freq_df(kept_concat_list, "concat")
stability_all <- rbind(metab_stab, taxa_stab, concat_stab)

# Always save stability outputs alongside the main run outputs
stab_dir <- file.path(results_base_dir, "StabilityReports")
dir.create(stab_dir, recursive = TRUE, showWarnings = FALSE)

stab_file_all <- file.path(stab_dir, paste0(base_study_name, "_feature_stability_all.csv"))
write.csv(stability_all, stab_file_all, row.names = FALSE)
message("Feature stability report saved to: ", stab_file_all)

# Aggregate fold feature counts across all runs
if (length(fold_logs_all) > 0) {
  all_fold_logs <- do.call(rbind, fold_logs_all)
  foldlog_dir <- file.path(results_base_dir, "FeatureCountLogs")
  dir.create(foldlog_dir, recursive = TRUE, showWarnings = FALSE)
  foldlog_all_file <- file.path(foldlog_dir, paste0(base_study_name, "_ALLRUNS_fold_feature_counts.csv"))
  write.csv(all_fold_logs, foldlog_all_file, row.names = FALSE)
  message("Aggregated fold feature count logs saved to: ", foldlog_all_file)
  
  message(sprintf(
    "Fold keep summary (avg across all folds/runs): metab=%.1f taxa=%.1f concat=%.1f",
    mean(all_fold_logs$n_metab_keep, na.rm = TRUE),
    mean(all_fold_logs$n_taxa_keep, na.rm = TRUE),
    mean(all_fold_logs$n_concat_keep, na.rm = TRUE)
  ))
}

# ---------------------------
# NEW: Coefficient/effect stability report
# - Only meaningful when methods produce effect sizes:
#   - metabolomics: limma (logFC)
#   - taxa: wilcox (median diff)
# ---------------------------
summarize_effects <- function(effects_list, view_name) {
  if (length(effects_list) == 0) {
    return(data.frame(
      View = character(0),
      Feature = character(0),
      N_Selected = integer(0),
      Mean_Effect = numeric(0),
      SD_Effect = numeric(0),
      Median_Effect = numeric(0),
      Prop_Positive = numeric(0),
      Mean_AdjP = numeric(0),
      Level1 = character(0),
      Level2 = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  all_df <- do.call(rbind, effects_list)
  
  # extra guard: list exists but rbind produced 0 rows (e.g., all were empty)
  if (is.null(all_df) || nrow(all_df) == 0) {
    return(data.frame(
      View = character(0),
      Feature = character(0),
      N_Selected = integer(0),
      Mean_Effect = numeric(0),
      SD_Effect = numeric(0),
      Median_Effect = numeric(0),
      Prop_Positive = numeric(0),
      Mean_AdjP = numeric(0),
      Level1 = character(0),
      Level2 = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  agg <- all_df %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(
      N_Selected = dplyr::n(),
      Mean_Effect = mean(Effect, na.rm = TRUE),
      SD_Effect = sd(Effect, na.rm = TRUE),
      Median_Effect = median(Effect, na.rm = TRUE),
      Prop_Positive = mean(Effect > 0, na.rm = TRUE),
      Mean_AdjP = mean(AdjP, na.rm = TRUE),
      Level1 = dplyr::first(Level1),
      Level2 = dplyr::first(Level2),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(N_Selected), dplyr::desc(abs(Mean_Effect)))
  
  out <- as.data.frame(agg)
  out$View <- view_name
  out <- out[, c("View", "Feature", "N_Selected", "Mean_Effect", "SD_Effect", "Median_Effect",
                 "Prop_Positive", "Mean_AdjP", "Level1", "Level2")]
  out
}

metab_eff_stab <- summarize_effects(effects_metab_list, "metabolomics")
taxa_eff_stab  <- summarize_effects(effects_taxa_list,  "taxa")

if (length(effects_metab_list) == 0) message("No metabolomics effects collected (metab_reduction not limma?).")
if (length(effects_taxa_list) == 0) message("No taxa effects collected (taxa_reduction not wilcox?).")

coef_stability_all <- rbind(metab_eff_stab, taxa_eff_stab)

coef_file_all <- file.path(stab_dir, paste0(base_study_name, "_coefficient_stability_all.csv"))
write.csv(coef_stability_all, coef_file_all, row.names = FALSE)
message("Coefficient/effect stability report saved to: ", coef_file_all)

# Helpful: merge frequency + coefficient stability for metab/taxa
freq_metab <- stability_all[stability_all$View == "metabolomics", c("Feature", "Count_Selected", "Proportion_Selected")]
freq_taxa  <- stability_all[stability_all$View == "taxa", c("Feature", "Count_Selected", "Proportion_Selected")]

metab_merged <- dplyr::left_join(freq_metab, metab_eff_stab, by = "Feature")
taxa_merged  <- dplyr::left_join(freq_taxa,  taxa_eff_stab,  by = "Feature")

write.csv(metab_merged, file.path(stab_dir, paste0(base_study_name, "_metabolomics_freq_plus_effects.csv")), row.names = FALSE)
write.csv(taxa_merged,  file.path(stab_dir, paste0(base_study_name, "_taxa_freq_plus_effects.csv")), row.names = FALSE)
message("Merged frequency+effects saved for metabolomics and taxa.")
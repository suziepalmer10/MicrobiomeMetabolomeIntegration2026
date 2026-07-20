assert_xgb_predictors_clean <- function(x, context = "XGBoost predictors") {
  predictor_names <- names(x)
  
  if (is.null(predictor_names)) {
    stop(context, ": predictor matrix has no column names.")
  }
  
  if (anyDuplicated(predictor_names)) {
    duplicated_names <- unique(predictor_names[duplicated(predictor_names)])
    
    stop(
      context,
      ": duplicate predictor names detected: ",
      paste(duplicated_names, collapse = ", ")
    )
  }
  
  if (".outcome" %in% predictor_names) {
    stop(
      context,
      ": .outcome is present among predictors. ",
      "This indicates true outcome bleed-over."
    )
  }
  
  bad_predictor_names <- predictor_names[
    !grepl("^(m__|t__)", predictor_names)
  ]
  
  if (length(bad_predictor_names) > 0) {
    stop(
      context,
      ": non-metabolite/non-taxon predictor columns detected: ",
      paste(bad_predictor_names, collapse = ", "),
      ". Expected all XGBoost predictors to start with m__ or t__."
    )
  }
  
  invisible(TRUE)
}


fit <- function(trainData, response_variable, response_type,
                stratify_variable = NULL, seed = 123,
                tune_grid = NULL, cv_folds = 5, positive_class = NULL,
                use_gpu = TRUE, gpu_device = "cuda", nthread = 16) {
  
  set.seed(seed)
  trainData <- as.data.frame(trainData)
  
  if (!response_variable %in% names(trainData)) {
    stop("response_variable '", response_variable, "' not found in trainData.")
  }
  
  drop_cols <- c(response_variable)
  
  if (!is.null(stratify_variable) && nzchar(stratify_variable) &&
      stratify_variable %in% names(trainData)) {
    drop_cols <- c(drop_cols, stratify_variable)
  }
  
  x <- trainData[, setdiff(names(trainData), drop_cols), drop = FALSE]
  y_raw <- trainData[[response_variable]]
  
  if (ncol(x) == 0) {
    stop("No predictor columns remain after dropping response/stratify columns.")
  }
  
  assert_xgb_predictors_clean(
    x,
    context = "XGBoost fit() predictor matrix"
  )
  
  if (anyNA(x)) {
    bad_cols <- names(which(colSums(is.na(x)) > 0))
    
    stop(
      "Predictor matrix contains NA values. Columns with NA: ",
      paste(bad_cols, collapse = ", ")
    )
  }
  
  non_numeric <- names(x)[!vapply(x, is.numeric, logical(1))]
  
  if (length(non_numeric) > 0) {
    stop(
      "All predictors must be numeric for XGBoost. Non-numeric columns: ",
      paste(non_numeric, collapse = ", ")
    )
  }
  
  x_mat <- as.matrix(x)
  feature_names <- colnames(x_mat)
  
  if (is.null(tune_grid)) {
    tune_grid <- expand.grid(
      nrounds = c(300, 600),
      max_depth = c(3, 6),
      eta = c(0.05, 0.10),
      gamma = 0,
      colsample_bytree = c(0.7, 1.0),
      min_child_weight = 1,
      subsample = c(0.7, 1.0)
    )
  }
  
  if (!is.data.frame(tune_grid)) {
    stop("tune_grid must be a data.frame or expand.grid result.")
  }
  
  required_cols <- c(
    "nrounds", "max_depth", "eta", "gamma",
    "colsample_bytree", "min_child_weight", "subsample"
  )
  
  missing_cols <- setdiff(required_cols, names(tune_grid))
  
  if (length(missing_cols) > 0) {
    stop(
      "tune_grid is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (response_type == "binary") {
    y_factor <- as.factor(y_raw)
    
    if (anyNA(y_factor)) {
      stop("Binary response contains NA values.")
    }
    
    if (length(levels(y_factor)) != 2) {
      stop("Binary response must have exactly 2 levels.")
    }
    
    levs <- levels(y_factor)
    
    if (is.null(positive_class) || !nzchar(positive_class)) {
      positive_class <- levs[2]
    }
    
    pos_match <- which(tolower(levs) == tolower(positive_class))
    
    if (length(pos_match) != 1) {
      stop(
        "positive_class '", positive_class, "' not found in outcome levels: ",
        paste(levs, collapse = ", ")
      )
    }
    
    positive_class <- levs[pos_match]
    negative_class <- setdiff(levs, positive_class)
    
    y_factor <- factor(y_factor, levels = c(positive_class, negative_class))
    y_num <- as.numeric(y_factor == positive_class)
    
    class_counts <- table(y_factor)
    min_class <- min(class_counts)
    
    if (min_class < 2) {
      stop(
        "Binary training fold has fewer than 2 observations in the minority class. ",
        "Inner cross-validation is not feasible. Class counts: ",
        paste(names(class_counts), as.integer(class_counts), collapse = ", ")
      )
    }
    
    cv_folds_used <- min(cv_folds, min_class)
    
    if (cv_folds_used < 2) {
      stop("cv_folds must be at least 2 for binary classification.")
    }
    
    if (cv_folds_used < cv_folds) {
      warning(
        "Reduced cv_folds from ", cv_folds, " to ", cv_folds_used,
        " due to small class count in this training fold."
      )
    }
    
    inner_index_train <- caret::createFolds(
      y_factor,
      k = cv_folds_used,
      returnTrain = TRUE
    )
    
    all_idx <- seq_len(length(y_num))
    
    xgb_folds <- lapply(
      inner_index_train,
      function(tr_idx) setdiff(all_idx, tr_idx)
    )
    
    objective <- "binary:logistic"
    eval_metric <- "auc"
    maximize_metric <- TRUE
    label_vec <- y_num
    
  } else if (response_type == "continuous") {
    y_num <- as.numeric(y_raw)
    
    if (anyNA(y_num)) {
      stop("Response contains NA values (continuous).")
    }
    
    if (length(y_num) < 3) {
      stop("Continuous training data has too few rows for CV.")
    }
    
    cv_folds_used <- min(cv_folds, nrow(x_mat))
    
    if (cv_folds_used < 2) {
      stop("cv_folds must be at least 2 for continuous outcomes.")
    }
    
    if (cv_folds_used < cv_folds) {
      warning(
        "Reduced cv_folds from ", cv_folds, " to ", cv_folds_used,
        " due to small training set size."
      )
    }
    
    set.seed(seed)
    
    shuffled <- sample(seq_len(nrow(x_mat)))
    fold_ids <- cut(seq_along(shuffled), breaks = cv_folds_used, labels = FALSE)
    xgb_folds <- split(shuffled, fold_ids)
    
    objective <- "reg:squarederror"
    eval_metric <- "rmse"
    maximize_metric <- FALSE
    label_vec <- y_num
    positive_class <- NULL
    
  } else {
    stop("Unsupported response_type: ", response_type)
  }
  
  dtrain <- xgboost::xgb.DMatrix(
    data = x_mat,
    label = label_vec
  )
  
  colnames(dtrain) <- feature_names
  
  base_params <- list(
    booster = "gbtree",
    objective = objective,
    eval_metric = eval_metric,
    tree_method = "hist",
    nthread = nthread,
    verbosity = 1
  )
  
  if (isTRUE(use_gpu)) {
    base_params$device <- gpu_device
  } else {
    base_params$device <- "cpu"
  }
  
  message("xgboost version: ", as.character(utils::packageVersion("xgboost")))
  message("use_gpu=", use_gpu, " gpu_device=", gpu_device, " nthread=", nthread)
  message("CUDA_VISIBLE_DEVICES=", Sys.getenv("CUDA_VISIBLE_DEVICES", unset = "unset"))
  message("base_params:")
  print(base_params)
  
  score_grid <- vector("list", nrow(tune_grid))
  
  for (i in seq_len(nrow(tune_grid))) {
    tg <- tune_grid[i, , drop = FALSE]
    
    params_i <- utils::modifyList(
      base_params,
      list(
        max_depth = as.integer(tg$max_depth),
        eta = as.numeric(tg$eta),
        gamma = as.numeric(tg$gamma),
        colsample_bytree = as.numeric(tg$colsample_bytree),
        min_child_weight = as.numeric(tg$min_child_weight),
        subsample = as.numeric(tg$subsample)
      )
    )
    
    message("Starting xgb.cv: ", Sys.time(), " device=", gpu_device, " use_gpu=", use_gpu)
    message("CUDA_VISIBLE_DEVICES=", Sys.getenv("CUDA_VISIBLE_DEVICES", unset = "unset"))
    message("params_i:")
    print(params_i)
    
    cv_fit <- xgboost::xgb.cv(
      params = params_i,
      data = dtrain,
      nrounds = as.integer(tg$nrounds),
      folds = xgb_folds,
      prediction = FALSE,
      verbose = 0
    )
    
    message("Finished xgb.cv: ", Sys.time())
    
    elog <- cv_fit$evaluation_log
    
    metric_col <- paste0("test_", eval_metric, "_mean")
    
    if (!metric_col %in% colnames(elog)) {
      stop(
        "Expected CV metric column '", metric_col, "' not found in evaluation_log. ",
        "Columns are: ", paste(colnames(elog), collapse = ", ")
      )
    }
    
    final_score <- elog[[metric_col]][nrow(elog)]
    
    score_grid[[i]] <- data.frame(
      tune_grid[i, , drop = FALSE],
      score = as.numeric(final_score),
      stringsAsFactors = FALSE
    )
  }
  
  score_df <- do.call(rbind, score_grid)
  
  if (maximize_metric) {
    best_idx <- which.max(score_df$score)
  } else {
    best_idx <- which.min(score_df$score)
  }
  
  best_tune <- tune_grid[best_idx, , drop = FALSE]
  
  final_params <- utils::modifyList(
    base_params,
    list(
      max_depth = as.integer(best_tune$max_depth),
      eta = as.numeric(best_tune$eta),
      gamma = as.numeric(best_tune$gamma),
      colsample_bytree = as.numeric(best_tune$colsample_bytree),
      min_child_weight = as.numeric(best_tune$min_child_weight),
      subsample = as.numeric(best_tune$subsample)
    )
  )
  
  message("Starting xgb.train: ", Sys.time(), " device=", gpu_device, " use_gpu=", use_gpu)
  message("CUDA_VISIBLE_DEVICES=", Sys.getenv("CUDA_VISIBLE_DEVICES", unset = "unset"))
  message("final_params:")
  print(final_params)
  
  final_model <- xgboost::xgb.train(
    params = final_params,
    data = dtrain,
    nrounds = as.integer(best_tune$nrounds),
    verbose = 0
  )
  
  message("Finished xgb.train: ", Sys.time())
  
  fit_obj <- list(
    finalModel = final_model,
    bestTune = best_tune,
    
    # Stored for metadata compatibility only.
    # Do not use trainingData to extract XGBoost feature names.
    trainingData = data.frame(.outcome = y_raw, x, check.names = FALSE),
    
    # Source of truth for XGBoost feature importance and prediction.
    feature_names = feature_names,
    
    response_type = response_type,
    positive_class = positive_class,
    use_gpu = isTRUE(use_gpu),
    gpu_device = if (isTRUE(use_gpu)) gpu_device else "cpu",
    nthread = nthread,
    metric = eval_metric,
    cv_results = score_df
  )
  
  class(fit_obj) <- c("xgb_train_wrapper", "list")
  
  return(fit_obj)
}


predict <- function(fit, newData, response_variable, response_type,
                    stratify_variable = NULL, positive_class = NULL) {
  
  newData <- as.data.frame(newData)
  
  if (!response_variable %in% names(newData)) {
    stop("response_variable '", response_variable, "' not found in newData.")
  }
  
  drop_cols <- c(response_variable)
  
  if (!is.null(stratify_variable) && nzchar(stratify_variable) &&
      stratify_variable %in% names(newData)) {
    drop_cols <- c(drop_cols, stratify_variable)
  }
  
  x_new <- newData[, setdiff(names(newData), drop_cols), drop = FALSE]
  y_true <- newData[[response_variable]]
  
  if (ncol(x_new) == 0) {
    stop("No predictor columns remain in newData after dropping response/stratify columns.")
  }
  
  if (anyNA(x_new)) {
    bad_cols <- names(which(colSums(is.na(x_new)) > 0))
    
    stop(
      "Prediction matrix contains NA values. Columns with NA: ",
      paste(bad_cols, collapse = ", ")
    )
  }
  
  non_numeric <- names(x_new)[!vapply(x_new, is.numeric, logical(1))]
  
  if (length(non_numeric) > 0) {
    stop(
      "All predictors in newData must be numeric for XGBoost. Non-numeric columns: ",
      paste(non_numeric, collapse = ", ")
    )
  }
  
  if (is.null(fit$finalModel) || is.null(fit$feature_names)) {
    stop("fit object does not contain expected XGBoost model metadata.")
  }
  
  missing_features <- setdiff(fit$feature_names, colnames(x_new))
  
  if (length(missing_features) > 0) {
    stop(
      "newData is missing required predictor columns: ",
      paste(missing_features, collapse = ", ")
    )
  }
  
  x_new <- x_new[, fit$feature_names, drop = FALSE]
  
  assert_xgb_predictors_clean(
    x_new,
    context = "XGBoost predict() aligned predictor matrix"
  )
  
  dnew <- xgboost::xgb.DMatrix(
    data = as.matrix(x_new)
  )
  
  colnames(dnew) <- fit$feature_names
  
  p <- stats::predict(
    fit$finalModel,
    newdata = dnew
  )
  
  if (response_type == "binary") {
    y_true <- as.factor(y_true)
    
    if (anyNA(y_true)) {
      stop("Binary truth contains NA values.")
    }
    
    if (length(levels(y_true)) != 2) {
      stop("Binary truth must have exactly 2 levels.")
    }
    
    levs <- levels(y_true)
    
    if (is.null(positive_class) || !nzchar(positive_class)) {
      if (!is.null(fit$positive_class) && nzchar(fit$positive_class)) {
        positive_class <- fit$positive_class
      } else {
        positive_class <- levs[2]
      }
    }
    
    pos_match <- which(tolower(levs) == tolower(positive_class))
    
    if (length(pos_match) != 1) {
      stop(
        "positive_class '", positive_class, "' not found in truth levels: ",
        paste(levs, collapse = ", ")
      )
    }
    
    positive_class_raw <- levs[pos_match]
    truth <- as.numeric(y_true == positive_class_raw)
    
    return(list(
      pred = as.numeric(p),
      truth = as.numeric(truth)
    ))
  }
  
  if (response_type == "continuous") {
    truth <- as.numeric(y_true)
    
    if (anyNA(truth)) {
      stop("Continuous truth contains NA values.")
    }
    
    return(list(
      pred = as.numeric(p),
      truth = truth
    ))
  }
  
  stop("Unsupported response_type: ", response_type)
}

fit <- function(trainData, response_variable, response_type,
                stratify_variable = NULL, seed = 123,
                tune_grid = NULL, cv_folds = 5, positive_class = NULL, 
                use_gpu = FALSE, gpu_device = "cuda", nthread = NULL) {
  
  set.seed(seed)
  trainData <- as.data.frame(trainData)
  
  drop_cols <- c(response_variable)
  if (!is.null(stratify_variable) && stratify_variable != "" &&
      stratify_variable %in% names(trainData)) {
    drop_cols <- c(drop_cols, stratify_variable)
  }
  
  x <- trainData[, setdiff(names(trainData), drop_cols), drop = FALSE]
  y <- trainData[[response_variable]]
  
  if (is.null(tune_grid)) {
    p <- ncol(x)
    mtry_vals <- unique(pmax(1, pmin(p, round(seq(1, sqrt(p), length.out = 10)))))
    tune_grid <- expand.grid(mtry = mtry_vals)
  }
  
  if (response_type == "binary") {
    y <- as.factor(y)
    if (length(levels(y)) != 2) stop("Binary response must have exactly 2 levels.")
    
    levs <- levels(y)
    
    # If not provided, fall back to old convention: 2nd level is positive
    if (is.null(positive_class) || !nzchar(positive_class)) {
      positive_class <- levs[2]
    }
    
    # Case-insensitive match to actual factor levels
    pos_match <- which(tolower(levs) == tolower(positive_class))
    if (length(pos_match) != 1) {
      stop("positive_class '", positive_class, "' not found in outcome levels: ",
           paste(levs, collapse = ", "))
    }
    positive_class <- levs[pos_match]
    negative_class <- setdiff(levs, positive_class)
    
    # IMPORTANT: caret twoClassSummary treats FIRST level as the "event"/positive
    y <- factor(y, levels = c(positive_class, negative_class))
    levs <- levels(y)  # event=levs[1] == positive_class
    
    # guard: can't have more folds than smallest class count
    min_class <- min(table(y))
    if (cv_folds > min_class) {
      cv_folds <- min_class
      warning("Reduced cv_folds to ", cv_folds, " due to small class count in training fold.")
    }
    
    # stratified inner CV folds
    inner_index <- caret::createFolds(y, k = cv_folds, returnTrain = TRUE)
    
    # optional diagnostics
    cat("\n[RF inner-CV diagnostics] seed=", seed,
        " cv_folds=", cv_folds,
        " positive_class=", positive_class,
        " levels=", paste(levs, collapse = ", "), "\n", sep = "")
    for (i in seq_along(inner_index)) {
      tr <- inner_index[[i]]
      te <- setdiff(seq_along(y), tr)
      cat("  Inner fold ", i, ":\n", sep = "")
      cat("    train:", paste(names(table(y[tr])), as.integer(table(y[tr])), collapse = " "), "\n")
      cat("    test :", paste(names(table(y[te])), as.integer(table(y[te])), collapse = " "), "\n")
    }
    
    ctrl <- caret::trainControl(
      method = "cv",
      number = cv_folds,
      index = inner_index,
      classProbs = TRUE,
      summaryFunction = caret::twoClassSummary,
      savePredictions = "final"
    )
    
    metric <- "ROC"
    
  } else if (response_type == "continuous") {
    
    y <- as.numeric(y)
    ctrl <- caret::trainControl(method = "cv", number = cv_folds)
    metric <- "RMSE"
    
  } else {
    stop("Unsupported response_type: ", response_type)
  }
  
  caret::train(
    x = x,
    y = y,
    method = "rf",
    tuneGrid = tune_grid,
    trControl = ctrl,
    metric = metric
  )
}


predict <- function(fit, newData, response_variable, response_type,
                    stratify_variable = NULL, positive_class = NULL) {
  
  newData <- as.data.frame(newData)
  
  drop_cols <- c(response_variable)
  if (!is.null(stratify_variable) && stratify_variable != "" &&
      stratify_variable %in% names(newData)) {
    drop_cols <- c(drop_cols, stratify_variable)
  }
  
  x_new <- newData[, setdiff(names(newData), drop_cols), drop = FALSE]
  y_true <- newData[[response_variable]]
  
  if (response_type == "binary") {
    y_true <- as.factor(y_true)
    if (length(levels(y_true)) != 2) stop("Binary truth must have exactly 2 levels.")
    levs <- levels(y_true)
    
    # If not provided, fall back to old convention: 2nd level is positive
    if (is.null(positive_class) || !nzchar(positive_class)) {
      positive_class <- levs[2]
    }
    
    # Case-insensitive match to actual truth levels
    pos_match <- which(tolower(levs) == tolower(positive_class))
    if (length(pos_match) != 1) {
      stop("positive_class '", positive_class, "' not found in truth levels: ",
           paste(levs, collapse = ", "))
    }
    positive_class <- levs[pos_match]
    
    truth <- as.numeric(y_true == positive_class)
    
    prob <- caret::predict.train(fit, newdata = x_new, type = "prob")
    if (!(positive_class %in% colnames(prob))) {
      stop("Probability column for positive_class='", positive_class, "' not found. ",
           "Prob cols: ", paste(colnames(prob), collapse = ", "))
    }
    
    p <- prob[, positive_class]
    return(list(pred = as.numeric(p), truth = truth))
  }
  
  if (response_type == "continuous") {
    truth <- as.numeric(y_true)
    p <- caret::predict.train(fit, newdata = x_new)
    return(list(pred = as.numeric(p), truth = truth))
  }
  
  stop("Unsupported response_type: ", response_type)
}

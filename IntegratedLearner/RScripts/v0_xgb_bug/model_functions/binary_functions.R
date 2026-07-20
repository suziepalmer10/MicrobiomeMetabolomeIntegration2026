#notes on inputs for these functions are working correctly
#Accuracy assumes that y_true and predictions are categorical with the same classes.
#AUROC assumes that y_true is a binary vector (e.g., 0/1) and predictions are probability
#scores or predicted probabilities.
#Kappa assumes that y_true and predictions are factors or categorical variables.


# Convert truth labels to 0/1 robustly.
# If y_true is factor/character with 2 levels, maps second level -> 1, first -> 0
truth_to_01 <- function(y_true) {
  if (is.factor(y_true)) {
    lv <- levels(y_true)
    if (length(lv) != 2) stop("Binary truth must have exactly 2 levels.")
    return(as.numeric(y_true == lv[2]))
  }
  if (is.character(y_true)) {
    lv <- sort(unique(y_true))
    if (length(lv) != 2) stop("Binary truth must have exactly 2 unique values.")
    return(as.numeric(y_true == lv[2]))
  }
  # numeric/integer/logical
  y <- as.numeric(y_true)
  # if it's already 0/1, keep
  u <- unique(y[is.finite(y)])
  if (length(u) == 2 && all(sort(u) == c(0, 1))) return(y)
  # otherwise, interpret "higher" value as positive
  if (length(u) == 2) return(as.numeric(y == max(u)))
  stop("Binary truth must be 0/1 or a 2-level factor/character.")
}

accuracy_calculation <- function(y_true, predictions, threshold = 0.5) {
  y <- truth_to_01(y_true)
  pred_class <- ifelse(predictions > threshold, 1, 0)
  mean(pred_class == y, na.rm = TRUE)
}

auroc_calculation <- function(y_true, predictions) {
  y <- truth_to_01(y_true)
  if (length(unique(y[is.finite(y)])) < 2) return(NA_real_)
  roc_obj <- pROC::roc(response = y, predictor = predictions, quiet = TRUE, direction = "<")
  as.numeric(pROC::auc(roc_obj))
}

kappa_calculation <- function(y_true, predictions, threshold = 0.5) {
  y <- truth_to_01(y_true)
  pred_class <- ifelse(predictions > threshold, 1, 0)
  pred_df <- data.frame(y_true_val = y, predictions_val = pred_class)
  irr::kappa2(pred_df, weight = "unweighted")$value
}

auprc_calculation <- function(y_true, predictions) {
  y <- truth_to_01(y_true)
  if (length(unique(y[is.finite(y)])) < 2) return(NA_real_)
  pos_scores <- predictions[y == 1]
  neg_scores <- predictions[y == 0]
  if (length(pos_scores) == 0 || length(neg_scores) == 0) return(NA_real_)
  PRROC::pr.curve(scores.class0 = pos_scores,
                  scores.class1 = neg_scores,
                  curve = FALSE)$auc.integral
}

sensitivity_calculation <- function(y_true, predictions, threshold = 0.5) {
  y <- truth_to_01(y_true)
  pred_class <- ifelse(predictions > threshold, 1, 0)
  TP <- sum(pred_class == 1 & y == 1)
  FN <- sum(pred_class == 0 & y == 1)
  if ((TP + FN) == 0) return(NA_real_)
  TP / (TP + FN)
}

specificity_calculation <- function(y_true, predictions, threshold = 0.5) {
  y <- truth_to_01(y_true)
  pred_class <- ifelse(predictions > threshold, 1, 0)
  TN <- sum(pred_class == 0 & y == 0)
  FP <- sum(pred_class == 1 & y == 0)
  if ((TN + FP) == 0) return(NA_real_)
  TN / (TN + FP)
}

precision_calculation <- function(y_true, predictions, threshold = 0.5) {
  y <- truth_to_01(y_true)
  pred_class <- ifelse(predictions > threshold, 1, 0)
  TP <- sum(pred_class == 1 & y == 1)
  FP <- sum(pred_class == 1 & y == 0)
  if ((TP + FP) == 0) return(NA_real_)
  TP / (TP + FP)
}

f1_calculation <- function(y_true, predictions, threshold = 0.5) {
  prec <- precision_calculation(y_true, predictions, threshold)
  rec  <- sensitivity_calculation(y_true, predictions, threshold)
  if (is.na(prec) || is.na(rec) || (prec + rec) == 0) return(NA_real_)
  2 * prec * rec / (prec + rec)
}

balanced_accuracy_calculation <- function(y_true, predictions, threshold = 0.5) {
  sens <- sensitivity_calculation(y_true, predictions, threshold)
  spec <- specificity_calculation(y_true, predictions, threshold)
  if (is.na(sens) || is.na(spec)) return(NA_real_)
  (sens + spec) / 2
}


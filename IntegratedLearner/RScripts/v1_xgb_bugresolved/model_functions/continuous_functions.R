# CALCULATING R^2 VALUES
r2_calculation <- function(y_true, predictions) {
  TSS <- sum((y_true - mean(y_true))^2)
  RSS <- sum((y_true - predictions)^2)
  R2 <- 1 - (RSS / TSS)
  return(R2)
}

# CALCULATING RMSE VALUES
rmse_calculation <- function(y_true, predictions) {
  MSE <- mean((y_true - predictions)^2)
  RMSE <- sqrt(MSE)
  return(RMSE)
}

# CALCULATING MAE VALUES
mae_calculation <- function(y_true, predictions) {
  MAE <- mean(abs(y_true - predictions))
  return(MAE)
}


# MAPE calculation
mape_calculation <- function(y_true, predictions) {
  y_true <- as.numeric(y_true)
  predictions <- as.numeric(predictions)
  # Avoid division by zero — remove or handle zeros in y_true
  nonzero_idx <- which(y_true != 0)
  if (length(nonzero_idx) == 0) {
    warning("All true values are zero; MAPE is undefined.")
    return(NA)
  }
  mape_val <- mean(abs((y_true[nonzero_idx] - predictions[nonzero_idx]) / y_true[nonzero_idx])) * 100
  return(mape_val)
}

# Median Absolute Error calculation
median_absolute_error_calculation <- function(y_true, predictions) {
  y_true <- as.numeric(y_true)
  predictions <- as.numeric(predictions)
  
  medae <- median(abs(y_true - predictions))
  return(medae)
}










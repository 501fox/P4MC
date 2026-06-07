#' Get Prediction for P4MC 
#'
#' This is an optimized version for multi-model evaluation. 
#' It skips the expensive MCMC sampling and directly accepts pre-computed samples.
#'
#' @param j The index of the test sample.
#' @param test_data The test data.
#' @param samples_y0 Pre-computed MCMC samples assuming class 0.
#' @param samples_y1 Pre-computed MCMC samples assuming class 1.
#' @param fit_model The trained classifier model.
#' @param classifier The type of classifier ('logistic', 'rf', 'svm', 'xgboost').
#' @param pi0 Prior probability for class 0. Default is 0.5.
#' @param pi1 Prior probability for class 1. Default is 0.5.
#' @return A numeric value representing the predicted probability of class 0.

get_prediction_P4MC_LPD <- function(j, test_data, samples_y0, samples_y1, fit_model, classifier = 'logistic', pi0 = 0.5, pi1 = 0.5) {
  
  # 1. Dynamically extract preoperative variables (X1) for patient j
  pre_cols <- grep("^X1_", names(test_data), value = TRUE)
  new_X1 <- test_data[j, pre_cols, drop = FALSE]
  
  # 2. Get number of MCMC samples (S)
  S <- nrow(samples_y0) 
  
  # 3. Construct prediction data frames (R_X2 + X1) by repeating patient data S times
  input0 <- data.frame('R_X2' = samples_y0$RX2_sampled)
  input0 <- cbind(input0, new_X1[rep(1, S), ])
  rownames(input0) <- NULL
  
  input1 <- data.frame('R_X2' = samples_y1$RX2_sampled)
  input1 <- cbind(input1, new_X1[rep(1, S), ])
  rownames(input1) <- NULL
  
  # 4. Helper function to extract class 1 probabilities for different models
  get_prob_class1 <- function(model, newdata, clf) {
    if (clf == 'logistic') {
      preds <- predict(model, newdata, type="response")
      if (any(preds < 0 | preds > 1, na.rm = TRUE)) {
        preds <- plogis(preds)
      }
      return(preds)
    } else if (clf == 'rf') {
      return(predict(model, newdata, type="prob")[, 2])
    } else if (clf == 'svm') {
      pred_obj <- predict(model, newdata, probability=TRUE)
      return(attr(pred_obj, "probabilities")[, "1"])
    } else if (clf == 'xgboost') {
      # Ensure column order matches XGBoost training data perfectly
      expected_cols <- model$feature_names
      newdata <- newdata[, expected_cols, drop = FALSE]
      newdata_x <- xgboost::xgb.DMatrix(data = data.matrix(newdata))
      return(predict(model, newdata_x))
    }
  }
  
  # 5. Batch predict probabilities for sampled points
  prob_1_under_input0 <- get_prob_class1(fit_model, input0, classifier)
  prob_1_under_input1 <- get_prob_class1(fit_model, input1, classifier)
  
  # 6. Bayesian integration inference using class priors
  p0 <- pi0 * mean(1 - prob_1_under_input0)
  p1 <- pi1 * mean(prob_1_under_input1)
  
  # 7. Return normalized probability for class 0
  return(p0 / (p0 + p1))
}
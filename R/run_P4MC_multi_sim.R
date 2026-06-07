#' Run P4MC Evaluation for Multiple Classifiers
#'
#' @param datasets A data.frame containing the full dataset.
#' @param k The number of folds. Default is 10.
#' @param S The iterations of MCMC sampling. Default is 200.
#' @param grid The number of grid points for KDE. Default is 25.
#' @param seed Random seed. Default is 123.
#' @param calculate_auc Logical. Whether to calculate AUC. Default is FALSE.
#' @param classifiers A vector of classifiers to evaluate.
#' @return A list containing the averaged summary table and raw predictions.
#' @export

run_P4MC_multi_sim <- function(datasets, k = 10, S = 200, grid = 25,
                               seed = 123,
                               calculate_auc = FALSE,
                               classifiers = c('logistic', 'svm', 'rf', 'xgboost')) {
  
  if (!requireNamespace("ks", quietly = TRUE)) stop("Package 'ks' required.")
  if (!requireNamespace("caret", quietly = TRUE)) stop("Package 'caret' required.")
  if (calculate_auc && !requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required when calculate_auc = TRUE.")
  }
  
  set.seed(seed)
  folds <- caret::createMultiFolds(y = datasets$Y, k = k, times = 1)
  
  if (length(folds) != k) {
    k <- length(folds)
    warning(paste("Adjusted k to", k))
  }
  
  metrics_store <- list()
  raw_predictions <- data.frame() # Store raw predictions for calibration curves
  
  for (clf in classifiers) {
    metrics_store[[clf]] <- list(
      ACC = numeric(k), F1_0 = numeric(k), F1_1 = numeric(k), 
      MACRO = numeric(k), SENS = numeric(k), SPEC = numeric(k), 
      BRIER = numeric(k), AUC = numeric(k)
    )
  }
  
  for (i in seq_len(k)) {
    fold_start_time <- Sys.time()
    cat(sprintf("  |-- [Rep: %d | Fold: %2d/%d] Starting feature extraction... ", seed, i, k))
    
    train_idx <- folds[[i]]
    train <- datasets[train_idx, ]
    test  <- datasets[-train_idx, ]
    
    pre_cols <- grep("^X1_", names(datasets), value = TRUE)
    
    train_means <- apply(train[, pre_cols], 2, mean, na.rm = TRUE)
    train_sds   <- apply(train[, pre_cols], 2, sd, na.rm = TRUE)
    
    test_scaled <- test
    test_scaled[, pre_cols] <- scale(test[, pre_cols], center = train_means, scale = train_sds)
    
    marg_base <- dimreduce_and_logit_sim(train = train, isscale = TRUE, classifier = classifiers[1])
    marg0 <- marg_base[[1]]
    marg1 <- marg_base[[2]]
    
    kde_cols <- c("R_X2", pre_cols[1:3]) 
    dat_y0 <- marg0[, kde_cols]
    dat_y1 <- marg1[, kde_cols]
    
    shrink_factor <- 1.0  
    
    H_y0 <- if(nrow(dat_y0) > ncol(dat_y0)) ks::Hpi(x = dat_y0) * shrink_factor else diag(apply(dat_y0, 2, var, na.rm=TRUE)) * shrink_factor
    H_y1 <- if(nrow(dat_y1) > ncol(dat_y1)) ks::Hpi(x = dat_y1) * shrink_factor else diag(apply(dat_y1, 2, var, na.rm=TRUE)) * shrink_factor
    
    kde_y0 <- ks::kde(dat_y0, H = H_y0, gridsize = rep(grid, 4)) 
    kde_y1 <- ks::kde(dat_y1, H = H_y1, gridsize = rep(grid, 4)) 
    
    cat("KDE built -> Sampling... ")
    
    mcmc_samples_y0 <- vector("list", nrow(test_scaled))
    mcmc_samples_y1 <- vector("list", nrow(test_scaled))
    
    for (j in 1:nrow(test_scaled)) {
      newX1 <- as.numeric(test_scaled[j, pre_cols]) 
      mcmc_samples_y0[[j]] <- p4mc_kde_sample_sim(newX1, target_y = 0, S, grid, kde_y0, kde_y1)
      mcmc_samples_y1[[j]] <- p4mc_kde_sample_sim(newX1, target_y = 1, S, grid, kde_y0, kde_y1)
    }
    
    cat("Fast predicting... ")
    
    for (clf in classifiers) {
      if (clf == classifiers[1]) {
        fit_model <- marg_base[[3]]
      } else {
        # Fixed: Changed to dimreduce_and_logit_sim
        marg_temp <- dimreduce_and_logit_sim(train = train, isscale = TRUE, classifier = clf)
        fit_model <- marg_temp[[3]]
      }
      
      P_score0 <- sapply(1:nrow(test_scaled), function(j) {
        get_prediction_P4MC_sim(j, test_scaled, mcmc_samples_y0[[j]], mcmc_samples_y1[[j]], fit_model, clf)
      })
      P_score1 <- 1 - P_score0
      
      auc_val <- NA
      if (calculate_auc) {
        tryCatch({
          roc_obj <- pROC::roc(response = test$Y, predictor = P_score1,
                               levels = c(0, 1), direction = "<", quiet = TRUE)
          auc_val <- as.numeric(roc_obj$auc)
        }, error = function(e) { auc_val <<- NA })
      }
      metrics_store[[clf]]$AUC[i] <- auc_val
      
      true_labels_numeric <- as.numeric(as.character(test$Y))
      metrics_store[[clf]]$BRIER[i] <- mean((P_score1 - true_labels_numeric)^2, na.rm = TRUE)
      
      # Store true labels and predicted probabilities per fold
      raw_predictions <- rbind(raw_predictions, data.frame(
        Model = paste0("P4MC_", clf),
        Fold = i,
        Observed = true_labels_numeric,
        Predicted = P_score1
      ))
      
      pred_class <- ifelse(P_score1 > P_score0, 1, 0)
      tab <- table(factor(pred_class, levels=c(0,1)), factor(test$Y, levels=c(0,1)))
      
      n00 <- tab[1, 1]; n01 <- tab[1, 2]
      n10 <- tab[2, 1]; n11 <- tab[2, 2]
      
      metrics_store[[clf]]$ACC[i] <- sum(diag(tab)) / sum(tab)
      
      rec1 <- if((n11 + n01) == 0) 0 else n11 / (n11 + n01)
      metrics_store[[clf]]$SENS[i] <- rec1
      prec1 <- if((n11 + n10) == 0) 0 else n11 / (n11 + n10)
      metrics_store[[clf]]$F1_1[i] <- if((prec1 + rec1) == 0) 0 else 2 * (prec1 * rec1) / (prec1 + rec1)
      
      rec0 <- if((n00 + n10) == 0) 0 else n00 / (n00 + n10)
      metrics_store[[clf]]$SPEC[i] <- rec0
      prec0 <- if((n00 + n01) == 0) 0 else n00 / (n00 + n01)
      metrics_store[[clf]]$F1_0[i] <- if((prec0 + rec0) == 0) 0 else 2 * (prec0 * rec0) / (prec0 + rec0)
      
      metrics_store[[clf]]$MACRO[i] <- (metrics_store[[clf]]$F1_1[i] + metrics_store[[clf]]$F1_0[i]) / 2
    } 
    
    fold_duration <- difftime(Sys.time(), fold_start_time, units = "secs")
    cat(sprintf("Done (%.1f s)\n", as.numeric(fold_duration)))
    
  } # End of folds
  
  final_results <- data.frame()
  for (clf in classifiers) {
    clf_df <- data.frame(
      Replication = seed,
      Model       = paste0("P4MC_", clf),
      Accuracy    = mean(metrics_store[[clf]]$ACC, na.rm=TRUE),
      F1_Class0   = mean(metrics_store[[clf]]$F1_0, na.rm=TRUE),
      F1_Class1   = mean(metrics_store[[clf]]$F1_1, na.rm=TRUE),
      Macro_F1    = mean(metrics_store[[clf]]$MACRO, na.rm=TRUE),
      Sensitivity = mean(metrics_store[[clf]]$SENS, na.rm=TRUE),
      Specificity = mean(metrics_store[[clf]]$SPEC, na.rm=TRUE),
      Brier_Score = mean(metrics_store[[clf]]$BRIER, na.rm=TRUE),
      AUC         = mean(metrics_store[[clf]]$AUC, na.rm=TRUE)
    )
    final_results <- rbind(final_results, clf_df)
  }
  
  # Return summary metrics and raw predictions
  return(list(SummaryTable = final_results, RawPreds = raw_predictions))
}
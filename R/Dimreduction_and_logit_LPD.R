#' Dimensionality Reduction and Classifier Training (with Latent Space SMOTE)
#'
#' @param train A data.frame with features and binary outcome `Y`.
#' @param isscale Whether to standardize X1 and X2. Default is TRUE.
#' @param classifier The type of classifier ('logistic', 'svm', 'rf', 'xgboost').
#' @return A list with two marginal distributions and the fitted model object.

dimreduce_and_logit_LPD <- function(train, isscale = TRUE, classifier = 'logistic') {
  stopifnot("Y" %in% colnames(train))
  
  # Dynamically extract feature columns
  pre_cols <- grep("^X1_", colnames(train), value = TRUE)   
  intra_cols <- grep("^X2_", colnames(train), value = TRUE) 
  
  dim_reduce <- function(data){
    X1 <- data[, pre_cols]   
    X1 <- apply(X1, 2, function(x) as.numeric(as.character(x))) 
    X2 <- data[, intra_cols] 
    
    if (isscale) {
      X1 <- scale(X1)
      X2 <- scale(X2)
    }
    
    # Kernel PCA feature extraction on X1
    sigma_est <- kernlab::sigest(X1)[2]
    kpc <- kernlab::kpca(X1, kernel = "rbfdot", kpar = list(sigma = as.numeric(sigma_est)), features = 1)
    RX1 <- as.matrix(kernlab::rotated(kpc)[, 1, drop = FALSE]) 
    
    # MRCE and residual space projection
    fit <- MRCE::mrce(X=as.matrix(RX1), Y=as.matrix(X2), lam1=0.1, lam2=10^(-4), method="single")
    X2h <- as.matrix(RX1) %*% fit$Bhat
    sig_fit <- crossprod(X2h) / nrow(X1)
    sig_hat <- crossprod(as.matrix(X2)) / nrow(X1)
    sig_res <- sig_hat - sig_fit
    
    tao_0 <- eigen(sig_res)$vectors
    beta_hat <- t(tao_0) %*% t(as.matrix(X2)) %*% as.matrix(RX1) %*% solve(crossprod(as.matrix(RX1)))
    
    # Final extracted intraoperative latent feature: R(X2)
    R_X2 <- as.matrix(X2) %*% tao_0 %*% beta_hat 
    return(R_X2 = R_X2)
  }
  
  # Extract manifold grouped by outcome Y
  R_X2_y0 <- dim_reduce(train[train$Y == 0, ])
  R_X2_y1 <- dim_reduce(train[train$Y == 1, ])
  
  # Combine extracted intraoperative features
  R_X2_all <- rbind(
    data.frame(V1 = R_X2_y0[, 1], index = as.numeric(rownames(R_X2_y0))),
    data.frame(V1 = R_X2_y1[, 1], index = as.numeric(rownames(R_X2_y1)))
  )
  R_X2_all <- R_X2_all[order(R_X2_all$index), ]
  
  # Combine purified intraoperative feature R(X2) with original preoperative features X1
  X1_selected <- train[, pre_cols]
  dat_logit <- cbind(R_X2_all$V1, X1_selected)
  colnames(dat_logit)[1] <- "R_X2" 
  
  if (isscale) {
    dat_logit <- as.data.frame(scale(dat_logit))
  }
  
  # =========================================================
  # Extract true marginal distributions 
  # (KDE must use real data, not SMOTE synthetic samples!)
  # =========================================================
  marg0 <- dat_logit[R_X2_all$index %in% rownames(R_X2_y0), ]
  marg1 <- dat_logit[R_X2_all$index %in% rownames(R_X2_y1), ]
  
  # Construct training data for the classifier
  dat_train <- dat_logit
  dat_train$Y <- as.factor(train$Y)
  
  # =========================================================
  # Apply SMOTE in the latent space 
  # (Generated synthetic samples are ONLY fed to the classifier)
  # =========================================================
  if (exists("apply_smote")) {
    dat_train <- apply_smote(dat_train, target_col = "Y", k_neighbors = 5)
  }
  
  # --- Train classifiers using the balanced dataset ---
  if (classifier == 'svm') {
    fit_model <- e1071::svm(Y~., data=dat_train, kernel="linear", cost=1, probability=TRUE)
  } else if (classifier == 'rf') {
    fit_model <- randomForest::randomForest(Y ~ ., data=dat_train, ntree=100)
  } else if (classifier == 'xgboost') {
    label_vec <- as.numeric(as.character(dat_train$Y))
    x_matrix <- data.matrix(dat_train[, setdiff(names(dat_train), "Y")])
    x_dmat <- xgboost::xgb.DMatrix(data = x_matrix, label = label_vec)
    param <- list(
      objective = "binary:logistic", 
      eval_metric = "logloss",
      max_depth = 4, 
      min_child_weight = 3, 
      subsample = 0.8,          
      colsample_bytree = 0.5, 
      eta = 0.1, 
      gamma = 0.5, 
      nthread = 1
    )
    fit_model <- xgboost::xgb.train(param, x_dmat, nrounds = 50, verbose=0)
  } else {
    fit_model <- glm(Y ~ ., data=dat_train, family=binomial(link = 'logit'))
  }
  
  list(marg0=marg0, marg1=marg1, fit_model=fit_model)
}
#' P4MC KDE Sampling Function
#'
#' @param pre_data The preoperative variable vector of the sampled patient (e.g., X1_1, X1_2, X1_3).
#' @param target_y The target class label (0 or 1), determining whether to use kde_y0 or kde_y1.
#' @param S The number of MCMC sampling iterations.
#' @param grid The number of grid points for KDE on the first dimension (R_X2).
#' @param kde_y0 The KDE object for class Y=0.
#' @param kde_y1 The KDE object for class Y=1.
#'
#' @return A data.frame with two columns: RX2_sampled and P.

p4mc_kde_sample_LPD <- function(pre_data, target_y, S, grid, kde_y0, kde_y1){
  
  # 1. Select the corresponding KDE object based on target class
  kde <- if (target_y == 0) kde_y0 else kde_y1
  x <- kde$eval.points[[1]]
  
  # 2. Match condition slice indices (Using the first three preoperative dimensions)
  ind1 <- which.min(abs(kde$eval.points[[2]] - pre_data[1]))
  ind2 <- which.min(abs(kde$eval.points[[3]] - pre_data[2]))
  ind3 <- which.min(abs(kde$eval.points[[4]] - pre_data[3]))
  
  # 3. Extract the conditional probability density y at the matched slice
  y <- kde$estimate[1:grid, ind1, ind2, ind3]
  
  # Fallback: If the slice is extremely sparse, use marginal density to prevent B-spline failure
  if (all(is.na(y)) || sum(y, na.rm = TRUE) <= 1e-8) {
    y <- apply(kde$estimate[1:grid, , , ], 1, mean, na.rm = TRUE)
  }
  
  y[is.na(y)] <- 0
  
  # ==========================================================
  # B-spline Regression Fit
  # ==========================================================
  err <- sapply(3:10, function(k){
    fit_temp <- lm(y ~ splines::bs(x, degree = k))
    sum(abs(y - predict(fit_temp, data.frame(x))), na.rm = TRUE)
  })
  
  best_df <- which.min(err) + 2
  fit <- lm(y ~ splines::bs(x, df = best_df, degree = 3))
  
  # ==========================================================
  # Metropolis-Hastings Sampling
  # ==========================================================
  theta <- x[which.max(y)]
  
  # Determine truncation range using the top 20 points with the highest probabilities
  top_x <- x[order(y, decreasing = TRUE)[1:min(20, length(x))]]
  min_val <- min(top_x, na.rm = TRUE)
  max_val <- max(top_x, na.rm = TRUE)
  
  # Fallback: Prevent truncation error if min equals max
  if (is.na(min_val) || is.na(max_val) || min_val >= max_val) {
    min_val <- min(x, na.rm = TRUE)
    max_val <- max(x, na.rm = TRUE)
  }
  
  RX2_sampled <- numeric(S)
  
  for (s in seq_len(S)) {
    
    # Proposal distribution: Truncated Normal
    theta.star <- EnvStats::rnormTrunc(
      n = 1, 
      mean = theta, 
      sd = 0.6,
      min = min_val, 
      max = max_val
    )
    
    # Predict probability densities for proposal and current points using B-spline
    p_theta_star <- as.numeric(predict(fit, data.frame(x = theta.star)))
    p_theta      <- as.numeric(predict(fit, data.frame(x = theta)))
    
    # Calculate M-H acceptance ratio
    if (is.na(p_theta_star) || p_theta_star <= 0) {
      r <- 0
    } else if (is.na(p_theta) || p_theta <= 0) {
      r <- 1 # Force accept if current point is invalid
    } else {
      q_forward <- EnvStats::dnormTrunc(
        theta.star, 
        mean = theta, 
        sd = 0.6,
        min = min_val, 
        max = max_val
      )
      
      q_backward <- EnvStats::dnormTrunc(
        theta, 
        mean = theta.star, 
        sd = 0.6,
        min = min_val, 
        max = max_val
      )
      
      if (is.na(q_forward) || q_forward <= 0 ||
          is.na(q_backward) || q_backward <= 0) {
        r <- 0
      } else {
        r <- (p_theta_star * q_backward) / (p_theta * q_forward)
      }
    }
    
    # Determine acceptance
    if (!is.na(r) && runif(1) < min(1, r)) {
      theta <- theta.star
    }
    
    RX2_sampled[s] <- theta
  }
  
  P <- as.numeric(predict(fit, data.frame(x = RX2_sampled)))
  P[is.na(P) | P < 0] <- 0
  
  return(data.frame(RX2_sampled = RX2_sampled, P = P))
}
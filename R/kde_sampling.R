#' Internal KDE Sampling Function
#'
#' This function performs posterior sampling of a latent variable based on conditional KDE slices.
#' It uses B-spline regression to approximate the conditional density and employs Metropolis-Hastings
#' sampling to draw samples from the approximated posterior.
#'
#' @param newX2 A numeric vector or 1-row matrix with 4 elements. The first three elements are used for
#' conditional slicing of the KDE, and the fourth determines which KDE object to use.
#' @param S The iterations of MCMC sampling.
#' @param grid The number of grid points to extract from the KDE in the first dimension.
#' @param kde_d2 KDE object used when \code{newX2[4] == 2}.
#' @param kde_d3 KDE object used when \code{newX2[4] == 3}.
#'
#' @return A \code{data.frame} with two columns:
#' \describe{
#'   \item{RX1i}{A numeric vector of sampled values (length S).}
#'   \item{P}{Predicted density values from the fitted B-spline model.}
#' }
#'
internal_kde_sample <- function(newX2, S, grid, kde_d2, kde_d3){
  kde <- if (newX2[4] == 2) kde_d2 else kde_d3
  x <- kde$eval.points[[1]]

  ind1 <- which.min(abs(kde$eval.points[[2]] - newX2[ ,1]))
  ind2 <- which.min(abs(kde$eval.points[[3]] - newX2[ ,2]))
  ind3 <- which.min(abs(kde$eval.points[[4]] - newX2[ ,3]))
  y <- kde$estimate[1:grid, ind1, ind2, ind3]

  err <- sapply(1:30, function(k){
    fit <- lm(y ~ splines::bs(x, degree=k))
    sum(abs(y - predict(fit, data.frame(x))))
  })
  k.best <- which.min(err)
  fit <- lm(y ~ splines::bs(x, degree=k.best))

  theta <- x[which.max(y)]
  min_val <- min(x[order(y, decreasing=TRUE)[1:20]])
  max_val <- max(x[order(y, decreasing=TRUE)[1:20]])

  RX1i <- numeric(S)
  set.seed(1111)
  for (s in seq_len(S)) {
    theta.star <- EnvStats::rnormTrunc(n = 1, mean = theta, sd = 0.6, min = min_val, max = max_val) # #proposal theta star

    p_theta_star <- predict(fit, data.frame(x=theta.star))
    p_theta <- predict(fit, data.frame(x=theta))
    if (is.na(p_theta_star) || p_theta_star <= 0) {
      r <- 0
    } else {
      r <- (p_theta_star * dnorm(theta, theta.star, 0.5)) /
        (p_theta * dnorm(theta.star, theta, 0.1))
    }
    if (!is.na(r) && runif(1) < min(1, r)) theta <- theta.star
    RX1i[s] <- theta
  }
  P <- predict(fit, data.frame(x=RX1i))
  data.frame(RX1i=RX1i, P=P)
}

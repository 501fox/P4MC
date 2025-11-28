#' Dimensionality Reduction and Logistic Regression
#'
#' @param train A data.frame with features and binary outcome `grading`.
#' @param isscale Whether to standardize X1 and X2. Default is TRUE.
#' @return A list with two marginal distributions and the fitted glm object.

dimreduce_and_logit <- function(train, isscale = TRUE) {
  stopifnot("grading" %in% colnames(train))

  pre_cols <- c('age','platelet','diameter','differentiation')
  intra_cols <- c('bleeding','numb.nodes.sweeping','numb.nodes.positive')

  dim_reduce <- function(data){
    X2 <- data[, pre_cols]   # preoperative
    X1 <- data[, intra_cols]   # intraoperative
    if (isscale) {
      X2 <- scale(X2)
      X1 <- scale(X1)
    }
    RX2 <- princomp(X2, cor=TRUE)$scores[, 1] # dim = 1
    fit <- MRCE::mrce(X=as.matrix(RX2), Y=as.matrix(X1), lam1=10^(-1.5), lam2=10^(-0.5), method="single")
    X1h <- as.matrix(RX2) %*% fit$Bhat
    sig_fit <- crossprod(X1h) / nrow(X2)
    sig_hat <- crossprod(as.matrix(X1)) / nrow(X2)
    sig_res <- sig_hat - sig_fit
    tao_0 <- eigen(sig_res)$vectors
    beta_hat <- t(tao_0) %*% t(as.matrix(X1)) %*% as.matrix(RX2) %*% solve(crossprod(as.matrix(RX2)))
    RX1_y <- as.matrix(X1) %*% tao_0 %*% beta_hat
    return(RX1_y=RX1_y)
  }

  RX1_y0 <- dim_reduce(train[train$grading == 0, ])
  RX1_y1 <- dim_reduce(train[train$grading == 1, ])

  RX1 <- rbind(
    data.frame(V1 = RX1_y0[, 1], index = as.numeric(rownames(RX1_y0))),
    data.frame(V1 = RX1_y1[, 1], index = as.numeric(rownames(RX1_y1)))
  )
  RX1 <- RX1[order(RX1$index), ]

  X2_selected <- train[, pre_cols]
  dat_logit <- cbind(RX1$V1, X2_selected)
  dat_logit$differentiation <- as.factor(dat_logit$differentiation)
  if (isscale) dat_logit[, 1:4] <- scale(dat_logit[, 1:4])

  logit_fit <- glm(as.factor(train$grading) ~ ., data=dat_logit, family=binomial(link = 'logit'))

  marg0 <- dat_logit[RX1$index %in% rownames(RX1_y0), ]
  marg1 <- dat_logit[RX1$index %in% rownames(RX1_y1), ]
  list(marg0=marg0, marg1=marg1, logit_fit=logit_fit)
}

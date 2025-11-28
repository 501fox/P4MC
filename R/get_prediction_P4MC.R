#' Get prediction for P4MC using KDE and logistic regression
#'
#' @param i The index of the test sample.
#' @param test_data The test data.
#' @param S The iterations of MCMC sampling.
#' @param grid The number of grid points to extract from the KDE in the first dimension.
#' @param kde.y0.d2 KDE estimates for class 2, differentiation = 2.
#' @param kde.y0.d3 KDE estimates for class 3, differentiation = 3.
#' @param kde.y1.d2 KDE estimates for class 2, differentiation = 2 for class 1.
#' @param kde.y1.d3 KDE estimates for class 3, differentiation = 3 for class 1.
#' @param logit_fit The logistic regression model.
#' @return A numeric value representing the predicted probability.

get_prediction_P4MC <- function(i, test_data, S, grid, kde.y0.d2, kde.y0.d3, kde.y1.d2, kde.y1.d3, logit_fit) {

  newX2 = test_data[i, setdiff(names(test_data), 'grading')]

  # caculate RX1
  new_RX1_P.y0 = internal_kde_sample(newX2, S, grid, kde.y0.d2, kde.y0.d3)
  new_RX1_P.y1 = internal_kde_sample(newX2, S, grid, kde.y1.d2, kde.y1.d3)

  # the list of variables in Logistic regression
  input0 = data.frame(
    'RX1$V1' = new_RX1_P.y0$RX1i,
    'age' = as.numeric(rep(newX2[1], S)),
    'platelet' = as.numeric(rep(newX2[2], S)),
    'diameter' = as.numeric(rep(newX2[3], S)),
    'differentiation' = rep(newX2[4], S)
  )

  input1 = data.frame(
    'RX1$V1' = new_RX1_P.y1$RX1i,
    'age' = as.numeric(rep(newX2[1], S)),
    'platelet' = as.numeric(rep(newX2[2], S)),
    'diameter' = as.numeric(rep(newX2[3], S)),
    'differentiation' = rep(newX2[4], S)
  )

  colnames(input0) = c('RX1$V1', 'age', 'platelet', 'diameter', 'differentiation')
  colnames(input1) = c('RX1$V1', 'age', 'platelet', 'diameter', 'differentiation')

  # the prediction results
  p0 = mean(1 - (1 / (exp(-predict(logit_fit, input0)) + 1)))
  p1 = mean(1 / (exp(-predict(logit_fit, input1)) + 1))

  return(p0 / (p0 + p1))
}


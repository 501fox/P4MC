#' Clinical LPD Dataset
#'
#' A clinical dataset containing preoperative and intraoperative features for patients.
#' This dataset is used to train and validate the P4MC model on real-world data.
#'
#' @format A data frame with 8 variables:
#' \describe{
#'   \item{bleeding}{Intraoperative blood loss (numeric).}
#'   \item{numb.nodes.sweeping}{Number of lymph nodes swept/harvested (numeric).}
#'   \item{numb.nodes.positive}{Number of positive lymph nodes (numeric).}
#'   \item{differentiation}{Tumor differentiation grade (e.g.,1, 2, 3) (numeric).}
#'   \item{age}{Patient age in years (numeric).}
#'   \item{platelet}{Platelet count (numeric).}
#'   \item{diameter}{Tumor diameter (numeric).}
#'   \item{grading}{Binary outcome variable (0 or 1).}
#' }
#' @usage data(data_LPD)
"data_LPD"

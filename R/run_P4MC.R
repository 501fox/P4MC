#' Run P4MC evaluation with Comprehensive F1 Scores + Optional AUC
#'
#' @param datasets A data.frame containing the full dataset.
#' @param k The Number of folds. Default is 10.
#' @param S The iterations of MCMC sampling. Default is 200.
#' @param grid The number of grid points. Default is 25.
#' @param seed Random seed. Default is 123.
#' @param csv_path The file path to save the result CSV.
#' @param calculate_auc Logical. Whether to calculate AUC. Default is FALSE.
#' @return A list with detailed scores and the final summary dataframe.
#' @export
run_P4MC <- function(datasets, k = 10, S = 200, grid = 25,
                                      seed = 123,
                                      csv_path = "P4MC_results_comprehensive.csv",
                                      calculate_auc = FALSE) { # 新增开关

  # 1. 检查依赖
  if (!requireNamespace("ks", quietly = TRUE)) stop("Package 'ks' required.")
  if (!requireNamespace("caret", quietly = TRUE)) stop("Package 'caret' required.")
  if (calculate_auc && !requireNamespace("pROC", quietly = TRUE)) {
    stop("Package 'pROC' is required when calculate_auc = TRUE.")
  }

  # 2. 准备数据划分
  cat("==========================================================\n")
  cat(sprintf(">>> Initializing P4MC Evaluation\n"))
  if(calculate_auc) cat(">>> AUC Calculation: ENABLED\n") else cat(">>> AUC Calculation: DISABLED\n")
  cat("==========================================================\n")

  set.seed(seed)
  # 注意：确保 datasets$grading 是 factor 类型，否则 createMultiFolds 可能会有不同行为
  folds <- caret::createMultiFolds(y=datasets$grading, k = k, times = 1)

  if (length(folds) != k) {
    k <- length(folds)
    warning(paste("Adjusted k to", k))
  }

  # 3. 初始化存储容器
  ACC_vec       <- numeric(k)
  F1_class0_vec <- numeric(k)
  F1_class1_vec <- numeric(k)
  F1_macro_vec  <- numeric(k)
  SENS_vec      <- numeric(k)
  SPEC_vec      <- numeric(k)
  AUC_vec       <- numeric(k) # 新增 AUC 存储容器 (默认全是 0)
  TAB           <- vector("list", k)

  # 如果不计算 AUC，填充为 NA 以便区分
  if (!calculate_auc) AUC_vec <- rep(NA, k)

  total_start_time <- Sys.time()

  # 4. 循环 K 折
  for (i in seq_len(k)) {

    # --- 进度提醒 ---
    fold_start_time <- Sys.time()
    cat(sprintf("[Fold %d/%d] Processing... ", i, k))

    # A. 数据切分
    train_idx <- folds[[i]]
    train <- datasets[train_idx, ]
    test  <- datasets[-train_idx, ]

    # B. 预处理
    train$differentiation <- as.numeric(as.character(train$differentiation))
    test$differentiation <- as.numeric(as.character(test$differentiation))

    # C. P4MCc 训练
    marg <- dimreduce_and_logit(train = train, isscale = TRUE)
    marg0 <- marg[[1]]; marg1 <- marg[[2]]; logit_fit <- marg[[3]]

    kde.y0.d2 <- ks::kde(marg0[marg0$differentiation == 2, 1:4], gridsize = rep(grid, 4))
    kde.y0.d3 <- ks::kde(marg0[marg0$differentiation == 3, 1:4], gridsize = rep(grid, 4))
    kde.y1.d2 <- ks::kde(marg1[marg1$differentiation == 2, 1:4], gridsize = rep(grid, 4))
    kde.y1.d3 <- ks::kde(marg1[marg1$differentiation == 3, 1:4], gridsize = rep(grid, 4))

    # D. Test 标准化
    test_scaled <- test[, c('age', 'platelet', 'diameter', 'differentiation', 'grading')]
    test_scaled[, c('age', 'platelet', 'diameter')] <- scale(test_scaled[, c('age', 'platelet', 'diameter')])
    test_scaled$differentiation <- as.factor(test_scaled$differentiation)

    # E. 预测
    P_score0 <- sapply(1:nrow(test_scaled), function(j) {
      get_prediction_P4MC(j, test_scaled, S, grid, kde.y0.d2, kde.y0.d3, kde.y1.d2, kde.y1.d3, logit_fit)
    })
    P_score1 <- 1 - P_score0

    # --- 新增：计算 AUC (带容错保护) ---
    auc_val <- NA
    if (calculate_auc) {
      # 使用 tryCatch 防止因某折缺少某一类样本而导致报错
      tryCatch({
        # 确保 response 是 factor 且 levels 正确
        roc_obj <- pROC::roc(response = test$grading, predictor = P_score1,
                             levels = c(0, 1), direction = "<", quiet = TRUE)
        auc_val <- as.numeric(roc_obj$auc)
      }, error = function(e) {
        # 如果报错（通常是因为只有一种类别），保持 NA
        auc_val <<- NA
      })
      AUC_vec[i] <- auc_val
    }

    # 类别判定
    pred_class <- ifelse(P_score1 > P_score0, 1, 0)

    # F. 生成混淆矩阵
    tab <- table(factor(pred_class, levels=c(0,1)),
                 factor(test$grading, levels=c(0,1)))
    TAB[[i]] <- tab

    # --- 详细指标计算 ---
    n00 <- tab[1, 1]; n01 <- tab[1, 2]
    n10 <- tab[2, 1]; n11 <- tab[2, 2]

    # 1. Accuracy
    ACC_vec[i] <- sum(diag(tab)) / sum(tab)

    # 2. Class 1 (传统正例)
    rec1 <- if((n11 + n01) == 0) 0 else n11 / (n11 + n01)
    SENS_vec[i] <- rec1
    prec1 <- if((n11 + n10) == 0) 0 else n11 / (n11 + n10)
    F1_class1_vec[i] <- if((prec1 + rec1) == 0) 0 else 2 * (prec1 * rec1) / (prec1 + rec1)

    # 3. Class 0 (以 0 为正例)
    rec0 <- if((n00 + n10) == 0) 0 else n00 / (n00 + n10)
    SPEC_vec[i] <- rec0
    prec0 <- if((n00 + n01) == 0) 0 else n00 / (n00 + n01)
    F1_class0_vec[i] <- if((prec0 + rec0) == 0) 0 else 2 * (prec0 * rec0) / (prec0 + rec0)

    # 4. Macro F1
    F1_macro_vec[i] <- (F1_class1_vec[i] + F1_class0_vec[i]) / 2

    # 打印单折耗时及核心指标
    fold_duration <- difftime(Sys.time(), fold_start_time, units = "secs")

    # 构建打印字符串
    log_str <- sprintf("Done %.1fs | Acc:%.3f | MacroF1:%.3f",
                       as.numeric(fold_duration), ACC_vec[i], F1_macro_vec[i])
    if (calculate_auc) {
      log_str <- paste0(log_str, sprintf(" | AUC:%.3f", ifelse(is.na(AUC_vec[i]), 0, AUC_vec[i])))
    }
    log_str <- paste0(log_str, "\n")
    cat(log_str)
  }

  # 5. 结果汇总与 CSV 生成
  cat("\n==========================================================\n")
  cat("Creating Final Report (Format: Mean ± SD)\n")

  fmt <- function(vec) sprintf("%.3f ± %.3f", mean(vec, na.rm=TRUE), sd(vec, na.rm=TRUE))

  # 构建最终的 DataFrame
  final_df <- data.frame(
    Model       = "P4MC",
    Accuracy    = fmt(ACC_vec),
    F1_Class0   = fmt(F1_class0_vec),
    F1_Class1   = fmt(F1_class1_vec),
    Macro_F1    = fmt(F1_macro_vec),
    Sensitivity = fmt(SENS_vec),
    Specificity = fmt(SPEC_vec),
    stringsAsFactors = FALSE
  )

  # 只有当 calculate_auc = TRUE 时，才添加 AUC 列
  if (calculate_auc) {
    final_df$AUC <- fmt(AUC_vec)
  }

  print(final_df)

  # 保存 CSV
  if (!is.null(csv_path)) {
    write.csv(final_df, csv_path, row.names = FALSE)
    cat(sprintf("\nResults successfully saved to: %s\n", csv_path))
  }

  total_duration <- difftime(Sys.time(), total_start_time, units = "mins")
  cat(sprintf("Total Runtime: %.2f mins\n", as.numeric(total_duration)))

  return(list(
    confusion_tables = TAB,
    scores = data.frame(
      Accuracy    = ACC_vec,
      F1_Class0   = F1_class0_vec,
      F1_Class1   = F1_class1_vec,
      Macro_F1    = F1_macro_vec,
      Sensitivity = SENS_vec,
      Specificity = SPEC_vec,
      AUC         = AUC_vec # 如果没算，这里就是 NA
    ),
    final_table = final_df
  ))
}

library(survival)
library(glmnet)

time_col <- "time"
status_col <- "status"

# 1. 严格过滤：只保留生存时间大于 0 的有效样本，并处理缺失值
df_analysis_filtered <- df_analysis[!is.na(df_analysis[[time_col]]) & df_analysis[[time_col]] > 0, ]

# 2. 自动挑选缺失率低于 30% 的数值型特征
all_numeric_cols <- names(df_analysis_filtered)[sapply(df_analysis_filtered, is.numeric)]
candidate_cols <- setdiff(all_numeric_cols, c(time_col, status_col))

missing_rates <- sapply(df_analysis_filtered[, candidate_cols], function(x) mean(is.na(x)))
reasonable_cols <- names(missing_rates[missing_rates < 0.3])

# 3. 提取数据并去除含 NA 的行
sub_df <- df_analysis_filtered[, c(time_col, status_col, reasonable_cols)]
df_clean <- sub_df[complete.cases(sub_df), ]

cat(sprintf("过滤无效时间和缺失值后的有效样本量: %d\n", nrow(df_clean)))

# 4. 构建特征矩阵 X，并剔除方差为 0 的列
X_raw <- as.matrix(df_clean[, reasonable_cols])
valid_cols <- c()
for(col in reasonable_cols) {
  x_vec <- X_raw[, col]
  if(isTRUE(all(is.finite(x_vec))) && isTRUE(var(x_vec) > 0)) {
    valid_cols <- c(valid_cols, col)
  }
}

X <- X_raw[, valid_cols, drop = FALSE]
Y <- Surv(time = df_clean[[time_col]], event = df_clean[[status_col]])

cat(sprintf("最终参与 Cox-Lasso 筛选的有效特征数: %d\n", ncol(X)))

# 5. 安全运行 Cox-Lasso 交叉验证（显式指定 cox.ties 以消除警告）
set.seed(42)
cv_fit <- cv.glmnet(X, Y, family = "cox", alpha = 1, standardize = TRUE, nfolds = 10, cox.ties = "breslow")

# 6. 提取最优 lambda 下的回归系数，并制作特征重要性排名表
best_lambda <- cv_fit$lambda.min
lasso_coefs <- coef(cv_fit, s = "lambda.min")
coef_matrix <- as.matrix(lasso_coefs)

selected_vars <- coef_matrix[coef_matrix[, 1] != 0, , drop = FALSE]

ranking_df <- data.frame(
  Feature = rownames(selected_vars),
  Coefficient = selected_vars[, 1],
  Abs_Coefficient = abs(selected_vars[, 1])
)

ranking_df <- ranking_df[order(ranking_df$Abs_Coefficient, decreasing = TRUE), ]
ranking_df$Rank <- 1:nrow(ranking_df)

# 7. 打印最终的特征排名榜
print("=== 机器学习 Cox-Lasso 特征筛选与排名结果 ===")
print(ranking_df[, c("Rank", "Feature", "Coefficient", "Abs_Coefficient")])

# 8. 自动检查 phenoage0 或 phenoage_advance0 是否上榜
target_vars <- c("phenoage0", "phenoage_advance0")
for(v in target_vars) {
  if(v %in% ranking_df$Feature) {
    p_rank <- ranking_df$Rank[ranking_df$Feature == v]
    cat(sprintf("\n👉 恭喜！'%s' 成功通过特征筛选，在所有变量中排名第 %d 位！\n", v, p_rank))
  } else {
    cat(sprintf("\n❌ 很遗憾，'%s' 的系数被压缩成了 0，在这次筛选中被剔除了。\n", v))
  }
}
library(survival)

# 1. 提取刚才 Lasso 筛选出来的特征名称（假设您的 ranking_df 还在内存中）
# 剔除掉生存时间和状态本身（如果它们在里面）
lasso_vars <- ranking_df$Feature
lasso_vars <- setdiff(lasso_vars, c("time", "status"))

# 2. 动态拼装标准多因素 Cox 回归的公式
# 格式为：Surv(time, status) ~ 变量1 + 变量2 + ...
cox_formula <- as.formula(paste("Surv(time, status) ~", paste(lasso_vars, collapse = " + ")))

# 3. 运行传统的多因素 Cox 回归 (Maximum Likelihood Estimation)
multivar_cox <- coxph(cox_formula, data = df_clean)

# 4. 查看详细的回归结果（包含系数、HR值、95%置信区间以及 p 值）
summary_cox <- summary(multivar_cox)
print(summary_cox)

# 5. 如果你想更直观地提取出各变量的 HR 值和 p 值，可以这样做：
coef_table <- as.data.frame(summary_cox$coefficients)
conf_int <- as.data.frame(summary_cox$conf.int)

result_df <- data.frame(
  Feature = rownames(coef_table),
  Coefficient = coef_table$coef,   HR = conf_int$`exp(coef)`,
  HR_lower = conf_int$`lower .95`,
  HR_upper = conf_int$`upper .95`,
  P_value = coef_table$`Pr(>|z|)`
)

print("=== 传统多因素 Cox 回归最终结果与效应量 (HR) ===")
print(result_df)

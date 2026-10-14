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
# 1. 加载必要的依赖包
library(mclust)       # 用于高斯混合模型 (GMM) 聚类
library(survival)     # 用于生存分析
library(survminer)    # 用于绘制高颜值生存曲线 (KM Curve)
library(dplyr)        # 用于数据处理
library(ranger)       # 用于高效计算特征重要性
library(ggplot2)      # 用于绘图

# -------------------------------------------------------------
# 阶段一：数据准备与核心特征锁定
# -------------------------------------------------------------
# 使用截图中确认存在的准确列名（包含 egfr, phenoage0, bmi 以及对应的 z 标准化轴）
cluster_vars <- c("egfr", "phenoage0", "bmi", "z_feat_alb", "z_feat_creat", "z_feat_glu", "z_feat_crp")
time_col     <- "time"
status_col   <- "status"

# 提取并严格过滤缺失值
required_cols <- c(time_col, status_col, cluster_vars)
df_sub <- df_analysis_filtered[, intersect(names(df_analysis_filtered), required_cols)]
df_cluster_clean <- na.omit(df_sub)

cat(sprintf("👉 成功锁定有效聚类样本量: %d\n", nrow(df_cluster_clean)))

# 矩阵标准化 (Z-score)
X_scaled <- scale(df_cluster_clean[, cluster_vars])


# -------------------------------------------------------------
# 阶段二：运行 GMM (Gaussian Mixture Model) 无监督聚类
# -------------------------------------------------------------
set.seed(2026)
# Mclust 会自动根据贝叶斯信息准则 (BIC) 寻找最优聚类数和模型结构
gmm_model <- Mclust(X_scaled)

# 打印聚类摘要（可查看最优分为几类）
print(summary(gmm_model))

# 将聚类标签安全赋回数据框并转换为因子
df_cluster_clean$Cluster <- as.factor(gmm_model$classification)

# 打印各亚群的绝对样本量与百分比占比
cat("\n=== 各聚类亚群样本量分布 ===\n")
print(table(df_cluster_clean$Cluster))
cat("\n=== 各聚类亚群百分比占比 ===\n")
print(prop.table(table(df_cluster_clean$Cluster)) * 100)


# -------------------------------------------------------------
# 阶段三：生存预后效能检验 (Validation)
# -------------------------------------------------------------
# 构建生存对象
surv_obj <- Surv(time = df_cluster_clean[[time_col]], event = df_cluster_clean[[status_col]])

# 1. 单/多因素 Cox 回归检验亚群间的死亡风险差异
cox_model <- coxph(surv_obj ~ Cluster, data = df_cluster_clean)
print("=== 聚类亚群的 Cox 预后回归结果 ===")
print(summary(cox_model))

# 2. 绘制高颜值 Kaplan-Meier 生存曲线
fit_km <- survfit(surv_obj ~ Cluster, data = df_cluster_clean)

p_km <- ggsurvplot(
  fit_km,
  data = df_cluster_clean,
  pval = TRUE,                  # 自动计算并展示 Log-rank p 值
  risk.table = TRUE,            # 展示下方风险表
  surv.median.line = "hv",      # 展示中位生存时间线
  palette = "jco",              # 采用高阶期刊配色 (JCO 风格)
  ggtheme = theme_classic(),
  title = "Kaplan-Meier Survival Curves by Data-Driven GMM Clusters",
  xlab = "Follow-up Time (Days/Months)",
  ylab = "Survival Probability"
)

print(p_km)


# -------------------------------------------------------------
# 阶段四：审稿人防御策略 —— 使用 ranger 评估特征对“区分亚群”的贡献度
# -------------------------------------------------------------
set.seed(2026)
ranger_model <- ranger(
  formula = Cluster ~ ., 
  data = df_cluster_clean[, c(cluster_vars, "Cluster")], 
  importance = "impurity",
  num.trees = 500
)

# 提取特征重要性得分
importance_scores <- ranger_model$variable.importance
importance_df <- data.frame(
  Feature = names(importance_scores),
  Importance = as.numeric(importance_scores)
) %>% arrange(desc(Importance))

print("=== 各特征在驱动聚类分群中的重要性排序 (Ranger) ===")
print(importance_df)

# 绘制特征重要性条形图
p_imp <- ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(
    title = "Feature Importance in Driving Cluster Separation",
    subtitle = "Evaluated via Ranger Random Forest (Impurity-based)",
    x = "Clinical & Biological Features",
    y = "Variable Importance (Decrease in Impurity)"
  )

print(p_imp)

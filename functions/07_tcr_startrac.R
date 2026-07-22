# AUTO-GENERATED FROM batch6.ipynb.
# Keep this file synchronized with records/function_provenance.csv.

# ===== Original batch6.ipynb cell 669 =====
library(dplyr)

#' 计算 STARTRAC 风格的 Expansion 指数
#' 
#' @param df 包含 clone_id 和分组信息的 Dataframe
#' @param group_col 亚群列名 (如 "cell_type")
#' @param clone_col 克隆列名 (如 "CTstrict" 或 "clone_id")
#' @param patient_col 患者列名 (可选，用于按病人分别统计)
#' 
#' @return 包含每个组 expansion 指数的统计表
calc_startrac_expansion <- function(df, group_col, clone_col, patient_col = NULL) {
  
  # 定义计算 expa 的内部核心公式
  calculate_expa <- function(clones) {
    N <- length(clones) # 总细胞数
    if (N <= 1) return(0) # 只有一个细胞谈不上扩增
    
    # 计算每个克隆的频率
    counts <- table(clones)
    p <- counts / N
    
    # 计算香农熵 (Shannon Entropy)
    entropy <- -sum(p * log(p))
    
    # STARTRAC Expansion 公式: 1 - (H / log(N))
    expa <- 1 - (entropy / log(N))
    return(expa)
  }

  # 按照指定的分组进行计算
  grouping_vars <- c(patient_col, group_col)
  
  res <- df %>%
    filter(!is.na(.data[[clone_col]])) %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(
      Total_Cells = n(),
      Unique_Clones = n_distinct(.data[[clone_col]]),
      STARTRAC_Expa = calculate_expa(.data[[clone_col]]),
      .groups = "drop"
    )
  
  return(res)
}

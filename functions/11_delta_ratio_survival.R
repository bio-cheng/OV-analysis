# AUTO-GENERATED FROM batch6.ipynb.
# Keep this file synchronized with records/function_provenance.csv.

# ===== Original batch6.ipynb cell 188 =====
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

#' 绘制亚群比例变化 (Delta) 与 PFS 相关性图（增加死亡事件标记）
#' 
#' @param meta 元数据框
#' @param target_cell 目标细胞类型
#' @param cell_type_col 细胞类型列名
#' @param time_col 时间点列名
#' @param patient_col 患者 ID 列名
#' @param pfs_col PFS 时间列名 (连续变量)
#' @param group_col 响应分组列名 (用于点上色)
#' @param death_col 死亡/终点事件状态列名 (新增参数)
plot_delta_pfs_correlation <- function(meta, 
                                       target_cell, 
                                       cell_type_col = "cell_type3", 
                                       time_col = "sample_type_rn", 
                                       patient_col = "patient_id",
                                       pfs_col = "pfs_time",
                                       group_col = "pfs_group",
                                       death_col = "is_death") { # <-- 新增参数
  
  # 1. 计算每个样本的比例，同时在 group_by 中保留 death_col
  prop_data <- meta %>%
    filter(!!sym(time_col) %in% c("D0", "D7")) %>%
    group_by(!!sym(patient_col), !!sym(group_col), !!sym(pfs_col), !!sym(death_col), !!sym(time_col)) %>% 
    summarise(
      Total = n(),
      Target = sum(!!sym(cell_type_col) == target_cell),
      Percentage = (Target / Total) * 100,
      .groups = "drop"
    )
  
  # 2. 计算 Delta 值 (D7 - D0) 并保留生存时间与死亡状态
  delta_data <- prop_data %>%
    pivot_wider(id_cols = c(!!sym(patient_col), !!sym(group_col), !!sym(pfs_col), !!sym(death_col)), 
                names_from = !!sym(time_col), 
                values_from = Percentage) %>%
    # 确保只保留配对齐全的患者
    filter(!is.na(D0) & !is.na(D7)) %>%
    mutate(Delta = D7 - D0)
  
  # 3. 绘图：x轴为比例差值，y轴为PFS时间
  p <- ggplot(delta_data, aes(x = Delta, y = !!sym(pfs_col))) +
    # 添加拟合线（线性回归）
    geom_smooth(method = "lm", color = "black", fill = "grey80", alpha = 0.2, linetype = "dashed") +
    
    # 叠加基础散点，颜色区分响应组
    geom_point(aes(color = !!sym(group_col)), size = 4, alpha = 0.8) +
    
    # -------------------------------------------------------------------
    # 【核心修改】关键步骤：为死亡/事件发生（Yes/1/TRUE）的点叠加一层黑色外圈
    # -------------------------------------------------------------------
    geom_point(data = filter(delta_data, !!sym(death_col) %in% c("Yes", "1", 1, TRUE)),
               shape = 1,          # 1号形状为空心圆圈
               size = 5.8,         # 尺寸设为 5.8（比基础点的 size=4 稍大），形成完美外轮廓
               color = "black",    # 轮廓圈颜色
               stroke = 1.2,       # 轮廓圈粗细
               show.legend = FALSE) + 
    
    # 添加相关系数 (Pearson 或 Spearman)
    stat_cor(method = "pearson", size = 7, label.x = 0.3, label.y.npc = "top") +
    
    # 添加患者标签防止点太挤（如需开启可取消注释）
    # geom_text_repel(aes(label = !!sym(patient_col)), size = 4, max.overlaps = 10) +
    
    scale_fill_manual(values = cell_type_colors_updated) +
    scale_color_manual(values = cell_type_colors_updated) +
    
    # 样式设置
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      strip.text = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8)
    ) +
    labs(
      title = paste(target_cell, "Delta vs OS"),
      subtitle = "Delta = % (D7) - % (D0) | Black circle outline indicates death/event",
      x = paste("Delta Percentage of", target_cell),
      y = "OS Time (Months)",
      color = "Group"
    )
  
  return(p)
}

# ===== Original batch6.ipynb cell 190 =====
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

#' 绘制亚群比例比值 (Ratio) 与 PFS 相关性图（增加死亡事件标记）
#'
#' @param meta 元数据框
#' @param target_cell 目标细胞类型
#' @param cell_type_col 细胞类型列名
#' @param time_col 时间点列名
#' @param patient_col 患者 ID 列名
#' @param pfs_col PFS 时间列名 (连续变量)
#' @param group_col 响应分组列名 (用于点上色)
#' @param death_col 死亡/终点事件状态列名
plot_ratio_pfs_correlation <- function(meta,
                                       target_cell,
                                       cell_type_col = "cell_type3",
                                       time_col = "sample_type_rn",
                                       patient_col = "patient_id",
                                       pfs_col = "pfs_time",
                                       group_col = "pfs_group",
                                       death_col = "is_death") {

  # 1. 计算每个样本的比例
  prop_data <- meta %>%
    filter(!!sym(time_col) %in% c("D0", "D7")) %>%
    group_by(!!sym(patient_col), !!sym(group_col), !!sym(pfs_col), !!sym(death_col), !!sym(time_col)) %>%
    summarise(
      Total = n(),
      Target = sum(!!sym(cell_type_col) == target_cell),
      Percentage = (Target / Total) * 100,
      .groups = "drop"
    )

  # 2. 计算 Ratio (D7 / D0) 并保留生存时间与死亡状态
  ratio_data <- prop_data %>%
    pivot_wider(
      id_cols = c(!!sym(patient_col), !!sym(group_col), !!sym(pfs_col), !!sym(death_col)),
      names_from = !!sym(time_col),
      values_from = Percentage
    ) %>%
    # 过滤掉 D0 缺失或 D0 为 0 的患者（避免比率无限大）
    filter(!is.na(D0), !is.na(D7), D0 > 0) %>%
    mutate(Ratio = D7 / D0)

  # 3. 绘图：x轴为 Ratio，y轴为 PFS 时间
  p <- ggplot(ratio_data, aes(x = Ratio, y = !!sym(pfs_col))) +
    # 拟合线（线性回归）
    geom_smooth(method = "lm", color = "black", fill = "grey80", alpha = 0.2, linetype = "dashed") +

    # 基础散点，颜色区分组别
    geom_point(aes(color = !!sym(group_col)), size = 4, alpha = 0.8) +

    # 为发生事件的患者添加黑色空心圆标记
    geom_point(
      data = filter(ratio_data, !!sym(death_col) %in% c("Yes", "1", 1, TRUE)),
      shape = 1, size = 5.8, color = "black", stroke = 1.2,
      show.legend = FALSE
    ) +

    # 相关系数
    stat_cor(method = "pearson", size = 7, label.x.npc = "left", label.y.npc = "top") +

    scale_fill_manual(values = cell_type_colors_updated) +
    scale_color_manual(values = cell_type_colors_updated) +

    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      strip.text = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8)
    ) +
    labs(
      title = paste(target_cell, "Ratio vs PFS"),
      subtitle = "Ratio = % (D7) / % (D0) | Black circle outline indicates event",
      x = paste("Ratio of", target_cell, "(D7 / D0)"),
      y = "PFS Time (Months)",
      color = "Group"
    )

  return(p)
}

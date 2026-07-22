# AUTO-GENERATED PROVENANCE FILE -- DO NOT EDIT
# Source: ../../batch6.ipynb
# Each block was copied verbatim; it may depend on interactive
# objects created outside this extraction. Use the source index and
# README before running individual blocks.


# ===== batch6.ipynb cell 180 (markdown) =====
# ## plot_delta_pfs_correlation

# ===== batch6.ipynb cell 181 (code) =====
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

#' 绘制亚群比例变化 (Delta) 与 PFS 相关性图
#' 
#' @param meta 元数据框
#' @param target_cell 目标细胞类型
#' @param cell_type_col 细胞类型列名
#' @param time_col 时间点列名
#' @param patient_col 患者 ID 列名
#' @param pfs_col PFS 时间列名 (连续变量)
#' @param group_col 响应分组列名 (用于点上色)
plot_delta_pfs_correlation <- function(meta, 
                                       target_cell, 
                                       cell_type_col = "cell_type3", 
                                       time_col = "sample_type_rn", 
                                       patient_col = "patient_id",
                                       pfs_col = "pfs_time",
                                       group_col = "pfs_group") {
  
  # 1. 计算每个样本的比例
  prop_data <- meta %>%
    filter(!!sym(time_col) %in% c("D0", "D7")) %>%
    group_by(!!sym(patient_col), !!sym(group_col), !!sym(pfs_col), !!sym(time_col)) %>%
    summarise(
      Total = n(),
      Target = sum(!!sym(cell_type_col) == target_cell),
      Percentage = (Target / Total) * 100,
      .groups = "drop"
    )
  
  # 2. 计算 Delta 值 (D7 - D0) 并保留 PFS 时间
  delta_data <- prop_data %>%
    pivot_wider(id_cols = c(!!sym(patient_col), !!sym(group_col), !!sym(pfs_col)), 
                names_from = !!sym(time_col), 
                values_from = Percentage) %>%
    # 确保只保留配对齐全的患者
    filter(!is.na(D0) & !is.na(D7)) %>%
    mutate(Delta = D7 - D0)
  
  # 3. 绘图：x轴为比例差值，y轴为PFS时间
  p <- ggplot(delta_data, aes(x = Delta, y = !!sym(pfs_col))) +
    # 添加拟合线（线性回归）
    geom_smooth(method = "lm", color = "black", fill = "grey80", alpha = 0.2, linetype = "dashed") +
    # 叠加散点，颜色区分响应组
    geom_point(aes(color = !!sym(group_col)), size = 4, alpha = 0.8) +
    # 添加相关系数 (Pearson 或 Spearman)
    stat_cor(method = "pearson", size = 7, label.x= 0.3, label.y.npc = "top") +
    # 添加患者标签防止点太挤
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
      subtitle = "Delta = % (D7) - % (D0)",
      x = paste("Delta Percentage of", target_cell),
      y = "OS Time (Months)",
      color = "Group"
    )
  
  return(p)
}

# ===== batch6.ipynb cell 182 (code) =====
options(repr.plot.width = 6, repr.plot.height =6)
# 假设你的 pfs 时间列名为 "pfs_time"
p <- plot_delta_pfs_correlation(group_col = "lm",
  meta = seurat_merge@meta.data, 
  target_cell = "CD4_Treg",
  pfs_col = "os_time"
)
print(p)

# ===== batch6.ipynb cell 183 (code) =====
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

# ===== batch6.ipynb cell 184 (code) =====
options(repr.plot.width = 6, repr.plot.height =6)
# 假设你的 pfs 时间列名为 "pfs_time"
p <- plot_delta_pfs_correlation(group_col = "lm",death_col = "is_death_e",
  meta = seurat_merge@meta.data, 
  target_cell = "CD4_LAG3",
  pfs_col = "os_time"
)
print(p)

# ===== batch6.ipynb cell 185 (code) =====
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

# ===== batch6.ipynb cell 186 (code) =====
options(repr.plot.width = 6, repr.plot.height =6)
plot_ratio_pfs_correlation(seurat_merge@meta.data,
                                       "CD4_Treg",
                                       cell_type_col = "cell_type3",
                                       time_col = "sample_type_rn",
                                       patient_col = "patient_id",
                                       pfs_col = "pfs_time",
                                       group_col = "lm",
                                       death_col = "is_death_e") + ggtitle("CD4_Treg Ratio vs PFS") + labs(y = "PFS Time (Months)")

# ===== batch6.ipynb cell 400 (code) =====
# ---------------- 统一数据准备函数 ----------------
get_clean_prop_data <- function(seurat_obj, group_var, cell_cluster, methods) {
  # 统一只保留计算比例必需的列，防止 distinct() 失效
  meta <- seurat_obj@meta.data %>%
    dplyr::select(orig.ident, patient_id, !!sym(group_var), !!sym(cell_cluster))
  
  if (methods == "part") {
    res <- meta %>% 
      group_by(orig.ident) %>% 
      mutate(sample_cells = n()) %>% 
      group_by(orig.ident, !!sym(cell_cluster)) %>%
      summarise(
        cluster_freq = n(),
        sample_cells = first(sample_cells),
        patient_id = first(patient_id),
        group_val = first(!!sym(group_var)),
        .groups = "drop"
      ) %>% 
      mutate(percent = cluster_freq / sample_cells * 100)
  } else {
    # 如果是 total 模式，假设 n_sample_cells 已经在 meta 里
    res <- seurat_obj@meta.data %>%
      group_by(orig.ident, !!sym(cell_cluster)) %>%
      summarise(
        cluster_freq = n(),
        n_sample_cells = first(n_sample_cells),
        patient_id = first(patient_id),
        group_val = first(!!sym(group_var)),
        .groups = "drop"
      ) %>% 
      mutate(percent = cluster_freq / n_sample_cells * 100)
  }
  colnames(res)[colnames(res) == "group_val"] <- group_var
  return(res)
}

calculate_cluster_diff <- function(seurat_obj, group_var = "new_group", 
                                  cell_cluster = "seurat_clusters",
                                  group1, group2, methods = "part", 
                                  paired = FALSE) {
  
  # 使用统一函数
  results <- get_clean_prop_data(seurat_obj, group_var, cell_cluster, methods)
  
  diff_stats <- results %>% #filter(!orig.ident %in% c("BZ035-02","BZ036-02")) %>% 
    filter(!!sym(group_var) %in% c(group1, group2)) %>%
    group_by(!!sym(cell_cluster)) %>%
    summarise(
      mean_group1 = mean(percent[!!sym(group_var) == group1], na.rm = TRUE),
      mean_group2 = mean(percent[!!sym(group_var) == group2], na.rm = TRUE),
      logFC = log2((mean_group1 + 0.01) / (mean_group2 + 0.01)),
      p_value = {
        d1 <- percent[!!sym(group_var) == group1]
        d2 <- percent[!!sym(group_var) == group2]
        # 核心：完全匹配 t.test 默认参数
        if(paired) {
            # 配对逻辑：按 patient_id 排序对齐
            tmp <- cur_data() %>% filter(!!sym(group_var) %in% c(group1, group2)) %>%
                   select(patient_id, !!sym(group_var), percent) %>%
                   tidyr::pivot_wider(names_from = !!sym(group_var), values_from = percent) %>%
                   drop_na()
            if(nrow(tmp) < 2) NA else t.test(tmp[[group1]], tmp[[group2]], paired = TRUE)$p.value
        } else {
            if(length(d1) < 2 || length(d2) < 2) NA else t.test(d1, d2, var.equal = FALSE)$p.value
        }
      }, .groups = "drop"
    ) %>%
    mutate(comparison = paste(group1, "vs", group2), .before = 1)
  
  return(diff_stats)
}

# ===== batch6.ipynb cell 403 (code) =====
# 确保加载所有必要的 R 包
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr) # stat_cor 需要此包

#' 绘制两个细胞亚群在不同患者水平的比例相关性散点图 (支持按组染色)
#'
#' @param seurat_obj Seurat 对象
#' @param cell_type1 要计算的第一个细胞类型的名称 (如 "CD4_Treg")
#' @param cell_type2 要计算的第二个细胞类型的名称 (如 "CD8_Tex")
#' @param cell_type_var 存储细胞类型的 Metadata 列名 (如 "cell_type")
#' @param patient_var 样本/患者层级的 Metadata 列名 (如 "orig.ident")
#' @param group_var 用于分组染色的 Metadata 列名 (如 "outcome_group")，默认 NULL
#' @param method 相关性计算方法，"pearson" (默认) 或 "spearman"
#'
#' @return 返回一个列表，包含: plot (ggplot对象) 和 data (宽格式数据框)
plot_celltype_correlation <- function(meta, 
                                      cell_type1, 
                                      cell_type2, 
                                      cell_type_var = "cell_type", 
                                      patient_var = "orig.ident", 
                                      group_var = NULL,     # <--- 新增分组参数
                                      method = "pearson") {
  
  # 1. 提取 metadata
  
  # 2. 检查输入的列名和细胞类型是否存在
  if (!cell_type_var %in% colnames(meta)) stop(paste("Column", cell_type_var, "not found in metadata!"))
  if (!patient_var %in% colnames(meta)) stop(paste("Column", patient_var, "not found in metadata!"))
  if (!is.null(group_var) && !group_var %in% colnames(meta)) stop(paste("Column", group_var, "not found!"))
  
  # 3. 确定需要保留的分组变量
  group_vars_to_keep <- patient_var
  if (!is.null(group_var)) {
    group_vars_to_keep <- c(patient_var, group_var)
  }
  
  # 4. 计算每个患者所有亚群的比例 (携带组别信息)
  prop_df <- meta %>%
    # 按患者 (和组别) 计算总细胞数
    group_by(across(all_of(group_vars_to_keep))) %>%
    mutate(Total_Patient_Cells = n()) %>%
    # 再按细胞亚群计算各自的数量
    group_by(across(all_of(c(group_vars_to_keep, cell_type_var, "Total_Patient_Cells")))) %>%
    summarise(Cell_Count = n(), .groups = "drop") %>%
    # 计算比例 (%)
    mutate(Proportion = (Cell_Count / Total_Patient_Cells) * 100)
  
  # 5. 提取并转置为宽矩阵
  wide_df <- prop_df %>%
    filter(.data[[cell_type_var]] %in% c(cell_type1, cell_type2)) %>%
    select(all_of(c(group_vars_to_keep, cell_type_var)), Proportion) %>%
    pivot_wider(names_from = all_of(cell_type_var), 
                values_from = Proportion, 
                values_fill = 0) 
  
  # 6. 绘图：根据是否有 group_var 决定是否映射颜色
  if (is.null(group_var)) {
    # 基础绘图 (无分组)
    p <- ggplot(wide_df, aes(x = .data[[cell_type1]], y = .data[[cell_type2]])) +
      geom_point(size = 3, alpha = 0.8, color = "#2c7bb6") +
      geom_smooth(method = "lm", color = "#d73027", fill = "grey80", alpha = 0.3) +
      stat_cor(method = method, size = 5)
  } else {
    # 分组绘图 (按组别映射 color 和 fill)
    p <- ggplot(wide_df, aes(x = .data[[cell_type1]], y = .data[[cell_type2]], 
                             color = .data[[group_var]], fill = .data[[group_var]])) +
      geom_point(size = 3, alpha = 0.8) +
      # geom_smooth 会自动为每个组单独拟合一条线
      geom_smooth(method = "lm", alpha = 0.2) +
      # stat_cor 会自动为每个组计算相关系数并根据颜色分层显示文本
      stat_cor(method = method, size = 8, show.legend = FALSE)
  }
  
  # 7. 添加统一的主题和标签
  subtitle_text <- paste0("At the ", patient_var, " level (", tools::toTitleCase(method), " correlation)")
  if (!is.null(group_var)) subtitle_text <- paste0(subtitle_text, " grouped by ", group_var)
  
  p <- p + theme_classic() +
    labs(
      title = "Cell Type Proportion Correlation",
      subtitle = subtitle_text,
      x = paste0("% of ", cell_type1, " (of total cells)"),
      y = paste0("% of ", cell_type2, " (of total cells)"),
      color = ifelse(is.null(group_var), "", group_var),
      fill = ifelse(is.null(group_var), "", group_var)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12, color = "black"),
      legend.text = element_text( size = 14),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey95"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
  
  return(list(plot = p, data = wide_df))
}

# ===== batch6.ipynb cell 404 (code) =====
options(repr.plot.width = 6, repr.plot.height =6)
cor_result <- plot_celltype_correlation(
  meta = meta %>% filter(sample_type_rn == "D7"),
  cell_type1 = "CD4_Treg",
  cell_type2 = "CD4_LAG3",
  cell_type_var = "cell_type3", # 你对象中实际包含亚群名字的 metadata 列名
  patient_var = "orig.ident",
#    group_var = "sample_type2",
  method = "spearman" # 如果数据明显不符合正态分布，可改为 "spearman"
)
   
# 打印出图
print(cor_result$plot + labs(x = "% of CD4_Treg "))

# ===== batch6.ipynb cell 415 (code) =====
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

plot_geneset_correlation <- function(seurat_obj, 
                                     gene_list,
                                     gene_set_name = "CustomScore",
                                     subset_col = "sample_type_rn", 
                                     subset_val = "02",
                                     prop_cell_type = "CD8_Tex", 
                                     feature_cell_type = "CD4_Treg", 
                                     cell_type_col = "cell_type3", 
                                     sample_col = "patient_id",
                                     color_col = NULL) {
  
  # 1. 数据预过滤
  message("正在进行子集过滤...")
  sub_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj@meta.data[[subset_col]] %in% subset_val])
  
  # 2. 计算基因集得分
  genes_to_use <- intersect(gene_list, rownames(sub_obj))
  if(length(genes_to_use) == 0) stop("指定的基因集中没有基因存在于对象中！")
  
  sub_obj <- AddModuleScore(sub_obj, features = list(genes_to_use), name = gene_set_name, ctrl = 100)
  score_column_name <- paste0(gene_set_name, "1")
  
  # 3. 计算【比例亚群】占比 (同时保留分组列)
  group_cols <- c(sample_col)
  if(!is.null(color_col)) group_cols <- c(group_cols, color_col)
  
  prop_data <- sub_obj@meta.data %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      Total_Cells = n(),
      Target_Cells = sum(.data[[cell_type_col]] %in% prop_cell_type),
      Proportion = (Target_Cells / Total_Cells) * 100,
      .groups = "drop"
    )
  
  # 4. 计算【特征亚群】平均得分
  feature_data <- sub_obj@meta.data %>%
    filter(.data[[cell_type_col]] %in% feature_cell_type) %>%
    group_by(.data[[sample_col]]) %>%
    summarise(
      Mean_Score = mean(.data[[score_column_name]], na.rm = TRUE),
      Cell_Count = n(),
      .groups = "drop"
    ) %>%
    filter(Cell_Count >= 5)
  
  # 5. 合并数据
  merged_data <- inner_join(prop_data, feature_data, by = sample_col)
  
  # 6. 绘图：核心逻辑修改点
  # 全局 aes 只包含 X 和 Y，这样 stat_cor 只会计算一个总的相关系数
  p <- ggplot(merged_data, aes(x = Proportion, y = Mean_Score)) +
    # 这里的回归线也是基于全局的（即所有点一条线）
    geom_smooth(method = "lm", color = "grey30", fill = "grey80", alpha = 0.2, linetype = "dashed") +
    # 仅在 geom_point 中映射颜色
    {if(!is.null(color_col)) 
      geom_point(aes(color = .data[[color_col]]), size = 5, alpha = 0.8) 
     else 
      geom_point(size = 5, alpha = 0.8, color = "midnightblue")
    } +
    theme_pubr() +
    # stat_cor 放在外面，不受 color 影响，只显示一个 R 和 P
    stat_cor(method = "spearman", size = 6) + 
    labs(
      title = paste("Global Correlation in", subset_val),
#      subtitle = paste0("X: % of ", prop_cell_type, "\nY: ", gene_set_name, " in ", feature_cell_type),
      x = paste(prop_cell_type, "Proportion (%)"),
      y = paste(gene_set_name, "Module Score"),
      color = color_col
    )
  
  if(!is.null(color_col)) {
    p <- p + scale_color_manual(values = cell_type_colors_updated)+
    theme(
      strip.background = element_rect(fill = "white"),       
      strip.text = element_text(size = 14, face = "bold"),   
      axis.text = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 14,color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right", 
    legend.title =  element_blank(),
       legend.text = element_text(size = 14, color = "black"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8) # ggplot2 新版推荐用 linewidth 替代 size
    )
  }
  
  return(list(plot = p, data = merged_data))
}

# ===== batch6.ipynb cell 416 (code) =====
# 1. 定义从图片中提取的最新特征基因集
cd8_features <- list(
  CytoxScore = c("KLRD1","KLRC1","FGFBP2",
  "CX3CR1","ZEB2","TBX21",
  "GZMH","PRF1"),
  
  ExhaustionScore = c(
  # 抗原经历
  "PDCD1","ENTPD1","TOX","LAG3","TIGIT",

  # 趋化/肿瘤相关
  "CXCR3","CXCR6","CXCR5",

  # 效应功能
  "GZMB","PRF1","IFNG",

  # 干性/前体
  "TCF7","IL7R","SLAMF6"
),
  
  ActivationScore = c("CD69", "CCR7", "CD27", "BTLA", "CD40LG", "IL2RA", "CD3E", 
                      "CD47", "EOMES", "GNLY", "GZMA", "GZMB", "PRF1", "IFNG", 
                      "CD8A", "CD8B", "FASLG", "LAMP1", "LAG3", "CTLA4", 
                      "HLA-DRA", "TNFRSF4", "ICOS", "TNFRSF9", "TNFRSF18")
)

# ===== batch6.ipynb cell 417 (code) =====


# ===== batch6.ipynb cell 418 (code) =====
unique(seurat_merge$cell_type3)

# ===== batch6.ipynb cell 419 (code) =====
immune_programs_targeted$CD8_all__cytotoxicity$genes

# ===== batch6.ipynb cell 420 (code) =====
immune_programs_targeted$CD8_all__cytotoxicity$genes

# ===== batch6.ipynb cell 421 (code) =====
# 定义你的基因集
my_genes <- c("GZMB", "GZMA",  "PRF1", "GNLY", "IFNG", "NKG7", "LAMP1", "ITGA1")

# 调用函数
result <- plot_geneset_correlation(
  seurat_obj = seurat_merge,
  gene_list = immune_programs_targeted$CD4_LAG3__checkpoint_exhaustion$genes,
  gene_set_name = "CD8 Cytotoxicity",
  subset_col = "sample_type_rn",
  subset_val = c("D7"),                 # 只看 02 阶段
  prop_cell_type = "CD4_LAG3",       # X轴：Treg 的占比
  feature_cell_type = immune_programs_targeted$CD4_LAG3__checkpoint_exhaustion$subtypes,     # Y轴：CD8_Tex 的杀伤得分
#feature_cell_type = c( "CD8_Tem"), 
  cell_type_col = "cell_type3",
  sample_col = "patient_id",
    color_col = "lm"
)

# ===== batch6.ipynb cell 422 (code) =====
options(repr.plot.width = 6.5, repr.plot.height = 5)
# 展示图表
print(result$plot + labs(x = "% of CD4_LAG3 in total cells"))

# ===== batch6.ipynb cell 423 (code) =====
ggsave("plot_v4/global_correlation in D7_CD4_LAG3_with CD8 Cytotoxicity.pdf", width = 6.5, height = 5)

# ===== batch6.ipynb cell 453 (code) =====
library(Seurat)
library(ggplot2)
library(dplyr)

# ---------------------------------------------------------
# 第一步：数据准备 (保持不变，但需确保时间点排序)
# ---------------------------------------------------------
meta_data <- seurat_merge@meta.data %>% 
  filter(sample_type_rn != "beforeICI")

# 强制设置时间点顺序，确保面积图横轴逻辑正确
# 请根据你的实际时间点名称修改 levels
time_levels <- c("D0", "D7", "postICI") 
meta_data$sample_type_rn <- factor(meta_data$sample_type_rn, levels = time_levels)

cell_counts <- meta_data %>%
  group_by(sample_type_rn, cell_type_major) %>%
  summarise(Count = n(), .groups = "drop")

cell_freq <- cell_counts %>%
  group_by(sample_type_rn) %>%
  mutate(Proportion = Count / sum(Count) * 100)

# ---------------------------------------------------------
# 第二步：绘制百分比堆积面积图 (参考 image_0d3a61.png)
# ---------------------------------------------------------
options(repr.plot.width = 5, repr.plot.height = 4)

ggplot(cell_freq, 
       aes(x = as.numeric(sample_type_rn), # 转换为数值以实现面积连续连接
           y = Proportion, 
           fill = cell_type_major, 
           group = cell_type_major)) +
  
  # 1. 绘制堆积面积图
  # position = "stack" 因为我们已经计算好了 Proportion (总和100)
  geom_area(alpha = 0.8, size = 0.3, color = "white") + 
  
  # 2. 映射颜色
  scale_fill_manual(values = cell_type_colors_updated) +
  
  # 3. 恢复 X 轴标签
  scale_x_continuous(breaks = 1:length(time_levels), 
                     labels = time_levels,
                     expand = c(0, 0)) + # 去除左右多余间隙
  
  # 4. 规范 Y 轴范围
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.1)) +
  
  # 5. 美化主题
  theme_classic() +
  labs(
    title = "Immune Cell Composition Over Time",
    x = "Time Point",
    y = "Proportion (%)",
    fill = "Cell Type"
  ) +
  theme(
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# ===== batch6.ipynb cell 454 (code) =====
ggsave("plot_v3/Immune Cell Composition Over Time.pdf", width = 5, height = 4)

# ===== batch6.ipynb cell 587 (markdown) =====
# # dotplot for single gene

# ===== batch6.ipynb cell 588 (code) =====
cd4_suppressive_markers <- unique(c(
  "FOXP3", "IL2RA", "IL7R", "CTLA4",
  "TIGIT", "LAG3", "PDCD1", "ICOS",
  "TNFRSF18", "TNFRSF4", "HLA-DRA",
  "ENTPD1", "NT5E",
  "TGFB1", "LRRC32",
  "IL10", "EBI3", "IL12A"
))

# ===== batch6.ipynb cell 589 (code) =====


# ===== batch6.ipynb cell 590 (code) =====


# ===== batch6.ipynb cell 591 (code) =====
library(Seurat)
library(ggplot2)
library(dplyr)

# --- 假设你的显著性结果表是 res_sig ---
# 包含列：gene, Comparison, Pval (或者是你之前算好的显著性标记)
# 我们只取 D7_vs_D0 且显著的部分
sig_genes <- res_sig %>% 
  filter(Comparison == "D7_vs_D0" & p_val < 0.05) %>% 
  pull(gene)

# --- 提取 Seurat DotPlot 数据 ---
# features: 你的基因 list
# idents: 设置为时间点列，如 "sample_type_rn"
p <- DotPlot(cd8_obj, features = cd8_features$ActivationScore, group.by = "sample_type_rn")
plot_data <- p$data # 提取 DotPlot 的计算结果


# ===== batch6.ipynb cell 592 (code) =====


# ===== batch6.ipynb cell 593 (code) =====
# 提取 D7 显著基因 (黑圈)
sig_D7 <- res_sig %>% filter(Comparison == "D7_vs_D0" & p_val < 0.05) %>% pull(gene) %>% unique()

# 提取 postICI 显著基因 (红圈)
sig_post <- res_sig %>% filter(Comparison == "postICI_vs_D0" & p_val < 0.05) %>% pull(gene) %>% unique()

# 在 plot_data 中打标签
plot_data <- plot_data %>%
  mutate(sig_mark = case_when(
    id == "D7" & features.plot %in% sig_D7 ~ "black",
    id == "postICI" & features.plot %in% sig_post ~ "red",
    TRUE ~ "none"
  ))

# ===== batch6.ipynb cell 594 (code) =====
# 设置表达量颜色映射
col_fun <- colorRampPalette(c("#4575B4", "white", "#D73027"))(100)

p_dot <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  # -------------------------------------------------------
  # 1. 基础图层：绘制所有基因的表达点
  # -------------------------------------------------------
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  
  # -------------------------------------------------------
  # 2. 黑色圈图层：标记 D7 显著基因
  # -------------------------------------------------------
  geom_point(data = filter(plot_data, sig_mark == "black"),
             aes(size = pct.exp), 
             shape = 21,      # 21号形状有边框
             color = "black", # 边框为黑色
             fill = NA,       # 内部透明，透出底层颜色
             stroke = 1.2) +  # 边框粗细
  
  # -------------------------------------------------------
  # 3. 红色圈图层：标记 postICI 显著基因
  # -------------------------------------------------------
  geom_point(data = filter(plot_data, sig_mark == "red"),
             aes(size = pct.exp), 
             shape = 21, 
             color = "red",   # 边框为红色
             fill = NA, 
             stroke = 1.2) +
  
  # -------------------------------------------------------
  # 4. 样式美化与坐标轴设置
  # -------------------------------------------------------
# 设置表达量颜色映射
  scale_color_gradientn(colors = col_fun) +
  scale_size_continuous(range = c(1, 8)) +
  theme_bw() +
  labs(title = "Differential Gene Expression Across Timepoints",
       subtitle = "Black circle: D7 vs D0 | Red circle: postICI vs D0",
       x = "Genes", y = "Timepoints", 
       color = "Avg Expression", size = "Percent Expressed") +
  theme(
    # 标题与副标题居中
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
    
    # 坐标轴文字
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(size = 14, angle = 35, hjust = 1, face = "italic", color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    
    # 图例设置
    legend.position = "bottom",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    
    # 面板设置
    panel.grid.major = element_line(color = "grey95"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.background = element_rect(fill = "whi
                                    
te"),       
    strip.text = element_text(size = 14, face = "bold")
  )

# 调整画布尺寸并打印
options(repr.plot.width = 10, repr.plot.height = 4.5)
print(p_dot)

# ===== batch6.ipynb cell 595 (code) =====
ggsave("plot_v4/Fig3d.pdf",width = 12, height = 5)

# ===== batch6.ipynb cell 664 (code) =====
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

# ===== batch6.ipynb cell 677 (code) =====
# 执行计算
startrac_results <- calc_startrac_expansion(
  df = tcr_meta, 
  group_col = "cell_type_major", 
  clone_col = "tcr_id",     # 替换为你实际的克隆 ID 列名
  patient_col = "orig.ident"  # 按病人计算，方便后续画箱线图做统计
)

# ===== batch6.ipynb cell 679 (code) =====
meta_info <- tcr_meta %>% 
  select(orig.ident, sample_type_rn) %>% 
  distinct()

startrac_results <- left_join(startrac_results, meta_info, by = "orig.ident")

#head(startrac_results)
options(repr.plot.width = 3.5, repr.plot.height = 4.5)
ggplot(startrac_results %>% filter(sample_type_rn != "beforeICI"), aes(x = sample_type_rn, y = STARTRAC_Expa, fill = sample_type_rn)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.5) +

  theme_bw() +
 # scale_fill_npg() +
  labs(
    title = "CD4 STARTRAC Expansion",
    y = "Expansion Index",
    x = ""
  ) + scale_fill_manual(values = cell_type_colors_updated) +
  theme(legend.text = element_text(size = 12, color = "black"),
      strip.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
   # legend.position = "none" # 隐藏图例，因为 X 轴已经有名字了
  ) + stat_compare_means( comparisons = list(c("D0", "D7"), c("D7", "postICI"), c("D0", "postICI")))+
    theme(
      strip.background = element_rect(fill = "white"),       
      strip.text = element_text(size = 14, face = "bold"),   
      axis.text = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 14, angle = 30, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none", 
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8) # ggplot2 新版推荐用 linewidth 替代 size
    )

# ===== batch6.ipynb cell 680 (code) =====
write.csv(startrac_results, "plot_v3/boxplot_CD4 STARTRAC Expansion.csv")
ggsave("plot_v3/boxplot_CD4 STARTRAC Expansion.pdf", width = 3.5, height = 4.5)

# ===== batch6.ipynb cell 693 (code) =====
meta_info <- tcr_meta %>% 
  select(orig.ident, group_rn) %>% 
  distinct()

startrac_results <- left_join(startrac_results, meta_info, by = "orig.ident")

#head(startrac_results)

options(repr.plot.width = 4, repr.plot.height = 5)
ggplot(startrac_results, aes(x = group_rn, y = STARTRAC_Expa, fill = group_rn)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.5) +
 # facet_wrap(~cell_type_major, scales = "free_y") + # 按细胞亚群分面
  theme_bw() +
  scale_fill_npg() +
  labs(
    title = "CD4 STARTRAC Expansion",
    y = "Expansion Index",
    x = ""
  ) +
  theme(strip.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
   # legend.position = "none" # 隐藏图例，因为 X 轴已经有名字了
  ) +
    theme(
      strip.background = element_rect(fill = "white"),       
      strip.text = element_text(size = 14, face = "bold"),   
      axis.text = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 14, angle = 30, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none", 
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8) # ggplot2 新版推荐用 linewidth 替代 size
    )+ stat_compare_means( method = "wilcox.test", method.args = list(alternative = "greater"), comparisons = list(c("PFS>6_D0", "PFS<6_D0"), c("PFS>6_D7", "PFS<6_D7"), c("PFS>6_postICI", "PFS<6_postICI")))

# ===== batch6.ipynb cell 694 (code) =====
# 执行计算
startrac_results <- calc_startrac_expansion(
  df = tcr_meta, 
  group_col = "cell_type_major", 
  clone_col = "tcr_id",     # 替换为你实际的克隆 ID 列名
  patient_col = "orig.ident"  # 按病人计算，方便后续画箱线图做统计
)

# ===== batch6.ipynb cell 699 (code) =====
meta_info <- tcr_meta %>% 
  select(orig.ident, sample_type_rn) %>% 
  distinct()
# group_response
startrac_results <- left_join(startrac_results, meta_info, by = "orig.ident")

#head(startrac_results)

# ===== batch6.ipynb cell 700 (code) =====
options(repr.plot.width = 8, repr.plot.height = 8)
ggplot(startrac_results %>% filter(), aes(x = sample_type_rn, y = STARTRAC_Expa, fill = sample_type_rn)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  facet_wrap(~cell_type3, scales = "free_y") + # 按细胞亚群分面
  theme_bw() +
  scale_fill_manual(values = cell_type_colors_updated)+
  labs(
    title = "CD4 subtype STARTRAC Expansion",
    y = "Expansion Index",
    x = ""
  ) +
  theme(strip.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
   # legend.position = "none" # 隐藏图例，因为 X 轴已经有名字了
  )  +
    theme(
      strip.background = element_rect(fill = "white"),       
      strip.text = element_text(size = 14, face = "bold"),   
      axis.text = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 14, angle = 30, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none", 
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8) # ggplot2 新版推荐用 linewidth 替代 size
    )+ stat_compare_means( comparisons = list(c("D0", "D7"), c("D7", "postICI"), c("D0", "postICI")))

# ===== batch6.ipynb cell 701 (code) =====


# ===== batch6.ipynb cell 702 (code) =====
options(repr.plot.width = 5.8, repr.plot.height =4.8)
ggplot(startrac_results %>% filter(cell_type3 %in% c("CD4_Treg", "CD4_LAG3")), aes(x = sample_type_rn, y = STARTRAC_Expa, color = sample_type_rn)) +
  geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 1) +
  facet_wrap(~cell_type3, scales = "free_y") + # 按细胞亚群分面
  theme_bw() +
  scale_color_manual(values = cell_type_colors_updated)+
  labs(
    title = "CD4 subtype STARTRAC Expansion",
    y = "Expansion Index",
    x = ""
  ) +
  theme(strip.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
   # legend.position = "none" # 隐藏图例，因为 X 轴已经有名字了
  )  +
    theme(
      strip.background = element_rect(fill = "white"),       
      strip.text = element_text(size = 14, face = "bold"),   
      axis.text = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 14, angle = 30, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none", 
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8) # ggplot2 新版推荐用 linewidth 替代 size
    )+ stat_compare_means( comparisons = list(c("D0", "D7"), c("D7", "postICI"), c("D0", "postICI")))

# ===== batch6.ipynb cell 703 (code) =====


# ===== batch6.ipynb cell 704 (code) =====
ggsave()

# ===== batch6.ipynb cell 705 (code) =====
write.csv(startrac_results, "plot_v4/Fig3f.csv")
ggsave("plot_v4/boxplot_group_by_pfs_CD4 subtype STARTRAC Expansion.pdf", width = 5.8, height = 4.8)

# ===== batch6.ipynb cell 706 (code) =====
require(ggsci)

# ===== batch6.ipynb cell 707 (code) =====
meta_info <- tcr_meta %>% 
  select(orig.ident,lm_group, lm) %>% 
  distinct()
# group_response
startrac_results <- left_join(startrac_results, meta_info, by = "orig.ident")

#head(startrac_results)
#head(startrac_results)

options(repr.plot.width = 7, repr.plot.height = 5)
ggplot(startrac_results %>% filter(cell_type3 %in% c("CD4_LAG3", "CD4_Treg")), aes(x = lm_group, y = STARTRAC_Expa, color = lm)) +
  geom_boxplot(outlier.shape = NA, alpha = 1) +
  geom_jitter(width = 0.2, alpha = 1) +
  facet_wrap(~cell_type3, scales = "free_y") + # 按细胞亚群分面
  theme_bw() +
  scale_fill_npg() +
  labs(
    title = "STARTRAC Expansion Index",
    y = "Expansion Index",
    x = ""
  ) +scale_color_manual(values = cell_type_colors_updated)+
  theme(strip.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
   # legend.position = "none" # 隐藏图例，因为 X 轴已经有名字了
  )  +
    theme(
      strip.background = element_rect(fill = "white"),       
      strip.text = element_text(size = 14, face = "bold"),   
      axis.text = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 14, angle = 30, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none", 
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8) # ggplot2 新版推荐用 linewidth 替代 size
    )+ stat_compare_means( method.args = list(alternative = "less"), comparisons = list(c("noLM_D0", "LM_D0"), c("noLM_D7", "LM_D7"), c("noLM_postICI", "LM_postICI")))

# ===== batch6.ipynb cell 708 (code) =====
write.csv(startrac_results, "plot_v4/FigS3e.csv")
ggsave("plot_v4/FigS3e.pdf", width = 7, height = 5)

# ===== batch6.ipynb cell 820 (code) =====
ggsave("plot_v4/IFNG_Dynamics_postICI vs D7.pdf", width = 7, height=5.5)

# ===== batch6.ipynb cell 822 (code) =====
PlotTargetCellCorrelationHeatmap <- function(
  seurat_obj,
  target_cell,
  target_stage = NULL,

  cell_col = "cell_type3",
  time_col = "sample_type_rn",
  patient_col = "patient_id",

  # NULL 表示和所有其他细胞类型计算相关性
  other_cells = NULL,

  # NULL 表示用 target_stage 下所有细胞作为分母
  # 如果想在某一大类细胞内算比例，可指定 denominator_cells
  denominator_cells = NULL,

  abundance_mode = c("percentage", "fraction", "count"),

  cor_method = c("spearman", "pearson"),
  p_adjust_method = "BH",
  p_cutoff = 0.05,
  fdr_cutoff = NULL,
  sig_by = c("p", "FDR"),

  min_patients = 5,
  min_total_cells_per_patient = 10,

  exclude_target_from_other = TRUE,

  show_cor_value = FALSE,
  cor_digits = 2,

  cluster_rows = FALSE,
  cluster_cols = FALSE,

  low_color = "#2166AC",
  mid_color = "white",
  high_color = "#B2182B",

  tile_color = "white",
  tile_linewidth = 0.35,

  axis_text_x_size = 11,
  axis_text_y_size = 12,
  sig_text_size = 6,
  cor_text_size = 3.5,

  title = NULL,
  output_prefix = NULL,
  output_pdf = FALSE,
  pdf_width = 10,
  pdf_height = 4.8
) {

  abundance_mode <- match.arg(abundance_mode)
  cor_method <- match.arg(cor_method)
  sig_by <- match.arg(sig_by)

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
  })

  target_cell <- as.character(target_cell)

  if (length(target_cell) < 1) {
    stop("target_cell 至少需要指定一个细胞类型。")
  }

  if (length(target_cell) != 2) {
    warning(
      "当前函数推荐 target_cell 指定两个细胞类型；",
      "你现在指定了 ",
      length(target_cell),
      " 个，函数仍会继续运行。"
    )
  }

  target_cell_label <- paste(target_cell, collapse = " + ")

  # ============================================================
  # 1. 提取并过滤 meta
  # ============================================================
  meta <- seurat_obj@meta.data

  required_cols <- c(
    cell_col,
    patient_col
  )

  if (!is.null(target_stage)) {
    required_cols <- c(required_cols, time_col)
  }

  missing_cols <- setdiff(required_cols, colnames(meta))

  if (length(missing_cols) > 0) {
    stop(
      "seurat_obj@meta.data 中缺少以下列：",
      paste(missing_cols, collapse = ", ")
    )
  }

  message(">>> 提取 Seurat meta.data...")
  message(">>> Target cell(s): ", target_cell_label)

  if (!is.null(target_stage)) {
    message(">>> Target stage: ", paste(target_stage, collapse = " + "))
  }

  meta_use <- meta %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::mutate(
      patient_id_internal = as.character(.data[[patient_col]]),
      cell_type_internal = as.character(.data[[cell_col]])
    ) %>%
    dplyr::filter(
      !is.na(patient_id_internal),
      !is.na(cell_type_internal)
    )

  if (!is.null(target_stage)) {
    meta_use <- meta_use %>%
      dplyr::mutate(
        stage_internal = as.character(.data[[time_col]])
      ) %>%
      dplyr::filter(
        !is.na(stage_internal),
        stage_internal %in% as.character(target_stage)
      )
  }

  if (nrow(meta_use) == 0) {
    stop("过滤 target_stage 后没有可用细胞。")
  }

  # ============================================================
  # 2. 确定 denominator
  # ============================================================
  if (!is.null(denominator_cells)) {
    denominator_cells <- as.character(denominator_cells)

    meta_denom <- meta_use %>%
      dplyr::filter(cell_type_internal %in% denominator_cells)

  } else {
    meta_denom <- meta_use
  }

  if (nrow(meta_denom) == 0) {
    stop("denominator_cells 过滤后没有细胞。")
  }

  # ============================================================
  # 3. 确定 other cells
  # ============================================================
  all_detected_cells <- sort(unique(meta_denom$cell_type_internal))

  missing_target <- setdiff(target_cell, all_detected_cells)

  if (length(missing_target) > 0) {
    warning(
      "以下 target_cell 在当前筛选数据中没有出现，将以 0 计数：",
      paste(missing_target, collapse = ", ")
    )
  }

  if (is.null(other_cells)) {
    other_cells <- all_detected_cells
  } else {
    other_cells <- as.character(other_cells)
  }

  if (isTRUE(exclude_target_from_other)) {
    other_cells <- setdiff(other_cells, target_cell)
  }

  other_cells <- intersect(other_cells, all_detected_cells)

  if (length(other_cells) == 0) {
    stop("没有可用于相关性分析的 other_cells。")
  }

  message(">>> Other cells 数量: ", length(other_cells))

  cells_to_calculate <- unique(c(target_cell, other_cells))

  # ============================================================
  # 4. 每个 patient 的总细胞数
  # ============================================================
  patient_total <- meta_denom %>%
    dplyr::group_by(patient_id_internal) %>%
    dplyr::summarise(
      Total_Cells = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      Total_Cells >= min_total_cells_per_patient
    )

  if (nrow(patient_total) < min_patients) {
    stop("满足 min_total_cells_per_patient 的患者数不足。")
  }

  # ============================================================
  # 5. 每个 patient × cell type 的丰度
  # ============================================================
  message(">>> 计算 patient-level cell abundance...")

  cell_count <- meta_denom %>%
    dplyr::filter(cell_type_internal %in% cells_to_calculate) %>%
    dplyr::count(
      patient_id_internal,
      cell_type_internal,
      name = "Target_Cells"
    )

  all_patient_cell_df <- tidyr::expand_grid(
    patient_id_internal = patient_total$patient_id_internal,
    cell_type_internal = cells_to_calculate
  )

  abundance_df <- all_patient_cell_df %>%
    dplyr::left_join(
      patient_total,
      by = "patient_id_internal"
    ) %>%
    dplyr::left_join(
      cell_count,
      by = c("patient_id_internal", "cell_type_internal")
    ) %>%
    dplyr::mutate(
      Target_Cells = dplyr::coalesce(Target_Cells, 0L),
      Fraction = Target_Cells / Total_Cells,
      Percentage = Fraction * 100,
      Cell_Value = dplyr::case_when(
        abundance_mode == "percentage" ~ Percentage,
        abundance_mode == "fraction" ~ Fraction,
        abundance_mode == "count" ~ as.numeric(Target_Cells),
        TRUE ~ Percentage
      )
    )

  abundance_wide <- abundance_df %>%
    dplyr::select(
      patient_id_internal,
      cell_type_internal,
      Cell_Value
    ) %>%
    tidyr::pivot_wider(
      names_from = cell_type_internal,
      values_from = Cell_Value
    )

  # ============================================================
  # 6. 相关性函数
  # ============================================================
  safe_cor_test <- function(x, y, method = "spearman") {

    ok <- !is.na(x) &
      !is.na(y) &
      is.finite(x) &
      is.finite(y)

    x <- x[ok]
    y <- y[ok]

    n <- length(x)

    if (
      n < min_patients ||
        stats::sd(x, na.rm = TRUE) == 0 ||
        stats::sd(y, na.rm = TRUE) == 0
    ) {
      return(data.frame(
        Cor = NA_real_,
        Pval = NA_real_,
        N = n
      ))
    }

    tryCatch(
      {
        ct <- suppressWarnings(
          stats::cor.test(
            x,
            y,
            method = method,
            exact = FALSE
          )
        )

        data.frame(
          Cor = unname(ct$estimate),
          Pval = ct$p.value,
          N = n
        )
      },
      error = function(e) {
        data.frame(
          Cor = NA_real_,
          Pval = NA_real_,
          N = n
        )
      }
    )
  }

  # ============================================================
  # 7. 计算 target cell vs other cells 相关性
  # ============================================================
  message(">>> 计算 target cell 与 other cell 的相关性...")

  cor_list <- list()

  pair_df <- tidyr::expand_grid(
    target_cell_name = target_cell,
    other_cell_name = other_cells
  )

  for (i in seq_len(nrow(pair_df))) {

    tc <- pair_df$target_cell_name[i]
    oc <- pair_df$other_cell_name[i]

    if (!(tc %in% colnames(abundance_wide))) {
      tmp <- data.frame(
        target_cell_name = tc,
        other_cell_name = oc,
        Cor = NA_real_,
        Pval = NA_real_,
        N = 0
      )
    } else if (!(oc %in% colnames(abundance_wide))) {
      tmp <- data.frame(
        target_cell_name = tc,
        other_cell_name = oc,
        Cor = NA_real_,
        Pval = NA_real_,
        N = 0
      )
    } else {
      res <- safe_cor_test(
        x = abundance_wide[[tc]],
        y = abundance_wide[[oc]],
        method = cor_method
      )

      tmp <- data.frame(
        target_cell_name = tc,
        other_cell_name = oc,
        Cor = res$Cor,
        Pval = res$Pval,
        N = res$N
      )
    }

    cor_list[[i]] <- tmp
  }

  cor_df <- dplyr::bind_rows(cor_list) %>%
    dplyr::mutate(
      FDR = stats::p.adjust(Pval, method = p_adjust_method)
    )

  # ============================================================
  # 8. 显著性标记
  # ============================================================
  if (sig_by == "FDR") {

    if (is.null(fdr_cutoff)) {
      fdr_cutoff <- p_cutoff
    }

    cor_df <- cor_df %>%
      dplyr::mutate(
        sig_value = FDR,
        sig_flag = !is.na(FDR) & FDR < fdr_cutoff
      )

  } else {

    cor_df <- cor_df %>%
      dplyr::mutate(
        sig_value = Pval,
        sig_flag = !is.na(Pval) & Pval < p_cutoff
      )
  }

  cor_df <- cor_df %>%
    dplyr::mutate(
      sig_label = dplyr::case_when(
        is.na(sig_value) ~ "",
        sig_value < 0.001 ~ "***",
        sig_value < 0.01 ~ "**",
        sig_value < 0.05 ~ "*",
        sig_value < 0.1 ~ "†",
        TRUE ~ ""
      ),
      cor_label = ifelse(
        is.na(Cor),
        "",
        sprintf(
          paste0("%.", cor_digits, "f"),
          Cor
        )
      ),
      tile_label = dplyr::case_when(
        show_cor_value & sig_label != "" ~ paste0(cor_label, "\n", sig_label),
        show_cor_value & sig_label == "" ~ cor_label,
        !show_cor_value ~ sig_label,
        TRUE ~ sig_label
      )
    )

  # ============================================================
  # 9. 排序
  # ============================================================
  target_cell_order <- target_cell

other_order_df <- cor_df %>%
  dplyr::group_by(other_cell_name) %>%
  dplyr::summarise(
    mean_cor = mean(Cor, na.rm = TRUE),
    mean_abs_cor = mean(abs(Cor), na.rm = TRUE),
    min_p = min(Pval, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    mean_cor = ifelse(is.nan(mean_cor), NA_real_, mean_cor),
    mean_abs_cor = ifelse(is.nan(mean_abs_cor), NA_real_, mean_abs_cor),
    min_p = ifelse(is.infinite(min_p), NA_real_, min_p)
  ) %>%
  dplyr::arrange(
    is.na(mean_cor),
    mean_cor,
    other_cell_name
  )

other_cell_order <- other_order_df$other_cell_name

  if (isTRUE(cluster_cols) && length(other_cells) >= 2) {

    mat_for_cluster <- cor_df %>%
      dplyr::select(
        target_cell_name,
        other_cell_name,
        Cor
      ) %>%
      tidyr::pivot_wider(
        names_from = other_cell_name,
        values_from = Cor
      ) %>%
      tibble::column_to_rownames("target_cell_name") %>%
      as.matrix()

    mat_for_cluster[is.na(mat_for_cluster)] <- 0

    if (ncol(mat_for_cluster) >= 2) {
      hc_cols <- stats::hclust(
        stats::dist(t(mat_for_cluster))
      )
      other_cell_order <- colnames(mat_for_cluster)[hc_cols$order]
    }
  }

  if (isTRUE(cluster_rows) && length(target_cell) >= 2) {

    mat_for_cluster <- cor_df %>%
      dplyr::select(
        target_cell_name,
        other_cell_name,
        Cor
      ) %>%
      tidyr::pivot_wider(
        names_from = other_cell_name,
        values_from = Cor
      ) %>%
      tibble::column_to_rownames("target_cell_name") %>%
      as.matrix()

    mat_for_cluster[is.na(mat_for_cluster)] <- 0

    if (nrow(mat_for_cluster) >= 2) {
      hc_rows <- stats::hclust(
        stats::dist(mat_for_cluster)
      )
      target_cell_order <- rownames(mat_for_cluster)[hc_rows$order]
    }
  }

  cor_df <- cor_df %>%
    dplyr::mutate(
      target_cell_name = factor(
        target_cell_name,
        levels = rev(target_cell_order)
      ),
      other_cell_name = factor(
        other_cell_name,
        levels = other_cell_order
      )
    )

  # ============================================================
  # 10. 画 heatmap
  # ============================================================
  if (is.null(title)) {
    title <- paste0(
      "Correlation between target cell types and other cell types"
    )
  }

  y_label_text <- dplyr::case_when(
    abundance_mode == "percentage" ~ "cell proportion (%)",
    abundance_mode == "fraction" ~ "cell fraction",
    abundance_mode == "count" ~ "cell count",
    TRUE ~ "cell abundance"
  )

  subtitle_text <- paste0(
    "Stage: ",
    ifelse(
      is.null(target_stage),
      "all",
      paste(target_stage, collapse = " + ")
    ),
    "; correlation: ",
    cor_method,
    "; value: patient-level ",
    y_label_text,
    "; significance by ",
    sig_by,
    ifelse(
      sig_by == "FDR",
      paste0(" < ", ifelse(is.null(fdr_cutoff), p_cutoff, fdr_cutoff)),
      paste0(" < ", p_cutoff)
    )
  )

  p <- ggplot(
    cor_df,
    aes(
      x = other_cell_name,
      y = target_cell_name,
      fill = Cor
    )
  ) +
    geom_tile(
      color = tile_color,
      linewidth = tile_linewidth
    ) +
    geom_text(
      aes(label = tile_label),
      size = ifelse(show_cor_value, cor_text_size, sig_text_size),
      fontface = "bold",
      color = "black",
      lineheight = 0.85
    ) +
    scale_fill_gradient2(
      low = low_color,
      mid = mid_color,
      high = high_color,
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "grey90",
      name = paste0("Correlation r\n", cor_method)
    ) +
    labs(
      title = title,
      subtitle = subtitle_text,
      x = "Other cell types",
      y = "Target cell types"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 15
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = 10
      ),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = axis_text_x_size,
        color = "black"
      ),
      axis.text.y = element_text(
        size = axis_text_y_size,
        face = "bold",
        color = "black"
      ),
      axis.title = element_text(
        size = 13,
        face = "bold"
      ),
      legend.title = element_text(
        size = 11,
        face = "bold"
      ),
      legend.text = element_text(
        size = 10
      ),
      panel.grid = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      )
    )

  # ============================================================
  # 11. 输出
  # ============================================================
  if (!is.null(output_prefix)) {

    utils::write.csv(
      cor_df,
      paste0(output_prefix, "_target_cell_vs_other_cell_correlation.csv"),
      row.names = FALSE
    )

    utils::write.csv(
      abundance_df,
      paste0(output_prefix, "_cell_abundance_by_patient.csv"),
      row.names = FALSE
    )
  }

  if (isTRUE(output_pdf)) {

    if (is.null(output_prefix)) {
      output_prefix <- "TargetCellCorrelationHeatmap"
    }

    grDevices::pdf(
      paste0(output_prefix, "_heatmap.pdf"),
      width = pdf_width,
      height = pdf_height
    )
    print(p)
    grDevices::dev.off()
  }

  message(">>> 分析完成！")

  return(list(
    plot = p,
    data = cor_df,
    abundance_data = abundance_df,
    abundance_wide = abundance_wide,
    target_cell = target_cell,
    other_cells = other_cells,
    target_stage = target_stage,
    abundance_mode = abundance_mode,
    cor_method = cor_method,
    sig_by = sig_by
  ))
}

# ===== batch6.ipynb cell 823 (code) =====
res <- PlotTargetCellCorrelationHeatmap(
  seurat_obj = seurat_merge,
  target_cell = c("CD4_Treg", "CD4_LAG3"),
  target_stage = "D7",
  show_cor_value = F
)
options(repr.plot.width = 15, repr.plot.height = 4)
print(res$plot)

# ===== batch6.ipynb cell 825 (code) =====
PlotTwoCellTypeProportionCorrelation <- function(
  seurat_obj,
  target_cell,
  target_stage = NULL,

  cell_col = "cell_type3",
  time_col = "sample_type_rn",
  patient_col = "patient_id",

  # 点的颜色分组列，可以是 meta.data 里的临床分组、治疗组、response 等
  color_col = NULL,
  color_values = NULL,

  # NULL 表示以当前 stage 的所有细胞为分母
  # 如果想在某一类细胞中计算比例，可以指定 denominator_cells
  denominator_cells = NULL,

  min_total_cells_per_patient = 10,
  min_patients = 5,

  cor_method = c("spearman", "pearson"),

  point_size = 5,
  point_alpha = 0.8,
  point_color = "midnightblue",

  smooth_method = "lm",
  smooth_color = "grey30",
  smooth_fill = "grey80",
  smooth_alpha = 0.2,
  smooth_linetype = "dashed",

  stat_cor_size = 6,

  title = NULL,
  subtitle = NULL,

  xlab = NULL,
  ylab = NULL,

  output_prefix = NULL,
  output_pdf = FALSE,
  pdf_width = 6.5,
  pdf_height = 5.8
) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(ggpubr)
  })

  cor_method <- match.arg(cor_method)

  target_cell <- as.character(target_cell)

  if (length(target_cell) != 2) {
    stop("target_cell 必须指定且只能指定两个细胞类型，例如 c('CD4_Treg', 'CD4_LAG3')。")
  }

  cell_x <- target_cell[1]
  cell_y <- target_cell[2]

  # ============================================================
  # 1. 提取 meta.data
  # ============================================================
  meta <- seurat_obj@meta.data

  required_cols <- c(
    cell_col,
    patient_col
  )

  if (!is.null(target_stage)) {
    required_cols <- c(required_cols, time_col)
  }

  if (!is.null(color_col)) {
    required_cols <- c(required_cols, color_col)
  }

  missing_cols <- setdiff(required_cols, colnames(meta))

  if (length(missing_cols) > 0) {
    stop(
      "seurat_obj@meta.data 中缺少以下列：",
      paste(missing_cols, collapse = ", ")
    )
  }

  message(">>> 提取 Seurat meta.data...")
  message(">>> Cell X: ", cell_x)
  message(">>> Cell Y: ", cell_y)

  if (!is.null(target_stage)) {
    message(">>> Target stage: ", paste(target_stage, collapse = " + "))
  }

  meta_use <- meta %>%
    tibble::rownames_to_column("cell_id") %>%
    dplyr::mutate(
      patient_id_internal = as.character(.data[[patient_col]]),
      cell_type_internal = as.character(.data[[cell_col]])
    ) %>%
    dplyr::filter(
      !is.na(patient_id_internal),
      !is.na(cell_type_internal)
    )

  if (!is.null(target_stage)) {
    meta_use <- meta_use %>%
      dplyr::mutate(
        stage_internal = as.character(.data[[time_col]])
      ) %>%
      dplyr::filter(
        !is.na(stage_internal),
        stage_internal %in% as.character(target_stage)
      )
  }

  if (nrow(meta_use) == 0) {
    stop("过滤 target_stage 后没有可用细胞。")
  }

  # ============================================================
  # 2. 分母细胞
  # ============================================================
  if (!is.null(denominator_cells)) {
    denominator_cells <- as.character(denominator_cells)

    meta_denom <- meta_use %>%
      dplyr::filter(cell_type_internal %in% denominator_cells)

  } else {
    meta_denom <- meta_use
  }

  if (nrow(meta_denom) == 0) {
    stop("denominator_cells 过滤后没有细胞。")
  }

  missing_target_cells <- setdiff(
    target_cell,
    unique(meta_denom$cell_type_internal)
  )

  if (length(missing_target_cells) > 0) {
    warning(
      "以下 target_cell 在当前筛选数据中没有出现，将以 0 计数：",
      paste(missing_target_cells, collapse = ", ")
    )
  }

  # ============================================================
  # 3. 每个 patient 的总细胞数
  # ============================================================
  patient_total <- meta_denom %>%
    dplyr::group_by(patient_id_internal) %>%
    dplyr::summarise(
      Total_Cells = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(
      Total_Cells >= min_total_cells_per_patient
    )

  if (nrow(patient_total) < min_patients) {
    stop("满足 min_total_cells_per_patient 的患者数不足。")
  }

  # ============================================================
  # 4. 每个 patient 中两个细胞类型的比例
  # ============================================================
  message(">>> 计算 patient-level cell proportion...")

  cell_count <- meta_denom %>%
    dplyr::filter(cell_type_internal %in% target_cell) %>%
    dplyr::count(
      patient_id_internal,
      cell_type_internal,
      name = "Target_Cells"
    )

  all_patient_cell_df <- tidyr::expand_grid(
    patient_id_internal = patient_total$patient_id_internal,
    cell_type_internal = target_cell
  )

  prop_long <- all_patient_cell_df %>%
    dplyr::left_join(
      patient_total,
      by = "patient_id_internal"
    ) %>%
    dplyr::left_join(
      cell_count,
      by = c("patient_id_internal", "cell_type_internal")
    ) %>%
    dplyr::mutate(
      Target_Cells = dplyr::coalesce(Target_Cells, 0L),
      Proportion = Target_Cells / Total_Cells * 100
    )

  prop_wide <- prop_long %>%
    dplyr::select(
      patient_id_internal,
      cell_type_internal,
      Proportion
    ) %>%
    tidyr::pivot_wider(
      names_from = cell_type_internal,
      values_from = Proportion
    )

  # ============================================================
  # 5. 提取 color_col
  # ============================================================
  if (!is.null(color_col)) {

    color_df <- meta_use %>%
      dplyr::select(
        patient_id_internal,
        color_value = dplyr::all_of(color_col)
      ) %>%
      dplyr::mutate(
        color_value = as.character(color_value)
      ) %>%
      dplyr::filter(!is.na(color_value)) %>%
      dplyr::group_by(patient_id_internal) %>%
      dplyr::summarise(
        color_value = dplyr::first(color_value),
        n_color_values = dplyr::n_distinct(color_value),
        .groups = "drop"
      )

    if (any(color_df$n_color_values > 1, na.rm = TRUE)) {
      warning(
        "部分 patient 存在多个 color_col 取值；当前使用 first(color_col)。"
      )
    }

    merged_data <- prop_wide %>%
      dplyr::left_join(
        color_df %>%
          dplyr::select(patient_id_internal, color_value),
        by = "patient_id_internal"
      )

  } else {

    merged_data <- prop_wide
  }

  merged_data <- merged_data %>%
    dplyr::filter(
      !is.na(.data[[cell_x]]),
      !is.na(.data[[cell_y]]),
      is.finite(.data[[cell_x]]),
      is.finite(.data[[cell_y]])
    )

  if (nrow(merged_data) < min_patients) {
    stop("合并后可用于相关性分析的 patient 数不足。")
  }

  # ============================================================
  # 6. 相关性计算
  # ============================================================
  cor_test <- tryCatch(
    {
      stats::cor.test(
        merged_data[[cell_x]],
        merged_data[[cell_y]],
        method = cor_method,
        exact = FALSE
      )
    },
    error = function(e) NULL
  )

  if (!is.null(cor_test)) {
    cor_r <- unname(cor_test$estimate)
    cor_p <- cor_test$p.value
  } else {
    cor_r <- NA_real_
    cor_p <- NA_real_
  }

  message(">>> Patients included: ", nrow(merged_data))
  message(">>> Correlation r: ", signif(cor_r, 4))
  message(">>> P value: ", signif(cor_p, 4))

  # ============================================================
  # 7. 作图
  # ============================================================
  if (is.null(title)) {
    title <- paste0(
      "Global Correlation",
      ifelse(
        is.null(target_stage),
        "",
        paste0(" in ", paste(target_stage, collapse = " + "))
      )
    )
  }

  if (is.null(subtitle)) {
    subtitle <- paste0(
      "X: % of ",
      cell_x,
      " | Y: % of ",
      cell_y,
      " | N = ",
      nrow(merged_data)
    )
  }

  if (is.null(xlab)) {
    xlab <- paste0("% of ", cell_x, " in total cells")
  }

  if (is.null(ylab)) {
    ylab <- paste0("% of ", cell_y, " in total cells")
  }

  p <- ggplot(
    merged_data,
    aes(
      x = .data[[cell_x]],
      y = .data[[cell_y]]
    )
  ) +
    geom_smooth(
      method = smooth_method,
      color = smooth_color,
      fill = smooth_fill,
      alpha = smooth_alpha,
      linetype = smooth_linetype
    ) +
    {
      if (!is.null(color_col)) {
        geom_point(
          aes(color = color_value),
          size = point_size,
          alpha = point_alpha
        )
      } else {
        geom_point(
          size = point_size,
          alpha = point_alpha,
          color = point_color
        )
      }
    } +
    ggpubr::theme_pubr() +
    ggpubr::stat_cor(
      method = cor_method,
      size = stat_cor_size
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab,
      color = color_col
    )

  if (!is.null(color_col) && !is.null(color_values)) {
    p <- p +
      scale_color_manual(
        values = cell_type_colors_updated
      )
  }

  p <- p +
    theme(
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 14, face = "bold"),

      axis.text = element_text(
        size = 14,
        color = "black"
      ),
      axis.text.x = element_text(
        size = 14,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 12,
        color = "black"
      ),
      axis.title = element_text(
        size = 16,
        face = "bold"
      ),

      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 15
      ),
      plot.subtitle = element_text(
        hjust = 0.5
      ),

      legend.position = ifelse(
        is.null(color_col),
        "none",
        "right"
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = 14,
        color = "black"
      ),

      panel.grid.major = element_line(
        color = "grey90"
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      )
    )

  # ============================================================
  # 8. 输出
  # ============================================================
  if (!is.null(output_prefix)) {
    utils::write.csv(
      merged_data,
      paste0(output_prefix, "_two_celltype_proportion_correlation_data.csv"),
      row.names = FALSE
    )

    utils::write.csv(
      prop_long,
      paste0(output_prefix, "_two_celltype_proportion_long.csv"),
      row.names = FALSE
    )
  }

  if (isTRUE(output_pdf)) {

    if (is.null(output_prefix)) {
      output_prefix <- paste0(
        "TwoCellTypeCorrelation_",
        cell_x,
        "_vs_",
        cell_y
      )
    }

    grDevices::pdf(
      paste0(output_prefix, "_scatter.pdf"),
      width = pdf_width,
      height = pdf_height
    )
    print(p)
    grDevices::dev.off()
  }

  message(">>> 分析完成！")

  return(list(
    plot = p,
    data = merged_data,
    proportion_long = prop_long,
    cor_test = cor_test,
    cor_r = cor_r,
    cor_p = cor_p,
    target_cell = target_cell,
    target_stage = target_stage,
    color_col = color_col
  ))
}

# ===== batch6.ipynb cell 826 (code) =====
res <- PlotTwoCellTypeProportionCorrelation(
  seurat_obj = seurat_merge,
  target_cell = c("CD4_LAG3", "CD8_NKT"),
  target_stage = "D7",
  color_col = "lm", color_values = cell_type_colors_updated
)
options(repr.plot.width = 6.5, repr.plot.height = 5.5)
print(res$plot)

# ===== batch6.ipynb cell 827 (code) =====
ggsave("plot_v4/D7 correlation CD4_LAG3 vs CD8_NKT.pdf", width = 6.5, height = 5.5)

# ===== batch6.ipynb cell 839 (code) =====
#' 计算 CD4 T 细胞功能评分
#' 
#' @param cd4_subset CD4 亚群的 Seurat 对象
#' @param ctrl 对照基因数量，默认为 100
#' @return 返回增加了 7 列评分（以 _score 结尾）的 Seurat 对象
score_cd4_subsets <- function(cd4_subset, ctrl = 100) {
  
  # 1. 定义基因集
  gene_sets <- list(
    Exhaustion     = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4", "TOX", "TOX2", 
                       "NR4A1", "NR4A2", "NR4A3", "BATF", "PRDM1", "LAYN", "CXCL13"),
    Helper         = c("IL2", "IFNG", "TNF", "CD40LG", "IL21", "CSF2"),
    Treg           = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "ENTPD1", "TNFRSF18", "CCR8", "LRRC32"),
    Tr1_like       = c("IL10", "LAG3", "ITGA2", "MAF", "PRDM1", "ICOS", "CTLA4"),
    Cytotoxic      = c("GZMB", "GZMK", "PRF1", "NKG7", "GNLY", "CST7", "IFNG", "CCL5", "EOMES", "RUNX3"),
    Acute_Activate = c("CD69", "IL2RA", "CD40LG", "ICOS", "TNFRSF4", "TNFRSF9", "JUN", "FOS", "NFKBIA", "DUSP1"),
    TypeI_IFN      = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2", "OAS1", "OAS2", "IFI6", "IFI44L", "RSAD2")
  )
  
  message("开始计算 CD4 功能评分...")
  
  # 2. 循环打分
  for(name in names(gene_sets)) {
    # 提取在当前对象中存在的基因
    valid_genes <- intersect(gene_sets[[name]], rownames(cd4_subset))
    
    if(length(valid_genes) < 3) {
      warning(paste("基因集", name, "在对象中找到的基因过少，可能影响准确性。"))
    }
    
    # 执行 Seurat 打分
    cd4_subset <- AddModuleScore(
      object = cd4_subset,
      features = list(valid_genes),
      name = name,
      ctrl = ctrl
    )
    
    # 清理列名：AddModuleScore 会生成 "Name1" 形式的列名
    # 我们将其重命名为 "Name_score" 并移除中间生成的列
    old_col <- paste0(name, "1")
    new_col <- paste0(name, "_score")
    
    cd4_subset@meta.data[[new_col]] <- cd4_subset@meta.data[[old_col]]
    cd4_subset@meta.data[[old_col]] <- NULL
    
    message(paste("✓ 完成评分:", name, "(匹配基因:", length(valid_genes), "个)"))
  }
  
  message("所有评分计算完成！")
  return(cd4_subset)
}

# ===== batch6.ipynb cell 843 (code) =====
cd8_obj <- score_cd4_subsets(cd8_obj)

# ===== batch6.ipynb cell 844 (code) =====
library(pheatmap)
library(dplyr)
library(tidyr)

# 1. 定义我们之前计算的 7 个得分列名
score_cols <- c("Exhaustion_score", "Helper_score", "Treg_score", 
                "Tr1_like_score", "Cytotoxic_score", "Acute_Activate_score", "TypeI_IFN_score")

# 2. 提取元数据并计算每个 cell_type3 的平均得分
plot_matrix <- cd8_obj@meta.data %>%
  group_by(cell_type3) %>%
  summarise(across(all_of(score_cols), \(x) mean(x, na.rm = TRUE))) %>%
  # 将 cell_type3 设为行名
  tibble::column_to_rownames("cell_type3")

# 3. 为了美观，去掉列名中的 "_score" 后缀
colnames(plot_matrix) <- gsub("_score", "", colnames(plot_matrix))

# ===== batch6.ipynb cell 845 (code) =====
options(repr.plot.width = 6, repr.plot.height = 5, warnings = F)
# 绘图设置
pheatmap(
  t(plot_matrix),            # 转置矩阵，让亚群作为列，功能模块作为行
  scale = "row",             # 按行（即功能模块）进行标准化，对比各亚群差异
  clustering_method = "ward.D2", 
  cluster_rows = TRUE,       # 对功能模块进行聚类
  cluster_cols = TRUE,       # 对细胞亚群进行聚类
  
  # 配色：使用经典的蓝-白-红配色 (Firebrick)
  color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  
  # 界面优化
  border_color = "white",    # 格子间隙为白色，更清爽
  main = "Functional Characteristics of CD4 Subsets",
  fontsize_row = 15,
  fontsize_col = 15,
  angle_col = 45,            # 亚群名称倾斜
  
  # 图例说明：展示的是 Z-score
  legend_labels = "Z-score"
)

# ===== batch6.ipynb cell 848 (code) =====
#' 计算 CD8 T 细胞功能评分
#' 
#' @param cd8_subset CD8 亚群的 Seurat 对象
#' @param ctrl 对照基因数量，默认为 100
#' @return 返回增加了 5 列评分（以 _score 结尾）的 Seurat 对象
score_cd8_subsets <- function(cd8_subset, ctrl = 100) {
  
  # 1. 定义 CD8 特异性基因集
  gene_sets <- list(
    Cytox      = c("PRF1", "IFNG", "GNLY", "NKG7", "GZMB", "GZMA", 
                   "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7"),
    
    Exhaustion = c("LAG3", "TIGIT", "PDCD1", "CTLA4", "HAVCR2"),
    
    Activation = c("CD69", "CCR7", "CD27", "BTLA", "CD40LG", "IL2RA", "CD3E", 
                   "CD47", "EOMES", "GNLY", "GZMA", "GZMB", "PRF1", "IFNG", 
                   "CD8A", "CD8B", "FASLG", "LAMP1", "LAG3", "CTLA4", 
                   "HLA-DRA", "TNFRSF4", "ICOS", "TNFRSF9", "TNFRSF18"),
                   
    Naive      = c("CCR7", "SELL", "TCF7", "LEF1"),
    
    Prolif     = c("MKI67", "TOP2A", "PCNA", "BIRC5", "MCM2")
  )
  
  message("开始计算 CD8 功能评分...")
  
  # 2. 循环执行 AddModuleScore
  for(name in names(gene_sets)) {
    # 提取在当前对象中存在的基因
    valid_genes <- intersect(gene_sets[[name]], rownames(cd8_subset))
    
    if(length(valid_genes) < 2) {
      warning(paste("基因集", name, "在对象中找到的基因过少（少于2个），跳过该评分。"))
      next
    }
    
    # 执行 Seurat 打分
    # AddModuleScore 默认会将结果存为 "name1"
    cd8_subset <- AddModuleScore(
      object = cd8_subset,
      features = list(valid_genes),
      name = name,
      ctrl = ctrl
    )
    
    # 3. 整理列名：将 "Name1" 转换为 "Name_score"
    old_col <- paste0(name, "1")
    new_col <- paste0(name, "_score")
    
    if(old_col %in% colnames(cd8_subset@meta.data)) {
      cd8_subset@meta.data[[new_col]] <- cd8_subset@meta.data[[old_col]]
      cd8_subset@meta.data[[old_col]] <- NULL # 移除原始生成的列
      message(paste("✓ 完成评分:", name, "(匹配基因:", length(valid_genes), "个)"))
    }
  }
  
  message("所有 CD8 评分计算完成！")
  return(cd8_subset)
}

# ===== batch6.ipynb cell 849 (code) =====
cd8_obj <- score_cd8_subsets(cd8_obj)

# ===== batch6.ipynb cell 850 (code) =====
library(pheatmap)
library(dplyr)
library(tidyr)

# 1. 定义 CD8 特有的 5 个得分列名
score_cols <- c("Cytox_score", "Exhaustion_score", "Activation_score", 
                "Naive_score", "Prolif_score")

# 2. 提取元数据并计算每个细胞亚群 (cell_type3) 的平均得分
# 确保此时的 cd8_obj 已经运行过 score_cd8_subsets 函数
plot_matrix <- cd8_obj@meta.data %>%
  group_by(cell_type3) %>%
  summarise(across(all_of(score_cols), \(x) mean(x, na.rm = TRUE))) %>%
  # 将 cell_type3 设为行名
  tibble::column_to_rownames("cell_type3")

# 3. 去掉列名中的 "_score" 后缀，使图形整洁
colnames(plot_matrix) <- gsub("_score", "", colnames(plot_matrix))

# 4. 绘图设置
options(repr.plot.width = 5.5, repr.plot.height = 5, warnings = FALSE)

pheatmap(
  t(plot_matrix),            # 转置矩阵：行是功能模块，列是细胞亚群
  scale = "row",             # 按行标准化 (Z-score)，突出各亚群在同一功能上的高低差异
  clustering_method = "ward.D2", 
  cluster_rows = TRUE,       # 对功能模块聚类
  cluster_cols = TRUE,       # 对细胞亚群聚类
  
  # 配色：红蓝经典配色
  color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
  
  # 界面美化
  border_color = "white",
  main = "Functional Characteristics of CD8 Subsets",
  fontsize_row = 14,
  fontsize_col = 14,
  angle_col = 45,            # X轴标签倾斜
  
  # 展示 Z-score 范围
  legend_labels = "Z-score"
)

# ===== batch6.ipynb cell 851 (code) =====


# ===== batch6.ipynb cell 852 (code) =====
ggsave("CD8_functional characteristics of CD8 subsets.pdf", width = 5.5, height = 5)

# ===== batch6.ipynb cell 889 (code) =====
plot_clusters_group <- function(meta_data, 
                                group_var = "new_group", 
                                color_var = "therapy_effect_adj",  
                                cell_cluster = "seurat_clusters",
                                show_clusters = NULL, # 新增参数：指定显示的亚群，例如 c("Cluster1", "Cluster2")
                                exclude_cell_type = NULL, # 新增参数：排除的细胞类型
                                levels = c("PR_BL", "PR_CT", "PR_DR", "SD_BL", "SD_CT", "SD_DR", "PD_BL", "PD_CT"),
                                comparision = list(c("PR_BL", "SD_BL"), c("PR_BL", "PR_DR")), 
                                methods = "part", 
                                line = FALSE) {
  
  require(ggplot2)
  require(dplyr)
  require(ggpubr)
  
  # 1. 计算频率和百分比
  if(methods == "part"){
    # 注意：计算百分比的基数（分母）通常应该基于过滤掉 exclude_cell_type 后的总数
    results <- meta_data %>% 
      filter(!(cell_type3 %in% exclude_cell_type), !is.na(cell_type3)) %>%  
      group_by(orig.ident) %>%  
      mutate(sample_cells = n()) %>%  
      ungroup() %>%  
      group_by(orig.ident, !!sym(cell_cluster)) %>%
      mutate(cluster_freq = n()) %>%  
      ungroup() %>%  
      dplyr::select(orig.ident, patient_id, !!sym(cell_cluster), sample_cells, cluster_freq, !!sym(color_var), !!sym(group_var)) %>%  
      distinct() %>%  
      mutate(percent = cluster_freq / sample_cells * 100)
  } else {
    results <- meta_data %>%  
      group_by(orig.ident, !!sym(cell_cluster)) %>%
      mutate(cluster_freq = n()) %>%  
      ungroup() %>%  
      # 注意：此处 n_sample_cells 需确保存在于 meta_data 中，否则需根据 orig.ident 计算
      dplyr::select(orig.ident, patient_id, !!sym(cell_cluster), cluster_freq, !!sym(color_var), !!sym(group_var)) %>%  
      group_by(orig.ident) %>%
      mutate(total_cells = sum(cluster_freq)) %>% # 如果没有外部预设总数，动态计算总数
      ungroup() %>%
      distinct() %>%  
      mutate(percent = cluster_freq / total_cells * 100)
  }

  # --- 修改部分：筛选指定的亚群 ---
  if(!is.null(show_clusters)) {
    results <- results %>% filter(!!sym(cell_cluster) %in% show_clusters)
    # 转换为 factor 以保证按照指定顺序显示
    results[[cell_cluster]] <- factor(results[[cell_cluster]], levels = show_clusters)
  }
  # ----------------------------

  # 2. 调整 group_var 的因子水平
  results[[group_var]] <- factor(results[[group_var]], levels = levels)

  # 3. 绘图准备
  cell_cluster_formula <- as.formula(paste0("~ ", cell_cluster))
  y_labs = ifelse(methods == "part", "Percentage in samples", "Percentage in samples")

  # 4. 构建基础图形
  p <- ggplot(results, aes(x = !!sym(group_var), y = percent)) +  
    geom_boxplot(outlier.shape = NA) 

  # 添加配对连线
  if(line) {
    p <- p + geom_line(aes(group = patient_id), color = "grey70", alpha = 0.8, linetype = "solid")
  }

  p <- p +  
    geom_jitter(aes(color = !!sym(color_var)), width = 0.1) + 
    facet_wrap(cell_cluster_formula, scales = "free_y") +  
    theme_classic() +
    theme(
      strip.text = element_text(size = 14),  
      axis.text = element_text(size = 14),
      axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14)
    ) + scale_color_manual(values = cell_type_colors_updated) +
    stat_compare_means(comparisons = comparision, method = "t.test", size = 5,  method.args = list(alternative = "less"), ) +  
    labs(y = y_labs)+
    theme(
      strip.background = element_rect(fill = "white"),       
      strip.text = element_text(size = 14, face = "bold"),   
      axis.text = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 14, angle = 30, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 16, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none", 
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8) # ggplot2 新版推荐用 linewidth 替代 size
    )
  
  return(p)
}

# ===== batch6.ipynb cell 893 (code) =====
options(repr.plot.width = 7.2, repr.plot.height = 5.5, warnings = F)
plot_clusters_group(seurat_merge@meta.data, 
                    "lm_group", 
                    show_clusters = c("CD4_LAG3", "CD4_Treg"), 
                    "lm", "cell_type3", line = F,
              levels =  sort(unique(seurat_merge$lm_group)),
             comparision = list(c("noLM_D0", "LM_D0"), c("noLM_D7", "LM_D7")))

# ===== batch6.ipynb cell 894 (code) =====
ggsave("plot_v4/Fig3g.pdf", width = 7.2, height = 5.5)

# ===== batch6.ipynb cell 933 (code) =====
plot_stage_cell_survival <- function(
  meta, 
  target_cell, 
  target_stage,
  cell_col = "cell_type3", 
  time_col = "sample_type_rn", 
  patient_col = "patient_id",
  surv_time_col = "os_time",
  surv_status_col = "is_death",

  start_at_zero = TRUE,
  break_time_by = NULL,
  x_axis_padding_frac = 0.05,

  handle_zero_time = c("keep", "epsilon", "drop"),
  epsilon_time = 0.001,

  # ============================================================
  # 字体和版式参数
  # ============================================================
  main_title_font_size = 18,
  subtitle_font_size = 15,
  axis_title_font_size = 16,
  axis_text_font_size = 14,
  legend_font_size = 14,

  logrank_label_font_size = 5.5,
  logrank_label_x_frac = 0.05,
  logrank_label_y = 0.15,

  risk_table_font_size = 6,
  risk_table_axis_font_size = 14,
  risk_table_title_font_size = 14,
  risk_table_height = 0.32,

  # TRUE: 保留 risk table 左侧 Low / High
  # FALSE: 更容易和主图横坐标严格对齐
  risk_table_show_y_text = TRUE
) {
  
  handle_zero_time <- match.arg(handle_zero_time)
  
  suppressPackageStartupMessages({
    library(survival)
    library(survminer)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  })

  message(sprintf(">>> 步骤 1: 提取阶段 [%s] 的数据并计算细胞比例...", target_stage))
  
  # ============================================================
  # 1. 检查必要列是否存在
  # ============================================================
  required_cols <- c(
    cell_col,
    time_col,
    patient_col,
    surv_time_col,
    surv_status_col
  )

  missing_cols <- setdiff(required_cols, colnames(meta))

  if (length(missing_cols) > 0) {
    stop(
      "meta 中缺少以下列：",
      paste(missing_cols, collapse = ", ")
    )
  }

  # ============================================================
  # 2. 过滤指定阶段，计算每个患者在该阶段的细胞比例
  # ============================================================
  prop_data <- meta %>%
    dplyr::filter(.data[[time_col]] == target_stage) %>%
    dplyr::group_by(.data[[patient_col]]) %>%
    dplyr::summarise(
      Total_Cells = dplyr::n(),
      Target_Cells = sum(.data[[cell_col]] == target_cell, na.rm = TRUE),
      Percentage = (Target_Cells / Total_Cells) * 100,
      Surv_Time = dplyr::first(.data[[surv_time_col]]),
      Raw_Status = dplyr::first(.data[[surv_status_col]]),
      n_surv_time_values = dplyr::n_distinct(.data[[surv_time_col]]),
      n_surv_status_values = dplyr::n_distinct(.data[[surv_status_col]]),
      .groups = "drop"
    ) %>%
    dplyr::filter(Total_Cells > 0)

  inconsistent_info <- prop_data %>%
    dplyr::filter(
      n_surv_time_values > 1 |
        n_surv_status_values > 1
    )

  if (nrow(inconsistent_info) > 0) {
    warning(
      "检测到部分 patient 在该阶段存在多个生存时间或状态值；",
      "当前函数使用 first()。"
    )
  }

  prop_data <- prop_data %>%
    dplyr::select(
      dplyr::all_of(patient_col),
      Total_Cells,
      Target_Cells,
      Percentage,
      Surv_Time,
      Raw_Status
    )
  
  if (nrow(prop_data) < 5) {
    stop("该阶段的有效样本量不足 N < 5，无法进行生存分析。")
  }

  message(">>> 步骤 2: 清洗生存状态并剔除缺失值...")
  
  # ============================================================
  # 3. 清洗生存状态
  # ============================================================
  surv_df <- prop_data %>%
    dplyr::mutate(
      Surv_Time = suppressWarnings(as.numeric(Surv_Time)),
      Raw_Status_chr = trimws(as.character(Raw_Status)),
      Surv_Status = dplyr::case_when(
        Raw_Status_chr %in% c(
          "Yes", "yes", "Y", "y",
          "Dead", "Death", "deceased", "Deceased",
          "1", "TRUE", "True", "true",
          "是", "死亡",
          "event", "Event",
          "PD", "progression", "Progression",
          "Progressed", "progressed",
          "Relapse", "relapse",
          "Recurrence", "recurrence"
        ) ~ 1,

        Raw_Status_chr %in% c(
          "No", "no", "N", "n",
          "Alive", "alive",
          "0", "FALSE", "False", "false",
          "否", "存活",
          "censored", "Censored",
          "未进展",
          "No progression", "no progression"
        ) ~ 0,

        TRUE ~ suppressWarnings(as.numeric(Raw_Status_chr))
      )
    ) %>%
    dplyr::filter(
      !is.na(Surv_Time),
      !is.na(Surv_Status),
      !is.na(Percentage),
      is.finite(Percentage)
    ) %>%
    dplyr::select(-Raw_Status_chr)

  invalid_status <- surv_df %>%
    dplyr::filter(!Surv_Status %in% c(0, 1))

  if (nrow(invalid_status) > 0) {
    warning("检测到 Surv_Status 中存在非 0/1 的值；这些样本将被剔除。")

    surv_df <- surv_df %>%
      dplyr::filter(Surv_Status %in% c(0, 1))
  }
  
  if (nrow(surv_df) < 5) {
    stop("具有完整生存信息的样本量不足。")
  }
  
  # ============================================================
  # 4. 处理 Surv_Time <= 0 的情况
  # ============================================================
  n_zero_time <- sum(surv_df$Surv_Time <= 0, na.rm = TRUE)
  
  if (n_zero_time > 0) {
    message(">>> 检测到 Surv_Time <= 0 的患者数: ", n_zero_time)
    
    if (handle_zero_time == "epsilon") {
      message(">>> 将 Surv_Time <= 0 替换为 epsilon_time = ", epsilon_time)
      surv_df <- surv_df %>%
        dplyr::mutate(
          Surv_Time = ifelse(Surv_Time <= 0, epsilon_time, Surv_Time)
        )
    }
    
    if (handle_zero_time == "drop") {
      message(">>> 删除 Surv_Time <= 0 的患者")
      surv_df <- surv_df %>%
        dplyr::filter(Surv_Time > 0)
    }
    
    if (handle_zero_time == "keep") {
      message(">>> 保留 Surv_Time <= 0 的患者")
    }
  }
  
  if (nrow(surv_df) < 5) {
    stop("处理 Surv_Time <= 0 后，样本量不足。")
  }

  message(">>> 步骤 3: 根据细胞比例中位数分组并拟合生存曲线...")
  
  # ============================================================
  # 5. Median cutoff 分组
  #    Percentage > median 为 High
  #    Percentage <= median 为 Low
  # ============================================================
  threshold <- median(surv_df$Percentage, na.rm = TRUE)

  surv_df <- surv_df %>%
    dplyr::mutate(
      Group = ifelse(Percentage >= threshold, "High", "Low"),
      Group = factor(Group, levels = c("Low", "High"))
    )

  if (length(unique(stats::na.omit(surv_df$Group))) < 2) {
    stop("分组失败：High/Low 只有一组。")
  }

  message(">>> 分组总人数:")
  print(table(surv_df$Group, useNA = "ifany"))
  
  message(">>> 各组事件数:")
  print(table(surv_df$Group, surv_df$Surv_Status, useNA = "ifany"))
  
  if (n_zero_time > 0) {
    message(">>> Surv_Time <= 0 的患者:")
    print(
      surv_df %>%
        dplyr::filter(Surv_Time <= 0) %>%
        dplyr::select(
          dplyr::all_of(patient_col),
          Surv_Time,
          Surv_Status,
          Percentage,
          Group
        )
    )
  }

  # ============================================================
  # 6. 拟合 KM
  # ============================================================
  fit <- survival::survfit(
    survival::Surv(Surv_Time, Surv_Status) ~ Group,
    data = surv_df
  )

  # ============================================================
  # 7. Log-rank P
  # ============================================================
  surv_diff <- survival::survdiff(
    survival::Surv(Surv_Time, Surv_Status) ~ Group,
    data = surv_df
  )

  logrank_p <- 1 - stats::pchisq(
    surv_diff$chisq,
    df = length(surv_diff$n) - 1
  )

  format_p <- function(p) {
    if (is.na(p)) {
      return("NA")
    } else if (p < 0.001) {
      return("<0.001")
    } else {
      return(sprintf("%.3f", p))
    }
  }

  stat_label <- paste0(
    "Log-rank P = ",
    format_p(logrank_p)
  )

  # ============================================================
  # 8. 设置 x 轴
  # ============================================================
  max_time <- max(surv_df$Surv_Time, na.rm = TRUE)
  min_time <- min(surv_df$Surv_Time, na.rm = TRUE)

  max_time_int <- max(1, ceiling(max_time))
  min_time_int <- floor(min_time)
  
  if (is.null(break_time_by)) {
    break_time_by <- max(1, ceiling(max_time_int / 5))
  } else {
    break_time_by <- ceiling(break_time_by)
    break_time_by <- max(1, break_time_by)
  }
  
  integer_breaks <- seq(
    0,
    max_time_int,
    by = break_time_by
  )

  x_padding <- max_time_int * x_axis_padding_frac

  if (start_at_zero) {
    xlim_use <- c(
      -x_padding,
      max_time_int + x_padding
    )
  } else {
    xlim_use <- c(
      min_time_int - x_padding,
      max_time_int + x_padding
    )
  }

  ylab_text <- ifelse(
    grepl("os", surv_time_col, ignore.case = TRUE),
    "Overall survival probability",
    "Progression-free survival probability"
  )

  # ============================================================
  # 9. 绘图
  # ============================================================
  p <- survminer::ggsurvplot(
    fit, 
    data = surv_df,

    pval = FALSE,
    conf.int = FALSE,

    risk.table = TRUE, 
    risk.table.col = "strata",
    risk.table.height = risk_table_height,

    # number at risk 字体大小
    fontsize = risk_table_font_size,

    risk.table.y.text = risk_table_show_y_text,
    risk.table.y.text.col = risk_table_show_y_text,

    palette = c("#377EB8", "#E41A1C"),

    title = paste0("Survival Analysis: ", target_cell, " at ", target_stage),
    subtitle = paste0("Cutoff (median) = ", round(threshold, 2), "%"),

    legend.title = "Proportion",
    legend.labs = c("Low", "High"),

    xlab = "Time",
    ylab = ylab_text,

    # 这里先用带 padding 的 xlim
    xlim = xlim_use,
    break.time.by = break_time_by,

    font.title = c(main_title_font_size, "bold"),
    font.subtitle = c(subtitle_font_size),
    font.x = c(axis_title_font_size, "bold"),
    font.y = c(axis_title_font_size, "bold"),
    font.tickslab = axis_text_font_size,
    font.legend = legend_font_size,

    ggtheme = ggplot2::theme_bw(base_size = axis_text_font_size) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = main_title_font_size
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5,
          size = subtitle_font_size
        ),
        axis.title.x = ggplot2::element_text(
          size = axis_title_font_size,
          face = "bold"
        ),
        axis.title.y = ggplot2::element_text(
          size = axis_title_font_size,
          face = "bold"
        ),
        axis.text.x = ggplot2::element_text(
          size = axis_text_font_size,
          color = "black"
        ),
        axis.text.y = ggplot2::element_text(
          size = axis_text_font_size,
          color = "black"
        ),
        legend.title = ggplot2::element_text(
          size = legend_font_size,
          face = "bold"
        ),
        legend.text = ggplot2::element_text(
          size = legend_font_size
        ),
        plot.margin = ggplot2::margin(
          t = 5,
          r = 10,
          b = 5,
          l = 10
        )
      ),

    # 保留 risk table 网格风格
    # 不使用 theme_cleantable()
    tables.theme = ggplot2::theme_bw(base_size = risk_table_axis_font_size)
  )

  # ============================================================
  # 10. 统一 x 轴 scale
  #     主图和 risk table 都用这个，确保一致
  # ============================================================
  common_x_scale <- function() {
    ggplot2::scale_x_continuous(
      limits = xlim_use,
      breaks = integer_breaks,
      labels = function(x) sprintf("%d", round(x)),
      expand = ggplot2::expansion(mult = c(0, 0))
    )
  }

  # ============================================================
  # 11. 主图：统一横坐标和字体
  # ============================================================
  p$plot <- p$plot +
    common_x_scale() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = main_title_font_size
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = subtitle_font_size
      ),
      axis.title.x = ggplot2::element_text(
        size = axis_title_font_size,
        face = "bold"
      ),
      axis.title.y = ggplot2::element_text(
        size = axis_title_font_size,
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        size = axis_text_font_size,
        color = "black"
      ),
      axis.text.y = ggplot2::element_text(
        size = axis_text_font_size,
        color = "black"
      ),
      legend.title = ggplot2::element_text(
        size = legend_font_size,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = legend_font_size
      ),
      plot.margin = ggplot2::margin(
        t = 5,
        r = 10,
        b = 5,
        l = 10
      )
    )

  # ============================================================
  # 12. 添加 Log-rank P 标注
  # ============================================================
  p$plot <- p$plot +
    ggplot2::annotate(
      "text",
      x = xlim_use[1] + diff(xlim_use) * logrank_label_x_frac,
      y = logrank_label_y,
      hjust = 0,
      size = logrank_label_font_size,
      label = stat_label
    )

  # ============================================================
  # 13. Risk table：统一横坐标，保留网格，只改字体
  # ============================================================
  p$table <- p$table +
    common_x_scale()

  # 关键：
  # 不再给 p$table + theme(...)
  # 而是直接修改 p$table$theme 中的文字元素
  # 这样可以保留 risk table 原有网格、边框、背景风格
  # 同时避免 element_blank / element_text merge 报错

  p$table$theme$plot.title <- ggplot2::element_text(
    size = risk_table_title_font_size,
    face = "bold",
    hjust = 0,
    color = "black"
  )

  p$table$theme$axis.title.x <- ggplot2::element_text(
    size = risk_table_axis_font_size,
    face = "bold",
    color = "black"
  )

  # y 轴标题容易撑开左侧宽度，影响和主图时间轴对齐
  p$table$theme$axis.title.y <- ggplot2::element_text(angle =90,
    size = risk_table_axis_font_size,
    face = "bold",
    color = "black"
  )


  p$table$theme$axis.text.x <- ggplot2::element_text(
    size = risk_table_axis_font_size,
    color = "black"
  )

  p$table$theme$axis.text.y <- if (risk_table_show_y_text) {
    ggplot2::element_text(
      size = risk_table_axis_font_size,
      color = "black"
    )
  } else {
    ggplot2::element_blank()
  }

  p$table$theme$strip.text <- ggplot2::element_text(
    size = risk_table_axis_font_size,
    face = "bold",
    color = "black"
  )

  p$table$theme$legend.position <- "none"

  p$table$theme$plot.margin <- ggplot2::margin(
    t = 0,
    r = 10,
    b = 5,
    l = 10
  )

  message(">>> 分析完成！")
  
  return(list(
    plot = p, 
    data = surv_df, 
    median_cutoff = threshold,
    logrank_p = logrank_p,
    logrank_label = stat_label,
    break_time_by = break_time_by,
    xlim = xlim_use,
    x_breaks = integer_breaks,
    fit = fit,
    survdiff = surv_diff
  ))
}

# ===== batch6.ipynb cell 935 (code) =====
# 假设你想分析 "D0" 阶段，"CD8+ T" 细胞的比例对预后的影响
res <- plot_stage_cell_survival(
  meta = seurat_merge@meta.data,
  target_cell = "Neu_IL1R1",
  target_stage = "D7",
  cell_col = "cell_type3",
  time_col = "sample_type_rn",
  patient_col = "orig.ident",
  surv_time_col = "pfs_time",
  surv_status_col = "pfs_status",
    start_at_zero = T,
    risk_table_height = 0.25,
)


# ===== batch6.ipynb cell 937 (code) =====
options(repr.plot.width = 6, repr.plot.height = 7)
# 显示图像
print(res$plot)

# ===== batch6.ipynb cell 938 (code) =====
write.csv(res$data, "plot_v4/figS3m_2.csv")

# ===== batch6.ipynb cell 940 (code) =====
pdf("plot_v4/Neu_IL1R1_D7_PFS.pdf", width = 6, height = 7)
print(res$plot)
dev.off()

# ===== batch6.ipynb cell 941 (code) =====
ggsave("plot_v4/CD8_D0_PFS.pdf", width = 6, height = 7)

# ===== batch6.ipynb cell 944 (code) =====
plot_stage_cell_survival <- function(
  meta, 
  target_cell, 
  target_stage,
  cell_col = "cell_type3", 
  time_col = "sample_type_rn", 
  patient_col = "patient_id",
  surv_time_col = "os_time",
  surv_status_col = "is_death",
  max_ci_upper_to_show = Inf
) {
  
  suppressPackageStartupMessages({
    library(survival)
    library(survminer)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
  })

  message(sprintf(">>> 步骤 1: 提取阶段 [%s] 的数据并计算细胞比例...", target_stage))
  
  prop_data <- meta %>%
    dplyr::filter(.data[[time_col]] == target_stage) %>%
    dplyr::group_by(.data[[patient_col]]) %>%
    dplyr::summarise(
      Total_Cells = dplyr::n(),
      Target_Cells = sum(.data[[cell_col]] == target_cell, na.rm = TRUE),
      Percentage = (Target_Cells / Total_Cells) * 100,
      Surv_Time = dplyr::first(.data[[surv_time_col]]),
      Raw_Status = dplyr::first(.data[[surv_status_col]]),
      .groups = "drop"
    ) %>%
    dplyr::filter(Total_Cells > 0)
  
  if (nrow(prop_data) < 5) {
    stop("该阶段的有效样本量不足 N < 5，无法进行生存分析。")
  }

  message(">>> 步骤 2: 清洗生存状态并剔除缺失值...")
  
  surv_df <- prop_data %>%
    dplyr::mutate(
      Surv_Time = as.numeric(Surv_Time),
      Surv_Status = dplyr::case_when(
        Raw_Status %in% c("Yes", "yes", "Dead", "Death", "1", 1, TRUE, "是", "死亡", "event", "Event") ~ 1,
        Raw_Status %in% c("No", "no", "Alive", "0", 0, FALSE, "否", "存活", "censored", "Censored") ~ 0,
        TRUE ~ suppressWarnings(as.numeric(Raw_Status))
      )
    ) %>%
    dplyr::filter(
      !is.na(Surv_Time),
      !is.na(Surv_Status),
      !is.na(Percentage),
      is.finite(Percentage)
    )

  if (nrow(surv_df) < 5) {
    stop("具有完整生存信息的样本量不足。")
  }

  message(">>> 步骤 3: 根据细胞比例的中位数分组并拟合生存曲线...")
  
  threshold <- median(surv_df$Percentage, na.rm = TRUE)

  surv_df$Group <- ifelse(
    surv_df$Percentage >= threshold,
    "High",
    "Low"
  )

  surv_df$Group <- factor(
    surv_df$Group,
    levels = c("Low", "High")
  )

  fit <- survival::survfit(
    survival::Surv(Surv_Time, Surv_Status) ~ Group,
    data = surv_df
  )

  # Cox model: High vs Low HR
  cox_group <- tryCatch(
    survival::coxph(
      survival::Surv(Surv_Time, Surv_Status) ~ Group,
      data = surv_df
    ),
    error = function(e) NULL
  )

  if (!is.null(cox_group)) {
    cox_sum <- summary(cox_group)

    cox_HR <- cox_sum$coefficients[1, "exp(coef)"]
    cox_CI_low <- cox_sum$conf.int[1, "lower .95"]
    cox_CI_high <- cox_sum$conf.int[1, "upper .95"]
  } else {
    cox_HR <- NA_real_
    cox_CI_low <- NA_real_
    cox_CI_high <- NA_real_
  }

  # Log-rank test
  surv_diff <- survival::survdiff(
    survival::Surv(Surv_Time, Surv_Status) ~ Group,
    data = surv_df
  )

  logrank_p <- 1 - stats::pchisq(
    surv_diff$chisq,
    df = length(surv_diff$n) - 1
  )

  if (
    !is.na(cox_HR) &&
    !is.na(cox_CI_low) &&
    !is.na(cox_CI_high) &&
    cox_CI_high <= max_ci_upper_to_show
  ) {
    stat_label <- paste0(
      "High vs Low HR = ", sprintf("%.2f", cox_HR),
      "\n95% CI ", sprintf("%.2f", cox_CI_low), "-", sprintf("%.2f", cox_CI_high),
      "\nLog-rank P = ", signif(logrank_p, 3)
    )
  } else {
    stat_label <- paste0(
      "HR not shown: unstable CI",
      "\nLog-rank P = ", signif(logrank_p, 3)
    )
  }

  ylab_text <- ifelse(
    grepl("os", surv_time_col, ignore.case = TRUE),
    "Overall survival probability",
    "Progression-free survival probability"
  )

  p <- survminer::ggsurvplot(
    fit, 
    data = surv_df,
    font.title = c(18, "bold"),
    font.x = 14,
    font.y = 14,
    font.tickslab = 12,
    font.legend = 12,
    pval = FALSE,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.col = "strata",
    palette = c("#377EB8", "#E41A1C"),
    title = paste0("Survival Analysis: ", target_cell, " at ", target_stage),
    subtitle = paste0("Cutoff (median) = ", round(threshold, 2), "%"),
    legend.title = "Proportion",
    legend.labs = c("Low", "High"),
    xlab = surv_time_col,
    ylab = ylab_text,
    ggtheme = theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 12)
      )
  )

  p$plot <- p$plot +
    ggplot2::annotate(
      "text",
      x = max(surv_df$Surv_Time, na.rm = TRUE) * 0.05,
      y = 0.15,
      hjust = 0,
      size = 5,
      label = stat_label
    )

  message(">>> 分析完成！")
  
  return(list(
    plot = p,
    data = surv_df,
    median_cutoff = threshold,
    cox_group = cox_group,
    HR = cox_HR,
    CI_low = cox_CI_low,
    CI_high = cox_CI_high,
    logrank_p = logrank_p
  ))
}

# ===== batch6.ipynb cell 953 (code) =====
immune_programs_targeted <- list(

  # -----------------------------
  # B lineage
  # -----------------------------
  B_cell__BCR_antigen_presentation = list(
    subtypes = c("B_cell", "B_Naive"),
    genes = c(
      "MS4A1", "CD79A", "CD79B", "CD74",
      "HLA-DRA", "HLA-DPA1", "HLA-DPB1",
      "BANK1", "PAX5", "CD40"
    )
  ),

  B_Naive__naive_B_program = list(
    subtypes = c("B_Naive"),
    genes = c(
      "MS4A1", "CD79A", "CD79B", "BANK1",
      "TCL1A", "IGHD", "IGHM", "PAX5", "CD22"
    )
  ),

  Plasma_cells__antibody_secreting = list(
    subtypes = c("Plasma cells"),
    genes = c(
      "PRDM1", "IRF4", "XBP1", "MZB1", "JCHAIN",
      "SDC1", "TNFRSF17", "CD38", "FKBP11",
      "DERL3", "IGHG1", "IGHG3", "IGHM"
    )
  ),


  # -----------------------------
  # CD4 T lineage
  # -----------------------------
  CD4_Naive__naive_memory = list(
    subtypes = c("CD4_Naive"),
    genes = c(
      "CCR7", "SELL", "IL7R", "TCF7", "LEF1",
      "LTB", "MAL", "S1PR1", "BACH2"
    )
  ),

  CD4_Tem__effector_memory = list(
    subtypes = c("CD4_Tem"),
    genes = c(
      "IL7R", "CCL5", "GZMK", "CXCR3",
      "ANXA1", "S100A4", "KLRB1", "IL32"
    )
  ),

  CD4_Th17__Th17_program = list(
    subtypes = c("CD4_Th17"),
    genes = c(
      "RORC", "IL17A", "IL17F", "CCR6",
      "KLRB1", "IL23R", "RORA", "IL22", "CCL20"
    )
  ),

  CD4_Treg__suppressive_program = list(
    subtypes = c("CD4_Treg"),
    genes = c(
      "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
      "TNFRSF18", "ENTPD1", "CCR8", "IL10",
      "TGFB1", "BATF", "PRDM1"
    )
  ),

  CD4_LAG3__checkpoint_exhaustion = list(
    subtypes = c("CD4_LAG3"),
    genes = c(
      "LAG3", "PDCD1", "HAVCR2", "TIGIT",
      "CTLA4", "TOX", "TOX2", "ENTPD1", "CXCL13"
    )
  ),

  CD4_all__Th1_antiviral_helper = list(
    subtypes = c(
      "CD4_LAG3", "CD4_Naive", "CD4_Other",
      "CD4_Tem", "CD4_Th17", "CD4_Treg"
    ),
    genes = c(
      "TBX21", "IFNG", "CXCR3", "IL12RB2",
      "STAT1", "CCL5", "CCR5", "CD40LG"
    )
  ),


  # -----------------------------
  # CD8 T lineage
  # -----------------------------
  CD8_Tem_TRM__cytolytic_effector = list(
    subtypes = c("CD8_Tem", "CD8_TRM", "CD8_NKT", "CD8_Prolif"),
    genes = c(
      "GZMA", "GZMB", "GZMH", "PRF1", "GNLY",
      "NKG7", "CTSW", "CST7", "IFNG"
    )
  ),

  CD8_TRM__tissue_resident_memory = list(
    subtypes = c("CD8_TRM"),
    genes = c(
      "ITGAE", "CXCR6", "ZNF683", "CD69",
      "XCL1", "XCL2", "GZMB", "CXCR3"
    )
  ),

  CD8_Naive__naive_memory = list(
    subtypes = c("CD8_Naive"),
    genes = c(
      "CCR7", "SELL", "IL7R", "TCF7", "LEF1",
      "LTB", "MAL", "S1PR1", "BACH2"
    )
  ),

  CD8_Prolif__proliferation = list(
    subtypes = c("CD8_Prolif"),
    genes = c(
      "MKI67", "TOP2A", "STMN1", "HMGB2",
      "TUBA1B", "TYMS", "PCNA", "MCM5"
    )
  ),

  # CD8_Stress__stress_response = list(
  #   subtypes = c("CD8_Stress"),
  #   genes = c(
  #     "HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1",
  #     "FOS", "JUN", "DUSP1", "IER2", "ATF3"
  #   )
  # ),

  CD8_all__checkpoint_exhaustion = list(
    subtypes = c(
      "CD8_Naive", "CD8_NKT", "CD8_Prolif",
      "CD8_Stress", "CD8_Tem", "CD8_TRM"
    ),
    genes = c(
      "PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4",
      "TOX", "TOX2", "ENTPD1", "CXCL13", "LAYN"
    )
  ),
    
   CD8_all__activation = list(
    subtypes = c(
      "CD8_Naive", "CD8_NKT", "CD8_Prolif",
      "CD8_Stress", "CD8_Tem", "CD8_TRM"
    ),
    genes = c("CD69", "CCR7", "CD27", "BTLA", "CD40LG", "IL2RA", "CD3E", 
                      "CD47", "EOMES", "GNLY", "GZMA", "GZMB", "PRF1", "IFNG", 
                      "CD8A", "CD8B", "FASLG", "LAMP1", "LAG3", "ICOS"
    )
  ),   

    
   CD8_all__cytotoxicity = list(
    subtypes = c(
      "CD8_Naive", "CD8_NKT", "CD8_Prolif",
      "CD8_Stress", "CD8_Tem", "CD8_TRM"
    ),
    genes = c("PRF1", "IFNG", "GNLY", "NKG7", "GZMB", "GZMA",  "CTSW", "CST7", "FASLG",
  "CD8A", "CD8B", "KLRG1", "KLRK1", "LAMP1")
  ),   

   
  CD8_all__chemokine_trafficking = list(
    subtypes = c(
      "CD8_Naive", "CD8_NKT", "CD8_Prolif",
      "CD8_Stress", "CD8_Tem", "CD8_TRM"
    ),
    genes = c(
      "CXCR3", "CCR5", "CCL5",  "ITGAL", "ITGB2", "SELPLG", "CXCR6"
    )
  ),


  # -----------------------------
  # MAIT
  # -----------------------------
  MAIT__innate_like_cytotoxicity = list(
    subtypes = c("MAIT"),
    genes = c(
      "TRAV1-2", "KLRB1", "SLC4A10", "NCR3",
      "CXCR6", "IL7R", "GZMK", "NKG7", "IFNG", "CCL5"
    )
  ),


  # -----------------------------
  # NK lineage
  # -----------------------------
  NK_C2__NK_cytotoxicity = list(
    subtypes = c("NK_C2", "NK_Prolif", "NK_XCL1"),
    genes = c(
      "NKG7", "GNLY", "PRF1", "GZMB", "GZMH",
      "CTSW", "CST7", "KLRD1", "KLRK1",
      "NCR1", "NCR3", "FCGR3A"
    )
  ),


  NK_Prolif__proliferation = list(
    subtypes = c("NK_Prolif"),
    genes = c(
      "MKI67", "TOP2A", "STMN1", "HMGB2",
      "TUBA1B", "TYMS", "PCNA", "MCM5"
    )
  ),

  NK_XCL1__DC_recruiting_NK = list(
    subtypes = c("NK_XCL1"),
    genes = c(
      "XCL1", "XCL2", "CCL5", "IFNG",
      "NKG7", "GNLY", "GZMB", "PRF1"
    )
  ),



  # -----------------------------
  # DC lineage
  # -----------------------------
  DC_cDC1__cross_presentation = list(
    subtypes = c("DC_cDC1"),
    genes = c(
      "CLEC9A", "XCR1", "BATF3", "IRF8", "CADM1",
      "THBD", "WDFY4", "TAP1", "TAP2", "B2M",
      "HLA-A", "HLA-B", "HLA-C", "PSMB8", "PSMB9"
    )
  ),

  DC_cDC2__MHCII_maturation = list(
    subtypes = c("DC_cDC2"),
    genes = c(
      "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
      "CD74", "CIITA", "CD80", "CD86", "CD40",
      "CD83", "CCR7", "LAMP3"
    )
  ),

  DC_all__antiviral_ISG = list(
    subtypes = c("DC_cDC1", "DC_cDC2"),
    genes = c(
      "STAT1", "STAT2", "IRF7", "ISG15",
      "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2",
      "OAS1", "OAS2", "OAS3", "RSAD2",
      "IFI6", "IFI27", "IFI44L"
    )
  ),

  DC_all__viral_sensing = list(
    subtypes = c("DC_cDC1", "DC_cDC2"),
    genes = c(
      "TLR3", "TLR7", "TLR8", "TLR9", "MYD88",
      "TICAM1", "DDX58", "IFIH1", "MAVS",
      "MB21D1", "TMEM173", "TBK1", "IRF3", "IRF7"
    )
  ),


  # -----------------------------
  # Monocyte lineage
  # -----------------------------
  Mono_Classical__inflammatory_monocyte = list(
    subtypes = c("Mono_Classical_1", "Mono_Classical_2"),
    genes = c(
      "S100A8", "S100A9", "S100A12", "LST1",
      "FCN1", "VCAN", "IL1B", "TNF", "CXCL8",
      "CCL2", "NLRP3", "NFKBIA"
    )
  ),

  Mono_IFN__ISG_antiviral = list(
    subtypes = c("Mono_IFN"),
    genes = c(
      "STAT1", "STAT2", "IRF7", "ISG15",
      "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2",
      "OAS1", "OAS2", "OAS3", "RSAD2",
      "IFI6", "IFI27", "IFI44L"
    )
  ),

  Mono_Intermediate__APC_inflammatory = list(
    subtypes = c("Mono_Intermediate"),
    genes = c(
      "CD14", "FCGR3A", "HLA-DRA", "HLA-DRB1",
      "CD74", "LST1", "MS4A7", "IL1B", "TNF"
    )
  ),

  Mono_NonClassical__patrolling_surveillance = list(
    subtypes = c("Mono_NonClassical"),
    genes = c(
      "FCGR3A", "MS4A7", "CX3CR1", "LST1",
      "IFITM3", "TYROBP", "AIF1", "LILRB1",
      "HLA-DRA", "CD74"
    )
  ),

  Mono_all__MDSC_like_positive = list(
    subtypes = c(
      "Mono_Classical_1", "Mono_Classical_2", "Mono_IFN",
      "Mono_Intermediate", "Mono_NonClassical", "Mono_Other"
    ),
    genes = c(
      "CD14", "ITGAM", "CD33", "S100A8", "S100A9",
      "VCAN", "FCN1", "IL1B", "ARG1", "PTGS2", "IL10"
    )
  ),

  Mono_all__HLA_DR_APC = list(
    subtypes = c(
      "Mono_Classical_1", "Mono_Classical_2", "Mono_IFN",
      "Mono_Intermediate", "Mono_NonClassical", "Mono_Other"
    ),
    genes = c("HLA-DRA", "HLA-DRB1", "CD74", "HLA-DPA1", "HLA-DPB1")
  ),


  # -----------------------------
  # Neutrophil lineage
  # -----------------------------
  Neu_CXCR2__mature_neutrophil_chemotaxis = list(
    subtypes = c("Neu_CXCR2"),
    genes = c(
      "CXCR2", "FCGR3B", "CSF3R", "S100A8",
      "S100A9", "MMP9", "LCN2", "CEACAM8"
    )
  ),

  Neu_IFN__ISG_antiviral = list(
    subtypes = c("Neu_IFN"),
    genes = c(
      "STAT1", "STAT2", "IRF7", "ISG15",
      "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2",
      "OAS1", "OAS2", "OAS3", "RSAD2",
      "IFI6", "IFI27", "IFI44L"
    )
  ),

  Neu_IL1R1__inflammatory_neutrophil = list(
    subtypes = c("Neu_IL1R1"),
    genes = c(
      "IL1R1", "IL1B", "CXCL8", "S100A8",
      "S100A9", "S100A12", "NFKBIA", "TNF", "NLRP3"
    )
  ),

  Neu_all__PMN_MDSC_suppressive = list(
    subtypes = c("Neu_CXCR2", "Neu_IFN", "Neu_IL1R1"),
    genes = c(
      "S100A8", "S100A9", "ARG1", "CXCL8",
      "MMP9", "LCN2", "CEACAM8", "FCGR3B",
      "CSF3R", "CYBB", "NCF1", "NCF2"
    )
  )
)

# ===== batch6.ipynb cell 954 (code) =====
calculate_targeted_sample_diffs <- function(
  seu,
  program_list,
  patient_col = "patient_id",
  timepoint_col = "timepoint",
  subtype_col = "celltype_sub",
  baseline_value = "D0",
  treatment_value = "D7",
  assay = NULL,
  slot = "data",
  min_cells_per_patient_program = 3, # 适当放低门槛，避免聚合后样本全被滤掉
  min_genes_present = 3,
  score_method = c("mean_z", "mean_expr")
) {
  score_method <- match.arg(score_method)

  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(dplyr)
    library(tidyr)
    library(tibble)
  })

  if (is.null(assay)) {
    assay <- DefaultAssay(seu)
  }

  # -----------------------------
  # 1. 提取并过滤 D0 和 D7 的细胞
  # -----------------------------
  meta <- seu@meta.data

  required_cols <- c(patient_col, timepoint_col, subtype_col)
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
  }

  valid_timepoints <- c(baseline_value, treatment_value)
  target_time_cells <- rownames(meta)[meta[[timepoint_col]] %in% valid_timepoints]

  if (length(target_time_cells) == 0) {
    stop("未找到任何属于指定时间点的细胞。")
  }

  seu_sub <- subset(seu, cells = target_time_cells)
  meta_sub <- seu_sub@meta.data
  expr <- GetAssayData(seu_sub, assay = assay, slot = slot)

  # -----------------------------
  # 2. 基因表达归一化
  # -----------------------------
  if (score_method == "mean_z") {
    gene_mean <- Matrix::rowMeans(expr)
    gene_sq_mean <- Matrix::rowMeans(expr ^ 2)
    gene_sd <- sqrt(gene_sq_mean - gene_mean ^ 2)
    gene_sd[gene_sd == 0 | is.na(gene_sd)] <- 1

    expr_z <- (expr - gene_mean) / gene_sd
    score_expr <- expr_z
  } else {
    score_expr <- expr
  }

  # -----------------------------
  # 3. 按 Program 计算单细胞得分
  # -----------------------------
  score_records <- list()

  for (program_name in names(program_list)) {
    target_subtypes <- program_list[[program_name]]$subtypes
    genes <- unique(program_list[[program_name]]$genes)
    genes_present <- intersect(genes, rownames(score_expr))

    if (length(genes_present) < min_genes_present) {
      next
    }

    target_cells <- rownames(meta_sub)[meta_sub[[subtype_col]] %in% target_subtypes]
    if (length(target_cells) == 0) next

    program_score <- Matrix::colMeans(score_expr[genes_present, target_cells, drop = FALSE])

    tmp <- data.frame(
      cell_id = target_cells,
      patient_id = meta_sub[target_cells, patient_col, drop = TRUE],
      timepoint = meta_sub[target_time_cells, timepoint_col, drop = TRUE][match(target_cells, target_time_cells)],
      program = program_name,
      score = as.numeric(program_score),
      stringsAsFactors = FALSE
    )
    score_records[[program_name]] <- tmp
  }

  cell_scores <- dplyr::bind_rows(score_records)
  if (nrow(cell_scores) == 0) {
    stop("没有计算出任何有效的特征得分。")
  }

  # -----------------------------
  # 4. 聚合至 Patient 层面
  # -----------------------------
  patient_scores_long <- cell_scores %>%
    dplyr::group_by(patient_id, timepoint, program) %>%
    dplyr::summarise(
      patient_score = mean(score, na.rm = TRUE),
      n_cells = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_cells >= min_cells_per_patient_program)

  # -----------------------------
  # 5. 核心修改：极其稳健的差异表格计算逻辑
  # -----------------------------
  message(">>> 正在进行 D7 vs D0 统计学检验...")
  
  # 转为标准对比宽表
  comparison_df <- patient_scores_long %>%
    dplyr::select(patient_id, timepoint, program, patient_score) %>%
    tidyr::pivot_wider(names_from = timepoint, values_from = patient_score)
  
  # 规避缺少时间列的极端情况
  if (!baseline_value %in% colnames(comparison_df)) comparison_df[[baseline_value]] <- NA_real_
  if (!treatment_value %in% colnames(comparison_df)) comparison_df[[treatment_value]] <- NA_real_

  # 获取所有待计算的 programs
  all_programs <- unique(comparison_df$program)
  res_list <- list()

  for (prog in all_programs) {
    # 提取当前特征的独立子表
    sub_dat <- comparison_df %>% dplyr::filter(program == prog)
    
    vec_base <- sub_dat[[baseline_value]]
    vec_treat <- sub_dat[[treatment_value]]
    
    # 过滤出成对的有效数据
    paired_idx <- which(!is.na(vec_base) & !is.na(vec_treat))
    v_base_paired <- vec_base[paired_idx]
    v_treat_paired <- vec_treat[paired_idx]
    
    # 初始化统计量
    p_unpaired_wilcox <- NA_real_
    p_paired_wilcox   <- NA_real_
    p_paired_ttest    <- NA_real_
    
    # 严格校验 1：非配对检验（剔除NA后必须有样本，且不能完全没有方差变化）
    v_b_clean <- na.omit(vec_base)
    v_t_clean <- na.omit(vec_treat)
    if (length(v_b_clean) >= 2 && length(v_t_clean) >= 2) {
      if (sd(c(v_b_clean, v_t_clean)) > 0) {
        p_unpaired_wilcox <- tryCatch({
          wilcox.test(v_t_clean, v_b_clean, paired = FALSE)$p.value
        }, error = function(e) NA_real_)
      }
    }
    
    # 严格校验 2：配对检验（配对样本量必须 >= 3，且差值不能全部为 0）
    if (length(paired_idx) >= 3) {
      diff_vec <- v_treat_paired - v_b_paired
      if (sd(diff_vec) > 0) {
        p_paired_wilcox <- tryCatch({
          wilcox.test(v_treat_paired, v_b_paired, paired = TRUE)$p.value
        }, error = function(e) NA_real_)
        
        p_paired_ttest <- tryCatch({
          t.test(v_treat_paired, v_b_paired, paired = TRUE)$p.value
        }, error = function(e) NA_real_)
      }
    }
    
    # 组装当前行的结果
    res_list[[prog]] <- data.frame(
      program = prog,
      N_with_Baseline  = length(v_b_clean),
      N_with_Treatment = length(v_t_clean),
      N_Paired_Pairs   = length(paired_idx),
      Mean_Baseline    = mean(v_b_clean, na.rm = TRUE),
      Mean_Treatment   = mean(v_t_clean, na.rm = TRUE),
      Mean_Diff        = mean(v_t_clean, na.rm = TRUE) - mean(v_b_clean, na.rm = TRUE),
      Pval_Unpaired_Wilcox = p_unpaired_wilcox,
      Pval_Paired_Wilcox   = p_paired_wilcox,
      Pval_Paired_Ttest    = p_paired_ttest,
      stringsAsFactors = FALSE
    )
  }

  diff_table <- dplyr::bind_rows(res_list)
  
  # 计算 FDR 校正（仅对非 NA 值的 P 值进行有效校正）
  if (nrow(diff_table) > 0) {
    diff_table <- diff_table %>%
      dplyr::mutate(
        FDR_Unpaired_Wilcox = p.adjust(Pval_Unpaired_Wilcox, method = "BH"),
        FDR_Paired_Wilcox   = p.adjust(Pval_Paired_Wilcox, method = "BH")
      ) %>%
      dplyr::arrange(Pval_Unpaired_Wilcox)
  }

  return(list(
    diff_table = diff_table,
    patient_scores_long = patient_scores_long,
    raw_comparison_df = comparison_df
  ))
}

# ===== batch6.ipynb cell 955 (code) =====
diff_program <- calculate_targeted_sample_diffs (
  seurat_merge,
  immune_programs_targeted,
  patient_col = "orig.ident",
  timepoint_col = "sample_type_rn",
  subtype_col = "cell_type3",
  baseline_value = "D0",
  treatment_value = "D7",      # 新增：对比的时间点
  assay = NULL,
  slot = "data",
  min_cells_per_patient_program = 5, # 适当降低阈值，因为分两个时间点后单点细胞数会变少
  min_genes_present = 3,
  score_method = c("mean_expr")
)

# ===== batch6.ipynb cell 956 (code) =====
plot_program_volcano <- function(
  diff_results, 
  title = "Differential Programs: D7 vs D0",
  diff_cutoff = 0.1,
  p_cutoff = 0.05,
  max_overlaps = 20,
  label_size = 5,
  label_celltypes = NULL,       # 新增：只显示指定细胞类型的 label，例如 c("CD4_LAG3", "CD8_Tem_TRM")
  label_sig_only = FALSE,       # 是否只标注显著点
  celltype_col = NULL           # 如果 diff_results 里有单独的 cell type 列，可指定；否则从 program 前缀提取
) {
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(ggrepel)
  })

  # ============================================================
  # 1. 检查必要列
  # ============================================================
  required_cols <- c("program", "Mean_Diff", "Pval_Unpaired_Wilcox")
  missing_cols <- setdiff(required_cols, colnames(diff_results))
  
  if (length(missing_cols) > 0) {
    stop("输入的数据框缺少必要的列: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.null(celltype_col) && !(celltype_col %in% colnames(diff_results))) {
    stop("指定的 celltype_col 不存在于 diff_results 中: ", celltype_col)
  }

  # ============================================================
  # 2. 生成标签格式
  #    CD4_LAG3__checkpoint_exhaustion
  #    -> CD4_LAG3 | checkpoint_exhaustion
  # ============================================================
  df <- diff_results %>%
    dplyr::mutate(
      program = as.character(program),
      program_label = gsub("__", " | ", program, fixed = TRUE)
    )

  # ============================================================
  # 3. 提取 cell type
  #    默认从 program 的 "__" 前面提取
  # ============================================================
  if (is.null(celltype_col)) {
    df <- df %>%
      dplyr::mutate(
        program_celltype = sub("__.*$", "", program)
      )
  } else {
    df <- df %>%
      dplyr::mutate(
        program_celltype = as.character(.data[[celltype_col]])
      )
  }

  # ============================================================
  # 4. 构建显著性分类
  # ============================================================
  up_label   <- paste0("p<", p_cutoff, " & Diff>", diff_cutoff)
  down_label <- paste0("p<", p_cutoff, " & Diff< -", diff_cutoff)
  abs_label  <- paste0("|Diff|>", diff_cutoff)
  ns_label   <- "Not significant"

  df <- df %>%
    dplyr::mutate(
      Pval_Unpaired_Wilcox = as.numeric(Pval_Unpaired_Wilcox),
      Mean_Diff = as.numeric(Mean_Diff),
      neg_log10_p = -log10(pmax(Pval_Unpaired_Wilcox, .Machine$double.xmin)),
      
      significance = dplyr::case_when(
        Pval_Unpaired_Wilcox < p_cutoff & Mean_Diff > diff_cutoff ~ up_label,
        Pval_Unpaired_Wilcox < p_cutoff & Mean_Diff < -diff_cutoff ~ down_label,
        abs(Mean_Diff) > diff_cutoff ~ abs_label,
        TRUE ~ ns_label
      )
    )

  # ============================================================
  # 5. 控制哪些点显示 label
  # ============================================================
  if (is.null(label_celltypes)) {
    df <- df %>%
      dplyr::mutate(label_flag = TRUE)
  } else {
    df <- df %>%
      dplyr::mutate(
        label_flag = program_celltype %in% label_celltypes
      )
  }

  if (isTRUE(label_sig_only)) {
    df <- df %>%
      dplyr::mutate(
        label_flag = label_flag &
          Pval_Unpaired_Wilcox < p_cutoff &
          abs(Mean_Diff) > diff_cutoff
      )
  }

  label_df <- df %>%
    dplyr::filter(label_flag)

  # ============================================================
  # 6. 颜色映射
  # ============================================================
  color_values <- c(
    "red", "blue", "orange", "gray"
  )
  
  names(color_values) <- c(
    up_label,
    down_label,
    abs_label,
    ns_label
  )

  # ============================================================
  # 7. 绘图
  # ============================================================
  p <- ggplot(
    df,
    aes(x = Mean_Diff, y = neg_log10_p)
  ) +
    geom_point(
      aes(color = significance, size = abs(Mean_Diff)),
      alpha = 0.8
    ) +
    geom_vline(
      xintercept = c(-diff_cutoff, diff_cutoff),
      linetype = "dashed",
      color = "gray50"
    ) +
    geom_hline(
      yintercept = -log10(p_cutoff),
      linetype = "dashed",
      color = "gray50"
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = program_label),
      size = label_size,
      box.padding = 0.5,
      point.padding = 0.1,
      max.overlaps = max_overlaps,
      segment.color = "grey50"
    ) +
    scale_color_manual(values = color_values) +
    labs(
      title = title,
      x = "Mean Score Difference (D7 - D0)",
      y = "-Log10(Unpaired Wilcoxon P-value)",
      color = "Significance",
      size = "|Mean Diff|"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 14),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 16),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )

  return(p)
}

# ===== batch6.ipynb cell 957 (code) =====
# 2. 设置画布尺寸（完美匹配你的 options 设定）
options(repr.plot.width = 6, repr.plot.height = 6)

# 3. 一键出图 (可通过 diff_cutoff 调节两边虚线的位置)
volcano_plot <- plot_program_volcano(
  diff_results = diff_program$diff_table, 
  title = "Immune Program Changes: D7 vs D0",
  diff_cutoff = 0.05,
  p_cutoff = 0.05,
  label_celltypes = c("CD4_LAG3", "CD8_Prolif", "NK_Prolif", "NK_XCL1", "Plamsa_cells", "CD4_Treg", "DC_cDC2", "Mono_all"),
  max_overlaps = 20,
  label_size = 5
)

# ===== batch6.ipynb cell 961 (code) =====
ggsave("plot_v4/immune Program changes D7 vs D0.pdf", width = 6, height = 6)

# ===== batch6.ipynb cell 963 (code) =====
calculate_targeted_patient_scores_all_samples <- function(
  seu,
  program_list,
  patient_col = "patient_id",
  timepoint_col = "sample_type_rn",
  subtype_col = "cell_type3",
  sample_col = NULL,              # 可选：如 "orig.ident"；如果不提供，则用 patient_id + timepoint 作为 sample_id
  assay = NULL,
  slot = "data",
  min_cells_per_patient_program = 10,
  min_genes_present = 3,
  score_method = c("AddModuleScore", "mean_z", "mean_expr")
) {

  score_method <- match.arg(score_method)

  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(dplyr)
    library(tidyr)
    library(tibble)
  })

  if (is.null(assay)) {
    assay <- DefaultAssay(seu)
  }

  # ============================================================
  # 0. helper: compatible with Seurat v4/v5
  # ============================================================
  get_assay_data_safe <- function(obj, assay, slot) {
    tryCatch(
      {
        Seurat::GetAssayData(obj, assay = assay, slot = slot)
      },
      error = function(e) {
        Seurat::GetAssayData(obj, assay = assay, layer = slot)
      }
    )
  }

  # ============================================================
  # 1. Check metadata
  # ============================================================
  meta <- seu@meta.data

  required_cols <- c(patient_col, timepoint_col, subtype_col)

  if (!is.null(sample_col)) {
    required_cols <- unique(c(required_cols, sample_col))
  }

  missing_cols <- setdiff(required_cols, colnames(meta))

  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
  }

  # ============================================================
  # 2. Use all cells, no timepoint subsetting
  # ============================================================
  meta_use <- meta %>%
    dplyr::mutate(
      patient_id_internal = as.character(.data[[patient_col]]),
      timepoint_internal = as.character(.data[[timepoint_col]]),
      subtype_internal = as.character(.data[[subtype_col]])
    )

  if (is.null(sample_col)) {
    meta_use$sample_id_internal <- paste(
      meta_use$patient_id_internal,
      meta_use$timepoint_internal,
      sep = "__"
    )
  } else {
    meta_use$sample_id_internal <- as.character(meta_use[[sample_col]])
  }

  valid_cells <- rownames(meta_use)[
    !is.na(meta_use$patient_id_internal) &
      !is.na(meta_use$timepoint_internal) &
      !is.na(meta_use$subtype_internal) &
      !is.na(meta_use$sample_id_internal)
  ]

  if (length(valid_cells) == 0) {
    stop("No valid cells found after checking metadata.")
  }

  seu_use <- subset(seu, cells = valid_cells)
  meta_use <- seu_use@meta.data

  meta_use$patient_id_internal <- as.character(meta_use[[patient_col]])
  meta_use$timepoint_internal <- as.character(meta_use[[timepoint_col]])
  meta_use$subtype_internal <- as.character(meta_use[[subtype_col]])

  if (is.null(sample_col)) {
    meta_use$sample_id_internal <- paste(
      meta_use$patient_id_internal,
      meta_use$timepoint_internal,
      sep = "__"
    )
  } else {
    meta_use$sample_id_internal <- as.character(meta_use[[sample_col]])
  }

  seu_use@meta.data <- meta_use

  message(">>> Cells used: ", ncol(seu_use))
  message(">>> Patients: ", length(unique(meta_use$patient_id_internal)))
  message(">>> Timepoints: ", paste(sort(unique(meta_use$timepoint_internal)), collapse = ", "))
  message(">>> Samples: ", length(unique(meta_use$sample_id_internal)))
  message(">>> Programs to score: ", length(program_list))

  # ============================================================
  # 3. Calculate program scores
  #    每个 program 单独 subset target cells 后计算 score
  # ============================================================
  score_records <- list()

  for (program_name in names(program_list)) {

    message(">>> Scoring program: ", program_name)

    target_subtypes <- program_list[[program_name]]$subtypes
    genes <- unique(program_list[[program_name]]$genes)

    # ------------------------------------------------------------
    # 3.1 target cells: all timepoints + compatible cell subtypes
    # ------------------------------------------------------------
    meta_prog <- seu_use@meta.data

    target_cells <- rownames(meta_prog)[
      !is.na(meta_prog$subtype_internal) &
        meta_prog$subtype_internal %in% target_subtypes
    ]

    if (length(target_cells) == 0) {
      message("Skip ", program_name, ": no target cells found.")
      next
    }

    # ------------------------------------------------------------
    # 3.2 subset object for this program
    # ------------------------------------------------------------
    sub_obj <- subset(seu_use, cells = target_cells)

    valid_genes <- intersect(genes, rownames(sub_obj))
    missing_genes <- setdiff(genes, valid_genes)

    if (length(valid_genes) < min_genes_present) {
      message(
        "Skip ", program_name,
        ": only ", length(valid_genes),
        " genes present. Required >= ", min_genes_present
      )
      next
    }

    if (length(missing_genes) > 0) {
      message(
        "Program ", program_name,
        ": missing genes removed: ",
        paste(missing_genes, collapse = ", ")
      )
    }

    # ------------------------------------------------------------
    # 3.3 Score calculation
    # ------------------------------------------------------------
    if (score_method == "AddModuleScore") {

      sub_obj <- Seurat::AddModuleScore(
        object = sub_obj,
        features = list(valid_genes),
        name = "ModuleScore",
        assay = assay
      )

      score_col_name <- "ModuleScore1"

      sc_data <- Seurat::FetchData(
        sub_obj,
        vars = c(
          score_col_name,
          patient_col,
          timepoint_col,
          subtype_col,
          "sample_id_internal"
        )
      )

      sc_data <- sc_data %>%
        tibble::rownames_to_column("cell_id") %>%
        dplyr::rename(
          score = dplyr::all_of(score_col_name),
          patient_id = dplyr::all_of(patient_col),
          timepoint = dplyr::all_of(timepoint_col),
          celltype_sub = dplyr::all_of(subtype_col),
          sample_id = sample_id_internal
        )
    }

    if (score_method %in% c("mean_z", "mean_expr")) {

      expr <- get_assay_data_safe(
        obj = sub_obj,
        assay = assay,
        slot = slot
      )

      expr_valid <- expr[valid_genes, , drop = FALSE]

      if (score_method == "mean_z") {

        gene_mean <- Matrix::rowMeans(expr_valid)
        gene_sq_mean <- Matrix::rowMeans(expr_valid ^ 2)
        gene_sd <- sqrt(gene_sq_mean - gene_mean ^ 2)
        gene_sd[gene_sd == 0 | is.na(gene_sd)] <- 1

        expr_use <- expr_valid
        expr_use <- expr_use - gene_mean
        expr_use <- expr_use / gene_sd

      } else {

        expr_use <- expr_valid
      }

      cell_score <- Matrix::colMeans(expr_use)

      sc_data <- Seurat::FetchData(
        sub_obj,
        vars = c(
          patient_col,
          timepoint_col,
          subtype_col,
          "sample_id_internal"
        )
      )

      sc_data$score <- as.numeric(cell_score[rownames(sc_data)])

      sc_data <- sc_data %>%
        tibble::rownames_to_column("cell_id") %>%
        dplyr::rename(
          patient_id = dplyr::all_of(patient_col),
          timepoint = dplyr::all_of(timepoint_col),
          celltype_sub = dplyr::all_of(subtype_col),
          sample_id = sample_id_internal
        ) %>%
        dplyr::select(
          cell_id,
          sample_id,
          patient_id,
          timepoint,
          celltype_sub,
          score
        )
    }

    sc_data <- sc_data %>%
      dplyr::mutate(
        program = program_name,
        n_genes = length(valid_genes),
        n_target_cells_total = length(target_cells),
        target_subtypes = paste(target_subtypes, collapse = ";"),
        valid_genes = paste(valid_genes, collapse = ";"),
        score_method = score_method
      )

    score_records[[program_name]] <- sc_data
  }

  cell_scores <- dplyr::bind_rows(score_records)

  if (nrow(cell_scores) == 0) {
    stop("No valid program scores were calculated.")
  }

  # ============================================================
  # 4. Aggregate to sample level
  #    sample_score = mean(Cell_Score)
  #    默认 sample_id = patient_id__timepoint
  # ============================================================
  patient_scores_long <- cell_scores %>%
    dplyr::group_by(sample_id, patient_id, timepoint, program) %>%
    dplyr::summarise(
      patient_score = mean(score, na.rm = TRUE),
      n_cells = dplyr::n(),
      n_genes = max(n_genes, na.rm = TRUE),
      target_subtypes = dplyr::first(target_subtypes),
      valid_genes = dplyr::first(valid_genes),
      score_method = dplyr::first(score_method),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_cells >= min_cells_per_patient_program)

  if (nrow(patient_scores_long) == 0) {
    stop(
      "No sample-level program scores remained after filtering by ",
      "min_cells_per_patient_program = ",
      min_cells_per_patient_program
    )
  }

  # ============================================================
  # 5. Wide sample × program matrix
  # ============================================================
  patient_program_matrix <- patient_scores_long %>%
    dplyr::select(sample_id, program, patient_score) %>%
    tidyr::pivot_wider(
      names_from = program,
      values_from = patient_score
    ) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  patient_program_matrix_t <- t(patient_program_matrix)

  # ============================================================
  # 6. Sample metadata table
  # ============================================================
  sample_metadata <- patient_scores_long %>%
    dplyr::select(sample_id, patient_id, timepoint) %>%
    dplyr::distinct()

  # ============================================================
  # 7. Summary table
  # ============================================================
  program_summary <- patient_scores_long %>%
    dplyr::group_by(program) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      n_patients = dplyr::n_distinct(patient_id),
      n_timepoints = dplyr::n_distinct(timepoint),
      median_cells_per_sample = median(n_cells, na.rm = TRUE),
      min_cells_per_sample = min(n_cells, na.rm = TRUE),
      max_cells_per_sample = max(n_cells, na.rm = TRUE),
      n_genes = max(n_genes, na.rm = TRUE),
      target_subtypes = dplyr::first(target_subtypes),
      valid_genes = dplyr::first(valid_genes),
      score_method = dplyr::first(score_method),
      .groups = "drop"
    )

  message(">>> Finished.")
  message(">>> Valid cell-level scores: ", nrow(cell_scores))
  message(">>> Valid sample-program scores: ", nrow(patient_scores_long))
  message(">>> Programs retained: ", ncol(patient_program_matrix))
  message(">>> Samples retained: ", nrow(patient_program_matrix))

  return(list(
    cell_scores = cell_scores,
    patient_scores_long = patient_scores_long,
    patient_program_matrix = patient_program_matrix,
    patient_program_matrix_t = patient_program_matrix_t,
    sample_metadata = sample_metadata,
    program_summary = program_summary,
    score_method = score_method,
    assay = assay,
    slot = slot
  ))
}

# ===== batch6.ipynb cell 965 (code) =====
score_res_all <- calculate_targeted_patient_scores_all_samples(
  seu = seurat_merge,
  program_list = immune_programs_targeted,
  patient_col = "patient_id",
  timepoint_col = "sample_type_rn",
  subtype_col = "cell_type3",
  sample_col = "orig.ident",
  assay = "RNA",
  slot = "data",
  min_cells_per_patient_program = 5,
  min_genes_present = 3,
  score_method = "AddModuleScore"
)

# ===== batch6.ipynb cell 966 (code) =====
plot_program_boxplot_from_score_res <- function(
  score_res,
  programs,
  group_col = "sample_type_rn",
  sample_col = "sample_id",
  patient_col = "patient_id",
  score_col = "patient_score",

  group_values = NULL,
  group_labels = NULL,
  comparisons = NULL,

  test_method = c("wilcox.test", "t.test"),
  alternative = c("two.sided", "less", "greater"),
  paired = FALSE,
  pair_col = "patient_id",
  line = FALSE,

  p_adjust_method = "BH",
  show_p_adj = FALSE,
p_label_type = c("p_exact", "p_plain", "p", "star_p", "star"),
p_digits = 4,

  show_points = TRUE,
  point_alpha = 0.85,
  point_size = 2.8,
  jitter_width = 0.10,
  box_width = 0.58,
  box_alpha = 0.85,
  bracket_size = 0.55,
  p_label_size = 4.2,

  facet_ncol = NULL,
  facet_nrow = NULL,
  scales = "free_y",

  group_colors = NULL,

  y_label = "Mean Module Score",
  x_label = "Treatment Stage",
  title = NULL,
  subtitle = NULL,

  legend_position = "bottom",
  strip_text_size = 13,
  axis_text_size = 12,
  axis_text_x_size = 12,
  axis_title_size = 14,
  legend_text_size = 12,
  legend_title_size = 12,
  title_size = 15,
  subtitle_size = 11,

  output_pdf = NULL,
  width = 9,
  height = 4.5
) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(ggpubr)
  })

  test_method <- match.arg(test_method)
  alternative <- match.arg(alternative)
  p_label_type <- match.arg(p_label_type)

  # ============================================================
  # 1. Extract score table
  # ============================================================
  if (is.list(score_res) && "patient_scores_long" %in% names(score_res)) {
    score_df <- score_res$patient_scores_long
  } else if (is.data.frame(score_res)) {
    score_df <- score_res
  } else {
    stop("score_res 必须是 score_res_all 结果，或包含 sample-level score 的 data.frame。")
  }

  # 兼容你之前函数输出的 timepoint 列
  if (!(group_col %in% colnames(score_df))) {
    if (group_col == "sample_type_rn" && "timepoint" %in% colnames(score_df)) {
      message("group_col = 'sample_type_rn' 不在 score 表中，自动使用 'timepoint' 列。")
      group_col_use <- "timepoint"
    } else {
      stop("score_res$patient_scores_long 中不存在 group_col: ", group_col)
    }
  } else {
    group_col_use <- group_col
  }

  required_cols <- c("program", score_col, group_col_use, patient_col)

  if (!is.null(sample_col) && sample_col %in% colnames(score_df)) {
    required_cols <- unique(c(required_cols, sample_col))
  }

  if (paired || line) {
    if (!(pair_col %in% colnames(score_df))) {
      stop("paired = TRUE 或 line = TRUE 时，score 表中必须包含 pair_col: ", pair_col)
    }
    required_cols <- unique(c(required_cols, pair_col))
  }

  missing_cols <- setdiff(required_cols, colnames(score_df))

  if (length(missing_cols) > 0) {
    stop("score 表中缺少以下列: ", paste(missing_cols, collapse = ", "))
  }

  programs <- as.character(programs)

  missing_programs <- setdiff(programs, unique(score_df$program))

  if (length(missing_programs) > 0) {
    stop(
      "以下 programs 不在 score_res_all$patient_scores_long 中: ",
      paste(missing_programs, collapse = ", ")
    )
  }

  # ============================================================
  # 2. Prepare plot data
  # ============================================================
  plot_df <- score_df %>%
    dplyr::filter(program %in% programs) %>%
    dplyr::mutate(
      group_raw = as.character(.data[[group_col_use]]),
      score = as.numeric(.data[[score_col]]),
      patient_id_plot = as.character(.data[[patient_col]])
    ) %>%
    dplyr::filter(
      !is.na(group_raw),
      !is.na(score),
      is.finite(score)
    )

  if (paired || line) {
    plot_df <- plot_df %>%
      dplyr::mutate(pair_id_plot = as.character(.data[[pair_col]]))
  } else {
    plot_df$pair_id_plot <- plot_df$patient_id_plot
  }

  if (!is.null(group_values)) {
    group_values <- as.character(group_values)

    plot_df <- plot_df %>%
      dplyr::filter(group_raw %in% group_values)
  } else {
    group_values <- unique(plot_df$group_raw)
  }

  if (is.null(group_labels)) {
    group_labels <- group_values
  }

  if (length(group_values) != length(group_labels)) {
    stop("group_values 和 group_labels 长度必须一致。")
  }

  plot_df <- plot_df %>%
    dplyr::mutate(
      group = factor(
        group_raw,
        levels = group_values,
        labels = group_labels
      ),
      program_label = gsub("__", " | ", program, fixed = TRUE),
      program_label = factor(
        program_label,
        levels = gsub("__", " | ", programs, fixed = TRUE)
      )
    ) %>%
    dplyr::filter(!is.na(group))

  if (nrow(plot_df) == 0) {
    stop("筛选 programs 和 group 后没有可用于绘图的数据。")
  }

  # ============================================================
  # 3. Comparisons
  # ============================================================
  if (is.null(comparisons)) {
    if (length(group_labels) == 2) {
      comparisons <- list(c(group_labels[1], group_labels[2]))
    } else {
      comparisons <- utils::combn(group_labels, 2, simplify = FALSE)
    }
  }

  comparisons <- lapply(comparisons, as.character)

  invalid_comparisons <- comparisons[
    !vapply(
      comparisons,
      function(x) all(x %in% group_labels),
      logical(1)
    )
  ]

  if (length(invalid_comparisons) > 0) {
    stop(
      "comparisons 中存在不属于 group_labels 的分组名。当前 group_labels 为: ",
      paste(group_labels, collapse = ", ")
    )
  }

  # ============================================================
  # 4. P value helper
  # ============================================================
format_p <- function(p, prefix = "P", digits = 4) {
  dplyr::case_when(
    is.na(p) ~ paste0(prefix, "=NA"),
    TRUE ~ paste0(
      prefix,
      "=",
      formatC(p, format = "fg", digits = digits)
    )
  )
}

format_p_plain <- function(p, digits = 4) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    TRUE ~ formatC(p, format = "fg", digits = digits)
  )
}

    
  format_sig <- function(p) {
    dplyr::case_when(
      is.na(p) ~ "NA",
      p < 0.001 ~ "***",
      p < 0.01 ~ "**",
      p < 0.05 ~ "*",
      p < 0.1 ~ "†",
      TRUE ~ "ns"
    )
  }

  run_one_test <- function(df, g1, g2) {

    df_sub <- df %>%
      dplyr::filter(as.character(group) %in% c(g1, g2)) %>%
      dplyr::mutate(group = droplevels(group))

    if (!paired) {

      x1 <- df_sub$score[as.character(df_sub$group) == g1]
      x2 <- df_sub$score[as.character(df_sub$group) == g2]

      if (length(x1) < 2 || length(x2) < 2) {
        return(NA_real_)
      }

      pval <- tryCatch(
        {
          if (test_method == "t.test") {
            stats::t.test(
              x2,
              x1,
              alternative = alternative,
              paired = FALSE
            )$p.value
          } else {
            stats::wilcox.test(
              x2,
              x1,
              alternative = alternative,
              paired = FALSE,
              exact = FALSE
            )$p.value
          }
        },
        error = function(e) NA_real_
      )

      return(pval)
    }

    if (paired) {

      wide_df <- df_sub %>%
        dplyr::select(pair_id_plot, group, score) %>%
        dplyr::mutate(group = as.character(group)) %>%
        dplyr::distinct(pair_id_plot, group, .keep_all = TRUE) %>%
        tidyr::pivot_wider(
          names_from = group,
          values_from = score
        )

      if (!(g1 %in% colnames(wide_df)) || !(g2 %in% colnames(wide_df))) {
        return(NA_real_)
      }

      x1 <- wide_df[[g1]]
      x2 <- wide_df[[g2]]

      keep <- !is.na(x1) & !is.na(x2)

      x1 <- x1[keep]
      x2 <- x2[keep]

      if (length(x1) < 2 || length(x2) < 2) {
        return(NA_real_)
      }

      pval <- tryCatch(
        {
          if (test_method == "t.test") {
            stats::t.test(
              x2,
              x1,
              alternative = alternative,
              paired = TRUE
            )$p.value
          } else {
            stats::wilcox.test(
              x2,
              x1,
              alternative = alternative,
              paired = TRUE,
              exact = FALSE
            )$p.value
          }
        },
        error = function(e) NA_real_
      )

      return(pval)
    }
  }

  # ============================================================
  # 5. P value calculation per program
  # ============================================================
  stat_df <- plot_df %>%
    dplyr::group_by(program, program_label) %>%
    dplyr::group_modify(function(df_prog, key) {

      out <- lapply(seq_along(comparisons), function(i) {

        comp <- comparisons[[i]]

        pval <- run_one_test(
          df = df_prog,
          g1 = comp[1],
          g2 = comp[2]
        )

        data.frame(
          group1 = comp[1],
          group2 = comp[2],
          comparison_id = i,
          Pval = pval,
          stringsAsFactors = FALSE
        )
      })

      dplyr::bind_rows(out)
    }) %>%
    dplyr::ungroup()

  stat_df <- stat_df %>%
    dplyr::mutate(
      P_adj = stats::p.adjust(Pval, method = p_adjust_method)
    )

  if (isTRUE(show_p_adj)) {
    stat_df <- stat_df %>%
      dplyr::mutate(
        P_show = P_adj,
        P_prefix = "FDR"
      )
  } else {
    stat_df <- stat_df %>%
      dplyr::mutate(
        P_show = Pval,
        P_prefix = "P"
      )
  }

  stat_df <- stat_df %>%
    dplyr::mutate(
    P_label = format_p(P_show, prefix = P_prefix, digits = p_digits),
    P_label_plain = format_p_plain(P_show, digits = p_digits),
    Sig_label = format_sig(P_show),
    label = dplyr::case_when(
      p_label_type == "p_exact" ~ P_label,
      p_label_type == "p_plain" ~ P_label_plain,
      p_label_type == "p" ~ P_label,
      p_label_type == "star_p" ~ paste0(Sig_label, "\n", P_label),
      p_label_type == "star" ~ Sig_label,
      TRUE ~ P_label
    ),
      program_label = factor(
        as.character(program_label),
        levels = levels(plot_df$program_label)
      )
    )

  # ============================================================
  # 6. y position
  # ============================================================
  y_pos_df <- plot_df %>%
    dplyr::group_by(program, program_label) %>%
    dplyr::summarise(
      y_min = min(score, na.rm = TRUE),
      y_max = max(score, na.rm = TRUE),
      y_range = y_max - y_min,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      y_range = ifelse(
        is.na(y_range) | y_range == 0,
        abs(y_max) * 0.2 + 1,
        y_range
      )
    )

  stat_df <- stat_df %>%
    dplyr::left_join(
      y_pos_df,
      by = c("program", "program_label")
    ) %>%
    dplyr::mutate(
      y.position = y_max + comparison_id * 0.18 * y_range
    )

  # ============================================================
  # 7. Colors
  # ============================================================
  if (is.null(group_colors)) {
    default_cols <- c(
      "D0" = "#E76F73",
      "D7" = "#4C78A8",
      "postICI" = "#45B7A5",
      "PostICI" = "#45B7A5",
      "postICB" = "#45B7A5",
      "PFS<6" = "#4C78A8",
      "PFS>6" = "#E76F73"
    )

    group_colors <- default_cols[group_labels]

    if (any(is.na(group_colors))) {
      fallback_cols <- c(
        "#E76F73", "#4C78A8", "#45B7A5",
        "#984EA3", "#FF7F00", "#7A7A7A"
      )
      group_colors <- fallback_cols[seq_along(group_labels)]
      names(group_colors) <- group_labels
    }
  }

  # ============================================================
  # 8. Plot
  # ============================================================
  if (is.null(title)) {
    title <- "Sample-level immune program score by treatment stage"
  }

  if (is.null(subtitle)) {
    subtitle <- NULL
  }

  p <- ggplot(
    plot_df,
    aes(
      x = group,
      y = score
    )
  ) +
    geom_boxplot(
      aes(color = group),
      outlier.shape = NA,
      width = box_width,
      alpha = box_alpha,
      linewidth = 0.75
    )

  if (isTRUE(line)) {
    p <- p +
      geom_line(
        aes(group = pair_id_plot),
        color = "grey70",
        alpha = 0.75,
        linewidth = 0.55
      )
  }

  if (isTRUE(show_points)) {
    p <- p +
      geom_jitter(
        aes(color = group),
        width = jitter_width,
        size = point_size,
        alpha = point_alpha
      )
  }

  p <- p +
    ggpubr::stat_pvalue_manual(
      stat_df,
      label = "label",
      xmin = "group1",
      xmax = "group2",
      y.position = "y.position",
      tip.length = 0.01,
      bracket.size = bracket_size,
      size = p_label_size,
      inherit.aes = FALSE,
      step.group.by = "program_label"
    ) +
    facet_wrap(
      ~ program_label,
      nrow = facet_nrow,
      ncol = facet_ncol,
      scales = scales
    ) +
    scale_color_manual(values = group_colors) +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_label,
      y = y_label,
      color = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = title_size
      ),
      plot.subtitle = element_text(
        hjust = 0.5,
        size = subtitle_size
      ),
      strip.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.7
      ),
      strip.text = element_text(
        size = strip_text_size,
        face = "bold",
        color = "black"
      ),
      axis.text = element_text(
        size = axis_text_size,
        color = "black"
      ),
axis.text.x = element_text(
  size = axis_text_x_size,
  angle = 30,
  hjust = 1,
  vjust = 1,
  color = "black"
),
      axis.title = element_text(
        size = axis_title_size,
        face = "bold",
        color = "black"
      ),
      legend.position = legend_position,
      legend.title = element_text(
        size = legend_title_size,
        face = "bold"
      ),
      legend.text = element_text(size = legend_text_size),
      panel.grid.major = element_line(
        color = "grey92",
        linewidth = 0.35
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      )
    )

  # ============================================================
  # 9. Output
  # ============================================================
  if (!is.null(output_pdf)) {
    grDevices::pdf(
      output_pdf,
      width = width,
      height = height
    )
    print(p)
    grDevices::dev.off()
  }

  return(list(
    plot = p,
    plot_data = plot_df,
    stat_data = stat_df,
    programs = programs,
    group_values = group_values,
    group_labels = group_labels,
    group_col = group_col_use,
    test_method = test_method,
    paired = paired,
    p_adjust_method = p_adjust_method,
    show_p_adj = show_p_adj
  ))
}

# ===== batch6.ipynb cell 976 (code) =====
options(repr.plot.width = 12.5, repr.plot.height = 5.8)

box_res <- plot_program_boxplot_from_score_res(
  score_res = score_res_all,
  programs = c(
    "CD8_all__cytotoxicity",
    "CD8_all__activation",
    "CD4_LAG3__checkpoint_exhaustion",
      "CD4_Treg__suppressive_program"
  ),
  group_col = "sample_type_rn",
  group_values = c("D0", "D7", "postICI"),
  group_labels = c("D0", "D7", "postICI"),
  comparisons = list(
    c("D0", "D7"),
    c("D0", "postICI"),
    c("D7", "postICI")
  ),
    box_alpha = 1,
    point_alpha = 1,
     alternative = c("greater"),
  test_method = "t.test",
  show_p_adj = FALSE,
    point_size = 1.5,
  p_label_type = "p",
  facet_ncol = 4,
  legend_position = "bottom",
  title = NULL,
  subtitle = NULL,
  y_label = "Mean Module Score",
  x_label = "Treatment Stage"
)

box_res$plot

# ===== batch6.ipynb cell 977 (code) =====
ggsave("plot_v4/Fig3b.pdf", width = 12.5, height = 5.8)

# ===== batch6.ipynb cell 978 (code) =====
calculate_targeted_patient_scores <- function(
  seu,
  program_list,
  patient_col = "patient_id",
  timepoint_col = "sample_type_rn",
  subtype_col = "cell_type3",
  baseline_value = "D7",
  assay = NULL,
  slot = "data",
  min_cells_per_patient_program = 10,
  min_genes_present = 3,
  score_method = c("AddModuleScore", "mean_z", "mean_expr")
) {

  score_method <- match.arg(score_method)

  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(dplyr)
    library(tidyr)
    library(tibble)
  })

  if (is.null(assay)) {
    assay <- DefaultAssay(seu)
  }

  # ============================================================
  # 0. helper: compatible with Seurat v4/v5
  # ============================================================
  get_assay_data_safe <- function(obj, assay, slot) {
    tryCatch(
      {
        Seurat::GetAssayData(obj, assay = assay, slot = slot)
      },
      error = function(e) {
        Seurat::GetAssayData(obj, assay = assay, layer = slot)
      }
    )
  }

  # ============================================================
  # 1. Check metadata
  # ============================================================
  meta <- seu@meta.data

  required_cols <- c(patient_col, timepoint_col, subtype_col)
  missing_cols <- setdiff(required_cols, colnames(meta))

  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
  }

  # ============================================================
  # 2. Subset baseline cells
  # ============================================================
  baseline_cells <- rownames(meta)[
    !is.na(meta[[timepoint_col]]) &
      meta[[timepoint_col]] == baseline_value
  ]

  if (length(baseline_cells) == 0) {
    stop("No cells found for baseline_value = ", baseline_value)
  }

  seu_base <- subset(seu, cells = baseline_cells)
  meta_base <- seu_base@meta.data

  message(">>> Baseline cells: ", ncol(seu_base))
  message(">>> Baseline patients: ", length(unique(meta_base[[patient_col]])))
  message(">>> Programs to score: ", length(program_list))

  # ============================================================
  # 3. Calculate program scores
  #    与 plot_geneset_pseudobulk_survival_v3 保持一致：
  #    每个 program 单独 subset target cells 后计算 score
  # ============================================================
  score_records <- list()

  for (program_name in names(program_list)) {

    message(">>> Scoring program: ", program_name)

    target_subtypes <- program_list[[program_name]]$subtypes
    genes <- unique(program_list[[program_name]]$genes)

    # ------------------------------------------------------------
    # 3.1 target cells: baseline + compatible cell subtypes
    # ------------------------------------------------------------
    target_cells <- rownames(meta_base)[
      !is.na(meta_base[[subtype_col]]) &
        meta_base[[subtype_col]] %in% target_subtypes
    ]

    if (length(target_cells) == 0) {
      message("Skip ", program_name, ": no target cells found.")
      next
    }

    # ------------------------------------------------------------
    # 3.2 subset object for this program
    # ------------------------------------------------------------
    sub_obj <- subset(seu_base, cells = target_cells)

    valid_genes <- intersect(genes, rownames(sub_obj))
    missing_genes <- setdiff(genes, valid_genes)

    if (length(valid_genes) < min_genes_present) {
      message(
        "Skip ", program_name,
        ": only ", length(valid_genes),
        " genes present. Required >= ", min_genes_present
      )
      next
    }

    if (length(missing_genes) > 0) {
      message(
        "Program ", program_name,
        ": missing genes removed: ",
        paste(missing_genes, collapse = ", ")
      )
    }

    # ------------------------------------------------------------
    # 3.3 Score calculation
    # ------------------------------------------------------------
    if (score_method == "AddModuleScore") {

      sub_obj <- Seurat::AddModuleScore(
        object = sub_obj,
        features = list(valid_genes),
        name = "ModuleScore",
        assay = assay
      )

      score_col_name <- "ModuleScore1"

      sc_data <- Seurat::FetchData(
        sub_obj,
        vars = c(score_col_name, patient_col, subtype_col)
      )

      sc_data <- sc_data %>%
        tibble::rownames_to_column("cell_id") %>%
        dplyr::rename(
          score = dplyr::all_of(score_col_name),
          patient_id = dplyr::all_of(patient_col),
          celltype_sub = dplyr::all_of(subtype_col)
        )
    }

    if (score_method %in% c("mean_z", "mean_expr")) {

      expr <- get_assay_data_safe(
        obj = sub_obj,
        assay = assay,
        slot = slot
      )

      # ----------------------------------------------------------
      # mean_z:
      # gene-wise z-score is calculated within this program's
      # target cells, matching plot_geneset_pseudobulk_survival_v3
      # ----------------------------------------------------------
      if (score_method == "mean_z") {

        gene_mean <- Matrix::rowMeans(expr)
        gene_sq_mean <- Matrix::rowMeans(expr ^ 2)
        gene_sd <- sqrt(gene_sq_mean - gene_mean ^ 2)
        gene_sd[gene_sd == 0 | is.na(gene_sd)] <- 1

        expr_use <- expr
        expr_use <- expr_use - gene_mean
        expr_use <- expr_use / gene_sd

      } else {

        expr_use <- expr
      }

      cell_score <- Matrix::colMeans(
        expr_use[valid_genes, , drop = FALSE]
      )

      sc_data <- Seurat::FetchData(
        sub_obj,
        vars = c(patient_col, subtype_col)
      )

      sc_data$score <- as.numeric(cell_score[rownames(sc_data)])

      sc_data <- sc_data %>%
        tibble::rownames_to_column("cell_id") %>%
        dplyr::rename(
          patient_id = dplyr::all_of(patient_col),
          celltype_sub = dplyr::all_of(subtype_col)
        ) %>%
        dplyr::select(
          cell_id,
          patient_id,
          celltype_sub,
          score
        )
    }

    sc_data <- sc_data %>%
      dplyr::mutate(
        program = program_name,
        n_genes = length(valid_genes),
        n_target_cells_total = length(target_cells),
        target_subtypes = paste(target_subtypes, collapse = ";"),
        valid_genes = paste(valid_genes, collapse = ";"),
        score_method = score_method
      )

    score_records[[program_name]] <- sc_data
  }

  cell_scores <- dplyr::bind_rows(score_records)

  if (nrow(cell_scores) == 0) {
    stop("No valid program scores were calculated.")
  }

  # ============================================================
  # 4. Aggregate to patient level
  #    与 plot_geneset_pseudobulk_survival_v3 一致：
  #    patient_score = mean(Cell_Score)
  # ============================================================
  patient_scores_long <- cell_scores %>%
    dplyr::group_by(patient_id, program) %>%
    dplyr::summarise(
      patient_score = mean(score, na.rm = TRUE),
      n_cells = dplyr::n(),
      n_genes = max(n_genes, na.rm = TRUE),
      target_subtypes = dplyr::first(target_subtypes),
      valid_genes = dplyr::first(valid_genes),
      score_method = dplyr::first(score_method),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_cells >= min_cells_per_patient_program)

  if (nrow(patient_scores_long) == 0) {
    stop(
      "No patient-level program scores remained after filtering by ",
      "min_cells_per_patient_program = ",
      min_cells_per_patient_program
    )
  }

  # ============================================================
  # 5. Wide patient × program matrix
  # ============================================================
  patient_program_matrix <- patient_scores_long %>%
    dplyr::select(patient_id, program, patient_score) %>%
    tidyr::pivot_wider(
      names_from = program,
      values_from = patient_score
    ) %>%
    tibble::column_to_rownames("patient_id") %>%
    as.matrix()

  patient_program_matrix_t <- t(patient_program_matrix)

  # ============================================================
  # 6. Summary table
  # ============================================================
  program_summary <- patient_scores_long %>%
    dplyr::group_by(program) %>%
    dplyr::summarise(
      n_patients = dplyr::n(),
      median_cells_per_patient = median(n_cells, na.rm = TRUE),
      min_cells_per_patient = min(n_cells, na.rm = TRUE),
      max_cells_per_patient = max(n_cells, na.rm = TRUE),
      n_genes = max(n_genes, na.rm = TRUE),
      target_subtypes = dplyr::first(target_subtypes),
      valid_genes = dplyr::first(valid_genes),
      score_method = dplyr::first(score_method),
      .groups = "drop"
    )

  message(">>> Finished.")
  message(">>> Valid cell-level scores: ", nrow(cell_scores))
  message(">>> Valid patient-program scores: ", nrow(patient_scores_long))
  message(">>> Programs retained: ", ncol(patient_program_matrix))
  message(">>> Patients retained: ", nrow(patient_program_matrix))

  return(list(
    cell_scores = cell_scores,
    patient_scores_long = patient_scores_long,
    patient_program_matrix = patient_program_matrix,
    patient_program_matrix_t = patient_program_matrix_t,
    program_summary = program_summary,
    score_method = score_method,
    baseline_value = baseline_value,
    assay = assay,
    slot = slot
  ))
}

# ===== batch6.ipynb cell 980 (code) =====
score_res <- calculate_targeted_patient_scores(
  seu = seurat_merge,
  program_list = immune_programs_targeted,
  patient_col = "patient_id",
  timepoint_col = "sample_type_rn",
  subtype_col = "cell_type3",
  baseline_value = "D7",
  assay = "RNA",
  slot = "data",
  min_cells_per_patient_program = 5,
  min_genes_present = 3,
  score_method = "AddModuleScore"
)

patient_score_mat <- score_res$patient_program_matrix
program_by_patient_mat <- score_res$patient_program_matrix_t

# ===== batch6.ipynb cell 982 (code) =====
calculate_targeted_patient_delta_scores <- function(
  seu,
  program_list,
  patient_col = "patient_id",
  timepoint_col = "sample_type_rn",
  subtype_col = "cell_type3",

  baseline_value = "D0",
  followup_value = "D7",

  assay = NULL,
  slot = "data",

  min_cells_per_patient_program = 10,
  min_genes_present = 3,

  score_method = c("AddModuleScore", "mean_z", "mean_expr"),

  addmodulescore_ctrl = 100
) {

  score_method <- match.arg(score_method)

  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(dplyr)
    library(tidyr)
    library(tibble)
  })

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seu)
  }

  # ============================================================
  # 0. helper: compatible with Seurat v4/v5
  # ============================================================
  get_assay_data_safe <- function(obj, assay, slot) {
    tryCatch(
      {
        Seurat::GetAssayData(obj, assay = assay, slot = slot)
      },
      error = function(e) {
        Seurat::GetAssayData(obj, assay = assay, layer = slot)
      }
    )
  }

  # ============================================================
  # 1. Check metadata
  # ============================================================
  meta <- seu@meta.data

  required_cols <- c(patient_col, timepoint_col, subtype_col)
  missing_cols <- setdiff(required_cols, colnames(meta))

  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
  }

  timepoint_values <- c(baseline_value, followup_value)

  # ============================================================
  # 2. Subset D0 + D7 cells
  # ============================================================
  selected_cells <- rownames(meta)[
    !is.na(meta[[timepoint_col]]) &
      meta[[timepoint_col]] %in% timepoint_values
  ]

  if (length(selected_cells) == 0) {
    stop(
      "No cells found for baseline_value/followup_value = ",
      paste(timepoint_values, collapse = ", ")
    )
  }

  seu_tp <- subset(seu, cells = selected_cells)
  meta_tp <- seu_tp@meta.data

  message(">>> Selected timepoints: ", baseline_value, " + ", followup_value)
  message(">>> Selected cells: ", ncol(seu_tp))
  message(">>> Selected patients: ", length(unique(meta_tp[[patient_col]])))
  message(">>> Programs to score: ", length(program_list))
  message(">>> Delta direction: ", followup_value, " - ", baseline_value)

  # ============================================================
  # 3. Calculate program scores
  #    每个 program 单独 subset target cells 后计算 score
  #    注意：
  #    对 mean_z 来说，z-score 在 D0+D7 的 target cells 合并后计算，
  #    这样 D0 和 D7 的 score 在同一尺度上，delta 才有意义。
  # ============================================================
  score_records <- list()

  for (program_name in names(program_list)) {

    message(">>> Scoring program: ", program_name)

    target_subtypes <- program_list[[program_name]]$subtypes
    genes <- unique(program_list[[program_name]]$genes)

    target_subtypes <- as.character(target_subtypes)
    genes <- as.character(genes)

    # ------------------------------------------------------------
    # 3.1 target cells: D0 + D7 + compatible cell subtypes
    # ------------------------------------------------------------
    target_cells <- rownames(meta_tp)[
      !is.na(meta_tp[[subtype_col]]) &
        meta_tp[[subtype_col]] %in% target_subtypes
    ]

    if (length(target_cells) == 0) {
      message("Skip ", program_name, ": no target cells found.")
      next
    }

    sub_obj <- subset(seu_tp, cells = target_cells)

    valid_genes <- intersect(genes, rownames(sub_obj))
    missing_genes <- setdiff(genes, valid_genes)

    if (length(valid_genes) < min_genes_present) {
      message(
        "Skip ", program_name,
        ": only ", length(valid_genes),
        " genes present. Required >= ", min_genes_present
      )
      next
    }

    if (length(missing_genes) > 0) {
      message(
        "Program ", program_name,
        ": missing genes removed: ",
        paste(missing_genes, collapse = ", ")
      )
    }

    # ------------------------------------------------------------
    # 3.2 Score calculation
    # ------------------------------------------------------------
    if (score_method == "AddModuleScore") {

      sub_obj <- Seurat::AddModuleScore(
        object = sub_obj,
        features = list(valid_genes),
        name = "ModuleScore",
        assay = assay,
        ctrl = addmodulescore_ctrl
      )

      score_col_name <- "ModuleScore1"

      sc_data <- Seurat::FetchData(
        sub_obj,
        vars = c(
          score_col_name,
          patient_col,
          timepoint_col,
          subtype_col
        )
      )

      sc_data <- sc_data %>%
        tibble::rownames_to_column("cell_id") %>%
        dplyr::rename(
          score = dplyr::all_of(score_col_name),
          patient_id = dplyr::all_of(patient_col),
          timepoint = dplyr::all_of(timepoint_col),
          celltype_sub = dplyr::all_of(subtype_col)
        )
    }

    if (score_method %in% c("mean_z", "mean_expr")) {

      expr <- get_assay_data_safe(
        obj = sub_obj,
        assay = assay,
        slot = slot
      )

      expr_valid <- expr[valid_genes, , drop = FALSE]

      # ----------------------------------------------------------
      # mean_z:
      # gene-wise z-score calculated within this program's
      # D0 + D7 target cells, so D7-D0 is comparable.
      # ----------------------------------------------------------
      if (score_method == "mean_z") {

        gene_mean <- Matrix::rowMeans(expr_valid)
        gene_sq_mean <- Matrix::rowMeans(expr_valid ^ 2)
        gene_sd <- sqrt(gene_sq_mean - gene_mean ^ 2)
        gene_sd[gene_sd == 0 | is.na(gene_sd)] <- 1

        expr_use <- as.matrix(expr_valid)
        expr_use <- sweep(expr_use, 1, gene_mean, "-")
        expr_use <- sweep(expr_use, 1, gene_sd, "/")

        cell_score <- colMeans(expr_use, na.rm = TRUE)

      } else {

        cell_score <- Matrix::colMeans(expr_valid)
      }

      sc_data <- Seurat::FetchData(
        sub_obj,
        vars = c(
          patient_col,
          timepoint_col,
          subtype_col
        )
      )

      sc_data$score <- as.numeric(cell_score[rownames(sc_data)])

      sc_data <- sc_data %>%
        tibble::rownames_to_column("cell_id") %>%
        dplyr::rename(
          patient_id = dplyr::all_of(patient_col),
          timepoint = dplyr::all_of(timepoint_col),
          celltype_sub = dplyr::all_of(subtype_col)
        ) %>%
        dplyr::select(
          cell_id,
          patient_id,
          timepoint,
          celltype_sub,
          score
        )
    }

    sc_data <- sc_data %>%
      dplyr::mutate(
        patient_id = as.character(patient_id),
        timepoint = as.character(timepoint),
        program = program_name,
        n_genes = length(valid_genes),
        n_target_cells_total = length(target_cells),
        target_subtypes = paste(target_subtypes, collapse = ";"),
        valid_genes = paste(valid_genes, collapse = ";"),
        score_method = score_method
      ) %>%
      dplyr::filter(
        timepoint %in% timepoint_values
      )

    score_records[[program_name]] <- sc_data
  }

  cell_scores <- dplyr::bind_rows(score_records)

  if (nrow(cell_scores) == 0) {
    stop("No valid program scores were calculated.")
  }

  # ============================================================
  # 4. Aggregate to patient × timepoint × program level
  #    patient_score = mean(Cell_Score)
  # ============================================================
  patient_scores_long <- cell_scores %>%
    dplyr::group_by(patient_id, timepoint, program) %>%
    dplyr::summarise(
      patient_score = mean(score, na.rm = TRUE),
      n_cells = dplyr::n(),
      n_genes = max(n_genes, na.rm = TRUE),
      target_subtypes = dplyr::first(target_subtypes),
      valid_genes = dplyr::first(valid_genes),
      score_method = dplyr::first(score_method),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_cells >= min_cells_per_patient_program)

  if (nrow(patient_scores_long) == 0) {
    stop(
      "No patient-level program scores remained after filtering by ",
      "min_cells_per_patient_program = ",
      min_cells_per_patient_program
    )
  }

  message(">>> Patient-program-timepoint scores: ", nrow(patient_scores_long))

  # ============================================================
  # 5. Calculate delta score: D7 - D0
  # ============================================================
  baseline_scores <- patient_scores_long %>%
    dplyr::filter(timepoint == baseline_value) %>%
    dplyr::select(
      patient_id,
      program,
      baseline_score = patient_score,
      baseline_n_cells = n_cells
    )

  followup_scores <- patient_scores_long %>%
    dplyr::filter(timepoint == followup_value) %>%
    dplyr::select(
      patient_id,
      program,
      followup_score = patient_score,
      followup_n_cells = n_cells,
      n_genes,
      target_subtypes,
      valid_genes,
      score_method
    )

  delta_scores_long <- followup_scores %>%
    dplyr::inner_join(
      baseline_scores,
      by = c("patient_id", "program")
    ) %>%
    dplyr::mutate(
      delta_score = followup_score - baseline_score,
      delta_label = paste0(followup_value, "-", baseline_value)
    ) %>%
    dplyr::select(
      patient_id,
      program,
      delta_score,
      followup_score,
      baseline_score,
      followup_n_cells,
      baseline_n_cells,
      n_genes,
      target_subtypes,
      valid_genes,
      score_method,
      delta_label
    )

  if (nrow(delta_scores_long) == 0) {
    stop(
      "No paired patient-program scores found for ",
      followup_value,
      " - ",
      baseline_value,
      ". Check whether the same patients have both timepoints."
    )
  }

  # ============================================================
  # 6. Wide patient × program delta matrix
  # ============================================================
  patient_delta_matrix <- delta_scores_long %>%
    dplyr::select(patient_id, program, delta_score) %>%
    tidyr::pivot_wider(
      names_from = program,
      values_from = delta_score
    ) %>%
    tibble::column_to_rownames("patient_id") %>%
    as.matrix()

  patient_delta_matrix_t <- t(patient_delta_matrix)

  # 带 patient_id 列的 data.frame，方便后续和 clinical_df merge
  patient_delta_df <- delta_scores_long %>%
    dplyr::select(patient_id, program, delta_score) %>%
    tidyr::pivot_wider(
      names_from = program,
      values_from = delta_score
    )

  # ============================================================
  # 7. Optional baseline/followup matrices from paired samples
  # ============================================================
  patient_baseline_matrix <- delta_scores_long %>%
    dplyr::select(patient_id, program, baseline_score) %>%
    tidyr::pivot_wider(
      names_from = program,
      values_from = baseline_score
    ) %>%
    tibble::column_to_rownames("patient_id") %>%
    as.matrix()

  patient_followup_matrix <- delta_scores_long %>%
    dplyr::select(patient_id, program, followup_score) %>%
    tidyr::pivot_wider(
      names_from = program,
      values_from = followup_score
    ) %>%
    tibble::column_to_rownames("patient_id") %>%
    as.matrix()

  # ============================================================
  # 8. Summary table
  # ============================================================
  program_summary <- delta_scores_long %>%
    dplyr::group_by(program) %>%
    dplyr::summarise(
      n_paired_patients = dplyr::n(),
      median_baseline_cells_per_patient = median(baseline_n_cells, na.rm = TRUE),
      median_followup_cells_per_patient = median(followup_n_cells, na.rm = TRUE),
      min_baseline_cells_per_patient = min(baseline_n_cells, na.rm = TRUE),
      min_followup_cells_per_patient = min(followup_n_cells, na.rm = TRUE),
      max_baseline_cells_per_patient = max(baseline_n_cells, na.rm = TRUE),
      max_followup_cells_per_patient = max(followup_n_cells, na.rm = TRUE),
      n_genes = max(n_genes, na.rm = TRUE),
      target_subtypes = dplyr::first(target_subtypes),
      valid_genes = dplyr::first(valid_genes),
      score_method = dplyr::first(score_method),
      delta_label = dplyr::first(delta_label),
      .groups = "drop"
    )

  timepoint_summary <- patient_scores_long %>%
    dplyr::group_by(program, timepoint) %>%
    dplyr::summarise(
      n_patients = dplyr::n(),
      median_cells_per_patient = median(n_cells, na.rm = TRUE),
      min_cells_per_patient = min(n_cells, na.rm = TRUE),
      max_cells_per_patient = max(n_cells, na.rm = TRUE),
      .groups = "drop"
    )

  message(">>> Finished.")
  message(">>> Valid cell-level scores: ", nrow(cell_scores))
  message(">>> Valid patient-timepoint-program scores: ", nrow(patient_scores_long))
  message(">>> Valid paired delta scores: ", nrow(delta_scores_long))
  message(">>> Programs retained in delta matrix: ", ncol(patient_delta_matrix))
  message(">>> Patients retained in delta matrix: ", nrow(patient_delta_matrix))
  message(">>> Delta matrix: ", followup_value, " - ", baseline_value)

  return(list(
    cell_scores = cell_scores,

    patient_scores_long = patient_scores_long,
    delta_scores_long = delta_scores_long,

    patient_delta_matrix = patient_delta_matrix,
    patient_delta_matrix_t = patient_delta_matrix_t,
    patient_delta_df = patient_delta_df,

    patient_baseline_matrix = patient_baseline_matrix,
    patient_followup_matrix = patient_followup_matrix,

    program_summary = program_summary,
    timepoint_summary = timepoint_summary,

    score_method = score_method,
    baseline_value = baseline_value,
    followup_value = followup_value,
    delta_label = paste0(followup_value, "-", baseline_value),

    assay = assay,
    slot = slot
  ))
}

# ===== batch6.ipynb cell 983 (code) =====
delta_res <- calculate_targeted_patient_delta_scores(
  seu = seurat_merge,
  program_list = immune_programs_targeted,
  patient_col = "patient_id",
  timepoint_col = "sample_type_rn",
  subtype_col = "cell_type3",
  baseline_value = "D0",
  followup_value = "D7",
  score_method = "mean_z",
  min_cells_per_patient_program = 10,
  min_genes_present = 3
)

# AUTO-GENERATED FROM batch6.ipynb.
# Keep this file synchronized with records/function_provenance.csv.

# ===== Original batch6.ipynb cell 883 =====
plot_clusters_group <- function(meta_data, group_var = "new_group", color_var = "therapy_effect_adj", 
                          cell_cluster = "seurat_clusters",
                          levels = c("PR_BL", "PR_CT", "PR_DR", "SD_BL", "SD_CT", "SD_DR", "PD_BL", "PD_CT"),
                          comparision = list(c("PR_BL", "SD_BL"), c("PR_BL", "PR_DR")), 
                          methods = "part", line = FALSE) {
  
  require(ggplot2)
  require(dplyr)
  require(ggpubr) # 用于 stat_compare_means
   
  # 1. 计算频率和百分比
  if(methods == "part"){
    results <- meta_data %>% filter(!(cell_type3 %in% exclude_cell_type), !is.na(cell_type3)) %>% 
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
      dplyr::select(orig.ident, patient_id, !!sym(cell_cluster), n_sample_cells, cluster_freq, !!sym(color_var), !!sym(group_var)) %>% 
      distinct() %>% 
      mutate(percent = cluster_freq / n_sample_cells * 100)
  }

  # 2. 调整 group_var 的因子水平
  results[[group_var]] <- factor(results[[group_var]], levels = levels)

  # 3. 绘图准备
  cell_cluster_formula <- as.formula(paste0("~ ", cell_cluster))
  y_labs = ifelse(methods == "part", "Percentage in Major Celltype", "Percentage in samples")

  # 4. 构建基础图形
  p <- ggplot(results, aes(x = !!sym(group_var), y = percent)) + 
    geom_boxplot(outlier.shape = NA) # 如果有连线，通常隐藏 boxplot 的离群点防止重叠

  # --- 修改部分：添加配对连线 ---
  if(line) {
    # 使用 geom_line，group 映射为 patient_id
    # alpha 设置透明度防止遮挡，color 可以固定也可以根据属性映射
    p <- p + geom_line(aes(group = patient_id), color = "grey70", alpha = 0.8, linetype = "solid")
  }
  # -------------------------

  p <- p + 
    geom_jitter(aes(color = !!sym(color_var)), width = 0.1) + # 缩小 jitter 宽度，方便连线对齐
    facet_wrap(cell_cluster_formula, scales = "free_y") + 
    theme_classic() +
    theme(
      strip.text = element_text(size = 14), 
      axis.text = element_text(size = 14),
      axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      plot.title = element_text(size = 14)
    ) + stat_compare_means(comparisons = comparision, method = "t.test") + 
    labs(y = y_labs)
  return(p)
}

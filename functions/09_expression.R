# AUTO-GENERATED FROM batch6.ipynb.
# Keep this file synchronized with records/function_provenance.csv.

# ===== Original batch6.ipynb cell 281 =====
library(Seurat)
library(dplyr)

#' 计算样本水平的差异表达基因 (支持指定两组比较)
#' 
#' @param seurat_obj Seurat 对象
#' @param cell_type_name 指定的细胞类型名称
#' @param group_col 包含分组信息的列名
#' @param comparisons 字符向量，长度为2，指定要比较的组 (例如 c("R", "NR"))
#' @param cell_type_col 细胞类型所在的列名
#' @param sample_col 样本 ID 列 (用于 Pseudo-bulk 聚合)
#' @param assay 使用的 Assay
#' 
#' @return 返回差异分析结果表格
get_sample_level_degs_specific <- function(seurat_obj, 
                                           cell_type_name, 
                                           group_col, 
                                           comparisons, 
                                           cell_type_col = "cell_type3", 
                                           sample_col = "patient_id",
                                           assay = "RNA") {
  
  # 1. 检查参数有效性
  if (length(comparisons) != 2) {
    stop("comparisons 参数必须是长度为 2 的字符向量，例如 c('GroupA', 'GroupB')")
  }
  
  # 2. 提取指定细胞类型且属于指定分组的子集
  message("正在提取子集: ", cell_type_name, " (", comparisons[1], " vs ", comparisons[2], ")")
  
  # 逻辑：先选细胞类型，再选分组
  sub_indices <- which(seurat_obj@meta.data[[cell_type_col]] == cell_type_name & 
                       seurat_obj@meta.data[[group_col]] %in% comparisons)
  
  if (length(sub_indices) == 0) stop("未找到匹配的细胞或分组，请检查名称是否正确。")
  
  sub_obj <- seurat_obj[, sub_indices]
  
  # 3. 聚合样本水平表达量 (Pseudo-bulk)
  # 注意：这里按 sample_col 聚合，确保每个点代表一个病人
  pb_matrix <- AverageExpression(sub_obj, 
                                 group.by = sample_col, 
                                 assays = assay, 
                                 slot = "data")[[assay]]
  
  # 4. 建立样本与组别的映射关系
  sample_info <- sub_obj@meta.data %>%
    select(!!sym(sample_col), !!sym(group_col)) %>%
    distinct() %>%
    filter(!!sym(group_col) %in% comparisons)
  
  # 确保矩阵列名与映射表一致
  common_samples <- intersect(colnames(pb_matrix), sample_info[[sample_col]])
  pb_matrix <- pb_matrix[, common_samples]
  sample_info <- sample_info[match(common_samples, sample_info[[sample_col]]), ]
  
  # 分组样本列表
  g1_samples <- sample_info[[sample_col]][sample_info[[group_col]] == comparisons[1]]
  g2_samples <- sample_info[[sample_col]][sample_info[[group_col]] == comparisons[2]]
  
  # 检查样本量
  if (length(g1_samples) < 2 || length(g2_samples) < 2) {
    warning("其中一组样本量少于 2，统计效力可能不足或无法计算 P 值。")
  }
  
  # 5. Wilcoxon 检验核心逻辑
  message("正在进行样本水平 Wilcoxon 检验...")
  
  results <- apply(pb_matrix, 1, function(gene_expr) {
    vec1 <- gene_expr[g1_samples]
    vec2 <- gene_expr[g2_samples]
    
    # 过滤极低表达基因
    if (mean(vec1) < 0.01 && mean(vec2) < 0.01) return(c(log2FC = 0, p_val = 1))
    
    # 计算 log2FC: log2( (mean_g1 + eps) / (mean_g2 + eps) )
    # 正值代表在 comparisons[1] 中高表达
    log2fc <- log2((mean(vec1) + 1e-6) / (mean(vec2) + 1e-6))
    
    # 统计检验 (捕获可能因样本太少导致的报错)
    p_val <- tryCatch(wilcox.test(vec1, vec2)$p.value, error = function(e) 1)
    
    return(c(log2FC = log2fc, p_val = p_val))
  })
  
  # 6. 整理结果
  deg_df <- as.data.frame(t(results)) %>%
    mutate(gene = rownames(pb_matrix),
           group_high = comparisons[1],
           group_low = comparisons[2]) %>%
    mutate(p_val_adj = p.adjust(p_val, method = "BH")) %>%
    select(gene, log2FC, p_val, p_val_adj, group_high, group_low) %>%
    arrange(p_val)
  
  message("计算完成！")
  return(deg_df)
}


# Fixed gene sets used for Figure 3d; copied from the confirmed source block.
fig3d_cd8_features <- list(
  CytoxScore = c("PRF1", "IFNG", "GNLY", "NKG7", "GZMB", "GZMA", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7"),
  ActivationScore = c("CD69", "CCR7", "CD27", "BTLA", "CD40LG", "IL2RA", "CD3E", "CD47", "EOMES", "GNLY", "GZMA", "GZMB", "PRF1", "IFNG", "CD8A", "CD8B", "FASLG", "LAMP1", "LAG3", "CTLA4", "HLA-DRA", "TNFRSF4", "ICOS", "TNFRSF9", "TNFRSF18")
)

make_fig3d_dotplot <- function(cd8_obj, res_sig, genes = fig3d_cd8_features$ActivationScore) {
  dot <- Seurat::DotPlot(cd8_obj, features = genes, group.by = "sample_type_rn")
  plot_data <- dot$data
  sig_d7 <- res_sig |>
    dplyr::filter(Comparison == "D7_vs_D0", p_val < 0.05) |>
    dplyr::pull(gene) |> unique()
  sig_post <- res_sig |>
    dplyr::filter(Comparison == "postICI_vs_D0", p_val < 0.05) |>
    dplyr::pull(gene) |> unique()
  plot_data <- plot_data |>
    dplyr::mutate(sig_mark = dplyr::case_when(
      id == "D7" & features.plot %in% sig_d7 ~ "black",
      id == "postICI" & features.plot %in% sig_post ~ "red",
      TRUE ~ "none"
    ))
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = features.plot, y = id)) +
    ggplot2::geom_point(ggplot2::aes(size = pct.exp, color = avg.exp.scaled)) +
    ggplot2::geom_point(data = dplyr::filter(plot_data, sig_mark == "black"), ggplot2::aes(size = pct.exp), shape = 21, color = "black", fill = NA, stroke = 1.2) +
    ggplot2::geom_point(data = dplyr::filter(plot_data, sig_mark == "red"), ggplot2::aes(size = pct.exp), shape = 21, color = "red", fill = NA, stroke = 1.2) +
    ggplot2::scale_color_gradientn(colors = grDevices::colorRampPalette(c("#4575B4", "white", "#D73027"))(100)) +
    ggplot2::scale_size_continuous(range = c(1, 8)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "Differential Gene Expression Across Timepoints", subtitle = "Black: D7 vs D0; red: post-ICI vs D0", x = "Genes", y = "Timepoints", color = "Average expression", size = "Percent expressed") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "bottom")
  list(plot = p, data = plot_data, significant_genes = res_sig)
}

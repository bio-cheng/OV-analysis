# AUTO-GENERATED FROM batch6.ipynb.
# Keep this file synchronized with records/function_provenance.csv.

# ===== Original batch6.ipynb cell 405 =====
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


# Standard manuscript-style volcano wrapper used by Fig. 2c and Fig. 3a.
plot_cell_frequency_volcano <- function(diff_results, title) {
  required <- c("cell_type3", "logFC", "p_value")
  absent <- setdiff(required, colnames(diff_results))
  if (length(absent)) stop("Missing frequency-result fields: ", paste(absent, collapse = ", "))
  df <- diff_results |>
    dplyr::mutate(significance = dplyr::case_when(
      p_value < 0.05 & logFC > 0 ~ "p<0.05 & logFC>0",
      p_value < 0.05 & logFC < 0 ~ "p<0.05 & logFC<0",
      TRUE ~ "Not sig"
    ))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = logFC, y = -log10(p_value), color = significance)) +
    ggplot2::geom_point(ggplot2::aes(size = abs(logFC)), alpha = 0.8) +
    ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
    ggrepel::geom_text_repel(ggplot2::aes(label = cell_type3), size = 4, max.overlaps = 30) +
    ggplot2::scale_color_manual(values = c("p<0.05 & logFC>0" = "red", "p<0.05 & logFC<0" = "blue", "Not sig" = "gray")) +
    ggplot2::labs(title = title, x = "Log2 fold change", y = "-Log10(P)", color = "Significance", size = "|logFC|") +
    ggplot2::theme_bw(base_size = 12) + ggplot2::theme(legend.position = "bottom", plot.title = ggplot2::element_text(hjust = 0.5))
  list(plot = p, data = df)
}

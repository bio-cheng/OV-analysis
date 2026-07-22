# AUTO-GENERATED PROVENANCE FILE -- DO NOT EDIT
# Source: ../../batch6.ipynb
# Each block was copied verbatim; it may depend on interactive
# objects created outside this extraction. Use the source index and
# README before running individual blocks.


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

# ===== batch6.ipynb cell 480 (code) =====
diff_results_bl <- calculate_cluster_diff(seurat_merge, 
                                  group_var = "group_rn",
                                  cell_cluster = "cell_type3",
                                  group1 = "PFS>6_D0", 
                                  group2 = "PFS<6_D0", methods = "part", 
                                  paired = FALSE
)

# ===== batch6.ipynb cell 482 (code) =====
# 加载必要包
library(ggplot2)
library(dplyr)
library(ggrepel)

# 添加显著性标记列
df <- diff_results_bl %>%
  mutate(
    significance = case_when(
      p_value < 0.05 & logFC > 0 ~ "p<0.05 & logFC>0",
      p_value < 0.05 & logFC < -0 ~ "p<0.05 & logFC < 0",
#      abs(logFC) > 1 ~ "|logFC|>1",
      TRUE ~ "Not sig"
    )
  )


options(repr.plot.width = 6, repr.plot.height = 6)
# 创建火山图
ggplot(df, aes(x = logFC, y = -log10(p_value), color = significance)) +
  geom_point(aes(color = significance, size = abs(logFC)), alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  geom_text_repel(
    aes(label = cell_type3),
    size = 5,
    box.padding = 0.5,
    point.padding = 0.1,
    max.overlaps = 20,
    segment.color = "grey50"
  ) +
  scale_color_manual(
    values = c(
      "p<0.05 & logFC>0" = "red",
      "p<0.05 & logFC < 0" = "blue",
      "Not sig" = "gray"
    )
  ) +
  labs(
    title = "Differential Cell Freq: PFS>6_D0 vs PFS<6_D0",
    x = "Log2 Fold Change (logFC)",
    y = "-Log10(p-value)",
    color = "Significance",
    size = "|logFC|"
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
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# ===== batch6.ipynb cell 484 (code) =====
write.csv(df, "plot_v4/scatter_cell_type_3_PFS>6_D0_vs_PFS<6_D0.csv")
ggsave("plot_v4/scatter_cell_type_3_PFS>6_D0_vs_PFS<6_D0.pdf", width = 6, height = 6)

# ===== batch6.ipynb cell 787 (code) =====
get_prop_matrix <- function(seurat_obj, sample_col = "orig.ident", 
                            celltype_col = "cell_type2", scale = TRUE) {
  
  # 1. 计算频数表
    meta <- seurat_obj@meta.data #%>% filter(!cell_type3 %in% exclude_cell_type)
  count_df <- as.data.frame(table(meta[[sample_col]], 
                                  meta[[celltype_col]]))
  colnames(count_df) <- c("Sample", "CellType", "Count")
  
  # 2. 转换为宽表 (比例)
  prop_df <- count_df %>%
    group_by(Sample) %>%
    mutate(Prop = Count / sum(Count)) %>%
    dplyr::select(-Count) %>%
    pivot_wider(names_from = CellType, values_from = Prop) %>%
    column_to_rownames("Sample")
  
  prop_matrix <- as.matrix(prop_df)
  
  # 3. 标准化 (按列/细胞类型进行 Z-score)
  if (scale) {
    prop_matrix <- scale(prop_matrix)
    # 处理可能的恒定比例导致的NaN
    prop_matrix[is.nan(prop_matrix)] <- 0
  }
  
  return(prop_matrix)
}

# ===== batch6.ipynb cell 788 (code) =====
get_prop_matrix <- function(seurat_obj, sample_col = "orig.ident", 
                            celltype_col = "cell_type2", scale = TRUE) {
  
  # 1. 计算频数表
  count_df <- as.data.frame(table(seurat_obj@meta.data[[sample_col]], 
                                  seurat_obj@meta.data[[celltype_col]]))
  colnames(count_df) <- c("Sample", "CellType", "Count")
  
  # 2. 转换为宽表 (比例)
  prop_df <- count_df %>%
    group_by(Sample) %>%
    mutate(Prop = Count / sum(Count)) %>%
    dplyr::select(-Count) %>%
    pivot_wider(names_from = CellType, values_from = Prop) %>%
    column_to_rownames("Sample")
  
  prop_matrix <- as.matrix(prop_df)
  
  # 3. 标准化 (按列/细胞类型进行 Z-score)
  if (scale) {
    prop_matrix <- scale(prop_matrix)
    # 处理可能的恒定比例导致的NaN
    prop_matrix[is.nan(prop_matrix)] <- 0
  }
  
  return(prop_matrix)
}

# ===== batch6.ipynb cell 789 (code) =====
plot_prop_heatmap <- function(prop_matrix, seurat_obj = NULL, 
                              annotation_cols = NULL, 
                              ann_colors = NULL, ...) {
  
  require(pheatmap)
  require(dplyr)
  require(tibble)
  
  mat_to_plot <- t(prop_matrix)
  
  # 1. 处理注释信息
  annotation_info <- list()
  if (!is.null(seurat_obj) && !is.null(annotation_cols)) {
    meta <- seurat_obj@meta.data %>%
      dplyr::select(orig.ident, all_of(annotation_cols)) %>%
      distinct() %>%
      tibble::remove_rownames() %>% 
      tibble::column_to_rownames("orig.ident")
    
    meta <- meta[colnames(mat_to_plot), , drop = FALSE]
    annotation_info$annotation_col <- meta
    
    # 2. 处理自定义颜色
    if (!is.null(ann_colors)) {
      annotation_info$annotation_colors <- ann_colors
    }
  }

  # 3. 基础绘图参数
  base_params <- list(
    mat = mat_to_plot,
    color = colorRampPalette(c("#4575b4", "white", "#d73027"))(100),
    border_color = NA,
    clustering_method = "ward.D2",
    main = "Cell Type Proportion Heatmap",
    fontsize = 12
  )

  # 4. 合并参数并绘图
  final_params <- modifyList(base_params, list(...))
  final_params <- modifyList(final_params, annotation_info)
  
  do.call(pheatmap, final_params)
}

# ===== batch6.ipynb cell 790 (code) =====
mat <- get_prop_matrix(seurat_merge, sample_col = "orig.ident", celltype_col = "cell_type3")

# ===== batch6.ipynb cell 791 (code) =====
options(repr.plot.width = 15, repr.plot.height = 10)
plot_prop_heatmap(mat)

# ===== batch6.ipynb cell 792 (code) =====


# ===== batch6.ipynb cell 793 (code) =====
library(tidyverse)
library(ggsci)

plot_cns_clinical_matrix <- function(seurat_obj) {
  
  # 1. 提取并清理元数据
  meta <- seurat_obj@meta.data %>% 
    filter(sample_type != "02.2") %>%
    select(patient_id, sample_type, response, pfs_6m) %>%
    distinct() %>%
    mutate(timepoint = sample_type) %>%
    filter(!is.na(timepoint))
  
  # 2. 采样状态
  sampling_status <- meta %>%
    complete(patient_id, timepoint) %>%
    mutate(Status = ifelse(!is.na(sample_type), "Sampled", "Missing")) %>%
    select(patient_id, timepoint, Status) %>%
    pivot_wider(names_from = timepoint, values_from = Status)
  
  # 3. 临床信息
  clinical_info <- meta %>%
    select(patient_id, response, pfs_6m) %>%
    distinct()
  
  # 4. 合并并转为绘图长表
  plot_data <- sampling_status %>%
    left_join(clinical_info, by = "patient_id") %>%
    arrange(patient_id, desc(response)) %>%
    mutate(patient_id = factor(patient_id, levels = unique(patient_id))) %>%
    pivot_longer(
      cols = -patient_id, 
      names_to = "Attribute", 
      values_to = "Value"
    ) %>%
    mutate(
      Attribute = factor(
        Attribute, 
        levels = rev(c("D0", "D7", "postICI", "response", "pfs_6m"))
      )
    )
  
  # 5. CNS 级别配色
  cns_colors <- c(
    "Sampled" = "#716db2", 
    "Missing" = "#E0E0E0",
    "R" = "#3CB371", 
    "NR" = "#CD5C5C",
    "FPS>=6 month" = "#4682B4", 
    "FPS<6 month" = "#F4A460"
  )
  
  # 6. 绘图：去掉 tile 边框和黑色网格线
  ggplot(plot_data, aes(x = patient_id, y = Attribute, fill = Value)) +
    geom_tile(color = NA) +
    scale_fill_manual(values = cns_colors, na.value = "white") +
    scale_x_discrete(position = "top") +
    labs(
      title = "Clinical Cohort Landscape", 
      subtitle = "Patient IDs (Columns) vs. Clinical Features (Rows)",
      x = NULL, 
      y = NULL
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 90, 
        hjust = -0.5, 
        vjust = -1, 
        size = 15, 
        family = "sans"
      ),
      axis.text.y = element_text(
        size = 12, 
        face = "bold", 
        color = "black"
      ),
      legend.position = "right",
      legend.text = element_text(
        size = 14, 
        face = "bold", 
        color = "black"
      ),
      plot.title = element_text(
        face = "bold", 
        size = 14
      ),
      strip.background = element_blank()
    )
}

# ===== batch6.ipynb cell 794 (code) =====
seurat_merge$sample_type <- seurat_merge$sample_type_rn

# ===== batch6.ipynb cell 795 (code) =====
unique(seurat_merge$sample_type_rn)

# ===== batch6.ipynb cell 796 (code) =====
seurat_merge$response[seurat_merge$patient_id %in% c("BZ005", "BZ039")] <- NA

# ===== batch6.ipynb cell 797 (code) =====
seurat_merge$pfs_6m <- ifelse(seurat_merge$pfs_time > 6, "FPS>=6 month", "FPS<6 month")

# ===== batch6.ipynb cell 798 (code) =====
p_cns <- plot_cns_clinical_matrix(seurat_merge)
options(repr.plot.width = 12, repr.plot.height = 3)
# 在 Jupyter/R 中展示
print(p_cns)

# ===== batch6.ipynb cell 799 (code) =====
ggsave("clinical_cohort_patient_annotation.pdf", width = 12, height = 3)

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


# ===== batch6.ipynb cell 936 (code) =====
write.csv(res$data, "plot_v4/fig2g_2.csv", row.names = F)

# ===== batch6.ipynb cell 947 (code) =====
# 假设你想分析 "D0" 阶段，"CD8+ T" 细胞的比例对预后的影响
res <- plot_stage_cell_survival(
  meta = seurat_merge@meta.data,
  target_cell = "CD8_TRM",
  target_stage = "D0",
  cell_col = "cell_type3",
  time_col = "sample_type_rn",
  patient_col = "patient_id",
  surv_time_col = "pfs_time",
  surv_status_col = "pfs_status"
)

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

# ===== batch6.ipynb cell 986 (code) =====
plot_patient_program_heatmap <- function(
  patient_score_mat,
  clinical_df = NULL,
  patient_col = "patient_id",
  order_by = c("PFS_time", "pfs_time", "response", "none"),
  scale_rows = TRUE,
  cluster_rows = TRUE,
  cluster_columns = FALSE,

  # clinical annotation 参数
  clinical_anno_cols = c(
    "pfs_group",
    "response",
    "PFS_status"
  ),
  show_clinical_annotation = TRUE,
  clinical_anno_colors = NULL,

  # function group 分组参数
  function_group_map = NULL,
  use_manual_group_map = TRUE,
  include_composite_axes = TRUE,
  include_unmapped = TRUE,
  show_function_group_bar = TRUE,

  function_group_order = c(
    "Cytotoxic / effector",
    "Antigen presentation / APC",
    "Antiviral IFN / sensing",
    "Suppressive / checkpoint",
    "Trafficking / chemotaxis",
    "Other"
  ),

  function_group_colors = NULL,

  row_label_wrap_width = 40,
  output_pdf = NULL,
  width = 13,
  height = 10,
  title = "Baseline immune programs by function group"
) {

  order_by <- match.arg(order_by)

  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(dplyr)
    library(tibble)
    library(grid)
  })

  # ============================================================
  # 0. 内置 manual function group map
  # ============================================================
  make_internal_manual_function_group_map <- function(
    candidate_features = NULL,
    include_composite_axes = TRUE,
    default_group = "Other"
  ) {

    manual_group_map <- c(
      # ================== Cytotoxic / effector ==================
      CD8_Tem_TRM__cytolytic_effector       = "Cytotoxic / effector",
      NK_C2__NK_cytotoxicity                = "Cytotoxic / effector",
      MAIT__innate_like_cytotoxicity        = "Cytotoxic / effector",
      CD8_all__activation                   = "Cytotoxic / effector",
      CD8_all__cytotoxicity                 = "Cytotoxic / effector",
      CD8_TRM__tissue_resident_memory       = "Cytotoxic / effector",
      CD4_Tem__effector_memory              = "Cytotoxic / effector",
      CD4_Th17__Th17_program                = "Cytotoxic / effector",
      CD8_Prolif__proliferation             = "Cytotoxic / effector",
      CD8_Stress__stress_response           = "Cytotoxic / effector",
      NK_Prolif__proliferation              = "Cytotoxic / effector",

      # ============== Antigen presentation / APC ================
      B_cell__BCR_antigen_presentation      = "Antigen presentation / APC",
      B_Naive__naive_B_program              = "Antigen presentation / APC",
      Plasma_cells__antibody_secreting      = "Antigen presentation / APC",
      DC_cDC1__cross_presentation           = "Antigen presentation / APC",
      DC_cDC2__MHCII_maturation             = "Antigen presentation / APC",
      Mono_all__HLA_DR_APC                  = "Antigen presentation / APC",
      Mono_Intermediate__APC_inflammatory   = "Antigen presentation / APC",
      Mono_Classical__inflammatory_monocyte = "Antigen presentation / APC",

      # ============== Antiviral IFN / sensing ===================
      DC_all__viral_sensing                 = "Antiviral IFN / sensing",
      DC_all__antiviral_ISG                 = "Antiviral IFN / sensing",
      Mono_IFN__ISG_antiviral               = "Antiviral IFN / sensing",
      Neu_IFN__ISG_antiviral                = "Antiviral IFN / sensing",
      CD4_all__Th1_antiviral_helper         = "Antiviral IFN / sensing",

      # ============== Suppressive / checkpoint ==================
      CD4_Treg__suppressive_program         = "Suppressive / checkpoint",
      CD4_LAG3__checkpoint_exhaustion       = "Suppressive / checkpoint",
      Mono_all__MDSC_like_positive          = "Suppressive / checkpoint",
      Neu_all__PMN_MDSC_suppressive         = "Suppressive / checkpoint",
      CD8_all__checkpoint_exhaustion        = "Suppressive / checkpoint",
      Neu_IL1R1__inflammatory_neutrophil    = "Suppressive / checkpoint",

      # ============== Trafficking / chemotaxis ===================
      CD8_all__chemokine_trafficking        = "Trafficking / chemotaxis",
      NK_XCL1__DC_recruiting_NK             = "Trafficking / chemotaxis",
      CD4_Naive__naive_memory               = "Trafficking / chemotaxis",
      CD8_Naive__naive_memory               = "Trafficking / chemotaxis",
      Mono_NonClassical__patrolling_surveillance =
        "Trafficking / chemotaxis",
      Neu_CXCR2__mature_neutrophil_chemotaxis =
        "Trafficking / chemotaxis"
    )

    if (isTRUE(include_composite_axes)) {
      composite_map <- c(
        Axis_Cytotoxic_Lymphocyte           = "Cytotoxic / effector",
        Axis_Antigen_Presentation           = "Antigen presentation / APC",
        Axis_Antiviral_Innate_Sensing       = "Antiviral IFN / sensing",
        Axis_Immunosuppressive              = "Suppressive / checkpoint",
        Axis_Trafficking                    = "Trafficking / chemotaxis",
        Ratio_Cytotoxic_to_Suppressive      = "Cytotoxic / effector",
        Ratio_APC_to_Myeloid_Suppressive    = "Antigen presentation / APC"
      )

      manual_group_map <- c(manual_group_map, composite_map)
    }

    if (is.null(candidate_features)) {
      return(manual_group_map)
    }

    candidate_features <- as.character(candidate_features)

    missing_features <- setdiff(
      candidate_features,
      names(manual_group_map)
    )

    if (length(missing_features) > 0) {
      missing_map <- stats::setNames(
        rep(default_group, length(missing_features)),
        missing_features
      )

      manual_group_map <- c(manual_group_map, missing_map)
    }

    manual_group_map[candidate_features]
  }

  infer_function_group <- function(x) {

    x_lower <- tolower(x)

    dplyr::case_when(
      grepl(
        "cytolytic|cytotoxic|nk_cytotoxicity|innate_like_cytotoxicity|activation|effector|proliferation|prolif|stress|th17",
        x_lower
      ) ~ "Cytotoxic / effector",

      grepl(
        "cross_presentation|mhcii|bcr_antigen_presentation|hla_dr_apc|apc_inflammatory|antigen_presentation|maturation|naive_b|antibody|plasma",
        x_lower
      ) ~ "Antigen presentation / APC",

      grepl(
        "viral_sensing|isg|ifn|antiviral|th1_antiviral",
        x_lower
      ) ~ "Antiviral IFN / sensing",

      grepl(
        "suppressive|mdsc|checkpoint|exhaustion|treg|lag3|regulatory|il1r1",
        x_lower
      ) ~ "Suppressive / checkpoint",

      grepl(
        "trafficking|chemotaxis|dc_recruiting|patrolling|surveillance|cxcr2|naive_memory",
        x_lower
      ) ~ "Trafficking / chemotaxis",

      TRUE ~ "Other"
    )
  }

  # ============================================================
  # 1. 准备 matrix：program x patient
  # ============================================================
  mat <- t(patient_score_mat)

  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]

  if (nrow(mat) == 0) {
    stop("没有可用于绘图的 immune programs。")
  }

  programs <- rownames(mat)

  # ============================================================
  # 2. function group map
  # ============================================================
  if (is.null(function_group_map)) {

    if (isTRUE(use_manual_group_map)) {
      function_group_map <- make_internal_manual_function_group_map(
        candidate_features = programs,
        include_composite_axes = include_composite_axes,
        default_group = "Other"
      )
    } else {
      function_group_map <- stats::setNames(
        infer_function_group(programs),
        programs
      )
    }

  } else {

    if (is.null(names(function_group_map))) {
      stop("function_group_map 必须是 named vector。")
    }

    missing_programs <- setdiff(programs, names(function_group_map))

    if (length(missing_programs) > 0) {
      missing_map <- stats::setNames(
        infer_function_group(missing_programs),
        missing_programs
      )

      function_group_map <- c(function_group_map, missing_map)
    }

    function_group_map <- function_group_map[programs]
  }

  row_function_group <- unname(function_group_map[programs])
  row_function_group[is.na(row_function_group)] <- "Other"

  if (!isTRUE(include_unmapped)) {
    keep_rows <- row_function_group != "Other"
    mat <- mat[keep_rows, , drop = FALSE]
    row_function_group <- row_function_group[keep_rows]
    programs <- rownames(mat)
  }

  present_groups <- unique(row_function_group)

  function_group_order_use <- c(
    intersect(function_group_order, present_groups),
    setdiff(present_groups, function_group_order)
  )

  row_function_group <- factor(
    row_function_group,
    levels = function_group_order_use
  )

  # ============================================================
  # 3. 行标准化
  # ============================================================
  if (scale_rows) {

    mat_plot <- t(scale(t(mat)))
    mat_plot[is.na(mat_plot)] <- 0

    legend_title <- "Row Z-score"

    col_fun <- circlize::colorRamp2(
      c(-2, 0, 2),
      c("#2166AC", "white", "#B2182B")
    )

  } else {

    mat_plot <- mat

    q <- stats::quantile(
      mat_plot,
      probs = c(0.05, 0.5, 0.95),
      na.rm = TRUE
    )

    col_fun <- circlize::colorRamp2(
      c(q[1], q[2], q[3]),
      c("#2166AC", "white", "#B2182B")
    )

    legend_title <- "Score"
  }

  # ============================================================
  # 4. 按 clinical 信息排序 patient columns
  # ============================================================
  if (!is.null(clinical_df) && order_by != "none") {

    clinical_sub <- clinical_df %>%
      dplyr::filter(.data[[patient_col]] %in% colnames(mat_plot)) %>%
      dplyr::distinct(.data[[patient_col]], .keep_all = TRUE)

    if (order_by %in% c("PFS_time", "pfs_time")) {

      pfs_col <- dplyr::case_when(
        "PFS_time" %in% colnames(clinical_sub) ~ "PFS_time",
        "pfs_time" %in% colnames(clinical_sub) ~ "pfs_time",
        TRUE ~ NA_character_
      )

      if (!is.na(pfs_col)) {
        clinical_sub <- clinical_sub %>%
          dplyr::mutate(
            .pfs_order = suppressWarnings(as.numeric(.data[[pfs_col]]))
          ) %>%
          dplyr::arrange(desc(.pfs_order))
      }
    }

    if (order_by == "response" && "response" %in% colnames(clinical_sub)) {
      clinical_sub <- clinical_sub %>%
        dplyr::mutate(
          response_order = dplyr::case_when(
            response %in% c("CR", "PR", "R") ~ 1,
            response %in% c("SD") ~ 2,
            response %in% c("PD", "NR") ~ 3,
            TRUE ~ 4
          )
        ) %>%
        dplyr::arrange(response_order)
    }

    patient_order <- clinical_sub[[patient_col]]

    patient_order <- unique(
      c(
        as.character(patient_order),
        setdiff(colnames(mat_plot), as.character(patient_order))
      )
    )

    mat_plot <- mat_plot[, patient_order, drop = FALSE]
  }

  # ============================================================
  # 5. column annotation：只展示指定临床信息
  # ============================================================
  top_anno <- NULL

  if (
    isTRUE(show_clinical_annotation) &&
    !is.null(clinical_df) &&
    !is.null(clinical_anno_cols) &&
    length(clinical_anno_cols) > 0
  ) {

    if (
      length(clinical_anno_cols) == 1 &&
      clinical_anno_cols %in% c("none", "None", "NONE")
    ) {
      clinical_anno_cols <- character(0)
    }

    if (length(clinical_anno_cols) > 0) {

      anno_df <- clinical_df %>%
        dplyr::filter(.data[[patient_col]] %in% colnames(mat_plot)) %>%
        dplyr::distinct(.data[[patient_col]], .keep_all = TRUE) %>%
        tibble::column_to_rownames(patient_col)

      anno_df <- anno_df[colnames(mat_plot), , drop = FALSE]

      anno_cols <- intersect(
        clinical_anno_cols,
        colnames(anno_df)
      )

      missing_anno_cols <- setdiff(
        clinical_anno_cols,
        colnames(anno_df)
      )

      if (length(missing_anno_cols) > 0) {
        message(
          ">>> 以下 clinical_anno_cols 在 clinical_df 中不存在，已跳过: ",
          paste(missing_anno_cols, collapse = ", ")
        )
      }

      if (length(anno_cols) > 0) {

        anno_use <- anno_df[, anno_cols, drop = FALSE]

        anno_use[] <- lapply(anno_use, function(x) {
          x <- as.character(x)
          x[is.na(x) | x == "" | x == "NA"] <- "Unknown"
          x
        })

        if (is.null(clinical_anno_colors)) {

          top_anno <- ComplexHeatmap::HeatmapAnnotation(
            df = anno_use,
            annotation_name_gp = grid::gpar(fontsize = 10),
            annotation_legend_param = list(
              title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
              labels_gp = grid::gpar(fontsize = 9)
            )
          )

        } else {

          clinical_anno_colors_use <- clinical_anno_colors[
            names(clinical_anno_colors) %in% anno_cols
          ]

          top_anno <- ComplexHeatmap::HeatmapAnnotation(
            df = anno_use,
            col = clinical_anno_colors_use,
            annotation_name_gp = grid::gpar(fontsize = 10),
            annotation_legend_param = list(
              title_gp = grid::gpar(fontsize = 10, fontface = "bold"),
              labels_gp = grid::gpar(fontsize = 9)
            )
          )
        }
      }
    }
  }

  # ============================================================
  # 6. row annotation：function group 色条
  # ============================================================
  if (is.null(function_group_colors)) {
    function_group_colors <- c(
      "Cytotoxic / effector"        = "#2A9D8F",
      "Antigen presentation / APC"  = "#4F7DC9",
      "Antiviral IFN / sensing"     = "#5E6ACB",
      "Suppressive / checkpoint"    = "#E76F51",
      "Trafficking / chemotaxis"    = "#7A7A7A",
      "Other"                       = "#BDBDBD"
    )
  }

  missing_group_colors <- setdiff(
    levels(row_function_group),
    names(function_group_colors)
  )

  if (length(missing_group_colors) > 0) {
    extra_cols <- stats::setNames(
      rep("#BDBDBD", length(missing_group_colors)),
      missing_group_colors
    )

    function_group_colors <- c(function_group_colors, extra_cols)
  }

  row_anno <- NULL

  if (isTRUE(show_function_group_bar)) {
    row_anno <- ComplexHeatmap::rowAnnotation(
      Function = row_function_group,
      col = list(
        Function = function_group_colors[levels(row_function_group)]
      ),
      annotation_name_gp = grid::gpar(
        fontsize = 9,
        fontface = "bold"
      ),
      annotation_legend_param = list(
        Function = list(
          title = "Function group",
          title_gp = grid::gpar(
            fontsize = 9,
            fontface = "bold"
          ),
          labels_gp = grid::gpar(fontsize = 8)
        )
      ),
      width = grid::unit(5, "mm")
    )
  }

  # ============================================================
  # 7. row labels
  # ============================================================
  wrap_label <- function(x, width = 40) {
    vapply(
      x,
      function(s) {
        paste(
          strwrap(
            gsub("__", " | ", s, fixed = TRUE),
            width = width
          ),
          collapse = "\n"
        )
      },
      character(1)
    )
  }

  row_labels <- wrap_label(
    rownames(mat_plot),
    width = row_label_wrap_width
  )

  # ============================================================
  # 8. heatmap
  # ============================================================
  ht <- ComplexHeatmap::Heatmap(
    mat_plot,
    name = legend_title,
    col = col_fun,

    top_annotation = top_anno,
    left_annotation = row_anno,

    row_split = row_function_group,
    cluster_rows = cluster_rows,
    cluster_row_slices = FALSE,
    row_gap = grid::unit(2, "mm"),

    cluster_columns = cluster_columns,

    show_row_names = TRUE,
    show_column_names = TRUE,

    row_labels = row_labels,

    row_names_gp = grid::gpar(fontsize = 11),
    column_names_gp = grid::gpar(fontsize = 11),

    row_title_gp = grid::gpar(
      fontsize = 12,
      fontface = "bold"
    ),
    row_title_rot = 0,

    column_title = title,
    column_title_gp = grid::gpar(
      fontsize = 13,
      fontface = "bold"
    ),

    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 9, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 8)
    ),

    na_col = "grey90"
  )

  # ============================================================
  # 9. 输出
  # ============================================================
  if (!is.null(output_pdf)) {

    grDevices::pdf(
      output_pdf,
      width = width,
      height = height
    )

    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = TRUE
    )

    grDevices::dev.off()

  } else {

    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = TRUE
    )
  }

  invisible(list(
    heatmap = ht,
    matrix_plot = mat_plot,
    function_group = row_function_group,
    function_group_map = function_group_map,
    clinical_annotation_columns = clinical_anno_cols
  ))
}

# ===== batch6.ipynb cell 987 (code) =====
ht_res <- plot_patient_program_heatmap(
  patient_score_mat = patient_score_mat,
  clinical_df = clinical_df,
  patient_col = "patient_id",
  order_by = "pfs_time",
  scale_rows = TRUE,
  cluster_rows = TRUE,
  cluster_columns = T,
    clinical_anno_cols = c(
    "response",
    "pfs_group"
  ),
  output_pdf = "plot_v4/baseline_patient_level_immune_program_heatmap.pdf",
  width = 15,
  height = 12,
  title = "Baseline patient-level immune programs"
)

# ===== batch6.ipynb cell 1014 (code) =====
plot_geneset_pseudobulk_survival_v3 <- function(
  seurat_obj,
  gene_set,
  target_cell,
  target_stage = "D0",
  score_name = "GeneSet_Score",
  cell_col = "cell_type3",
  time_col = "sample_type_rn",
  patient_col = "patient_id",
  surv_time_col = "pfs_time",
  surv_status_col = "pfs_statu",
  min_cells = 5,
  score_method = c("AddModuleScore", "mean_z", "mean_expr"),

  # 为了和 plot_stage_cell_survival 一致，默认使用 median_gt
  # median_gt: Mean_Score > median 为 High
  # median_ge: Mean_Score >= median 为 High
  cutoff_rule = c("median_gt", "median_ge"),

  start_at_zero = TRUE,
  break_time_by = NULL,
  x_axis_padding_frac = 0.05,

  handle_zero_time = c("keep", "epsilon", "drop"),
  epsilon_time = 0.001,

  # ============================================================
  # 字体和版式参数，与 plot_stage_cell_survival 风格对齐
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

  risk_table_show_y_text = TRUE
) {

  score_method <- match.arg(score_method)
  cutoff_rule <- match.arg(cutoff_rule)
  handle_zero_time <- match.arg(handle_zero_time)

  suppressPackageStartupMessages({
    library(Seurat)
    library(survival)
    library(survminer)
    library(dplyr)
    library(Matrix)
    library(ggplot2)
    library(tidyr)
  })

  target_cell_label <- paste(target_cell, collapse = " + ")

  message(sprintf(
    ">>> 步骤 1: subset [%s] 细胞在 [%s] 阶段...",
    target_cell_label,
    target_stage
  ))

  meta <- seurat_obj@meta.data

  required_meta_cols <- c(
    cell_col,
    time_col,
    patient_col,
    surv_time_col,
    surv_status_col
  )

  missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))

  if (length(missing_meta_cols) > 0) {
    stop(
      "meta.data 中缺少以下列：",
      paste(missing_meta_cols, collapse = ", ")
    )
  }

  cells_to_keep <- rownames(meta)[
    !is.na(meta[[time_col]]) &
      !is.na(meta[[cell_col]]) &
      meta[[time_col]] == target_stage &
      meta[[cell_col]] %in% target_cell
  ]

  if (length(cells_to_keep) < 10) {
    stop("符合条件的细胞数量极少，无法进行分析。")
  }

  sub_obj <- subset(seurat_obj, cells = cells_to_keep)

  valid_genes <- intersect(gene_set, rownames(sub_obj))
  missing_genes <- setdiff(gene_set, valid_genes)

  if (length(valid_genes) == 0) {
    stop("提供的基因集中，没有任何基因存在于 Seurat 对象中。")
  }

  if (length(missing_genes) > 0) {
    warning("以下基因未找到并被剔除: ", paste(missing_genes, collapse = ", "))
  }

  message(">>> 有效基因: ", paste(valid_genes, collapse = ", "))

  # ============================================================
  # 2. 计算基因集得分
  # ============================================================
  if (score_method == "AddModuleScore") {

    sub_obj <- Seurat::AddModuleScore(
      object = sub_obj,
      features = list(valid_genes),
      name = "ModuleScore"
    )

    score_col_name <- "ModuleScore1"

    sc_data <- Seurat::FetchData(
      sub_obj,
      vars = c(
        score_col_name,
        patient_col,
        surv_time_col,
        surv_status_col
      )
    )

    colnames(sc_data) <- c(
      "Cell_Score",
      "Patient_ID",
      "Surv_Time",
      "Raw_Status"
    )
  }

  if (score_method %in% c("mean_z", "mean_expr")) {

    expr <- Seurat::GetAssayData(
      sub_obj,
      slot = "data"
    )

    expr_valid <- expr[valid_genes, , drop = FALSE]

    if (score_method == "mean_z") {

      gene_mean <- Matrix::rowMeans(expr_valid)
      gene_sq_mean <- Matrix::rowMeans(expr_valid ^ 2)
      gene_sd <- sqrt(gene_sq_mean - gene_mean ^ 2)
      gene_sd[gene_sd == 0 | is.na(gene_sd)] <- 1

      expr_use <- as.matrix(expr_valid)
      expr_use <- sweep(expr_use, 1, gene_mean, "-")
      expr_use <- sweep(expr_use, 1, gene_sd, "/")

      cell_score <- colMeans(expr_use, na.rm = TRUE)
    }

    if (score_method == "mean_expr") {
      cell_score <- Matrix::colMeans(expr_valid)
    }

    sc_data <- Seurat::FetchData(
      sub_obj,
      vars = c(
        patient_col,
        surv_time_col,
        surv_status_col
      )
    )

    sc_data$Cell_Score <- as.numeric(cell_score[rownames(sc_data)])

    sc_data <- sc_data %>%
      dplyr::select(
        Cell_Score,
        Patient_ID = dplyr::all_of(patient_col),
        Surv_Time = dplyr::all_of(surv_time_col),
        Raw_Status = dplyr::all_of(surv_status_col)
      )
  }

  # ============================================================
  # 3. Pseudobulk 聚合到患者层面
  # ============================================================
  message(">>> 步骤 2: pseudobulk 聚合到患者层面...")

  pb_data <- sc_data %>%
    dplyr::group_by(Patient_ID) %>%
    dplyr::summarise(
      Cell_Count = dplyr::n(),
      Mean_Score = mean(Cell_Score, na.rm = TRUE),
      Surv_Time = dplyr::first(Surv_Time),
      Raw_Status = dplyr::first(Raw_Status),
      n_surv_time_values = dplyr::n_distinct(Surv_Time),
      n_surv_status_values = dplyr::n_distinct(Raw_Status),
      .groups = "drop"
    ) %>%
    dplyr::filter(Cell_Count >= min_cells)

  inconsistent_info <- pb_data %>%
    dplyr::filter(
      n_surv_time_values > 1 |
        n_surv_status_values > 1
    )

  if (nrow(inconsistent_info) > 0) {
    warning(
      "检测到部分 patient 在该阶段/细胞亚群中存在多个生存时间或状态值；",
      "当前函数使用 first()。"
    )
  }

  pb_data <- pb_data %>%
    dplyr::select(
      Patient_ID,
      Cell_Count,
      Mean_Score,
      Surv_Time,
      Raw_Status
    )

  if (nrow(pb_data) < 5) {
    stop("细胞数达标的患者样本量不足，无法进行生存分析。")
  }

  # ============================================================
  # 4. 清洗生存数据
  # ============================================================
  message(">>> 步骤 3: 清洗生存信息...")

  surv_df <- pb_data %>%
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
      !is.na(Mean_Score),
      is.finite(Mean_Score)
    ) %>%
    dplyr::filter(
      Surv_Status %in% c(0, 1)
    ) %>%
    dplyr::select(-Raw_Status_chr)

  if (nrow(surv_df) < 5) {
    stop("具有完整生存信息的样本量不足。")
  }

  # ============================================================
  # 5. 处理 Surv_Time <= 0
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

  # ============================================================
  # 6. 分组逻辑
  #    默认 median_gt，与 plot_stage_cell_survival 一致：
  #    Mean_Score > median 为 High
  #    Mean_Score <= median 为 Low
  # ============================================================
  message(">>> 步骤 4: 根据基因集 pseudobulk score 中位数分组...")

  threshold <- median(surv_df$Mean_Score, na.rm = TRUE)

  if (cutoff_rule == "median_gt") {
    surv_df$Group <- ifelse(
      surv_df$Mean_Score >= threshold,
      "High",
      "Low"
    )
  }

  if (cutoff_rule == "median_ge") {
    surv_df$Group <- ifelse(
      surv_df$Mean_Score >= threshold,
      "High",
      "Low"
    )
  }

  surv_df$Group <- factor(
    surv_df$Group,
    levels = c("Low", "High")
  )

  if (length(unique(stats::na.omit(surv_df$Group))) < 2) {
    stop("分组失败：High/Low 只有一组。")
  }

  message(">>> 纳入患者数: ", nrow(surv_df))
  message(">>> 事件数 Events: ", sum(surv_df$Surv_Status == 1, na.rm = TRUE))
  message(">>> 删失数 Censored: ", sum(surv_df$Surv_Status == 0, na.rm = TRUE))

  message(">>> 分组总人数:")
  print(table(surv_df$Group, useNA = "ifany"))

  message(">>> 各组事件数:")
  print(table(surv_df$Group, surv_df$Surv_Status, useNA = "ifany"))

  # ============================================================
  # 7. 拟合 KM
  # ============================================================
  fit <- survival::survfit(
    survival::Surv(Surv_Time, Surv_Status) ~ Group,
    data = surv_df
  )

  # ============================================================
  # 8. Log-rank P
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
  # 9. x 轴设置，与 plot_stage_cell_survival 一致
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

  common_x_scale <- function() {
    ggplot2::scale_x_continuous(
      limits = xlim_use,
      breaks = integer_breaks,
      labels = function(x) sprintf("%d", round(x)),
      expand = ggplot2::expansion(mult = c(0, 0))
    )
  }

  # ============================================================
  # 10. 绘图，与 plot_stage_cell_survival 风格对齐
  # ============================================================
  p <- survminer::ggsurvplot(
    fit,
    data = surv_df,

    pval = FALSE,
    conf.int = FALSE,

    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.height = risk_table_height,

    # risk table 中 number at risk 的字体
    fontsize = risk_table_font_size,

    risk.table.y.text = risk_table_show_y_text,
    risk.table.y.text.col = risk_table_show_y_text,

    palette = c("#377EB8", "#E41A1C"),

    title = paste0(
      "Survival Analysis: ",
      score_name,
      " in ",
      target_cell_label,
      " at ",
      target_stage
    ),

    subtitle = paste0(
      "Cutoff (median) = ",
      signif(threshold, 3),
      " | Valid genes = ",
      length(valid_genes),
      " | N = ",
      nrow(surv_df)
    ),

    legend.title = "Score",
    legend.labs = c("Low", "High"),

    xlab = "Time",
    ylab = ylab_text,

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
    tables.theme = ggplot2::theme_bw(base_size = risk_table_axis_font_size)
  )

  # ============================================================
  # 11. 主图：统一 x 轴、字体、Log-rank P 标注
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
    ) +
    ggplot2::annotate(
      "text",
      x = xlim_use[1] + diff(xlim_use) * logrank_label_x_frac,
      y = logrank_label_y,
      hjust = 0,
      size = logrank_label_font_size,
      label = stat_label
    )

  # ============================================================
  # 12. Risk table：统一 x 轴，保留网格，只改字体
  # ============================================================
  p$table <- p$table +
    common_x_scale()

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

  p$table$theme$axis.title.y <- ggplot2::element_text(angle = 90,
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
    valid_genes = valid_genes,
    missing_genes = missing_genes,
    median_cutoff = threshold,
    score_method = score_method,
    cutoff_rule = cutoff_rule,
    target_cell = target_cell,
    target_cell_label = target_cell_label,
    target_stage = target_stage,
    logrank_p = logrank_p,
    logrank_label = stat_label,
    fit = fit,
    survdiff = surv_diff,
    break_time_by = break_time_by,
    xlim = xlim_use,
    x_breaks = integer_breaks
  ))
}

# ===== batch6.ipynb cell 1017 (code) =====
res_add <- plot_geneset_pseudobulk_survival_v3(
  seurat_obj = seurat_merge,
  gene_set = immune_programs_targeted$CD4_LAG3__checkpoint_exhaustion$genes,
  target_cell = immune_programs_targeted$CD4_LAG3__checkpoint_exhaustion$subtypes,
  cell_col = "cell_type3",
  target_stage = "D0",
  score_name = "MHCII_maturation",
  surv_time_col = "pfs_time",
  surv_status_col = "pfs_status",
  score_method = "AddModuleScore",risk_table_height = 0.25
)

# ===== batch6.ipynb cell 1019 (code) =====
options(repr.plot.width = 6, repr.plot.height = 7)
res_add$plot + ggtitle("MHCII_maturation in DC_cDC2 at D0")

# ===== batch6.ipynb cell 1020 (code) =====
pdf("plot_v4/MHCII_maturation in DC_cDC2 at D0.pdf", width = 6, height = 7)
print(res_add$plot + ggtitle("MHCII_maturation in DC_cDC2 at D0"))
dev.off()

# ===== batch6.ipynb cell 1021 (code) =====
write.csv(res_add$data, "plot_v4/fig2g_3.csv")

# ===== batch6.ipynb cell 1022 (code) =====
ggsave("plot_v4/Checkpoint exhaustion in CD4_LAG3 at D0_PFS.pdf", width = 6, height = 7)

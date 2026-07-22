# AUTO-GENERATED FROM batch6.ipynb.
# Keep this file synchronized with records/function_provenance.csv.

# ===== Original batch6.ipynb cell 408 =====
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

# ===== Original batch6.ipynb cell 420 =====
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

# ===== Original batch6.ipynb cell 830 =====
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

# AUTO-GENERATED FROM batch6.ipynb.
# Keep this file synchronized with records/function_provenance.csv.

# ===== Original batch6.ipynb cell 938 =====
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

# ===== Original batch6.ipynb cell 1020 =====
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

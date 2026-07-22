# AUTO-GENERATED FROM batch6.ipynb.
# Keep this file synchronized with records/function_provenance.csv.

# ===== Original batch6.ipynb cell 959 =====
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

# ===== Original batch6.ipynb cell 961 =====
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

# ===== Original batch6.ipynb cell 968 =====
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

# ===== Original batch6.ipynb cell 971 =====
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

# ===== Original batch6.ipynb cell 983 =====
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

# ===== Original batch6.ipynb cell 991 =====
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

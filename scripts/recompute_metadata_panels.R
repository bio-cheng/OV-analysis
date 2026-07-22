#!/usr/bin/env Rscript
# Clean, metadata-only recomputation for panels whose inputs do not require the 4-GB Seurat object.
# These are analytical checks, not a claim of pixel-identical regeneration of the manuscript artwork.

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
if (is.na(script_file) || !nzchar(script_file)) stop("Run with Rscript analysis/recompute_metadata_panels.R")
script_dir <- dirname(normalizePath(script_file))
package_dir <- dirname(script_dir)
project_dir <- normalizePath(file.path(package_dir, "..", ".."))
output_dir <- file.path(package_dir, "outputs", "recomputed_metadata")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (pkg in c("dplyr", "tidyr", "ggplot2", "ggrepel")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing package: ", pkg)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
})

meta_file <- file.path(project_dir, "batch_6_merge_meta_filted_02.2.rds")
meta <- readRDS(meta_file)
required <- c("orig.ident", "patient_id", "cell_type3", "group_rn", "sample_type_rn")
missing <- setdiff(required, colnames(meta))
if (length(missing)) stop("Metadata columns missing: ", paste(missing, collapse = ", "))

calc_frequency <- function(data, group_col, group_a, group_b) {
  use <- data %>%
    transmute(sample_id = as.character(orig.ident), patient_id = as.character(patient_id),
              cell_type = as.character(cell_type3), group = as.character(.data[[group_col]])) %>%
    filter(!is.na(sample_id), !is.na(cell_type), group %in% c(group_a, group_b))
  sample_info <- use %>% distinct(sample_id, patient_id, group)
  totals <- use %>% count(sample_id, name = "total_cells")
  all_types <- sort(unique(use$cell_type))
  grid <- tidyr::crossing(sample_info, cell_type = all_types)
  counts <- use %>% count(sample_id, cell_type, name = "cell_count")
  prop <- grid %>% left_join(totals, by = "sample_id") %>% left_join(counts, by = c("sample_id", "cell_type")) %>%
    mutate(cell_count = coalesce(cell_count, 0L), proportion = cell_count / total_cells)
  stats <- prop %>% group_by(cell_type) %>% group_modify(~ {
    x <- .x$proportion[.x$group == group_a]
    y <- .x$proportion[.x$group == group_b]
    p <- tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
    data.frame(mean_group_a = mean(x), mean_group_b = mean(y), n_group_a = length(x), n_group_b = length(y), p_value = p)
  }) %>% ungroup() %>%
    mutate(log2FC = log2((mean_group_a + 1e-6) / (mean_group_b + 1e-6)), p_adj_BH = p.adjust(p_value, method = "BH"))
  list(sample_proportions = prop, statistics = stats)
}

plot_volcano <- function(stats, title, filename) {
  d <- stats %>% mutate(significance = case_when(
    p_value < 0.05 & log2FC > 0 ~ "Higher in first group",
    p_value < 0.05 & log2FC < 0 ~ "Higher in second group",
    TRUE ~ "Not significant"
  ))
  p <- ggplot(d, aes(log2FC, -log10(p_value), color = significance)) +
    geom_point(aes(size = abs(log2FC)), alpha = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    ggrepel::geom_text_repel(data = filter(d, p_value < 0.05), aes(label = cell_type), max.overlaps = Inf, size = 3) +
    scale_color_manual(values = c("Higher in first group" = "#C62828", "Higher in second group" = "#1565C0", "Not significant" = "grey60")) +
    labs(title = title, x = "log2 fold change", y = "-log10(Wilcoxon P)", color = NULL, size = "|log2FC|") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
  ggsave(filename, p, width = 6, height = 6)
}

fig2c <- calc_frequency(meta, "group_rn", "PFS>6_D0", "PFS<6_D0")
write.csv(fig2c$statistics, file.path(output_dir, "Fig2c_baseline_cell_frequency_recomputed.csv"), row.names = FALSE)
plot_volcano(fig2c$statistics, "Fig. 2c - baseline cell-frequency change (PFS >6 vs <6 months)", file.path(output_dir, "Fig2c_baseline_cell_frequency_recomputed.pdf"))

meta_d0_d7 <- meta %>% mutate(stage_compare = as.character(sample_type_rn))
fig3a <- calc_frequency(meta_d0_d7, "stage_compare", "D7", "D0")
write.csv(fig3a$statistics, file.path(output_dir, "Fig3a_D7_vs_D0_cell_frequency_recomputed.csv"), row.names = FALSE)
plot_volcano(fig3a$statistics, "Fig. 3a - cell-frequency change (D7 vs D0)", file.path(output_dir, "Fig3a_D7_vs_D0_cell_frequency_recomputed.pdf"))

patient_proportions <- function(data, stage, cells) {
  use <- data %>% filter(sample_type_rn == stage, !is.na(patient_id), !is.na(cell_type3)) %>%
    transmute(patient_id = as.character(patient_id), cell_type = as.character(cell_type3))
  totals <- use %>% count(patient_id, name = "total_cells")
  grid <- tidyr::crossing(patient_id = totals$patient_id, cell_type = cells)
  grid %>% left_join(totals, by = "patient_id") %>%
    left_join(use %>% filter(cell_type %in% cells) %>% count(patient_id, cell_type, name = "cell_count"), by = c("patient_id", "cell_type")) %>%
    mutate(cell_count = coalesce(cell_count, 0L), proportion = cell_count / total_cells) %>%
    select(patient_id, cell_type, proportion) %>%
    pivot_wider(names_from = cell_type, values_from = proportion, values_fill = 0)
}

plot_correlation <- function(df, x, y, label, filename) {
  test <- cor.test(df[[x]], df[[y]], method = "spearman", exact = FALSE)
  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(size = 2.5, color = "#2A6F97") + geom_smooth(method = "lm", se = TRUE, color = "grey35") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.4,
             label = sprintf("Spearman r = %.2f\nP = %.3g\nn = %d", unname(test$estimate), test$p.value, nrow(df)), size = 4) +
    labs(title = label, x = paste0(x, " proportion"), y = paste0(y, " proportion")) + theme_bw(base_size = 12)
  ggsave(filename, p, width = 6, height = 5)
  data.frame(x = x, y = y, n = nrow(df), spearman_r = unname(test$estimate), p_value = test$p.value)
}

cor_cells <- c("CD4_Treg", "CD4_LAG3", "CD8_TRM", "CD8_NKT")
d7_prop <- patient_proportions(meta, "D7", cor_cells)
write.csv(d7_prop, file.path(output_dir, "Fig3_D7_patient_proportions.csv"), row.names = FALSE)
cor_results <- bind_rows(
  plot_correlation(d7_prop, "CD4_Treg", "CD4_LAG3", "Fig. 3e - D7 CD4_Treg versus CD4_LAG3", file.path(output_dir, "Fig3e_Treg_vs_LAG3_recomputed.pdf")),
  plot_correlation(d7_prop, "CD4_Treg", "CD8_TRM", "Fig. 3h - D7 CD4_Treg versus CD8_TRM", file.path(output_dir, "Fig3h_Treg_vs_CD8_TRM_recomputed.pdf")),
  plot_correlation(d7_prop, "CD4_LAG3", "CD8_NKT", "Fig. 3h - D7 CD4_LAG3 versus CD8_NKT", file.path(output_dir, "Fig3h_LAG3_vs_CD8_NKT_recomputed.pdf"))
)
write.csv(cor_results, file.path(output_dir, "Fig3e_h_correlation_statistics_recomputed.csv"), row.names = FALSE)

message("Metadata-only recomputation completed: ", output_dir)

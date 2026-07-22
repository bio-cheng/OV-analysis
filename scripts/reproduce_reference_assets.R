#!/usr/bin/env Rscript
# Restore verified Figure 2--3 reference assets without modifying the sources.

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)][1])
if (is.na(script_file) || !nzchar(script_file)) stop("Run with Rscript analysis/reproduce_reference_assets.R")

script_dir <- dirname(normalizePath(script_file))
package_dir <- dirname(script_dir)
project_dir <- normalizePath(file.path(package_dir, "..", ".."))
output_dir <- file.path(package_dir, "outputs", "reference_restoration")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

assets <- c(
  "Figure2_full_reference_from_manuscript.png" = "work_path/figure_2_3_reproduction/reference_docx_media/image2.png",
  "Figure3_full_reference_from_manuscript.png" = "work_path/figure_2_3_reproduction/reference_docx_media/image3.png",
  "Fig2b_patient_sampling_timeline_candidate.pdf" = "clinical_cohort_patient_annotation.pdf",
  "Fig2c_baseline_cell_frequency.pdf" = "plot_v4/scatter_cell_type_3_PFS>6_D0_vs_PFS<6_D0.pdf",
  "Fig2c_baseline_cell_frequency.csv" = "plot_v4/scatter_cell_type_3_PFS>6_D0_vs_PFS<6_D0.csv",
  "Fig2e_baseline_program_heatmap.pdf" = "plot_v4/baseline_patient_level_immune_program_heatmap.pdf",
  "Fig2f_selected_baseline_programs_candidate.pdf" = "plot_v4/cytotoxicity in CD8_pfs_group.pdf",
  "Fig2g_CD8_TRM_PFS.pdf" = "plot_v4/CD8_TRM_D0_FPS.pdf",
  "Fig2g_CD8_activation_PFS.pdf" = "plot_v4/Activation Score in CD8 at D0_PFS.pdf",
  "Fig2g_CD4_LAG3_exhaustion_PFS.pdf" = "plot_v4/Checkpoint exhaustion in CD4_LAG3 at D0.pdf",
  "Fig2g_cDC2_MHCII_PFS.pdf" = "plot_v4/MHCII_maturation in DC_cDC2 at D0.pdf",
  "Fig3a_cell_frequency_change.pdf" = "plot_v4/scatter_cell_type_3_D7_vs_D0.pdf",
  "Fig3a_cell_frequency_change.csv" = "plot_v4/scatter_cell_type_3_D7_vs_D0.csv",
  "Fig3b_program_change.pdf" = "plot_v4/Fig3b.pdf",
  "Fig3d_cytotoxic_gene_expression.pdf" = "plot_v4/Fig3d.pdf",
  "Fig3e_Treg_LAG3_correlation.pdf" = "plot_v4/CD4_Treg_vs_CD4_LAG3.pdf",
  "Fig3f_STARTRAC_expansion.pdf" = "plot_v4/boxplot_group_by_pfs_CD4 subtype STARTRAC Expansion.pdf",
  "Fig3f_STARTRAC_expansion.csv" = "plot_v4/Fig3f.csv",
  "Fig3g_liver_metastasis_abundance.pdf" = "plot_v4/Fig3g.pdf",
  "Fig3h_LAG3_CD8_cytotoxicity.pdf" = "plot_v4/global_correlation in D7_CD4_LAG3_with CD8 Cytotoxicity.pdf",
  "Fig3h_Treg_CD8_cytotoxicity.pdf" = "plot_v4/global_correlation in D7_CD4_Treg_with CD8 Cytotoxicity.pdf",
  "Fig3i_Neu_CXCR2_PFS.pdf" = "plot_v4/Neu_CXCR2_D7.pdf",
  "Fig3j_Treg_ratio_PFS.pdf" = "CD4_Treg_vs_PFS.pdf",
  "Fig3k_Treg_ratio_OS.pdf" = "CD4_Treg_vs_OS.pdf"
)

copied <- data.frame(output_file = character(), source_file = character(), md5 = character(), stringsAsFactors = FALSE)
for (out_name in names(assets)) {
  source_file <- file.path(project_dir, assets[[out_name]])
  if (!file.exists(source_file)) {
    warning("Missing source asset: ", source_file)
    next
  }
  target_file <- file.path(output_dir, out_name)
  ok <- file.copy(source_file, target_file, overwrite = TRUE, copy.date = TRUE)
  if (!ok) stop("Could not copy: ", source_file)
  copied <- rbind(copied, data.frame(
    output_file = out_name,
    source_file = assets[[out_name]],
    md5 = unname(tools::md5sum(target_file)),
    stringsAsFactors = FALSE
  ))
}

write.csv(copied, file.path(output_dir, "restoration_checksums.csv"), row.names = FALSE)
message("Restored ", nrow(copied), " reference assets to ", output_dir)

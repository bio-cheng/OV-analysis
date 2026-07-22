# Shared package paths, data loaders and resource-table helpers.
# This file is intentionally small; scientific functions live in numbered files.

# In Jupyter, source() does not reliably expose the sourced file path.  The
# notebooks set this option explicitly; the fallback supports direct R sourcing.
figure_package_dir <- getOption('figure_package_dir', NA_character_)
if (is.na(figure_package_dir) || !dir.exists(figure_package_dir)) {
  this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
  if (!is.na(this_file) && file.exists(this_file)) figure_package_dir <- dirname(dirname(this_file))
}
if (is.na(figure_package_dir) || !dir.exists(figure_package_dir)) {
  stop("Cannot determine figure_2_3_reproduction package directory. Set options(figure_package_dir = '...').")
}
figure_package_dir <- normalizePath(figure_package_dir)
figure_data_dir <- file.path(figure_package_dir, "data")
figure_resource_dir <- file.path(figure_package_dir, "resources")

required_packages <- c(
  "Seurat", "dplyr", "tidyr", "tibble", "ggplot2", "ggpubr", "ggrepel",
  "Matrix", "survival", "survminer", "ComplexHeatmap", "circlize", "patchwork",
  "pheatmap", "ggsci"
)

check_figure_packages <- function(packages = required_packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
  invisible(TRUE)
}

data_path <- function(name) {
  path <- file.path(figure_data_dir, name)
  if (!file.exists(path)) stop("Missing linked input: ", path)
  path
}

load_final_metadata <- function() readRDS(data_path("metadata_final.rds"))

load_seurat_final <- function(attach_final_metadata = TRUE) {
  seu <- readRDS(data_path("seurat_filtered.rds"))
  if (attach_final_metadata) {
    meta <- load_final_metadata()
    common <- intersect(colnames(seu), rownames(meta))
    if (length(common) != ncol(seu)) {
      stop("The Seurat object and final metadata do not have identical cell IDs.")
    }
    seu@meta.data <- meta[colnames(seu), , drop = FALSE]
  }
  seu
}

load_tcr_annotation <- function() readRDS(data_path("tcr_annotation.rds"))

make_clinical_table <- function(meta) {
  required <- c("patient_id", "pfs_time", "pfs_status", "pfs_group", "response")
  absent <- setdiff(required, colnames(meta))
  if (length(absent)) stop("Final metadata lacks: ", paste(absent, collapse = ", "))
  meta |>
    dplyr::select(dplyr::any_of(c(required, "os_time", "is_death", "is_death_e", "lm"))) |>
    dplyr::distinct() |>
    dplyr::group_by(patient_id) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), dplyr::first), .groups = "drop")
}

prepare_figure2b_metadata <- function(seu) {
  seu$sample_type <- as.character(seu$sample_type_rn)
  seu$pfs_6m <- ifelse(seu$pfs_time >= 6, "FPS>=6 month", "FPS<6 month")
  seu
}

make_resource_dir <- function(figure_name) {
  out <- file.path(figure_resource_dir, figure_name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}

write_resource_table <- function(x, filename, figure_name) {
  out <- file.path(make_resource_dir(figure_name), filename)
  utils::write.csv(as.data.frame(x), out, row.names = FALSE)
  message("Resource table: ", out)
  invisible(out)
}

save_ggplot_pdf <- function(plot, filename, figure_name, width, height) {
  out <- file.path(make_resource_dir(figure_name), filename)
  ggplot2::ggsave(out, plot = plot, width = width, height = height)
  message("Plot: ", out)
  invisible(out)
}

save_survival_pdf <- function(surv_result, filename, figure_name, width = 6, height = 7) {
  out <- file.path(make_resource_dir(figure_name), filename)
  grDevices::pdf(out, width = width, height = height)
  print(surv_result$plot)
  grDevices::dev.off()
  message("Survival plot: ", out)
  invisible(out)
}

program_resource_table <- function(program_list) {
  dplyr::bind_rows(lapply(names(program_list), function(nm) {
    x <- program_list[[nm]]
    data.frame(program = nm, subtypes = paste(x$subtypes, collapse = ";"),
               genes = paste(x$genes, collapse = ";"), n_genes = length(x$genes),
               stringsAsFactors = FALSE)
  }))
}

make_d7_d0_ratio_table <- function(meta, target_cell, outcome_col, group_col = "lm",
                                   death_col = "is_death_e", cell_type_col = "cell_type3") {
  meta |>
    dplyr::filter(sample_type_rn %in% c("D0", "D7")) |>
    dplyr::group_by(patient_id, .data[[group_col]], .data[[outcome_col]], .data[[death_col]], sample_type_rn) |>
    dplyr::summarise(total_cells = dplyr::n(), target_cells = sum(.data[[cell_type_col]] == target_cell),
                     percentage = 100 * target_cells / total_cells, .groups = "drop") |>
    tidyr::pivot_wider(names_from = sample_type_rn, values_from = percentage) |>
    dplyr::filter(!is.na(D0), !is.na(D7), D0 > 0) |>
    dplyr::mutate(ratio_D7_D0 = D7 / D0)
}

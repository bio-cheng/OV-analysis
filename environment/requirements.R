# R package requirements for the Figure 2/3 notebooks.
required_packages <- c(
  "Seurat", "dplyr", "tidyr", "tibble", "ggplot2", "ggpubr", "ggrepel",
  "Matrix", "survival", "survminer", "ComplexHeatmap", "circlize", "patchwork",
  "pheatmap", "ggsci", "readxl", "openxlsx"
)

check_required_packages <- function(packages = required_packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Install missing R packages: ", paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

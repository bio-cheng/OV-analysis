# OV-analysis

Reproducible analysis code for Figures 2 and 3 of the oncolytic-virus peripheral-blood study. The original exploratory notebook (`batch6.ipynb`) was stateful and mixed data loading, calculations, plotting, and revised clinical annotations. This repository reorganizes the figure-specific workflow into reusable R functions and two ordered Jupyter notebooks.

> 本仓库仅包含代码、数据规范和图—代码对应记录；不包含受控临床数据、Seurat 对象、原始论文或生成图。

## Repository layout

| Path | Purpose |
| --- | --- |
| `notebooks/Figure2.ipynb` | Figure 2a–g, ordered by panel. |
| `notebooks/Figure3.ipynb` | Figure 3a–l, ordered by panel. |
| `functions/` | Shared R functions for cell frequencies, program scoring, survival, TCR/STARTRAC, correlations, expression, and metastasis analyses. |
| `data/` | Empty by design. Add local symlinks/copies following [`data/README.md`](data/README.md). |
| `records/panel_code_mapping.csv` | Panel-to-notebook/function/input/resource-table index. |
| `records/manual_review.md` | User-confirmed decisions used to resolve ambiguous source calls. |
| `provenance/` | Relevant source cells extracted from the original notebook, retaining source cell indices. |
| `scripts/` | Helper scripts for provenance extraction and independent restoration/metadata checks. |
| `environment/requirements.R` | Required R packages. |

## Data policy

Large Seurat objects, clinical spreadsheets, TCR annotations, and the manuscript are not versioned. They may contain controlled or sensitive data. Place local **symlinks** (recommended) or local copies in `data/` using the exact filenames listed in [`data/README.md`](data/README.md). The `.gitignore` prevents those files from being committed.

## Environment

Use the R Jupyter kernel named `seurat_fresh` (R 4.3.2 in the source environment), with packages listed in `environment/requirements.R`.

```r
source("environment/requirements.R")
check_required_packages()
```

## Run order

1. Configure local input links in `data/`.
2. Open `notebooks/Figure2.ipynb` with the `seurat_fresh` kernel and run from top to bottom.
3. Open `notebooks/Figure3.ipynb` with the same kernel and run from top to bottom.

The notebooks create `resources/Figure2/` and `resources/Figure3/` locally. Each analytical panel exports its PDF and/or resource table there; these generated outputs are intentionally ignored by Git.


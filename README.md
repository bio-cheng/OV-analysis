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

## Confirmed analysis choices

- Survival event field: `pfs_status` in the final metadata object.
- Figure 2e–f scoring: `AddModuleScore` at baseline (`D0`).
- Figure 2b: `plot_cns_clinical_matrix`.
- Figure 2g checkpoint-exhaustion survival: CD4_LAG3 gene set with `pfs_status`.
- Figure 3c: four confirmed programs, `AddModuleScore`, D0/D7/postICI comparison.
- Figure 3f: CD4_Treg/CD4_LAG3 STARTRAC expansion plot.
- Figure 3h: D7 CD4_Treg proportion versus CD8 cytotoxicity score in CD8 cytotoxicity subtypes.
- Figure 3k: CD4_Treg D7/D0 ratio using `os_time`; the supplied source title remains `CD4_Treg Ratio vs PFS` while the y-axis is OS time.

Details and provenance are recorded in [`records/manual_review.md`](records/manual_review.md) and [`records/panel_code_mapping.csv`](records/panel_code_mapping.csv).

## Validation status

- All extracted R function files parse and source successfully in `seurat_fresh`.
- Both notebooks parse as R and target the `seurat_fresh` kernel.
- Figure 2 was executed against the full Seurat object: panels 2b–e and 2g produced PDFs/resource tables. The Figure 2f failure was traced to comparison labels and fixed in the committed notebook (`PFS <6 months` vs `PFS >6 months`); rerun it after configuring data locally.
- Figure 3 has syntax/function checks but has not been fully executed in this repository snapshot.

## Notes

Panels 2a, 2d, and 3l are manuscript schematics rather than direct computational outputs. The code instead exports their associated program-definition resource table where applicable.

# Local input data

This directory is intentionally empty in Git. Create local symlinks (or copies) with the names below before running the notebooks. The filenames are used by `functions/00_setup.R`.

```bash
ln -s /path/to/batch_6_merge_seurat_obj_filted_02.2.rds data/seurat_filtered.rds
ln -s /path/to/batch_6_merge_seurat_obj.rds data/seurat_full.rds
ln -s /path/to/batch_6_merge_meta_filted_02.2.rds data/metadata_final.rds
ln -s /path/to/tmp_seuart_batch6_meta.RDS data/metadata_intermediate.rds
ln -s /path/to/merged_tcr_annotation.RDS data/tcr_annotation.rds
ln -s /path/to/IDOV溶瘤病毒外周血标本-zhc-20260330.xlsx data/clinical_20260330.xlsx
ln -s '/path/to/IDOV溶瘤病毒外周血标本-zhc-20260420-更新(1).xlsx' data/clinical_20260420.xlsx
```

`seurat_filtered.rds`, `metadata_final.rds`, and `tcr_annotation.rds` are required for the Figure 2/3 notebooks. The remaining inputs are retained for provenance and auxiliary scripts.

Do not commit these links or any underlying data files.

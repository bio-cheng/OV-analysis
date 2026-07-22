# Manual-review record

## User-confirmed mappings applied on 2026-07-18

- Final survival-event field: `pfs_status` in `batch_6_merge_meta_filted_02.2.rds`.
- Figure 2e–f scores: `AddModuleScore`; Figure 2 uses baseline D0.
- Figure 2b: `plot_cns_clinical_matrix`.
- Figure 2g checkpoint-exhaustion panel: `plot_geneset_pseudobulk_survival_v3` with the CD4_LAG3 checkpoint gene set and `pfs_status`.
- Figure 3c: `plot_program_boxplot_from_score_res` using the supplied four programs.
- Figure 3f: `calc_startrac_expansion` followed by the supplied CD4_Treg/CD4_LAG3 faceted box plot.
- Figure 2c/Figure 3a: `calculate_cluster_diff` with the confirmed group parameters.
- Figure 3h: at D7, use CD4_Treg proportion against the CD8 cytotoxicity gene set (`CD8_all__cytotoxicity`) scored in its defined CD8 subtypes.
- Figure 3k: use `plot_ratio_pfs_correlation(..., pfs_col = 'os_time', group_col = 'lm', death_col = 'is_death_e')`. The supplied original title text is retained as `CD4_Treg Ratio vs PFS`; the y-axis is `OS Time (Months)`.

## No remaining analytical ambiguity

| Item | Handling |
| --- | --- |
| Fig. 2a/2d and Fig. 3l | These are manuscript artwork, not direct code output; the verified full-figure images are retained in `outputs/reference_restoration/`. |

The two previous Fig. 3h/3k review items are now resolved and are reflected in `notebooks/Figure3.ipynb` and `records/panel_code_mapping.csv`.

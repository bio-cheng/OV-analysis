# AUTO-GENERATED PROVENANCE FILE -- DO NOT EDIT
# Source: ../../batch6.ipynb
# Each block was copied verbatim; it may depend on interactive
# objects created outside this extraction. Use the source index and
# README before running individual blocks.


# ===== batch6.ipynb cell 93 (code) =====
seurat_merge <- readRDS("batch_6_merge_seurat_obj_filted_02.2.rds")

# ===== batch6.ipynb cell 94 (code) =====
meta <- readRDS("batch_6_merge_meta_filted_02.2.rds")

# ===== batch6.ipynb cell 95 (code) =====
length(unique(meta$orig.ident))

# ===== batch6.ipynb cell 96 (code) =====
saveRDS(seurat_merge@meta.data, "batch_6_merge_meta_filted_02.2.rds")

# ===== batch6.ipynb cell 97 (code) =====
saveRDS(seurat_merge@meta.data, "tmp_seuart_batch6_meta.RDS")

# ===== batch6.ipynb cell 98 (code) =====
1

# ===== batch6.ipynb cell 99 (code) =====


# ===== batch6.ipynb cell 100 (code) =====
saveRDS(seurat_merge, "batch_6_merge_seurat_obj.rds")

# ===== batch6.ipynb cell 101 (code) =====
seurat_merge <- readRDS("batch_6_merge_seurat_obj.rds")

# ===== batch6.ipynb cell 102 (code) =====
table(seurat_merge$cell_type3)

# ===== batch6.ipynb cell 103 (code) =====


# ===== batch6.ipynb cell 104 (code) =====
meta <- readRDS("tmp_seuart_batch6_meta.RDS")

# ===== batch6.ipynb cell 105 (code) =====
seurat_merge@meta.data <- meta

# ===== batch6.ipynb cell 202 (code) =====
tb <- read_xlsx("IDOV溶瘤病毒外周血标本-zhc-20260420-更新(1).xlsx", skip = 1)

# ===== batch6.ipynb cell 203 (code) =====
tb2 <- read_xlsx("IDOV溶瘤病毒外周血标本-zhc.xlsx")

# ===== batch6.ipynb cell 204 (code) =====
tb2$patient_id <- gsub("-01$", "", tb2$`样本编号1`)
tb2$patient_id <- gsub("BZ-", "BZ", tb2$patient_id )

# ===== batch6.ipynb cell 205 (code) =====
p_delta$response <- tb2$`response`[match(p_delta$patient_id, tb$patient_id)]

# ===== batch6.ipynb cell 206 (code) =====
p_delta

# ===== batch6.ipynb cell 207 (code) =====
seurat_merge$os_time <- tb$`OS(m)`[match(seurat_merge$patient_id, tb$patient_id)]
seurat_merge$is_death <- tb$`是否死亡`[match(seurat_merge$patient_id, tb$patient_id)]

# ===== batch6.ipynb cell 208 (code) =====
seurat_merge$is_death_e <- ifelse(seurat_merge$is_death == "是", "Yes", "No")

# ===== batch6.ipynb cell 209 (code) =====
seurat_merge$os_status  <- ifelse(seurat_merge$is_death_e == "Yes", 1, 0)

# ===== batch6.ipynb cell 210 (code) =====
meta$os_time <- tb$`OS(m)`[match(meta$patient_id, tb$patient_id)]
meta$is_death <- tb$`是否死亡`[match(meta$patient_id, tb$patient_id)]

# ===== batch6.ipynb cell 211 (code) =====
meta$os_statu <- ifelse(meta$is_death == "是", 1, 0)

# ===== batch6.ipynb cell 212 (code) =====


# ===== batch6.ipynb cell 213 (code) =====
table(seurat_merge$is_death_e)

# ===== batch6.ipynb cell 214 (code) =====


# ===== batch6.ipynb cell 215 (code) =====
tb$patient_id <- gsub("-01$", "", tb$`样本编号1`)
tb$patient_id <- gsub("BZ-", "BZ", tb$patient_id )

# ===== batch6.ipynb cell 216 (code) =====
p_delta$os_time <- tb$`OS(m)`[match(p_delta$patient_id, tb$patient_id)]
p_delta$is_death <- tb$`是否死亡`[match(p_delta$patient_id, tb$patient_id)]

# ===== batch6.ipynb cell 217 (code) =====
p_delta$pfs_time <- tb$`PFS(m)`[match(p_delta$patient_id, tb$patient_id)]

# ===== batch6.ipynb cell 218 (code) =====
p_delta$is_death_e <- ifelse(p_delta$is_death == "是", "Yes", "No")

# ===== batch6.ipynb cell 219 (code) =====
unique(meta$)

# ===== batch6.ipynb cell 220 (code) =====


# ===== batch6.ipynb cell 221 (code) =====
head(meta)

# ===== batch6.ipynb cell 222 (code) =====


# ===== batch6.ipynb cell 223 (code) =====


# ===== batch6.ipynb cell 224 (code) =====
p_delta$lm <- meta$lm[match(p_delta$patient_id, meta$patient_id)]

# ===== batch6.ipynb cell 225 (code) =====
p_delta$response <- meta$response[match(p_delta$patient_id, meta$patient_id)]

# ===== batch6.ipynb cell 428 (code) =====
lm_patient <- c("002", "003", "004", "007", "013", "017", "026", "027", "033", "040")

# ===== batch6.ipynb cell 429 (code) =====
lm_patient <- paste0("BZ", lm_patient)

# ===== batch6.ipynb cell 430 (code) =====
seurat_merge$lm <- ifelse(seurat_merge$patient_id %in% lm_patient, "LM", "noLM")

# ===== batch6.ipynb cell 431 (code) =====
seurat_merge$lm_group <- paste0(seurat_merge$lm, "_", seurat_merge$sample_type_rn)

# ===== batch6.ipynb cell 432 (code) =====
unique(seurat_merge$lm_group)

# ===== batch6.ipynb cell 433 (code) =====


# ===== batch6.ipynb cell 434 (code) =====


# ===== batch6.ipynb cell 435 (code) =====


# ===== batch6.ipynb cell 436 (code) =====
seurat_merge$sample_type_rn[seurat_merge$sample_type2 == "02.2"] <- "beforeICI"

# ===== batch6.ipynb cell 437 (code) =====
table(seurat_merge$sample_type_rn, seurat_merge$sample_type2)

# ===== batch6.ipynb cell 438 (code) =====
seurat_merge

# ===== batch6.ipynb cell 439 (code) =====
min(seurat_merge$nFeature_RNA)

# ===== batch6.ipynb cell 440 (code) =====
pkgs <- c("Seurat", "harmony", "DoubletFinder")

data.frame(
  package = pkgs,
  version = sapply(pkgs, function(x) {
    if (requireNamespace(x, quietly = TRUE)) {
      as.character(packageVersion(x))
    } else {
      NA
    }
  })
)

# ===== batch6.ipynb cell 441 (code) =====


# ===== batch6.ipynb cell 442 (code) =====
seurat_merge$group_rn <- paste0(seurat_merge$pfs_group, "_", seurat_merge$sample_type_rn)

# ===== batch6.ipynb cell 443 (code) =====
seurat_merge$pfs_group <- ifelse(seurat_merge$pfs_time < 6, "PFS<6", "PFS>6")

# ===== batch6.ipynb cell 630 (code) =====
tcr_seq <- readRDS("merged_tcr_annotation.RDS")

# ===== batch6.ipynb cell 631 (code) =====
seurat_merge$cell_type_major <- unlist(lapply(strsplit(seurat_merge$cell_type3, "_"), "[", 1))

# ===== batch6.ipynb cell 632 (code) =====
cd8_obj <- seurat_merge[,seurat_merge$cell_type_major == "CD8"]

# ===== batch6.ipynb cell 633 (code) =====
unique(cd8_obj$cell_type3)

# ===== batch6.ipynb cell 634 (code) =====
cd8_obj <- seurat_merge[,seurat_merge$cell_type_major == "CD4"]

# ===== batch6.ipynb cell 635 (code) =====
rm(seurat_merge)

# ===== batch6.ipynb cell 636 (code) =====
gc()

# ===== batch6.ipynb cell 637 (code) =====
tcr_seq

# ===== batch6.ipynb cell 638 (code) =====


# ===== batch6.ipynb cell 639 (code) =====
cd8_obj$tcr_id <- tcr_seq$combine_seq[match(colnames(cd8_obj), tcr_seq$barcode_id)]
cd8_obj$n_clonal <- tcr_seq$n_clonal[match(colnames(cd8_obj), tcr_seq$barcode_id)]
cd8_obj$n_clonal_s <- tcr_seq$n_clonal_s[match(colnames(cd8_obj), tcr_seq$barcode_id)]
cd8_obj@meta.data$n_clonal_s[is.na(cd8_obj@meta.data$n_clonal_s)] <- 0
cd8_obj@meta.data$n_clonal[is.na(cd8_obj@meta.data$n_clonal)] <- 0

# ===== batch6.ipynb cell 640 (code) =====


# ===== batch6.ipynb cell 641 (code) =====
cd8_obj$cell_type3[cd8_obj$cell_type3 == "CD4_Activated"] <- "CD4_LAG3"

# ===== batch6.ipynb cell 642 (code) =====


# ===== batch6.ipynb cell 643 (code) =====
cd8_obj@meta.data$expand <- ifelse(cd8_obj$n_clonal_s > 1, "Exp", "NoExp")

# ===== batch6.ipynb cell 925 (code) =====
survival <- read_xlsx("IDOV溶瘤病毒外周血标本-zhc-20260330.xlsx")

# ===== batch6.ipynb cell 926 (code) =====
unique(seurat_merge$patient_id)

# ===== batch6.ipynb cell 927 (code) =====
seurat_merge$pfs_time <- 

# ===== batch6.ipynb cell 928 (code) =====
survival$patient_id <- gsub("-01$", "", survival$`样本编号1`)
survival$patient_id <- gsub("-", "", survival$patient_id)

# ===== batch6.ipynb cell 929 (code) =====
seurat_merge$pfs_time <- survival$`PFS(m)`[match(seurat_merge$patient_id, survival$patient_id)]

# ===== batch6.ipynb cell 930 (code) =====
seurat_merge$pfs_statu <- 1

# ===== batch6.ipynb cell 1003 (code) =====
pd <- c("BZ005", "BZ017", "BZ021")

# ===== batch6.ipynb cell 1004 (code) =====
table(seurat_merge$is_death )

# ===== batch6.ipynb cell 1005 (code) =====
seurat_merge$pfs_status <- ifelse(seurat_merge$patient_id %in% pd | seurat_merge$is_death == "是", 1, 0)

# ===== batch6.ipynb cell 1006 (code) =====
unique(seurat_merge$patient_id )

# ===== batch6.ipynb cell 1007 (code) =====
seurat_merge$pfs_status[seurat_merge$patient_id == "BZ024"] <- 1
seurat_merge$pfs_status[seurat_merge$patient_id == "BZ029"] <- 0

# ===== batch6.ipynb cell 1008 (code) =====
table(seurat_merge$pfs_status)

# ===== batch6.ipynb cell 1009 (code) =====
seurat_merge$pfs_status[is.na(seurat_merge$pfs_status)] <- 0

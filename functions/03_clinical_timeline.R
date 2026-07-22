# AUTO-GENERATED FROM batch6.ipynb.
# Keep this file synchronized with records/function_provenance.csv.

# ===== Original batch6.ipynb cell 798 =====
library(tidyverse)
library(ggsci)

plot_cns_clinical_matrix <- function(seurat_obj) {
  
  # 1. 提取并清理元数据
  meta <- seurat_obj@meta.data %>% 
    filter(sample_type != "02.2") %>%
    select(patient_id, sample_type, response, pfs_6m) %>%
    distinct() %>%
    mutate(timepoint = sample_type) %>%
    filter(!is.na(timepoint))
  
  # 2. 采样状态
  sampling_status <- meta %>%
    complete(patient_id, timepoint) %>%
    mutate(Status = ifelse(!is.na(sample_type), "Sampled", "Missing")) %>%
    select(patient_id, timepoint, Status) %>%
    pivot_wider(names_from = timepoint, values_from = Status)
  
  # 3. 临床信息
  clinical_info <- meta %>%
    select(patient_id, response, pfs_6m) %>%
    distinct()
  
  # 4. 合并并转为绘图长表
  plot_data <- sampling_status %>%
    left_join(clinical_info, by = "patient_id") %>%
    arrange(patient_id, desc(response)) %>%
    mutate(patient_id = factor(patient_id, levels = unique(patient_id))) %>%
    pivot_longer(
      cols = -patient_id, 
      names_to = "Attribute", 
      values_to = "Value"
    ) %>%
    mutate(
      Attribute = factor(
        Attribute, 
        levels = rev(c("D0", "D7", "postICI", "response", "pfs_6m"))
      )
    )
  
  # 5. CNS 级别配色
  cns_colors <- c(
    "Sampled" = "#716db2", 
    "Missing" = "#E0E0E0",
    "R" = "#3CB371", 
    "NR" = "#CD5C5C",
    "FPS>=6 month" = "#4682B4", 
    "FPS<6 month" = "#F4A460"
  )
  
  # 6. 绘图：去掉 tile 边框和黑色网格线
  ggplot(plot_data, aes(x = patient_id, y = Attribute, fill = Value)) +
    geom_tile(color = NA) +
    scale_fill_manual(values = cns_colors, na.value = "white") +
    scale_x_discrete(position = "top") +
    labs(
      title = "Clinical Cohort Landscape", 
      subtitle = "Patient IDs (Columns) vs. Clinical Features (Rows)",
      x = NULL, 
      y = NULL
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(
        angle = 90, 
        hjust = -0.5, 
        vjust = -1, 
        size = 15, 
        family = "sans"
      ),
      axis.text.y = element_text(
        size = 12, 
        face = "bold", 
        color = "black"
      ),
      legend.position = "right",
      legend.text = element_text(
        size = 14, 
        face = "bold", 
        color = "black"
      ),
      plot.title = element_text(
        face = "bold", 
        size = 14
      ),
      strip.background = element_blank()
    )
}

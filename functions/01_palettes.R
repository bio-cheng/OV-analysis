# AUTO-GENERATED FROM batch6.ipynb.
# Keep this file synchronized with records/function_provenance.csv.

# ===== Original batch6.ipynb cell 403 =====
# 定义优化后的高对比度颜色向量 (结合了新版谱系色调逻辑)
cell_type_colors_updated <- c(
  
  # --- B & Plasma 谱系 (重构为紫/深紫色系，避免与 CD4 蓝色冲突) ---
  "B_cell"       = "#6a51a3", # 中紫色
  "B_Naive"      = "#9e9ac8", # 浅灰紫
  "Plasma cells" = "#3f007d", # 极暗紫
  
  # --- CD4 谱系 (蓝-青-深蓝切换) ---
  "CD4_Naive"      = "#084594", # 深蓝
  "CD4_Th17"       = "#41b6c4", # 蓝绿 (Teal)
  "CD4_Tem"        = "#225ea8", # 皇家蓝
  "CD4_Treg"       = "#7fcdbb", # 浅青色
  "CD4_Activated"  = "#0c2c84", # 午夜蓝
  "CD4_Other"      = "#c6dbef", # 极浅蓝
  
  # --- CD8 谱系 (红-橙-紫-粉切换) ---
  "CD8_TRM"        = "#f768a1", # 亮粉色
  "CD8_Tem"        = "#e31a1c", # 正红
  "CD8_Naive"      = "#feb24c", # 橙黄
  "CD8_Prolif"     = "#800026", # 深血红 (用户备注亮粉，保留用户Hex码)
  "CD8_NKT"        = "#8c6bb1", # 中紫色
  "CD8_Stress"     = "#b10026", # 铁锈红
  "CD8_MAIT"       = "#ec7014", # 亮橙/棕 (映射自参考 MAIT)
  "gdT"            = "#cc4c02", # 赭石色 (映射自参考 T_ActEarly，保持在泛T细胞暖色系)
  
  # --- NK 谱系 (绿-黄-浅绿切换) ---
  "NK_C2"          = "#238b45", # 深绿 (映射自参考 NK)
  "NK_Prolif"      = "#74c476", # 嫩草绿
  "NK_XCL1"        = "#addd8e", # 浅黄绿 (映射自参考 Prolif_NK)
  "NK_PTGDS"       = "#006d2c", # 墨绿色 (补充未提及的 NK 亚群)
  
  # --- Mono 谱系 (保持原有的高对比棕/金/橙，与 T 细胞完美区分) ---
  "Mono_Classical_1" = "#8B4513", # 马鞍棕
  "Mono_Classical_2" = "#FFD700", # 亮金色
  "Mono_IFN"         = "#D2691E", # 巧克力色
  "Mono_Intermediate"= "#F5DEB3", # 浅麦色
  "Mono_NonClassical"= "#A0522D", # 黄褐色
  "Mono_Other"       = "#DAA520", # 暗金菊
  
  # --- Neu 谱系 (保持青/水鸭色系) ---

    
  "Neu_CXCR2" = "#008080", # 深水鸭青
  "Neu_IFN"   = "#40E0D0", # 绿松石
  "Neu_IL1R1" = "#4682B4", # 钢蓝
  
  # --- DC 谱系 (保持灰/黑色系) ---
  "DC_DC1" = "#2F4F4F", # 深石板灰
  "DC_DC2" = "#A9A9A9", # 暗灰
  
  # --- 特殊分类 ---
  "New Clones" = "#E0E0E0",  # 极浅灰 (背景色)

  # --- 特殊分类 ---
  "PFS>6" = "#C96A4A",  # 极浅灰 (背景色)
   "PFS<6" = "#2F6F8F",  # 极浅灰 (背景色)

      # --- 特殊分类 ---
  "R" = "#C96A4A",  # 极浅灰 (背景色)
   "NR" = "#2F6F8F",  # 极浅灰 (背景色)

      # --- 特殊分类 ---
  "noLM" = "#C96A4A",  # 极浅灰 (背景色)
   "LM" = "#2F6F8F",  # 极浅灰 (背景色)
    
      # --- 特殊分类 ---
  "D0" = "#EF767A",  # 
   "D7" = "#456990",  # 
   "postICI" = "#48C0AA",
      "CD4" = "#5CA7D3",
  "CD8" = "#225EA8",
  "NK" = "#1D91C0",
  "B" = "#FEB24C",
 # "Plasma cells" = "#F03B20",
  "Mono" = "#7FCDBB",
  "DC" = "#238B45",
  "Neu" = "#BCBDDC"
)

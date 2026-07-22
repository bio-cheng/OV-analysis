# AUTO-GENERATED FROM batch6.ipynb.
# Keep this file synchronized with records/function_provenance.csv.

# ===== Original batch6.ipynb cell 958 =====
immune_programs_targeted <- list(

  # -----------------------------
  # B lineage
  # -----------------------------
  B_cell__BCR_antigen_presentation = list(
    subtypes = c("B_cell", "B_Naive"),
    genes = c(
      "MS4A1", "CD79A", "CD79B", "CD74",
      "HLA-DRA", "HLA-DPA1", "HLA-DPB1",
      "BANK1", "PAX5", "CD40"
    )
  ),

  B_Naive__naive_B_program = list(
    subtypes = c("B_Naive"),
    genes = c(
      "MS4A1", "CD79A", "CD79B", "BANK1",
      "TCL1A", "IGHD", "IGHM", "PAX5", "CD22"
    )
  ),

  Plasma_cells__antibody_secreting = list(
    subtypes = c("Plasma cells"),
    genes = c(
      "PRDM1", "IRF4", "XBP1", "MZB1", "JCHAIN",
      "SDC1", "TNFRSF17", "CD38", "FKBP11",
      "DERL3", "IGHG1", "IGHG3", "IGHM"
    )
  ),


  # -----------------------------
  # CD4 T lineage
  # -----------------------------
  CD4_Naive__naive_memory = list(
    subtypes = c("CD4_Naive"),
    genes = c(
      "CCR7", "SELL", "IL7R", "TCF7", "LEF1",
      "LTB", "MAL", "S1PR1", "BACH2"
    )
  ),

  CD4_Tem__effector_memory = list(
    subtypes = c("CD4_Tem"),
    genes = c(
      "IL7R", "CCL5", "GZMK", "CXCR3",
      "ANXA1", "S100A4", "KLRB1", "IL32"
    )
  ),

  CD4_Th17__Th17_program = list(
    subtypes = c("CD4_Th17"),
    genes = c(
      "RORC", "IL17A", "IL17F", "CCR6",
      "KLRB1", "IL23R", "RORA", "IL22", "CCL20"
    )
  ),

  CD4_Treg__suppressive_program = list(
    subtypes = c("CD4_Treg"),
    genes = c(
      "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
      "TNFRSF18", "ENTPD1", "CCR8", "IL10",
      "TGFB1", "BATF", "PRDM1"
    )
  ),

  CD4_LAG3__checkpoint_exhaustion = list(
    subtypes = c("CD4_LAG3"),
    genes = c(
      "LAG3", "PDCD1", "HAVCR2", "TIGIT",
      "CTLA4", "TOX", "TOX2", "ENTPD1", "CXCL13"
    )
  ),

  CD4_all__Th1_antiviral_helper = list(
    subtypes = c(
      "CD4_LAG3", "CD4_Naive", "CD4_Other",
      "CD4_Tem", "CD4_Th17", "CD4_Treg"
    ),
    genes = c(
      "TBX21", "IFNG", "CXCR3", "IL12RB2",
      "STAT1", "CCL5", "CCR5", "CD40LG"
    )
  ),


  # -----------------------------
  # CD8 T lineage
  # -----------------------------
  CD8_Tem_TRM__cytolytic_effector = list(
    subtypes = c("CD8_Tem", "CD8_TRM", "CD8_NKT", "CD8_Prolif"),
    genes = c(
      "GZMA", "GZMB", "GZMH", "PRF1", "GNLY",
      "NKG7", "CTSW", "CST7", "IFNG"
    )
  ),

  CD8_TRM__tissue_resident_memory = list(
    subtypes = c("CD8_TRM"),
    genes = c(
      "ITGAE", "CXCR6", "ZNF683", "CD69",
      "XCL1", "XCL2", "GZMB", "CXCR3"
    )
  ),

  CD8_Naive__naive_memory = list(
    subtypes = c("CD8_Naive"),
    genes = c(
      "CCR7", "SELL", "IL7R", "TCF7", "LEF1",
      "LTB", "MAL", "S1PR1", "BACH2"
    )
  ),

  CD8_Prolif__proliferation = list(
    subtypes = c("CD8_Prolif"),
    genes = c(
      "MKI67", "TOP2A", "STMN1", "HMGB2",
      "TUBA1B", "TYMS", "PCNA", "MCM5"
    )
  ),

  # CD8_Stress__stress_response = list(
  #   subtypes = c("CD8_Stress"),
  #   genes = c(
  #     "HSPA1A", "HSPA1B", "HSP90AA1", "DNAJB1",
  #     "FOS", "JUN", "DUSP1", "IER2", "ATF3"
  #   )
  # ),

  CD8_all__checkpoint_exhaustion = list(
    subtypes = c(
      "CD8_Naive", "CD8_NKT", "CD8_Prolif",
      "CD8_Stress", "CD8_Tem", "CD8_TRM"
    ),
    genes = c(
      "PDCD1", "LAG3", "HAVCR2", "TIGIT", "CTLA4",
      "TOX", "TOX2", "ENTPD1", "CXCL13", "LAYN"
    )
  ),
    
   CD8_all__activation = list(
    subtypes = c(
      "CD8_Naive", "CD8_NKT", "CD8_Prolif",
      "CD8_Stress", "CD8_Tem", "CD8_TRM"
    ),
    genes = c("CD69", "CCR7", "CD27", "BTLA", "CD40LG", "IL2RA", "CD3E", 
                      "CD47", "EOMES", "GNLY", "GZMA", "GZMB", "PRF1", "IFNG", 
                      "CD8A", "CD8B", "FASLG", "LAMP1", "LAG3", "ICOS"
    )
  ),   

    
   CD8_all__cytotoxicity = list(
    subtypes = c(
      "CD8_Naive", "CD8_NKT", "CD8_Prolif",
      "CD8_Stress", "CD8_Tem", "CD8_TRM"
    ),
    genes = c("PRF1", "IFNG", "GNLY", "NKG7", "GZMB", "GZMA",  "CTSW", "CST7", "FASLG",
  "CD8A", "CD8B", "KLRG1", "KLRK1", "LAMP1")
  ),   

   
  CD8_all__chemokine_trafficking = list(
    subtypes = c(
      "CD8_Naive", "CD8_NKT", "CD8_Prolif",
      "CD8_Stress", "CD8_Tem", "CD8_TRM"
    ),
    genes = c(
      "CXCR3", "CCR5", "CCL5",  "ITGAL", "ITGB2", "SELPLG", "CXCR6"
    )
  ),


  # -----------------------------
  # MAIT
  # -----------------------------
  MAIT__innate_like_cytotoxicity = list(
    subtypes = c("MAIT"),
    genes = c(
      "TRAV1-2", "KLRB1", "SLC4A10", "NCR3",
      "CXCR6", "IL7R", "GZMK", "NKG7", "IFNG", "CCL5"
    )
  ),


  # -----------------------------
  # NK lineage
  # -----------------------------
  NK_C2__NK_cytotoxicity = list(
    subtypes = c("NK_C2", "NK_Prolif", "NK_XCL1"),
    genes = c(
      "NKG7", "GNLY", "PRF1", "GZMB", "GZMH",
      "CTSW", "CST7", "KLRD1", "KLRK1",
      "NCR1", "NCR3", "FCGR3A"
    )
  ),


  NK_Prolif__proliferation = list(
    subtypes = c("NK_Prolif"),
    genes = c(
      "MKI67", "TOP2A", "STMN1", "HMGB2",
      "TUBA1B", "TYMS", "PCNA", "MCM5"
    )
  ),

  NK_XCL1__DC_recruiting_NK = list(
    subtypes = c("NK_XCL1"),
    genes = c(
      "XCL1", "XCL2", "CCL5", "IFNG",
      "NKG7", "GNLY", "GZMB", "PRF1"
    )
  ),



  # -----------------------------
  # DC lineage
  # -----------------------------
  DC_cDC1__cross_presentation = list(
    subtypes = c("DC_cDC1"),
    genes = c(
      "CLEC9A", "XCR1", "BATF3", "IRF8", "CADM1",
      "THBD", "WDFY4", "TAP1", "TAP2", "B2M",
      "HLA-A", "HLA-B", "HLA-C", "PSMB8", "PSMB9"
    )
  ),

  DC_cDC2__MHCII_maturation = list(
    subtypes = c("DC_cDC2"),
    genes = c(
      "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
      "CD74", "CIITA", "CD80", "CD86", "CD40",
      "CD83", "CCR7", "LAMP3"
    )
  ),

  DC_all__antiviral_ISG = list(
    subtypes = c("DC_cDC1", "DC_cDC2"),
    genes = c(
      "STAT1", "STAT2", "IRF7", "ISG15",
      "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2",
      "OAS1", "OAS2", "OAS3", "RSAD2",
      "IFI6", "IFI27", "IFI44L"
    )
  ),

  DC_all__viral_sensing = list(
    subtypes = c("DC_cDC1", "DC_cDC2"),
    genes = c(
      "TLR3", "TLR7", "TLR8", "TLR9", "MYD88",
      "TICAM1", "DDX58", "IFIH1", "MAVS",
      "MB21D1", "TMEM173", "TBK1", "IRF3", "IRF7"
    )
  ),


  # -----------------------------
  # Monocyte lineage
  # -----------------------------
  Mono_Classical__inflammatory_monocyte = list(
    subtypes = c("Mono_Classical_1", "Mono_Classical_2"),
    genes = c(
      "S100A8", "S100A9", "S100A12", "LST1",
      "FCN1", "VCAN", "IL1B", "TNF", "CXCL8",
      "CCL2", "NLRP3", "NFKBIA"
    )
  ),

  Mono_IFN__ISG_antiviral = list(
    subtypes = c("Mono_IFN"),
    genes = c(
      "STAT1", "STAT2", "IRF7", "ISG15",
      "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2",
      "OAS1", "OAS2", "OAS3", "RSAD2",
      "IFI6", "IFI27", "IFI44L"
    )
  ),

  Mono_Intermediate__APC_inflammatory = list(
    subtypes = c("Mono_Intermediate"),
    genes = c(
      "CD14", "FCGR3A", "HLA-DRA", "HLA-DRB1",
      "CD74", "LST1", "MS4A7", "IL1B", "TNF"
    )
  ),

  Mono_NonClassical__patrolling_surveillance = list(
    subtypes = c("Mono_NonClassical"),
    genes = c(
      "FCGR3A", "MS4A7", "CX3CR1", "LST1",
      "IFITM3", "TYROBP", "AIF1", "LILRB1",
      "HLA-DRA", "CD74"
    )
  ),

  Mono_all__MDSC_like_positive = list(
    subtypes = c(
      "Mono_Classical_1", "Mono_Classical_2", "Mono_IFN",
      "Mono_Intermediate", "Mono_NonClassical", "Mono_Other"
    ),
    genes = c(
      "CD14", "ITGAM", "CD33", "S100A8", "S100A9",
      "VCAN", "FCN1", "IL1B", "ARG1", "PTGS2", "IL10"
    )
  ),

  Mono_all__HLA_DR_APC = list(
    subtypes = c(
      "Mono_Classical_1", "Mono_Classical_2", "Mono_IFN",
      "Mono_Intermediate", "Mono_NonClassical", "Mono_Other"
    ),
    genes = c("HLA-DRA", "HLA-DRB1", "CD74", "HLA-DPA1", "HLA-DPB1")
  ),


  # -----------------------------
  # Neutrophil lineage
  # -----------------------------
  Neu_CXCR2__mature_neutrophil_chemotaxis = list(
    subtypes = c("Neu_CXCR2"),
    genes = c(
      "CXCR2", "FCGR3B", "CSF3R", "S100A8",
      "S100A9", "MMP9", "LCN2", "CEACAM8"
    )
  ),

  Neu_IFN__ISG_antiviral = list(
    subtypes = c("Neu_IFN"),
    genes = c(
      "STAT1", "STAT2", "IRF7", "ISG15",
      "IFIT1", "IFIT2", "IFIT3", "MX1", "MX2",
      "OAS1", "OAS2", "OAS3", "RSAD2",
      "IFI6", "IFI27", "IFI44L"
    )
  ),

  Neu_IL1R1__inflammatory_neutrophil = list(
    subtypes = c("Neu_IL1R1"),
    genes = c(
      "IL1R1", "IL1B", "CXCL8", "S100A8",
      "S100A9", "S100A12", "NFKBIA", "TNF", "NLRP3"
    )
  ),

  Neu_all__PMN_MDSC_suppressive = list(
    subtypes = c("Neu_CXCR2", "Neu_IFN", "Neu_IL1R1"),
    genes = c(
      "S100A8", "S100A9", "ARG1", "CXCL8",
      "MMP9", "LCN2", "CEACAM8", "FCGR3B",
      "CSF3R", "CYBB", "NCF1", "NCF2"
    )
  )
)

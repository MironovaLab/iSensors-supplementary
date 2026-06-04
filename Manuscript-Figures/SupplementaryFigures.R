library(Seurat)
library(iSensors)
library(tidyverse)
library(RColorBrewer)
library(scales)

iSensors_obj <- readRDS(file = "02-single-cell-data-analysis/in/iSensors-Martin-Arevalillo-auxin-root2025_mean.rds")


#For supplementary Figure A
p <- DimPlot(
  iSensors_obj,
  reduction = "umap",
  split.by = "orig.ident2",
  combine = TRUE,
  #  label = TRUE,
  label.size = 2.2,        # smaller cluster labels
  repel = TRUE
) +
  scale_color_manual(
    values = Seurat::DiscretePalette(
      length(unique(Idents(iSensors_obj))),
      palette = "polychrome"   # colorblind-friendly, high contrast
    ) 
  )+  guides(color = guide_legend(override.aes = list(size = 2)))

p
p+theme(
  legend.text  = element_text(size = 8),   # smaller cluster names
  axis.line        = element_blank(),
  axis.ticks       = element_blank(),
  axis.text        = element_blank(),
  axis.title       = element_blank(),
  legend.key.height = unit(0.4, "cm"),  # 🔑 reduces vertical spacing
  
)
ggsave("Manuscript-Figures/out/Supplementary_AuxinTreated_DimPlot_splitted.pdf", last_plot(), width = 8, height = 4.0, dpi = 300)

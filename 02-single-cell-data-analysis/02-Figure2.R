library(Seurat)
library(iSensors)
library(tidyverse)
library(RColorBrewer)
library(scales)

seurat_obj <- readRDS(file = "02-single-cell-data-analysis/in/MartinArevalillo2025.rds")
isensors_obj <- readRDS(file = "02-single-cell-data-analysis/in/iSensors-Martin-Arevalillo-auxin-root2025.rds")

#For Figure A
p <- DimPlot(
  seurat_obj,
  reduction = "umap",
  combine = TRUE,
  #  label = TRUE,
  label.size = 2.2,        # smaller cluster labels
  repel = TRUE
) +
  scale_color_manual(
    values = Seurat::DiscretePalette(
      length(unique(Idents(seurat_obj))),
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
ggsave("02-single-cell-data-analysis/out/AuxinTreated_DimPlot_integrated.pdf", last_plot(), width = 6, height = 4.0, dpi = 300)


#For supplementary Figure A
p <- DimPlot(
  seurat_obj,
  reduction = "umap",
  split.by = "orig.ident2",
  combine = TRUE,
  #  label = TRUE,
  label.size = 2.2,        # smaller cluster labels
  repel = TRUE
) +
  scale_color_manual(
    values = Seurat::DiscretePalette(
      length(unique(Idents(seurat_obj))),
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
ggsave("02-single-cell-data-analysis/out/AuxinTreated_DimPlot_splitted.pdf", last_plot(), width = 8, height = 4.0, dpi = 300)


#For Figure B
#Publication ready featureplots
feat = "ATH-aux-trans-ARF"
feat = "ATH-aux-trans-PolarAuxinTransport"

vals <- FetchData(iSensors_obj, vars = feat)[, 1]
lim  <- max(abs(quantile(vals, probs = c(0.02, 0.98), na.rm = TRUE)))
lims <- c(0, lim)


p <- FeaturePlot(
  object = iSensors_obj,
  features = feat,
  split.by = "orig.ident2",
  combine = TRUE
)

p <- p &
  scale_color_gradientn(
    colors = rdylbu5,
    limits = lims,
    oob = scales::squish
  ) &
  guides(
    color = guide_colorbar(
      #      title = "iSensor activity",
      #      title.position = "top",
      title.hjust = 0.4,
      barwidth = unit(3, "cm"),
      barheight = unit(0.4, "cm"),
      ticks = FALSE
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    axis.line  = element_blank(),
    axis.ticks = element_blank(),
    axis.text  = element_blank(),
    axis.title = element_blank()
  )


p
ggsave("out/auxin/AuxinTreated_ARF_FeaturePlot.pdf", last_plot(), width = 4, height = 4.0, dpi = 300)
ggsave("out/auxin/AuxinTreated_PAT_FeaturePlot.pdf", last_plot(), width = 4, height = 4.0, dpi = 300)

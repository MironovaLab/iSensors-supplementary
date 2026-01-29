library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(ggrepel)

output_dir <- "02-single-cell-data-analysis/out/"

seurat_obj <- readRDS("D:/MyPAPERS/iSensors/Data/iSensor_Shahan_output_v0.4_6meta_panels.rds")

seurat_obj@assays
DimPlot(seurat_obj)
meta <- seurat_obj@meta.data

DefaultAssay(seurat_obj_spearman) <- "RNA"


clusters_to_keep <- c("QC", "CSC", "LRP", "LRCEI", "SI", "CEI", 
                      "Columella", "CSCD", "Young LRC", 
                      "Protoxylem_m1", "Protoxylem_m2", "LRC", "Dying LRC", 
                      "Procambium_m1", "Procambium_m2",
                      "Endodermis_m1", "Endodermis_m2",
                      "PPP_m1","PPP_m2", "XPP_m1", "XPP_m2",
                      "Cortex_m1", "Cortex_m2", 
                      "Atrichoblast_m1", "Trichoblast_m1", "Atrichoblast_m2", "Trichoblast_m2")


seurat_obj_spearman <- subset(
  seurat_obj, 
  idents = clusters_to_keep
)
DimPlot(seurat_obj_spearman)

seurat_obj_spearman <- RenameIdents(
  seurat_obj_spearman,
  `PPP_m1` = "Pericycle_m1",
  `PPP_m2` = "Pericycle_m2",
  `XPP_m1` = "Pericycle_m1",
  `XPP_m2` = "Pericycle_m2",
  `LRCEI` = "Initials",
  `SI` = "Initials",
  `CEI` = "Initials"
)

seurat_obj_spearman$seurat_obj_spearman <- Idents(seurat_obj_spearman)

DefaultAssay(seurat_obj_spearman) <- "RNA"
seurat_obj_spearman <- NormalizeData(seurat_obj_spearman)
seurat_obj_spearman <- FindVariableFeatures(seurat_obj_spearman, selection.method = "vst", nfeatures = 3000)
seurat_obj_spearman <- ScaleData(seurat_obj_spearman, vars.to.regress = NULL)  # e.g., c("percent.mt","nCount_RNA")
seurat_obj_spearman <- RunPCA(seurat_obj_spearman, npcs = 50, verbose = FALSE)
dims_use <- 1:30 

seurat_obj_spearman <- FindNeighbors(seurat_obj_spearman, dims = dims_use)
seurat_obj_spearman <- FindClusters(seurat_obj_spearman, resolution = 0.4)  # try a sweep below
seurat_obj_spearman <- RunUMAP(seurat_obj_spearman, dims = dims_use)        # or RunTSNE(seurat_obj, dims = dims_use)

# Quick look
DimPlot(seurat_obj_spearman, reduction = "umap", label = TRUE) + NoLegend()

Idents(seurat_obj_spearman) <- "seurat_obj_spearman"

saveRDS(seurat_obj_spearman, file = "02-single-cell-data-analysis/in/shahan_obj_spearman_iSensors.rds")

#Generating plots ----


paired12 <- c(
  "#e31a1c", "#1f78b4", "#b2df8a", "#33a02c",
  "#fb9a99", "#a6cee3", "#fdbf6f", "#ff7f00",
  "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"
)

if (length(levels(Idents(seurat_obj_spearman))) > length(paired12)) {
  times <- ceiling(length(levels(Idents(seurat_obj_spearman))) / length(paired12))
  paired12 <- rep(paired12, times = times)
}

new_order <- c("Initials", "QC", "CSC", "CSCD", "Columella", "Young LRC", "LRC","Dying LRC", 
               "Protoxylem_m1", "Protoxylem_m2", 
               "Procambium_m1", "Procambium_m2", 
               "Pericycle_m1","Pericycle_m2", 
               "Cortex_m1", "Cortex_m2", 
               "Endodermis_m1", "Endodermis_m2",
               "Atrichoblast_m1", "Atrichoblast_m2", 
               "Trichoblast_m1", "Trichoblast_m2", "LRP")

seurat_obj_spearman$seurat_obj_spearman <- factor(
  seurat_obj_spearman$seurat_obj_spearman,
  levels = new_order
)

p_old <- DimPlot(
  seurat_obj_spearman,
  reduction = "umap",
  group.by = "seurat_obj_spearman",
  cols = paired12[seq_along(levels(Idents(seurat_obj_spearman)))],
  label = FALSE) + theme(plot.title = element_blank()) +  
  guides(
    colour = guide_legend(
      override.aes = list(
        shape = 15,    # square
        size  = 3      # bigger size
      ),
      ncol = 1         # single column legend
    ))

print(p_old)

file_path <- file.path(output_dir, "shahan_dataset_spearman_dimplot.svg")
ggsave(file_path, plot = p_old,
       width = 6, height = 6, dpi = 300)

#FeaturePlot
DefaultAssay(seurat_obj_spearman) <- "iSensor_mean_normed"

head(seurat_obj_spearman@assays$iSensor_mean_normed@counts@Dimnames)
feat <- "AT-aux-trans-A-ARF"
p <- FeaturePlot(
  object = seurat_obj_spearman,
  features = feat,
#  split.by = "orig.ident2",
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
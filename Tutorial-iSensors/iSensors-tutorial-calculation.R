#This script uses the seurat object loaded in iSensors-tutorial-test-load.R script
#run it after the script or upload the saved version:
#seu <- readRDS("path/to/your_scRNAseq_seurat.rds")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(iSensors)
})


DimPlot(
  seu,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  split.by = "orig.ident2"
) + ggtitle("Stem cell niche: control vs auxin")

#Loading iSensors

panel_set <- LoadSensors(setName = 'AuxinPanel', 
                         species = 'AT', 
                         hormone = 'aux', 
                         randomInfo = list('n' = 3, 'sizes' = c(100, 200, 300), majortrend = TRUE))
summary(panel_set$panels)
AT_aux_trans_Transport <- panel_set$panels$AT_aux_trans_Transport

iSensors_obj <- CalcSensors(seu,
                      seurLayer = "data",
                      panelSet = panel_set,
                      signals = c("mean_normed", "median"))
summary(iSensors_obj@assays)

DefaultAssay(iSensors_obj) <- "iSensors_mean_normed"
DimPlot(iSensors_obj, split.by = "orig.ident2",   group.by = "cluster_annot")

DimPlot(iSensors_obj)

FeaturePlot(iSensors_obj, features = "AT-aux-trans-Transport")

library(ggplot2)

# Define a gray–yellow–red palette
gry_yel_red <- c("grey85", "#FEE08B", "#D73027")

p_feat <- FeaturePlot(
  iSensors_obj,
  features = "AT-aux-trans-Transport",
  reduction = "umap",
  cols = gry_yel_red,
  min.cutoff = "q05",
  max.cutoff = "q5",
) &
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    plot.title       = element_blank()
  ) &
  coord_equal()

p_feat

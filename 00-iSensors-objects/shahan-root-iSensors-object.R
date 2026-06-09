#Calculation iSensors for small shahan dataset ----

install.packages("devtools")
devtools::install_github("MironovaLab/iSensors")
library(iSensors)
library(Seurat)
library(ggplot2)
library(scales)
library(Seurat)
library(patchwork)
library(cowplot)
set.seed(100)

#calculating iSensors for the shahan gound truth data ----
seurat_obj <- readRDS("D:/!GitHub/DigitalSensor-Toolbox/iSensors-supplementary/00-iSensors-objects/data/02_ra_integrated_Shahan_20250114.rds")
Assays(seurat_obj)
seurat_obj[["integrated"]] <- NULL
seurat_obj[["ALT"]] <- NULL

clusters_to_keep <- c("QC", "CSC", "LRP", "LRCEI", "SI", "CEI", 
                      "Columella", "CSCD", "Young LRC", "LRC", "Dying LRC", 
                      "Protoxylem_m1", "Protoxylem_m2",   
                      "Procambium_m1", "Procambium_m2", 
                      "Endodermis_m1", "Endodermis_m2",
                      "PPP_m1","PPP_m2", "XPP_m1", "XPP_m2", 
                      "Cortex_m1", "Cortex_m2", 
                      "Atrichoblast_m1", "Trichoblast_m1", 
                      "Atrichoblast_m2", "Trichoblast_m2")

seurat_obj <- subset(
  seurat_obj, 
  idents = clusters_to_keep
)
DimPlot(seurat_obj)


seurat_obj <- RenameIdents(
  seurat_obj,
  `PPP_m1` = "Pericycle_m1",
  `PPP_m2` = "Pericycle_m2",
  `XPP_m1` = "Pericycle_m1",
  `XPP_m2` = "Pericycle_m2",
  `LRCEI` = "Initials",
  `SI` = "Initials",
  `CEI` = "Initials"
)
seurat_obj$root_tip <- Idents(seurat_obj)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = NULL)  # e.g., c("percent.mt","nCount_RNA")
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
dims_use <- 1:30 

seurat_obj <- FindNeighbors(seurat_obj, dims = dims_use)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)  # try a sweep below
seurat_obj <- RunUMAP(seurat_obj, dims = dims_use)        # or RunTSNE(seurat_obj, dims = dims_use)



# Quick look
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + NoLegend()

 
 AuxinPanel <- LoadSensors(setName = 'Auxin', species = 'ATH', hormone = 'aux',
                           customPanels = FALSE,
                           randomInfo = list('n' = 3, 'sizes' = c(100, 200, 300), 
                                             majortrend = TRUE))
 
 iSensors_obj <- CalcSensors(
   seurat_obj,
   seurLayer = 'data',
   panelSet = AuxinPanel,
   signals = c("mean"))
 
 Assays(iSensors_obj)
 
 Idents(iSensors_obj) <- "root_tip"
 
 saveRDS(iSensors_obj, file = "D:/!GitHub/DigitalSensor-Toolbox/iSensors-supplementary/00-iSensors-objects/data/shahan_iSensors_obj_groundtruth.rds")

 
 
# shahan epidermis
epidermis_obj <- readRDS(
   "00-iSensors-objects/data/shahan_epidermis.rds"
 )
 DimPlot(epidermis_obj)
 AuxinPanel <- LoadSensors(setName    = "Auxin",
                             species    = "ATH",
                             hormone    = "aux",
                             customPanels = FALSE,
                             randomInfo = list(n = 3, sizes = c(100, 200, 300),
                                               majortrend = TRUE))
 epidermis_obj <- CalcSensors(
     epidermis_obj,
     seurLayer = "data",
     panelSet  = AuxinPanel,
     signals   = "mean"
   )
 
 Assays(epidermis_obj) 
 epidermis_obj[["integrated"]] <- NULL
 epidermis_obj[["ALT"]] <- NULL
 saveRDS(epidermis_obj, file = "D:/!GitHub/DigitalSensor-Toolbox/iSensors-supplementary/00-iSensors-objects/data/shahan_epidermis_iSensors_obj.rds")
 
#subsetting and normalizing shahan root tip data ----

seurat_obj <- readRDS("D:/!GitHub/DigitalSensor-Toolbox/iSensors-supplementary/00-iSensors-objects/data/02_ra_integrated_Shahan_20250114.rds")
Assays(seurat_obj)
# 
seurat_obj[["integrated"]] <- NULL
seurat_obj[["ALT"]] <- NULL

meta <- seurat_obj@meta.data
levels(meta$Atlas_new)


DimPlot(seurat_obj)
#reducing the object for Spearmann analysis
clusters_to_keep <- c("QC", "CSC", "LRP", "LRCEI", "SI", "CEI", 
                      "Columella", "CSCD", "Young LRC", "C1", "LRC", "Dying LRC", 
                      "Protoxylem_m1", "Protoxylem_m2", "Protoxylem_t",
                      "PSE_m1", "PSE_m2", "PSE_t",
                      "PCC_m1", "PCC_m2", "PCC_t",
                      "Metaxylem_m1", "Metaxylem_m2", "Metaxylem_t", 
                      "Procambium_m1", "Procambium_m2", "Procambium_t",
                      "Endodermis_m1", "Endodermis_m2", "Endodermis_t",
                      "PPP_m1","PPP_m2", "XPP_m1", "XPP_m2", "XPP_t", "PPP_t",
                      "Cortex_m1", "Cortex_m2", "Cortex_t",
                      "Atrichoblast_m1", "Trichoblast_m1", "Trichoblast_t", 
                      "Atrichoblast_m2", "Trichoblast_m2", "Atrichoblast_t")

seurat_obj <- subset(
  seurat_obj, 
  idents = clusters_to_keep
)
DimPlot(seurat_obj)


seurat_obj$root_tip <- Idents(seurat_obj)

DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = NULL)  # e.g., c("percent.mt","nCount_RNA")
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
dims_use <- 1:30 

seurat_obj <- FindNeighbors(seurat_obj, dims = dims_use)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)  # try a sweep below
seurat_obj <- RunUMAP(seurat_obj, dims = dims_use)        # or RunTSNE(seurat_obj, dims = dims_use)

# Quick look
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + NoLegend()


Idents(seurat_obj) <- "root_tip"

saveRDS(seurat_obj, file = "D:/!GitHub/DigitalSensor-Toolbox/iSensors-supplementary/00-iSensors-objects/data/shahan_seurat_obj_root_tip.rds")
Assays(seurat_obj)

## calculating iSensors object

seurat_obj <- readRDS("D:/!GitHub/DigitalSensor-Toolbox/iSensors-supplementary/00-iSensors-objects/data/shahan_seurat_obj_root_tip.rds")

AuxinPanel <- LoadSensors(setName = 'Auxin', species = 'ATH', hormone = 'aux',
                                customPanels = FALSE,
                                randomInfo = list('n' = 3, 'sizes' = c(100, 200, 300), 
                                                  majortrend = TRUE))

iSensors_obj <- CalcSensors(
  seurat_obj,
  seurLayer = 'data',
  panelSet = AuxinPanel,
  signals = c("mean"))

Assays(iSensors_obj)




saveRDS(iSensors_obj, file = "D:/!GitHub/DigitalSensor-Toolbox/iSensors-supplementary/00-iSensors-objects/data/shahan-roottip-iSensors-obj.rds")

# test iSensors object ----

iSensors_obj <-readRDS(file = "D:/!GitHub/DigitalSensor-Toolbox/iSensors-supplementary/00-iSensors-objects/data/shahan-roottip-iSensors-obj.rds")


DefaultAssay(iSensors_obj) <- "iSensors_mean"



feat = "ATH-aux-trans-ARF"
feat = "ATH-aux-trans-PolarAuxinTransport"

vals <- FetchData(iSensors_obj, vars = feat)[, 1]
lim  <- max(abs(quantile(vals, probs = c(0.02, 0.98), na.rm = TRUE)))
lims <- c(0, lim)

p <- FeaturePlot(
  object = iSensors_obj,
  features = feat,
  combine = TRUE
)
library(RColorBrewer)
rdylbu5 <- rev(brewer.pal(n = 5, name = "RdYlBu"))

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


library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

DimPlot(iSensors_obj)
assay_use <- "iSensors_mean"

DefaultAssay(iSensors_obj) <- assay_use

# Все iSensors из assay
isensors_all <- rownames(iSensors_obj[[assay_use]])

# Среднее значение iSensor по кластерам
avg_list <- AverageExpression(
  iSensors_obj,
  assays = assay_use,
  features = isensors_all,
#  group.by = cluster_col,
  slot = "data",
  verbose = FALSE
)

avg_mat <- avg_list[[assay_use]]

# Опционально: z-score по каждому iSensor, чтобы сравнивать паттерны между кластерами
avg_z <- t(scale(t(avg_mat)))
avg_z[is.na(avg_z)] <- 0
avg_z <- avg_mat

# Long format
hm_df <- as.data.frame(avg_z) %>%
  tibble::rownames_to_column("iSensor") %>%
  pivot_longer(
    -iSensor,
    names_to = "cluster",
    values_to = "value_z"
  )

# Heatmap
p_heat <- ggplot(hm_df, aes(x = cluster, y = iSensor, fill = value_z)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    oob = squish,
    name = "iSensors value"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 6),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    legend.key.height = unit(0.35, "cm")
  )

p_heat



range(avg_mat, na.rm = TRUE)
summary(as.numeric(avg_mat))

#Calculation iSensors for small shahan dataset ----

install.packages("devtools")
devtools::install_github("MironovaLab/iSensors")
devtools::install_github("MironovaLab/iSensors", ref = "iSensors-v1.2.0-dev")
library(iSensors)
library(Seurat)
library(ggplot2)
library(scales)
library(Seurat)
library(patchwork)
library(cowplot)
set.seed(100)

#subsetting and normalizing shahan data ----

seurat_obj <- readRDS("D:/MyPAPERS/iSensors/Data/02_ra_integrated_Shahan_20250114.rds")

meta <- seurat_obj@meta.data

DimPlot(seurat_obj)

#reducing the object for Spearmann analysis
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

saveRDS(seurat_obj_spearman, file = "D:/MyPAPERS/iSensors/Data/shahan_obj_spearman.rds")

## calculating iSensors object

seurat_obj <- readRDS("D:/MyPAPERS/iSensors/Data/shahan_obj_spearman.rds")
seurat_obj[["integrated"]] <- NULL
seurat_obj[["ALT"]] <- NULL

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




saveRDS(iSensors_obj, file = "D:/!GitHub/iSensors-supplementary/00-iSensors-objects/data/shahan-iSensors-obj.rds")

# test iSensors object ----

iSensors_obj <-readRDS(file = "D:/!GitHub/iSensors-supplementary/00-iSensors-objects/data/shahan-iSensors-obj.rds")


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

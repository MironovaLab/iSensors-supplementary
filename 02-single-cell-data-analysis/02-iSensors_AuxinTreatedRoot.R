#THis script is to explore the original Auxin dataset and to run iSEnsors on it
library(Seurat)
library(iSensors)
library(tidyverse)
library(RColorBrewer)
library(scales)


rdylbu5 <- rev(brewer.pal(n = 5, name = "RdYlBu"))


#The datasets is found here: https:/www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241573

load("C:/!Victoria/!Projects/DigitalSensors/Data/Martin-Arevalillo-2025/GSE241573_Seurat_Object_Ath_Root_Auxin_treatment_DR5_DR15_IR8_ER13.RData")
ls()
load("E:/Users/vmironova/OneDrive - Radboud Universiteit/Datasets/scTranscriptomics/Martin-Arevalillo-2025/GSE241573_Seurat_Object_Ath_Root_Auxin_treatment_DR5_DR15_IR8_ER13.RData")

seurat_obj <-
  UpdateSeuratObject(Seurat_Object_Ath_Root_Auxin_treatment_DR5_DR15_IR8_ER13)

saveRDS(seurat_obj, file = "02-single-cell-data-analysis/in/MartinArevalillo2025.rds")

#seurat_obj <- readRDS(file = "02-single-cell-data-analysis/in/MartinArevalillo2025.rds")

DimPlot(seurat_obj)

seurat_obj$orig.ident2 <- recode(
  seurat_obj$orig.ident2,
  "Ctr" = "Control",
  "Aux" = "Auxin"
)

# Optional: enforce order
seurat_obj$orig.ident2 <- factor(
  seurat_obj$orig.ident2,
  levels = c("Control", "Auxin")
)


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

#This is the publication ready figure for Auxin data
ggsave("out/auxin/AuxinTreated_DimPlot_split.pdf", last_plot(), width = 8, height = 4.0, dpi = 300)

meta <- seurat_obj@meta.data
set.seed(100)
#meta <- iSensors_obj@meta.data

#Create iSensors object ----

#Loading custom panels
customTransPanel <- LoadSensors(setName = 'Auxin', species = 'ATH', hormone = 'aux',
                                customPanels = FALSE,
                                randomInfo = list('n' = 3, 'sizes' = c(100, 200, 300), 
                                majortrend = TRUE))

#With meta-panels
# customTransPanel <- LoadSensors(setName = 'Auxin', species = 'AT', hormone = 'aux',
#                                 customPanels = FALSE,
#                                 randomInfo = list('n' = 3, 'sizes' = c(100, 200, 300), 
#                                                   majortrend = TRUE),
#                                 metaPanels = list('Response1' = list('srcPanels' = c('AT_aux_trans_A_ARF', 'AT_aux_trans_IAA'),
#                                                                      'rule' = prod),
#                                                   'Response2' = list('srcPanels' = c('AT_aux_trans_ARF', 'AT_aux_trans_IAA'),
#                                                                      'rule' = prod),
#                                                   'Response3' = list('srcPanels' = c('AT_aux_trans_A_ARF', 'AT_aux_trans_IAA'),
#                                                                      'rule' = mean),
#                                                   'Response4' = list('srcPanels' = c('AT_aux_trans_ARF', 'AT_aux_trans_IAA'),
#                                                                      'rule' = mean))
# )
iSensors_obj <- CalcSensors(
  seurat_obj,
  seurLayer = 'data',
  panelSet = customTransPanel,
  signals = c("mean", "mean_normed")
)

DimPlot(
  iSensors_obj,
  reduction = "umap",
  split.by = "orig.ident",
  combine = TRUE
)

DefaultAssay(iSensors_obj) <-"iSensors_mean"
#DefaultAssay(iSensors_obj) <-"iSensors_mean_normed"

isensors_list <- unlist(iSensors_obj@assays$iSensors_mean@counts@Dimnames[1])
isensors_list

saveRDS(iSensors_obj, file = "02-single-cell-data-analysis/in/iSensors-Martin-Arevalillo-auxin-root2025_mean.rds")

isensors_list
#Summary all the plots
pdf("02-single-cell-data-analysis/out/iSensors_Auxin_FeaturePlots_split_mean.pdf", width = 12, height = 6)

getwd()

for (feat in isensors_list) {
  
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
      oob = squish
    ) &
    guides(color = guide_colorbar(title = "iSensor")) &
    theme(legend.position = "right")
  
  
  print(p)
  
}

dev.off()

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



# Compute unified limits from ALL cells (both conditions)

saveRDS(iSensors_obj, file = "C:/!Victoria/!Projects/DigitalSensors/Data/Martin-Arevalillo-2025/iSensors-Martin-Arevalillo-auxin-root2025.rds")
saveRDS(iSensors_obj, file = "D:/!GitHub/iSensors-supplementary/02-single-cell-data-analysis/iSensors-Martin-Arevalillo-auxin-root2025.rds")

# check if the data contains mTQ reporter expression


genes<- unlist(iSensors_obj@assays$RNA@counts@Dimnames[1])
genes

any(!grepl("^AT[1-5]G[0-9]{5}$", genes))

sensors<- genes[(!grepl("^AT[1-5]G[0-9]{5}$", genes))]
sensors
#Seurat object does not contain reporter genes

library(Matrix)

mat <- readMM("C:/!Victoria/!Projects/DigitalSensors/Data/Martin-Arevalillo-2025/GSE241573_Expression_Matrix_Normed_Counts_Ath_Root_Auxin_treatment_DR5_DR15_IR8_ER13.mtx")
# it alsodoes not contain mTQ

#IAA30

DefaultAssay(iSensors_obj) <-"RNA"

p<- FeaturePlot(
  object = iSensors_obj,
  features = "AT3G62100",
  split.by = "orig.ident2",
  combine = TRUE)
p <- p &
  scale_color_gradientn(
    colors = rdylbu5,
    limits = lims,
    oob = squish
  ) &
  guides(color = guide_colorbar(title = "iSensor")) &
  theme(legend.position = "right")
p
ggsave(plot = get_last_plot(), file = "C:/!Victoria/!Projects/DigitalSensors/Data/Martin-Arevalillo-2025/IAA30.jpg")

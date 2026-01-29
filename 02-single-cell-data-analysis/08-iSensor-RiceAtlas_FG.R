#Description: Rice atlas of high quality nuclei from eight organs (crouw root, culm, leaf, flag leaf, shoot apex, tiller bud, panicle and seed)
#Paper: https://www.nature.com/articles/s41586-025-09251-0
#Data (rds file): https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232863 

#Loading/installing required packages
library(Seurat)
library(Matrix)
library(dplyr)
#BiocManager::install('Rsamtools')
#install.packages("remotes")
#remotes::install_github("timoast/signac")
library(Signac)
#devtools::install_github("MironovaLab/iSensors")

#Inspecting dataset
Rice_Wang <- readRDS("X:/scRNA-seq datasets/Wang et al. 2025 Rice Atlas/GSE232863_scRNA_omics.Rds")
#Rice_Wang <- readRDS("//mironova.z.science.ru.nl/mironova/scRNA-seq datasets/Wang et al. 2025 Rice Atlas/GSE232863_scRNA_omics.Rds")
unique(Rice_Wang@meta.data$tissue_cluster)
DimPlot(Rice_Wang)
unique(Rice_Wang@meta.data$tissue)

#Running iSensor v1.2.0 on Rice data
setwd("//mironova.z.science.ru.nl/mironova/Fabienne/iSensor/DigitalSensor-Toolbox")
library(iSensors)
packageVersion("iSensors")
library(Seurat)

#Rice_Wang <- readRDS("//mironova.z.science.ru.nl/mironova/scRNA-seq datasets/Wang et al. 2025 Rice Atlas/GSE232863_scRNA_omics.Rds")
seurat_obj <- Rice_Wang
DimPlot(seurat_obj) #Check if data was loaded properly

#iSensorsTransPanelCreate 
library(stringr)
library(magrittr)
library(dplyr)
library(Biostrings)

#Make custom panels --> Already included in newest version iSensors
#iSensorsTransPanelCreate(panel_name = "OS_aux_trans_ARF", gene_list = "OSgenes_PLAZA/OS_aux_trans_ARF.txt", species = "OS", 
#                         panel_description = "ARF Trans pannel rice")
#iSensorsTransPanelCreate(panel_name = "OS_aux_trans_IAA", gene_list = "OSgenes_PLAZA/OS_aux_trans_IAA.txt", species = "OS", 
#                         panel_description = "IAA Trans pannel rice")
#iSensorsTransPanelCreate(panel_name = "OS_aux_trans_PAT", gene_list = "OSgenes_PLAZA/OS_aux_trans_PAT.txt", species = "OS", 
#                         panel_description = "PAT Trans pannel rice")
#iSensorsTransPanelCreate(panel_name = "OS_aux_trans_receptors", gene_list = "OSgenes_PLAZA/OS_aux_trans_receptors.txt", species = #"OS", 
#                        panel_description = "Receptors Trans pannel rice")
#iSensorsTransPanelCreate(panel_name = "OS_aux_trans_synthesis", gene_list = "OSgenes_PLAZA/OS_aux_trans_synthesis.txt", species = "OS",
#                         panel_description = "Synthesis Trans pannel rice")
#iSensorsTransPanelCreate(panel_name = "OS_aux_trans_transport", gene_list = "OSgenes_PLAZA/OS_aux_trans_transport.txt", species = "OS", 
#                         panel_description = "Transport Trans pannel rice")

#Copied conjugation pannel from develop branch in iSensors folder from develop branch

RicePanel <- LoadSensors(setName = 'RicePanelSet',
                         species = "OSA",
                         customPanels = TRUE, #Already loaded in newest iSensor package
                         randomInfo = list('n' = 3, 'sizes' = c(100, 200, 300), majortrend = TRUE),
                         metaPanels = list(
                           'meta1' = list('srcPanels' = c("AT_aux_cis_DR5_ARF1", "AT_aux_cistrans_DR5_ARF5_2_up"), rule = mean),
                           'meta2' = list('srcPanels' = c("AT_aux_cis_DR5_ARF1", "AT_aux_cistrans_DR5_ARF5_2_up"), rule = prod))
)
RicePanel #Check before running

result <- CalcSensors(seurat_obj,
                      panelSet = RicePanel,
                      signals = c("mean_normed", "median"))
saveRDS(result, file = "RiceiSensorResult_Conjugation_290126.rds")

#Average expression Heatmap
library(pheatmap)
Idents(result) <- "tissue_cluster"
avg_expr <- AverageExpression(object = result, assays = "iSensors_mean_normed", layer = "counts")
dense_mat <- as.matrix(avg_expr$"iSensors_mean_normed")
dense_mat[is.nan(dense_mat)] <- 0 #Remove NAN values
colnames(dense_mat)
write.csv(dense_mat, file = "Rice_iSensors_mean_normed_290126.csv", row.names = TRUE)
#saveRDS(dense_mat, file = "Rice_iSensor_mean_normed_matrix.rds")

#type_order <- c("Epidermis", "epidermis-near-root-hair", "Exodermis", "Sclerenchyma", "Cortex", "Endodermis", "Vascular-cylinder", "Meristem", "Root-cap")
#dense_mat <- dense_mat[,type_order]

pheatmap(dense_mat)
pheatmap(dense_mat, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE)

#Running iSensor via custom panels (Old pipeline)
#Panels made using PLAZA Homologous gene family, based on panels OSA
library(iSensors)
packageVersion('iSensors') 

seurat_obj <- Rice_Wang
DimPlot(seurat_obj) #Check if data was loaded properly

RicePanelSet <- read_genePanels(panelsDir = 'OSgenes_PLAZA', species = 'OS', hormone = 'aux', type = 'trans', format = 'txt')
names(RicePanelSet)

#seurat_obj <- iSensor_pipeline(seuratObject = seurat_obj, usePanelPreset = RicePanelSet, presetPanelName = 'RicePanelSet', metaPanels = NULL)  
seurat_obj <- iSensor_pipeline(seuratObject = seurat_obj, usePanelPreset = RicePanelSet, presetPanelName = 'RicePanelSet', 
                               randPanels = 1, randSize = c(5), 
                               metaPanels = list('Response2' = list('srcPanels' = c('OS_aux_trans_ARF', 'OS_aux_trans_IAA'), 'rule' = prod),
                                                 'Response4' = list('srcPanels' = c('OS_aux_trans_ARF', 'OS_aux_trans_IAA'), 'rule' = mean)))
#Could not reproduce response and response3 Carmen, as there is no separate A-ARF panel for rice

seurat_obj@assays$iSensor_mean_normed@counts
saveRDS(seurat_obj, file = "RiceAtlas_Wang_iSensor_metapanels.rds")

#seurat_obj <- readRDS("//mironova.z.science.ru.nl/mironova/Fabienne/iSensor/DigitalSensor-Toolbox/RiceAtlas_Wang_iSensor.rds")
DefaultAssay(seurat_obj) <- 'iSensor_mean_normed'
FeaturePlot(object = seurat_obj, features = 'OS-aux-trans-synthesis')
all_panels <- rownames(seurat_obj@assays$iSensor_mean_normed)
DoHeatmap(object = seurat_obj, features = all_panels, slot = 'counts')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Subset plots 
Roots <- subset(Rice_Wang, subset = tissue == "Root")
Roots@meta.data
#library(future)
#options(future.globals.maxSize = 4 *  1024^3)
#plan("multicore", workers = 4)
Roots <- SCTransform(Roots)
Roots <- RunPCA(Roots,verbose = FALSE)
Roots <- RunUMAP(Roots, dims = 1:50, verbose = FALSE)
Roots <- FindNeighbors(Roots, dims = 1:50, verbose = FALSE)
Roots <- FindClusters(Roots, resolution = 1,verbose = FALSE)
Idents(Roots) <- "cluster_names"
DimPlot(Roots)
#saveRDS(Roots, "RiceAtlas_Root.rds")

Roots <- readRDS("RiceAtlas_Root.rds")
Roots_subset <- subset(Roots, subset = cluster_names %in% c("Cortex", "Endodermis", "Epidermis", "Meristem", "Root_cap", "Vascular_cylinder"))
Idents(Roots_subset) <- "cluster_names"


Roots_subset <- RenameIdents(Roots_subset,
                             'Cortex' = "cortex",
                             "Endodermis" = "endodermis", 
                             "Epidermis" = "epidermis", 
                             "Meristem" = "meristem",
                             "Root_cap" = "root cap",
                             "Vascular_cylinder" = "vascular"
)
Roots_subset$CellTypes<-Idents(Roots_subset)

cluster_colors <- c(
  "cortex" = "#a6cee3",
  "endodermis" = "#1f78b4", 
  "epidermis" = "#b2df8a", 
  "meristem" = "#33a02c",
  "root cap" = "#fb9a99",
  "vascular" = "#e31a1c"
)
DimPlot(Roots_subset, cols = cluster_colors, label = FALSE, label.size = 5) 

#Leaf
Leaf <- subset(Rice_Wang, subset = tissue == "leaf")
Leaf@meta.data
#library(future)
#options(future.globals.maxSize = 4 *  1024^3)
#plan("multicore", workers = 4)
Leaf <- SCTransform(Leaf)
Leaf <- RunPCA(Leaf,verbose = FALSE)
Leaf <- RunUMAP(Leaf, dims = 1:50, verbose = FALSE)
Leaf <- FindNeighbors(Leaf, dims = 1:50, verbose = FALSE)
Leaf <- FindClusters(Leaf, resolution = 1,verbose = FALSE)
Idents(Leaf) <- "cluster_names"
DimPlot(Leaf)
#saveRDS(Leaf, "RiceAtlas_Leaf.rds")

Leaf <- readRDS("RiceAtlas_Leaf.rds")
Leaf_subset <- subset(Leaf, subset = cluster_names %in% c("Epidermis", "Mesophyll", "Guard_cell", "Vascular_cylinder"))
Idents(Leaf_subset) <- "cluster_names"

Leaf_subset <- RenameIdents(Leaf_subset,
                            "Epidermis" = "epidermis", 
                            "Mesophyll" = "mesophyll",
                            "Guard_cell" = "guard cell",
                            "Vascular_cylinder" = "vascular"
)
Leaf_subset$CellTypes<-Idents(Leaf_subset)

cluster_colors <- c(
  "cortex" = "#a6cee3",
  "endodermis" = "#1f78b4", 
  "epidermis" = "#b2df8a", 
  "meristem" = "#33a02c",
  "root cap" = "#fb9a99",
  "vascular" = "#e31a1c",
  "mesophyll" = "#fdbf6f", 
  "guard cell" = "#ff7f00"
)
DimPlot(Leaf_subset, cols = cluster_colors, label = FALSE, label.size = 5) 


#Stem
stem <- subset(Rice_Wang, subset = tissue == "ST")
stem@meta.data
#library(future)
#options(future.globals.maxSize = 4 *  1024^3)
#plan("multicore", workers = 4)
stem <- SCTransform(stem)
stem <- RunPCA(stem,verbose = FALSE)
stem <- RunUMAP(stem, dims = 1:50, verbose = FALSE)
stem <- FindNeighbors(stem, dims = 1:50, verbose = FALSE)
stem <- FindClusters(stem, resolution = 1,verbose = FALSE)
Idents(stem) <- "cluster_names"
DimPlot(stem)
#saveRDS(stem, "RiceAtlas_Stem.rds")

stem <- readRDS("RiceAtlas_Stem.rds")
stem_subset <- subset(stem, subset = cluster_names %in% c("Epidermis", "Cortex", "Vascular_cylinder"))
Idents(stem_subset) <- "cluster_names"
stem_subset <- RenameIdents(stem_subset,
                            "Epidermis" = "epidermis", 
                            "Cortex" = "cortex",
                            "Vascular_cylinder" = "vascular"
)
stem_subset$CellTypes<-Idents(stem_subset)

cluster_colors <- c(
  "cortex" = "#a6cee3",
  "endodermis" = "#1f78b4", 
  "epidermis" = "#b2df8a", 
  "meristem" = "#33a02c",
  "root cap" = "#fb9a99",
  "vascular" = "#e31a1c",
  "mesophyll" = "#fdbf6f", 
  "guard cell" = "#ff7f00"
)
DimPlot(stem_subset, cols = cluster_colors, label = FALSE, label.size = 5) 

#Flower
flower <- subset(Rice_Wang, subset = tissue == "SP")
flower@meta.data
#library(future)
#options(future.globals.maxSize = 4 *  1024^3)
#plan("multicore", workers = 4)
flower <- SCTransform(flower)
flower <- RunPCA(flower,verbose = FALSE)
flower <- RunUMAP(flower, dims = 1:50, verbose = FALSE)
flower <- FindNeighbors(flower, dims = 1:50, verbose = FALSE)
flower <- FindClusters(flower, resolution = 1,verbose = FALSE)
Idents(flower) <- "cluster_names"
DimPlot(flower)
#saveRDS(flower, "RiceAtlas_Flower.rds")

flower <- readRDS("RiceAtlas_Flower.rds")
flower_subset <- subset(flower, subset = cluster_names %in% c("Epidermis", "Tapetum", "Vascular_cylinder", "Pollen_mother_cell", "Ovary_wall", "Ovule", "Stigmata"))
Idents(flower_subset) <- "cluster_names"

flower_subset <- RenameIdents(flower_subset,
                              "Epidermis" = "epidermis", 
                              "Tapetum" = "tapetum",
                              "Vascular_cylinder" = "vascular", 
                              "Pollen_mother_cell" = "pollen",
                              "Ovary_wall" = "pistil",
                              "Ovule" = "pistil",
                              "Stigmata" = "pistil"
)
flower_subset$CellTypes<-Idents(flower_subset)

cluster_colors <- c(
  "cortex" = "#a6cee3",
  "endodermis" = "#1f78b4", 
  "epidermis" = "#b2df8a", 
  "meristem" = "#33a02c",
  "root cap" = "#fb9a99",
  "vascular" = "#e31a1c",
  "mesophyll" = "#fdbf6f", 
  "guard cell" = "#ff7f00",
  "tapetum" = "#cab2d6",
  "pollen" = "#6a3d9a",
  "pistil" = "#ffff9a"
)
DimPlot(flower_subset, cols = cluster_colors, label = FALSE, label.size = 5) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Root subset
library(iSensor)
packageVersion('iSensor') 

#seurat_obj <- readRDS("~/R PhD/Data/Wang 2025/GSE232863_scRNA_omics.Rds")
seurat_obj <- Rice_Wang
DimPlot(seurat_obj) #Check if data was loaded properly

RicePanelSet <- read_genePanels(panelsDir = 'OSgenes_PLAZA', species = 'OS', hormone = 'aux', type = 'trans', format = 'txt')
names(RicePanelSet)

seurat_obj_roots <- subset(seurat_obj, subset = tissue == "Root")
DimPlot(seurat_obj_roots, group.by = "cluster_names")

#seurat_obj_roots <- iSensor_pipeline(seuratObject = seurat_obj_roots, usePanelPreset = RicePanelSet, presetPanelName = 'RicePanelSet', randPanels = 1, randSize = c(5),
#                               metaPanels = list('R2D2' = list('srcPanels' = c('OS_aux_trans_ARF', 'OS_aux_trans_IAA'), 'rule' = prod)))
seurat_obj_roots <- readRDS("//mironova.z.science.ru.nl/mironova/Fabienne/iSensor/DigitalSensor-Toolbox/RiceAtlas_roots_Wang_iSensor.rds")

seurat_obj_roots@assays$iSensor_mean_normed@counts
all_panels <- rownames(seurat_obj_roots@assays$iSensor_mean_normed)
saveRDS(seurat_obj_roots, file = "RiceAtlas_roots_Wang_iSensor.rds")

DefaultAssay(seurat_obj_roots) <- 'iSensor_mean_normed'
FeaturePlot(object = seurat_obj_roots, features = 'OS-aux-trans-synthesis')
DoHeatmap(object = seurat_obj_roots, features = all_panels, slot = 'counts')

#Average expression Heatmap
library(pheatmap)
Idents(seurat_obj_roots) <- "cluster_names"
avg_expr <- AverageExpression(object = seurat_obj_roots, assays = "iSensor_mean_normed", layer = "counts")
dense_mat <- as.matrix(avg_expr$"iSensor_mean_normed")
dense_mat[is.nan(dense_mat)] <- 0 #Remove NAN values
colnames(dense_mat)

type_order <- c("Epidermis", "epidermis-near-root-hair", "Exodermis", "Sclerenchyma", "Cortex", "Endodermis", "Vascular-cylinder", "Meristem", "Root-cap")
dense_mat <- dense_mat[,type_order]

pheatmap(dense_mat)
pheatmap(dense_mat, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Combined visualisation Arabidopsis (Gui et al. 2025) and Rice (Wang et al. 2025) iSensor 
library(dplyr)
library(tidyr)
library(tibble)
library(pheatmap)

Rice <- readRDS("//mironova.z.science.ru.nl/mironova/Fabienne/iSensor/DigitalSensor-Toolbox/Rice_iSensor_mean_normed_matrix.rds")
pheatmap(Rice, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE)


Arabidopsis <- read.csv("//mironova.z.science.ru.nl/mironova/Fabienne/iSensor/DigitalSensor-Toolbox/Guo_results_select_panels.csv")
Arabidopsis <- Arabidopsis %>%
  mutate(tissue_organ = paste0(tissue, "_", organ))
Arabidopsis_pivot <- Arabidopsis %>%
  select(-organ, -tissue) %>%                     # remove unwanted columns
  pivot_longer(cols = -tissue_organ,              # all except 'tissue_organ'
               names_to = "panel", 
               values_to = "value") %>%           # long format
  pivot_wider(names_from = tissue_organ, 
              values_from = value)  
Arabidopsis_matrix <- Arabidopsis_pivot %>%
  column_to_rownames(var = "panel") %>%
  as.matrix()
pheatmap(Arabidopsis_matrix, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE)
 

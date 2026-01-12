##### Installs ComplexHeatmap #####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

##### Installs other packages #####

install.packages("devtools")
devtools::install_github("MironovaLab/iSensors")
install.packages("Seurat")
install.packages("circlize")

##### Setting the environment #####

library(Seurat)
library(iSensors)
#library(RColorBrewer)
#library(gplots)
library(ComplexHeatmap)
library(circlize)
library(forcats)

##### Creating iSensors signal matrix #####

in_path <- "C://Users/lenna/Documents/iSensors_test/"
seurat_obj <- readRDS(paste0(in_path, "seurat_auxin_microarrays_no_cb.rds"))
seurat_obj <- iSensor_pipeline(seuratObject=seurat_obj, species = 'AT', hormone = 'aux')
#seurat_obj <- iSensor_pipeline(seuratObject=seurat_obj, species = 'AT', hormone = 'aux', type = c('cis','trans'))
#seurat_obj <- iSensor_pipeline(seuratObject=seurat_obj, species = 'AT', hormone = 'aux', type = c('cistrans'))
#print(names(seurat_obj@assays))   # Check the assays
#seurat_obj@assays$iSensor_mean_normed@counts[1:10, 1:3]   #Check the data
matrix_data <- LayerData(seurat_obj, assay = "iSensor_mean_normed", layer = "counts")
matrix_data_new <- as.matrix(matrix_data)   # Transform data.frame to matrix

#write.table(matrix_data_new, file="matrix_data_new.txt") # saves in file, keeps the rownames
matrix_data_new <- as.matrix(read.table("matrix_data_new.txt",header=TRUE,row.names=1)) # says first column are rownames

##### Define breaks and colors #####

#mycols <- colorRamp2(breaks = c(min(matrix_data_new), 1, max(matrix_data_new)), 
#                     colors = c("#4575b4", "white", "#d73027"))
mycols <- colorRamp2(breaks = c(min(matrix_data_new), 1, max(matrix_data_new)), 
                     colors = c("#4575b4", "#ffffbf", "#d73027"))

##### Plot (advanced)

# Makes clustering and defines the order of rows and columns

hclust_obj <- Heatmap(matrix_data_new,
                      column_km = 2,    #split columns into 2
                      name = "S value", #title of legend
                      col = mycols,
                      column_title = "Microarrays", row_title = "Gene panels",
                      row_names_gp = gpar(fontsize = 7), # Text size for row names
                      column_names_gp = gpar(fontsize = 7) # Text size for row names
)
hclust_obj
column_order_list <- column_order(hclust_obj)
cluster1_columns <- colnames(matrix_data_new)[column_order_list[[1]]]
cluster2_columns <- colnames(matrix_data_new)[column_order_list[[2]]]
col_order_vector <- c(cluster2_columns, cluster1_columns)
row_order_vector <- row_order(hclust_obj)
row_order_names <- rownames(matrix_data_new)[row_order_vector]

# Prepares annotation 1 (control and auxin)

treatment <- substr(colnames(matrix_data_new), start = 1, stop = 5)

df_cols <- data.frame(
  microarray = colnames(matrix_data_new),
  treatment = treatment
)

col_1 = list(auxin = c("contr" = "green", "treat" = "darkred"))   # Define colors for each levels of qualitative variables
ha_1 <- HeatmapAnnotation(
  auxin = df_cols$treatment,
  col = col_1
#  height = unit(10, "cm") # Set the height of the entire column annotation (doesn't work in my hands)
)

# Removes undesired rows (ARFsyn)

rows_to_remove <- c(
  "AT-aux-cis-DR5-ARF5syn", 
  "AT-aux-cis-IR8-ARF5syn", 
  "AT-aux-cistrans-DR5-ARF5syn-up", 
  "AT-aux-cistrans-DR5-ARF5syn-down",
  "AT-aux-cistrans-IR8-ARF5syn-up",
  "AT-aux-cistrans-IR8-ARF5syn-down"
)
matrix_data_short <- matrix_data_new[!rownames(matrix_data_new) %in% rows_to_remove, ]
row_order_names_short <- row_order_names[!row_order_names %in% rows_to_remove]

# Prepares annotation 2 (panel type)

Panel_name <- rownames(matrix_data_short)
Panel_type <- c()
i = 1
for (i in (1:nrow(matrix_data_short))) {
  if (grepl("-cis-", Panel_name[i], ignore.case = TRUE)) {
    Panel_type <- c(Panel_type, 'Cis')
  } else {
    if (grepl("-trans-", Panel_name[i], ignore.case = TRUE)) {
      Panel_type <- c(Panel_type, 'Trans')
    } else {
      if (grepl("cistrans", Panel_name[i], ignore.case = TRUE)) {
        Panel_type <- c(Panel_type, 'CisTrans')
      } else {
        Panel_type <- c(Panel_type, 'Neg control')
      }
    }
  }
  i = i + 1
}

col_2 = list(Panel_type = c("Cis" = "#fdbf6f", "Trans" = "#ff7f00", "CisTrans" = "#b15928", "Neg control" = "gray"))   # Define colors for each levels of qualitative variables
ha_2 <- rowAnnotation(
  Panel_type = Panel_type,
  col = col_2
)

# Plot

png(
  filename = "heatmap_short.png",
  width = 21, 
  height = 30, 
  units = "cm",
  res = 300   # Resolution (the size should be in physical units)
)
Heatmap(matrix_data_short,
        top_annotation = ha_1,
        right_annotation = ha_2,
        column_order = col_order_vector,
        row_order = row_order_names_short,
        name = "S value", #title of legend
        col = mycols,
        column_title = "Microarrays", row_title = "Gene panels",
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        column_names_gp = gpar(fontsize = 7) # Text size for row names
)
dev.off()


#This script uses the seurat object loaded in iSensors-tutorial-test-load.R script
#run it after the script or upload the saved version:
#seu <- readRDS("path/to/your_scRNAseq_seurat.rds")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(iSensors)
  library(pheatmap)
  library(RColorBrewer)
  library(tibble)
  
})


DimPlot(
  seu,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  split.by = "orig.ident2"
) + ggtitle("EarlyMeristem: control vs auxin")

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

p_feat <- FeaturePlot(
  iSensors_obj,
  features = c("AT-aux-trans-Transport", "AT-aux-trans-PolarAuxinTransport", "AT-aux-trans-ARF", "AT-aux-trans-Synthesis", "AT-aux-trans-IAA", "AT-aux-trans-ConjugationDeconjugation"),
label = TRUE)
p_feat

p_feat <- FeaturePlot(
  seu,
  features = c("AT3G11260", "AT5G66770"), split.by = "orig.ident2",
  label = TRUE)
p_feat



# OrRd (n=5) palette
expr_cols <- brewer.pal(n = 3, name = "OrRd")

p <- FeaturePlot(
  iSensors_obj,
  features = c("AT-aux-trans-Transport","AT-aux-trans-ARF", "AT-aux-trans-Synthesis","AT-aux-trans-ConjugationDeconjugation"),
  reduction = "umap",
  cols = expr_cols,
  min.cutoff = 0,
  max.cutoff = "q97", 
  ncol = 2, 
  order = TRUE,
#  label = TRUE
) &
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank()
  ) &
  coord_equal()
p

# Save at 300 dpi
ggsave(
  filename = "Tutorial-iSensors/out/FeaturePlot_iSensors_4panels_pattern.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)

p <- FeaturePlot(
  iSensors_obj,
  features = c("AT-aux-trans-ARF", "AT-aux-trans-Synthesis"),
  reduction = "umap",
  cols = expr_cols,
  min.cutoff = 0,
  max.cutoff = "q97", 
  ncol = 2, 
  order = TRUE,
  split.by = "orig.ident2"
  #  label = TRUE
) &
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank()
  ) &
  coord_equal()
p

# Save at 300 dpi
ggsave(
  filename = "Tutorial-iSensors/out/FeaturePlot_iSensors_2panels_auxinresponse.png",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300
)


sensors <- c(
  "AT-aux-trans-Transport",
  "AT-aux-trans-ARF",
  "AT-aux-trans-Synthesis",
  "AT-aux-trans-ConjugationDeconjugation"
)


# Fetch sensor scores + metadata
df <- FetchData(
  iSensors_obj,
  vars = c(sensors, "cluster_annot", "orig.ident2")
)

# Average per cell type per condition
avg_df <- df %>%
  group_by(orig.ident2, cluster_annot) %>%
  summarise(across(all_of(sensors), mean, na.rm = TRUE), .groups = "drop")

# Control matrix
mat_ctr <- avg_df %>%
  dplyr::filter(orig.ident2 == "Ctr") %>%
  dplyr::select(cluster_annot, all_of(sensors)) %>%
  as.data.frame() %>%                     # <- important
  column_to_rownames("cluster_annot") %>%
  t()

# Auxin matrix
mat_aux <- avg_df %>%
  dplyr::filter(orig.ident2 == "Aux") %>%
  dplyr::select(cluster_annot, all_of(sensors)) %>%
  as.data.frame() %>%                     # <- important
  column_to_rownames("cluster_annot") %>%
  t()

dim(mat_ctr)
dim(mat_aux)

rownames(mat_ctr)  # should be sensors
colnames(mat_ctr)  # should be cluster_annot

celltype_order <- c("Daughters_LRCEI", "Daughters_CEI", "TransitInitials", "Initials")

mat_ctr <- mat_ctr[, intersect(celltype_order, colnames(mat_ctr)), drop = FALSE]
mat_aux <- mat_aux[, intersect(celltype_order, colnames(mat_aux)), drop = FALSE]


# Ensure same column order
common_ct <- intersect(colnames(mat_ctr), colnames(mat_aux))

# Rename columns to be unique and encode condition
colnames(mat_ctr) <- paste0(common_ct, "__Ctr")
colnames(mat_aux) <- paste0(common_ct, "__Aux")

# Combine (Ctr | Aux)
mat <- cbind(mat_ctr, mat_aux)

# Column annotation
annotation_col <- data.frame(
  Condition = rep(c("Ctr", "Aux"), each = length(common_ct)),
  row.names = colnames(mat)
)

rownames(annotation_col) <- colnames(mat)
colnames(mat)

cols <- brewer.pal(n = 5, name = "OrRd")

labels_col <- sub("__Ctr$|__Aux$", "", colnames(mat))





mat <- cbind(mat_ctr, mat_aux)
mat_z <- row_z(mat)
annotation_col <- data.frame(
  Condition = rep(c("Ctr", "Aux"), each = length(common_ct)),
  row.names = colnames(mat_z)
)


labels_col <- sub("__Ctr$|__Aux$", "", colnames(mat_z))


ann_colors <- list(
  Condition = c(
    Ctr = "#8DA0A8",   # muted blue-gray
    Aux = "#E6A57E"    # soft orange
  )
)



png(
  filename = "Tutorial-iSensors/out/Pheatmap_iSensors_4panels_aux_vs_ctr.png",
  width = 1800,
  height = 900,
  res = 300,
  bg = "white"   # <-- critical
)

pheatmap(
  mat_z,
  color = colorRampPalette(brewer.pal(5, "OrRd"))(100),
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  annotation_legend = FALSE,   # remove condition legend
  labels_col = labels_col,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_row = 6,
  fontsize_col = 6,
  main = "iSensors activity (ctr vs aux)"
)

dev.off()


#Differentially expressed iSensors

DefaultAssay(iSensors_obj) <- "iSensors_mean_normed"

DEG1 <- FindAllMarkers(iSensors_obj,group.by = "orig.ident2")
DEG3 <- FindAllMarkers(iSensors_obj,group.by = c("cluster_annot", "orig.ident2"))



sensors <- c(
  "AT-aux-trans-Transport",
  "AT-aux-trans-ARF",
  "AT-aux-trans-Synthesis",
  "AT-aux-trans-ConjugationDeconjugation"
)

res_ct <- lapply(sort(unique(iSensors_obj$cluster_annot)), function(ct) {
  obj_ct <- subset(iSensors_obj, subset = cluster_annot == ct)
  
  FindMarkers(
    object = obj_ct,
    ident.1 = "Aux",
    ident.2 = "Ctr",
    group.by = "orig.ident2",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0
  ) %>%
    rownames_to_column("iSensor") %>%
    mutate(celltype = ct)
}) %>% bind_rows()

res_ct
write.csv(res_ct, file =  "Tutorial-iSensors/out/iSensors-diffexpression.csv")

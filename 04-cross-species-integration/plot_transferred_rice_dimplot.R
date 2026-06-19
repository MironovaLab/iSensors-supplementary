#!/usr/bin/env Rscript

# Build RNA/ortholog-based UMAPs for the Wang rice root object annotated by
# Arabidopsis Atlas_new label transfer.
#
# Run from iSensors-supplementary/:
# Rscript 04-cross-species-integration/plot_transferred_rice_dimplot.R

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(100)

input_rds <- paste0(
  "00-iSensors-objects/data/arabidopsis_rice_transfer/",
  "wang_rice_root_label_transferred.rds"
)
output_dir <- "00-iSensors-objects/data/arabidopsis_rice_transfer/dimplots"
output_rds <- paste0(
  "00-iSensors-objects/data/arabidopsis_rice_transfer/",
  "wang_rice_root_label_transferred_umap.rds"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading transferred rice object")
rice <- readRDS(input_rds)
DefaultAssay(rice) <- "RNA"

message("Calculating RNA PCA and UMAP")
rice <- NormalizeData(rice, verbose = FALSE)
rice <- FindVariableFeatures(
  rice,
  selection.method = "vst",
  nfeatures = min(3000L, nrow(rice)),
  verbose = FALSE
)
rice <- ScaleData(
  rice,
  features = VariableFeatures(rice),
  verbose = FALSE
)
rice <- RunPCA(
  rice,
  features = VariableFeatures(rice),
  npcs = 30,
  verbose = FALSE
)
rice <- RunUMAP(
  rice,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap.ortholog",
  reduction.key = "orthUMAP_",
  seed.use = 100,
  verbose = FALSE
)

theme_dim <- theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.height = grid::unit(0.32, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

fine_labels <- sort(unique(as.character(rice$predicted.Atlas_new)))
fine_non_uncertain <- setdiff(fine_labels, "Uncertain")
fine_colors <- setNames(hue_pal(l = 62, c = 100)(length(fine_non_uncertain)),
                        fine_non_uncertain)
fine_colors <- c(fine_colors, Uncertain = "grey85")

broad_labels <- sort(unique(as.character(rice$predicted.Atlas_broad)))
broad_colors <- setNames(hue_pal(l = 58, c = 90)(length(broad_labels)),
                         broad_labels)

rice_labels <- sort(unique(as.character(rice$rice_annotation_original)))
rice_colors <- setNames(hue_pal(l = 60, c = 95)(length(rice_labels)),
                        rice_labels)

p_fine <- DimPlot(
  rice,
  reduction = "umap.ortholog",
  group.by = "predicted.Atlas_new",
  cols = fine_colors,
  raster = TRUE,
  raster.dpi = c(300, 300),
  pt.size = 0.12,
  shuffle = TRUE,
  seed = 100
) +
  ggtitle("Transferred Arabidopsis subtypes") +
  guides(colour = guide_legend(
    title = "Atlas_new",
    override.aes = list(size = 2, alpha = 1),
    ncol = 2
  )) +
  theme_dim

p_broad <- DimPlot(
  rice,
  reduction = "umap.ortholog",
  group.by = "predicted.Atlas_broad",
  cols = broad_colors,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  raster = TRUE,
  raster.dpi = c(300, 300),
  pt.size = 0.15,
  shuffle = TRUE,
  seed = 100
) +
  ggtitle("Transferred broad lineage") +
  guides(colour = guide_legend(
    title = "Broad lineage",
    override.aes = list(size = 2, alpha = 1)
  )) +
  theme_dim

p_original <- DimPlot(
  rice,
  reduction = "umap.ortholog",
  group.by = "rice_annotation_original",
  cols = rice_colors,
  label = TRUE,
  repel = TRUE,
  label.size = 3,
  raster = TRUE,
  raster.dpi = c(300, 300),
  pt.size = 0.15,
  shuffle = TRUE,
  seed = 100
) +
  ggtitle("Original Wang annotation") +
  guides(colour = guide_legend(
    title = "Wang cell type",
    override.aes = list(size = 2, alpha = 1)
  )) +
  theme_dim

p_comparison <- p_original + p_broad +
  plot_annotation(title = "Wang rice root: original and transferred annotations")

ggsave(
  file.path(output_dir, "rice_transferred_Atlas_new_dimplot.png"),
  p_fine,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(output_dir, "rice_transferred_Atlas_new_dimplot.pdf"),
  p_fine,
  width = 12,
  height = 8,
  device = cairo_pdf
)
ggsave(
  file.path(output_dir, "rice_transferred_broad_dimplot.png"),
  p_broad,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(output_dir, "rice_original_vs_transferred_dimplot.png"),
  p_comparison,
  width = 16,
  height = 7,
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(output_dir, "rice_original_vs_transferred_dimplot.pdf"),
  p_comparison,
  width = 16,
  height = 7,
  device = cairo_pdf
)

# Keep normalized data and reductions for immediate downstream plotting, but
# omit scale.data because it is straightforward to regenerate and adds size.
rice <- DietSeurat(
  rice,
  assays = "RNA",
  dimreducs = c("pca", "umap.ortholog"),
  graphs = NULL,
  layers = c("counts", "data"),
  misc = TRUE
)
saveRDS(rice, output_rds, compress = FALSE)

message("Saved plots in: ", output_dir)
message("Saved UMAP-enriched object: ", output_rds)

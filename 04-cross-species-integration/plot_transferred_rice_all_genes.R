#!/usr/bin/env Rscript

# Rebuild the Wang rice root object with its complete RNA gene set, calculate
# PCA/UMAP and de novo clusters, then visualize and summarize the transferred
# Arabidopsis annotations.
#
# Run from iSensors-supplementary/:
# Rscript 04-cross-species-integration/plot_transferred_rice_all_genes.R

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

set.seed(100)

source_rds <- "00-iSensors-objects/data/Wang_rice.Rds"
predictions_tsv <- paste0(
  "00-iSensors-objects/data/arabidopsis_rice_transfer/",
  "wang_rice_root_predictions.tsv"
)
output_dir <- paste0(
  "00-iSensors-objects/data/arabidopsis_rice_transfer/",
  "all_gene_umap"
)
output_rds <- paste0(
  "00-iSensors-objects/data/arabidopsis_rice_transfer/",
  "wang_rice_root_label_transferred_all_genes_umap.rds"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

log_step <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

rna_counts <- function(obj) {
  assay <- obj[["RNA"]]
  if (is.null(assay)) stop("Rice object has no RNA assay")
  layers <- Layers(assay)
  count_layers <- grep("^counts($|\\.)", layers, value = TRUE)
  if (!length(count_layers)) {
    stop("RNA assay has no counts layer: ", paste(layers, collapse = ", "))
  }
  if (length(count_layers) > 1L) {
    assay <- JoinLayers(assay, layers = count_layers, new = "counts")
  }
  LayerData(assay, layer = "counts")
}

named_palette <- function(labels, palette = "Dark 3") {
  labels <- sort(unique(as.character(labels)))
  setNames(grDevices::hcl.colors(length(labels), palette = palette), labels)
}

annotation_centroids <- function(df, label_col, min_cells = 20L) {
  dt <- as.data.table(df)
  dt[, .(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    cells = .N
  ), by = label_col][cells >= min_cells]
}

normalized_mutual_information <- function(x, y) {
  tab <- table(x, y)
  pxy <- tab / sum(tab)
  px <- rowSums(pxy)
  py <- colSums(pxy)
  nz <- which(pxy > 0, arr.ind = TRUE)
  mi <- sum(vapply(seq_len(nrow(nz)), function(i) {
    r <- nz[i, 1L]
    c <- nz[i, 2L]
    pxy[r, c] * log(pxy[r, c] / (px[r] * py[c]))
  }, numeric(1)))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 || hy == 0) return(NA_real_)
  mi / sqrt(hx * hy)
}

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  expected <- sum_rows * sum_cols / choose2(n)
  maximum <- (sum_rows + sum_cols) / 2
  if (maximum == expected) return(NA_real_)
  (sum_cells - expected) / (maximum - expected)
}

log_step("Loading full Wang rice atlas")
source <- readRDS(source_rds)
if (!"tissue" %in% colnames(source[[]])) stop("Missing metadata: tissue")

root_cells <- rownames(source[[]])[
  tolower(as.character(source$tissue)) == "root"
]
pred <- fread(predictions_tsv)
if (!setequal(root_cells, pred$cell)) {
  stop("Prediction cells do not exactly match rice root cells")
}

log_step("Extracting all RNA genes for ", length(root_cells), " root cells")
counts <- rna_counts(source)[, root_cells, drop = FALSE]
meta <- source[[]][root_cells, , drop = FALSE]
rm(source)
invisible(gc())

rice <- CreateSeuratObject(
  counts = counts,
  assay = "RNA",
  meta.data = meta,
  project = "Wang_rice_root_all_genes"
)
rm(counts, meta)
invisible(gc())

pred_meta <- as.data.frame(pred)
rownames(pred_meta) <- pred_meta$cell
pred_meta$cell <- NULL
rice <- AddMetaData(rice, pred_meta)

log_step("Running all-gene RNA preprocessing")
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
  npcs = 50,
  verbose = FALSE
)
rice <- RunUMAP(
  rice,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap.allgenes",
  reduction.key = "allGeneUMAP_",
  n.neighbors = 30,
  min.dist = 0.15,
  metric = "cosine",
  seed.use = 100,
  verbose = FALSE
)
rice <- FindNeighbors(rice, reduction = "pca", dims = 1:30, verbose = FALSE)
rice <- FindClusters(
  rice,
  resolution = 0.6,
  cluster.name = "all_gene_cluster",
  random.seed = 100,
  verbose = FALSE
)

embedding <- as.data.table(Embeddings(rice, "umap.allgenes"), keep.rownames = "cell")
setnames(embedding, c("allGeneUMAP_1", "allGeneUMAP_2"), c("UMAP_1", "UMAP_2"))
plot_data <- merge(
  embedding,
  as.data.table(rice[[]], keep.rownames = "cell"),
  by = "cell",
  sort = FALSE
)

old_colors <- c(
  Cortex = "#D55E00",
  Endodermis = "#E69F00",
  Epidermis = "#009E73",
  epidermis_near_root_hair = "#56B4E9",
  Exodermis = "#0072B2",
  Meristem = "#CC79A7",
  Root_cap = "#F0E442",
  Sclerenchyma = "#999999",
  Vascular_cylinder = "#6A3D9A"
)
broad_colors <- c(
  Cortex = "#D55E00",
  Endodermis = "#E69F00",
  Epidermis = "#009E73",
  Meristem_initials = "#CC79A7",
  Pericycle = "#56B4E9",
  Root_cap = "#F0E442",
  Vascular = "#6A3D9A"
)
cluster_colors <- named_palette(plot_data$all_gene_cluster, "Dark 3")
fine_colors <- named_palette(
  plot_data[predicted.Atlas_new != "Uncertain", predicted.Atlas_new],
  "Dynamic"
)

base_theme <- theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.height = grid::unit(0.35, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

point_layer <- function(data, color_col, size = 0.28, alpha = 0.82) {
  geom_point(
    data = data,
    aes(x = UMAP_1, y = UMAP_2, color = .data[[color_col]]),
    size = size,
    alpha = alpha,
    shape = 16,
    stroke = 0
  )
}

make_labeled_plot <- function(color_col, colors, title, legend_title) {
  centers <- annotation_centroids(plot_data, color_col)
  ggplot(plot_data, aes(UMAP_1, UMAP_2)) +
    point_layer(plot_data, color_col) +
    geom_label_repel(
      data = centers,
      aes(label = .data[[color_col]]),
      size = 3,
      label.size = 0.15,
      fill = alpha("white", 0.85),
      seed = 100,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_color_manual(values = colors, breaks = names(colors), drop = FALSE) +
    labs(title = title, color = legend_title) +
    guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    base_theme
}

p_clusters <- make_labeled_plot(
  "all_gene_cluster",
  cluster_colors,
  "De novo clusters from all rice genes",
  "RNA cluster"
)
p_old <- make_labeled_plot(
  "rice_annotation_original",
  old_colors,
  "Original Wang annotation",
  "Wang cell type"
)
p_broad <- make_labeled_plot(
  "predicted.Atlas_broad",
  broad_colors,
  "Transferred broad lineage",
  "Transferred lineage"
)

uncertain <- plot_data[predicted.Atlas_new == "Uncertain"]
confident <- plot_data[predicted.Atlas_new != "Uncertain"]
p_fine <- ggplot(plot_data, aes(UMAP_1, UMAP_2)) +
  geom_point(
    data = uncertain,
    color = "grey88",
    size = 0.14,
    alpha = 0.3,
    shape = 16,
    stroke = 0
  ) +
  point_layer(confident, "predicted.Atlas_new", size = 0.34, alpha = 0.95) +
  scale_color_manual(
    values = fine_colors,
    breaks = names(fine_colors),
    drop = FALSE
  ) +
  labs(
    title = "Confident transferred Arabidopsis subtypes",
    subtitle = "Uncertain cells are shown underneath in grey",
    color = "Atlas_new"
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 2.5, alpha = 1),
    ncol = 2
  )) +
  base_theme

p_fine_faceted <- p_fine +
  facet_wrap(~rice_broad_original, ncol = 3, scales = "free") +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

p_structure <- p_clusters + p_old +
  plot_annotation(
    title = "Wang rice root: all-gene RNA structure and original annotation"
  )
p_transfer <- p_broad + p_fine +
  plot_annotation(
    title = "Wang rice root: transferred Arabidopsis annotations"
  )

save_plot <- function(name, plot, width, height) {
  ggsave(
    file.path(output_dir, paste0(name, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 350,
    bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(name, ".pdf")),
    plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
}

save_plot("rice_all_genes_clusters_vs_original", p_structure, 16, 7)
save_plot("rice_all_genes_broad_vs_fine_transfer", p_transfer, 18, 8)
save_plot("rice_all_genes_transferred_fine_faceted", p_fine_faceted, 16, 13)
save_plot("rice_all_genes_de_novo_clusters", p_clusters, 9, 7)
save_plot("rice_all_genes_original_annotation", p_old, 9, 7)
save_plot("rice_all_genes_transferred_broad", p_broad, 9, 7)

# Summary statistics ---------------------------------------------------------
stats <- as.data.table(rice[[]], keep.rownames = "cell")
stats[, accepted := predicted.Atlas_new != "Uncertain"]

coverage_old <- stats[, .(
  cells = .N,
  accepted_cells = sum(accepted),
  accepted_percent = 100 * mean(accepted),
  median_score = median(prediction.score.max),
  median_margin = median(prediction.score.margin),
  raw_subtypes = uniqueN(predicted.Atlas_new_raw),
  accepted_subtypes = uniqueN(predicted.Atlas_new[accepted])
), by = rice_annotation_original][order(-cells)]

coverage_broad <- stats[, .(
  cells = .N,
  accepted_cells = sum(accepted),
  accepted_percent = 100 * mean(accepted),
  median_score = median(prediction.score.max),
  median_margin = median(prediction.score.margin)
), by = rice_broad_original][order(-cells)]

fine_summary <- stats[, .(
  raw_cells = .N,
  accepted_cells = sum(accepted),
  accepted_percent = 100 * mean(accepted),
  median_score = median(prediction.score.max),
  median_margin = median(prediction.score.margin)
), by = .(
  rice_annotation_original,
  predicted.Atlas_new_raw
)][order(rice_annotation_original, -raw_cells)]

old_by_new <- dcast(
  stats,
  rice_annotation_original ~ predicted.Atlas_new,
  value.var = "cell",
  fun.aggregate = length
)
old_by_raw <- dcast(
  stats,
  rice_annotation_original ~ predicted.Atlas_new_raw,
  value.var = "cell",
  fun.aggregate = length
)
old_broad_by_new_broad <- dcast(
  stats,
  rice_broad_original ~ predicted.Atlas_broad,
  value.var = "cell",
  fun.aggregate = length
)
cluster_by_old <- dcast(
  stats,
  all_gene_cluster ~ rice_annotation_original,
  value.var = "cell",
  fun.aggregate = length
)
cluster_by_new <- dcast(
  stats,
  all_gene_cluster ~ predicted.Atlas_new,
  value.var = "cell",
  fun.aggregate = length
)

accepted_stats <- stats[accepted == TRUE]
agreement_metrics <- data.table(
  comparison = c(
    "All-gene cluster vs original Wang annotation",
    "All-gene cluster vs accepted transferred fine annotation",
    "Original Wang annotation vs accepted transferred fine annotation"
  ),
  cells = c(nrow(stats), nrow(accepted_stats), nrow(accepted_stats)),
  adjusted_rand_index = c(
    adjusted_rand_index(stats$all_gene_cluster, stats$rice_annotation_original),
    adjusted_rand_index(
      accepted_stats$all_gene_cluster,
      accepted_stats$predicted.Atlas_new
    ),
    adjusted_rand_index(
      accepted_stats$rice_annotation_original,
      accepted_stats$predicted.Atlas_new
    )
  ),
  normalized_mutual_information = c(
    normalized_mutual_information(
      stats$all_gene_cluster,
      stats$rice_annotation_original
    ),
    normalized_mutual_information(
      accepted_stats$all_gene_cluster,
      accepted_stats$predicted.Atlas_new
    ),
    normalized_mutual_information(
      accepted_stats$rice_annotation_original,
      accepted_stats$predicted.Atlas_new
    )
  )
)

fwrite(coverage_old, file.path(output_dir, "annotation_coverage_by_Wang_type.tsv"),
       sep = "\t")
fwrite(coverage_broad, file.path(output_dir, "annotation_coverage_by_broad_type.tsv"),
       sep = "\t")
fwrite(fine_summary, file.path(output_dir, "fine_label_summary_by_Wang_type.tsv"),
       sep = "\t")
fwrite(old_by_new, file.path(output_dir, "Wang_vs_accepted_transferred_counts.tsv"),
       sep = "\t")
fwrite(old_by_raw, file.path(output_dir, "Wang_vs_raw_transferred_counts.tsv"),
       sep = "\t")
fwrite(
  old_broad_by_new_broad,
  file.path(output_dir, "old_vs_transferred_broad_counts.tsv"),
  sep = "\t"
)
fwrite(cluster_by_old, file.path(output_dir, "RNA_cluster_vs_Wang_counts.tsv"),
       sep = "\t")
fwrite(cluster_by_new, file.path(output_dir, "RNA_cluster_vs_transferred_counts.tsv"),
       sep = "\t")
fwrite(agreement_metrics, file.path(output_dir, "annotation_agreement_metrics.tsv"),
       sep = "\t")

summary_lines <- c(
  paste("Rice root cells:", ncol(rice)),
  paste("RNA genes:", nrow(rice)),
  paste("Variable genes used for PCA:", length(VariableFeatures(rice))),
  paste("All-gene RNA clusters:", uniqueN(stats$all_gene_cluster)),
  sprintf(
    "Confident fine-label assignments: %d (%.1f%%)",
    sum(stats$accepted),
    100 * mean(stats$accepted)
  ),
  "",
  "Important interpretation:",
  paste(
    "Transferred fine labels were inferred within Wang-derived broad",
    "lineages, so broad old/new agreement is partly imposed by design."
  ),
  paste(
    "ARI/NMI therefore describe association and cluster alignment,",
    "not independent validation of biological truth."
  )
)
writeLines(summary_lines, file.path(output_dir, "summary.txt"))

# Retain all genes, normalized data, PCA and UMAP, but remove scale.data to
# reduce the saved object size.
rice <- DietSeurat(
  rice,
  assays = "RNA",
  dimreducs = c("pca", "umap.allgenes"),
  graphs = NULL,
  layers = c("counts", "data"),
  misc = TRUE
)
saveRDS(rice, output_rds, compress = FALSE)

log_step("Saved plots and statistics to ", output_dir)
log_step("Saved all-gene UMAP object to ", output_rds)

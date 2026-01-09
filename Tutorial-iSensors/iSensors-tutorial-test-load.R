# Load minimal tutorial dataset and rebuild Seurat object
# Input file was created in iSensors-tutorial-testdata-praparation-clean.R:

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
})

# ----------------------------
# 1) Load the minimal object
# ----------------------------

x <- readRDS("Tutorial-iSensors/data/stemcellniche_minimal_github.rds")

# ----------------------------
# 2) Sanity checks and alignment
# ----------------------------

counts <- x$counts
meta   <- x$meta
umap   <- x$umap

# Ensure sparse matrix
if (!inherits(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

# Ensure metadata rownames are cell IDs
if (is.null(rownames(meta))) {
  stop("meta must have rownames corresponding to cell IDs")
}

# Ensure consistent cell ordering
cells <- intersect(colnames(counts), rownames(meta))
stopifnot(length(cells) > 0)

counts <- counts[, cells, drop = FALSE]
meta   <- meta[cells, , drop = FALSE]

# ----------------------------
# 3) Rebuild Seurat object
# ----------------------------

seu <- CreateSeuratObject(
  counts = counts,
  assay = "RNA",
  meta.data = meta
)

DefaultAssay(seu) <- "RNA"

# ----------------------------
# 4) Add UMAP embedding
# ----------------------------
umap <- umap[cells, , drop = FALSE]
colnames(umap) <- c("UMAP_1", "UMAP_2")

seu[["umap"]] <- CreateDimReducObject(
  embeddings = umap,
  key = "UMAP_",
  assay = "RNA"
)

# ----------------------------
# 5) Minimal preprocessing for plotting
# ----------------------------

# Required for FeaturePlot / DotPlot
seu <- NormalizeData(seu, verbose = FALSE)

# Optional (safe to skip if not needed)
# seu <- FindVariableFeatures(seu, verbose = FALSE)

# ----------------------------
# 6) Example plot (sanity check)
# ----------------------------

DimPlot(
  seu,
  reduction = "umap",
  group.by = "cluster_annot",
  label = TRUE,
  repel = TRUE,
  split.by = "orig.ident2"
) + ggtitle("Stem cell niche: control vs auxin")

#!/usr/bin/env Rscript

# Create a memory-conscious, balanced full-root Shahan reference while
# retaining every gene from the RNA counts assay.
#
# Run from iSensors-supplementary/:
# Rscript 00-iSensors-objects/create-balanced-shahan-full-genes.R

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

options(Seurat.object.assay.version = "v5")

input_rds <- "00-iSensors-objects/data/02_ra_integrated_Shahan_20250114.rds"
output_rds <- "00-iSensors-objects/data/shahan_balanced_full_root_full_genes.rds"
sampling_tsv <- "00-iSensors-objects/data/shahan_balanced_full_root_sampling.tsv"

cells_per_label <- 500L
seed <- 100L
set.seed(seed)

log_step <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

rna_counts <- function(obj) {
  assay <- obj[["RNA"]]
  if (is.null(assay)) stop("Input object has no RNA assay")

  layers <- Layers(assay)
  count_layers <- grep("^counts($|\\.)", layers, value = TRUE)
  if (!length(count_layers)) {
    stop("RNA assay has no counts layer. Available layers: ",
         paste(layers, collapse = ", "))
  }
  if (length(count_layers) > 1L) {
    log_step("Joining ", length(count_layers), " RNA count layers")
    assay <- JoinLayers(assay, layers = count_layers, new = "counts")
  }
  LayerData(assay, layer = "counts")
}

log_step("Loading full Shahan object: ", input_rds)
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object")
if (!"Atlas_new" %in% colnames(obj[[]])) {
  stop("Input metadata does not contain Atlas_new")
}

labels <- as.character(obj$Atlas_new)
eligible <- which(!is.na(labels) & nzchar(labels))
cells_by_label <- split(colnames(obj)[eligible], labels[eligible])

sampled_cells <- unlist(lapply(cells_by_label, function(cells) {
  if (length(cells) <= cells_per_label) cells
  else sample(cells, cells_per_label)
}), use.names = FALSE)

sampling <- data.frame(
  Atlas_new = names(cells_by_label),
  available = lengths(cells_by_label),
  sampled = pmin(lengths(cells_by_label), cells_per_label),
  stringsAsFactors = FALSE
)
write.table(
  sampling,
  sampling_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

log_step(
  "Selected ", length(sampled_cells), " cells across ",
  length(cells_by_label), " Atlas_new labels"
)

counts <- rna_counts(obj)
counts <- counts[, sampled_cells, drop = FALSE]
meta <- obj[[]][sampled_cells, , drop = FALSE]

# Rebuilding from counts avoids carrying integrated/ALT assays, graphs,
# reductions, scaled matrices, and command histories into the compact object.
small <- CreateSeuratObject(
  counts = counts,
  assay = "RNA",
  meta.data = meta,
  project = "Shahan_balanced_full_root_full_genes"
)
Idents(small) <- "Atlas_new"

small@misc$balanced_reference <- list(
  source_file = input_rds,
  annotation = "Atlas_new",
  cells_per_label = cells_per_label,
  seed = seed,
  created = as.character(Sys.time())
)

saveRDS(small, output_rds, compress = FALSE)

log_step(
  "Saved ", ncol(small), " cells x ", nrow(small),
  " genes to ", output_rds
)


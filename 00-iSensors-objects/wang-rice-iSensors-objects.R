# Generate iSensors object and AverageExpression table for Wang et al. 2025 rice atlas
# Run with working directory = iSensors-supplementary/
#
# Input:  00-iSensors-objects/data/Wang_rice_iSensor.RDS
#           Seurat object with RNA assay; rice gene names (Os-prefix)
#           Key metadata: tissue_cluster_names = "Organ--CellType"
#                         cluster_names        = cell type names
#                         tissue               = organ names
# Output:
#   00-iSensors-objects/data/Wang_rice_iSensors.rds
#     — Seurat object with RNA + iSensors_mean (OSA panel) assays
#   00-iSensors-objects/data/wang_rice_avg_exp_iSensors_mean.rds
#     — average expression matrix (sensors × tissue_cluster groups)
#       column names are "Organ--CellType"

library(iSensors)
library(Seurat)
library(tidyverse)

data_dir  <- "00-iSensors-objects/data"
input_rds <- file.path(data_dir, "Wang_rice_iSensor.RDS")

# ── Load object ───────────────────────────────────────────────────────────────
rice_obj <- readRDS(input_rds)
message("Loaded: ", ncol(rice_obj), " cells; assays: ",
        paste(Assays(rice_obj), collapse = ", "))

# ── Keep only RNA assay ───────────────────────────────────────────────────────
assays_to_remove <- setdiff(Assays(rice_obj), "RNA")
for (a in assays_to_remove) rice_obj[[a]] <- NULL
DefaultAssay(rice_obj) <- "RNA"
message("Assays after cleanup: ", paste(Assays(rice_obj), collapse = ", "))
message("RNA genes: ", nrow(rice_obj))
message("Sample RNA gene names: ", paste(head(rownames(rice_obj), 8), collapse = ", "))
gc()

# ── Load rice (OSA) Auxin panel ───────────────────────────────────────────────
AuxinPanel <- LoadSensors(
  setName      = "Auxin",
  species      = "OSA",
  hormone      = "aux",
  customPanels = FALSE,
  randomInfo   = list(n = 3, sizes = c(100, 200, 300), majortrend = TRUE)
)
message("OSA Auxin panel loaded: ", length(AuxinPanel), " sensors")

# ── Subset RNA to panel genes (memory-efficient, same as Guo pipeline) ────────
panel_genes <- unique(unlist(AuxinPanel))
genes_keep  <- intersect(panel_genes, rownames(rice_obj))
message("Panel genes found in RNA: ", length(genes_keep), " of ", length(panel_genes))

if (length(genes_keep) == 0) stop("No panel genes found — check species name or gene ID format")

obj_for_calc <- subset(rice_obj, features = genes_keep)
gc()

# ── Calculate iSensors ────────────────────────────────────────────────────────
message("Running CalcSensors on ", ncol(obj_for_calc), " cells x ",
        nrow(obj_for_calc), " genes ...")
obj_for_calc <- CalcSensors(
  obj_for_calc,
  seurLayer = "data",
  panelSet  = AuxinPanel,
  signals   = "mean"
)
message("iSensors_mean assay created: ",
        nrow(obj_for_calc[["iSensors_mean"]]), " sensors")

# ── Transfer iSensors_mean back to the full object ────────────────────────────
rice_obj[["iSensors_mean"]] <- obj_for_calc[["iSensors_mean"]]
rm(obj_for_calc); gc()

# ── Save updated Seurat object ────────────────────────────────────────────────
out_rds <- file.path(data_dir, "Wang_rice_iSensors.rds")
saveRDS(rice_obj, out_rds)
message("Saved Seurat object: ", out_rds)

# ── AverageExpression per Organ--CellType group ───────────────────────────────
# group.by = "tissue_cluster_names" groups cells by "Organ--CellType" values
# matching the column naming used throughout the Figure 7 pipeline
DefaultAssay(rice_obj) <- "iSensors_mean"
avg_mat <- AverageExpression(
  rice_obj,
  assays   = "iSensors_mean",
  slot     = "data",
  group.by = "tissue_cluster_names",
  verbose  = FALSE
)[["iSensors_mean"]]

message("AverageExpression: ", nrow(avg_mat), " sensors x ", ncol(avg_mat), " groups")
message("Sample group names: ", paste(head(colnames(avg_mat), 5), collapse = ", "))
message("Sensors: ", paste(rownames(avg_mat), collapse = ", "))

# Quick sanity check
valid_vals <- sum(!is.na(avg_mat) & avg_mat != 0)
message("Valid (non-NA, non-zero) entries: ", valid_vals,
        " of ", prod(dim(avg_mat)))

# ── Save averaged expression matrix ──────────────────────────────────────────
out_avg <- file.path(data_dir, "wang_rice_avg_exp_iSensors_mean.rds")
saveRDS(avg_mat, out_avg)
message("Saved avg_exp matrix: ", out_avg)
